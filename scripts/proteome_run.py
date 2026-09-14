"""
proteome_run.py

Phase C of the proteome-scale batch pipeline: drives the existing analysis
pipeline across every protein in the local UniProt proteome index
(scripts/local_store.py) instead of one protein at a time from the Shiny app.

Deliberately does NOT reimplement per-task dispatch: it builds the same
{"id", "analysis", "args", "output"} task dicts modules/progress_server.py
uses and calls scripts/task_runners.py's existing TASK_RUNNERS /
afdb_presearch() directly, in-process. Those runners already resolve their
backend via scripts/providers.py, so whichever mode batch.config.json's
"backends" block selects (remote or local/bulk) is picked up automatically —
this driver has no backend-selection logic of its own.

Concurrency is a bounded worker pool across *proteins* (each protein's own
tasks run sequentially within its worker) — sized from batch.config.json's
batch_run.max_workers (default: os.cpu_count()). This avoids
run_tag_sites_from_json.py's per-stderr-line status-file rewrite (the O(L^2)
issue noted in the plan): progress here is a single line appended to a JSONL
file per finished protein, not a live per-task status file, so resuming a
killed run only means skipping accessions already marked "success" in that
file — no partial in-run state to reconcile.

What this does NOT do (yet): consolidate the ~100k loose per-protein output
files into one queryable store (Parquet/SQLite) — each protein's outputs are
left as the same flat files utils/results.py already reads for a single run.
Also does not implement per-task remote fallback when a local backend can't
resolve a protein (domains_bulk.py/structure_bulk.py degrade to an empty/
"not found" result; genewise_bulk.py and conservation_local.py raise) —
doing that safely at proteome scale needs its own rate-limited dispatch
(batch.config.json's batch_run.max_concurrent_network_calls anticipates this)
and is left for a follow-up; failures are recorded per protein/task in the
JSONL log instead of silently falling back to a queued EBI job.

Usage
-----
    python scripts/proteome_run.py --out-dir data/runs/proteome_v1 --limit 50
    python scripts/proteome_run.py --out-dir data/runs/proteome_v1 --workers 8
    python scripts/proteome_run.py --out-dir data/runs/proteome_v1 --tasks domains,plddt
"""

import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))

from local_store import open_index
from providers import _load_config as _load_batch_config
from task_runners import TASK_RUNNERS, afdb_presearch

DEFAULT_TASKS = ["domains", "plddt", "modifications", "uniprot", "scores", "blast"]

_TABLES = _REPO_ROOT / "tables"


def _batch_run_config():
    cfg = _load_batch_config().get("batch_run", {})
    return {
        "max_workers": cfg.get("max_workers") or os.cpu_count() or 4,
        "genomic_flank_bp": cfg.get("genomic_flank_bp", 2000),
    }


def build_tasks_for_protein(accession, seq, working_dir, run_name, task_types,
                            email="", taxid="6239", include_reagents=False):
    """Build the task list for one protein, in the same {"id", "analysis",
    "args", "output"} shape modules/progress_server.py / task_runners.py use.

    email is passed through to every task's args (some remote-backend fallback
    paths still expect the key to exist) but is never required to be a real
    address when every configured backend is local — none of the four local/
    bulk backends use it.
    """
    fasta_path = os.path.join(working_dir, f"{run_name}.fa")
    with open(fasta_path, "w") as f:
        f.write(f">{run_name}\n{seq}\n")

    common = {"input_file": fasta_path, "fasta": fasta_path, "email": email,
              "working_dir": working_dir, "run_name": run_name, "taxid": taxid}

    tasks = []

    if "domains" in task_types:
        out = f"{working_dir}/{run_name}_domains.txt"
        tasks.append({"id": "domains", "analysis": "domains", "output": out,
                      "args": {**common, "output": out}})

    if "plddt" in task_types:
        out = f"{working_dir}/{run_name}_plddt.txt"
        tasks.append({"id": "plddt", "analysis": "plddt", "output": out,
                      "args": {**common, "output": out, "pdb": "", "existing_AF2": 1}})

    if "modifications" in task_types:
        out = f"{working_dir}/{run_name}_mods.txt"
        tasks.append({"id": "modifications", "analysis": "modifications", "output": out,
                      "args": {**common, "output": out,
                               "sites_file": str(_TABLES / "modification_sites.txt")}})

    if "uniprot" in task_types:
        out = f"{working_dir}/{run_name}_uniprotfeat.txt"
        tasks.append({"id": "uniprot", "analysis": "uniprot", "output": out,
                      "args": {**common, "output": out, "accession": accession,
                               "features_file": str(_TABLES / "uniprot_features.txt")}})

    if "scores" in task_types:
        out = f"{working_dir}/{run_name}_scores.tsv"
        tasks.append({"id": "scores", "analysis": "scores", "output": out,
                      "args": {**common, "output": out, "window": 9,
                               "scores_file": str(_TABLES / "hydrophobicity_kyte-doolittle.tsv")}})

    if "blast" in task_types:
        out = f"{working_dir}/{run_name}_conservation.jsd"
        tasks.append({"id": "blast", "analysis": "blast", "output": out,
                      "args": {**common, "output": out, "evalue": "1e-10",
                               "max_hits": "20", "db": "uniprotkb"}})

    if include_reagents and "plddt" in task_types:
        # genomic_fasta only needs to be non-empty to trigger the genewise
        # pre-step in task_runners.run_reagents(); the bulk backend resolves
        # the actual genomic region from the accession, ignoring its content
        # (see genewise_bulk.py) — the remote backend, if configured instead,
        # would need a real genomic FASTA here, which this driver doesn't build
        out = f"{working_dir}/{run_name}_reagents.tsv"
        tasks.append({"id": "reagents", "analysis": "reagents", "output": out,
                      "args": {**common, "output": out, "genomic_fasta": fasta_path}})

    return tasks


def run_protein(accession, seq, out_dir, task_types, include_reagents=False):
    """Run every configured task for one protein; returns a result dict
    {"accession", "status", "tasks": {analysis: "ok"|"<error message>"}}.
    Never raises — a task's exception is caught and recorded so one protein's
    failure can't take down the whole batch.
    """
    working_dir = os.path.join(out_dir, accession)
    os.makedirs(working_dir, exist_ok=True)
    run_name = accession

    tasks = build_tasks_for_protein(accession, seq, working_dir, run_name, task_types,
                                    include_reagents=include_reagents)

    afdb_presearch(tasks)

    task_results = {}
    for task in tasks:
        runner = TASK_RUNNERS.get(task["analysis"])
        try:
            result = runner(task["args"], report=None, job_id_cb=None, resume_job_ids=None)
            if isinstance(result, dict) and result.get("ebi_status") in ("pending", "expired"):
                task_results[task["analysis"]] = f"ebi_status:{result['ebi_status']}"
            elif task["output"] and not os.path.exists(task["output"]):
                task_results[task["analysis"]] = "no output file produced"
            else:
                task_results[task["analysis"]] = "ok"
        except Exception as exc:
            task_results[task["analysis"]] = f"error: {exc}"

    overall = "success" if all(v == "ok" for v in task_results.values()) else "partial"
    if all(v != "ok" for v in task_results.values()):
        overall = "failed"
    return {"accession": accession, "status": overall, "tasks": task_results}


def _load_completed(status_path):
    """Read the JSONL status log and return the set of accessions already
    marked "success" — proteome_run.main()'s resume mechanism.
    """
    completed = set()
    if not os.path.exists(status_path):
        return completed
    with open(status_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue
            if entry.get("status") == "success":
                completed.add(entry.get("accession"))
    return completed


def main(out_dir, task_types=None, limit=None, accessions=None, workers=None,
         include_reagents=False, force=False):
    """Run the configured tasks for every protein in local_store.py's proteins
    table (or `accessions`, if given), skipping accessions already marked
    "success" in {out_dir}/_status.jsonl unless force=True.
    """
    task_types = task_types or DEFAULT_TASKS
    os.makedirs(out_dir, exist_ok=True)
    status_path = os.path.join(out_dir, "_status.jsonl")

    conn = open_index()
    try:
        if accessions:
            placeholders = ",".join("?" * len(accessions))
            rows = conn.execute(
                f"SELECT accession, sequence FROM proteins WHERE accession IN ({placeholders})",
                accessions,
            ).fetchall()
        else:
            query = "SELECT accession, sequence FROM proteins ORDER BY accession"
            if limit:
                query += f" LIMIT {int(limit)}"
            rows = conn.execute(query).fetchall()
    finally:
        conn.close()

    completed = set() if force else _load_completed(status_path)
    todo = [(acc, seq) for acc, seq in rows if acc not in completed]

    cfg = _batch_run_config()
    n_workers = workers or cfg["max_workers"]
    print(f"[proteome_run] {len(rows)} proteins total, {len(completed)} already "
          f"complete, {len(todo)} to run, {n_workers} workers, tasks={task_types}")

    n_success = n_partial = n_failed = 0
    with open(status_path, "a") as status_f, ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(run_protein, acc, seq, out_dir, task_types, include_reagents): acc
            for acc, seq in todo
        }
        for i, future in enumerate(as_completed(futures), 1):
            acc = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover — run_protein itself never raises
                result = {"accession": acc, "status": "failed", "tasks": {"_driver": f"error: {exc}"}}

            result["timestamp"] = time.time()
            status_f.write(json.dumps(result) + "\n")
            status_f.flush()

            if result["status"] == "success":
                n_success += 1
            elif result["status"] == "partial":
                n_partial += 1
            else:
                n_failed += 1

            if i % 100 == 0 or i == len(todo):
                print(f"[proteome_run] {i}/{len(todo)} done "
                      f"(success={n_success} partial={n_partial} failed={n_failed})")

    print(f"[proteome_run] finished: success={n_success} partial={n_partial} "
          f"failed={n_failed} -> {status_path}")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, help="output directory (one subdir per protein)")
    parser.add_argument("--tasks", type=str, default=None,
                        help=f"comma-separated task types (default: {','.join(DEFAULT_TASKS)})")
    parser.add_argument("--limit", type=int, default=None, help="only process the first N proteins")
    parser.add_argument("--accessions-file", type=str, default=None,
                        help="only process accessions listed in this file (one per line)")
    parser.add_argument("--workers", type=int, default=None,
                        help="worker pool size (default: batch.config.json's batch_run.max_workers, or os.cpu_count())")
    parser.add_argument("--reagents", action="store_true",
                        help="also run CRISPR reagent design (requires genewise)")
    parser.add_argument("--force", action="store_true",
                        help="reprocess accessions even if already marked complete")
    args = parser.parse_args()

    accessions = None
    if args.accessions_file:
        with open(args.accessions_file) as f:
            accessions = [ln.strip() for ln in f if ln.strip()]

    task_types = args.tasks.split(",") if args.tasks else None

    main(args.out_dir, task_types=task_types, limit=args.limit, accessions=accessions,
         workers=args.workers, include_reagents=args.reagents, force=args.force)
