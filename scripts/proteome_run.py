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

Concurrency is a bounded pool of worker *processes* across proteins (each protein's
own tasks run sequentially within its worker; processes, not threads, because JSD
scoring is GIL-bound Python/numpy) — sized from batch.config.json's
batch_run.max_workers (default: os.cpu_count()). This avoids
run_tag_sites_from_json.py's per-stderr-line status-file rewrite (the O(L^2)
issue noted in the plan): progress here is a single line appended to a JSONL
file per finished protein, not a live per-task status file, so resuming a
killed run only means skipping accessions already marked "success" in that
file — no partial in-run state to reconcile.

Each protein's outputs are left as the same flat files utils/results.py already reads
for a single run; scripts/build_proteome_db.py consolidates them into one SQLite file
(see its docstring, and run its --estimate first).

What this does NOT do (yet): implement per-task remote fallback when a local backend can't
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

With backends.conservation = "local" and blast among the tasks, a batched DIAMOND
pass (scripts/conservation_presearch.py) runs first and each protein's task then only
does MAFFT + JSD; --no-presearch-conservation turns that off.
"""

import gzip
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))

from local_store import _reference_dir, open_index
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

    if include_reagents or "reagents" in task_types:
        # genomic_fasta only needs to be non-empty to trigger the genewise
        # pre-step in task_runners.run_reagents(); the bulk backend resolves
        # the actual genomic region from the accession, ignoring its content
        # (see genewise_bulk.py) — the remote backend, if configured instead,
        # would need a real genomic FASTA here, which this driver doesn't build
        out = f"{working_dir}/{run_name}_reagents.tsv"
        tasks.append({"id": "reagents", "analysis": "reagents", "output": out,
                      "args": {**common, "output": out, "genomic_fasta": fasta_path}})

    return tasks


def _empty_result_reason(task):
    """Reason string when a task legitimately produced no output file, else None."""
    if task["analysis"] == "blast":
        # a batch-searched sequence with zero DIAMOND hits has nothing to align
        from conservation_local import _load_cached_hits

        with open(task["args"]["input_file"]) as f:
            seq = "".join(ln.strip() for ln in f if not ln.startswith(">"))
        if _load_cached_hits(seq) == []:
            return "skipped: no homologs found"
    return None


def run_protein(accession, seq, out_dir, task_types, include_reagents=False):
    """Run every configured task for one protein; returns a result dict
    {"accession", "status", "tasks": {analysis: "ok"|"<error message>"}}.
    Never raises — a task's exception is caught and recorded so one protein's
    failure can't take down the whole batch.
    """
    t_start = time.time()
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
                task_results[task["analysis"]] = _empty_result_reason(task) or \
                    "no output file produced"
            else:
                task_results[task["analysis"]] = "ok"
        except Exception as exc:
            # a protein with no AlphaFold model is a legitimate empty result, not a failure
            if task["analysis"] == "plddt" and "No PDB path set" in str(exc):
                task_results["plddt"] = "skipped: no AFDB model"
            else:
                task_results[task["analysis"]] = f"error: {exc}"

    # "skipped: ..." marks a legitimate empty result, so it counts as done for resume
    done = [v == "ok" or v.startswith("skipped:") for v in task_results.values()]
    overall = "success" if all(done) else "partial"
    if not any(done):
        overall = "failed"
    return {"accession": accession, "status": overall, "tasks": task_results,
            "seconds": round(time.time() - t_start, 1)}


def _append_status(status_path, result):
    """Append one result line, reopening the file each time (a long-lived handle on a
    network mount can go stale and silently wedge the recorder) and retrying on errors.
    """
    line = json.dumps(result) + "\n"
    for attempt in range(6):
        try:
            with open(status_path, "a") as f:
                f.write(line)
            return
        except OSError as exc:
            print(f"[proteome_run] warning: status write failed ({exc}); retry {attempt + 1}",
                  flush=True)
            time.sleep(2 ** attempt)
    # outputs are on disk either way; a missing status line only means this protein is redone
    print(f"[proteome_run] warning: gave up recording {result['accession']}", flush=True)


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


def _should_presearch(task_types, presearch):
    """True when the batched conservation search should run: forced on/off by
    `presearch`, else automatic when blast is a task and the conservation
    backend is local.
    """
    if presearch is not None:
        return presearch and "blast" in task_types
    backend = _load_batch_config().get("backends", {}).get("conservation", "remote")
    return "blast" in task_types and backend == "local"


def _read_fasta_gz(path):
    """Read a (gzipped) FASTA into {first header token: sequence}."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name = {}, None
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = ""
            else:
                records[name] += line.strip()
    return records


def isoform_rows(table_seqs):
    """Return [(id, seq)] for sequences the proteins table lacks: UniProt isoform records
    (id = accession-N) first, then WormBase-only sequences (id = WormBase transcript name).

    local_store.py is built from the canonical-only UniProt JSON, so isoform records that
    are only in the proteome FASTA (and WormBase isoforms matching nothing in UniProt) never
    reach the run. Sequences already in the table, or repeated here, are skipped.
    """
    ref_dir, proteome_id = _reference_dir()
    ref_cfg = _load_batch_config().get("reference_data", {})
    wb_path = Path(ref_cfg["protein_fasta"])
    if not wb_path.is_absolute():
        wb_path = _REPO_ROOT / wb_path

    seen = set(table_seqs)
    rows = []
    # UniProt FASTA headers look like sp|A0A0K3AUE4-10|SEA2_CAEEL; the middle field is the id
    uniprot = _read_fasta_gz(ref_dir / f"{proteome_id}.fasta.gz")
    for header, seq in uniprot.items():
        acc = header.split("|")[1] if "|" in header else header
        if "-" in acc and seq not in seen:
            seen.add(seq)
            rows.append((acc, seq))
    uniprot_seqs = set(uniprot.values())
    for name, seq in _read_fasta_gz(wb_path).items():
        # WormBase-only: matches neither the table nor any UniProt FASTA record
        if seq not in seen and seq not in uniprot_seqs:
            seen.add(seq)
            rows.append((name, seq))
    return rows


def main(out_dir, task_types=None, limit=None, accessions=None, workers=None,
         include_reagents=False, force=False, presearch=None, isoforms=False, rows=None):
    """Run the configured tasks for every protein in local_store.py's proteins
    table (or `accessions`, if given), skipping accessions already marked
    "success" in {out_dir}/_status.jsonl unless force=True. presearch=None
    batch-searches conservation hits up front when that backend is local
    (scripts/conservation_presearch.py); True/False forces it on/off.
    isoforms=True runs isoform_rows() instead of the table, without the uniprot task.
    rows=[(id, sequence)] runs exactly those proteins (ids need not be in the table), which
    is how scripts/build_proteome_db.py drives the reagent stage; "reagents" is then an
    ordinary task type.
    """
    task_types = task_types or DEFAULT_TASKS
    if isoforms:
        # curated UniProt features are canonical-numbered, so they do not map onto isoforms
        task_types = [t for t in task_types if t != "uniprot"]
    os.makedirs(out_dir, exist_ok=True)
    status_path = os.path.join(out_dir, "_status.jsonl")

    conn = open_index()
    try:
        if rows is not None:
            rows = list(rows)
        elif isoforms:
            rows = isoform_rows(r[0] for r in conn.execute("SELECT sequence FROM proteins"))
            if limit:
                rows = rows[:int(limit)]
        elif accessions:
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

    # one batched DIAMOND pass first, so each protein's conservation task only has to
    # run MAFFT + JSD (see conservation_presearch.py)
    if _should_presearch(task_types, presearch):
        from conservation_presearch import run_presearch

        run_presearch(todo)

    # worker processes (not threads): JSD scoring is Python + numpy and would serialize on
    # the GIL; single-threaded BLAS keeps the pool within the configured thread cap
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"

    n_success = n_partial = n_failed = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(run_protein, acc, seq, out_dir, task_types, include_reagents): acc
            for acc, seq in todo
        }
        try:
            for i, future in enumerate(as_completed(futures), 1):
                acc = futures[future]
                try:
                    result = future.result()
                except Exception as exc:  # pragma: no cover — run_protein itself never raises
                    result = {"accession": acc, "status": "failed",
                              "tasks": {"_driver": f"error: {exc}"}}

                result["timestamp"] = time.time()
                _append_status(status_path, result)

                if result["status"] == "success":
                    n_success += 1
                elif result["status"] == "partial":
                    n_partial += 1
                else:
                    n_failed += 1

                if i % 100 == 0 or i == len(todo):
                    print(f"[proteome_run] {i}/{len(todo)} done "
                          f"(success={n_success} partial={n_partial} failed={n_failed})",
                          flush=True)
        except BaseException:
            # without this, leaving the with-block waits for every queued protein to finish
            # while nothing is being recorded
            pool.shutdown(wait=False, cancel_futures=True)
            raise

    print(f"[proteome_run] finished: success={n_success} partial={n_partial} "
          f"failed={n_failed} -> {status_path}")


if __name__ == "__main__":
    from argparse import ArgumentParser, BooleanOptionalAction

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
    parser.add_argument("--presearch-conservation", action=BooleanOptionalAction, default=None,
                        help="batch-search conservation hits up front (default: automatic when "
                             "backends.conservation is local and blast is a task)")
    parser.add_argument("--isoforms", action="store_true",
                        help="run UniProt isoform records and WormBase-only sequences missing "
                             "from the proteins table (no uniprot task)")
    args = parser.parse_args()

    accessions = None
    if args.accessions_file:
        with open(args.accessions_file) as f:
            accessions = [ln.strip() for ln in f if ln.strip()]

    task_types = args.tasks.split(",") if args.tasks else None

    main(args.out_dir, task_types=task_types, limit=args.limit, accessions=accessions,
         workers=args.workers, include_reagents=args.reagents, force=args.force,
         presearch=args.presearch_conservation, isoforms=args.isoforms)
