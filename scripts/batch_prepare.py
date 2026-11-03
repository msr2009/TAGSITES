"""
batch_prepare.py

The batch app's Progress step: materialize every protein of a batch as a run folder and warm
the Results load cache.

Database proteins go through proteome_export (prepare_run, export_analysis, then
export_reagents with arms skipped: review never shows arms, and regenerating them for every
residue x guide costs ~1 s per protein; batch_output builds arms only for accepted sites).
Proteins missing from the database run the ordinary EBI pipeline when the batch opted in
(run_tag_sites_from_json, at most batch_app.ebi_max_concurrent at once), using
params/worm_default.json as the analysis template.

Resumable: a protein whose run folder is already complete is skipped, and the EBI runner
resumes pending EBI jobs from {rn}.status.json instead of resubmitting.

    python scripts/batch_prepare.py <batch.json> [--workers 6]

Matt Rich, 2026
"""

import argparse
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

_SCRIPTS = Path(__file__).parent
_ROOT = _SCRIPTS.parent
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_ROOT))

import batch_state as BS  # noqa: E402
import run_status  # noqa: E402

EBI_TEMPLATE = _ROOT / "params" / "worm_default.json"
_lock = threading.Lock()


def _settle(path, batch, key, **fields):
    """Thread-safe update of one protein entry, saved to the batch file."""
    with _lock:
        BS.update_protein(path, batch, key, **fields)


def run_is_complete(run_json):
    """True when a run folder's analyses (and its reagents, when present) finished."""
    if not run_json or not os.path.exists(run_json):
        return False
    with open(run_json) as fh:
        j = json.load(fh)
    gb = j["global"]
    st = run_status.load_status(gb["working_dir"], gb["run_name"])
    if gb.get("source") == "proteome_db":
        return all(st.get(t, {}).get("status") in ("success", "failed")
                   for t in ("Analyses_db", "Reagents_db"))
    return bool(st) and all(t.get("status") == "success" for t in st.values())


def export_db_protein(batch_path, protein):
    """Export one database protein (no arms); return the fields to store on its entry."""
    import proteome_db as P
    import proteome_export as X
    out_dir = str(BS.protein_dir(batch_path, protein["key"])) + "/"
    rn = BS.safe_name(protein["key"])
    t0 = time.time()
    conn = P.open_db()
    try:
        run_json = X.prepare_run(protein["pid"], out_dir, rn, conn)
        X.export_analysis(protein["pid"], out_dir, rn, conn, render_alignment=False)
        reagent_log = X.export_reagents(protein["pid"], out_dir, rn, conn, with_arms=False)
    finally:
        conn.close()
    fields = {"run_json": run_json, "seconds": round(time.time() - t0, 2)}
    # a product whose batch reagent design failed can be viewed but yields no reagents
    if protein.get("reagent_status") != "success":
        fields.update(status="failed", error=(reagent_log or ["no reagents"])[0])
    else:
        fields.update(status="ready", error="")
    return fields


def _ebi_sequence(protein):
    """(sequence, gene) for an EBI protein: its FASTA record, else its UniProt entry."""
    if protein.get("sequence"):
        return protein["sequence"], protein.get("gene") or ""
    import uniprot_api
    acc = uniprot_api.parse_accession(protein["input"])
    if not acc:
        raise ValueError("not in the database and not a UniProt accession or FASTA record")
    fasta = uniprot_api.fetch_fasta(acc)
    seq = "".join(ln.strip() for ln in fasta.splitlines() if not ln.startswith(">"))
    entry = uniprot_api.fetch_entry(acc)
    genes = entry.get("genes") or [{}]
    gene = (genes[0].get("geneName") or {}).get("value", "")
    return seq, gene


def run_ebi_protein(batch_path, protein, settings):
    """Run the EBI pipeline for one protein missing from the database."""
    import batch_run_tag_sites as BR
    from providers import _load_config as _load_batch_config
    seq, gene = _ebi_sequence(protein)
    out_dir = str(BS.protein_dir(batch_path, protein["key"])) + "/"
    rn = BS.safe_name(protein["key"])
    os.makedirs(out_dir, exist_ok=True)
    fasta = "{}{}.fa".format(out_dir, rn)
    with open(fasta, "w") as fh:
        fh.write(">{}\n{}\n".format(rn, seq))
    with open(EBI_TEMPLATE) as fh:
        template = json.load(fh)
    cfg = _load_batch_config().get("batch_run", {})
    taxid = str(settings.get("taxid") or cfg.get("taxid") or "6239")
    row = {"run_name": rn, "input_file": fasta, "gene": gene, "working_dir": out_dir}
    run_json, run = BR.build_protein_run_json(row, template, settings.get("email", ""),
                                              str(Path(batch_path).parent), taxid=taxid,
                                              flank_bp=int(cfg.get("genomic_flank_bp", 2000)))
    with open(run_json, "w") as fh:
        json.dump(run, fh, indent=4)
    t0 = time.time()
    name, ok, message = BR.run_one(run_json, force=False)
    has_reagents = "REAGENTS_reagents" in run.get("tasks", {})
    fields = {"run_json": run_json, "gene": protein.get("gene") or gene or None,
              "seconds": round(time.time() - t0, 1),
              "reagent_status": "success" if ok and has_reagents else "failed"}
    if ok and has_reagents:
        fields.update(status="ready", error="")
    else:
        fields.update(status="failed", error=message if not ok else "no genomic region found")
    return fields


def _job(batch_path, batch, protein, settings):
    """Materialize one protein; record success or the error on its entry."""
    key = protein["key"]
    _settle(batch_path, batch, key, status="running", error="")
    try:
        if protein.get("source") == "db":
            fields = export_db_protein(batch_path, protein)
        else:
            fields = run_ebi_protein(batch_path, protein, settings)
    except Exception as e:
        fields = {"status": "failed", "error": "{}: {}".format(type(e).__name__, e)}
    _settle(batch_path, batch, key, **fields)
    return key, fields


def todo(batch):
    """Proteins still to materialize: not done on disk and not skipped."""
    settings = batch.get("settings", {})
    out = []
    for p in batch["proteins"]:
        if p.get("duplicate_of"):
            continue
        if p.get("source") == "missing" and not settings.get("ebi_for_missing"):
            continue
        # a folder that finished is kept as is, whatever the entry says (resume)
        if p.get("status") in ("ready", "in_progress", "accepted", "skipped", "failed") \
                and run_is_complete(p.get("run_json")):
            continue
        out.append(p)
    return out


def prepare(batch_path, batch, workers=6, ebi_workers=3, on_done=None):
    """Materialize every protein still to do; DB exports and EBI runs use separate pools."""
    settings = batch.get("settings", {})
    jobs = todo(batch)
    for p in jobs:
        if p.get("source") == "missing":
            p["source"] = "ebi"
    db_jobs = [p for p in jobs if p.get("source") == "db"]
    ebi_jobs = [p for p in jobs if p.get("source") == "ebi"]
    with ThreadPoolExecutor(max_workers=max(1, workers)) as dbp, \
            ThreadPoolExecutor(max_workers=max(1, ebi_workers)) as ebp:
        futs = [dbp.submit(_job, batch_path, batch, p, settings) for p in db_jobs]
        futs += [ebp.submit(_job, batch_path, batch, p, settings) for p in ebi_jobs]
        # report each protein as it lands, so Results can open on the first one
        for f in as_completed(futs):
            key, fields = f.result()
            if on_done:
                on_done(key, fields)
    return len(jobs)


def render_in_subprocess(batch_path):
    """Start the deferred PDF rendering in its own low-priority process; return the Popen.

    A separate process, not a thread: rendering is CPU-bound Python, and in the app's process
    it would hold the GIL and stall the Shiny event loop (protein switches slowed ~5x).
    """
    import subprocess
    cmd = [sys.executable, str(Path(__file__).resolve()), str(batch_path), "--render-alignments"]
    return subprocess.Popen(cmd, cwd=str(_ROOT), stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL, preexec_fn=lambda: os.nice(10))


def render_pending_alignments(batch, stop=None):
    """Render the deferred alignment PDFs one at a time; return how many were rendered.

    stop() returning True ends it early. The app runs this via render_in_subprocess.
    """
    import proteome_export as X
    from task_registry import companion_path
    n = 0
    for p in batch["proteins"]:
        if stop and stop():
            break
        path = p.get("run_json")
        if p.get("source") != "db" or not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            tasks = json.load(fh).get("tasks", {})
        for t in tasks.values():
            if t.get("type") != "blast":
                continue
            aln = companion_path(t["args"]["output"], "blast", "alignment")
            if os.path.exists(aln) and not os.path.exists(aln.removesuffix(".aln") + ".pdf"):
                X.render_alignment_pdf(aln)
                n += 1
    return n


def preload(batch, keys=None):
    """Fill the Results load cache for the batch's reviewable proteins; return seconds taken."""
    from config import RESULTS_TYPE_DICT
    from utils.results import (cached_load, load_data_from_json, load_reagents_df,
                               load_run_metadata)
    t0 = time.time()
    for p in batch["proteins"]:
        if keys is not None and p["key"] not in keys:
            continue
        path = p.get("run_json")
        if not path or not os.path.exists(path):
            continue
        try:
            with open(path) as fh:
                j = json.load(fh)
            cached_load("data", j, load_data_from_json, RESULTS_TYPE_DICT)
            cached_load("meta", j, load_run_metadata)
            cached_load("reagents", j, load_reagents_df)
        except Exception as e:   # a bad folder must not stop the rest; Results reports it
            print("preload {}: {}".format(p["key"], e), file=sys.stderr)
    return time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Materialize a batch's run folders.")
    ap.add_argument("batch_json")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--render-alignments", action="store_true",
                    help="only render the deferred alignment PDFs")
    a = ap.parse_args()
    b = BS.load_batch(a.batch_json)
    if a.render_alignments:
        print("{} alignment PDF(s) rendered".format(render_pending_alignments(b)))
        sys.exit(0)
    cfg = BS.app_config()
    t = time.time()
    n = prepare(a.batch_json, b, a.workers or int(cfg.get("workers", 6)),
                int(cfg.get("ebi_max_concurrent", 3)),
                on_done=lambda k, f: print("{:<16} {:<8} {}".format(k, f.get("status"),
                                                                   f.get("error", ""))))
    print("{} protein(s) in {:.1f} s; {}".format(n, time.time() - t, json.dumps(BS.counts(b))))
