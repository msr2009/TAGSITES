"""
proteome_export.py

Materialize one protein of the proteome database (scripts/build_proteome_db.py) as an ordinary
run folder: {rn}.run.json, the per-task result files the Results tab reads, and the reagents
TSV the Reagents tab reads. Every existing loader, download and bundle path then works on it
unchanged. The files invert the ingest in build_proteome_db.parse_protein (scaled integers
back to floats, NULL back to the -1000 gap sentinel, modification stop back to exclusive).

Three steps, each safe to run in a worker thread; the Progress tab shows the last two as the
tasks "Analyses" and "Reagents" (status ids Analyses_db / Reagents_db in {rn}.status.json):

    prepare_run     instant   {rn}.fa, {rn}.run.json, status seeded as pending
    export_analysis <0.5 s    track files, isoforms, alignment, structure
    export_reagents 3-10 s    reagents TSV with regenerated arms, genotyping, gene model

    python scripts/proteome_db.py export trxr-1 /path/to/out_dir

Matt Rich, 2026
"""

import json
import os
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(_REPO_ROOT))

import proteome_db as P  # noqa: E402
import run_status  # noqa: E402
from build_proteome_db import db_config  # noqa: E402
from task_registry import TASK_DEFS, companion_path, task_output_suffix  # noqa: E402

LABEL = "DB"
TASK_ANALYSES, TASK_REAGENTS = "Analyses_db", "Reagents_db"
# range-feature tasks: task type -> (features.track, stop shift back to the raw file)
RANGE_TASKS = {"domains": ("domain", 0), "modifications": ("modification", 1),
               "topology": ("topology", 0), "uniprot": ("uniprot", 0)}
# task type -> run_status column of the batch run, for the log
BATCH_COL = {"blast": "blast", "plddt": "plddt", "scores": "scores", "domains": "domains",
             "modifications": "modifications", "uniprot": "uniprot"}


def _write(path, text):
    """Write text to path through a temp file so a reader never sees a partial file."""
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def _tid(task):
    """The run JSON task id for a task type."""
    return "{}_{}".format(LABEL, task)


def _out(wd, rn, task):
    """The main output path of a task type for this run (reagents follow the app's naming)."""
    if task == "reagents":
        return "{}{}.reagents{}".format(wd, rn, task_output_suffix(task))
    return "{}{}.{}{}".format(wd, rn, _tid(task), task_output_suffix(task))


def _status(wd, rn, tid, **fields):
    """Update one task entry in {rn}.status.json."""
    run_status.update_task(wd, rn, tid, **fields)


def _protein(pid, conn):
    """The proteins row as a dict (with the sequence) plus the reagent status."""
    row = conn.execute("SELECT p.id, p.kind, p.accession, p.gene_name, p.seq_length, "
                       "p.sequence, p.has_structure, r.status, r.error FROM proteins p "
                       "LEFT JOIN reagent_status r ON r.pid = p.pid WHERE p.pid = ?",
                       (pid,)).fetchone()
    if row is None:
        raise KeyError("no protein with pid={}".format(pid))
    keys = ("id", "kind", "accession", "gene_name", "seq_length", "sequence", "has_structure",
            "reagent_status", "reagent_error")
    return dict(zip(keys, row))


def _tracks(pid, conn):
    """Task types that have data for this protein, in the order they are written."""
    n_cons, n_hydro, n_plddt = conn.execute(
        "SELECT COUNT(conservation), COUNT(hydrophobicity), COUNT(plddt) FROM residues "
        "WHERE pid = ?", (pid,)).fetchone()
    tracks = []
    if n_cons:
        tracks.append("blast")
    if n_plddt:
        tracks.append("plddt")
    if n_hydro:
        tracks.append("scores")
    counts = dict(conn.execute("SELECT track, COUNT(*) FROM features WHERE pid = ? "
                               "GROUP BY track", (pid,)).fetchall())
    for task, (track, _) in RANGE_TASKS.items():
        if counts.get(track):
            tracks.append(task)
    return tracks


def _meta(conn):
    """The database's build stamp for provenance."""
    return dict(conn.execute("SELECT key, value FROM meta WHERE key IN ('built', 'git_commit')"))


def prepare_run(pid, out_dir, run_name=None, conn=None):
    """Write {rn}.fa and {rn}.run.json and seed the two tasks as pending; return the JSON path."""
    from modules.setup_logic import build_global_block, build_reagents_entry, build_task_entry

    own = conn is None
    conn = conn or P.open_db()
    try:
        prot = _protein(pid, conn)
        tracks = _tracks(pid, conn)
        meta = _meta(conn)
    finally:
        if own:
            conn.close()
    rn = run_name or prot["id"]
    wd = str(out_dir).rstrip("/") + "/"
    Path(wd).mkdir(parents=True, exist_ok=True)
    fasta = "{}{}.fa".format(wd, rn)
    _write(fasta, ">{}\n{}\n".format(rn, prot["sequence"]))
    pdb = "{}{}.pdb".format(wd, rn) if prot["has_structure"] else ""
    gb = build_global_block("", rn, wd, fasta, pdb=pdb, taxid="6239")
    gb["source"] = "proteome_db"
    gb["proteome_db"] = {"path": str(P.default_path()), "pid": pid, "id": prot["id"],
                         "accession": prot["accession"], "gene_name": prot["gene_name"],
                         "kind": prot["kind"], "built": meta.get("built"),
                         "git_commit": meta.get("git_commit")}
    entries = []
    for task in tracks:
        args = {"pdb": pdb, "existing_AF2": 0} if task == "plddt" and pdb else {}
        entries.append(build_task_entry(_tid(task), task, args, task_output_suffix(task), wd, rn))
    run = {"global": gb, "tasks": {}}
    for e in entries:
        run["tasks"].update(e)
    # a reagents entry only where the batch designed some; the Reagents tab polls for its TSV
    reagents_ok = prot["reagent_status"] == "success"
    if reagents_ok:
        defaults = {k: v.get("default") for k, v in TASK_DEFS["reagents"]["params"].items()}
        run["tasks"]["REAGENTS_reagents"] = build_reagents_entry(
            defaults, "{}{}_genewise.genewise_genomic.fa".format(wd, rn),
            task_output_suffix("reagents"), wd, rn, offtarget=False, offtarget_spacer=False)
    path = "{}{}.run.json".format(wd, rn)
    _write(path, json.dumps(run, indent=4))
    for tid in (TASK_ANALYSES, TASK_REAGENTS):
        _status(wd, rn, tid, status="pending", job_ids=[], message="", output=path, log="")
    return path


def _lines(rows):
    """Tab-join rows of values into file text."""
    return "".join("\t".join(str(v) for v in r) + "\n" for r in rows)


def _dense(res, col, fmt, extra=None):
    """Rows pos-ordered 1..last with -1000 for gaps, as the loaders expect (row i = pos i+1)."""
    vals = res[col].dropna()
    if vals.empty:
        return ""
    rows = []
    for pos in range(1, int(vals.index.max()) + 1):
        v = vals.get(pos)
        cells = [pos, fmt % v if v is not None else -1000]
        if extra:
            cells.append(extra(pos))
        rows.append(cells)
    return _lines(rows)


def _sparse(res, col, fmt, seq):
    """Rows pos, value, aa for the positions that have a value."""
    vals = res[col].dropna()
    return _lines([(p, fmt % v, seq[p - 1] if 0 < p <= len(seq) else "-")
                   for p, v in vals.items()])


def _isoforms_json(pid, prot, conn):
    """The .isoforms.json dict rebuilt from the isoforms and isoform_spans tables, or None."""
    isos, spans = P.get_isoforms(pid, conn)
    if isos.empty:
        return None
    out = []
    for r in isos.itertuples():
        mine = spans[spans["iso_idx"] == r.iso_idx]
        out.append({
            "accession": r.accession, "name": r.name, "length": int(r.length),
            "is_query": bool(r.is_query),
            "present": [[int(s), int(e)] for s, e in zip(*[mine[mine.kind == "present"][c]
                                                            for c in ("start", "stop")])],
            "skipped": [[int(s), int(e)] for s, e in zip(*[mine[mine.kind == "skipped"][c]
                                                            for c in ("start", "stop")])],
            "inserts": [{"after": int(a), "length": int(n)} for a, n in
                        zip(mine[mine.kind == "insert"]["start"],
                            mine[mine.kind == "insert"]["length"])]})
    return {"query_len": int(prot["seq_length"]), "query_accession": prot["accession"],
            "source": isos.iloc[0]["source"], "isoforms": out}


def _batch_status(pid, conn):
    """{task column: 'ok'/'failed'/...} from the batch run's last status for this protein."""
    row = conn.execute("SELECT domains, plddt, modifications, uniprot, scores, blast "
                       "FROM run_status WHERE pid = ?", (pid,)).fetchone()
    return dict(zip(("domains", "plddt", "modifications", "uniprot", "scores", "blast"),
                    row)) if row else {}


def export_analysis(pid, out_dir, run_name=None, conn=None):
    """Write the track files, isoforms, alignment and structure; log each to the Analyses task."""
    own = conn is None
    conn = conn or P.open_db()
    prot = _protein(pid, conn)
    rn = run_name or prot["id"]
    wd = str(out_dir).rstrip("/") + "/"
    log = []
    _status(wd, rn, TASK_ANALYSES, status="running", stage="writing", started_at=time.time(),
            log="", message="")
    try:
        seq, batch = prot["sequence"], _batch_status(pid, conn)
        res = P.get_residues(pid, conn).set_index("pos", drop=False)
        feats = P.get_features(pid, conn=conn)
        have = set(_tracks(pid, conn))

        def note(task, text, batch_task=None):
            """Log one track with the batch's own status for it."""
            b = batch.get(BATCH_COL.get(batch_task or task, ""), None)
            log.append("{}: {}{}".format(task, text, " (batch: {})".format(b) if b else ""))
            _status(wd, rn, TASK_ANALYSES, log="\n".join(log))

        if "blast" in have:
            out = _out(wd, rn, "blast")
            _write(out, "# precomputed by {}; row i = residue i\n".format(Path(P.default_path()).name)
                   + _dense(res, "conservation", "%.5f",
                            lambda p: seq[p - 1] if p <= len(seq) else "-"))
            n_cons = int(res["conservation"].notna().sum())
            iso = _isoforms_json(pid, prot, conn)
            if iso:
                _write(companion_path(out, "blast", "isoforms"), json.dumps(iso))
            note("conservation", "{} positions, {} isoforms".format(
                n_cons, len(iso["isoforms"]) if iso else 0), "blast")
            note("alignment", _copy_alignment(pid, prot, out, conn))
        elif batch.get("blast") not in (None, "ok"):
            note("conservation", "no data in database", "blast")
        if "plddt" in have:
            out = _out(wd, rn, "plddt")
            _write(out, _dense(res, "plddt", "%.2f"))
            _write(companion_path(out, "plddt", "sasa"), _sparse(res, "rsasa", "%.4f", seq))
            _write(companion_path(out, "plddt", "hydrophobic"),
                   _sparse(res, "struct_hydro", "%.4f", seq))
            patches = feats[feats.track == "hydrophobic_patch"]
            _write(companion_path(out, "plddt", "patches"), _lines(
                [(r.source, r.start, r.stop, r.description) for r in patches.itertuples()]))
            note("plddt", "{} residues, {} patch residues".format(
                int(res["plddt"].notna().sum()), len(patches)))
            note("structure", _copy_structure(prot, wd, rn, run_dir=db_config()["run_dir"]))
        if "scores" in have:
            _write(_out(wd, rn, "scores"), _dense(res, "hydrophobicity", "%.6f",
                                                  lambda p: seq[p - 1] if p <= len(seq) else "-"))
            note("scores", "{} positions".format(int(res["hydrophobicity"].notna().sum())))
        for task, (track, shift) in RANGE_TASKS.items():
            sub = feats[feats.track == track]
            if task in have:
                _write(_out(wd, rn, task), _lines(
                    [(r.source, r.start, r.stop + shift, r.description) for r in sub.itertuples()]))
                note(task, "{} features".format(len(sub)))
            elif task != "topology":
                # the batch ingests a task's file only when its last status was ok
                note(task, "no data in database")
        _status(wd, rn, TASK_ANALYSES, status="success", stage="done", message="",
                log="\n".join(log))
    except Exception as e:
        log.append("failed: {}: {}".format(type(e).__name__, e))
        _status(wd, rn, TASK_ANALYSES, status="failed", message=str(e), log="\n".join(log))
        raise
    finally:
        if own:
            conn.close()
    return log


def _copy_alignment(pid, prot, out, conn):
    """Copy the batch's .aln beside the jsd and render its PDF; return a log phrase."""
    row = conn.execute("SELECT aln_path FROM conservation_summary WHERE pid = ?",
                       (pid,)).fetchone()
    src = db_config()["run_dir"] / row[0] if row and row[0] else None
    if not src or not src.exists():
        return "not in the batch folder"
    dest = companion_path(out, "blast", "alignment")
    shutil.copyfile(src, dest)
    try:
        from build_heatmap_reportlab import plot_alignment_reportlab
        plot_alignment_reportlab(dest)
    except Exception as e:   # the alignment image is optional; the Results tab says so
        return "copied; PDF failed ({}: {})".format(type(e).__name__, e)
    return "copied and rendered"


def _copy_structure(prot, wd, rn, run_dir):
    """Copy the batch's AlphaFold PDB to {rn}.pdb; return a log phrase."""
    src = Path(run_dir) / prot["id"] / (prot["id"] + ".AF.pdb")
    if not src.exists():
        return "model not in the batch folder"
    shutil.copyfile(src, "{}{}.pdb".format(wd, rn))
    return "copied"


# columns of {run}.reagents.tsv, in the order design_tag_reagents writes them
REAGENT_COLUMNS = [
    "residue_index", "amino_acid", "insert_pos", "exon_index", "is_split_codon",
    "dist_to_5p_splice", "dist_to_3p_splice", "guide_strand", "spacer", "pam_seq",
    "pam_fwd_start", "cut_pos", "distance", "pam_in_arm", "recut_block_method",
    "mutation_desc", "left_arm", "right_arm", "left_arm_wt", "right_arm_wt", "rs3_score",
    "rs3_percentile", "offtarget_count", "offtarget_identical", "offtarget_detail",
    "offtarget_status"]


def _reagent_frame(df, arms):
    """The reagent rows with arms, shaped like the TSV the design step writes."""
    df = pd.concat([df, arms], axis=1)
    df["is_split_codon"] = df["is_split_codon"].astype(bool)
    for col in ("offtarget_count", "offtarget_identical"):
        df[col] = df[col].astype("Int64")
    for col in ("mutation_desc", "offtarget_detail", "rs3_score", "rs3_percentile"):
        df[col] = df[col].astype(object).where(df[col].notna(), "")
    df["offtarget_status"] = df["offtarget_status"].fillna("")
    return df[REAGENT_COLUMNS]


def export_reagents(pid, out_dir, run_name=None, conn=None):
    """Write the reagents TSV (arms regenerated), genotyping TSV and gene-model files."""
    from cds_check import compare_to_query, translate_model, write_sidecar
    from parse_genewise import write_genewise_gff
    from site_selection_util import save_fasta

    own = conn is None
    conn = conn or P.open_db()
    prot = _protein(pid, conn)
    rn = run_name or prot["id"]
    wd = str(out_dir).rstrip("/") + "/"
    log = []
    _status(wd, rn, TASK_REAGENTS, status="running", stage="reagents", started_at=time.time(),
            log="", message="")
    try:
        if prot["reagent_status"] != "success":
            msg = "batch design {}: {}".format(prot["reagent_status"] or "not run",
                                               prot["reagent_error"] or "no reason recorded")
            _status(wd, rn, TASK_REAGENTS, status="failed", message=msg, log=msg)
            return [msg]
        t0 = time.time()
        region = P.load_region(pid, conn)
        df = P.get_reagents(pid, conn=conn)
        arms = P.regenerate_all_arms(pid, df, region, conn)
        log.append("{} residue x guide rows, {} guides, arms regenerated in {:.1f} s".format(
            len(df), df["gid"].nunique(), time.time() - t0))
        _status(wd, rn, TASK_REAGENTS, log="\n".join(log))
        # gene model first (the GenBank export reads it), reagents TSV last: the Reagents tab
        # loads the TSV the moment it appears, so everything it needs must already exist
        prefix = "{}{}_genewise".format(wd, rn)
        save_fasta("{}_genomic".format(region["transcript_id"]), region["dna"],
                   prefix + ".genewise_genomic.fa")
        write_genewise_gff(prefix + ".genewise.out.txt", region["cds_df"], region["chrom"],
                           "proteome_db (GFF3-derived)")
        check = compare_to_query(translate_model(region["cds_df"], region["dna"]),
                                 prot["sequence"])
        write_sidecar("{}{}.reagents.model_check.json".format(wd, rn), check)
        log.append("gene model {}: {}".format(check["status"], check.get("message") or "exact"))
        geno = P.get_genotyping(pid, conn).drop(columns=["pid"])
        out = _out(wd, rn, "reagents")
        _write(companion_path(out, "reagents", "genotyping"), geno.to_csv(sep="\t", index=False))
        log.append("{} genotyping primer rows".format(len(geno)))
        _write(out, _reagent_frame(df, arms).to_csv(sep="\t", index=False))
        _status(wd, rn, TASK_REAGENTS, status="success", stage="done", message="",
                log="\n".join(log))
    except Exception as e:
        log.append("failed: {}: {}".format(type(e).__name__, e))
        _status(wd, rn, TASK_REAGENTS, status="failed", message=str(e), log="\n".join(log))
        raise
    finally:
        if own:
            conn.close()
    return log


def export_run(pid, out_dir, run_name=None, conn=None):
    """All three steps in order; return the run JSON path (reagent failures do not raise)."""
    path = prepare_run(pid, out_dir, run_name, conn)
    export_analysis(pid, out_dir, run_name, conn)
    try:
        export_reagents(pid, out_dir, run_name, conn)
    except Exception as e:   # recorded in the status file; the analysis files are still good
        print("reagents failed: {}: {}".format(type(e).__name__, e), file=sys.stderr)
    return path
