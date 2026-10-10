"""
build_proteome_db.py

Consolidates the proteome batch run (scripts/proteome_run.py's ~28,600 per-protein
folders, plus the Pfam / DeepTMHMM caches in data/reference/) into ONE SQLite file for
simple lookup: scripts/proteome_db.py reads it. Tables:

  ingest   proteins, wb_transcripts, features, residues, conservation_summary,
           isoforms, isoform_spans, run_status, meta
  scores   score_configs, site_scores, tag_sites  (the app's own scorer, so values match
           the app exactly; suggested sites via utils.scoring.pick_suggested_sites)
  reagents regions, insertion_sites, guides, site_guides, edits, genotyping_primers,
           reagent_status  (normalised so the four ~1 kb arm columns are never stored;
           arms are regenerated from the genome, see proteome_db.regenerate_arms)

ALWAYS RUN `--estimate` FIRST. It builds a small random sample, measures every table with
SQLite's dbstat, projects the full size and run time, and writes nothing large. A full
build refuses to start without a fresh estimate and `--yes`, and refuses outright when the
projected size does not fit the free space on the build directory.

The database is assembled on local disk (build_dir) and then copied to its final path:
SQLite locking and journalling are unreliable on network volumes.

Conventions, all recorded in the `meta` table:
  * proteins.sequence is the sequence the batch actually analysed ({id}.fa in its folder),
    not local_store's: 28 proteins differ (another UniProt release) and are flagged by
    sequence_differs, since per-residue values index into the analysed sequence
  * a task's file is ingested only if that task's LAST status is "ok" (folders contain stale
    files from earlier failed attempts); an empty ok file means "ran, found nothing"
  * numbers are stored as scaled integers (SCALE) to keep ~16M residue rows small
  * modification ranges are stored with an INCLUSIVE stop: regex_sites.py writes
    match.end()+1, one past the last residue. The app's scorer still sees the raw file, so
    site_scores match the app, which treats those ranges as one residue too long
  * residues.hydrophobicity is the Kyte-Doolittle track exactly as the app positions it:
    row i of {id}_scores.tsv -> pos i+1 (the window START, not its centre)
  * residues.conservation: row i of the .jsd -> pos i+1; the -1000 gap sentinel is NULL

Matt Rich, 2026
"""

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import os
import random
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(_REPO_ROOT))

# value x scale, rounded to int
SCALE = {"conservation": 100000, "hydrophobicity": 10000, "plddt": 100, "rsasa": 10000,
         "struct_hydro": 10000, "score": 100}

# per-protein files written by proteome_run.build_tasks_for_protein
SUFFIX = {"domains": "_domains.txt", "modifications": "_mods.txt",
          "uniprot": "_uniprotfeat.txt", "scores": "_scores.tsv", "plddt": "_plddt.txt",
          "blast": "_conservation.jsd"}
TASKS = list(SUFFIX)

# table -> (what its size scales with in the estimate)
SCALES_WITH = {
    "proteins": "proteins", "conservation_summary": "proteins", "isoforms": "proteins",
    "isoform_spans": "proteins", "run_status": "proteins", "tag_sites": "proteins",
    "curated_sites": "fixed",
    "regions": "proteins", "guides": "proteins", "edits": "proteins",
    "reagent_status": "proteins",
    "features": "residues", "residues": "residues", "site_scores": "residues",
    "insertion_sites": "residues", "site_guides": "residues",
    "genotyping_primers": "residues",
    "wb_transcripts": "fixed", "meta": "fixed", "score_configs": "fixed",
}
STAGE_TABLES = {
    "ingest": ["proteins", "wb_transcripts", "features", "residues", "conservation_summary",
               "isoforms", "isoform_spans", "run_status", "meta"],
    "scores": ["score_configs", "site_scores", "tag_sites", "curated_sites"],
    "reagents": ["regions", "insertion_sites", "guides", "site_guides", "edits",
                 "genotyping_primers", "reagent_status"],
}

SCHEMA = {
    "proteins": """CREATE TABLE proteins (
        pid INTEGER PRIMARY KEY, id TEXT UNIQUE NOT NULL, kind TEXT, accession TEXT,
        gene_name TEXT, wormbase_gene TEXT, wormbase_transcript TEXT, wb_gene TEXT,
        wb_transcript TEXT, crc64 TEXT, seq_length INTEGER, sequence TEXT,
        has_structure INTEGER, structure_source TEXT, sequence_differs INTEGER DEFAULT 0)""",
    "wb_transcripts": """CREATE TABLE wb_transcripts (
        transcript TEXT PRIMARY KEY, wb_gene TEXT, locus TEXT, uniprot_tag TEXT,
        crc64 TEXT, pid INTEGER)""",
    "features": """CREATE TABLE features (
        pid INTEGER, track TEXT, source TEXT, start INTEGER, stop INTEGER, description TEXT)""",
    "residues": """CREATE TABLE residues (
        pid INTEGER, pos INTEGER, conservation INTEGER, hydrophobicity INTEGER,
        plddt INTEGER, rsasa INTEGER, struct_hydro INTEGER, PRIMARY KEY (pid, pos))
        WITHOUT ROWID""",
    "conservation_summary": """CREATE TABLE conservation_summary (
        pid INTEGER PRIMARY KEY, n_hits INTEGER, n_species INTEGER, best_evalue REAL,
        n_swissprot INTEGER, n_trembl INTEGER, mean_jsd REAL, aln_path TEXT, hits_path TEXT)""",
    "isoforms": """CREATE TABLE isoforms (
        pid INTEGER, iso_idx INTEGER, accession TEXT, name TEXT, length INTEGER,
        is_query INTEGER, source TEXT, PRIMARY KEY (pid, iso_idx)) WITHOUT ROWID""",
    "isoform_spans": """CREATE TABLE isoform_spans (
        pid INTEGER, iso_idx INTEGER, kind TEXT, start INTEGER, stop INTEGER, length INTEGER)""",
    "run_status": """CREATE TABLE run_status (
        pid INTEGER PRIMARY KEY, status TEXT, domains TEXT, plddt TEXT, modifications TEXT,
        uniprot TEXT, scores TEXT, blast TEXT, timestamp REAL)""",
    "meta": "CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)",
    "score_configs": """CREATE TABLE score_configs (
        config_id INTEGER PRIMARY KEY, hash TEXT UNIQUE, json TEXT, criteria TEXT,
        max_score REAL)""",
    "site_scores": """CREATE TABLE site_scores (
        config_id INTEGER, pid INTEGER, pos INTEGER, score INTEGER, masked INTEGER,
        bits INTEGER, missing INTEGER, PRIMARY KEY (config_id, pid, pos)) WITHOUT ROWID""",
    "tag_sites": """CREATE TABLE tag_sites (
        pid INTEGER, pos INTEGER, source TEXT, rank INTEGER, config_id INTEGER)""",
    "curated_sites": """CREATE TABLE curated_sites (
        pid INTEGER, wb_gene TEXT, gene TEXT, allele TEXT, pos INTEGER, verified INTEGER,
        note TEXT, source TEXT)""",
    "regions": """CREATE TABLE regions (
        pid INTEGER PRIMARY KEY, transcript_id TEXT, chrom TEXT, start INTEGER, stop INTEGER,
        strand TEXT, flank_5p INTEGER, flank_3p INTEGER)""",
    "insertion_sites": """CREATE TABLE insertion_sites (
        pid INTEGER, residue_index INTEGER, amino_acid TEXT, insert_pos INTEGER,
        exon_index INTEGER, is_split_codon INTEGER, dist_to_5p_splice INTEGER,
        dist_to_3p_splice INTEGER, PRIMARY KEY (pid, residue_index)) WITHOUT ROWID""",
    "guides": """CREATE TABLE guides (
        gid INTEGER PRIMARY KEY, pid INTEGER, strand TEXT, spacer TEXT, pam_seq TEXT,
        pam_fwd_start INTEGER, cut_pos INTEGER, rs3_score REAL, rs3_percentile REAL,
        offtarget_count INTEGER, offtarget_identical INTEGER, offtarget_status TEXT,
        offtarget_detail TEXT)""",
    "site_guides": """CREATE TABLE site_guides (
        pid INTEGER, residue_index INTEGER, gid INTEGER, distance INTEGER, pam_in_arm TEXT,
        recut_block_method TEXT, edit_id INTEGER,
        PRIMARY KEY (pid, residue_index, gid)) WITHOUT ROWID""",
    "edits": "CREATE TABLE edits (edit_id INTEGER PRIMARY KEY, text TEXT UNIQUE)",
    "genotyping_primers": """CREATE TABLE genotyping_primers (
        pid INTEGER, residue_index INTEGER, amplicon_type TEXT, fwd_seq TEXT, fwd_tm REAL,
        rev_seq TEXT, rev_tm REAL, product_size INTEGER, offtarget_amplicons INTEGER,
        offtarget_detail TEXT, PRIMARY KEY (pid, residue_index, amplicon_type))
        WITHOUT ROWID""",
    "reagent_status": """CREATE TABLE reagent_status (
        pid INTEGER PRIMARY KEY, status TEXT, seconds REAL, error TEXT)""",
}
INDEXES = [
    "CREATE INDEX idx_proteins_gene ON proteins(gene_name)",
    "CREATE INDEX idx_proteins_wbgene ON proteins(wormbase_gene)",
    "CREATE INDEX idx_proteins_wbtx ON proteins(wormbase_transcript)",
    "CREATE INDEX idx_proteins_wbname ON proteins(wb_transcript)",
    "CREATE INDEX idx_proteins_acc ON proteins(accession)",
    "CREATE INDEX idx_wbtx_gene ON wb_transcripts(wb_gene)",
    "CREATE INDEX idx_wbtx_locus ON wb_transcripts(locus)",
    "CREATE INDEX idx_wbtx_pid ON wb_transcripts(pid)",
    "CREATE INDEX idx_features_pid ON features(pid, track)",
    "CREATE INDEX idx_features_desc ON features(track, description)",
    "CREATE INDEX idx_isospans_pid ON isoform_spans(pid, iso_idx)",
    "CREATE INDEX idx_tag_sites_pid ON tag_sites(pid)",
    "CREATE INDEX idx_guides_pid ON guides(pid)",
]
# per-protein child tables, for delete-then-insert on --update
CHILD_TABLES = {
    # run_status is not listed: it is rewritten by INSERT OR REPLACE before parsing starts,
    # so deleting it here would remove the row that was just written
    "ingest": ["features", "residues", "conservation_summary", "isoforms", "isoform_spans"],
    "scores": ["site_scores", "tag_sites"],
    "reagents": ["regions", "insertion_sites", "guides", "site_guides",
                 "genotyping_primers", "reagent_status"],
}


# ── configuration ─────────────────────────────────────────────────────────────

def _abs(path):
    """Resolve a config path against the repo root."""
    p = Path(path)
    return p if p.is_absolute() else _REPO_ROOT / p


def db_config():
    """The proteome_db block of batch.config.json with defaults and absolute paths."""
    from providers import _load_config

    full = _load_config()
    cfg = full.get("proteome_db", {})
    run_dir = _abs(cfg.get("run_dir") or "data/runs/proteome_v1")
    build_dir = Path(cfg["build_dir"]) if cfg.get("build_dir") else Path(tempfile.gettempdir())
    ref = _abs(full.get("reference_data", {}).get("out_dir", "data/reference"))
    return {
        "path": _abs(cfg.get("path") or run_dir / "proteome.sqlite3"),
        "run_dir": run_dir,
        "build_dir": build_dir,
        "scratch_dir": Path(cfg["scratch_dir"]) if cfg.get("scratch_dir") else build_dir / "proteome_scratch",
        "chunk_size": int(cfg.get("chunk_size", 200)),
        "reagent_chunk": int(cfg.get("reagent_chunk", 500)),
        "max_workers": int(full.get("batch_run", {}).get("max_workers") or os.cpu_count() or 4),
        "ref_dir": ref,
        "wb_fasta": _abs(full.get("reference_data", {}).get("protein_fasta", "")),
        "curated_sites": [str(_abs(x)) for x in
                          cfg.get("curated_sites", ["data/internal_batch_sites.json"])],
    }


# ── proteins, status, cache names ─────────────────────────────────────────────

def read_wb_fasta(path):
    """{transcript: (header attributes dict, sequence)} from the WormBase protein FASTA."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    out, name, attrs, chunks = {}, None, {}, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out[name] = (attrs, "".join(chunks))
                head = line[1:].split(None, 1)
                name = head[0]
                attrs = dict(re.findall(r"(\w+)=(\S+)", head[1])) if len(head) > 1 else {}
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        out[name] = (attrs, "".join(chunks))
    return out


def load_proteins(cfg):
    """All proteome ids as dicts, plus the wb_transcripts rows.

    Canonical entries come from local_store; UniProt isoform ids and WormBase-only ids from
    proteome_run.isoform_rows() (the same list the batch --isoforms run used).
    """
    from Bio.SeqUtils.CheckSum import crc64

    import local_store
    import proteome_run

    conn = local_store.open_index()
    canon = conn.execute("SELECT accession, sequence, gene_name, wormbase_gene, "
                         "wormbase_transcript FROM proteins ORDER BY accession").fetchall()
    conn.close()
    by_acc = {r[0]: r for r in canon}
    wb = read_wb_fasta(cfg["wb_fasta"])
    wb_by_crc = {}
    for name, (_, seq) in wb.items():
        wb_by_crc.setdefault(crc64(seq), []).append(name)

    proteins = []
    for acc, seq, gene, wbg, wbt in canon:
        proteins.append({"id": acc, "kind": "uniprot", "accession": acc, "gene_name": gene,
                         "wormbase_gene": wbg, "wormbase_transcript": wbt, "sequence": seq})
    for pid_, seq in proteome_run.isoform_rows([p["sequence"] for p in proteins]):
        base = pid_.rsplit("-", 1)[0] if pid_.rsplit("-", 1)[-1].isdigit() else None
        if base and base in by_acc:
            r = by_acc[base]
            proteins.append({"id": pid_, "kind": "uniprot_isoform", "accession": base,
                             "gene_name": r[2], "wormbase_gene": r[3],
                             "wormbase_transcript": r[4], "sequence": seq})
        else:
            attrs = wb.get(pid_, ({}, ""))[0]
            proteins.append({"id": pid_, "kind": "wormbase_only", "accession": None,
                             "gene_name": attrs.get("locus"), "wormbase_gene": None,
                             "wormbase_transcript": pid_, "sequence": seq})
    # Every per-residue value indexes into the sequence the batch actually analysed, which is
    # the {id}.fa in its folder. It can differ from this machine's local_store when the two
    # were built from different UniProt releases (28 of 28,626 proteins did), so the folder's
    # sequence wins and the difference is flagged.
    from concurrent.futures import ThreadPoolExecutor

    def analysed(p):
        try:
            with open(Path(cfg["run_dir"]) / p["id"] / (p["id"] + ".fa")) as f:
                return "".join(ln.strip() for ln in f if not ln.startswith(">"))
        except OSError:
            return None

    with ThreadPoolExecutor(32) as ex:
        for p, seq in zip(proteins, ex.map(analysed, proteins)):
            p["sequence_differs"] = int(bool(seq) and seq != p["sequence"])
            if seq:
                p["sequence"] = seq
    for p in proteins:
        p["crc64"] = crc64(p["sequence"])
        names = sorted(wb_by_crc.get(p["crc64"], []))
        p["wb_names"] = names
        p["wb_transcript"] = p["id"] if p["id"] in wb else (names[0] if names else None)
        p["wb_gene"] = wb[p["wb_transcript"]][0].get("gene") if p["wb_transcript"] in wb else None
    return proteins, wb


def read_last_status(run_dir):
    """{id: status record}: top-level fields from the last line per accession, with each
    task's status taken from the last line that ran it (a task-subset re-run adds to the
    record instead of replacing it)."""
    last = {}
    with open(Path(run_dir) / "_status.jsonl") as f:
        for line in f:
            if line.strip():
                rec = json.loads(line)
                prev = last.get(rec["accession"])
                # carry over tasks the newer line did not run
                if prev:
                    rec["tasks"] = {**prev.get("tasks", {}), **rec.get("tasks", {})}
                last[rec["accession"]] = rec
    return last


# ── per-protein parsing (runs in worker processes) ────────────────────────────

def _rows(path, maxsplit=-1):
    """Rows of a headerless TSV as string lists, skipping blanks and # comments."""
    out = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                out.append(line.rstrip("\n").split("\t", maxsplit))
    return out


def _scaled(text, scale):
    """Round float(text) x scale to an int; None for blank, NaN or the -1000 gap sentinel."""
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or v == -1000:
        return None
    return int(round(v * scale))


def _range_features(path, track, stop_shift=0):
    """(track, source, start, stop, description) rows from a source/start/stop/description file."""
    out = []
    for r in _rows(path, maxsplit=3):
        if len(r) >= 4:
            out.append((track, r[0], int(r[1]), int(r[2]) + stop_shift, r[3]))
    return out


def _conservation_summary(hits_path, mean_jsd, aln_rel, hits_rel):
    """Counts and best e-value from a .json.json DIAMOND hit list."""
    with open(hits_path) as f:
        data = json.load(f)
    hits = data.get("hits", []) if isinstance(data, dict) else data
    best, species = None, set()
    n_sp = n_tr = 0
    for h in hits:
        species.add(h.get("hit_os"))
        n_sp += h.get("source_db") == "swissprot"
        n_tr += h.get("source_db") == "rhabditida_trembl"
        for hsp in h.get("hit_hsps", []):
            try:
                e = float(hsp.get("hsp_expect"))
            except (TypeError, ValueError):
                continue
            best = e if best is None else min(best, e)
    return (len(hits), len(species), best, n_sp, n_tr, mean_jsd, aln_rel, hits_rel)


def _isoform_rows(path):
    """(isoforms rows, isoform_spans rows) from a .isoforms.json (pid left for the writer)."""
    with open(path) as f:
        data = json.load(f)
    isos, spans = [], []
    for i, iso in enumerate(data.get("isoforms", [])):
        isos.append((i, iso.get("accession"), iso.get("name"), iso.get("length"),
                     int(bool(iso.get("is_query"))), data.get("source")))
        for kind in ("present", "skipped"):
            for s, e in iso.get(kind, []):
                spans.append((i, kind, s, e, e - s + 1))
        for ins in iso.get("inserts", []):
            spans.append((i, "insert", ins["after"], ins["after"], ins["length"]))
    return isos, spans


def _reagents_frame(job):
    """The per-residue reagent columns the scorer reads, or None when none were designed."""
    import pandas as pd

    rows = job.get("reagents")
    if not rows:
        return None   # no design: the guide/splice criteria stay "missing", as in the app
    return pd.DataFrame(rows, columns=["residue_index", "distance", "dist_to_5p_splice",
                                       "dist_to_3p_splice"])


def _score_protein(job, files):
    """Tag-site scores and suggested sites for one protein via the app's own scorer."""
    from config import RESULTS_TYPE_DICT
    from utils.results import load_data_from_json
    from utils.scoring import pick_suggested_sites, score_tag_sites

    tasks = {}
    for task in ("blast", "plddt", "scores", "domains", "modifications", "uniprot"):
        if files.get(task):
            tasks[task] = {"type": task, "args": {"output": files[task]}}
    if files.get("topology"):
        tasks["topology"] = {"type": "topology", "args": {"output": files["topology"]}}
    # the loader prints a line for every missing file; the batch does not need them
    with contextlib.redirect_stdout(io.StringIO()):
        aa_df, range_df, _, _ = load_data_from_json({"tasks": tasks}, RESULTS_TYPE_DICT)
        scores = score_tag_sites(aa_df, range_df, job["sequence"], _reagents_frame(job),
                                 job["scores_config"])
    keys = [c["key"] for c in job["scores_config"]["criteria"]]
    rows = []
    for pos, rec in scores.iterrows():
        bits = missing = 0
        for i, k in enumerate(keys):
            v = rec[k]
            if v is None or (isinstance(v, float) and math.isnan(v)):
                missing |= 1 << i
            elif bool(v):
                bits |= 1 << i
        rows.append((int(pos), int(round(float(rec["score"]) * SCALE["score"])),
                     int(bool(rec["masked"])), bits, missing))
    suggested = [(p, r + 1) for r, p in enumerate(pick_suggested_sites(scores))]
    return rows, suggested


def parse_protein(job):
    """Read one protein's files; return every row for the requested stages.

    A task's file is read only if its last status is "ok". Never raises: a problem is
    returned in "errors" so one bad folder cannot stop the build.
    """
    t_start = time.time()
    out = {"pid": job["pid"], "errors": [], "features": [], "residues": [], "summary": None,
           "isoforms": [], "spans": [], "site_scores": [], "suggested": [], "seconds": 0.0}
    folder, pid_ = Path(job["folder"]), job["id"]
    ok = {t: job["tasks"].get(t) == "ok" for t in TASKS}
    path = {t: str(folder / (pid_ + SUFFIX[t])) for t in TASKS}
    values = {}   # pos -> [conservation, hydrophobicity, plddt, rsasa, struct_hydro]

    def put(pos, col, val):
        if val is not None:
            values.setdefault(pos, [None] * 5)[col] = val

    def attempt(label, fn):
        try:
            fn()
        except Exception as e:   # one unreadable file must not lose the whole protein
            out["errors"].append("{}: {}: {}".format(label, type(e).__name__, e))

    if "ingest" in job["stages"]:
        for task, track, shift in (("domains", "domain", 0),
                                   ("modifications", "modification", -1),
                                   ("uniprot", "uniprot", 0)):
            if ok[task]:
                attempt(task, lambda t=task, k=track, s=shift: out["features"].extend(
                    _range_features(path[t], k, s)))
        if job.get("topology_file"):
            attempt("topology", lambda: out["features"].extend(
                _range_features(job["topology_file"], "topology")))
        if ok["plddt"]:
            def structure():
                for r in _rows(path["plddt"]):
                    put(int(r[0]), 2, _scaled(r[1], SCALE["plddt"]))
                base = path["plddt"][:-len(".txt")]
                for suffix, col, scale in ((".sasa.txt", 3, SCALE["rsasa"]),
                                           (".hydro.txt", 4, SCALE["struct_hydro"])):
                    if os.path.exists(base + suffix):
                        for r in _rows(base + suffix):
                            put(int(r[0]), col, _scaled(r[1], scale))
                if os.path.exists(base + ".patches.txt"):
                    out["features"].extend(_range_features(base + ".patches.txt",
                                                           "hydrophobic_patch"))
            attempt("plddt", structure)
        if ok["scores"]:
            def kd():
                for i, r in enumerate(_rows(path["scores"])):
                    put(i + 1, 1, _scaled(r[1], SCALE["hydrophobicity"]))
            attempt("scores", kd)
        mean_jsd = None
        if ok["blast"]:
            def jsd():
                nonlocal mean_jsd
                vals = []
                for i, r in enumerate(_rows(path["blast"])):
                    v = _scaled(r[1], SCALE["conservation"])
                    put(i + 1, 0, v)
                    if v is not None:
                        vals.append(v)
                if vals:
                    mean_jsd = sum(vals) / len(vals) / SCALE["conservation"]
            attempt("blast", jsd)
        blast_status = job["tasks"].get("blast", "")
        if ok["blast"] or blast_status.startswith("skipped"):
            hits_file = path["blast"].replace(".jsd", ".json.json")
            rel = lambda p: os.path.relpath(p, job["run_dir"])
            if os.path.exists(hits_file):
                attempt("hits", lambda: out.update(summary=_conservation_summary(
                    hits_file, mean_jsd, rel(path["blast"].replace(".jsd", ".aln")),
                    rel(hits_file))))
            else:
                out["summary"] = (0, 0, None, 0, 0, mean_jsd, None, None)
            iso_file = path["blast"].replace(".jsd", ".isoforms.json")
            if ok["blast"] and os.path.exists(iso_file) and os.path.getsize(iso_file) > 0:
                def isos():
                    out["isoforms"], out["spans"] = _isoform_rows(iso_file)
                attempt("isoforms", isos)
        out["residues"] = [(pos, *vals) for pos, vals in sorted(values.items())]

    if "scores" in job["stages"]:
        files = {t: path[t] for t in TASKS if ok[t]}
        files["topology"] = job.get("topology_file")
        try:
            out["site_scores"], out["suggested"] = _score_protein(job, files)
        except Exception as e:
            out["errors"].append("scoring: {}: {}".format(type(e).__name__, e))
    out["seconds"] = time.time() - t_start
    return out


def parse_chunk(jobs):
    """Worker entry point: parse a list of jobs."""
    return [parse_protein(j) for j in jobs]


# ── writing ───────────────────────────────────────────────────────────────────

def create_schema(conn, tables=None):
    """Create the listed tables (default all) if they do not exist."""
    for name, ddl in SCHEMA.items():
        if tables is None or name in tables:
            conn.execute(ddl.replace("CREATE TABLE", "CREATE TABLE IF NOT EXISTS", 1))


def upgrade_schema(conn):
    """Add columns introduced after a database was first built (so --update can reuse it)."""
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if cols and "sequence_differs" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN sequence_differs INTEGER DEFAULT 0")


def create_indexes(conn):
    """Create the lookup indexes after bulk inserts (much faster than maintaining them)."""
    for ddl in INDEXES:
        conn.execute(ddl.replace("CREATE INDEX", "CREATE INDEX IF NOT EXISTS", 1))


def table_in(conn, name):
    """True when the table exists in this database."""
    return conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                        (name,)).fetchone() is not None


def delete_children(conn, pid, stage):
    """Remove one protein's rows from the stage's per-protein tables (for --update)."""
    for t in CHILD_TABLES[stage]:
        if table_in(conn, t):
            conn.execute("DELETE FROM {} WHERE pid=?".format(t), (pid,))


def insert_proteins(conn, proteins, wb, status):
    """Insert/refresh the proteins and wb_transcripts tables; assigns pid to each dict."""
    existing = dict(conn.execute("SELECT id, pid FROM proteins"))
    next_pid = max(existing.values(), default=0) + 1
    rows = []
    for p in proteins:
        if p["id"] in existing:
            p["pid"] = existing[p["id"]]
        else:
            p["pid"], next_pid = next_pid, next_pid + 1
        st = status.get(p["id"], {}).get("tasks", {})
        has = int(st.get("plddt") == "ok")
        rows.append((p["pid"], p["id"], p["kind"], p["accession"], p["gene_name"],
                     p["wormbase_gene"], p["wormbase_transcript"], p["wb_gene"],
                     p["wb_transcript"], p["crc64"], len(p["sequence"]), p["sequence"], has,
                     "afdb" if has else None, p.get("sequence_differs", 0)))
    conn.executemany("INSERT OR REPLACE INTO proteins VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    pid_by_crc = {}
    for p in proteins:
        pid_by_crc.setdefault(p["crc64"], p["pid"])
    from Bio.SeqUtils.CheckSum import crc64
    conn.executemany("INSERT OR REPLACE INTO wb_transcripts VALUES (?,?,?,?,?,?)", [
        (name, a.get("gene"), a.get("locus"), a.get("uniprot"), crc64(seq),
         pid_by_crc.get(crc64(seq))) for name, (a, seq) in wb.items()])


def write_result(conn, res, stage_set, config_id):
    """Insert one parsed protein's rows."""
    pid = res["pid"]
    if "ingest" in stage_set:
        conn.executemany("INSERT INTO features VALUES (?,?,?,?,?,?)",
                         [(pid, *r) for r in res["features"]])
        conn.executemany("INSERT INTO residues VALUES (?,?,?,?,?,?,?)",
                         [(pid, *r) for r in res["residues"]])
        if res["summary"]:
            conn.execute("INSERT OR REPLACE INTO conservation_summary VALUES (?,?,?,?,?,?,?,?,?)",
                         (pid, *res["summary"]))
        conn.executemany("INSERT INTO isoforms VALUES (?,?,?,?,?,?,?)",
                         [(pid, *r) for r in res["isoforms"]])
        conn.executemany("INSERT INTO isoform_spans VALUES (?,?,?,?,?,?)",
                         [(pid, *r) for r in res["spans"]])
    if "scores" in stage_set:
        conn.executemany("INSERT INTO site_scores VALUES (?,?,?,?,?,?,?)",
                         [(config_id, pid, *r) for r in res["site_scores"]])
        conn.executemany("INSERT INTO tag_sites VALUES (?,?,?,?,?)",
                         [(pid, pos, "suggested", rank, config_id) for pos, rank in res["suggested"]])


def register_scores_config(conn):
    """Insert the current scores.config.json (once, by hash); returns (config_id, config)."""
    from utils.scoring import load_scoring_config, score_max

    config = load_scoring_config()
    text = json.dumps(config, sort_keys=True)
    digest = hashlib.sha1(text.encode()).hexdigest()[:12]
    row = conn.execute("SELECT config_id FROM score_configs WHERE hash=?", (digest,)).fetchone()
    if row:
        return row[0], config
    cur = conn.execute("INSERT INTO score_configs (hash, json, criteria, max_score) "
                       "VALUES (?,?,?,?)",
                       (digest, text, json.dumps([c["key"] for c in config["criteria"]]),
                        score_max(config)))
    return cur.lastrowid, config


def load_curated_sites(conn, proteins, paths):
    """Curated tag sites from JSON lists (gene, allele, wormbase_id, site, flank_n, flank_c).

    Each entry attaches to the protein of that WormBase gene whose sequence contains
    flank_n + flank_c with the tag between them AT the stated residue, so a curated site
    only lands on an isoform it genuinely fits. Entries that fit nothing are kept with
    pid NULL and verified=0 rather than dropped.
    """
    by_gene = {}
    for p in proteins:
        if p.get("wb_gene"):
            by_gene.setdefault(p["wb_gene"], []).append(p)
    conn.execute("DELETE FROM curated_sites")
    conn.execute("DELETE FROM tag_sites WHERE source = 'curated'")
    n_ok = n_all = 0
    for path in paths:
        if not Path(path).exists():
            continue
        for e in json.load(open(path)):
            n_all += 1
            needle = (e.get("flank_n") or "") + (e.get("flank_c") or "")
            site, hit = int(e["site"]), None
            for p in by_gene.get(e.get("wormbase_id"), []):
                i = p["sequence"].find(needle) if needle else -1
                # tag goes after the last N-flank residue, i.e. after residue i + len(flank_n)
                if i >= 0 and i + len(e.get("flank_n") or "") == site:
                    hit = p
                    break
            n_ok += hit is not None
            pid = hit["pid"] if hit else None
            conn.execute("INSERT INTO curated_sites VALUES (?,?,?,?,?,?,?,?)",
                         (pid, e.get("wormbase_id"), e.get("gene"), e.get("allele"), site,
                          int(hit is not None), e.get("isoform"), Path(path).name))
            if hit:
                conn.execute("INSERT INTO tag_sites VALUES (?,?,?,?,?)",
                             (pid, site, "curated", None, None))
    return n_ok, n_all


def write_status(conn, proteins, status):
    """run_status rows from the last status line per accession."""
    rows = []
    for p in proteins:
        rec = status.get(p["id"])
        if rec:
            t = rec["tasks"]
            rows.append((p["pid"], rec["status"], t.get("domains"), t.get("plddt"),
                         t.get("modifications"), t.get("uniprot"), t.get("scores"),
                         t.get("blast"), rec.get("timestamp")))
    conn.executemany("INSERT OR REPLACE INTO run_status VALUES (?,?,?,?,?,?,?,?,?)", rows)


def merged_stages(conn, stages):
    """Stage names recorded in meta plus `stages`, in pipeline order."""
    row = conn.execute("SELECT value FROM meta WHERE key = 'stages'").fetchone() \
        if table_in(conn, "meta") else None
    have = set((row[0] if row else "").split(",")) | set(stages)
    return ",".join(s for s in STAGE_TABLES if s in have)


def write_meta(conn, cfg, stages):
    """Record build provenance and every convention applied."""
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=_REPO_ROOT,
                                capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = ""
    meta = {
        "built": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": commit,
        "stages": merged_stages(conn, stages), "run_dir": str(cfg["run_dir"]),
        "scale": json.dumps(SCALE),
        "modification_stop": "inclusive (raw file stop - 1); the app scores the raw file",
        "hydrophobicity_pos": "row i of _scores.tsv -> pos i+1 (window start), as in the app",
        "conservation_pos": "row i of .jsd -> pos i+1; -1000 gap sentinel -> NULL",
        "ingest_rule": "task file read only if its last status is ok",
        "sequence": "proteins.sequence is the batch's own {id}.fa; sequence_differs=1 marks "
                    "proteins whose analysed sequence differs from local_store (a different "
                    "UniProt release), and those have no wb_transcript link",
    }
    for name, path in (("pfam_scan_cache", "pfam_scan_cache/_meta.json"),
                       ("conservation_hits", "conservation_hits/_meta.json"),
                       ("topology_cache", "topology_cache/_meta.json")):
        f = cfg["ref_dir"] / path
        if f.exists():
            meta[name] = json.dumps(json.load(open(f)))[:2000]
    conn.executemany("INSERT OR REPLACE INTO meta VALUES (?,?)", list(meta.items()))


# ── stages ────────────────────────────────────────────────────────────────────

def make_jobs(proteins, status, cfg, stages, scores_config, topo_dir):
    """One job dict per protein for the worker pool."""
    jobs = []
    for p in proteins:
        topo = None
        for name in p["wb_names"]:
            f = topo_dir / "{}_topology.txt".format(name)
            if f.exists():
                topo = str(f)
                break
        jobs.append({"pid": p["pid"], "id": p["id"], "sequence": p["sequence"],
                     "folder": str(cfg["run_dir"] / p["id"]), "run_dir": str(cfg["run_dir"]),
                     "tasks": status.get(p["id"], {}).get("tasks", {}), "stages": stages,
                     "topology_file": topo, "scores_config": scores_config})
    return jobs


def fetch_reagent_rows(conn, jobs):
    """Attach each job's per-residue (nearest guide distance, splice distances) rows."""
    by_pid = {j["pid"]: [] for j in jobs}
    for batch in range(0, len(jobs), 500):
        pids = [j["pid"] for j in jobs[batch:batch + 500]]
        # only residues with a guide, as in the reagents TSV; negatives are sentinels
        for pid, res, dist, d5, d3 in conn.execute(
                "SELECT i.pid, i.residue_index, MIN(CASE WHEN s.distance >= 0 THEN s.distance END),"
                " i.dist_to_5p_splice, i.dist_to_3p_splice FROM insertion_sites i"
                " JOIN site_guides s ON s.pid = i.pid AND s.residue_index = i.residue_index"
                " WHERE i.pid IN ({}) GROUP BY i.pid, i.residue_index".format(
                    ",".join("?" * len(pids))), pids):
            by_pid[pid].append((res, dist, d5, d3))
    for j in jobs:
        j["reagents"] = by_pid[j["pid"]]


def run_parse(conn, jobs, stages, workers, chunk, config_id, update=False, progress=True):
    """Parse jobs in worker processes and insert their rows; returns the error list."""
    stage_set = set(stages)
    chunks = [jobs[i:i + chunk] for i in range(0, len(jobs), chunk)]
    # scoring sees the designed reagents when the database has any
    with_reagents = "scores" in stage_set and \
        conn.execute("SELECT 1 FROM site_guides LIMIT 1").fetchone() is not None
    errors, done, t0, busy = [], 0, time.time(), 0.0
    pending, queue = set(), iter(chunks)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        # a bounded window keeps only a few chunks' reagent rows in memory at once
        while True:
            while len(pending) < workers * 2:
                c = next(queue, None)
                if c is None:
                    break
                if with_reagents:
                    fetch_reagent_rows(conn, c)
                pending.add(pool.submit(parse_chunk, c))
            if not pending:
                break
            finished = next(as_completed(pending))
            pending.discard(finished)
            for res in finished.result():
                if update:
                    for s in stage_set:
                        delete_children(conn, res["pid"], s)
                write_result(conn, res, stage_set, config_id)
                errors.extend("{}: {}".format(res["pid"], e) for e in res["errors"])
                busy += res["seconds"]
                done += 1
            conn.commit()
            if progress and (done % (chunk * 5) < chunk or done == len(jobs)):
                print("[build_proteome_db] {}/{} proteins ({:.0f}s)".format(
                    done, len(jobs), time.time() - t0), flush=True)
    run_parse.busy_seconds = busy
    return errors


def build_core(conn, proteins, wb, status, cfg, stages, workers, update=False, everyone=None):
    """Ingest/scores stages into an open connection; returns (errors, seconds).

    `proteins` are the ones to parse; `everyone` (default: the same list) is every protein,
    used for the proteins/wb_transcripts tables so a subset update cannot break their links.
    """
    t0 = time.time()
    create_schema(conn)
    upgrade_schema(conn)
    insert_proteins(conn, everyone or proteins, wb, status)
    write_status(conn, proteins, status)
    config_id, scores_config = (None, None)
    if "scores" in stages:
        config_id, scores_config = register_scores_config(conn)
        n_ok, n_all = load_curated_sites(conn, everyone or proteins, cfg["curated_sites"])
        print("[build_proteome_db] curated tag sites: {}/{} verified against sequence".format(
            n_ok, n_all), flush=True)
    conn.commit()
    topo_dir = cfg["ref_dir"] / "topology_cache"
    jobs = make_jobs(proteins, status, cfg, stages, scores_config, topo_dir)
    # never fewer chunks than workers, or a small run leaves most of the pool idle
    chunk = max(1, min(cfg["chunk_size"], -(-len(jobs) // workers)))
    errors = run_parse(conn, jobs, stages, workers, chunk, config_id, update)
    return errors, time.time() - t0


def connect_build(path):
    """A fast-and-unsafe build connection: the file is rebuilt from scratch if interrupted."""
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA cache_size=-200000")
    return conn


# ── reagent stage ─────────────────────────────────────────────────────────────

def _int(text):
    """int(text), or None for blank/NaN."""
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _float(text):
    """float(text), or None for blank/NaN."""
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def parse_reagent_folder(job):
    """Compact rows from one protein's reagent files (arms are deliberately dropped).

    Reads {id}_reagents.tsv, .genotyping.tsv and {id}_genewise.region.json. A guide that
    serves many insertion sites becomes one `guides` row plus one small `site_guides`
    row per site; repeated mutation_desc text is deduplicated by the writer.
    """
    folder, id_ = Path(job["folder"]), job["id"]
    out = {"pid": job["pid"], "ok": False, "error": None, "region": None, "sites": [],
           "guides": [], "site_guides": [], "genotyping": [], "tsv_bytes": 0}
    try:
        out["region"] = json.load(open(folder / (id_ + "_genewise.region.json")))
        reagents = folder / (id_ + "_reagents.tsv")
        csv.field_size_limit(1 << 28)
        guide_idx, seen_sites = {}, set()
        with open(reagents, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                res = int(r["residue_index"])
                if res not in seen_sites:
                    seen_sites.add(res)
                    out["sites"].append((res, r["amino_acid"], int(r["insert_pos"]),
                                         int(r["exon_index"]), int(r["is_split_codon"] == "True"),
                                         int(r["dist_to_5p_splice"]), int(r["dist_to_3p_splice"])))
                key = (r["guide_strand"], int(r["pam_fwd_start"]), r["spacer"])
                if key not in guide_idx:
                    guide_idx[key] = len(out["guides"])
                    out["guides"].append((r["guide_strand"], r["spacer"], r["pam_seq"],
                                          int(r["pam_fwd_start"]), int(r["cut_pos"]),
                                          _float(r["rs3_score"]), _float(r["rs3_percentile"]),
                                          _int(r["offtarget_count"]), _int(r["offtarget_identical"]),
                                          r["offtarget_status"], r["offtarget_detail"] or None))
                out["site_guides"].append((res, guide_idx[key], int(r["distance"]),
                                           r["pam_in_arm"], r["recut_block_method"],
                                           r["mutation_desc"] or None))
        out["tsv_bytes"] = reagents.stat().st_size
        geno = folder / (id_ + "_reagents.genotyping.tsv")
        if geno.exists():
            with open(geno, newline="") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    out["genotyping"].append((
                        int(r["residue_index"]), r["amplicon_type"], r["fwd_seq"],
                        _float(r["fwd_tm"]), r["rev_seq"], _float(r["rev_tm"]),
                        _int(r["product_size"]), _int(r.get("offtarget_amplicons")),
                        r.get("offtarget_detail") or None))
            out["tsv_bytes"] += geno.stat().st_size
        out["ok"] = True
    except Exception as e:   # a failed design is recorded, never fatal
        out["error"] = "{}: {}".format(type(e).__name__, e)
    return out


def parse_reagent_chunk(jobs):
    """Worker entry point for the reagent files."""
    return [parse_reagent_folder(j) for j in jobs]


def write_reagent_result(conn, res, status_rec, caches):
    """Insert one protein's reagent rows; caches = {"edits": {text: id}, "next_gid": [int]}."""
    pid = res["pid"]
    delete_children(conn, pid, "reagents")
    seconds = status_rec.get("seconds") if status_rec else None
    if not res["ok"]:
        # the design step's own message (e.g. "no GFF3 transcript ...") beats the follow-on
        # "file not found" from parsing outputs that were never written
        task_msg = (status_rec or {}).get("tasks", {}).get("reagents")
        err = task_msg if task_msg and task_msg != "ok" else (res["error"] or "failed")
        conn.execute("INSERT OR REPLACE INTO reagent_status VALUES (?,?,?,?)",
                     (pid, "failed", seconds, str(err)[:300]))
        return
    r = res["region"]
    conn.execute("INSERT OR REPLACE INTO regions VALUES (?,?,?,?,?,?,?,?)",
                 (pid, r["transcript_id"], r["chrom"], r["start"], r["stop"], r["strand"],
                  r["flank_5p"], r["flank_3p"]))
    conn.executemany("INSERT INTO insertion_sites VALUES (?,?,?,?,?,?,?,?)",
                     [(pid, *row) for row in res["sites"]])
    gids = []
    for g in res["guides"]:
        gid = caches["next_gid"][0]
        caches["next_gid"][0] += 1
        gids.append(gid)
        conn.execute("INSERT INTO guides VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", (gid, pid, *g))
    sg = []
    for res_idx, gi, dist, in_arm, method, desc in res["site_guides"]:
        edit_id = None
        if desc:
            edit_id = caches["edits"].get(desc)
            if edit_id is None:
                edit_id = conn.execute("INSERT INTO edits (text) VALUES (?)", (desc,)).lastrowid
                caches["edits"][desc] = edit_id
        sg.append((pid, res_idx, gids[gi], dist, in_arm, method, edit_id))
    conn.executemany("INSERT OR REPLACE INTO site_guides VALUES (?,?,?,?,?,?,?)", sg)
    conn.executemany("INSERT INTO genotyping_primers VALUES (?,?,?,?,?,?,?,?,?,?)",
                     [(pid, *row) for row in res["genotyping"]])
    conn.execute("INSERT OR REPLACE INTO reagent_status VALUES (?,?,?,?)",
                 (pid, "success", seconds, None))


def reagent_caches(conn):
    """Writer-side counters: the next guide id and the mutation_desc -> edit_id map."""
    nxt = (conn.execute("SELECT MAX(gid) FROM guides").fetchone()[0] or 0) + 1
    edits = dict((t, i) for i, t in conn.execute("SELECT edit_id, text FROM edits"))
    return {"next_gid": [nxt], "edits": edits}


def check_reagent_backends():
    """The reagent stage needs the local genewise backend; fail early with the fix."""
    import providers

    if providers.backend_mode("genewise", default="remote") != "bulk":
        sys.exit("The reagent stage needs backends.genewise = \"bulk\". Run it with "
                 "TAGSITES_BATCH_CONFIG=batch.config.local.json (the config that holds the "
                 "local backends).")


def run_reagent_chunk(conn, items, cfg, workers, caches, keep_files=False):
    """Design, ingest and delete one chunk of proteins; returns stats.

    items = [(pid, id, sequence)]. Designs run through proteome_run.main into a scratch
    folder, are parsed into compact rows, and the folder is removed, which bounds transient
    disk to one chunk instead of ~20 kB of TSV per residue for the whole proteome.
    """
    import proteome_run

    scratch = Path(cfg["scratch_dir"]) / "chunk_{}".format(items[0][0])
    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True, exist_ok=True)
    os.environ["TAGSITES_BLAST_THREADS"] = "1"   # the pool already has one process per core
    t0 = time.time()
    proteome_run.main(str(scratch), task_types=["reagents"], workers=workers,
                      rows=[(i, seq) for _, i, seq in items], presearch=False, force=True)
    wall = time.time() - t0
    status = read_last_status(scratch)
    jobs = [{"pid": pid, "id": id_, "folder": str(scratch / id_)} for pid, id_, _ in items]
    chunks = [jobs[i:i + 25] for i in range(0, len(jobs), 25)]
    stats = {"n": len(items), "ok": 0, "wall": wall, "tsv_bytes": 0, "seconds": 0.0,
             "ok_residues": 0, "ok_seconds": 0.0}
    seq_len = {pid: len(seq) for pid, _, seq in items}
    id_of = {pid: id_ for pid, id_, _ in items}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for results in pool.map(parse_reagent_chunk, chunks):
            for res in results:
                rec = status.get(id_of[res["pid"]])
                write_reagent_result(conn, res, rec, caches)
                stats["tsv_bytes"] += res["tsv_bytes"]
                stats["seconds"] += (rec or {}).get("seconds") or 0
                if res["ok"]:
                    stats["ok"] += 1
                    stats["ok_residues"] += seq_len[res["pid"]]
                    stats["ok_seconds"] += (rec or {}).get("seconds") or 0
    conn.commit()
    if not keep_files:
        shutil.rmtree(scratch, ignore_errors=True)
    return stats


def reagent_todo(conn):
    """(pid, id, sequence) still needing reagents; every protein is designed on its own."""
    done = {r[0] for r in conn.execute("SELECT pid FROM reagent_status WHERE status='success'")}
    # identical sequences are not merged: paralogs sit at different genomic loci
    return [(pid, id_, seq) for pid, id_, seq in conn.execute(
        "SELECT pid, id, sequence FROM proteins ORDER BY pid") if pid not in done]


def build_reagents(args, cfg):
    """Resumable reagent stage: works on the local build copy, committing every chunk."""
    check_reagent_backends()
    if not (args.limit or args.accessions):
        gate(cfg, args)
    work = cfg["build_dir"] / "proteome_building.sqlite3"
    if not work.exists():
        if not cfg["path"].exists():
            sys.exit("Build the core database first (`--stage all`); the reagent stage adds to it.")
        shutil.copy(cfg["path"], work)
    conn = sqlite3.connect(work)
    conn.execute("PRAGMA synchronous=NORMAL")
    create_schema(conn)
    caches = reagent_caches(conn)
    todo = reagent_todo(conn)
    if args.accessions:
        wanted = set(args.accessions.split(","))
        todo = [t for t in todo if t[1] in wanted]
    elif args.limit:
        todo = todo[: args.limit]
    workers = args.workers or cfg["max_workers"]
    print("[reagents] {} proteins to design, chunks of {}, {} workers".format(
        len(todo), cfg["reagent_chunk"], workers), flush=True)
    t0, done = time.time(), 0
    for i in range(0, len(todo), cfg["reagent_chunk"]):
        stats = run_reagent_chunk(conn, todo[i:i + cfg["reagent_chunk"]], cfg, workers, caches)
        done += stats["n"]
        print("[reagents] {}/{} proteins ({:.0f} min, chunk ok {}/{})".format(
            done, len(todo), (time.time() - t0) / 60, stats["ok"], stats["n"]), flush=True)
    create_indexes(conn)
    write_meta(conn, cfg, ["reagents"])
    conn.commit()
    conn.close()
    tmp = cfg["path"].with_suffix(".sqlite3.tmp")
    shutil.copy(work, tmp)
    os.replace(tmp, cfg["path"])
    print("reagents stage written to", cfg["path"])


# ── size estimate ─────────────────────────────────────────────────────────────

def table_sizes(conn):
    """{table: bytes} including each table's indexes, from SQLite's dbstat."""
    rows = conn.execute("SELECT m.tbl_name, SUM(d.pgsize) FROM dbstat d "
                        "JOIN sqlite_master m ON m.name = d.name GROUP BY m.tbl_name").fetchall()
    return {t: int(b) for t, b in rows}


def pick_sample(proteins, status, n, seed):
    """A seeded random sample that always includes structure-less and isoform proteins."""
    rng = random.Random(seed)
    strata = {}
    for p in proteins:
        has = status.get(p["id"], {}).get("tasks", {}).get("plddt") == "ok"
        strata.setdefault((p["kind"], has), []).append(p)
    chosen = {}
    for members in strata.values():
        for p in rng.sample(members, min(len(members), 5)):
            chosen[p["id"]] = p
    pool = [p for p in proteins if p["id"] not in chosen]
    for p in rng.sample(pool, min(len(pool), max(0, n - len(chosen)))):
        chosen[p["id"]] = p
    return list(chosen.values())


def estimate(args, cfg, proteins, wb, status):
    """Build a sample database, project every table to the full proteome, and report."""
    stages = ["ingest", "scores"]
    sample = pick_sample(proteins, status, args.sample, args.seed)
    total_res = sum(len(p["sequence"]) for p in proteins)
    sample_res = sum(len(p["sequence"]) for p in sample)
    cfg["build_dir"].mkdir(parents=True, exist_ok=True)
    path = cfg["build_dir"] / "proteome_estimate.sqlite3"
    path.unlink(missing_ok=True)
    conn = connect_build(path)
    errors, secs = build_core(conn, sample, wb, status, cfg, stages, args.workers or cfg["max_workers"])
    create_indexes(conn)
    conn.commit()
    sizes = table_sizes(conn)
    conn_rows = {t: conn.execute("SELECT COUNT(*) FROM {}".format(t)).fetchone()[0]
                 for t in sizes}
    conn.close()

    proj, rows_out = {}, []
    for t, b in sorted(sizes.items(), key=lambda kv: -kv[1]):
        if conn_rows.get(t, 0) == 0:
            continue   # an empty table is just a minimum page; not worth projecting
        how = SCALES_WITH.get(t, "proteins")
        factor = {"proteins": len(proteins) / len(sample), "residues": total_res / sample_res,
                  "fixed": 1.0}[how]
        proj[t] = b * factor
        rows_out.append((t, how, b, proj[t]))
    reagent = estimate_reagents(args, cfg, proteins, wb, status, total_res) \
        if args.reagent_sample else {}
    for t, b in reagent.get("tables", {}).items():
        proj[t] = b
    total = sum(proj.values())

    print("\n== size estimate ({} sample proteins of {}, {:,} of {:,} residues) ==".format(
        len(sample), len(proteins), sample_res, total_res))
    print("{:<22}{:<10}{:>14}{:>16}".format("table", "scales", "sample", "projected"))
    for t, how, b, p in rows_out:
        print("{:<22}{:<10}{:>12.2f}MB{:>14.1f}MB".format(t, how, b / 1e6, p / 1e6))
    for t, b in reagent.get("tables", {}).items():
        print("{:<22}{:<10}{:>14}{:>14.1f}MB".format(t, "reagents", "(see below)", b / 1e6))
    print("{:<32}{:>28.2f} GB".format("TOTAL projected database", total / 1e9))
    workers = args.workers or cfg["max_workers"]
    busy = run_parse.busy_seconds / len(sample)
    print("ingest+scores: {:.2f} worker-seconds/protein (sample wall {:.0f}s) -> ~{:.0f} min "
          "for all at {} workers".format(busy, secs, busy * len(proteins) / workers / 60,
                                         workers))
    if reagent:
        print("reagents: {:.1f} s/protein/worker -> ~{:.1f} h for {} proteins at {} workers".format(
            reagent["seconds_per_protein"], reagent["hours"], reagent["n_eligible"],
            args.workers or cfg["max_workers"]))
        print("          transient TSV disk per chunk of {}: ~{:.1f} GB".format(
            cfg["reagent_chunk"], reagent["tsv_gb_per_chunk"]))
    for d in {str(cfg["build_dir"]), str(cfg["path"].parent)}:
        free = shutil.disk_usage(d).free
        print("free space at {}: {:.1f} GB".format(d, free / 1e9))
    if errors:
        print("{} parse warning(s) in the sample, e.g. {}".format(len(errors), errors[:3]))
    out = {"made": time.time(), "projected_bytes": total, "per_table": proj,
           "n_proteins": len(proteins), "sample": len(sample),
           "reagents": {k: v for k, v in reagent.items() if k != "tables"}}
    (cfg["build_dir"] / "proteome_db.estimate.json").write_text(json.dumps(out, indent=1))
    path.unlink(missing_ok=True)
    return out


def estimate_reagents(args, cfg, proteins, wb, status, total_res):
    """Design reagents for a small random sample into a temp DB; project size and time."""
    check_reagent_backends()
    rng = random.Random(args.seed + 1)
    pool = [p for p in proteins if len(p["sequence"]) <= args.reagent_max_len]
    chosen = rng.sample(pool, min(args.reagent_sample, len(pool)))
    path = cfg["build_dir"] / "proteome_estimate_reagents.sqlite3"
    path.unlink(missing_ok=True)
    conn = connect_build(path)
    create_schema(conn)
    insert_proteins(conn, chosen, wb, status)
    caches = reagent_caches(conn)
    workers = min(args.workers or cfg["max_workers"], len(chosen))
    stats = run_reagent_chunk(conn, [(p["pid"], p["id"], p["sequence"]) for p in chosen],
                              cfg, workers, caches, keep_files=False)
    create_indexes(conn)
    conn.commit()
    sizes = table_sizes(conn)
    failures = conn.execute("SELECT p.id, p.seq_length, r.error FROM reagent_status r "
                            "JOIN proteins p USING (pid) WHERE r.status != 'success'").fetchall()
    conn.close()
    path.unlink(missing_ok=True)
    for id_, length, err in failures:
        print("reagent sample failure: {} ({} aa): {}".format(id_, length, (err or "").split("\n")[0][:200]))
    if not stats["ok"]:
        print("reagent sample: no protein succeeded; cannot project reagent size")
        return {}
    ok_rate = stats["ok"] / stats["n"]
    reagent_bytes = sum(b for t, b in sizes.items() if t in STAGE_TABLES["reagents"])
    per_res = reagent_bytes / stats["ok_residues"]
    sec_per_res = stats["ok_seconds"] / stats["ok_residues"]
    n_eligible = int(len(proteins) * ok_rate)
    tables = {t: b / stats["ok_residues"] * total_res * ok_rate
              for t, b in sizes.items() if t in STAGE_TABLES["reagents"]}
    longest = sum(1 for p in proteins if len(p["sequence"]) > args.reagent_max_len)
    print("\nreagent sample: {}/{} designed ok, {:,} residues, {:.1f}s worker-time per protein, "
          "{:.3f}s per residue; sample capped at {} aa ({} longer proteins not sampled)".format(
              stats["ok"], stats["n"], stats["ok_residues"],
              stats["ok_seconds"] / stats["ok"], sec_per_res, args.reagent_max_len, longest))
    return {"tables": tables, "seconds_per_protein": stats["ok_seconds"] / stats["ok"],
            "hours": sec_per_res * total_res * ok_rate / (args.workers or cfg["max_workers"]) / 3600,
            "n_eligible": n_eligible, "ok_rate": ok_rate,
            "tsv_gb_per_chunk": stats["tsv_bytes"] / stats["n"] * cfg["reagent_chunk"] / 1e9,
            "bytes_per_residue": per_res}


# ── full build ────────────────────────────────────────────────────────────────

def gate(cfg, args):
    """Refuse a full build without a fresh estimate, --yes, and enough free space."""
    est_file = cfg["build_dir"] / "proteome_db.estimate.json"
    if not est_file.exists():
        sys.exit("No size estimate found. Run `build_proteome_db.py --estimate` first.")
    est = json.loads(est_file.read_text())
    if time.time() - est["made"] > 7 * 86400:
        sys.exit("The size estimate is over a week old; re-run `--estimate`.")
    free = shutil.disk_usage(cfg["build_dir"]).free
    need = est["projected_bytes"] * 1.3
    print("projected {:.2f} GB (x1.3 headroom = {:.2f} GB); free on build dir {:.1f} GB".format(
        est["projected_bytes"] / 1e9, need / 1e9, free / 1e9))
    if need > free:
        sys.exit("Not enough free space on {} for this build.".format(cfg["build_dir"]))
    if not args.yes:
        sys.exit("Re-run with --yes to build.")


def build(args, cfg, proteins, wb, status):
    """Build or update the database for the requested stages."""
    if args.stage == "reagents":
        return build_reagents(args, cfg)
    subset = bool(args.limit or args.accessions)
    everyone = proteins
    if args.accessions:
        wanted = set(args.accessions.split(","))
        proteins = [p for p in proteins if p["id"] in wanted]
    elif args.limit:
        proteins = proteins[: args.limit]
    if not subset:
        gate(cfg, args)
    stages = ["ingest", "scores"] if args.stage in ("all", "core") else [args.stage]
    cfg["build_dir"].mkdir(parents=True, exist_ok=True)
    work = cfg["build_dir"] / "proteome_building.sqlite3"
    if args.update and cfg["path"].exists():
        shutil.copy(cfg["path"], work)       # continue from the finished file
    else:
        work.unlink(missing_ok=True)
    conn = connect_build(work)
    errors, secs = build_core(conn, proteins, wb, status, cfg, stages,
                              args.workers or cfg["max_workers"], update=args.update,
                              everyone=everyone)
    create_indexes(conn)
    write_meta(conn, cfg, stages)
    conn.commit()
    conn.close()
    cfg["path"].parent.mkdir(parents=True, exist_ok=True)
    tmp = cfg["path"].with_suffix(".sqlite3.tmp")
    shutil.copy(work, tmp)
    os.replace(tmp, cfg["path"])
    work.unlink(missing_ok=True)
    print("built {} ({:.2f} GB) in {:.0f}s; {} parse warning(s)".format(
        cfg["path"], cfg["path"].stat().st_size / 1e9, secs, len(errors)))
    for e in errors[:10]:
        print("  warning:", e)
    return errors


def main(argv=None):
    """CLI entry point."""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--estimate", action="store_true",
                   help="sample, measure and project the database size; build nothing large")
    p.add_argument("--stage", default="all", choices=["all", "core", "ingest", "scores",
                                                      "reagents"])
    p.add_argument("--update", action="store_true",
                   help="refresh only the selected proteins in the existing database")
    p.add_argument("--yes", action="store_true", help="confirm a full build after --estimate")
    p.add_argument("--accessions", default=None, help="comma-separated ids (subset build)")
    p.add_argument("--limit", type=int, default=None, help="first N proteins (subset build)")
    p.add_argument("--sample", type=int, default=300, help="--estimate sample size")
    p.add_argument("--reagent-sample", type=int, default=0,
                   help="--estimate: also design reagents for this many proteins")
    p.add_argument("--reagent-max-len", type=int, default=2000,
                   help="--estimate: longest protein (aa) drawn for the reagent sample")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--out", default=None, help="write the database here instead of the config path")
    args = p.parse_args(argv)

    cfg = db_config()
    if args.out:
        cfg["path"] = Path(args.out)
    proteins, wb = load_proteins(cfg)
    status = read_last_status(cfg["run_dir"])
    if args.estimate:
        estimate(args, cfg, proteins, wb, status)
    else:
        build(args, cfg, proteins, wb, status)


if __name__ == "__main__":
    main()
