"""
proteome_db.py

Read-only lookup over the SQLite file built by scripts/build_proteome_db.py: find a protein
by accession, gene name or WormBase name, then pull its features, per-residue tracks,
tag-site scores, suggested sites, isoforms, run status and reagents as DataFrames, or run
any SQL with query(). The file is opened read-only and immutable because SQLite locking is
unreliable on network volumes.

Numbers in the database are scaled integers (build_proteome_db.SCALE); the get_* functions
return them unscaled. query() returns the raw stored values.

    python scripts/proteome_db.py lookup trxr-1
    python scripts/proteome_db.py sql "SELECT track, COUNT(*) FROM features GROUP BY track"

Matt Rich, 2026
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(_REPO_ROOT))

from build_proteome_db import SCALE  # noqa: E402

# the reagent parameters run_reagents / the batch stage designed with
ARM_LENGTH, PAM, SEED_LEN = 1000, "NGG", 15


def default_path():
    """The database path from batch.config.json's proteome_db block."""
    from providers import _load_config

    cfg = _load_config().get("proteome_db", {})
    p = Path(cfg.get("path") or "data/runs/proteome_v1/proteome.sqlite3")
    return p if p.is_absolute() else _REPO_ROOT / p


def open_db(path=None):
    """A read-only, immutable connection (no locking, safe on a network volume)."""
    path = Path(path) if path else default_path()
    if not path.exists():
        raise FileNotFoundError("{} not found; build it with scripts/build_proteome_db.py".format(path))
    return sqlite3.connect("file:{}?mode=ro&immutable=1".format(path), uri=True)


def query(sql, params=(), conn=None):
    """Run any SELECT and return a DataFrame of the raw stored values."""
    own = conn is None
    conn = conn or open_db()
    try:
        return pd.read_sql_query(sql, conn, params=params)
    finally:
        if own:
            conn.close()


def find_protein(term, conn=None):
    """Proteins matching a term: accession, id, gene name, WormBase gene / transcript / locus.

    Case-insensitive and exact (so `trxr-1` finds both trxr-1 proteins); returns one row per
    protein without the sequence.
    """
    sql = """SELECT DISTINCT p.pid, p.id, p.kind, p.accession, p.gene_name, p.wormbase_gene,
                    p.wormbase_transcript, p.wb_transcript, p.wb_gene, p.seq_length,
                    p.has_structure, p.sequence_differs
             FROM proteins p LEFT JOIN wb_transcripts w ON w.pid = p.pid
             WHERE p.id = :t COLLATE NOCASE OR p.accession = :t COLLATE NOCASE
                OR p.gene_name = :t COLLATE NOCASE OR p.wormbase_gene = :t COLLATE NOCASE
                OR p.wormbase_transcript = :t COLLATE NOCASE OR p.wb_transcript = :t COLLATE NOCASE
                OR p.wb_gene = :t COLLATE NOCASE OR w.locus = :t COLLATE NOCASE
             ORDER BY p.pid"""
    own = conn is None
    conn = conn or open_db()
    try:
        return pd.read_sql_query(sql, conn, params={"t": term})
    finally:
        if own:
            conn.close()


def get_features(pid, tracks=None, conn=None):
    """Range features (domain / modification / topology / uniprot / hydrophobic_patch)."""
    sql = "SELECT track, source, start, stop, description FROM features WHERE pid = ?"
    params = [pid]
    if tracks:
        sql += " AND track IN ({})".format(",".join("?" * len(tracks)))
        params += list(tracks)
    return query(sql + " ORDER BY start", params, conn)


def get_residues(pid, conn=None):
    """Per-residue tracks, unscaled, with the amino acid: pos, aa, conservation, ..."""
    own = conn is None
    conn = conn or open_db()
    try:
        seq = conn.execute("SELECT sequence FROM proteins WHERE pid = ?", (pid,)).fetchone()
        df = pd.read_sql_query("SELECT pos, conservation, hydrophobicity, plddt, rsasa, "
                               "struct_hydro FROM residues WHERE pid = ? ORDER BY pos",
                               conn, params=(pid,))
    finally:
        if own:
            conn.close()
    for col in ("conservation", "hydrophobicity", "plddt", "rsasa", "struct_hydro"):
        df[col] = df[col] / SCALE[col]
    if seq:
        df.insert(1, "aa", [seq[0][p - 1] if 0 < p <= len(seq[0]) else None for p in df["pos"]])
    return df


def get_scores(pid, config_id=1, conn=None):
    """Tag-site scores per residue with each criterion decoded to True/False/None."""
    own = conn is None
    conn = conn or open_db()
    try:
        crit = conn.execute("SELECT criteria, max_score FROM score_configs WHERE config_id = ?",
                            (config_id,)).fetchone()
        df = pd.read_sql_query("SELECT pos, score, masked, bits, missing FROM site_scores "
                               "WHERE pid = ? AND config_id = ? ORDER BY pos",
                               conn, params=(pid, config_id))
    finally:
        if own:
            conn.close()
    if crit is None:
        return df
    for i, key in enumerate(json.loads(crit[0])):
        df[key] = [None if (m >> i) & 1 else bool((b >> i) & 1)
                   for b, m in zip(df["bits"], df["missing"])]
    df["score"] = df["score"] / SCALE["score"]
    df["masked"] = df["masked"].astype(bool)
    return df.drop(columns=["bits", "missing"])


def get_sites(pid, conn=None):
    """Suggested (and curated) tag sites for a protein."""
    return query("SELECT pos, source, rank, config_id FROM tag_sites WHERE pid = ? "
                 "ORDER BY source, rank", (pid,), conn)


def get_isoforms(pid, conn=None):
    """Isoforms with their present / skipped / insert spans."""
    isos = query("SELECT * FROM isoforms WHERE pid = ? ORDER BY iso_idx", (pid,), conn)
    spans = query("SELECT iso_idx, kind, start, stop, length FROM isoform_spans WHERE pid = ?",
                  (pid,), conn)
    return isos, spans


def get_status(pid, conn=None):
    """The last batch-run status of each task for a protein."""
    return query("SELECT * FROM run_status WHERE pid = ?", (pid,), conn)


def get_reagents(pid, residue_index=None, conn=None):
    """Reagent rows (one per residue x guide) without arms; see regenerate_arms()."""
    sql = """SELECT s.residue_index, i.amino_acid, i.insert_pos, i.exon_index, i.is_split_codon,
                    i.dist_to_5p_splice, i.dist_to_3p_splice, g.gid, g.strand AS guide_strand,
                    g.spacer, g.pam_seq, g.pam_fwd_start, g.cut_pos, s.distance, s.pam_in_arm,
                    s.recut_block_method, e.text AS mutation_desc, g.rs3_score,
                    g.rs3_percentile, g.offtarget_count, g.offtarget_identical,
                    g.offtarget_status, g.offtarget_detail
             FROM site_guides s
             JOIN guides g ON g.gid = s.gid
             JOIN insertion_sites i ON i.pid = s.pid AND i.residue_index = s.residue_index
             LEFT JOIN edits e ON e.edit_id = s.edit_id
             WHERE s.pid = ?"""
    params = [pid]
    if residue_index is not None:
        sql += " AND s.residue_index = ?"
        params.append(residue_index)
    return query(sql + " ORDER BY s.residue_index, s.distance", params, conn)


def get_genotyping(pid, conn=None):
    """Genotyping primer pairs, one row per residue x amplicon_type."""
    return query("SELECT * FROM genotyping_primers WHERE pid = ? ORDER BY residue_index",
                 (pid,), conn)


def regenerate_arms(pid, residue_index, gid, conn=None, arm_length=ARM_LENGTH, pam=PAM):
    """Rebuild left_arm / right_arm / left_arm_wt / right_arm_wt for one residue x guide.

    The database stores no arms. They are recomputed from the genome with the same
    deterministic functions design_reagents uses: the gene region (transcript + stored
    flank), then disrupt_pam for the stored guide, then exon/intron casing.
    """
    from crispr_util import build_frame_lookup, disrupt_pam
    from design_tag_reagents import _case_arm
    from genome_regions import get_transcript_region

    own = conn is None
    conn = conn or open_db()
    try:
        reg = conn.execute("SELECT transcript_id, flank_5p, flank_3p FROM regions WHERE pid = ?",
                           (pid,)).fetchone()
        site = conn.execute("SELECT insert_pos FROM insertion_sites WHERE pid = ? AND "
                            "residue_index = ?", (pid, residue_index)).fetchone()
        guide = conn.execute("SELECT strand, pam_fwd_start FROM guides WHERE gid = ?",
                             (gid,)).fetchone()
        method = conn.execute("SELECT recut_block_method FROM site_guides WHERE pid = ? AND "
                              "residue_index = ? AND gid = ?",
                              (pid, residue_index, gid)).fetchone()
    finally:
        if own:
            conn.close()
    if not (reg and site and guide and method):
        raise KeyError("no reagent for pid={} residue={} gid={}".format(pid, residue_index, gid))
    # flank is clipped at chromosome ends, so max() of the stored sides reproduces the region
    region = get_transcript_region(reg[0], flank=max(reg[1], reg[2]))
    dna, cds_df = region["dna"], region["cds_df"]
    insert_pos = site[0]
    frame_lookup = build_frame_lookup(cds_df, dna)
    left_start = max(0, insert_pos - arm_length)
    right_end = min(len(dna), insert_pos + arm_length)
    left_raw, right_raw = dna[left_start:insert_pos], dna[insert_pos:right_end]
    left, right = left_raw, right_raw
    # an insertion inside the seed already blocks re-cutting; every other method edits the PAM
    if method[0] != "insertion":
        mutated, _, _ = disrupt_pam(dna, pam, guide[1], guide[0], frame_lookup, seed_len=SEED_LEN)
        left, right = mutated[left_start:insert_pos], mutated[insert_pos:right_end]
    return {"left_arm": _case_arm(left, left_start, frame_lookup),
            "right_arm": _case_arm(right, insert_pos, frame_lookup),
            "left_arm_wt": _case_arm(left_raw, left_start, frame_lookup),
            "right_arm_wt": _case_arm(right_raw, insert_pos, frame_lookup)}


def _summary(pid, conn):
    """A printable one-protein summary for the lookup command."""
    feats = get_features(pid, conn=conn)
    lines = []
    if not feats.empty:
        counts = feats.groupby("track").size().to_dict()
        lines.append("  features: " + ", ".join("{} {}".format(k, v) for k, v in counts.items()))
    sites = get_sites(pid, conn=conn)
    if not sites.empty:
        lines.append("  suggested tag sites: " + ", ".join(
            "{} (#{})".format(r.pos, r.rank) for r in sites.itertuples()))
    status = get_status(pid, conn=conn)
    if not status.empty:
        lines.append("  last run: {}".format(status.iloc[0]["status"]))
    n_res = conn.execute("SELECT COUNT(*) FROM site_guides WHERE pid = ?", (pid,)).fetchone()[0]
    if n_res:
        lines.append("  reagents: {} residue x guide rows".format(n_res))
    return lines


def main(argv=None):
    """CLI: `lookup <term>` or `sql "<query>"`."""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", default=None, help="database path (default: batch.config.json)")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("lookup", help="find proteins by name and summarise them").add_argument("term")
    sub.add_parser("sql", help="run a SELECT").add_argument("statement")
    args = p.parse_args(argv)
    conn = open_db(args.db)
    if args.cmd == "sql":
        print(query(args.statement, conn=conn).to_string(index=False))
        return
    found = find_protein(args.term, conn=conn)
    if found.empty:
        print("no protein matches {!r}".format(args.term))
        return
    for r in found.itertuples():
        print("{} [{}] {} {} aa{}".format(r.id, r.kind, r.gene_name or "", r.seq_length,
                                         "" if r.has_structure else " (no structure)"))
        print("\n".join(_summary(r.pid, conn)))


if __name__ == "__main__":
    main()
