"""
local_store.py

CRC64-keyed local cache/index over the bulk reference data downloaded by
scripts/reference_data.py (Phase A). Backs the local/bulk analysis backends:
looking up a protein by accession or by sequence checksum should be a single
indexed SQLite query instead of a network round-trip.

Schema (one row per canonical UniProt entry in the proteome):
  proteins(accession PK, crc64, sequence, seq_length, gene_name,
           wormbase_gene, wormbase_transcript)
  index on crc64 (mirrors uniprot_api.checksum_lookup()'s CRC64-of-sequence key,
  so identical sequences — including sibling isoforms — collapse to one lookup)

wormbase_transcript/wormbase_gene come straight out of the UniProt JSON's
genes[].orfNames — this is the accession<->WormBase-name cross-reference
scripts/genome_regions.py (Phase A2) needs to look CDS exons up in the GFF3,
with no separate mapping file or per-protein lookup required.

Usage
-----
    python scripts/local_store.py --build          # (re)build the index
    python scripts/local_store.py --stats          # print row counts
"""

import gzip
import json
import sqlite3
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from providers import _load_config as _load_batch_config  # reuses the same batch.config.json


def _reference_dir(cfg=None):
    cfg = cfg or _load_batch_config().get("reference_data", {})
    out_dir = cfg.get("out_dir", "data/reference")
    d = Path(out_dir)
    if not d.is_absolute():
        d = _REPO_ROOT / d
    return d, cfg.get("uniprot_proteome_id", "UP000001940")


def _db_path(cfg=None):
    ref_dir, _ = _reference_dir(cfg)
    return ref_dir / "local_store.sqlite3"


def _first_wormbase_orf_name(entry):
    """Pull the first WormBase-sourced orfNames value (gene) and its evidence id
    (transcript-ish sequence name), e.g. ("C10C5.1g", "C10C5.1") for pezo-1.
    Returns (wormbase_gene_evidence_id, wormbase_transcript) — either may be None.
    """
    for gene in entry.get("genes", []):
        for orf in gene.get("orfNames", []):
            transcript = orf.get("value")
            for ev in orf.get("evidences", []):
                if ev.get("source") == "WormBase":
                    return ev.get("id"), transcript
            if transcript:
                return None, transcript
    return None, None


def _gene_name(entry):
    for gene in entry.get("genes", []):
        name = gene.get("geneName", {}).get("value")
        if name:
            return name
    return None


def build_index(cfg=None, force=False):
    """(Re)build the SQLite index from the downloaded UniProt proteome JSON.

    Reads <out_dir>/<proteome_id>.json (one JSON array of full UniProtKB
    entries, from reference_data.fetch_uniprot()). Raises FileNotFoundError
    with a pointer to that step if it hasn't been run yet.
    """
    ref_dir, proteome_id = _reference_dir(cfg)
    json_path = ref_dir / f"{proteome_id}.json.gz"
    if not json_path.exists():
        raise FileNotFoundError(
            f"{json_path} not found — run `python scripts/reference_data.py "
            "--only uniprot` first."
        )

    db_path = _db_path(cfg)
    if db_path.exists():
        if not force:
            print(f"[skip] {db_path} already present; pass force=True to rebuild")
            return db_path
        db_path.unlink()

    print(f"[local_store] loading {json_path} …")
    with gzip.open(json_path, "rt") as f:
        payload = json.load(f)
    entries = payload["results"] if isinstance(payload, dict) else payload

    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE proteins (
            accession           TEXT PRIMARY KEY,
            crc64               TEXT,
            sequence            TEXT,
            seq_length          INTEGER,
            gene_name           TEXT,
            wormbase_gene       TEXT,
            wormbase_transcript TEXT
        )
    """)
    conn.execute("CREATE INDEX idx_proteins_crc64 ON proteins(crc64)")

    rows = []
    for entry in entries:
        acc = entry.get("primaryAccession")
        seq = entry.get("sequence", {})
        wb_gene, wb_transcript = _first_wormbase_orf_name(entry)
        rows.append((
            acc,
            seq.get("crc64"),
            seq.get("value"),
            seq.get("length"),
            _gene_name(entry),
            wb_gene,
            wb_transcript,
        ))

    conn.executemany(
        "INSERT OR REPLACE INTO proteins VALUES (?, ?, ?, ?, ?, ?, ?)", rows
    )
    conn.commit()
    print(f"[local_store] indexed {len(rows):,} proteins -> {db_path}")
    conn.close()
    return db_path


def open_index(cfg=None):
    """Open a read-only connection to the built index; raises if it hasn't been built."""
    db_path = _db_path(cfg)
    if not db_path.exists():
        raise FileNotFoundError(
            f"{db_path} not found — run `python scripts/local_store.py --build` first."
        )
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)


def lookup_by_accession(accession, conn=None):
    """Return the indexed row for `accession` as a dict, or None if absent."""
    own_conn = conn is None
    conn = conn or open_index()
    try:
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT * FROM proteins WHERE accession = ?", (accession,)
        ).fetchone()
        return dict(row) if row else None
    finally:
        if own_conn:
            conn.close()


def lookup_by_crc64(crc64, conn=None):
    """Return all indexed rows matching a CRC64 checksum (usually 0 or 1; more
    than one indicates two accessions with byte-identical sequences).
    """
    own_conn = conn is None
    conn = conn or open_index()
    try:
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT * FROM proteins WHERE crc64 = ?", (crc64,)
        ).fetchall()
        return [dict(r) for r in rows]
    finally:
        if own_conn:
            conn.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="(re)build the index")
    parser.add_argument("--force", action="store_true", help="rebuild even if the index already exists")
    parser.add_argument("--stats", action="store_true", help="print row counts for the existing index")
    args = parser.parse_args()

    if args.build:
        build_index(force=args.force)
    if args.stats:
        conn = open_index()
        n = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        n_wb = conn.execute("SELECT COUNT(*) FROM proteins WHERE wormbase_transcript IS NOT NULL").fetchone()[0]
        print(f"proteins: {n:,}  (with WormBase transcript name: {n_wb:,})")
    if not args.build and not args.stats:
        parser.error("specify --build and/or --stats")
