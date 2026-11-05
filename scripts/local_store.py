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

  domains(accession, ipr_id, description, external_db_match_id, start, stop)
  index on accession — one row per protein2ipr.dat match (see
  reference_data.fetch_interpro()), backing scripts/domains_bulk.py

  sequence_aliases(crc64, accession, sequence, source)
  every accession known to carry a sequence: the proteome's own entries plus every AFDB
  tarball model (source "afdb_model"). UniProt often holds one sequence under two entries
  (unc-25: Q7K7A3 in the reference proteome, G5EDB7 in AFDB), and services keyed by
  accession may know only one of them; aliases_for_sequence() returns them all

wormbase_transcript/wormbase_gene come straight out of the UniProt JSON's
genes[].orfNames — this is the accession<->WormBase-name cross-reference
scripts/genome_regions.py (Phase A2) needs to look CDS exons up in the GFF3,
with no separate mapping file or per-protein lookup required.

Usage
-----
    python scripts/local_store.py --build          # (re)build the index
    python scripts/local_store.py --build-aliases  # (re)build sequence_aliases
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


def build_domains_index(cfg=None, force=False):
    """Load protein2ipr.filtered.tsv (see reference_data.fetch_interpro()) into
    the `domains` table of the same local_store.sqlite3 database, indexed by
    accession. Raises FileNotFoundError pointing at that step if it hasn't run.
    """
    ref_dir, _ = _reference_dir(cfg)
    tsv_path = ref_dir / "protein2ipr.filtered.tsv"
    if not tsv_path.exists():
        raise FileNotFoundError(
            f"{tsv_path} not found — run `python scripts/reference_data.py "
            "--only interpro` first."
        )

    db_path = _db_path(cfg)
    conn = sqlite3.connect(db_path)
    existing = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='domains'"
    ).fetchone()
    if existing:
        if not force:
            print(f"[skip] domains table already present in {db_path}; pass force=True to rebuild")
            conn.close()
            return db_path
        conn.execute("DROP TABLE domains")

    conn.execute("""
        CREATE TABLE domains (
            accession             TEXT,
            ipr_id                TEXT,
            description           TEXT,
            external_db_match_id  TEXT,
            start                 INTEGER,
            stop                  INTEGER
        )
    """)
    conn.execute("CREATE INDEX idx_domains_accession ON domains(accession)")

    print(f"[local_store] loading {tsv_path} …")
    rows = []
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue
            acc, ipr_id, desc, ext_id, start, stop = parts
            rows.append((acc, ipr_id, desc, ext_id, int(start), int(stop)))

    conn.executemany("INSERT INTO domains VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    print(f"[local_store] indexed {len(rows):,} domain matches -> {db_path}")
    conn.close()
    return db_path


def afdb_model_path(accession, cfg=None):
    """Path of an accession's AFDB model in the extracted tarball (it may not exist)."""
    ref_dir, _ = _reference_dir(cfg)
    return ref_dir / "afdb_pdbs" / f"AF-{accession}-F1-model_v6.pdb.gz"


def _pdb_sequence(path):
    """One-letter sequence of a gzipped AFDB PDB from its SEQRES records (CA atoms if none)."""
    from Bio.SeqUtils import seq1

    seqres, ca = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("SEQRES"):
                seqres.extend(line[19:].split())
            elif line.startswith("ATOM") and line[12:16] == " CA ":
                # SEQRES precedes the coordinates, so reaching atoms with it means it is complete
                if seqres:
                    break
                ca.append(line[17:20])
    return seq1("".join(seqres or ca))


def build_aliases_index(cfg=None, force=False):
    """(Re)build sequence_aliases from the proteins table and every AFDB tarball model."""
    from Bio.SeqUtils.CheckSum import crc64

    db_path = _db_path(cfg)
    if not db_path.exists():
        raise FileNotFoundError(f"{db_path} not found — run `python scripts/local_store.py "
                                "--build` first.")
    conn = sqlite3.connect(db_path)
    if conn.execute("SELECT name FROM sqlite_master WHERE type='table' "
                    "AND name='sequence_aliases'").fetchone():
        if not force:
            print(f"[skip] sequence_aliases already present in {db_path}; pass force=True")
            conn.close()
            return db_path
        conn.execute("DROP TABLE sequence_aliases")
    conn.execute("""
        CREATE TABLE sequence_aliases (
            crc64     TEXT,
            accession TEXT,
            sequence  TEXT,
            source    TEXT,
            has_model INTEGER,
            PRIMARY KEY (crc64, accession)
        )
    """)
    conn.executemany("INSERT OR IGNORE INTO sequence_aliases VALUES "
                     "(?, ?, ?, 'reference_proteome', 0)",
                     conn.execute("SELECT crc64, accession, sequence FROM proteins").fetchall())

    pdb_dir = afdb_model_path("X", cfg).parent
    rows = []
    # every model file, including ones structure_bulk fetched and cached after extraction
    for i, path in enumerate(sorted(pdb_dir.glob("AF-*-F1-model_v6.pdb.gz"))):
        acc = path.name[3:-len("-F1-model_v6.pdb.gz")]
        try:
            seq = _pdb_sequence(path)
        except (OSError, EOFError) as e:
            print(f"[local_store] unreadable model {path.name}: {e}")
            continue
        rows.append((crc64(seq).replace("CRC-", ""), acc, seq))
        if (i + 1) % 5000 == 0:
            print(f"[local_store] read {i + 1:,} models …")
    # has_model = the accession's model carries exactly this sequence; a model built on an
    # older sequence version lands under its own checksum, leaving the proteome row at 0
    conn.executemany("INSERT INTO sequence_aliases VALUES (?, ?, ?, 'afdb_model', 1) "
                     "ON CONFLICT (crc64, accession) DO UPDATE SET has_model = 1", rows)
    conn.commit()
    n, n_crc = conn.execute("SELECT COUNT(*), COUNT(DISTINCT crc64) FROM sequence_aliases").fetchone()
    shared = conn.execute("SELECT COUNT(*) FROM (SELECT crc64 FROM sequence_aliases GROUP BY crc64 "
                          "HAVING COUNT(*) > 1)").fetchone()[0]
    print(f"[local_store] {n:,} accessions over {n_crc:,} sequences ({shared:,} sequences "
          f"under more than one accession) -> {db_path}")
    conn.close()
    return db_path


def aliases_for_sequence(seq, conn=None):
    """Accessions whose indexed sequence equals seq, proteome entries first; [] without the table."""
    return [acc for acc, _, _ in alias_rows(seq, conn)]


def alias_rows(seq, conn=None):
    """(accession, source, has_model) whose indexed sequence equals seq, proteome entries first.

    has_model = AFDB has a model of exactly this sequence under that accession.

    Selenocysteine is compared as cysteine, since AFDB models carry it as CYS.
    """
    from Bio.SeqUtils.CheckSum import crc64

    seq = "".join(str(seq).split()).upper().rstrip("*")
    if not seq:
        return []
    own_conn = conn is None
    conn = conn or open_index()
    try:
        if not conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                            "AND name='sequence_aliases'").fetchone():
            return []
        # a model's sequence has C where the proteome's has U, so both checksums are looked up
        keys = {crc64(seq).replace("CRC-", ""), crc64(seq.replace("U", "C")).replace("CRC-", "")}
        rows = conn.execute(
            "SELECT accession, sequence, source, has_model FROM sequence_aliases WHERE crc64 IN ({}) "
            "ORDER BY source != 'reference_proteome', accession".format(",".join("?" * len(keys))),
            list(keys)).fetchall()
    finally:
        if own_conn:
            conn.close()
    # the checksum only narrows; the sequence itself decides
    want = seq.replace("U", "C")
    out, seen = [], set()
    for acc, s, source, has_model in rows:
        if s.replace("U", "C") == want and acc not in seen:
            seen.add(acc)
            out.append((acc, source, has_model))
    return out


def lookup_domains_by_accession(accession, conn=None):
    """Return all InterPro match rows for `accession` as a list of dicts."""
    own_conn = conn is None
    conn = conn or open_index()
    try:
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT * FROM domains WHERE accession = ?", (accession,)
        ).fetchall()
        return [dict(r) for r in rows]
    finally:
        if own_conn:
            conn.close()


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
    parser.add_argument("--build", action="store_true", help="(re)build the proteins index")
    parser.add_argument("--build-domains", action="store_true", help="(re)build the domains index")
    parser.add_argument("--build-aliases", action="store_true",
                        help="(re)build the sequence_aliases index (proteome + AFDB models)")
    parser.add_argument("--force", action="store_true", help="rebuild even if the index already exists")
    parser.add_argument("--stats", action="store_true", help="print row counts for the existing index")
    args = parser.parse_args()

    if args.build:
        build_index(force=args.force)
    if args.build_domains:
        build_domains_index(force=args.force)
    if args.build_aliases:
        build_aliases_index(force=args.force)
    if args.stats:
        conn = open_index()
        n = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        n_wb = conn.execute("SELECT COUNT(*) FROM proteins WHERE wormbase_transcript IS NOT NULL").fetchone()[0]
        print(f"proteins: {n:,}  (with WormBase transcript name: {n_wb:,})")
        try:
            n_dom = conn.execute("SELECT COUNT(*) FROM domains").fetchone()[0]
            print(f"domain matches: {n_dom:,}")
        except sqlite3.OperationalError:
            pass
        try:
            n_al = conn.execute("SELECT COUNT(*) FROM sequence_aliases").fetchone()[0]
            print(f"sequence aliases: {n_al:,}")
        except sqlite3.OperationalError:
            pass
    if not (args.build or args.build_domains or args.build_aliases or args.stats):
        parser.error("specify --build, --build-domains, --build-aliases, and/or --stats")
