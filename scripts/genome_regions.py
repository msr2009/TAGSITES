"""
genome_regions.py

Phase A2 of the proteome-scale batch pipeline: parses the WormBase GFF3
annotation (downloaded via scripts/reference_data.py) exactly once into a
SQLite index, and extracts genomic sequence directly from the local genome
FASTA — replacing both run_genewise.py's per-protein EBI Genewise jobs (see
scripts/genewise_remote.py) and design_guides_across_region.py's per-gene
`grep` over the whole GFF3 file (O(genes x filesize); at proteome scale that
dominates wall clock on its own).

No bedtools dependency: sequence extraction uses Biopython (already a
project dependency) instead of shelling out to `bedtools getfasta`.

Key output: get_cds_dataframe(transcript_id) returns the same column shape
as scripts/parse_genewise.py's parse_genewise() — {name, source, type,
start (0-indexed), stop (0-indexed), score, strand, frame, note} — so
scripts/genewise_bulk.py (not yet implemented) can hand its result straight
to enumerate_insertion_sites() exactly as the remote Genewise backend does.

Matching a UniProt accession to a specific GFF3 transcript is NOT a plain
string match: UniProt's cross-reference (see local_store.py's
wormbase_transcript) is the gene-level name (e.g. "C10C5.1"), while GFF3 CDS
rows carry the isoform-lettered protein_id (e.g. "C10C5.1g.1") — the isoform
letter has to come from the UniProt entry's ALTERNATIVE PRODUCTS block. That
resolution isn't implemented here yet; get_cds_dataframe() takes the GFF3
transcript_id (protein_id) directly.

Usage
-----
    python scripts/genome_regions.py --build
    python scripts/genome_regions.py --stats
"""

import gzip
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

sys.path.insert(0, str(Path(__file__).parent))
from providers import _load_config as _load_batch_config

_REPO_ROOT = Path(__file__).parent.parent

_ATTR_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^;]*)")


def _parse_attributes(attr_field):
    """Parse a GFF3 column-9 attribute string into a dict; values are left as-is
    (not URL-decoded — none of the fields this module reads need it)."""
    return dict(_ATTR_RE.findall(attr_field))


def _reference_dir(cfg=None):
    cfg = cfg or _load_batch_config().get("reference_data", {})
    out_dir = Path(cfg.get("out_dir", "data/reference"))
    if not out_dir.is_absolute():
        out_dir = _REPO_ROOT / out_dir
    return out_dir, cfg


def _gff3_path(cfg=None):
    ref_dir, cfg = _reference_dir(cfg)
    species = cfg.get("ensembl_species", "caenorhabditis_elegans")
    assembly = cfg.get("ensembl_assembly", "WBcel235")
    release = cfg.get("ensembl_release", "114")
    return ref_dir / f"{species.capitalize()}.{assembly}.{release}.gff3.gz"


def _genome_fasta_path(cfg=None):
    ref_dir, cfg = _reference_dir(cfg)
    species = cfg.get("ensembl_species", "caenorhabditis_elegans")
    assembly = cfg.get("ensembl_assembly", "WBcel235")
    return ref_dir / f"{species.capitalize()}.{assembly}.dna_sm.toplevel.fa.gz"


def _db_path(cfg=None):
    ref_dir, _ = _reference_dir(cfg)
    return ref_dir / "genome_regions.sqlite3"


# ── GFF3 -> SQLite (parsed once, not re-scanned per gene) ────────────────────

def build_index(cfg=None, force=False):
    """Parse the GFF3 exactly once into two tables:

      cds_exons(transcript_id, chrom, start, stop, strand, phase, protein_id)
        — one row per CDS row in the file, 1-based inclusive coordinates as
          given (converted to 0-based only in get_cds_dataframe(), to match
          parse_genewise.parse_genewise()'s convention exactly).
      transcript_spans(transcript_id, chrom, start, stop, strand)
        — the full mRNA/transcript span (min/max over its CDS rows would
          exclude UTRs; this uses the GFF3 "mRNA" feature row directly).

    Indexed on transcript_id (both tables) so a lookup is O(log n), not the
    O(genes x filesize) cost of grep-ing the whole file per gene.
    """
    gff_path = _gff3_path(cfg)
    if not gff_path.exists():
        raise FileNotFoundError(
            f"{gff_path} not found — run `python scripts/reference_data.py "
            "--only gff3` first."
        )

    db_path = _db_path(cfg)
    if db_path.exists():
        if not force:
            print(f"[skip] {db_path} already present; pass force=True to rebuild")
            return db_path
        db_path.unlink()

    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE cds_exons (
            transcript_id TEXT,
            chrom         TEXT,
            start         INTEGER,
            stop          INTEGER,
            strand        TEXT,
            phase         INTEGER,
            row_order     INTEGER
        )
    """)
    conn.execute("CREATE INDEX idx_cds_transcript ON cds_exons(transcript_id)")
    conn.execute("""
        CREATE TABLE transcript_spans (
            transcript_id TEXT PRIMARY KEY,
            chrom         TEXT,
            start         INTEGER,
            stop          INTEGER,
            strand        TEXT
        )
    """)

    print(f"[genome_regions] parsing {gff_path} …")
    cds_rows = []
    span_rows = []
    row_order = 0
    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, source, ftype, start, stop, score, strand, phase, attrs = parts
            if ftype == "CDS":
                attr = _parse_attributes(attrs)
                transcript_id = attr.get("protein_id") or attr.get("Parent", "").split(":")[-1]
                if not transcript_id:
                    continue
                row_order += 1
                cds_rows.append((
                    transcript_id, chrom, int(start), int(stop), strand,
                    0 if phase in (".", "") else int(phase), row_order,
                ))
            elif ftype in ("mRNA", "transcript"):
                attr = _parse_attributes(attrs)
                transcript_id = attr.get("transcript_id") or attr.get("ID", "").split(":")[-1]
                if transcript_id:
                    span_rows.append((transcript_id, chrom, int(start), int(stop), strand))

    conn.executemany(
        "INSERT INTO cds_exons VALUES (?, ?, ?, ?, ?, ?, ?)", cds_rows
    )
    conn.executemany(
        "INSERT OR REPLACE INTO transcript_spans VALUES (?, ?, ?, ?, ?)", span_rows
    )
    conn.commit()
    n_transcripts = len({r[0] for r in cds_rows})
    print(f"[genome_regions] indexed {len(cds_rows):,} CDS rows across "
          f"{n_transcripts:,} transcripts, {len(span_rows):,} transcript spans -> {db_path}")
    conn.close()
    return db_path


def open_index(cfg=None):
    db_path = _db_path(cfg)
    if not db_path.exists():
        raise FileNotFoundError(
            f"{db_path} not found — run `python scripts/genome_regions.py --build` first."
        )
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)


def get_cds_dataframe(transcript_id, conn=None):
    """Return a DataFrame in exactly parse_genewise.parse_genewise()'s shape
    for `transcript_id`'s CDS rows: columns name/source/type/start(0-indexed)/
    stop(0-indexed)/score/strand/frame/note, sorted by start.

    Raises ValueError if the transcript isn't in the index (mirrors
    parse_genewise() raising when a Genewise .out.txt has no CDS rows).
    """
    own_conn = conn is None
    conn = conn or open_index()
    try:
        rows = conn.execute(
            "SELECT chrom, start, stop, strand, phase FROM cds_exons "
            "WHERE transcript_id = ? ORDER BY row_order",
            (transcript_id,),
        ).fetchall()
    finally:
        if own_conn:
            conn.close()

    if not rows:
        raise ValueError(f"transcript_id {transcript_id!r} not found in genome_regions index")

    df = pd.DataFrame(rows, columns=["name", "start", "stop", "strand", "frame"])
    df["source"] = "WormBase"
    df["type"] = "cds"
    df["score"] = 0.0
    df["note"] = ""
    # 1-based inclusive (GFF3) -> 0-based, matching parse_genewise.parse_genewise()
    df["start"] = df["start"] - 1
    df["stop"] = df["stop"] - 1
    df = df[["name", "source", "type", "start", "stop", "score", "strand", "frame", "note"]]
    return df.sort_values("start").reset_index(drop=True)


def get_transcript_span(transcript_id, conn=None):
    """Return (chrom, start, stop, strand) for a transcript's full mRNA span
    (1-based inclusive, as in the GFF3), or None if not indexed.
    """
    own_conn = conn is None
    conn = conn or open_index()
    try:
        row = conn.execute(
            "SELECT chrom, start, stop, strand FROM transcript_spans WHERE transcript_id = ?",
            (transcript_id,),
        ).fetchone()
    finally:
        if own_conn:
            conn.close()
    return row


# ── Genome FASTA -> in-memory sequence lookup ─────────────────────────────────

_genome_cache = None  # lazy singleton: {chrom: sequence str}


def _load_genome(cfg=None):
    global _genome_cache
    if _genome_cache is not None:
        return _genome_cache

    fasta_path = _genome_fasta_path(cfg)
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"{fasta_path} not found — run `python scripts/reference_data.py "
            "--only genome` first."
        )

    print(f"[genome_regions] loading genome FASTA {fasta_path} into memory …")
    genome = {}
    with gzip.open(fasta_path, "rt") as f:
        chrom = None
        chunks = []
        for line in f:
            if line.startswith(">"):
                if chrom is not None:
                    genome[chrom] = "".join(chunks)
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if chrom is not None:
            genome[chrom] = "".join(chunks)
    _genome_cache = genome
    return genome


def extract_sequence(chrom, start, stop, strand="+", cfg=None):
    """Extract genomic sequence [start, stop] (1-based inclusive, GFF3-style)
    from the local genome FASTA, reverse-complemented if strand == "-" so the
    returned sequence is always in the coding (5'->3') orientation — this is
    what removes the need for run_genewise.py's separate reverse-complement
    submission once a gene's strand is known from the annotation.
    """
    genome = _load_genome(cfg)
    if chrom not in genome:
        raise KeyError(f"chromosome/scaffold {chrom!r} not found in genome FASTA")
    seq = genome[chrom][start - 1:stop]
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


def get_transcript_region(transcript_id, conn=None, cfg=None):
    """One-call convenience for a bulk Genewise-replacement backend: return
    {"dna", "cds_df", "chrom", "start", "stop", "strand"} for `transcript_id`,
    where "dna" is the transcript's full genomic span extracted and oriented
    to the coding strand, and "cds_df" is get_cds_dataframe()'s exon table
    with start/stop shifted to be local offsets into "dna" — i.e. exactly the
    (cds_df, dna) pair enumerate_insertion_sites() expects, with no separate
    reverse-complement submission needed since strand is already known here.
    """
    own_conn = conn is None
    conn = conn or open_index()
    try:
        df = get_cds_dataframe(transcript_id, conn=conn)
        span = get_transcript_span(transcript_id, conn=conn)
    finally:
        if own_conn:
            conn.close()

    if span is None:
        raise ValueError(f"transcript_id {transcript_id!r} has no indexed transcript span")
    chrom, start, stop, strand = span

    dna = extract_sequence(chrom, start, stop, strand, cfg=cfg)

    local_df = df.copy()
    local_df["start"] -= (start - 1)
    local_df["stop"] -= (start - 1)
    if strand == "-":
        # exon coordinates were chromosome-absolute in the + orientation;
        # extract_sequence already reverse-complemented "dna", so exon offsets
        # must be flipped into that same reversed frame
        span_len = stop - start + 1
        new_start = span_len - 1 - local_df["stop"]
        new_stop = span_len - 1 - local_df["start"]
        local_df["start"], local_df["stop"] = new_start, new_stop
        local_df = local_df.sort_values("start").reset_index(drop=True)

    return {
        "dna": dna, "cds_df": local_df,
        "chrom": chrom, "start": start, "stop": stop, "strand": strand,
    }


def _translate_cds(cds_df, dna):
    """Translate a get_transcript_region()-style local CDS table against its
    DNA into a protein string (no trailing stop codon), or None if the CDS
    length isn't a multiple of 3 or contains an untranslatable codon.
    """
    cds_seq = "".join(dna[row.start:row.stop + 1] for row in cds_df.itertuples())
    if len(cds_seq) == 0 or len(cds_seq) % 3 != 0:
        return None
    try:
        protein = str(Seq(cds_seq).translate())
    except TranslationError:
        return None
    return protein[:-1] if protein.endswith("*") else protein


def resolve_transcript_for_accession(accession, wormbase_gene, expected_sequence, conn=None, cfg=None):
    """Find the GFF3 transcript_id whose CDS translates to `expected_sequence`
    (a UniProt protein sequence), starting from candidates whose transcript_id
    has `wormbase_gene` as a locus prefix (e.g. wormbase_gene "C10C5.1g" ->
    candidates "C10C5.1g.1", "C10C5.1g.2", ...).

    Matching by prefix alone isn't sufficient: WormBase sometimes has several
    numbered transcript versions under the same isoform-lettered locus name
    (splice/UTR variants) with no further hint of which one the FASTA/JSON
    sequence in local_store.py corresponds to — verified in a 500-accession
    sample, ~8% of prefix matches were ambiguous (multiple candidates) and
    ~1% had none. Translating each candidate's CDS and comparing to the known
    protein sequence resolves the ambiguity outright, and also catches the
    rare true mismatch a prefix match alone can't detect.

    Returns the resolved transcript_id, or None if no candidate's translation
    matches (the caller should fall back to genewise_remote.py in that case).
    """
    own_conn = conn is None
    conn = conn or open_index()
    try:
        candidates = [
            r[0] for r in conn.execute(
                "SELECT DISTINCT transcript_id FROM transcript_spans WHERE transcript_id LIKE ?",
                (wormbase_gene + ".%",),
            ).fetchall()
        ]
        for transcript_id in candidates:
            try:
                region = get_transcript_region(transcript_id, conn=conn, cfg=cfg)
            except (ValueError, KeyError):
                continue
            protein = _translate_cds(region["cds_df"], region["dna"])
            if protein == expected_sequence:
                return transcript_id
        return None
    finally:
        if own_conn:
            conn.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="(re)build the GFF3 index")
    parser.add_argument("--force", action="store_true", help="rebuild even if the index already exists")
    parser.add_argument("--stats", action="store_true", help="print row counts for the existing index")
    args = parser.parse_args()

    if args.build:
        build_index(force=args.force)
    if args.stats:
        conn = open_index()
        n_cds = conn.execute("SELECT COUNT(*) FROM cds_exons").fetchone()[0]
        n_transcripts = conn.execute("SELECT COUNT(DISTINCT transcript_id) FROM cds_exons").fetchone()[0]
        n_spans = conn.execute("SELECT COUNT(*) FROM transcript_spans").fetchone()[0]
        print(f"CDS rows: {n_cds:,}  transcripts: {n_transcripts:,}  transcript spans: {n_spans:,}")
    if not args.build and not args.stats:
        parser.error("specify --build and/or --stats")
