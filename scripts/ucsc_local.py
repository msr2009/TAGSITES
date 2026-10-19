"""ucsc_local.py — offline copies of the UCSC data the gene-model and phastCons steps use.

Bulk mode (batch.config.json backends.genewise == "bulk") must not make per-protein REST
calls, so `python scripts/reference_data.py --only ucsc` downloads, for the assembly named by
reference_data.ucsc_genome (default ce11):

  {genome}.unipAliSwissprot.bb / .unipAliTrembl.bb   UniProt-to-genome alignments (bigPsl)
  {genome}.{phastcons track}.bw                       per-base phastCons scores (bigWig)

and build_index() dumps both bigBeds into {genome}.unipali.sqlite3, keyed by accession and by
the aligned protein's CRC64 (pyBigWig cannot search a bigBed by name). Rows have the same
field names as the REST getData/track items, so ucsc_gene_model.py treats both alike.

UCSC chromosome names are mapped onto the local genome's (Ensembl) names by
reference_data.ucsc_chrom_alias (explicit pairs, default chrM -> MtDNA) and otherwise by
dropping a leading "chr". build_index() checks that by translating a sample of Swiss-Prot
alignments from the local genome: an assembly or chromosome-name mismatch fails the build
rather than every later lookup. A sampled alignment counts as consistent when its translation
is within cds_check's limits of the aligned protein or a part of it: UniProt alignments can
be partial (paralog hits) or a codon short at the ends, which says nothing about the mapping.

Matt Rich, 2026
"""

import json
import os
import random
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

TRACKS = ("unipAliSwissprot", "unipAliTrembl")
DEFAULT_GENOME = "ce11"
DEFAULT_PHASTCONS = {"ce11": "phastCons26way"}
DEFAULT_ALIAS = {"chrM": "MtDNA"}
VALIDATE_SAMPLE = 50  # Swiss-Prot alignments translated from the local genome at build time
VALIDATE_MIN = 0.8  # fraction that must translate consistently with their protein


def settings():
    """(genome, phastcons_track, alias dict, folder) from batch.config.json's reference_data."""
    import reference_data

    cfg = reference_data._load_config()
    genome = cfg.get("ucsc_genome", DEFAULT_GENOME)
    track = cfg.get("ucsc_phastcons_track", DEFAULT_PHASTCONS.get(genome, ""))
    alias = dict(DEFAULT_ALIAS, **(cfg.get("ucsc_chrom_alias") or {}))
    return genome, track, alias, Path(cfg["out_dir"]) / "ucsc"


def paths():
    """Local file paths for the configured assembly (they may not exist yet)."""
    genome, track, _, folder = settings()
    return {
        "unipAliSwissprot": folder / f"{genome}.unipAliSwissprot.bb",
        "unipAliTrembl": folder / f"{genome}.unipAliTrembl.bb",
        "phastcons": folder / f"{genome}.{track}.bw",
        "index": folder / f"{genome}.unipali.sqlite3",
    }


def available():
    """True when the accession index has been built."""
    return paths()["index"].exists()


def local_chrom(ucsc_chrom):
    """UCSC chromosome name -> the local genome FASTA's name."""
    _, _, alias, _ = settings()
    if ucsc_chrom in alias:
        return alias[ucsc_chrom]
    return ucsc_chrom.removeprefix("chr")


def _bigbed_fields(bb):
    """Column names after chrom/chromStart/chromEnd, from the bigBed's autoSql."""
    names = []
    for line in bb.SQL().decode().splitlines():
        line = line.strip()
        # field lines look like:  type name;  "comment"
        if ";" in line and not line.startswith(("table", "(", ")", '"')):
            names.append(line.split(";")[0].split()[-1])
    return names[3:]


def build_index(force=False):
    """Dump the two UniProt bigBeds into the SQLite accession index, then validate it."""
    import pyBigWig
    from Bio.SeqUtils.CheckSum import crc64

    p = paths()
    if p["index"].exists() and not force:
        print(f"[skip] {p['index']} already present")
        return p["index"]
    tmp = p["index"].with_suffix(".tmp")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    conn.execute("""CREATE TABLE items (name TEXT, acc TEXT, is_main INTEGER, track TEXT,
                    status TEXT, chrom TEXT, chromStart INTEGER, chromEnd INTEGER, strand TEXT,
                    blockSizes TEXT, chromStarts TEXT, crc64 TEXT, json TEXT)""")
    n = 0
    for track in TRACKS:
        bb = pyBigWig.open(str(p[track]))
        fields = _bigbed_fields(bb)
        # every chromosome, every item
        for chrom, size in bb.chroms().items():
            for start, end, rest in bb.entries(chrom, 0, size) or []:
                row = dict(zip(fields, rest.split("\t")))
                seq = row.get("oSequence", "")
                item = {"chrom": chrom, "chromStart": start, "chromEnd": end,
                        "name": row.get("name", ""), "strand": row.get("strand", "+"),
                        "blockSizes": row.get("blockSizes", ""),
                        "chromStarts": row.get("chromStarts", ""),
                        "acc": row.get("acc", ""), "isMain": row.get("isMain", ""),
                        "status": row.get("status", ""), "geneName": row.get("geneName", "")}
                conn.execute("INSERT INTO items VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    item["name"], item["acc"], int(bool(item["isMain"])), track, item["status"],
                    chrom, start, end, item["strand"], item["blockSizes"], item["chromStarts"],
                    crc64(seq).replace("CRC-", "") if seq else "", json.dumps(item)))
                n += 1
        bb.close()
    conn.execute("CREATE INDEX i_name ON items(name)")
    conn.execute("CREATE INDEX i_acc ON items(acc)")
    conn.execute("CREATE INDEX i_crc ON items(crc64)")
    conn.commit()
    conn.close()
    _validate(tmp)
    os.replace(tmp, p["index"])
    print(f"[ucsc] indexed {n:,} UniProt alignments -> {p['index']}")
    return p["index"]


def _validate(index_path):
    """Translate sampled Swiss-Prot alignments from the local genome; raise on a mismatch."""
    import pyBigWig
    from cds_check import compare_to_query, translate_model
    from ucsc_gene_model import item_to_model

    p = paths()
    conn = sqlite3.connect(index_path)
    rows = conn.execute("SELECT json FROM items WHERE track='unipAliSwissprot'").fetchall()
    conn.close()
    random.seed(0)
    sample = [json.loads(r[0]) for r in random.sample(rows, min(VALIDATE_SAMPLE, len(rows)))]
    bb = pyBigWig.open(str(p["unipAliSwissprot"]))
    fields = _bigbed_fields(bb)
    ok = 0
    for item in sample:
        # the aligned protein itself, to compare against
        entries = bb.entries(item["chrom"], item["chromStart"], item["chromEnd"]) or []
        seqs = [dict(zip(fields, rest.split("\t"))) for s, e, rest in entries
                if s == item["chromStart"] and e == item["chromEnd"]]
        prot = next((r.get("oSequence", "") for r in seqs if r.get("name") == item["name"]), "")
        try:
            model = item_to_model(item, flank_bp=0, backend="local")
        except Exception:
            continue
        tr = translate_model(model["cds_df"], model["dna"])
        # partial alignments translate to a piece of the protein; others to (nearly) all of it
        if prot and (tr in prot or compare_to_query(tr, prot)["status"] != "fail"):
            ok += 1
    bb.close()
    frac = ok / max(len(sample), 1)
    print(f"[ucsc] validation: {ok}/{len(sample)} sampled alignments translate consistently "
          "from the local genome")
    if frac < VALIDATE_MIN:
        raise RuntimeError("UniProt alignments do not translate from the local genome "
                           f"({frac:.0%}); check reference_data.ucsc_genome / ucsc_chrom_alias "
                           "against the Ensembl genome")


def lookup(accession="", crc=""):
    """Index rows (REST-shaped item dicts plus '_track') for an accession and/or sequence CRC64.

    A bare accession also returns its isoforms, matching ucsc_rest.search_accession.
    """
    conn = sqlite3.connect(f"file:{paths()['index']}?mode=ro", uri=True)
    try:
        rows = []
        if accession:
            base = accession.split("-")[0]
            rows += conn.execute("SELECT json, track FROM items WHERE name = ? OR name LIKE ? "
                                 "OR acc = ?", (accession, base + "-%", base)).fetchall()
        if crc:
            rows += conn.execute("SELECT json, track FROM items WHERE crc64 = ?",
                                 (crc,)).fetchall()
    finally:
        conn.close()
    out, seen = [], set()
    for js, track in rows:
        item = dict(json.loads(js), _track=track)
        key = (item["name"], item["chrom"], item["chromStart"], track)
        if key not in seen:
            seen.add(key)
            out.append(item)
    return out


def phastcons_values(chrom, start0, end):
    """{0-based position: score} from the local bigWig over [start0, end); unscored omitted."""
    import math

    import pyBigWig

    bw = pyBigWig.open(str(paths()["phastcons"]))
    try:
        vals = bw.values(chrom, start0, end)
    finally:
        bw.close()
    return {start0 + i: v for i, v in enumerate(vals) if not math.isnan(v)}


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Build or query the local UCSC UniProt index.")
    parser.add_argument("--build", action="store_true", help="build the accession index")
    parser.add_argument("--force", action="store_true", help="rebuild even if present")
    parser.add_argument("--lookup", default="", help="print index rows for an accession")
    args = parser.parse_args()
    if args.build:
        build_index(force=args.force)
    if args.lookup:
        for row in lookup(args.lookup):
            print(row)
