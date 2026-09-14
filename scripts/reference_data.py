"""
reference_data.py

Phase A of the proteome-scale batch pipeline: one-time bulk downloads that
replace the per-protein EBI/Ensembl calls the interactive pipeline makes.
See ~/.claude/plans/i-m-considering-a-large-humble-sun.md.

Sources (verified reachable at write time; WormBase's own downloads portal is
Cloudflare-gated against non-browser clients, so genome/annotation comes from
Ensembl's plain FTP mirror instead — same WBcel235 assembly and gene models):

  uniprot   UniProt proteome UP000001940 (canonical + isoforms), FASTA + JSON
  afdb      AlphaFold DB proteome tarball for UP000001940 (~2.8 GB)
  genome    Ensembl C. elegans genome FASTA (soft-masked, toplevel)
  gff3      Ensembl C. elegans GFF3 annotation (CDS exons, gene/transcript spans)
  interpro  protein2ipr.dat.gz (~13 GB compressed, all of UniProt), streamed
            and filtered down to the accessions in the UniProt proteome FASTA
            without ever writing the full decompressed file to disk

All downloads are resumable (curl -C -) and skipped if the target file already
exists with a non-zero size — re-running this script after an interrupted
download picks up where it left off; use --force to redo a step anyway.

Paths and URLs are read from batch.config.json (repo root) under a
"reference_data" key when present; see DEFAULTS below for the fallback used
when that file/key is absent.

Usage
-----
    python scripts/reference_data.py --all
    python scripts/reference_data.py --only uniprot,genome,gff3
    python scripts/reference_data.py --only interpro --force
"""

import gzip
import json
import subprocess
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent
_CONFIG_PATH = _REPO_ROOT / "batch.config.json"

DEFAULTS = {
    "out_dir": str(_REPO_ROOT / "data" / "reference"),
    "uniprot_proteome_id": "UP000001940",
    "afdb_tarball_url": (
        "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"
        "UP000001940_6239_CAEEL_v6.tar"
    ),
    "ensembl_release": "114",
    "ensembl_species": "caenorhabditis_elegans",
    "ensembl_assembly": "WBcel235",
    "interpro_release": "110.0",
}


def _load_config():
    """Merge batch.config.json's "reference_data" block (if any) over DEFAULTS."""
    cfg = dict(DEFAULTS)
    if _CONFIG_PATH.exists():
        with open(_CONFIG_PATH) as f:
            user_cfg = json.load(f).get("reference_data", {})
        cfg.update(user_cfg)
    return cfg


def _out_dir(cfg):
    d = Path(cfg["out_dir"])
    d.mkdir(parents=True, exist_ok=True)
    return d


def _curl_download(url, dest, force=False, extra_args=None, resumable=True):
    """Download via curl; skip if dest already has non-zero size and force isn't
    set. Returns True if a download ran, False if skipped.

    resumable=True adds -C - (byte-range resume) for static files that support
    HTTP Range (AFDB/Ensembl FTP). Set resumable=False for compute-on-the-fly
    endpoints (e.g. UniProt's /stream) that don't reliably support Range —
    resuming against one of those risks silently splicing two different result
    orderings together instead of actually continuing the same download.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {dest} already present ({dest.stat().st_size:,} bytes)")
        return False
    print(f"[fetch] {url} -> {dest}")
    cmd = ["curl", "-fL", "--retry", "5", "--retry-delay", "10"]
    if resumable:
        cmd += ["-C", "-"]
    cmd += ["-o", str(dest), url]
    if extra_args:
        cmd[1:1] = extra_args
    subprocess.run(cmd, check=True)
    return True


# ── UniProt proteome (sequences, CRC64s, curated features) ───────────────────

def fetch_uniprot(cfg, force=False):
    """Download the proteome's canonical+isoform FASTA and per-entry JSON stream.

    The FASTA gives sequences/CRC64s for local_store's checksum-keyed cache;
    the JSON stream gives curated features (uniprot_features.py's network
    calls) and WormBase/Ensembl cross-references (the accession<->GFF-name map
    genome_regions.py needs — see the plan's Phase A2).
    """
    out = _out_dir(cfg)
    pid = cfg["uniprot_proteome_id"]

    fasta_url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=proteome:{pid}&format=fasta&includeIsoform=true"
    )
    _curl_download(fasta_url, out / f"{pid}.fasta", force=force, resumable=False)

    # JSON stream: one array of full entry records (features, xrefs, CRC64 checksum)
    json_url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=proteome:{pid}&format=json"
    )
    _curl_download(json_url, out / f"{pid}.json", force=force, resumable=False)


# ── AlphaFold DB structures ───────────────────────────────────────────────────

def fetch_afdb(cfg, force=False):
    """Download the AFDB proteome tarball (one PDB per canonical UniProt entry)."""
    out = _out_dir(cfg)
    url = cfg["afdb_tarball_url"]
    dest = out / Path(url).name
    _curl_download(url, dest, force=force)


# ── Ensembl genome + annotation ───────────────────────────────────────────────

def fetch_genome(cfg, force=False):
    """Download the soft-masked toplevel genome FASTA (for bedtools/pyfaidx region extraction)."""
    out = _out_dir(cfg)
    release, species, assembly = cfg["ensembl_release"], cfg["ensembl_species"], cfg["ensembl_assembly"]
    fname = f"{species.capitalize()}.{assembly}.dna_sm.toplevel.fa.gz"
    url = f"https://ftp.ensembl.org/pub/release-{release}/fasta/{species}/dna/{fname}"
    _curl_download(url, out / fname, force=force)


def fetch_gff3(cfg, force=False):
    """Download the GFF3 annotation (CDS exons per transcript — the genewise_bulk.py input)."""
    out = _out_dir(cfg)
    release, species, assembly = cfg["ensembl_release"], cfg["ensembl_species"], cfg["ensembl_assembly"]
    fname = f"{species.capitalize()}.{assembly}.{release}.gff3.gz"
    url = f"https://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{fname}"
    _curl_download(url, out / fname, force=force)


# ── InterPro domain matches (filtered to this proteome's accessions) ─────────

def _load_proteome_accessions(cfg):
    """Read accessions out of the UniProt proteome FASTA already downloaded by fetch_uniprot()."""
    out = _out_dir(cfg)
    fasta_path = out / f"{cfg['uniprot_proteome_id']}.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"{fasta_path} not found — run fetch_uniprot() (or --only uniprot) first; "
            "the InterPro filter needs the proteome's accession list."
        )
    accessions = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # UniProt FASTA header: >sp|ACCESSION|... or >tr|ACCESSION|...
                parts = line[1:].split("|")
                if len(parts) >= 2:
                    accessions.add(parts[1])
    return accessions


def fetch_interpro(cfg, force=False):
    """Stream protein2ipr.dat.gz (all of UniProt, ~13 GB compressed) through a
    gunzip pipe and keep only rows for this proteome's accessions, so the full
    decompressed file (tens of GB) is never written to disk. Output is a plain
    TSV with the same protein2ipr.dat columns: accession, ipr_id, description,
    external_db_match_id, start, stop.
    """
    out = _out_dir(cfg)
    dest = out / "protein2ipr.filtered.tsv"
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {dest} already present ({dest.stat().st_size:,} bytes)")
        return

    accessions = _load_proteome_accessions(cfg)
    print(f"[interpro] filtering to {len(accessions):,} proteome accessions")

    url = (
        "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/"
        f"{cfg['interpro_release']}/protein2ipr.dat.gz"
    )
    print(f"[fetch+filter] {url} -> {dest}")

    # curl streams compressed bytes to stdout; gunzip decompresses on the fly;
    # the accession set is checked in this process rather than shelling out to
    # grep/awk, since a plain substring/set membership test is both simpler and
    # correct here (protein2ipr rows are already one-accession-per-line, so
    # there's no risk of a grep pattern matching the wrong column).
    curl = subprocess.Popen(
        ["curl", "-fsL", url],
        stdout=subprocess.PIPE,
    )
    matched = 0
    total = 0
    with gzip.GzipFile(fileobj=curl.stdout) as gz, open(dest, "w") as out_f:
        for raw_line in gz:
            total += 1
            line = raw_line.decode("utf-8", errors="replace")
            acc = line.split("\t", 1)[0]
            if acc in accessions:
                out_f.write(line)
                matched += 1
            if total % 5_000_000 == 0:
                print(f"[interpro] scanned {total:,} rows, matched {matched:,}…")
    ret = curl.wait()
    if ret != 0:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"curl exited {ret} while streaming {url}")
    print(f"[interpro] done: {matched:,} matching rows out of {total:,} scanned -> {dest}")


STEPS = {
    "uniprot":  fetch_uniprot,
    "afdb":     fetch_afdb,
    "genome":   fetch_genome,
    "gff3":     fetch_gff3,
    "interpro": fetch_interpro,  # depends on "uniprot" having run first
}


def main(steps, force=False):
    cfg = _load_config()
    _out_dir(cfg)
    for step in steps:
        if step not in STEPS:
            raise ValueError(f"unknown step {step!r}; choose from {sorted(STEPS)}")
    for step in steps:
        print(f"=== {step} ===")
        STEPS[step](cfg, force=force)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true",
                        help="run every step (uniprot, afdb, genome, gff3, interpro)")
    parser.add_argument("--only", type=str, default=None,
                        help="comma-separated subset of steps to run, e.g. uniprot,genome,gff3")
    parser.add_argument("--force", action="store_true",
                        help="redo a step even if its output already exists")
    args = parser.parse_args()

    if args.all:
        selected = list(STEPS.keys())
    elif args.only:
        selected = [s.strip() for s in args.only.split(",") if s.strip()]
    else:
        parser.error("specify --all or --only <steps>")

    main(selected, force=args.force)
