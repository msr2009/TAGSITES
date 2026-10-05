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
  swissprot        all of reviewed UniProt (Swiss-Prot, ~575k sequences),
                   indexed with `diamond makedb` — one of two local search
                   targets for scripts/conservation_local.py (see its own
                   docstring and batch.config.json's conservation_local key)
  rhabditida_trembl  all unreviewed UniProt entries under the Rhabditida
                   order (taxid 6236, ~2.08M sequences) — the other local
                   search target; a large one-time download, comparable in
                   scale to the interpro step above. Rhabditida (not the
                   full Nematoda phylum, taxid 6231, ~2.4M sequences) is the
                   order that actually contains C. elegans, giving more
                   taxonomically relevant hits per sequence searched.
  pfam      Pfam-A.hmm + Pfam-A.clans.tsv (current_release) — the local Pfam
            domain-scanning model library for scripts/pfam_scan.py, used by
            scripts/build_pfam_cache.py (bulk pre-scan) and scripts/domains_scan.py
            (per-protein backend, on-demand fallback for cache misses)

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
import subprocess
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.parent

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
    # Rhabditida order taxid, for the rhabditida_trembl ortholog reference database
    "rhabditida_taxid": 6236,
    # Source protein FASTA for scripts/build_pfam_cache.py's bulk Pfam scan —
    # organism-agnostic key so switching organisms only means pointing this
    # (and out_dir) at a different config, no code changes.
    "protein_fasta": str(
        _REPO_ROOT / "data" / "reference"
        / "c_elegans.canonical_bioproject.current.protein.fa.gz"
    ),
}


def _load_config():
    """Merge the active config's "reference_data" block (if any) over DEFAULTS."""
    # providers resolves TAGSITES_BATCH_CONFIG, so overrides apply here too
    from providers import _load_config as _load_batch_config

    cfg = dict(DEFAULTS)
    cfg.update(_load_batch_config().get("reference_data", {}))
    return cfg


def _out_dir(cfg):
    d = Path(cfg["out_dir"])
    d.mkdir(parents=True, exist_ok=True)
    return d


def _verify_gzip(path):
    """Raise if `path` isn't a valid, uncorrupted gzip stream.

    UniProt's /stream endpoint has been observed to close the connection early
    on a large uncompressed response without curl detecting it as an error (no
    Content-Length on a streamed response, so curl can't tell truncation from
    a normal end-of-stream) — the result is a silently truncated JSON file
    that "skip if already present" would then trust forever. Requesting the
    compressed form and verifying it here catches that: a truncated gzip
    member fails integrity checking, where a truncated plain-text stream
    would not have raised anything at all.
    """
    result = subprocess.run(["gzip", "-t", str(path)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"downloaded file {path} failed gzip integrity check (likely "
            f"truncated mid-transfer): {result.stderr.strip()}"
        )


def _curl_download(url, dest, force=False, extra_args=None, resumable=True, verify_gzip=False,
                   stream_retries=3):
    """Download via curl; skip if dest already has non-zero size and force isn't
    set. Returns True if a download ran, False if skipped.

    resumable=True adds -C - (byte-range resume) for static files that support
    HTTP Range (AFDB/Ensembl FTP). Set resumable=False for compute-on-the-fly
    endpoints (e.g. UniProt's /stream) that don't reliably support Range —
    resuming against one of those risks silently splicing two different result
    orderings together instead of actually continuing the same download.

    Non-resumable downloads use HTTP/1.1 (--http1.1) and retry the *whole*
    transfer up to stream_retries times on failure: curl's own --retry only
    covers a fresh connection attempt before any data arrives, not a stream
    that fails partway through — observed in practice as "curl: (92) HTTP/2
    stream ... INTERNAL_ERROR" after several minutes and tens of MB into a
    large UniProt /stream response. Each retry deletes whatever partial file
    is on disk first, since there's nothing to resume from.

    verify_gzip=True runs _verify_gzip() after a fresh download and deletes+
    raises on failure, so a truncated transfer can't masquerade as a
    successfully cached file on the next run.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {dest} already present ({dest.stat().st_size:,} bytes)")
        return False
    print(f"[fetch] {url} -> {dest}")
    cmd = ["curl", "-fL", "--retry", "5", "--retry-delay", "10"]
    if resumable:
        cmd += ["-C", "-"]
    else:
        cmd += ["--http1.1"]
    cmd += ["-o", str(dest), url]
    if extra_args:
        cmd[1:1] = extra_args

    attempts = stream_retries if not resumable else 1
    for attempt in range(1, attempts + 1):
        try:
            subprocess.run(cmd, check=True)
            break
        except subprocess.CalledProcessError:
            dest.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            print(f"[retry] transfer failed (attempt {attempt}/{attempts}); restarting {url}")

    if verify_gzip:
        try:
            _verify_gzip(dest)
        except RuntimeError:
            dest.unlink(missing_ok=True)
            raise
    return True


# ── UniProt proteome (sequences, CRC64s, curated features) ───────────────────

def fetch_uniprot(cfg, force=False):
    """Download the proteome's canonical+isoform FASTA and per-entry JSON stream.

    The FASTA gives sequences/CRC64s for local_store's checksum-keyed cache;
    the JSON stream gives curated features (uniprot_features.py's network
    calls) and WormBase/Ensembl cross-references (the accession<->GFF-name map
    genome_regions.py needs — see the plan's Phase A2).

    Both are requested compressed (format=...&compressed=true, saved as .gz)
    rather than plain text: UniProt's /stream endpoint has been observed to
    silently truncate a large plain response mid-transfer without curl
    reporting an error, whereas a truncated gzip stream fails integrity
    checking — see _verify_gzip()'s docstring.
    """
    out = _out_dir(cfg)
    pid = cfg["uniprot_proteome_id"]

    fasta_url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=proteome:{pid}&format=fasta&includeIsoform=true&compressed=true"
    )
    _curl_download(fasta_url, out / f"{pid}.fasta.gz", force=force,
                   resumable=False, verify_gzip=True)

    # JSON stream: one array of full entry records (features, xrefs, CRC64 checksum)
    json_url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=proteome:{pid}&format=json&compressed=true"
    )
    _curl_download(json_url, out / f"{pid}.json.gz", force=force,
                   resumable=False, verify_gzip=True)


# ── AlphaFold DB structures ───────────────────────────────────────────────────

def fetch_afdb(cfg, force=False):
    """Download the AFDB proteome tarball, then extract just its .pdb.gz
    members (skipping the .cif.gz half of each pair — the pipeline only ever
    parses PDB) into out_dir/afdb_pdbs/, named AF-<accession>-F1-model_v6.pdb.gz,
    so scripts/structure_bulk.py can open a specific accession's structure by
    deterministic filename in O(1) rather than scanning the tar (which has no
    index — a plain tarfile.extractfile(name) lookup is O(n) per call unless
    the whole member list has already been built once).
    """
    import tarfile

    out = _out_dir(cfg)
    url = cfg["afdb_tarball_url"]
    dest = out / Path(url).name
    downloaded = _curl_download(url, dest, force=force)

    pdb_dir = out / "afdb_pdbs"
    marker = pdb_dir / ".extracted"
    if marker.exists() and not force and not downloaded:
        print(f"[skip] {pdb_dir} already extracted")
        return

    pdb_dir.mkdir(parents=True, exist_ok=True)
    print(f"[extract] {dest} .pdb.gz members -> {pdb_dir}")
    n = 0
    with tarfile.open(dest, "r") as tar:
        for member in tar:
            if not member.name.endswith(".pdb.gz"):
                continue
            with tar.extractfile(member) as src, open(pdb_dir / member.name, "wb") as dst:
                dst.write(src.read())
            n += 1
            if n % 5000 == 0:
                print(f"[extract] {n:,} structures so far…")
    marker.write_text(f"{n}\n")
    print(f"[extract] done: {n:,} structures -> {pdb_dir}")


# ── Ensembl genome + annotation ───────────────────────────────────────────────

def fetch_genome(cfg, force=False):
    """Download the soft-masked toplevel genome FASTA (for bedtools/pyfaidx region extraction)."""
    out = _out_dir(cfg)
    release, species, assembly = cfg["ensembl_release"], cfg["ensembl_species"], cfg["ensembl_assembly"]
    fname = f"{species.capitalize()}.{assembly}.dna_sm.toplevel.fa.gz"
    url = f"https://ftp.ensembl.org/pub/release-{release}/fasta/{species}/dna/{fname}"
    _curl_download(url, out / fname, force=force, verify_gzip=True)


def fetch_gff3(cfg, force=False):
    """Download the GFF3 annotation (CDS exons per transcript — the genewise_bulk.py input)."""
    out = _out_dir(cfg)
    release, species, assembly = cfg["ensembl_release"], cfg["ensembl_species"], cfg["ensembl_assembly"]
    fname = f"{species.capitalize()}.{assembly}.{release}.gff3.gz"
    url = f"https://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{fname}"
    _curl_download(url, out / fname, force=force, verify_gzip=True)


# ── InterPro domain matches (filtered to this proteome's accessions) ─────────

def _load_proteome_accessions(cfg):
    """Read accessions out of the UniProt proteome FASTA already downloaded by fetch_uniprot()."""
    out = _out_dir(cfg)
    fasta_path = out / f"{cfg['uniprot_proteome_id']}.fasta.gz"
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"{fasta_path} not found — run fetch_uniprot() (or --only uniprot) first; "
            "the InterPro filter needs the proteome's accession list."
        )
    accessions = set()
    with gzip.open(fasta_path, "rt") as f:
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


def fetch_pfam(cfg, force=False):
    """Download Pfam-A's HMM library and clan table for scripts/pfam_scan.py's
    local de novo domain scanning (scripts/build_pfam_cache.py's bulk pass and
    scripts/domains_scan.py's on-demand fallback). Uses Pfam's "current_release"
    alias so this always pulls the latest models; the actual version string is
    recorded separately (see _record_pfam_version below) so a cache built
    against an older release can be detected and invalidated later.
    """
    out = _out_dir(cfg)
    base = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"

    hmm_gz = out / "Pfam-A.hmm.gz"
    hmm_path = out / "Pfam-A.hmm"
    clans_path = out / "Pfam-A.clans.tsv"

    _curl_download(f"{base}/Pfam-A.hmm.gz", hmm_gz, force=force, resumable=True)
    if not hmm_path.exists() or force:
        print(f"[gunzip] {hmm_gz} -> {hmm_path}")
        with gzip.open(hmm_gz, "rb") as src, open(hmm_path, "wb") as dst:
            dst.write(src.read())

    clans_gz = out / "Pfam-A.clans.tsv.gz"
    _curl_download(f"{base}/Pfam-A.clans.tsv.gz", clans_gz, force=force, resumable=True)
    if not clans_path.exists() or force:
        print(f"[gunzip] {clans_gz} -> {clans_path}")
        with gzip.open(clans_gz, "rb") as src, open(clans_path, "wb") as dst:
            dst.write(src.read())

    _record_pfam_version(base, out, force=force)


def _record_pfam_version(base_url, out_dir, force=False):
    """Fetch Pfam.version.gz and write out/Pfam.version.txt (plain text) — the
    release string scripts/build_pfam_cache.py stamps into its cache's
    _meta.json so a later Pfam-A.hmm.gz update invalidates stale hit caches.
    """
    dest = out_dir / "Pfam.version.txt"
    if dest.exists() and not force:
        return
    result = subprocess.run(
        ["curl", "-fsL", f"{base_url}/Pfam.version.gz"], check=True, capture_output=True,
    )
    dest.write_bytes(gzip.decompress(result.stdout))
    print(f"[pfam] recorded version info -> {dest}")


# ── Local ortholog reference DBs (conservation_local.py's DIAMOND search targets) ──

def _build_diamond_db(fasta_gz_path, db_path, force=False):
    """Run `diamond makedb` directly against a gzipped FASTA (DIAMOND reads
    gzip input natively — verified against the installed diamond v2.2.6 — so
    there's no need to decompress to disk first).
    """
    if db_path.exists() and not force:
        print(f"[skip] {db_path} already present")
        return
    print(f"[makedb] building DIAMOND database from {fasta_gz_path} -> {db_path}")
    subprocess.run(
        ["diamond", "makedb", "--in", str(fasta_gz_path), "-d", str(db_path.with_suffix(""))],
        check=True,
    )
    print(f"[makedb] done -> {db_path}")


def fetch_swissprot(cfg, force=False):
    """Download all of reviewed UniProt (Swiss-Prot, ~575k sequences) and
    build a DIAMOND database from it — one of two local search targets for
    scripts/conservation_local.py (see its docstring and batch.config.json's
    conservation_local.search_databases list), used in place of the EBI BLAST
    job scripts/conservation_remote.py submits.
    """
    out = _out_dir(cfg)
    fasta_path = out / "swissprot.fasta.gz"
    db_path = out / "swissprot.dmnd"
    url = "https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=fasta&compressed=true"
    _curl_download(url, fasta_path, force=force, resumable=False, verify_gzip=True)
    _build_diamond_db(fasta_path, db_path, force=force)


def fetch_rhabditida_trembl(cfg, force=False):
    """Download all unreviewed UniProt entries under the Rhabditida order
    (taxid 6236, ~2.08M sequences — a large one-time download, comparable in
    scale to the interpro step) and build a DIAMOND database from it — the
    other local search target for scripts/conservation_local.py. Restricted
    to reviewed:false so there's no overlap with fetch_swissprot()'s output
    (a reviewed nematode entry lives in Swiss-Prot, not TrEMBL).
    """
    out = _out_dir(cfg)
    taxid = cfg["rhabditida_taxid"]
    fasta_path = out / "rhabditida_trembl.fasta.gz"
    db_path = out / "rhabditida_trembl.dmnd"
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=taxonomy_id:{taxid}+AND+reviewed:false&format=fasta&compressed=true"
    )
    _curl_download(url, fasta_path, force=force, resumable=False, verify_gzip=True)
    _build_diamond_db(fasta_path, db_path, force=force)


STEPS = {
    "uniprot":           fetch_uniprot,
    "afdb":              fetch_afdb,
    "genome":            fetch_genome,
    "gff3":              fetch_gff3,
    "interpro":          fetch_interpro,   # depends on "uniprot" having run first
    "swissprot":         fetch_swissprot,
    "rhabditida_trembl": fetch_rhabditida_trembl,
    "pfam":              fetch_pfam,
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
