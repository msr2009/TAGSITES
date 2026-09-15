"""
build_pfam_cache.py

One-time bulk pre-scan: runs scripts/pfam_scan.py's hmmsearch over every
sequence in the configured protein FASTA (reference_data.protein_fasta) and
writes one flat domains file per sequence into
{reference_data.out_dir}/pfam_scan_cache/{seq_name}_domains.txt — same
headerless 4-column (source, start, stop, description) format every other
domains backend writes. scripts/domains_scan.py reads this cache directory
as its fast path, falling back to an on-demand scan only for sequences not
covered here.

Organism-agnostic by design: both the source FASTA and the cache directory
come from config (reference_data.protein_fasta / reference_data.out_dir),
so running this for a different organism is a config change, not a code
change. See ~/.claude/plans/i-m-considering-a-large-humble-sun.md.

Usage
-----
    python scripts/build_pfam_cache.py
    python scripts/build_pfam_cache.py --force
"""

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pfam_scan import _load_scan_config, scan_sequences
from reference_data import _load_config as _load_reference_data_config


def _read_fasta_records(path):
    """Yield (name, sequence) pairs from a (optionally gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def _cache_dir(ref_cfg):
    d = Path(ref_cfg["out_dir"]) / "pfam_scan_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_domains_file(path, hits):
    """Write one sequence's Pfam hits as a headerless 4-col TSV; an empty
    file is the legal "scanned, zero hits" result (matches domains_bulk.py's
    convention for empty-but-resolved output).
    """
    with open(path, "w") as f:
        for pfam_acc, description, start, stop, _score in hits:
            print(f"Pfam\t{start}\t{stop}\t{description}", file=f)


def _write_meta(cache_dir, force=False):
    """Record the Pfam release + scan config this cache was built with, so a
    later Pfam-A.hmm.gz update or threshold change can be detected instead of
    silently serving stale hits.
    """
    meta_path = cache_dir / "_meta.json"
    if meta_path.exists() and not force:
        return
    ref_cfg = _load_reference_data_config()
    version_path = Path(ref_cfg["out_dir"]) / "Pfam.version.txt"
    pfam_release = version_path.read_text().splitlines()[0] if version_path.exists() else "unknown"
    scan_cfg = _load_scan_config()
    meta_path.write_text(json.dumps({"pfam_release": pfam_release, **scan_cfg}, indent=2))


def main(force=False):
    """Scan every sequence in the configured protein FASTA once, skipping any
    sequence whose cache file already exists unless force=True — cheap to
    re-run after new isoforms are added, since only new/missing ids get scanned.
    """
    ref_cfg = _load_reference_data_config()
    fasta_path = Path(ref_cfg["protein_fasta"])
    if not fasta_path.exists():
        raise FileNotFoundError(f"{fasta_path} not found (reference_data.protein_fasta)")

    cache_dir = _cache_dir(ref_cfg)
    _write_meta(cache_dir, force=force)

    all_records = list(_read_fasta_records(fasta_path))
    # single directory listing instead of one .exists() stat per isoform —
    # out_dir may be a network-mounted volume, where 28k+ individual stat()
    # calls in a loop is drastically slower than one listdir.
    already_cached = {p.name for p in cache_dir.iterdir()} if not force else set()
    todo = [(name, seq) for name, seq in all_records
            if force or f"{name}_domains.txt" not in already_cached]

    print(f"[build_pfam_cache] {len(all_records)} sequences total, "
          f"{len(all_records) - len(todo)} already cached, {len(todo)} to scan")
    if not todo:
        print("[build_pfam_cache] nothing to do")
        return

    results = scan_sequences(todo)
    for name, _seq in todo:
        _write_domains_file(cache_dir / f"{name}_domains.txt", results.get(name, []))

    n_with_hits = sum(1 for hits in results.values() if hits)
    print(f"[build_pfam_cache] done: {len(todo)} scanned, "
          f"{n_with_hits} with >=1 hit -> {cache_dir}")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="rescan every sequence even if already cached")
    args = parser.parse_args()
    main(force=args.force)
