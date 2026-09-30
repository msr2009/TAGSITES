"""
build_topology_cache.py

Bulk pre-scan / benchmark driver: runs scripts/deeptmhmm_runner.py's
predict_topology() over sequences in the configured protein FASTA
(reference_data.protein_fasta) and writes one flat topology file per sequence
into {reference_data.out_dir}/topology_cache/{seq_name}_topology.txt — same
headerless 4-column (source, start, stop, description) format every other
range task writes. scripts/topology_deeptmhmm.py reads this cache directory
as its fast path, falling back to an on-demand prediction only for sequences
not covered here.

Unlike build_pfam_cache.py (one hmmsearch call over the whole proteome),
DeepTMHMM predictions are batched: each predict.py subprocess invocation
loads its own ESM1b model and writes one embedding file per sequence into a
temp dir that is deleted after the batch completes (see
deeptmhmm_runner.py's docstring). Batching bounds that scratch to one batch
at a time rather than the whole proteome at once, at the cost of paying the
model-load overhead once per batch instead of once total.

Organism-agnostic by design: both the source FASTA and the cache directory
come from config (reference_data.protein_fasta / reference_data.out_dir), and
DeepTMHMM's device selection (CPU vs. CUDA) is entirely internal to
predict.py — this driver never mentions a device. See
~/.claude/plans/i-m-considering-a-large-humble-sun.md.

Usage
-----
    python scripts/build_topology_cache.py                  # full proteome
    python scripts/build_topology_cache.py --limit 200       # benchmark
    python scripts/build_topology_cache.py --batch-size 100
    python scripts/build_topology_cache.py --force
"""

import gzip
import json
import platform
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from deeptmhmm_runner import predict_topology, _deeptmhmm_config
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
    d = Path(ref_cfg["out_dir"]) / "topology_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_regions_file(path, regions):
    """Write one sequence's topology regions as a headerless 4-col TSV; an
    empty file is the legal "predicted, zero regions" result (matches
    domains_bulk.py's convention for empty-but-resolved output).
    """
    with open(path, "w") as f:
        for region, start, stop in regions:
            print(f"DeepTMHMM\t{start}\t{stop}\t{region}", file=f)


def _write_meta(cache_dir, force=False):
    """Record which machine/device produced this cache and the deeptmhmm config
    used, so a mixed CPU/CUDA cache (e.g. Mac benchmark + Linux/3090 full run)
    is auditable rather than silently indistinguishable.
    """
    meta_path = cache_dir / "_meta.json"
    if meta_path.exists() and not force:
        return
    meta_path.write_text(json.dumps({
        "deeptmhmm_config": _deeptmhmm_config(),
        "host": platform.node(),
        "machine": platform.machine(),
    }, indent=2))


def main(force=False, limit=None, batch_size=200):
    """Predict topology for every sequence in the configured protein FASTA,
    skipping any sequence whose cache file already exists unless force=True.
    Processes shortest-first so a --limit benchmark run returns a meaningful
    per-sequence rate quickly, rather than stalling on the longest isoforms
    first.
    """
    ref_cfg = _load_reference_data_config()
    fasta_path = Path(ref_cfg["protein_fasta"])
    if not fasta_path.exists():
        raise FileNotFoundError(f"{fasta_path} not found (reference_data.protein_fasta)")

    cache_dir = _cache_dir(ref_cfg)
    _write_meta(cache_dir, force=force)

    all_records = list(_read_fasta_records(fasta_path))
    # single directory listing instead of one .exists() stat per isoform —
    # out_dir may be a network-mounted volume, where thousands of individual
    # stat() calls in a loop is drastically slower than one listdir (see
    # build_pfam_cache.py's history for the real bug this avoids).
    already_cached = {p.name for p in cache_dir.iterdir()} if not force else set()
    todo = [(name, seq) for name, seq in all_records
            if force or f"{name}_topology.txt" not in already_cached]
    todo.sort(key=lambda ns: len(ns[1]))

    if limit is not None:
        todo = todo[:limit]

    print(f"[build_topology_cache] {len(all_records)} sequences total, "
          f"{len(all_records) - len(todo) if limit is None else '?'} already cached, "
          f"{len(todo)} to predict")
    if not todo:
        print("[build_topology_cache] nothing to do")
        return

    n_done, n_with_hits = 0, 0
    start_time = time.time()
    for batch_start in range(0, len(todo), batch_size):
        batch = todo[batch_start:batch_start + batch_size]
        batch_time = time.time()
        results = predict_topology(batch)
        for name, _seq in batch:
            regions = results.get(name, [])
            _write_regions_file(cache_dir / f"{name}_topology.txt", regions)
            if regions:
                n_with_hits += 1
        n_done += len(batch)
        elapsed = time.time() - batch_time
        print(f"[build_topology_cache] batch {batch_start}-{batch_start + len(batch)}: "
              f"{elapsed:.1f}s ({elapsed / len(batch):.2f}s/seq)")

    total_elapsed = time.time() - start_time
    print(f"[build_topology_cache] done: {n_done} predicted, {n_with_hits} with >=1 region "
          f"-> {cache_dir} ({total_elapsed:.1f}s total, {total_elapsed / n_done:.2f}s/seq avg)")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="repredict every sequence even if already cached")
    parser.add_argument("--limit", type=int, default=None,
                        help="only predict the N shortest not-yet-cached sequences (benchmarking)")
    parser.add_argument("--batch-size", type=int, default=200,
                        help="sequences per predict.py subprocess invocation (default 200)")
    args = parser.parse_args()
    main(force=args.force, limit=args.limit, batch_size=args.batch_size)
