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
model-load overhead once per batch instead of once total. A batch closes at
batch_size sequences or batch_residues total residues, whichever comes first,
so the long isoforms processed last do not blow up scratch or VRAM.

Robustness for unattended runs:
  - a failing batch is bisected down to single sequences, so one poison
    sequence costs only itself; it is logged to _failed.tsv and the run goes on
  - a sequence DeepTMHMM left out of its output is logged to _failed.tsv, never
    cached as "zero regions" (an empty cache file means "predicted, none found")
  - cache files are written atomically (.tmp then rename), so a killed run
    never leaves a truncated file that the resume check counts as done
  - sequences longer than max_length are not predicted and are logged to
    _skipped.tsv (rewritten on every run that sets a cap)
  - _meta.json keeps a "runs" list (host, device, counts, timing) so a cache
    mixing Mac CPU and Linux CUDA results stays auditable

Organism-agnostic by design: both the source FASTA and the cache directory
come from config (reference_data.protein_fasta / reference_data.out_dir), and
DeepTMHMM's device selection (CPU vs. CUDA) is entirely internal to
predict.py — this driver only *reports* the device in _meta.json.
Knobs live in the batch config's "deeptmhmm" block (batch_size,
batch_residues, max_length); CLI flags override them.

Running on the Linux / RTX 3090 box
-----------------------------------
    1. Env (see ~/.claude/plans/i-m-considering-a-large-humble-sun.md): python=3.11,
       pytorch<2.6 built for CUDA, fair-esm==0.4.0, PeptideBuilder, biopython,
       h5py, numpy, matplotlib, tqdm. DeepTMHMM's academic-license checkout stays
       outside the repo and must not be redistributed.
    2. Write a config override with that machine's deeptmhmm.install_dir/python and
       a local reference_data.out_dir / protein_fasta (only the protein FASTA is
       needed there), then:
           TAGSITES_BATCH_CONFIG=/path/to/override.json \\
               python scripts/build_topology_cache.py
    3. Copy results back — per-isoform files merge by plain copy, no merge step:
           rsync -a remote:OUT/topology_cache/*_topology.txt data/reference/topology_cache/
       and merge the remote _meta.json "runs" entries and _skipped.tsv/_failed.tsv.

Usage
-----
    python scripts/build_topology_cache.py                  # full proteome
    python scripts/build_topology_cache.py --limit 200       # benchmark
    python scripts/build_topology_cache.py --batch-size 100
    python scripts/build_topology_cache.py --max-residues 30000
    python scripts/build_topology_cache.py --max-length 5000
    python scripts/build_topology_cache.py --force
"""

import gzip
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from deeptmhmm_runner import _deeptmhmm_config, predict_topology
from reference_data import _load_config as _load_reference_data_config

DEFAULT_BATCH_SIZE = 200
DEFAULT_BATCH_RESIDUES = 60000


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
    """Return the topology cache directory, creating it and clearing any stale
    .tmp files left behind by a killed run.
    """
    d = Path(ref_cfg["out_dir"]) / "topology_cache"
    d.mkdir(parents=True, exist_ok=True)
    for stale in d.glob("*.tmp"):
        stale.unlink()
    return d


def _write_regions_file(path, regions):
    """Atomically write one sequence's topology regions as a headerless 4-col TSV;
    an empty file is the legal "predicted, zero regions" result (matches
    domains_bulk.py's convention for empty-but-resolved output).
    """
    tmp_path = Path(str(path) + ".tmp")
    with open(tmp_path, "w") as f:
        for region, start, stop in regions:
            print(f"DeepTMHMM\t{start}\t{stop}\t{region}", file=f)
    os.replace(tmp_path, path)


def _detect_device(python_bin):
    """Report "cuda" or "cpu" by asking the DeepTMHMM env's torch, the same
    check predict.py makes; "unknown" if torch can't be queried.
    """
    try:
        out = subprocess.run(
            [python_bin, "-c", "import torch; print(torch.cuda.is_available())"],
            capture_output=True,
            text=True,
            timeout=120,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    if out == "True":
        return "cuda"
    return "cpu" if out == "False" else "unknown"


def _read_meta(meta_path):
    """Load _meta.json as {"runs": [...]}, migrating the original single-run
    object (deeptmhmm_config/host/machine) into runs[0].
    """
    if not meta_path.exists():
        return {"runs": []}
    meta = json.loads(meta_path.read_text())
    if "runs" in meta:
        return meta
    # legacy single-run format from the first benchmark pass
    return {"runs": [{**meta, "device": "cpu", "note": "migrated from single-run _meta.json"}]}


def _append_meta_run(cache_dir, run):
    """Append one invocation's record to _meta.json's "runs" list."""
    meta_path = cache_dir / "_meta.json"
    meta = _read_meta(meta_path)
    meta["runs"].append(run)
    meta_path.write_text(json.dumps(meta, indent=2))


def _append_log(path, header, rows):
    """Append tab-separated rows to a log file, writing the header if new."""
    if not rows:
        return
    is_new = not path.exists()
    with open(path, "a") as f:
        if is_new:
            print(header, file=f)
        for row in rows:
            print("\t".join(str(x) for x in row), file=f)


def _make_batches(records, max_seqs, max_residues):
    """Split records into consecutive batches closing at max_seqs sequences or
    max_residues total residues; a single over-budget sequence gets its own batch.
    """
    batches, current, residues = [], [], 0
    for name, seq in records:
        # close the running batch before it would exceed either budget
        if current and (len(current) >= max_seqs or residues + len(seq) > max_residues):
            batches.append(current)
            current, residues = [], 0
        current.append((name, seq))
        residues += len(seq)
    if current:
        batches.append(current)
    return batches


def _predict_isolated(batch):
    """Predict a batch, bisecting on failure down to single sequences. Returns
    ({name: regions}, [(name, length, reason), ...]) where the failure list holds
    sequences that errored alone or that DeepTMHMM left out of its output.
    """
    try:
        results = predict_topology(batch)
    except RuntimeError as e:
        # a single sequence that still fails is the poison one: log it, move on
        if len(batch) == 1:
            name, seq = batch[0]
            return {}, [(name, len(seq), str(e)[-300:].replace("\n", " | "))]
        mid = len(batch) // 2
        left_ok, left_bad = _predict_isolated(batch[:mid])
        right_ok, right_bad = _predict_isolated(batch[mid:])
        return {**left_ok, **right_ok}, left_bad + right_bad
    ok = {name: results[name] for name, _ in batch if name in results}
    # a missing key means "not predicted", not "zero regions"
    missing = [(name, len(seq), "not in TMRs.gff3") for name, seq in batch if name not in results]
    return ok, missing


def main(force=False, limit=None, batch_size=None, max_residues=None, max_length=None):
    """Predict topology for every sequence in the configured protein FASTA,
    skipping any sequence whose cache file already exists unless force=True.
    Processes shortest-first so a --limit benchmark run returns a meaningful
    per-sequence rate quickly, rather than stalling on the longest isoforms
    first. Knob defaults come from the batch config's "deeptmhmm" block.
    """
    ref_cfg = _load_reference_data_config()
    fasta_path = Path(ref_cfg["protein_fasta"])
    if not fasta_path.exists():
        raise FileNotFoundError(f"{fasta_path} not found (reference_data.protein_fasta)")

    dtm_cfg = _deeptmhmm_config()
    batch_size = batch_size or dtm_cfg.get("batch_size") or DEFAULT_BATCH_SIZE
    max_residues = max_residues or dtm_cfg.get("batch_residues") or DEFAULT_BATCH_RESIDUES
    max_length = max_length or dtm_cfg.get("max_length")

    cache_dir = _cache_dir(ref_cfg)

    all_records = list(_read_fasta_records(fasta_path))
    # single directory listing instead of one .exists() stat per isoform —
    # out_dir may be a network-mounted volume, where thousands of individual
    # stat() calls in a loop is drastically slower than one listdir (see
    # build_pfam_cache.py's history for the real bug this avoids).
    already_cached = {p.name for p in cache_dir.iterdir()} if not force else set()
    todo = [
        (name, seq)
        for name, seq in all_records
        if force or f"{name}_topology.txt" not in already_cached
    ]
    n_already_cached = len(all_records) - len(todo)
    todo.sort(key=lambda ns: len(ns[1]))

    # over-cap sequences are recorded, never silently dropped
    skipped = []
    if max_length:
        skipped = [
            (name, len(seq), f"longer than max_length={max_length}")
            for name, seq in todo
            if len(seq) > max_length
        ]
        todo = [(name, seq) for name, seq in todo if len(seq) <= max_length]
        skipped_path = cache_dir / "_skipped.tsv"
        skipped_path.unlink(missing_ok=True)
        _append_log(skipped_path, "name\tlength\treason", skipped)

    if limit is not None:
        todo = todo[:limit]

    print(
        f"[build_topology_cache] {len(all_records)} sequences total, "
        f"{n_already_cached} already cached, {len(skipped)} skipped (> max_length), "
        f"{len(todo)} to predict"
    )
    if not todo:
        print("[build_topology_cache] nothing to do")
        return

    batches = _make_batches(todo, batch_size, max_residues)
    run = {
        "host": platform.node(),
        "machine": platform.machine(),
        "device": _detect_device(dtm_cfg["python"]),
        "started": datetime.now().isoformat(timespec="seconds"),
        "deeptmhmm_config": dtm_cfg,
        "max_length": max_length,
    }
    n_done, n_with_hits, n_failed = 0, 0, 0
    start_time = time.time()
    try:
        for i, batch in enumerate(batches):
            batch_time = time.time()
            ok, failed = _predict_isolated(batch)
            for name, regions in ok.items():
                _write_regions_file(cache_dir / f"{name}_topology.txt", regions)
                if regions:
                    n_with_hits += 1
            _append_log(cache_dir / "_failed.tsv", "name\tlength\treason", failed)
            n_done += len(ok)
            n_failed += len(failed)
            elapsed = time.time() - batch_time
            print(
                f"[build_topology_cache] batch {i + 1}/{len(batches)} "
                f"({len(batch)} seqs, {sum(len(s) for _, s in batch)} aa, "
                f"{len(failed)} failed): {elapsed:.1f}s ({elapsed / len(batch):.2f}s/seq)"
            )
    finally:
        # recorded even on Ctrl-C so the partial run is still attributed to a device
        run.update(
            {
                "finished": datetime.now().isoformat(timespec="seconds"),
                "n_predicted": n_done,
                "n_with_regions": n_with_hits,
                "n_failed": n_failed,
                "n_skipped": len(skipped),
                "wall_seconds": round(time.time() - start_time, 1),
            }
        )
        _append_meta_run(cache_dir, run)

    total_elapsed = time.time() - start_time
    print(
        f"[build_topology_cache] done: {n_done} predicted, {n_with_hits} with >=1 region, "
        f"{n_failed} failed, {len(skipped)} skipped -> {cache_dir} "
        f"({total_elapsed:.1f}s total, {total_elapsed / max(n_done, 1):.2f}s/seq avg)"
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true", help="repredict every sequence even if already cached"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="only predict the N shortest not-yet-cached sequences (benchmarking)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="max sequences per predict.py invocation "
        "(default: deeptmhmm.batch_size in config, else 200)",
    )
    parser.add_argument(
        "--max-residues",
        type=int,
        default=None,
        help="max total residues per predict.py invocation "
        "(default: deeptmhmm.batch_residues in config, else 60000)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="skip (and log to _skipped.tsv) sequences longer than this "
        "(default: deeptmhmm.max_length in config, else no cap)",
    )
    args = parser.parse_args()
    main(
        force=args.force,
        limit=args.limit,
        batch_size=args.batch_size,
        max_residues=args.max_residues,
        max_length=args.max_length,
    )
