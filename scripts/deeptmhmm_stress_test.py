"""
deeptmhmm_stress_test.py

Standalone long-sequence stress test for DeepTMHMM, meant to be copied to the
Linux / RTX 3090 box and run there before the full-proteome topology pass.
Stdlib only: it needs no TAGSITES checkout, just this file, the protein FASTA,
and DeepTMHMM's checkout plus its conda env's python.

For each target isoform it runs DeepTMHMM's predict.py alone, polls nvidia-smi
for peak VRAM, and records wall time, peak host RSS, exit status and whether
TMRs.gff3 was produced (out-of-memory and other failures are caught, not fatal).
Results go to a JSON file to copy back, and a summary prints at the end,
including the max_length to set in batch.config.json's "deeptmhmm" block if the
largest target fails.

Targets, in order of precedence:
  --ladder          probe ascending lengths (nearest isoform in the FASTA to each
                    of 2000/4000/6000/8000/10000/13000/15188 aa), stopping at the
                    first failure; the largest pass is the recommended cap
  --names A,B       specific isoform names
  (default)         W06H8.8f, the 15,188 aa titin-like isoform

Usage
-----
    python deeptmhmm_stress_test.py \\
        --install-dir /path/to/DeepTMHMM-Academic-License-v1.0 \\
        --python /path/to/envs/deeptmhmm/bin/python \\
        --fasta c_elegans.canonical_bioproject.current.protein.fa.gz

    # if the default target fails, find the largest length that works
    python deeptmhmm_stress_test.py ... --ladder

    # CPU-only sanity check of the harness itself
    python deeptmhmm_stress_test.py ... --names T10H9.4
"""

import gzip
import json
import platform
import resource
import shutil
import subprocess
import tempfile
import threading
import time
from argparse import ArgumentParser
from pathlib import Path

DEFAULT_TARGET = "W06H8.8f"
LADDER_LENGTHS = [2000, 4000, 6000, 8000, 10000, 13000, 15188]
OOM_MARKERS = ("out of memory", "cuda error", "cublas", "killed")


def read_fasta_records(path):
    """Return [(name, sequence), ...] from a (optionally gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name, chunks = [], None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def torch_device_info(python_bin):
    """Ask the DeepTMHMM env's torch whether CUDA is available and on which GPU."""
    code = (
        "import torch; ok = torch.cuda.is_available(); "
        "print(ok, torch.__version__, torch.cuda.get_device_name(0) if ok else '-')"
    )
    out = subprocess.run([python_bin, "-c", code], capture_output=True, text=True, check=False)
    parts = out.stdout.split(None, 2)
    if len(parts) < 3:
        return {"cuda": False, "torch": "unknown", "gpu": "-", "error": out.stderr[-300:]}
    return {"cuda": parts[0] == "True", "torch": parts[1], "gpu": parts[2].strip()}


def gpu_memory_used_mib():
    """Current GPU 0 memory use in MiB via nvidia-smi, or None if unavailable."""
    if not shutil.which("nvidia-smi"):
        return None
    out = subprocess.run(
        ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits", "-i", "0"],
        capture_output=True,
        text=True,
        check=False,
    )
    try:
        return int(out.stdout.strip().splitlines()[0])
    except (ValueError, IndexError):
        return None


def poll_peak_vram(stop_event, peak, interval=0.5):
    """Track the maximum GPU memory use seen until stop_event is set."""
    while not stop_event.is_set():
        used = gpu_memory_used_mib()
        if used is not None and used > peak["mib"]:
            peak["mib"] = used
        stop_event.wait(interval)


def run_one(name, seq, install_dir, python_bin):
    """Run predict.py on one sequence; return a result dict (never raises on failure)."""
    tmp_root = tempfile.mkdtemp(prefix="dtm_stress_")
    result = {"name": name, "length": len(seq)}
    try:
        fasta_path = Path(tmp_root) / "input.fasta"
        fasta_path.write_text(f">{name}\n{seq}\n")
        out_dir = Path(tmp_root) / "out"

        baseline = gpu_memory_used_mib()
        peak = {"mib": baseline or 0}
        stop_event = threading.Event()
        poller = threading.Thread(target=poll_peak_vram, args=(stop_event, peak), daemon=True)
        poller.start()

        rss_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        start = time.time()
        proc = subprocess.run(
            [python_bin, "predict.py", "--fasta", str(fasta_path), "--output-dir", str(out_dir)],
            cwd=str(install_dir),
            capture_output=True,
            text=True,
            check=False,
        )
        result["wall_seconds"] = round(time.time() - start, 1)
        stop_event.set()
        poller.join()

        # ru_maxrss is KB on Linux, bytes on macOS; children max only rises, so a
        # value equal to the pre-run max means this run did not set a new high
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        divisor = 1024 if platform.system() == "Linux" else 1024 * 1024
        result["peak_host_rss_mib"] = round(rss / divisor) if rss > rss_before else None

        gff3 = out_dir / "TMRs.gff3"
        result["exit_code"] = proc.returncode
        result["tmrs_gff3"] = gff3.exists()
        result["peak_vram_mib"] = peak["mib"] if baseline is not None else None
        result["vram_baseline_mib"] = baseline
        stderr_tail = proc.stderr[-1500:]
        result["stderr_tail"] = stderr_tail
        result["oom_suspected"] = any(m in stderr_tail.lower() for m in OOM_MARKERS)
        # predict.py can crash in its plotting step after writing TMRs.gff3, so the
        # file's presence, not the exit code, decides pass/fail (see the plan notes)
        result["passed"] = gff3.exists()
        if gff3.exists():
            result["n_tmhelix"] = sum(
                1 for line in gff3.read_text().splitlines() if "\tTMhelix\t" in line
            )
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
    return result


def pick_targets(records, names, ladder):
    """Resolve targets to [(name, seq)], ascending by length for the ladder."""
    by_name = dict(records)
    if ladder:
        picked = {}
        for length in LADDER_LENGTHS:
            # nearest isoform by length to each rung; dedupe so rungs don't repeat
            name, seq = min(records, key=lambda r: abs(len(r[1]) - length))
            picked[name] = seq
        return sorted(picked.items(), key=lambda ns: len(ns[1]))
    missing = [n for n in names if n not in by_name]
    if missing:
        raise SystemExit(f"isoform(s) not in FASTA: {', '.join(missing)}")
    return [(n, by_name[n]) for n in names]


def format_result(r):
    """One-line human summary of a run_one result."""
    status = "PASS" if r["passed"] else ("OOM?" if r["oom_suspected"] else "FAIL")
    vram = f"{r['peak_vram_mib']} MiB" if r["peak_vram_mib"] is not None else "n/a"
    return (
        f"{status:5s} {r['name']:12s} {r['length']:6d} aa  {r['wall_seconds']:8.1f}s  "
        f"peak VRAM {vram}  peak host RSS {r['peak_host_rss_mib']} MiB"
    )


def main(install_dir, python_bin, fasta, names, ladder, out_json):
    """Run the stress test and write/print results."""
    install_dir = Path(install_dir)
    if not (install_dir / "predict.py").exists():
        raise SystemExit(f"predict.py not found in {install_dir}")

    device = torch_device_info(python_bin)
    print(f"host={platform.node()} device={device}")
    if not device["cuda"]:
        print("WARNING: CUDA not available to torch; this run measures CPU, not the 3090.")

    records = read_fasta_records(fasta)
    targets = pick_targets(records, names, ladder)
    print(f"targets: {', '.join(f'{n} ({len(s)} aa)' for n, s in targets)}\n")

    results = []
    for name, seq in targets:
        print(f"running {name} ({len(seq)} aa) ...", flush=True)
        r = run_one(name, seq, install_dir, python_bin)
        results.append(r)
        print("  " + format_result(r), flush=True)
        # the ladder stops at the first failure: longer rungs cannot pass either
        if ladder and not r["passed"]:
            break

    passed = [r for r in results if r["passed"]]
    failed = [r for r in results if not r["passed"]]
    summary = {
        "host": platform.node(),
        "device": device,
        "results": results,
        "largest_passing_length": max((r["length"] for r in passed), default=None),
        "smallest_failing_length": min((r["length"] for r in failed), default=None),
    }
    Path(out_json).write_text(json.dumps(summary, indent=2))

    print("\n== summary ==")
    for r in results:
        print(format_result(r))
    if failed:
        print(
            f"\nFAILED. Set deeptmhmm.max_length below {summary['smallest_failing_length']}"
            f" (largest passing: {summary['largest_passing_length']})."
        )
        print("stderr tail of first failure:\n" + failed[0]["stderr_tail"])
    else:
        print("\nAll targets passed: no max_length cap needed for these lengths.")
    print(f"\nresults written to {out_json}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--install-dir", required=True, help="DeepTMHMM checkout (has predict.py)")
    parser.add_argument("--python", required=True, help="python of the DeepTMHMM conda env")
    parser.add_argument("--fasta", required=True, help="protein FASTA (.gz ok)")
    parser.add_argument("--names", default=DEFAULT_TARGET, help="comma-separated isoform names")
    parser.add_argument("--ladder", action="store_true", help="probe ascending lengths")
    parser.add_argument("--out-json", default="deeptmhmm_stress_results.json")
    args = parser.parse_args()
    main(
        args.install_dir,
        args.python,
        args.fasta,
        [n.strip() for n in args.names.split(",") if n.strip()],
        args.ladder,
        args.out_json,
    )
