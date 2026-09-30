"""
topology_deeptmhmm.py

Local backend for the topology task: transmembrane helix / signal peptide /
inside-outside orientation prediction via scripts/deeptmhmm_runner.py
(DeepTMHMM 1.0, academic license), run in a dedicated conda env as a
subprocess.

Fast path: scripts/build_topology_cache.py's bulk pre-scan already wrote this
sequence's result to {out_dir}/topology_cache/{seq_name}_topology.txt — just
copy it through. Slow path (cache miss, e.g. a brand-new pasted sequence):
predict it on demand via deeptmhmm_runner.predict_topology() and write the
result to both outputfile and the cache path, so it's warm for next time.
Rows are (source, start, stop, description) same as every other range task,
with source always "DeepTMHMM" and description DeepTMHMM's own region word
(TMhelix/signal/inside/outside/...) — not translated into Phobius's
vocabulary. See ~/.claude/plans/i-m-considering-a-large-humble-sun.md.
"""

import sys
from pathlib import Path

from site_selection_util import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from deeptmhmm_runner import predict_topology
from reference_data import _load_config as _load_reference_data_config
from progress import report as _report, resolve_reporter


def _cache_path(seq_name):
    ref_cfg = _load_reference_data_config()
    cache_dir = Path(ref_cfg["out_dir"]) / "topology_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{seq_name}_topology.txt"


def _parse_cached_file(path):
    """Read a previously-written 4-col topology TSV back into row tuples."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(tuple(line.split("\t")))
    return rows


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Write fasta_in's DeepTMHMM topology regions to outputfile: reuse a
    cached bulk-scan result if one exists for this sequence name, else predict
    on demand. email/job_id_cb/resume_job_ids are accepted but unused (no
    network job is submitted) — kept so providers.resolve("topology") can call
    any backend with the identical signature.
    """
    reporter = resolve_reporter(report)
    seq_name, seq = read_fasta(fasta_in)
    cache_path = _cache_path(seq_name)

    if cache_path.exists():
        region_rows = _parse_cached_file(cache_path)
        source_note = "cached bulk scan"
    else:
        _report(reporter, f"No cached DeepTMHMM prediction for {seq_name} — predicting on demand.",
                stage="topology_deeptmhmm")
        regions = predict_topology([(seq_name, str(seq))]).get(seq_name, [])
        region_rows = [("DeepTMHMM", str(start), str(stop), region)
                       for region, start, stop in regions]
        with open(cache_path, "w") as f:
            for row in region_rows:
                print("\t".join(row), file=f)
        source_note = "on-demand prediction"

    with open(outputfile, "w") as f_out:
        for row in region_rows:
            print("\t".join(row), file=f_out)

    descriptions = sorted({row[3] for row in region_rows if len(row) > 3 and row[3]})
    if region_rows:
        shown = ", ".join(descriptions)
        summary = f"Found {len(region_rows)} topology region(s) for {seq_name} ({source_note}): {shown}"
    else:
        summary = f"Found 0 topology regions for {seq_name} ({source_note})"
    _report(reporter, summary, stage="done")
