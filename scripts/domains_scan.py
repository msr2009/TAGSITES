"""
domains_scan.py

Local backend for domain annotation: genuine de novo Pfam-A scanning via
scripts/pfam_scan.py (pyhmmer), instead of the accession-keyed precomputed
lookup domains_bulk.py does. Every sequence gets a real answer regardless of
UniProt membership.

Fast path: scripts/build_pfam_cache.py's bulk pre-scan already wrote this
sequence's result to {out_dir}/pfam_scan_cache/{seq_name}_domains.txt — just
copy it through. Slow path (cache miss, e.g. a brand-new pasted sequence):
scan it on demand via pfam_scan.scan_sequences() and write the result to
both outputfile and the cache path, so it's warm for next time. Either way
this produces the same (source, start, stop, description) TSV domains_remote.py
does, with source always "Pfam" (see pfam_scan.py's module docstring for why
only Pfam is scanned).
"""

import sys
from pathlib import Path

from site_selection_util import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from pfam_scan import scan_sequences
from reference_data import _load_config as _load_reference_data_config
from progress import report as _report, resolve_reporter


def _cache_path(seq_name):
    ref_cfg = _load_reference_data_config()
    cache_dir = Path(ref_cfg["out_dir"]) / "pfam_scan_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{seq_name}_domains.txt"


def _parse_cached_file(path):
    """Read a previously-written 4-col domains TSV back into row tuples."""
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
    """Write fasta_in's Pfam domain hits to outputfile: reuse a cached bulk-scan
    result if one exists for this sequence name, else scan on demand.
    email/job_id_cb/resume_job_ids are accepted but unused (no network job is
    submitted) — kept so providers.resolve("domains") can call any backend
    with the identical signature.
    """
    reporter = resolve_reporter(report)
    seq_name, seq = read_fasta(fasta_in)
    cache_path = _cache_path(seq_name)

    if cache_path.exists():
        domain_rows = _parse_cached_file(cache_path)
        source_note = "cached bulk scan"
    else:
        _report(reporter, f"No cached Pfam scan for {seq_name} — scanning on demand.",
                stage="domains_scan")
        hits = scan_sequences([(seq_name, str(seq))]).get(seq_name, [])
        domain_rows = [("Pfam", str(start), str(stop), description)
                       for pfam_acc, description, start, stop, score in hits]
        with open(cache_path, "w") as f:
            for row in domain_rows:
                print("\t".join(row), file=f)
        source_note = "on-demand scan"

    with open(outputfile, "w") as f_out:
        for row in domain_rows:
            print("\t".join(row), file=f_out)

    descriptions = sorted({row[3] for row in domain_rows if len(row) > 3 and row[3]})
    if domain_rows:
        shown = ", ".join(descriptions[:5])
        if len(descriptions) > 5:
            shown += ", ..."
        summary = f"Found {len(domain_rows)} domain hit(s) for {seq_name} ({source_note}): {shown}"
    else:
        summary = f"Found 0 domain hits for {seq_name} ({source_note})"
    _report(reporter, summary, stage="done")
