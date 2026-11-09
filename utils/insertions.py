"""Query insertions: residues the query has and most aligned orthologs lack.

The conservation scorer writes -1000 for a query residue whose alignment column is mostly
gaps (scripts/score_conservation_py3.py), which the loaders read as NaN. This module recovers
how many orthologs lack the residue so scoring can credit deep insertions without crediting
shallow alignments, and turns the result into `query_insertion` range rows.
"""

import os
import re
import sys

import numpy as np

INSERTION_SOURCE = "query_insertion"
_DESC = "absent from {absent} spanning orthologs ({pct}% of spanning)"
_DESC_RE = re.compile(r"absent from (\d+) spanning orthologs \((\d+)% of spanning\)")


def _read_alignment(aln_path):
    """(names, rows) of an alignment file, parsed by the same readers the scorer used."""
    scripts_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import score_conservation_py3 as sc

    names, rows = sc.read_clustal_alignment(aln_path)
    if not names:
        names, rows = sc.read_fasta_alignment(aln_path)
    return names, rows


def _reference_name(jsd_path):
    """The `# reference sequence:` name from a .jsd header, or None."""
    with open(jsd_path) as fh:
        for line in fh:
            if line.startswith("# reference sequence:"):
                return line.split(":", 1)[1].strip()
            if not line.startswith("#"):
                break
    return None


def _jsd_rows(jsd_path):
    """Number of data rows (one per query residue) in a .jsd file."""
    with open(jsd_path) as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith("#"))


def _query_row(jsd_path, names, rows):
    """Index of the query in the alignment, or None when it cannot be told unambiguously.

    A .jsd the scorer wrote names its reference. A database export has no header, so the row
    is the alignment sequence whose residue count equals the .jsd's row count (the first on a tie).
    """
    ref = _reference_name(jsd_path)
    if ref is not None:
        return names.index(ref) if ref in names else None
    n = _jsd_rows(jsd_path)
    hits = [i for i, r in enumerate(rows) if len(r) - r.count("-") == n]
    # equal-length orthologs can tie; the query leads the alignment, so row 0 settles it
    if len(hits) == 1 or (hits and hits[0] == 0):
        return hits[0]
    return None


def absent_counts(jsd_path, aln_path):
    """Per query residue (jsd row order) the (n_absent, n_span) ortholog counts, or None.

    An ortholog spans a column only between its first and last aligned residue, so a
    fragment that ends before the column is not counted as lacking the residue.
    """
    try:
        names, rows = _read_alignment(aln_path)
        q = _query_row(jsd_path, names, rows)
        if q is None:
            return None
        grid = np.array([list(r) for r in rows])
    except Exception:
        return None
    present = grid != "-"
    others = np.delete(present, q, axis=0)
    if others.shape[0] == 0:
        return None
    # first / last aligned column of each ortholog
    first = np.where(others.any(axis=1), others.argmax(axis=1), grid.shape[1])
    last = np.where(others.any(axis=1),
                    grid.shape[1] - 1 - others[:, ::-1].argmax(axis=1), -1)
    cols = np.arange(grid.shape[1])
    spans = (cols[None, :] >= first[:, None]) & (cols[None, :] <= last[:, None])
    n_span = spans.sum(axis=0)
    n_absent = (spans & ~others).sum(axis=0)
    # keep only the columns where the query has a residue, in order = jsd rows
    keep = present[q]
    return n_absent[keep], n_span[keep]


def insertion_rows(jsd_path, aln_path, missing_pos):
    """Range rows (source, start, stop, description) for runs of insertion residues.

    `missing_pos` is the set of 1-based positions whose conservation is NaN. Each run keeps
    its weakest residue's counts, so a run qualifies for credit only if every residue does.
    """
    counts = absent_counts(jsd_path, aln_path)
    if counts is None:
        return []
    n_absent, n_span = counts
    rows, run = [], []
    # walk residues in order, closing a run at the first residue that is not an insertion
    for pos in range(1, len(n_absent) + 1):
        if pos in missing_pos and n_absent[pos - 1] > 0 and n_span[pos - 1] > 0:
            run.append(pos)
            continue
        if run:
            rows.append(_row(run, n_absent, n_span))
            run = []
    if run:
        rows.append(_row(run, n_absent, n_span))
    return rows


def _row(run, n_absent, n_span):
    """One range row for a run of consecutive residues."""
    absent = min(int(n_absent[p - 1]) for p in run)
    pct = min(int(100 * n_absent[p - 1] / n_span[p - 1]) for p in run)
    return [INSERTION_SOURCE, run[0], run[-1], _DESC.format(absent=absent, pct=pct)]


def parse_description(desc):
    """(n_absent, fraction) from a query_insertion description, or None."""
    m = _DESC_RE.search(str(desc))
    return (int(m.group(1)), int(m.group(2)) / 100.0) if m else None
