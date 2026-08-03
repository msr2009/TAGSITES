"""isoform_align.py — map a UniProt isoform sequence onto query residue coordinates.

Isoform pairs are near-identical (same gene, alternative splicing), so a global
alignment with a stiff gap-open penalty cleanly separates shared blocks from
splice-in/splice-out differences — no need for the heavier machinery used for
cross-species orthologs (Clustal Omega via EBI).
"""

from Bio.Align import PairwiseAligner, substitution_matrices


def _spans_from_positions(positions):
    """Collapse a sorted list of 1-based positions into inclusive [start, stop] spans."""
    if not positions:
        return []
    spans = []
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            spans.append([start, prev])
            start = prev = pos
    spans.append([start, prev])
    return spans


def _make_aligner():
    """Build a PairwiseAligner tuned for near-identical isoform pairs."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def align_isoform_to_query(query_seq, iso_seq):
    """Align an isoform sequence to the query and map it onto query coordinates.

    Returns a dict: {"present": [[start, stop], ...], "skipped": [[start, stop], ...],
    "inserts": [{"after": pos, "length": n, "seq": str}, ...]}. Positions are 1-based,
    inclusive, in query coordinates. "after" is the query position immediately before
    the insertion point (0 if the insert precedes the query's first residue).
    """
    aligner = _make_aligner()
    alignment = aligner.align(query_seq, iso_seq)[0]  # best-scoring alignment
    aligned_query, aligned_iso = alignment[0], alignment[1]

    present_positions = []
    skipped_positions = []
    inserts = []
    qpos = 0          # 1-based query position, incremented as we consume query residues
    cur_insert_len = 0

    def _flush_insert():
        nonlocal cur_insert_len
        if cur_insert_len:
            inserts.append({"after": qpos, "length": cur_insert_len})
            cur_insert_len = 0

    for q, h in zip(aligned_query, aligned_iso):
        if q != "-":
            _flush_insert()
            qpos += 1
            if h != "-":
                present_positions.append(qpos)
            else:
                skipped_positions.append(qpos)
        else:
            # isoform residue with no corresponding query position: an insertion
            if h != "-":
                cur_insert_len += 1
    _flush_insert()

    return {
        "present": _spans_from_positions(present_positions),
        "skipped": _spans_from_positions(skipped_positions),
        "inserts": inserts,
    }
