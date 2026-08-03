"""isoform_align.py — map a UniProt isoform sequence onto query residue coordinates.

Isoform pairs are near-identical (same gene, alternative splicing): a shared core plus
splice-in/splice-out differences and, sometimes, divergent unique N/C-termini from
alternative first/last exons. Local alignment (not global) is used deliberately — a
global alignment with a gap-open penalty will happily substitute a same-length unique
terminus residue-by-residue rather than pay for opening a gap, hiding real divergence
as a false match. Local alignment naturally stops extending once a region stops paying
for itself, so divergent termini simply fall outside the aligned block.
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
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def align_isoform_to_query(query_seq, iso_seq):
    """Align an isoform sequence to the query and map it onto query coordinates.

    Returns a dict: {"present": [[start, stop], ...], "skipped": [[start, stop], ...],
    "inserts": [{"after": pos, "length": n}, ...]}. Positions are 1-based, inclusive,
    in query coordinates. "after" is the query position immediately before the
    insertion point (0 if the insert precedes the query's first residue).
    """
    aligner = _make_aligner()
    alignment = aligner.align(query_seq, iso_seq)[0]  # best-scoring local alignment
    aligned_query, aligned_iso = alignment[0], alignment[1]
    q_blocks, i_blocks = alignment.aligned
    q_start, q_end = int(q_blocks[0][0]), int(q_blocks[-1][1])  # 0-based, half-open
    i_start, i_end = int(i_blocks[0][0]), int(i_blocks[-1][1])

    present_positions = []
    skipped_positions = list(range(1, q_start + 1))  # query residues before the aligned block
    inserts = []
    if i_start:  # isoform residues before the aligned block: a unique N-terminus/extension
        inserts.append({"after": q_start, "length": i_start})

    qpos = q_start  # 1-based query position, incremented as we consume query residues
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
            if h == q:
                present_positions.append(qpos)
            else:
                skipped_positions.append(qpos)
        else:
            # isoform residue with no corresponding query position: an insertion
            if h != "-":
                cur_insert_len += 1
    _flush_insert()

    skipped_positions.extend(range(q_end + 1, len(query_seq) + 1))  # query residues after the block
    if len(iso_seq) - i_end:  # isoform residues after the aligned block: a unique C-terminus
        inserts.append({"after": q_end, "length": len(iso_seq) - i_end})

    return {
        "present": _spans_from_positions(sorted(present_positions)),
        "skipped": _spans_from_positions(sorted(set(skipped_positions))),
        "inserts": inserts,
    }
