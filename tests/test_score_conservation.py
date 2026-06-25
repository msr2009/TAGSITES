"""
test_score_conservation.py — offline unit tests for scripts/score_conservation_py3.py

Tests cover:
  read_fasta_alignment   — parse a FASTA-format multiple alignment
  js_divergence          — Jensen-Shannon divergence score per column
  weighted_freq_count_pseudocount — frequency counting helper
  gap_percentage         — gap fraction in a column
"""

import math
import pytest

from score_conservation_py3 import (
    read_fasta_alignment,
    js_divergence,
    weighted_freq_count_pseudocount,
    gap_percentage,
    PSEUDOCOUNT,
)

# Flat uniform background distribution over 20 amino acids (no gap)
_BG_FLAT20 = [1.0 / 20.0] * 20


# ── read_fasta_alignment ──────────────────────────────────────────────────────

class TestReadFastaAlignment:

    def test_returns_names_and_alignment(self, DATA):
        names, aln = read_fasta_alignment(str(DATA / "mini_alignment.fa"))
        assert isinstance(names, list)
        assert isinstance(aln, list)

    def test_correct_sequence_count(self, DATA):
        names, aln = read_fasta_alignment(str(DATA / "mini_alignment.fa"))
        assert len(names) == 3
        assert len(aln) == 3

    def test_sequences_same_length(self, DATA):
        _, aln = read_fasta_alignment(str(DATA / "mini_alignment.fa"))
        lengths = {len(s) for s in aln}
        assert len(lengths) == 1, f"unequal sequence lengths: {lengths}"

    def test_names_match_headers(self, DATA):
        names, _ = read_fasta_alignment(str(DATA / "mini_alignment.fa"))
        assert "seq1" in names
        assert "seq2" in names

    def test_sequences_are_uppercase(self, DATA):
        _, aln = read_fasta_alignment(str(DATA / "mini_alignment.fa"))
        for seq in aln:
            assert seq == seq.upper(), "alignment contains lowercase"


# ── js_divergence ─────────────────────────────────────────────────────────────

class TestJsDivergence:

    def _weights(self, n):
        return [1.0] * n

    def test_fully_conserved_column_high_score(self):
        """All 'A' → column is highly conserved → JSD should be > 0."""
        col = ["A", "A", "A", "A", "A"]
        score = js_divergence(col, None, _BG_FLAT20, self._weights(len(col)))
        assert score > 0.0

    def test_fully_conserved_greater_than_diverse(self):
        """A conserved column should score higher than a diverse one."""
        col_conserved = ["A", "A", "A", "A"]
        col_diverse = ["A", "R", "N", "D"]
        s_conserved = js_divergence(col_conserved, None, _BG_FLAT20, self._weights(4))
        s_diverse = js_divergence(col_diverse, None, _BG_FLAT20, self._weights(4))
        assert s_conserved > s_diverse

    def test_returns_float(self):
        col = ["A", "A", "V"]
        score = js_divergence(col, None, _BG_FLAT20, self._weights(3))
        assert isinstance(score, float)

    def test_score_nonnegative(self):
        col = ["A", "R", "N", "D", "C"]
        score = js_divergence(col, None, _BG_FLAT20, self._weights(5))
        assert score >= 0.0

    def test_gap_penalizes_score(self):
        """A column with gaps should score ≤ the same column without gaps."""
        col_no_gap = ["A", "A", "A"]
        col_with_gap = ["A", "A", "-"]
        s_no = js_divergence(col_no_gap, None, _BG_FLAT20, [1.0] * 3)
        s_gap = js_divergence(col_with_gap, None, _BG_FLAT20, [1.0] * 3)
        assert s_no >= s_gap


# ── weighted_freq_count_pseudocount ──────────────────────────────────────────

class TestWeightedFreqCount:

    def test_sums_to_one(self):
        col = ["A", "R", "N"]
        weights = [1.0, 1.0, 1.0]
        fc = weighted_freq_count_pseudocount(col, weights, PSEUDOCOUNT)
        assert abs(sum(fc) - 1.0) < 1e-6

    def test_length_is_21(self):
        col = ["A"]
        fc = weighted_freq_count_pseudocount(col, [1.0], PSEUDOCOUNT)
        assert len(fc) == 21  # 20 AAs + gap


# ── gap_percentage ────────────────────────────────────────────────────────────

class TestGapPercentage:

    def test_no_gaps(self):
        assert gap_percentage(["A", "R", "N"]) == pytest.approx(0.0)

    def test_all_gaps(self):
        assert gap_percentage(["-", "-"]) == pytest.approx(1.0)

    def test_half_gaps(self):
        assert gap_percentage(["A", "-"]) == pytest.approx(0.5)
