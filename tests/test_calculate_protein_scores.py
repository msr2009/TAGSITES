"""
test_calculate_protein_scores.py — offline unit tests for scripts/calculate_protein_scores.py

Verifies:
  - load_scores loads all rows (not just the first)
  - calculate_property_sliding_window produces correct per-window averages
"""

import pytest
from pathlib import Path

import sys
from pathlib import Path as _P
# conftest already adds scripts/ to sys.path
from calculate_protein_scores import load_scores, calculate_property_sliding_window

REPO_ROOT = _P(__file__).parent.parent


# ── load_scores ───────────────────────────────────────────────────────────────

class TestLoadScores:

    def test_mini_all_rows_loaded(self, DATA):
        """Regression test: scores_mini.tsv has 3 rows; all must be in the dict."""
        scores = load_scores(DATA / "scores_mini.tsv")
        assert len(scores) == 3, f"expected 3 entries, got {len(scores)}"

    def test_mini_correct_values(self, DATA):
        scores = load_scores(DATA / "scores_mini.tsv")
        assert scores["A"] == pytest.approx(1.8)
        assert scores["C"] == pytest.approx(2.5)
        assert scores["D"] == pytest.approx(-3.5)

    def test_full_table_all_standard_aas(self):
        """The real hydrophobicity table must contain all 20 standard AAs."""
        scores = load_scores(REPO_ROOT / "tables" / "hydrophobicity_kyle-doolittle.tsv")
        standard = set("ACDEFGHIKLMNPQRSTVWY")
        missing = standard - set(scores.keys())
        assert not missing, f"Missing AAs in score table: {missing}"

    def test_full_table_values_are_floats(self):
        scores = load_scores(REPO_ROOT / "tables" / "hydrophobicity_kyle-doolittle.tsv")
        for aa, val in scores.items():
            assert isinstance(val, float), f"value for {aa} is not float: {val}"


# ── calculate_property_sliding_window ────────────────────────────────────────

class TestCalculatePropertySlidingWindow:

    def _simple_scores(self):
        return {"A": 1.0, "V": 2.0, "L": 3.0}

    def test_window1_equals_per_residue(self):
        """With window=1 each element equals the residue's own score."""
        scores = self._simple_scores()
        result = calculate_property_sliding_window("AVL", 1, scores)
        assert result == pytest.approx([1.0, 2.0, 3.0])

    def test_window3_average(self):
        """With window=3 the single result is the mean of all three values."""
        scores = self._simple_scores()
        result = calculate_property_sliding_window("AVL", 3, scores)
        assert result == pytest.approx([(1.0 + 2.0 + 3.0) / 3])

    def test_output_length_window1(self):
        scores = self._simple_scores()
        seq = "AVLAVL"
        result = calculate_property_sliding_window(seq, 1, scores)
        assert len(result) == len(seq)

    def test_output_length_window3(self):
        scores = self._simple_scores()
        seq = "AVLAVL"
        result = calculate_property_sliding_window(seq, 3, scores)
        assert len(result) == len(seq) - 3 + 1  # == 4

    def test_unknown_aa_scored_as_zero(self):
        """Amino acids absent from the score dict contribute 0 to the window."""
        scores = {"A": 1.0}  # X not in scores
        result = calculate_property_sliding_window("AX", 2, scores)
        # window = [A, X], A contributes 1.0, X contributes 0 → mean = 0.5
        assert result == pytest.approx([0.5])

    def test_real_sequence_with_full_table(self):
        """Smoke-test: full hydrophobicity table on SNB-1 sequence."""
        scores = load_scores(REPO_ROOT / "tables" / "hydrophobicity_kyle-doolittle.tsv")
        seq = "MDAQGDAGAQGGSQGGPRPSNKRLQQTQAQVDEVVGIMKVNVEKVLERDQKLSQ"
        result = calculate_property_sliding_window(seq, 5, scores)
        assert len(result) == len(seq) - 5 + 1
        assert all(isinstance(v, float) for v in result)
