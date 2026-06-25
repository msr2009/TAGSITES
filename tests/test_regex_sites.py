"""
test_regex_sites.py — offline unit tests for scripts/regex_sites.py

Uses the snb-1.fa fixture and a purpose-built sites_mini.txt.
Output format: modification<TAB>start+1<TAB>end+1<TAB>description (1-based).
"""

import io
import pytest

from regex_sites import main as regex_main


# ── helpers ───────────────────────────────────────────────────────────────────

def _run_regex(fasta_path, sites_path):
    """Run regex_main and return output rows as list of split lists."""
    buf = io.StringIO()
    regex_main(str(fasta_path), str(sites_path), buf)
    buf.seek(0)
    rows = [line.strip().split("\t") for line in buf if line.strip()]
    return rows


# ── tests ─────────────────────────────────────────────────────────────────────

class TestRegexSitesMain:

    def test_n_terminal_met_one_row(self, DATA):
        """^M matches only once at the N-terminus."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        met_rows = [r for r in rows if r[3] == "n_terminal_met"]
        assert len(met_rows) == 1

    def test_n_terminal_met_coordinates(self, DATA):
        """^M at position 1 → start=1, end=2 (both 1-based, end exclusive)."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        met_row = [r for r in rows if r[3] == "n_terminal_met"][0]
        assert met_row[1] == "1"   # start+1 = 1
        assert met_row[2] == "2"   # end+1 = 2

    def test_type_column_is_modification(self, DATA):
        """First column is always 'modification'."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        for row in rows:
            assert row[0] == "modification"

    def test_lysine_count(self, DATA):
        """SNB-1 contains 10 K residues; 'lysine' pattern K should yield 10 rows."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        lys_rows = [r for r in rows if r[3] == "lysine"]
        assert len(lys_rows) == 10

    def test_lysine_first_position(self, DATA):
        """First K in SNB-1 is at 1-based position 22."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        lys_rows = sorted([r for r in rows if r[3] == "lysine"], key=lambda r: int(r[1]))
        assert lys_rows[0][1] == "22"

    def test_no_match_produces_no_rows_for_absent_pattern(self, DATA, tmp_path):
        """A pattern that matches nothing produces no output rows for that description."""
        sites = tmp_path / "nomatch.txt"
        sites.write_text("impossible\tZZZZZ\n")
        rows = _run_regex(DATA / "snb-1.fa", sites)
        assert rows == []

    def test_columns_have_four_fields(self, DATA):
        """Every output row should have exactly 4 tab-separated fields."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        for row in rows:
            assert len(row) == 4, f"Row has {len(row)} fields: {row}"

    def test_coordinates_are_numeric(self, DATA):
        """start and end columns should be parseable as ints."""
        rows = _run_regex(DATA / "snb-1.fa", DATA / "sites_mini.txt")
        for row in rows:
            int(row[1])
            int(row[2])
