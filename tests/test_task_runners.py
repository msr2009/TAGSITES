"""
test_task_runners.py — offline unit tests for scripts/task_runners.py and scripts/progress.py

Covers: _str, _float, _int, _bool coercion helpers (task_runners.py);
        parse_line (progress.py).
"""

import pytest
from unittest.mock import patch

from task_runners import _str, _float, _int, _bool, afdb_presearch
from progress import parse_line


# ── _str ──────────────────────────────────────────────────────────────────────

class TestStr:

    def test_converts_int_to_str(self):
        assert _str(42) == "42"

    def test_converts_float_to_str(self):
        assert _str(3.14) == "3.14"

    def test_none_returns_default(self):
        assert _str(None) == ""

    def test_empty_string_returns_default(self):
        assert _str("") == ""

    def test_custom_default_on_none(self):
        assert _str(None, default="x") == "x"

    def test_non_empty_string_passed_through(self):
        assert _str("hello") == "hello"


# ── _float ────────────────────────────────────────────────────────────────────

class TestFloat:

    def test_converts_string_float(self):
        assert _float("3.14") == pytest.approx(3.14)

    def test_converts_int_to_float(self):
        assert _float(5) == pytest.approx(5.0)

    def test_none_returns_default(self):
        assert _float(None) == pytest.approx(0.0)

    def test_bad_string_returns_default(self):
        assert _float("not-a-number") == pytest.approx(0.0)

    def test_custom_default(self):
        assert _float(None, default=-1.0) == pytest.approx(-1.0)

    def test_scientific_notation(self):
        assert _float("1e-10") == pytest.approx(1e-10)


# ── _int ──────────────────────────────────────────────────────────────────────

class TestInt:

    def test_converts_string_int(self):
        assert _int("7") == 7

    def test_converts_float_to_int(self):
        assert _int(3.9) == 3

    def test_none_returns_default(self):
        assert _int(None) == 0

    def test_bad_string_returns_default(self):
        assert _int("nope") == 0

    def test_custom_default(self):
        assert _int(None, default=99) == 99


# ── _bool ─────────────────────────────────────────────────────────────────────

class TestBool:

    def test_true_bool_passthrough(self):
        assert _bool(True) is True

    def test_false_bool_passthrough(self):
        assert _bool(False) is False

    def test_string_true_case_insensitive(self):
        assert _bool("True") is True
        assert _bool("true") is True
        assert _bool("TRUE") is True

    def test_string_false(self):
        assert _bool("False") is False
        assert _bool("false") is False

    def test_string_zero(self):
        assert _bool("0") is False

    def test_empty_string_is_false(self):
        assert _bool("") is False

    def test_none_returns_default(self):
        assert _bool(None, default=False) is False

    def test_none_custom_default_true(self):
        assert _bool(None, default=True) is True

    def test_nonzero_int_is_true(self):
        assert _bool(1) is True

    def test_zero_int_is_false(self):
        assert _bool(0) is False

    def test_any_non_false_string_is_true(self):
        assert _bool("yes") is True
        assert _bool("1") is True
        assert _bool("anything") is True


# ── afdb_presearch ────────────────────────────────────────────────────────────

class TestAfdbPresearch:
    """afdb_presearch skips AFDB lookup when the plddt task already has a pdb path."""

    def _task(self, pdb=""):
        return {"id": "X_plddt", "analysis": "plddt", "args": {
            "pdb": pdb, "working_dir": "/wd", "run_name": "run",
            "input_file": "/wd/run.fa", "email": "t@example.com",
        }}

    def test_skip_when_pdb_provided(self):
        """No AFDB call when plddt task has a non-empty pdb path."""
        with patch("existing_AF_model.search_AFDB") as mock_search:
            result = afdb_presearch([self._task(pdb="/wd/run.AF.pdb")])
        mock_search.assert_not_called()
        assert result is False

    def test_skip_when_no_plddt_task(self):
        """No AFDB call when there is no plddt task in the list."""
        blast = {"id": "X_blast", "analysis": "blast", "args": {}}
        with patch("existing_AF_model.search_AFDB") as mock_search:
            result = afdb_presearch([blast])
        mock_search.assert_not_called()
        assert result is False

    def test_triggers_when_pdb_empty(self):
        """AFDB search runs when plddt task has an empty pdb (no file supplied)."""
        with patch("existing_AF_model.search_AFDB", return_value=1) as mock_search:
            result = afdb_presearch([self._task(pdb="")])
        mock_search.assert_called_once()
        assert result is False  # return_value=1 means no AFDB hit found


# ── parse_line ────────────────────────────────────────────────────────────────

class TestParseLine:

    def test_tagged_info_line(self):
        stage, msg, level = parse_line("[blast_submit] Submitting job to EBI")
        assert stage == "blast_submit"
        assert msg == "Submitting job to EBI"
        assert level == "info"

    def test_tagged_warning_line(self):
        stage, msg, level = parse_line("[ebi_poll] WARNING: job queued, retrying")
        assert stage == "ebi_poll"
        assert msg == "job queued, retrying"
        assert level == "warning"

    def test_tagged_error_line(self):
        stage, msg, level = parse_line("[blast_submit] ERROR: bad request")
        assert stage == "blast_submit"
        assert msg == "bad request"
        assert level == "error"

    def test_untagged_line_passes_through(self):
        stage, msg, level = parse_line("Some library noise")
        assert stage == ""
        assert msg == "Some library noise"
        assert level == "info"

    def test_timestamp_prefix_stripped(self):
        stage, msg, level = parse_line("12:34:56 [blast_submit] done")
        assert stage == "blast_submit"
        assert msg == "done"
        assert level == "info"

    def test_no_closing_bracket_treated_as_untagged(self):
        stage, msg, level = parse_line("[unclosed message")
        assert stage == ""

    def test_trailing_newline_stripped(self):
        stage, msg, level = parse_line("[stage] message\n")
        assert msg == "message"

    def test_empty_string(self):
        stage, msg, level = parse_line("")
        assert stage == ""
        assert level == "info"
