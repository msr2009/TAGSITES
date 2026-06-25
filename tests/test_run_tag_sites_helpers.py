"""
test_run_tag_sites_helpers.py — offline unit tests for pure helpers in
scripts/run_tag_sites_from_json.py

Tests cover:
  build_task_args_string  — CLI argument string builder
  searchAFDB_required     — detects plddt task needing AFDB search; mutates JSON
  genewise_required       — detects reagents task needing Genewise run; read-only
"""

import copy
import pytest

from run_tag_sites_from_json import (
    build_task_args_string,
    searchAFDB_required,
    genewise_required,
)


# ── build_task_args_string ────────────────────────────────────────────────────

class TestBuildTaskArgsString:

    def test_basic(self):
        result = build_task_args_string({"a": 1, "b": "", "c": "x"})
        assert result == " --a 1 --c x"

    def test_empty_string_values_skipped(self):
        result = build_task_args_string({"x": "", "y": ""})
        assert result == ""

    def test_empty_dict(self):
        assert build_task_args_string({}) == ""

    def test_exclude_removes_key(self):
        result = build_task_args_string({"a": 1, "b": 2, "c": 3}, EXCLUDE=["b"])
        assert "--b" not in result
        assert "--a 1" in result
        assert "--c 3" in result

    def test_integer_zero_is_kept(self):
        """Only literal empty string '' is skipped; 0 must be included."""
        result = build_task_args_string({"flag": 0})
        assert "--flag 0" in result

    def test_false_is_kept(self):
        result = build_task_args_string({"flag": False})
        assert "--flag False" in result

    def test_insertion_order_preserved(self):
        """Args appear in dict insertion order."""
        result = build_task_args_string({"z": 1, "a": 2, "m": 3})
        assert result.index("--z") < result.index("--a") < result.index("--m")

    def test_leading_space(self):
        """Result for non-empty input always starts with a space (ready to append)."""
        result = build_task_args_string({"a": 1})
        assert result.startswith(" ")


# ── searchAFDB_required ───────────────────────────────────────────────────────

def _make_plddt_json(pdb="", species_taxid=None):
    """Build a minimal JSON dict with one plddt task."""
    args = {
        "pdb": pdb,
        "output": "test.plddt.txt",
        "existing_AF2": 0,
        "input_file": "test.fa",
    }
    if species_taxid is not None:
        args["species_taxid"] = species_taxid
    return {
        "PLDDT_plddt": {
            "analysis": "plddt",
            "args": args,
        }
    }


class TestSearchAFDBRequired:

    def test_returns_args_string_when_pdb_empty(self):
        j = _make_plddt_json(pdb="", species_taxid=6239)
        result = searchAFDB_required(j)
        assert result is not None
        assert isinstance(result, str)

    def test_returns_none_when_pdb_set(self):
        j = _make_plddt_json(pdb="existing.pdb")
        result = searchAFDB_required(j)
        assert result is None

    def test_returns_none_no_plddt_task(self):
        j = {"OTHER_scores": {"analysis": "scores", "args": {"pdb": ""}}}
        result = searchAFDB_required(j)
        assert result is None

    def test_mutates_taxid_from_species_taxid(self):
        j = _make_plddt_json(pdb="", species_taxid=6239)
        searchAFDB_required(j)
        assert j["PLDDT_plddt"]["args"]["taxid"] == 6239

    def test_taxid_defaults_to_1_when_missing(self):
        j = _make_plddt_json(pdb="")  # no species_taxid
        searchAFDB_required(j)
        assert j["PLDDT_plddt"]["args"]["taxid"] == 1

    def test_pdb_output_existing_af2_excluded_from_result(self):
        j = _make_plddt_json(pdb="", species_taxid=1)
        result = searchAFDB_required(j)
        # These keys must be excluded from the returned args string
        assert "--pdb" not in result
        assert "--output" not in result
        assert "--existing_AF2" not in result

    def test_input_file_included_in_result(self):
        j = _make_plddt_json(pdb="", species_taxid=1)
        result = searchAFDB_required(j)
        assert "--input_file" in result


# ── genewise_required ─────────────────────────────────────────────────────────

def _make_reagents_json(genewise="", genomic_fasta="test_genomic.fa"):
    """Build a minimal JSON dict with one reagents task."""
    return {
        "REAGENTS_reagents": {
            "analysis": "reagents",
            "args": {
                "genewise": genewise,
                "genomic_fasta": genomic_fasta,
                "n_guides": 3,
            }
        }
    }


class TestGenewiseRequired:

    def test_returns_key_when_genewise_empty(self):
        j = _make_reagents_json(genewise="")
        key = genewise_required(j, {}, {}, "scripts/", "python")
        assert key == "REAGENTS_reagents"

    def test_returns_key_when_genewise_key_absent(self):
        j = {"REAGENTS_reagents": {"analysis": "reagents", "args": {"n_guides": 3}}}
        key = genewise_required(j, {}, {}, "scripts/", "python")
        assert key == "REAGENTS_reagents"

    def test_returns_none_when_genewise_set(self):
        j = _make_reagents_json(genewise="run/snb-1.genewise.out.txt")
        key = genewise_required(j, {}, {}, "scripts/", "python")
        assert key is None

    def test_returns_none_no_reagents_task(self):
        j = {"OTHER_scores": {"analysis": "scores", "args": {}}}
        key = genewise_required(j, {}, {}, "scripts/", "python")
        assert key is None

    def test_does_not_mutate_json(self):
        j = _make_reagents_json(genewise="")
        j_copy = copy.deepcopy(j)
        genewise_required(j, {}, {}, "scripts/", "python")
        assert j == j_copy
