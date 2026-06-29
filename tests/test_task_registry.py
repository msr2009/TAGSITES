"""
test_task_registry.py — offline unit tests for scripts/task_registry.py

Covers: task_defaults, task_choices, task_tooltips, task_hidden,
        task_script, task_output_suffix, task_companions,
        companion_path, result_type, GLOBAL_KEYS, AVAILABLE_TASKS.
All tests run against the real task_definitions.json.
"""

import pytest

from task_registry import (
    task_defaults,
    task_choices,
    task_tooltips,
    task_hidden,
    task_script,
    task_output_suffix,
    task_companions,
    companion_path,
    result_type,
    GLOBAL_KEYS,
    AVAILABLE_TASKS,
    SELECTABLE_TASKS,
)


# ── GLOBAL_KEYS ───────────────────────────────────────────────────────────────

def test_global_keys_is_set():
    assert isinstance(GLOBAL_KEYS, set)


def test_global_keys_contains_canonical_fields():
    for k in ("email", "working_dir", "run_name", "input_file",
               "pdb", "scripts_folder", "selected_sites"):
        assert k in GLOBAL_KEYS, f"expected '{k}' in GLOBAL_KEYS"


# ── AVAILABLE_TASKS / SELECTABLE_TASKS ───────────────────────────────────────

def test_available_tasks_includes_core_types():
    for t in ("blast", "plddt", "modifications", "domains", "scores", "reagents"):
        assert t in AVAILABLE_TASKS, f"'{t}' missing from AVAILABLE_TASKS"


def test_reagents_not_selectable():
    # reagents is auto-added; users should not be able to select it manually
    assert "reagents" not in SELECTABLE_TASKS


def test_core_analyses_selectable():
    for t in ("blast", "plddt", "modifications", "domains", "scores"):
        assert t in SELECTABLE_TASKS, f"'{t}' should be selectable"


# ── task_defaults ──────────────────────────────────────────────────────────────

class TestTaskDefaults:

    def test_returns_dict(self):
        assert isinstance(task_defaults("blast"), dict)

    def test_blast_has_evalue(self):
        d = task_defaults("blast")
        assert "evalue" in d

    def test_plddt_existing_af2_default_is_one(self):
        d = task_defaults("plddt")
        assert d["existing_AF2"] == 1

    def test_plddt_pdb_default_is_empty(self):
        d = task_defaults("plddt")
        assert d["pdb"] == ""


# ── task_choices ──────────────────────────────────────────────────────────────

class TestTaskChoices:

    def test_returns_dict(self):
        assert isinstance(task_choices("blast"), dict)

    def test_blast_evalue_has_choices(self):
        ch = task_choices("blast")
        assert "evalue" in ch
        assert isinstance(ch["evalue"], list)
        assert len(ch["evalue"]) > 0

    def test_blast_db_has_choices(self):
        ch = task_choices("blast")
        assert "db" in ch
        assert "uniprotkb" in ch["db"]

    def test_params_without_choices_absent(self):
        # blast's taxid has no choices dropdown
        ch = task_choices("blast")
        assert "taxid" not in ch


# ── task_tooltips ─────────────────────────────────────────────────────────────

class TestTaskTooltips:

    def test_returns_dict(self):
        assert isinstance(task_tooltips("blast"), dict)

    def test_blast_evalue_has_tooltip(self):
        tips = task_tooltips("blast")
        assert "evalue" in tips
        assert isinstance(tips["evalue"], str) and len(tips["evalue"]) > 0

    def test_params_without_tooltip_absent(self):
        # modifications has sites_file with a tooltip; domains has no params → empty dict
        tips = task_tooltips("domains")
        assert isinstance(tips, dict)


# ── task_hidden ───────────────────────────────────────────────────────────────

class TestTaskHidden:

    def test_returns_set(self):
        assert isinstance(task_hidden("plddt"), set)

    def test_plddt_pdb_and_existing_af2_hidden(self):
        h = task_hidden("plddt")
        assert "pdb" in h
        assert "existing_AF2" in h

    def test_reagents_genewise_hidden(self):
        h = task_hidden("reagents")
        assert "genewise" in h

    def test_blast_no_hidden_params(self):
        h = task_hidden("blast")
        assert len(h) == 0


# ── task_script / task_output_suffix ─────────────────────────────────────────

class TestTaskScriptAndSuffix:

    def test_blast_script(self):
        assert task_script("blast") == "blast_orthologs.py"

    def test_plddt_script(self):
        assert task_script("plddt") == "extract_from_pdb.py"

    def test_blast_output_suffix(self):
        assert task_output_suffix("blast") == ".jsd"

    def test_plddt_output_suffix(self):
        assert task_output_suffix("plddt") == ".txt"

    def test_reagents_output_suffix(self):
        assert task_output_suffix("reagents") == ".tsv"


# ── task_companions ───────────────────────────────────────────────────────────

class TestTaskCompanions:

    def test_blast_has_alignment_companion(self):
        c = task_companions("blast")
        assert "alignment" in c

    def test_plddt_has_sasa_companion(self):
        c = task_companions("plddt")
        assert "sasa" in c

    def test_modifications_has_no_companions(self):
        c = task_companions("modifications")
        assert c == {}


# ── companion_path ─────────────────────────────────────────────────────────────

class TestCompanionPath:

    def test_plddt_sasa_path(self):
        out = "/wd/run1.PLDDT_plddt.txt"
        result = companion_path(out, "plddt", "sasa")
        assert result == "/wd/run1.PLDDT_plddt.sasa.txt"

    def test_blast_alignment_path(self):
        out = "/wd/run1.WORM_blast.jsd"
        result = companion_path(out, "blast", "alignment")
        assert result == "/wd/run1.WORM_blast.aln"

    def test_unknown_companion_key_returns_empty(self):
        out = "/wd/run1.WORM_blast.jsd"
        result = companion_path(out, "blast", "nonexistent_key")
        assert result == ""

    def test_wrong_suffix_returns_empty(self):
        # output path doesn't end in the task's own suffix
        out = "/wd/run1.wrong_suffix.txt"
        result = companion_path(out, "blast", "alignment")
        assert result == ""

    def test_task_with_no_companions_returns_empty(self):
        out = "/wd/run1.DOMAINS_domains.txt"
        result = companion_path(out, "domains", "anything")
        assert result == ""


# ── result_type ────────────────────────────────────────────────────────────────

class TestResultType:

    def test_blast_is_continuous(self):
        assert result_type("blast") == "continuous"

    def test_plddt_is_continuous(self):
        assert result_type("plddt") == "continuous"

    def test_modifications_is_range(self):
        assert result_type("modifications") == "range"

    def test_domains_is_range(self):
        assert result_type("domains") == "range"

    def test_reagents_is_none(self):
        assert result_type("reagents") == "none"
