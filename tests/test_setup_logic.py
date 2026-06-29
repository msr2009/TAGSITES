"""
test_setup_logic.py — offline unit tests for modules/setup_logic.py

Covers: make_task, build_global_block, build_task_entry, build_defaults_entry,
        build_reagents_entry, build_run_json.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "modules"))

from setup_logic import (
    make_task,
    build_global_block,
    build_task_entry,
    build_defaults_entry,
    build_reagents_entry,
    build_run_json,
)


# ── make_task ─────────────────────────────────────────────────────────────────

class TestMakeTask:

    def test_id_is_label_underscore_type(self):
        t = make_task("blast", "WORM", {"evalue": "1e-10"}, {})
        assert t["id"] == "WORM_blast"

    def test_name_is_task_type(self):
        t = make_task("plddt", "HUM", {}, {})
        assert t["name"] == "plddt"

    def test_params_copied(self):
        cfg = {"evalue": "1e-5", "max_hits": 50}
        t = make_task("blast", "A", cfg, {})
        assert t["params"] == {"evalue": "1e-5", "max_hits": 50}

    def test_params_independent_of_cfg(self):
        # mutating the returned params should not affect the original cfg
        cfg = {"evalue": "1e-5"}
        t = make_task("blast", "A", cfg, {})
        t["params"]["evalue"] = "1e-100"
        assert cfg["evalue"] == "1e-5"

    def test_tooltips_copied(self):
        tips = {"evalue": "the e-value"}
        t = make_task("blast", "A", {}, tips)
        assert t["tooltips"] == {"evalue": "the e-value"}

    def test_choices_optional_defaults_empty(self):
        t = make_task("blast", "A", {}, {})
        assert t["choices"] == {}

    def test_choices_stored_when_provided(self):
        ch = {"db": ["a", "b"]}
        t = make_task("blast", "A", {}, {}, choices_cfg=ch)
        assert t["choices"] == {"db": ["a", "b"]}


# ── build_global_block ────────────────────────────────────────────────────────

class TestBuildGlobalBlock:

    def _make(self, **kw):
        return build_global_block(
            email="t@test.com", run_name="run1", working_dir="/wd/",
            input_file="/wd/run1.fa", **kw
        )

    def test_required_keys_present(self):
        g = self._make()
        for k in ("email", "run_name", "working_dir", "input_file",
                   "pdb", "scripts_folder", "selected_sites"):
            assert k in g, f"missing key: {k}"

    def test_email_stored(self):
        g = self._make()
        assert g["email"] == "t@test.com"

    def test_pdb_defaults_empty(self):
        g = self._make()
        assert g["pdb"] == ""

    def test_pdb_stored_when_set(self):
        g = self._make(pdb="/wd/run1.pdb")
        assert g["pdb"] == "/wd/run1.pdb"

    def test_selected_sites_defaults_to_empty_list(self):
        g = self._make()
        assert g["selected_sites"] == []

    def test_genomic_file_absent_by_default(self):
        g = self._make()
        assert "genomic_file" not in g

    def test_genomic_file_included_when_passed(self):
        g = self._make(genomic_file="/wd/run1.genomic.fa")
        assert g["genomic_file"] == "/wd/run1.genomic.fa"

    def test_scripts_folder_default(self):
        g = self._make()
        assert g["scripts_folder"] == "./scripts/"


# ── build_task_entry ──────────────────────────────────────────────────────────

class TestBuildTaskEntry:

    def _make(self, task_id="WORM_blast", task_name="blast", args=None):
        collected = args or {"evalue": "1e-10"}
        return build_task_entry(task_id, task_name, collected, ".jsd", "/wd/", "run1")

    def test_top_level_key_is_task_id(self):
        entry = self._make()
        assert "WORM_blast" in entry

    def test_type_field_is_task_name(self):
        entry = self._make()
        assert entry["WORM_blast"]["type"] == "blast"

    def test_args_present(self):
        entry = self._make()
        assert "args" in entry["WORM_blast"]

    def test_output_path_constructed(self):
        entry = self._make()
        assert entry["WORM_blast"]["args"]["output"] == "/wd/run1.WORM_blast.jsd"

    def test_collected_args_preserved(self):
        entry = self._make(args={"evalue": "1e-5", "max_hits": 100})
        assert entry["WORM_blast"]["args"]["evalue"] == "1e-5"
        assert entry["WORM_blast"]["args"]["max_hits"] == 100

    def test_no_global_keys_in_args(self):
        # global fields should not bleed into task args
        global_keys = {"email", "working_dir", "run_name", "input_file",
                        "pdb", "scripts_folder", "genomic_file", "selected_sites"}
        entry = self._make()
        for k in global_keys:
            assert k not in entry["WORM_blast"]["args"], f"global key '{k}' leaked into args"


# ── build_defaults_entry ──────────────────────────────────────────────────────

class TestBuildDefaultsEntry:

    def test_shape_matches_task_entry(self):
        """build_defaults_entry and build_task_entry must produce identical structure."""
        args = {"evalue": "1e-10"}
        task  = build_task_entry("A_blast", "blast", args, ".jsd", "/wd/", "run")
        dflt  = build_defaults_entry("A_blast", "blast", args, ".jsd", "/wd/", "run")
        assert task == dflt

    def test_output_path_constructed(self):
        entry = build_defaults_entry("X_plddt", "plddt", {}, ".txt", "/d/", "r")
        assert entry["X_plddt"]["args"]["output"] == "/d/r.X_plddt.txt"


# ── build_reagents_entry ──────────────────────────────────────────────────────

class TestBuildReagentsEntry:

    def _make(self):
        return build_reagents_entry(
            defaults={"n_guides": 5, "arm_length": 1000},
            genomic_path="/wd/run1.genomic.fa",
            out_suffix=".tsv",
            working_dir="/wd/",
            run_name="run1",
        )

    def test_type_is_reagents(self):
        entry = self._make()
        assert entry["type"] == "reagents"

    def test_genewise_starts_empty(self):
        entry = self._make()
        assert entry["args"]["genewise"] == ""

    def test_genomic_fasta_stored(self):
        entry = self._make()
        assert entry["args"]["genomic_fasta"] == "/wd/run1.genomic.fa"

    def test_output_constructed(self):
        entry = self._make()
        assert entry["args"]["output"] == "/wd/run1.reagents.tsv"

    def test_defaults_preserved(self):
        entry = self._make()
        assert entry["args"]["n_guides"] == 5
        assert entry["args"]["arm_length"] == 1000


# ── build_run_json ────────────────────────────────────────────────────────────

class TestBuildRunJson:

    def _global(self):
        return build_global_block(
            email="t@test.com", run_name="run1",
            working_dir="/wd/", input_file="/wd/run1.fa",
        )

    def test_top_level_keys(self):
        j = build_run_json(self._global(), [])
        assert set(j.keys()) == {"global", "tasks"}

    def test_global_block_stored(self):
        g = self._global()
        j = build_run_json(g, [])
        assert j["global"] is g

    def test_task_entries_merged_into_tasks(self):
        entry1 = build_task_entry("A_blast", "blast", {}, ".jsd", "/wd/", "run1")
        entry2 = build_task_entry("B_plddt", "plddt", {}, ".txt", "/wd/", "run1")
        j = build_run_json(self._global(), [entry1, entry2])
        assert "A_blast" in j["tasks"]
        assert "B_plddt" in j["tasks"]

    def test_reagents_absent_when_none(self):
        j = build_run_json(self._global(), [])
        assert "REAGENTS_reagents" not in j["tasks"]

    def test_reagents_included_under_fixed_key(self):
        reagents = build_reagents_entry({}, "/wd/r.fa", ".tsv", "/wd/", "run1")
        j = build_run_json(self._global(), [], reagents_entry=reagents)
        assert "REAGENTS_reagents" in j["tasks"]
        assert j["tasks"]["REAGENTS_reagents"] is reagents
