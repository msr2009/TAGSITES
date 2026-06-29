"""
test_json_roundtrip.py — end-to-end contract test for the JSON rewrite.

Builds run-JSON using setup_logic producers, then reads it back with
progress_logic.parse_run and utils/results.load_run_metadata.
Verifies that global keys stay in the global block (not duplicated per task),
that type/args/output survive the trip, and that global fields are merged into
each task's args by the consumer.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "modules"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from setup_logic import (
    build_global_block,
    build_task_entry,
    build_reagents_entry,
    build_run_json,
)
from progress_logic import parse_run
from utils.results import load_run_metadata


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def sample_run_json(tmp_path):
    """Produce a realistic run-JSON dict using the official builder functions."""
    fa = tmp_path / "run1.fa"
    fa.write_text(">run1\nMSTGKKVL\n")

    g = build_global_block(
        email="t@test.com",
        run_name="run1",
        working_dir=str(tmp_path) + "/",
        input_file=str(fa),
        pdb="",
        genomic_file="",
    )
    blast_entry  = build_task_entry("WORM_blast",  "blast",  {"evalue": "1e-10"}, ".jsd", str(tmp_path) + "/", "run1")
    plddt_entry  = build_task_entry("HUM_plddt",   "plddt",  {"existing_AF2": 1}, ".txt", str(tmp_path) + "/", "run1")
    scores_entry = build_task_entry("ANY_scores",  "scores", {"window": 5},       ".tsv", str(tmp_path) + "/", "run1")
    return build_run_json(g, [blast_entry, plddt_entry, scores_entry])


# ── JSON shape tests ──────────────────────────────────────────────────────────

class TestRunJsonShape:

    def test_top_level_keys(self, sample_run_json):
        assert set(sample_run_json.keys()) == {"global", "tasks"}

    def test_global_contains_required_fields(self, sample_run_json):
        g = sample_run_json["global"]
        for k in ("email", "run_name", "working_dir", "input_file", "pdb", "scripts_folder"):
            assert k in g, f"missing global key: {k}"

    def test_tasks_contains_expected_ids(self, sample_run_json):
        tasks = sample_run_json["tasks"]
        assert "WORM_blast" in tasks
        assert "HUM_plddt" in tasks
        assert "ANY_scores" in tasks

    def test_each_task_has_type_and_args(self, sample_run_json):
        for tid, v in sample_run_json["tasks"].items():
            assert "type" in v, f"task {tid} missing 'type'"
            assert "args" in v, f"task {tid} missing 'args'"

    def test_type_values_correct(self, sample_run_json):
        tasks = sample_run_json["tasks"]
        assert tasks["WORM_blast"]["type"] == "blast"
        assert tasks["HUM_plddt"]["type"] == "plddt"
        assert tasks["ANY_scores"]["type"] == "scores"

    def test_global_keys_not_duplicated_per_task(self, sample_run_json):
        """Global fields must live in 'global' only, not repeated in each task's args."""
        global_keys = {"email", "working_dir", "run_name", "input_file",
                        "pdb", "scripts_folder", "selected_sites"}
        for tid, v in sample_run_json["tasks"].items():
            for k in global_keys:
                assert k not in v["args"], (
                    f"global key '{k}' leaked into task '{tid}' args"
                )

    def test_output_path_in_args(self, sample_run_json):
        assert "output" in sample_run_json["tasks"]["WORM_blast"]["args"]

    def test_json_serializable(self, sample_run_json):
        """The dict must round-trip through json.dumps/loads without error."""
        s = json.dumps(sample_run_json)
        restored = json.loads(s)
        assert set(restored.keys()) == {"global", "tasks"}


# ── parse_run round-trip tests ────────────────────────────────────────────────

class TestParseRunRoundTrip:

    def test_global_block_returned(self, sample_run_json):
        g, _ = parse_run(sample_run_json)
        assert g["email"] == "t@test.com"
        assert g["run_name"] == "run1"

    def test_task_list_length(self, sample_run_json):
        _, tasks = parse_run(sample_run_json)
        assert len(tasks) == 3

    def test_task_analysis_matches_type(self, sample_run_json):
        _, tasks = parse_run(sample_run_json)
        by_id = {t["id"]: t for t in tasks}
        assert by_id["WORM_blast"]["analysis"] == "blast"
        assert by_id["HUM_plddt"]["analysis"] == "plddt"

    def test_global_merged_into_task_args(self, sample_run_json):
        """After parse_run each task's args must contain global fields."""
        _, tasks = parse_run(sample_run_json)
        for t in tasks:
            assert "email" in t["args"], f"email missing from task '{t['id']}' args"
            assert "run_name" in t["args"], f"run_name missing from task '{t['id']}' args"

    def test_task_level_arg_preserved(self, sample_run_json):
        _, tasks = parse_run(sample_run_json)
        blast = next(t for t in tasks if t["id"] == "WORM_blast")
        assert blast["args"]["evalue"] == "1e-10"

    def test_output_in_task_entry(self, sample_run_json):
        _, tasks = parse_run(sample_run_json)
        for t in tasks:
            assert t["output"] != "", f"task '{t['id']}' has empty output"

    def test_roundtrip_preserves_task_ids(self, sample_run_json):
        _, tasks = parse_run(sample_run_json)
        ids = {t["id"] for t in tasks}
        assert ids == {"WORM_blast", "HUM_plddt", "ANY_scores"}


# ── load_run_metadata round-trip ──────────────────────────────────────────────

class TestLoadRunMetadataRoundTrip:

    def test_query_seq_read_from_fasta(self, sample_run_json):
        meta = load_run_metadata(sample_run_json)
        assert meta["query_seq"] == "MSTGKKVL"

    def test_seq_len_matches(self, sample_run_json):
        meta = load_run_metadata(sample_run_json)
        assert meta["seq_len"] == len("MSTGKKVL")

    def test_pdb_path_empty_when_none(self, sample_run_json):
        meta = load_run_metadata(sample_run_json)
        assert meta["pdb_path"] == ""

    def test_accepts_json_serialized_form(self, sample_run_json):
        """Metadata loading must work when given the JSON-serialized form."""
        meta = load_run_metadata(json.dumps(sample_run_json))
        assert meta["query_seq"] == "MSTGKKVL"


# ── reagents integration ──────────────────────────────────────────────────────

class TestReagentsIntegration:

    def test_reagents_entry_included_in_tasks(self, tmp_path):
        fa = tmp_path / "run1.fa"
        fa.write_text(">run1\nMSTGK\n")
        g = build_global_block("t@t.com", "run1", str(tmp_path) + "/", str(fa))
        r = build_reagents_entry({"n_guides": 5}, str(tmp_path / "run1.genomic.fa"), ".tsv",
                                  str(tmp_path) + "/", "run1")
        j = build_run_json(g, [], reagents_entry=r)
        assert "REAGENTS_reagents" in j["tasks"]
        assert j["tasks"]["REAGENTS_reagents"]["type"] == "reagents"

    def test_reagents_genewise_empty_survives_roundtrip(self, tmp_path):
        fa = tmp_path / "run1.fa"
        fa.write_text(">run1\nMSTGK\n")
        g = build_global_block("t@t.com", "run1", str(tmp_path) + "/", str(fa))
        r = build_reagents_entry({}, str(tmp_path / "gen.fa"), ".tsv", str(tmp_path) + "/", "run1")
        j = build_run_json(g, [], reagents_entry=r)
        _, tasks = parse_run(j)
        reagents = next(t for t in tasks if t["id"] == "REAGENTS_reagents")
        assert reagents["args"]["genewise"] == ""
