"""
test_progress_logic.py — offline unit tests for modules/progress_logic.py

Covers: parse_run (global merge, type→analysis, empty cases),
        display_params (key filtering, runner-default substitution),
        status_label.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "modules"))

from progress_logic import parse_run, display_params, status_label


# ── helpers ───────────────────────────────────────────────────────────────────

def _run_json(tasks=None, global_block=None):
    """Build a minimal run-JSON dict."""
    g = global_block or {
        "email": "t@test.com", "run_name": "run1",
        "working_dir": "/wd/", "input_file": "/wd/run1.fa",
        "pdb": "", "scripts_folder": "./scripts/", "selected_sites": [],
    }
    return {"global": g, "tasks": tasks or {}}


# ── parse_run ─────────────────────────────────────────────────────────────────

class TestParseRun:

    def test_returns_tuple_of_global_and_list(self):
        g, tasks = parse_run(_run_json())
        assert isinstance(g, dict)
        assert isinstance(tasks, list)

    def test_global_block_returned_unchanged(self):
        j = _run_json()
        g, _ = parse_run(j)
        assert g == j["global"]

    def test_empty_tasks_gives_empty_list(self):
        _, tasks = parse_run(_run_json())
        assert tasks == []

    def test_task_id_preserved(self):
        j = _run_json(tasks={"WORM_blast": {"type": "blast", "args": {"output": "/o.jsd"}}})
        _, tasks = parse_run(j)
        assert tasks[0]["id"] == "WORM_blast"

    def test_type_becomes_analysis(self):
        j = _run_json(tasks={"X_plddt": {"type": "plddt", "args": {"output": "/o.txt"}}})
        _, tasks = parse_run(j)
        assert tasks[0]["analysis"] == "plddt"

    def test_output_extracted_correctly(self):
        j = _run_json(tasks={"A_scores": {"type": "scores", "args": {"output": "/wd/run1.A_scores.tsv"}}})
        _, tasks = parse_run(j)
        assert tasks[0]["output"] == "/wd/run1.A_scores.tsv"

    def test_global_merged_into_args(self):
        """Global keys should appear in each task's merged args."""
        j = _run_json(tasks={"A_blast": {"type": "blast", "args": {"evalue": "1e-5", "output": ""}}})
        _, tasks = parse_run(j)
        args = tasks[0]["args"]
        assert args["email"] == "t@test.com"
        assert args["run_name"] == "run1"

    def test_task_args_override_global(self):
        """Task-level args with same key as global should win."""
        j = _run_json(
            global_block={"email": "g@g.com", "run_name": "run1",
                           "working_dir": "/wd/", "input_file": "/wd/run1.fa",
                           "pdb": "", "scripts_folder": "./scripts/", "selected_sites": []},
            tasks={"T": {"type": "blast", "args": {"email": "task@task.com", "output": ""}}},
        )
        _, tasks = parse_run(j)
        assert tasks[0]["args"]["email"] == "task@task.com"

    def test_multiple_tasks_returned(self):
        j = _run_json(tasks={
            "A_blast":  {"type": "blast",  "args": {"output": ""}},
            "B_plddt":  {"type": "plddt",  "args": {"output": ""}},
        })
        _, tasks = parse_run(j)
        assert len(tasks) == 2

    def test_missing_global_does_not_crash(self):
        """parse_run should handle a JSON with no global key gracefully."""
        j = {"tasks": {"A": {"type": "blast", "args": {"output": ""}}}}
        g, tasks = parse_run(j)
        assert g == {}
        assert len(tasks) == 1


# ── display_params ────────────────────────────────────────────────────────────

class TestDisplayParams:

    def _args(self, **extra):
        base = {
            "email": "t@t.com", "working_dir": "/wd/", "run_name": "run1",
            "input_file": "/wd/run1.fa", "pdb": "", "scripts_folder": "./",
            "selected_sites": [], "output": "/wd/run1.out",
            "evalue": "1e-5", "max_hits": 100,
        }
        base.update(extra)
        return base

    def test_global_keys_dropped(self):
        args = self._args()
        result = display_params(args, {"email": "t@t.com"}, analysis="blast")
        for k in ("email", "working_dir", "run_name", "input_file",
                   "pdb", "scripts_folder", "selected_sites"):
            assert k not in result, f"global key '{k}' should be dropped"

    def test_output_dropped(self):
        args = self._args()
        result = display_params(args, {}, analysis="blast")
        assert "output" not in result

    def test_task_params_kept(self):
        args = self._args()
        result = display_params(args, {}, analysis="blast")
        assert "evalue" in result
        assert "max_hits" in result

    def test_blank_evalue_substituted_with_runner_default_for_blast(self):
        args = self._args(evalue="")
        result = display_params(args, {}, analysis="blast")
        assert result["evalue"] == "1e-10"

    def test_non_blank_evalue_kept_as_is(self):
        args = self._args(evalue="1e-5")
        result = display_params(args, {}, analysis="blast")
        assert result["evalue"] == "1e-5"

    def test_hidden_plddt_params_dropped(self):
        # pdb and existing_AF2 are hidden for plddt
        args = self._args(pdb="", existing_AF2=1)
        result = display_params(args, {}, analysis="plddt")
        assert "existing_AF2" not in result

    def test_no_analysis_keeps_all_non_global_non_output(self):
        args = self._args()
        result = display_params(args, {}, analysis=None)
        assert "evalue" in result
        assert "output" not in result


# ── status_label ──────────────────────────────────────────────────────────────

class TestStatusLabel:

    def test_returns_status_when_present(self):
        assert status_label({"status": "success"}) == "success"

    def test_defaults_to_pending_when_absent(self):
        assert status_label({}) == "pending"

    def test_other_statuses_passed_through(self):
        for s in ("running", "failed", "pending"):
            assert status_label({"status": s}) == s
