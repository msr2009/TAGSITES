"""
test_run_status.py — offline unit tests for scripts/run_status.py

Covers: status_path, load_status, save_status, update_task, is_complete.
All filesystem operations use pytest's tmp_path fixture.
"""

import json
import os

import pytest

from run_status import (
    status_path,
    load_status,
    save_status,
    update_task,
    is_complete,
)


# ── status_path ───────────────────────────────────────────────────────────────

class TestStatusPath:

    def test_combines_dir_and_name(self, tmp_path):
        p = status_path(str(tmp_path), "myrun")
        assert p == str(tmp_path / "myrun.status.json")

    def test_uses_os_sep(self, tmp_path):
        p = status_path(str(tmp_path), "myrun")
        assert os.path.basename(p) == "myrun.status.json"


# ── load_status ───────────────────────────────────────────────────────────────

class TestLoadStatus:

    def test_returns_empty_dict_when_file_missing(self, tmp_path):
        result = load_status(str(tmp_path), "nofile")
        assert result == {}

    def test_returns_dict_when_file_exists(self, tmp_path):
        data = {"WORM_blast": {"status": "success", "output": "/out.jsd"}}
        path = tmp_path / "run1.status.json"
        path.write_text(json.dumps(data))
        result = load_status(str(tmp_path), "run1")
        assert result == data

    def test_preserves_all_fields(self, tmp_path):
        data = {"T": {"status": "running", "job_id": "abc123", "log": "line1\nline2"}}
        (tmp_path / "r.status.json").write_text(json.dumps(data))
        result = load_status(str(tmp_path), "r")
        assert result["T"]["job_id"] == "abc123"
        assert result["T"]["log"] == "line1\nline2"


# ── save_status ───────────────────────────────────────────────────────────────

class TestSaveStatus:

    def test_creates_file(self, tmp_path):
        save_status(str(tmp_path), "run1", {"A": {"status": "pending"}})
        assert (tmp_path / "run1.status.json").exists()

    def test_round_trips_data(self, tmp_path):
        data = {"A": {"status": "success", "output": "/a.jsd"}}
        save_status(str(tmp_path), "run1", data)
        result = load_status(str(tmp_path), "run1")
        assert result == data

    def test_creates_working_dir_if_absent(self, tmp_path):
        wd = str(tmp_path / "nested" / "dir")
        save_status(wd, "run1", {})
        assert os.path.exists(os.path.join(wd, "run1.status.json"))

    def test_overwrites_previous_content(self, tmp_path):
        save_status(str(tmp_path), "run1", {"A": {"status": "pending"}})
        save_status(str(tmp_path), "run1", {"A": {"status": "success"}})
        result = load_status(str(tmp_path), "run1")
        assert result["A"]["status"] == "success"


# ── update_task ───────────────────────────────────────────────────────────────

class TestUpdateTask:

    def test_creates_task_entry_if_absent(self, tmp_path):
        update_task(str(tmp_path), "run1", "A_blast", status="running")
        result = load_status(str(tmp_path), "run1")
        assert "A_blast" in result
        assert result["A_blast"]["status"] == "running"

    def test_merges_fields_into_existing_entry(self, tmp_path):
        save_status(str(tmp_path), "run1", {"A": {"status": "pending", "job_id": "x"}})
        update_task(str(tmp_path), "run1", "A", status="running", log="started")
        result = load_status(str(tmp_path), "run1")
        assert result["A"]["status"] == "running"
        assert result["A"]["job_id"] == "x"      # untouched field preserved
        assert result["A"]["log"] == "started"

    def test_does_not_clobber_other_tasks(self, tmp_path):
        save_status(str(tmp_path), "run1", {
            "A": {"status": "success"},
            "B": {"status": "pending"},
        })
        update_task(str(tmp_path), "run1", "A", status="failed")
        result = load_status(str(tmp_path), "run1")
        assert result["B"]["status"] == "pending"


# ── is_complete ───────────────────────────────────────────────────────────────

class TestIsComplete:

    def test_true_when_success_and_output_exists(self, tmp_path):
        out = tmp_path / "out.jsd"
        out.write_text("data")
        entry = {"status": "success", "output": str(out)}
        assert is_complete(entry) is True

    def test_false_when_output_missing(self, tmp_path):
        entry = {"status": "success", "output": str(tmp_path / "missing.jsd")}
        assert is_complete(entry) is False

    def test_false_when_status_not_success(self, tmp_path):
        out = tmp_path / "out.jsd"
        out.write_text("data")
        entry = {"status": "running", "output": str(out)}
        assert is_complete(entry) is False

    def test_false_when_output_empty_string(self):
        entry = {"status": "success", "output": ""}
        assert is_complete(entry) is False

    def test_output_path_arg_overrides_entry(self, tmp_path):
        out = tmp_path / "real.jsd"
        out.write_text("data")
        entry = {"status": "success", "output": str(tmp_path / "fake.jsd")}
        assert is_complete(entry, output_path=str(out)) is True

    def test_false_when_entry_empty(self):
        assert is_complete({}) is False
