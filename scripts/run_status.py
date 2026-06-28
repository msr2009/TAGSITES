"""run_status.py — helpers for per-run task status files.

Status file path: {working_dir}/{run_name}.status.json
Schema per task entry:
    {
        "status":  "pending | running | success | failed",
        "job_id":  "ebi-jobId-or-empty",
        "message": "error message if failed",
        "output":  "/path/to/expected/output/file",
        "stage":   "current high-level step token, e.g. 'blast_submit' (optional)",
        "log":     "full progress log as a newline-joined string (optional)"
    }

A task is "complete" when status == "success" AND its output file exists on disk.
Both the CLI orchestrator (run_tag_sites_from_json.py) and the Shiny app read/write
this file so they share the same resume/skip-completed logic.
"""

import json
import os
import threading

_lock = threading.Lock()


def status_path(working_dir, run_name):
    """Return the canonical path to the status file for a run."""
    return os.path.join(working_dir, f"{run_name}.status.json")


def load_status(working_dir, run_name):
    """Load the status file; return empty dict if the file doesn't exist yet."""
    path = status_path(working_dir, run_name)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def save_status(working_dir, run_name, status_dict):
    """Write the full status dict atomically (tmp → rename) to avoid partial reads."""
    path = status_path(working_dir, run_name)
    os.makedirs(working_dir, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(status_dict, f, indent=2)
    os.replace(tmp, path)


def is_complete(task_entry, output_path=None):
    """True when the task succeeded and its expected output file exists on disk.

    output_path overrides task_entry["output"] if provided.
    """
    op = output_path or task_entry.get("output", "")
    return (
        task_entry.get("status") == "success"
        and bool(op)
        and os.path.exists(op)
    )


def update_task(working_dir, run_name, task_id, **fields):
    """Read-modify-write a single task entry; lock ensures concurrent threads don't clobber each other."""
    with _lock:
        status = load_status(working_dir, run_name)
        if task_id not in status:
            status[task_id] = {}
        status[task_id].update(fields)
        save_status(working_dir, run_name, status)
