"""progress_logic.py — Pure functions for the Progress tab.

No Shiny imports; all functions accept/return plain Python dicts and are
independently testable without running the app.
"""

from config import GLOBAL_KEYS, task_hidden

# effective defaults applied by the runners when a param is blank — shown in
# the progress pane so users see the value that will actually be used
_RUNNER_DEFAULTS = {
    "blast": {"evalue": "1e-10", "max_hits": 100, "db": "uniprotkb", "length": 0, "taxid": "1"},
    "plddt": {"evalue": "1e-10", "percent_id": 99, "taxid": "1"},
}


def parse_run(run_json):
    """Split a run JSON dict into (global_block, [task_entry, ...]).

    Each task entry is a dict with keys: id, analysis, args, output.
    Global fields are merged into each task's args so runners and other
    consumers keep working unchanged.
    """
    global_block = run_json.get("global", {})
    tasks = []
    for k, v in run_json.get("tasks", {}).items():
        # merge global into args at read time — not stored per-task in the file
        merged_args = {**global_block, **v.get("args", {})}
        tasks.append({
            "id":       k,
            "analysis": v.get("type", ""),
            "args":     merged_args,
            "output":   v.get("args", {}).get("output", ""),
        })
    return global_block, tasks


def display_params(args, global_block, analysis=None):
    """Return an ordered dict of per-job params suitable for display.

    Drops: global keys (shared / not job-specific), hidden params for the
    analysis type, and 'output' (internal path). Substitutes runner defaults
    for any param that is blank so the displayed value matches runtime behavior.
    """
    hidden = task_hidden(analysis) if analysis else set()
    skip = GLOBAL_KEYS | hidden | {"output"}
    defaults = _RUNNER_DEFAULTS.get(analysis, {})
    result = {}
    for k, v in args.items():
        if k in skip:
            continue
        # show the effective default when the stored value is blank
        if (v == "" or v is None) and k in defaults:
            v = defaults[k]
        result[k] = v
    return result


def status_label(entry):
    """Return the status string from a status-file entry, defaulting to 'pending'."""
    return entry.get("status", "pending")


def all_terminal(status_dict, tasks):
    """True once every task has reached a terminal status ('success' or 'failed').

    Unlike watching the batch-runner's own extended-task status, this checks
    each task individually, so a batch that "finished" but left some tasks
    'queued' (still waiting on a stalled EBI job) does not count as done.
    """
    if not tasks:
        return False
    return all(status_dict.get(t["id"], {}).get("status") in ("success", "failed")
              for t in tasks)
