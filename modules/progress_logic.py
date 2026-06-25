"""progress_logic.py — Pure functions for the Progress tab.

No Shiny imports; all functions accept/return plain Python dicts and are
independently testable without running the app.
"""

from config import EXCLUDE_ARGS

# keys that are part of the global block and not meaningful per-job params
_GLOBAL_KEYS = {"email", "working_dir", "run_name", "input_file", "pdb",
                "scripts_folder", "genomic_file"}

# effective defaults applied by the runners when a param is blank — shown in
# the progress pane so users see the value that will actually be used
_RUNNER_DEFAULTS = {
    "blast": {"evalue": "1e-10", "max_hits": 100, "db": "uniprotkb", "length": 0, "taxid": "1"},
    "plddt": {"evalue": "1e-10", "percent_id": 99, "taxid": "1"},
}


def parse_run(run_json):
    """Split a run JSON dict into (global_block, [task_entry, ...]).

    Each task entry is a dict with keys: id, analysis, args, output.
    The 'scripts' and 'global' top-level keys are consumed; everything else
    is treated as a task.
    """
    global_block = run_json.get("global", {})
    tasks = []
    for k, v in run_json.items():
        if k in ("scripts", "global"):
            continue
        tasks.append({
            "id":       k,
            "analysis": v.get("analysis", ""),
            "args":     v.get("args", {}),
            "output":   v.get("args", {}).get("output", ""),
        })
    return global_block, tasks


def display_params(args, global_block, analysis=None):
    """Return an ordered dict of per-job params suitable for display.

    Drops: keys in the global block (shared / not job-specific), keys in
    EXCLUDE_ARGS (internal pipeline flags), and 'output' (internal path).
    Substitutes runner defaults for any param that is blank so the displayed
    value matches what will actually be used at runtime.
    """
    skip = set(global_block.keys()) | set(EXCLUDE_ARGS) | {"output"}
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
