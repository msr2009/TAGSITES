"""task_registry.py — load task_definitions.json and expose accessor helpers.

No Shiny imports; safe to use from both the app (via config.py) and standalone
CLI scripts (run_tag_sites_from_json.py, task_runners.py).
"""
import json
from pathlib import Path

# task_definitions.json lives one level above scripts/
_ROOT = Path(__file__).parent.parent
REGISTRY_PATH = str(_ROOT / "task_definitions.json")

try:
    _REG = json.load(open(REGISTRY_PATH))
except FileNotFoundError:
    raise RuntimeError(f"Could not find task_definitions.json at {REGISTRY_PATH}")

# top-level sections
TASK_DEFS         = _REG["tasks"]
GLOBAL_DEFAULTS   = _REG.get("global_defaults", {})
GLOBAL_TOOLTIPS   = _REG.get("global_tooltips", {})

# derived lists
AVAILABLE_TASKS   = list(TASK_DEFS)
SELECTABLE_TASKS  = [t for t, d in TASK_DEFS.items() if d.get("selectable", True)]

# canonical set of keys that belong in the global block, not per-task args
GLOBAL_KEYS = {
    "email", "working_dir", "run_name", "input_file",
    "pdb", "scripts_folder", "genomic_file", "selected_sites",
}

# --- per-type accessors ---

def task_defaults(task_type):
    """Return {param: default_value} for all params of task_type."""
    return {p: cfg["default"] for p, cfg in TASK_DEFS[task_type]["params"].items()}


def task_choices(task_type):
    """Return {param: [choices]} for params that have a dropdown; others absent."""
    return {
        p: cfg["choices"]
        for p, cfg in TASK_DEFS[task_type]["params"].items()
        if "choices" in cfg
    }


def task_tooltips(task_type):
    """Return {param: tooltip_string} for params that have a tooltip."""
    return {
        p: cfg["tooltip"]
        for p, cfg in TASK_DEFS[task_type]["params"].items()
        if cfg.get("tooltip")
    }


def task_hidden(task_type):
    """Return set of param names that should not appear in the UI form."""
    return {p for p, cfg in TASK_DEFS[task_type]["params"].items() if cfg.get("hidden")}


def task_script(task_type):
    """Return the script filename (basename only) for task_type."""
    return TASK_DEFS[task_type]["script"]


def task_output_suffix(task_type):
    """Return the output file suffix (e.g. '.jsd') for task_type."""
    return TASK_DEFS[task_type]["output_suffix"]


def task_companions(task_type):
    """Return {companion_key: suffix} for sidecar output files, or {}."""
    return TASK_DEFS[task_type].get("companions", {})


def companion_path(output_path, task_type, companion_key):
    """Derive a companion file path by swapping the task's output_suffix for companion_key's suffix."""
    base_suffix   = task_output_suffix(task_type)
    comp_suffix   = task_companions(task_type).get(companion_key, "")
    if not comp_suffix or not output_path.endswith(base_suffix):
        return ""
    return output_path[: -len(base_suffix)] + comp_suffix


def result_type(task_type):
    """Return 'continuous', 'range', or 'none' for task_type."""
    return TASK_DEFS[task_type].get("result_type", "none")
