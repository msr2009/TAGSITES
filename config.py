from pathlib import Path
import sys

# add scripts/ to path so task_registry is importable without installing
_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT / "scripts"))

from task_registry import (
    TASK_DEFS, GLOBAL_DEFAULTS, GLOBAL_TOOLTIPS,
    AVAILABLE_TASKS, SELECTABLE_TASKS, GLOBAL_KEYS,
    task_defaults, task_choices, task_tooltips, task_hidden,
    task_script, task_output_suffix, task_companions, companion_path, result_type,
)

# species → NCBI taxid; "Other (search...)" is a sentinel for dynamic search
DEFAULT_SPECIES = {
    "Homo sapiens":             9606,
    "Mus musculus":             10090,
    "Rattus norvegicus":        10116,
    "Danio rerio":              7955,
    "Drosophila melanogaster":  7227,
    "Caenorhabditis elegans":   6239,
    "Saccharomyces cerevisiae": 559292,
    "Escherichia coli":         562,
    "Other (search...)":        None,
}

# derive result-type dicts from the registry instead of maintaining separate lists
RESULTS_TYPE_DICT = {
    "CONTINUOUS": [t for t in TASK_DEFS if TASK_DEFS[t].get("result_type") == "continuous"],
    "RANGE":      [t for t in TASK_DEFS if TASK_DEFS[t].get("result_type") == "range"],
}

# derive per-analysis color palettes from the registry
ANALYSIS_COLORS = {t: d["colors"] for t, d in TASK_DEFS.items() if "colors" in d}

# stable colors for domain/annotation sources (not task-type-keyed, stays in config)
DOMAIN_SOURCE_COLORS = {
    "Phobius":      "#9467bd",   # purple
    "Pfam":         "#17becf",   # teal
    "modification": "#d62728",   # red
}

# isoform classification colors (constitutive → green, intermediate → amber, unique → red)
ISOFORM_CLASS_COLORS = {
    "constitutive": "#00c020",
    "intermediate": "#7bafd4",
    "unique":       "#e81010",
}

# EXCLUDE_ARGS: union of all hidden params across all task types + "output".
# Kept for backward compatibility with any call sites not yet migrated to task_hidden().
EXCLUDE_ARGS = list(
    {"output"} | {p for t in TASK_DEFS for p in task_hidden(t)}
)

# convenience: the full in-memory task parameter map (used by setup_server for now)
# maps task_type -> {"params": {name: {default, choices?, tooltip?, hidden?}}, ...}
TASK_PARAMETERS = TASK_DEFS

# --- backward-compat shims (removed once all modules are migrated) ---

# path to the registry file (was DEFAULT_JSON pointing at default_json.json)
DEFAULT_JSON = str(_ROOT / "task_definitions.json")

# flat global defaults dict (was INPUT_JSON["global"]; setup_server/setup_ui read email default from here)
_GLOBAL_DEFAULTS_WITH_EMPTY = {"email": "", "working_dir": "", "run_name": "",
                                "input_file": "", "pdb": "", "scripts_folder": "./scripts/"}

# INPUT_JSON shim: expose the same dict shape the old modules expected
INPUT_JSON = {
    "global":  _GLOBAL_DEFAULTS_WITH_EMPTY,
    "scripts": {t: task_script(t) for t in TASK_DEFS},
    "global_tooltips": GLOBAL_TOOLTIPS,
    "analyses": {
        t: {
            "name": "",
            "type": t,
            "args": {
                **{p: (cfg.get("choices", [cfg["default"]])[0] if "choices" in cfg else cfg["default"])
                   for p, cfg in d["params"].items()},
                "output": task_output_suffix(t),
            },
            "tooltips": task_tooltips(t),
        }
        for t, d in TASK_DEFS.items()
    },
}
