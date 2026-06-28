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

