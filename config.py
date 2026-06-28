from pathlib import Path
import json

# project root is the directory containing this file
_ROOT = Path(__file__).parent

# edit the filename below for your defaults
DEFAULT_JSON = str(_ROOT / "default_json.json")

# you can change these if you'd like (except for Other)
DEFAULT_SPECIES = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Rattus norvegicus": 10116,
    "Danio rerio": 7955,
    "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239,
    "Saccharomyces cerevisiae": 559292,
    "Escherichia coli": 562,
    "Other (search...)": None  # marker for dynamic search
}

# don't touch this!
try:
    INPUT_JSON = json.load(open(DEFAULT_JSON, "r"))
except FileNotFoundError:
    raise RuntimeError(f"Could not find default_json.json at {DEFAULT_JSON}")
TASK_PARAMETERS = INPUT_JSON["analyses"]
AVAILABLE_TASKS = list(INPUT_JSON["analyses"].keys())
# reagents is auto-injected at save time when a genomic FASTA is present; not user-selectable
SELECTABLE_TASKS = [t for t in AVAILABLE_TASKS if t != "reagents"]
EXCLUDE_ARGS = ["fasta", "input", "output", "pdb"]
GLOBAL_TOOLTIPS = INPUT_JSON.get("global_tooltips", {})

# analysis result types — add new analysis types here
RESULTS_TYPE_DICT = {
    "CONTINUOUS": ["blast", "plddt", "scores"],
    "RANGE": ["domains", "modifications", "isoforms"]
}

# three-color palettes per analysis type — cycled when multiple tracks of same type exist
# colors are drawn from the Tableau-10 palette so they are maximally distinct
ANALYSIS_COLORS = {
    "blast":  ["#2e7d32", "#00838f", "#8bc34a"],   # deep green, teal, lime
    "plddt":  ["#1565c0", "#6a1b9a", "#ad1457"],   # deep blue, deep purple, deep pink
    "scores": ["#e65100", "#b71c1c", "#795548"],   # deep orange, deep red, brown
}

# stable colors for each domain/annotation source (plot boxes + 3D domain coloring)
DOMAIN_SOURCE_COLORS = {
    "Phobius":      "#9467bd",   # purple
    "Pfam":         "#17becf",   # teal
    "modification": "#d62728",   # red
}

# isoform classification colors (constitutive → green, intermediate → amber, unique → red)
ISOFORM_CLASS_COLORS = {
    "constitutive":  "#00c020",   # bright green
    "intermediate":  "#7bafd4",   # muted blue
    "unique":        "#e81010",   # bright red
}
