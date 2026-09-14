"""providers.py — resolves which backend implementation an analysis task uses.

Batch-scale proteome runs can swap a network-backed analysis for a local/bulk
equivalent by setting a mode in batch.config.json (repo root). The interactive
Shiny app never sets that file, so every analysis always resolves to its
"remote" backend — today's behavior, unchanged, by construction rather than
by convention.

Backend modules are flat siblings of this file, named <analysis>_<mode>.py
(e.g. domains_remote.py / domains_bulk.py, conservation_remote.py /
conservation_local.py), each exposing a main(...) with the same signature as
the analysis script's own main(). Import is lazy — only the module for the
resolved mode is ever loaded, so e.g. a machine without DIAMOND installed
never imports conservation_local.py.
"""

import importlib
import json
import sys
from pathlib import Path

_CONFIG_PATH = Path(__file__).parent.parent / "batch.config.json"
_config_cache = None


def _load_config():
    """Load batch.config.json once and cache it; {} if the file doesn't exist."""
    global _config_cache
    if _config_cache is None:
        if _CONFIG_PATH.exists():
            with open(_CONFIG_PATH) as f:
                _config_cache = json.load(f)
        else:
            _config_cache = {}
    return _config_cache


def _reset_config_cache():
    """Test/tooling hook: force the next _load_config() call to re-read the file."""
    global _config_cache
    _config_cache = None


def backend_mode(analysis, default="remote"):
    """Return the configured backend mode for `analysis` (e.g. "domains",
    "conservation", "structure", "genewise"). Falls back to `default` — always
    "remote" unless a caller overrides it — when batch.config.json is absent or
    has no entry for this analysis.
    """
    config = _load_config()
    return config.get("backends", {}).get(analysis, default)


def resolve(analysis, default="remote"):
    """Import <analysis>_<mode>.py for the configured (or default) mode and
    return its main() function.
    """
    mode = backend_mode(analysis, default=default)
    module_name = f"{analysis}_{mode}"
    scripts_dir = str(Path(__file__).parent)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    module = importlib.import_module(module_name)
    return module.main
