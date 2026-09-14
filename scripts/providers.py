"""providers.py — resolves which backend implementation an analysis task uses.

Batch-scale proteome runs can swap a network-backed analysis for a local/bulk
equivalent by setting a mode in a config file — batch.config.json (repo root)
by default. The interactive Shiny app never sets the TAGSITES_BATCH_CONFIG
environment variable and the checked-in batch.config.json ships with every
backend set to "remote", so the app always resolves to its "remote" backend —
today's behavior, unchanged, by construction rather than by convention.

Batch runs (scripts/proteome_run.py) that want local/bulk backends must NOT
just edit batch.config.json in place: that file is read by this same fixed
path regardless of caller, so doing that would silently switch the
interactive app's backends too. Instead, point TAGSITES_BATCH_CONFIG at a
separate file (e.g. a gitignored batch.config.local.json) before running the
batch driver; the app's own config is untouched.

Backend modules are flat siblings of this file, named <analysis>_<mode>.py
(e.g. domains_remote.py / domains_bulk.py, conservation_remote.py /
conservation_local.py), each exposing a main(...) with the same signature as
the analysis script's own main(). Import is lazy — only the module for the
resolved mode is ever loaded, so e.g. a machine without DIAMOND installed
never imports conservation_local.py.
"""

import importlib
import json
import os
import sys
from pathlib import Path

_DEFAULT_CONFIG_PATH = Path(__file__).parent.parent / "batch.config.json"
_config_cache = None
_config_cache_path = None


def _config_path():
    """Path to the active config file: TAGSITES_BATCH_CONFIG if set (relative
    paths resolve against the repo root), else the checked-in batch.config.json.
    """
    override = os.environ.get("TAGSITES_BATCH_CONFIG")
    if not override:
        return _DEFAULT_CONFIG_PATH
    p = Path(override)
    return p if p.is_absolute() else _DEFAULT_CONFIG_PATH.parent / p


def _load_config():
    """Load the active config file once per path and cache it; {} if it doesn't
    exist. Re-reads if TAGSITES_BATCH_CONFIG changes the active path (e.g. a
    test sets the env var mid-process) rather than serving a stale cache.
    """
    global _config_cache, _config_cache_path
    path = _config_path()
    if _config_cache is None or _config_cache_path != path:
        if path.exists():
            with open(path) as f:
                _config_cache = json.load(f)
        else:
            _config_cache = {}
        _config_cache_path = path
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
