"""TAGSITES batch: review tag sites for many proteins at once and output their reagents.

Local only (needs the proteome database). Runs with batch.config.worm_local.json unless
TAGSITES_BATCH_CONFIG names another config; port 8001, so it can run beside app_modular.py.

    conda activate tagsites
    python app_batch.py
"""

import os
import sys
from pathlib import Path

# the config must be chosen before any module reads it at import time
os.environ.setdefault("TAGSITES_BATCH_CONFIG", "batch.config.worm_local.json")

from shiny import App  # noqa: E402

# same module name the batch modules import, so they share one batch-file lock
sys.path.insert(0, str(Path(__file__).parent / "scripts"))
from batch_state import app_config  # noqa: E402
from utils.results import enable_load_cache  # noqa: E402

# protein switches in Review read parsed results from memory (filled by Prepare)
enable_load_cache(int(app_config().get("cache_entries", 500)))

from server_batch import app_server  # noqa: E402
from ui_batch import app_ui  # noqa: E402

app = App(app_ui, app_server, static_assets=Path(__file__).parent / "www")

if __name__ == "__main__":
    app.run(port=int(app_config().get("port", 8001)))
