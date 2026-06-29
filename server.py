import io
import zipfile
from pathlib import Path

from shiny import reactive, render

from modules.setup_server import setup_server
from modules.progress_server import progress_server
from modules.results_server import results_server
from modules.reagents_server import reagents_server

_ROOT = Path(__file__).parent


def app_server(input, output, session):
    # per-session shared state: the path to the active run JSON
    # created inside app_server so each browser session gets its own value
    shared_values = reactive.Value("")

    # chosen tag sites — list of residue positions (ints) from the Results page
    shared_sites = reactive.Value([])

    # incremented by progress_server each time a run completes so results_server
    # reloads even when the JSON path hasn't changed
    shared_results_trigger = reactive.Value(0)

    setup_server("setup", shared_json=shared_values)
    progress_server("progress", shared_json=shared_values,
                    shared_results_trigger=shared_results_trigger)
    results_server("results", shared_json=shared_values, shared_sites=shared_sites,
                   shared_results_trigger=shared_results_trigger)
    reagents_server("reagents", shared_json=shared_values, shared_sites=shared_sites)

    @render.download(filename="tagsites-app.zip", media_type="application/zip")
    async def download_app():
        """Bundle app source files into a zip for local deployment."""
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for arcname, filepath in _collect_app_files(_ROOT):
                zf.write(filepath, arcname)
        yield buf.getvalue()


def _collect_app_files(root):
    """Yield (arcname, filepath) pairs for app source files to include in the download zip."""
    for pattern in ("*.py", "*.json", "*.yml", "*.yaml", "*.txt"):
        for f in sorted(root.glob(pattern)):
            yield f.name, f
    readme = root / "README.md"
    if readme.exists():
        yield "README.md", readme
    for subdir in ("modules", "scripts", "utils", "tables", "www", "params"):
        d = root / subdir
        if not d.is_dir():
            continue
        for f in sorted(d.rglob("*")):
            if f.is_file() and "__pycache__" not in str(f) and not f.name.startswith("."):
                yield str(f.relative_to(root)), f
