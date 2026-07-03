import asyncio
import io
import json
import os
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

from shiny import reactive, render, ui

from modules.setup_server import setup_server
from modules.progress_server import progress_server
from modules.results_server import results_server
from modules.reagents_server import reagents_server
from utils.bundle import is_bundle_zip, extract_and_rebase

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

    # session-scoped temp dir for extracted bundle files; created lazily on first
    # bundle upload and deleted when the browser session ends
    _session_temp = [""]

    def _get_session_dir():
        """Return (creating if needed) the session-scoped temp dir."""
        if not _session_temp[0]:
            _session_temp[0] = tempfile.mkdtemp(prefix="tagsites_")
        return _session_temp[0]

    @session.on_ended
    def _cleanup_session():
        """Remove the session temp dir when the browser disconnects."""
        d = _session_temp[0]
        if d and os.path.isdir(d):
            shutil.rmtree(d, ignore_errors=True)

    @reactive.effect
    @reactive.event(shared_values)
    async def _handle_bundle_upload():
        """Intercept bundle ZIP uploads: extract, rebase, and re-run reagents if missing.

        When any tab sets shared_values to a bundle ZIP path, this effect:
        1. Extracts the bundle into a session-scoped temp dir and rewrites all paths.
        2. If the reagents TSV is absent but genewise output files are present, patches
           the reagents args to point at them so the runner skips the slow EBI step.
        3. Runs the reagents task in a thread if the TSV still doesn't exist.
        Non-ZIP uploads (bare run JSON files) pass through unchanged.
        The loop cannot repeat: the rewritten JSON path is not a ZIP.
        """
        path = shared_values.get()
        if not path or not is_bundle_zip(path):
            return
        try:
            new_path = extract_and_rebase(path, _get_session_dir())
        except Exception as exc:
            ui.notification_show(f"Could not extract bundle: {exc}",
                                 type="error", duration=8)
            return

        # check if the reagents TSV was included; if not, re-run the task
        try:
            j = json.load(open(new_path))
            reagents_task = j.get("tasks", {}).get("REAGENTS_reagents")
            if reagents_task:
                args    = reagents_task.get("args", {})
                tsv_out = args.get("output", "")
                if tsv_out and not os.path.exists(tsv_out):
                    # patch genewise paths if available so the runner skips the EBI call
                    wd = j["global"].get("working_dir", "")
                    rn = j["global"].get("run_name", "")
                    gw_out = os.path.join(wd, f"{rn}_genewise.genewise.out.txt")
                    gw_fa  = os.path.join(wd, f"{rn}_genewise.genewise_genomic.fa")
                    if os.path.exists(gw_out) and os.path.exists(gw_fa):
                        args = dict(args, genewise=gw_out, genomic_fasta=gw_fa)

                    ui.notification_show("Reagents TSV not in bundle — recomputing…",
                                         type="message", duration=5)
                    sys.path.insert(0, str(Path(__file__).parent / "scripts"))
                    from task_runners import run_reagents
                    loop = asyncio.get_running_loop()
                    await loop.run_in_executor(None, lambda: run_reagents(args))
                    if os.path.exists(tsv_out):
                        ui.notification_show("Reagents recomputed.", type="message", duration=4)
                        shared_results_trigger.set(shared_results_trigger.get() + 1)
                    else:
                        ui.notification_show(
                            "Reagents recompute failed — run the Reagents task manually.",
                            type="warning", duration=8,
                        )
        except Exception as exc:
            ui.notification_show(f"Reagents check failed: {exc}",
                                 type="warning", duration=6)

        shared_values.set(new_path)

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
