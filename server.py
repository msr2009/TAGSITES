import asyncio
import io
import json
import os
import shutil
import sys
import tempfile
import time
import zipfile
from pathlib import Path

from shiny import reactive, render, ui

from modules.setup_server import setup_server
from modules.progress_server import progress_server
from modules.results_server import results_server
from modules.reagents_server import reagents_server
from utils.bundle import is_bundle_zip, extract_and_rebase

_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT / "scripts"))
import ebi_health
from config import EBI_HEALTH_BANNER


def app_server(input, output, session):
    # per-session shared state: the path to the active run JSON
    # created inside app_server so each browser session gets its own value
    shared_values = reactive.Value("")

    # chosen tag sites — list of residue positions (ints) from the Results page
    shared_sites = reactive.Value([])

    # incremented by progress_server each time a run completes so results_server
    # reloads even when the JSON path hasn't changed
    shared_results_trigger = reactive.Value(0)

    # bumped by setup_server after Save Analysis is pressed, so progress_server
    # auto-starts all tasks once the new run JSON has been parsed
    shared_autostart = reactive.Value(0)

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
            await session.send_custom_message("tagsites_bundle_done", {})
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

        await session.send_custom_message("tagsites_bundle_done", {})
        shared_values.set(new_path)

    setup_server("setup", shared_json=shared_values, shared_autostart=shared_autostart)
    progress_server("progress", shared_json=shared_values,
                    shared_results_trigger=shared_results_trigger,
                    shared_autostart=shared_autostart)
    results_server("results", shared_json=shared_values, shared_sites=shared_sites,
                   shared_results_trigger=shared_results_trigger)
    reagents_server("reagents", shared_json=shared_values, shared_sites=shared_sites)

    # ── EBI health banner ──────────────────────────────────────────────────────
    # gated by config.EBI_HEALTH_BANNER; when off, nothing here runs or probes EBI
    if EBI_HEALTH_BANNER:
        # latest ebi_health verdict; set only when level or problems change, so the
        # banner does not re-render (and un-dismiss) on every poll
        ebi_state = reactive.Value(None)

        # the Setup module's inputs ("setup-email" is not a valid id at the root scope)
        _setup_input = session.make_scope("setup").input

        def _setup_email():
            """The Setup tab's Email field (also filled from a loaded run JSON), or ''."""
            with reactive.isolate():
                return _setup_input["email"]() if "email" in _setup_input else ""

        async def _refresh_ebi(force=False):
            """Fetch the (process-wide cached) EBI verdict off the event loop."""
            try:
                res = await asyncio.to_thread(ebi_health.check, force, _setup_email())
            except Exception:
                return
            old = ebi_state.get()
            # skip no-op updates: same level, problems and submission result
            key = ("level", "problems", "accepting", "email")
            if old is None or [old.get(k) for k in key] != [res.get(k) for k in key]:
                ebi_state.set(res)

        @reactive.effect
        async def _poll_ebi():
            """Re-read the EBI verdict every minute (network only when the cache expires)."""
            reactive.invalidate_later(60)
            await _refresh_ebi()

        @reactive.effect
        async def _ebi_email_changed():
            """Test submissions once a real email is entered (debounced against typing)."""
            email = _setup_input["email"]() if "email" in _setup_input else ""
            # placeholder or half-typed: nothing to test yet
            if not ebi_health._usable_email(email) or "." not in email.split("@")[-1]:
                return
            await asyncio.sleep(1.5)
            # still typing: a later run of this effect handles the final address
            if _setup_email() != email:
                return
            await _refresh_ebi()

        @reactive.effect
        @reactive.event(input.ebi_recheck)
        async def _recheck_ebi():
            """Force a fresh EBI probe from the banner's Check-again button."""
            await _refresh_ebi(force=True)
            res = ebi_state.get()
            if res and res["level"] == "ok":
                ui.notification_show("EBI services are working normally.", duration=4)

        @render.ui
        def ebi_banner():
            """Warning strip above the tabs while EBI is slow, impaired, blocking or down."""
            res = ebi_state.get()
            if not res:
                return ui.span()
            # green strip for a healthy EBI (the banner is opt-in, so showing OK is wanted)
            if res["level"] == "ok":
                lat = res.get("latency_s")
                checked = time.strftime("%H:%M:%S", time.localtime(res["checked_at"]))
                return ui.div(
                    ui.tags.strong("✓ " + ebi_health.summary(res)),
                    (f" Accepting jobs for {res['email']}: "
                     f"{({True: 'yes', False: 'no'}).get(res.get('accepting'), 'unknown')}"
                     if res.get("email") else
                     " Job submission not tested yet (enter an email in Setup).")
                    + (f" · slowest service {lat:.2f} s" if lat is not None else "")
                    + f" · checked {checked} ",
                    ui.input_action_button("ebi_recheck", "Check again",
                                           class_="btn-sm btn-outline-secondary"),
                    class_="alert alert-success",
                    style="margin:6px 0;padding:8px 12px;font-size:0.9em;",
                )
            cls = "alert-info" if res["level"] == "slow" else (
                "alert-warning" if res["level"] == "degraded" else "alert-danger")
            return ui.div(
                ui.tags.strong("⚠ " + ebi_health.summary(res)),
                ui.tags.ul(*[ui.tags.li(p) for p in res["problems"]],
                           style="margin:4px 0 4px 0;"),
                ui.input_action_button("ebi_recheck", "Check again",
                                       class_="btn-sm btn-outline-secondary"),
                ui.tags.button(type="button", class_="btn-close",
                               **{"data-bs-dismiss": "alert", "aria-label": "Close"}),
                class_=f"alert {cls} alert-dismissible",
                style="margin:6px 0;padding:8px 40px 8px 12px;font-size:0.9em;",
            )

    @render.download_button(filename="tagsites-app.zip", media_type="application/zip")
    def download_app():
        """Bundle app source files into a zip on disk and return its path (sends Content-Length)."""
        # strict_timestamps=False: deployed files have pre-1980 mtimes, which zipfile rejects
        zip_path = os.path.join(_get_session_dir(), "tagsites-app.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, strict_timestamps=False) as zf:
            for arcname, filepath in _collect_app_files(_ROOT):
                zf.write(filepath, arcname)
        return zip_path


def _collect_app_files(root):
    """Yield (arcname, filepath) pairs for app source files to include in the download zip."""
    for pattern in ("*.py", "*.json", "*.yml", "*.yaml", "*.txt"):
        for f in sorted(root.glob(pattern)):
            # *.local.json holds secrets (UCSC BLAT API key) and machine-specific overrides
            if f.name.endswith(".local.json"):
                continue
            yield f.name, f
    readme = root / "README.md"
    if readme.exists():
        yield "README.md", readme
    for subdir in ("modules", "scripts", "utils", "tables", "taxid_lists", "docs", "www", "params"):
        d = root / subdir
        if not d.is_dir():
            continue
        for f in sorted(d.rglob("*")):
            if f.is_file() and "__pycache__" not in str(f) and not f.name.startswith("."):
                yield str(f.relative_to(root)), f
