"""progress_server.py — Shiny module server for the Progress tab.

Responsibilities:
- Display the current run JSON (or allow uploading one)
- Launch all analysis tasks concurrently via reactive.extended_task
- Show per-task status cards (gray → amber → green/red) that update while
  tasks run by polling the run-status JSON file every 2 seconds
- Persist task status + EBI jobIds to {working_dir}/{run_name}.status.json
  so a crashed session can be diagnosed and the CLI can resume
- Offer a download button to zip completed output files
"""

import asyncio
import io
import json
import os
import sys
import zipfile
from pathlib import Path

from shiny import reactive, render, ui, module

# ensure scripts/ is importable from the app context
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import run_status as _run_status
from task_runners import TASK_RUNNERS

# card background colors by status
_STATUS_COLORS = {
    "pending": "#cccccc",
    "running": "#fd7e14",
    "success": "#28a745",
    "failed":  "#dc3545",
}


@module.server
def progress_server(input, output, session, shared_json):

    # per-session state derived from the loaded JSON
    working_dir = reactive.Value("")
    run_name    = reactive.Value("")
    task_list   = reactive.Value([])  # [{id, analysis, output}, ...]

    # ── parse JSON whenever the shared path changes ────────────────────────────

    @reactive.effect
    @reactive.event(shared_json)
    def _parse_json():
        """Extract working_dir, run_name, and task list from the loaded JSON."""
        path = shared_json.get()
        if not path:
            return
        try:
            with open(path) as f:
                j = json.load(f)
        except Exception:
            return
        working_dir.set(j["global"].get("working_dir", ""))
        run_name.set(j["global"].get("run_name", ""))
        tasks = [
            {"id": k, "analysis": v["analysis"], "output": v["args"].get("output", "")}
            for k, v in j.items()
            if k not in ("scripts", "global")
        ]
        task_list.set(tasks)
        # initialise any not-yet-seen tasks as "pending" in the status file
        wd, rn = working_dir.get(), run_name.get()
        if wd and rn:
            for t in tasks:
                entry = _run_status.load_status(wd, rn).get(t["id"], {})
                if not entry:
                    _run_status.update_task(wd, rn, t["id"],
                                            status="pending", job_id="",
                                            message="", output=t["output"])

    # ── upload/display JSON ────────────────────────────────────────────────────

    @render.ui
    def show_or_upload_json():
        """Display the loaded JSON config or prompt to upload one."""
        elems = [ui.input_file("upload_json", "Upload / Reload JSON File", accept=[".json"])]
        path = shared_json.get()
        if path:
            try:
                with open(path) as f:
                    content = json.load(f)
                elems.append(ui.pre(json.dumps(content, indent=4)))
            except FileNotFoundError:
                elems.append(ui.div("Error: saved JSON file not found."))
        else:
            elems.append(ui.div("No JSON file loaded. Upload one or save an analysis in Setup."))
        return ui.div(*elems)

    @reactive.effect
    @reactive.event(input.upload_json)
    def _handle_upload():
        """Store an uploaded JSON path in shared state."""
        info = input.upload_json()
        if info:
            shared_json.set(info[0]["datapath"])

    # ── run-button state ───────────────────────────────────────────────────────

    @reactive.effect
    def _update_run_button():
        """Enable Run when JSON is loaded and no tasks are currently running."""
        has_json  = bool(shared_json.get())
        is_running = run_all_tasks.status() == "running"
        ui.update_action_button("run_analysis", disabled=not has_json or is_running)

    # ── concurrent task execution ──────────────────────────────────────────────

    @reactive.extended_task
    async def run_all_tasks(json_path, wd, rn):
        """Run all analysis tasks concurrently in a thread pool.

        Each task writes its status to the status file as it progresses
        so the polling UI can show live per-task updates.
        """
        with open(json_path) as f:
            j = json.load(f)
        tasks = {k: v for k, v in j.items() if k not in ("scripts", "global")}

        def _run_one(task_id, config):
            """Run a single task synchronously, updating status file on state changes."""
            analysis = config["analysis"]
            args     = config["args"]
            output   = args.get("output", "")
            runner   = TASK_RUNNERS.get(analysis)
            if runner is None:
                _run_status.update_task(wd, rn, task_id, status="failed",
                                        message=f"Unknown analysis type '{analysis}'")
                return

            # check if already complete (resume: skip finished tasks)
            entry = _run_status.load_status(wd, rn).get(task_id, {})
            if _run_status.is_complete(entry, output):
                return  # already done on a previous run

            _run_status.update_task(wd, rn, task_id, status="running", job_id="")

            # progress_cb captures EBI jobIds into the status file
            def _cb(job_id, status_str):
                _run_status.update_task(wd, rn, task_id, job_id=job_id)

            try:
                runner(args, progress_cb=_cb)
                _run_status.update_task(wd, rn, task_id, status="success", output=output)
            except Exception as exc:
                _run_status.update_task(wd, rn, task_id, status="failed",
                                        message=str(exc))

        loop = asyncio.get_running_loop()
        coros = [
            loop.run_in_executor(None, _run_one, tid, cfg)
            for tid, cfg in tasks.items()
        ]
        await asyncio.gather(*coros)
        return "done"

    @reactive.effect
    @reactive.event(input.run_analysis)
    def _launch():
        """Launch all tasks when the user clicks Run."""
        run_all_tasks.invoke(shared_json.get(), working_dir.get(), run_name.get())

    # ── live status polling ────────────────────────────────────────────────────

    @reactive.calc
    def _current_status():
        """Re-read the status file; auto-refreshes every 2 s while tasks are running."""
        if run_all_tasks.status() == "running":
            reactive.invalidate_later(2)
        wd, rn = working_dir.get(), run_name.get()
        if not wd or not rn:
            return {}
        return _run_status.load_status(wd, rn)

    # ── task cards ─────────────────────────────────────────────────────────────

    @render.ui
    def task_cards():
        """Render one status card per task; colors reflect current status."""
        tasks  = task_list.get()
        status = _current_status()

        if not tasks:
            return ui.div("No tasks loaded — save an analysis in Setup or upload a JSON.")

        cards = []
        for t in tasks:
            entry   = status.get(t["id"], {})
            st      = entry.get("status", "pending")
            job_id  = entry.get("job_id", "")
            message = entry.get("message", "")
            color   = _STATUS_COLORS.get(st, "#cccccc")

            body = [
                ui.tags.strong(t["id"]),
                ui.p(f"Analysis: {t['analysis']}", style="margin:0;"),
                ui.p(f"Status: {st}", style="margin:0;"),
            ]
            if job_id:
                body.append(ui.p(f"Job ID: {job_id}",
                                 style="margin:0; font-size:0.8em; font-family:monospace;"))
            if message:
                body.append(ui.p(f"Error: {message}",
                                 style="margin:0; font-size:0.85em; color:#dc3545;"))

            cards.append(
                ui.card(
                    ui.card_body(*body),
                    style=f"border-left: 6px solid {color}; margin-bottom: 8px;",
                )
            )
        return ui.div(*cards)

    # ── download completed results ─────────────────────────────────────────────

    @render.download(filename=lambda: f"{run_name.get() or 'results'}.zip")
    def download_results():
        """Zip all completed output files and stream to the browser."""
        tasks  = task_list.get()
        status = _current_status()

        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for t in tasks:
                entry = status.get(t["id"], {})
                if entry.get("status") == "success":
                    out = t.get("output", "")
                    if out and os.path.exists(out):
                        zf.write(out, os.path.basename(out))
        buf.seek(0)
        yield buf.read()
