"""progress_server.py — Shiny module server for the Progress tab.

Thin reactive layer: loads the run JSON via shared_json, polls the status
file every 2 s while tasks are running, and renders a collapsible accordion
with one panel per job (status badge always visible; params hidden by default).
All data-shaping is delegated to progress_logic.py.
"""

import asyncio
import io
import json
import os
import sys
import zipfile
from pathlib import Path

from shiny import module, reactive, render, ui

# ensure scripts/ is importable from the app context
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import run_status as _run_status
from task_runners import TASK_RUNNERS, afdb_presearch
from progress import make_status_reporter, report as _progress_report

from modules.progress_logic import display_params, parse_run, status_label
from modules.setup_ui import label_with_tip
from modules.json_card import json_upload_card as _json_upload_card


# ── helpers ───────────────────────────────────────────────────────────────────

def _status_badge(state):
    """Return a colored status badge element for an accordion title."""
    return ui.tags.span(state, class_=f"status-badge status-{state}")


def _stage_chip(stage):
    """Return a small muted stage chip, or empty span when no stage is set."""
    if not stage or stage in ("done", "failed"):
        return ui.span()
    return ui.tags.span(stage, class_="stage-chip")


def _param_grid(args, global_block, analysis=None):
    """Build a read-only 2-column grid of per-job parameters."""
    params = display_params(args, global_block, analysis=analysis)
    if not params:
        return ui.p("No configurable parameters.", class_="section-hint")
    cells = []
    for k, v in params.items():
        cells.append(ui.div(k.replace("_", " "), class_="param-label"))
        cells.append(ui.div(str(v) if v != "" else "—", class_="param-value"))
    return ui.div(*cells, class_="param-grid")


def _build_body(task, global_block, entry):
    """Build the full accordion panel body: param grid + log textarea + re-run button."""
    children = [_param_grid(task["args"], global_block, analysis=task.get("analysis"))]
    log = entry.get("log", "")
    children.append(
        ui.tags.textarea(
            log,
            readonly=True,
            class_="task-log",
            rows=6,
            id=f"task-log-{task['id']}",
        )
    )
    children.append(
        ui.input_action_button(
            f"rerun_{task['id']}", "↻ Re-run",
            class_="btn-outline-warning btn-sm mt-2",
        )
    )
    return ui.div(*children)


# ── module server ─────────────────────────────────────────────────────────────

@module.server
def progress_server(input, output, session, shared_json, shared_results_trigger):

    # per-session state derived from the loaded JSON
    global_block = reactive.Value({})
    working_dir  = reactive.Value("")
    run_name     = reactive.Value(None)
    task_list    = reactive.Value([])   # [{"id","analysis","args","output"}, ...]

    # re-run queue: tracks which tasks the user has asked to re-run
    rerun_queue    = reactive.Value(frozenset())  # reactive mirror of _click_baseline
    _click_baseline = {}   # {tid: click_count captured at last Run} — plain dict, not reactive
    _baseline_tick  = reactive.Value(0)           # bumped on each Run to trigger re-sync

    # ── parse JSON whenever the shared path changes ────────────────────────────

    @reactive.effect
    @reactive.event(shared_json)
    def _parse_json():
        """Extract global block and task list from the loaded run JSON."""
        path = shared_json.get()
        if not path:
            return
        try:
            with open(path) as f:
                j = json.load(f)
        except Exception:
            return

        gb, tasks = parse_run(j)
        global_block.set(gb)
        working_dir.set(gb.get("working_dir", ""))
        run_name.set(gb.get("run_name", ""))
        task_list.set(tasks)

        # seed any not-yet-seen tasks as "pending" in the status file
        wd, rn = gb.get("working_dir", ""), gb.get("run_name", "")
        if wd and rn:
            for t in tasks:
                entry = _run_status.load_status(wd, rn).get(t["id"], {})
                if not entry:
                    _run_status.update_task(wd, rn, t["id"],
                                            status="pending", job_id="",
                                            message="", output=t["output"])

    # ── upload handler ─────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.upload_json)
    def _handle_upload():
        """Store an uploaded JSON path in shared state."""
        info = input.upload_json()
        if info:
            shared_json.set(info[0]["datapath"])

    # ── re-run queue: track per-task re-run button clicks ────────────────────

    @reactive.effect
    def _sync_rerun_queue():
        """Auto-launch any task whose re-run button was just clicked.

        If a run is already in progress the clicked tasks are queued instead
        (same as before) so they can be included in the next manual Run.
        """
        _baseline_tick.get()      # re-evaluate after each Run resets the baseline
        tasks = task_list.get()   # re-evaluate when the task list changes
        newly_clicked = set()
        for t in tasks:
            tid = t["id"]
            try:
                clicks = input[f"rerun_{tid}"]() or 0
            except Exception:
                clicks = 0
            if clicks > _click_baseline.get(tid, 0):
                newly_clicked.add(tid)
        rerun_queue.set(frozenset(newly_clicked))
        # auto-launch immediately; fall back to queue when a run is in progress
        if newly_clicked and run_all_tasks.status() != "running":
            _do_launch(frozenset(newly_clicked))

    # ── run-button state ───────────────────────────────────────────────────────

    @reactive.effect
    def _toggle_run():
        """Enable Run/Re-run All when JSON is loaded and no tasks are currently running."""
        has_json   = bool(shared_json.get())
        is_running = run_all_tasks.status() == "running"
        disabled   = not has_json or is_running
        ui.update_action_button("run_analysis", disabled=disabled)
        ui.update_action_button("rerun_all",    disabled=disabled)

    # ── run header ─────────────────────────────────────────────────────────────

    @render.ui
    def json_upload_card():
        """Render the JSON upload card; collapsed and warning-free when a run is loaded."""
        return _json_upload_card(
            "upload_json", "ts-prog-json-body", "Upload run JSON",
            run_name.get() is not None,
        )

    @render.ui
    def run_header():
        """Show a one-line run summary when a JSON is loaded; nothing otherwise."""
        if run_name.get() is None:
            return ui.span()
        rn = run_name.get()
        wd = working_dir.get()
        n  = len(task_list.get())
        return ui.p(
            ui.tags.strong(rn or "(unnamed)"),
            f" · {wd or '—'} · {n} task{'s' if n != 1 else ''}",
            class_="run-header",
        )

    # ── concurrent task execution ──────────────────────────────────────────────

    @reactive.extended_task
    async def run_all_tasks(json_path, wd, rn, rerun_set):
        """Run all analysis tasks concurrently in a thread pool.

        Each task writes its status to the status file as it progresses
        so the polling UI can show live per-task updates.
        rerun_set is a frozenset of task IDs to re-run even if already complete.
        """
        with open(json_path) as f:
            j = json.load(f)
        _, tasks = parse_run(j)

        def _run_one(task):
            """Run a single task synchronously, updating the status file on each state change."""
            tid      = task["id"]
            analysis = task["analysis"]
            args     = task["args"]
            output   = task["output"]
            runner   = TASK_RUNNERS.get(analysis)
            if runner is None:
                _run_status.update_task(wd, rn, tid, status="failed",
                                        message=f"Unknown analysis type '{analysis}'",
                                        log=f"ERROR: Unknown analysis type '{analysis}'")
                return

            # skip if already complete and not explicitly queued for re-run
            entry = _run_status.load_status(wd, rn).get(tid, {})
            if _run_status.is_complete(entry, output) and tid not in rerun_set:
                return

            log = []
            _run_status.update_task(wd, rn, tid, status="running", job_id="", log="", stage="")
            reporter = make_status_reporter(wd, rn, tid, log)

            try:
                runner(args, report=reporter)
                _progress_report(reporter, "Done.", stage="done")
                _run_status.update_task(wd, rn, tid, status="success", output=output)
            except Exception as exc:
                import traceback
                log.append(traceback.format_exc())
                _run_status.update_task(wd, rn, tid, status="failed", message=str(exc),
                                        log="\n".join(log), stage="failed")

        loop = asyncio.get_running_loop()

        # AFDB pre-step: must finish before any task starts so all analyses share
        # the same reference sequence (AF sequence replaces user .fa if a model is found)
        plddt_task = next(
            (t for t in tasks if t["analysis"] == "plddt" and not t["args"].get("pdb")),
            None,
        )
        if plddt_task:
            afdb_log = []
            _run_status.update_task(wd, rn, plddt_task["id"],
                                    status="running", stage="afdb_lookup")
            afdb_reporter = make_status_reporter(wd, rn, plddt_task["id"], afdb_log)
            await loop.run_in_executor(None, lambda: afdb_presearch(tasks, report=afdb_reporter))

        await asyncio.gather(*[
            loop.run_in_executor(None, _run_one, t) for t in tasks
        ])
        return "done"

    def _do_launch(rerun_set):
        """Shared launch logic: stamp queued tasks, advance baseline, invoke runner."""
        wd, rn = working_dir.get(), run_name.get()
        if wd and rn:
            for t in task_list.get():
                if t["id"] in rerun_set:
                    _run_status.update_task(wd, rn, t["id"], status="queued", log="", stage="")
        for t in task_list.get():
            try:
                _click_baseline[t["id"]] = input[f"rerun_{t['id']}"]() or 0
            except Exception:
                pass
        _baseline_tick.set(_baseline_tick.get() + 1)
        run_all_tasks.invoke(shared_json.get(), wd, rn, rerun_set)

    # plain mutable guard — avoids re-signaling on the same completion
    _completion_count = [0]
    _completion_notified = [False]

    @reactive.effect
    def _on_tasks_complete():
        """Increment shared_results_trigger once when tasks finish so results reload."""
        status = run_all_tasks.status()
        if status == "running":
            _completion_notified[0] = False
        elif status == "success" and not _completion_notified[0]:
            _completion_notified[0] = True
            _completion_count[0] += 1
            shared_results_trigger.set(_completion_count[0])

    @reactive.effect
    @reactive.event(input.run_analysis)
    def _launch():
        """Launch all tasks; any task queued for re-run bypasses the is_complete guard."""
        _do_launch(frozenset(rerun_queue.get()))

    @reactive.effect
    @reactive.event(input.rerun_all)
    def _launch_all():
        """Re-run every task unconditionally."""
        _do_launch(frozenset(t["id"] for t in task_list.get()))

    # ── live status polling ────────────────────────────────────────────────────

    @reactive.calc
    def _current_status():
        """Re-read the status file; auto-refreshes every 2 s while tasks are running."""
        if run_all_tasks.status() == "running":
            reactive.invalidate_later(2)
        wd, rn = working_dir.get(), run_name.get()
        if not wd or not rn:
            return {}
        status = _run_status.load_status(wd, rn)
        # Also poll while tasks are "queued": covers the race between _do_launch
        # writing "queued" and the extended_task transitioning to "running".
        if any(e.get("status") == "queued" for e in status.values()):
            reactive.invalidate_later(2)
        return status

    # ── task accordion ─────────────────────────────────────────────────────────
    #
    # Two-step render to prevent the 2-second status poll from collapsing open
    # panels:
    #   1. task_cards — builds the accordion DOM; fires only when task_list
    #      changes (new JSON loaded), so open/closed state is never reset.
    #   2. _update_status_badges — fires every poll cycle and surgically updates
    #      each panel's title and body via update_accordion_panel, leaving the
    #      open/closed state untouched.

    _ACCORDION_ID = "task_accordion"

    @render.ui
    @reactive.event(task_list, global_block)
    def task_cards():
        """Render accordion skeleton; re-renders only when the task list changes."""
        tasks = task_list.get()
        gb    = global_block.get()

        if not tasks:
            return ui.span()

        # read status directly (not via reactive calc) to avoid a polling dependency
        wd, rn = working_dir.get(), run_name.get()
        status = _run_status.load_status(wd, rn) if (wd and rn) else {}

        panels = []
        for t in tasks:
            tid   = t["id"]
            label = tid.rsplit("_", 1)[0]
            entry = status.get(tid, {})
            state = status_label(entry)
            stage = entry.get("stage", "")
            panels.append(
                ui.accordion_panel(
                    ui.span(label,
                            ui.tags.span(t["analysis"], class_="task-type-badge"),
                            _status_badge(state),
                            _stage_chip(stage)),
                    _build_body(t, gb, entry),
                    value=tid,
                )
            )

        return ui.div(
            ui.p("Each panel shows live job status; expand to see parameters.",
                 class_="section-hint"),
            ui.accordion(*panels, open=False, multiple=True,
                         class_="task-accordion", id=_ACCORDION_ID),
        )

    @reactive.effect
    async def _update_status_badges():
        """Update accordion titles + push log text in-place (no body re-render)."""
        tasks  = task_list.get()
        status = _current_status()   # provides the 2 s invalidation while running
        queued = rerun_queue.get()   # tasks the user has queued for re-run
        log_updates = []
        for t in tasks:
            tid   = t["id"]
            label = tid.rsplit("_", 1)[0]
            entry = status.get(tid, {})
            # show "queued" badge when task is pending re-run; otherwise real status
            state = "queued" if tid in queued else status_label(entry)
            stage = entry.get("stage", "")
            # update title only — omitting body preserves textarea scroll/size
            ui.update_accordion_panel(
                _ACCORDION_ID, tid,
                title=ui.span(label,
                              ui.tags.span(t["analysis"], class_="task-type-badge"),
                              _status_badge(state),
                              _stage_chip(stage)),
            )
            log_updates.append({"id": f"task-log-{tid}", "log": entry.get("log", "")})
        # push all log updates to JS in one message
        await session.send_custom_message("tagsites_update_logs", {"updates": log_updates})

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
