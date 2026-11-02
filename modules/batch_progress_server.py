"""batch_progress_server.py — Progress tab of the batch app.

Runs scripts/batch_prepare.py in a background thread: database proteins are exported to run
folders (seconds for a whole batch), opted-in EBI proteins are analysed (minutes each), the
Results load cache is filled, and the deferred alignment PDFs are rendered last. The table
polls the in-memory batch while the thread runs. Results opens once the exports finish;
EBI proteins join the list as they complete.
"""

import sys
import threading
import time
from pathlib import Path

from shiny import module, reactive, render, ui

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import batch_prepare as BP  # noqa: E402
import batch_state as BS  # noqa: E402

_STATUS_CLASS = {"accepted": "text-success", "ready": "", "in_progress": "text-primary",
                 "failed": "text-danger", "missing": "text-muted", "running": "text-warning",
                 "pending": "text-muted"}


@module.server
def batch_progress_server(input, output, session, bstate):
    """bstate: the shared batch state (see batch_setup_server)."""
    job = {"thread": None, "phase": "", "t0": 0.0, "error": "", "exported": False,
           "opened": False, "render": None}

    def _run(path, data):
        """Background thread: materialize, preload, then render the alignment PDFs."""
        cfg = BS.app_config()
        try:
            job["phase"] = "exporting"
            BP.prepare(path, data, int(cfg.get("workers", 6)),
                       int(cfg.get("ebi_max_concurrent", 3)))
            job["phase"] = "loading into memory"
            BP.preload(data)
            job["exported"] = True
            # alignment images render in a separate process; Review never waits on them
            job["phase"] = "done (alignment images rendering in the background)"
            job["render"] = BP.render_in_subprocess(path)
            job["render"].wait()
            job["phase"] = "done"
        except Exception as e:
            job["error"] = "{}: {}".format(type(e).__name__, e)
            job["phase"] = "failed"
            job["exported"] = True

    @session.on_ended
    def _stop():
        """End the PDF-rendering process with the session (it is this app's own child)."""
        proc = job.get("render")
        if proc is not None and proc.poll() is None:
            proc.terminate()

    @reactive.effect
    @reactive.event(bstate["start"])
    def _start():
        """Start (or resume) the batch whenever Setup hands one over."""
        if bstate["data"] is None:
            return
        if job["thread"] is not None and job["thread"].is_alive():
            ui.notification_show("A batch is still preparing.", type="warning")
            return
        job.update(phase="starting", t0=time.time(), error="", exported=False, opened=False)
        job["thread"] = threading.Thread(target=_run, args=(bstate["path"], bstate["data"]),
                                         daemon=True)
        job["thread"].start()

    def _running():
        """True while the background thread is alive."""
        return job["thread"] is not None and job["thread"].is_alive()

    tick = reactive.Value(0)     # this tab's own refresh; bstate rev moves only on change
    _last_sig = [None]

    def _signature():
        """Statuses of every protein, to tell whether the shared state changed."""
        data = bstate["data"]
        return tuple(p.get("status") for p in data["proteins"]) if data else None

    @reactive.effect
    def _poll():
        """While running, refresh this tab every second; open Results when exports land."""
        bstate["start"].get()
        running = _running()
        if running:
            reactive.invalidate_later(1)
        with reactive.isolate():
            tick.set(tick.get() + 1)
            # other tabs (the Results protein list) re-render only when a status moved
            sig = _signature()
            if sig != _last_sig[0]:
                _last_sig[0] = sig
                bstate["rev"].set(bstate["rev"].get() + 1)
        if job["exported"] and not job["opened"]:
            job["opened"] = True
            ui.update_navset("main_tabs", selected="results", session=session.root_scope())

    @reactive.effect
    @reactive.event(input.review)
    def _review():
        """Go to Results."""
        ui.update_navset("main_tabs", selected="results", session=session.root_scope())

    @render.ui
    def summary():
        """Batch name, counts, phase and elapsed time."""
        bstate["rev"].get()
        tick.get()
        data = bstate["data"]
        if data is None:
            return ui.p("Set up or open a batch first.", class_="text-muted")
        c = BS.counts(data)
        parts = ["{} {}".format(v, k.replace("_", " ")) for k, v in c.items() if v]
        phase = job["phase"] or "idle"
        if job["t0"] and not phase.startswith(("done", "failed", "idle")):
            phase += " · {:.0f} s".format(time.time() - job["t0"])
        reviewable = c["ready"] + c["in_progress"] + c["accepted"]
        return ui.card(
            ui.h4(data["name"]),
            ui.p(" · ".join(parts)),
            ui.p("Phase: " + phase, class_="text-muted"),
            ui.p(job["error"], class_="text-danger") if job["error"] else ui.div(),
            ui.input_action_button("review", "Review sites →", class_="btn-success",
                                   disabled=reviewable == 0),
        )

    @render.ui
    def protein_table():
        """One row per protein: input, product, source, status, seconds, error."""
        bstate["rev"].get()
        tick.get()
        data = bstate["data"]
        if data is None:
            return ui.div()
        rows = []
        for p in data["proteins"]:
            st = p.get("status", "pending")
            rows.append(ui.tags.tr(
                ui.tags.td(p["input"]), ui.tags.td(p.get("product_id") or "—"),
                ui.tags.td(p.get("gene") or ""), ui.tags.td(p.get("source") or ""),
                ui.tags.td(st.replace("_", " "), class_=_STATUS_CLASS.get(st, "")),
                ui.tags.td(p.get("seconds", "")),
                ui.tags.td(p.get("error", ""), class_="small text-muted")))
        return ui.tags.table(
            ui.tags.thead(ui.tags.tr(*[ui.tags.th(h) for h in (
                "Input", "Product", "Gene", "Source", "Status", "s", "Note")])),
            ui.tags.tbody(*rows), class_="table table-sm")
