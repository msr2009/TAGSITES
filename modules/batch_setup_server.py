"""batch_setup_server.py — Setup tab of the batch app.

Resolves a pasted or uploaded protein list against the proteome database
(scripts/batch_inputs.py), lets the user pick among a gene's products, and writes the batch
file (scripts/batch_state.py). Opening a saved batch skips all of that. Either way the batch
lands in the shared state and the Progress tab takes over.
"""

import os
import sys
import time
from pathlib import Path

from shiny import module, reactive, render, ui

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import batch_bundle as BB  # noqa: E402
import batch_inputs as BI  # noqa: E402
import batch_state as BS  # noqa: E402
import reagent_bundle as RB  # noqa: E402

from modules.setup_server import _off_loop  # noqa: E402

_TAGS = RB.load_tag_table()


def _entry_note(e):
    """Short status text for one resolved entry."""
    if e.get("duplicate_of"):
        return "duplicate of {} (dropped)".format(e["duplicate_of"])
    if e["source"] == "missing":
        return "not in database"
    if e.get("reagent_status") != "success":
        return "in database · reagent design failed"
    return "in database"


@module.server
def batch_setup_server(input, output, session, bstate):
    """bstate: shared {"path", "data"} dict plus reactive "rev" (refresh) and "start" (run)."""
    entries = reactive.Value([])        # resolved entries, before the batch is saved
    listing = reactive.Value(0)         # bumped to refresh the saved-batch list

    # ── defaults ─────────────────────────────────────────────────────────────
    @reactive.effect
    def _default_name():
        """Fill the batch name with batch_YYYYMMDD once."""
        with reactive.isolate():
            if not input.batch_name():
                ui.update_text("batch_name", value=time.strftime("batch_%Y%m%d"))

    @reactive.effect
    @reactive.event(input.repair_type)
    def _arm_default():
        """Arm length follows the repair type's default, as in the Reagents tab."""
        ui.update_numeric("arm_length", value=RB.ARM_DEFAULTS.get(input.repair_type(), 500))

    @render.ui
    def insert_select():
        """Insert library multi-select (tables/tags.tsv plus Custom)."""
        choices = list(_TAGS) + ["Custom"]
        return ui.input_selectize("inserts", "Insert(s)", choices=choices,
                                  selected=[c for c in choices if c == "FLAG"][:1],
                                  multiple=True, width="100%")

    # ── open a saved batch ───────────────────────────────────────────────────
    @render.ui
    def open_batch_card():
        """Select of saved batches under the batch folder, newest first."""
        listing.get()
        rows = BS.list_batches()
        upload = ui.input_file("batch_zip", "…or load a batch ZIP (⬇ Batch in Review)",
                               accept=[".zip"], width="100%")
        if not rows:
            return ui.div(ui.p("No saved batches in {}".format(BS.batch_root()),
                               class_="text-muted small"), upload)
        choices = {r["path"]: "{} · {}/{} accepted · {}".format(
            r["name"], r["accepted"], r["n"], r["created"]) for r in rows}
        return ui.div(
            ui.input_select("open_path", None, choices=choices, width="100%"),
            ui.input_action_button("open_batch", "Open batch", class_="btn-primary"),
            ui.p("Settings and product choices come from the saved file; finished "
                 "proteins are not exported again.", class_="text-muted small mt-2"),
            upload,
        )

    def _hand_over(path):
        """Put a batch file in the shared state and resume it in Progress."""
        bstate["path"], bstate["data"] = path, BS.load_batch(path)
        bstate["rev"].set(bstate["rev"].get() + 1)
        bstate["start"].set(bstate["start"].get() + 1)
        ui.update_navset("main_tabs", selected="progress", session=session.root_scope())

    @reactive.effect
    @reactive.event(input.batch_zip)
    async def _load_zip():
        """Restore an uploaded batch ZIP under the batch folder, then open it."""
        info = input.batch_zip()
        if not info:
            return
        zpath = info[0]["datapath"]
        if not BB.is_batch_zip(zpath):
            ui.notification_show("Not a batch ZIP (no BATCH_MANIFEST.json). A single-protein "
                                 "bundle opens in the single-protein app.", type="error",
                                 duration=8)
            return
        ui.notification_show("Restoring batch…", id="batch_unpack", duration=None)
        try:
            path = await _off_loop(BB.unpack, zpath)
        except Exception as e:
            ui.notification_show("Could not restore the batch: {}".format(e), type="error",
                                 duration=10)
            return
        finally:
            ui.notification_remove("batch_unpack")
        listing.set(listing.get() + 1)
        _hand_over(path)

    @reactive.effect
    @reactive.event(input.open_batch)
    def _open():
        """Load a saved batch into the shared state and resume it in Progress."""
        path = input.open_path()
        if not path or not os.path.exists(path):
            return
        _hand_over(path)

    # ── resolve a new list ───────────────────────────────────────────────────
    @reactive.effect
    @reactive.event(input.list_file)
    def _read_upload():
        """Copy an uploaded list into the text box, where it can still be edited."""
        info = input.list_file()
        if info:
            ui.update_text_area("protein_list", value=Path(info[0]["datapath"]).read_text())

    @reactive.effect
    @reactive.event(input.resolve)
    async def _resolve():
        """Look every term up in the database off the event loop."""
        parsed = BI.parse_list(input.protein_list(), _TAGS)
        if not parsed:
            ui.notification_show("Paste or upload at least one protein.", type="warning")
            return
        try:
            res = await _off_loop(BI.resolve, parsed, budget=120)
        except Exception as e:
            ui.notification_show("Lookup failed: {}".format(e), type="error", duration=10)
            return
        entries.set(res)

    @render.ui
    def resolved_card():
        """Resolution table: input, chosen product (dropdown when several), status, insert."""
        ents = entries.get()
        if not ents:
            return ui.div()
        rows = []
        for i, e in enumerate(ents):
            # a gene with several products gets a dropdown, preselected per batch_inputs
            if len(e["choices"]) > 1 and not e.get("duplicate_of"):
                pick = ui.input_select("pick_{}".format(i), None, width="100%",
                                       choices={c["id"]: BI.choice_label(c)
                                                for c in e["choices"]},
                                       selected=e["product_id"])
            elif e["choices"]:
                pick = BI.choice_label(e["choices"][0])
            else:
                pick = "—"
            note = _entry_note(e)
            cls = ("text-danger" if e["source"] == "missing" or "failed" in note
                   else "text-muted" if e.get("duplicate_of") else "")
            rows.append(ui.tags.tr(
                ui.tags.td(e["input"]), ui.tags.td(e.get("gene") or ""),
                ui.tags.td(pick, style="min-width:360px"),
                ui.tags.td(note, class_=cls), ui.tags.td(e.get("insert") or "batch default")))
        n_db = sum(e["source"] == "db" and not e.get("duplicate_of") for e in ents)
        n_miss = sum(e["source"] == "missing" for e in ents)
        return ui.card(
            ui.card_header("{} in database · {} not found · {} duplicate".format(
                n_db, n_miss, sum(bool(e.get("duplicate_of")) for e in ents))),
            ui.tags.table(
                ui.tags.thead(ui.tags.tr(*[ui.tags.th(h) for h in
                                           ("Input", "Gene", "Product", "Status", "Insert")])),
                ui.tags.tbody(*rows), class_="table table-sm"),
        )

    def _chosen_entries():
        """Resolved entries with each dropdown's product applied, duplicates re-marked."""
        out = []
        for i, e in enumerate(entries.get()):
            if len(e["choices"]) > 1 and not e.get("duplicate_of"):
                try:
                    pid = input["pick_{}".format(i)]()
                except Exception:
                    pid = e["product_id"]
                if pid and pid != e["product_id"]:
                    e = BI.choose(e, pid)
            out.append(e)
        return BI.mark_duplicates(out)

    # ── save & continue ──────────────────────────────────────────────────────
    @reactive.effect
    @reactive.event(input.save_continue)
    def _save():
        """Write the batch file, then hand over to Progress."""
        ents = _chosen_entries()
        if not ents:
            ui.notification_show("Look up a protein list first.", type="warning")
            return
        name = BS.safe_name(input.batch_name() or time.strftime("batch_%Y%m%d"))
        if BS.batch_file(name).exists():
            ui.notification_show("A batch named '{}' exists; open it or choose another "
                                 "name.".format(name), type="error", duration=8)
            return
        inserts = list(input.inserts() or [])
        if "Custom" in inserts and not input.custom_insert().strip():
            ui.notification_show("Enter the custom insert DNA.", type="warning")
            return
        if input.ebi_for_missing() and not input.email().strip():
            ui.notification_show("EBI analyses need an email address.", type="warning")
            return
        settings = {"repair_type": input.repair_type(), "arm_length": int(input.arm_length()),
                    "inserts": inserts, "custom_insert": input.custom_insert().strip(),
                    "ebi_for_missing": bool(input.ebi_for_missing()),
                    "email": input.email().strip(),
                    "site_filters": {
                        "n_suggested": max(1, int(input.n_suggested() or 2)),
                        "no_nterm": bool(input.filt_no_nterm()),
                        "no_cterm": bool(input.filt_no_cterm()),
                        "guide_within": int(input.filt_guide_bp() or 0) if input.filt_guide()
                        else 0,
                        "constitutive_only": bool(input.filt_constitutive())}}
        path = BS.new_batch(name, settings, ents)
        bstate["path"], bstate["data"] = path, BS.load_batch(path)
        bstate["rev"].set(bstate["rev"].get() + 1)
        bstate["start"].set(bstate["start"].get() + 1)
        listing.set(listing.get() + 1)
        ui.update_navset("main_tabs", selected="progress", session=session.root_scope())
