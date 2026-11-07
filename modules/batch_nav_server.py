"""batch_nav_server.py — moves the Results tab through a batch and records the review.

Switching protein points shared_json at that protein's run folder, so the unchanged Results
module loads it (from the in-memory cache the Progress tab filled). Selected sites mirror
into the batch file; each site's guide defaults to the nearest one and can be changed here.
Accept stores the protein's sites and guides and moves to the next protein still to review.
Output Reagents builds the ZIP of every accepted protein (scripts/batch_output.py).

Statuses: opening a ready protein marks it in_progress; editing an accepted protein's sites
returns it to in_progress until it is accepted again. Skip marks a protein skipped (red, italic
in the list; written to {batch}.skipped.txt) and moves on; Unskip restores its earlier status.
Editing a skipped protein's sites leaves it skipped; accepting it un-skips it.
Opening a protein that is not accepted always selects freshly computed suggestions; only
accepted proteins keep their stored sites (server_batch._kept_sites).
"""

import asyncio
import os
import sys
import time
from pathlib import Path

from shiny import module, reactive, render, ui

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import batch_bundle as BB  # noqa: E402
import batch_output as BO  # noqa: E402
import batch_state as BS  # noqa: E402
import reagent_bundle as RB  # noqa: E402

_NAVIGABLE = ("ready", "in_progress", "accepted", "skipped", "failed")
_MARK = {"accepted": "✓", "in_progress": "●", "ready": "○", "failed": "✕", "skipped": "⤼"}


def _same(a, b):
    """True when two paths name the same file."""
    return bool(a and b) and os.path.normpath(a) == os.path.normpath(b)


def _guide_text(row):
    """Guide choice label: distance, strand, spacer; RS3 / off-target as badges only."""
    bits = ["{} bp".format(int(row["distance"])), "strand {}".format(row["guide_strand"]),
            str(row["spacer"])]
    rs3 = row.get("rs3_score", "")
    if rs3 == rs3 and rs3 not in ("", None):
        bits.append("RS3 {:.2f}".format(float(rs3)))
    ot = row.get("offtarget_count", "")
    if ot == ot and ot not in ("", None):
        bits.append("off-target {}".format(int(float(ot))))
    return " · ".join(bits)


@module.server
def batch_nav_server(input, output, session, bstate, shared_json, shared_sites, loaded_path,
                     shared_regions=None):
    """bstate: shared batch state; loaded_path: set by results_server when a run is loaded."""
    current = reactive.Value(None)       # key of the protein on screen
    _reagents = {}                       # key -> arm-less reagents frame (guide choices)

    def _data():
        """The batch dict, or None."""
        return bstate["data"]

    def _save():
        """Persist the batch file."""
        if bstate["path"] and bstate["data"] is not None:
            BS.save_batch(bstate["path"], bstate["data"])

    def _navigable():
        """Proteins that have a run folder to show, in batch order."""
        data = _data()
        if data is None:
            return []
        return [p for p in data["proteins"]
                if p.get("status") in _NAVIGABLE and p.get("run_json")]

    def _protein(key):
        """The entry for key, or None."""
        data = _data()
        return BS.protein(data, key) if data is not None and key else None

    def _reagents_df(p):
        """The protein's reagents frame (arm-less is enough for guide choices), cached."""
        if p["key"] not in _reagents:
            try:
                _, tsv, *_ = RB.run_paths(p["run_json"])
                _reagents[p["key"]] = RB.load_reagents_tsv(tsv) if tsv and os.path.exists(
                    tsv) else None
            except Exception:
                _reagents[p["key"]] = None
        return _reagents[p["key"]]

    # ── protein list ─────────────────────────────────────────────────────────
    @reactive.effect
    @reactive.event(bstate["start"])
    def _batch_changed():
        """Reset the nav when another batch is opened; nothing survives from the last one.

        Without this the key of a protein from the previous batch stayed current. It is not
        in the new batch, so the counter read '0 / 3' and Results went on showing that
        protein until the first Next.
        """
        _reagents.clear()
        with reactive.isolate():
            # the sites are deliberately left alone: results_server writes them back to
            # whatever run it last loaded, so clearing them here would empty the previous
            # batch's selection on disk. Loading this batch's protein replaces them
            shared_json.set("")
            # back to None so _fill_list picks this batch's own starting protein, and so the
            # move onto it is a change even when both batches hold the same key
            current.set(None)
            bstate["rev"].set(bstate["rev"].get() + 1)

    @reactive.effect
    def _fill_list():
        """Refresh the protein dropdown (labels carry the status for colouring)."""
        bstate["rev"].get()
        items = _navigable()
        with reactive.isolate():
            cur = current.get()
        choices = {p["key"]: "{}|{} {}  {}".format(
            p["status"], _MARK.get(p["status"], ""), BS.name_prefix(p),
            "({} sites)".format(len(p.get("sites") or [])) if p.get("sites") else "")
            for p in items}
        if cur is None and items:
            data = _data()
            cur = data.get("current_key") if data.get("current_key") in choices else None
            cur = cur or BS.next_unaccepted(data) or items[0]["key"]
        ui.update_selectize("protein", choices=choices, selected=cur)
        if cur != current.get():
            current.set(cur)

    @reactive.effect
    @reactive.event(input.protein)
    def _picked():
        """Dropdown choice -> current protein."""
        key = input.protein()
        if key and key != current.get():
            current.set(key)

    @reactive.effect
    def _show():
        """Load the current protein into Results and mark it in progress."""
        key = current.get()
        p = _protein(key)
        if p is None:
            return
        with reactive.isolate():
            data = _data()
            data["current_key"] = key
            if p.get("status") == "ready":
                p["status"] = "in_progress"
                bstate["rev"].set(bstate["rev"].get() + 1)
            _save()
            if input.protein() != key:
                ui.update_selectize("protein", selected=key)
            if not _same(shared_json.get(), p["run_json"]):
                shared_json.set(p["run_json"])

    def _step(delta):
        """Move delta proteins through the navigable list (wrapping)."""
        items = [p["key"] for p in _navigable()]
        if not items:
            return
        cur = current.get()
        i = items.index(cur) if cur in items else -1
        current.set(items[(i + delta) % len(items)])

    @reactive.effect
    @reactive.event(input.prev)
    def _prev():
        """Previous protein."""
        _step(-1)

    @reactive.effect
    @reactive.event(input.next)
    def _next():
        """Next protein."""
        _step(1)

    # ── sites and guides mirror into the batch ───────────────────────────────
    @reactive.effect
    def _mirror_sites():
        """Copy the selected sites of the loaded protein into its batch entry."""
        sites = sorted(int(s) for s in shared_sites.get())
        path = loaded_path.get()
        # a dependency, not an isolated read: the regions are published when the plot is
        # built, which can be after the sites arrive, and the first protein of a batch is
        # where that order shows up most often
        published = (shared_regions.get() or {}) if shared_regions is not None else {}
        with reactive.isolate():
            p = _protein(current.get())
            # only once Results has loaded this protein (sites still show the last one before)
            if p is None or not _same(path, p["run_json"]):
                return
            guides = {str(k): v for k, v in (p.get("guides") or {}).items()
                      if int(k) in sites}
            df = _reagents_df(p)
            for r, g in RB.default_guides(df, [s for s in sites if str(s) not in guides]).items():
                guides[str(r)] = g
            # where each site came from, as the reagent name spells it ("r1", "Cterm");
            # hand-added sites have none. The labels count only while they carry this
            # protein's own run path, so a switch cannot stamp the previous protein's onto it
            ranks = (published.get("ranks") or {}) if _same(published.get("path"),
                                                            p["run_json"]) else {}
            site_ranks = {str(s): str(ranks[s]) for s in sites if s in ranks}
            chose = sites != sorted(p.get("sites") or []) or guides != (p.get("guides") or {})
            if not chose and site_ranks == (p.get("site_ranks") or {}):
                return
            p["sites"], p["guides"], p["site_ranks"] = sites, guides, site_ranks
            # ranks arriving late are bookkeeping, not an edit: only a changed choice of
            # site or guide sends an accepted protein back for another look
            if chose and p.get("status") == "accepted":
                p["status"] = "in_progress"
            _save()
            bstate["rev"].set(bstate["rev"].get() + 1)

    @render.ui
    def guides():
        """One guide select per selected site, nearest first; badges never filter."""
        bstate["rev"].get()
        p = _protein(current.get())
        if p is None:
            return ui.div()
        sites = sorted(int(s) for s in (p.get("sites") or []))
        df = _reagents_df(p)
        if df is None:
            return ui.p("No reagents for this protein: " + (p.get("error") or "not designed"),
                        class_="text-danger small")
        if not sites:
            return ui.p("Select sites on the plot; each gets its nearest guide.",
                        class_="text-muted small")
        rows = []
        for r in sites:
            sub = df[df["residue_index"] == r].sort_values("distance", kind="stable")
            if sub.empty:
                rows.append(ui.div(ui.tags.label("{}".format(r), class_="ts-guide-site"),
                                   ui.span("no guide within the arms of this site",
                                           class_="text-danger small"), class_="ts-guide-row"))
                continue
            choices = {str(g): _guide_text(row) for g, row in
                       zip(sub["guide_id"], sub.to_dict("records"))}
            sel = (p.get("guides") or {}).get(str(r)) or next(iter(choices))
            rows.append(ui.div(
                ui.tags.label("{}{}".format(sub.iloc[0]["amino_acid"], r),
                              class_="ts-guide-site"),
                ui.input_select("guide_{}".format(r), None, choices=choices, selected=sel,
                                width="520px"),
                class_="ts-guide-row"))
        return ui.div(*rows)

    @reactive.effect
    def _guide_picks():
        """Store a changed guide choice on the protein's entry."""
        p = _protein(current.get())
        if p is None:
            return
        df = _reagents_df(p)
        if df is None:
            return
        picks = {}
        for r in (p.get("sites") or []):
            try:
                v = input["guide_{}".format(r)]()
            except Exception:
                v = None
            # a select left from the previous protein can share this id; keep only guides
            # that belong to this residue of this protein
            valid = set(df.loc[df["residue_index"] == int(r), "guide_id"])
            if v and v in valid:
                picks[str(r)] = v
        with reactive.isolate():
            guides = dict(p.get("guides") or {})
            new = {**guides, **picks}
            if new != guides:
                p["guides"] = new
                if p.get("status") == "accepted":
                    p["status"] = "in_progress"
                    bstate["rev"].set(bstate["rev"].get() + 1)
                _save()

    # ── counter, accept, save ────────────────────────────────────────────────
    @render.ui
    def counter():
        """'position / total · n accepted'."""
        bstate["rev"].get()
        items = [p["key"] for p in _navigable()]
        data = _data()
        if data is None:
            return ui.span()
        cur = current.get()
        pos = items.index(cur) + 1 if cur in items else 0
        c = BS.counts(data)
        text = "{} / {} · {} accepted".format(pos, len(items), c["accepted"])
        if c["skipped"]:
            text += " · {} skipped".format(c["skipped"])
        return ui.span(text, class_="ts-batch-count")

    @render.ui
    def accept_button():
        """Accept Sites → Next, disabled for a protein without reagents."""
        bstate["rev"].get()
        p = _protein(current.get())
        no_reagents = p is None or p.get("status") == "failed"
        label = "No reagents" if no_reagents and p is not None else "Accept sites → Next"
        return ui.input_action_button("accept", label, class_="btn-sm btn-success",
                                      disabled=no_reagents)

    @render.ui
    def skip_button():
        """Skip → Next (red); Unskip on a protein already skipped."""
        bstate["rev"].get()
        p = _protein(current.get())
        skipped = p is not None and p.get("status") == "skipped"
        return ui.input_action_button("skip", "Unskip" if skipped else "Skip → Next",
                                      class_="btn-sm " + ("btn-outline-danger" if skipped
                                                          else "btn-danger"),
                                      disabled=p is None)

    @render.ui
    def protein_note():
        """Why a protein is red, or its gene-model warning."""
        bstate["rev"].get()
        p = _protein(current.get())
        if p is None:
            return ui.div()
        if p.get("status") == "failed":
            return ui.p(p.get("error") or "No reagents for this protein.",
                        class_="text-danger small")
        return ui.div()

    @reactive.effect
    @reactive.event(input.accept)
    def _accept():
        """Store sites and guides as accepted, then go to the next protein to review."""
        p = _protein(current.get())
        if p is None:
            return
        if not p.get("sites"):
            ui.notification_show("Select at least one site before accepting.", type="warning")
            return
        missing = [s for s in p["sites"] if str(s) not in (p.get("guides") or {})]
        if missing:
            ui.notification_show("No guide for site(s) {}; remove them or pick another "
                                 "site.".format(", ".join(map(str, missing))), type="warning",
                                 duration=8)
            return
        p.pop("skipped_from", None)
        p.update(status="accepted", accepted_at=time.strftime("%Y-%m-%d %H:%M:%S"))
        _save()
        bstate["rev"].set(bstate["rev"].get() + 1)
        _go_next(p["key"])

    def _go_next(after):
        """Move to the next protein still to review, or say there is none."""
        nxt = BS.next_unaccepted(_data(), after=after)
        if nxt is None:
            ui.notification_show("Every protein is accepted or skipped. Output reagents when "
                                 "ready.", type="message", duration=6)
            return
        current.set(nxt)

    @reactive.effect
    @reactive.event(input.skip)
    def _skip():
        """Skip the protein and go to the next one; on a skipped protein, restore its status."""
        p = _protein(current.get())
        if p is None:
            return
        # Unskip: back to whatever it was before the skip
        if p.get("status") == "skipped":
            p["status"] = p.pop("skipped_from", None) or "in_progress"
            p["skipped_at"] = None
            _save()
            bstate["rev"].set(bstate["rev"].get() + 1)
            return
        p.update(skipped_from=p.get("status"), status="skipped",
                 skipped_at=time.strftime("%Y-%m-%d %H:%M:%S"))
        _save()
        bstate["rev"].set(bstate["rev"].get() + 1)
        _go_next(p["key"])

    @reactive.effect
    @reactive.event(input.save)
    def _checkpoint():
        """Save, plus a timestamped copy of the batch file."""
        if bstate["path"] is None:
            return
        copy = BS.checkpoint(bstate["path"], _data())
        ui.notification_show("Saved ({})".format(os.path.basename(copy)), duration=4)

    # ── output ───────────────────────────────────────────────────────────────
    def _zip_name():
        """Download file name."""
        data = _data()
        return "{}_reagents.zip".format(data["name"] if data else "batch")

    @render.download_button(filename=_zip_name, media_type="application/zip")
    async def output_reagents():
        """ZIP of every accepted protein's reagents plus batch-wide files (built off-loop)."""
        data = _data()
        if data is None or not BS.counts(data)["accepted"]:
            ui.notification_show("Accept at least one protein first.", type="warning")
            return
        ui.notification_show("Building reagents for {} protein(s)…".format(
            BS.counts(data)["accepted"]), id="batch_zip", duration=None)
        try:
            payload = await asyncio.to_thread(BO.zip_bytes, data)
        finally:
            ui.notification_remove("batch_zip")
        yield payload

    def _batch_zip_name():
        """Batch ZIP download file name."""
        data = _data()
        return "{}_batch.zip".format(data["name"] if data else "batch")

    @render.download_button(filename=_batch_zip_name, media_type="application/zip")
    async def download_batch():
        """The whole batch as one ZIP of single-app run bundles (scripts/batch_bundle.py)."""
        data = _data()
        if data is None:
            return
        _save()
        ui.notification_show("Packing the batch…", id="batch_pack", duration=None)
        try:
            payload = await asyncio.to_thread(BB.pack_bytes, data)
        finally:
            ui.notification_remove("batch_pack")
        yield payload
