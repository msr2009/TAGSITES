from shiny import reactive, ui, render, module
import base64, json, os
from modules.json_card import json_upload_card as _json_upload_card

from utils.results import (
    load_data_from_json,
    load_run_metadata,
    build_plot_payload,
    assign_task_colors,
    residue_colors_for_track,
    residue_colors_gradient,
    residue_colors_for_annotations,
    residue_colors_for_phobius,
    residue_colors_for_isoforms,
    residue_colors_jet,
    _guess_analysis_type,
    _pick_gradient_cmap,
    _VIRIDIS, _PLASMA, _COOL, _PLDDT_GRADIENT,
)
from config import RESULTS_TYPE_DICT, DOMAIN_SOURCE_COLORS

# pLDDT 4-band legend items (values stored 0-1 after /100 at load time)
_PLDDT_LEGEND = [
    {"color": "#0053D6", "label": "Very high (≥0.90)"},
    {"color": "#65CBF3", "label": "Confident (0.70–0.90)"},
    {"color": "#FFDB13", "label": "Low (0.50–0.70)"},
    {"color": "#FF7D45", "label": "Very low (<0.50)"},
]


def _build_colors_and_legend(track, task_name, scheme, aa_df, range_df, seq_len, hex_color):
    """Single entry point: compute per-residue structure colors AND the matching legend.

    Returns (colors, legend) where legend is a dict ready for JSON serialization,
    or (None, None) when the track cannot be resolved.
    """
    # ── Annotation coloring (domains / modifications / Phobius) ─────────────────
    if track == "__domains__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_annotations(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    if track == "__phobius__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_phobius(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    if track == "__isoforms__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_isoforms(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    # ── Continuous track coloring ────────────────────────────────────────────────
    if aa_df is None or task_name not in aa_df.columns:
        return None, None

    if scheme == "continuous":
        colors = residue_colors_gradient(aa_df, task_name, hex_color)
    else:
        colors = residue_colors_for_track(aa_df, task_name, hex_color)

    # pLDDT categorical (4-band) — excludes SASA and forced-continuous variants
    if (_guess_analysis_type(task_name) == "plddt"
            and not task_name.endswith("_sasa")
            and scheme != "continuous"):
        legend = {"type": "categorical", "items": _PLDDT_LEGEND}
    else:
        _CMAP_NAMES = {
            id(_VIRIDIS): "viridis", id(_PLASMA): "plasma",
            id(_COOL): "cool",       id(_PLDDT_GRADIENT): "plddt",
        }
        cmap_name = _CMAP_NAMES.get(id(_pick_gradient_cmap(task_name)), "viridis")
        legend = {"type": "gradient", "label": task_name, "vmin": 0, "vmax": 1,
                  "cmap": cmap_name}

    return colors, legend


@module.server
def results_server(input, output, session, shared_json, shared_sites, shared_results_trigger=None):

    # ── Per-run data ────────────────────────────────────────────────────────────
    aa_data     = reactive.Value()     # continuous scores DataFrame
    range_data  = reactive.Value()     # range annotations DataFrame
    aln_meta    = reactive.Value([])   # list of (path, task_name, params)
    run_name    = reactive.Value(None)
    run_meta    = reactive.Value({})   # {query_seq, pdb_path, seq_len}
    task_colors = reactive.Value({})   # task_name → hex color (stable across renders)

    # ── Color-by choices (drives the button row above the structure) ────────────
    color_by_choices = reactive.Value({})  # val → label dict, populated on load

    # canonical on-disk path (not the Shiny upload temp path) — used for saves
    json_path = reactive.Value(None)

    # ── Selection state ─────────────────────────────────────────────────────────
    pending_sites = reactive.Value(set())  # amber — highlighted but not committed

    # ── Load results ────────────────────────────────────────────────────────────

    async def _do_load(json_content):
        """Parse a loaded JSON dict and push all results to the UI."""
        aa_df, range_df, alns = load_data_from_json(json_content, RESULTS_TYPE_DICT)
        meta = load_run_metadata(json_content)
        aa_data.set(aa_df)
        range_data.set(range_df)
        aln_meta.set(alns)
        run_name.set(json_content["global"]["run_name"])
        run_meta.set(meta)
        task_colors.set(assign_task_colors(aa_df) if aa_df is not None else {})

        # build color-by button choices
        # pLDDT tracks get two buttons: categorical (4-band) and continuous (gradient)
        choices = {"(none)": "N→C"}
        if aa_df is not None:
            for col in aa_df.columns[1:]:
                if _guess_analysis_type(col) == "plddt" and not col.endswith("_sasa"):
                    choices[col + ":categorical"] = f"{col} (4-band)"
                    choices[col + ":continuous"]  = f"{col} (gradient)"
                elif col.endswith("_sasa"):
                    choices[col] = "Solv Access"
                else:
                    choices[col] = col
        if range_df is not None and not range_df.empty:
            if not range_df[range_df["source"].isin({"Pfam", "modification"})].empty:
                choices["__domains__"] = "Domains"
            if not range_df[range_df["source"] == "Phobius"].empty:
                choices["__phobius__"] = "Phobius"
            if not range_df[range_df["source"] == "isoforms"].empty:
                choices["__isoforms__"] = "Isoforms"
        color_by_choices.set(choices)
        ui.update_select("color_by", choices=choices, selected="(none)")

        await _send_plot(aa_df, range_df, meta, json_content["global"]["run_name"])
        # always reset structure colors on load — don't rely on color_by reactive firing
        # (it won't fire if color_by was already "(none)" before this load)
        jet_colors = residue_colors_jet(meta.get("seq_len", 0))
        await session.send_custom_message("tagsites_set_colors",
                                          {"colors": jet_colors, "legend": {"type": "rainbow"}})
        if meta.get("pdb_path"):
            await _send_struct(meta["pdb_path"])
        # record the canonical on-disk path so _save_sites can write back
        wd = json_content["global"].get("working_dir", "")
        rn_str = json_content["global"].get("run_name", "")
        canon = os.path.join(wd, rn_str + ".json") if wd and rn_str else None
        json_path.set(canon)

        # restore previously saved site selection (or start fresh)
        saved = json_content["global"].get("selected_sites", [])
        pending_sites.set(set())
        shared_sites.set(sorted(int(s) for s in saved))

    @reactive.effect
    @reactive.event(input.json_file_input)
    def load_results():
        """Propagate a Results-tab JSON upload to shared state.

        Validates the file first, then sets shared_json so that auto_load_results
        (and the Progress / Reagents tabs) all react to the same path — same path
        as when JSON is loaded through any other tab.
        """
        info = input.json_file_input()
        if not info:
            return
        path = info[0]["datapath"]
        try:
            with open(path) as f:
                json.load(f)
        except FileNotFoundError:
            ui.notification_show("Uploaded file not found.", type="error", duration=6)
            return
        except json.JSONDecodeError:
            ui.notification_show("Could not parse JSON — file may be malformed.",
                                 type="error", duration=6)
            return
        shared_json.set(path)

    @reactive.effect
    async def auto_load_results():
        """Auto-plot results whenever the JSON path changes or a run completes."""
        path = shared_json.get()
        # take a dependency on the completion trigger so we reload after tasks finish
        # even when the JSON path hasn't changed
        if shared_results_trigger is not None:
            shared_results_trigger.get()
        if not path:
            return
        try:
            with open(path) as f:
                json_content = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return
        await _do_load(json_content)

    @reactive.effect
    def _save_sites():
        """Write selected_sites back to the run JSON whenever they change."""
        sites = shared_sites.get()
        path  = json_path.get()
        if not path or not os.path.exists(path):
            return
        try:
            with open(path) as f:
                j = json.load(f)
            j.setdefault("global", {})["selected_sites"] = sites
            with open(path, "w") as f:
                json.dump(j, f, indent=2)
        except Exception:
            pass

    def _input_names():
        """Build the dict of namespaced Shiny input IDs to send to JS."""
        return {
            "residue_click":    session.ns("residue_click"),
            "residue_dblclick": session.ns("residue_dblclick"),
            "struct_click":     session.ns("struct_click"),
            "struct_dblclick":  session.ns("struct_dblclick"),
            "remove_site":      session.ns("remove_site"),
        }

    async def _send_plot(aa_df, range_df, meta, title):
        """Build and send all plot data to the native canvas renderer."""
        payload = build_plot_payload(aa_df, range_df, title=title)
        payload["seq"] = meta.get("query_seq", "")
        payload["inputs"] = _input_names()
        await session.send_custom_message("tagsites_set_plot", payload)

    async def _send_struct(pdb_path):
        """Read PDB file and send to 3Dmol viewer."""
        try:
            with open(pdb_path) as f:
                pdb_str = f.read()
        except OSError:
            return
        await session.send_custom_message("tagsites_init_struct", {
            "pdb": pdb_str,
            "inputs": _input_names(),
        })

    # ── Click/selection handlers ─────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.residue_click)
    def on_residue_click():
        """Toggle a residue in/out of the pending set."""
        pos = input.residue_click()
        if pos is None:
            return
        pos = int(pos)
        p = set(pending_sites.get())
        if pos in p:
            p.discard(pos)
        else:
            if pos not in set(shared_sites.get()):
                p.add(pos)
        pending_sites.set(p)

    @reactive.effect
    @reactive.event(input.residue_dblclick)
    def on_residue_dblclick():
        """Commit a residue on double-click; remove it if already committed."""
        pos = input.residue_dblclick()
        if pos is None:
            return
        pos = int(pos)
        sites = list(shared_sites.get())
        if pos in sites:
            sites.remove(pos)
            shared_sites.set(sites)
        else:
            p = set(pending_sites.get())
            p.discard(pos)
            pending_sites.set(p)
            sites.append(pos)
            sites.sort()
            shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.struct_click)
    def on_struct_click():
        """Toggle a residue from 3D structure click."""
        pos = input.struct_click()
        if pos is None:
            return
        pos = int(pos)
        p = set(pending_sites.get())
        if pos in p:
            p.discard(pos)
        else:
            if pos not in set(shared_sites.get()):
                p.add(pos)
        pending_sites.set(p)

    @reactive.effect
    @reactive.event(input.struct_dblclick)
    def on_struct_dblclick():
        """Commit a residue from 3D double-click; remove it if already committed."""
        pos = input.struct_dblclick()
        if pos is None:
            return
        pos = int(pos)
        sites = list(shared_sites.get())
        if pos in sites:
            sites.remove(pos)
            shared_sites.set(sites)
        else:
            p = set(pending_sites.get())
            p.discard(pos)
            pending_sites.set(p)
            sites.append(pos)
            sites.sort()
            shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.add_highlighted_button)
    def on_add_highlighted():
        """Commit all pending residues at once."""
        new_sites = set(pending_sites.get())
        sites = sorted(set(shared_sites.get()) | new_sites)
        shared_sites.set(sites)
        pending_sites.set(set())

    @reactive.effect
    @reactive.event(input.add_nterm_button)
    def on_add_nterm():
        """Add position 1 (N-terminus) to committed sites."""
        sites = sorted(set(shared_sites.get()) | {1})
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.add_cterm_button)
    def on_add_cterm():
        """Add the last residue (C-terminus) to committed sites."""
        meta = run_meta.get() or {}
        seq_len = meta.get("seq_len", 0)
        if not seq_len:
            return
        sites = sorted(set(shared_sites.get()) | {seq_len})
        shared_sites.set(sites)

    # TODO: implement suggested-site algorithm (score-based auto-selection)
    # Button is disabled in results_ui.py until this is ready.
    # When implementing: enable the button, fill this handler with the selection logic,
    # and update the button class from btn-secondary back to btn-outline-success.
    @reactive.effect
    @reactive.event(input.add_suggested_button)
    def on_add_suggested():
        pass

    @reactive.effect
    @reactive.event(input.clear_highlights_button)
    def on_clear():
        """Clear all committed and pending tag sites."""
        shared_sites.set([])
        pending_sites.set(set())

    @reactive.effect
    @reactive.event(input.remove_site)
    def on_remove_site():
        """Remove a committed site via its chip × button."""
        pos = input.remove_site()
        if pos is None:
            return
        pos = int(pos)
        sites = [s for s in shared_sites.get() if s != pos]
        shared_sites.set(sites)

    # ── Sync JS highlight state whenever Python state changes ───────────────────

    @reactive.effect
    async def sync_js_states():
        """Push pending + committed state to both JS surfaces."""
        await session.send_custom_message("tagsites_set_states", {
            "pending":   sorted(pending_sites.get()),
            "committed": sorted(shared_sites.get()),
        })

    # ── Color-by handler ─────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.color_by)
    async def sync_js_colors():
        """Compute per-residue colors for the selected track and send to 3D viewer."""
        track = input.color_by()
        if not track or track == "(none)":
            meta = run_meta.get()
            seq_len = meta.get("seq_len", 0) if meta else 0
            colors = residue_colors_jet(seq_len) if seq_len else []
            await session.send_custom_message("tagsites_set_colors",
                                              {"colors": colors,
                                               "legend": {"type": "rainbow"}})
            return

        task_name, scheme = track.rsplit(":", 1) if ":" in track else (track, "categorical")
        meta = run_meta.get()
        seq_len = meta.get("seq_len", 0) if meta else 0
        hex_color = task_colors.get().get(task_name, "#888888")

        colors, legend = _build_colors_and_legend(
            track, task_name, scheme,
            aa_data.get(), range_data.get(), seq_len, hex_color,
        )
        if colors is None:
            return

        await session.send_custom_message("tagsites_set_colors",
                                          {"colors": colors, "legend": legend})

    # ── Outputs ──────────────────────────────────────────────────────────────────

    @render.ui
    def json_upload_card():
        """Render the JSON upload card; collapsed and warning-free when data is loaded."""
        return _json_upload_card(
            "json_file_input", "ts-json-body", "Upload results JSON",
            run_name.get() is not None,
        )

    @render.ui
    def color_buttons_ui():
        """Render color-by buttons above the structure viewer."""
        choices = color_by_choices.get()
        if not choices:
            return ui.div()
        color_by_input_id = session.ns("color_by")
        _isoforms_tooltip = (
            "Isoform regions inferred from same-organism BLAST hits (UniProt isoforms). "
            "Segments reflect protein sequence coverage, not genomic exon boundaries."
        )
        buttons = [
            ui.tags.button(
                label,
                onclick=f"tsSetColorBy(this, '{val}', '{color_by_input_id}')",
                class_="btn btn-sm ts-colorby-btn" + (" ts-colorby-active" if val == "(none)" else ""),
                title=_isoforms_tooltip if val == "__isoforms__" else None,
            )
            for val, label in choices.items()
        ]
        return ui.div(*buttons, class_="ts-colorby-row")

    @render.ui
    def chosen_sites_display():
        """Render chosen-site chips with × remove buttons."""
        sites = shared_sites.get()
        if not sites:
            return ui.div(
                ui.span("No sites chosen yet — click residues in the plot or structure to select.",
                        class_="ts-empty-hint"),
                id="ts-sites-box",
            )
        meta = run_meta.get() or {}
        seq = meta.get("query_seq", "")
        remove_id = session.ns("remove_site")
        chips = []
        for pos in sites:
            aa = seq[pos - 1] if seq and 0 < pos <= len(seq) else "?"
            label = f"{aa}{pos}"
            chip = ui.span(
                label,
                ui.tags.button(
                    "×",
                    class_="ts-chip-remove",
                    onclick=f"tsRemoveSite({pos}, '{remove_id}')",
                    title=f"Remove {label}",
                ),
                class_="ts-site-chip",
            )
            chips.append(chip)
        return ui.div(*chips, id="ts-sites-box")

    @render.ui
    def alignments_container():
        """Render blast alignment PNGs as collapsible accordion panes."""
        alns = aln_meta.get()
        if not alns:
            return ui.div()

        panels = []
        for aln_path, task_name, params in alns:
            png_path = aln_path.removesuffix(".aln") + ".png"

            param_rows = [
                ui.tags.tr(ui.tags.td(k), ui.tags.td(str(v)))
                for k, v in params.items()
            ]
            param_table = ui.tags.table(
                ui.tags.tbody(*param_rows),
                class_="ts-aln-params",
            ) if param_rows else ui.div()

            if os.path.exists(png_path):
                with open(png_path, "rb") as f:
                    b64 = base64.b64encode(f.read()).decode("ascii")
                img = ui.tags.img(
                    src=f"data:image/png;base64,{b64}",
                    style="height:auto; display:block;",
                )
                body = ui.div(
                    param_table,
                    ui.div(img, class_="ts-aln-svg-wrap",
                           style="overflow-x:auto; overflow-y:hidden;"),
                )
            else:
                body = ui.div(
                    param_table,
                    ui.p(f"Alignment image not found: {os.path.basename(png_path)}",
                         style="color:#c00;"),
                )

            panels.append(ui.accordion_panel(task_name, body))

        if not panels:
            return ui.div()
        return ui.div(
            ui.h4("Alignments"),
            ui.accordion(*panels, id="ts_aln_accordion", open=False),
        )
