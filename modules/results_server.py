from shiny import reactive, ui, render, module
import json, os

from utils.results import (
    load_data_from_json,
    load_run_metadata,
    build_plot_payload,
    assign_task_colors,
    residue_colors_for_track,
    residue_colors_gradient,
    residue_colors_for_annotations,
    _guess_analysis_type,
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
        legend = {"type": "gradient", "label": task_name, "vmin": 0, "vmax": 1}

    return colors, legend


@module.server
def results_server(input, output, session, shared_json, shared_sites):

    # ── Per-run data ────────────────────────────────────────────────────────────
    aa_data     = reactive.Value()     # continuous scores DataFrame
    range_data  = reactive.Value()     # range annotations DataFrame
    aln_meta    = reactive.Value([])   # list of (path, task_name, params)
    run_name    = reactive.Value(None)
    run_meta    = reactive.Value({})   # {query_seq, pdb_path, seq_len}
    task_colors = reactive.Value({})   # task_name → hex color (stable across renders)

    # ── Color-by choices (drives the button row above the structure) ────────────
    color_by_choices = reactive.Value({})  # val → label dict, populated on load

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
        choices = {"(none)": "None"}
        if aa_df is not None:
            for col in aa_df.columns[1:]:
                if _guess_analysis_type(col) == "plddt" and not col.endswith("_sasa"):
                    choices[col + ":categorical"] = f"{col} (4-band)"
                    choices[col + ":continuous"]  = f"{col} (gradient)"
                else:
                    choices[col] = col
        if range_df is not None and not range_df.empty:
            choices["__domains__"] = "Domains"
        color_by_choices.set(choices)
        ui.update_select("color_by", choices=choices, selected="(none)")

        await _send_plot(aa_df, range_df, meta, json_content["global"]["run_name"])
        # always reset structure colors on load — don't rely on color_by reactive firing
        # (it won't fire if color_by was already "(none)" before this load)
        await session.send_custom_message("tagsites_set_colors", {"colors": [], "legend": None})
        if meta.get("pdb_path"):
            await _send_struct(meta["pdb_path"])
        pending_sites.set(set())
        shared_sites.set([])

    @reactive.effect
    @reactive.event(input.json_file_input)
    async def load_results():
        """Auto-plot when the user selects an uploaded JSON file."""
        if not input.json_file_input():
            return
        try:
            with open(input.json_file_input()[0]["datapath"]) as f:
                json_content = json.load(f)
        except FileNotFoundError:
            ui.notification_show("Uploaded file not found.", type="error", duration=6)
            return
        except json.JSONDecodeError:
            ui.notification_show("Could not parse JSON — file may be malformed.",
                                 type="error", duration=6)
            return
        await _do_load(json_content)

    @reactive.effect
    async def auto_load_results():
        """Auto-plot results whenever the pipeline writes a new shared JSON path."""
        path = shared_json.get()
        if not path:
            return
        try:
            with open(path) as f:
                json_content = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return
        await _do_load(json_content)

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
    @reactive.event(input.add_suggested_button)
    def on_add_suggested():
        """Placeholder — suggested site logic to be implemented."""
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
            await session.send_custom_message("tagsites_set_colors",
                                              {"colors": [], "legend": None})
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
    def json_card_warning():
        """Show a warning badge in the upload card header when no data is loaded."""
        if run_name.get() is not None:
            return ui.span()
        return ui.span("⚠ no current JSON",
                       style="color:#dc3545; font-size:0.8em; font-style:italic;")

    @render.ui
    def color_buttons_ui():
        """Render color-by buttons above the structure viewer."""
        choices = color_by_choices.get()
        if not choices:
            return ui.div()
        color_by_input_id = session.ns("color_by")
        buttons = [
            ui.tags.button(
                label,
                onclick=f"tsSetColorBy(this, '{val}', '{color_by_input_id}')",
                class_="btn btn-sm ts-colorby-btn" + (" ts-colorby-active" if val == "(none)" else ""),
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
        """Render blast alignment SVGs as collapsible accordion panes."""
        alns = aln_meta.get()
        if not alns:
            return ui.div()

        panels = []
        for aln_path, task_name, params in alns:
            svg_path = aln_path.removesuffix("aln") + "svg"

            param_rows = [
                ui.tags.tr(ui.tags.td(k), ui.tags.td(str(v)))
                for k, v in params.items()
            ]
            param_table = ui.tags.table(
                ui.tags.tbody(*param_rows),
                class_="ts-aln-params",
            ) if param_rows else ui.div()

            if os.path.exists(svg_path):
                with open(svg_path) as f:
                    svg_content = ui.HTML(f.read())
                body = ui.div(
                    param_table,
                    ui.div(svg_content, class_="ts-aln-svg-wrap"),
                )
            else:
                body = ui.div(
                    param_table,
                    ui.p(f"Alignment image not found: {os.path.basename(svg_path)}",
                         style="color:#c00;"),
                )

            panels.append(ui.accordion_panel(task_name, body))

        if not panels:
            return ui.div()
        return ui.div(
            ui.h4("Alignments"),
            ui.accordion(*panels, id="ts_aln_accordion", open=False),
        )
