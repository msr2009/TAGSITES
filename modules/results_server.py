from shiny import reactive, ui, render, module
import json, os

from utils.results import (
    load_data_from_json,
    load_run_metadata,
    build_plot_payload,
    assign_task_colors,
    residue_colors_for_track,
    residue_colors_gradient,
    residue_colors_for_domains,
    _guess_analysis_type,
)
from config import RESULTS_TYPE_DICT


@module.server
def results_server(input, output, session, shared_json, shared_sites):

    # ── Per-run data ────────────────────────────────────────────────────────────
    aa_data     = reactive.Value()     # continuous scores DataFrame
    range_data  = reactive.Value()     # range annotations DataFrame
    aln_meta    = reactive.Value([])   # list of (path, task_name, params)
    run_name    = reactive.Value()
    run_meta    = reactive.Value({})   # {query_seq, pdb_path, seq_len}
    task_colors = reactive.Value({})   # task_name → hex color (stable across renders)

    # ── Color-by choices (drives the button row above the structure) ────────────
    color_by_choices = reactive.Value({})  # val → label dict, populated on load

    # ── Selection state ─────────────────────────────────────────────────────────
    pending_sites = reactive.Value(set())  # amber — highlighted but not committed

    # ── Load results ────────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.plot_results_button)
    async def load_results():
        """Load analysis results from disk when the user clicks Plot Results."""
        json_content = {}
        try:
            if shared_json.get():
                with open(shared_json.get(), "r") as f:
                    json_content = json.load(f)
            elif input.json_file_input():
                with open(input.json_file_input()[0]["datapath"], "r") as f:
                    json_content = json.load(f)
            else:
                ui.notification_show("No JSON file loaded.", type="warning", duration=4)
                return
        except FileNotFoundError:
            ui.notification_show("JSON file not found. Has the analysis been saved?",
                                 type="error", duration=6)
            return
        except json.JSONDecodeError:
            ui.notification_show("Could not parse JSON file — file may be malformed.",
                                 type="error", duration=6)
            return

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

        # send plot data + sequence to the native canvas renderer
        await _send_plot(aa_df, range_df, meta, json_content["global"]["run_name"])

        # send PDB to 3D viewer if available
        if meta.get("pdb_path"):
            await _send_struct(meta["pdb_path"])

        # reset selection state on new load
        pending_sites.set(set())
        shared_sites.set([])

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
    @reactive.event(input.clear_highlights_button)
    def on_clear():
        """Clear all pending highlights without committing."""
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
        aa_df = aa_data.get()
        range_df = range_data.get()
        meta = run_meta.get()
        seq_len = meta.get("seq_len", 0) if meta else 0
        if not track or track == "(none)":
            await session.send_custom_message("tagsites_set_colors", {"colors": []})
            return
        # parse optional scheme suffix added for pLDDT tracks (e.g. "AF2_plddt:continuous")
        if ":" in track:
            task_name, scheme = track.rsplit(":", 1)
        else:
            task_name, scheme = track, "categorical"

        if track == "__domains__" and range_df is not None and seq_len:
            colors = residue_colors_for_domains(range_df, seq_len)
        elif aa_df is not None and task_name in aa_df.columns:
            hex_color = task_colors.get().get(task_name, "#888888")
            if scheme == "continuous":
                colors = residue_colors_gradient(aa_df, task_name, hex_color)
            else:
                colors = residue_colors_for_track(aa_df, task_name, hex_color)
        else:
            return
        await session.send_custom_message("tagsites_set_colors", {"colors": colors})

    # ── Outputs ──────────────────────────────────────────────────────────────────

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
