from shiny import ui, module


@module.ui
def results_ui():
    return ui.page_fluid(
        # tagsites assets (3Dmol must load before tagsites.js)
        ui.tags.script(src="3Dmol-min.js"),
        ui.tags.link(rel="stylesheet", href="tagsites.css"),
        ui.tags.script(src="tagsites.js"),

        ui.output_ui("json_upload_card"),

        ui.h4("Analysis Results"),

        # ── Main content: plot (left) + structure (right) ─────────────────────
        ui.layout_columns(
            # left column: toolbar + native canvas plot
            ui.div(
                # pan / set-region / zoom-in / zoom-out / reset toolbar
                ui.div(
                    ui.tags.button("✥", onclick="tsSetMode('pan')",
                                   class_="btn btn-sm ts-mode-btn", title="Pan"),
                    ui.tags.button("⬚", onclick="tsSetMode('zoom')",
                                   class_="btn btn-sm ts-mode-btn ts-mode-active", title="Set region"),
                    ui.tags.button("⊕", onclick="tsZoomIn()",
                                   class_="btn btn-sm ts-mode-btn", title="Zoom in"),
                    ui.tags.button("⊖", onclick="tsZoomOut()",
                                   class_="btn btn-sm ts-mode-btn", title="Zoom out"),
                    ui.tags.button("⟳", onclick="tsResetZoom()",
                                   class_="btn btn-sm btn-outline-secondary", title="Reset zoom"),
                    class_="ts-plot-toolbar",
                ),
                ui.div(
                    ui.tags.canvas(id="ts_plot_div"),
                    id="ts-plot-wrap",
                ),
                class_="ts-main-left",
            ),
            # right column: 3D structure pane (hidden until PDB available)
            ui.div(
                ui.div(
                    # hidden select registers the color_by input with Shiny;
                    # the visible button row drives it via Shiny.setInputValue
                    ui.div(
                        ui.input_select("color_by", "", choices={"(none)": "N→C"}),
                        style="display:none; position:absolute; pointer-events:none;",
                    ),
                    # color-by buttons above the viewer
                    ui.output_ui("color_buttons_ui"),
                    # the viewer itself
                    ui.div(id="ts-viewer-container"),
                    # color legend (populated by JS on color change)
                    ui.div(id="ts-struct-legend"),
                    id="ts-struct-pane",
                ),
                class_="ts-struct-col",
            ),
            col_widths=(7, 5),
        ),

        # ── Chosen tag sites: chips | Add highlighted | N-term | C-term | Add suggested | Clear ──
        ui.div(
            ui.output_ui("chosen_sites_display"),
            ui.input_action_button("add_highlighted_button",
                                   "Add highlighted",
                                   class_="btn btn-success btn-sm ts-sites-btn"),
            ui.input_action_button("add_nterm_button",
                                   "Add N-term",
                                   class_="btn btn-outline-primary btn-sm ts-sites-btn"),
            ui.input_action_button("add_cterm_button",
                                   "Add C-term",
                                   class_="btn btn-outline-danger btn-sm ts-sites-btn"),
            ui.input_action_button("add_suggested_button",
                                   "Add suggested",
                                   class_="btn btn-outline-success btn-sm ts-sites-btn"),
            ui.input_action_button("clear_highlights_button",
                                   "Clear",
                                   class_="btn btn-secondary btn-sm ts-sites-btn"),
            class_="ts-sites-row",
        ),
        ui.br(),

        # ── Alignment accordion ───────────────────────────────────────────────
        ui.output_ui("alignments_container"),
    )
