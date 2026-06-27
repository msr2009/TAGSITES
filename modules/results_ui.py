from shiny import ui, module


@module.ui
def results_ui():
    return ui.page_fluid(
        # tagsites assets (3Dmol must load before tagsites.js)
        ui.tags.script(src="3Dmol-min.js"),
        ui.tags.link(rel="stylesheet", href="tagsites.css"),
        ui.tags.script(src="tagsites.js"),

        ui.h2("Results"),

        # file upload + plot button — compact inline row
        ui.div(
            ui.input_file("json_file_input", "Upload results JSON", accept=[".json"]),
            ui.input_action_button("plot_results_button", "Plot Results",
                                   class_="btn-primary"),
            class_="ts-controls-row",
        ),
        ui.hr(),

        # ── Main content: plot (left) + structure (right) ─────────────────────
        ui.layout_columns(
            # left column: toolbar + native canvas plot
            ui.div(
                # zoom / pan / reset toolbar
                ui.div(
                    ui.tags.button("⊕ Zoom", onclick="tsSetMode('zoom')",
                                   class_="btn btn-sm ts-mode-btn ts-mode-active"),
                    ui.tags.button("✥ Pan",  onclick="tsSetMode('pan')",
                                   class_="btn btn-sm ts-mode-btn"),
                    ui.tags.button("⟳ Reset", onclick="tsResetZoom()",
                                   class_="btn btn-sm btn-outline-secondary"),
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
                        ui.input_select("color_by", "", choices={"(none)": "None"}),
                        style="display:none; position:absolute; pointer-events:none;",
                    ),
                    # color-by buttons above the viewer
                    ui.output_ui("color_buttons_ui"),
                    # the viewer itself
                    ui.div(id="ts-viewer-container"),
                    # background selector + download below the viewer
                    ui.div(
                        ui.tags.button("White",
                                       onclick="tsSetBackground('white')",
                                       class_="btn btn-sm ts-bg-btn"),
                        ui.tags.button("Black",
                                       onclick="tsSetBackground('black')",
                                       class_="btn btn-sm ts-bg-btn"),
                        ui.tags.button("Transparent",
                                       onclick="tsSetBackground('transparent')",
                                       class_="btn btn-sm ts-bg-btn ts-bg-btn-active",
                                       id="ts-bg-transparent"),
                        ui.tags.button("⬇ Download",
                                       onclick="tsDownloadStructure()",
                                       class_="btn btn-sm btn-outline-primary ts-dl-btn"),
                        class_="ts-struct-bottom",
                    ),
                    id="ts-struct-pane",
                ),
            ),
            col_widths=(8, 4),
        ),

        # ── Chosen tag sites: chips | Add highlighted | Clear ─────────────────
        ui.div(
            ui.output_ui("chosen_sites_display"),
            ui.input_action_button("add_highlighted_button",
                                   "Add highlighted",
                                   class_="btn btn-success btn-sm ts-sites-btn"),
            ui.input_action_button("clear_highlights_button",
                                   "Clear",
                                   class_="btn btn-secondary btn-sm ts-sites-btn"),
            class_="ts-sites-row",
        ),
        ui.br(),

        # ── Alignment accordion ───────────────────────────────────────────────
        ui.output_ui("alignments_container"),
    )
