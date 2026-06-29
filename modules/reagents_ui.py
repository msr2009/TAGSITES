"""reagents_ui.py — Shiny UI module for the Design Reagents tab."""

from shiny import ui, module
from modules.ui_helpers import COMPACT_FILE_CSS


_STYLE = """
    .tab-pane { font-size: 0.875rem; }
    .form-label { font-size: 0.8rem; font-weight: 500; margin-bottom: 0.15rem; color: #343a40; }
    .form-control, .form-select { font-size: 0.8rem; padding: 0.2rem 0.45rem; }
    .shiny-input-container { margin-bottom: 0.35rem !important; }

    /* param grid reused from progress tab */
    .param-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 0.3rem 0.75rem;
                  margin-bottom: 0.5rem; }
    .param-label { font-size: 0.75rem; font-weight: 500; color: #6c757d;
                   text-transform: capitalize; }
    .param-value { font-size: 0.8rem; color: #212529; word-break: break-all; }

    /* run header */
    .run-header { font-size: 0.8rem; color: #495057; margin-bottom: 0.3rem; }
    .run-header strong { color: #212529; }

    /* section hints */
    .section-hint { font-size: 0.8rem; color: #6c757d; margin-bottom: 0.3rem; margin-top: 0; }

    /* per-site accordion */
    .accordion-button { font-size: 0.875rem; font-weight: 600; padding: 0.45rem 0.75rem; }
    .accordion-body   { padding: 0.5rem 0.75rem; }
    .card-body        { padding: 0.45rem 0.75rem; }
    .card             { margin-bottom: 0.4rem; }
    .reagents-accordion .accordion-button { background: #f0f4ff; color: #212529; }
    .reagents-accordion .accordion-button:not(.collapsed) { background: #dde4f8; }

    /* guide sub-panels */
    .guide-panel { border: 1px solid #adb5bd; border-radius: 4px;
                   padding: 0.5rem 0.75rem; margin-bottom: 0.5rem; background: #fff; }
    .guide-panel.guide-best { border-color: #4a90d9; background: #f5f9ff; }
    .guide-header { font-size: 0.8rem; font-weight: 600; margin-bottom: 0.4rem;
                    display: flex; align-items: center; gap: 0.5rem; }
    .dist-badge { font-size: 0.7rem; background: #e9ecef; border-radius: 3px;
                  padding: 1px 6px; color: #495057; }

    /* ASCII diagram */
    .guide-diagram pre {
        font-size: 0.72rem; background: #f8f9fa; border: 1px solid #adb5bd;
        border-radius: 3px; padding: 0.4rem 0.6rem; overflow-x: auto;
        white-space: pre; margin-top: 0.4rem; margin-bottom: 0.4rem;
    }

    /* warning / error notices */
    .ts-warn { font-size: 0.8rem; color: #856404; background: #fff3cd;
               border: 1px solid #ffc107; border-radius: 4px; padding: 0.3rem 0.6rem;
               margin-top: 0.3rem; }
    .ts-error { font-size: 0.8rem; color: #842029; background: #f8d7da;
                border: 1px solid #f5c2c7; border-radius: 4px; padding: 0.3rem 0.6rem;
                margin-top: 0.3rem; }

    /* other-guides toggle button */
    .guide-toggle-btn { font-size: 0.78rem; color: #4a90d9; text-decoration: none;
                        margin: 0.3rem 0; display: inline-block; }
    .guide-toggle-btn:hover { color: #2265b0; }

    /* download row */
    .download-row { display: flex; gap: 0.5rem; margin: 0.5rem 0; flex-wrap: wrap; }

"""


@module.ui
def reagents_ui():
    return ui.page_fluid(

        ui.tags.style(_STYLE + COMPACT_FILE_CSS),

        # ── JSON upload card ──────────────────────────────────────────────────
        ui.output_ui("json_upload_card"),

        # ── Run summary ───────────────────────────────────────────────────────
        ui.output_ui("run_header"),

        # ── Output parameters card ────────────────────────────────────────────
        ui.div(
            ui.div(
                ui.div(
                    ui.span("Output parameters"),
                    class_="card-header fw-semibold",
                ),
                ui.div(
                    # Tag / insert sequence — select is static to preserve selection;
                    # server populates choices via ui.update_select on session start
                    ui.input_select("tag_name", "Tag / insert sequence",
                                    choices={"Custom": "Custom…"}, selected="Custom"),
                    ui.output_ui("tag_custom_area"),

                    ui.hr(style="margin: 0.5rem 0;"),

                    # Repair type
                    ui.input_radio_buttons(
                        "repair_type", "Repair type",
                        choices={
                            "ha":       "Homology arm only",
                            "ssodn":    "ssODN",
                            "amplicon": "Amplicon (PCR)",
                            "plasmid":  "Plasmid (coming soon)",
                        },
                        selected="ha",
                        inline=True,
                    ),

                    # Arm length — default set by server when repair_type changes
                    ui.input_numeric(
                        "arm_length", "Homology arm length (bp)",
                        value=500, min=10, max=2000,
                    ),
                    ui.output_ui("arm_length_warning"),

                    # Repair-type-specific options rendered from server
                    ui.output_ui("repair_options"),

                    class_="card-body",
                ),
                class_="card mb-3",
            ),
        ),

        # ── Download row (top) ────────────────────────────────────────────────
        ui.div(
            ui.output_ui("download_row_top"),
            class_="download-row",
        ),

        # ── Per-site cards ────────────────────────────────────────────────────
        ui.output_ui("site_cards"),

        # ── Download row (bottom) ─────────────────────────────────────────────
        ui.div(
            ui.output_ui("download_row_bottom"),
            class_="download-row",
        ),
    )
