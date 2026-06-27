"""progress_ui.py — Shiny UI module for the Progress tab."""

from shiny import ui, module


_STYLE = """
    /* inherit base sizing from setup tab */
    .tab-pane { font-size: 0.875rem; }
    .form-label { font-size: 0.8rem; font-weight: 500; margin-bottom: 0.15rem; color: #343a40; }
    .form-control, .form-select { font-size: 0.8rem; padding: 0.2rem 0.45rem; }
    .shiny-input-container { margin-bottom: 0.35rem !important; }

    /* accordion — same palette as setup */
    .accordion-button { font-size: 0.875rem; font-weight: 600; padding: 0.45rem 0.75rem; }
    .accordion-body   { padding: 0.5rem 0.75rem; }
    .card-body        { padding: 0.45rem 0.75rem; }
    .card             { margin-bottom: 0.4rem; }

    /* task accordion */
    .task-accordion .accordion-button { background: #f8f9fa; color: #212529; }
    .task-accordion .accordion-button:not(.collapsed) { background: #e9ecef; }

    /* badge chips in accordion titles */
    .task-type-badge {
        font-size: 0.7rem; font-weight: 400; color: #6c757d;
        background: #e9ecef; border-radius: 3px; padding: 1px 5px;
        margin-left: 0.4rem; vertical-align: middle;
    }

    /* colored status badge */
    .status-badge {
        font-size: 0.68rem; font-weight: 600; border-radius: 3px;
        padding: 1px 6px; margin-left: 0.5rem; vertical-align: middle;
        letter-spacing: 0.02em; text-transform: uppercase; color: #fff;
    }
    .status-pending { background: #adb5bd; }
    .status-running { background: #fd7e14; animation: pulse 1.2s ease-in-out infinite; }
    .status-success { background: #28a745; }
    .status-failed  { background: #dc3545; }
    @keyframes pulse {
        0%, 100% { opacity: 1; }
        50%       { opacity: 0.6; }
    }

    /* param grid inside accordion body */
    .param-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 0.3rem 0.75rem;
                  margin-bottom: 0.5rem; }
    .param-label { font-size: 0.75rem; font-weight: 500; color: #6c757d;
                   text-transform: capitalize; }
    .param-value { font-size: 0.8rem; color: #212529; word-break: break-all; }

    /* job-id / error lines */
    .job-id-line  { font-size: 0.75rem; font-family: monospace; color: #495057; margin-top: 0.3rem; }
    .error-line   { font-size: 0.78rem; color: #dc3545; margin-top: 0.3rem; }

    /* run header summary */
    .run-header { font-size: 0.8rem; color: #495057; margin-bottom: 0.3rem; }
    .run-header strong { color: #212529; }

    /* section hints */
    .section-hint { font-size: 0.8rem; color: #6c757d; margin-bottom: 0.3rem; margin-top: 0; }

    /* tooltip icon */
    .tip-icon { cursor: help; color: #adb5bd; font-size: 0.75em; margin-left: 3px; }

    /* sticky action bar */
    .save-bar {
        position: sticky; bottom: 0; background: #fff;
        border-top: 2px solid #dee2e6; padding: 0.5rem 1rem;
        margin-top: 0.75rem; z-index: 100;
        display: flex; align-items: center; gap: 0.75rem;
    }

    /* file input sizing */
    .shiny-input-file .btn { font-size: 0.8rem !important; padding: 0.2rem 0.5rem !important; }
    .shiny-input-file .form-control { font-size: 0.8rem; }

    /* current-stage chip in accordion title */
    .stage-chip {
        font-size: 0.68rem; font-weight: 400; color: #fff;
        background: #6c757d; border-radius: 3px; padding: 1px 5px;
        margin-left: 0.3rem; vertical-align: middle; font-style: italic;
    }

    /* per-task log output */
    .task-log {
        width: 100%; margin-top: 0.5rem;
        font-family: monospace; font-size: 0.72rem;
        background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 3px;
        padding: 0.3rem 0.5rem; resize: vertical; color: #212529;
    }
"""


@module.ui
def progress_ui():
    return ui.page_fluid(

        ui.tags.style(_STYLE),

        ui.output_ui("json_upload_card"),

        ui.output_ui("run_header"),

        # ── task accordion ──────────────────────────────────────────────────────
        ui.output_ui("task_cards"),

        # ── sticky action bar ───────────────────────────────────────────────────
        ui.div(
            ui.input_action_button("run_analysis", "▶ Run Analysis",
                                   disabled=True, class_="btn-primary btn-sm"),
            ui.download_button("download_results", "⬇ Download Results",
                               class_="btn-secondary btn-sm"),
            class_="save-bar",
        ),
    )
