"""batch_nav_ui.py — the batch navigation bar shown above the Results panel."""

from shiny import module, ui

# selectize renderer that colours each protein by its review status. Option labels arrive as
# "status|text" from batch_nav_server; the status picks the colour and is not displayed.
_RENDER = ui.js_eval("""{
  option: function(item, escape) {
    var p = item.label.split('|'); var st = p[0];
    return '<div class="ts-batch-opt ts-batch-' + st + '">' + escape(p.slice(1).join('|'))
           + '</div>';
  },
  item: function(item, escape) {
    var p = item.label.split('|'); var st = p[0];
    return '<div class="ts-batch-opt ts-batch-' + st + '">' + escape(p.slice(1).join('|'))
           + '</div>';
  }
}""")

_CSS = """
.ts-batch-bar { display:flex; flex-wrap:wrap; gap:8px; align-items:center;
  padding:8px 12px; background:#f1f5fb; border:1px solid #cfe0f5; border-radius:6px;
  margin-bottom:10px; }
.ts-batch-bar .form-group { margin-bottom:0; }
.ts-batch-bar .selectize-control { min-width:340px; }
.ts-batch-opt.ts-batch-accepted { color:#198754; font-weight:600; }
.ts-batch-opt.ts-batch-in_progress { color:#0d6efd; }
.ts-batch-opt.ts-batch-ready { color:#495057; }
.ts-batch-opt.ts-batch-failed { color:#dc3545; }
.ts-batch-opt.ts-batch-pending, .ts-batch-opt.ts-batch-running,
.ts-batch-opt.ts-batch-missing { color:#adb5bd; }
.ts-batch-count { font-weight:600; margin:0 6px; }
.ts-guide-row { display:flex; gap:8px; align-items:center; margin:2px 0; }
.ts-guide-row .form-group { margin-bottom:0; }
.ts-guide-row label.ts-guide-site { min-width:70px; font-weight:600; }
#results-json_upload_card { display:none; }
"""


@module.ui
def batch_nav_ui():
    """Prev · protein list · Next · counter · Accept → Next · Save · Batch ZIP · Output."""
    return ui.div(
        ui.tags.style(_CSS),
        ui.div(
            ui.input_action_button("prev", "◀ Prev", class_="btn-sm btn-outline-secondary"),
            ui.input_selectize("protein", None, choices={}, width="360px",
                               options={"render": _RENDER}),
            ui.input_action_button("next", "Next ▶", class_="btn-sm btn-outline-secondary"),
            ui.output_ui("counter", inline=True),
            ui.output_ui("accept_button", inline=True),
            ui.input_action_button("save", "Save batch", class_="btn-sm btn-outline-primary"),
            ui.download_button("download_batch", "⬇ Batch (.zip)",
                               class_="btn-sm btn-outline-primary"),
            ui.download_button("output_reagents", "⬇ Output reagents (.zip)",
                               class_="btn-sm btn-primary"),
            class_="ts-batch-bar",
        ),
        ui.output_ui("protein_note"),
    )


@module.ui
def batch_guides_ui():
    """Per-site guide picker; ui_batch.py places it inside the Results tab (results_ui's
    `between` slot), below the result panes and above Advanced / Rescoring."""
    return ui.accordion(
        ui.accordion_panel("Guide per site (nearest by default)", ui.output_ui("guides")),
        id="guide_accordion", open=True,
    )
