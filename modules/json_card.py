"""json_card.py — shared JSON-upload card used by Progress, Results, and Reagents tabs."""

from shiny import ui


def json_upload_card(input_id, body_id, title, has_json):
    """Return a collapsible Bootstrap card containing a JSON file-upload input."""
    warning = ui.span() if has_json else ui.span(
        "⚠ no current JSON",
        style="color:#dc3545; font-size:0.8em; font-style:italic;",
    )
    return ui.div(
        ui.div(
            ui.div(
                ui.span(title),
                warning,
                class_="card-header d-flex justify-content-between align-items-center",
                style="cursor:pointer;",
                **{"data-bs-toggle": "collapse",
                   "data-bs-target": "#{}".format(body_id),
                   "aria-expanded": "false" if has_json else "true"},
            ),
            ui.div(
                ui.div(
                    ui.input_file(input_id, None, accept=[".json"]),
                    class_="card-body py-2",
                ),
                id=body_id,
                class_="collapse" if has_json else "collapse show",
            ),
            class_="card",
        ),
        class_="mb-2",
    )
