"""batch_progress_ui.py — Progress tab of the batch app: export, EBI runs and preload."""

from shiny import module, ui


@module.ui
def batch_progress_ui():
    """Summary line, a Review button, and the per-protein status table."""
    return ui.page_fluid(
        ui.output_ui("summary"),
        ui.output_ui("protein_table"),
    )
