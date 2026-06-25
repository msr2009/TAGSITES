from shiny import reactive, render, ui, module
import json


@module.server
def reagents_server(input, output, session, shared_json):

    run_name = reactive.Value()
    subj_fasta = reactive.Value()

    @render.ui
    def show_or_upload_json_reagents():
        """Placeholder output for the reagents panel (to be built out in Phase 3)."""
        return ui.div("Reagents design coming soon.")
