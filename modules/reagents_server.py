from shiny import reactive, render, ui, module
import json


@module.server
def reagents_server(input, output, session, shared_json, shared_sites):

    run_name = reactive.Value()
    subj_fasta = reactive.Value()

    @render.ui
    def show_or_upload_json_reagents():
        """Show chosen tag sites from the Results page; reagent design coming in Phase 3."""
        sites = shared_sites.get()
        if not sites:
            return ui.div(
                ui.p("No tag sites selected yet."),
                ui.p("Go to the Results tab, explore the data, "
                     "and click residues to choose insertion sites.",
                     style="color:#555;"),
            )
        site_list = ", ".join(str(s) for s in sites)
        return ui.div(
            ui.p(ui.strong(f"Chosen tag sites ({len(sites)}): "), site_list),
            ui.hr(),
            ui.p("Reagent design coming soon.", style="color:#888; font-style:italic;"),
        )
