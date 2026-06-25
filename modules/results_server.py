from shiny import reactive, ui, render, module
import json, os

from shinywidgets import render_plotly

from utils.results import load_data_from_json, plot_results
from config import RESULTS_TYPE_DICT


@module.server
def results_server(input, output, session, shared_json):

    aa_data = reactive.Value()
    range_data = reactive.Value()
    aln_files = reactive.Value([])
    run_name = reactive.Value()
    subj_fasta = reactive.Value()

    @reactive.effect
    @reactive.event(input.plot_results_button)
    def load_results():
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

        aa_data_df, range_data_df, aln_file_list = load_data_from_json(json_content, RESULTS_TYPE_DICT)
        aa_data.set(aa_data_df)
        range_data.set(range_data_df)
        aln_files.set(aln_file_list)
        run_name.set(json_content["global"]["run_name"])
        subj_fasta.set(json_content["global"]["input_file"])

    @reactive.calc
    def results_figure():
        """Build the Plotly figure — cached until data changes, avoiding repeated disk reads."""
        if run_name.get() is None:
            return None
        return plot_results(aa_data.get(), subj_fasta.get(), range_data.get(),
                            title=run_name.get())

    @output
    @render_plotly
    def results_plot():
        """Render the cached results figure."""
        return results_figure()

    @render.ui
    def alignments_container():
        """Render alignment SVGs as inline HTML."""
        if not aln_files.get():
            return ui.div()
        images = []
        for a in aln_files.get():
            svg_path = a.removesuffix("aln") + "svg"
            if os.path.exists(svg_path):
                with open(svg_path) as f:
                    images.append(ui.HTML(f.read()))
            else:
                images.append(ui.p(f"Alignment image not found: {os.path.basename(svg_path)}"))
        return ui.div(*images)
