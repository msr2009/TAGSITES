from shiny import reactive, ui, render

import json
import pandas as pd
import plotly.graph_objs as go
from shinywidgets import render_widget, render_plotly


from utils.results import load_data_from_json, plot_results
from utils.results import plot_alignment_matrix, make_alignment_subplots
from config import RESULTS_TYPE_DICT

def results_server(input, output, session, shared_json):

	aa_data = reactive.Value()
	range_data = reactive.Value()
	aln_files = reactive.Value([])
	run_name = reactive.Value()
	subj_fasta = reactive.Value()
	
	@reactive.Effect
	@reactive.event(input.plot_results_button)
	def load_results():
		json_content = {}

		#first check if shared_json
		if shared_json.get():
			json_content = json.load(open(shared_json.get(), "r"))
		else:
			#we need to query if there's an uploaded json in the upload box
			if input.json_file_input():
				json_content = json.load(open(input.json_file_input()[0]["datapath"], "r"))
			else:
				print("upload json in order to plot results")

		aa_data_df, range_data_df, aln_file_list = load_data_from_json(json_content, RESULTS_TYPE_DICT)
		aa_data.set(aa_data_df)
		range_data.set(range_data_df)
		aln_files.set(aln_file_list)
		run_name.set(json_content["global"]["run_name"])
		subj_fasta.set(json_content["global"]["input_file"])

    # Output Plotly plot
	@output
	@render_plotly
	def results_plot():
		# Check if data is available to plot
		if aa_data.get() is not None:
			# Assuming `data.get()` returns a DataFrame, use `plot_results` to generate the plot
			return plot_results(aa_data.get(), subj_fasta.get(), range_data.get(), title=run_name.get())
		else:
			return ui.p("No data to display yet.")

	@output
	@render_plotly
	def alignment_plot():
		#if aln_files() returns a file, then use "plot_alignment_matrix"
		if aln_files.get() is not None:
			# make an alignment plot 
			return make_alignment_subplots(aln_files.get())
		else:
			return ui.p("No alignment data yet.")

