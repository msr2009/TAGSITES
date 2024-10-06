from shiny import reactive, render, ui
import json, subprocess
import pandas as pd
import plotly.graph_objs as go
from shinywidgets import render_widget, render_plotly

from utils.results import load_data_from_json, plot_results
from utils.results import plot_alignment_matrix, make_alignment_subplots
from config import RESULTS_TYPE_DICT

def reagents_server(input, output, session, shared_json):

	run_name = reactive.Value()
	subj_fasta = reactive.Value()
	genewise_data = reactive.Value()

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

		run_name.set(json_content["global"]["run_name"])
		subj_fasta.set(json_content["global"]["input_file"])
