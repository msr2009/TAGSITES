from shiny import ui
from shinywidgets import output_widget

#Page 3: Results
def results_ui():
	return ui.page_fluid(
	
		ui.h2("Results"),
	
		ui.input_file("json_file_input", "Upload a JSON file for results"),
		ui.input_action_button("plot_results_button", "Plot Results"),
		
		ui.hr(),
		#output plotly plot
		output_widget("results_plot"),
		ui.hr(),
		output_widget("alignment_plot")
	)
