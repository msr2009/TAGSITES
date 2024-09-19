from shiny import ui

#Page 2: Progress
def progress_ui():
	return ui.page_fluid(
		ui.h2("Progress"),
		ui.hr(),

		#load analyses json
		ui.h3("Load Analysis"),
		#load json button
#		ui.input_file("analysis_json", "Load Analysis JSON file",
#				accept=['.txt']),

		ui.output_ui("show_or_upload_json"),
		ui.hr(),

		#maybe something here to show the steps of the analysis?

		#run analysis button
		ui.input_action_button("run_analysis", "Run Analysis", disabled=True)
	)

