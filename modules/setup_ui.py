from shiny import ui

from config import INPUT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS

##SETUP TAB UI
def setup_ui():

	#check for output directory in json, 
	#if there isn't one, then use os.get_cwd()

	return ui.page_fluid(
		
		# Global parameters input
		ui.input_text("email", "Email Address (required)",
					   value=INPUT_JSON["global"]["email"]),
		
		ui.hr(),

		# Input Sequence
		ui.h3("Sequence Input"),
		ui.row(
			ui.column(2, ui.input_file("input_file", "Choose Input File (.fa, .fasta, or .pdb)", accept=[".fa", ".fasta",".pdb"])),
#			ui.column(2, ui.input_radio_buttons("input_type", "Input Type",
#					{"fasta": "FASTA", "pdb": "PDB"}, selected=None))
		),
		ui.hr(),

		# Task selector and task management
		ui.h3("Analysis Steps"),
		ui.row(
			ui.column(1, ui.input_text("run_name", "Analysis Name")),
			ui.column(3, ui.input_text("working_dir", "Output Directory", width="100%")),
		),
		ui.hr(),

		ui.row(
			ui.column(1, ui.input_select("task_selector", "Select Analysis", AVAILABLE_TASKS)),
			ui.column(1, ui.input_text("task_desc_name", "Name (required)")),
		),
		ui.input_action_button("add_task", "Add Task", disabled=True),
		ui.hr(),

		# Container for dynamically added tasks
		ui.output_ui("task_params_container"),	

		ui.hr(),
		ui.input_action_button("save_analysis", "Save Analysis", disabled=True)
	)

