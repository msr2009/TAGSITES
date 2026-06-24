from shiny import ui

from config import INPUT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS, UNIPROT_SPECIES, DEFAULT_SPECIES

from utils.helpers import load_taxonomic_mapping

def setup_ui():

	return ui.page_fluid(
		
		# Global parameters input
		ui.input_text("email", "Email Address (required)",
					   value=INPUT_JSON["global"]["email"]),
		
		ui.hr(),

		# Input Sequence
		ui.h3("Sequence Input"),
		ui.row(
			ui.column(2, ui.input_file("input_file", "Choose Input File (.fa, .fasta, or .pdb)", accept=[".fa", ".fasta",".pdb"])),
			ui.column(2, ui.input_file("input_genomic", "FASTA file of genomic region", accept=[".fasta", ".fa"])),
		),
#		ui.input_selectize("species_search", "Species", choices=[], multiple=False, selected=None,
#						options={'create': False}, width="700px"),

#		ui.input_selectize("species_search", "Species",
#				choices=list(DEFAULT_SPECIES.keys()), multiple=False,
#				selected=None, options={'create': False}, width="700px"),
#		ui.output_ui("species_search_ui"),

		ui.hr(),

		# Task selector and task management
		ui.h3("Setup Analyses"),
		ui.row(
			ui.column(1, ui.input_text("run_name", "Analysis Name")),
			ui.column(3, ui.input_text("working_dir", "Output Directory", width="100%")),
		),
		
		ui.row(
		ui.column(2, ui.input_selectize("load_default_tasks", "Load Default Tasks?", choices=[], multiple=False, selected=None,
				options={'create': False}, width="400px")),
		ui.column(1, ui.input_action_button("load_defaults_button", "Load", disabled=True)),
		),
		ui.hr(),

		ui.row(
			ui.column(1, ui.input_select("task_selector", "Select analysis task", AVAILABLE_TASKS)),
			ui.column(1, ui.input_text("task_desc_name", "Task name (required)")),
		),
		ui.input_action_button("add_task", "Add Task", disabled=True),
		ui.hr(),

		# Container for dynamically added tasks
		ui.output_ui("task_params_container"),	

		ui.hr(),
		
		ui.row(
			ui.column(1, ui.input_action_button("save_analysis", "Save Analysis", disabled=True)),
			ui.column(4, ui.output_text_verbatim("tip_save_analysis")),
		),

		ui.hr(),
		
		ui.row(
			ui.column(1, ui.input_action_button("save_default_tasks_button", "Save as Default?", disabled=True)),
			ui.column(2, ui.input_text("new_default_name", "", ""),
		)
	)
)
