from shiny import ui, reactive, render, Inputs, Outputs, Session

from config import INPUT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS, EXCLUDE_ARGS
from config import UNIPROT_SPECIES

from scripts.site_selection_util import get_sequence, save_fasta 
from utils.helpers import load_taxonomic_mapping, update_shared_dict

import json, os, shutil, pandas
from pathlib import Path


def setup_server(input: Inputs, output: Outputs, session: Session, shared_json):
	

		# Dynamically populate the selectize input with species name
#	@reactive.event(input.species_search, ignore_none=False)
	def populate_species_selectize():
		# Send the list of species as choices to the selectize box
		species_list = list(taxonomic_mapping.keys())
		ui.update_selectize("species_search", choices=species_list)
		ui.update_selectize("species_search", selected="")

	#load uniprot species file and update dropdown
	taxonomic_mapping = load_taxonomic_mapping(UNIPROT_SPECIES)
	session.on_flush(populate_species_selectize, once=True)

	@reactive.effect
	@reactive.event(input.species_search)
	def selected_species():	
		selected_species = input.species_search()
		if selected_species != "":
			INPUT_JSON["global"]["species"] = selected_species
			INPUT_JSON["global"]["taxid"] = taxonomic_mapping[selected_species]
			print(taxonomic_mapping[selected_species])

	#for updating working directory once run_name is set
	@reactive.effect
	@reactive.event(input.run_name)
	def update_working_dir():
		
		out_dir = os.getcwd()	
		if INPUT_JSON["global"]["working_dir"] != "":
			out_dir = INPUT_JSON["global"]["working_dir"]

#		out_dir = INPUT_JSON["global"]["working_dir"]

		if input.run_name():
			ui.update_text("working_dir", value=f"{out_dir}/{input.run_name()}/")
			INPUT_JSON["global"]["working_dir"] = input.working_dir()

	# checks that a name was added for a task before allowing button push
	@reactive.effect
	@reactive.event(input.task_desc_name)
	def update_task_name_button_state():
		if input.task_desc_name():
			ui.update_action_button("add_task", disabled=False)
		else:
			ui.update_action_button("add_task", disabled=True)
	
	@reactive.effect
	def update_save_analysis_button_state():
		if input.email() and input.input_file() and input.run_name() and input.working_dir():
			ui.update_action_button('save_analysis', disabled=False)
		else:
			ui.update_action_button('save_analysis', disabled=True)
		

	selected_tasks = reactive.Value([])  # To store added tasks
	task_values = reactive.Value({})  # to maintain inputted values

	# builds the row for parameters for each task
	def build_task_params_input(task):
		param_list = []
		for p in task["params"]:
			if p in EXCLUDE_ARGS: continue

			#if parameters are a list, then make it a dropdown
			if isinstance(task["params"][p], list) == True:
				#retrieve the current value for a dropdown
				current_value = task_values().get(task['id'], {}).get(p)
				#use the currect value if it exists, otherwise use default
				selected_value = current_value if current_value is not None else task["params"][p][0]

				print("{}:{} current: {}, selected: {}".format(task['id'], p, current_value, selected_value))

				param_list.append(ui.input_select(f"{task['id']}_{p}", 
						label=p,
						selected=selected_value,
						choices=task["params"][p], size=1))

			#otherwise it's free text input
			else:
				param_list.append(ui.input_text(f"{task['id']}_{p}",
						value=task_values().get(task['id'], {}).get(p, task["params"][p]),			
						label=p))

		#add remove task button
#		param_list.append(ui.input_action_button("remove_task", "Remove"))
		return param_list

	#build list of tasks (this gets updated dynamically)
	@reactive.effect
	@reactive.event(input.add_task)
	def add_task():

		task_name = input.task_selector()
		if task_name:

			#before making a new task, save the existing data in the old tasks
			for task in selected_tasks():
				task_params = {}
				for param in task['params'].keys():
					input_id = f"{task['id']}_{param}"
					if input_id in input:
						task_params[param] = input[input_id]()
				task_values()[task['id']] = task_params

			task_id = f"{input.task_desc_name()}_{task_name}"
			
			#grab appropriate arguments for each task
			default_task_params = INPUT_JSON["analyses"][task_name]["args"]
			task_params = {p: default_task_params[p] for p in default_task_params if p not in EXCLUDE_ARGS}
			
			#grab tooltips for each task
			default_task_tips = INPUT_JSON["analyses"][task_name]["tooltips"]
			task_tooltips = {p: default_task_tips[p] for p in default_task_tips if p not in EXCLUDE_ARGS}

			# Append the task with default parameters
			selected_tasks.set(selected_tasks() + [{"id": task_id, "name": task_name, "params": task_params, "tooltips": task_tooltips}])
			print(selected_tasks())
			ui.update_text("task_desc_name", value="")

	@reactive.effect
	@reactive.event(input.remove_task)
	def remove_task():
		return

	#render task parameters
	@output
	@render.ui
	def task_params_container():
		task_ui_list = []
		for task in selected_tasks():
#			print(task)
			task_id = task['id']
			task_params = task['params']

			task_ui = ui.div(
				ui.h3(f"{'_'.join(task['id'].split('_')[0:-1])} ({task['name']}) parameters"),
				ui.row(
					ui.row(build_task_params_input(task)),
					ui.input_action_button(f"{task['id']}_remove", "Remove")
				),
				ui.hr()
			)

			task_ui_list.append(task_ui)

		return ui.div(*task_ui_list)
	

	#write json entry (including global args)
	def write_task_json(task, global_params):
		task_json = {
			task['id']: {
				"analysis": task['name'],
				"args": {}
				}
			}
		#add all the inputted parameters (and defaults if not changed)
		for p in TASK_PARAMETERS[task["name"]]["args"]:
			input_id = f"{task['id']}_{p}"
			input_set = input[input_id].is_set()
			param_value = ""
			if not input_set:
				param_value = TASK_PARAMETERS[task["name"]]["args"][p]
			else:
				param_value = input[input_id]()  # Add parentheses to call the reactive value

			print(p,param_value)
			
			task_json[task['id']]['args'][p] = param_value
		#and add all the globals -- we can sort out the requirements later
		task_json[task['id']]['args'] = {**task_json[task['id']]['args'], **global_params}
		#and update output name
		task_json[task['id']]['args']['output'] = f"{input.working_dir()}{input.run_name()}.{task['id']}{TASK_PARAMETERS[task['name']]['args']['output']}"
		return task_json

	#save analyses and parameters to json
	@reactive.effect
	@reactive.event(input.save_analysis)
	def save_analysis():

		#first, make working directory if it doesn't exist
		try:
			Path(input.working_dir()).mkdir()
		except FileExistsError:
			pass

		#update global parameters in INPUT_JSON
		INPUT_JSON["global"]["email"] = input.email()
		INPUT_JSON["global"]["run_name"] = input.run_name()
		
		#save uploaded input file with better name
		tmp_input = Path(input.input_file()[0]["datapath"])
		new_input_filename = input.working_dir() + input.run_name() + tmp_input.suffix
		INPUT_JSON["global"]["input_file"] = new_input_filename

		#now, copy the temp input file into the working directory
		shutil.copy(tmp_input, new_input_filename)

		#if input is pdb, then we need to extract the fasta
		if tmp_input.suffix == ".pdb":
			INPUT_JSON["global"]["pdb"] = INPUT_JSON["global"]["input_file"]
			#extract and save FASTA
			new_fasta_name = INPUT_JSON["global"]["pdb"].replace(".pdb", ".fa")
			save_fasta(Path(INPUT_JSON["global"]["pdb"]).stem, 
					   get_sequence(INPUT_JSON["global"]["pdb"]),
					   new_fasta_name)
			#and update input_file in json
			INPUT_JSON["global"]["input_file"] = new_fasta_name
	
		#open a file to dump json into
		out_json_name = input.working_dir() + "/" + input.run_name() + ".json"
		out_json_file = open(out_json_name, "w")
		out_json = {}
		out_json["scripts"] = INPUT_JSON["scripts"]
		out_json["global"] = INPUT_JSON["global"]

		#print all tasks to json file
		task_list = []

		for task in selected_tasks():
			#we need to check plddt and input type
			if task["name"] == "plddt":
				#if there isn't a pdb file in the input
				if INPUT_JSON["global"]["pdb"] == "" or task["params"]["pdb"] == "":
					task["params"]["existing_AF2"] = 1
	
#			print(f"Task ID: {task['id']}, Task Name: {task['name']}")  # Debug: Print task details
			out_json = {**out_json, **write_task_json(task, INPUT_JSON["global"])}

			# Collect the task parameters (assuming they exist)
#			task_params = {param: input[f"{task['id']}_{param}"]() for param in task['params'].keys()}
#			task_list.append({"id": task['id'], "name": task['name'], "params": task_params})

		#write task jsons to new json file in working directory
		print(json.dumps(out_json, indent=4), file=out_json_file)
		out_json_file.close()
		shared_json.set(out_json_name)
