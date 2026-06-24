from shiny import ui, reactive, render, Inputs, Outputs, Session

#from config import INPUT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS, EXCLUDE_ARGS
from config import DEFAULT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS, EXCLUDE_ARGS

from config import UNIPROT_SPECIES, DEFAULT_SPECIES

from scripts.site_selection_util import get_sequence, save_fasta 
from utils.helpers import load_taxonomic_mapping, update_shared_dict

import json, os, shutil, pandas, copy, pickle, time
from pathlib import Path

#def setup_server(input: Inputs, output: Outputs, session: Session, shared_json):
def setup_server(input, output, session, shared_json):
	
	selected_tasks = reactive.Value([])  # To store added tasks
	task_values = reactive.Value({})  # to maintain inputted values
#	taxonomic_mapping = reactive.Value(None) # to hold taxonomy
#	use_search = reactive.Value(False) # if species search is required
#	search_mode = reactive.Value(False) # user is using the search bar 


	#doing this to reset INPUT_JSON upon refreshing app
	INPUT_JSON = json.load(open(DEFAULT_JSON,"r"))

#	def load_taxonomy():
#		taxonomic_mapping.set(load_taxonomic_mapping(UNIPROT_SPECIES))
#		taxonomic_mapping.set(pickle.load(open("uniprot_species.pkl", "rb")))
#		populate_species_selectize()

	# Dynamically populate the selectize input with species name
#	@reactive.effect
#	def populate_species_selectize():
#		# Send the list of species as choices to the selectize box
#		mapping = taxonomic_mapping()
#		if mapping is None:
#			return
#		
#		species_list = list(taxonomic_mapping().keys())
#		ui.update_selectize("species_search", choices=species_list)
#		ui.update_selectize("species_search", selected="")

#	session.on_flush(load_taxonomy, once=True)

#	with open("uniprot_species.pkl", "rb") as f:
#		full_taxonomic_mapping = pickle.load(f)
#	taxonomic_mapping = reactive.Value(full_taxonomic_mapping)

	def populate_default_params():
		# load names of params json files in ./params/
		params_dir = Path("./params/")
		params_files = [j.name.removesuffix(".json") for j in list(params_dir.glob("*.json"))]
		ui.update_selectize("load_default_tasks", choices=params_files) 
		ui.update_selectize("load_default_tasks", selected="")

	def populate_tables_list():
		# load names of files in ./tables/
		tables_dir = Path("./tables/")
		tables_files = list(tables_dir.glob(""))
		return tables_files

	#load uniprot species file and update dropdown
#	taxonomic_mapping = load_taxonomic_mapping(UNIPROT_SPECIES)
#	session.on_flush(populate_species_selectize, once=True)
	session.on_flush(populate_default_params, once=True)

### this is all maybe getting cut

#	@reactive.effect
#	@reactive.event(input.species_search)
#	def handle_primary_species():
#		selected_species = input.species_search()

#		if selected_species == "Other (search...)":
#			#show search box
#			use_search.set(True)
#			INPUT_JSON["global"]["species"] = None
#			INPUT_JSON["global"]["species_taxid"] = None
		
#		elif use_search():
#			use_search.set(False)
#			INPUT_JSON["global"]["species"] = selected_species
#			INPUT_JSON["global"]["species_taxid"] = DEFAULT_SPECIES.get(selected_species, None)
		
#		else:
#			#one of the defaults was selected
#			INPUT_JSON["global"]["species"] = selected_species
#			INPUT_JSON["global"]["species_taxid"] = DEFAULT_SPECIES.get(selected_species, None)
#		print("selected species taxid: ", INPUT_JSON["global"]["species_taxid"])

#	@output
#	@render.ui
#	def species_search_ui():
#		print("use_search: ", use_search())
#		if not use_search(): #we aren't searching
#			return None
#		else:
#			#otherwise we need the search textbox
#			return ui.input_text( 
#				"species_query",
#				"Search species",
#				placeholder = "Start typing..."
#				)

#	@reactive.effect
#	@reactive.event(input.species_search)
#	def handle_search_selection():
#		#only update if dynamic search is active
#		if use_search() and input.species_search() and input.species_search() != "Other (search...)":
#			selected_species = input.species_search()
##			INPUT_JSON["global"]["species"] = selected_species
#			INPUT_JSON["global"]["species_taxid"] = taxonomic_mapping()[selected_species]


#	@reactive.effect
#	@reactive.event(input.species_query)
#	def filter_species():	
#		query = input.species_query()
#		if not query or len(query) < 3:
#			return

#		mapping = taxonomic_mapping()
#		matches = [name for name in mapping.keys() if query.lower() in name.lower()][:50]

#		ui.update_selectize(
#			"species_search",
#			choices=matches,
#			selected=None
#		)

#	@reactive.effect
#	@reactive.event(input.species_search)
#	def selected_species():	
#		selected_species = input.species_search()
#		if selected_species != "":
#			INPUT_JSON["global"]["species"] = selected_species
#			INPUT_JSON["global"]["species_taxid"] = taxonomic_mapping[selected_species]
#			print(taxonomic_mapping[selected_species])

	#for updating working directory once run_name is set
	@reactive.effect
	@reactive.event(input.run_name)
	def update_working_dir():
		out_dir = os.getcwd()
		if INPUT_JSON["global"]["working_dir"] != "":
			out_dir = INPUT_JSON["global"]["working_dir"]
		if input.run_name():
			ui.update_text("working_dir", value=f"{out_dir}/{input.run_name()}/")

	#can only put defaults button when you've chosen one set of defaults
	@reactive.effect
	def update_load_defaults_button_state():
		if input.load_default_tasks():
			ui.update_action_button("load_defaults_button", disabled=False)
		else:
			ui.update_action_button("load_defaults_button", disabled=True)

	# checks that a name was added for a task before allowing button push
	@reactive.effect
	@reactive.event(input.task_desc_name)
	def update_task_name_button_state():
		if input.task_desc_name():
			ui.update_action_button("add_task", disabled=False)
		else:
			ui.update_action_button("add_task", disabled=True)
	
	@reactive.Calc
	def requirements_filled():
		print("checking requirements")
		return bool(input.email()) and bool(input.input_file) and bool(input.working_dir()) and bool(len(selected_tasks())>0)

	@reactive.effect
	@reactive.event(requirements_filled)
	def update_save_analysis_button_state():
#		if input.email() and input.input_file() and input.run_name() and input.working_dir() and len(selected_tasks()) > 0:
		if requirements_filled():
			ui.update_action_button('save_analysis', disabled=False)
		else:
			ui.update_action_button('save_analysis', disabled=True)
	
	@render.text
	def tip_save_analysis():
		if not requirements_filled():
			return "Requires inputs, analysis name, and at least one task."
		else:
			if input.save_analysis() == 0:
				return "Click to save analysis"
			else:
				return f"Saved to {input.working_dir()}{input.run_name()}.json"

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

#				print("{}:{} current: {}, selected: {}".format(task['id'], p, current_value, selected_value))

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


	#push button and add default tasks to list
	@reactive.effect
	@reactive.event(input.load_defaults_button)
	def add_tasks_from_default():
		default_params_file = input.load_default_tasks()
		if default_params_file != "":
			default_task_json = json.load(open("./params/{}.json".format(input.load_default_tasks()),"r"))
			for task_id in default_task_json:
				print(task_id)
				
				selected_tasks.set(selected_tasks() + [{"id": task_id, 
													    "name": default_task_json[task_id]["analysis"], 
														"params": default_task_json[task_id]["args"]}])
#			selected_tasks.set(selected_tasks() + [default_task_json])
			save_task_values()
		print(selected_tasks())

	#update save_default_tasks_button
	@reactive.effect
	def update_save_default_tasks_button():
		if input.new_default_name() != "" and len(selected_tasks()) > 0:
			ui.update_action_button("save_default_tasks_button", disabled=False)	
		else:
			ui.update_action_button("save_default_tasks_button", disabled=True)	

	#save current tasks as a default to ./params/
	@reactive.effect
	@reactive.event(input.save_default_tasks_button)
	def make_new_default():
		save_task_values()
		new_default_filename = f"./params/{input.new_default_name()}.json"
		new_default_json_file = open(new_default_filename, "w")
		def_json = {}
		for task in selected_tasks():
			def_json = {**def_json, **write_task_json(task)}
		json.dump(def_json, new_default_json_file, indent=4)
		new_default_json_file.close()

	def save_task_values():
		for task in selected_tasks():
			task_params = {}
			for param in task['params'].keys():
				input_id = f"{task['id']}_{param}"
				if input_id in input:
					task_params[param] = input[input_id]()				
			task_values()[task['id']] = task_params
		print(task_values())

	#build list of tasks (this gets updated dynamically)
	@reactive.effect
	@reactive.event(input.add_task)
	def add_task(custom_params=None):
		task_name = input.task_selector()

		if task_name:
			#before making a new task, save the existing data in the old tasks
			save_task_values()

			#then, make a new task
			task_id = f"{input.task_desc_name()}_{task_name}"
			
			#grab appropriate arguments for each task
			default_task_params = INPUT_JSON["analyses"][task_name]["args"]
			task_params = {p: default_task_params[p] for p in default_task_params if p not in EXCLUDE_ARGS}
			
			#grab tooltips for each task
			default_task_tips = INPUT_JSON["analyses"][task_name]["tooltips"]
			task_tooltips = {p: default_task_tips[p] for p in default_task_tips if p not in EXCLUDE_ARGS}

			# Append the task with default parameters
			selected_tasks.set(selected_tasks() + [{"id": task_id, "name": task_name, "params": task_params, "tooltips": task_tooltips}])
#			print(selected_tasks())
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
	def write_task_json(task, global_params=None):
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
#			print(input_id, input_set)
			param_value = ""
			if not input_set:
				param_value = TASK_PARAMETERS[task["name"]]["args"][p]
			else:
				param_value = input[input_id]()  # Add parentheses to call the reactive value

#			print(p,param_value)
			
			task_json[task['id']]['args'][p] = param_value
#			print(task_json)a
		#and add all the globals -- we can sort out the requirements later
		if global_params != None:
			task_json[task['id']]['args'] = {**task_json[task['id']]['args'], **global_params}
		#and update output name
		task_json[task['id']]['args']['output'] = f"{input.working_dir()}{input.run_name()}.{task['id']}{TASK_PARAMETERS[task['name']]['args']['output']}"
		return task_json

	#save analyses and parameters to json
	@reactive.effect
	@reactive.event(input.save_analysis)
	def save_analysis():
	
#		print("SAVING")

		#first, make working directory if it doesn't exist
		try:
			Path(input.working_dir()).mkdir()
		except FileExistsError:
			pass

		#update global parameters in INPUT_JSON
		INPUT_JSON["global"]["email"] = input.email()
		INPUT_JSON["global"]["run_name"] = input.run_name()
		INPUT_JSON["global"]["working_dir"] = input.working_dir()
		
#		print("BEFORE INPUTS", INPUT_JSON)

		#save uploaded input file with better name
		tmp_input = Path(input.input_file()[0]["datapath"])
		new_input_filename = input.working_dir() + input.run_name() 
		if tmp_input.suffix == ".pdb":
				new_input_filename += ".pdb"
		else:
			new_input_filename += ".fa" #change .fasta to .fa

		
		INPUT_JSON["global"]["input_file"] = new_input_filename
		#copy into working directory
		shutil.copy(tmp_input, new_input_filename)

		#save uploaded genomic fasta with better name
		tmp_genomic = Path(input.input_genomic()[0]["datapath"])
		new_genomic_filename = input.working_dir() + input.run_name() + ".genomic" + tmp_genomic.suffix
		INPUT_JSON["global"]["genomic_file"] = new_genomic_filename
		#copy into working directory
		shutil.copy(tmp_genomic, new_genomic_filename)

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
#				print("PARAMS", task["params"])
#				print("INPUT_GLOBAL", INPUT_JSON["global"])
				#if there isn't a pdb file in the input
				if INPUT_JSON["global"]["pdb"] != "":	#############I DONT THINK THIS WORKS
					task["params"]["existing_AF2"] = 0
				else:
					task["params"]["existing_AF2"] = 1
#				print("EXISTINGAF2", task["params"])
	
			##############
			#if we have a genomic sequence, then also add a genewise task
			###############

#			print(f"Task ID: {task['id']}, Task Name: {task['name']}")  # Debug: Print task details
			out_json = {**out_json, **write_task_json(task, INPUT_JSON["global"])}

			# Collect the task parameters (assuming they exist)
#			task_params = {param: input[f"{task['id']}_{param}"]() for param in task['params'].keys()}
#			task_list.append({"id": task['id'], "name": task['name'], "params": task_params})

		#write task jsons to new json file in working directory
		json.dump(out_json, out_json_file, indent=4)
		out_json_file.close()
		shared_json.set(out_json_name)

	################
	####TOOLTIPS####
	################
	
	#not really a tooltip, but some text next to the save_analysis button
	#change to path to json once there's a task and required names/dirs







