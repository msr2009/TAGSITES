from shiny import ui, reactive, render, Inputs, Outputs, Session

from config import INPUT_JSON, TASK_PARAMETERS, AVAILABLE_TASKS, EXCLUDE_ARGS

import json, os, shutil
from pathlib import Path

def setup_server(input: Inputs, output: Outputs, session: Session):
		
	#for updating working directory once run_name is set
	@reactive.effect
	@reactive.event(input.run_name)
	def update_working_dir():
	
		out_dir = os.getcwd()
		if INPUT_JSON["global"]["directory"] != "":
			out_dir = INPUT_JSON["global"]["directory"]

		if input.run_name():
			ui.update_text("working_dir", value=f"{out_dir}/{input.run_name()}")


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
				param_list.append(ui.input_select(f"{task['id']}_{p}", 
						label=p,
						selected=task_values().get(task['id'], {}).get(p, task["params"][p][0]),
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
			#task_params = TASK_PARAMETERS["analyses"][task["name"]]["args"]

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
		for p in task['params'].keys():
			input_id = f"{task['id']}_{p}"
			param_value = input[input_id]()  # Add parentheses to call the reactive value
			task_json[task['id']]['args'][p] = param_value
		#and add all the globals -- we can sort out the requirements later
		task_json[task['id']]['args'] = {**task_json[task['id']]['args'], **global_params}
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

		#global parameters
		global_parameters = {
						"email": input.email(),
						"working_dir": input.working_dir(),
						"run_name": input.run_name(),
						}
		
		#save uploaded input file with better name
		tmp_input = Path(input.input_file()[0]["datapath"])
		new_input_filename = input.run_name() + tmp_input.suffix
		global_parameters["input_file"] = new_input_filename

		#now, copy the temp input file into the working directory
		shutil.copy(tmp_input, input.working_dir()+"/"+new_input_filename)
		
		#open a file to dump json into
		out_json_file = open(input.working_dir() + "/" + input.run_name() + ".json", "w")
		out_json = {}
		out_json["scripts"] = INPUT_JSON["scripts"]

		#print all tasks to json file
		task_list = []

		for task in selected_tasks():
			print(f"Task ID: {task['id']}, Task Name: {task['name']}")  # Debug: Print task details
			out_json = {**out_json, **write_task_json(task, global_parameters)}

			# Collect the task parameters (assuming they exist)
			task_params = {param: input[f"{task['id']}_{param}"]() for param in task['params'].keys()}
			task_list.append({"id": task['id'], "name": task['name'], "params": task_params})

		#write task jsons to new json file in working directory
		print(out_json)
		print(json.dumps(out_json, indent=4), file=out_json_file)
		out_json_file.close()

