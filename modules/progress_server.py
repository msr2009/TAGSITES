from shiny import reactive, render, ui
import json, subprocess

from utils.helpers import update_shared_dict

def progress_server(input, output, session, shared_json):

	# Display JSON content if it exists
	@output
	@render.ui
	def show_or_upload_json():
		ui_elements = []

		# make a persistent file upload button
		ui_elements.append(ui.input_file("upload_json", "Upload/Reload JSON File"))

		# Check if session.user_data has the "setup_json_out" saved
		if shared_json.get():
			# Display the JSON content
			filename = shared_json.get()
			try:
				with open(filename, "r") as f:
					json_content = json.load(f)
				# Pretty print JSON in the UI
				ui_elements.append(ui.pre(json.dumps(json_content, indent=4)))
			except FileNotFoundError:
				ui_elements.append(ui.div("Error: Saved JSON file not found."))
		else:
			# Prompt the user to upload a new JSON file
			ui_elements.append(ui.div("No JSON file loaded. Please upload a JSON file."))

		#return ui elements
		return ui.div(*ui_elements)

	# Handle the file upload
	@reactive.Effect
	@reactive.event(input.upload_json)
	def handle_file_upload():
		file_info = input.upload_json()
		if file_info:
			# Save the uploaded file path to session.user_data["setup_json_out"]
			uploaded_file_path = file_info[0]["datapath"]  # Get the file path of the uploaded file
#			shared_values.set("setup_json_out") = uploaded_file_path
			shared_json.set(uploaded_file_path)

			# Optionally read and print the file to the console (for testing)
			with open(uploaded_file_path, "r") as f:
				json_data = json.load(f)
				print("Uploaded JSON content:", json.dumps(json_data, indent=4))

	# upon loading a json file, you can run the analysis
	@reactive.effect
	def update_run_analysis_button_state():
		if shared_json.get():
			ui.update_action_button('run_analysis', disabled=False)
		else:
			ui.update_action_button('run_analysis', disabled=True)

	@reactive.effect
	@reactive.event(input.run_analysis)
	def run_analysis():
		#load json file
		with open(shared_json.get()) as f:
			INPUT_JSON = json.load(f)
		subprocess.call('python {}run_tag_sites_from_json.py -i {}'.format(INPUT_JSON["global"]["scripts-folder"], shared_json.get()), shell=True)
