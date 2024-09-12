"""
run_tag_sites_from_json.py

takes json file with analyses to be run 
(e.g., from shiny app)

submits asynchronous jobs using subprocess

Matt Rich, 09/24
"""

import subprocess, os, json

def main(json_in):
#	print(json_in) 
	
	task_scripts = json_in["scripts"]
	
	commands = []

	#loop through tasks and make calls
	for task in json_in:
		#skip scripts
		if task == "scripts": continue #this is just a mapping dict
		#build command for analysis calls
		interpreter = "python"
		scripts_folder = "scripts/"
		script = "{} {}".format(interpreter, scripts_folder) + task_scripts[json_in[task]["analysis"]]
		task_args = json_in[task]["args"]
		arguments = [" --{} {}".format(arg, task_args[arg]) for arg in task_args if task_args[arg] != ""]

		commands.append(" ".join([script] + arguments))
	print("\n".join(commands))
	print("#############\n############\n#############")

	#run the commands
	procs = [ subprocess.Popen(i, shell=True) for i in commands ]
	for p in procs:
		p.wait()

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-i', '--input', '--json', action = 'store', type = str,
					dest = 'JSON_IN', help = "HELP")
	args = parser.parse_args()
	
	try:
		main(json.load(open(args.JSON_IN, "r")))	
	except json.decoder.JSONDecodeError:
		print("Incorrectly formatted JSON input!")


