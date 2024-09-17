"""
run_tag_sites_from_json.py

takes json file with analyses to be run 
(e.g., from shiny app)

submits asynchronous jobs using subprocess

Matt Rich, 09/24
"""

import subprocess, os, json

#function to determine whether we need to search AFDB
def searchAFDB_required(j):
	for key in j:
		if j[key]["analysis"] == "plddt": #if pLDDT in analyses
			if j[key]["args"]["existing_AF2"]: #if reqd
				#return args
				return build_task_args_string(j[key]["args"], ["pdb", "existing_AF2", "output"])
	return None

#function to create argument string from json entry
def build_task_args_string(targs, EXCLUDE=[]):
	return "".join([" --{} {}".format(arg, targs[arg]) for arg in targs if targs[arg] != "" and arg not in EXCLUDE])

def main(json_in):
	
	task_scripts = json_in.pop("scripts", None)
	#json_in is now just tasks
	
	#print(json_in)
	scripts_folder = "scripts/"
	interpreter = "python"
	
	#check if afdb search required
	search_afdb = searchAFDB_required(json_in)
	print(search_afdb)

	if search_afdb != None:
		#do AFDB search
		af2_search_call = "python {}existing_AF_model.py {}".format(scripts_folder, search_afdb)
		print("SEARCHING AFDB: " + af2_search_call)
		subprocess.call(af2_search_call, shell=True)
		#if successful, then extract fasta from pdb
		#then update json to save old inputs (user_[input]) 
		#correct inputs and outputs for those in found pdb
		#set existing_AF2=0, 

	return
	commands = []

	#loop through tasks and make calls
	for task in json_in:
		#skip scripts
		if task == "scripts": continue #this is just a mapping dict
		#build command for analysis calls
		script = "{} {}".format(interpreter, scripts_folder) + task_scripts[json_in[task]["analysis"]]
#		task_args = json_in[task]["args"]
#		arguments = [" --{} {}".format(arg, task_args[arg]) for arg in task_args if task_args[arg] != ""]
		arguments = build_task_args_string(json_in[task]["args"])	

		commands.append(script + arguments)


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


