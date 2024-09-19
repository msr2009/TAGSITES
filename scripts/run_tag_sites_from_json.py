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

def main(json_input_file):
	
	#load json
	try:
		json_in = json.load(open(json_input_file, "r"))
	except json.decoder.JSONDecodeError:
		raise IOError("Incorrectly formatted JSON file!")

	task_scripts = json_in.pop("scripts", None)
	global_args = json_in.pop("global", None)
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
		af_proc = subprocess.call(af2_search_call, shell=True)
#		af_exit = af_proc.returncode

		#if existing_AF_model.py finds new pdb, then it doesn't return 1
		#it will write new fasta as [run_name].AF.fa
		if af_proc != 1:
			#then we need to rename original fasta and update json
			os.rename("{}/{}.fa".format(global_args["working_dir"], global_args["run_name"]),
					  "{}/user_{}.fa".format(global_args["working_dir"], global_args["run_name"]))
				
		#then update json to save old inputs (user_[input]) 
		#correct inputs and outputs for those in found pdb
		#set existing_AF2=0, 
		
		for task in json_in: 
			#for each task, move the current input_file (.fa) to user_input_file
			#make an old_args to store these
			json_in[task]["old_args"] = {}
			json_in[task]["old_args"]["user_input_file"] = "user" + json_in[task]["args"]["input_file"]
			#and rename the input_file to the new .AF.fa
			json_in[task]["args"]["input_file"] = json_in[task]["args"]["input_file"].replace(".fa", ".AF.fa")
			#if plddt task, then set existing_AF2=0
			if json_in[task]["analysis"] == "plddt":
				json_in[task]["args"]["existing_AF2"] = 0
				json_in[task]["args"]["pdb"] = "{}/{}.AF.pdb".format(global_args["working_dir"], global_args["run_name"])
		
		#also append global and scripts for new json
		
	
		#write a new json file for posterity
		new_json_out = {**json_in, **global_args, **task_scripts}
		new_json_out_file = open(json_input_file.replace(".json", ".AF.json"), "w")
		print(json.dumps(json_in, indent=4), file=new_json_out_file) 
		new_json_out_file.close()
	
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
					dest = 'JSON_INPUT', help = "HELP")
	args = parser.parse_args()
	
	main(args.JSON_INPUT)	


