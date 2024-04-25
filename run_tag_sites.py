"""
call_interpro.py

A script to perform an InterPro search.

Uses EBI's API (iprscan5). 

Parameters:
	- fasta_in (str): name of fasta file containing seq

Returns: 
	- interpro_out (str): BED format str with for interpro results

Matt Rich, 4/2024
"""

import subprocess, sys, os
from shutil import copy2
from datetime import datetime
from site_selection_util import *

def main(infile, email, workingdir, name, inputtype, clients_folder, run_af2, sites_file):
	
	#make workingdir
	if os.path.exists(workingdir) == False:
		os.mkdir(workingdir)

	#list to store commands
	commands = []

	#new file names
	fastafile = "{}/{}.fa".format(workingdir, name)
	pdbfile = "{}/{}.pdb".format(workingdir, name)
	
	#determine input type
	if inputtype == None:
		inputtype = check_input_type(infile)
	
	if inputtype == "fasta":
		seq = get_sequence(infile)
		save_fasta(name, seq, fastafile)
	elif inputtype == "err":
		sys.exit("can't figure out input type. use --input-type argument.")
	
	#if input is fasta, then we need to find/make a pdb
	#we want to do this first, because there might not an existing one, 
	#and should fail if we can't run AF2 locally.
	if inputtype == "fasta" or "uniprot":
		
		#check if there's one we can download
		existing_call = "python existing_AF_model.py \
							-f {} --dir {} \
							--email {} --name {} \
							--clients-folder {}".format(infile, workingdir,
														email, name, clients_folder)
		try:
			subprocess.run(existing_call, shell=True, check=True)
		except subprocess.CalledProcessError:
			# no existing model, so we can make one if possible
			if run_af2:
				local_af2_call = "python local_AF2_prediction.py -f {} --dir {}".format(fastafile, 
																						workingdir)
				commands.append(local_af2_call) 
			else:
				sys.exit()("no existing alphafold prediction for input and \
								I can't do it myself. please make one!")
	
	if inputtype == "pdb":
		copy2(infile, pdbfile)

	#after dealing with making or getting a pdb, we can compile our list of jobs			
	#interpro
	interpro_call = "python call_interpro.py -f {} --dir {} --email {} --clients-folder {}".format(fastafile, 
																									workingdir, 
																									email, 
																									clients_folder)
	commands.append(interpro_call)
	
	#regex_sites (functional sites)
	sites_call = "python regex_sites.py -f {} --sites {} -o {}/{}.sites.txt".format(fastafile, 
																					sites_file, 
																					workingdir, 
																					name)
	commands.append(sites_call)
	
	#pLDDT
	plddt_call = "python pLDDT_minima.py -p {} --sites_out {} --smoothed_plddt_out {}".format(pdbfile, 
																"{}/{}.minima.txt".format(workingdir,name), 
																"{}/{}.smoothed.txt".format(workingdir,name)) 
																						
	commands.append(plddt_call)

	print("\n".join(commands))

	#run the commands
	procs = [ subprocess.Popen(i, shell=True) for i in commands ]
	for p in procs:
		p.wait()
	

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-i', '--input', action='store', type=str, dest='INFILE', 
		help = "input file. can be fasta or pdb (determined automatically).", required=True)
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('-n', '--name', action='store', type=str, dest='NAME', 
		help = "job name for files, fasta header, etc...", required=True) 

	parser.add_argument('--input-type', action='store', type=str, dest='INPUTTYPE', 
		help = "input type (fasta, pdb).", default=None)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER',
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")

	parser.add_argument('--run-af2', action='store_true', dest='RUN_AF2', 
		help = "run alphafold locally. requires localcolabfold (default=False).", default=False)
	parser.add_argument('--sites_file', action='store', type=str, dest='SITES_FILE', 
		help = "tab-delimited file containing regular expressions for modification sites (default=modifications.txt)", 
		default="modifications.txt")
	
	args = parser.parse_args()
	
	main(args.INFILE, args.EMAIL, args.WORKINGDIR, args.NAME, args.INPUTTYPE,
					args.CLIENTS_FOLDER+"/", args.RUN_AF2, args.SITES_FILE)	

