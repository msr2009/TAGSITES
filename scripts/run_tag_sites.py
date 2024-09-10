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

def main(infile, email, workingdir, name, inputtype, clients_folder, 
				run_af2, sites_file, close_db, close_nhits, close_eval,
				close_hitlength, close_align_full):
	
	#make workingdir
	if os.path.exists(workingdir) == False:
		print("creating working directory: {}".format(workingdir), file=sys.stderr)
		os.mkdir(workingdir)
	else:
		print("working directory exists: {}".format(workingdir), file=sys.stderr)

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
		
		print("input is {}, checking if there's an existing AF model".format(inputtype), file=sys.stderr)

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
		print("input is {}, assuming this is an AF2 model w/ pLDDT scores".format(inputtype), file=stderr)
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

	#conservation calls
	_close_align_full = ""
	if close_align_full:
		_close_align_full = "--align-full-sequence"
	blast_close_call = "python blast_close_orthologs.py -f {} --email {} \
														--dir {} \
														--clients-folder {} \
														--name {} \
														--evalue {} \
														-n {}\
														--db {} \
														{}".format(
														fastafile, email,
														workingdir,
														clients_folder, name,
														1e-50, close_nhits,
														close_db, 
														close_taxids, _close_align_full)

	commands.append(blast_close_call)




	#protein properties calls
	hydro_call = "python calculate_protein_scores.py -f {}/{}.fa -s kyle-doolittle.tsv -o {}/{}.hydrophob.tsv".format(workingdir, name, workingdir, name)
	commands.append(hydro_call)

	print("#################\nRUNNING ANALYSES\n#################")
	print("\n".join(commands))
	print("######################################################")

	#run the commands
	procs = [ subprocess.Popen(i, shell=True) for i in commands ]
	for p in procs:
		p.wait()
	

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	req_args = parser.add_argument_group("required arguments")
	req_args.add_argument('-i', '--input', action='store', type=str, dest='INFILE', 
		help = "input file. can be fasta or pdb (determined automatically).", required=True)
	req_args.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	req_args.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	req_args.add_argument('-n', '--name', action='store', type=str, dest='NAME', 
		help = "job name for files, fasta header, etc...", required=True) 

	other_admin = parser.add_argument_group("other administrative arguments")
	other_admin.add_argument('--input-type', action='store', type=str, dest='INPUTTYPE', 
		help = "input type (fasta, pdb).", default=None, choices=["fasta", "pdb"])
	other_admin.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER',
		help = "path to EBI webservice clients (default=current directory)",
		default = "./")

	af_args = parser.add_argument_group("alphafold arguments")
	af_args.add_argument('--run-af2', action='store_true', dest='RUN_AF2', 
		help = "run alphafold locally. requires localcolabfold (default=False).", default=False)

	cons_args = parser.add_argument_group("conservation analysis arguments")
	cons_args.add_argument('--db', action='store', type=str, dest='BLAST_DB', 
		help = "database to blast for orthologs (default=uniprotkb)",
		default="uniprotkb")
	cons_args.add_argument('--n-hits', action='store', type=int, dest='N_HITS', 
		help = "maximum number of blast hits to return (default=100)",
		default=100)
	cons_args.add_argument('-e', '--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for keeping hits in analysis (1e-10)",
		default=1e-10)
	cons_args.add_argument('-l', '--length', action='store', type=float, dest='LENGTH', 
		help = "minimum length of BLAST matches, as percent of input (0.7)",
		default=0.7)
	cons_args.add_argument('--align-full-sequence', action='store_true', dest='ALIGN_FULL', 
		help = "fetch and align full protein sequence? by default will use BLAST alignment.", default=False)
	cons_args.add_argument('--taxids', action='store', dest='TAXIDS', 
		help = "file containing taxonomic IDs for uniprot to search for BLAST hits.", 
		default = None)
	

	pred_args = parser.add_argument_group("functional site prediction arguments")
	pred_args.add_argument('--sites-file', action='store', type=str, dest='SITES_FILE',
		help = "tab-delimited file containing regular expressions for modification sites (default=modifications.txt)", 
		default="modifications.txt")
	
	args = parser.parse_args()
	
	main(args.INFILE, args.EMAIL, args.WORKINGDIR, args.NAME, args.INPUTTYPE,
					args.CLIENTS_FOLDER+"/", args.RUN_AF2, args.SITES_FILE, 
					args.BLAST_DB, args.N_HITS, args.EVALUE, args.LENGTH,
					args.ALIGN_FULL)	

