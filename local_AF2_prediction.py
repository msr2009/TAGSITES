"""
local_AF2_prediction.py

A wrapper for calling colabfold for AlphaFold predictions.

REQUIRES: local installation of localcolabfold
			https://github.com/YoshitakaMo/localcolabfold

Parameters:
	- fasta_in (str): name of fasta file containing seq

Returns: 
	- localcolabfold output: folder containing lots of files, including relaxed
	  PDBs.

Matt Rich, 4/2024
"""

import subprocess, sys
from site_selection_util import read_fasta

def check_localcolabfold_install():
	check = subprocess.run("colabfold_batch -h", shell=True, 
							capture_output=True, text=True).stdout
	if check.startswith("usage: colabfold_batch") != True:
		sys.exit("LocalColabfold installation not found. Exiting.")

def main(fasta_in, workingdir):
	
	#first, check if localcolabfold exists
	check_localcolabfold_install()

	#read sequence from fasta
	name, seq = read_fasta(fasta_in)

	#make localcolabfold command
	results_folder = "{}/{}_AF2/".format(workingdir, name)
	lcf_call = "colabfold_batch --num-models 5 \
								--recycle-early-stop-tolerance 0.1 \
								--amber \
								--use-gpu-relax \
								{} \
								{}".format(fasta_in, results_folder)

	#call localcolabfold command
	subprocess.run(lcf_call, shell=True)

	#move best PDB model to workingdir
	subprocess.run("cp {}/*_relaxed_rank_001*.pdb {}".format(results_folder, workingdir), shell=True)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)

	args = parser.parse_args()
	
	main(args.FASTA_IN, args.WORKINGDIR)	

