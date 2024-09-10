"""
existing_AF2_model.py

A script to check whether the AlphaFold database has a predicted structure
for your protein sequence. 

Uses EBI's API (ncbiblast and dbfetch). 

1) BLASTs sequence against Uniprot to find accession. Will take best match above thresholds.
2) Tries to fetch PDB file from AF DB for that accession.

Parameters:
	- fasta_in (str): name of fasta file containing seq
	- evalue (float): evalue threshold for match
	- percentid (float): %ID threshold for match

Returns: 
	- pdb_out (str): raw string of PDB file (if found) or "pdb not found"

Matt Rich, 4/2024
"""

import subprocess, re
from site_selection_util import get_sequence, uniprot_accession_regex, save_fasta

def main(fasta_in, email, workingdir, name, evalue, percentid, clients_folder):
	
	match_accession = ""
	match_eval = 1e-200
	match_id = 100.0
	move_fasta = False

	#if the input is an accession
	if uniprot_accession_regex(fasta_in) != None:
		match_accession = fasta_in
		move_fasta = True
	#otherwise, we need to find a match in the AF2 db
	else:
		#read sequence from fasta
		seq = get_sequence(fasta_in)

		#make ncbiblast command
		ncbi_call = "python {}ncbiblast.py \
						--email {} \
						--program blastp \
						--stype protein \
						--sequence {} \
						--database afdb \
						--outformat tsv \
						--outfile {}/{}.ncbiblast".format(clients_folder, email, seq, workingdir, name)
		#call ncbiblast command
		subprocess.run(ncbi_call, shell=True)

		#open ncbiblast output
		ncbi_out = open("{}/{}.ncbiblast.tsv.tsv".format(workingdir, name), "r")
		#we only care about the best match (the second line)
		l = ncbi_out.readlines()[1].strip().split("\t")
		match_id = float(l[7])
		match_eval = float(l[9])
		match_accession = l[2]

	#download AF2 model if it exists
	if match_id >= percentid and match_eval <= evalue and match_accession != "":
		#if we have a close match, then see if there's an AF2 model
		#this uses dbfetch from EBI
		dbfetch_call = "python {}dbfetch.py fetchData \
						afdb:{} pdb raw".format(clients_folder, match_accession)
	
		#call dbfetch command
		pdb_out = subprocess.run(dbfetch_call, shell=True, capture_output=True, text=True).stdout

		#if this call returns an error, there isn't a predicted structure
		no_af_error = "ERROR 11 Unable to connect to database [afdb]."
		if pdb_out != no_af_error:
			#otherwise, we found a good pdb
			f_out = open("{}/{}.pdb".format(workingdir, name), "w")
			#save it in the working dir
			print(pdb_out, file=f_out)
			print("saving {} pdb file to {}/{}.pdb".format(match_accession, workingdir, name))
			f_out.close()
			
			#also, make a fasta from the pdb if we went straight from accession
			if move_fasta:
				save_fasta(name, 
							get_sequence("{}/{}.pdb".format(workingdir, name)), 
							"{}/{}.fa".format(workingdir, name))
			print("saving fasta from {} pdb to {}/{}.fa".format(match_accession, workingdir, name))

			return 0
		else:
			print("pdb not found")
			return 1

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--name', action='store', type=str, dest='NAME', 
		help = "name for output", required=True)
	parser.add_argument('--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for BLAST hit (1e-20)", default=1e-50)
	parser.add_argument('--percentid', action='store', type=float, dest='PERCENTID', 
		help = "Identity threshold for BLAST hit (95)", default=99)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")

	args = parser.parse_args()
	
	main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME, args.EVALUE,
					args.PERCENTID, args.CLIENTS_FOLDER)	

