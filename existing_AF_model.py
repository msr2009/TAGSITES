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

import subprocess
from site_selection_util import read_fasta

#EBI_CLIENTS = "./ebi_api_clients/"

def main(fasta_in, email, workingdir, evalue, percentid, clients_folder):
	
	#read sequence from fasta
	name, seq = read_fasta(fasta_in)

	#make ncbiblast command
	ncbi_call = "python {}ncbiblast.py \
					--email {} \
					--program blastp \
					--stype protein \
					--sequence {} \
					--database uniprotkb \
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

	#check if there's a close match
	if match_id >= percentid and match_eval <= evalue:
		#if we have a close match, then see if there's an AF2 model
		#this uses dbfetch from EBI
		dbfetch_call = "python {}dbfetch.py fetchData \
							afdb:{} pdb raw".format(clients_folder, match_accession)

		#call dbfetch command
		pdb_out = subprocess.run(dbfetch_call, shell=True, capture_output=True, text=True).stdout

		#if this call returns an error, there isn't a predicted structure
		no_af_error = "ERROR 11 Unable to connect to database [afdb]."
		if pdb_out != no_af_error:
			f_out = open("{}/{}.pdb".format(workingdir, name), "w")
			print(pdb_out, file=f_out)
			print("saving pdb file to {}/{}.pdb".format(workingdir, name))
			f_out.close()
		else:
			print("pdb not found")

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for BLAST hit (1e-20)", default=1e-20)
	parser.add_argument('--percentid', action='store', type=float, dest='PERCENTID', 
		help = "Identity threshold for BLAST hit (95)", default=95)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")

	args = parser.parse_args()
	
	main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.EVALUE,
					args.PERCENTID, args.CLIENTS_FOLDER)	

