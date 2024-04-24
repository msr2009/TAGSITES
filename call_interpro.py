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

import subprocess
from site_selection_util import read_fasta

def main(fasta_in, email, workingdir, clients_folder):
	
	#read sequence from fasta
	name, seq = read_fasta(fasta_in)

	#make interpro command
	interpro_call = "python {}iprscan5.py \
					--email {} \
					--stype p \
					--sequence {} \
					--outformat tsv \
					--outfile {}/{}.interpro".format(clients_folder, email, seq, workingdir, name)
	#call interpro command
	subprocess.run(interpro_call, shell=True)

	#process interpro output
	for line in open("{}/{}.interpro.tsv.tsv".format(workingdir, name), "r").readlines():
		l = line.strip().split("\t")
		#print output
		f_out = open("{}/{}.interpro.out", "w")
		print("\t".join([l[3], l[6], l[7], l[5]]), file=f_out)
		f_out.close()

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER',
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")

	args = parser.parse_args()
	
	main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.CLIENTS_FOLDER)	

