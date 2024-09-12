import re, sys

"""
regex_sites.py

uses regular expressions to identify sites in a protein sequence that may be
functional

Matt Rich, 4/2024
"""

from site_selection_util import read_fasta

def main(fasta_in, site_res, fout):

	#read sequence from fasta
	name, seq = read_fasta(fasta_in)
	
	site_ranges = []

	for line in open(site_res, "r"):
		l = line.strip().split("\t")
		description = l[0] 
		regex = l[1]
		for match in re.finditer(regex, str(seq)):
			site_ranges = [description, match.start(), match.end()]
			print("\t".join(str(x) for x in site_ranges), file=fout)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', '--input_file', action = 'store', type = str, dest = 'FASTA', 
		help = "fasta file containing protein sequence", required=True)
	parser.add_argument('--sites_file', action = 'store', type = str, dest = 'SITES_FILE', 
		help = "tab-delimited file containing name and regular expression for site", required=True)
	parser.add_argument('--output', action = 'store', type = str, dest = 'OUTPUT_FILE', 
		help = "file to store output, default=STDOUT", default = None)

	
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', '--working_dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--name', '--run_name', action='store', type=str, dest='NAME', 
		help = "prefix name for output", required=True)

	

	args = parser.parse_args()
	
	if args.OUTPUT_FILE != None:
		fout = open(args.OUTPUT_FILE, "w")		
		main(args.FASTA, args.SITES_FILE, fout)
		fout.close()
	else:
		main(args.FASTA, args.SITES_FILE, sys.stdout)
		

