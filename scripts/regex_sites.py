import re, sys

"""
regex_sites.py

uses regular expressions to identify sites in a protein sequence that may be
functional

Matt Rich, 4/2024
"""

from site_selection_util import read_fasta
from progress import report as _report, resolve_reporter

def main(fasta_in, site_res, fout, report=None):

	reporter = resolve_reporter(report)

	#read sequence from fasta
	name, seq = read_fasta(fasta_in)

	site_ranges = []
	counts = {}

	for line in open(site_res, "r"):
		l = line.strip().split("\t")
		description = l[0]
		regex = l[1]
		for match in re.finditer(regex, str(seq)):
			site_ranges = ["modification", match.start()+1, match.end()+1, description]
			print("\t".join(str(x) for x in site_ranges), file=fout)
			counts[description] = counts.get(description, 0) + 1

	total = sum(counts.values())
	if total:
		breakdown = ", ".join(f"{n} {desc}" for desc, n in counts.items())
		summary = f"Found {total} modification site(s): {breakdown}"
	else:
		summary = "Found 0 modification sites"
	_report(reporter, summary, stage="done")

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', '--input_file', action = 'store', type = str, dest = 'FASTA', 
		help = "fasta file containing protein sequence", required=True)
	parser.add_argument('--sites_file', action = 'store', type = str, dest = 'SITES_FILE', 
		help = "tab-delimited file containing name and regular expression for site", required=True)
	parser.add_argument('--output', action = 'store', type = str, dest = 'OUTPUT_FILE', 
		help = "file to store output, default=STDOUT", default = None)

	args, unknowns = parser.parse_known_args()
	
	if args.OUTPUT_FILE != None:
		fout = open(args.OUTPUT_FILE, "w")		
		main(args.FASTA, args.SITES_FILE, fout)
		fout.close()
	else:
		main(args.FASTA, args.SITES_FILE, sys.stdout)
		

