import re, sys

"""
regex_sites.py

uses regular expressions to identify sites in a protein sequence that may be
functional

Matt Rich, 4/2024
"""

def main(seq, site_res, fout):
	site_ranges = []

	for line in open(site_res, "r"):
		l = line.strip().split("\t")
		desc = l[0] 
		regex = l[1]
		for match in re.finditer(regex, seq):
			site_ranges = [desc, match.start(), match.end()]
			print("\t".join(str(x) for x in site_ranges), file=fout)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-s', '--seq', action = 'store', type = str, dest = 'SEQ', 
		help = "protein sequence")
	parser.add_argument('-f', action = 'store', type = str, dest = 'SITE_FILE', 
		help = "tab-delimited file containing name and regular expression for site")
	parser.add_argument('-o', action = 'store', type = str, dest = 'OUTPUT_FILE', 
		help = "file to store output, default=STDOUT", default = None)
	args = parser.parse_args()
	
	if args.OUTPUT_FILE != None:
		fout = open(args.OUTPUT_FILE, "w")		
		main(args.SEQ, args.SITE_FILE, args.OUTPUT_FILE)
		fout.close()
	else:
		main(args.SEQ, args.SITE_FILE, sys.stdout)
		

