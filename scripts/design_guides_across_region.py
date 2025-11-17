"""
design_guides_across_region.py

script that designs CRISPR guide rnas in region and will filter based on blast 

requires bedtools. optionally: local blast installation

Matt Rich, 2023
"""

import subprocess, re, sys
from Bio import SeqIO

def get_regions_by_gene(gene, extra):
	#call grep to get regions
	#this is a bit hacky -- it will only work when gene names are at the end of the line
	#this is true for wormbase gff3's though. 
	grep_call =	"grep WormBase {} | grep CDS | grep {}$".format(args.GFF, gene)
	grep_output = subprocess.run(grep_call, shell=True, capture_output=True, text=True)
	regions = [x.split("\t") for x in grep_output.stdout.split("\n")][0:-1]
	
	return [x[0] + ":" + str(int(x[3])-extra) + "-" + str(int(x[4])+extra) for x in regions]

def merge_regions(regions):
	merged_regions = []
	#make bed file for intersecting
	merge_bed = open("tmp_merge.bed","w")
	for r in sorted(regions):
		print(convert_region_to_BED(r), file=merge_bed)
	merge_bed.close()

	merge_call = "bedtools merge -i tmp_merge.bed"
	merged_output = subprocess.run(merge_call, shell=True, capture_output=True, text=True).stdout
	
	for m in [x.strip().split("\t") for x in merged_output.strip().split("\n")]:
		merged_regions.append("{}:{}-{}".format(m[0], m[1], m[2]))
	return merged_regions
		
def call_bedtools_getfasta(region, fasta, strand="+"):
	bed_in = "\t".join([convert_region_to_BED(region), "foo", "0", strand])
#	print("BED_IN: " + bed_in)
	gf_call = 'bedtools getfasta -s -tab -fi {} -bed <(echo "{}")'.format(fasta, bed_in)
#	print("GF_CALL: " + gf_call)
	gf_output = subprocess.run(gf_call, shell=True, capture_output=True, 
								text=True, executable='/bin/bash')
	return(gf_output.stdout.strip().split("\t"))

def convert_region_to_BED(region):
	return region.strip("\n").replace(":", "\t").replace("-", "\t")

#rewrite to maintain region info for sorting
def find_targets(seq, seqstart, seqend, rc=False):
	
	if rc:
		seq = reverse_complement(seq)
		return [ [seq[s.start()-args.GUIDE_LENGTH:s.start()],
				  seqend-s.start()+args.GUIDE_LENGTH, 
				  seqend-s.start()+1] for s in re.finditer(r'(?=({}))'.format(iupac_to_regex(args.PAM)), seq) ]	
	else:
		return [ [seq[s.start()-args.GUIDE_LENGTH:s.start()],
				  seqstart+s.start()-args.GUIDE_LENGTH,
				  seqstart+s.start()-1] for s in re.finditer(r'(?=({}))'.format(iupac_to_regex(args.PAM)), seq) ]	

def filter_target(target):
	filters = []
	notes = []
	#no 3' C
	if target[-1] == "C":
		filters.append("3'C")
	#if blast, then do blast search, filter for any with more than one hit
	if args.USE_BLAST:
		blast_results = blast_guide("tmp", target, args.BLASTDB, args.EVAL, args.PAM, args.MINLENGTH, args.MISMATCH, args.FASTA)
		if len(blast_results) > 1:
			filters.append("OFFTARGETS")
			notes.append(",".join(blast_results))
	if len(filters) == 0:
		filters.append("PASS")
	return ";".join(filters), ";".join(notes)

def blast_guide(seqname, seq, db, e_value, pam, min_length, mismatch, fasta):
	#need to make a fasta file for input
	tmp_fasta = open("tmp.fasta", "w")
	print(">{}\n{}".format(seqname, seq), file=tmp_fasta)
	tmp_fasta.close()

	#wormbase parameters
	word_size = 5
	gapopen = 5
	gapextend = 2
	penalty = -3	

	#write blast call
	blast_call = "blastn -query tmp.fasta -db {} -evalue {} -word_size {} -gapopen {} -gapextend {} -penalty {} -outfmt 6".format(db, e_value, word_size, gapopen, gapextend, penalty)
	blast_output = subprocess.run(blast_call, shell=True, text=True,
									capture_output=True).stdout.strip().split("\n")

#	print(blast_output)

	#if there are multiple blast hits, see if any of them have PAMs
	ontargets = []
	if len(blast_output) != 1:
		for h in blast_output:
			l = h.strip().split("\t")
			if int(l[7]) == len(seq) and int(l[7])-int(l[6])+1 >= min_length and int(l[4]) <= mismatch:
				pam_check = False
				if int(l[8]) < int(l[9]):
					pam_check = check_pam(l[1], int(l[8])-1, int(l[9])+3, "+")
				else:
					pam_check = check_pam(l[1], int(l[9])-4, int(l[8]), "-")
				if pam_check:
					ontargets.append("{}:{}-{}".format(l[1], l[8], l[9]))
#	subprocess.run("rm tmp.fasta")
	return ontargets

def check_pam(chrom, start, end, strand):
	#print("{}:{}-{}".format(chrom, start, end))
	hit_fa = call_bedtools_getfasta("{}:{}-{}".format(chrom, start, end), args.FASTA, strand)
#	print(hit_fa)
	if re.search("{}$".format(iupac_to_regex(args.PAM)), hit_fa[1]) != None:
		return True
	else:
		return False

def determine_input_type(s):
	if s.endswith(".fa") or s.endswith(".fasta"):
		return "fasta"
	elif re.match(r".+:\d+-\d+", s): 
		return "region"
	elif re.match(r"^[actgACTG]+$", s):
		return "seq"
	elif re.match(r"^[A-Za-z\d\-.]+$", s):
		return "gene"
	else:
		raise ValueError("Unknown input type. Please define with --input-type")

def iupac_to_regex(iupac_nucleotides):
	iupac_dict = {
		'A': 'A',
		'C': 'C',
		'G': 'G',
		'T': 'T',
		'U': 'U',  # RNA Thymine
		'R': '[AG]',  # Purine (A or G)
		'Y': '[CT]',  # Pyrimidine (C or T)
		'S': '[GC]',  # Strong (G or C)
		'W': '[AT]',  # Weak (A or T)
		'K': '[GT]',  # Keto (G or T)
		'M': '[AC]',  # Amino (A or C)
		'B': '[CGT]',  # Not A (C or G or T)
		'D': '[AGT]',  # Not C (A or G or T)
		'H': '[ACT]',  # Not G (A or C or T)
		'V': '[ACG]',  # Not T (A or C or G)
		'N': '[ACGT]'  # Any nucleotide
	}
	regex_patterns = [iupac_dict.get(nucleotide, '') for nucleotide in iupac_nucleotides]
	return ''.join(regex_patterns)

def reverse_complement(seq):
	dna_dict = {
		"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"
	}
	return "".join([dna_dict[x.upper()] for x in seq[::-1]])

#a little error for BLASTdbs not being real
class BLASTdbError(Exception):
	pass

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	#required arguements
	parser.add_argument('-i', '--input', action = 'store', type = str, dest = "INPUT",
		help = "input. Can be sequence, gene name, or region (e.g., I:123-456)",
		required=True)
	parser.add_argument("--input-type", action = 'store', type = str, dest = "INPUTTYPE", 
		help = "input type. must be 'seq', 'region', 'gene', or 'fasta'", 
		choices=["seq", "region", "gene", "fasta"])
	parser.add_argument('--gff', action = 'store', type = str, dest = 'GFF',
		help = "GFF3 file containing exons and gene names as locus field")
	parser.add_argument('-f', '--fasta', action = 'store', type = str, dest = "FASTA", 
		help = "fasta file")
	
	#blast args
	parser.add_argument('--blast', action = 'store_true', dest = "USE_BLAST", 
		help = "use BLAST to search against genome for multiple guide hits? requires --blastdb")
	parser.add_argument('--blastdb', action = 'store', type = str, dest = "BLASTDB",
		help = "location of blastdb to search against")
	parser.add_argument('--minlength', action = 'store', type = int, dest = "MINLENGTH",
		help = "minimum length count BLAST match with PAM as hit. default=17", default = 17)
	parser.add_argument('--mismatch', action = 'store', type = int, dest = "MISMATCH", 
		help = "maxmimum mismatches allowed in potential guide hit. default=0", default = 0)
	parser.add_argument('--e_value', action = 'store', type = float, dest = "EVAL",
		help = "BLAST E-value cutoff", default = 1)
	
	#guide args
	parser.add_argument('--PAM', action = 'store', type = str, dest = "PAM", 
		help = "PAM sequence (default = NGG)", default = "NGG")
	parser.add_argument('--guide_length', action = 'store', type = int, dest = "GUIDE_LENGTH", 
		help = "length of guide target, default = 20", default=20)
	parser.add_argument('--match_regex', action = 'store', type = str, dest = "MATCH_REGEX",
		help = "check if guide matches regular expression, e.g., for base-editing.", default=None)

	#sequence args
	parser.add_argument('-b', '--prefix', action = 'store', type = str, dest = "PREFIX", 
		help = "prefix sequence to add to guide (e.g., for cloning)", default="")
	parser.add_argument('-a', '--suffix', action = 'store', type = str, dest = "SUFFIX", 
		help = "suffix sequence to add to guide (e.g., for cloning)", default="")

	args = parser.parse_args()
	
	
	#try to figure out what the input is
	input_type = ""
	if args.INPUTTYPE != None:
		input_type = args.INPUTTYPE
	else:
		input_type = determine_input_type(args.INPUT)

	print("input type detected as {}".format(input_type), file=sys.stderr)

	if input_type not in ["seq", "fasta"]:
		if args.FASTA == None:
			raise IOError("genomic FASTA file required for 'gene' or 'region' input!")

	if args.USE_BLAST:
		#check if BLAST db exists, because otherwise this'll just run and not
		#check anything....
		check_blast_call = "blastdbcmd -db {} -info".format(args.BLASTDB)
		check_blast = subprocess.run(check_blast_call, shell=True, text=True,
						capture_output=True)
		if check_blast.stdout.startswith("Database:"):
			print("BLAST db found. ({})".format(args.BLASTDB), file=sys.stderr)
		else:
			raise BLASTdbError("BLAST database not found! Either run without --blast or confirm location.")

	reg = []
	if input_type == "gene":
		all_gene_reg = get_regions_by_gene(args.INPUT,args.GUIDE_LENGTH+len(args.PAM))
		#call bedtools merge to get intersection of all regions
		reg = merge_regions(all_gene_reg)
	elif input_type == "region":
		reg = [args.INPUT]
	
	region_seqs = []
	if input_type == "seq":
		region_seqs = [["SEQ:1-{}()".format(args.GUIDE_LENGTH), args.INPUT]]
	elif input_type == "fasta":
		#import fasta sequence using biopython
		record = [r for r in SeqIO.parse(open(args.INPUT,"r"), "fasta")][0]
		region_seqs = [["{}:1-{}()".format(record.id, len(record.seq)), str(record.seq)]]
	else:
		region_seqs = [call_bedtools_getfasta(r, args.FASTA) for r in reg] 

	#print header
	print("\t".join(["chrom", "start", "stop", "spacer", "fullseq", "filter", "notes"]))

	for r in region_seqs:
		potential_targets = []
		chrom = r[0].split(":")[0]
		start, stop = r[0].split(":")[1].split("-")			
		stop = stop.split("(")[0]
	
		if input_type == "seq":
			potential_targets = [[r[1], "1", "20"]]
		else:	
			potential_targets = find_targets(r[1], int(start), int(stop)) + \
								find_targets(r[1], int(start), int(stop), rc=True)
		
		for pt in potential_targets:
			if len(pt[0]) == 0:
					continue

			filter_str, filter_notes = filter_target(pt[0])

			#check for match to regular expression
			if args.MATCH_REGEX != None:
#				print(args.MATCH_REGEX, re.search(pt[0], args.MATCH_REGEX))
				if re.search(args.MATCH_REGEX, pt[0]) != None:
					filter_notes += "RE_MATCH:{}".format(args.MATCH_REGEX)
			
			#print output for each guide
			print("\t".join([chrom, str(pt[1]), str(pt[2]), pt[0],
					args.PREFIX+pt[0]+args.SUFFIX, filter_str,
					";".join([filter_notes])]), file=sys.stdout)

