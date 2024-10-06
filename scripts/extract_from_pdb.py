from site_selection_util import save_fasta, check_input_type
from site_selection_util import extract_bfactors_from_pdb, extract_seq_from_pdb

def main(input_file, output_file):
	
	bf_out = open(output_file, "w")
	bf = []

	#we extract the pLDDT values from bfactors
	bf = extract_bfactors_from_pdb(input_file)
	
	for aa in range(len(bf)):
		print("\t".join([str(aa+1), str(bf[aa])]), file=bf_out)
	
	return 0
		
if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	#required arguments
	parser.add_argument('--pdb', action = 'store', type = str, dest = "INPUT_PDB", 
					help = "AlphaFold2 PDB file. B-factors must contain pLDDTs for each residue.")
	parser.add_argument('--output', action = 'store', type=str, dest="OUTPUT",
					help = "output filename")
	parser.add_argument('--working_dir', action = 'store', type = str, dest = 'WORKINGDIR', 
					help = "working directory")

	args, unknowns = parser.parse_known_args() 

#	outputfilename = args.WORKINGDIR + "/" + args.OUTPUT

	if check_input_type(args.INPUT_PDB) != "pdb":
		raise IOError("File type must be .pdb!")

	main(args.INPUT_PDB, args.OUTPUT)
	
