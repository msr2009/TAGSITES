from site_selection_util import save_fasta, check_input_type
from site_selection_util import extract_bfactors_from_pdb, extract_seq_from_pdb

from existing_AF2_model import search_AFDB

def main(input_file, outputfile, evalue, percentid, email, name, workingdir, clients_folder, writefasta, searchAF2):
	
	bf_out = open(outputfile, "w")
	bf = []

	#we extract the pLDDT values from bfactors
	bf = extract_bfactors_from_pdb(input_file)
	
	for aa in range(len(bf)):
		print("\t".join([str(aa+1), bf[aa]]), bf_out)
		
	if writefasta:
		save_fasta(inputfile.split(".")[0:-1], 
					extract_seq_from_pdb(input_file),
					inputfile.replace(".pdb", ".fa"))

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	#required arguments
	parser.add_argument('--input_file', action = 'store', type = str, dest = "INPUT_FILE", 
					help = "AlphaFold2 PDB file. B-factors must contain pLDDTs for each residue.")
	parser.add_argument('--output', action = 'store', type=str, dest="OUTPUT",
					help = "output filename")
	parser.add_argument('--write_fasta', action = 'store_true', dest="WRITE_FASTA", 
					help="output FASTA of protein sequence from PDB file",
					default=False)
	parser.add_argument('--existing_AF2', action = 'store_true', dest="args.SEARCHAFDB", 
					help="search for matching protein in AFDB (req'd if input is fasta)",
					default=False)

	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--working_dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--run_name', action='store', type=str, dest='NAME', 
		help = "name for output", required=True)
	parser.add_argument('--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for BLAST hit (1e-100)", default=1e-100)
	parser.add_argument('--percentid', action='store', type=float, dest='PERCENTID', 
		help = "Identity threshold for BLAST hit (99)", default=99)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients (default=./scripts/)",
		default = "./scripts/")

	outputfilename = args.WORKINGDIR + "/" args.OUTPUT

	if check_input_type(args.INPUT_FILE) != "pdb":
		raise IOError("File type must be .pdb!")
		args.SEARCHAFDB = True

	args = parser.parse_args(args.INPUT_FILE, args.outputfilename, 
							 args.EVALUE, args.PERCENTID, 
							 args.EMAIL, args.RUN_NAME, args.WORKINGDIR, args.CLIENTS_FOLDER, 
							 args.WRITE_FASTA)

	
