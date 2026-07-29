import os
import sys
from site_selection_util import save_fasta, check_input_type, three_to_one
from site_selection_util import extract_bfactors_from_pdb, calc_sasa_shrakerupley, calc_hydrophobic_patches
from calculate_protein_scores import load_scores

#default hydrophobicity table, resolved relative to the repo so the CLI doesn't need a new required arg
DEFAULT_HYDRO_TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables",
									"hydrophobicity_kyte-doolittle.tsv")

def main(input_file, output_file, hydro_table=None, dilation_radius=4.0, min_seed_size=3, min_total_hydro=10.0):

	#output name = "xxx.x_plddt.txt"

	#we extract the pLDDT values from bfactors
	bf_out = open(output_file, "w")
#	bf = []
	bf = extract_bfactors_from_pdb(input_file)
	for aa in range(len(bf)):
		print("\t".join([str(aa+1), str(bf[aa])]), file=bf_out)
	bf_out.close()

	#now we extract SASA values; name matches companion_path convention (.txt → .sasa.txt)
	sasa_out = open(output_file.replace(".txt", ".sasa.txt"), "w")
	sasa = calc_sasa_shrakerupley(input_file)
	for s in sasa:
		print("\t".join([s[0], s[1], three_to_one(s[2])]), file=sasa_out)
	sasa_out.close()

	#now we identify hydrophobic surface patches (see calc_hydrophobic_patches docstring
	#for method reference); companions follow the same .txt-suffix-swap convention
	hydro_scores = load_scores(hydro_table or DEFAULT_HYDRO_TABLE)
	scores, patches = calc_hydrophobic_patches(input_file, hydro_scores,
												dilation_radius=dilation_radius,
												min_seed_size=min_seed_size)

	hydro_out = open(output_file.replace(".txt", ".hydro.txt"), "w")
	seq = [three_to_one(s[2]) for s in sasa]
	for i, score in enumerate(scores):
		print("\t".join([str(i+1), str(score), seq[i]]), file=hydro_out)
	hydro_out.close()

	patches_out = open(output_file.replace(".txt", ".patches.txt"), "w")
	n = 0
	for patch in patches:
		#skip weakly hydrophobic patches; total_hydro is the summed per-residue
		#Kyte-Doolittle score over seed residues, not an area
		if patch["total_hydro"] <= min_total_hydro:
			continue
		description = "patch{}|{:.2f}|area={:.1f}".format(n, patch["total_hydro"], patch["total_area_A2"])
		for resnum in patch["members"]:
			print("\t".join(["hydrophobic_patch", str(resnum), str(resnum), description]), file=patches_out)
		n += 1
	patches_out.close()

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
	parser.add_argument('--hydro_table', action = 'store', type = str, dest = 'HYDRO_TABLE', default = None,
					help = "tab-delimited per-amino-acid hydrophobicity table; defaults to the bundled Kyte-Doolittle table")
	parser.add_argument('--dilation_radius', action = 'store', type = float, dest = 'DILATION_RADIUS', default = 4.0,
					help = "distance (Angstroms) to grow each hydrophobic patch by, pulling in nearby residues regardless of their own exposure/hydrophobicity; 0 = no dilation")
	parser.add_argument('--min_seed_size', action = 'store', type = int, dest = 'MIN_SEED_SIZE', default = 3,
					help = "minimum hydrophobic surface residue count for a patch core to be reported")
	parser.add_argument('--min_total_hydro', action = 'store', type = float, dest = 'MIN_TOTAL_HYDRO', default = 10.0,
					help = "minimum summed Kyte-Doolittle hydrophobicity (dimensionless) over a patch's seed residues for it to be reported")

	args, unknowns = parser.parse_known_args()

#	outputfilename = args.WORKINGDIR + "/" + args.OUTPUT

	# no PDB means the AFDB lookup found no structure (or none was supplied);
	# this task's pLDDT/SASA/hydrophobicity outputs simply aren't available —
	# warn clearly and exit instead of crashing, so the rest of the pipeline's
	# independent tasks (blast, domains, modifications, ...) still complete.
	if not args.INPUT_PDB:
		print("No AlphaFold structure available: no PDB was provided and none was "
			  "found in the AlphaFold DB for this protein. Skipping pLDDT/SASA/"
			  "hydrophobicity analysis; other analyses are unaffected. To include "
			  "this analysis, supply a PDB file manually with --pdb.", file=sys.stderr)
		sys.exit(1)

	if check_input_type(args.INPUT_PDB) != "pdb":
		raise IOError("File type must be .pdb!")

	main(args.INPUT_PDB, args.OUTPUT, args.HYDRO_TABLE, args.DILATION_RADIUS, args.MIN_SEED_SIZE, args.MIN_TOTAL_HYDRO)

