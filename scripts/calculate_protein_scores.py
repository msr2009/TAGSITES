"""
calculate_protein_scores.py

Matt Rich, 4/2024
"""

from site_selection_util import *
import sys

def load_scores(scores_file):
	"""
	Load hydrophobicity scores from a tab-delimited file into a dictionary.

	Args:
		scores_file (str): Path to the tab-delimited file containing amino acid scores.

	Returns:
		dict: A dictionary containing amino acid scores.
	"""
	scores_dict = {}
	with open(scores_file, 'r') as file:
		for line in file:
			line = line.strip()
			if line:
				aa, score = line.split('\t')
				scores_dict[aa] = float(score)
		return scores_dict


def calculate_property_sliding_window(protein_sequence, window_size, scores):
    """
    Calculate the average hydrophobicity of a protein sequence using a sliding window approach.

    Args:
        protein_sequence (str): The protein sequence as a string of amino acid symbols.
        window_size (int): Size of the sliding window.
		scores (dict): amino acid scores to use in calculation

    Returns:
        list: A list containing average hydrophobicity values for each sliding window.
    """

    # Initialize list to store average hydrophobicity values for each window
    average_values = []

    # Iterate through the protein sequence with a sliding window
    for i in range(len(protein_sequence) - window_size + 1):
        window_sequence = protein_sequence[i:i+window_size]
        window_value = 0.0

        # Calculate hydrophobicity for amino acids in the current window
        for aa in window_sequence:
            if aa.upper() in scores:
                window_value += scores[aa.upper()]

        # Calculate the average hydrophobicity for the current window
        average_window_value = window_value / window_size
        average_values.append(average_window_value)

    return average_values

def main(fasta_in, fout, window, scores_file):

	#read fasta
	name, seq = read_fasta(fasta_in)

	scores_dict = load_scores(scores_file)

	#calculate scores
	h_scores = calculate_property_sliding_window(seq, window, scores_dict)

	for i in range(len(h_scores)):
		print("\t".join([str(i), str(h_scores[i])]), file=fout)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "fasta file containing sequence", required=True)
	parser.add_argument('-s', '--scores', action='store', type=str, dest='SCORES_FILE', 
		help = "tab-delimited file containing property scores for each amino acid", required=True)
	parser.add_argument('-w', '--window', action='store', type=int, dest='WINDOW', 
		help = "sliding window width, default = 11", default = 11)
	parser.add_argument('-o', '--output', action='store', type=str, dest='OUTFILE', 
		help = "output file, default=stdout", default = sys.stdout)
	args = parser.parse_args()
	
	if args.OUTFILE != sys.stdout:
		args.OUTFILE = open(args.OUTFILE, "w")

	main(args.FASTA_IN, args.OUTFILE, args.WINDOW, args.SCORES_FILE)	

	args.OUTFILE.close()

