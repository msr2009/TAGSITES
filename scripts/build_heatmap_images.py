"""
build_alignment_images.py

takes fasta alignment file (output of clustal) and 
creates matplotlib heatmap for alignment.

this is an alternative to using plotly to make an interactive heatmap, 
which runs really slow in the tag sites app.

Matt Rich, 11/2025
"""

import matplotlib
matplotlib.use("agg")  # non-interactive backend; required when called from worker threads
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# Define a basic ClustalX-style color scheme for amino acids
CLUSTALX_COLOR_SCHEME = {
        'A': 'green',  'I': 'green',  'L': 'green',  'M': 'green',  'V': 'green',  # Hydrophobic
        'F': 'orange', 'Y': 'orange', 'W': 'orange',  # Aromatic
        'H': 'blue',   'K': 'blue',   'R': 'blue',  # Positively charged
        'D': 'red',    'E': 'red',    # Negatively charged
        'S': 'magenta', 'T': 'magenta',  # Polar, uncharged
        'N': 'cyan',   'Q': 'cyan',    # Polar, uncharged
        'G': 'yellow', 'P': 'yellow', 'C': 'yellow', # Special cases
        "-": 'white', 'X': 'white' # Gap
}

AA_VAL_MAP = { 'A': 0, 'I': 7,  'L': 9,  'M': 10,  'V': 17,  # Hydrophobic
        'F': 4, 'Y': 19, 'W': 18,  # Aromatic
        'H': 6, 'K': 8, 'R': 14,  # Positively charged
        'D': 2, 'E': 3,    # Negatively charged
        'S': 15, 'T': 16,  # Polar, uncharged
        'N': 11,   'Q': 13,    # Polar, uncharged
        'G': 5, 'P': 12, 'C': 1, # Special cases
        '-': 20, 'X': 21 # Gap
}

AA_NUM_MAP = { v:k for k,v in AA_VAL_MAP.items() } 

CLUSTAL_CMAP = ListedColormap([CLUSTALX_COLOR_SCHEME[AA_NUM_MAP[x]] for x in sorted(AA_NUM_MAP.keys())])

def parse_fasta_alignment(fasta_file):
	sequences = []
	headers = []
	
	with open(fasta_file, 'r') as handle:
		for record in SeqIO.parse(handle, "fasta"):
				sequences.append(list(record.seq))
				headers.append(record.id)
	return headers, sequences

#function to convert seq matrix and aa:color mapping dict
#we need to do this this way because we're using a heatmap 
#to plot the alignment
def color_seq_for_heatmap_bydict(seq_mat, color_mapping, nan_color="white"):
    seq_mat_colors = []
    for seq in seq_mat:
        seq_colors = []
        for res in seq:
            seq_colors.append(color_mapping[res])
        seq_mat_colors.append(seq_colors)

    return seq_mat_colors

#function to make mapping for x ticks based on subject sequence
def x_ticks_by_subject(fasta_file, increment=10):
    headers, sequences = parse_fasta_alignment(fasta_file)
    seq = sequences[0]

    ticks_loc = []
    tick_labs = []
    ticks_counter = 0
    for i in range(len(seq)):
        if seq[i] != "-":
            ticks_counter += 1

        if ticks_counter % increment == 0 and ticks_counter != 0:
            ticks_loc.append(i)
            tick_labs.append(ticks_counter)

    return (ticks_loc, tick_labs)

def plot_alignment_matrix_matplotlib(fasta_file, outfmt, outfile):
    #import fasta
	headers, sequences = parse_fasta_alignment(fasta_file)

	#update output filename
	if outfile == "":
		outfile = fasta_file.removesuffix("aln")+outfmt
	print(outfile)


	#convert sequences into numbers ("A"=0..."-"=20, etc)
	dat_colors = color_seq_for_heatmap_bydict(sequences, AA_VAL_MAP)
	num_sequences = len(sequences)
	alignment_length = len(sequences[0])

	#fig size based on length and number of seqs
	fig, ax = plt.subplots(figsize=(alignment_length/8, num_sequences))

	fs = 6 #fontsize

	ax.imshow(dat_colors, cmap=CLUSTAL_CMAP, aspect="equal")
	ax.set_yticks(range(len(headers)), labels=headers, fontsize=fs)
	
	#set x-ticks based on position in first (top) sequence
	xticks = x_ticks_by_subject(fasta_file, 20)
	ax.set_xticks(xticks[0], xticks[1])

	ax.set_title(" ".join(fasta_file.removesuffix(".aln").split(".")))
	ax.set_xlabel("position (aa)")
    
	#annotate text for each aa
	for i in range(len(sequences)):
		for j in range(len(sequences[i])):
			text = ax.text(j, i, sequences[i][j], ha="center", va="center", color="k", fontsize=fs)
    
	plt.tight_layout()
	plt.savefig(outfile, bbox_inches="tight")


if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-a', '--aln', action = 'store', type = str, dest = 'ALN', 
		help = "fasta file containing alignment", required=True)
	parser.add_argument('-f', '--format', action = 'store', type = str, dest = 'FMT', 
		help = "format for saving figure (svg)", default="svg")
	parser.add_argument('-o', '--output', action = 'store', type = str, dest = 'OUT', 
		help = "output file name. overrides automatic output naming.", default="")
	args = parser.parse_args()
	
	plot_alignment_matrix_matplotlib(args.ALN, args.FMT, args.OUT)	

