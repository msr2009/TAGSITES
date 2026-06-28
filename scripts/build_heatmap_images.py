"""
build_heatmap_images.py — ClustalX-style alignment heatmap renderer.
Output is PNG at 96 DPI with per-cell amino acid letter annotations.
Matt Rich, 11/2025
"""

import matplotlib
matplotlib.use("agg")  # non-interactive backend; required when called from worker threads
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# ClustalX-style color scheme for amino acids
CLUSTALX_COLOR_SCHEME = {
    'A': 'green',  'I': 'green',  'L': 'green',  'M': 'green',  'V': 'green',
    'F': 'orange', 'Y': 'orange', 'W': 'orange',
    'H': 'blue',   'K': 'blue',   'R': 'blue',
    'D': 'red',    'E': 'red',
    'S': 'magenta', 'T': 'magenta',
    'N': 'cyan',   'Q': 'cyan',
    'G': 'yellow', 'P': 'yellow', 'C': 'yellow',
    '-': 'white',  'X': 'white',
}

AA_VAL_MAP = {
    'A': 0,  'I': 7,  'L': 9,  'M': 10, 'V': 17,
    'F': 4,  'Y': 19, 'W': 18,
    'H': 6,  'K': 8,  'R': 14,
    'D': 2,  'E': 3,
    'S': 15, 'T': 16,
    'N': 11, 'Q': 13,
    'G': 5,  'P': 12, 'C': 1,
    '-': 20, 'X': 21,
}

AA_NUM_MAP   = {v: k for k, v in AA_VAL_MAP.items()}
CLUSTAL_CMAP = ListedColormap([CLUSTALX_COLOR_SCHEME[AA_NUM_MAP[x]] for x in sorted(AA_NUM_MAP.keys())])


def parse_fasta_alignment(fasta_file):
    """Return (headers, sequences) from a FASTA alignment file."""
    headers, sequences = [], []
    with open(fasta_file) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            headers.append(record.id)
            sequences.append(list(record.seq))
    return headers, sequences


def _seqs_to_color_array(sequences):
    """Convert sequence list-of-lists to an int8 numpy array for imshow."""
    # build ASCII → int lookup once; missing chars default to 21 (gap)
    lookup = np.full(128, 21, dtype=np.int8)
    for aa, val in AA_VAL_MAP.items():
        lookup[ord(aa)] = val
    flat = np.frombuffer(
        "".join("".join(seq) for seq in sequences).encode("ascii"),
        dtype=np.uint8,
    )
    return lookup[flat].reshape(len(sequences), len(sequences[0]))


def x_ticks_by_subject(sequences, increment=20):
    """Build (tick_positions, tick_labels) from the first sequence, skipping gaps."""
    seq = sequences[0]
    ticks_loc, tick_labs, counter = [], [], 0
    for i, res in enumerate(seq):
        if res != "-":
            counter += 1
            if counter % increment == 0:
                ticks_loc.append(i)
                tick_labs.append(counter)
    return ticks_loc, tick_labs


def plot_alignment_matrix_matplotlib(fasta_file, outfmt=None, outfile=""):
    """Render a ClustalX-style alignment heatmap; always writes PNG at 300 DPI."""
    headers, sequences = parse_fasta_alignment(fasta_file)

    outfile = outfile or fasta_file.removesuffix(".aln") + ".png"
    # ensure .png extension regardless of outfmt arg
    if not outfile.endswith(".png"):
        outfile = outfile.rsplit(".", 1)[0] + ".png"

    dat_colors = _seqs_to_color_array(sequences)
    num_seqs   = len(sequences)
    aln_len    = len(sequences[0])

    # target ~14 px per cell at screen DPI (96) → 14/96 ≈ 0.146 inches per cell
    _DPI      = 96
    _CELL_W   = 14 / _DPI   # inches per column
    _CELL_H   = 18 / _DPI   # inches per row (taller than wide for readability)
    fig_w = max(aln_len * _CELL_W, 4)
    fig_h = max(num_seqs * _CELL_H, 1)
    fs    = 7

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.imshow(dat_colors, cmap=CLUSTAL_CMAP, aspect="equal", vmin=0, vmax=len(AA_NUM_MAP) - 1)
    ax.set_yticks(range(num_seqs), labels=headers, fontsize=fs)

    ticks_loc, tick_labs = x_ticks_by_subject(sequences, increment=20)
    ax.set_xticks(ticks_loc, tick_labs, fontsize=fs)

    ax.set_title(" ".join(fasta_file.removesuffix(".aln").split(".")))
    ax.set_xlabel("position (aa)")

    for i in range(num_seqs):
        for j in range(aln_len):
            ax.text(j, i, sequences[i][j],
                    ha="center", va="center", color="k", fontsize=fs)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight", dpi=_DPI)
    plt.close(fig)
    return outfile


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-a", "--aln", dest="ALN", required=True,
                        help="FASTA alignment file")
    parser.add_argument("-o", "--output", dest="OUT", default="",
                        help="output PNG path (default: <aln>.png)")
    args = parser.parse_args()

    out = plot_alignment_matrix_matplotlib(args.ALN, outfile=args.OUT)
    print(out)
