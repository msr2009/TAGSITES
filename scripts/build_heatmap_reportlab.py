"""
build_heatmap_reportlab.py — ClustalX-style alignment heatmap renderer using ReportLab.
Produces a true vector PDF; faster and smaller than the Pillow PNG approach.
Matt Rich, 2026
"""

import os
from pathlib import Path
from Bio import SeqIO
from reportlab.pdfgen import canvas


# ClustalX color scheme as 0–1 RGB tuples for ReportLab
CLUSTALX_RGB = {
    'A': (0, 0.71, 0),    'I': (0, 0.71, 0),   'L': (0, 0.71, 0),
    'M': (0, 0.71, 0),    'V': (0, 0.71, 0),
    'F': (1, 0.65, 0),    'Y': (1, 0.65, 0),   'W': (1, 0.65, 0),
    'H': (0.12, 0.39, 1), 'K': (0.12, 0.39, 1),'R': (0.12, 0.39, 1),
    'D': (0.86, 0, 0),    'E': (0.86, 0, 0),
    'S': (0.78, 0, 0.78), 'T': (0.78, 0, 0.78),
    'N': (0, 0.78, 0.78), 'Q': (0, 0.78, 0.78),
    'G': (0.82, 0.82, 0), 'P': (0.82, 0.82, 0),'C': (0.82, 0.82, 0),
    '-': (0.94, 0.94, 0.94), 'X': (0.94, 0.94, 0.94),
}
DEFAULT_RGB  = (0.78, 0.78, 0.78)

# layout constants (points — same numeric values as the pixel layout in build_heatmap_pillow)
CELL_W      = 20
CELL_H      = 26
LABEL_PAD   = 10
TICK_EVERY  = 20
TITLE_H     = 36
TICK_BAND_H = 28
MARGIN_BOT  = 32
FONT_SIZE   = 11   # cell letters — sized to fit within CELL_W=20
LABEL_FONT  = 9
TICK_FONT   = 11
LABEL_MAXLEN = 30


def parse_fasta_alignment(fasta_file):
    """Return (headers, sequences) from a FASTA alignment file."""
    headers, sequences = [], []
    with open(fasta_file) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            headers.append(record.id)
            sequences.append(str(record.seq))
    return headers, sequences



def plot_alignment_reportlab(fasta_file, outfile=""):
    """Render a ClustalX-style alignment heatmap as a vector PDF."""
    headers, sequences = parse_fasta_alignment(fasta_file)
    outfile = outfile or str(fasta_file).removesuffix(".aln") + ".pdf"

    headers = [h[:LABEL_MAXLEN] for h in headers]
    num_seqs = len(sequences)
    aln_len  = len(sequences[0]) if sequences else 0

    # estimate label column width from character count (Helvetica ~0.55× em)
    label_w  = max(len(h) for h in headers) * LABEL_FONT * 0.55 + LABEL_PAD
    grid_w   = aln_len * CELL_W
    grid_h   = num_seqs * CELL_H
    page_w   = label_w + grid_w
    page_h   = TITLE_H + TICK_BAND_H + grid_h + MARGIN_BOT
    grid_top = TITLE_H + TICK_BAND_H   # distance from top of page to top of grid

    # ReportLab origin is bottom-left; fy() converts from top-down coords
    def fy(y): return page_h - y

    c = canvas.Canvas(outfile, pagesize=(page_w, page_h))

    # title
    title = Path(str(fasta_file).removesuffix(".aln")).name.replace(".", " ")
    c.setFont("Helvetica-Bold", FONT_SIZE + 2)
    c.setFillColorRGB(0, 0, 0)
    c.drawCentredString(page_w / 2, fy(TITLE_H - 8), title)

    # tick positions: count non-gap residues in reference (first) sequence
    ref_seq = sequences[0] if sequences else ""
    tick_positions, pos_counter = [], 0
    for j, aa in enumerate(ref_seq):
        if aa != "-":
            pos_counter += 1
            if pos_counter % TICK_EVERY == 0:
                tick_positions.append((j, pos_counter))

    def draw_ticks(y_grid_edge, above):
        """Draw tick marks and position labels at the top or bottom of the grid."""
        c.setFont("Helvetica", TICK_FONT)
        tick_len   = 5 if above else -5
        label_dy   = tick_len + (2 if above else -TICK_FONT - 2)
        for col_j, pos in tick_positions:
            x = label_w + col_j * CELL_W + CELL_W / 2
            y = fy(y_grid_edge)
            c.setStrokeColorRGB(0.3, 0.3, 0.3)
            c.line(x, y, x, y + tick_len)
            c.setFillColorRGB(0.24, 0.24, 0.24)
            c.drawCentredString(x, y + label_dy, str(pos))

    draw_ticks(grid_top, above=True)

    # cells
    for i, seq in enumerate(sequences):
        y0 = grid_top + i * CELL_H

        # row label
        c.setFont("Helvetica", LABEL_FONT)
        c.setFillColorRGB(0, 0, 0)
        # vertical centre of the row
        c.drawRightString(label_w - LABEL_PAD, fy(y0 + CELL_H * 0.35), headers[i])

        c.setFont("Helvetica-Bold", FONT_SIZE)
        for j, aa in enumerate(seq):
            x0  = label_w + j * CELL_W
            bg  = CLUSTALX_RGB.get(aa.upper(), DEFAULT_RGB)
            c.setFillColorRGB(*bg)
            c.rect(x0, fy(y0 + CELL_H), CELL_W, CELL_H, fill=1, stroke=0)
            c.setFillColorRGB(0, 0, 0)
            # baseline offset ≈ 30% of cell height keeps text visually centred
            c.drawCentredString(x0 + CELL_W / 2, fy(y0 + CELL_H) + CELL_H * 0.28, aa)

    draw_ticks(grid_top + grid_h, above=False)

    c.save()
    return outfile


if __name__ == "__main__":
    from argparse import ArgumentParser
    import time

    parser = ArgumentParser(description="ReportLab vector PDF alignment heatmap renderer")
    parser.add_argument("-a", "--aln", dest="ALN", required=True, help="FASTA alignment file")
    parser.add_argument("-o", "--output", dest="OUT", default="", help="output PDF path")
    args = parser.parse_args()

    t0  = time.time()
    out = plot_alignment_reportlab(args.ALN, outfile=args.OUT)
    print(f"Written: {out}  ({time.time() - t0:.2f}s)")
