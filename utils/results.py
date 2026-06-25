"""utility functions for loading and visualizing analysis results"""

import os
import pandas as pd
import json
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from Bio import SeqIO
from numpy import linspace


def load_data_from_json(json_in, type_dict):
    """Load per-task result TSVs referenced in the run JSON.

    Returns (aa_df, range_df, aln_file_list). Tasks whose output files don't
    exist yet are skipped with a warning rather than crashing.
    """
    df = pd.DataFrame(columns=["pos"])
    range_cols = ["source", "start", "stop", "description"]
    dat = pd.DataFrame(columns=range_cols)
    alns = []

    # accept dict, JSON string, or file path
    j = {}
    if isinstance(json_in, dict):
        j = json_in
    elif isinstance(json_in, str):
        try:
            j = json.loads(json_in)
        except ValueError:
            try:
                with open(json_in, "r") as f:
                    j = json.load(f)
            except (ValueError, FileNotFoundError):
                raise IOError(f"Cannot read JSON from: {json_in}")

    for task in j:
        if task in ["scripts", "global"]:
            continue

        analysis = j[task]["analysis"]
        output_path = j[task]["args"].get("output", "")

        # skip tasks whose output files haven't been written yet
        if not output_path or not os.path.exists(output_path):
            print(f"Output not found for task '{task}' ({analysis}): {output_path}")
            continue

        if analysis in type_dict["CONTINUOUS"]:
            try:
                read_df = pd.read_csv(output_path, sep="\t",
                                      header=None, comment="#", na_values=[-1000])
                tmp_df = pd.DataFrame()
                tmp_df["pos"] = [x + 1 for x in range(read_df.shape[0] + 1)]
                divisor = 100.0 if analysis == "plddt" else 1.0
                tmp_df[task] = read_df.iloc[:, 1].div(divisor)
                df = pd.merge(df, tmp_df, how="outer", on="pos")
                # blast task also has an alignment file to display
                if analysis == "blast":
                    alns.append(output_path.replace(".jsd", ".aln"))
            except Exception as e:
                print(f"Error loading continuous data for task '{task}': {e}")
                continue

        elif analysis in type_dict["RANGE"]:
            try:
                tmp_dat = pd.read_csv(output_path, sep="\t", names=range_cols)
                dat = pd.concat([dat, tmp_dat], ignore_index=True)
            except Exception as e:
                print(f"Error loading range data for task '{task}': {e}")
                continue

        else:
            print(f"Unknown analysis type '{analysis}' for task '{task}'. Skipping.")
            continue

    return df, dat, alns


###FOR PLOTTING RESULTS DATA

def plot_results(aa_df, fasta, range_df, title="Results Plot"):
    """Build the main three-row Plotly results figure."""
    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        row_heights=[.7, .4, .1],
    )

    # row 1: continuous per-position scores (JSD, pLDDT, etc.)
    for c in range(1, len(aa_df.columns)):
        fig.add_trace(
            go.Scatter(x=aa_df.iloc[:, 0], y=aa_df.iloc[:, c],
                       mode="lines", name=aa_df.columns[c]),
            row=1, col=1
        )
    fig.update_xaxes(title_text="", row=1, col=1)
    fig.update_yaxes(title_text="Score", row=1, col=1)

    # row 2: categorical range data (domains, PTM sites)
    # only Phobius, Pfam, and modification sources are plotted
    height = 2
    preds = {"Phobius": 2, "Pfam": 4, "modification": 6}
    for i in range(range_df.shape[0]):
        source = range_df.iloc[i]["source"]
        start = range_df.iloc[i]["start"]
        stop = range_df.iloc[i]["stop"]
        desc = range_df.iloc[i]["description"]
        if source in preds:
            y0 = preds[source] - height * .3
            y1 = preds[source] + height * .3
            fig.add_trace(
                go.Scatter(
                    x=[start, stop, stop, start, start],
                    y=[y0, y0, y1, y1, y0],
                    mode="lines", fill="toself",
                    name=f"{desc}: {start}-{stop}",
                    showlegend=False
                ),
                row=2, col=1
            )
    fig.update_yaxes(title_text="", row=2, col=1)
    fig.update_layout(yaxis2=dict(
        tickmode="array",
        tickvals=list(preds.values()),
        ticktext=list(preds.keys())
    ))

    # row 3: single-row subject sequence heatmap
    fig.add_trace(make_alignment_subplots([fasta]).data[0], row=3, col=1)
    fig.update_xaxes(title_text="Position (aa)", row=3, col=1)

    fig.update_layout(
        title=title,
        showlegend=True,
        xaxis=dict(fixedrange=False),
        yaxis=dict(fixedrange=True),
        clickmode="event",
        hovermode="closest",
    )
    return go.FigureWidget(fig)


####FOR PLOTTING ALIGNMENTS

# ClustalX-style amino acid color scheme
CLUSTALX_COLOR_SCHEME = {
    "A": "green",  "I": "green",  "L": "green",  "M": "green",  "V": "green",
    "F": "orange", "Y": "orange", "W": "orange",
    "H": "blue",   "K": "blue",   "R": "blue",
    "D": "red",    "E": "red",
    "S": "magenta", "T": "magenta",
    "N": "cyan",   "Q": "cyan",
    "G": "yellow", "P": "yellow", "C": "yellow",
    "-": "white",  "X": "white",
}


def parse_fasta_alignment(fasta_file):
    """Parse a FASTA alignment file into (headers, sequences) lists."""
    sequences = []
    headers = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(list(record.seq))
            headers.append(record.id)
    return headers, sequences


def color_seq_for_heatmap(seq_mat, color_mapping=None, nan_color="black"):
    """Convert a sequence matrix to numeric z-values for Plotly heatmap rendering."""
    colors = sorted(set(color_mapping.values()))
    n_colors = len(colors)

    edge_color_tuples = []
    color_z_dict = {}
    bin_edges = linspace(0, 1, n_colors + 1)
    for i in range(len(bin_edges) - 1):
        edge_color_tuples.append((bin_edges[i], colors[i]))
        edge_color_tuples.append((bin_edges[i + 1], colors[i]))
        color_z_dict[colors[i]] = (bin_edges[i] + bin_edges[i + 1]) / 2.0

    seq_mat_colors = [
        [color_z_dict[color_mapping[res]] for res in seq]
        for seq in seq_mat
    ]
    return seq_mat_colors, edge_color_tuples


def plot_alignment_matrix(fasta_file, title="alignment title", colormap=None):
    """Build a Plotly Heatmap trace for one alignment FASTA file."""
    headers, sequences = parse_fasta_alignment(fasta_file)
    if colormap is None:
        colormap = CLUSTALX_COLOR_SCHEME
    dat_colors, scale_tuples = color_seq_for_heatmap(sequences, colormap)

    return go.Heatmap(
        z=dat_colors[::-1],
        y=headers[::-1],
        text=sequences[::-1],
        texttemplate="%{text}",
        colorscale=scale_tuples,
        autocolorscale=False,
        showlegend=False,
        hoverinfo="skip",
        x0=1, dx=1,
    )


def make_alignment_subplots(fasta_file_list, title=""):
    """Combine multiple alignment FASTA files into a stacked Plotly figure."""
    fa_lengths = []
    aln_heatmaps = []
    for fasta_file in fasta_file_list:
        fa_lengths.append(len(parse_fasta_alignment(fasta_file)[0]))
        aln_heatmaps.append(plot_alignment_matrix(fasta_file, title))

    fig = make_subplots(
        rows=len(fasta_file_list), cols=1,
        row_heights=[5 * x for x in fa_lengths]
    )
    for i, hm in enumerate(aln_heatmaps):
        fig.add_trace(hm, row=i + 1, col=1)

    fig.update_traces(showscale=False)
    fig.update_xaxes({"fixedrange": False})
    fig.update_yaxes({"tickfont_size": 8, "fixedrange": True})
    fig.update_layout(title=title, showlegend=False)
    return fig
