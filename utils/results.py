"""utility functions for loading and visualizing analysis results"""

import os
import re
import pandas as pd
import json
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from Bio import SeqIO
from numpy import linspace

from config import ANALYSIS_COLORS, DOMAIN_SOURCE_COLORS


# keyword → short label mapping for Phobius region descriptions from EBI InterProScan5
_PHOBIUS_KEYWORDS = [
    ("cytoplasm",              "Cytoplasmic"),
    ("extracellular region",   "Extracellular"),
    ("embedded in the membrane", "Transmembrane"),
    ("signal peptide",         "Signal peptide"),
    ("signal sequence",        "Signal peptide"),
]


def _translate_phobius_desc(desc):
    """Map verbose Phobius region description to a short topology label."""
    d = desc.lower()
    for keyword, label in _PHOBIUS_KEYWORDS:
        if keyword in d:
            return label
    return desc


def _merge_phobius_intervals(phobius_df):
    """Merge overlapping/contiguous Phobius rows that share the same translated label.

    Handles cases like signal peptides, which InterProScan5 emits as several
    overlapping sub-region rows (N-terminal, hydrophobic, C-terminal, whole span).
    Returns a new DataFrame with the merged rows.
    """
    merged_rows = []
    for desc, group in phobius_df.groupby("description", sort=False):
        intervals = sorted(zip(group["start"].astype(int), group["stop"].astype(int)))
        cur_start, cur_stop = intervals[0]
        for start, stop in intervals[1:]:
            # contiguous or overlapping: extend current interval
            if start <= cur_stop + 1:
                cur_stop = max(cur_stop, stop)
            else:
                merged_rows.append({"source": "Phobius", "start": cur_start,
                                    "stop": cur_stop, "description": desc})
                cur_start, cur_stop = start, stop
        merged_rows.append({"source": "Phobius", "start": cur_start,
                            "stop": cur_stop, "description": desc})
    return pd.DataFrame(merged_rows, columns=["source", "start", "stop", "description"])


def load_data_from_json(json_in, type_dict):
    """Load per-task result TSVs referenced in the run JSON.

    Returns (aa_df, range_df, aln_meta_list).
    aln_meta_list: list of (aln_path, task_name, params_dict) for blast tasks.
    Tasks whose output files don't exist yet are skipped with a warning.
    """
    df = pd.DataFrame(columns=["pos"])
    range_cols = ["source", "start", "stop", "description"]
    dat = pd.DataFrame(columns=range_cols)
    alns = []  # list of (aln_path, task_name, params_dict)

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
                # positions are 1-based, one row per residue in the output file
                tmp_df["pos"] = [x + 1 for x in range(read_df.shape[0])]
                divisor = 100.0 if analysis == "plddt" else 1.0
                tmp_df[task] = read_df.iloc[:, 1].div(divisor)
                df = pd.merge(df, tmp_df, how="outer", on="pos")
                # pLDDT tasks also produce a companion SASA file
                if analysis == "plddt":
                    sasa_path = output_path.replace("plddt.txt", "sasa.txt")
                    if os.path.exists(sasa_path):
                        try:
                            sasa_df = pd.read_csv(sasa_path, sep="\t", header=None,
                                                  comment="#", na_values=[-1000])
                            stmp = pd.DataFrame()
                            stmp["pos"] = sasa_df.iloc[:, 0].astype(int)
                            sasa_vals = sasa_df.iloc[:, 1].rolling(window=3, center=True, min_periods=1).mean().clip(upper=1.0)
                            stmp[task + "_sasa"] = sasa_vals
                            df = pd.merge(df, stmp, how="outer", on="pos")
                        except Exception as e:
                            print(f"Error loading SASA for task '{task}': {e}")
                    else:
                        print(f"SASA file not found for task '{task}': {sasa_path}")
                # blast tasks also produce an alignment file
                if analysis == "blast":
                    aln_path = output_path.replace(".jsd", ".aln")
                    # collect user-visible params for the alignment pane
                    params = {
                        k: v for k, v in j[task]["args"].items()
                        if k not in ("output", "email", "working_dir", "run_name",
                                     "input_file", "pdb", "scripts_folder",
                                     "genomic_file", "fasta")
                        and v not in ("", None)
                    }
                    alns.append((aln_path, task, params))
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

    # translate verbose Phobius descriptions and merge overlapping intervals per label
    if not dat.empty:
        mask = dat["source"] == "Phobius"
        dat.loc[mask, "description"] = dat.loc[mask, "description"].map(_translate_phobius_desc)
        phobius_merged = _merge_phobius_intervals(dat[mask])
        dat = pd.concat([dat[~mask], phobius_merged], ignore_index=True)

    return df, dat, alns


def load_run_metadata(json_in):
    """Extract run-level metadata needed by the results UI (sequence, PDB path, task params).

    Returns a dict with keys:
      query_seq  — the full protein sequence as a string (or "" on error)
      pdb_path   — absolute path to a PDB file if one exists, else ""
      seq_len    — integer sequence length
    """
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
                return {"query_seq": "", "pdb_path": "", "seq_len": 0}

    g = j.get("global", {})
    query_seq = ""
    fasta_path = g.get("input_file", "")
    if fasta_path and os.path.exists(fasta_path):
        try:
            with open(fasta_path) as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    query_seq = str(rec.seq)
                    break
        except Exception:
            pass

    # PDB: first check explicit global pdb field, then look for AFDB-downloaded file
    pdb_path = g.get("pdb", "")
    if not pdb_path or not os.path.exists(pdb_path):
        working_dir = g.get("working_dir", "")
        run_name = g.get("run_name", "")
        af_pdb = os.path.join(working_dir, f"{run_name}.AF.pdb")
        pdb_path = af_pdb if os.path.exists(af_pdb) else ""

    return {
        "query_seq": query_seq,
        "pdb_path": pdb_path,
        "seq_len": len(query_seq),
    }


###FOR PLOTTING RESULTS DATA

def assign_task_colors(aa_df):
    """Assign a stable hex color to each task column, cycling through the palette per analysis type."""
    colors = {}
    counters = {}
    for col_name in aa_df.columns[1:]:
        analysis = _guess_analysis_type(col_name)
        palette = ANALYSIS_COLORS.get(analysis, ["#888888", "#aaaaaa", "#666666"])
        idx = counters.get(analysis, 0)
        colors[col_name] = palette[idx % len(palette)]
        counters[analysis] = idx + 1
    return colors


def _guess_analysis_type(task_name):
    """Return the analysis type keyword (blast/plddt/scores/…) matching a task name."""
    t = task_name.lower()
    # sasa is co-produced by the pLDDT task; group it there for color palette purposes
    if t.endswith("_sasa"):
        return "plddt"
    for analysis in ANALYSIS_COLORS:
        if analysis in t:
            return analysis
    return "unknown"


def plot_results(aa_df, range_df, title="Results Plot"):
    """Build the two-row Plotly results figure (continuous scores + range annotations).

    Colors are stable: driven by ANALYSIS_COLORS and DOMAIN_SOURCE_COLORS from config.
    The row-3 sequence heatmap has been removed; the sequence strip is a separate
    interactive JS component.
    """
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.65, 0.35],
    )

    task_colors = assign_task_colors(aa_df)

    # row 1: continuous per-position scores — one line per task column
    for col_name in aa_df.columns[1:]:
        color = task_colors.get(col_name, "#888888")
        fig.add_trace(
            go.Scatter(
                x=aa_df["pos"], y=aa_df[col_name],
                mode="lines", name=col_name,
                line=dict(color=color),
            ),
            row=1, col=1
        )
    fig.update_yaxes(title_text="Score", row=1, col=1, fixedrange=True)

    # row 2: range/categorical data — filled rectangles per annotation
    height = 2
    y_positions = {"Phobius": 2, "Pfam": 4, "modification": 6}
    for _, row in range_df.iterrows():
        source = row["source"]
        if source not in y_positions:
            continue
        start, stop, desc = row["start"], row["stop"], row["description"]
        y0 = y_positions[source] - height * 0.3
        y1 = y_positions[source] + height * 0.3
        color = DOMAIN_SOURCE_COLORS.get(source, "#888888")
        fig.add_trace(
            go.Scatter(
                x=[start, stop, stop, start, start],
                y=[y0, y0, y1, y1, y0],
                mode="lines", fill="toself",
                fillcolor=color, line=dict(color=color),
                name=f"{desc}: {start}-{stop}",
                showlegend=False,
                hovertemplate=f"{source}: {desc} ({start}-{stop})<extra></extra>",
            ),
            row=2, col=1
        )
    fig.update_layout(yaxis2=dict(
        tickmode="array",
        tickvals=list(y_positions.values()),
        ticktext=list(y_positions.keys()),
        fixedrange=True,
    ))
    fig.update_xaxes(title_text="Position (aa)", row=2, col=1)

    fig.update_layout(
        title=title,
        showlegend=True,
        xaxis=dict(fixedrange=False),
        hovermode="closest",
        height=500,
        margin=dict(l=60, r=20, t=40, b=40),
    )
    return go.Figure(fig)


def _task_color(task_name):
    """Return the first palette color for a task's analysis type (used as gradient endpoint)."""
    analysis = _guess_analysis_type(task_name)
    palette = ANALYSIS_COLORS.get(analysis, ["#888888"])
    return palette[0]


def residue_colors_for_track(aa_df, task_name, hex_color):
    """Compute per-residue hex colors for a continuous track.

    pLDDT tracks use the AlphaFold 4-band categorical scheme; SASA and all
    other tracks use a white → hex_color gradient.
    Use residue_colors_gradient() directly to force gradient mode on pLDDT.
    """
    if task_name not in aa_df.columns:
        return []
    if _guess_analysis_type(task_name) == "plddt" and not task_name.endswith("_sasa"):
        return _residue_colors_plddt(aa_df, task_name)
    return residue_colors_gradient(aa_df, task_name, hex_color)


def residue_colors_gradient(aa_df, task_name, hex_color=None):
    """Viridis gradient for any continuous track (hex_color ignored, kept for API compat)."""
    # viridis sampled at 9 stops (low → high)
    _VIRIDIS = [
        (68,  1,  84), (72,  40, 120), (62,  74, 137), (49, 104, 142),
        (38, 130, 142), (31, 158, 137), (53, 183, 121), (110, 206,  88),
        (181, 222,  43), (253, 231,  37),
    ]

    if task_name not in aa_df.columns:
        return []
    series = aa_df.set_index("pos")[task_name]
    # normalize to fixed [0, 1] so color position is absolute, not relative to data range
    normed = series

    def _viridis_hex(t):
        # linearly interpolate between the two nearest stops
        t = max(0.0, min(1.0, t))
        idx = t * (len(_VIRIDIS) - 1)
        lo, hi = int(idx), min(int(idx) + 1, len(_VIRIDIS) - 1)
        frac = idx - lo
        r = int(_VIRIDIS[lo][0] + (_VIRIDIS[hi][0] - _VIRIDIS[lo][0]) * frac)
        g = int(_VIRIDIS[lo][1] + (_VIRIDIS[hi][1] - _VIRIDIS[lo][1]) * frac)
        b = int(_VIRIDIS[lo][2] + (_VIRIDIS[hi][2] - _VIRIDIS[lo][2]) * frac)
        return f"#{r:02x}{g:02x}{b:02x}"

    colors = []
    for pos in sorted(series.index):
        v = normed.get(pos, float("nan"))
        colors.append("#aaaaaa" if pd.isna(v) else _viridis_hex(v))
    return colors


def _residue_colors_plddt(aa_df, task_name):
    """Apply the standard AlphaFold pLDDT 4-band color scheme.

    Values are stored 0–1 (divided by 100 at load time):
      ≥ 0.90  very high confidence → dark blue  #0053D6
      ≥ 0.70  confident            → light blue #65CBF3
      ≥ 0.50  low confidence       → yellow     #FFDB13
      <  0.50 very low             → orange     #FF7D45
    """
    series = aa_df.set_index("pos")[task_name]
    colors = []
    for pos in sorted(series.index):
        v = series.get(pos, float("nan"))
        if pd.isna(v):
            colors.append("#ffffff")
        elif v >= 0.90:
            colors.append("#0053D6")
        elif v >= 0.70:
            colors.append("#65CBF3")
        elif v >= 0.50:
            colors.append("#FFDB13")
        else:
            colors.append("#FF7D45")
    return colors


def residue_colors_for_domains(range_df, seq_len):
    """Compute per-residue hex colors for domain annotations.

    Each residue gets the color of the first domain covering it; overlapping
    domains are colored by their source priority. Unannotated residues get a
    light gray.
    Returns a list of seq_len hex strings (index 0 = residue 1).
    """
    colors = ["#e8e8e8"] * seq_len
    # apply by source priority: modifications on top
    priority = ["Pfam", "Phobius", "modification"]
    for source in priority:
        subset = range_df[range_df["source"] == source]
        color = DOMAIN_SOURCE_COLORS.get(source, "#888888")
        for _, row in subset.iterrows():
            start, stop = int(row["start"]), int(row["stop"])
            for pos in range(start - 1, min(stop, seq_len)):
                colors[pos] = color
    return colors


# Qualitative palette for unique-per-annotation coloring — saturated, maximally distinct
_ANNOTATION_PALETTE = [
    "#e41a1c", "#1565c0", "#2e7d32", "#e65100", "#6a1b9a",
    "#00838f", "#ad1457", "#9e9d24", "#283593", "#00695c",
    "#b71c1c", "#0277bd", "#33691e", "#bf360c", "#4a148c",
    "#006064", "#880e4f", "#f9a825", "#1a237e", "#004d40",
]

# Sources shown in the 2D feature panel
_FEATURE_SOURCES = {"Pfam", "Phobius", "modification"}

# Sources painted in the __domains__ structure colorscheme (Phobius has its own scheme)
_DOMAIN_STRUCT_SOURCES = {"Pfam", "modification"}

# paint priority: last painted wins; modifications end up on top
_FEATURE_PAINT_ORDER = ["Pfam", "Phobius", "modification"]


def _annotation_color_map(range_df):
    """Assign unique palette colors to each visible (source, description) pair.

    Only considers rows whose source is in _FEATURE_SOURCES so non-displayed
    InterPro annotations don't consume palette slots.
    Returns an ordered dict {(source, desc): hex_color}.
    """
    color_map = {}
    if range_df is None:
        return color_map
    for _, row in range_df.iterrows():
        if row["source"] not in _FEATURE_SOURCES:
            continue
        key = (str(row["source"]), str(row["description"]))
        if key not in color_map:
            color_map[key] = _ANNOTATION_PALETTE[len(color_map) % len(_ANNOTATION_PALETTE)]
    return color_map


def residue_colors_for_annotations(range_df, seq_len):
    """Unique color per Pfam/modification annotation for the __domains__ structure scheme.

    Phobius is intentionally excluded here — it has its own __phobius__ colorscheme.
    Returns (per_residue_colors, legend_items).
    """
    color_map = _annotation_color_map(range_df)
    colors = ["#e8e8e8"] * seq_len

    # paint only domain-struct sources; modifications land on top
    for source in _FEATURE_PAINT_ORDER:
        if source not in _DOMAIN_STRUCT_SOURCES:
            continue
        for _, row in range_df[range_df["source"] == source].iterrows():
            key = (str(row["source"]), str(row["description"]))
            color = color_map.get(key)
            if color is None:
                continue
            start, stop = int(row["start"]), int(row["stop"])
            for pos in range(start - 1, min(stop, seq_len)):
                colors[pos] = color

    legend_items = [
        {"color": color, "label": f"{desc} ({src})"}
        for (src, desc), color in color_map.items()
        if src in _DOMAIN_STRUCT_SOURCES
    ]
    return colors, legend_items


def residue_colors_for_phobius(range_df, seq_len):
    """Color structure by Phobius topology region using the shared annotation palette.

    Returns (per_residue_colors, legend_items).
    """
    color_map = _annotation_color_map(range_df)
    colors = ["#e8e8e8"] * seq_len
    for _, row in range_df[range_df["source"] == "Phobius"].iterrows():
        key = ("Phobius", str(row["description"]))
        color = color_map.get(key)
        if color is None:
            continue
        start, stop = int(row["start"]), int(row["stop"])
        for pos in range(start - 1, min(stop, seq_len)):
            colors[pos] = color
    legend_items = [
        {"color": color, "label": desc}
        for (src, desc), color in color_map.items()
        if src == "Phobius"
    ]
    return colors, legend_items


# 5-stop jet colormap matching the CSS legend gradient (blue→cyan→green→yellow→red)
_JET_STOPS = [
    (0.00, (0,   0,   255)),
    (0.25, (0,   255, 255)),
    (0.50, (0,   255, 0  )),
    (0.75, (255, 255, 0  )),
    (1.00, (255, 0,   0  )),
]


def _jet_hex(t):
    """Return hex color for t ∈ [0,1] using the 5-stop jet gradient."""
    t = max(0.0, min(1.0, t))
    for i in range(len(_JET_STOPS) - 1):
        t0, c0 = _JET_STOPS[i]
        t1, c1 = _JET_STOPS[i + 1]
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            r, g, b = (int(c0[j] + f * (c1[j] - c0[j])) for j in range(3))
            return f"#{r:02x}{g:02x}{b:02x}"
    return "#ff0000"


def residue_colors_jet(seq_len):
    """N→C jet gradient per residue, matching the structure legend CSS gradient."""
    if seq_len <= 1:
        return ["#0000ff"] * max(seq_len, 0)
    return [_jet_hex(i / (seq_len - 1)) for i in range(seq_len)]


def build_plot_payload(aa_df, range_df, title="Results"):
    """Serialize plot data for the native canvas renderer (tagsites_set_plot message).

    Returns a dict with lineTracks, rangeFeatures, and title, ready for JSON serialization.
    NaN values in aa_df are converted to None (→ null in JSON).
    """
    y_positions = {"Phobius": 2, "Pfam": 4, "modification": 6}

    line_tracks = []
    data_max = -float("inf")
    if aa_df is not None:
        colors = assign_task_colors(aa_df)
        for col in aa_df.columns[1:]:  # skip 'pos'
            vals = [None if pd.isna(v) else float(v) for v in aa_df[col]]
            line_tracks.append({"name": col, "color": colors.get(col, "#888888"), "values": vals})
            col_max = aa_df[col].dropna().max()
            if pd.notna(col_max):
                data_max = max(data_max, float(col_max))
    y_max = max(1.1, data_max) if data_max > -float("inf") else 1.1

    range_features = []
    if range_df is not None and not range_df.empty:
        color_map = _annotation_color_map(range_df)
        for _, row in range_df.iterrows():
            src = row["source"]
            if src not in y_positions:
                continue
            key = (str(src), str(row["description"]))
            range_features.append({
                "source": src,
                "start": int(row["start"]),
                "stop": int(row["stop"]),
                "desc": str(row["description"]),
                "color": color_map.get(key, "#888888"),
                "yRow": y_positions[src],
            })

    return {"title": title, "lineTracks": line_tracks, "rangeFeatures": range_features,
            "yMax": y_max}


def _hex_to_rgb(hex_color):
    """Convert '#rrggbb' to (r, g, b) integers."""
    h = hex_color.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


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
