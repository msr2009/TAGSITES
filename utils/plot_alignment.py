# utils/plot_alignment.py

import plotly.graph_objects as go
from Bio import SeqIO

# Define a basic ClustalX-style color scheme for amino acids
CLUSTALX_COLOR_SCHEME = {
        'A': 'green',  'I': 'green',  'L': 'green',  'M': 'green',  'V': 'green',  # Hydrophobic
        'F': 'orange', 'Y': 'orange', 'W': 'orange',  # Aromatic
        'H': 'blue',   'K': 'blue',   'R': 'blue',  # Positively charged
        'D': 'red',    'E': 'red',    # Negatively charged
        'S': 'magenta', 'T': 'magenta',  # Polar, uncharged
        'N': 'cyan',   'Q': 'cyan',    # Polar, uncharged
        'G': 'yellow', 'P': 'yellow',  # Special cases
        '-': 'white'  # Gap
}

def parse_fasta_alignment(fasta_file):
        sequences = []
        headers = []
    
        with open(fasta_file, 'r') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                        sequences.append(str(record.seq))
                        headers.append(record.id)
    
        return headers, sequences

def create_color_matrix(sequences):
        color_matrix = []
    
        for seq in sequences:
                row_colors = [CLUSTALX_COLOR_SCHEME.get(aa, 'white') for aa in seq]
                color_matrix.append(row_colors)
    
        return color_matrix

def plot_alignment_matrix(headers, sequences, color_matrix):
        num_sequences = len(sequences)
        alignment_length = len(sequences[0])
    
        fig = go.Figure()

        # Add subject sequence (first row)
        for j, color in enumerate(color_matrix[0]):
                fig.add_trace(go.Scatter(
                        x=[j], y=[0],
                        mode="markers",
                        marker=dict(color=color, size=20),
                        hoverinfo="text",
                        text=f"Subject<br>Position: {j+1}<br>{sequences[0][j]}",
                        showlegend=False,
                        customdata=[j+1],
                        name='subject'
                ))

        fig.update_yaxes(range=[-0.5, num_sequences + 0.5])

        # Add non-interactive rows for other sequences
        for i in range(1, num_sequences):
                for j, color in enumerate(color_matrix[i]):
                        fig.add_trace(go.Scatter(
                                x=[j], y=[i + 1],
                                mode="markers",
                                marker=dict(color=color, size=20),
                                hoverinfo="none",
                                showlegend=False
                        ))

        fig.update_layout(
                title="Protein Sequence Alignment",
                xaxis=dict(
                        tickmode="array",
                        tickvals=list(range(alignment_length)),
                        ticktext=[str(i+1) for i in range(alignment_length)],
                        title="Position",
                        showgrid=False
                ),
                yaxis=dict(
                        tickmode="array",
                        tickvals=list(range(num_sequences + 1)),
                        ticktext=["Subject"] + headers[1:],
                        title="Sequences",
                        showgrid=False
                ),
                plot_bgcolor='white',
                xaxis_showgrid=False,
                yaxis_showgrid=False,
                clickmode='event+select'
        )

        return fig
