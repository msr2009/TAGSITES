"""
utility functions for showing results
"""

import pandas as pd
import json, sys
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from Bio import SeqIO
from numpy import linspace

"""
	input: json file containing analyses and output locations
	output: pandas dataframe containing continuous data (e.g., conservation)
			pandas dataframe containing discrete data (e.g., domains)
"""
def load_data_from_json(json_in, type_dict):
	#initialize output data
	df = pd.DataFrame(columns=["pos"])
	range_cols = ["source", "start", "stop", "description"]
	dat = pd.DataFrame(columns=range_cols)
	alns = []

	j = {}

	#read from a bunch of different ways (dict -> json_str -> json_file)
	if isinstance(json_in, dict):
		j = json_in
	if isinstance(json_in, str):
		try:
			j = json.loads(json_in)
		except ValueError:
			print('input is str, but not json. trying loading from file')
			try:
				with open(json_in, "r") as f:
					j = json.load(f)
			except ValueError:
				raise IOError("Can't read json from this file")

	#loop through json file to get tasks
	for task in j:
		if task in ["scripts", "global"]: continue #this is just a mapping dict
		
		#for each task, determine whether it's a 
		# -- continuous variable (and goes in a dataframe) or
		# -- range (and goes in a dictionary)
		if j[task]["analysis"] in type_dict["CONTINUOUS"]:
			tmp_df = pd.DataFrame()
			#add to dataframe
			#print(f"{task}: df") #add to dict
			read_df = pd.read_csv(j[task]["args"]["output"], sep='\t',
								  comment="#", na_values=[-1000])
			tmp_df["pos"] = range(1, read_df.shape[0])
			divisor = 1.0
			if j[task]["analysis"] == "plddt":
				divisor=100.0
			tmp_df[task] = read_df.iloc[:, 1].div(divisor)
			df = pd.merge(df, tmp_df, how="outer", on="pos")

			#also, if the analysis is "blast", then add path to alns list
			if j[task]["analysis"] == "blast":
				alns.append(j[task]["args"]["output"].replace(".jsd", ".aln"))

		elif j[task]["analysis"] in type_dict["RANGE"]: 
			#print(f"{task}: dict") #add to dict
			#input is [source, start, end, description]
			#output dict should be {source: {desc_start_end: [description, start, end]}}
			tmp_dat = pd.read_csv(j[task]["args"]["output"], sep="\t", names = range_cols)
			dat = pd.concat([dat, tmp_dat], ignore_index=True)
		else:
			raise TypeError(f"Can't determine the analysis type of {task}. Skipping.")
			continue

	return (df, dat, alns)


###FOR PLOTTING RESULTS DATA

def plot_results(aa_df: pd.DataFrame, fasta: str, range_df: pd.DataFrame, title: str = "Results Plot"):
	# make subplots for the various data 
	fig = make_subplots(
		rows=3, cols=1,
		shared_xaxes=True,
		row_heights=[.7, .4, .1],
#		vertical_spacing=0.1
		)

	# in the first subplot we're putting aa data (e.g., jsd, plddt)
	# Add traces based on columns in the dataframe (assuming x, y format)
	for c in range(1, len(aa_df.columns)):
		fig.add_trace(
			go.Scatter(x=aa_df.iloc[:, 0], y=aa_df.iloc[:, c], mode='lines', name=aa_df.columns[c]),
			row=1, col=1
		)
	
	# update the axes for the top subplot
	fig.update_xaxes(title_text="", row=1, col=1)
	fig.update_yaxes(title_text="Score", row=1, col=1)


	#in the second subplot, we put categorical range-based data
	#look for phobius, pfam, and regex in ranges df -- these are the 
	#datasets we want to plot
	height = 2
	preds={"Phobius": 2, "Pfam": 4, "modification": 6}
	for i in range(range_df.shape[0]):
		source = range_df.iloc[i]["source"]
		start = range_df.iloc[i]["start"]
		stop = range_df.iloc[i]["stop"]
		desc = range_df.iloc[i]["description"]
		if source in preds:
			#add rectangles
			y0 = preds[source]-height*.3
			y1 = preds[source]+height*.3
			fig.add_trace(
				go.Scatter(
					x=[start, stop, stop, start, start],
					y=[y0, y0, y1, y1, y0],
					mode="lines",
					fill="toself",
					name=f"{desc}: {start}-{stop}",
					showlegend=False
				),	
				row=2, col=1
			)			
		
	# update the axes for the bottom plot	
	fig.update_yaxes(title_text="", row=2, col=1)
	fig.update_layout(yaxis2 = dict(tickmode="array", 
									tickvals=list(preds.values()), 
									ticktext=list(preds.keys())))
	

	#the last subplot is a single-row heatmap
	#of the subject sequence
#	aln_fig_widget = go.FigureWidget(make_alignment_subplots([fasta])).data[0]
#	fig.add_trace(aln_fig_widget, row=3, col=1)
	fig.add_trace(make_alignment_subplots([fasta]).data[0], row=3, col=1)
		
	# update axes for bottom axis
	fig.update_xaxes(title_text="Position (aa)", row=3, col=1)

	# Customize layout for both subplots
	fig.update_layout(
		title=title,
		showlegend=True,
		xaxis=dict(fixedrange=False),
		yaxis=dict(fixedrange=True),
		clickmode="event",
		hovermode="closest"
		)
	
#	return fig
	fig_widget = go.FigureWidget(fig)
	print(fig_widget)
	return fig_widget

####FOR PLOTTING ALIGNMENTS
# utils/plot_alignment.py

# Define a basic ClustalX-style color scheme for amino acids
CLUSTALX_COLOR_SCHEME = {
	'A': 'green',  'I': 'green',  'L': 'green',  'M': 'green',  'V': 'green',  # Hydrophobic
	'F': 'orange', 'Y': 'orange', 'W': 'orange',  # Aromatic
	'H': 'blue',   'K': 'blue',   'R': 'blue',  # Positively charged
	'D': 'red',    'E': 'red',    # Negatively charged
	'S': 'magenta', 'T': 'magenta',  # Polar, uncharged
	'N': 'cyan',   'Q': 'cyan',    # Polar, uncharged
	'G': 'yellow', 'P': 'yellow', 'C': 'yellow', # Special cases
	'-': 'white', 'X': 'white' # Gap
}

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
def color_seq_for_heatmap(seq_mat, color_mapping=None, nan_color="black"):
    #how many colors
    n_colors = 21 #max aa colors + stop

	#first, exclude any values that aren't in your sequences
	#sometimes there is a bug concerning 

    colors = sorted(list(set(color_mapping.values())))
    if color_mapping:
        n_colors = len(colors)
    
    aa_z_dict = {} ###we'll populate this
    edge_color_tuples = [] ###we'll populate this too
    color_z_dict = {}
    
    #split [0,1] into bins
    bin_edges = linspace(0,1,n_colors+1)
    #loop through edges and populate with colors and aas
    for i in range(len(bin_edges)-1):
        edge_color_tuples.append((bin_edges[i], colors[i]))
        edge_color_tuples.append((bin_edges[i+1], colors[i]))
        
        color_z_dict[colors[i]] = (bin_edges[i] + bin_edges[i+1])/2.0

    seq_mat_colors = []
    for seq in seq_mat:
        seq_colors = []
        for res in seq:
            seq_colors.append(color_z_dict[color_mapping[res]])
        seq_mat_colors.append(seq_colors)
    
    return seq_mat_colors, edge_color_tuples

def plot_alignment_matrix(fasta_file, title="alignment title", colormap=None):
	headers, sequences = parse_fasta_alignment(fasta_file)
	if colormap == None:
		colormap = CLUSTALX_COLOR_SCHEME
	(dat_colors, scale_tuples) = color_seq_for_heatmap(sequences, colormap)

	num_sequences = len(sequences)
	alignment_length = len(sequences[0])

	aln_hm = go.Heatmap(
			z=dat_colors[::-1],
			y=headers[::-1],
			text=sequences[::-1],
			texttemplate='%{text}',
			colorscale=scale_tuples,
			autocolorscale=False,
			showlegend=False)    	

	return aln_hm

def make_alignment_subplots(fasta_file_list, title=""):
	
#	fig = make_subplots(rows=len(fasta_file_list), 
#						cols=1)
	fa_lengths = []
	aln_heatmaps = []
	for i in range(len(fasta_file_list)):
		#read fasta once to get the lengths
		fa_lengths.append(len(parse_fasta_alignment(fasta_file_list[i])[0]))
#		fig.add_trace(plot_alignment_matrix(fasta_file_list[i], title),
#				 row=i+1, col=1)
		aln_heatmaps.append(plot_alignment_matrix(fasta_file_list[i], title))

	fig = make_subplots(rows=len(fasta_file_list), 
						cols=1, 
						row_heights=fa_lengths)
	
	for hm in range(len(aln_heatmaps)):
		fig.add_trace(aln_heatmaps[hm], row=hm+1, col=1)

	fig.update_traces(showscale=False)
	fig.update_yaxes({"tickfont_size":8})
	fig.update_layout(
		title=title,
		showlegend=False,
		xaxis=dict(fixedrange=False),
		yaxis=dict(fixedrange=True),
		#width=1000, height=500
	)

	return fig

