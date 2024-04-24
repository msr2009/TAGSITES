#from site_selection_util import *

import numpy as np
from scipy.signal import find_peaks
import re, string, sys
from Bio import SeqIO, PDB 

def smooth_iterative(data, windowsize):
    """
    function to (potentially) iteratively smooth data. Calls numpy convolve
    for smoothing.
    
    Parameters:
        data (list): data to be smoothed
        windowsize (int or list): if windowsize is a list, smoothing is 
				performed iteratively using window sizes in list.
        
    Returns:
        smoothed data (list)
    """
    
    # If windowsize is an int, change it into a list
    if not isinstance(windowsize, list):
        windowsize = list(windowsize)
     
    # Perform iterative smoothing
    tmp_smoothed = data
    for w in windowsize:
        window = np.ones(w) / w
        tmp_smoothed = np.convolve(tmp_smoothed, window, mode='valid')
        
    return tmp_smoothed 
        
def split_minima(minima, min_dist):
    """
    Recursively splits a set of local minima into subsets that are separated
    by at least `min_dist` distance between adjacent elements.
    
    Parameters:
        minima (list): List of indices of local minima.
        min_dist (int): Minimum distance allowed between adjacent elements.
    
    Returns:
        List of mean indices for each subset.
    """
    
    # Check if minima is a single value
    if len(minima) == 1:
        return [minima]
    
    # Compute the differences between adjacent elements in the `minima` list.
    diffs = np.diff(minima)
    
    # Find the index of the largest difference between adjacent elements.
    max_diff_idx = np.argmax(diffs)
    
    # If the largest difference is less than the minimum distance, return
    # the mean index of the `minima` list.
    if diffs[max_diff_idx] < min_dist:
        return [np.mean(minima)]
    
    # Recursively split the set of local minima at the index with the
    # largest difference.
    left_subset = minima[:max_diff_idx+1]
    right_subset = minima[max_diff_idx+1:]
    subsets = []
    if len(left_subset) > 0:
        subsets += split_minima(left_subset, min_dist)
    if len(right_subset) > 0:
        subsets += split_minima(right_subset, min_dist)
    
    # Compute the mean index for each subset and return the list of means.
    return [np.floor(np.mean(subset)) for subset in subsets]


def find_insertion_sites(data, window_size, num_minima, min_dist, exclude_list=[], max_score=100):
    """
    Uses pLDDT scores from AlphaFold2 structural predictions to 
    identify sites that may be permissive for tag insertion. 
    Specifically, finds local minima in the pLDDT scores. These 
    minima are then filtered by various methods (they must be far
    enough apart, they don't overlap functional sites, etc...).
    
    Parameters:
        data (list): pLDDT scores for prediction
        window_size (int or list): if windowsize is a list, smoothing is performed
            iteratively using window sizes in list.
        num_minima (int): number of sites in return
        min_dist (int): minimum distance between 'clusters' of sites
        exclude_list (list): list of indices to exclude from returned sites
        max_score (int): pLDDT threshold for returning sites

    Returns:
        List of potentially permissive sites for tag insertion
    """    
    
    # Apply sliding window smoothing and identify local minima
    smoothed = smooth_iterative(data, window_size)
    minima_indices, _ = find_peaks(-smoothed, prominence=1)
    
    # Split minima that are too close together
    if len(minima_indices) > 1:
        minima_indices = np.sort(minima_indices)
        minima_indices = split_minima(minima_indices, min_dist)
    
    # Exclude minima in exclude_list
    # these could be functional sites, etc...
    indices = []
    for idx in minima_indices:
        if idx not in exclude_list:
            indices.append(idx)
    
    # Calculate scores for each local minimum and sort
    scores = []
    indices_to_keep = []
    for idx in indices:
        score = smoothed[int(idx)]
        if score < max_score:
            scores.append(smoothed[int(idx)])
            indices_to_keep.append(int(idx))
    sorted_indices = [x for _, x in sorted(zip(scores, indices_to_keep))]
    
    # Return the top num_minima minima
    return sorted_indices[:num_minima]

def pLDDT_tag_sites(pdb_file, max_score, window, passes, nsites, min_dist):
	"""
	Identify potential internal tagging sites using pLDDT from alphafold.

	Parameters:
	- pdb_file (str): name of PDB file from AF2 prediction
	- max_score (int): score threshold to call site
	- window (int): window size for smoothing pLDDT
	- passes (int): number of smoothing passes
	- nsites (int): number of minima to report
	- min_dist (int): minimum distance (in residues) between minima

	Returns: 
	- list: indices of potential tag sites
	- array (2d): smoothed pLDDT, for plotting
	"""

	pdb_parser = PDB.PDBParser()
	structure = pdb_parser.get_structure("foo", pdb_file)

	window_size = passes*[window]

	bf = [float(r["CA"].get_bfactor()) for r in PDB.Selection.unfold_entities(structure, "R")]
	
	minima = find_insertion_sites(np.array(bf), window_size, nsites, min_dist, [], max_score)
	smoothed = smooth_iterative(bf, window_size)

	return minima, smoothed

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	#required arguments
	parser.add_argument('-p', '--pdb', action = 'store', type = str, dest = "PDB", 
					help = "AlphaFold2 PDB file. B-factors must contain pLDDTs for each residue.")
	
	parser.add_argument('-s', '--sites_out', action = 'store', type = str, dest = "SITES_OUTPUT", 
					help = "prints sites to file, default=STDOUT", default = None)
	parser.add_argument('--smoothed_plddt_out', action = 'store', type = str, dest = 'SMOOTHED_OUT', 
					help = "print smoothed pLDDT to file", default = None)

	#optional arguments
	parser.add_argument('--score', action = 'store', type = int, dest = "MAX_SCORE", 
					help = "maximum smoothed pLDDT for calling a tag site (60)", default = 60)
	parser.add_argument('--window', action = 'store', type = int, dest = "SMOOTHING_WINDOW_SIZE", 
					help = "width of window for smoothing pLDDT (20)", default = 20)
	parser.add_argument('--passes', action = 'store', type = int, dest = "SMOOTHING_PASSES", 
					help = "number of pLDDT smoothing passes (2)", default = 2)
	parser.add_argument('--min-dist', action = 'store', type = int, dest = "SITE_DISTANCE", 
					help = "minimum distance between reported tagging sites (20)", 
					default = 20)
	parser.add_argument('--n-sites', action = 'store', type = int, dest = "N_SITES", 
					help = "number of potential tagging sites to report (100)", 
					default = 100) 

	args = parser.parse_args()

	potential_sites, smoothed_plddt = pLDDT_tag_sites(args.PDB, args.MAX_SCORE, 
													  args.SMOOTHING_WINDOW_SIZE, 
													  args.SMOOTHING_PASSES, 
													  args.N_SITES, 
													  args.SITE_DISTANCE)

	#print minimum sites
	if args.SITES_OUTPUT != None:
		f_out = open(args.SITES_OUTPUT, "w")
		print("\n".join([str(x) for x in potential_sites]), file=f_out)
		f_out.close()
	else:
		print("\n".join([str(x) for x in potential_sites]), file=sys.stdout)
		
	#print full smoothed output
	if args.SMOOTHED_OUT != None:
		f_out = open(args.SMOOTHED_OUT, "w")
		print("\n".join(["\t".join([str(x[0]), str(x[1])]) for x in zip(range(len(smoothed_plddt)), smoothed_plddt)]), file=f_out)
		f_out.close()
