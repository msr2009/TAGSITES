import numpy as np
from scipy.signal import find_peaks
import re, string, sys
from Bio import SeqIO, PDB 

def filter_functional_sites(seq, sites_to_filter):
    """
    function for finding known functional sites in sequence
    e.g., myristoylation = N-MGxxxS
    
    Parameters:
        seq (str): sequence of protein
        sites_to_filter (dict): dictionary of "function: [regex]"
    
    Returns:
        list of sites to be excluded. 
    """
    
    filtered_indices = []
    
    for r in sites_to_filter:
        for match in re.finditer(sites_to_filter[r], seq):
            for i in range(match.start(), match.end()):
                filtered_indices.append([r, i])
    
    return list(zip(*filtered_indices))

def smooth_iterative(data, windowsize):
    """
    function to (potentially) iteratively smooth data. Calls numpy convolve
    for smoothing.
    
    Parameters:
        data (list): data to be smoothed
        windowsize (int or list): if windowsize is a list, smoothing is performed
            iteratively using window sizes in list.
        
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

def insert_tag(target_seq, tag_seq, site):
    """
    Function for inserting tag sequences at specific sites
    in a target protein sequence.
    
    Parameters:
        target_seq (str): FASTA string of protein getting tagged
        tag_seq (str): FASTA string of tag being inserted
        site (int or list): site or list of sites for tag insertion

    Both sequences can be either protein or DNA. If DNA, coding
    sequence must be UPPERCASE and introns (if any) must be 
    lowercase.
    
    Returns: FASTA string of tagged targets. 
    """   
    return

def get_sequence(infile):
    """
    Wrapper for reading sequences either from 
        (a) a FASTA file or 
        (b) a PDB file.
        
    Uses Biopython for both. Translates if necessary,
    after removing lowercase noncoding sequence.
    
    Parameters: 
        infile (str)
    
    Returns: protein sequence (str).
    """
    
    # if it's a PDB, we can read the 
    # protein seq from the PDB file
    if infile.endswith(".pdb"):
        for s in SeqIO.parse(infile, "pdb-seqres"):
            return str(s.seq)
            break
            
    # if it's a fasta file, we need to determine 
    # (1) if it's DNA or protein
    # (2) if there are lowercase letters (only DNA)
    elif infile.endswith(".fasta") or infile.endswith(".fa"):
        for s in SeqIO.parse(infile, "fasta"): 
            # return the seq if it's protein
            if sequence_type(s.seq) == "protein":
                return str(s.seq)
            # if DNA, then first remove lowercase 
            elif sequence_type(s.seq) == "dna":
                #then translate uppercase to protein
                return translate_seq(remove_lowercase(str(s.seq)))
            elif sequence_type(s.seq) == "unknown":
                raise ValueError("not a valid sequence!")
            

def sequence_type(sequence):
    """
    Determines whether a sequence is a DNA, protein, or 
    contains unknown bases. From ChatGPT.

    Parameters:
        sequence (str): The sequence to check.

    Returns:
        str: "DNA" if the sequence is a DNA sequence, 
             "protein" if the sequence is a protein sequence, or 
             "unknown" if the sequence contains unknown bases.
    """
    dna_bases = set('ACTGN')
    protein_bases = set('ACDEFGHIKLMNPQRSTVWYX')

    # Convert the sequence to uppercase to make it case-insensitive
    sequence = sequence.upper()

    # Check whether the sequence contains only DNA bases
    if set(sequence).issubset(dna_bases):
        return "dna"

    # Check whether the sequence contains protein bases
    elif set(sequence).issubset(protein_bases):
        return "protein"

    # If the sequence contains unknown bases, return "unknown"
    else:
        return "unknown"
    
def remove_lowercase(seq):
    """
    Function to remove lowercase letters from a sequence. 
    
    From stackoverflow: https://stackoverflow.com/questions/15437589/fast-way-to-remove-lowercase-substrings-from-string
    """
    table = str.maketrans('', '', string.ascii_lowercase)
    return seq.translate(table)

def translate_seq(dna_seq):
    """
    Translates a DNA sequence into a protein sequence.

    Parameters:
        dna_seq (str): The DNA sequence to be translated.

    Returns: the translated protein sequence (str).
    """
    # Define the genetic code as a dictionary
    codons = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # Make sure the input sequence length is a multiple of 3
    if len(dna_seq) % 3 != 0:
        raise ValueError("Input sequence length must be a multiple of 3.")

    # Translate the sequence
    translated_seq = ""
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if codon in codons:
            translated_seq += codons[codon]
        else:
            translated_seq += "X"

    return translated_seq

def three_to_one(three_letter_code):
    """
    Convert three-letter amino acid codes to one-letter amino acid codes.

    Parameters:
    - three_letter_code (str): The three-letter amino acid code to be converted.

    Returns:
    - str: The corresponding one-letter amino acid code, or 'X' if the input is not recognized.
    """
    # Define a dictionary mapping three-letter codes to one-letter codes
    aa_mapping = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    # Convert the input to uppercase for case-insensitivity
    three_letter_code = three_letter_code.upper()

    # Use dictionary lookup to get the corresponding one-letter code, or 'X' if not recognized
    one_letter_code = aa_mapping.get(three_letter_code, 'X')

    return one_letter_code

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

def read_fasta(fasta_in):
	"""
	Wrapper function for reading a FASTA file using Biopython.
	Simple -- fasta file must have only one record.

	Parameters: 
	- fasta_in (str): name of fasta file

	Returns: 
	- dict: {seqname: sequence}
	"""

	from Bio import SeqIO

	for record in SeqIO.parse(fasta_in, "fasta"):
		return [record.id.strip(","), record.seq]
