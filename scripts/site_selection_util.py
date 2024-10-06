import numpy as np
from scipy.signal import find_peaks
import re, string, sys
from Bio import SeqIO, PDB 

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
	contains unknown bases. From ChatGPT and edited by me.

	If unknown bases, then use RE to see if it's a 
	uniprot accession.

	Parameters:
		sequence (str): The sequence to check.

	Returns:
		str:	"DNA" if the sequence is a DNA sequence, 
				"protein" if the sequence is a protein sequence, or 
				"unknown" if the sequence contains unknown bases.
	"""
	from re import match
   
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

def check_input_type(inputfile):
	"""
	determines whether input is PDB or FASTA.
	currently very dumb -- looks at file extension lol
	
	Parameters: 
	- inputfile (str): name of input file

	Returns:
	- inputtype (str): fasta or pdb
	"""
	
	file_extension = inputfile.split(".")[-1]

	if uniprot_accession_regex(inputfile) != None:
		return("uniprot")
	elif file_extension == "pdb":
		return("pdb")
	elif file_extension in ["fasta", "fsa", "fa"]:
		return("fasta")
	else:
		return("err")
	
def save_fasta(name, seq, outfile):
	"""
	writes fasta very simply.

	Parameters: 
	- name (str): name for header
	- seq (str): sequence
	- outfile (str): name for output file
	"""

	f_out = open(outfile, "w")
	print(">{}".format(name), file=f_out)
	print(seq, file=f_out)
	f_out.close()
	
	return

def uniprot_accession_regex(_str):
	"""
	according to uniprot (https://www.uniprot.org/help/accession_numbers)
	this regex will match all accessions

	Parameters:
		- _str (str): string to be matched by regex
	
	Returns:
		- None if no match
		- match object if match
	"""
	from re import match
	return match("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", _str)

def ncbiblast_call(script_dir, email, seq, db, scores, evalue, taxid, outfile):
	"""
	writes ebi ncbiblast call using common parameters
	"""
	
	if taxid == 1:
		ncbi_call = "python {}ncbiblast.py --email {} --program blastp --stype protein\
					--sequence {} --database {} --outformat json --alignments 0\
					--scores {} --exp {} --outfile {} --pollFreq 5".format(script_dir, email, seq, db, scores, evalue, outfile)
	else:
		ncbi_call = "python {}ncbiblast.py --email {} --program blastp --stype protein\
					--sequence {} --database {} --outformat json --alignments 0\
					--scores {} --exp {} --outfile {} --pollFreq 5 --taxid {}".format(script_dir, email, seq, db, scores, evalue, outfile, taxid)
	return ncbi_call


def extract_seq_from_pdb(pdb_file):
	"""
	use biopython to extract the amino acid sequence from a pdb file
	"""
	pdb_parser = PDB.PDBParser()
	structure = pdb_parser.get_structure("foo", pdb_file)
	seq = "".join([three_to_one(r.get_resname()) for r in PDB.Selection.unfold_entities(structure, "R")])
	return seq

def extract_bfactors_from_pdb(pdb_file):
	"""
	use biopython to extract b-factor values from a pdb file
	"""
	pdb_parser = PDB.PDBParser()
	structure = pdb_parser.get_structure("foo", pdb_file)
	bf = [float(r["CA"].get_bfactor()) for r in PDB.Selection.unfold_entities(structure, "R")]
	return bf

