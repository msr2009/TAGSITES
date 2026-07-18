import numpy as np
from scipy.signal import find_peaks
import re, string, sys
from Bio import SeqIO, PDB 

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley


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
	sequence = sequence.upper().strip("*")

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
	# Local import to avoid circular dependency (pLDDT_minima imports three_to_one from here)
	from pLDDT_minima import find_insertion_sites, smooth_iterative

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
	Removes trailing asterisks from sequence.
	Simple -- fasta file must have only one record.

	Parameters: 
	- fasta_in (str): name of fasta file

	Returns: 
	- dict: {seqname: sequence}
	"""

	from Bio import SeqIO

	for record in SeqIO.parse(fasta_in, "fasta"):
		return [record.id.strip(","), record.seq.strip("*")]

def resolve_taxids(taxid, taxid_file):
	"""
	Merge a manually-entered taxid (comma-separated) with taxids loaded from
	taxid_file (one per line, '#' comments allowed). Dedupes, preserves order,
	drops the "any species" sentinel ("", "1", "1.0").

	Parameters:
	- taxid (str): manually-entered taxid or comma-separated list
	- taxid_file (str): path to a file of taxids, or None/"" if unused

	Returns:
	- str: comma-joined taxid list, or "" if none given
	"""

	parts = [p.strip() for p in str(taxid).split(",") if p.strip() not in ("", "1", "1.0")]

	if taxid_file:
		with open(taxid_file) as f:
			parts += [ln.split("#")[0].strip() for ln in f if ln.split("#")[0].strip()]

	seen = set()
	deduped = []
	for p in parts:
		if p not in seen:
			seen.add(p)
			deduped.append(p)

	return ",".join(deduped)

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


def _compute_rsasa(struct, max_sasa_dict=None):
	"""
	shared helper: compute relative SASA (0-1) per residue using Shrake-Rupley.
	returns list of [resnum (int), rsasa (float), res_name (str)], in chain/residue order.
	"""

	#these are from Tien 2013
	if not max_sasa_dict:
		max_sasa_dict = {
			'ALA': 121.0, 'ARG': 265.0, 'ASN': 187.0, 'ASP': 187.0, 'CYS': 148.0,
			'GLN': 214.0, 'GLU': 214.0, 'GLY': 97.0, 'HIS': 216.0, 'ILE': 195.0,
			'LEU': 191.0, 'LYS': 230.0, 'MET': 203.0, 'PHE': 228.0, 'PRO': 154.0,
			'SER': 143.0, 'THR': 163.0, 'TRP': 264.0, 'TYR': 255.0, 'VAL': 165.0
		}

	sr = ShrakeRupley()
	sr.compute(struct, level="R")

	rsasa_list = []
	for chain in struct.get_chains():
		for residue in chain:
			res_id = residue.get_id()
			res_name = residue.get_resname()
			rsasa_list.append([res_id[1], residue.sasa/max_sasa_dict[res_name], res_name])

	return rsasa_list


def calc_sasa_shrakerupley(pdb_file, max_sasa_dict=None):

	parser = PDBParser()
	struct = parser.get_structure("foo", pdb_file)
	rsasa = _compute_rsasa(struct, max_sasa_dict)

	# preserve original string-typed return signature
	return [[str(resnum), str(val), res_name] for resnum, val, res_name in rsasa]


def _dilate_and_merge_patches(seed_patches, rep_atoms, atom_to_idx, ns, dilation_radius):
	"""
	grow each seed-based patch to include any spatially nearby residue (no exposure
	or hydrophobicity filter — the peptide backbone is part of the interface even
	when a sidechain points inward), then merge patches whose footprints now overlap.

	Parameters:
	- seed_patches (list): [{"seed_members": [idx, ...]}, ...] (residue indices)
	- rep_atoms (list): representative atom per residue index, or None
	- atom_to_idx (dict): Bio.PDB Atom -> residue index
	- ns (Bio.PDB.NeighborSearch): search built over atom_to_idx's atoms
	- dilation_radius (float): distance (Angstroms) to grow the footprint by; 0 = no dilation

	Returns:
	- list: [{"seed_members": set(idx), "members": set(idx)}, ...], overlapping patches merged
	"""
	dilated = []
	for p in seed_patches:
		seed_members = set(p["seed_members"])
		members = set(seed_members)
		if dilation_radius > 0:
			for i in seed_members:
				atom = rep_atoms[i]
				if atom is None:
					continue
				neighbors = ns.search(atom.coord, dilation_radius, level="A")
				members.update(atom_to_idx[a] for a in neighbors)
		dilated.append({"seed_members": seed_members, "members": members})

	# merge any patches whose dilated footprints now overlap
	merged_again = True
	while merged_again:
		merged_again = False
		for a in range(len(dilated)):
			for b in range(a + 1, len(dilated)):
				if dilated[a]["members"] & dilated[b]["members"]:
					dilated[a]["members"]      |= dilated[b]["members"]
					dilated[a]["seed_members"] |= dilated[b]["seed_members"]
					del dilated[b]
					merged_again = True
					break
			if merged_again:
				break

	return dilated


def calc_hydrophobic_patches(pdb_file, hydro_scores, rsasa_cutoff=0.20, dist_cutoff=8.0,
							  dilation_radius=8.0, min_seed_size=3, max_sasa_dict=None):
	"""
	identify hydrophobic patches on the solvent-exposed surface of a structure.

	method: Lijnzaad, Berendsen & Argos (1996), "A method for detecting hydrophobic
	patches on protein surfaces", Proteins 26:192-203. Surface-exposed hydrophobic
	residues ("seeds") are treated as nodes in a spatial graph (edges = CB-CB, or CA
	for Gly, within dist_cutoff); connected components of that graph are patch cores.
	Cores with fewer than min_seed_size seed residues are dropped as noise. Surviving
	cores are then dilated: any residue within dilation_radius of a seed member joins
	the patch's reported footprint regardless of its own exposure/hydrophobicity, since
	a residue with an interior-pointing sidechain can still have its backbone at the
	interface. Patches whose dilated footprints overlap are merged.

	Parameters:
	- pdb_file (str): PDB file (e.g. AlphaFold model)
	- hydro_scores (dict): per-amino-acid hydrophobicity scores, one-letter keys
		(e.g. loaded from tables/hydrophobicity_kyte-doolittle.tsv via
		calculate_protein_scores.load_scores)
	- rsasa_cutoff (float): minimum relative SASA to call a residue "surface" (default 0.20)
	- dist_cutoff (float): distance (Angstroms) for seed-to-seed core edges (default 8.0)
	- dilation_radius (float): distance (Angstroms) to grow each core's footprint by;
		0 = no dilation, i.e. patch == core (default 8.0)
	- min_seed_size (int): minimum seed-residue count for a core to be reported (default 3)
	- max_sasa_dict (dict): optional override for per-residue-type max SASA (Tien 2013 default)

	Returns:
	- list: per-residue continuous "surface hydrophobic exposure" score, min-max normalized
		to [0, 1], aligned to residue order
	- list: patches, each a dict {"members": [resnum, ...] (dilated footprint),
		"seed_members": [resnum, ...] (core only), "total_hydro": float,
		"total_area_A2": float (absolute exposed area of seed members)}
	"""

	parser = PDBParser()
	struct = parser.get_structure("foo", pdb_file)
	rsasa_list = _compute_rsasa(struct, max_sasa_dict)

	residues = PDB.Selection.unfold_entities(struct, "R")
	resnums   = [r[0] for r in rsasa_list]
	res_names = [r[2] for r in rsasa_list]
	rsasa = np.array([r[1] for r in rsasa_list])
	hydro = np.array([hydro_scores.get(three_to_one(name), 0.0) for name in res_names])

	if not max_sasa_dict:
		max_sasa_dict = {
			'ALA': 121.0, 'ARG': 265.0, 'ASN': 187.0, 'ASP': 187.0, 'CYS': 148.0,
			'GLN': 214.0, 'GLU': 214.0, 'GLY': 97.0, 'HIS': 216.0, 'ILE': 195.0,
			'LEU': 191.0, 'LYS': 230.0, 'MET': 203.0, 'PHE': 228.0, 'PRO': 154.0,
			'SER': 143.0, 'THR': 163.0, 'TRP': 264.0, 'TYR': 255.0, 'VAL': 165.0
		}
	abs_area = np.array([rsasa[i] * max_sasa_dict[res_names[i]] for i in range(len(residues))])

	# representative atom per residue for distance calcs: CB, else CA (e.g. Gly)
	rep_atoms = []
	for residue in residues:
		if "CB" in residue:
			rep_atoms.append(residue["CB"])
		elif "CA" in residue:
			rep_atoms.append(residue["CA"])
		else:
			rep_atoms.append(None)

	atom_to_idx = {atom: i for i, atom in enumerate(rep_atoms) if atom is not None}
	ns = PDB.NeighborSearch(list(atom_to_idx.keys()))

	# continuous score: sum of hydro*rsasa over neighbors within dist_cutoff, per residue
	raw_scores = np.zeros(len(residues))
	for i, atom in enumerate(rep_atoms):
		if atom is None:
			continue
		neighbors = ns.search(atom.coord, dist_cutoff, level="A")
		neighbor_idx = [atom_to_idx[a] for a in neighbors]
		raw_scores[i] = sum(hydro[j] * rsasa[j] for j in neighbor_idx)

	# min-max normalize to [0, 1]
	score_range = raw_scores.max() - raw_scores.min()
	if score_range > 0:
		norm_scores = (raw_scores - raw_scores.min()) / score_range
	else:
		norm_scores = np.zeros(len(raw_scores))

	# patch cores: connected components among surface hydrophobic ("seed") residues
	surface_hydrophobic = [i for i in range(len(residues))
							if rsasa[i] >= rsasa_cutoff and hydro[i] > 0]

	adjacency = {i: [] for i in surface_hydrophobic}
	surface_set = set(surface_hydrophobic)
	for i in surface_hydrophobic:
		atom = rep_atoms[i]
		if atom is None:
			continue
		neighbors = ns.search(atom.coord, dist_cutoff, level="A")
		for a in neighbors:
			j = atom_to_idx.get(a)
			if j is not None and j != i and j in surface_set:
				adjacency[i].append(j)

	# BFS connected components
	visited = set()
	seed_patches = []
	for start in surface_hydrophobic:
		if start in visited:
			continue
		component = []
		queue = [start]
		visited.add(start)
		while queue:
			node = queue.pop()
			component.append(node)
			for neighbor in adjacency[node]:
				if neighbor not in visited:
					visited.add(neighbor)
					queue.append(neighbor)
		if len(component) >= min_seed_size:
			seed_patches.append({"seed_members": component})

	dilated_patches = _dilate_and_merge_patches(seed_patches, rep_atoms, atom_to_idx, ns, dilation_radius)

	patches = []
	for p in dilated_patches:
		seed_idx = sorted(p["seed_members"])
		patches.append({
			"members":      [resnums[i] for i in sorted(p["members"])],
			"seed_members": [resnums[i] for i in seed_idx],
			"total_hydro":  sum(hydro[i] * rsasa[i] for i in seed_idx),
			"total_area_A2": sum(abs_area[i] for i in seed_idx),
		})

	return list(norm_scores), patches
