"""
score_conservation_py3.py

reimplementing Capra and Singh's JS divergence code in python3

Matt Rich, 2024
"""

import math, sys

PSEUDOCOUNT = .0000001

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"] 

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {}
for i, aa in enumerate(amino_acids):
    aa_to_index[aa] = i

################################################################################
# Frequency Count and Gap Penalty
################################################################################

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
	""" Return the weighted frequency count for a column--with pseudocount."""

	# if the weights do not match, use equal weight
	if len(seq_weights) != len(col):	
		seq_weights = [1.] * len(col)

	aa_num = 0
	freq_counts = len(amino_acids)*[pc_amount] # in order defined by amino_acids

	for aa in amino_acids:
		for j in range(len(col)):
			if col[j] == aa:
				freq_counts[aa_num] += 1 * seq_weights[j]
		aa_num += 1

	for j in range(len(freq_counts)):
		freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

	return freq_counts

def weighted_gap_penalty(col, seq_weights):
	""" Calculate the simple gap penalty multiplier for the column. If the 
	sequences are weighted, the gaps, when penalized, are weighted 
	accordingly. """

	# if the weights do not match, use equal weight
	if len(seq_weights) != len(col):
		seq_weights = [1.] * len(col)
    
	gap_sum = 0.
	for i in range(len(col)):
		if col[i] == '-':
			gap_sum += seq_weights[i]

	return 1 - (gap_sum / sum(seq_weights))


def gap_percentage(col):
	"""Return the percentage of gaps in col."""
	num_gaps = 0.

	for aa in col:
		if aa == '-': num_gaps += 1

	return num_gaps / len(col)

################################################################################
# Shannon Entropy
################################################################################

def shannon_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	"""Calculates the Shannon entropy of the column col. sim_matrix  and 
	bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The 
	entropy will be between zero and one because of its base. See p.13 of 
	Valdar 02 for details. The information score 1 - h is returned for the sake 
	of consistency with other scores."""

	fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

	h = 0. 
	for i in range(len(fc)):
		if fc[i] != 0:
			h += fc[i] * math.log(fc[i])

#	h /= math.log(len(fc))
	h /= math.log(min(len(fc), len(col)))

	inf_score = 1 - (-1 * h)

	if gap_penalty == 1: 
		return inf_score * weighted_gap_penalty(col, seq_weights)
	else: 
		return inf_score

################################################################################
# Property Entropy
################################################################################

def property_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	"""Calculate the entropy of a column col relative to a partition of the 
	amino acids. Similar to Mirny '99. sim_matrix and bg_distr are ignored, but 
	could be used to define the sets. """

	# Mirny and Shakn. '99
	property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

	# Williamson '95
	# property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C'], ['-']]

	fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

	# sum the aa frequencies to get the property frequencie
	prop_fc = [0.] * len(property_partition)
	for p in range(len(property_partition)):
		for aa in property_partition[p]:
			prop_fc[p] += fc[aa_to_index[aa]]

	h = 0. 
	for i in range(len(prop_fc)):
		if prop_fc[i] != 0:
			h += prop_fc[i] * math.log(prop_fc[i])

	h /= math.log(min(len(property_partition), len(col)))

	inf_score = 1 - (-1 * h)

	if gap_penalty == 1: 
		return inf_score * weighted_gap_penalty(col, seq_weights)
	else: 
		return inf_score


################################################################################
# Property Relative Entropy
################################################################################

def property_relative_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	"""Calculate the relative entropy of a column col relative to a
	partition of the amino acids. Similar to Williamson '95. sim_matrix is
	ignored, but could be used to define the sets. See shannon_entropy()
	for more general info. """

	# Mirny and Shakn. '99
	#property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

	# Williamson '95
	property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']]

	prop_bg_freq = []
	if len(bg_distr) == len(property_partition):
		prop_bg_freq = bg_distr
	else:
		prop_bg_freq = [0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000] # from BL62

	#fc = weighted_freq_count_ignore_gaps(col, seq_weights)
	fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

	# sum the aa frequencies to get the property frequencies
	prop_fc = [0.] * len(property_partition)
	for p in range(len(property_partition)):
		for aa in property_partition[p]:
			prop_fc[p] += fc[aa_to_index[aa]]

	d = 0. 
	for i in range(len(prop_fc)):
		if prop_fc[i] != 0 and prop_bg_freq[i] != 0:
			d += prop_fc[i] * math.log(prop_fc[i] / prop_bg_freq[i], 2)


	if gap_penalty == 1: 
		return d * weighted_gap_penalty(col, seq_weights)
	else: 
		return d

################################################################################
# von Neumann Entropy
################################################################################

def vn_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	""" Calculate the von Neuman Entropy as described in Caffrey et al. 04.
	This code was adapted from the implementation found in the PFAAT project 
	available on SourceForge. bg_distr is ignored."""

	aa_counts = [0.] * 20
	for aa in col:
		if aa != '-': aa_counts[aa_to_index[aa]] += 1

	dm_size = 0
	dm_aas = []
	for i in range(len(aa_counts)):
		if aa_counts[i] != 0:
			dm_aas.append(i)
			dm_size += 1

	if dm_size == 0: return 0.0

	row_i = 0
	col_i = 0
	dm = zeros((dm_size, dm_size), Float32)
	for i in range(dm_size):
		row_i = dm_aas[i]
		for j in range(dm_size):
			col_i = dm_aas[j]
			dm[i][j] = aa_counts[row_i] * sim_matrix[row_i][col_i]

	ev = la.eigenvalues(dm).real

	temp = 0.
	for e in ev:
		temp += e

	if temp != 0:
		for i in range(len(ev)):
			ev[i] = ev[i] / temp

	vne = 0.0
	for e in ev:
		if e > (10**-10):
			vne -= e * math.log(e) / math.log(20)

	if gap_penalty == 1: 
	#return (1-vne) * weighted_gap_penalty(col, seq_weights)
		return (1-vne) * weighted_gap_penalty(col, [1.] * len(col))
	else: 
		return 1 - vne

################################################################################
# Relative Entropy 
################################################################################

def relative_entropy(col, sim_matix, bg_distr, seq_weights, gap_penalty=1):
	"""Calculate the relative entropy of the column distribution with a 
	background distribution specified in bg_distr. This is similar to the 
	approach proposed in Wang and Samudrala 06. sim_matrix is ignored."""

	distr = bg_distr[:]

	fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

	# remove gap count
	if len(distr) == 20: 
		new_fc = fc[:-1]
		s = sum(new_fc)
		for i in range(len(new_fc)):
			new_fc[i] = new_fc[i] / s
		fc = new_fc

	if len(fc) != len(distr): return -1

	d = 0.
	for i in range(len(fc)):
		if distr[i] != 0.0:
			d += fc[i] * math.log(fc[i]/distr[i])

	d /= math.log(len(fc))

	if gap_penalty == 1: 
		return d * weighted_gap_penalty(col, seq_weights)
	else: 
		return d

################################################################################
# Jensen-Shannon Divergence
################################################################################

def js_divergence(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	""" Return the Jensen-Shannon Divergence for the column with the background
	distribution bg_distr. sim_matrix is ignored. JSD is the default method."""

	distr = bg_distr[:]

	fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

	# if background distrubtion lacks a gap count, remove fc gap count
	if len(distr) == 20: 
		new_fc = fc[:-1]
		s = sum(new_fc)
		for i in range(len(new_fc)):
			new_fc[i] = new_fc[i] / s
		fc = new_fc

	if len(fc) != len(distr): return -1

	# make r distriubtion
	r = []
	for i in range(len(fc)):
		r.append(.5 * fc[i] + .5 * distr[i])

	d = 0.
	for i in range(len(fc)):
		if r[i] != 0.0:
			if fc[i] == 0.0:
				d += distr[i] * math.log(distr[i]/r[i], 2)
			elif distr[i] == 0.0:
				d += fc[i] * math.log(fc[i]/r[i], 2) 
			else:
				d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)
		
	# d /= 2 * math.log(len(fc))
	d /= 2

	if gap_penalty == 1: 
		return d * weighted_gap_penalty(col, seq_weights)
	else: 
		return d

################################################################################
# Mutation Weighted Pairwise Match
################################################################################

def sum_of_pairs(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
	""" Sum the similarity matrix values for all pairs in the column. 
	This method is similar to those proposed in Valdar 02. bg_distr is ignored."""

	xsum = 0.
	max_sum = 0.

	for i in range(len(col)):
		for j in range(i):
			if col[i] != '-' and col[j] != '-':
				max_sum += seq_weights[i] * seq_weights[j]
				xsum += seq_weights[i] * seq_weights[j] * sim_matrix[aa_to_index[col[i]]][aa_to_index[col[j]]]

	if max_sum != 0: 
		xsum /= max_sum
	else:
		xsum = 0.

	if gap_penalty == 1: 
		return xsum * weighted_gap_penalty(col, seq_weights)
	else:
		return xsum



################################################################################
# Window Score
################################################################################

def window_score(scores, window_len, lam=.5):
	""" This function takes a list of scores and a length and transforms them 
	so that each position is a weighted average of the surrounding positions. 
	Positions with scores less than zero are not changed and are ignored in the 
	calculation. Here window_len is interpreted to mean window_len residues on 
	either side of the current residue. """

	w_scores = scores[:]

	for i in range(window_len, len(scores) - window_len):
		if scores[i] < 0: 
			continue

		xsum = 0.
		num_terms = 0.
		for j in range(i - window_len, i + window_len + 1):
			if i != j and scores[j] >= 0:
				num_terms += 1
				xsum += scores[j]

		if num_terms > 0:
			w_scores[i] = (1 - lam) * (xsum / num_terms) + lam * scores[i]

	return w_scores


def calc_z_scores(scores, score_cutoff):
	"""Calculates the z-scores for a set of scores. Scores below
	score_cutoff are not included."""

	average = 0.
	std_dev = 0.
	z_scores = []
	num_scores = 0

	for s in scores:
		if s > score_cutoff:
			average += s
			num_scores += 1
	if num_scores != 0:
		average /= num_scores

	for s in scores:
		if s > score_cutoff:
			std_dev += ((s - average)**2) / num_scores
	std_dev = math.sqrt(std_dev)

	for s in scores:
		if s > score_cutoff and std_dev != 0:
			z_scores.append((s-average)/std_dev)
		else:
			z_scores.append(-1000.0)

	return z_scores


################################################################################
################################################################################
################################################################################
#  END CONSERVATION SCORES
################################################################################
################################################################################
################################################################################

def read_scoring_matrix(sm_file):
	""" Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
	as an array. """
	aa_index = 0
	first_line = 1
	row = []
	list_sm = [] # hold the matrix in list form

	try:
		matrix_file = open(sm_file, 'r')

		for line in matrix_file:

			if line[0] != '#' and first_line:
				first_line = 0
				if len(amino_acids) == 0:
					for c in line.split():
						aa_to_index[string.lower(c)] = aa_index
						amino_acids.append(string.lower(c))
						aa_index += 1

			elif line[0] != '#' and first_line == 0:
				if len(line) > 1:
					row = line.split()
					list_sm.append(row)

	except IOError:
		print("Could not load similarity matrix: %s. Using identity matrix..." % sm_file)
		return identity(20)
	
	# if matrix is stored in lower tri form, copy to upper
	if len(list_sm[0]) < 20:
		for i in range(0,19):
			for j in range(i+1, 20):
				list_sm[i].append(list_sm[j][i])

	for i in range(len(list_sm)):
		for j in range(len(list_sm[i])):
			list_sm[i][j] = float(list_sm[i][j])

	return list_sm
	#sim_matrix = array(list_sm, type=Float32)
	#return sim_matrix

def calculate_sequence_weights(msa):
	""" Calculate the sequence weights using the Henikoff '94 method
	for the given msa. """

	seq_weights = [0.] * len(msa)
	for i in range(len(msa[0])):
		freq_counts = [0] * len(amino_acids)

		col = []
		for j in range(len(msa)):
			if msa[j][i] != '-': # ignore gaps
				freq_counts[aa_to_index[msa[j][i]]] += 1

		num_observed_types = 0
		for j in range(len(freq_counts)):
			if freq_counts[j] > 0: num_observed_types +=1

		for j in range(len(msa)):
			d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
			if d > 0:
				seq_weights[j] += 1. / d

	for w in range(len(seq_weights)):
		seq_weights[w] /= len(msa[0])

	return seq_weights


def load_sequence_weights(fname):
	"""Read in a sequence weight file f and create sequence weight list. 
	The weights are in the same order as the sequences each on a new line. """
	seq_weights = []

	try:
		f = open(fname)

		for line in f:
			l = line.split()
			if line[0] != '#' and len(l) == 2:
				seq_weights.append(float(l[1]))

	except IOError:
		pass
		#print "No sequence weights. Can't find: ", fname

	return seq_weights

def get_column(col_num, alignment):
	"""Return the col_num column of alignment as a list."""
	col = []
	for seq in alignment:
		if col_num < len(seq): col.append(seq[col_num])

	return col

def get_distribution_from_file(fname):
	""" Read an amino acid distribution from a file. The probabilities should
	be on a single line separated by whitespace in alphabetical order as in 
	amino_acids above. # is the comment character."""

	distribution = []
	try:
		f = open(fname)
		for line in f:
			if line[0] == '#': continue
			line = line[:-1]
			distribution = line.split()
			distribution = map(float, distribution)
		
	except IOError:
		print("Using default (BLOSUM62) background.")
		return []
	
	# use a range to be flexible about round off
	if .997 > sum(distribution) or sum(distribution) > 1.003:
		print("Distribution does not sum to 1. Using default (BLOSUM62) background.")
		print(sum(distribution))
		return []

	return distribution

def read_fasta_alignment(filename):
	""" Read in the alignment stored in the FASTA file, filename. Return two
	lists: the identifiers and sequences. """

	f = open(filename)

	names = []
	alignment = []
	cur_seq = ''

	for line in f:
		line = line[:-1]
		if len(line) == 0: continue

		if line[0] == ';': continue
		if line[0] == '>':
			names.append(line[1:].replace('\r', ''))

			if cur_seq != '':
				cur_seq = cur_seq.upper()
				for i, aa in enumerate(cur_seq):
					if aa not in iupac_alphabet:
						cur_seq = cur_seq.replace(aa, '-')
				alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
				cur_seq = ''
		elif line[0] in iupac_alphabet:
			cur_seq += line.replace('\r', '')

	# add the last sequence
	cur_seq = cur_seq.upper()
	for i, aa in enumerate(cur_seq):
		if aa not in iupac_alphabet:
			cur_seq = cur_seq.replace(aa, '-')
	alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

	return names, alignment

def read_clustal_alignment(filename):
	""" Read in the alignment stored in the CLUSTAL file, filename. Return
	two lists: the names and sequences. """

	names = []
	alignment = []

	f = open(filename)

	for line in f:
		line = line[:-1]
		if len(line) == 0: continue
		if '*' in line: continue

		if 'CLUSTAL' in line: continue

		t = line.split()

		if len(t) == 2 and t[1][0] in iupac_alphabet:
			if t[0] not in names:
				names.append(t[0])
				alignment.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
			else:
				alignment[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X','-').replace('\r', '')
		   
	return names, alignment

def main(align_file, window_size, win_lam, outfile_name, 
				s_matrix_file, bg_distribution, scoring_function, use_seq_weights, 
				gap_cutoff, use_gap_penalty, seq_specific_output,
				normalize_scores):
	

	# read scoring matrix
	s_matrix = read_scoring_matrix(s_matrix_file)

	names = []
	alignment = []
	seq_weights = []

	align_suffix = align_file.split('.')[-1]

	# read input alignment
	try: 
		names, alignment = read_clustal_alignment(align_file)
		if names == []:
			names, alignment = read_fasta_alignment(align_file)
	except IOError:
		print("Could not find %s. Exiting...".format(align_file))
		sys.exit(1)

	if len(alignment) != len(names) or alignment == []:
		print("Unable to parse alignment.")
		sys.exit(1)

	# make sure seq lengths make sense
	seq_len = len(alignment[0])
	for i, seq in enumerate(alignment):
		if len(seq) != seq_len:
			print("ERROR: Sequences of different lengths: %s (%d) != %s (%d).\n" % (names[0], seq_len, names[i], len(seq)))
			sys.exit(1)

	# deal with sequence weighting
	if use_seq_weights:	
		seq_weights = load_sequence_weights(align_file.replace('.%s' % align_suffix, '.weights'))
		if seq_weights == []:
			seq_weights = calculate_sequence_weights(alignment)

	if len(seq_weights) != len(alignment): seq_weights = [1.] * len(alignment)

	# handle print of output relative to specific sequence
	ref_seq_num = None
	if seq_specific_output and seq_specific_output not in names:
		print("Sequence %s not found in alignment. Using default output format...\n" % seq_specific_output)
		seq_specific_output = 0
	elif seq_specific_output in names:
		ref_seq_num = names.index(seq_specific_output)

	# calculate scores
	scores = []
	for i in range(len(alignment[0])):
		col = get_column(i, alignment)

		if len(col) == len(alignment):
			if gap_percentage(col) <= gap_cutoff:
				scores.append(scoring_function(col, s_matrix, bg_distribution, seq_weights, use_gap_penalty))
			else:
				scores.append(-1000.)
	
	if window_size > 0:
		scores = window_score(scores, window_size, win_lam)

	if normalize_scores:
		scores = calc_z_scores(scores, -999)

	#print to file or stdout
	#headers -- doesn't currently output background
	print("# {} -- {} - \
			window_size: {} - \
			window lambda: {} - \
			background: {} - \
			seq. weighting: {} - \
			gap penalty: {} - \
			normalized: {}".format(align_file, scoring_function.__name__,
								   window_size, win_lam, "background",
								   use_seq_weights, use_gap_penalty,
								   normalize_scores), file=outfile_name)
	if seq_specific_output:
		print("# reference sequence: {}".format(seq_specific_output), file=outfile_name)
		print("# {}".format("\t".join(["align_column_number", "score", "amino acid"])), file=outfile_name)
	else:
		print("# {}".format("\t".join(["align_column_number", "score", "column"])), file=outfile_name)
	
	#print scores
	for i, score in enumerate(scores):
		if seq_specific_output:
			cur_aa = get_column(i, alignment)[ref_seq_num]
			if cur_aa == "-": continue
			print("{}\t{:.5f}\t{}".format(i, score, cur_aa), file=outfile_name)
		else:
			print("{}\t{:.5f}\t{}".format(i, score, "".join(get_column(i, alignment))), file = outfile_name)


if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	
	parser.add_argument('-i', '--input', action = 'store', type = str, dest = 'INPUT', 
		help = "input alignment in either FASTA or CLUSTAL format [filename]",
		required = True)

	parser.add_argument('-a', action = 'store', type = str, dest = 'REF', 
		help = "reference sequence. Print scores in reference to a specific sequence (ignoring gaps). Default prints the entire column. [sequence name]")
	
	parser.add_argument('-b', action = 'store', type = float, dest = "LAMBDA",
		help = "lambda for window heuristic linear combination. Default=.5 [real in [0,1]]", default=0.5)

	parser.add_argument('-d', action = 'store', type = str, dest = 'BACKGROUND', 
		help = "background distribution file, e.g., swissprot.distribution. Default=BLOSUM62 background [filename]")

	parser.add_argument('-g', action = 'store', type = float, dest = 'GAPCUTOFF', 
		help = "gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps. Default=.3 [real in [0, 1)]", default = 0.3)

	parser.add_argument('-l', action = 'store_false', dest = 'SEQWEIGHT', 
		help = "use sequence weighting. Default=True [True|False]", default=True)

	parser.add_argument('-m', action = 'store', type = str, dest = 'SIMMATRIX', 
		help = "similarity matrix file, e.g., matrix/blosum62.bla or .qij. Default=identity matrix [filename]",
		default = "./matrix/blosum62.bla")

	parser.add_argument('-n', action = 'store_true', dest = 'NORMALIZE', 
		help = "normalize scores. Print the z-score (over the alignment) of each column raw score. Default=False",
		default=False)

	parser.add_argument('-o', action = 'store', type = str, dest = 'OUTFILE', 
		help = "name of output file. Default=STDOUT [filename]")

	parser.add_argument('-p', action = 'store_false', dest = 'GAPPENALTY', 
		help = "use gap penalty. Lower the score of columns that contain gaps. Default=True [True|False]", 
		default=True)

	parser.add_argument('-s', action = 'store', type = str, dest = 'SCORE', 
		choices=["shannon_entropy", "property_entropy",
				"property_relative_entropy", "vn_entropy", "relative_entropy",
				"js_divergence", "sum_of_pairs"], 
		help = "conservation estimation method.", default="js_divergence")

	parser.add_argument('-w', action = 'store', type = int, dest = 'WINDOW', 
		help = "window size. Number of residues on either side included in the window. Default=3 [int]", 
		default = 3)

	args = parser.parse_args()
	
	if args.OUTFILE != None:
		outfile = open(args.OUTFILE, "w")
	else:
		outfile = sys.stdout
	
	bg_dist = []
	bg_name = []
	
	blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
	
	if args.BACKGROUND != None:
		bg_dist = get_distribution_from_file(args.BACKGROUND)
		if bg_dist == []:
			print("can't read background distribution file, using BLOSUM62 instead")
			bg_dist = blosum_background_distr
	else:
		bg_dist = blosum_background_distr

	
	main(args.INPUT, args.WINDOW, args.LAMBDA, outfile, 
		 args.SIMMATRIX, bg_dist, eval(args.SCORE), args.SEQWEIGHT, 
		 args.GAPCUTOFF, args.GAPPENALTY, args.REF, args.NORMALIZE)
	
	outfile.close()
