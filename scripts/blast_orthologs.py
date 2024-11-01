"""
blast_orthologs.py

Matt Rich, 9/2024
"""

import subprocess, json
from site_selection_util import read_fasta, ncbiblast_call

from pathlib import Path

def hit_to_dict(j):
	"""
	converts json blast output to something a little easier for us to deal with.

	parameters
		- j (dict): json dictionary
	returns 
		- hit_dict (dict)
	"""
	hit_dict = {}
	hit_dict["acc"] = j["hit_acc"]
	hit_dict["species"] = j["hit_os"]
	hit_dict["evalue"] = float(j["hit_hsps"][0]["hsp_expect"])
	hit_dict["identity"] = float(j["hit_hsps"][0]["hsp_identity"])
	hit_dict["length"] = int(j["hit_hsps"][0]["hsp_hit_to"]) - int(j["hit_hsps"][0]["hsp_hit_from"])
	hit_dict["hitseq"] = j["hit_hsps"][0]["hsp_hseq"].replace("-", "")

	return hit_dict

def main(fasta_in, email, workingdir, name, output, 
				n, evalue, db, length_percent,
				align_full_seqs, taxid, clients_folder, exclude_paralogs):

	name, seq = read_fasta(fasta_in)
	seq_len = float(len(seq))
	
	out_prefix = Path(output).with_suffix('')

	############################
	# BLAST TO FIND ORTHOLOGS
	###########################

	#blast seq against database
	ncbi_call = ncbiblast_call(clients_folder, email, seq, db, n, evalue, taxid, out_prefix)

	print(ncbi_call)
	#call ncbiblast command
	subprocess.run(ncbi_call, shell=True)

	###########################
	# PROCESS BLAST OUTPUT
	###########################

	#processing ncbiblast output
	#read in from json output
	with open("{}.json.json".format(out_prefix)) as f:
		blast_output = json.load(f)

	#relevant details:
	qlen = blast_output["query_len"]
	#hits in blast_output["hits"]

	#save query species
	query_species = blast_output["hits"][0]["hit_os"]

	blast_hits = {}
	blast_hits[query_species] = []

	for h in blast_output["hits"]:
		
		d_hit = hit_to_dict(h)
		
		#skip hits that don't meet evalue and match length criteria
		if d_hit["evalue"] > evalue:
			continue
		if d_hit["length"]/float(len(seq)) < length_percent or d_hit["length"]/float(len(seq)) > 1/length_percent:
			continue
		
		if exclude_paralogs:	#do shenanigans to filter for best hits
			#store all isoforms/paralogs from target species
			if d_hit["species"] == query_species:
				blast_hits[query_species].append(d_hit)
			#otherwise, save teh best hit from other species
			else:
				#if we haven't seen that species yet
				if d_hit["species"] not in blast_hits:
					blast_hits[d_hit["species"]] = [d_hit]
				#otherwise, we need to see if the new hit is better
				###I think we can just assume that this is true?
				elif blast_hits[d_hit["species"]][0]["evalue"] > d_hit["evalue"]:
					blast_hits[d_hit["species"]] = [d_hit]
		else:
			if d_hit["species"] in blast_hits:
				blast_hits[d_hit["species"]].append(d_hit)
			else:
				blast_hits[d_hit["species"]] = [d_hit]
		
		#check each time if there are _n_ hits in the blast_hits dict, 
		#stop iterating once we get enough.
		if len(blast_hits) >= n: 
			break

	
	###########################
	# GET FULL SEQS OF HITS
	###########################

	#now we have blast hits and can make a fasta file to submit to an aligner. 
	fasta_str = ""
	for s in blast_hits:
		for h in blast_hits[s]:
			if align_full_seqs: 
				## then we need to dbfetch the accessions we want
				dbfetch_call = "python {}/dbfetch.py fetchData UNIPROT:{} fasta raw".format(clients_folder, h["acc"])
#				print(dbfetch_call)
				
				acc_fasta = subprocess.run(dbfetch_call, shell=True,capture_output=True, text=True).stdout
				
				#here, filter for length of full seq
				hit_len = float(len("".join(acc_fasta.split("\n")[1:])))
				if hit_len/seq_len < length_percent or hit_len/seq_len > 1/length_percent:
						continue

				#truncate the name for each seq	
				acc_fasta_name = acc_fasta.split("\n")[0].split()[0]
				acc_fasta_seq = "\n".join(acc_fasta.split("\n")[1:])
				
				#add to fasta_str
				fasta_str += acc_fasta_name+"\n"
				fasta_str += acc_fasta_seq.rstrip('\n')+"\n"
#				print(fasta_str)
			else:
				#we can make a fasta str (and save a file) for the alignment
				fasta_str += ">{}\n{}\n".format(h["acc"], h["hitseq"])

	
	###########################
	# ALIGN WITH CLUSTALO
	###########################

	#output fasta of seqs to align
	fasta_out = open("{}.fasta".format(out_prefix), "w")
	print(fasta_str, file=fasta_out)
	fasta_out.close()

	#call clustalomega to realign ortholog seqs
	# I can't figure out how to format these for submission to EBI, so 
	# I think it's better to simply run clustal locally. (below)
#	ebi_clustal_call = "python {}/clustalo.py --email {} --sequence {} --stype protein --outfile {}".format(clients_folder, email, fasta_str, out_prefix)
#	subprocess.run(ebi_clustal_call, shell=True)

	local_clustal_call = "clustalo -i {}.fasta --outfmt fa -o {}.aln --force".format(out_prefix, out_prefix)
#	local_clustal_call = "clustalo -i {}/{}.blast_close.fasta".format(workingdir, name)
	print(local_clustal_call)
	subprocess.run(local_clustal_call, shell=True)
	
	###########################
	# CALCULATE JSD
	###########################

	#calculate jensen-shannon divergence (from Capra and Singh)
	#we need the name of the best hit, which will be our input seq, 
	#this is the first seq in the fasta file
	best_hit_name = open("{}.aln".format(out_prefix), "r").readline()[1:]

	js_call = 'python {}score_conservation_py3.py -i {}.aln -a "{}" -o {}.jsd -m {}matrix/blosum62.bla'.format(clients_folder, out_prefix,
																					  best_hit_name.rstrip(),
																					  out_prefix, clients_folder)

	print(js_call)
	subprocess.run(js_call, shell=True)
	return 0

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	parser.add_argument('-f', '--fasta', '--input_file', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)	

	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', '--working_dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--name', '--run_name', action='store', type=str, dest='NAME', 
		help = "prefix name for output", required=True)
	parser.add_argument('--output', action='store', type=str, dest='OUTPUT',
		help = "user-supplied output filename", default=None)
		
	parser.add_argument('--taxid', action='store', type=str, dest='TAXID',
		help = "taxid to use for blast search", default = 1)
	parser.add_argument('--taxid_file', action='store', type=str, dest='TAX_FILE', 
		help = "file containing taxids for BLAST, one per line", default=None)
	parser.add_argument('-e', '--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for keeping hits in analysis (1e-10)", default=1e-10)
	parser.add_argument('-n', '--max_hits', action='store', type=int, dest='MAX_HITS', 
		help = "max number of hits to keep (100)", default=100)
	parser.add_argument('--db', action='store', type=str, dest='DB', 
		help = "Uniprot database to search (uniprotkb)", default="uniprotkb")
	parser.add_argument('-l', '--length', action='store', type=float, dest='LENGTH', 
		help = "minimum length of matches, as percent of input (.7)", default=0.7)
	parser.add_argument('--align-blast-sequence', action='store_false', dest='FULLSEQS', 
		help = "only align BLAST hit sequences", default=True)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients", default = "./scripts/")	

	parser.add_argument("--exclude-paralogs", action="store_true", dest="EXCLUDE_PARALOGS", 
		help="filter BLAST results, return only best match per non-subject species", 
		default=False)

	args, unknowns = parser.parse_known_args()

	parsed_taxids = []
	if args.TAX_FILE != None:
		parsed_taxids = [x for x in open(args.TAX_FILE, "r").readlines().strip()]

	main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME, args.OUTPUT,
		 args.MAX_HITS, args.EVALUE, args.DB, args.LENGTH, args.FULLSEQS,
		 args.TAXID, args.CLIENTS_FOLDER, args.EXCLUDE_PARALOGS)	

