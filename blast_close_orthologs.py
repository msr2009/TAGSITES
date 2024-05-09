"""
blast_close_orthologs.py



Matt Rich, 4/2024
"""

import subprocess, json
from site_selection_util import read_fasta


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

def main(fasta_in, email, workingdir, name, n, evalue, db, length_percent, align_full_seqs, clients_folder):

	name, seq = read_fasta(fasta_in)

	out_prefix = "{}/{}.blast_close".format(workingdir, name)

	############################
	# BLAST TO FIND ORTHOLOGS
	###########################

	#blast seq against database
	#make ncbiblast command
	ncbi_call = "python {}ncbiblast.py \
					--email {} \
					--program blastp \
					--stype protein \
					--sequence {} \
					--database {} \
					--outformat out,json,accs,ids \
					--alignments 0 --scores {} --exp {} \
					--outfile {} --pollFreq 10".format(clients_folder, email, seq, db, 5*n, evalue, out_prefix)

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
		if d_hit["evalue"] > evalue and d_hit["length"]/float(len(seq)) < length_percent:
			continue
		if d_hit["length"]/float(len(seq)) < length_percent or d_hit["length"]/float(len(seq)) > 1/length_percent:
			continue

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
				print(dbfetch_call)
				acc_fasta = subprocess.run(dbfetch_call, shell=True,capture_output=True, text=True).stdout
				fasta_str += acc_fasta.rstrip('\n')+"\n"
			else:
				#we can make a fasta str (and save a file) for the alignment
				fasta_str += ">{}\n{}\n".format(h["acc"], h["hitseq"])

	
	###########################
	# ALIGN WITH CLUSTALO
	###########################

	#output fasta of seqs to align
	fasta_out = open("{}/{}.blast_close.fasta".format(workingdir, name), "w")
	print(fasta_str, file=fasta_out)
	fasta_out.close()

	#call clustalomega to realign ortholog seqs
	# I can't figure out how to format these for submission to EBI, so 
	# I think it's better to simply run clustal locally. (below)
#	ebi_clustal_call = "python {}/clustalo.py --email {} --sequence {} --stype protein --outfile {}".format(clients_folder, email, fasta_str, out_prefix)
#	subprocess.run(ebi_clustal_call, shell=True)

	local_clustal_call = "clustalo -i {}/{}.blast_close.fasta -o {}/{}.blast_close.aln".format(workingdir, name, workingdir, name)
#	local_clustal_call = "clustalo -i {}/{}.blast_close.fasta".format(workingdir, name)
	print(local_clustal_call)
	subprocess.run(local_clustal_call, shell=True)
	
	###########################
	# CALCULATE JSD
	###########################

	#calculate jensen-shannon divergence (from Capra and Singh)
	#we need the name of the best hit, which will be our input seq, 
	#this is the first seq in the fasta file
	best_hit_name = open("{}/{}.blast_close.fasta".format(workingdir, name), "r").readline()[1:]

	js_call = "python {}/score_conservation_py3.py -i {}/{}.blast_close.aln -a {} -o {}/{}.blast_close.jsd".format(clients_folder, 
																											workingdir, name, 
																											best_hit_name,
																											workingdir, name)

	print(js_call)
	subprocess.run(js_call, shell=True)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()

	parser.add_argument('-f', '--fasta', action='store', type=str, dest='FASTA_IN', 
		help = "name of fasta file containing seq.", required=True)	

	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--dir', action='store', type=str, dest='WORKINGDIR', 
		help = "working directory for output", required=True)
	parser.add_argument('--name', action='store', type=str, dest='NAME', 
		help = "name for output", required=True)

	parser.add_argument('-e', '--evalue', action='store', type=float, dest='EVALUE', 
		help = "evalue threshold for keeping hits in analysis (1e-10)", default=1e-10)
	parser.add_argument('-n', '--maxhits', action='store', type=int, dest='MAX_HITS', 
		help = "max number of hits to keep (100)", default=100)
	parser.add_argument('--db', action='store', type=str, dest='DB', 
		help = "Uniprot database to search (uniprotkb)", default="uniprotkb")
	parser.add_argument('-l', '--length', action='store', type=float, dest='LENGTH', 
		help = "minimum length of matches, as percent of input (.7)", default=0.7)
	parser.add_argument('--align-full-sequence', action='store_true', dest='FULLSEQS', 
		help = "fetch and align full protein sequences (by default, uses BLAST alignment)")
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")	

	args = parser.parse_args()
	
	main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME, 
		 args.MAX_HITS, args.EVALUE, args.DB, args.LENGTH, args.FULLSEQS,
		 args.CLIENTS_FOLDER+"/")	

