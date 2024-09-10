"""
get_species_taxonomy.py

A script to determine uniprot taxonomy of a given species
e.g., the host of a BLAST hit for a given sequence

Uses EBI's API (dbfetch). 

1) Recursively checks for "PARENT ID" fields in the results of a 
   dbfetch query of the taxonomy database, and stores the results
   in a dict.

Parameters:
	- taxid (string): taxonomic id of hit
	- taxlevel (string): level to output

Returns: 
	- by default, the taxonomy w/ taxids 
	- if taxlevel is defined, then just that taxid

Matt Rich, 9/2024
"""

import subprocess, re
from site_selection_util import get_sequence, uniprot_accession_regex, save_fasta

def dbfetch_call(taxid, email, clients_folder):
	dbfetch_call = "python {}/dbfetch.py fetchData \
					taxonomy:{}".format(clients_folder, taxid)
	
	#call dbfetch command
	tax_out = subprocess.run(dbfetch_call, shell=True, capture_output=True, text=True).stdout
	return tax_out

def tax_to_dict(tax_output):
		return {y[0].strip(): y[1].strip() for y in [x.split(":") for x in tax_output.strip("//\n\n").split("\n")]}
					
def main(taxid, taxlevel, email, clients_folder):

	#initialize with first taxid
	tax_dicts = [tax_to_dict(dbfetch_call(taxid, email, clients_folder))]
	while tax_dicts[-1]["PARENT ID"] != "1":
		parent = tax_dicts[-1]["PARENT ID"]
		tax_dicts.append(tax_to_dict(dbfetch_call(parent, email, clients_folder)))

	if taxlevel != None:
		for t in tax_dicts:
			if t["RANK"] == taxlevel:
				return t["ID"]
				break
	else:
		out = ""
		for t in tax_dicts:
			out += "\t".join([t["ID"], t["PARENT ID"], t["RANK"], t["SCIENTIFIC NAME"]])
			out += "\n"
		return out
	
if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-t', '--taxid', action='store', type=str, dest='TAXID', 
		help = "uniprot taxonomic ID to query", required=True)
	parser.add_argument('-l', '--taxlevel', action='store', type=str, dest='TAXLEVEL', 
		choices=["species", "genus", "subfamily", "family", "superfamily", "order", "class", "phylum", "kingdom", "superkingdom"],
		help = "taxonomic level to output", default=None)
	parser.add_argument('--email', action='store', type=str, dest='EMAIL', 
		help = "email address, required by EBI job submission.", required=True)
	parser.add_argument('--clients-folder', action='store', type=str, dest='CLIENTS_FOLDER', 
		help = "path to EBI webservice clients (default=./ebi_api_clients/)",
		default = "./ebi_api_clients/")

	args = parser.parse_args()
	
	output = main(args.TAXID, args.TAXLEVEL, args.EMAIL, args.CLIENTS_FOLDER)	
	print(output)
	
