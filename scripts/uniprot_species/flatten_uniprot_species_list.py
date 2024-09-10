"""
flatten_uniprot_species_list.py

script to download species list from uniprot, parse, 
and output a "flattened" dataset.

written with some help from chatgpt, but it was wrong in a few places...

Matt Rich, 07/2024
"""

import re, sys, requests

def fetch_file(url):
	response = requests.get(url)
	response.raise_for_status()  # Check if the request was successful
	return response.text.splitlines()

def parse_speclist(lines):
    
    species_data = []
    species_entry = {}
    parsing_species = False

    for line in lines:
        if not parsing_species:
            if line.startswith('_____'):
                parsing_species = True
            continue
            
        if parsing_species: 
            if re.match(r'^[A-Z]{2,5}', line):
                if species_entry:
                    #species_data.append(species_entry)
                    print_flattened_species(species_entry)
                    species_entry = {}
        
                parts = line.split()
                code = parts[0]
                classifier = parts[1]
                taxon_node = int(parts[2][:-1])
                scientific_name = ' '.join(parts[3:])[2:]
                species_entry = {
                    'code': code,
                    'classifier': classifier,
                    'taxon_node': taxon_node,
                    'scientific_name': scientific_name,
                    'common_name': '',
                    'synonyms': ''
                    }

            elif line.strip().startswith('C='):
                species_entry['common_name'] = line.strip()[2:].strip()

            elif line.strip().startswith('S='):
                species_entry['synonyms'] = line.strip()[2:].strip()
        
    return species_data

def print_flattened_species(entry, f_out=sys.stdout):
    out_line = "\t".join([entry["code"], entry["classifier"], str(entry['taxon_node']), 
                          entry['scientific_name'], entry['common_name'], entry['synonyms']])
    print(out_line, file=f_out)

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--url', action = 'store', type = str, dest = 'SPEC_URL', 
		help = "url for uniprot species list",
		default = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt")
	parser.add_argument('--file', action = 'store', type = str, dest = 'SPEC_FILE',
		help = "uniprot species list file. will use this instead of url.",
		default = None)
	args = parser.parse_args()
	
	spec_lines = []
	if args.SPEC_FILE != None:
		spec_lines = open(args.SPEC_FILE, "r").readlines()
	else:
		spec_lines = fetch_file(args.SPEC_URL)
	parse_speclist(spec_lines)

