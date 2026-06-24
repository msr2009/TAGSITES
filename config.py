#file containing flattened list of uniprot species and taxids
UNIPROT_SPECIES = "./uniprot_species.flat.txt"

#edit the filename below for your defaults
DEFAULT_JSON = "./default_json.json"

#you can change these if you'd like (except for Other)
DEFAULT_SPECIES = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Rattus norvegicus": 10116,
    "Danio rerio": 7955,
    "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239,
    "Saccharomyces cerevisiae": 559292,
    "Escherichia coli": 562,
    "Other (search...)": None  # marker for dynamic search
}

#don't touch this!
import json
INPUT_JSON = json.load(open(DEFAULT_JSON, "r"))
TASK_PARAMETERS = INPUT_JSON["analyses"]
AVAILABLE_TASKS = list(INPUT_JSON["analyses"].keys())
EXCLUDE_ARGS = ["fasta", "input", "output", "pdb"]

#analysis result types
#add new analysis types here
RESULTS_TYPE_DICT = {
	"CONTINUOUS": ["blast", "plddt", "scores"], 
	"RANGE": ["domains", "modifications"]
}
