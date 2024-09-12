#edit the filename below for your defaults
DEFAULT_JSON = "./default_json.txt"

#don't touch this!
import json
INPUT_JSON = json.load(open(DEFAULT_JSON, "r"))
TASK_PARAMETERS = INPUT_JSON["analyses"]
AVAILABLE_TASKS = list(INPUT_JSON["analyses"].keys())
EXCLUDE_ARGS = ["fasta", "input", "output"]
