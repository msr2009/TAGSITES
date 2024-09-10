import json



INPUT_JSON = json.load(open("./default_json.txt", "r"))
TASK_PARAMETERS = INPUT_JSON["analyses"]
AVAILABLE_TASKS = list(INPUT_JSON["analyses"].keys())
EXCLUDE_ARGS = ["fasta", "input", "output"]
