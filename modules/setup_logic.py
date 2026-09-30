"""setup_logic.py — pure functions for building Setup data structures.

No Shiny imports; all functions accept/return plain Python dicts and are
independently testable without running the app.
"""


def make_task(task_type, label, params_cfg, tooltips_cfg, choices_cfg=None):
    """Build an in-memory task dict from user selection + config defaults."""
    return {
        "id":       f"{label}_{task_type}",
        "name":     task_type,
        "params":   dict(params_cfg),    # {param: current_value} — visible params only
        "tooltips": dict(tooltips_cfg),
        "choices":  dict(choices_cfg or {}),  # {param: [options]} for dropdown params
    }


def build_global_block(email, run_name, working_dir, input_file,
                       pdb="", genomic_file="", scripts_folder="./scripts/", taxid=""):
    """Assemble the 'global' block that goes into the run JSON."""
    g = {
        "email":       email,
        "run_name":    run_name,
        "working_dir": working_dir,
        "input_file":  input_file,
        "pdb":         pdb,
        "scripts_folder": scripts_folder,
        "selected_sites": [],
    }
    # only include genomic_file when a genomic FASTA was uploaded
    if genomic_file:
        g["genomic_file"] = genomic_file
    # species taxid, when known — progress_logic.parse_run() merges the global block
    # into every task's args, so this reaches any runner that wants it (the reagents
    # off-target screen scopes its blastn by it). A task-level taxid still wins.
    if taxid:
        g["taxid"] = str(taxid)
    return g


def build_task_entry(task_id, task_name, collected_args, out_suffix, working_dir, run_name):
    """Build a single task JSON entry — task args only, no global keys.

    Global fields (email, working_dir, etc.) are stored once in the run JSON's
    global block and merged at run time, not copied into every task.
    """
    args = dict(collected_args)
    args["output"] = f"{working_dir}{run_name}.{task_id}{out_suffix}"
    return {task_id: {"type": task_name, "args": args}}


def build_defaults_entry(task_id, task_name, collected_args, out_suffix, working_dir, run_name):
    """Build one entry for a saved-defaults file — task args only, no global block."""
    args = dict(collected_args)
    args["output"] = f"{working_dir}{run_name}.{task_id}{out_suffix}"
    return {task_id: {"type": task_name, "args": args}}


def build_reagents_entry(defaults, genomic_path, out_suffix, working_dir, run_name,
                         offtarget=None, offtarget_spacer=None):
    """Build the auto-injected REAGENTS_reagents entry for CRISPR precompute.

    Global fields are omitted — they're merged at run time from the global block.
    offtarget=False skips the EBI BLAST off-target jobs, the slowest part of a run;
    offtarget_spacer=False keeps the region screen (primer co-amplification) but
    skips the slower per-guide search. None leaves task_definitions.json's defaults
    in place.
    """
    args = dict(defaults)
    args["genomic_fasta"] = genomic_path
    args["genewise"] = ""          # runner calls Genewise automatically
    args["output"] = f"{working_dir}{run_name}.reagents{out_suffix}"
    if offtarget is not None:
        args["offtarget"] = bool(offtarget)
    if offtarget_spacer is not None:
        args["offtarget_spacer"] = bool(offtarget_spacer)
    return {"type": "reagents", "args": args}


def build_run_json(global_block, task_entries, reagents_entry=None):
    """Assemble the full top-level run JSON dict from its parts."""
    tasks = {}
    for entry in task_entries:
        tasks.update(entry)
    # only include reagents when a genomic FASTA was uploaded
    if reagents_entry is not None:
        tasks["REAGENTS_reagents"] = reagents_entry
    return {"global": global_block, "tasks": tasks}
