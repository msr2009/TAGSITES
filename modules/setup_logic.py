"""setup_logic.py — pure functions for building Setup data structures.

No Shiny imports; all functions accept/return plain Python dicts and are
independently testable without running the app.
"""


def make_task(task_type, label, params_cfg, tooltips_cfg):
    """Build an in-memory task dict from user selection + config defaults."""
    return {
        "id": f"{label}_{task_type}",
        "name": task_type,
        "params": dict(params_cfg),
        "tooltips": dict(tooltips_cfg),
    }


def build_global_block(base_global, email, run_name, working_dir,
                       input_file, pdb="", genomic_file=""):
    """Assemble the 'global' block that goes into the run JSON."""
    g = dict(base_global)
    g["email"] = email
    g["run_name"] = run_name
    g["working_dir"] = working_dir
    g["input_file"] = input_file
    g["pdb"] = pdb
    # only include genomic_file when a genomic FASTA was uploaded
    if genomic_file:
        g["genomic_file"] = genomic_file
    return g


def build_task_entry(task_id, task_name, collected_args, global_block, out_suffix, working_dir, run_name):
    """Merge task args over global block and return a single-task JSON entry.

    Global block keys override task args (matching original write_task_json behavior),
    so shared fields like email/working_dir/input_file always come from the run context.
    Output path is always rewritten from working_dir + run_name + task_id + suffix.
    """
    args = {**collected_args, **global_block}
    args["output"] = f"{working_dir}{run_name}.{task_id}{out_suffix}"
    return {task_id: {"analysis": task_name, "args": args}}


def build_defaults_entry(task_id, task_name, collected_args, out_suffix, working_dir, run_name):
    """Build one entry for a saved-defaults file — task args only, no global block."""
    args = dict(collected_args)
    args["output"] = f"{working_dir}{run_name}.{task_id}{out_suffix}"
    return {task_id: {"analysis": task_name, "args": args}}


def build_reagents_entry(defaults_cfg, genomic_path, email, working_dir, run_name, input_file):
    """Build the auto-injected REAGENTS_reagents entry for CRISPR precompute."""
    args = dict(defaults_cfg)
    args["genomic_fasta"] = genomic_path
    args["genewise"] = ""           # runner calls Genewise automatically
    args["email"] = email
    args["working_dir"] = working_dir
    args["run_name"] = run_name
    args["input_file"] = input_file
    out_suffix = defaults_cfg.get("output", ".tsv")
    args["output"] = f"{working_dir}{run_name}.reagents{out_suffix}"
    return {"analysis": "reagents", "args": args}


def build_run_json(scripts, global_block, task_entries, reagents_entry=None):
    """Assemble the full top-level run JSON dict from its parts."""
    result = {"scripts": scripts, "global": global_block}
    for entry in task_entries:
        result.update(entry)
    # only include reagents when a genomic FASTA was uploaded
    if reagents_entry is not None:
        result["REAGENTS_reagents"] = reagents_entry
    return result
