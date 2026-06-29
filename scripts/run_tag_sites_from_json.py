"""
run_tag_sites_from_json.py

CLI orchestrator: reads a run JSON and submits each analysis as a
concurrent subprocess.  Idempotent: completed tasks (status == success AND
output file exists) are skipped; failed/pending tasks are re-run.
Add --force to re-run everything regardless of prior status.

Matt Rich, 09/24 / updated 2026
"""

import subprocess, os, json, sys, threading
from pathlib import Path

# ensure scripts/ is on sys.path for sibling imports
_SCRIPTS = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS))
import run_status
import progress
from task_registry import task_script, task_hidden, GLOBAL_KEYS


def _stderr_reader(proc, task_id, working_dir, run_name):
    """Drain proc.stderr line-by-line, parsing [stage] events into the status file."""
    log_lines = []
    stage = ""
    for line in proc.stderr:
        s, msg, level = progress.parse_line(line)
        if s:
            stage = s
        log_lines.append(line.rstrip("\n"))
        run_status.update_task(working_dir, run_name, task_id,
                               stage=stage, log="\n".join(log_lines))


def genewise_required(tasks):
    """Return the key of a 'reagents' task that still needs Genewise, or None."""
    for key, v in tasks.items():
        if v["type"] == "reagents" and v["args"].get("genewise", "") == "":
            return key
    return None


def searchAFDB_required(tasks, global_args):
    """Return AFDB search args string if there's a plddt task with no pdb path, else None."""
    # existing_AF_model.py accepts: fasta/input_file, email, working_dir, run_name,
    # taxid, evalue, percent_id — exclude everything else
    _AFDB_EXCLUDE = {"pdb", "output", "existing_AF2", "scripts_folder",
                     "selected_sites", "genomic_file", "taxid_file"}
    for key, v in tasks.items():
        if (v["type"] == "plddt" and v["args"].get("pdb", "") == ""
                and v["args"].get("existing_AF2", 1) != 0):
            merged = {**global_args, **v["args"]}
            merged.setdefault("taxid", merged.get("species_taxid", 1))
            return build_task_args_string(merged, EXCLUDE=list(_AFDB_EXCLUDE))
    return None


def build_task_args_string(targs, EXCLUDE=None):
    """Build a CLI argument string from a merged args dict, skipping empty and excluded keys."""
    if EXCLUDE is None:
        EXCLUDE = []
    return "".join(
        f" --{arg} {targs[arg]}"
        for arg in targs
        if targs[arg] != "" and arg not in EXCLUDE
    )


def _write_json(path, global_args, tasks):
    """Persist the current run state back to disk."""
    with open(path, "w") as f:
        json.dump({"global": global_args, "tasks": tasks}, f, indent=4)


def main(json_input_file, force=False):
    """Run all tasks in the JSON; skip completed tasks unless force=True."""
    json_in = {}
    try:
        with open(json_input_file) as f:
            json_in = json.load(f)
    except json.decoder.JSONDecodeError:
        raise IOError("Incorrectly formatted JSON file!")

    global_args = json_in["global"]
    tasks       = json_in["tasks"]

    # use absolute path to scripts/ so this script is runnable from any cwd
    scripts_folder = str(_SCRIPTS) + "/"
    interpreter = "python"

    working_dir = global_args.get("working_dir", ".")
    run_name    = global_args.get("run_name", "run")

    # load existing status file (empty dict if first run)
    status = run_status.load_status(working_dir, run_name)

    ###############################################
    # PRE-STEP: Genewise (for reagents tasks)
    ###############################################
    reagents_task_key = genewise_required(tasks)
    if reagents_task_key is not None:
        task_args   = tasks[reagents_task_key]["args"]
        genomic_fa  = task_args.get("genomic_fasta", "")
        gw_out_prefix = os.path.join(working_dir, run_name + "_genewise")
        gw_out_file   = gw_out_prefix + ".genewise.out.txt"

        # skip genewise pre-step if output already exists
        if not force and os.path.exists(gw_out_file):
            print(f"Skipping Genewise pre-step: {gw_out_file} already exists")
            tasks[reagents_task_key]["args"]["genewise"]      = gw_out_file
            tasks[reagents_task_key]["args"]["genomic_fasta"] = gw_out_prefix + ".genewise_genomic.fa"
        elif not genomic_fa:
            print("ERROR: reagents task has empty 'genewise' and 'genomic_fasta' — cannot run Genewise.",
                  file=sys.stderr)
        else:
            protein_fa  = global_args.get("input_file", "")
            email       = global_args.get("email", "")
            gw_call = (
                f"{interpreter} {scripts_folder}run_genewise.py "
                f"--protein_fasta {protein_fa} --genomic_fasta {genomic_fa} "
                f"--email {email} --outprefix {gw_out_prefix}"
            )
            print("RUNNING GENEWISE (both orientations)\n\n" + gw_call)
            ret = subprocess.call(gw_call, shell=True)
            if ret == 0:
                tasks[reagents_task_key]["args"]["genewise"]      = gw_out_file
                tasks[reagents_task_key]["args"]["genomic_fasta"] = gw_out_prefix + ".genewise_genomic.fa"
                _write_json(json_input_file, global_args, tasks)
            else:
                print(f"WARNING: run_genewise.py exited {ret}; reagents task may fail.",
                      file=sys.stderr)

    ###############################################
    # PRE-STEP: AFDB search (for plddt tasks)
    ###############################################
    search_afdb = searchAFDB_required(tasks, global_args)
    if search_afdb is not None:
        af2_pdb = os.path.join(working_dir, run_name + ".AF.pdb")
        if not force and os.path.exists(af2_pdb):
            print(f"Skipping AFDB search: {af2_pdb} already exists")
        else:
            af2_search_call = f"{interpreter} {scripts_folder}existing_AF_model.py {search_afdb}"
            print("SEARCHING AFDB\n\n" + af2_search_call)
            af_proc = subprocess.call(af2_search_call, shell=True)

            if af_proc != 1:
                # AFDB found a model — rename files and update plddt task args
                orig_fa = os.path.join(working_dir, run_name + ".fa")
                user_fa = os.path.join(working_dir, "user_" + run_name + ".fa")
                af_fa   = os.path.join(working_dir, run_name + ".AF.fa")
                if os.path.exists(orig_fa):
                    os.rename(orig_fa, user_fa)
                if os.path.exists(af_fa):
                    os.rename(af_fa, orig_fa)

                # record original input_file in global before overwriting
                global_args["user_input_file"] = global_args.get("input_file", "")
                global_args["input_file"] = orig_fa

                for key, v in tasks.items():
                    if v["type"] == "plddt":
                        v["args"]["existing_AF2"] = 0
                        v["args"]["pdb"] = af2_pdb

                orig_json_path = os.path.join(working_dir, run_name + ".json")
                user_json_path = os.path.join(working_dir, "user_" + run_name + ".json")
                if os.path.exists(orig_json_path):
                    os.rename(orig_json_path, user_json_path)
                _write_json(json_input_file, global_args, tasks)

    ###############################################
    # MAIN TASK LOOP — re-read JSON (may have been rewritten above)
    ###############################################
    with open(json_input_file) as f:
        json_in = json.load(f)
    global_args = json_in["global"]
    tasks       = json_in["tasks"]

    # build the EXCLUDE list: global keys + hidden params from registry + output
    hidden_params = {p for ttype in tasks.values()
                       for p in task_hidden(ttype["type"])}
    EXCLUDE_CLI = list(GLOBAL_KEYS | hidden_params | {"output"})

    task_commands = []
    for task_id, v in tasks.items():
        output_path = v["args"].get("output", "")
        task_entry  = status.get(task_id, {})

        if not force and run_status.is_complete(task_entry, output_path):
            print(f"Skipping {task_id}: already complete ({output_path})")
            continue

        # merge global into args before building the CLI string
        merged_args = {**global_args, **v["args"]}
        script      = f"{interpreter} {scripts_folder}" + task_script(v["type"])
        arguments   = build_task_args_string(merged_args, EXCLUDE=EXCLUDE_CLI)
        task_commands.append((task_id, script + arguments, output_path))

        # mark as pending before launch so a crash shows "pending" not stale "success"
        run_status.update_task(working_dir, run_name, task_id,
                               status="pending", output=output_path, job_id="", message="")

    if not task_commands:
        print("All tasks already complete. Use --force to re-run.")
        return

    print("#############\n#############\n#############")
    for _, cmd, _ in task_commands:
        print(cmd)
    print("#############\n#############\n#############")

    # launch all remaining tasks concurrently; reader threads drain stderr live
    procs = []
    for task_id, cmd, output_path in task_commands:
        run_status.update_task(working_dir, run_name, task_id, status="running", stage="", log="")
        proc = subprocess.Popen(cmd, shell=True,
                                stderr=subprocess.PIPE, text=True, bufsize=1)
        reader = threading.Thread(target=_stderr_reader,
                                  args=(proc, task_id, working_dir, run_name),
                                  daemon=True)
        reader.start()
        procs.append((task_id, proc, output_path, reader))

    for task_id, proc, output_path, reader in procs:
        ret = proc.wait()
        reader.join()   # drain any remaining stderr before writing final status
        if ret == 0 and (not output_path or os.path.exists(output_path)):
            run_status.update_task(working_dir, run_name, task_id, status="success", stage="done")
        else:
            msg = f"exit code {ret}" if ret != 0 else "output file not found"
            run_status.update_task(working_dir, run_name, task_id,
                                   status="failed", message=msg, stage="failed")

    print("####DONE WITH ANALYSIS####")


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-i", "--input", "--json", action="store", type=str,
                        dest="JSON_INPUT", help="path to run JSON file")
    parser.add_argument("--force", action="store_true", dest="FORCE",
                        help="re-run all tasks even if already complete", default=False)
    args = parser.parse_args()

    main(args.JSON_INPUT, force=args.FORCE)
