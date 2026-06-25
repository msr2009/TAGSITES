"""
run_tag_sites_from_json.py

CLI orchestrator: reads a run JSON and submits each analysis as a
concurrent subprocess.  Idempotent: completed tasks (status == success AND
output file exists) are skipped; failed/pending tasks are re-run.
Add --force to re-run everything regardless of prior status.

Matt Rich, 09/24 / updated 2026
"""

import subprocess, os, json, sys
from pathlib import Path

# ensure scripts/ is on sys.path for sibling imports
_SCRIPTS = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS))
import run_status


def genewise_required(j, global_args, task_scripts, scripts_folder, interpreter):
    """Detect a 'reagents' task that still needs Genewise (empty 'genewise' arg).

    Runs run_genewise.py synchronously for both strand orientations,
    picks the winner, then rewrites the JSON with resolved paths.
    """
    for key in j:
        if j[key]["analysis"] == "reagents":
            if j[key]["args"].get("genewise", "") == "":
                return key
    return None


def searchAFDB_required(j):
    """Return AFDB search args if there's a plddt task with no pdb path."""
    for key in j:
        if j[key]["analysis"] == "plddt":
            if j[key]["args"]["pdb"] == "":
                try:
                    j[key]["args"]["taxid"] = j[key]["args"]["species_taxid"]
                except KeyError:
                    j[key]["args"]["taxid"] = 1
                return build_task_args_string(j[key]["args"], ["pdb", "output", "existing_AF2"])
    return None


def build_task_args_string(targs, EXCLUDE=[]):
    """Build a CLI argument string from a task args dict."""
    return "".join(
        [f" --{arg} {targs[arg]}" for arg in targs if targs[arg] != "" and arg not in EXCLUDE]
    )


def main(json_input_file, force=False):
    """Run all tasks in the JSON; skip completed tasks unless force=True."""
    json_in = {}
    try:
        with open(json_input_file) as f:
            json_in = json.load(f)
    except json.decoder.JSONDecodeError:
        raise IOError("Incorrectly formatted JSON file!")

    task_scripts = json_in.pop("scripts", None)
    global_args  = json_in.pop("global", None)
    # json_in is now just task entries

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
    reagents_task_key = genewise_required(json_in, global_args, task_scripts, scripts_folder, interpreter)
    if reagents_task_key is not None:
        task_args   = json_in[reagents_task_key]["args"]
        genomic_fa  = task_args.get("genomic_fasta", "")
        gw_out_prefix = os.path.join(working_dir, run_name + "_genewise")
        gw_out_file   = gw_out_prefix + ".genewise.out.txt"

        # skip genewise pre-step if output already exists
        if not force and os.path.exists(gw_out_file):
            print(f"Skipping Genewise pre-step: {gw_out_file} already exists")
            json_in[reagents_task_key]["args"]["genewise"]      = gw_out_file
            json_in[reagents_task_key]["args"]["genomic_fasta"] = gw_out_prefix + ".genewise_genomic.fa"
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
                json_in[reagents_task_key]["args"]["genewise"]      = gw_out_file
                json_in[reagents_task_key]["args"]["genomic_fasta"] = gw_out_prefix + ".genewise_genomic.fa"
                # persist resolved paths
                new_json = {**json_in, "global": global_args, "scripts": task_scripts}
                with open(json_input_file, "w") as f:
                    json.dump(new_json, f, indent=4)
            else:
                print(f"WARNING: run_genewise.py exited {ret}; reagents task may fail.",
                      file=sys.stderr)

    ###############################################
    # PRE-STEP: AFDB search (for plddt tasks)
    ###############################################
    search_afdb = searchAFDB_required(json_in)
    if search_afdb is not None:
        af2_pdb = os.path.join(working_dir, run_name + ".AF.pdb")
        if not force and os.path.exists(af2_pdb):
            print(f"Skipping AFDB search: {af2_pdb} already exists")
        else:
            af2_search_call = f"{interpreter} {scripts_folder}existing_AF_model.py {search_afdb}"
            print("SEARCHING AFDB\n\n" + af2_search_call)
            af_proc = subprocess.call(af2_search_call, shell=True)

            if af_proc != 1:
                # AFDB found a model — rename files and update task args
                orig_fa = os.path.join(working_dir, run_name + ".fa")
                user_fa = os.path.join(working_dir, "user_" + run_name + ".fa")
                af_fa   = os.path.join(working_dir, run_name + ".AF.fa")
                if os.path.exists(orig_fa):
                    os.rename(orig_fa, user_fa)
                if os.path.exists(af_fa):
                    os.rename(af_fa, orig_fa)

                for task in json_in:
                    json_in[task]["old_args"] = {}
                    json_in[task]["old_args"]["user_input_file"] = "user" + json_in[task]["args"]["input_file"]
                    if json_in[task]["analysis"] == "plddt":
                        json_in[task]["args"]["existing_AF2"] = 0
                        json_in[task]["args"]["pdb"] = af2_pdb

                orig_json_path = os.path.join(working_dir, run_name + ".json")
                user_json_path = os.path.join(working_dir, "user_" + run_name + ".json")
                if os.path.exists(orig_json_path):
                    os.rename(orig_json_path, user_json_path)
                new_json = {**json_in, "global": global_args, "scripts": task_scripts}
                with open(json_input_file, "w") as f:
                    json.dump(new_json, f, indent=4)

    ###############################################
    # MAIN TASK LOOP — re-read JSON (may have been rewritten above)
    ###############################################
    with open(json_input_file) as f:
        json_in = json.load(f)

    task_commands = []
    for task in json_in:
        if task in ("scripts", "global"):
            continue

        output_path = json_in[task]["args"].get("output", "")
        task_entry  = status.get(task, {})

        if not force and run_status.is_complete(task_entry, output_path):
            print(f"Skipping {task}: already complete ({output_path})")
            continue

        script    = f"{interpreter} {scripts_folder}" + json_in["scripts"][json_in[task]["analysis"]]
        arguments = build_task_args_string(json_in[task]["args"])
        task_commands.append((task, script + arguments, output_path))

        # mark as pending before launch so a crash shows "pending" not stale "success"
        run_status.update_task(working_dir, run_name, task,
                               status="pending", output=output_path, job_id="", message="")

    if not task_commands:
        print("All tasks already complete. Use --force to re-run.")
        return

    print("#############\n#############\n#############")
    for _, cmd, _ in task_commands:
        print(cmd)
    print("#############\n#############\n#############")

    # launch all remaining tasks concurrently, then collect exit codes
    procs = []
    for task_id, cmd, output_path in task_commands:
        run_status.update_task(working_dir, run_name, task_id, status="running")
        procs.append((task_id, subprocess.Popen(cmd, shell=True), output_path))

    for task_id, proc, output_path in procs:
        ret = proc.wait()
        if ret == 0 and (not output_path or os.path.exists(output_path)):
            run_status.update_task(working_dir, run_name, task_id, status="success")
        else:
            msg = f"exit code {ret}" if ret != 0 else "output file not found"
            run_status.update_task(working_dir, run_name, task_id,
                                   status="failed", message=msg)

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
