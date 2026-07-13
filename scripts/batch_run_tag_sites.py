"""
batch_run_tag_sites.py

CLI batch driver: runs run_tag_sites_from_json.main() for many proteins at
once. Reads a manifest (one row per protein) and a shared task-template JSON
(same shape as params/worm_default.json), builds a per-protein run JSON via
modules/setup_logic.py, then dispatches them with bounded concurrency so we
don't overwhelm EBI's REST API with simultaneous submissions.

Matt Rich, 2026
"""

import csv
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

_SCRIPTS = Path(__file__).parent
_ROOT = _SCRIPTS.parent
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_ROOT / "modules"))

import run_status
import run_tag_sites_from_json
from task_registry import task_defaults, task_output_suffix
from setup_logic import build_global_block, build_task_entry, build_reagents_entry, build_run_json


def _sniff_delimiter(manifest_path):
    """Detect comma vs tab delimiter from the manifest's header line."""
    with open(manifest_path) as f:
        header = f.readline()
    return "\t" if "\t" in header else ","


def read_manifest(manifest_path):
    """Read the batch manifest into a list of row dicts, keyed by column name."""
    delimiter = _sniff_delimiter(manifest_path)
    with open(manifest_path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def build_protein_run_json(row, template, default_email, output_root):
    """Build the full run JSON dict + its destination path for one manifest row."""
    run_name = row["run_name"]
    working_dir = row.get("working_dir") or os.path.join(output_root, run_name)
    if not working_dir.endswith("/"):
        working_dir += "/"
    os.makedirs(working_dir, exist_ok=True)

    global_block = build_global_block(
        email=row.get("email") or default_email,
        run_name=run_name,
        working_dir=working_dir,
        input_file=row["input_file"],
        pdb=row.get("pdb", ""),
        genomic_file=row.get("genomic_file", ""),
    )

    task_entries = []
    reagents_template = None
    for task_id, entry in template.items():
        if entry["type"] == "reagents":
            reagents_template = entry
            continue
        task_entries.append(
            build_task_entry(
                task_id,
                entry["type"],
                entry["args"],
                task_output_suffix(entry["type"]),
                working_dir,
                run_name,
            )
        )

    reagents_entry = None
    genomic_file = row.get("genomic_file", "")
    if genomic_file:
        defaults = reagents_template["args"] if reagents_template else task_defaults("reagents")
        reagents_entry = build_reagents_entry(
            defaults=defaults,
            genomic_path=genomic_file,
            out_suffix=task_output_suffix("reagents"),
            working_dir=working_dir,
            run_name=run_name,
        )
    elif reagents_template:
        print(
            f"WARNING: skipping reagents task for {run_name}: no genomic_file in manifest",
            file=sys.stderr,
        )

    run_json = build_run_json(global_block, task_entries, reagents_entry)
    run_json_path = os.path.join(working_dir, f"{run_name}.run.json")
    return run_json_path, run_json


def run_one(run_json_path, force):
    """Run one protein's pipeline and return (run_name, ok, message)."""
    run_name = Path(run_json_path).name.removesuffix(".run.json")
    try:
        run_tag_sites_from_json.main(run_json_path, force=force)
    except Exception as e:
        return run_name, False, f"exception: {e}"

    with open(run_json_path) as f:
        working_dir = json.load(f)["global"]["working_dir"]
    status = run_status.load_status(working_dir, run_name)
    failed = [tid for tid, t in status.items() if t.get("status") != "success"]
    if failed:
        return run_name, False, f"incomplete tasks: {', '.join(failed)}"
    return run_name, True, "ok"


def main(manifest_path, params_path, email, output_root, max_concurrent=3, force=False):
    """Build and run one pipeline per manifest row, capped at max_concurrent concurrent proteins."""
    with open(params_path) as f:
        template = json.load(f)
    rows = read_manifest(manifest_path)
    if not rows:
        print("Manifest is empty; nothing to do.", file=sys.stderr)
        return

    run_json_paths = []
    for row in rows:
        run_json_path, run_json = build_protein_run_json(row, template, email, output_root)
        with open(run_json_path, "w") as f:
            json.dump(run_json, f, indent=4)
        run_json_paths.append(run_json_path)

    results = []
    with ThreadPoolExecutor(max_workers=max_concurrent) as pool:
        futures = {pool.submit(run_one, p, force): p for p in run_json_paths}
        for future in as_completed(futures):
            results.append(future.result())

    print("\n#############\nBATCH SUMMARY\n#############")
    n_failed = 0
    for run_name, ok, message in sorted(results, key=lambda r: r[0]):
        status_word = "OK" if ok else "FAILED"
        print(f"{run_name}: {status_word} ({message})")
        n_failed += not ok

    if n_failed:
        print(f"\n{n_failed}/{len(results)} protein(s) failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "--manifest",
        required=True,
        help="CSV/TSV file, one row per protein (run_name, input_file, "
        "and optional pdb, genomic_file, working_dir, email columns)",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="task-template JSON shared across the batch (same shape as params/worm_default.json)",
    )
    parser.add_argument(
        "--email", default="", help="default EBI e-mail; overridden per-row by an 'email' column"
    )
    parser.add_argument(
        "--output-root",
        default="./batch_runs",
        help="base directory for per-protein working dirs when a row "
        "omits 'working_dir' (default: ./batch_runs)",
    )
    parser.add_argument(
        "--max-concurrent",
        type=int,
        default=3,
        help="max proteins to run simultaneously (default: 3, to avoid "
        "overwhelming EBI's REST API)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="re-run all tasks even if already complete",
    )
    args = parser.parse_args()

    main(
        args.manifest,
        args.params,
        args.email,
        args.output_root,
        max_concurrent=args.max_concurrent,
        force=args.force,
    )
