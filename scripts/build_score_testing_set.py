"""
build_score_testing_set.py

One-time (repeatable) consolidation: collect completed runs scattered across
batch folders into a single data/SCORE_TESTING/genes/ tree with a fresh
manifest (run_json, site, gene, allele), so scripts/batch_rescore_tagsites.py
has one stable dataset to iterate scoring schemes against without needing to
re-run the analysis pipeline.

Only the files needed for re-scoring are copied (run JSON, input FASTA,
AF2/pLDDT + companions, BLAST .jsd, domains/modifications/UniProt tracks,
reagents TSV) — never regenerable intermediates (BLAST .aln/.pdf/.json,
genewise fwd/rc, interpro raw TSV). The set of files to copy is derived from
the run JSON itself plus task_registry's companion map, not hardcoded
suffixes, so it stays correct if task_definitions.json changes.

Paths embedded in each run's <run>.run.json are rewritten from the old
working_dir to the new one; shared/repo-relative paths (./scripts/,
./tables/, data/internal_batch/inputs/, ...) are left untouched since they're
outside the run's own working_dir and aren't copied.

Matt Rich, 2026
"""

import argparse
import copy
import json
import os
import shutil
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).parent
_ROOT = _SCRIPTS.parent
sys.path.insert(0, str(_SCRIPTS))

from batch_run_tag_sites import read_manifest
from task_registry import task_companions, companion_path


def parse_run_name(run_name):
    """INTERNAL_<gene>_<allele>_<token> -> (gene, allele, token), or (None, None, None) if malformed."""
    parts = run_name.split("_")
    if len(parts) != 4 or parts[0] != "INTERNAL":
        return None, None, None
    return parts[1], parts[2], parts[3]


def load_sites(sites_manifest_paths):
    """Merge one or more run_json/site manifests into {run_name: site}.

    Later paths override earlier ones, so a small supplementary manifest (e.g.
    for a run added after the main batch, like syd-2) can be layered on top of
    an older batch's plot_manifest.tsv without editing it.
    """
    sites = {}
    for path in sites_manifest_paths:
        for row in read_manifest(path):
            run_name = Path(row["run_json"]).name.removesuffix(".run.json")
            sites[run_name] = row["site"]
    return sites


def is_run_complete(run_json_path):
    """True if run_json_path loads and every task's declared output file exists on disk."""
    try:
        with open(run_json_path) as f:
            run_json = json.load(f)
    except (OSError, json.JSONDecodeError):
        return False

    input_file = run_json.get("global", {}).get("input_file", "")
    if input_file and not os.path.exists(input_file):
        return False

    for entry in run_json.get("tasks", {}).values():
        output = entry.get("args", {}).get("output", "")
        if output and not os.path.exists(output):
            return False
    return True


def resolve_run_dir(run_name, source_roots):
    """Return the path to the first complete <run_name>.run.json found across source_roots, or None.

    Checked in the given root order; a root's copy only counts if every task
    output it declares actually exists (guards against stray partial dirs).
    """
    for root in source_roots:
        candidate = Path(root) / run_name / f"{run_name}.run.json"
        if candidate.exists() and is_run_complete(candidate):
            return candidate
    return None


def files_for_run(run_json_path):
    """Return (set of on-disk paths to copy, parsed run_json dict) for one run.

    Paths are read directly out of the run JSON (not assumed from run_json_path's
    location) — a run's declared working_dir can, in principle, point elsewhere
    than the directory the .run.json itself was found in.
    """
    with open(run_json_path) as f:
        run_json = json.load(f)

    paths = {str(run_json_path)}
    g = run_json.get("global", {})
    for key in ("input_file", "pdb"):
        v = g.get(key)
        if v:
            paths.add(v)

    for entry in run_json.get("tasks", {}).values():
        ttype = entry["type"]
        args = entry.get("args", {})
        output = args.get("output", "")
        if output:
            paths.add(output)
            # only pLDDT's companions (sasa/hydrophobic/patches) feed scoring or
            # plotting; blast's (alignment/isoforms) are alignment-view-only and
            # regenerable, so skip them here even though task_companions() has them
            if ttype == "plddt":
                for companion_key in task_companions(ttype):
                    cp = companion_path(output, ttype, companion_key)
                    if cp and os.path.exists(cp):
                        paths.add(cp)
        # extra path-valued args that live inside the run's own working_dir
        if ttype == "plddt" and args.get("pdb"):
            paths.add(args["pdb"])
        if ttype == "reagents":
            for k in ("genomic_fasta", "genewise"):
                if args.get(k):
                    paths.add(args[k])

    return paths, run_json


def rewrite_run_json(run_json, old_wd, new_wd):
    """Return a deep copy of run_json with every string value prefixed by old_wd rewritten to new_wd.

    Values outside old_wd (./scripts/, ./tables/, data/internal_batch/inputs/,
    ...) are left exactly as-is — they're shared/repo-relative, not part of
    this run's own working_dir, and are never copied.
    """
    def _rewrite(v):
        if isinstance(v, str) and old_wd and v.startswith(old_wd):
            return new_wd + v[len(old_wd):]
        return v

    out = copy.deepcopy(run_json)
    g = out.setdefault("global", {})
    for k in list(g.keys()):
        g[k] = _rewrite(g[k])
    g["working_dir"] = new_wd

    for entry in out.get("tasks", {}).values():
        args = entry.get("args", {})
        for k in list(args.keys()):
            args[k] = _rewrite(args[k])

    return out


def main(manifest_path, sites_manifest_path, source_roots, dest_root, out_manifest_path):
    """Consolidate every run named in manifest_path into dest_root/genes/, writing out_manifest_path."""
    rows = read_manifest(manifest_path)
    sites = load_sites(sites_manifest_path)

    dest_genes_rel = f"{dest_root.rstrip('/')}/genes"
    dest_genes = Path(dest_genes_rel)
    dest_genes.mkdir(parents=True, exist_ok=True)

    out_rows = []
    for row in rows:
        run_name = row["run_name"]
        gene, allele, _token = parse_run_name(run_name)
        if gene is None:
            print(f"WARNING: {run_name} doesn't parse as INTERNAL_<gene>_<allele>_<token>; skipping",
                  file=sys.stderr)
            continue

        site = sites.get(run_name)
        if site is None:
            print(f"WARNING: no site found for {run_name} in {sites_manifest_path}; skipping",
                  file=sys.stderr)
            continue

        src_run_json = resolve_run_dir(run_name, source_roots)
        if src_run_json is None:
            print(f"WARNING: no complete run found for {run_name} in {source_roots}; skipping",
                  file=sys.stderr)
            continue

        paths, run_json = files_for_run(src_run_json)
        old_wd = run_json.get("global", {}).get("working_dir", "")
        new_wd = f"{dest_genes_rel}/{run_name}/"
        dest_dir = dest_genes / run_name
        dest_dir.mkdir(parents=True, exist_ok=True)

        copied, skipped = 0, 0
        for p in paths:
            if not p.startswith(old_wd):
                # shared/repo-relative input outside this run's own working_dir; leave in place
                continue
            rel = p[len(old_wd):]
            dest_path = dest_dir / rel
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            if os.path.exists(p):
                shutil.copy2(p, dest_path)
                copied += 1
            else:
                print(f"WARNING: {run_name}: expected file missing, not copied: {p}", file=sys.stderr)
                skipped += 1

        new_run_json = rewrite_run_json(run_json, old_wd, new_wd)
        new_run_json_path = dest_dir / f"{run_name}.run.json"
        with open(new_run_json_path, "w") as f:
            json.dump(new_run_json, f, indent=4)

        out_rows.append({
            "run_json": str(new_run_json_path),
            "site": site,
            "gene": gene,
            "allele": allele,
        })
        print(f"{run_name}: copied {copied} file(s)" + (f", {skipped} missing" if skipped else ""))

    with open(out_manifest_path, "w", newline="") as f:
        f.write("run_json\tsite\tgene\tallele\n")
        for r in out_rows:
            f.write(f"{r['run_json']}\t{r['site']}\t{r['gene']}\t{r['allele']}\n")

    print(f"\nWrote {len(out_rows)}/{len(rows)} run(s) to {dest_genes_rel}/")
    print(f"Wrote manifest: {out_manifest_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="data/20260728_batch_manifest.tsv",
        help="manifest of runs to consolidate (column: run_name)")
    parser.add_argument("--sites-manifest", nargs="+",
        default=["data/20260727_batch/plot_manifest.tsv", "data/20260728_batch/plot_manifest.tsv"],
        help="one or more manifests providing each run's tagged site (columns: "
             "run_json, site); later manifests override earlier ones for the same run")
    parser.add_argument("--source-roots", nargs="+",
        default=["data/20260727_batch", "data/20260728_batch"],
        help="batch folders to search for each run, in preference order")
    parser.add_argument("--dest-root", default="data/SCORE_TESTING",
        help="destination root; runs are written to <dest-root>/genes/<run_name>/")
    parser.add_argument("--output-manifest", default="data/SCORE_TESTING/score_testing_manifest.tsv",
        help="path to write the fresh consolidated manifest")
    args = parser.parse_args()

    main(args.manifest, args.sites_manifest, args.source_roots, args.dest_root, args.output_manifest)
