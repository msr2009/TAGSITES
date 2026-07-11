"""bundle.py — portable session bundles for TAGSITES.

A bundle is a flat ZIP holding everything needed to restore a run on any
instance (local or shinyapps.io): the run JSON, its status JSON, the query
FASTA/PDB, every task output plus companions (.aln/.pdf/.isoforms/.sasa),
and the reagents TSV.

Restore works by extracting the bundle into a directory and rewriting the
run JSON's path fields to point there, so the existing readers (which build
paths from working_dir / args.output) resolve correctly with no other changes.

No Shiny imports — safe to use from the app and from CLI scripts.
"""

import io
import json
import os
import zipfile

from scripts.task_registry import task_companions, companion_path

MANIFEST_NAME = "MANIFEST.json"

# run-JSON global fields that hold file paths and must be re-based on extract
_GLOBAL_PATH_KEYS = ("input_file", "pdb", "genomic_file")
# reagents-task arg fields that hold file paths (present only when reagents ran)
_REAGENT_PATH_KEYS = ("genomic_fasta", "genewise")
# task arg fields that hold file paths and must be re-based on extract
_TASK_PATH_KEYS = ("pdb",) + _REAGENT_PATH_KEYS


# ── shared JSON loader ──────────────────────────────────────────────────────

def _load_json(path):
    """Load a JSON file into a dict; raise on failure."""
    with open(path) as f:
        return json.load(f)


# ── download side: enumerate + zip ──────────────────────────────────────────

def bundle_file_map(run_json_path):
    """Return {arcname: abspath} of every existing file needed to restore a run.

    Paths are enumerated exactly as the readers derive them (task outputs,
    companions via companion_path, alignment PNGs, globals), so the bundle is
    complete by construction. Files that don't exist are skipped.
    """
    j = _load_json(run_json_path)
    g = j.get("global", {})
    wd = g.get("working_dir", "")
    rn = g.get("run_name", "")

    files = {}

    def _add(path):
        """Add a file to the map under its basename if it exists."""
        if path and os.path.exists(path):
            files[os.path.basename(path)] = path

    # run JSON itself + status JSON
    _add(run_json_path)
    if wd and rn:
        _add(os.path.join(wd, f"{rn}.status.json"))

    # global inputs (query FASTA, structure, genomic)
    for key in _GLOBAL_PATH_KEYS:
        _add(g.get(key, ""))
    # AFDB-downloaded PDB is not in a global field — reconstruct it
    if wd and rn:
        _add(os.path.join(wd, f"{rn}.AF.pdb"))

    # per-task outputs + companions
    for _, v in j.get("tasks", {}).items():
        analysis = v.get("type", "")
        output = v.get("args", {}).get("output", "")
        if not output:
            continue
        _add(output)
        for key in task_companions(analysis):
            _add(companion_path(output, analysis, key))
        # blast alignment image is derived from the .aln path, not a companion
        if analysis == "blast":
            aln = companion_path(output, analysis, "alignment")
            if aln.endswith(".aln"):
                _add(aln[:-len(".aln")] + ".pdf")

    # genewise outputs: allow fast reagents re-run after bundle restore
    # (stored at {wd}{rn}_genewise.* regardless of what the reagents args say)
    if wd and rn:
        pfx = os.path.join(wd, f"{rn}_genewise")
        for suffix in (".genewise.out.txt", ".genewise_genomic.fa"):
            _add(pfx + suffix)

    return files


def make_bundle(run_json_path):
    """Return ZIP bytes: a flat bundle of all restore files plus a MANIFEST."""
    j = _load_json(run_json_path)
    rn = j.get("global", {}).get("run_name", "")
    files = bundle_file_map(run_json_path)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        manifest = {"tagsites_bundle": 1, "run_name": rn,
                    "run_json": os.path.basename(run_json_path),
                    "files": sorted(files)}
        zf.writestr(MANIFEST_NAME, json.dumps(manifest, indent=2))
        for arcname, path in files.items():
            zf.write(path, arcname)
    buf.seek(0)
    return buf.read()


# ── upload side: detect + extract + rebase ──────────────────────────────────

def _is_junk_entry(name):
    """True for macOS Finder zip cruft: __MACOSX/ AppleDouble sidecar entries.

    Zips created via Finder's "Compress" action embed a ._<name> resource-fork
    shadow file next to every real entry (under a top-level __MACOSX/ folder).
    These shadow files match naming heuristics like "*.run.json" but hold
    binary data, so they must be excluded from bundle detection and extraction.
    """
    base = os.path.basename(name)
    return name.startswith("__MACOSX/") or base.startswith("._")


def is_bundle_zip(path):
    """True if path is a ZIP holding a MANIFEST or exactly one *.run.json."""
    if not path or not zipfile.is_zipfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as zf:
            names = zf.namelist()
    except Exception:
        return False
    names = [n for n in names if not _is_junk_entry(n)]
    if MANIFEST_NAME in names:
        return True
    return sum(1 for n in names if n.endswith(".run.json")) == 1


def _rebase(path, dest_dir):
    """Point a stored file path at dest_dir, keeping its basename."""
    if not path:
        return path
    return os.path.join(dest_dir, os.path.basename(path))


def extract_and_rebase(zip_path, dest_dir):
    """Extract a bundle into dest_dir and rewrite its run JSON's paths.

    Returns the path to the rewritten run JSON (written into dest_dir).
    Members are written by basename only to avoid zip path traversal.
    """
    os.makedirs(dest_dir, exist_ok=True)
    run_json_member = None
    manifest = {}
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        # MANIFEST names the run JSON explicitly (naming may be .json or .run.json)
        if MANIFEST_NAME in names:
            try:
                manifest = json.loads(zf.read(MANIFEST_NAME))
            except Exception:
                manifest = {}
        for name in names:
            base = os.path.basename(name)
            if not base or name == MANIFEST_NAME or _is_junk_entry(name):
                continue   # skip directory entries, the manifest, and Finder zip cruft
            with zf.open(name) as src, open(os.path.join(dest_dir, base), "wb") as dst:
                dst.write(src.read())
            if base.endswith(".run.json"):
                run_json_member = base

    # prefer the manifest's declared run JSON, else the *.run.json we saw
    run_json_member = manifest.get("run_json") or run_json_member
    if run_json_member is None:
        raise ValueError("Bundle contains no run JSON")

    run_json_path = os.path.join(dest_dir, run_json_member)
    j = _load_json(run_json_path)

    # rebase the global block: working_dir becomes dest_dir, file fields keep basename
    g = j.get("global", {})
    g["working_dir"] = os.path.join(dest_dir, "")   # trailing sep for wd+rn+suffix
    for key in _GLOBAL_PATH_KEYS:
        if g.get(key):
            g[key] = _rebase(g[key], dest_dir)

    # rebase each task's output + path args (pdb, genomic_fasta, genewise)
    for _, v in j.get("tasks", {}).items():
        args = v.get("args", {})
        if args.get("output"):
            args["output"] = _rebase(args["output"], dest_dir)
        for key in _TASK_PATH_KEYS:
            if args.get(key):
                args[key] = _rebase(args[key], dest_dir)

    with open(run_json_path, "w") as f:
        json.dump(j, f, indent=4)
    return run_json_path
