"""
batch_bundle.py

A saved batch as one portable ZIP, built from the single-protein app's bundles:

    BATCH_MANIFEST.json        {"tagsites_batch": 1, "name", "proteins": {key: member}}
    {name}.batch.json          settings, statuses, sites, guides
    proteins/{key}.zip         one ordinary run bundle per protein (utils/bundle.py format)

Each protein ZIP is exactly what the single-protein app downloads and uploads, so any one of
them can also be opened there. Database proteins leave out their reagents TSV: the batch keeps
it without arms (see batch_prepare), which the single app could not use. Uploaded to the
single app, such a bundle recomputes its reagents (its existing "TSV not in bundle" path);
loaded back into the batch app, the arm-less TSV is re-exported from the database (~0.1 s).

The manifest is not named MANIFEST.json, so the single app never mistakes the batch ZIP for
a run bundle.

    python scripts/batch_bundle.py pack <batch.json> [--out batch.zip]
    python scripts/batch_bundle.py unpack <batch.zip>

Matt Rich, 2026
"""

import argparse
import io
import json
import os
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(_ROOT))

import batch_state as BS  # noqa: E402
from utils.bundle import MANIFEST_NAME, bundle_file_map, extract_and_rebase  # noqa: E402

BATCH_MANIFEST = "BATCH_MANIFEST.json"


def protein_bundle(protein):
    """One protein's run bundle as ZIP bytes (single-app format; see module docstring)."""
    run_json = protein["run_json"]
    files = bundle_file_map(run_json)
    # a database protein's TSV has blank arms; leave it out so loaders rebuild it
    if protein.get("source") == "db":
        with open(run_json) as fh:
            out = json.load(fh).get("tasks", {}).get("REAGENTS_reagents", {}).get(
                "args", {}).get("output", "")
        if out:
            files.pop(os.path.basename(out), None)
    with open(run_json) as fh:
        rn = json.load(fh).get("global", {}).get("run_name", "")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(MANIFEST_NAME, json.dumps(
            {"tagsites_bundle": 1, "run_name": rn, "run_json": os.path.basename(run_json),
             "files": sorted(files)}, indent=2))
        for arcname, path in files.items():
            zf.write(path, arcname)
    return buf.getvalue()


def pack(batch, dest):
    """Write the batch ZIP to dest (path or buffer); return the number of protein bundles."""
    members = {}
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_STORED) as zf:
        # protein bundles are already deflated; storing them avoids compressing twice
        for p in batch["proteins"]:
            if not p.get("run_json") or not os.path.exists(p["run_json"]):
                continue
            member = "proteins/{}.zip".format(BS.safe_name(p["key"]))
            zf.writestr(member, protein_bundle(p))
            members[p["key"]] = member
        zf.writestr("{}.batch.json".format(batch["name"]), json.dumps(batch, indent=1,
                                                                       default=BS._plain))
        zf.writestr(BATCH_MANIFEST, json.dumps(
            {"tagsites_batch": 1, "name": batch["name"], "proteins": members}, indent=1))
    return len(members)


def pack_bytes(batch):
    """The batch ZIP as bytes."""
    buf = io.BytesIO()
    pack(batch, buf)
    return buf.getvalue()


def is_batch_zip(path):
    """True for a ZIP written by pack()."""
    try:
        with zipfile.ZipFile(path) as zf:
            return BATCH_MANIFEST in zf.namelist()
    except Exception:
        return False


def _free_name(name, root):
    """name, or name_2, name_3 ... when a batch of that name already exists under root."""
    cand, i = name, 2
    while BS.batch_file(cand, root).exists():
        cand, i = "{}_{}".format(name, i), i + 1
    return cand


def unpack(zip_path, root=None):
    """Restore a batch ZIP under the batch folder; return the new batch file path.

    Each protein bundle goes through the single app's extract_and_rebase, then database
    proteins get their arm-less reagents TSV back from the database.
    """
    root = Path(root) if root else BS.batch_root()
    with zipfile.ZipFile(zip_path) as zf:
        manifest = json.loads(zf.read(BATCH_MANIFEST))
        batch = json.loads(zf.read("{}.batch.json".format(manifest["name"])))
        name = _free_name(BS.safe_name(manifest["name"]), root)
        batch["name"] = name
        path = BS.batch_file(name, root)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = tempfile.mkdtemp(prefix="tagsites_batchzip_")
        try:
            for p in batch["proteins"]:
                member = manifest["proteins"].get(p["key"])
                if not member:
                    p["run_json"] = None
                    continue
                # member names come from our own manifest, but stay inside tmp regardless
                local = os.path.join(tmp, os.path.basename(member))
                with zf.open(member) as src, open(local, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                p["run_json"] = extract_and_rebase(local, str(BS.protein_dir(path, p["key"])))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    BS.save_batch(path, batch)
    _restore_reagents(batch)
    BS.save_batch(path, batch)
    return str(path)


def _restore_reagents(batch):
    """Re-export the arm-less reagents TSV of database proteins left out of the bundle."""
    import proteome_export as X
    for p in batch["proteins"]:
        if p.get("source") != "db" or not p.get("run_json") or p.get("reagent_status") != "success":
            continue
        with open(p["run_json"]) as fh:
            gb = json.load(fh)["global"]
        try:
            X.export_reagents(p["pid"], gb["working_dir"], gb["run_name"], with_arms=False)
        except Exception as e:   # recorded on the entry; the analyses still load
            p.update(status="failed", error="reagents not restored: {}".format(e))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Pack or unpack a batch ZIP.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("pack")
    a1.add_argument("batch_json")
    a1.add_argument("--out", default=None)
    a2 = sub.add_parser("unpack")
    a2.add_argument("zip_path")
    a = ap.parse_args()
    if a.cmd == "pack":
        b = BS.load_batch(a.batch_json)
        out = a.out or os.path.join(os.path.dirname(a.batch_json), b["name"] + "_batch.zip")
        print("wrote {} ({} proteins)".format(out, pack(b, out)))
    else:
        print("restored to", unpack(a.zip_path))
