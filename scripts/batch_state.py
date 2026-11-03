"""
batch_state.py

The batch app's one persistent record: {batch_dir}/{name}/{name}.batch.json. It holds the
batch settings (repair type, inserts, arm length) and, per protein, the resolved product, its
run folder, review status, accepted sites and the guide chosen for each site. Every change is
written atomically (temp file then replace, as run_status does), so a crash loses at most the
last click; checkpoint() also keeps timestamped copies.

Per-protein status: pending -> ready (exported) -> in_progress (sites edited) -> accepted;
failed carries the error; skipped (set from Review, any reviewable status) is listed, one
protein per line, in {name}.skipped.txt beside the batch file, rewritten on every save. Each protein's run folder is {batch_dir}/{name}/{key}/, an ordinary
run folder every app loader reads unchanged.

    python scripts/batch_state.py list
    python scripts/batch_state.py show <batch.json>

Matt Rich, 2026
"""

import argparse
import glob
import json
import os
import re
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import _load_config as _load_batch_config  # noqa: E402

STATUSES = ("pending", "running", "ready", "in_progress", "accepted", "skipped", "failed",
            "missing")
CHECKPOINTS_KEPT = 5
_DEFAULT_DIR = "~/tagsites_local/batches"
_lock = threading.Lock()


def app_config():
    """The batch_app block of the active batch config ({} when absent)."""
    return _load_batch_config().get("batch_app", {})


def batch_root():
    """Directory holding every batch folder (batch_app.batch_dir, local disk)."""
    return Path(os.path.expanduser(app_config().get("batch_dir") or _DEFAULT_DIR))


def safe_name(text):
    """A file-system-safe token (letters, digits, . _ -)."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_") or "x"


def batch_file(name, root=None):
    """Path of a batch's JSON file."""
    root = Path(root) if root else batch_root()
    return root / name / "{}.batch.json".format(name)


def protein_dir(batch_path, key):
    """Run folder of one protein of a batch."""
    return Path(batch_path).parent / safe_name(key)


def name_prefix(protein):
    """Reagent-name prefix {gene}_{uniprot_id}; the product id alone when there is no gene."""
    gene = protein.get("gene") or ""
    pid = protein.get("product_id") or protein.get("input") or "protein"
    # a gene named after its own id (WormBase-only rows) would only repeat it
    if gene and gene != pid:
        return "{}_{}".format(safe_name(gene), safe_name(pid))
    return safe_name(pid)


def new_batch(name, settings, proteins, root=None):
    """Create a batch folder and file; return the file path."""
    path = batch_file(name, root)
    path.parent.mkdir(parents=True, exist_ok=True)
    # an input resolving to a product another input already uses adds nothing
    proteins = [p for p in proteins if not p.get("duplicate_of")]
    batch = {"name": name, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
             "settings": settings, "proteins": proteins, "current_key": None}
    for p in batch["proteins"]:
        p.setdefault("status", "pending" if p.get("pid") is not None else "missing")
        p.setdefault("sites", [])
        p.setdefault("guides", {})
        p.setdefault("error", "")
        p.setdefault("accepted_at", None)
    save_batch(path, batch)
    return str(path)


def load_batch(path):
    """Read a batch file."""
    with open(path) as fh:
        return json.load(fh)


def _plain(o):
    """JSON fallback for numpy scalars (sites and pids often arrive as numpy ints)."""
    if hasattr(o, "item"):
        return o.item()
    raise TypeError("not JSON serializable: {}".format(type(o).__name__))


def skipped_file(path):
    """{name}.skipped.txt beside a batch file."""
    return Path(str(path)[:-len(".batch.json")] + ".skipped.txt")


def _write_atomic(path, text):
    """Write text to path through a temp file."""
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def write_skipped(path, batch):
    """Rewrite {name}.skipped.txt: one skipped protein per line (absent while none ever were)."""
    out = skipped_file(path)
    rows = [p for p in batch["proteins"] if p.get("status") == "skipped"]
    # no file for a batch that never skipped anything; an emptied list keeps its header
    if not rows and not out.exists():
        return
    lines = ["# key\tname\tinput\tskipped_at"]
    lines += ["\t".join([p["key"], name_prefix(p), str(p.get("input") or ""),
                         str(p.get("skipped_at") or "")]) for p in rows]
    _write_atomic(out, "\n".join(lines) + "\n")


def save_batch(path, batch):
    """Write the batch atomically (tmp then replace), and its skipped list."""
    with _lock:
        _write_atomic(path, json.dumps(batch, indent=1, default=_plain))
        # checkpoint copies (.batch.<stamp>.json) have no skipped list of their own
        if str(path).endswith(".batch.json"):
            write_skipped(path, batch)


def checkpoint(path, batch):
    """Save, plus a timestamped copy; keep the newest CHECKPOINTS_KEPT copies."""
    save_batch(path, batch)
    stem = str(path)[:-len(".json")]
    copy = "{}.{}.json".format(stem, time.strftime("%Y%m%d-%H%M%S"))
    save_batch(copy, batch)
    old = sorted(glob.glob(stem + ".2*.json"))
    # oldest first; drop all but the newest few
    for extra in old[:-CHECKPOINTS_KEPT]:
        os.remove(extra)
    return copy


def protein(batch, key):
    """The protein entry with this key, or None."""
    return next((p for p in batch["proteins"] if p["key"] == key), None)


def update_protein(path, batch, key, **fields):
    """Set fields on one protein entry and save the batch."""
    p = protein(batch, key)
    if p is None:
        raise KeyError(key)
    p.update(fields)
    save_batch(path, batch)
    return p


def counts(batch):
    """{status: n} over the batch's proteins."""
    out = {s: 0 for s in STATUSES}
    for p in batch["proteins"]:
        out[p.get("status", "pending")] = out.get(p.get("status", "pending"), 0) + 1
    return out


def next_unaccepted(batch, after=None):
    """Key of the first reviewable, not-yet-accepted protein (after `after` when given)."""
    keys = [p["key"] for p in batch["proteins"]
            if p.get("status") in ("ready", "in_progress")]
    if not keys:
        return None
    if after in keys:
        later = keys[keys.index(after) + 1:]
        return later[0] if later else keys[0]
    if after is not None:
        order = [p["key"] for p in batch["proteins"]]
        pos = order.index(after) if after in order else -1
        later = [k for k in keys if order.index(k) > pos]
        return later[0] if later else keys[0]
    return keys[0]


def list_batches(root=None):
    """Summaries of every batch under the root, newest first."""
    root = Path(root) if root else batch_root()
    out = []
    # checkpoints are {name}.batch.{stamp}.json, so this glob sees only live files
    for path in glob.glob(str(root / "*" / "*.batch.json")):
        try:
            b = load_batch(path)
        except Exception:
            continue
        c = counts(b)
        out.append({"path": path, "name": b.get("name"), "created": b.get("created", ""),
                    "n": len(b.get("proteins", [])), "accepted": c.get("accepted", 0)})
    return sorted(out, key=lambda r: r["created"], reverse=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Inspect batch-app state files.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    s = sub.add_parser("show")
    s.add_argument("path")
    a = ap.parse_args()
    if a.cmd == "list":
        for r in list_batches():
            print("{created}  {name}  {accepted}/{n} accepted  {path}".format(**r))
    else:
        b = load_batch(a.path)
        print(json.dumps(counts(b)))
        for p in b["proteins"]:
            print("{:<16} {:<12} {:<20} sites={}".format(
                p["key"], p.get("status"), name_prefix(p), p.get("sites")))
