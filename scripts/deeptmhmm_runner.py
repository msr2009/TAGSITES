"""
deeptmhmm_runner.py

Shared subprocess bridge to DeepTMHMM 1.0 (academic license), used by
scripts/topology_deeptmhmm.py (per-protein backend) and
scripts/build_topology_cache.py (bulk pre-scan). DeepTMHMM's own predict.py
is not importable — it runs argparse at module scope, resolves its model
files relative to cwd, and exit(1)s if the output directory already exists —
so it is driven as a subprocess in a dedicated conda env, with vendor code
left unmodified (required by its academic license).

DeepTMHMM writes {output_dir}/TMRs.gff3 with rows
`prot_id \t region \t start \t end`, region one of inside/outside/TMhelix/
signal/periplasm/Beta sheet. Those region words are kept verbatim here (not
translated into Phobius's Transmembrane/Cytoplasmic/... vocabulary) — see
~/.claude/plans/i-m-considering-a-large-humble-sun.md.

generate_esm_embeddings() writes one embedding file per sequence into
{output_dir}/embeddings and never deletes them; since that directory lives
inside the temp output dir this module creates per call, deleting the whole
temp dir after parsing reclaims that scratch.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import _load_config as _load_batch_config

_KNOWN_REGIONS = {"inside", "outside", "TMhelix", "signal", "periplasm", "Beta sheet"}


def _deeptmhmm_config():
    """Read the deeptmhmm config block (install_dir, python, ...) via the same
    active-config path providers.py already resolves, so TAGSITES_BATCH_CONFIG
    selects it consistently with every other backend.
    """
    cfg = _load_batch_config().get("deeptmhmm", {})
    if "install_dir" not in cfg or "python" not in cfg:
        raise RuntimeError(
            "deeptmhmm config block missing install_dir/python — set batch.config.json's "
            '"deeptmhmm": {"install_dir": ..., "python": ...} to the DeepTMHMM checkout '
            "and its conda env's python interpreter."
        )
    return cfg


def _parse_tmrs_gff3(path):
    """Read TMRs.gff3 into {seq_id: [(region, start, stop), ...]}, one entry per
    sequence including sequences with zero regions (an empty list is a legal,
    meaningful "predicted, no regions" result, not "not predicted").
    """
    results = {}
    current_id = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line == "//" or line.startswith("##"):
                # a "# <id> Length: ..." comment line announces the next record
                if line.startswith("# ") and "Length:" in line:
                    current_id = line[2:].split(" Length:")[0].strip()
                    results.setdefault(current_id, [])
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            seq_id, region, start, stop = fields[0], fields[1], fields[2], fields[3]
            if region not in _KNOWN_REGIONS:
                print(f"[deeptmhmm_runner] warning: unrecognized region '{region}' for {seq_id}",
                      file=sys.stderr)
            results.setdefault(seq_id, []).append((region, int(start), int(stop)))
    return results


def predict_topology(records):
    """Run DeepTMHMM over `records` ([(name, seq), ...]) and return
    {name: [(region, start, stop), ...]}, including an empty list for any
    sequence DeepTMHMM assigned zero regions.

    Sequences are deduplicated before invoking: predict.py keys its internal
    prediction maps by sequence content, not id, so two identical sequences
    under different ids would otherwise silently collapse to one output
    record. The result is fanned back out to every original name.
    """
    cfg = _deeptmhmm_config()
    install_dir = Path(cfg["install_dir"])
    python_bin = cfg["python"]

    # dedupe by sequence, keep first name seen as the representative id
    seq_to_names = {}
    for name, seq in records:
        seq_to_names.setdefault(seq, []).append(name)
    unique_records = [(names[0], seq) for seq, names in seq_to_names.items()]

    tmp_root = tempfile.mkdtemp(prefix="deeptmhmm_")
    try:
        fasta_path = Path(tmp_root) / "input.fasta"
        with open(fasta_path, "w") as f:
            for name, seq in unique_records:
                f.write(f">{name}\n{seq}\n")

        out_dir = Path(tmp_root) / "out"
        result = subprocess.run(
            [python_bin, "predict.py", "--fasta", str(fasta_path), "--output-dir", str(out_dir)],
            cwd=str(install_dir), capture_output=True, text=True,
        )
        gff3_path = out_dir / "TMRs.gff3"
        if not gff3_path.exists():
            raise RuntimeError(
                f"DeepTMHMM did not produce TMRs.gff3 (exit {result.returncode}). "
                f"stderr tail:\n{result.stderr[-2000:]}"
            )
        by_representative = _parse_tmrs_gff3(gff3_path)
    finally:
        # reclaims the ESM embeddings written under out_dir/embeddings
        shutil.rmtree(tmp_root, ignore_errors=True)

    # fan the deduped result back out to every original name
    results = {}
    for seq, names in seq_to_names.items():
        representative = names[0]
        regions = by_representative.get(representative, [])
        for name in names:
            results[name] = regions
    return results
