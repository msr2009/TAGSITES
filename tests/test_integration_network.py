"""
test_integration_network.py — network-dependent integration tests

All tests in this file are marked @pytest.mark.network and are SKIPPED by
default.  To run them:

    TAGSITES_EMAIL=your@email.com pytest -m network --run-network

These tests call the EBI REST API (~10–120 s each) and require a valid
registered e-mail address.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent


# ── Genewise both-strand selection ────────────────────────────────────────────

@pytest.mark.network
def test_run_genewise_selects_forward_for_src1(DATA, email, tmp_path):
    """
    run_genewise.py on src-1 should select the + orientation (score ~1138 bits)
    and write the expected output files.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    outprefix = tmp_path / "src-1_gw"
    ret = subprocess.call(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "run_genewise.py"),
            "--protein_fasta", str(DATA / "src-1.fa"),
            "--genomic_fasta", str(DATA / "src-1_genomic.fa"),
            "--email", email,
            "--outprefix", str(outprefix),
        ]
    )
    assert ret == 0, "run_genewise.py exited with non-zero status"

    winner = Path(str(outprefix) + ".genewise.out.txt")
    orientation_file = Path(str(outprefix) + ".genewise_orientation.txt")
    assert winner.exists(), "winner .genewise.out.txt not written"
    assert orientation_file.exists(), "orientation file not written"

    orientation = orientation_file.read_text().strip()
    assert orientation == "+", f"expected '+', got '{orientation}'"

    # Score of the winner should be high
    from parse_genewise import parse_genewise_score
    score = parse_genewise_score(winner)
    assert score > 500, f"winner score unexpectedly low: {score}"


# ── End-to-end pipeline via JSON ─────────────────────────────────────────────

@pytest.mark.network
def test_full_pipeline_reagents_tsv(DATA, email, tmp_path):
    """
    Run run_tag_sites_from_json.py with a JSON that has an empty 'genewise' arg
    (triggering the genewise_required pre-step) and assert the reagents TSV is
    produced with non-zero rows.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    working_dir = tmp_path / "run"
    working_dir.mkdir()

    # Copy fixtures into a working dir (pipeline writes outputs relative to working_dir)
    import shutil
    prot_fa  = working_dir / "snb-1.fa"
    genomic_fa = working_dir / "snb-1_genomic.fa"
    shutil.copy(DATA / "snb-1.fa", prot_fa)
    shutil.copy(DATA / "snb-1_genomic.fa", genomic_fa)

    output_tsv = working_dir / "snb-1.reagents.tsv"

    run_json = {
        "scripts": {
            "reagents": "design_tag_reagents.py"
        },
        "global": {
            "email": email,
            "working_dir": str(working_dir),
            "run_name": "snb-1",
            "input_file": str(prot_fa),
            "pdb": "",
            "scripts_folder": str(REPO_ROOT / "scripts") + "/"
        },
        "REAGENTS_reagents": {
            "analysis": "reagents",
            "args": {
                "genomic_fasta": str(genomic_fa),
                "genewise": "",   # empty → triggers genewise_required pre-step
                "n_guides": 3,
                "arm_length": 200,
                "PAM": "NGG",
                "guide_length": 20,
                "cut_offset": 3,
                "output": str(output_tsv)
            }
        }
    }

    json_path = working_dir / "snb-1.json"
    json_path.write_text(json.dumps(run_json, indent=4))

    ret = subprocess.call(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "run_tag_sites_from_json.py"),
            "-i", str(json_path),
        ]
    )
    assert ret == 0, "run_tag_sites_from_json.py exited with non-zero status"
    assert output_tsv.exists(), "reagents TSV was not written"

    import pandas as pd
    df = pd.read_csv(output_tsv, sep="\t")
    assert len(df) > 0, "reagents TSV is empty"
    assert "spacer" in df.columns
    assert "left_arm" in df.columns
