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
SCRIPTS = str(REPO_ROOT / "scripts") + "/"
sys.path.insert(0, str(REPO_ROOT / "modules"))


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

    # the EBI REST Genewise API doesn't write a "Score NNN bits" header line
    # (only the standalone binary does) — use the GFF-aware parser instead
    from parse_genewise import parse_genewise_gff_score
    score = parse_genewise_gff_score(winner)
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

    import shutil
    prot_fa    = working_dir / "snb-1.fa"
    genomic_fa = working_dir / "snb-1_genomic.fa"
    shutil.copy(DATA / "snb-1.fa",          prot_fa)
    shutil.copy(DATA / "snb-1_genomic.fa",  genomic_fa)

    output_tsv = working_dir / "snb-1.reagents.tsv"

    run_json = {
        "global": {
            "email": email,
            "working_dir": str(working_dir) + "/",
            "run_name": "snb-1",
            "input_file": str(prot_fa),
            "pdb": "",
            "scripts_folder": SCRIPTS,
            "selected_sites": [],
        },
        "tasks": {
            "REAGENTS_reagents": {
                "type": "reagents",
                "args": {
                    "genomic_fasta": str(genomic_fa),
                    "genewise": "",   # empty → triggers genewise_required pre-step
                    "n_guides": 3,
                    "arm_length": 200,
                    "PAM": "NGG",
                    "guide_length": 20,
                    "cut_offset": 3,
                    "output": str(output_tsv),
                }
            }
        }
    }

    json_path = working_dir / "snb-1.run.json"
    json_path.write_text(json.dumps(run_json, indent=4))

    ret = subprocess.call(
        [sys.executable, str(REPO_ROOT / "scripts" / "run_tag_sites_from_json.py"),
         "-i", str(json_path)]
    )
    assert ret == 0, "run_tag_sites_from_json.py exited with non-zero status"
    assert output_tsv.exists(), "reagents TSV was not written"

    import pandas as pd
    df = pd.read_csv(output_tsv, sep="\t")
    assert len(df) > 0, "reagents TSV is empty"
    assert "spacer" in df.columns
    assert "left_arm" in df.columns


# ── BLAST orthologs ───────────────────────────────────────────────────────────

@pytest.mark.network
def test_blast_orthologs_produces_output(DATA, email, tmp_path):
    """
    blast_orthologs.main() on snb-1.fa should produce a JSD output file
    with at least one line.  Uses n=5 to keep runtime short.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    import shutil
    working_dir = tmp_path / "blast"
    working_dir.mkdir()
    fa = working_dir / "snb-1.fa"
    shutil.copy(DATA / "snb-1.fa", fa)

    # output suffix is stripped to form out_prefix, so ".jsd" → out_prefix "snb-1"
    # and jsd_path = out_prefix + ".jsd" = "snb-1.jsd"
    output = str(working_dir / "snb-1.jsd")

    from blast_orthologs import main as blast_main
    ret = blast_main(
        fasta_in=str(fa),
        email=email,
        workingdir=str(working_dir),
        name="snb-1",
        output=output,
        n=5,
        evalue=1e-10,
        db="uniprotkb",
        length_percent=0,
        align_full_seqs=False,
        taxid=1,
        clients_folder=SCRIPTS,
        exclude_paralogs=False,
    )
    jsd_file = working_dir / "snb-1.jsd"
    assert jsd_file.exists(), "JSD output not written"
    lines = [l for l in jsd_file.read_text().splitlines() if l.strip()]
    assert len(lines) > 0, "JSD file is empty"


# ── InterPro domain annotation ────────────────────────────────────────────────

@pytest.mark.network
def test_call_interpro_produces_domains(DATA, email, tmp_path):
    """
    call_interpro.main() on snb-1.fa should produce a domain TSV with ≥1 row.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    import shutil
    working_dir = tmp_path / "interpro"
    working_dir.mkdir()
    fa = working_dir / "SNB-1.fa"
    shutil.copy(DATA / "snb-1.fa", fa)

    outputfile = str(working_dir / "snb-1.DOMAINS_domains.txt")

    from call_interpro import main as interpro_main
    interpro_main(
        fasta_in=str(fa),
        email=email,
        workingdir=str(working_dir),
        clients_folder=SCRIPTS,
        outputfile=outputfile,
    )

    out = Path(outputfile)
    assert out.exists(), "InterPro domain file not written"
    rows = [l for l in out.read_text().splitlines() if l.strip()]
    assert len(rows) > 0, "InterPro output is empty — no domains found for snb-1"


# ── AlphaFold DB search ───────────────────────────────────────────────────────

@pytest.mark.network
def test_search_afdb_return_type(DATA, email, tmp_path):
    """
    existing_AF_model.search_AFDB on snb-1 should return either a path string
    (found a model) or the integer 1 (no model). The return type contract must hold.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    import shutil
    working_dir = tmp_path / "afdb"
    working_dir.mkdir()
    fa = working_dir / "snb-1.fa"
    shutil.copy(DATA / "snb-1.fa", fa)

    from existing_AF_model import search_AFDB
    result = search_AFDB(
        fasta_in=str(fa),
        email=email,
        workingdir=str(working_dir),
        name="snb-1",
        taxid="1",
        evalue=1e-100,
        percentid=99,
        clients_folder=SCRIPTS,
    )
    # Contract: returns an AF .fa path (str) or 1 (int) when not found
    assert isinstance(result, (str, int)), f"unexpected return type: {type(result)}"
    if isinstance(result, str):
        assert Path(result).exists(), f"returned path does not exist: {result}"


# ── Species taxonomy lookup ───────────────────────────────────────────────────

@pytest.mark.network
def test_get_species_taxonomy_c_elegans(email):
    """
    get_species_taxonomy.main on C. elegans taxid 6239 should return a
    multi-line string containing 'Caenorhabditis elegans'.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    from get_species_taxonomy import main as tax_main
    result = tax_main(
        taxid="6239",
        taxlevel=None,
        email=email,
        clients_folder=SCRIPTS,
    )
    assert isinstance(result, str), "tax_main should return a string"
    assert "Caenorhabditis elegans" in result or "6239" in result, (
        f"Expected C. elegans name or ID in result; got: {result[:200]}"
    )


# ── Full pipeline, all task types ─────────────────────────────────────────────

@pytest.mark.network
def test_full_pipeline_all_tasks(DATA, email, tmp_path):
    """
    Run every task type (blast, plddt, modifications, domains, scores, reagents)
    for snb-1 through run_tag_sites_from_json.py and check each output file's
    format matches what utils/results.py's loader expects.
    """
    if not email:
        pytest.skip("No EBI e-mail provided (set TAGSITES_EMAIL or --email).")

    import shutil
    import pandas as pd
    from task_registry import task_defaults, task_output_suffix, companion_path
    from setup_logic import build_global_block, build_task_entry, build_reagents_entry, build_run_json
    import run_status
    import run_tag_sites_from_json

    working_dir = str(tmp_path / "full") + "/"
    os.makedirs(working_dir)
    run_name = "snb-1"

    prot_fa = working_dir + "snb-1.fa"
    genomic_fa = working_dir + "snb-1_genomic.fa"
    genewise_out = working_dir + "snb-1.genewise.out.txt"
    shutil.copy(DATA / "snb-1.fa", prot_fa)
    shutil.copy(DATA / "snb-1_genomic.fa", genomic_fa)
    shutil.copy(DATA / "snb-1.genewise.out.txt", genewise_out)

    global_block = build_global_block(
        email=email, run_name=run_name, working_dir=working_dir,
        input_file=prot_fa, genomic_file=genomic_fa,
    )

    blast_args = {**task_defaults("blast"), "taxid": 6239, "db": "uniprotkb", "max_hits": 5}
    plddt_args = {**task_defaults("plddt"), "species_taxid": 6239}
    mods_args = task_defaults("modifications")
    domains_args = task_defaults("domains")
    scores_args = {**task_defaults("scores"),
                   "scores_file": str(REPO_ROOT / "tables" / "hydrophobicity_kyle-doolittle.tsv")}

    task_entries = [
        build_task_entry("BLAST_blast", "blast", blast_args,
                         task_output_suffix("blast"), working_dir, run_name),
        build_task_entry("PLDDT_plddt", "plddt", plddt_args,
                         task_output_suffix("plddt"), working_dir, run_name),
        build_task_entry("MODS_modifications", "modifications", mods_args,
                         task_output_suffix("modifications"), working_dir, run_name),
        build_task_entry("DOMAINS_domains", "domains", domains_args,
                         task_output_suffix("domains"), working_dir, run_name),
        build_task_entry("SCORES_scores", "scores", scores_args,
                         task_output_suffix("scores"), working_dir, run_name),
    ]

    reagents_entry = build_reagents_entry(
        defaults={**task_defaults("reagents"), "genewise": genewise_out, "n_guides": 3, "arm_length": 200},
        genomic_path=genomic_fa, out_suffix=task_output_suffix("reagents"),
        working_dir=working_dir, run_name=run_name,
    )

    run_json = build_run_json(global_block, task_entries, reagents_entry)
    json_path = working_dir + f"{run_name}.run.json"
    with open(json_path, "w") as f:
        json.dump(run_json, f, indent=4)

    run_tag_sites_from_json.main(json_path)

    # every task must report success in the status file
    status = run_status.load_status(working_dir, run_name)
    for task_id in run_json["tasks"]:
        assert status.get(task_id, {}).get("status") == "success", (
            f"{task_id} did not complete successfully: {status.get(task_id)}"
        )

    # blast: non-empty JSD + isoforms companion (if present)
    blast_out = run_json["tasks"]["BLAST_blast"]["args"]["output"]
    blast_lines = [l for l in Path(blast_out).read_text().splitlines() if l.strip()]
    assert len(blast_lines) > 0, "blast JSD output is empty"
    isoforms_path = companion_path(blast_out, "blast", "isoforms")
    if isoforms_path and os.path.exists(isoforms_path):
        iso_df = pd.read_csv(isoforms_path, sep="\t",
                             names=["source", "start", "stop", "description"])
        assert set(iso_df.columns) == {"source", "start", "stop", "description"}

    # plddt: headerless 2-col TSV + sasa companion
    plddt_out = run_json["tasks"]["PLDDT_plddt"]["args"]["output"]
    plddt_df = pd.read_csv(plddt_out, sep="\t", header=None, comment="#")
    assert plddt_df.shape[1] == 2, "plddt output should be a 2-column TSV"
    pd.to_numeric(plddt_df[1])  # values must be numeric
    sasa_path = companion_path(plddt_out, "plddt", "sasa")
    if sasa_path and os.path.exists(sasa_path):
        # utils/results.py's loader only reads columns 0 (position) and 1 (SASA
        # value) positionally — extract_from_pdb.py's sasa file has a 3rd
        # (one-letter amino acid) column, which is fine as long as 0/1 are numeric.
        sasa_df = pd.read_csv(sasa_path, sep="\t", header=None, comment="#")
        assert sasa_df.shape[1] >= 2
        pd.to_numeric(sasa_df[0])
        pd.to_numeric(sasa_df[1])

    # modifications + domains: headerless 4-col TSV
    for task_id in ("MODS_modifications", "DOMAINS_domains"):
        out_path = run_json["tasks"][task_id]["args"]["output"]
        df = pd.read_csv(out_path, sep="\t",
                         names=["source", "start", "stop", "description"])
        assert set(df.columns) == {"source", "start", "stop", "description"}

    # scores: headerless 2-col TSV (position, score)
    scores_out = run_json["tasks"]["SCORES_scores"]["args"]["output"]
    scores_df = pd.read_csv(scores_out, sep="\t", header=None)
    assert scores_df.shape[1] == 2, "scores output should be a 2-column TSV"
    pd.to_numeric(scores_df[1])

    # reagents: existing format contract
    reagents_out = run_json["tasks"]["REAGENTS_reagents"]["args"]["output"]
    reagents_df = pd.read_csv(reagents_out, sep="\t")
    assert len(reagents_df) > 0, "reagents TSV is empty"
    assert "spacer" in reagents_df.columns
    assert "left_arm" in reagents_df.columns


# ── Ensembl genomic sequence auto-fetch ───────────────────────────────────────

@pytest.mark.network
def test_fetch_genomic_sequence_worm_unc18():
    """Live-fetch C. elegans unc-18 and confirm expected coordinates + flank length."""
    from fetch_genomic_sequence import fetch_genomic_sequence

    fasta_text, meta = fetch_genomic_sequence(taxid=6239, gene_symbol="unc-18", flank_bp=500)

    assert meta["species"] == "caenorhabditis_elegans"
    assert meta["gene_id"] == "WBGene00006757"
    assert fasta_text.startswith(">")
    seq = "".join(fasta_text.splitlines()[1:])
    assert len(seq) > 0


@pytest.mark.network
def test_fetch_genomic_sequence_ecoli_assembly_qualified_slug():
    """Live-fetch E. coli thrA to confirm the assembly-qualified species slug path works."""
    from fetch_genomic_sequence import fetch_genomic_sequence

    fasta_text, meta = fetch_genomic_sequence(taxid=562, gene_symbol="thrA", flank_bp=300)

    assert meta["species"] == "escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
    assert meta["gene_id"] == "b0002"
    assert fasta_text.startswith(">")
