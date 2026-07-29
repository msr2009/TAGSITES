"""
test_uniprot_network.py — network-dependent contract tests for the UniProt REST helpers.

All tests are marked @pytest.mark.network and are SKIPPED by default. To run:

    pytest -m network --run-network tests/test_uniprot_network.py

These guard against UniProt changing its response shape or feature-type
vocabulary out from under scripts/uniprot_api.py and scripts/uniprot_features.py
-- a silent change there would degrade to an empty feature track with no error.
Each call is a single sub-second GET; no job queue, no email required.
"""

import json
from pathlib import Path

import pytest

from uniprot_api import checksum_lookup, fetch_entry
from uniprot_features import load_feature_config, entry_to_feature_rows, main as uniprot_features_main

REPO_ROOT = Path(__file__).parent.parent
FEATURES_FILE = REPO_ROOT / "tables" / "uniprot_features.txt"
FIXTURE = REPO_ROOT / "tests" / "data" / "Q22038.uniprot.json"

RHO1_ACCESSION = "Q22038"
RHO1_TAXID = 6239


@pytest.mark.network
def test_fetch_entry_response_shape():
    """fetch_entry() still returns the fields we depend on."""
    entry = fetch_entry(RHO1_ACCESSION)
    assert entry["primaryAccession"] == RHO1_ACCESSION
    assert len(entry["sequence"]["value"]) == 192
    assert entry["features"]


@pytest.mark.network
def test_live_feature_types_include_tracked_blocking_type():
    """UniProt hasn't renamed 'Lipidation' or moved it off Cys189 for rho-1."""
    entry = fetch_entry(RHO1_ACCESSION)
    live_types = {f["type"] for f in entry["features"]}
    feature_config = load_feature_config(str(FEATURES_FILE))
    blocking_types = {t for t, (sev, _) in feature_config.items() if sev == "blocking"}
    assert live_types & blocking_types

    lipidation = [f for f in entry["features"] if f["type"] == "Lipidation"]
    assert lipidation
    assert lipidation[0]["location"]["start"]["value"] == 189


@pytest.mark.network
def test_checksum_lookup_resolves_rho1():
    """checksum: query syntax and response parsing still work for a known sequence."""
    entry = json.loads(FIXTURE.read_text())
    seq = entry["sequence"]["value"]
    assert checksum_lookup(seq, taxid=RHO1_TAXID) == RHO1_ACCESSION


@pytest.mark.network
def test_isoform_suffixed_accession_resolves_to_canonical_features(tmp_path):
    """An isoform-suffixed accession (e.g. 'G5EE56-1') must resolve to the base
    accession's full feature set, not the isoform entry's near-empty one.

    Regression test: UniProt's uniprotkb/{acc}.json endpoint accepts an
    isoform-suffixed accession (200 OK) but silently returns an incomplete/empty
    feature list for it -- src-1 (G5EE56-1) reported 0 UniProt features before
    this was fixed to strip the isoform suffix before fetching.
    """
    out = tmp_path / "out.txt"
    fasta = tmp_path / "in.fa"
    fasta.write_text(">G5EE56-1\nM\n")

    uniprot_features_main(str(fasta), str(out), str(FEATURES_FILE))

    rows = out.read_text().splitlines()
    assert any("Lipidation" in r and r.startswith("UniProt\t2\t2") for r in rows), rows


@pytest.mark.network
def test_committed_fixture_matches_live_blocking_features():
    """The committed Q22038 fixture hasn't drifted from live UniProt.

    If this fails, refresh tests/data/Q22038.uniprot.json from the live entry.
    """
    feature_config = load_feature_config(str(FEATURES_FILE))
    live_entry = fetch_entry(RHO1_ACCESSION)
    fixture_entry = json.loads(FIXTURE.read_text())

    live_blocking = {r for r in entry_to_feature_rows(live_entry, feature_config) if r[0] == "UniProt"}
    fixture_blocking = {r for r in entry_to_feature_rows(fixture_entry, feature_config) if r[0] == "UniProt"}
    assert live_blocking == fixture_blocking
