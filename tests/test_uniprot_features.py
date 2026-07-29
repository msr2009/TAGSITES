"""
test_uniprot_features.py — offline unit tests for scripts/uniprot_features.py

Tests cover load_feature_config and entry_to_feature_rows against a committed
UniProt entry fixture (tests/data/Q22038.uniprot.json — C. elegans rho-1), so
no network access is required. See tests/test_uniprot_network.py for the
live-API contract guards that this fixture must stay in sync with.
"""

import json
from pathlib import Path

import pytest

from site_selection_util import uniprot_accession_regex
import uniprot_features
from uniprot_features import load_feature_config, entry_to_feature_rows, main as uniprot_features_main

REPO_ROOT = Path(__file__).parent.parent
FEATURES_FILE = REPO_ROOT / "tables" / "uniprot_features.txt"
FIXTURE = REPO_ROOT / "tests" / "data" / "Q22038.uniprot.json"


@pytest.fixture(scope="module")
def feature_config():
    return load_feature_config(str(FEATURES_FILE))


@pytest.fixture(scope="module")
def q22038():
    return json.loads(FIXTURE.read_text())


class TestLoadFeatureConfig:

    def test_known_types_present(self, feature_config):
        assert feature_config["Lipidation"] == ("blocking", "Lipidation")
        assert feature_config["Binding site"] == ("site", "Binding site")

    def test_blank_lines_and_comments_ignored(self, tmp_path):
        f = tmp_path / "features.txt"
        f.write_text("# comment\n\nLipidation\tblocking\tLipidation\n")
        config = load_feature_config(str(f))
        assert config == {"Lipidation": ("blocking", "Lipidation")}

    def test_short_lines_skipped(self, tmp_path):
        f = tmp_path / "features.txt"
        f.write_text("Lipidation\tblocking\n")  # missing label field
        assert load_feature_config(str(f)) == {}


class TestEntryToFeatureRows:

    def test_lipidation_is_blocking(self, q22038, feature_config):
        rows = entry_to_feature_rows(q22038, feature_config)
        # evidence is ECO:0000250/UniProtKB (similarity-inferred), not a PubMed citation
        assert ("UniProt", 189, 189, "Lipidation: S-geranylgeranyl cysteine (predicted)") in rows

    def test_propeptide_is_blocking(self, q22038, feature_config):
        rows = entry_to_feature_rows(q22038, feature_config)
        assert ("UniProt", 190, 192, "Propeptide: Removed in mature form (predicted)") in rows

    def test_binding_site_and_motif_are_informational(self, q22038, feature_config):
        rows = entry_to_feature_rows(q22038, feature_config)
        sources = {(desc.split(":")[0], source) for source, _, _, desc in rows}
        assert ("Binding site", "UniProt_site") in sources
        assert ("Motif", "UniProt_site") in sources
        assert not any(source == "UniProt" and desc.startswith(("Binding site", "Motif"))
                       for source, _, _, desc in rows)

    def test_binding_site_uses_ligand_name_when_description_blank(self, q22038, feature_config):
        rows = entry_to_feature_rows(q22038, feature_config)
        assert any(source == "UniProt_site" and start == 15 and "GTP" in desc
                   for source, start, _, desc in rows)

    def test_pubmed_evidence_appears_in_suffix(self, feature_config):
        entry = {"features": [{"type": "Modified residue",
                                "location": {"start": {"value": 10}, "end": {"value": 10}},
                                "description": "Phosphotyrosine",
                                "evidences": [{"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "123"}]}]}
        rows = entry_to_feature_rows(entry, feature_config)
        assert rows == [("UniProt_site", 10, 10, "Modified residue: Phosphotyrosine (PubMed:123)")]

    def test_no_evidence_at_all_has_no_suffix(self, feature_config):
        entry = {"features": [{"type": "Motif",
                                "location": {"start": {"value": 10}, "end": {"value": 10}},
                                "description": "Foo"}]}
        rows = entry_to_feature_rows(entry, feature_config)
        assert rows == [("UniProt_site", 10, 10, "Motif: Foo")]

    def test_unlisted_types_dropped(self, q22038, feature_config):
        rows = entry_to_feature_rows(q22038, feature_config)
        assert not any("Chain" in desc for _, _, _, desc in rows)

    def test_missing_description_falls_back_to_type(self, feature_config):
        entry = {"features": [{"type": "Lipidation",
                                "location": {"start": {"value": 5}, "end": {"value": 5}}}]}
        rows = entry_to_feature_rows(entry, feature_config)
        assert rows == [("UniProt", 5, 5, "Lipidation: Lipidation")]

    def test_empty_features_list_yields_no_rows(self, feature_config):
        assert entry_to_feature_rows({"features": []}, feature_config) == []

    def test_missing_features_key_yields_no_rows(self, feature_config):
        assert entry_to_feature_rows({}, feature_config) == []


class TestMainDegradesGracefullyOnFetchFailure:
    """A bad/rejected accession must not crash the task -- it should degrade to
    an empty feature file like the no-accession-found case does, not propagate
    fetch_entry's HTTPError uncaught. Regression for a bad --accession (e.g. a
    typo) or a stale/deleted accession returned by checksum lookup.
    """

    def test_fetch_entry_exception_writes_empty_file_and_does_not_raise(self, tmp_path, monkeypatch):
        def _boom(accession):
            raise Exception("simulated 400 Bad Request")
        monkeypatch.setattr(uniprot_features, "fetch_entry", _boom)

        fasta = tmp_path / "in.fa"
        fasta.write_text(">x\nM\n")
        out = tmp_path / "out.txt"

        uniprot_features_main(str(fasta), str(out), str(FEATURES_FILE), accession="NOTREAL123")

        assert out.exists()
        assert out.read_text() == ""


class TestAccessionRegexStripsIsoformSuffix:
    """uniprot_accession_regex uses re.match, which only anchors the start --
    callers must use match.group(0), not the original string, or an
    isoform-suffixed accession like 'G5EE56-1' leaks through unstripped.
    See tests/test_uniprot_network.py::test_isoform_suffixed_accession_resolves_to_canonical_features
    for the live regression this guards against.
    """

    def test_isoform_suffix_not_included_in_match(self):
        m = uniprot_accession_regex("G5EE56-1")
        assert m is not None
        assert m.group(0) == "G5EE56"

    def test_plain_accession_matches_in_full(self):
        m = uniprot_accession_regex("Q22038")
        assert m.group(0) == "Q22038"
