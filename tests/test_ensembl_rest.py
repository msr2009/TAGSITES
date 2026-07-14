"""
test_ensembl_rest.py — offline unit tests for scripts/ensembl_rest.py

Covers: xref_symbol/lookup_id/fetch_region_fasta/list_species (requests
mocked), first_gene_id (pure), resolve_species_slug (fast-path + dynamic
division-search fallback with caching).
"""

from unittest.mock import Mock

import pytest

import ensembl_rest


def _fake_response(json_data=None, text=None, status_ok=True):
    resp = Mock()
    resp.text = text
    if json_data is not None:
        resp.json = Mock(return_value=json_data)
    if status_ok:
        resp.raise_for_status = Mock()
    else:
        resp.raise_for_status = Mock(side_effect=ensembl_rest.requests.HTTPError("bad status"))
    return resp


@pytest.fixture(autouse=True)
def _clear_division_cache():
    """Each test gets a clean division cache so caching behavior is testable in isolation."""
    ensembl_rest._division_cache.clear()
    yield
    ensembl_rest._division_cache.clear()


# ── xref_symbol / lookup_id / fetch_region_fasta / list_species ───────────────

class TestXrefSymbol:

    def test_builds_correct_url_and_returns_json(self, monkeypatch):
        captured = {}

        def fake_get(url, params=None, timeout=None):
            captured["url"] = url
            captured["params"] = params
            return _fake_response(json_data=[{"type": "gene", "id": "WBGene00006757"}])

        monkeypatch.setattr(ensembl_rest.requests, "get", fake_get)
        result = ensembl_rest.xref_symbol("caenorhabditis_elegans", "unc-18")
        assert captured["url"] == "https://rest.ensembl.org/xrefs/symbol/caenorhabditis_elegans/unc-18"
        assert result == [{"type": "gene", "id": "WBGene00006757"}]

    def test_raises_on_bad_status(self, monkeypatch):
        monkeypatch.setattr(ensembl_rest.requests, "get",
                            lambda *a, **k: _fake_response(json_data=[], status_ok=False))
        with pytest.raises(Exception):
            ensembl_rest.xref_symbol("homo_sapiens", "NOPE")


class TestLookupId:

    def test_builds_correct_url(self, monkeypatch):
        captured = {}

        def fake_get(url, params=None, timeout=None):
            captured["url"] = url
            return _fake_response(json_data={"seq_region_name": "X", "start": 1, "end": 100,
                                             "strand": 1, "assembly_name": "WBcel235"})

        monkeypatch.setattr(ensembl_rest.requests, "get", fake_get)
        result = ensembl_rest.lookup_id("WBGene00006757")
        assert captured["url"] == "https://rest.ensembl.org/lookup/id/WBGene00006757"
        assert result["assembly_name"] == "WBcel235"


class TestFetchRegionFasta:

    def test_builds_correct_url_and_params(self, monkeypatch):
        captured = {}

        def fake_get(url, params=None, timeout=None):
            captured["url"] = url
            captured["params"] = params
            return _fake_response(text=">seq\nACGT\n")

        monkeypatch.setattr(ensembl_rest.requests, "get", fake_get)
        result = ensembl_rest.fetch_region_fasta(
            "caenorhabditis_elegans", "X", 100, 200, 1, expand_5prime=50, expand_3prime=60
        )
        assert captured["url"] == "https://rest.ensembl.org/sequence/region/caenorhabditis_elegans/X:100-200:1"
        assert captured["params"]["expand_5prime"] == 50
        assert captured["params"]["expand_3prime"] == 60
        assert result == ">seq\nACGT\n"


class TestListSpecies:

    def test_no_division_filter(self, monkeypatch):
        captured = {}

        def fake_get(url, params=None, timeout=None):
            captured["url"] = url
            captured["params"] = params
            return _fake_response(json_data={"species": [{"name": "homo_sapiens", "taxon_id": "9606"}]})

        monkeypatch.setattr(ensembl_rest.requests, "get", fake_get)
        result = ensembl_rest.list_species()
        assert "division" not in captured["params"]
        assert result == [{"name": "homo_sapiens", "taxon_id": "9606"}]

    def test_division_filter_passed_through(self, monkeypatch):
        captured = {}

        def fake_get(url, params=None, timeout=None):
            captured["params"] = params
            return _fake_response(json_data={"species": []})

        monkeypatch.setattr(ensembl_rest.requests, "get", fake_get)
        ensembl_rest.list_species(division="EnsemblBacteria")
        assert captured["params"]["division"] == "EnsemblBacteria"

    def test_raises_on_bad_status(self, monkeypatch):
        monkeypatch.setattr(ensembl_rest.requests, "get",
                            lambda *a, **k: _fake_response(json_data={}, status_ok=False))
        with pytest.raises(Exception):
            ensembl_rest.list_species()


# ── first_gene_id ──────────────────────────────────────────────────────────────

class TestFirstGeneId:

    def test_returns_first_gene_entry(self):
        xrefs = [{"type": "gene", "id": "ENSG00000012048"}, {"type": "gene", "id": "LRG_292"}]
        assert ensembl_rest.first_gene_id(xrefs) == "ENSG00000012048"

    def test_skips_non_gene_entries(self):
        xrefs = [{"type": "transcript", "id": "T1"}, {"type": "gene", "id": "G1"}]
        assert ensembl_rest.first_gene_id(xrefs) == "G1"

    def test_returns_none_when_empty(self):
        assert ensembl_rest.first_gene_id([]) is None

    def test_returns_none_when_no_gene_type(self):
        assert ensembl_rest.first_gene_id([{"type": "transcript", "id": "T1"}]) is None


# ── gene_candidates ─────────────────────────────────────────────────────────────

class TestGeneCandidates:

    def test_filters_out_lrg_duplicates(self):
        xrefs = [{"type": "gene", "id": "ENSG00000012048"}, {"type": "gene", "id": "LRG_292"}]
        assert ensembl_rest.gene_candidates(xrefs) == ["ENSG00000012048"]

    def test_keeps_multiple_non_lrg_candidates(self):
        xrefs = [{"type": "gene", "id": "ENSDARG00000035559"},
                 {"type": "gene", "id": "ENSDARG00000115148"}]
        assert ensembl_rest.gene_candidates(xrefs) == ["ENSDARG00000035559", "ENSDARG00000115148"]

    def test_falls_back_to_lrg_only_when_nothing_else_matches(self):
        """If every match happens to be an LRG_ id, don't discard everything."""
        xrefs = [{"type": "gene", "id": "LRG_292"}]
        assert ensembl_rest.gene_candidates(xrefs) == ["LRG_292"]

    def test_returns_empty_list_when_no_gene_type(self):
        assert ensembl_rest.gene_candidates([{"type": "transcript", "id": "T1"}]) == []


# ── resolve_species_slug ────────────────────────────────────────────────────────

class TestResolveSpeciesSlug:

    def test_fast_path_makes_no_network_call(self, monkeypatch):
        called = []
        monkeypatch.setattr(ensembl_rest.requests, "get", lambda *a, **k: called.append(1))
        result = ensembl_rest.resolve_species_slug(6239)
        assert result == "caenorhabditis_elegans"
        assert called == []

    def test_fast_path_accepts_string_taxid(self, monkeypatch):
        called = []
        monkeypatch.setattr(ensembl_rest.requests, "get", lambda *a, **k: called.append(1))
        assert ensembl_rest.resolve_species_slug("9606") == "homo_sapiens"
        assert called == []

    def test_dynamic_fallback_searches_divisions_in_order(self, monkeypatch):
        seen_divisions = []

        def fake_list_species(division=None):
            seen_divisions.append(division)
            if division == "EnsemblMetazoa":
                return [{"name": "some_metazoan", "taxon_id": "99999"}]
            return []

        monkeypatch.setattr(ensembl_rest, "list_species", fake_list_species)
        result = ensembl_rest.resolve_species_slug(99999)
        assert result == "some_metazoan"
        assert seen_divisions == ["EnsemblVertebrates", "EnsemblMetazoa"]

    def test_dynamic_fallback_caches_division_listing(self, monkeypatch):
        call_count = {"n": 0}

        def fake_list_species(division=None):
            call_count["n"] += 1
            return [{"name": "some_vertebrate", "taxon_id": "12345"}]

        monkeypatch.setattr(ensembl_rest, "list_species", fake_list_species)
        ensembl_rest.resolve_species_slug(12345)
        ensembl_rest.resolve_species_slug(12345)
        assert call_count["n"] == 1

    def test_returns_none_when_not_found_in_any_division(self, monkeypatch):
        monkeypatch.setattr(ensembl_rest, "list_species", lambda division=None: [])
        assert ensembl_rest.resolve_species_slug(424242) is None

    def test_multiple_matches_returns_sorted_first(self, monkeypatch):
        def fake_list_species(division=None):
            if division == "EnsemblVertebrates":
                return [{"name": "zzz_species", "taxon_id": "555"},
                        {"name": "aaa_species", "taxon_id": "555"}]
            return []

        monkeypatch.setattr(ensembl_rest, "list_species", fake_list_species)
        assert ensembl_rest.resolve_species_slug(555) == "aaa_species"
