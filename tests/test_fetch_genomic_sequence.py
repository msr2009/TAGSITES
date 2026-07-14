"""
test_fetch_genomic_sequence.py — offline unit tests for scripts/fetch_genomic_sequence.py

Mocks ensembl_rest's functions directly (no network) to drive the
organism-not-found / gene-not-found / success paths.
"""

import pytest

import ensembl_rest
from fetch_genomic_sequence import fetch_genomic_sequence


def test_raises_when_organism_not_on_ensembl(monkeypatch):
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: None)
    with pytest.raises(ValueError, match="not found on Ensembl"):
        fetch_genomic_sequence(taxid=1234567, gene_symbol="foo")


def test_raises_when_gene_symbol_not_found(monkeypatch):
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug",
                        lambda taxid: "caenorhabditis_elegans")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [])
    with pytest.raises(ValueError, match="No gene found"):
        fetch_genomic_sequence(taxid=6239, gene_symbol="nonexistent-gene")


def test_success_path_returns_fasta_and_meta(monkeypatch):
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug",
                        lambda taxid: "caenorhabditis_elegans")
    monkeypatch.setattr(ensembl_rest, "xref_symbol",
                        lambda species, symbol: [{"type": "gene", "id": "WBGene00006757"}])
    monkeypatch.setattr(ensembl_rest, "lookup_id", lambda gene_id: {
        "seq_region_name": "X", "start": 7682896, "end": 7686037,
        "strand": 1, "assembly_name": "WBcel235",
    })

    captured = {}

    def fake_fetch_region_fasta(species, seq_region, start, end, strand,
                                expand_5prime=0, expand_3prime=0):
        captured.update(locals())
        return ">chromosome:WBcel235:X:7681896:7687037:1\nACGT\n"

    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", fake_fetch_region_fasta)

    fasta_text, meta = fetch_genomic_sequence(taxid=6239, gene_symbol="unc-18", flank_bp=1000)

    assert fasta_text == ">chromosome:WBcel235:X:7681896:7687037:1\nACGT\n"
    assert meta == {
        "species": "caenorhabditis_elegans",
        "gene_id": "WBGene00006757",
        "assembly": "WBcel235",
        "region": "X:7682896-7686037:1",
        "ambiguous": False,
        "candidates": [],
    }
    assert captured["expand_5prime"] == 1000
    assert captured["expand_3prime"] == 1000


def test_success_path_ecoli_assembly_qualified_slug(monkeypatch):
    """The E. coli fast-path slug is assembly-qualified; confirm it flows through unchanged."""
    ecoli_slug = "escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: ecoli_slug)
    monkeypatch.setattr(ensembl_rest, "xref_symbol",
                        lambda species, symbol: [{"type": "gene", "id": "b0002"}])
    monkeypatch.setattr(ensembl_rest, "lookup_id", lambda gene_id: {
        "seq_region_name": "Chromosome", "start": 337, "end": 2799,
        "strand": 1, "assembly_name": "ASM584v2",
    })
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta",
                        lambda *a, **k: ">Chromosome\nATG...\n")

    fasta_text, meta = fetch_genomic_sequence(taxid=562, gene_symbol="thrA")

    assert meta["species"] == ecoli_slug
    assert meta["gene_id"] == "b0002"
    assert meta["assembly"] == "ASM584v2"
    assert meta["ambiguous"] is False


def test_lrg_duplicate_filtered_without_extra_lookup_call(monkeypatch):
    """The common human case (ENSG + LRG_) resolves via gene_candidates alone —
    _pick_most_complete_gene's span comparison should never kick in."""
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: "homo_sapiens")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [
        {"type": "gene", "id": "ENSG00000012048"}, {"type": "gene", "id": "LRG_292"},
    ])

    lookup_calls = []

    def fake_lookup_id(gene_id):
        lookup_calls.append(gene_id)
        return {"seq_region_name": "17", "start": 43044292, "end": 43170245,
                "strand": -1, "assembly_name": "GRCh38"}

    monkeypatch.setattr(ensembl_rest, "lookup_id", fake_lookup_id)
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", lambda *a, **k: ">seq\nACGT\n")

    fasta_text, meta = fetch_genomic_sequence(taxid=9606, gene_symbol="BRCA1")

    assert meta["gene_id"] == "ENSG00000012048"
    assert meta["ambiguous"] is False
    assert lookup_calls == ["ENSG00000012048"]   # only the chosen candidate, once


def test_duplicate_annotation_prefers_longest_genomic_span(monkeypatch):
    """Two non-LRG gene ids (e.g. zebrafish tp53 primary + a shorter duplicate
    annotation) — the longer (more complete) span should be chosen, and the
    ambiguity surfaced to the caller."""
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: "danio_rerio")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [
        {"type": "gene", "id": "ENSDARG00000035559"},
        {"type": "gene", "id": "ENSDARG00000115148"},
    ])

    def fake_lookup_id(gene_id):
        if gene_id == "ENSDARG00000035559":
            # real span: 11,579 bp (the complete gene model)
            return {"seq_region_name": "5", "start": 24086227, "end": 24097805,
                    "strand": 1, "assembly_name": "GRCz11"}
        # observed real duplicate span: 8,702 bp (shorter, incomplete)
        return {"seq_region_name": "ALT_CTG5_1_14", "start": 93172, "end": 101873,
                "strand": 1, "assembly_name": "GRCz11"}

    monkeypatch.setattr(ensembl_rest, "lookup_id", fake_lookup_id)
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", lambda *a, **k: ">seq\nACGT\n")

    fasta_text, meta = fetch_genomic_sequence(taxid=7955, gene_symbol="tp53")

    assert meta["gene_id"] == "ENSDARG00000035559"
    assert meta["ambiguous"] is True
    assert meta["candidates"] == ["ENSDARG00000035559", "ENSDARG00000115148"]


def test_equal_spans_falls_back_to_first_candidate(monkeypatch):
    """If every candidate has the same span, fall back to the first rather than picking arbitrarily."""
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: "danio_rerio")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [
        {"type": "gene", "id": "GENE_A"}, {"type": "gene", "id": "GENE_B"},
    ])
    monkeypatch.setattr(ensembl_rest, "lookup_id", lambda gene_id: {
        "seq_region_name": "5", "start": 1, "end": 100,
        "strand": 1, "assembly_name": "GRCz11",
    })
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", lambda *a, **k: ">seq\nACGT\n")

    fasta_text, meta = fetch_genomic_sequence(taxid=7955, gene_symbol="ambiguous-gene")

    assert meta["gene_id"] == "GENE_A"
    assert meta["ambiguous"] is True


def test_disambiguation_reuses_lookup_result_no_refetch(monkeypatch):
    """The winning candidate's coords (fetched while comparing spans) should be
    reused, not re-fetched via a second lookup_id call."""
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: "danio_rerio")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [
        {"type": "gene", "id": "GENE_A"}, {"type": "gene", "id": "GENE_B"},
    ])

    lookup_calls = []

    def fake_lookup_id(gene_id):
        lookup_calls.append(gene_id)
        span = {"GENE_A": (1, 1000), "GENE_B": (1, 100)}[gene_id]
        return {"seq_region_name": "5", "start": span[0], "end": span[1],
                "strand": 1, "assembly_name": "GRCz11"}

    monkeypatch.setattr(ensembl_rest, "lookup_id", fake_lookup_id)
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", lambda *a, **k: ">seq\nACGT\n")

    fetch_genomic_sequence(taxid=7955, gene_symbol="ambiguous-gene")

    assert lookup_calls == ["GENE_A", "GENE_B"]   # one lookup per candidate, no re-fetch of the winner


def test_disambiguation_falls_back_when_all_lookups_fail(monkeypatch):
    """If every candidate's lookup_id call fails, fall back to the first candidate
    and fetch its coords fresh rather than raising."""
    monkeypatch.setattr(ensembl_rest, "resolve_species_slug", lambda taxid: "danio_rerio")
    monkeypatch.setattr(ensembl_rest, "xref_symbol", lambda species, symbol: [
        {"type": "gene", "id": "GENE_A"}, {"type": "gene", "id": "GENE_B"},
    ])

    calls = {"n": 0}

    def fake_lookup_id(gene_id):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise Exception("boom")
        return {"seq_region_name": "5", "start": 1, "end": 100,
                "strand": 1, "assembly_name": "GRCz11"}

    monkeypatch.setattr(ensembl_rest, "lookup_id", fake_lookup_id)
    monkeypatch.setattr(ensembl_rest, "fetch_region_fasta", lambda *a, **k: ">seq\nACGT\n")

    fasta_text, meta = fetch_genomic_sequence(taxid=7955, gene_symbol="ambiguous-gene")

    assert meta["gene_id"] == "GENE_A"
