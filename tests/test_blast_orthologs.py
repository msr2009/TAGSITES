"""
test_blast_orthologs.py — offline unit tests for scripts/blast_orthologs.py

Tests cover hit_to_dict (pure JSON->dict converter), group_hits_by_species
(evalue/length filtering + per-species grouping), parse_dbfetch_fasta_record
(raw dbfetch FASTA -> header/seq), and ensure_query_in_alignment_set (the
Clustal Omega input-preparation step) — all offline, no network required.
"""

import pytest
from blast_orthologs import (
    hit_to_dict,
    group_hits_by_species,
    parse_dbfetch_fasta_record,
    ensure_query_in_alignment_set,
)


# ── fixture hit dict ──────────────────────────────────────────────────────────

def _make_hit(acc="P12345", species="Caenorhabditis elegans",
              expect="1e-50", identity="98.0",
              hit_from="10", hit_to="210", hseq="MKVL---ACDE"):
    return {
        "hit_acc": acc,
        "hit_os": species,
        "hit_hsps": [{
            "hsp_expect": expect,
            "hsp_identity": identity,
            "hsp_hit_from": hit_from,
            "hsp_hit_to": hit_to,
            "hsp_hseq": hseq,
        }]
    }


# ── hit_to_dict ───────────────────────────────────────────────────────────────

class TestHitToDict:

    def test_acc_passthrough(self):
        d = hit_to_dict(_make_hit(acc="Q8I4D9"))
        assert d["acc"] == "Q8I4D9"

    def test_species_passthrough(self):
        d = hit_to_dict(_make_hit(species="Homo sapiens"))
        assert d["species"] == "Homo sapiens"

    def test_evalue_is_float(self):
        d = hit_to_dict(_make_hit(expect="1e-50"))
        assert isinstance(d["evalue"], float)
        assert d["evalue"] == pytest.approx(1e-50)

    def test_identity_is_float(self):
        d = hit_to_dict(_make_hit(identity="95.5"))
        assert isinstance(d["identity"], float)
        assert d["identity"] == pytest.approx(95.5)

    def test_length_is_raw_difference(self):
        """length = hit_to - hit_from (no +1)."""
        d = hit_to_dict(_make_hit(hit_from="10", hit_to="210"))
        assert d["length"] == 200

    def test_hitseq_gaps_removed(self):
        """Dashes in hsp_hseq are stripped from hitseq."""
        d = hit_to_dict(_make_hit(hseq="MKVL---ACDE"))
        assert d["hitseq"] == "MKVLACDE"
        assert "-" not in d["hitseq"]

    def test_hitseq_no_gaps_unchanged(self):
        d = hit_to_dict(_make_hit(hseq="MKVLACDE"))
        assert d["hitseq"] == "MKVLACDE"

    def test_returns_dict(self):
        d = hit_to_dict(_make_hit())
        assert isinstance(d, dict)

    def test_expected_keys_present(self):
        d = hit_to_dict(_make_hit())
        for key in ["acc", "species", "evalue", "identity", "length", "hitseq"]:
            assert key in d, f"key {key!r} missing from hit_to_dict output"


# ── group_hits_by_species ───────────────────────────────────────────────────────

def _raw_hit(acc, species, expect="1e-50", hit_from="1", hit_to="101"):
    return {
        "hit_acc": acc,
        "hit_os": species,
        "hit_hsps": [{
            "hsp_expect": expect,
            "hsp_identity": "90.0",
            "hsp_hit_from": hit_from,
            "hsp_hit_to": hit_to,
            "hsp_hseq": "M" * (int(hit_to) - int(hit_from)),
        }],
    }


class TestGroupHitsBySpecies:

    def test_groups_by_species(self):
        hits = [_raw_hit("A1", "worm"), _raw_hit("A2", "human"), _raw_hit("A3", "worm")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-10,
                                        length_percent=0, exclude_paralogs=False, n=10)
        assert len(grouped["worm"]) == 2
        assert len(grouped["human"]) == 1

    def test_query_species_key_always_present(self):
        grouped = group_hits_by_species([], "worm", 100, evalue=1e-10,
                                        length_percent=0, exclude_paralogs=False, n=10)
        assert grouped == {"worm": []}

    def test_evalue_filter_excludes_worse_hits(self):
        hits = [_raw_hit("A1", "human", expect="1e-5")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-10,
                                        length_percent=0.5, exclude_paralogs=False, n=10)
        assert "human" not in grouped

    def test_length_filter_excludes_short_hits(self):
        # hit length 20/100 = 0.2, below length_percent=0.5 floor
        hits = [_raw_hit("A1", "human", hit_from="0", hit_to="20")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0.5, exclude_paralogs=False, n=10)
        assert "human" not in grouped

    def test_length_filter_disabled_when_zero(self):
        hits = [_raw_hit("A1", "human", hit_from="0", hit_to="20")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0, exclude_paralogs=False, n=10)
        assert "human" in grouped

    def test_exclude_paralogs_keeps_best_hit_per_species(self):
        hits = [
            _raw_hit("A1", "human", expect="1e-20"),
            _raw_hit("A2", "human", expect="1e-50"),  # better (lower) evalue
        ]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0, exclude_paralogs=True, n=10)
        assert len(grouped["human"]) == 1
        assert grouped["human"][0]["acc"] == "A2"

    def test_exclude_paralogs_worse_hit_does_not_replace_best(self):
        hits = [
            _raw_hit("A1", "human", expect="1e-50"),
            _raw_hit("A2", "human", expect="1e-20"),  # worse; should not replace A1
        ]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0, exclude_paralogs=True, n=10)
        assert grouped["human"][0]["acc"] == "A1"

    def test_exclude_paralogs_keeps_all_query_species_hits(self):
        hits = [_raw_hit("A1", "worm"), _raw_hit("A2", "worm")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0, exclude_paralogs=True, n=10)
        assert len(grouped["worm"]) == 2

    def test_stops_once_n_species_reached(self):
        hits = [_raw_hit("A1", "human"), _raw_hit("A2", "mouse"), _raw_hit("A3", "fly")]
        grouped = group_hits_by_species(hits, "worm", 100, evalue=1e-1,
                                        length_percent=0, exclude_paralogs=False, n=2)
        assert len(grouped) == 2


# ── parse_dbfetch_fasta_record ────────────────────────────────────────────────

class TestParseDbfetchFastaRecord:

    def test_extracts_header_word_and_sequence(self):
        raw = ">sp|P12345|SNB1_CAEEL Synaptobrevin\nMDAQGDAGAQ\nGGSQGGPRPS\n"
        name, seq = parse_dbfetch_fasta_record(raw)
        assert name == ">sp|P12345|SNB1_CAEEL"
        assert seq == "MDAQGDAGAQGGSQGGPRPS"

    def test_single_line_sequence(self):
        name, seq = parse_dbfetch_fasta_record(">acc1 desc\nMKVL\n")
        assert name == ">acc1"
        assert seq == "MKVL"


# ── ensure_query_in_alignment_set ─────────────────────────────────────────────

class TestEnsureQueryInAlignmentSet:

    def test_inserts_query_when_no_match(self):
        fasta_list = [">hit1\nAAAA\n", ">hit2\nBBBB\n"]
        result_list, input_match = ensure_query_in_alignment_set(
            fasta_list, "", "myquery", "MKVL")
        assert result_list[0] == ">myquery\nMKVL\n"
        assert input_match == "myquery"

    def test_replaces_first_element_only(self):
        fasta_list = [">hit1\nAAAA\n", ">hit2\nBBBB\n"]
        result_list, _ = ensure_query_in_alignment_set(fasta_list, "", "myquery", "MKVL")
        assert result_list[1] == ">hit2\nBBBB\n"
        assert len(result_list) == 2

    def test_no_change_when_match_already_found(self):
        fasta_list = [">hit1\nAAAA\n", ">hit2\nBBBB\n"]
        result_list, input_match = ensure_query_in_alignment_set(
            fasta_list, "hit1", "myquery", "MKVL")
        assert result_list == fasta_list
        assert input_match == "hit1"

    def test_empty_list_with_no_match(self):
        result_list, input_match = ensure_query_in_alignment_set([], "", "myquery", "MKVL")
        assert result_list == [">myquery\nMKVL\n"]
        assert input_match == "myquery"
