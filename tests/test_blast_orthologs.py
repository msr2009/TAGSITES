"""
test_blast_orthologs.py — offline unit tests for scripts/blast_orthologs.py

Tests cover hit_to_dict, the pure JSON→dict converter that is the only
offline-testable function in this script.
"""

import pytest
from blast_orthologs import hit_to_dict


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
