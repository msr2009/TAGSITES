"""
test_get_species_taxonomy.py — offline unit tests for scripts/get_species_taxonomy.py

Tests cover tax_to_dict, the pure text→dict parser.
"""

import pytest
from get_species_taxonomy import tax_to_dict


# A sample dbfetch taxonomy response for C. elegans
_ELEGANS_TAX = (
    "ID: 6239\n"
    "PARENT ID: 6237\n"
    "RANK: species\n"
    "SCIENTIFIC NAME: Caenorhabditis elegans\n"
    "//\n"
)

# A sample with extra whitespace around values
_WHITESPACE_TAX = (
    "ID:   9606  \n"
    "PARENT ID:   9605  \n"
    "RANK:   species  \n"
    "SCIENTIFIC NAME:   Homo sapiens  \n"
    "//\n"
)


class TestTaxToDict:

    def test_id_field(self):
        d = tax_to_dict(_ELEGANS_TAX)
        assert d["ID"] == "6239"

    def test_parent_id_field(self):
        d = tax_to_dict(_ELEGANS_TAX)
        assert d["PARENT ID"] == "6237"

    def test_rank_field(self):
        d = tax_to_dict(_ELEGANS_TAX)
        assert d["RANK"] == "species"

    def test_scientific_name_field(self):
        d = tax_to_dict(_ELEGANS_TAX)
        assert d["SCIENTIFIC NAME"] == "Caenorhabditis elegans"

    def test_returns_dict(self):
        assert isinstance(tax_to_dict(_ELEGANS_TAX), dict)

    def test_trailing_separator_stripped(self):
        """The '//' dbfetch record separator should not appear as a key."""
        d = tax_to_dict(_ELEGANS_TAX)
        assert "//" not in d

    def test_whitespace_trimmed_in_keys_and_values(self):
        d = tax_to_dict(_WHITESPACE_TAX)
        assert d["ID"] == "9606"
        assert d["SCIENTIFIC NAME"] == "Homo sapiens"

    def test_homo_sapiens_parent(self):
        d = tax_to_dict(_WHITESPACE_TAX)
        assert d["PARENT ID"] == "9605"
