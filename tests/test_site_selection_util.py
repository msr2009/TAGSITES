"""
test_site_selection_util.py — offline unit tests for scripts/site_selection_util.py

Tests cover the pure/deterministic functions only (no network, no external binaries):
  sequence_type, translate_seq, three_to_one, remove_lowercase,
  check_input_type, uniprot_accession_regex, ncbiblast_call, read_fasta,
  resolve_taxids.
"""

import pytest
from Bio import SeqIO

from site_selection_util import (
    sequence_type,
    translate_seq,
    three_to_one,
    remove_lowercase,
    check_input_type,
    uniprot_accession_regex,
    ncbiblast_call,
    read_fasta,
    resolve_taxids,
)


# ── sequence_type ─────────────────────────────────────────────────────────────

class TestSequenceType:

    def test_dna_uppercase(self):
        assert sequence_type("ACTGNN") == "dna"

    def test_dna_lowercase(self):
        # sequence_type uppercases internally
        assert sequence_type("actg") == "dna"

    def test_protein(self):
        assert sequence_type("MAEQWRKLHDS") == "protein"

    def test_protein_with_stop(self):
        # trailing '*' is stripped internally
        assert sequence_type("MKVL*") == "protein"

    def test_unknown(self):
        assert sequence_type("123!@#") == "unknown"

    def test_shared_bases_are_dna(self):
        # ACGT alone → dna (subset of DNA bases)
        assert sequence_type("ACGT") == "dna"

    def test_protein_unique_aa(self):
        # E,F,H not in DNA set → protein
        assert sequence_type("AEFHIKLM") == "protein"


# ── translate_seq ─────────────────────────────────────────────────────────────

class TestTranslateSeq:

    def test_atg_met(self):
        assert translate_seq("ATG") == "M"

    def test_stop_underscore(self):
        # stop codons → '_'
        assert translate_seq("TAA") == "_"
        assert translate_seq("TAG") == "_"
        assert translate_seq("TGA") == "_"

    def test_unknown_codon(self):
        # NNN is not in the codon table → 'X'
        assert translate_seq("NNN") == "X"

    def test_full_codon_chain(self):
        # ATG-GGT-CCG → MGP
        assert translate_seq("ATGGGTCCG") == "MGP"

    def test_length_not_multiple_of_3(self):
        with pytest.raises(ValueError):
            translate_seq("ATGC")

    def test_empty_string(self):
        assert translate_seq("") == ""


# ── three_to_one ──────────────────────────────────────────────────────────────

class TestThreeToOne:

    def test_ala(self):
        assert three_to_one("ALA") == "A"

    def test_met(self):
        assert three_to_one("MET") == "M"

    def test_case_insensitive(self):
        assert three_to_one("ala") == "A"
        assert three_to_one("Ala") == "A"

    def test_unknown(self):
        assert three_to_one("XYZ") == "X"

    def test_all_standard_aas(self):
        mapping = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        for three, one in mapping.items():
            assert three_to_one(three) == one, f"three_to_one({three!r}) != {one!r}"


# ── remove_lowercase ──────────────────────────────────────────────────────────

class TestRemoveLowercase:

    def test_removes_lowercase(self):
        assert remove_lowercase("ATGGTaCCg") == "ATGGTATCC".replace("a", "").replace("g", "") \
            or remove_lowercase("ATGGTaCCg") == "ATGGTCC"

    def test_all_uppercase_unchanged(self):
        assert remove_lowercase("ATGGT") == "ATGGT"

    def test_all_lowercase_becomes_empty(self):
        assert remove_lowercase("actg") == ""

    def test_mixed(self):
        result = remove_lowercase("AaTtGg")
        assert result == "ATG"


# ── check_input_type ──────────────────────────────────────────────────────────

class TestCheckInputType:

    def test_pdb_extension(self):
        assert check_input_type("structure.pdb") == "pdb"

    def test_fa_extension(self):
        assert check_input_type("protein.fa") == "fasta"

    def test_fasta_extension(self):
        assert check_input_type("protein.fasta") == "fasta"

    def test_fsa_extension(self):
        assert check_input_type("protein.fsa") == "fasta"

    def test_uniprot_accession_o(self):
        # UniProt accession starting with O/P/Q
        assert check_input_type("Q8I4D9") == "uniprot"

    def test_uniprot_accession_p(self):
        assert check_input_type("P12345") == "uniprot"

    def test_unknown_extension(self):
        assert check_input_type("data.csv") == "err"

    def test_uniprot_wins_over_extension(self):
        # A uniprot accession that happens to have no extension
        # should be identified as uniprot, not err
        assert check_input_type("Q8I4D9") == "uniprot"


# ── uniprot_accession_regex ───────────────────────────────────────────────────

class TestUniprotAccessionRegex:

    def test_q_accession(self):
        assert uniprot_accession_regex("Q8I4D9") is not None

    def test_p_accession(self):
        assert uniprot_accession_regex("P12345") is not None

    def test_o_accession(self):
        assert uniprot_accession_regex("O15350") is not None

    def test_a_accession(self):
        assert uniprot_accession_regex("A0A000ABC1") is not None

    def test_non_accession_rejected(self):
        assert uniprot_accession_regex("notanacc") is None

    def test_filename_rejected(self):
        assert uniprot_accession_regex("protein.fa") is None


# ── ncbiblast_call ────────────────────────────────────────────────────────────

class TestNcbiblastCall:

    def test_taxid_1_omits_taxid_flag(self):
        cmd = ncbiblast_call("scripts/", "test@example.com", "MKVL",
                             "uniprotkb", 100, 1e-10, 1, "out")
        assert "--taxid" not in cmd

    def test_other_taxid_includes_taxid_flag(self):
        cmd = ncbiblast_call("scripts/", "test@example.com", "MKVL",
                             "uniprotkb", 100, 1e-10, 6239, "out")
        assert "--taxid 6239" in cmd

    def test_returns_string(self):
        cmd = ncbiblast_call("scripts/", "test@example.com", "MKVL",
                             "uniprotkb", 100, 1e-10, 1, "out")
        assert isinstance(cmd, str)

    def test_email_in_command(self):
        cmd = ncbiblast_call("scripts/", "me@example.com", "MKVL",
                             "uniprotkb", 100, 1e-10, 1, "out")
        assert "me@example.com" in cmd


# ── read_fasta ────────────────────────────────────────────────────────────────

class TestReadFasta:

    def test_returns_id_and_seq(self, DATA):
        name, seq = read_fasta(DATA / "snb-1.fa")
        assert name == "SNB-1"

    def test_seq_length(self, DATA):
        name, seq = read_fasta(DATA / "snb-1.fa")
        assert len(seq) == 109

    def test_seq_starts_with_M(self, DATA):
        name, seq = read_fasta(DATA / "snb-1.fa")
        assert str(seq)[0] == "M"

    def test_trailing_stop_stripped(self, DATA):
        # snb-1.fa does not have a trailing '*', but if it did it should be stripped
        name, seq = read_fasta(DATA / "snb-1.fa")
        assert not str(seq).endswith("*")


# ── resolve_taxids ────────────────────────────────────────────────────────────

class TestResolveTaxids:

    def test_manual_taxid_only(self):
        assert resolve_taxids("9606", None) == "9606"

    def test_manual_taxid_list(self):
        assert resolve_taxids("9606,10090", None) == "9606,10090"

    def test_sentinel_any_species_dropped(self):
        assert resolve_taxids("1", None) == ""
        assert resolve_taxids("1.0", None) == ""
        assert resolve_taxids("", None) == ""

    def test_no_taxid_no_file_returns_empty(self):
        assert resolve_taxids("", None) == ""

    def test_file_only(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("9606\n10090\n7955\n")
        assert resolve_taxids("", str(f)) == "9606,10090,7955"

    def test_file_ignores_blank_lines(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("9606\n\n10090\n\n")
        assert resolve_taxids("", str(f)) == "9606,10090"

    def test_file_strips_comments(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("9606\t# Homo sapiens\n10090  # Mus musculus\n# full-line comment\n")
        assert resolve_taxids("", str(f)) == "9606,10090"

    def test_manual_and_file_merged(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("10090\n7955\n")
        assert resolve_taxids("9606", str(f)) == "9606,10090,7955"

    def test_duplicates_deduped_preserving_first_occurrence(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("9606\n10090\n")
        assert resolve_taxids("10090,9606", str(f)) == "10090,9606"

    def test_manual_sentinel_dropped_but_file_kept(self, tmp_path):
        f = tmp_path / "taxids.txt"
        f.write_text("9606\n")
        assert resolve_taxids("1", str(f)) == "9606"

    def test_no_taxid_file_arg(self):
        # taxid_file="" (falsy) behaves like None -- no file read attempted
        assert resolve_taxids("9606", "") == "9606"
