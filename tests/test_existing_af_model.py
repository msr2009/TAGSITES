"""
test_existing_af_model.py — offline unit tests for scripts/existing_AF_model.py

Tests cover checksum_lookup_response_to_accession, blast_tsv_line_to_afdb_hit,
and is_afdb_not_found — the pure "raw EBI/UniProt response -> internal shape"
parsers used by search_AFDB.
"""

from existing_AF_model import (
    checksum_lookup_response_to_accession,
    blast_tsv_line_to_afdb_hit,
    is_afdb_not_found,
)


# ── checksum_lookup_response_to_accession ────────────────────────────────────

class TestChecksumLookupResponseToAccession:

    def test_returns_first_accession_when_present(self):
        resp = {"results": [{"primaryAccession": "P12345"}, {"primaryAccession": "Q99999"}]}
        assert checksum_lookup_response_to_accession(resp) == "P12345"

    def test_returns_none_when_no_results(self):
        assert checksum_lookup_response_to_accession({"results": []}) is None

    def test_returns_none_when_results_key_absent(self):
        assert checksum_lookup_response_to_accession({}) is None


# ── blast_tsv_line_to_afdb_hit ────────────────────────────────────────────────

def _blast_tsv_line(accession="P12345", percent_id="98.5", evalue="1e-100"):
    # NCBIBLAST TSV columns (0-indexed): 0=query, 1=id, 2=acc, ... 7=%id ... 9=evalue
    fields = ["query", "id"] + [accession] + ["_"] * 4 + [percent_id, "_", evalue]
    return "\t".join(fields) + "\n"


class TestBlastTsvLineToAfdbHit:

    def test_parses_accession_percent_id_evalue(self):
        hit = blast_tsv_line_to_afdb_hit(_blast_tsv_line("Q8I4D9", "95.5", "1e-50"))
        assert hit["accession"] == "Q8I4D9"
        assert hit["percent_id"] == 95.5
        assert hit["evalue"] == 1e-50

    def test_percent_id_and_evalue_are_floats(self):
        hit = blast_tsv_line_to_afdb_hit(_blast_tsv_line())
        assert isinstance(hit["percent_id"], float)
        assert isinstance(hit["evalue"], float)

    def test_trailing_newline_stripped(self):
        hit = blast_tsv_line_to_afdb_hit(_blast_tsv_line(accession="P00000") + "\n")
        assert hit["accession"] == "P00000"


# ── is_afdb_not_found ──────────────────────────────────────────────────────────

class TestIsAfdbNotFound:

    def test_error_message_detected(self):
        assert is_afdb_not_found("ERROR 12 No entries found.") is True

    def test_valid_pdb_text_not_flagged(self):
        pdb_text = "HEADER    STRUCTURE\nATOM      1  N   MET A   1\n"
        assert is_afdb_not_found(pdb_text) is False

    def test_error_only_checked_in_first_50_chars(self):
        long_prefix = "X" * 60
        assert is_afdb_not_found(long_prefix + "ERROR") is False
