"""
test_call_interpro.py — offline unit tests for scripts/call_interpro.py

Tests cover iprscan_tsv_to_domains, the pure InterProScan5-TSV parser.
"""

from call_interpro import iprscan_tsv_to_domains


# A realistic InterProScan5 TSV row (11-column schema):
# acc, md5, len, source, sig_acc, description, start, stop, score, status, date
def _row(source="Pfam", sig_acc="PF00856", description="SNARE domain",
         start="10", stop="90", score="1.2e-30"):
    return "\t".join([
        "SNB-1", "abc123md5", "109", source, sig_acc, description,
        start, stop, score, "T", "01-01-2026",
    ])


class TestIprscanTsvToDomains:

    def test_single_row_parsed(self):
        rows = iprscan_tsv_to_domains(_row())
        assert rows == [("Pfam", "10", "90", "SNARE domain")]

    def test_multiple_rows_preserve_order(self):
        text = "\n".join([
            _row(source="Pfam", start="10", stop="90"),
            _row(source="PANTHER", start="1", stop="109"),
        ])
        rows = iprscan_tsv_to_domains(text)
        assert rows == [
            ("Pfam", "10", "90", "SNARE domain"),
            ("PANTHER", "1", "109", "SNARE domain"),
        ]

    def test_short_rows_skipped(self):
        text = "\n".join(["a\tb\tc", _row()])
        rows = iprscan_tsv_to_domains(text)
        assert rows == [("Pfam", "10", "90", "SNARE domain")]

    def test_empty_input_returns_empty_list(self):
        assert iprscan_tsv_to_domains("") == []

    def test_blank_lines_skipped(self):
        text = "\n\n" + _row() + "\n\n"
        assert iprscan_tsv_to_domains(text) == [("Pfam", "10", "90", "SNARE domain")]

    def test_extra_trailing_columns_ignored(self):
        row_with_go = _row() + "\tGO:0006886\tREACT_123"
        rows = iprscan_tsv_to_domains(row_with_go)
        assert rows == [("Pfam", "10", "90", "SNARE domain")]
