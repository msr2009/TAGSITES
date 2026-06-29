"""
test_extract_from_pdb.py — offline tests for scripts/extract_from_pdb.py

Uses the committed UNC-25 AlphaFold PDB fixture (433 residues).
Tests that extract_from_pdb.main() writes:
  - a pLDDT file  (residue_index  b-factor)
  - a SASA file   (residue_index  relative_sasa  one-letter-aa)
with expected row counts and value ranges.
"""

import pytest
from pathlib import Path

from extract_from_pdb import main as pdb_main

EXPECTED_ROWS = 433  # UNC-25 AF model has 433 resolved residues


@pytest.fixture(scope="module")
def pdb_outputs(DATA, tmp_path_factory):
    """Run extract_from_pdb.main() once per test-session and return output paths."""
    tmp = tmp_path_factory.mktemp("pdb_out")
    plddt_out = str(tmp / "unc25_plddt.txt")
    pdb_main(str(DATA / "unc-25.pdb"), plddt_out)
    sasa_out = plddt_out.replace(".txt", ".sasa.txt")
    return Path(plddt_out), Path(sasa_out)


# ── pLDDT output ──────────────────────────────────────────────────────────────

class TestPlddt:

    def test_plddt_file_created(self, pdb_outputs):
        plddt, _ = pdb_outputs
        assert plddt.exists(), "pLDDT output file was not written"

    def test_plddt_row_count(self, pdb_outputs):
        plddt, _ = pdb_outputs
        rows = [l for l in plddt.read_text().splitlines() if l.strip()]
        assert len(rows) == EXPECTED_ROWS, f"expected {EXPECTED_ROWS} rows, got {len(rows)}"

    def test_plddt_first_index_is_1(self, pdb_outputs):
        plddt, _ = pdb_outputs
        first_line = plddt.read_text().splitlines()[0]
        idx = int(first_line.split("\t")[0])
        assert idx == 1

    def test_plddt_values_in_range(self, pdb_outputs):
        plddt, _ = pdb_outputs
        for line in plddt.read_text().splitlines():
            if not line.strip():
                continue
            _, val = line.split("\t")
            v = float(val)
            assert 0.0 <= v <= 100.0, f"pLDDT out of range: {v}"

    def test_plddt_indices_are_sequential(self, pdb_outputs):
        plddt, _ = pdb_outputs
        indices = [int(l.split("\t")[0]) for l in plddt.read_text().splitlines() if l.strip()]
        assert indices == list(range(1, EXPECTED_ROWS + 1))


# ── SASA output ───────────────────────────────────────────────────────────────

class TestSasa:

    def test_sasa_file_created(self, pdb_outputs):
        _, sasa = pdb_outputs
        assert sasa.exists(), "SASA output file was not written"

    def test_sasa_row_count(self, pdb_outputs):
        _, sasa = pdb_outputs
        rows = [l for l in sasa.read_text().splitlines() if l.strip()]
        assert len(rows) == EXPECTED_ROWS

    def test_sasa_three_columns(self, pdb_outputs):
        _, sasa = pdb_outputs
        for line in sasa.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            assert len(parts) == 3, f"SASA row has {len(parts)} cols: {line!r}"

    def test_sasa_relative_values_numeric(self, pdb_outputs):
        _, sasa = pdb_outputs
        for line in sasa.read_text().splitlines():
            if not line.strip():
                continue
            _, rel_sasa, _ = line.split("\t")
            float(rel_sasa)  # must be parseable

    def test_sasa_aa_column_is_one_letter(self, pdb_outputs):
        _, sasa = pdb_outputs
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX")
        for line in sasa.read_text().splitlines():
            if not line.strip():
                continue
            _, _, aa = line.split("\t")
            assert aa in valid_aa, f"unexpected AA code: {aa!r}"
