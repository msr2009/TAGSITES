"""
test_parse_genewise.py — offline unit tests for scripts/parse_genewise.py

Uses committed fixtures in tests/data/:
  snb-1  — 109 residues, 2 exons, no split codons
  src-1  — 533 residues, 6 exons, 5 split codons; forward score 1138 bits
  src-1 rc — same gene, RC genomic input; garbage score 1.39 bits
"""

import pytest
from Bio import SeqIO

from parse_genewise import (
    parse_genewise,
    enumerate_insertion_sites,
    parse_genewise_score,
    cds_coverage,
)


# ── parse_genewise ────────────────────────────────────────────────────────────

class TestParseGenewise:

    def test_snb1_exon_count(self, DATA):
        cds = parse_genewise(DATA / "snb-1.genewise.out.txt")
        assert len(cds) == 2

    def test_snb1_cds_only(self, DATA):
        cds = parse_genewise(DATA / "snb-1.genewise.out.txt")
        assert (cds["type"] == "cds").all()

    def test_snb1_0indexed(self, DATA):
        cds = parse_genewise(DATA / "snb-1.genewise.out.txt")
        # The first exon starts at position 0 in the snb-1 coding-only genomic FASTA
        assert cds.iloc[0]["start"] == 0

    def test_src1_exon_count(self, DATA):
        cds = parse_genewise(DATA / "src-1.genewise.out.txt")
        assert len(cds) == 6

    def test_src1_sorted_by_start(self, DATA):
        cds = parse_genewise(DATA / "src-1.genewise.out.txt")
        starts = list(cds["start"])
        assert starts == sorted(starts)


# ── enumerate_insertion_sites ─────────────────────────────────────────────────

class TestEnumerateInsertionSites:

    def _snb1(self, DATA):
        cds = parse_genewise(DATA / "snb-1.genewise.out.txt")
        dna = str(list(SeqIO.parse(DATA / "snb-1_genomic.fa", "fasta"))[0].seq)
        return enumerate_insertion_sites(cds, dna)

    def _src1(self, DATA):
        cds = parse_genewise(DATA / "src-1.genewise.out.txt")
        dna = str(list(SeqIO.parse(DATA / "src-1_genomic.fa", "fasta"))[0].seq)
        return enumerate_insertion_sites(cds, dna)

    def test_snb1_residue_count(self, DATA):
        sites = self._snb1(DATA)
        assert len(sites) == 109

    def test_snb1_first_residue(self, DATA):
        sites = self._snb1(DATA)
        assert sites.iloc[0]["amino_acid"] == "M"
        assert sites.iloc[0]["insert_pos"] == 3

    def test_snb1_no_split_codons(self, DATA):
        sites = self._snb1(DATA)
        assert sites["is_split_codon"].sum() == 0

    def test_snb1_last_exon_sentinel(self, DATA):
        sites = self._snb1(DATA)
        # snb-1 only has 2 exons; last residue is on the last exon → sentinel
        assert sites.iloc[-1]["dist_to_3p_splice"] == 99999

    def test_src1_residue_count(self, DATA):
        sites = self._src1(DATA)
        assert len(sites) == 533

    def test_src1_first_residue(self, DATA):
        sites = self._src1(DATA)
        assert sites.iloc[0]["amino_acid"] == "M"
        assert sites.iloc[0]["insert_pos"] == 510

    def test_src1_split_codons(self, DATA):
        sites = self._src1(DATA)
        assert sites["is_split_codon"].sum() == 5

    def test_src1_last_exon_sentinel(self, DATA):
        sites = self._src1(DATA)
        assert sites.iloc[-1]["dist_to_3p_splice"] == 99999

    def test_src1_residue_indices_sequential(self, DATA):
        sites = self._src1(DATA)
        assert list(sites["residue_index"]) == list(range(1, 534))

    def test_snb1_translation_matches_fasta(self, DATA):
        """Translated residues from enumerated sites must equal the protein FASTA.
        Uses snb-1 (109 aa, no split codons, FASTA exactly matches Genewise input)."""
        from crispr_util import CODONS

        sites = self._snb1(DATA)
        protein_seq = str(list(SeqIO.parse(DATA / "snb-1.fa", "fasta"))[0].seq).rstrip("*")
        dna = str(list(SeqIO.parse(DATA / "snb-1_genomic.fa", "fasta"))[0].seq)

        mismatches = 0
        for _, row in sites.iterrows():
            pos = int(row["insert_pos"])
            codon = dna[pos - 3: pos].upper()
            aa_from_codon = CODONS.get(codon, "?")
            aa_from_fasta = protein_seq[int(row["residue_index"]) - 1]
            if aa_from_codon != aa_from_fasta:
                mismatches += 1

        assert mismatches == 0, f"{mismatches} residues do not translate to the expected AA"


# ── parse_genewise_score ──────────────────────────────────────────────────────

class TestParseGenewiseScore:

    def test_src1_fwd_score(self, DATA):
        score = parse_genewise_score(DATA / "src-1.genewise.out.txt")
        assert score == pytest.approx(1138.03, abs=0.1)

    def test_src1_rc_score(self, DATA):
        score = parse_genewise_score(DATA / "src-1.genewise_rc.out.txt")
        assert score == pytest.approx(1.39, abs=0.1)

    def test_returns_float(self, DATA):
        score = parse_genewise_score(DATA / "src-1.genewise.out.txt")
        assert isinstance(score, float)


# ── cds_coverage ─────────────────────────────────────────────────────────────

class TestCdsCoverage:

    def test_src1_fwd_coverage_near_1(self, DATA):
        cds = parse_genewise(DATA / "src-1.genewise.out.txt")
        cover = cds_coverage(cds, protein_length=536)
        assert cover == pytest.approx(1.0, abs=0.02)

    def test_src1_rc_coverage_low(self, DATA):
        cds = parse_genewise(DATA / "src-1.genewise_rc.out.txt")
        cover = cds_coverage(cds, protein_length=536)
        assert cover < 0.10
