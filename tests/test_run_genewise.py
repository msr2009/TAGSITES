"""
test_run_genewise.py — offline unit tests for scripts/run_genewise.py

Tests the two network-free functions:
  select_orientation — picks the higher-scoring result, fires warning when low
  write_rc_fasta     — produces exact reverse-complement of input
"""

import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO

from run_genewise import select_orientation, write_rc_fasta
from crispr_util import reverse_complement


# ── select_orientation ────────────────────────────────────────────────────────

class TestSelectOrientation:
    """src-1 fwd=1138 bits, rc=1.39 bits (fixtures in tests/data/)."""

    def test_picks_forward_when_fwd_wins(self, DATA):
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise_rc.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",  # rc_fa path; content doesn't matter here
            protein_length=536,
        )
        assert result["orientation"] == "+"
        assert result["fwd_score"] == pytest.approx(1138.03, abs=0.1)
        assert result["rc_score"] == pytest.approx(1.39, abs=0.1)

    def test_no_warning_for_high_score_winner(self, DATA):
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise_rc.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=536,
        )
        assert result["warning"] is None

    def test_picks_rc_when_rc_wins(self, DATA):
        """Swap fwd/rc: the high-score result is now in the rc slot."""
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise_rc.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=536,
        )
        assert result["orientation"] == "rc"
        assert result["winner_score"] == pytest.approx(1138.03, abs=0.1)

    def test_warning_fires_for_low_score_winner(self, DATA):
        """If both inputs are the garbage RC result, winner is 1.39 → warning."""
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise_rc.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise_rc.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=536,
        )
        assert result["warning"] is not None
        assert "WARNING" in result["warning"]

    def test_result_keys_present(self, DATA):
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise_rc.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",
        )
        for key in ("orientation", "out_txt", "genomic_fa",
                    "fwd_score", "rc_score", "winner_score", "warning"):
            assert key in result


# ── write_rc_fasta ────────────────────────────────────────────────────────────

class TestWriteRcFasta:

    def test_rc_sequence_is_correct(self, DATA, tmp_path):
        out = tmp_path / "rc.fa"
        write_rc_fasta(DATA / "snb-1_genomic.fa", out)
        original = str(list(SeqIO.parse(DATA / "snb-1_genomic.fa", "fasta"))[0].seq)
        rc_seq = str(list(SeqIO.parse(out, "fasta"))[0].seq)
        assert rc_seq == reverse_complement(original)

    def test_rc_header_contains_rc(self, DATA, tmp_path):
        out = tmp_path / "rc.fa"
        write_rc_fasta(DATA / "snb-1_genomic.fa", out)
        rec = list(SeqIO.parse(out, "fasta"))[0]
        assert "rc" in rec.id.lower() or "rc" in rec.description.lower()

    def test_rc_length_preserved(self, DATA, tmp_path):
        out = tmp_path / "rc.fa"
        write_rc_fasta(DATA / "snb-1_genomic.fa", out)
        original = str(list(SeqIO.parse(DATA / "snb-1_genomic.fa", "fasta"))[0].seq)
        rc_seq = str(list(SeqIO.parse(out, "fasta"))[0].seq)
        assert len(rc_seq) == len(original)
