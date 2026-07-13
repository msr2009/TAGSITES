"""
test_run_genewise.py — offline unit tests for scripts/run_genewise.py

Tests the two network-free functions:
  select_orientation — picks the higher-scoring result, fires warning when low
  write_rc_fasta     — produces exact reverse-complement of input

Also covers bad-input scenarios via UNC119 fixtures:
  unc119_truncated.genewise.out.txt  — score OK but CDS covers only 2 of 5 exons (~18%)
  unc119_wrong_seq.genewise.out.txt  — score 1.52 bits (garbage alignment, wrong sequence)
"""

import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO

from run_genewise import select_orientation, write_rc_fasta
from crispr_util import reverse_complement


# ── select_orientation ────────────────────────────────────────────────────────

class TestSelectOrientation:
    """src-1 fwd=1177.09 bits (GFF score), rc=-17.62 bits (fixtures in tests/data/)."""

    def test_picks_forward_when_fwd_wins(self, DATA):
        result = select_orientation(
            fwd_out_txt=DATA / "src-1.genewise.out.txt",
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=DATA / "src-1.genewise_rc.out.txt",
            rc_genomic_fa=DATA / "src-1_genomic.fa",  # rc_fa path; content doesn't matter here
            protein_length=536,
        )
        assert result["orientation"] == "+"
        assert result["fwd_score"] == pytest.approx(1177.09, abs=0.1)
        assert result["rc_score"] == pytest.approx(-17.62, abs=0.1)

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
        assert result["winner_score"] == pytest.approx(1177.09, abs=0.1)

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


# ── Bad-input scenarios (UNC119 fixtures) ─────────────────────────────────────

# UNC119 protein is 244 aa; full alignment covers 5 CDS exons = 732 nt = 100% coverage.
UNC119_PROTEIN_LEN = 244


class TestTruncatedGenomicSequence:
    """
    Fixture: unc119_truncated.genewise.out.txt
    Simulates a user providing a genomic region that is too short — only the
    first two of five UNC119 exons fall inside it.  Genewise score is
    reasonable (~150 bits, above the 50-bit floor) but CDS coverage is only
    ~18%, triggering the LOW_COVER_WARN path (threshold is 90%).
    """

    def test_warning_fires_for_low_coverage(self, DATA):
        trunc = DATA / "unc119_truncated.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=trunc,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",   # content unused by select_orientation
            rc_out_txt=trunc,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=UNC119_PROTEIN_LEN,
        )
        assert result["warning"] is not None, "expected a warning for low-coverage alignment"

    def test_warning_mentions_coverage(self, DATA):
        trunc = DATA / "unc119_truncated.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=trunc,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=trunc,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=UNC119_PROTEIN_LEN,
        )
        # Should mention coverage, not score — score is high enough
        assert "coverage" in result["warning"].lower() or "%" in result["warning"]

    def test_score_is_not_the_trigger(self, DATA):
        """Score 150 bits is above LOW_SCORE_WARN=50; warning is about coverage."""
        trunc = DATA / "unc119_truncated.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=trunc,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=trunc,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=UNC119_PROTEIN_LEN,
        )
        assert result["winner_score"] == pytest.approx(150.0, abs=0.1)

    def test_no_warning_without_protein_length(self, DATA):
        """Without protein_length the coverage check is skipped; 150 bits → no warning."""
        trunc = DATA / "unc119_truncated.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=trunc,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=trunc,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
        )
        assert result["warning"] is None


class TestWrongGenomicSequence:
    """
    Fixture: unc119_wrong_seq.genewise.out.txt
    Simulates a user providing a completely wrong genomic sequence.  Genewise
    produces a near-random 1.52-bit alignment — well below LOW_SCORE_WARN=50 —
    triggering the score-based warning regardless of protein_length.
    """

    def test_warning_fires_for_low_score(self, DATA):
        wrong = DATA / "unc119_wrong_seq.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=wrong,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=wrong,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=UNC119_PROTEIN_LEN,
        )
        assert result["warning"] is not None

    def test_warning_mentions_score(self, DATA):
        wrong = DATA / "unc119_wrong_seq.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=wrong,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=wrong,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
            protein_length=UNC119_PROTEIN_LEN,
        )
        assert "WARNING" in result["warning"]
        # Score-based warning always includes the score value
        assert "1.5" in result["warning"] or "bits" in result["warning"]

    def test_warning_fires_without_protein_length(self, DATA):
        """Score check does not require protein_length — bad sequence is caught unconditionally."""
        wrong = DATA / "unc119_wrong_seq.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=wrong,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=wrong,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
        )
        assert result["warning"] is not None

    def test_winner_score_is_very_low(self, DATA):
        wrong = DATA / "unc119_wrong_seq.genewise.out.txt"
        result = select_orientation(
            fwd_out_txt=wrong,
            fwd_genomic_fa=DATA / "src-1_genomic.fa",
            rc_out_txt=wrong,
            rc_genomic_fa=DATA / "src-1_genomic.fa",
        )
        assert result["winner_score"] < 50.0
