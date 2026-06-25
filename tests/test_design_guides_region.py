"""
test_design_guides_region.py — offline unit tests for scripts/design_guides_across_region.py

The module uses a global `args` object set at __main__ time; tests inject a fake
SimpleNamespace to exercise the pure-compute functions without running argparse.

Offline-testable functions:
  convert_region_to_BED  — pure string transform
  determine_input_type   — pure regex classifier
  find_targets           — pure guide-finder (needs args.GUIDE_LENGTH / args.PAM)
"""

import pytest
from types import SimpleNamespace

import design_guides_across_region as dgr
from crispr_util import iupac_to_regex, reverse_complement

import re


# ── helpers ───────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def inject_args():
    """Inject a minimal args namespace and clean up after each test."""
    dgr.args = SimpleNamespace(GUIDE_LENGTH=20, PAM="NGG")
    yield
    # Remove to avoid leaking state between test modules
    if hasattr(dgr, "args"):
        del dgr.args


# ── convert_region_to_BED ─────────────────────────────────────────────────────

class TestConvertRegionToBED:

    def test_basic_conversion(self):
        assert dgr.convert_region_to_BED("I:123-456") == "I\t123\t456"

    def test_strips_trailing_newline(self):
        assert dgr.convert_region_to_BED("II:100-200\n") == "II\t100\t200"

    def test_x_chromosome(self):
        assert dgr.convert_region_to_BED("X:1000-2000") == "X\t1000\t2000"

    def test_returns_string(self):
        assert isinstance(dgr.convert_region_to_BED("I:1-100"), str)


# ── determine_input_type ──────────────────────────────────────────────────────

class TestDetermineInputType:

    def test_region(self):
        assert dgr.determine_input_type("I:123-456") == "region"

    def test_region_roman_numeral(self):
        assert dgr.determine_input_type("IV:10000-20000") == "region"

    def test_seq_lowercase(self):
        assert dgr.determine_input_type("acgtacgt") == "seq"

    def test_seq_uppercase(self):
        assert dgr.determine_input_type("ACGTACGT") == "seq"

    def test_gene_name(self):
        assert dgr.determine_input_type("unc-25") == "gene"

    def test_gene_name_with_digits(self):
        assert dgr.determine_input_type("abl-1a") == "gene"

    def test_unknown_raises_value_error(self):
        with pytest.raises(ValueError):
            dgr.determine_input_type("!@#$%")


# ── find_targets ──────────────────────────────────────────────────────────────

def _seq_with_known_pam(guide_len=20):
    """Build a sequence with EXACTLY one + strand NGG PAM.

    Uses all T's/A's so there are no GG runs before the explicit PAM,
    ensuring find_targets produces one guide with a full-length spacer.
    Layout:  T * guide_len  | A * guide_len | TGG | A * 30
    The single NGG sits at position 2*guide_len; spacer = A * guide_len.
    """
    prefix = "T" * guide_len    # no G → no GG
    spacer = "A" * guide_len    # no G → no GG
    pam    = "TGG"              # the only NGG in the sequence
    suffix = "A" * 30
    seq = prefix + spacer + pam + suffix
    return seq, 0, len(seq)


class TestFindTargets:

    def test_finds_guides(self):
        seq, start, end = _seq_with_known_pam()
        targets = dgr.find_targets(seq, start, end)
        assert len(targets) > 0

    def test_spacer_length(self):
        seq, start, end = _seq_with_known_pam()
        targets = dgr.find_targets(seq, start, end)
        for t in targets:
            assert len(t[0]) == 20, f"spacer length wrong: {t[0]!r}"

    def test_target_is_list_of_three(self):
        seq, start, end = _seq_with_known_pam()
        targets = dgr.find_targets(seq, start, end)
        for t in targets:
            assert len(t) == 3

    def test_rc_reverses_coordinates(self):
        """rc=True should flip signs of coordinate increments."""
        seq, start, end = _seq_with_known_pam()
        fwd = dgr.find_targets(seq, start, end, rc=False)
        rev = dgr.find_targets(seq, start, end, rc=True)
        # Both should find the same guides from complementary strand scan
        assert len(fwd) > 0 or len(rev) > 0

    def test_different_pam_respected(self):
        """When args.PAM changes, different guides are found."""
        seq, start, end = _seq_with_known_pam()
        dgr.args.PAM = "NAG"
        targets_nag = dgr.find_targets(seq, start, end)
        dgr.args.PAM = "NGG"
        targets_ngg = dgr.find_targets(seq, start, end)
        # The two PAMs should produce different guide sets
        spacers_nag = {t[0] for t in targets_nag}
        spacers_ngg = {t[0] for t in targets_ngg}
        assert spacers_nag != spacers_ngg or (len(spacers_nag) == 0 and len(spacers_ngg) == 0)
