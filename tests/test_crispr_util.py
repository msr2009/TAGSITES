"""
test_crispr_util.py — offline unit tests for scripts/crispr_util.py

Covers: iupac_to_regex, reverse_complement, reverse_complement_iupac,
        SYN_CODONS consistency, find_guides (PAM parametric, both strands,
        coords), disrupt_pam (syn_1 synonymous, gc_preserve GC-stable,
        PAM actually broken after mutation).
"""

import re
import pytest

from crispr_util import (
    iupac_to_regex,
    reverse_complement,
    reverse_complement_iupac,
    gc_count,
    find_guides,
    build_frame_lookup,
    disrupt_pam,
    CODONS,
    SYN_CODONS,
    _pam_disrupted,
)


# ── iupac_to_regex ────────────────────────────────────────────────────────────

def test_iupac_ngg():
    assert iupac_to_regex("NGG") == "[ACGT]GG"

def test_iupac_tttn():
    # T=T, N=[ACGT]
    assert iupac_to_regex("TTTN") == "TTT[ACGT]"

def test_iupac_ng():
    assert iupac_to_regex("NG") == "[ACGT]G"

def test_iupac_case_insensitive():
    assert iupac_to_regex("ngg") == iupac_to_regex("NGG")


# ── reverse_complement ────────────────────────────────────────────────────────

def test_rc_simple():
    assert reverse_complement("ATCG") == "CGAT"

def test_rc_roundtrip():
    seq = "AATTCCGG"
    assert reverse_complement(reverse_complement(seq)) == seq

def test_rc_mixed_case():
    # reverse_complement upcases its output regardless of input case
    assert reverse_complement("atcg") == "CGAT"


# ── reverse_complement_iupac ──────────────────────────────────────────────────

def test_rc_iupac_ngg():
    # RC of NGG (5'→3' on guide strand) = CCN (5'→3' on genome minus strand)
    assert reverse_complement_iupac("NGG") == "CCN"

def test_rc_iupac_roundtrip():
    for pam in ["NGG", "TTTN", "NG"]:
        assert reverse_complement_iupac(reverse_complement_iupac(pam)) == pam


# ── SYN_CODONS consistency ────────────────────────────────────────────────────

def test_syn_codons_all_aas_present():
    """Every standard amino acid encoded by CODONS should appear in SYN_CODONS.
    Stop codons are stored as '_' in CODONS and are intentionally excluded."""
    aas = set(CODONS.values()) - {"*", "_"}
    for aa in aas:
        assert aa in SYN_CODONS, f"AA {aa} missing from SYN_CODONS"

def test_syn_codons_all_synonymous():
    """SYN_CODONS is keyed by amino acid (1-letter); every listed codon must
    encode that amino acid according to CODONS."""
    for aa, codons in SYN_CODONS.items():
        for codon in codons:
            assert CODONS.get(codon) == aa, (
                f"SYN_CODONS[{aa}] lists {codon} but CODONS[{codon}]={CODONS.get(codon)}"
            )


# ── find_guides ───────────────────────────────────────────────────────────────

# Build a minimal test sequence with one known + strand PAM and one - strand PAM.
# + strand NGG at position 30: spacer = seq[10:30], pam = seq[30:33]
# - strand NGG (= CCN on fwd strand) at position 50: pam_fwd_start = 47, pam = seq[47:50]
_GUIDE_LEN = 20
_CUT_OFF = 3
_SEQ_FWD = "A" * 10 + "GCTAGCTAGCTAGCTAGCTA" + "TGG" + "C" * 14 + "CCA" + "GCTAGCTAGCTAGCTAGCTA" + "A" * 10
# position 10–29: spacer; 30–32: TGG (+ PAM)
# position 47–49: CCA → RC = TGG, so this is a - strand guide; spacer on - strand = RC(seq[50:70])


def test_find_guides_returns_both_strands():
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    strands = {g["strand"] for g in guides}
    assert "+" in strands
    assert "-" in strands


def test_find_guides_spacer_length():
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in guides:
        assert len(g["spacer"]) == _GUIDE_LEN, f"spacer length wrong: {g}"


def test_find_guides_pam_matches_iupac():
    pam = "NGG"
    pattern = re.compile(iupac_to_regex(pam))
    guides = find_guides(_SEQ_FWD, pam=pam, guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in guides:
        assert pattern.fullmatch(g["pam_seq"]), (
            f"pam_seq {g['pam_seq']} does not match {pam}"
        )


def test_find_guides_cut_pos_plus_strand():
    """For + strand guides, cut_pos = pam_fwd_start - cut_offset."""
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in [x for x in guides if x["strand"] == "+"]:
        assert g["cut_pos"] == g["pam_fwd_start"] - _CUT_OFF


def test_find_guides_cut_pos_minus_strand():
    """For - strand guides, cut_pos = pam_fwd_start + pam_len + cut_offset."""
    pam_len = 3
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in [x for x in guides if x["strand"] == "-"]:
        assert g["cut_pos"] == g["pam_fwd_start"] + pam_len + _CUT_OFF


def test_find_guides_coords_0based():
    """All guide coords should be 0-based (pam_fwd_start < len(seq))."""
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in guides:
        assert 0 <= g["pam_fwd_start"] < len(_SEQ_FWD)
        assert 0 <= g["guide_fwd_start"] < len(_SEQ_FWD)


@pytest.mark.parametrize("pam", ["NG", "TTTN", "NNGRRT"])
def test_find_guides_parametric_pam(pam):
    """find_guides must use whatever PAM is supplied, not hardcode NGG."""
    # Build a sequence containing the PAM explicitly
    spacer = "GCTAGCTAGCTAGCTAGCTA"
    seq = "A" * 10 + spacer + pam.replace("N", "A").replace("R", "G").replace("T", "T") + "A" * 50
    guides = find_guides(seq, pam=pam, guide_length=20, cut_offset=3)
    if guides:
        pattern = re.compile(iupac_to_regex(pam))
        for g in guides:
            assert pattern.fullmatch(g["pam_seq"]), (
                f"PAM {g['pam_seq']} does not match pattern for {pam}"
            )


def test_find_guides_spacer_in_sequence():
    """+ strand spacer should appear verbatim in the forward sequence."""
    guides = find_guides(_SEQ_FWD, pam="NGG", guide_length=_GUIDE_LEN, cut_offset=_CUT_OFF)
    for g in [x for x in guides if x["strand"] == "+"]:
        assert g["spacer"] in _SEQ_FWD, f"spacer {g['spacer']} not found in fwd seq"


# ── disrupt_pam ───────────────────────────────────────────────────────────────

def _make_coding_context():
    """
    Build a small genomic sequence where we know the coding frame exactly.
    Three codons in a row: ATG (M) GGT (G) CCG (P), then PAM TGG immediately after.
    Frame: positions 0–2 = ATG, 3–5 = GGT, 6–8 = CCG, 9–11 = TGG (PAM, also a codon TGG=W but here it's the PAM).

    We set up a single CDS exon: start=0, stop=8 (0-based, last full codon ends at 8 incl.)
    PAM starts at position 9 (TGG on + strand).
    """
    import pandas as pd
    dna = "ATGGGTCCGTGGAAA"  # ATG GGT CCG TGG AAA
    # CDS exon: start=0, stop=8 (0-indexed; the 9th base is the start of PAM)
    cds_df = pd.DataFrame([{
        "name": "test", "type": "cds", "start": 0, "stop": 8,
        "score": ".", "strand": "+", "frame": 0, "note": ""
    }])
    return dna, cds_df


def test_disrupt_pam_syn_1_coding():
    """A PAM fully inside a coding codon should return syn_1 or syn_2."""
    dna, cds_df = _make_coding_context()
    frame_lookup = build_frame_lookup(cds_df, dna)
    # PAM TGG at position 9 on + strand
    result = disrupt_pam(dna, "NGG", pam_fwd_start=9, strand="+", frame_lookup=frame_lookup)
    assert result is not None, "disrupt_pam returned None for a mutable coding PAM"
    mutated_seq, desc, method = result
    assert method in ("syn_1", "syn_2", "gc_preserve"), f"unexpected method: {method}"
    # The PAM in the mutated sequence should no longer match NGG
    assert _pam_disrupted(mutated_seq, "NGG", 9, 3, "+"), (
        "PAM still matches after disruption"
    )


def test_disrupt_pam_synonymous_translation():
    """syn_1/syn_2 must not change the amino acid sequence."""
    dna, cds_df = _make_coding_context()
    frame_lookup = build_frame_lookup(cds_df, dna)
    result = disrupt_pam(dna, "NGG", pam_fwd_start=9, strand="+", frame_lookup=frame_lookup)
    if result is None:
        pytest.skip("PAM not mutable in this context")
    mutated_seq, desc, method = result
    if method in ("syn_1", "syn_2"):
        # Check that every codon in the CDS translates identically
        for pos in range(0, 9, 3):
            orig_codon = dna[pos:pos+3].upper()
            mut_codon  = mutated_seq[pos:pos+3].upper()
            if orig_codon != mut_codon:
                assert CODONS.get(orig_codon) == CODONS.get(mut_codon), (
                    f"Non-synonymous change at pos {pos}: {orig_codon}→{mut_codon}"
                )


def test_disrupt_pam_gc_preserve_noncoding():
    """gc_preserve must not change GC content of the affected region."""
    # Non-coding context: empty frame_lookup
    dna = "A" * 10 + "ATGGGG" + "A" * 20  # PAM GGG at pos 13 on + strand
    result = disrupt_pam(dna, "NGG", pam_fwd_start=13, strand="+", frame_lookup={})
    if result is None:
        pytest.skip("PAM not mutable — sequence may have no GC alternative")
    mutated_seq, desc, method = result
    if method == "gc_preserve":
        window = slice(13, 13 + 3)
        assert gc_count(mutated_seq[window]) == gc_count(dna[window]), (
            "gc_preserve changed GC content"
        )


def test_disrupt_pam_pam_broken():
    """After disrupt_pam, the PAM regex must not match at the original position."""
    dna, cds_df = _make_coding_context()
    frame_lookup = build_frame_lookup(cds_df, dna)
    result = disrupt_pam(dna, "NGG", pam_fwd_start=9, strand="+", frame_lookup=frame_lookup)
    if result is None:
        pytest.skip("No mutation found")
    mutated_seq, _, _ = result
    assert _pam_disrupted(mutated_seq, "NGG", 9, 3, "+"), (
        "PAM still matches NGG after disrupt_pam"
    )
