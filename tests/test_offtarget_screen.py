"""
Unit tests for scripts/offtarget_screen.py — the network-free core of the
off-target / primer-specificity screen.

Every HSP here is hand-built, so these tests pin the alignment arithmetic
(mismatch counting, the 3'-end rule, gap handling, PAM extraction, block
mapping) without touching EBI.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import offtarget_screen as ots  # noqa: E402


CFG = ots.load_config()


def mk(qseq, hseq, q_from=1, h_from=1001, h_to=None, acc="X1", desc="", identity=95.0,
       bits=200.0, klass="duplicate"):
    """Build one HSP dict the way parse_blast_json would."""
    ungapped = len(qseq.replace("-", ""))
    return {
        "acc": acc, "desc": desc, "os": "Caenorhabditis elegans",
        "q_from": q_from, "q_to": q_from + ungapped - 1,
        "h_from": h_from,
        "h_to": h_to if h_to is not None else h_from + len(hseq.replace("-", "")) - 1,
        "h_strand": "+" if (h_to is None or h_to >= h_from) else "-",
        "qseq": qseq, "hseq": hseq,
        "identity": identity, "bits": bits, "evalue": 1e-20,
        "align_len": len(qseq), "klass": klass,
    }


# ── score_window: mismatches, 3'-end rule, gaps, coverage ────────────────────

def test_perfect_window_is_flagged_and_perfect():
    q = "ACGTACGTACGTACGTACGT"
    h = ots.score_window(mk(q, q), 0, 20, "right", CFG)
    assert h["mismatches"] == 0 and h["gaps"] == 0
    assert h["perfect"] and h["three_prime_ok"] and h["flagged"]


def test_three_prime_terminal_mismatch_is_not_flagged():
    q = "ACGTACGTACGTACGTACGT"
    h = mk(q, q[:-1] + "A")          # last base differs
    r = ots.score_window(h, 0, 20, "right", CFG)
    assert r["mismatches"] == 1
    assert r["three_prime_ok"] is False
    assert r["flagged"] is False, "a 3'-terminal mismatch must not be reported"


def test_five_prime_mismatch_is_flagged():
    q = "ACGTACGTACGTACGTACGT"
    h = mk(q, "T" + q[1:])           # first base differs
    r = ots.score_window(h, 0, 20, "right", CFG)
    assert r["mismatches"] == 1 and r["three_prime_ok"] and r["flagged"]
    assert not r["perfect"]


def test_three_prime_side_left_flips_which_end_matters():
    q = "ACGTACGTACGTACGTACGT"
    first_differs = mk(q, "T" + q[1:])
    # With the 3' end on the left, a mismatch at the first base is now terminal
    assert ots.score_window(first_differs, 0, 20, "left", CFG)["flagged"] is False
    last_differs = mk(q, q[:-1] + "A")
    assert ots.score_window(last_differs, 0, 20, "left", CFG)["flagged"] is True


def test_too_many_mismatches_not_flagged():
    q = "ACGTACGTACGTACGTACGT"
    h = mk(q, "TGCA" + q[4:])        # complement of ACGT -> 4 mismatches, all 5'
    r = ots.score_window(h, 0, 20, "right", CFG)
    assert r["mismatches"] == 4 and r["three_prime_ok"]
    assert r["flagged"] is False, "max_mismatch default is 3"


def test_gap_in_subject_counts_and_blocks_perfect():
    q = "ACGTACGTACGTACGTACGT"
    h = mk(q, q[:5] + "-" + q[6:])
    r = ots.score_window(h, 0, 20, "right", CFG)
    assert r["gaps"] == 1 and not r["perfect"] and r["flagged"]


def test_partial_coverage_returns_none():
    q = "ACGTACGTAC"
    # HSP covers region 0-9 only; asking for 0-19 cannot be scored
    assert ots.score_window(mk(q, q), 0, 20, "right", CFG) is None


def test_window_offset_by_q_from():
    q = "ACGTACGTACGTACGTACGT"
    h = mk(q, q, q_from=101)         # HSP starts at region coord 100
    assert ots.score_window(h, 100, 120, "right", CFG)["perfect"]
    assert ots.score_window(h, 0, 20, "right", CFG) is None


# ── guide footprints and PAM extraction ──────────────────────────────────────

def test_guide_footprint_plus_strand():
    spacer, pam, side = ots.guide_footprint({"pam_fwd_start": 100, "guide_strand": "+"})
    assert spacer == (80, 100) and pam == (100, 103) and side == "right"


def test_guide_footprint_minus_strand():
    spacer, pam, side = ots.guide_footprint({"pam_fwd_start": 100, "guide_strand": "-"})
    assert spacer == (103, 123) and pam == (100, 103) and side == "left"


def test_offtarget_pam_plus_strand_read_directly():
    q = "A" * 20 + "CGG"
    h = mk(q, q)
    assert ots.offtarget_pam(h, (20, 23), "+") == "CGG"
    assert ots.pam_matches("CGG")


def test_offtarget_pam_minus_strand_is_reverse_complemented():
    # forward bases CCG read on the guide strand are CGG
    q = "CCG" + "A" * 20
    h = mk(q, q)
    assert ots.offtarget_pam(h, (0, 3), "-") == "CGG"


def test_offtarget_pam_absent_is_detected():
    q = "A" * 20 + "CTA"
    assert ots.pam_matches(ots.offtarget_pam(mk(q, q), (20, 23), "+")) is False


def test_offtarget_pam_uncovered_returns_none():
    q = "A" * 20
    assert ots.offtarget_pam(mk(q, q), (20, 23), "+") is None


def test_offtarget_pam_gap_returns_none():
    q = "A" * 20 + "CGG"
    h = mk(q, "A" * 20 + "C-G")
    assert ots.offtarget_pam(h, (20, 23), "+") is None


def test_screen_guide_requires_pam():
    spacer = "ACGTACGTACGTACGTACGT"
    row = {"pam_fwd_start": 20, "guide_strand": "+"}
    with_pam = ots.screen_guide([mk(spacer + "CGG", spacer + "CGG")], row, cfg=CFG)
    assert with_pam["n_total"] == 1 and with_pam["n_identical"] == 1
    no_pam = ots.screen_guide([mk(spacer + "CTA", spacer + "CTA")], row, cfg=CFG)
    assert no_pam["n_total"] == 0, "a near-match without a PAM cannot cut"


def test_screen_guide_reports_pam_unverified_when_not_covered():
    spacer = "ACGTACGTACGTACGTACGT"
    row = {"pam_fwd_start": 20, "guide_strand": "+"}
    res = ots.screen_guide([mk(spacer, spacer)], row, cfg=CFG)
    assert res["n_total"] == 1 and res["n_pam_unverified"] == 1
    assert res["n_identical"] == 0, "unverified PAM must not count as identical"


# ── classification ───────────────────────────────────────────────────────────

def test_self_hit_classified_by_identity_and_coverage():
    q = "ACGT" * 50
    h = mk(q, q, identity=100.0)
    ots.classify_hsps([h], query_len=200, exons=[], cfg=CFG)
    assert h["klass"] == "self"


def test_partial_high_identity_hit_is_a_duplicate_not_self():
    q = "ACGT" * 25                      # 100 bp of a 1000 bp region
    h = mk(q, q, identity=100.0)
    ots.classify_hsps([h], query_len=1000, exons=[], cfg=CFG)
    assert h["klass"] == "duplicate"


def test_transcript_detected_by_keyword():
    q = "ACGT" * 25
    h = mk(q, q, identity=99.9, desc="Caenorhabditis elegans snt-1 mRNA, complete cds")
    ots.classify_hsps([h], query_len=1000, exons=[], cfg=CFG)
    assert h["klass"] == "transcript"


def _no_keywords():
    """Config with the keyword shortcut off, to exercise the structural test alone."""
    cfg = dict(CFG)
    cfg["transcripts"] = dict(CFG["transcripts"], use_keywords=False)
    return cfg


def test_transcript_detected_structurally_without_keywords():
    """An mRNA is split into one HSP per exon and skips the intron between them."""
    exons = [(0, 100), (200, 300)]
    e1 = mk("ACGT" * 25, "ACGT" * 25, q_from=1, identity=100.0, acc="TX", desc="")
    e2 = mk("TGCA" * 25, "TGCA" * 25, q_from=201, identity=100.0, acc="TX", desc="")
    ots.classify_hsps([e1, e2], query_len=1000, exons=exons, cfg=_no_keywords())
    assert e1["klass"] == "transcript" and e2["klass"] == "transcript"


def test_single_exon_confined_hit_is_a_duplicate_not_a_transcript():
    """The gene-family case: a paralog's coding exon is also exon-confined.

    Classifying it as a transcript would silently discard exactly the hits this
    screen exists to find, so one exon must never be enough on its own.
    """
    q = "ACGT" * 25
    h = mk(q, q, q_from=1, identity=100.0, acc="DUP", desc="")
    ots.classify_hsps([h], query_len=1000, exons=[(0, 100), (200, 300)],
                      cfg=_no_keywords())
    assert h["klass"] == "duplicate"


def test_exon_spanning_hit_that_covers_the_intron_is_genomic():
    """Covering the intron proves the hit is genomic, not a transcript."""
    q = "ACGT" * 75                      # 300 bp, spanning both exons and the intron
    h = mk(q, q, q_from=1, identity=100.0, acc="GEN", desc="")
    ots.classify_hsps([h], query_len=1000, exons=[(0, 100), (200, 300)],
                      cfg=_no_keywords())
    assert h["klass"] == "duplicate"


def test_diverged_paralog_is_never_a_transcript():
    """A paralog at 85% identity must survive even when it looks exon-shaped."""
    e1 = mk("ACGT" * 25, "ACGT" * 25, q_from=1, identity=85.0, acc="PAR", desc="")
    e2 = mk("TGCA" * 25, "TGCA" * 25, q_from=201, identity=85.0, acc="PAR", desc="")
    ots.classify_hsps([e1, e2], query_len=1000, exons=[(0, 100), (200, 300)],
                      cfg=_no_keywords())
    assert e1["klass"] == "duplicate"


def test_genomic_hit_spanning_introns_is_not_a_transcript():
    q = "ACGT" * 25
    h = mk(q, q, identity=95.0, desc="")
    # Only a fifth of the footprint is exonic — intronic sequence aligned too
    ots.classify_hsps([h], query_len=1000, exons=[(0, 20)], cfg=_no_keywords())
    assert h["klass"] == "duplicate"


def test_classification_is_per_accession():
    """Two HSPs of one record share a verdict; a different accession is judged alone."""
    a1 = mk("ACGT" * 25, "ACGT" * 25, q_from=1, identity=100.0, acc="TX", desc="")
    a2 = mk("TGCA" * 25, "TGCA" * 25, q_from=201, identity=100.0, acc="TX", desc="")
    b1 = mk("ACGT" * 25, "ACGT" * 25, q_from=1, identity=100.0, acc="DUP", desc="")
    ots.classify_hsps([a1, a2, b1], query_len=1000, exons=[(0, 100), (200, 300)],
                      cfg=_no_keywords())
    assert a1["klass"] == a2["klass"] == "transcript"
    assert b1["klass"] == "duplicate"


# ── post filter ──────────────────────────────────────────────────────────────

def test_post_filter_drops_weak_duplicates_only():
    strong = mk("ACGT" * 25, "ACGT" * 25, identity=95.0, bits=300.0, acc="strong")
    weak = mk("ACGT" * 5, "ACGT" * 5, identity=70.0, bits=10.0, acc="weak")
    selfish = mk("ACGT" * 5, "ACGT" * 5, identity=70.0, bits=10.0, acc="self",
                 klass="self")
    kept, dropped = ots.apply_post_filter([strong, weak, selfish], CFG)
    assert dropped == 1
    assert {h["acc"] for h in kept} == {"strong", "self"}


def test_permissive_evalue_keeps_short_perfect_segment_until_post_filter():
    """A short perfect duplicate is exactly what a strict E-value would discard."""
    short_perfect = mk("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                       "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                       identity=100.0, bits=95.0)
    ots.classify_hsps([short_perfect], query_len=5000, exons=[], cfg=CFG)
    assert short_perfect["klass"] == "duplicate"
    kept, dropped = ots.apply_post_filter([short_perfect], CFG)
    assert dropped == 0 and len(kept) == 1


# ── concatenated spacer query (Screen B) ─────────────────────────────────────

def test_build_spacer_query_layout():
    seq, block = ots.build_spacer_query(["A" * 23, "C" * 23], separator_len=25)
    assert block == 23
    assert len(seq) == 23 + 25 + 23
    assert seq[23:48] == "N" * 25


def test_block_of_maps_coordinates_and_rejects_separators():
    block_len, sep, n = 23, 25, 3
    assert ots.block_of(1, block_len, sep, n) == (0, 0)
    assert ots.block_of(23, block_len, sep, n) == (0, 22)
    assert ots.block_of(24, block_len, sep, n) is None, "inside a separator"
    assert ots.block_of(49, block_len, sep, n) == (1, 0)
    assert ots.block_of(10_000, block_len, sep, n) is None, "past the last block"


def test_screen_spacer_hits_maps_back_to_the_right_spacer():
    s0, s1 = "ACGTACGTACGTACGTACGT", "TTTTGGGGCCCCAAAATTTT"
    spacers = [s0 + "CGG", s1 + "CGG"]
    seq, block = ots.build_spacer_query(spacers, 25)
    # An HSP over the second block only
    start = block + 25
    h = mk(spacers[1], spacers[1], q_from=start + 1)
    out, n_self = ots.screen_spacer_hits([h], spacers, block, 25, cfg=CFG)
    assert list(out) == [1] and out[1][0]["perfect"]
    assert out[1][0]["screen"] == "spacer" and n_self == 0


def test_screen_spacer_hits_ignores_hsp_straddling_a_separator():
    spacers = ["A" * 20 + "CGG", "C" * 20 + "CGG"]
    seq, block = ots.build_spacer_query(spacers, 25)
    straddle = mk("A" * 10, "A" * 10, q_from=block - 4)   # runs into the N-run
    hits, _ = ots.screen_spacer_hits([straddle], spacers, block, 25, cfg=CFG)
    assert hits == {}


def test_spacer_screen_excludes_the_on_target_site():
    """Without this every guide reports a perfect hit per copy of the source genome.

    ENA carries six independent C. elegans genome submissions, so the on-target
    alone produced a median of 7 spurious hits per guide before this exclusion.
    """
    spacers = ["ACGTACGTACGTACGTACGT" + "CGG"]
    _, block = ots.build_spacer_query(spacers, 25)
    on_target = mk(spacers[0], spacers[0], q_from=1, acc="CHR2", h_from=6828987)
    spans = {"CHR2": [(6827016, 6836486)]}
    hits, n_self = ots.screen_spacer_hits([on_target], spacers, block, 25,
                                          cfg=CFG, exclude_spans=spans)
    assert hits == {} and n_self == 1
    # the same hit elsewhere on the same chromosome is a genuine off-target
    elsewhere = mk(spacers[0], spacers[0], q_from=1, acc="CHR2", h_from=1_000_000)
    hits, n_self = ots.screen_spacer_hits([elsewhere], spacers, block, 25,
                                          cfg=CFG, exclude_spans=spans)
    assert list(hits) == [0] and n_self == 0


def test_self_spans_collects_per_accession():
    a = mk("ACGT" * 50, "ACGT" * 50, acc="CHR2", h_from=1000, klass="self")
    b = mk("ACGT" * 50, "ACGT" * 50, acc="CHR2", h_from=1150, klass="self")
    c = mk("ACGT" * 50, "ACGT" * 50, acc="OTHER", h_from=50, klass="duplicate")
    spans = ots.self_spans([a, b, c])
    assert set(spans) == {"CHR2"}, "only self-classified hits define the on-target locus"
    assert spans["CHR2"] == [(1000, 1349)], "overlapping self spans merge"
    assert ots.in_self_span("CHR2", (1100, 1120), spans)
    assert not ots.in_self_span("CHR2", (9000, 9020), spans)
    assert not ots.in_self_span("OTHER", (60, 80), spans)


# ── predicted amplicons ──────────────────────────────────────────────────────

def _pair_hsp():
    """An HSP covering region 0-199 where both primer windows sit."""
    q = "ACGTTGCA" * 25          # 200 bp
    return mk(q, q, q_from=1, h_from=5001)


def test_both_primers_in_one_duplicate_predicts_an_amplicon():
    h = _pair_hsp()
    amps = ots.predict_amplicons([h], (0, 20), (180, 200), CFG)
    assert len(amps) == 1
    assert amps[0]["product_size"] == 200 and amps[0]["perfect"]


def test_single_primer_in_a_duplicate_predicts_nothing():
    q = "ACGTTGCA" * 5           # only 40 bp — the reverse window is uncovered
    amps = ots.predict_amplicons([mk(q, q)], (0, 20), (180, 200), CFG)
    assert amps == [], "a lone primer in a duplicate is harmless"


def test_amplicon_beyond_max_product_is_dropped():
    cfg = dict(CFG)
    cfg["scoring"] = dict(CFG["scoring"], max_product=50)
    amps = ots.predict_amplicons([_pair_hsp()], (0, 20), (180, 200), cfg)
    assert amps == []


def test_amplicon_requires_intact_three_prime_ends():
    q = "ACGTTGCA" * 25
    # break the forward primer's 3'-terminal base at region coord 19
    hseq = q[:19] + ("A" if q[19] != "A" else "C") + q[20:]
    amps = ots.predict_amplicons([mk(q, hseq)], (0, 20), (180, 200), CFG)
    assert amps == []


# ── region summary ───────────────────────────────────────────────────────────

def test_identical_segments_needs_full_identity_and_length():
    long_perfect = mk("ACGT" * 50, "ACGT" * 50, identity=100.0)
    long_near = mk("ACGT" * 50, "ACGT" * 50, identity=98.0)
    short_perfect = mk("ACGT" * 5, "ACGT" * 5, identity=100.0)
    got = ots.identical_segments([long_perfect, long_near, short_perfect], CFG)
    assert got == [long_perfect]


@pytest.mark.parametrize("args,expected", [
    ((0, 0, 0), "no sites found"),
    ((1, 1, 0), "1 site (1 identical)"),
    ((3, 1, 0), "3 sites (1 identical, 2 near)"),
    ((2, 0, 2), "2 sites (2 PAM unverified)"),
])
def test_summarise(args, expected):
    assert ots.summarise(*args) == expected


# ── PAM resolution by fetching the subject flank ──────────────────────────────

def test_pam_fetch_span_plus_strand_sits_above_the_match():
    site = {"subject_span": (1000, 1019), "h_strand": "+"}
    assert ots.pam_fetch_span(site, 3) == (1020, 1022, False)


def test_pam_fetch_span_minus_strand_sits_below_and_needs_revcomp():
    """The query runs 5'->3', so on a minus hit the PAM is at lower coordinates."""
    site = {"subject_span": (1000, 1019), "h_strand": "-"}
    assert ots.pam_fetch_span(site, 3) == (997, 999, True)


def test_pam_fetch_span_guards_record_start_and_missing_span():
    assert ots.pam_fetch_span({"subject_span": (2, 21), "h_strand": "-"}, 3) is None
    assert ots.pam_fetch_span({"subject_span": None, "h_strand": "+"}, 3) is None


def test_apply_fetched_pam_settles_a_site():
    site = {"pam": None, "pam_ok": False, "pam_unverified": True}
    assert ots.apply_fetched_pam(site, "cgg") is True
    assert site["pam"] == "CGG" and site["pam_ok"] and not site["pam_unverified"]
    assert site["pam_source"] == "fetched"


def test_apply_fetched_pam_records_absence_of_a_pam():
    site = {"pam": None, "pam_ok": False, "pam_unverified": True}
    assert ots.apply_fetched_pam(site, "GTA") is True
    assert site["pam_ok"] is False and site["pam_unverified"] is False


def test_apply_fetched_pam_leaves_site_unverified_on_a_failed_fetch():
    """A failed fetch must stay unknown, never be read as "no PAM"."""
    for bases in ("", None, "NN", "XYZ"):
        site = {"pam": None, "pam_ok": False, "pam_unverified": True}
        assert ots.apply_fetched_pam(site, bases) is False
        assert site["pam_unverified"] is True
