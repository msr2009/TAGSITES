"""
reagent_sequences.py

Pure sequence-assembly logic for CRISPR editing reagents: truncate precomputed
homology arms, build ssODN / PCR primers, calculate Tm, and render ASCII diagrams.

All functions are stateless — no Shiny imports.  Designed for dual use:
in-process calls from reagents_server.py and standalone CLI.

Matt Rich, 2025
"""

import csv
import os
import re
import sys
from html import escape as _esc

from Bio.SeqUtils import MeltingTemp as _mt
import primer3

sys.path.insert(0, os.path.dirname(__file__))
from crispr_util import (
    iupac_to_regex,
    reverse_complement,
    reverse_complement_iupac,
    CODONS,
)


def _rc_case(seq):
    """Reverse complement preserving uppercase/lowercase case encoding."""
    _comp = {'A':'T','C':'G','G':'C','T':'A','N':'N',
             'a':'t','c':'g','g':'c','t':'a','n':'n'}
    return ''.join(_comp.get(c, c) for c in seq[::-1])


# ── Tag library ───────────────────────────────────────────────────────────────

def load_tags(tags_path):
    """Load tables/tags.tsv → OrderedDict {name: dna_sequence}."""
    tags = {}
    with open(tags_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            name = row.get('name', '').strip()
            seq  = row.get('dna_sequence', '').strip()
            if name and seq:
                tags[name] = seq
    return tags


# ── Arm handling ──────────────────────────────────────────────────────────────

def _farthest_mutation_reach(arm, arm_wt, side):
    """bp from the junction-adjacent boundary to the farthest-out base that
    differs between arm and arm_wt (case-insensitive — case encodes exon/intron,
    not edits). Returns 0 if arm_wt is not given, or the two are identical.

    side='left': junction is at the arm's end, so reach = distance from the
    leftmost (farthest) mismatch to the end.
    side='right': junction is at the arm's start, so reach = distance from the
    start to the rightmost (farthest) mismatch.
    """
    if not arm_wt or len(arm_wt) != len(arm):
        return 0
    diffs = [i for i in range(len(arm)) if arm[i].upper() != arm_wt[i].upper()]
    if not diffs:
        return 0
    return (len(arm) - min(diffs)) if side == 'left' else (max(diffs) + 1)


def truncate_arms(left_arm, right_arm, arm_length, left_arm_wt=None, right_arm_wt=None):
    """Truncate precomputed arms to arm_length bp from the insertion junction.

    Takes from the 3' end of left_arm and 5' end of right_arm so that the bases
    closest to the insertion site (including PAM-disrupting mutations) are retained.

    A recut-blocking mutation is not always within arm_length of the junction
    (e.g. the tag insert sits right at the cut site while the PAM it needs to
    disrupt lies farther out). When left_arm_wt/right_arm_wt are supplied, the
    arm on that side is extended out to the mutation instead of truncating it
    away — the homology arm then starts at the mutation, not at arm_length, so
    the returned arm (and the resulting repair template / primer) can be
    longer than arm_length on that side. Without wt arms, this floor can't be
    computed and truncation is exact, as before.

    Raises ValueError when the required length (arm_length, or farther if a
    mutation demands it) exceeds the stored arm length.
    """
    left_needed  = max(arm_length, _farthest_mutation_reach(left_arm, left_arm_wt, 'left'))
    right_needed = max(arm_length, _farthest_mutation_reach(right_arm, right_arm_wt, 'right'))
    if left_needed > len(left_arm) or right_needed > len(right_arm):
        raise ValueError(
            'Requested arm_length {} (extended to {}/{} to reach a recut-blocking '
            'mutation) exceeds stored arm length {} / {} — rerun the reagents '
            'pipeline with a larger arm_length.'.format(
                arm_length, left_needed, right_needed, len(left_arm), len(right_arm))
        )
    return left_arm[-left_needed:], right_arm[:right_needed]


def truncate_arms_with_tolerance(left_arm, right_arm, arm_length, tolerance=0.2,
                                 left_arm_wt=None, right_arm_wt=None):
    """Like truncate_arms, but each arm keeps extra headroom on its far/outer
    side (up to tolerance * arm_length, bounded by the stored arm length).

    The insertion-adjacent boundary never moves — left_arm still ends exactly
    at the junction, right_arm still starts exactly at the junction — only the
    far end gets extra bases so a primer3 search (see
    plasmid_assembly.sapi_amplification_primers) has room to pick a
    better-quality primer without ever touching the near/junction end.

    As in truncate_arms, a recut-blocking mutation farther than arm_length from
    the junction (detected via left_arm_wt/right_arm_wt) pushes that side's
    length out to the mutation, independent of the tolerance headroom.
    """
    left_needed  = max(arm_length, _farthest_mutation_reach(left_arm, left_arm_wt, 'left'))
    right_needed = max(arm_length, _farthest_mutation_reach(right_arm, right_arm_wt, 'right'))
    if left_needed > len(left_arm) or right_needed > len(right_arm):
        raise ValueError(
            'Requested arm_length {} (extended to {}/{} to reach a recut-blocking '
            'mutation) exceeds stored arm length {} / {} — rerun the reagents '
            'pipeline with a larger arm_length.'.format(
                arm_length, left_needed, right_needed, len(left_arm), len(right_arm))
        )
    base_widened = int(arm_length * (1 + tolerance))
    widened_left  = min(max(base_widened, left_needed), len(left_arm))
    widened_right = min(max(base_widened, right_needed), len(right_arm))
    return left_arm[-widened_left:], right_arm[:widened_right]


# ── Tm calculation ────────────────────────────────────────────────────────────

def calc_tm(seq):
    """Nearest-neighbour Tm (°C) for a DNA sequence; returns 0.0 if too short."""
    seq = seq.upper().replace('U', 'T')
    if len(seq) < 8:
        return 0.0
    try:
        return float(_mt.Tm_NN(seq))
    except Exception:
        return 0.0


# ── WT arm reconstruction (legacy-TSV fallback) ──────────────────────────────

_SYN_DESC_RE = re.compile(r'^([ACGTacgt]{3})>[ACGTacgt]{3} at genomic pos (\d+)$')
_MUT_DESC_RE = re.compile(r'^pos (\d+) ([ACGTacgt])>[ACGTacgt]')


def reconstruct_wt_arms(left_arm, right_arm, insert_pos, method, mutation_desc):
    """Rebuild the pre-mutation (WT) arms by reverting the recorded PAM edit.

    Fallback for reagent TSVs written before the left_arm_wt/right_arm_wt columns
    existed. The mutated arms plus the mutation_desc string fully determine the WT
    bases: syn_1 records the original codon and its genomic position; mut_1 records
    a single genomic position and original base. Genomic→arm mapping uses insert_pos
    (right_arm starts at insert_pos; left_arm ends there).

    Returns (left_wt, right_wt). When there is no mutation to revert the arms are
    already WT and returned unchanged; returns None only if the description can't
    be parsed (caller then keeps legacy PAM-only behaviour).
    """
    if method in ('none', 'insertion', ''):
        return left_arm, right_arm   # no PAM edit → arms are already WT
    if method not in ('syn_1', 'mut_1'):
        # multi-change methods (e.g. syn_seed) aren't reconstructable from the
        # description here — but they only appear in newer TSVs that carry the
        # left_arm_wt/right_arm_wt columns, so this path isn't reached for them.
        return None

    reverts = []   # (genomic_pos, original_base)
    m = _SYN_DESC_RE.match(mutation_desc or '')
    if m:
        orig_codon, cs = m.group(1).upper(), int(m.group(2))
        reverts = [(cs + k, orig_codon[k]) for k in range(3)]
    else:
        m = _MUT_DESC_RE.match(mutation_desc or '')
        if not m:
            return None
        reverts = [(int(m.group(1)), m.group(2).upper())]

    left_gstart = insert_pos - len(left_arm)   # genomic coord of left_arm[0]
    lchars, rchars = list(left_arm), list(right_arm)
    for gpos, orig in reverts:
        li = gpos - left_gstart
        if 0 <= li < len(lchars):
            lchars[li] = orig.lower() if lchars[li].islower() else orig.upper()
        ri = gpos - insert_pos
        if 0 <= ri < len(rchars):
            rchars[ri] = orig.lower() if rchars[ri].islower() else orig.upper()
    return ''.join(lchars), ''.join(rchars)


# ── Recut risk check ─────────────────────────────────────────────────────────

def check_recut_risk(left_arm, right_arm, insert_seq, spacer, pam_seq, guide_strand):
    """Return True if inserting insert_seq reconstitutes the guide target at the repair junction.

    Builds a window spanning both arm/insert junctions and searches for the
    original spacer+PAM sequence.  The left_arm and right_arm should already
    contain any PAM-disrupting mutations — this check catches the separate case
    where the insert itself provides bases that complete the PAM.
    """
    spacer = spacer.upper()
    pam_seq = pam_seq.upper()
    insert_seq = insert_seq.upper()

    # Window must span spacer + PAM across both junctions
    n = len(spacer) + len(pam_seq)
    junction = left_arm[-(n + len(insert_seq)):].upper() + insert_seq + right_arm[:(n + len(insert_seq))].upper()

    # For + strand: look for spacer immediately followed by PAM on the fwd strand
    # For - strand: look for RC(PAM) + RC(spacer) on the fwd strand
    if guide_strand == '+':
        pattern = re.escape(spacer) + iupac_to_regex(pam_seq)
    else:
        pattern = iupac_to_regex(reverse_complement_iupac(pam_seq)) + re.escape(reverse_complement(spacer))

    return bool(re.search(pattern, junction))


# ── ssODN assembly ────────────────────────────────────────────────────────────

def build_ssodn(left_arm, right_arm, insert_seq, cut_pos, insert_pos, consider_strand):
    """Build a single-stranded oligo donor from truncated arms and an insert.

    Strand selection (when consider_strand=True):
      cut 5' of insert (cut_pos < insert_pos) → reverse-complement the ssODN
      cut at or 3' of insert                  → keep forward strand

    Returns (sequence, total_length, strand_used '+'/'-').
    """
    oligo = left_arm + insert_seq + right_arm
    if consider_strand and cut_pos < insert_pos:
        oligo = _rc_case(oligo)   # preserves exonic/intronic case encoding
        strand = '-'
    else:
        strand = '+'
    return oligo, len(oligo), strand


# ── PCR primer design ─────────────────────────────────────────────────────────

def _tm_growth_bind_regions(template, tm_target):
    """Hand-rolled Tm-growth fallback: grow each binding region from its
    anchored end until Tm >= tm_target (or the template is exhausted).
    Returns (fwd_bind, rev_bind_rc) — same as the pre-primer3 implementation.
    """
    k_fwd = 12
    while k_fwd < len(template) and calc_tm(template[:k_fwd]) < tm_target:
        k_fwd += 1
    k_rev = 12
    while k_rev < len(template) and calc_tm(template[-k_rev:]) < tm_target:
        k_rev += 1
    return template[:k_fwd], reverse_complement(template[-k_rev:])


def _primer3_anchored_single(template, force_pos, pick_left, primer_opt_tm, max_len=36):
    """Design one primer whose start is pinned exactly at force_pos, letting
    primer3 choose the other end/length for Tm and primer quality (hairpins,
    GC clamp) instead of a naive Tm-growth loop.

    Uses primer3-py's true single-sided PRIMER_PICK_LEFT/RIGHT_PRIMER=0 mode
    (measured ~100x slower per call than an ordinary pair in earlier
    profiling) rather than the pair-then-discard trick used elsewhere in this
    module — acceptable here since this is called per user-selected reagent
    row (a handful of calls), not per genome-wide site.

    Returns (sequence, tm) or None if primer3 found no valid design.
    """
    seq_args = {'SEQUENCE_ID': 'design', 'SEQUENCE_TEMPLATE': template}
    if pick_left:
        seq_args['SEQUENCE_FORCE_LEFT_START'] = force_pos
    else:
        seq_args['SEQUENCE_FORCE_RIGHT_START'] = force_pos
    global_args = {
        'PRIMER_OPT_TM': primer_opt_tm,
        # Wider than the ±3 window used elsewhere in this module — the anchor
        # constraint already removes most of primer3's search freedom (only
        # length can vary, not position), so a tight Tm window on top of that
        # made the fallback fire far more often than warranted (profiled: ~30%
        # fallback rate at ±3 vs. ~15-20% at ±8 on random test sequences).
        'PRIMER_MIN_TM': primer_opt_tm - 8,
        'PRIMER_MAX_TM': primer_opt_tm + 8,
        'PRIMER_MAX_SIZE': max_len,
        'PRIMER_PICK_LEFT_PRIMER': 1 if pick_left else 0,
        'PRIMER_PICK_RIGHT_PRIMER': 0 if pick_left else 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_NUM_RETURN': 1,
    }
    res = primer3.bindings.design_primers(seq_args, global_args)
    key = 'PRIMER_LEFT_0_SEQUENCE' if pick_left else 'PRIMER_RIGHT_0_SEQUENCE'
    seq = res.get(key)
    if not seq:
        return None
    tm_key = 'PRIMER_LEFT_0_TM' if pick_left else 'PRIMER_RIGHT_0_TM'
    return seq, res.get(tm_key)


def design_pcr_primers(left_arm, right_arm, template, tm_target, phos,
                       cut_pos, insert_pos, max_bind_len=36):
    """Design PCR primers that add homology arms to a tag template.

    Forward primer = left_arm + a primer3-designed binding region pinned to
    start exactly at template[0]. Reverse primer = RC(right_arm) + a
    primer3-designed binding region pinned to end exactly at the last base of
    template. Both anchors match the arm/template boundary exactly — primer3
    only picks the binding region's length/Tm/quality (hairpins, GC clamp),
    replacing the previous naive Tm-growth loop.

    Falls back to the old hand-rolled Tm-growth loop (unconditionally, no
    anchoring) if primer3 finds no valid design for either side — e.g. a
    template too short for primer3's minimum primer size.

    Phosphorylation (for lambda-exo ssDNA): same strand logic as ssODN.
      cut 5' of insert → phosphorylate forward primer (to produce antisense ssDNA)
      cut at/3' of insert → phosphorylate reverse primer (to produce sense ssDNA)
    /5Phos/ IDT modifier is prepended to the sequence when phos=True.

    Returns ((fwd_suffix, fwd_seq), (rev_suffix, rev_seq), meta).
    fwd/rev_suffix are '_F' / '_R' (already include phos tag in seq, not suffix).
    meta = {"used_fallback": bool}.
    """
    template = template.upper()
    used_fallback = False

    fwd_result = _primer3_anchored_single(template, 0, True, tm_target, max_bind_len)
    rev_result = _primer3_anchored_single(template, len(template) - 1, False, tm_target, max_bind_len)

    if fwd_result and rev_result:
        fwd_bind = fwd_result[0]
        rev_bind = reverse_complement(rev_result[0])
    else:
        used_fallback = True
        fwd_bind, rev_bind = _tm_growth_bind_regions(template, tm_target)

    fwd_seq = left_arm + fwd_bind
    rev_seq = _rc_case(right_arm) + rev_bind   # preserve exonic/intronic case in arm overhang

    # Phosphorylate the primer that produces the desired repair strand
    if phos:
        if cut_pos < insert_pos:
            # want antisense ssDNA after lambda exo → phos fwd primer
            fwd_seq = '/5Phos/' + fwd_seq
        else:
            # want sense ssDNA → phos rev primer
            rev_seq = '/5Phos/' + rev_seq

    return ('_F', fwd_seq), ('_R', rev_seq), {'used_fallback': used_fallback}


# ── Genotyping (screening) primer design ─────────────────────────────────────

def _primer3_pair(template, target_start, target_len, excluded_regions,
                  primer_opt_tm, product_opt_size, product_size_range):
    """Design one ordinary primer3 pair on template; returns dict or None if none found."""
    seq_args = {
        'SEQUENCE_ID': 'design',
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_TARGET': [target_start, max(target_len, 1)],
        'SEQUENCE_EXCLUDED_REGION': excluded_regions,
    }
    global_args = {
        'PRIMER_OPT_TM': primer_opt_tm,
        'PRIMER_MIN_TM': primer_opt_tm - 3,
        'PRIMER_MAX_TM': primer_opt_tm + 3,
        'PRIMER_PRODUCT_OPT_SIZE': product_opt_size,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
        'PRIMER_NUM_RETURN': 1,
    }
    try:
        res = primer3.bindings.design_primers(seq_args, global_args)
    except OSError:
        # region too small/degenerate for primer3's geometry (e.g. a site near
        # the edge of the genomic record) — equivalent to "no valid design found"
        return None
    if res.get('PRIMER_PAIR_NUM_RETURNED', 0) < 1:
        return None
    return {
        'fwd_seq': res['PRIMER_LEFT_0_SEQUENCE'],
        'fwd_tm': res['PRIMER_LEFT_0_TM'],
        'rev_seq': res['PRIMER_RIGHT_0_SEQUENCE'],
        'rev_tm': res['PRIMER_RIGHT_0_TM'],
        'product_size': res['PRIMER_PAIR_0_PRODUCT_SIZE'],
    }


def _primer3_single(template, region_start, region_len, excluded_region, pick_left,
                    primer_opt_tm, product_size_range=None):
    """Design a single primer (left-only or right-only) on template; returns sequence or None.

    Designed as an ordinary primer3 *pair* (fast — primer3-py's single-sided
    PRIMER_PICK_LEFT/RIGHT_PRIMER=0 mode is ~100x slower per call, confirmed by
    profiling) and only the requested side of the pair is kept; the other side
    is discarded — so the *pair's* product size is irrelevant. product_size_range
    defaults to an effectively-unconstrained range so a small SEQUENCE_INCLUDED_REGION
    (smaller than primer3's built-in default product-size options) doesn't spuriously
    fail; pass an explicit range only if the discarded side's placement actually matters.
    """
    seq_args = {
        'SEQUENCE_ID': 'design',
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_INCLUDED_REGION': [region_start, region_len],
        'SEQUENCE_EXCLUDED_REGION': [excluded_region],
    }
    if product_size_range is None:
        # product size must be >= primer3's default max primer size (27);
        # the pair's product size is otherwise irrelevant since we discard
        # the unused side.
        lo = min(27, len(template))
        product_size_range = [lo, max(lo, len(template))]
    global_args = {
        'PRIMER_OPT_TM': primer_opt_tm,
        'PRIMER_MIN_TM': primer_opt_tm - 3,
        'PRIMER_MAX_TM': primer_opt_tm + 3,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
        'PRIMER_NUM_RETURN': 1,
    }
    try:
        res = primer3.bindings.design_primers(seq_args, global_args)
    except OSError:
        return None
    if res.get('PRIMER_PAIR_NUM_RETURNED', 0) < 1:
        return None
    key = 'PRIMER_LEFT_0_SEQUENCE' if pick_left else 'PRIMER_RIGHT_0_SEQUENCE'
    return res.get(key)


def _primer3_paired_to_fixed(edited_seq, fixed_seq, fixed_is_forward,
                             primer_opt_tm, product_opt_size, product_size_range):
    """Design the true partner primer for a fixed (already-chosen) primer, as a real pair.

    fixed_seq is pinned via SEQUENCE_PRIMER (forward-fixed) or SEQUENCE_PRIMER_REVCOMP
    (reverse-fixed) so primer3 evaluates real pair metrics (Tm balance, complementarity,
    product size) against it, rather than combining two independently-designed primers.
    """
    seq_args = {'SEQUENCE_ID': 'design', 'SEQUENCE_TEMPLATE': edited_seq}
    if fixed_is_forward:
        seq_args['SEQUENCE_PRIMER'] = fixed_seq
    else:
        seq_args['SEQUENCE_PRIMER_REVCOMP'] = fixed_seq
    global_args = {
        'PRIMER_OPT_TM': primer_opt_tm,
        'PRIMER_MIN_TM': primer_opt_tm - 10,
        'PRIMER_MAX_TM': primer_opt_tm + 10,
        'PRIMER_PRODUCT_OPT_SIZE': product_opt_size,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
        'PRIMER_NUM_RETURN': 1,
    }
    try:
        res = primer3.bindings.design_primers(seq_args, global_args)
    except OSError:
        return None
    if res.get('PRIMER_PAIR_NUM_RETURNED', 0) < 1:
        return None
    return {
        'fwd_seq': res['PRIMER_LEFT_0_SEQUENCE'],
        'fwd_tm': res['PRIMER_LEFT_0_TM'],
        'rev_seq': res['PRIMER_RIGHT_0_SEQUENCE'],
        'rev_tm': res['PRIMER_RIGHT_0_TM'],
        'product_size': res['PRIMER_PAIR_0_PRODUCT_SIZE'],
    }


def design_genotyping_primers(left_flank, right_flank, insert_sequence,
                              internal_threshold=500, primer_opt_tm=60.0,
                              product_opt_size=200, flank_min=50, flank_max=150):
    """Design genotyping (screening) PCR primers flanking a CRISPR insertion site.

    left_flank / right_flank: genomic sequence on each side of the insertion site
    (at least flank_max + ~60bp margin recommended).

    No insert sequence: one external pair spanning the insertion site, with both
    primers held flank_min-flank_max bp away — verifies the WT locus / detects any
    size change at the site.

    Insert sequence given: the external pair above, plus a 5' and 3' junction pair,
    attempted regardless of insert length. Each junction pair is designed as a
    genuine primer3 pair — the genomic-side primer is designed first (single-primer
    mode), then pinned via SEQUENCE_PRIMER/SEQUENCE_PRIMER_REVCOMP so primer3 designs
    a true, Tm-/dimer-compatible internal partner on the in-silico edited allele
    sequence (left_flank + insert_sequence + right_flank) — not matched after the fact.
    A junction pair is simply absent from the result if the insert is too short for
    primer3 to find a valid internal primer.

    internal_threshold is accepted for backward compatibility but no longer gates
    junction-primer design (kept to avoid changing the CLI/task interface).

    Returns {amplicon_type: {fwd_seq, fwd_tm, rev_seq, rev_tm, product_size}}.
    Amplicon types missing from the dict mean primer3 found no valid design.
    """
    left_flank = left_flank.upper()
    right_flank = right_flank.upper()
    insert_sequence = insert_sequence.upper()
    insert_len = len(insert_sequence)
    edited_seq = left_flank + insert_sequence + right_flank

    results = {}

    # ── External pair: spans the whole insert, primers held off the insert ──
    # Excluded-region lengths are clamped to the available flank on each side —
    # primer3 raises OSError if a region extends past the template's end, which
    # happens for sites near the edge of the genomic record (short left/right flank).
    insert_start = len(left_flank)
    left_excl_start = max(insert_start - flank_min, 0)
    left_excl_len   = min(flank_min, len(edited_seq) - left_excl_start)
    right_excl_start = insert_start + insert_len
    right_excl_len   = min(flank_min, len(edited_seq) - right_excl_start)
    excluded_regions = []
    if left_excl_len > 0:
        excluded_regions.append([left_excl_start, left_excl_len])
    if right_excl_len > 0:
        excluded_regions.append([right_excl_start, right_excl_len])

    ext = _primer3_pair(
        edited_seq, insert_start, max(insert_len, 1),
        excluded_regions=excluded_regions,
        primer_opt_tm=primer_opt_tm,
        product_opt_size=product_opt_size,
        product_size_range=[insert_len + 2 * flank_min, insert_len + 2 * flank_max],
    )
    if ext:
        results['external'] = ext

    if not insert_sequence:
        return results

    # ── 5' junction: fixed genomic forward primer, primer3 designs internal reverse ──
    fwd_region_len = min(len(left_flank), flank_max + 60)
    fwd_region_start = len(left_flank) - fwd_region_len
    fwd_fixed = _primer3_single(
        left_flank, fwd_region_start, fwd_region_len,
        excluded_region=[len(left_flank) - flank_min, flank_min],
        pick_left=True, primer_opt_tm=primer_opt_tm,
    )
    if fwd_fixed:
        j5 = _primer3_paired_to_fixed(
            edited_seq, fwd_fixed, fixed_is_forward=True,
            primer_opt_tm=primer_opt_tm, product_opt_size=product_opt_size,
            product_size_range=[flank_min, flank_max + 400],
        )
        if j5:
            results['5p_junction'] = j5

    # ── 3' junction: fixed genomic reverse primer, primer3 designs internal forward ──
    rev_region_len = min(len(right_flank), flank_max + 60)
    rev_fixed = _primer3_single(
        right_flank, 0, rev_region_len,
        excluded_region=[0, flank_min],
        pick_left=False, primer_opt_tm=primer_opt_tm,
    )
    if rev_fixed:
        j3 = _primer3_paired_to_fixed(
            edited_seq, rev_fixed, fixed_is_forward=False,
            primer_opt_tm=primer_opt_tm, product_opt_size=product_opt_size,
            product_size_range=[flank_min, flank_max + 400],
        )
        if j3:
            results['3p_junction'] = j3

    return results


# ── Translation helper ────────────────────────────────────────────────────────

def _build_aa_line(seq_raw, local_ins):
    """Return a string of same length as seq_raw with AAs centered in coding triplets.

    Uppercase = exonic (coding); lowercase = intronic (skipped).
    Codons are read outward from local_ins (the insert/junction point).
    """
    ann = [' '] * len(seq_raw)
    seq_raw = list(seq_raw)  # mutable so stop codon bases can be forced uppercase

    def _annotate(cod_pos, seq_raw):
        codon = ''.join(seq_raw[k].upper() for k in cod_pos)
        aa = CODONS.get(codon, 'X')
        if aa == '_':
            aa = '*'
            for k in cod_pos:   # stop codons always display uppercase
                seq_raw[k] = seq_raw[k].upper()
        ann[cod_pos[1]] = aa

    # Right side: codons read forward from local_ins
    i = local_ins
    while i < len(seq_raw):
        cod_pos = []
        j = i
        while j < len(seq_raw) and len(cod_pos) < 3:
            ch = seq_raw[j]
            if ch.isupper():
                cod_pos.append(j)
                j += 1
            elif ch.islower():
                j += 1   # skip intronic
            else:
                break
        if len(cod_pos) == 3:
            _annotate(cod_pos, seq_raw)
            i = j
        else:
            break

    # After the coding region, check if the very next codon (ignoring case) is a stop.
    # Only do this when no more uppercase bases exist to the right — i.e. we are at the
    # end of the last exon, not in an intron between exons (where TAA/TAG/TGA can appear
    # at splice sites).
    no_more_exon = not any(seq_raw[k].isupper() for k in range(i, len(seq_raw)))
    if no_more_exon:
        stop_pos = []
        k = i
        while k < len(seq_raw) and len(stop_pos) < 3:
            if seq_raw[k].upper() in 'ACGT':
                stop_pos.append(k)
            k += 1
        if len(stop_pos) == 3:
            codon = ''.join(seq_raw[p].upper() for p in stop_pos)
            if CODONS.get(codon) == '_':
                ann[stop_pos[1]] = '*'
                for p in stop_pos:
                    seq_raw[p] = seq_raw[p].upper()

    # Left side: codons read backward from local_ins
    i = local_ins - 1
    while i >= 0:
        cod_pos = []
        j = i
        while j >= 0 and len(cod_pos) < 3:
            ch = seq_raw[j]
            if ch.isupper():
                cod_pos.insert(0, j)   # prepend to keep left-to-right order
                j -= 1
            elif ch.islower():
                j -= 1
            else:
                break
        if len(cod_pos) == 3:
            _annotate(cod_pos, seq_raw)
            i = j
        else:
            break

    return ''.join(ann), ''.join(seq_raw)


# ── ASCII diagram ─────────────────────────────────────────────────────────────

def ascii_diagram(guide_row, left_arm, right_arm, left_wt=None, right_wt=None,
                  insert_seq=None, wt_is_true=False):
    """Return an HTML string for embedding in a <pre> block.

    Rows (top to bottom):
      1. WT genomic sequence (PAM restored, blue; extends into wt context)
      2. Translation (single-letter AA centered in each coding triplet)
      3. Repair template (arm with mutations red; unmutated PAM blue; context dimmed)
      4. Annotation (^ insert site green+bold, | cut site)
      5. Guide (>>>> PAM for + strand; PAM <<<< for - strand)

    left_wt / right_wt: WT (pre-mutation) arms trimmed to arm_length+10 for the WT
    display row and context.  If omitted, same as left_arm / right_arm.

    wt_is_true: True when left_wt/right_wt are the genuine unmutated arms (the
    left_arm_wt/right_arm_wt TSV columns).  In that case the repair row highlights
    every base that differs from WT — catching synonymous PAM-disrupting codon
    changes that alter a base outside the 3 bp PAM.  When False (older TSVs with no
    WT columns), it falls back to reddening the PAM positions when a mutation was
    recorded, and the WT row is only accurate at the PAM (legacy behaviour).

    insert_seq: ssODN tag/insert sequence. When given, the literal insert is
    spliced into the repair-template row (uppercase) between the two arms,
    and the arm bases in that row are shown lowercase — except any
    mutation within the arm, which stays uppercase and red.
    """
    insert_pos    = int(guide_row['insert_pos'])
    cut_pos       = int(guide_row['cut_pos'])
    pam_fwd_start = int(guide_row['pam_fwd_start'])
    pam_seq       = str(guide_row['pam_seq'])
    pam_len       = len(pam_seq)
    spacer_len    = len(str(guide_row['spacer']))
    strand        = str(guide_row['guide_strand'])
    recut_method  = str(guide_row.get('recut_block_method', 'none') or 'none')
    mutation_desc = str(guide_row.get('mutation_desc', '') or '')
    # a PAM mutation exists when method is syn_1 or mut_1 (not insertion or none)
    mutated = recut_method in ('syn_1', 'mut_1')
    dist    = abs(cut_pos - insert_pos)

    insert_seq  = insert_seq.upper() if insert_seq else ''
    ins_len     = len(insert_seq)
    ssodn_mode  = bool(insert_seq)

    # Display window is based on the WT arms (wider than truncated repair arms)
    l_wt = left_wt if left_wt is not None else left_arm
    r_wt = right_wt if right_wt is not None else right_arm
    MAX  = 100   # hard cap per side
    l_show = l_wt[-MAX:]
    r_show = r_wt[:MAX]

    win_len   = len(l_show) + len(r_show)
    local_ins = len(l_show)

    gds             = insert_pos - local_ins
    local_cut       = cut_pos       - gds
    local_pam_start = pam_fwd_start - gds

    def _in_pam(i):
        return local_pam_start <= i < local_pam_start + pam_len

    # WT fwd-strand bases at PAM positions (left-to-right in genomic coords).
    # pam_seq is 5'→3' on the guide strand.
    # For + strand: guide strand = fwd strand, so pam_seq IS the fwd bases.
    # For - strand: guide strand = rev strand, so fwd bases = RC(pam_seq).
    wt_pam_fwd = pam_seq.upper() if strand == '+' else reverse_complement(pam_seq).upper()

    # Arm offset within the display window
    left_ctx  = local_ins - len(left_arm)    # chars of WT context to the left of the arm
    right_ctx = len(r_show) - len(right_arm) # chars of WT context to the right of the arm

    pfx = '   '   # 3-char left prefix (space or '...')

    # ── Row 2: translation — computed first so stop codons can be forced uppercase
    wt_raw_base = l_show + r_show
    aa_raw, wt_raw = _build_aa_line(wt_raw_base, local_ins)

    # ── Row 1: WT sequence ────────────────────────────────────────────────────
    # wt_raw may differ from wt_raw_base at stop codon positions (forced uppercase)
    wt_parts = [pfx]
    for i, ch in enumerate(wt_raw):
        if _in_pam(i):
            k = i - local_pam_start
            base = wt_pam_fwd[k] if k < len(wt_pam_fwd) else ch.upper()
            wt_parts.append('<span style="color:#4a90d9">{}</span>'.format(_esc(base)))
        else:
            wt_parts.append(_esc(ch))

    # ── Row 3: repair template ────────────────────────────────────────────────
    # Context bases (outside the arm) shown in dim gray; arm bases shown normally
    # with mutated PAM positions in red, unmutated PAM in blue.
    def _arm_ch(i):
        """Character from the repair arm at display position i."""
        if local_ins > i >= left_ctx:
            return left_arm[i - left_ctx]
        if local_ins <= i < local_ins + len(right_arm):
            return right_arm[i - local_ins]
        return wt_raw[i]   # context: same as WT

    rep_parts = [pfx]
    for i, ch in enumerate(wt_raw):
        in_ctx = i < left_ctx or i >= local_ins + len(right_arm)
        arm_ch = _arm_ch(i)
        # a base is "mutated" if it differs from WT (case-insensitive: case encodes
        # exon/intron, not edits).  With genuine WT arms this catches synonymous
        # changes outside the PAM; otherwise fall back to the recorded PAM mutation.
        if wt_is_true:
            is_mut = (not in_ctx) and arm_ch.upper() != wt_raw[i].upper()
        else:
            is_mut = (not in_ctx) and _in_pam(i) and mutated
        if in_ctx:
            rep_parts.append('<span style="color:#bbb">{}</span>'.format(_esc(arm_ch)))
        elif is_mut:
            disp = arm_ch.upper() if ssodn_mode else arm_ch
            rep_parts.append('<span style="color:#c0392b">{}</span>'.format(_esc(disp)))
        elif _in_pam(i):
            disp = arm_ch.lower() if ssodn_mode else arm_ch
            rep_parts.append('<span style="color:#4a90d9">{}</span>'.format(_esc(disp)))
        else:
            disp = arm_ch.lower() if ssodn_mode else arm_ch
            rep_parts.append(_esc(disp))

    # ── Row 4: annotation (^ insert, | cut) ──────────────────────────────────
    # Both markers are shifted left by 0.5ch so they appear between two bases
    # rather than under the right-hand one (monospace cells are whole-width only).
    _HALF = 'position:relative;left:-0.5ch'
    ann_parts = [pfx]
    for i in range(win_len):
        if i == local_ins and 0 <= local_ins < win_len:
            ann_parts.append(
                '<span style="color:#27ae60;font-weight:bold;{}">' \
                '^</span>'.format(_HALF)
            )
        elif i == local_cut and 0 <= local_cut < win_len:
            ann_parts.append('<span style="{}">|</span>'.format(_HALF))
        else:
            ann_parts.append(' ')

    # ── Row 5: guide line ─────────────────────────────────────────────────────
    guide_slots = [''] * win_len
    if strand == '+':
        # + strand: >>>>PAM  (PAM label = pam_seq left-to-right in fwd coords)
        pam_label = pam_seq.upper()
        for k in range(spacer_len):
            p = local_pam_start - spacer_len + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#888">&gt;</span>'
        for k, ch in enumerate(pam_label):
            p = local_pam_start + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#4a90d9">{}</span>'.format(_esc(ch))
    else:
        # - strand: PAM<<<<  (PAM label reversed so letters match fwd-strand positions)
        # pam_seq is 5'→3' on the minus strand (right-to-left in fwd coords).
        # Reversing gives the same letters left-to-right in fwd display order.
        pam_label = pam_seq[::-1].lower()
        for k, ch in enumerate(pam_label):
            p = local_pam_start + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#4a90d9">{}</span>'.format(_esc(ch))
        for k in range(spacer_len):
            p = local_pam_start + pam_len + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#888">&lt;</span>'

    # ── Splice the literal insert into the display, between the two arms ────────
    # Inserting fresh tokens before index local_ins in every row keeps all five
    # rows column-aligned: each row gains exactly ins_len display positions.
    if ins_len:
        aa_raw = aa_raw[:local_ins] + (' ' * ins_len) + aa_raw[local_ins:]
        wt_parts[1 + local_ins:1 + local_ins] = [
            '<span style="color:#bbb">-</span>'
        ] * ins_len
        rep_parts[1 + local_ins:1 + local_ins] = [_esc(ch) for ch in insert_seq]
        ann_parts[1 + local_ins:1 + local_ins] = [' '] * ins_len
        guide_slots[local_ins:local_ins] = [''] * ins_len

    aa_line   = pfx + _esc(aa_raw)
    wt_line   = ''.join(wt_parts)
    rep_line  = ''.join(rep_parts)
    ann_line  = ''.join(ann_parts)
    guide_line = pfx + ''.join(s if s else ' ' for s in guide_slots)

    # ── Footer ────────────────────────────────────────────────────────────────
    if mutation_desc and recut_method not in ('none', ''):
        recut_str = '{}: {}'.format(_esc(recut_method), _esc(mutation_desc))
    elif recut_method == 'none':
        recut_str = 'no PAM disruption'
    else:
        recut_str = _esc(recut_method)

    footer = 'cut&#8596;insert: {} bp  PAM: {}  strand: {}  recut: {}'.format(
        dist, _esc(pam_seq.upper()), strand, recut_str
    )

    return '\n'.join([wt_line, aa_line, rep_line, ann_line, guide_line, footer])


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == '__main__':
    from argparse import ArgumentParser
    import pandas as pd

    parser = ArgumentParser(
        description='Assemble editing reagents from a design_tag_reagents TSV.'
    )
    parser.add_argument('--reagents_tsv',  required=True, help='Path to .reagents.tsv')
    parser.add_argument('--repair_type',   default='ssodn',
                        choices=['ha', 'ssodn', 'amplicon'],
                        help='Repair strategy (default: ssodn)')
    parser.add_argument('--arm_length',    type=int, default=60,
                        help='Homology arm length in bp (default: 60)')
    parser.add_argument('--insert_seq',    default='', help='Tag/insert DNA sequence')
    parser.add_argument('--consider_strand', action='store_true', default=True,
                        help='Apply repair-strand logic for ssODN (default: on)')
    parser.add_argument('--pcr_template',  default='', help='Template sequence for amplicon')
    parser.add_argument('--pcr_tm',        type=float, default=60.0,
                        help='Primer Tm target (°C) for amplicon (default: 60)')
    parser.add_argument('--pcr_phos',      action='store_true',
                        help='Add 5-phosphate to one PCR primer for lambda-exo')
    parser.add_argument('--output',        required=True, help='Output TSV path')
    args = parser.parse_args()

    df   = pd.read_csv(args.reagents_tsv, sep='\t')
    rows = []
    has_wt = 'left_arm_wt' in df.columns and 'right_arm_wt' in df.columns
    for _, row in df.iterrows():
        try:
            left, right = truncate_arms(
                row['left_arm'], row['right_arm'], args.arm_length,
                left_arm_wt=row['left_arm_wt'] if has_wt else None,
                right_arm_wt=row['right_arm_wt'] if has_wt else None,
            )
        except ValueError as e:
            print('WARNING: {}'.format(e), file=sys.stderr)
            continue

        site  = '{}{}'.format(row['amino_acid'], row['residue_index'])
        crrna = str(row['spacer'])
        cut   = int(row['cut_pos'])
        ins   = int(row['insert_pos'])

        if args.repair_type == 'ssodn':
            seq, total, strand = build_ssodn(
                left, right, args.insert_seq.upper(), cut, ins, args.consider_strand
            )
            rows.append({'name': '{}_ssODN'.format(site), 'sequence': seq,
                         'total_len': total, 'strand': strand})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

        elif args.repair_type == 'amplicon' and args.pcr_template:
            (_, fwd), (_, rev), meta = design_pcr_primers(
                left, right, args.pcr_template.upper(), args.pcr_tm,
                args.pcr_phos, cut, ins
            )
            if meta['used_fallback']:
                print('WARNING: {} primer3 design failed; fell back to Tm-growth '
                      'primers'.format(site), file=sys.stderr)
            rows.append({'name': '{}_F'.format(site), 'sequence': fwd})
            rows.append({'name': '{}_R'.format(site), 'sequence': rev})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

        else:  # ha only
            rows.append({'name': '{}_5HA'.format(site), 'sequence': left})
            rows.append({'name': '{}_3HA'.format(site), 'sequence': right})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep='\t', index=False)
    print('Wrote {} rows to {}'.format(len(out_df), args.output), file=sys.stderr)
