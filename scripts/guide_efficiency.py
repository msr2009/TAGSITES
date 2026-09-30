"""
guide_efficiency.py

Rule Set 3 (RS3) on-target efficiency scoring for CRISPR guides.

RS3 (DeWeirdt et al., Nat Commun 2022) predicts SpCas9 *cutting* activity from a
30-nt sequence context. Scores are approximately z-scored log-fold-changes from
pooled human/mouse screens: centered near 0, unbounded, and only meaningful as a
relative ranking between guides.

This module is display-only by design — nothing here selects, filters, or reorders
guides. It must also import cleanly when the optional `rs3` package is absent, so
the rs3 import happens lazily inside load_rs3() and every failure degrades to "no
score" rather than raising.

RS3 context layout (guide orientation, 30 nt total):
    [0:4]   4 nt upstream of the spacer
    [4:24]  20 nt protospacer
    [24:27] 3 nt PAM
    [27:30] 3 nt downstream

Matt Rich, 2025
"""

import sys
from bisect import bisect_left, bisect_right

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from crispr_util import reverse_complement

RS3_TRACRS = ('Hsu2013', 'Chen2013')
RS3_CONTEXT_LEN = 30
_UPSTREAM = 4       # nt of context 5' of the spacer
_DOWNSTREAM = 3     # nt of context 3' of the PAM


def rs3_supported(pam, guide_length):
    """True when RS3's 30-mer layout is valid for these guide settings (SpCas9 NGG, 20 nt)."""
    return str(pam).upper() == 'NGG' and int(guide_length) == 20


def guide_context_30mer(seq, guide):
    """Return the 30-nt RS3 context for one find_guides() dict, or None if out of bounds."""
    pam_len = len(guide['pam_seq'])
    # Slice in forward coords, then flip minus-strand guides into guide orientation.
    if guide['strand'] == '+':
        start = guide['guide_fwd_start'] - _UPSTREAM
        end = guide['pam_fwd_start'] + pam_len + _DOWNSTREAM
    else:
        start = guide['pam_fwd_start'] - _DOWNSTREAM
        end = guide['guide_fwd_end'] + _UPSTREAM
    # Never pad a truncated context — a guide too close to the region edge gets no score
    if start < 0 or end > len(seq):
        return None
    # Genomic FASTA may be soft-masked, so uppercase before anything else
    context = seq[start:end].upper()
    # Ambiguity codes (N and friends) make the context unscoreable, and would also
    # blow up reverse_complement, so drop those guides rather than guessing a base
    if set(context) - set('ACGT'):
        return None
    if guide['strand'] == '-':
        context = reverse_complement(context)
    return context if len(context) == RS3_CONTEXT_LEN else None


def guide_key(guide):
    """Stable identity for a guide, matching the guide_id the reagents UI synthesizes."""
    return '{}_{}_{}'.format(guide['strand'], guide['pam_fwd_start'], guide['spacer'])


def load_rs3():
    """Import rs3's predict_seq, returning (predict_fn, None) or (None, reason_string)."""
    try:
        from rs3.seq import predict_seq
    except ImportError as e:
        return None, 'RS3 unavailable ({}); guides will not be scored'.format(e)
    except Exception as e:
        # lightgbm/libomp dlopen and sklearn unpickle errors surface here, not as ImportError
        return None, 'RS3 unavailable ({}: {}); guides will not be scored'.format(
            type(e).__name__, e)
    return predict_seq, None


def score_guides(seq, guides, pam='NGG', guide_length=20, tracr='Hsu2013'):
    """Score every guide with RS3 in one batched call; returns (scores_by_key, note)."""
    if not rs3_supported(pam, guide_length):
        return {}, 'RS3 skipped: requires an NGG PAM and 20 nt spacers (got {}, {} nt)'.format(
            pam, guide_length)
    if not guides:
        return {}, None
    if tracr not in RS3_TRACRS:
        tracr = 'Hsu2013'

    # Build contexts first so a missing rs3 install costs nothing beyond the import attempt
    contexts_by_key = {}
    for g in guides:
        context = guide_context_30mer(seq, g)
        if context:
            contexts_by_key[guide_key(g)] = context
    if not contexts_by_key:
        return {}, 'RS3 skipped: no guide had a full 30 nt context'

    predict_seq, reason = load_rs3()
    if predict_seq is None:
        return {}, reason

    # One predict_seq call for the whole region: dedupe contexts, score, then map back
    unique = sorted(set(contexts_by_key.values()))
    try:
        predicted = predict_seq(unique, sequence_tracr=tracr)
    except Exception as e:
        return {}, 'RS3 scoring failed ({}: {}); guides will not be scored'.format(
            type(e).__name__, e)

    by_context = {ctx: round(float(score), 3) for ctx, score in zip(unique, predicted)}
    scores_by_key = {key: by_context[ctx] for key, ctx in contexts_by_key.items()}
    skipped = len(guides) - len(scores_by_key)
    note = None
    if skipped:
        note = '{} guide(s) near the region edge had no 30 nt context and were not scored'.format(
            skipped)
    return scores_by_key, note


def rs3_band(score):
    """Map an RS3 score to 'low'/'medium'/'high' using fixed absolute cutoffs."""
    if score is None:
        return None
    # Fixed (not percentile) cutoffs so a band means the same thing across different genes;
    # +/-0.5 is roughly half a standard deviation of the z-scored training data
    if score < -0.5:
        return 'low'
    if score <= 0.5:
        return 'medium'
    return 'high'


def rs3_percentiles(scores_by_key):
    """Rank-percentile (0-100) of each guide's score among all scored guides in the region."""
    if not scores_by_key:
        return {}
    values = sorted(scores_by_key.values())
    n = len(values)
    # Midrank percentile: guides sharing a score share one percentile instead of straddling
    # the jump between them
    percentiles = {}
    for key, score in scores_by_key.items():
        lower = bisect_left(values, score)
        ties = bisect_right(values, score) - lower
        percentiles[key] = round(100.0 * (lower + 0.5 * ties) / n, 1)
    return percentiles
