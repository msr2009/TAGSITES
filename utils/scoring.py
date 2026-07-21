"""utils/scoring.py — tag-site suggestion scoring (see GitHub issue #32)

Combines the per-residue analysis tracks already loaded for the Results tab
into a simple additive score: each of 8 criteria contributes +1 when
satisfied, giving an integer score in [0, 8] per position. A criterion whose
underlying data isn't available (e.g. no Reagents run) contributes 0 rather
than penalizing the position, so the score degrades gracefully to [0, 6].
"""

import pandas as pd

from utils.results import _guess_analysis_type

# amino acids treated as "small" for the flanking-residue criterion
_SMALL_AA = set("GSA")

# criteria are named for the boolean columns returned in the score DataFrame
_CRITERIA = [
    "not_in_domain",
    "not_in_modification",
    "not_in_hydrophobic_patch",
    "low_plddt",
    "low_conservation",
    "near_guide",
    "not_near_splice",
    "small_aa_flanked",
]

# short label for each criterion's *failed* state — shown as a tooltip flag
# when a position didn't earn that point
CRITERION_LABELS = {
    "not_in_domain": "in domain",
    "not_in_modification": "in modification",
    "not_in_hydrophobic_patch": "in hydrophobic patch",
    "low_plddt": "pLDDT>0.5",
    "low_conservation": "SJD>0.5",
    "near_guide": "no guide<5bp",
    "not_near_splice": "near splice site",
    "small_aa_flanked": "few small flanking",
}


def failed_flags(scores_df):
    """For each position, list the CRITERION_LABELS of criteria that were NOT satisfied.

    These are the "reasons points were lost" shown in the Results-tab tooltip,
    e.g. H28: 5 [pLDDT>0.5, SJD>0.5, few small flanking].
    """
    flags = []
    for _, row in scores_df.iterrows():
        flags.append([CRITERION_LABELS[c] for c in _CRITERIA if not row[c]])
    return flags


def score_green_hex(frac):
    """Interpolate white (#ffffff) to green (#28a745) at frac ∈ [0, 1]."""
    frac = max(0.0, min(1.0, frac))
    r0, g0, b0 = 0xff, 0xff, 0xff
    r1, g1, b1 = 0x28, 0xa7, 0x45
    r = int(r0 + frac * (r1 - r0))
    g = int(g0 + frac * (g1 - g0))
    b = int(b0 + frac * (b1 - b0))
    return f"#{r:02x}{g:02x}{b:02x}"


# Phobius topology labels (see utils.results._PHOBIUS_KEYWORDS) that describe bulk
# membrane orientation rather than a structurally distinct domain — these shouldn't
# count against a position for the "not_in_domain" criterion
_PHOBIUS_NON_DOMAIN_LABELS = {"Cytoplasmic", "Extracellular"}


def _range_mask(range_df, sources, seq_len, exclude_descriptions=None):
    """Return a set of 1-based positions covered by any interval whose source is in `sources`.

    Rows whose description is in `exclude_descriptions` are skipped even if their
    source matches.
    """
    positions = set()
    if range_df is None or range_df.empty:
        return positions
    subset = range_df[range_df["source"].isin(sources)]
    if exclude_descriptions:
        subset = subset[~subset["description"].isin(exclude_descriptions)]
    for _, row in subset.iterrows():
        start, stop = int(row["start"]), int(row["stop"])
        positions.update(range(max(1, start), min(stop, seq_len) + 1))
    return positions


def _find_column(aa_df, analysis_type):
    """Return the first aa_df column classified as `analysis_type` by _guess_analysis_type, or None."""
    if aa_df is None:
        return None
    for col in aa_df.columns[1:]:
        if _guess_analysis_type(col) == analysis_type:
            return col
    return None


def _small_aa_fraction(query_seq, pos, window, small_aa):
    """Fraction of residues in query_seq[pos-1-window : pos+window] (1-based pos) that are small AAs."""
    if not query_seq:
        return float("nan")
    lo = max(0, pos - 1 - window)
    hi = min(len(query_seq), pos + window)
    flank = query_seq[lo:hi]
    if not flank:
        return float("nan")
    return sum(1 for aa in flank if aa in small_aa) / len(flank)


def _reagents_min_distance(reagents_df, pos, column):
    """Minimum value of `column` across all reagents-TSV rows for this residue, or NaN if absent."""
    if reagents_df is None or column not in reagents_df.columns:
        return float("nan")
    rows = reagents_df[reagents_df["residue_index"] == pos]
    if rows.empty:
        return float("nan")
    vals = rows[column]
    vals = vals[vals >= 0]  # drop sentinel negative values (e.g. dist_to_3p_splice == -1)
    if vals.empty:
        return float("nan")
    return vals.min()


def score_tag_sites(aa_df, range_df, query_seq, reagents_df=None,
                     plddt_thresh=0.5, cons_thresh=0.5, guide_bp=5,
                     splice_bp=5, small_aa=None, small_frac=0.70, window=5):
    """Score every residue position for suitability as a tag insertion site.

    Each of 8 criteria (see issue #32) adds +1 to the position's score when
    satisfied; missing data contributes 0 rather than penalizing the position.

    Returns a DataFrame indexed by 1-based `pos` with one bool column per
    criterion plus an integer "score" column (sum of the criteria columns).
    """
    small_aa = set(small_aa) if small_aa else _SMALL_AA
    seq_len = len(query_seq) if query_seq else 0
    if aa_df is not None and not aa_df.empty:
        seq_len = max(seq_len, int(aa_df["pos"].max()))

    positions = list(range(1, seq_len + 1))
    out = pd.DataFrame(index=pd.Index(positions, name="pos"))

    domain_mask = _range_mask(range_df, {"Pfam", "Phobius"}, seq_len,
                               exclude_descriptions=_PHOBIUS_NON_DOMAIN_LABELS)
    mod_mask = _range_mask(range_df, {"modification"}, seq_len)
    patch_mask = _range_mask(range_df, {"hydrophobic_patch"}, seq_len)
    out["not_in_domain"] = [p not in domain_mask for p in positions]
    out["not_in_modification"] = [p not in mod_mask for p in positions]
    out["not_in_hydrophobic_patch"] = [p not in patch_mask for p in positions]

    plddt_col = _find_column(aa_df, "plddt")
    cons_col = _find_column(aa_df, "blast")
    if plddt_col is not None:
        plddt = aa_df.set_index("pos")[plddt_col]
        out["low_plddt"] = [plddt.get(p, float("nan")) < plddt_thresh for p in positions]
    else:
        out["low_plddt"] = False
    if cons_col is not None:
        cons = aa_df.set_index("pos")[cons_col]
        out["low_conservation"] = [cons.get(p, float("nan")) < cons_thresh for p in positions]
    else:
        out["low_conservation"] = False

    if reagents_df is not None:
        out["near_guide"] = [
            _reagents_min_distance(reagents_df, p, "distance") < guide_bp for p in positions
        ]
        # nearest of the two splice distances; NaN (unknown) is excluded rather than
        # compared directly, since min() with NaN operands is order-dependent
        out["not_near_splice"] = [
            min([d for d in (
                _reagents_min_distance(reagents_df, p, "dist_to_5p_splice"),
                _reagents_min_distance(reagents_df, p, "dist_to_3p_splice"),
            ) if pd.notna(d)], default=float("nan")) >= splice_bp
            for p in positions
        ]
    else:
        out["near_guide"] = False
        out["not_near_splice"] = False

    out["small_aa_flanked"] = [
        _small_aa_fraction(query_seq, p, window, small_aa) > small_frac for p in positions
    ]

    # NaN comparisons above evaluate False, which is the desired "no penalty" behavior
    out["score"] = out[_CRITERIA].sum(axis=1).astype(int)
    return out
