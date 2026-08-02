"""utils/scoring.py — tag-site suggestion scoring (see GitHub issues #32, #38)

Combines the per-residue analysis tracks already loaded for the Results tab
into a weighted score: each criterion in scores.config.json contributes its
`weight` when satisfied, giving a score in [0, sum(weights)] per position. A
criterion whose underlying data isn't available (e.g. no Reagents run)
contributes 0 rather than penalizing the position. All scoring knobs (weights,
thresholds, labels, the "small" amino-acid set, etc.) live in scores.config.json
— see that file's __doc__ for how to edit or extend it.

A criterion may also set `"mask": true` (weight 0). Mask criteria never add to
the weighted sum; instead, if any of them fires, the position's score is forced
to 0 regardless of the other criteria — used for categorical disqualifiers
(transmembrane/signal/propeptide spans, high-confidence modifications) that
should never be tagged (issue #38).
"""

import json
from pathlib import Path

import pandas as pd

from utils.results import _guess_analysis_type, _interp_cmap, _VIRIDIS

# default location of the editable scoring config, next to task_definitions.json
SCORING_CONFIG_PATH = Path(__file__).parent.parent / "scores.config.json"

# Phobius topology labels (see utils.results._PHOBIUS_KEYWORDS) that describe bulk
# membrane orientation rather than a structurally distinct domain — handled via
# the "exclude_descriptions" param on a range_absent criterion


def load_scoring_config(path=None):
    """Load scores.config.json (or a given path), returning its parsed criteria list."""
    path = path or SCORING_CONFIG_PATH
    with open(path) as f:
        return json.load(f)


def failed_flags(scores_df, config=None):
    """For each position, list the labels of criteria that were NOT satisfied (weight > 0 only).

    These are the "reasons points were lost" shown in the Results-tab tooltip,
    e.g. H28: 5 [pLDDT>0.5, SJD>0.5, few small flanking]. Mask criteria have
    weight 0 so they're excluded here; see mask_reasons() for those.
    """
    config = config or load_scoring_config()
    active = [c for c in config["criteria"] if c["weight"] > 0]
    flags = []
    for _, row in scores_df.iterrows():
        flags.append([c["label"] for c in active if not row[c["key"]]])
    return flags


def mask_reasons(scores_df, config=None):
    """For each position, list the labels of mask criteria that fired (empty if unmasked)."""
    config = config or load_scoring_config()
    masks = [c for c in config["criteria"] if c.get("mask")]
    reasons = []
    for _, row in scores_df.iterrows():
        reasons.append([c["label"] for c in masks if row[c["key"]]])
    return reasons


def score_max(config=None):
    """Sum of weights over non-mask criteria — the normalization ceiling for the score track."""
    config = config or load_scoring_config()
    return sum(c["weight"] for c in config["criteria"] if not c.get("mask"))


def score_cmap_hex(frac):
    """Map frac ∈ [0, 1] through the viridis colormap (dark purple = low, yellow = high)."""
    return _interp_cmap(frac, _VIRIDIS)


def _range_mask(range_df, sources, seq_len, exclude_descriptions=None, include_descriptions=None):
    """Return a set of 1-based positions covered by any interval whose source is in `sources`.

    Rows whose description is in `exclude_descriptions` are skipped even if their
    source matches. If `include_descriptions` is given, only rows whose description
    contains (case-insensitive substring) one of those strings are kept — descriptions
    are often composite (e.g. "Transmembrane: Helical (PubMed:123)"), so substring
    matching is used rather than equality.
    """
    positions = set()
    if range_df is None or range_df.empty:
        return positions
    subset = range_df[range_df["source"].isin(sources)]
    if exclude_descriptions:
        subset = subset[~subset["description"].isin(exclude_descriptions)]
    if include_descriptions:
        needles = [d.lower() for d in include_descriptions]
        subset = subset[subset["description"].str.lower().apply(
            lambda desc: any(n in desc for n in needles))]
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


def _reagents_min_by_residue(reagents_df, column):
    """Precompute, once per (reagents_df, column), a residue_index -> min(column) Series.

    Scoring every position of a run (rather than just the tagged site) turns the
    naive per-position `reagents_df[reagents_df["residue_index"] == pos]` scan into
    an O(positions x rows) cost against a 6-24MB DataFrame — the dominant cost of a
    full-protein rescore. groupby-once + O(1) lookup preserves the exact prior
    semantics (drop negative sentinels before the min; NaN when a residue has no
    rows) while making it O(rows + positions).
    """
    if reagents_df is None or column not in reagents_df.columns:
        return None
    non_negative = reagents_df[reagents_df[column] >= 0]
    return non_negative.groupby("residue_index")[column].min()


def _reagents_min_distance(reagents_df, pos, column, _cache=None):
    """Minimum value of `column` across all reagents-TSV rows for this residue, or NaN if absent.

    `_cache` is an optional {column: Series} map from `_reagents_min_by_residue`,
    used by the handlers below to avoid rescanning reagents_df per position.
    """
    if reagents_df is None or column not in reagents_df.columns:
        return float("nan")
    if _cache is not None:
        series = _cache.get(column)
        return series.get(pos, float("nan")) if series is not None else float("nan")
    rows = reagents_df[reagents_df["residue_index"] == pos]
    if rows.empty:
        return float("nan")
    vals = rows[column]
    vals = vals[vals >= 0]  # drop sentinel negative values (e.g. dist_to_3p_splice == -1)
    if vals.empty:
        return float("nan")
    return vals.min()


# ── criterion handlers ───────────────────────────────────────────────────────
# each handler takes (positions, params, ctx) and returns a list[bool] aligned
# with `positions`; ctx bundles the data already loaded for scoring. Adding a
# criterion that reuses one of these types is config-only (scores.config.json);
# a genuinely new kind of test needs a new handler registered in SCORE_HANDLERS.

def _handle_range_absent(positions, params, ctx):
    """True where the position is NOT covered by any range interval from the given sources."""
    mask = _range_mask(ctx["range_df"], set(params["sources"]), ctx["seq_len"],
                        exclude_descriptions=params.get("exclude_descriptions"),
                        include_descriptions=params.get("include_descriptions"))
    return [p not in mask for p in positions]


def _handle_range_present(positions, params, ctx):
    """True where the position IS covered by any range interval from the given sources."""
    mask = _range_mask(ctx["range_df"], set(params["sources"]), ctx["seq_len"],
                        exclude_descriptions=params.get("exclude_descriptions"),
                        include_descriptions=params.get("include_descriptions"))
    return [p in mask for p in positions]


def _handle_column_below(positions, params, ctx):
    """True where the aa_df column for `analysis` is below `threshold` (None if column missing)."""
    col = _find_column(ctx["aa_df"], params["analysis"])
    if col is None:
        return [None] * len(positions)
    series = ctx["aa_df"].set_index("pos")[col]
    threshold = params["threshold"]
    return [series.get(p, float("nan")) < threshold for p in positions]


def _reagents_cache_for(ctx, column):
    """Get (building + memoizing on ctx if needed) the per-residue min Series for `column`."""
    cache = ctx.setdefault("_reagents_cache", {})
    if column not in cache:
        cache[column] = _reagents_min_by_residue(ctx["reagents_df"], column)
    return cache


def _handle_reagent_min_below(positions, params, ctx):
    """True where the minimum reagents-TSV `column` value for this position is below `threshold`."""
    if ctx["reagents_df"] is None:
        return [None] * len(positions)
    column, threshold = params["column"], params["threshold"]
    cache = _reagents_cache_for(ctx, column)
    return [_reagents_min_distance(ctx["reagents_df"], p, column, cache) < threshold for p in positions]


def _handle_reagent_min_above(positions, params, ctx):
    """True where the nearest of several reagents-TSV distance columns is at least `threshold`."""
    if ctx["reagents_df"] is None:
        return [None] * len(positions)
    columns, threshold = params["columns"], params["threshold"]
    for c in columns:
        _reagents_cache_for(ctx, c)
    cache = ctx["_reagents_cache"]
    out = []
    for p in positions:
        # nearest of the distance columns; NaN (unknown) is excluded rather than
        # compared directly, since min() with NaN operands is order-dependent
        dists = [d for d in (_reagents_min_distance(ctx["reagents_df"], p, c, cache) for c in columns)
                 if pd.notna(d)]
        out.append(min(dists, default=float("nan")) >= threshold)
    return out


def _handle_flank_small_fraction(positions, params, ctx):
    """True where more than `fraction` of the +/-`window` flank are small amino acids."""
    small_aa = set(params["small_aa"])
    window, fraction = params["window"], params["fraction"]
    return [
        _small_aa_fraction(ctx["query_seq"], p, window, small_aa) > fraction
        for p in positions
    ]


SCORE_HANDLERS = {
    "range_absent":         _handle_range_absent,
	"range_present":		_handle_range_present,
    "column_below":         _handle_column_below,
    "reagent_min_below":    _handle_reagent_min_below,
    "reagent_min_above":    _handle_reagent_min_above,
    "flank_small_fraction": _handle_flank_small_fraction,
}


def score_tag_sites(aa_df, range_df, query_seq, reagents_df=None, config=None):
    """Score every residue position for suitability as a tag insertion site.

    Each criterion in `config` (scores.config.json by default) adds its `weight`
    to the position's score when satisfied; missing data contributes 0 rather
    than penalizing the position.

    Returns a DataFrame indexed by 1-based `pos` with one column per criterion
    (True/False, or None where the underlying data wasn't available — e.g. no
    Reagents run) plus a numeric "score" column (weighted sum, None counting
    as 0/not-satisfied).
    """
    config = config or load_scoring_config()
    criteria = config["criteria"]

    seq_len = len(query_seq) if query_seq else 0
    if aa_df is not None and not aa_df.empty:
        seq_len = max(seq_len, int(aa_df["pos"].max()))

    positions = list(range(1, seq_len + 1))
    out = pd.DataFrame(index=pd.Index(positions, name="pos"))

    ctx = {
        "aa_df": aa_df, "range_df": range_df, "query_seq": query_seq,
        "reagents_df": reagents_df, "seq_len": seq_len,
    }

    for c in criteria:
        handler = SCORE_HANDLERS[c["type"]]
        out[c["key"]] = handler(positions, c["params"], ctx)

    mask_criteria = [c for c in criteria if c.get("mask")]
    score_criteria = [c for c in criteria if not c.get("mask")]

    # NaN comparisons inside handlers evaluate False, and None (missing data)
    # is treated as not-satisfied — both are the desired "no penalty" behavior
    out["score"] = sum(
        c["weight"] * out[c["key"]].apply(lambda v: bool(v)) for c in score_criteria
    )

    # a fired mask criterion (weight 0, never a penalty for missing data) forces
    # the position's score to 0 regardless of how many other criteria it satisfied
    if mask_criteria:
        out["masked"] = pd.concat(
            [out[c["key"]].apply(lambda v: bool(v)) for c in mask_criteria], axis=1
        ).any(axis=1)
    else:
        out["masked"] = False
    out.loc[out["masked"], "score"] = 0

    return out


def write_run_scores(run_json_path, config=None, output_path=None):
    """Compute tag-site scores for a completed run and write them to a TSV.

    Loads the run JSON's analysis tracks (aa_df/range_df/reagents_df), scores
    every position, and writes pos + one bool column per criterion + score to
    `output_path` (default: <working_dir>/<run_name>.scores.tsv next to the run
    JSON). Returns the output path.
    """
    # local imports avoid a hard circular dependency at module load time
    from config import RESULTS_TYPE_DICT
    from utils.results import load_data_from_json, load_run_metadata, load_reagents_df

    with open(run_json_path) as f:
        run_json = json.load(f)

    aa_df, range_df, _ = load_data_from_json(run_json, RESULTS_TYPE_DICT)
    meta = load_run_metadata(run_json)
    reagents_df = load_reagents_df(run_json)

    config = config or load_scoring_config()
    scores_df = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""), reagents_df, config)

    if output_path is None:
        g = run_json.get("global", {})
        working_dir = g.get("working_dir", "")
        run_name = g.get("run_name", "")
        output_path = str(Path(working_dir) / f"{run_name}.scores.tsv")

    scores_df.to_csv(output_path, sep="\t")
    return output_path
