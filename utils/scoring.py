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

Keys starting with "user_" are reserved for runtime mask criteria synthesized
from user choices (isoform/topology restrictions — see utils/tag_filters.py and
position_mask_criterion() below) rather than scores.config.json; don't use that
prefix for a config-authored criterion key.
"""

import json
import os
from pathlib import Path

import pandas as pd

from utils.results import _guess_analysis_type, _interp_cmap, _VIRIDIS

# default location of the editable scoring config, next to task_definitions.json
SCORING_CONFIG_PATH = Path(__file__).parent.parent / "scores.config.json"

# Phobius topology labels (see utils.results._PHOBIUS_KEYWORDS) that describe bulk
# membrane orientation rather than a structurally distinct domain — handled via
# the "exclude_descriptions" param on a range_absent criterion


_CONFIG_CACHE = {}   # path -> (mtime, parsed config)

# a graded criterion below this contributes too little to count as satisfied, for the
# pass/fail views that cannot show a fraction (tooltips, the database's criterion bitmask)
PARTIAL_CREDIT_FLOOR = 0.5


def satisfaction(value):
    """A criterion's contribution as a 0-1 fraction of its weight.

    True/False and missing data (None/NaN) behave exactly as before — missing never
    penalizes. A graded handler returns a float, which passes through clipped.
    """
    if value is None:
        return 0.0
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        v = float(value)
        return 0.0 if v != v else min(max(v, 0.0), 1.0)
    return 1.0 if value else 0.0


def load_scoring_config(path=None):
    """Load scores.config.json (or a given path), returning its parsed criteria list.

    Parsed once per file version (keyed on mtime, so edits still take effect without a
    restart); each caller gets its own deep copy, since callers extend the criteria.
    """
    import copy
    path = str(path or SCORING_CONFIG_PATH)
    mtime = os.path.getmtime(path)
    hit = _CONFIG_CACHE.get(path)
    if hit is None or hit[0] != mtime:
        with open(path) as f:
            hit = (mtime, json.load(f))
        _CONFIG_CACHE[path] = hit
    return copy.deepcopy(hit[1])


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
        reasons = []
        for c in active:
            got = satisfaction(row[c["key"]])
            # a graded criterion is a reason only when it earned less than half its weight,
            # and it says how much it did earn rather than reading as a flat failure
            if got < PARTIAL_CREDIT_FLOOR:
                pts = c["weight"] * got
                reasons.append(c["label"] if pts == 0 else
                               "{} ({:.1f} pt)".format(c["label"], pts))
        flags.append(reasons)
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


def position_mask_criterion(masked_positions, key, label):
    """Build a runtime mask criterion (weight 0) blocking an explicit set of positions.

    Used to fold a user-driven restriction (see utils/tag_filters.py) into the same
    config-criteria pipeline as scores.config.json's mask criteria, so masked/score/
    mask_reasons all pick it up with no other code changes.
    """
    return {"key": key, "label": label, "weight": 0, "mask": True,
            "type": "position_masked", "params": {"positions": set(masked_positions)}}


def with_extra_criteria(config, extra):
    """Shallow copy of a scoring config with extra runtime criteria appended."""
    if not extra:
        return config
    out = dict(config)
    out["criteria"] = list(config["criteria"]) + list(extra)
    return out


# topology labels that scores.config.json's hard "in_membrane_or_cleaved" mask also covers;
# each group lifts together because the config matches descriptions as case-insensitive
# substrings ("signal" also matches "Signal peptide")
_HARD_MASK_TOPOLOGY_GROUPS = [
    {"Signal peptide", "signal"},
    {"Transmembrane", "TMhelix"},
]


def unmask_checked_topology(config, checked_labels):
    """Copy of config whose hard signal/transmembrane mask skips groups the user has checked."""
    lifted = set()
    for group in _HARD_MASK_TOPOLOGY_GROUPS:
        if group & set(checked_labels):
            lifted |= group
    return lift_hard_masks(config, lifted)


def lift_hard_masks(config, descriptions):
    """Copy of config whose hard "in_membrane_or_cleaved" mask no longer matches descriptions."""
    if not descriptions:
        return config
    out = dict(config)
    out["criteria"] = []
    for c in config["criteria"]:
        if c.get("key") == "in_membrane_or_cleaved":
            params = dict(c["params"])
            params["include_descriptions"] = [d for d in params["include_descriptions"]
                                              if d not in descriptions]
            c = {**c, "params": params}
        out["criteria"].append(c)
    return out


def range_positions(range_df, sources, seq_len):
    """1-based positions inside any annotation interval from the given sources."""
    return _range_mask(range_df, set(sources), seq_len)


# fallbacks for a config written before the site_selection block existed
_SITE_SELECTION_DEFAULTS = {
    "max_sites": 5,
    "min_spacing": 10,
    "technical_criteria": ["near_guide", "not_near_splice"],
    "cut_near_bp": 5,
    "cut_far_bp": 20,
    "weight_proximity": 1.0,
    "weight_margin": 0.5,
}

_SUGGESTED_MAX_SITES = _SITE_SELECTION_DEFAULTS["max_sites"]
_SUGGESTED_MIN_SPACING = _SITE_SELECTION_DEFAULTS["min_spacing"]


def site_selection_params(config=None):
    """The site_selection block of the scoring config, with defaults for any missing key."""
    config = config or load_scoring_config()
    out = dict(_SITE_SELECTION_DEFAULTS)
    out.update({k: v for k, v in (config.get("site_selection") or {}).items()
                if k != "__doc__"})
    return out


def _criteria_signature(scores, config):
    """Per-position tuple of which biological criteria are satisfied.

    Reagent-availability criteria are excluded, so a guide appearing or disappearing never
    splits a region. This reproduces the database's criterion bitmask from the scores frame,
    which the app has no bitmask for: a criterion counts as satisfied at PARTIAL_CREDIT_FLOOR,
    exactly as scripts/build_proteome_db.py packs its bits.
    """
    tech = set(site_selection_params(config)["technical_criteria"])
    keys = [c["key"] for c in config["criteria"]
            if not c.get("mask") and c["key"] not in tech and c["key"] in scores.columns]
    if not keys:
        return pd.Series([()] * len(scores), index=scores.index)
    flags = pd.concat(
        [scores[k].apply(lambda v: satisfaction(v) >= PARTIAL_CREDIT_FLOOR) for k in keys],
        axis=1)
    return pd.Series(list(map(tuple, flags.to_numpy())), index=scores.index)


def find_regions(scores, config=None):
    """Contiguous runs of unmasked residues tied at the best score with identical biology.

    Every residue in a region is an equivalent place to tag, so the choice among them is
    technical rather than biological (see pick_in_region). Masks and any change of satisfied
    biological criteria split a region. Returns dicts with start/stop/width/mid, widest first.
    """
    config = config or load_scoring_config()
    if scores is None or scores.empty or "score" not in scores.columns:
        return []
    masked = (scores["masked"].astype(bool) if "masked" in scores.columns
              else pd.Series(False, index=scores.index))
    free = scores[(scores["score"] > 0) & (~masked)]
    if free.empty:
        return []
    at_best = free[free["score"] >= free["score"].max()]
    sig = _criteria_signature(scores, config)
    pos = list(at_best.index)
    regions, start = [], pos[0]
    # a gap of more than one residue, or a change of criteria, ends the current region
    for a, b in zip(pos, pos[1:] + [None]):
        if b is None or b - a > 1 or sig.loc[b] != sig.loc[a]:
            regions.append({"start": int(start), "stop": int(a), "width": int(a - start + 1),
                            "mid": int((start + a) // 2)})
            start = b
    return sorted(regions, key=lambda g: (-g["width"], g["start"]))


def pick_in_region(region, scores, config=None):
    """The residue in a region closest to a guide cut and furthest from the region's edges.

    Every residue here scores the same biologically, so this is purely a technical tiebreak.
    With no guide anywhere in the region the proximity term is 0 throughout and the rule
    degenerates to the region's centre, which is the right default; the region is still offered.
    """
    p = site_selection_params(config)
    start, stop = region["start"], region["stop"]
    centre = (start + stop) // 2
    if "guide_distance" not in scores.columns:
        return centre
    sub = scores.loc[start:stop]
    dist = pd.to_numeric(sub["guide_distance"], errors="coerce")
    if dist.isna().all():
        return centre
    prox = dist.apply(lambda d: _ramp(d, p["cut_near_bp"], p["cut_far_bp"])).fillna(0.0)
    # depth into the region: 1 at the centre, 0 at either edge, and 1 for a width-1 region
    half = max((stop - start) // 2, 1)
    idx = sub.index.to_series()
    margin = (idx - start).combine(stop - idx, min) / half
    total = p["weight_proximity"] * prox + p["weight_margin"] * margin
    return int(total.idxmax())


def pick_suggested_regions(scores, max_sites=None, min_spacing=None, config=None):
    """The offered regions, each as start/stop plus the residue picked inside it, by position.

    The count is data-driven: a protein whose best score forms two regions offers two,
    never padded to max_sites with lower-scoring residues. Within each region the residue is
    chosen on technical grounds only (see pick_in_region). The Results plot draws the spans;
    everything else wants only the residue and calls pick_suggested_sites.
    """
    config = config or load_scoring_config()
    p = site_selection_params(config)
    max_sites = p["max_sites"] if max_sites is None else max_sites
    min_spacing = p["min_spacing"] if min_spacing is None else min_spacing
    picked = []
    for region in find_regions(scores, config):
        if len(picked) >= max_sites:
            break
        pos = pick_in_region(region, scores, config)
        # two regions can be adjacent when they satisfy different criteria at the same score;
        # offering both would spend the cap on what reads as one location, so skip to the
        # next region instead (the cap still fills if the protein has somewhere else)
        if all(abs(pos - q["pos"]) >= min_spacing for q in picked):
            picked.append({"start": int(region["start"]), "stop": int(region["stop"]),
                           "pos": int(pos)})
    return sorted(picked, key=lambda q: q["pos"])


def pick_suggested_sites(scores, max_sites=None, min_spacing=None, config=None):
    """One suggested site per equivalent-scoring region, as 1-based positions.

    Used by the Results tab (issue #32), scripts/suggest_sites.py and the proteome
    database's scores stage.
    """
    return [q["pos"] for q in pick_suggested_regions(scores, max_sites=max_sites,
                                                     min_spacing=min_spacing, config=config)]


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


def _nearest_guide_by_residue(reagents_df):
    """Per-residue cut-to-insert distance of the nearest guide, as a Series indexed by residue.

    Negative distances are sentinels and are dropped, matching _reagents_min_distance.
    """
    if reagents_df is None or "distance" not in reagents_df.columns:
        return None
    rows = reagents_df[reagents_df["distance"] >= 0]
    if rows.empty:
        return None
    nearest = rows.groupby("residue_index")["distance"].min()
    nearest.index.name = "pos"
    return nearest.rename("guide_distance")


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


def _ramp(value, low, high):
    """0-1 response: full credit at or past `low`, none at or past `high`, linear between."""
    if value != value:  # NaN: missing data, handled as not-satisfied without penalty
        return float("nan")
    if high == low:
        return 1.0 if value < high else 0.0
    frac = (high - value) / (high - low)
    return min(max(frac, 0.0), 1.0)


def _handle_column_ramp(positions, params, ctx):
    """Graded version of column_below: `low` earns the full weight, `high` earns none."""
    col = _find_column(ctx["aa_df"], params["analysis"])
    if col is None:
        return [None] * len(positions)
    series = ctx["aa_df"].set_index("pos")[col]
    low, high = params["low"], params["high"]
    return [_ramp(series.get(p, float("nan")), low, high) for p in positions]


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


def _handle_position_masked(positions, params, ctx):
    """True where the position is in params['positions'] — a runtime user-supplied mask."""
    blocked = params["positions"]
    return [p in blocked for p in positions]


SCORE_HANDLERS = {
    "range_absent":         _handle_range_absent,
	"range_present":		_handle_range_present,
    "column_below":         _handle_column_below,
    "column_ramp":          _handle_column_ramp,
    "reagent_min_below":    _handle_reagent_min_below,
    "reagent_min_above":    _handle_reagent_min_above,
    "flank_small_fraction": _handle_flank_small_fraction,
    "position_masked":      _handle_position_masked,
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
        c["weight"] * out[c["key"]].apply(satisfaction) for c in score_criteria
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

    # raw nearest-guide distance for the within-region pick (pick_in_region). This is not a
    # criterion and never contributes to the score; it rides along so pick_suggested_sites
    # needs nothing but this frame.
    nearest = _nearest_guide_by_residue(reagents_df)
    if nearest is None:
        out["guide_distance"] = float("nan")
    else:
        out = out.join(nearest)

    return out


def build_user_criteria(run_json, isoform_mode="none", isoform_names=None, topology_labels=None):
    """Build extra runtime mask criteria for CLI isoform/topology restriction flags.

    isoform_mode: "none" (no restriction, the default — returns [] unchanged), "constitutive"
    (mask anything not present in every isoform), or "tagged" (mask anything not present in
    exactly the isoforms named by isoform_names, matched against accession or name). Used by
    scripts/rescore.py and scripts/suggest_sites.py to let large-scale batch runs design
    constitutive and isoform-unique sites without going through the Results-tab UI.
    """
    # local imports avoid a hard circular dependency at module load time
    from config import RESULTS_TYPE_DICT
    from utils.results import load_data_from_json, load_run_metadata
    from utils.tag_filters import isoform_key, isoform_allowed_positions, topology_allowed_positions

    if isoform_mode == "none" and not topology_labels:
        return []

    aa_df, range_df, _, iso_result = load_data_from_json(run_json, RESULTS_TYPE_DICT)
    meta = load_run_metadata(run_json)
    seq_len = meta.get("seq_len", 0)
    isoforms = (iso_result or {}).get("isoforms", [])
    all_positions = set(range(1, seq_len + 1))
    extra = []

    if isoform_mode != "none":
        keys = [isoform_key(iso, i) for i, iso in enumerate(isoforms)]
        visible_keys = set(keys)
        if isoform_mode == "constitutive":
            allowed = isoform_allowed_positions(isoforms, visible_keys, set(), True, seq_len)
            label = "not constitutive"
        else:  # "tagged"
            names = set(isoform_names or [])
            tagged_keys = {k for k in keys if k in names}
            allowed = isoform_allowed_positions(isoforms, visible_keys, tagged_keys, False, seq_len)
            label = "not in selected isoforms"
        if allowed is not None:
            extra.append(position_mask_criterion(all_positions - allowed,
                                                 "user_isoform_mask", label))

    if topology_labels:
        topo_allowed = topology_allowed_positions(range_df, topology_labels, seq_len)
        if topo_allowed is not None:
            extra.append(position_mask_criterion(all_positions - topo_allowed,
                                                 "user_topology_mask",
                                                 "outside " + ", ".join(topology_labels)))
    return extra


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

    aa_df, range_df, _, _ = load_data_from_json(run_json, RESULTS_TYPE_DICT)
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
