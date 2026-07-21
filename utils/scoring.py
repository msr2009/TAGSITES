"""utils/scoring.py — tag-site suggestion scoring (see GitHub issue #32)

Combines the per-residue analysis tracks already loaded for the Results tab
into a weighted score: each criterion in scores.config.json contributes its
`weight` when satisfied, giving a score in [0, sum(weights)] per position. A
criterion whose underlying data isn't available (e.g. no Reagents run)
contributes 0 rather than penalizing the position. All scoring knobs (weights,
thresholds, labels, the "small" amino-acid set, etc.) live in scores.config.json
— see that file's __doc__ for how to edit or extend it.
"""

import json
from pathlib import Path

import pandas as pd

from utils.results import _guess_analysis_type

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
    e.g. H28: 5 [pLDDT>0.5, SJD>0.5, few small flanking].
    """
    config = config or load_scoring_config()
    active = [c for c in config["criteria"] if c["weight"] > 0]
    flags = []
    for _, row in scores_df.iterrows():
        flags.append([c["label"] for c in active if not row[c["key"]]])
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


# ── criterion handlers ───────────────────────────────────────────────────────
# each handler takes (positions, params, ctx) and returns a list[bool] aligned
# with `positions`; ctx bundles the data already loaded for scoring. Adding a
# criterion that reuses one of these types is config-only (scores.config.json);
# a genuinely new kind of test needs a new handler registered in SCORE_HANDLERS.

def _handle_range_absent(positions, params, ctx):
    """True where the position is NOT covered by any range interval from the given sources."""
    mask = _range_mask(ctx["range_df"], set(params["sources"]), ctx["seq_len"],
                        exclude_descriptions=params.get("exclude_descriptions"))
    return [p not in mask for p in positions]


def _handle_column_below(positions, params, ctx):
    """True where the aa_df column for `analysis` is below `threshold` (False if column missing)."""
    col = _find_column(ctx["aa_df"], params["analysis"])
    if col is None:
        return [False] * len(positions)
    series = ctx["aa_df"].set_index("pos")[col]
    threshold = params["threshold"]
    return [series.get(p, float("nan")) < threshold for p in positions]


def _handle_reagent_min_below(positions, params, ctx):
    """True where the minimum reagents-TSV `column` value for this position is below `threshold`."""
    if ctx["reagents_df"] is None:
        return [False] * len(positions)
    column, threshold = params["column"], params["threshold"]
    return [_reagents_min_distance(ctx["reagents_df"], p, column) < threshold for p in positions]


def _handle_reagent_min_above(positions, params, ctx):
    """True where the nearest of several reagents-TSV distance columns is at least `threshold`."""
    if ctx["reagents_df"] is None:
        return [False] * len(positions)
    columns, threshold = params["columns"], params["threshold"]
    out = []
    for p in positions:
        # nearest of the distance columns; NaN (unknown) is excluded rather than
        # compared directly, since min() with NaN operands is order-dependent
        dists = [d for d in (_reagents_min_distance(ctx["reagents_df"], p, c) for c in columns)
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

    Returns a DataFrame indexed by 1-based `pos` with one bool column per
    criterion plus a numeric "score" column (weighted sum of the criteria columns).
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

    # NaN comparisons inside handlers evaluate False, which is the desired
    # "no penalty" behavior
    out["score"] = sum(c["weight"] * out[c["key"]] for c in criteria)
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
