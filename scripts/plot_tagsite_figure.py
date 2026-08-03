"""
plot_tagsite_figure.py

Renders a run's full results plot (utils.results.plot_results — both the
continuous-score row and the categorical/range row) with a bold red vertical
line marking a known tag insertion site, annotated with the tag-site score at
that position, and exports it as a PNG. Used to visually check whether
pipeline signals (conservation, pLDDT, domains) line up with real
CRISPR-validated internal tag sites.

Parameters:
    - run_json (str): path to a completed run's <run_name>.run.json
    - site (int): residue position (1-based) to mark
    - outfile (str): PNG output path
    - approximate (bool): if True, label the site as approximate in the title

Returns:
    - (outfile, score, max_score) where score is the tag-site score at `site`
      (NaN if the position has no computed score) and max_score is the
      highest score achievable under the active scoring config

Matt Rich, 7/2026
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import RESULTS_TYPE_DICT
from utils.results import load_data_from_json, load_run_metadata, load_reagents_df, plot_results
from utils.scoring import score_tag_sites, load_scoring_config


def score_at_site(run_json, aa_df, range_df, site, config):
    """Compute the tag-site score and per-criterion pass/fail at `site` (NaN score if out of range).

    Returns (score, criteria) where criteria is a list of (label, passed) pairs,
    one per criterion in `config`, in config order. `passed` is True, False, or
    None if the underlying data wasn't available for that criterion (e.g. no
    Reagents run) — not the same as a criterion actually failing.
    """
    meta = load_run_metadata(run_json)
    reagents_df = load_reagents_df(run_json)
    scores_df = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""), reagents_df, config)
    if site not in scores_df.index:
        return float("nan"), [(c["label"], None) for c in config["criteria"]]
    row = scores_df.loc[site]
    criteria = [(c["label"], row[c["key"]]) for c in config["criteria"]]
    return row["score"], criteria


def main(run_json, site, outfile, approximate=False, config=None):
    config = config or load_scoring_config()
    max_score = sum(c["weight"] for c in config["criteria"])

    aa_df, range_df, _alns, _iso = load_data_from_json(run_json, RESULTS_TYPE_DICT)
    site_score, criteria = score_at_site(run_json, aa_df, range_df, site, config)

    run_name = Path(run_json).name.removesuffix(".run.json")
    title = f"{run_name} — tag site {site}" + (" (approximate)" if approximate else "")

    # SASA and hydrophobicity companion tracks clutter this figure; drop them,
    # keeping the underlying aa_df (used for scoring) untouched
    plot_df = aa_df[[c for c in aa_df.columns if not c.endswith(("_sasa", "_hydro"))]]

    fig = plot_results(plot_df, range_df, title=title)
    fig.add_vline(x=site, line_color="red", line_width=3, row="all", col="all")
    score_label = "score n/a" if site_score != site_score else f"score {site_score:.2f}/{max_score}"
    fig.add_annotation(
        x=site, y=1, yref="y domain", yanchor="bottom", row=1, col=1,
        text=score_label, showarrow=False,
        font=dict(color="red", size=14),
    )

    # criteria breakdown, printed as a subtitle line (a fixed location that
    # doesn't collide with the legend, which grows with how many domains/
    # patches a given protein has). None means the underlying data wasn't
    # available (e.g. no Reagents run) — shown as "n/a", not a failed check.
    def _mark(passed):
        if passed is None:
            return "–"
        return "✓" if passed else "✗"

    breakdown_text = "    ".join(f"{_mark(passed)} {label}" for label, passed in criteria)
    full_title = f"{title}<br><span style='font-size:12px'>{breakdown_text}</span>"
    fig.update_layout(title=dict(text=full_title), margin=dict(t=70))

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    fig.write_image(outfile, width=1400, height=700, scale=2)
    return outfile, site_score, max_score


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--run-json", action="store", type=str, dest="RUN_JSON",
        help="path to a completed run's <run_name>.run.json", required=True)
    parser.add_argument("--site", action="store", type=int, dest="SITE",
        help="residue position (1-based) to mark with a red line", required=True)
    parser.add_argument("--output", action="store", type=str, dest="OUTFILE",
        help="PNG output path", required=True)
    parser.add_argument("--approximate", action="store_true", default=False,
        help="label the site as approximate in the plot title")
    parser.add_argument("--config", action="store", type=str, dest="CONFIG", default=None,
        help="path to a scores.config.json variant (default: repo scores.config.json)")

    args = parser.parse_args()

    main(args.RUN_JSON, args.SITE, args.OUTFILE, args.approximate,
         config=load_scoring_config(args.CONFIG) if args.CONFIG else None)
