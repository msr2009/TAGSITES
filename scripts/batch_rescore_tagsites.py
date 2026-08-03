"""
batch_rescore_tagsites.py

CLI batch driver for iterating on tag-site SCORING SCHEMES: rescore and
replot a whole manifest of runs against one scores.config.json variant,
without re-running the analysis pipeline (everything is read from each run's
already-written task outputs, same as scripts/rescore.py does for one run).

Every position of every run is scored (not just the tagged site), so the
tagged site can be judged against the rest of its own protein. For each run
this writes a full per-position scores TSV and a per-protein annotated plot
(scripts/plot_tagsite_figure.py); across the whole batch it writes one
summary TSV (one row per allele, with each criterion's weighted contribution
at the tagged site plus tagged-vs-background stats) and one histogram
overlaying tagged-site scores on the pooled all-positions distribution.

Matt Rich, 2026
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_SCRIPTS = Path(__file__).parent
_ROOT = _SCRIPTS.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_SCRIPTS))

from config import RESULTS_TYPE_DICT
from utils.results import load_data_from_json, load_run_metadata, load_reagents_df
from utils.scoring import score_tag_sites, load_scoring_config

import plot_tagsite_figure
from batch_run_tag_sites import read_manifest


def _histogram_bins(max_score, n_default=20):
    """Bin edges spanning [0, max_score].

    scores.config.json weights are integers today (one bin per point, like
    batch_plot_tagsites.py's plot_score_histogram), but a scheme variant could
    use fractional weights — range(max_score + 2) would then break, so fall
    back to a fixed bin count whenever max_score isn't a (small) whole number.
    """
    if max_score <= 0:
        return [0, 1]
    if abs(max_score - round(max_score)) < 1e-9 and max_score <= 50:
        return [x - 0.5 for x in range(int(round(max_score)) + 2)]
    return list(np.linspace(0, max_score, n_default + 1))


def plot_tagged_vs_background(tagged_scores, all_scores, max_score, outfile):
    """Save a density-normalized histogram overlaying tagged-site scores on all-positions scores.

    all_scores pools every scored position across every run in the batch, so a
    scoring scheme that fails to separate real tag sites from background is
    visible as substantial overlap between the two distributions.
    """
    bins = _histogram_bins(max_score)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(all_scores, bins=bins, density=True, alpha=0.5, color="#999999",
            edgecolor="black", label=f"all positions (n={len(all_scores)})")
    ax.hist(tagged_scores, bins=bins, density=True, alpha=0.6, color="#c44e52",
            edgecolor="black", label=f"tagged sites (n={len(tagged_scores)})")
    ax.set_xlim(0, max_score)
    ax.set_xlabel(f"Tag-site score (0-{max_score})")
    ax.set_ylabel("Density")
    ax.set_title("Tagged-site scores vs. all-positions background")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    return outfile


def score_one_run(run_json_path, config):
    """Load one run's analysis tracks and score every position.

    Returns (scores_df, run_name) where scores_df is indexed by 1-based pos
    with one bool/None column per criterion plus "score" (see
    utils.scoring.score_tag_sites).
    """
    aa_df, range_df, _alns, _iso = load_data_from_json(run_json_path, RESULTS_TYPE_DICT)
    meta = load_run_metadata(run_json_path)
    reagents_df = load_reagents_df(run_json_path)
    scores_df = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""), reagents_df, config)
    run_name = Path(run_json_path).name.removesuffix(".run.json")
    return scores_df, run_name


def write_full_scores(scores_df, site, outfile):
    """Write one run's full per-position scores TSV, with a "tagged" column marking `site`."""
    df = scores_df.copy()
    df["tagged"] = False
    if site in df.index:
        df.loc[site, "tagged"] = True
    df.to_csv(outfile, sep="\t")


def summary_row(scores_df, site, gene, allele, run_name, config):
    """Build one summary-TSV row for a run: per-criterion contribution at `site`, plus totals.

    Criterion columns hold `weight` if the criterion passed at `site`, `0` if
    it failed, and "" (NA) if the handler returned None (data unavailable) —
    the "score" column alone can't tell missing-data apart from failure.
    """
    max_score = sum(c["weight"] for c in config["criteria"])
    row = {"gene": gene, "allele": allele, "site": site, "run_name": run_name}

    if site not in scores_df.index:
        print(f"WARNING: {run_name}: site {site} is outside the scored range "
              f"(1-{len(scores_df)}); recording as unscored", file=sys.stderr)
        for c in config["criteria"]:
            row[c["key"]] = ""
        row.update(score=float("nan"), max_score=max_score, n_positions=len(scores_df),
                   median_all=scores_df["score"].median() if len(scores_df) else float("nan"),
                   max_observed=scores_df["score"].max() if len(scores_df) else float("nan"),
                   pct_below="")
        return row, float("nan")

    site_row = scores_df.loc[site]
    for c in config["criteria"]:
        v = site_row[c["key"]]
        row[c["key"]] = "" if v is None else (c["weight"] if bool(v) else 0)

    site_score = site_row["score"]
    all_scores = scores_df["score"]
    row.update(
        score=site_score,
        max_score=max_score,
        n_positions=len(scores_df),
        median_all=all_scores.median(),
        max_observed=all_scores.max(),
        pct_below=(all_scores < site_score).mean(),
    )
    return row, site_score


def main(manifest_path, output_dir, prefix, config_path=None, make_plots=True):
    """Rescore and (optionally) replot every run in manifest_path under one scoring config."""
    rows = read_manifest(manifest_path)
    if not rows:
        print("Manifest is empty; nothing to do.", file=sys.stderr)
        return

    config = load_scoring_config(config_path)
    max_score = sum(c["weight"] for c in config["criteria"])
    criterion_keys = [c["key"] for c in config["criteria"]]

    output_dir = Path(output_dir)
    scores_dir = output_dir / "scores"
    plots_dir = output_dir / "plots"
    scores_dir.mkdir(parents=True, exist_ok=True)
    if make_plots:
        plots_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    tagged_scores = []
    all_scores_pooled = []

    for row in rows:
        run_json = row["run_json"]
        site = int(row["site"])
        gene = row.get("gene", "")
        allele = row.get("allele", "")
        approximate = row.get("approximate", "").strip().lower() in ("1", "true", "yes")

        try:
            scores_df, run_name = score_one_run(run_json, config)
        except Exception as e:
            print(f"WARNING: skipping {run_json}: {e}", file=sys.stderr)
            continue

        if not gene and not allele:
            # fall back to parsing INTERNAL_<gene>_<allele>_<token> from run_name
            parts = run_name.split("_")
            if len(parts) == 4 and parts[0] == "INTERNAL":
                gene, allele = parts[1], parts[2]

        write_full_scores(scores_df, site, scores_dir / f"{run_name}.scores.tsv")

        row_summary, site_score = summary_row(scores_df, site, gene, allele, run_name, config)
        summary_rows.append(row_summary)

        if site_score == site_score:  # skip NaN (out-of-range site) from the histogram
            tagged_scores.append(site_score)
        all_scores_pooled.extend(scores_df["score"].tolist())

        if make_plots:
            outfile = str(plots_dir / f"{run_name}.tagsite.png")
            try:
                plot_tagsite_figure.main(run_json, site, outfile, approximate=approximate, config=config)
            except Exception as e:
                print(f"WARNING: {run_name}: plot failed: {e}", file=sys.stderr)

    if not summary_rows:
        print("No runs scored successfully; nothing to write.", file=sys.stderr)
        return

    summary_path = output_dir / f"{prefix}.summary.tsv"
    fieldnames = ["gene", "allele", "site", "run_name", *criterion_keys,
                  "score", "max_score", "n_positions", "median_all", "max_observed", "pct_below"]
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote summary: {summary_path} ({len(summary_rows)} run(s))")

    if tagged_scores:
        hist_path = plot_tagged_vs_background(
            tagged_scores, all_scores_pooled, max_score,
            str(output_dir / f"{prefix}.score_histogram.png"),
        )
        print(f"Wrote tag-site score histogram: {hist_path}")
    else:
        print("No in-range tagged-site scores collected; skipping histogram.", file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True,
        help="CSV/TSV file, one row per run (run_json, site, and optional gene, "
             "allele, approximate columns; gene/allele fall back to parsing run_name)")
    parser.add_argument("--output-dir", required=True,
        help="directory for per-run scores TSVs, plots, the batch summary TSV, and histogram")
    parser.add_argument("--prefix", required=True,
        help="filename prefix for the summary TSV and histogram (e.g. <prefix>.summary.tsv)")
    parser.add_argument("--config", default=None,
        help="path to a scores.config.json variant (default: repo scores.config.json)")
    parser.add_argument("--no-plots", action="store_true", default=False,
        help="skip per-run plots; only write the scores TSVs, summary, and histogram")
    args = parser.parse_args()

    main(args.manifest, args.output_dir, args.prefix,
         config_path=args.config, make_plots=not args.no_plots)
