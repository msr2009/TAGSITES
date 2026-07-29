"""
batch_plot_tagsites.py

CLI batch driver: renders a per-protein tag-site figure (plot_tagsite_figure)
for every row in a manifest, then collects the tag-site score from each run
into a matplotlib histogram summarizing the whole batch.

Matt Rich, 2026
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

_SCRIPTS = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS))

import plot_tagsite_figure
from batch_run_tag_sites import read_manifest


def plot_score_histogram(scores, max_score, outfile):
    """Save a matplotlib histogram of tag-site scores across the batch, binned over [0, max_score]."""
    bins = [x - 0.5 for x in range(max_score + 2)]  # one bin per integer score, 0..max_score
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scores, bins=bins, range=(0, max_score), color="#4c72b0", edgecolor="black")
    ax.set_xlim(0, max_score)
    ax.set_xlabel(f"Tag-site score (0-{max_score})")
    ax.set_ylabel("Count")
    ax.set_title(f"Tag-site score distribution (n={len(scores)})")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    return outfile


def main(manifest_path, output_dir, histogram_name="tagsite_score_histogram.png"):
    """Render one annotated figure per manifest row and a batch-wide score histogram."""
    rows = read_manifest(manifest_path)
    if not rows:
        print("Manifest is empty; nothing to do.", file=sys.stderr)
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scores = []
    max_score = None
    for row in rows:
        run_json = row["run_json"]
        site = int(row["site"])
        approximate = row.get("approximate", "").strip().lower() in ("1", "true", "yes")
        run_name = Path(run_json).name.removesuffix(".run.json")
        outfile = row.get("outfile") or str(output_dir / f"{run_name}.tagsite.png")

        try:
            _, score, max_score = plot_tagsite_figure.main(run_json, site, outfile, approximate=approximate)
        except Exception as e:
            print(f"WARNING: skipping {run_name}: {e}", file=sys.stderr)
            continue

        if score == score:  # skip NaN scores in the histogram
            scores.append(score)
        else:
            print(f"WARNING: no score computed for {run_name} at site {site}", file=sys.stderr)

    if scores:
        hist_path = plot_score_histogram(scores, max_score, str(output_dir / histogram_name))
        print(f"Wrote tag-site score histogram ({len(scores)} sites): {hist_path}")
    else:
        print("No tag-site scores collected; skipping histogram.", file=sys.stderr)


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "--manifest",
        required=True,
        help="CSV/TSV file, one row per protein (run_json, site, and optional "
        "outfile, approximate columns)",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="directory for per-protein PNGs and the batch score histogram",
    )
    args = parser.parse_args()

    main(args.manifest, args.output_dir)
