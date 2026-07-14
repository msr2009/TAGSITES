"""
plot_tagsite_figure.py

Renders a run's results plot (utils.results.plot_results) with a bold red
vertical line marking a known tag insertion site, and exports it as a PNG.
Used to visually check whether pipeline signals (conservation, pLDDT,
domains) line up with real CRISPR-validated internal tag sites.

Parameters:
    - run_json (str): path to a completed run's <run_name>.run.json
    - site (int): residue position (1-based) to mark
    - outfile (str): PNG output path
    - approximate (bool): if True, label the site as approximate in the title

Returns:
    - PNG written to outfile

Matt Rich, 7/2026
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import RESULTS_TYPE_DICT
from utils.results import load_data_from_json, plot_results


def main(run_json, site, outfile, approximate=False):
    aa_df, range_df, _alns = load_data_from_json(run_json, RESULTS_TYPE_DICT)

    run_name = Path(run_json).name.removesuffix(".run.json")
    title = f"{run_name} — tag site {site}" + (" (approximate)" if approximate else "")

    fig = plot_results(aa_df, range_df, title=title)
    fig.add_vline(x=site, line_color="red", line_width=3, row="all", col="all")
    fig.update_layout(title=title)

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    fig.write_image(outfile, width=1400, height=700, scale=2)
    return outfile


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

    args = parser.parse_args()

    main(args.RUN_JSON, args.SITE, args.OUTFILE, args.approximate)
