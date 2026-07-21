"""
rescore.py

Standalone CLI: recompute tag-site scores for a completed run using a given
scoring config, and rewrite the run's <run>.scores.tsv (see utils/scoring.py,
issue #32). Useful for iterating on weights/thresholds in scores.config.json
without re-running the analysis pipeline.

Matt Rich, 2026
"""

import sys
from pathlib import Path

# ensure the repo root is on sys.path so `utils`/`config` are importable
# whether this script is run directly or from another directory
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

from utils.scoring import write_run_scores, load_scoring_config, SCORING_CONFIG_PATH


def main(run_json_path, config_path=None, output_path=None):
    """Load a scoring config, rescore a run, and print where the TSV was written."""
    config = load_scoring_config(config_path)
    out_path = write_run_scores(run_json_path, config=config, output_path=output_path)
    print(f"Wrote tag-site scores: {out_path}")
    return out_path


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Rescore a completed TAGSITES run with a given scoring config.")
    parser.add_argument("-i", "--input", required=True, dest="RUN_JSON",
                        help="path to the run's .run.json file")
    parser.add_argument("-c", "--config", default=str(SCORING_CONFIG_PATH), dest="CONFIG",
                        help="path to a scores.config.json (default: repo scores.config.json)")
    parser.add_argument("-o", "--output", default=None, dest="OUTPUT",
                        help="output TSV path (default: <working_dir>/<run_name>.scores.tsv)")
    args = parser.parse_args()

    main(args.RUN_JSON, config_path=args.CONFIG, output_path=args.OUTPUT)
