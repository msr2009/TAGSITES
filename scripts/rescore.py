"""
rescore.py

Standalone CLI: recompute tag-site scores for a completed run using a given
scoring config, and rewrite the run's <run>.scores.tsv (see utils/scoring.py,
issue #32). Useful for iterating on weights/thresholds in scores.config.json
without re-running the analysis pipeline.

Matt Rich, 2026
"""

import json
import sys
from pathlib import Path

# ensure the repo root is on sys.path so `utils`/`config` are importable
# whether this script is run directly or from another directory
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

from utils.scoring import (
    write_run_scores, load_scoring_config, with_extra_criteria, build_user_criteria,
    SCORING_CONFIG_PATH,
)


def main(run_json_path, config_path=None, output_path=None,
         isoform_mode="none", isoforms=None, topology=None):
    """Load a scoring config, rescore a run, and print where the TSV was written.

    isoform_mode/isoforms/topology default to no restriction, so the emitted TSV is
    unchanged from before these flags existed unless explicitly requested (see
    utils.scoring.build_user_criteria — used for large-scale constitutive/isoform-unique
    site design, see scripts/suggest_sites.py).
    """
    config = load_scoring_config(config_path)
    if isoform_mode != "none" or topology:
        with open(run_json_path) as f:
            run_json = json.load(f)
        extra = build_user_criteria(run_json, isoform_mode=isoform_mode,
                                    isoform_names=isoforms, topology_labels=topology)
        config = with_extra_criteria(config, extra)
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
    parser.add_argument("--isoform-mode", choices=["none", "constitutive", "tagged"],
                        default="none", dest="ISOFORM_MODE",
                        help="restrict scoring to constitutive sequence, or to specific "
                             "isoforms via --isoforms (default: no isoform restriction)")
    parser.add_argument("--isoforms", default=None, dest="ISOFORMS",
                        help="comma-separated isoform accessions/names, used with "
                             "--isoform-mode tagged")
    parser.add_argument("--topology", default=None, dest="TOPOLOGY",
                        help="comma-separated Phobius topology labels to restrict to, "
                             "e.g. Cytoplasmic (default: no topology restriction)")
    args = parser.parse_args()

    isoforms = [s.strip() for s in args.ISOFORMS.split(",")] if args.ISOFORMS else None
    topology = [s.strip() for s in args.TOPOLOGY.split(",")] if args.TOPOLOGY else None

    main(args.RUN_JSON, config_path=args.CONFIG, output_path=args.OUTPUT,
        isoform_mode=args.ISOFORM_MODE, isoforms=isoforms, topology=topology)
