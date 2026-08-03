"""
suggest_sites.py

Standalone CLI: rank and emit top-scoring tag sites for one or more completed runs,
optionally restricted to constitutive sequence or specific isoforms (see
utils.tag_filters.py, utils.scoring.build_user_criteria). Lets a large-scale batch
job design constitutive and isoform-unique sites for a whole gene set at once,
without going through the Results-tab Advanced/Rescoring accordion.

Matt Rich, 2026
"""

import csv
import json
import sys
from pathlib import Path

# ensure the repo root is on sys.path so `utils`/`config` are importable
# whether this script is run directly or from another directory
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

from config import RESULTS_TYPE_DICT
from utils.results import load_data_from_json, load_run_metadata, load_reagents_df
from utils.scoring import (
    score_tag_sites, load_scoring_config, pick_suggested_sites,
    with_extra_criteria, build_user_criteria, score_max, SCORING_CONFIG_PATH,
)
from utils.tag_filters import position_isoform_labels, describe_position_isoforms


def suggest_sites_for_run(run_json_path, config=None, n_sites=3, min_spacing=10,
                          isoform_mode="none", isoforms=None, topology=None):
    """Rank and return the top tag sites for one run.

    Returns a list of dicts with keys run_name, pos, amino_acid, score, max_score,
    isoform_specificity — one row per suggested site.
    """
    with open(run_json_path) as f:
        run_json = json.load(f)

    config = config or load_scoring_config()
    extra = build_user_criteria(run_json, isoform_mode=isoform_mode,
                                isoform_names=isoforms, topology_labels=topology)
    config = with_extra_criteria(config, extra)

    aa_df, range_df, _alns, iso_result = load_data_from_json(run_json, RESULTS_TYPE_DICT)
    meta = load_run_metadata(run_json)
    reagents_df = load_reagents_df(run_json)
    scores_df = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""), reagents_df, config)

    picked = pick_suggested_sites(scores_df, max_sites=n_sites, min_spacing=min_spacing)
    max_score = score_max(config)

    # isoform specificity is annotated for context only — never used to filter which
    # sites are suggested (see [[reagents-not-filtered-by-scoring-masks]] equivalent
    # reasoning: the isoform_mode restriction above already governs suggestion)
    isoforms_list = (iso_result or {}).get("isoforms", []) if iso_result else []
    labels_by_pos = position_isoform_labels(isoforms_list, meta.get("seq_len", 0))
    all_labels = [iso.get("name") or iso.get("accession") or "" for iso in isoforms_list]

    run_name = Path(run_json_path).name.removesuffix(".run.json")
    rows = []
    for pos in picked:
        row = scores_df.loc[pos]
        rows.append({
            "run_name": run_name,
            "pos": pos,
            "amino_acid": meta.get("query_seq", "")[pos - 1] if meta.get("query_seq") else "",
            "score": row["score"],
            "max_score": max_score,
            "isoform_specificity": describe_position_isoforms(pos, labels_by_pos, all_labels),
        })
    return rows


def main(run_json_paths, output_path=None, config_path=None, n_sites=3, min_spacing=10,
         isoform_mode="none", isoforms=None, topology=None):
    """Suggest tag sites across one or more runs and write (or print) a TSV."""
    config = load_scoring_config(config_path)
    rows = []
    for run_json_path in run_json_paths:
        rows.extend(suggest_sites_for_run(
            run_json_path, config=config, n_sites=n_sites, min_spacing=min_spacing,
            isoform_mode=isoform_mode, isoforms=isoforms, topology=topology,
        ))

    fieldnames = ["run_name", "pos", "amino_acid", "score", "max_score", "isoform_specificity"]
    out = open(output_path, "w", newline="") if output_path else sys.stdout
    try:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if output_path:
            out.close()
    return rows


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Suggest top-scoring tag sites for one or more runs, optionally "
                    "restricted to constitutive sequence or specific isoforms.")
    parser.add_argument("-i", "--input", required=True, dest="RUN_JSON", nargs="+",
                        help="path(s) to one or more run .run.json files")
    parser.add_argument("-o", "--output", default=None, dest="OUTPUT",
                        help="output TSV path (default: stdout)")
    parser.add_argument("-c", "--config", default=str(SCORING_CONFIG_PATH), dest="CONFIG",
                        help="path to a scores.config.json (default: repo scores.config.json)")
    parser.add_argument("--n-sites", type=int, default=3, dest="N_SITES",
                        help="maximum number of sites to suggest per run (default: 3)")
    parser.add_argument("--min-spacing", type=int, default=10, dest="MIN_SPACING",
                        help="minimum residue gap between suggested sites (default: 10)")
    parser.add_argument("--isoform-mode", choices=["none", "constitutive", "tagged"],
                        default="none", dest="ISOFORM_MODE",
                        help="restrict to constitutive sequence, or to specific isoforms "
                             "via --isoforms (default: no isoform restriction)")
    parser.add_argument("--isoforms", default=None, dest="ISOFORMS",
                        help="comma-separated isoform accessions/names, used with "
                             "--isoform-mode tagged")
    parser.add_argument("--topology", default=None, dest="TOPOLOGY",
                        help="comma-separated Phobius topology labels to restrict to, "
                             "e.g. Cytoplasmic (default: no topology restriction)")
    args = parser.parse_args()

    isoforms = [s.strip() for s in args.ISOFORMS.split(",")] if args.ISOFORMS else None
    topology = [s.strip() for s in args.TOPOLOGY.split(",")] if args.TOPOLOGY else None

    main(args.RUN_JSON, output_path=args.OUTPUT, config_path=args.CONFIG,
        n_sites=args.N_SITES, min_spacing=args.MIN_SPACING,
        isoform_mode=args.ISOFORM_MODE, isoforms=isoforms, topology=topology)
