"""uniprot_features.py

Fetches a UniProt entry and emits a curated subset of its feature table (PTMs,
lipidation, propeptides, binding sites, etc.) as range annotations, following
the same conventions as regex_sites.py / call_interpro.py.

Two output sources, distinguished by severity in tables/uniprot_features.txt:
  UniProt       - "blocking" features (lipidation, propeptide, signal peptide, ...):
                  tagging here is expected to disrupt processing/targeting/folding.
  UniProt_site  - "site" features (binding site, motif, PTM, ...): informational only.

Matt Rich, 2026
"""

import sys

from site_selection_util import get_sequence, read_fasta, uniprot_accession_regex
from uniprot_api import checksum_lookup, fetch_entry
from progress import report as _report, resolve_reporter


def load_feature_config(path):
    """Parse tables/uniprot_features.txt into {feature_type: (severity, label)}.

    Blank lines and lines starting with '#' are skipped. Fields are tab-delimited;
    consecutive tabs used for visual alignment are collapsed.
    """
    config = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [f for f in line.split("\t") if f != ""]
        if len(fields) < 3:
            continue
        feature_type, severity, label = fields[0], fields[1], fields[2]
        config[feature_type] = (severity, label)
    return config


def _feature_span(feature):
    """Return (start, stop) 1-based inclusive from a UniProt feature's location block."""
    loc = feature.get("location", {})
    start = loc.get("start", {}).get("value")
    stop = loc.get("end", {}).get("value")
    if start is None:
        return None
    if stop is None:
        stop = start
    return start, stop


def _feature_detail(feature, ftype):
    """Best available human-readable detail for a feature.

    UniProt often leaves 'description' blank for Binding site features and puts
    the actual info (e.g. the ligand name) in a separate 'ligand' field instead.
    Falls back to the feature type name if nothing else is available.
    """
    detail = feature.get("description") or ""
    if not detail:
        ligand = feature.get("ligand") or {}
        detail = ligand.get("name", "")
    return detail or ftype


def _evidence_suffix(feature):
    """' (PubMed:id1, id2)' for experimentally-cited features, ' (predicted)' for
    rule-based/inferred ones (ECO codes with no PubMed source), '' if unevidenced.
    """
    evidences = feature.get("evidences") or []
    pubmed_ids = []
    for e in evidences:
        if e.get("source") == "PubMed" and e.get("id") and e["id"] not in pubmed_ids:
            pubmed_ids.append(str(e["id"]))
    if pubmed_ids:
        return f" (PubMed:{', '.join(pubmed_ids)})"
    if evidences:
        return " (predicted)"
    return ""


def entry_to_feature_rows(entry_json, feature_config):
    """Convert a UniProt entry JSON into a list of (source, start, stop, description) rows.

    Only feature types present in feature_config are kept; everything else is dropped.
    Description includes the ligand name (for binding sites whose own description is
    blank) and an evidence suffix distinguishing experimentally-cited PubMed support
    from rule-based/predicted annotations.
    """
    rows = []
    for feature in entry_json.get("features", []) or []:
        ftype = feature.get("type")
        if ftype not in feature_config:
            continue
        span = _feature_span(feature)
        if span is None:
            continue
        start, stop = span
        severity, label = feature_config[ftype]
        source = "UniProt" if severity == "blocking" else "UniProt_site"
        detail = _feature_detail(feature, ftype) + _evidence_suffix(feature)
        rows.append((source, start, stop, f"{label}: {detail}"))
    return rows


def main(fasta_in, output, features_file, accession="", taxid="", report=None):
    """Fetch UniProt features for a protein and write them as a range TSV.

    Resolves an accession (explicit arg > accession-shaped FASTA header > checksum
    lookup), fetches the entry, and writes matching features. If no accession can
    be resolved, writes an empty file rather than raising -- an empty range file
    is handled fine by the rest of the pipeline.
    """
    reporter = resolve_reporter(report)
    feature_config = load_feature_config(features_file)

    acc = accession.strip() if accession else ""
    if not acc:
        name, _ = read_fasta(fasta_in)
        m = uniprot_accession_regex(name)
        if m is not None:
            acc = m.group(0)
    if not acc:
        _report(reporter, "No accession given; looking up UniProt by sequence checksum…", stage="uniprot")
        seq = get_sequence(fasta_in)
        try:
            acc = checksum_lookup(seq, taxid) or ""
        except Exception as e:
            _report(reporter, f"checksum lookup failed ({e})", stage="uniprot", level="warning")
            acc = ""
    else:
        # re.match() only anchors the start, so an isoform-suffixed accession
        # (e.g. "G5EE56-1") matches on its base prefix but was previously used
        # whole -- the UniProt entry endpoint accepts isoform accessions but
        # silently returns an incomplete/empty feature list for them, so we
        # always resolve to the canonical (base) accession.
        m = uniprot_accession_regex(acc)
        if m is not None:
            acc = m.group(0)

    if not acc:
        _report(reporter, "No UniProt accession found; writing empty feature table.", stage="uniprot", level="warning")
        open(output, "w").close()
        return

    _report(reporter, f"Fetching UniProt entry {acc}…", stage="uniprot")
    try:
        entry = fetch_entry(acc)
    except Exception as e:
        _report(reporter, f"fetching UniProt entry {acc} failed ({e}); writing empty feature table.",
                stage="uniprot", level="warning")
        open(output, "w").close()
        return
    rows = entry_to_feature_rows(entry, feature_config)

    with open(output, "w") as fout:
        for source, start, stop, desc in rows:
            print("\t".join(str(x) for x in (source, start, stop, desc)), file=fout)
    _report(reporter, f"Wrote {len(rows)} UniProt feature(s).", stage="uniprot")


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', '--input_file', action='store', type=str, dest='FASTA',
        help="fasta file containing protein sequence", required=True)
    parser.add_argument('--features_file', action='store', type=str, dest='FEATURES_FILE',
        help="tab-delimited file of UniProt feature types to track", default="./tables/uniprot_features.txt")
    parser.add_argument('--accession', action='store', type=str, dest='ACCESSION',
        help="UniProt accession to use directly, skipping lookup", default="")
    parser.add_argument('--taxid', action='store', type=str, dest='TAXID',
        help="taxonomy ID to constrain the checksum lookup", default="")
    parser.add_argument('--output', action='store', type=str, dest='OUTPUT_FILE',
        help="file to store output", required=True)

    args, unknowns = parser.parse_known_args()

    main(args.FASTA, args.OUTPUT_FILE, args.FEATURES_FILE, args.ACCESSION, args.TAXID)
