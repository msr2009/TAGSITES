"""derive_isoforms.py — identify protein isoforms and map them onto query coordinates.

Two sourcing strategies, tried in order:
1. UniProt: pull the ALTERNATIVE PRODUCTS isoform list for the query accession, plus any
   "computationally mapped potential isoform sequences" (Ensembl-predicted transcripts
   UniProt imports into TrEMBL as separate entries sharing the same GeneID — see
   scripts/uniprot_api.py:fetch_computationally_mapped_isoforms), and align each isoform
   sequence to the query (scripts/isoform_align.py). Grounded in curated/computed
   UniProtKB annotation rather than local sequence inference.
2. BLAST fallback: infer isoforms from same-organism, same-gene, near-100%-identity BLAST
   hits already fetched by blast_orthologs.py. Caveat: these are alignment-derived
   segments, not genomic exon boundaries, and a point mutation can look identical to a
   splice difference — used only when UniProt has no isoform data for the query at all.

Output is a JSON dict (see main()) rather than the flat presence-class TSV this script
used to emit; every isoform is a distinct record so callers can render one row each,
rather than a single collapsed constitutive/intermediate/unique track.
"""

import json
import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from isoform_align import align_isoform_to_query, _spans_from_positions
from uniprot_api import fetch_isoform_sequences


def _hseq_sig(hsp):
    """Gap-stripped hit sequence — used to deduplicate isoform hits."""
    return hsp.get("hsp_hseq", "").replace("-", "")


def _check_identity(hsp, min_perfect_len):
    """Return True if the HSP has a consecutive perfect-match run >= min_perfect_len.

    Organism + gene-name filtering already excludes paralogs; this guards against
    spurious same-gene hits with no substantial shared block.
    """
    qseq = hsp.get("hsp_qseq", "")
    hseq = hsp.get("hsp_hseq", "")
    if not qseq or not hseq:
        return True  # no alignment strings; accept if it passed organism/gene filters
    max_run = cur_run = 0
    for q, h in zip(qseq, hseq):
        if q != "-" and h != "-" and q == h:
            cur_run += 1
            if cur_run > max_run:
                max_run = cur_run
        else:
            cur_run = 0
    return max_run >= min_perfect_len


def _spans_from_hsp(hsp, query_len):
    """Map one BLAST HSP's aligned strings onto query coordinates.

    Positions outside the HSP's query span are marked skipped (the hit sequence
    wasn't observed there), matching what align_isoform_to_query would report for
    a fully-aligned pair.
    """
    qseq = hsp.get("hsp_qseq", "")
    hseq = hsp.get("hsp_hseq", "")
    qfrom = int(hsp.get("hsp_query_from", 1))
    qto = int(hsp.get("hsp_query_to", query_len))

    present, skipped, inserts = [], [], []
    qpos = qfrom - 1
    cur_insert_len = 0

    def _flush_insert():
        nonlocal cur_insert_len
        if cur_insert_len:
            inserts.append({"after": qpos, "length": cur_insert_len})
            cur_insert_len = 0

    for q, h in zip(qseq, hseq):
        if q != "-":
            _flush_insert()
            qpos += 1
            (present if h == q else skipped).append(qpos)
        elif h != "-":
            cur_insert_len += 1
    _flush_insert()

    skipped.extend(range(1, qfrom))
    skipped.extend(range(qto + 1, query_len + 1))

    return {
        "present": _spans_from_positions(sorted(set(present))),
        "skipped": _spans_from_positions(sorted(set(skipped))),
        "inserts": inserts,
    }


def _blast_fallback_isoforms(blast_output, min_perfect_len=40):
    """Infer per-isoform records from same-organism, same-gene BLAST hits.

    Returns a list of isoform records (see main() schema), or [] when fewer than
    2 distinct isoforms are detected (single-isoform gene → no track).
    """
    hits = blast_output.get("hits", [])
    if not hits:
        return []
    query_len = int(blast_output.get("query_len", 0))
    if query_len == 0:
        return []

    # the top hit is the self-hit; use it to establish query organism and gene
    query_ox = hits[0].get("hit_uni_ox", "")
    query_gn = hits[0].get("hit_uni_gn", "")

    # filter hits to same-organism, same-gene candidates with a shared conserved block
    candidates = []
    for h in hits:
        hsp = h["hit_hsps"][0]
        hit_ox = h.get("hit_uni_ox", "")
        hit_gn = h.get("hit_uni_gn", "")
        if not hit_ox or hit_ox != query_ox:
            continue
        # gene name filter: skip if both names are known but differ
        if query_gn and hit_gn and hit_gn != query_gn:
            continue
        if not _check_identity(hsp, min_perfect_len):
            continue
        candidates.append(h)

    # deduplicate by gap-stripped hit sequence so the canonical UniProt entry and its
    # -1 isoform variant (identical aligned region) don't double-count
    seen = set()
    unique_hits = []
    for h in candidates:
        sig = _hseq_sig(h["hit_hsps"][0])
        if sig not in seen:
            seen.add(sig)
            unique_hits.append(h)

    if len(unique_hits) <= 1:
        return []  # single isoform → nothing to show

    records = []
    for h in unique_hits:
        hsp = h["hit_hsps"][0]
        spans = _spans_from_hsp(hsp, query_len)
        is_query = h is hits[0]
        records.append({
            "accession": h.get("hit_acc", ""),
            "name": "query" if is_query else h.get("hit_acc", ""),
            "length": len(hsp.get("hsp_hseq", "").replace("-", "")) or query_len,
            "is_query": is_query,
            **spans,
        })
    return records


def derive_isoforms(blast_output, min_perfect_len=40):
    """Derive per-isoform records for a query, preferring UniProt over BLAST inference.

    Returns a dict: {"query_len", "query_accession", "source", "isoforms", "caveat"?}.
    "isoforms" is [] when only one isoform is detected (nothing to show).
    """
    hits = blast_output.get("hits", [])
    query_len = int(blast_output.get("query_len", 0))
    query_acc = hits[0].get("hit_acc", "") if hits else ""

    if query_acc:
        uniprot_records = fetch_isoform_sequences(query_acc)
        if len(uniprot_records) >= 2:
            query_seq = None
            isoforms = []
            for rec in uniprot_records:
                if rec["is_canonical"]:
                    query_seq = rec["sequence"]
            for rec in uniprot_records:
                if rec["is_canonical"]:
                    spans = {"present": [[1, len(rec["sequence"])]], "skipped": [], "inserts": []}
                else:
                    spans = (align_isoform_to_query(query_seq, rec["sequence"])
                             if query_seq else {"present": [], "skipped": [], "inserts": []})
                isoforms.append({
                    "accession": rec["accession"],
                    "name": rec["isoform_name"] or rec["accession"],
                    "length": len(rec["sequence"]),
                    "is_query": rec["is_canonical"],
                    **spans,
                })
            return {
                "query_len": query_len,
                "query_accession": query_acc,
                "source": "uniprot",
                "isoforms": isoforms,
            }

    # fallback: infer from BLAST hits
    fallback = _blast_fallback_isoforms(blast_output, min_perfect_len=min_perfect_len)
    result = {
        "query_len": query_len,
        "query_accession": query_acc,
        "source": "blast",
        "isoforms": fallback,
    }
    if fallback:
        result["caveat"] = ("Isoforms inferred from BLAST identity, not UniProt annotation — "
                             "a point mutation can appear identical to a splice difference.")
    return result


def main(blast_json_path, output_path, min_perfect_len=40):
    """Derive isoform records from a BLAST JSON file and write an isoforms JSON file."""
    with open(blast_json_path) as f:
        blast_output = json.load(f)
    result = derive_isoforms(blast_output, min_perfect_len=int(min_perfect_len))
    with open(output_path, "w") as f:
        json.dump(result, f)


if __name__ == "__main__":
    parser = ArgumentParser(description="Derive per-isoform coordinate mappings from BLAST hits")
    parser.add_argument("-i", "--input", "--blast_json", dest="BLAST_JSON",
                        required=True, help="path to *.json.json BLAST output")
    parser.add_argument("-o", "--output", dest="OUTPUT", required=True,
                        help="path to write isoforms JSON")
    parser.add_argument("--min_perfect_len", dest="MIN_PERFECT_LEN", default=40, type=int,
                        help="minimum consecutive perfect-match run length (default 40)")
    args, _ = parser.parse_known_args()
    main(args.BLAST_JSON, args.OUTPUT, args.MIN_PERFECT_LEN)
