"""derive_isoforms.py — identify protein isoforms from BLAST hits.

Detects same-organism, near-100%-identity BLAST hits (UniProt isoforms) and
classifies each query residue as constitutive, intermediate, or unique based on
what fraction of isoforms contain it identically.

Caveats: these are isoform-coverage segments inferred from protein alignment, not
genomic exon boundaries. Detection depends on isoforms being present in UniProtKB
as separate accessions.
"""

import json
import sys
from argparse import ArgumentParser
from pathlib import Path


def _hseq_sig(hsp):
    """Gap-stripped hit sequence — used to deduplicate isoform hits."""
    return hsp.get("hsp_hseq", "").replace("-", "")


def _check_identity(hsp, min_identity, min_perfect_len):
    """Return True if the HSP meets identity and perfect-run-length thresholds."""
    qseq = hsp.get("hsp_qseq", "")
    hseq = hsp.get("hsp_hseq", "")
    if not qseq or not hseq:
        # no alignment strings — fall back to summary identity percentage
        return float(hsp.get("hsp_identity", 0)) / 100.0 >= min_identity
    # fraction of non-gap query positions that match hit
    total = sum(1 for c in qseq if c != "-")
    if total == 0:
        return False
    identical = sum(1 for q, h in zip(qseq, hseq) if q != "-" and h != "-" and q == h)
    if identical / total < min_identity:
        return False
    # require at least one consecutive perfect-match run >= min_perfect_len
    # (guards against same-organism hits that are deeply gapped paralogs)
    max_run = cur_run = 0
    for q, h in zip(qseq, hseq):
        if q != "-" and h != "-" and q == h:
            cur_run += 1
            if cur_run > max_run:
                max_run = cur_run
        else:
            cur_run = 0
    return max_run >= min_perfect_len


def _mark_presence(hsp, presence, query_len):
    """Increment presence[pos] for each query position covered identically by this HSP."""
    qseq = hsp.get("hsp_qseq", "")
    hseq = hsp.get("hsp_hseq", "")
    qfrom = int(hsp.get("hsp_query_from", 1))
    if not qseq or not hseq:
        # fallback when alignment strings are absent: mark the full span
        qto = int(hsp.get("hsp_query_to", query_len))
        for pos in range(qfrom, min(qto, query_len) + 1):
            presence[pos] += 1
        return
    qpos = qfrom
    for q, h in zip(qseq, hseq):
        if q == "-":
            continue  # query insertion — no query position consumed
        # q is a real query residue; count it if hit has same residue at this column
        if 1 <= qpos <= query_len and h != "-" and q == h:
            presence[qpos] += 1
        qpos += 1


def derive_isoform_segments(blast_output, min_perfect_len=40, min_identity=0.98):
    """Classify query residues by isoform coverage from a BLAST result dict.

    Returns a list of (start, stop, description) tuples (1-based, inclusive),
    or [] when <=1 distinct isoform is detected (single-isoform gene → no track).
    Description is one of 'constitutive (N/N isoforms)', 'intermediate (k/N isoforms)',
    or 'unique (1/N isoforms)'.
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

    # filter hits to same-organism, same-gene, near-100%-identity candidates
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
        if not _check_identity(hsp, min_identity, min_perfect_len):
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

    n_iso = len(unique_hits)
    if n_iso <= 1:
        return []  # single isoform → nothing to show

    # per-residue presence counts (1-based; index 0 unused)
    presence = [0] * (query_len + 1)
    for h in unique_hits:
        _mark_presence(h["hit_hsps"][0], presence, query_len)

    # classify each position, then collapse consecutive same-class runs into segments
    def classify(count):
        if count >= n_iso:
            return f"constitutive ({n_iso}/{n_iso} isoforms)"
        elif count <= 1:
            return f"unique (1/{n_iso} isoforms)"
        else:
            return f"intermediate ({count}/{n_iso} isoforms)"

    segments = []
    cur_start = 1
    cur_cls = classify(presence[1])
    for pos in range(2, query_len + 1):
        cls = classify(presence[pos])
        if cls != cur_cls:
            segments.append((cur_start, pos - 1, cur_cls))
            cur_start = pos
            cur_cls = cls
    segments.append((cur_start, query_len, cur_cls))

    return segments


def main(blast_json_path, output_path, min_perfect_len=40, min_identity=0.98):
    """Derive isoform segments from a BLAST JSON file and write a range TSV."""
    with open(blast_json_path) as f:
        blast_output = json.load(f)
    segments = derive_isoform_segments(
        blast_output,
        min_perfect_len=int(min_perfect_len),
        min_identity=float(min_identity),
    )
    with open(output_path, "w") as f:
        for start, stop, desc in segments:
            f.write(f"isoforms\t{start}\t{stop}\t{desc}\n")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    parser = ArgumentParser(description="Derive isoform coverage segments from BLAST hits")
    parser.add_argument("-i", "--input", "--blast_json", dest="BLAST_JSON",
                        required=True, help="path to *.json.json BLAST output")
    parser.add_argument("-o", "--output", dest="OUTPUT", required=True,
                        help="path to write isoforms TSV")
    parser.add_argument("--min_perfect_len", dest="MIN_PERFECT_LEN", default=40, type=int,
                        help="minimum consecutive perfect-match run length (default 40)")
    parser.add_argument("--min_identity", dest="MIN_IDENTITY", default=0.98, type=float,
                        help="minimum fractional identity (default 0.98)")
    args, _ = parser.parse_known_args()
    main(args.BLAST_JSON, args.OUTPUT, args.MIN_PERFECT_LEN, args.MIN_IDENTITY)
