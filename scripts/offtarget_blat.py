"""
offtarget_blat.py

UCSC BLAT backend for Screen A of the off-target / primer-specificity screen.

Drop-in alternative to offtarget_remote.run_region_screen(): same signature, same
return shape, same EBI-shaped HSP dicts, so offtarget_screen.py is untouched.
Named *_blat.py to match the scripts/providers.py <analysis>_<mode> convention.

Why this exists: the EBI region query is the whole cost of the off-target screen
(173 s for snt-1, 764 s for col-103, both measured). BLAT answers the same
question against the same assembly in 0.6-2.3 s, which is what makes it possible
to leave the screen on by default.

Screen A only. BLAT cannot do Screen B: a guaranteed hit needs 2*stepSize +
tileSize - 1 = 20 perfect bases, and nothing shorter than stepSize + tileSize =
16 is reported at all, so a 20 nt spacer carrying 3 mismatches sits below the
tiling floor. Screen B stays on offtarget_remote.

Three properties of the UCSC API that shape this module, none of them documented
upstream and all verified against the live service:

  * GET is unusable. A 9.4 kb query returns 414 Request-URI Too Long, so the
    documented 75,000-base limit is unreachable that way. POST works.
  * No output format carries aligned sequence. psl, text and hgblat all return
    PSL columns only, so qseq/hseq have to be rebuilt from target sequence.
  * A PSL record is NOT a BLAST HSP. BLAT chains blocks across large target
    gaps — measured up to 16 blocks spanning 395 kb for 127 aligned bases. Fed
    in whole, such a record would produce phantom amplicons hundreds of kb wide.
    Every record is therefore split into one pseudo-HSP per ungapped block.

Matt Rich, 2026
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import http_retry
import offtarget_screen as ots
from crispr_util import reverse_complement
from progress import report as _report, resolve_reporter

BLAT_URL = "https://api.genome.ucsc.edu/blat/dna"
SEQUENCE_URL = "https://api.genome.ucsc.edu/getData/sequence"
GENOMES_URL = "https://api.genome.ucsc.edu/list/ucscGenomes"

KEY_ENV = "TAGSITES_UCSC_API_KEY"
KEY_FILE = Path(__file__).parent.parent / "ucsc.local.json"

# Karlin-Altschul parameters for ungapped blastn at reward +1 / penalty -2, the
# effective scoring of an exact PSL block. PSL carries no bit score, but
# post_filter.min_bitscore is a real threshold in offtarget_screen.py (which must
# stay byte-unchanged), so a score has to be supplied rather than omitted. These
# reproduce blastn's own statistics for the same alignment closely enough to keep
# the threshold meaning what it meant; see _bitscore().
_LAMBDA = 1.28
_K = 0.46


def load_api_key(cfg=None):
    """UCSC apiKey from the environment, else from the gitignored ucsc.local.json.

    shinyapps.io has no environment-variable or secrets mechanism, so the deployed
    app reads the file (which ships inside the private bundle) while CLI and local
    runs can use the env var instead.
    """
    key = os.environ.get(KEY_ENV, "").strip()
    if key:
        return key
    try:
        with open(KEY_FILE) as f:
            return str(json.load(f).get("api_key", "")).strip()
    except (OSError, ValueError):
        return ""


def assembly_for_taxid(taxid, cfg):
    """UCSC assembly name for a taxid, or '' when none is configured."""
    table = (cfg.get("blat") or {}).get("assemblies") or {}
    return str(table.get(str(taxid).strip(), "")).strip()


def available(taxid, cfg):
    """True when this taxid can be screened by BLAT (assembly known, key present)."""
    return bool(assembly_for_taxid(taxid, cfg)) and bool(load_api_key(cfg))


def fetch_assembly_table(deadline=None):
    """taxid -> assembly for every UCSC genome, for refreshing the config map.

    Not called at runtime; offtarget.config.json carries a static map so the app
    never depends on this endpoint. Kept here so the map can be regenerated.
    """
    resp = http_retry.request_with_retries("get", GENOMES_URL, timeout=60,
                                           deadline=deadline)
    out = {}
    for name, g in (resp.json().get("ucscGenomes") or {}).items():
        tax = str(g.get("taxId") or "").strip()
        # several assemblies share a taxid; keep the first and let the config pin it
        if tax and tax not in out:
            out[tax] = name
    return out


# ── BLAT ──────────────────────────────────────────────────────────────────────

def blat_dna(sequence, genome, api_key, deadline=None, max_items=None):
    """POST one BLAT DNA query; returns the list of PSL records.

    POST rather than GET: a multi-kb query overflows the server's URI limit.
    """
    body = {"genome": genome, "userSeq": str(sequence), "apiKey": api_key}
    if max_items:
        body["maxItemsOutput"] = int(max_items)
    resp = http_retry.request_with_retries("post", BLAT_URL, data=body,
                                           timeout=(10, 300), deadline=deadline)
    payload = resp.json()
    if payload.get("error"):
        raise RuntimeError("UCSC BLAT: {}".format(payload["error"]))
    return payload.get("blat", []) or []


def _ints(value):
    """PSL comma-list (or already-split list) to a list of ints."""
    if isinstance(value, str):
        return [int(x) for x in value.split(",") if x.strip() != ""]
    return [int(x) for x in value]


def split_psl_blocks(rec):
    """One ungapped block per entry: (q_start, t_start, size), 0-based, in target order.

    Coordinates are left in PSL's own frame — for a minus-strand record qStarts are
    on the reverse-complemented query — because that is the frame in which both
    query and target coordinates increase together, which is what group_blocks needs
    to measure gaps. Conversion to forward-query coordinates happens in psl_to_hsps.
    """
    return list(zip(_ints(rec["qStarts"]), _ints(rec["tStarts"]),
                    _ints(rec["blockSizes"])))


def group_blocks(blocks, max_gap=20):
    """Group consecutive blocks into local alignments, splitting at large gaps.

    BLAT chains blocks across enormous target gaps — measured up to 395 kb for 127
    aligned bases — so a whole record cannot be treated as one local alignment.
    Splitting only where a gap exceeds max_gap keeps ordinary gapped alignments
    whole while still breaking the chimeric chains apart. Keep max_gap small: a
    large budget merges distant blocks into a long alignment padded with gap
    columns, whose identity then falls under post_filter.min_identity, so it is
    dropped entirely. Measured on col-103, guide windows flagged fall from 210 at
    a 20 bp budget to 129 at 1000 bp. See offtarget.config.json's blat block.
    """
    groups = []
    for blk in blocks:
        if groups:
            pq, pt, ps = groups[-1][-1]
            q_gap, t_gap = blk[0] - (pq + ps), blk[1] - (pt + ps)
            # a negative gap means the blocks overlap or run backwards on the query,
            # which only happens across a chain boundary, so split there too
            if 0 <= q_gap <= max_gap and 0 <= t_gap <= max_gap:
                groups[-1].append(blk)
                continue
        groups.append([blk])
    return groups


# ── Target sequence ───────────────────────────────────────────────────────────

def _merge_ranges(ranges, max_gap=200000):
    """Merge sorted (start, end) ranges separated by less than max_gap.

    max_gap is deliberately large. Each fetch costs a round trip plus a politeness
    pause, while the sequence itself is cheap, so pulling one 200 kb window beats
    pulling twenty scattered 60 bp ones. Measured on col-103: 45.5 s to 4.6 s.
    """
    out = []
    for s, e in sorted(ranges):
        if out and s - out[-1][1] <= max_gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def fetch_sequence(genome, chrom, start, end, deadline=None):
    """Target sequence for a 0-based half-open range; needs no apiKey."""
    resp = http_retry.request_with_retries(
        "get", SEQUENCE_URL, timeout=(10, 120), deadline=deadline,
        params={"genome": genome, "chrom": chrom, "start": int(start), "end": int(end)})
    payload = resp.json()
    if payload.get("error"):
        raise RuntimeError("UCSC getData/sequence: {}".format(payload["error"]))
    return str(payload.get("dna", "")).upper()


def _target_cache(genome, needed, deadline=None, pause=0.0, max_gap=200000):
    """Fetch merged target ranges per chromosome; returns {chrom: [(start, end, seq)]}.

    Blocks are merged before fetching because adjacent blocks of one chained PSL
    record sit within a few kb of each other, so one request usually covers many.
    """
    cache = {}
    for chrom, ranges in sorted(needed.items()):
        spans = []
        for start, end in _merge_ranges(ranges, max_gap=max_gap):
            spans.append((start, end, fetch_sequence(genome, chrom, start, end,
                                                     deadline=deadline)))
            if pause:
                time.sleep(pause)
        cache[chrom] = spans
    return cache


def _slice(cache, chrom, start, end):
    """Target bases for a 0-based half-open range out of the fetch cache."""
    for s, e, seq in cache.get(chrom, []):
        if start >= s and end <= e:
            return seq[start - s:end - s]
    return ""


# ── PSL -> EBI-shaped HSPs ────────────────────────────────────────────────────

def _bitscore(matches, mismatches):
    """Bit score for an ungapped block under blastn's reward +1 / penalty -2."""
    raw = matches - 2 * mismatches
    return max(0.0, (_LAMBDA * raw - math.log(_K)) / math.log(2))


def _evalue(bits, query_len, target_len):
    """E-value from a bit score and the search space."""
    try:
        return float(query_len) * float(target_len) * math.pow(2.0, -bits)
    except (OverflowError, ValueError):
        return 0.0


def _revcomp_gapped(seq):
    """Reverse complement an aligned string, leaving gap characters as gaps."""
    return "".join("-" if c == "-" else reverse_complement(c) for c in reversed(seq))


def psl_to_hsps(records, query_seq, genome, min_align_len=0, deadline=None, pause=0.0,
                max_gap=200000, group_gap=20):
    """Convert PSL records into EBI-shaped HSP dicts, one per ungapped block.

    Blocks shorter than min_align_len are dropped BEFORE any sequence is fetched —
    that is what keeps reconstruction cheap, since the post-filter would discard
    them anyway (93 of 796 blocks survive for col-103, 22 of 241 for snap29).
    """
    query_seq = str(query_seq).upper()
    rc_query = reverse_complement(query_seq)
    groups = []
    needed = {}
    for rec in records:
        chrom = rec["tName"]
        minus = str(rec["strand"]).startswith("-")
        for blocks in group_blocks(split_psl_blocks(rec), max_gap=group_gap):
            span = sum(s for _, _, s in blocks)
            if span < min_align_len:
                continue
            t0, t1 = blocks[0][1], blocks[-1][1] + blocks[-1][2]
            groups.append((rec, chrom, minus, blocks))
            needed.setdefault(chrom, []).append((t0, t1))
    if not groups:
        return []

    cache = _target_cache(genome, needed, deadline=deadline, pause=pause,
                          max_gap=max_gap)

    hsps = []
    for rec, chrom, minus, blocks in groups:
        q_src = rc_query if minus else query_seq
        qcols, hcols = [], []
        prev = None
        ok = True
        for qs, ts, size in blocks:
            if prev is not None:
                # gaps between blocks: unmatched bases on one side, '-' on the other
                pq, pt, ps = prev
                q_gap, t_gap = qs - (pq + ps), ts - (pt + ps)
                qcols.append(q_src[pq + ps:qs] + "-" * t_gap)
                hcols.append("-" * q_gap + _slice(cache, chrom, pt + ps, ts))
            tseq = _slice(cache, chrom, ts, ts + size)
            qseg = q_src[qs:qs + size]
            if len(tseq) != size or len(qseg) != size:
                ok = False
                break
            qcols.append(qseg)
            hcols.append(tseq)
            prev = (qs, ts, size)
        if not ok:
            continue
        qseq, hseq = "".join(qcols), "".join(hcols)
        if len(qseq) != len(hseq):
            continue
        # a minus-strand alignment is built in the reverse-complemented query frame;
        # flipping both strings puts it back on the forward query, which is the sense
        # EBI's qseq/hseq pair already uses
        if minus:
            qseq, hseq = _revcomp_gapped(qseq), _revcomp_gapped(hseq)
        matches = sum(1 for a, b in zip(qseq, hseq) if a == b and a != "-")
        align_len = len(qseq)
        bits = _bitscore(matches, align_len - matches)
        q_lo = min(b[0] for b in blocks)
        q_hi = max(b[0] + b[2] for b in blocks)
        if minus:
            q_lo, q_hi = len(query_seq) - q_hi, len(query_seq) - q_lo
        t0 = blocks[0][1]
        t1 = blocks[-1][1] + blocks[-1][2]
        # One acc per local alignment, NOT per chromosome. classify_hsps groups by acc
        # and calls an accession "self" when its hits cover most of the query at high
        # identity — right for ENA (one accession is one record), wrong for UCSC (one
        # chromosome holds every locus on it). Keyed by chromosome, every hit on the
        # query's own chromosome inherited its self-classification and was dropped.
        hsps.append({
            "acc":       "{}:{}-{}".format(chrom, t0, t1),
            "desc":      "{} {}:{}-{} ({})".format(genome, chrom, t0 + 1, t1,
                                                   "-" if minus else "+"),
            "os":        genome,
            "q_from":    q_lo + 1,
            "q_to":      q_hi,
            # subject coordinates stay in the target's own forward numbering, with
            # h_from > h_to on the minus strand, matching BLAST's convention
            "h_from":    t1 if minus else (t0 + 1),
            "h_to":      (t0 + 1) if minus else t1,
            "h_strand":  "-" if minus else "+",
            "qseq":      qseq,
            "hseq":      hseq,
            "identity":  100.0 * matches / max(1, align_len),
            "bits":      bits,
            "evalue":    _evalue(bits, len(query_seq), int(rec.get("tSize", 0) or 0)),
            "align_len": align_len,
        })
    return hsps


# ── Screen A ──────────────────────────────────────────────────────────────────

def run_region_screen(region_seq, email, taxid, exons=None, cfg=None, report=None,
                      job_id_cb=None, resume_job_ids=None):
    """Screen A against a UCSC assembly; same contract as offtarget_remote's.

    email, job_id_cb and resume_job_ids are accepted and ignored: BLAT needs no
    submitter address and returns synchronously, so there is no job to resume.
    They stay in the signature so this module is interchangeable with the EBI one.
    """
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    blat_cfg = cfg.get("blat") or {}
    genome = assembly_for_taxid(taxid, cfg)
    if not genome:
        raise RuntimeError("no UCSC assembly configured for taxid {}".format(taxid))
    api_key = load_api_key(cfg)
    if not api_key:
        raise RuntimeError("no UCSC API key ({} or {})".format(KEY_ENV, KEY_FILE.name))

    deadline = http_retry.deadline_from(blat_cfg.get("total_timeout", 300))
    _report(reporter, "Submitting region off-target BLAT ({} bp vs {}, taxid {})".format(
        len(region_seq), genome, taxid or "unscoped"), stage="offtarget_region")

    t0 = time.perf_counter()
    records = blat_dna(region_seq, genome, api_key, deadline=deadline,
                       max_items=blat_cfg.get("max_items"))
    hsps = psl_to_hsps(records, region_seq, genome,
                       min_align_len=int(cfg["post_filter"]["min_align_len"]),
                       deadline=deadline,
                       pause=float(blat_cfg.get("request_pause", 0.0)),
                       max_gap=int(blat_cfg.get("fetch_merge_gap", 200000)),
                       group_gap=int(blat_cfg.get("max_alignment_gap", 20)))
    elapsed = time.perf_counter() - t0

    ots.classify_hsps(hsps, len(region_seq), exons or [], cfg)
    hsps, n_dropped = ots.apply_post_filter(hsps, cfg)
    dups = ots.duplicates(hsps)
    identical = ots.identical_segments(hsps, cfg)
    n_self = sum(1 for h in hsps if h["klass"] == "self")
    n_tx = sum(1 for h in hsps if h["klass"] == "transcript")
    _report(reporter, "Region screen ({} BLAT hits in {:.1f}s): {} duplicated segment(s) "
                      "({} identical), {} self, {} transcript, {} below thresholds".format(
                          len(records), elapsed, len(dups), len(identical), n_self,
                          n_tx, n_dropped),
            stage="offtarget_region")
    if identical:
        _report(reporter, "{} identical duplicated segment(s) >= {} bp — every oligo in "
                          "those windows will act on both copies".format(
                              len(identical), cfg["identical"]["identical_min_len"]),
                stage="offtarget_region", level="warning")

    return {
        "database":     genome,
        "is_fallback":  False,
        "duplicates":   dups,
        "identical":    identical,
        "excluded":     [h for h in hsps if h["klass"] in ("self", "transcript")],
        # Screen B needs these to drop each guide's own on-target site. They are
        # keyed by UCSC chromosome, so they only line up with a Screen B that
        # searched the same assembly — see design_tag_reagents._run_spacer_screen.
        "self_spans":   ots.self_spans(hsps),
        "own_accessions": sorted(ots.own_locus_accessions(hsps)),
        "n_below_threshold": n_dropped,
    }


def main(genomic_fasta, taxid, config=None, report=None):
    """CLI entry point: screen one genomic FASTA and print a summary."""
    cfg = ots.load_config(config)
    seq = ""
    with open(genomic_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq += line.strip()
    result = run_region_screen(seq, "", taxid, cfg=cfg, report=report)
    print("assembly:        {}".format(result["database"]))
    print("duplicates:      {}".format(len(result["duplicates"])))
    print("identical:       {}".format(len(result["identical"])))
    print("excluded:        {}".format(len(result["excluded"])))
    print("below threshold: {}".format(result["n_below_threshold"]))
    return result


def _parse_args(argv=None):
    """Command-line arguments for standalone use."""
    p = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    p.add_argument("--genomic-fasta", required=True, help="region FASTA (one record)")
    p.add_argument("--taxid", required=True, help="NCBI taxid of the source species")
    p.add_argument("--config", default=None, help="alternate offtarget.config.json")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()
    main(args.genomic_fasta, args.taxid, config=args.config)
