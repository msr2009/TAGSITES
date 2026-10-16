"""ucsc_rest.py — gene symbol -> flanked genomic FASTA via the UCSC Genome Browser REST API.

Second source in fetch_genomic_sequence's fallback chain. Keyless (only BLAT needs
the apiKey), and covers the assemblies in offtarget.config.json's blat.assemblies
map, which is every preset organism except E. coli.

Three calls per lookup, all through http_retry with the caller's deadline:
  /search          symbol -> candidate positions (1-based inclusive "chrom:start-end")
  /getData/track   the gene track over the chosen span, for its strand (search omits it)
  /getData/sequence  the flanked region (0-based half-open), via offtarget_blat

/search is a fuzzy full-text search, so only matches whose posName equals the
symbol (case-insensitive) on a gene track are kept. UCSC asks for roughly one
request per second; REQUEST_PAUSE_S spaces the three calls.

Matt Rich, 2026
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import http_retry
import offtarget_blat
from crispr_util import reverse_complement

SEARCH_URL = "https://api.genome.ucsc.edu/search"
TRACK_URL = "https://api.genome.ucsc.edu/getData/track"

# gene tracks in preference order: refGene carries the symbol as posName on every
# preset assembly; ensGene is the only one that does on sacCer3
GENE_TRACKS = ("refGene", "ncbiRefSeq", "ensGene", "sgdGene", "knownGene")

REQUEST_PAUSE_S = 1.0


def _get(url, params, cfg, deadline):
    """GET a UCSC endpoint through the shared retry/deadline wrapper; return parsed JSON."""
    t = cfg["timeouts"]
    ret = cfg["retries"]
    resp = http_retry.request_with_retries(
        "get", url, retries=ret["retries"], retry_wait=ret["retry_wait_s"], deadline=deadline,
        params=params, timeout=(t["connect_timeout_s"], t["read_timeout_s"]),
    )
    payload = resp.json()
    if payload.get("error"):
        raise RuntimeError(f"UCSC {url.rsplit('/', 1)[-1]}: {payload['error']}")
    return payload


def _pause(deadline):
    """Sleep REQUEST_PAUSE_S between UCSC calls, never past the deadline."""
    left = http_retry.remaining(deadline)
    time.sleep(REQUEST_PAUSE_S if left is None else min(REQUEST_PAUSE_S, left))


def _parse_position(position):
    """'chrX:7682896-7686200' -> ('chrX', 7682896, 7686200), 1-based inclusive."""
    chrom, rng = position.rsplit(":", 1)
    start, end = rng.replace(",", "").split("-")
    return chrom, int(start), int(end)


def _cluster_spans(spans):
    """Merge overlapping (chrom, start, end) transcript spans into gene-level spans."""
    out = []
    # sort so overlapping transcripts on one chromosome are adjacent
    for chrom, start, end in sorted(spans):
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            out[-1][2] = max(out[-1][2], end)
        else:
            out.append([chrom, start, end])
    return [tuple(s) for s in out]


def search_gene(genome, symbol, cfg, deadline=None):
    """Return (track, [gene spans]) for exact symbol matches, or (None, []) when none."""
    payload = _get(SEARCH_URL, {"search": symbol, "genome": genome}, cfg, deadline)
    by_track = {}
    for group in payload.get("positionMatches") or []:
        track = group.get("trackName")
        if track not in GENE_TRACKS:
            continue
        for m in group.get("matches") or []:
            if str(m.get("posName", "")).lower() == symbol.lower() and m.get("position"):
                by_track.setdefault(track, []).append(_parse_position(m["position"]))
    # first preferred track with any exact hit
    for track in GENE_TRACKS:
        if by_track.get(track):
            return track, _cluster_spans(by_track[track])
    return None, []


def gene_strand(genome, track, symbol, chrom, start, end, cfg, deadline=None):
    """Strand ('+'/'-') of the symbol's gene-track items in a 1-based span, or None if unclear."""
    payload = _get(TRACK_URL, {"genome": genome, "track": track, "chrom": chrom,
                               "start": start - 1, "end": end}, cfg, deadline)
    items = payload.get(track) or []
    # a dict keyed by chrom is returned for some tracks; flatten it
    if isinstance(items, dict):
        items = [it for chunk in items.values() for it in chunk]
    inside = [it for it in items
              if it.get("txStart", -1) >= start - 1 and it.get("txEnd", 0) <= end]
    # prefer items named for the symbol (refGene/ncbiRefSeq name2), so a nested
    # gene on the other strand can't make the call ambiguous; ensGene's name2 is
    # an id, so fall back to every item contained in the span
    named = [it for it in inside if str(it.get("name2", "")).lower() == symbol.lower()]
    strands = {it.get("strand") for it in (named or inside)}
    strands.discard(None)
    if len(strands) == 1:
        return strands.pop()
    return None


def fetch_region(genome, chrom, start, end, strand, flank, deadline=None):
    """Flanked region as a sequence string, in the gene's orientation; returns (seq, start, end).

    start/end are 1-based inclusive gene coordinates; the returned pair is the
    flanked region in the same convention, clamped at the chromosome start.
    """
    _pause(deadline)
    lo = max(1, start - flank)
    hi = end + flank
    seq = offtarget_blat.fetch_sequence(genome, chrom, lo - 1, hi, deadline=deadline)
    if not seq:
        raise RuntimeError(f"UCSC returned no sequence for {genome} {chrom}:{lo}-{hi}")
    if strand == "-":
        seq = reverse_complement(seq)
    return seq, lo, lo + len(seq) - 1


def find_gene(genome, symbol, cfg, deadline=None):
    """Resolve a symbol to the longest exact-match gene span; returns a coords dict or None.

    coords: chrom, start, end (1-based inclusive), strand ('+'/'-'), track,
    candidates (all spans as 'chrom:start-end', longest first).
    """
    track, spans = search_gene(genome, symbol, cfg, deadline=deadline)
    if not spans:
        return None
    # longest genomic span wins, as in fetch_genomic_sequence._pick_most_complete_gene
    spans.sort(key=lambda s: (-(s[2] - s[1]), s))
    chrom, start, end = spans[0]
    _pause(deadline)
    strand = gene_strand(genome, track, symbol, chrom, start, end, cfg, deadline=deadline)
    if strand is None:
        raise RuntimeError(f"UCSC {track}: could not determine strand for "
                           f"{symbol} at {chrom}:{start}-{end}")
    return {
        "chrom": chrom, "start": start, "end": end, "strand": strand, "track": track,
        "candidates": [f"{c}:{s}-{e}" for c, s, e in spans],
    }
