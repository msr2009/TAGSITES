"""
offtarget_screen.py

Pure, network-free core of the off-target / primer-specificity screen.

Two BLAST screens feed this module (submitted by offtarget_remote.py):

  Screen A — the whole genomic region as one query, finding *duplicated segments*
             elsewhere in the same species. Catches primer co-amplification (the
             collagen / gene-family case) and guides sitting in repeated sequence.
  Screen B — every spacer+PAM concatenated into one query separated by N-runs,
             finding *scattered* guide near-matches that have no regional homology.

Everything here is annotation only: nothing in this module filters or reorders
guides. The key trick is that a BLAST HSP carries the aligned subject sequence, so
an oligo's mismatches AND an off-target's PAM can both be read straight out of the
alignment — no per-hit sequence refetch is needed.

Coordinate conventions:
  - region coordinates are 0-based, half-open [start, end), matching crispr_util
  - BLAST query/subject coordinates from EBI JSON are 1-based inclusive
  - an oligo's 3' end is named explicitly ('right' or 'left' region coordinate)
    rather than derived from a strand, because primers and minus-strand guides
    flip that sense independently and deriving it twice invites sign errors

Matt Rich, 2025
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from crispr_util import iupac_to_regex, reverse_complement

CONFIG_PATH = Path(__file__).parent.parent / "offtarget.config.json"


def load_config(path=None):
    """Load offtarget.config.json (or a given path)."""
    with open(path or CONFIG_PATH) as f:
        return json.load(f)


# ── BLAST JSON parsing ────────────────────────────────────────────────────────

def _hit_strand(hsp, h_from, h_to):
    """Subject strand of an HSP, from hsp_strand when EBI supplies it."""
    raw = str(hsp.get("hsp_strand", "") or "")
    if "/" in raw:
        return "-" if raw.split("/")[-1].strip().lower().startswith("minus") else "+"
    return "-" if h_to < h_from else "+"


def parse_blast_json(payload):
    """Flatten EBI ncbiblast JSON into one dict per HSP (see module docstring for coords)."""
    if isinstance(payload, (bytes, bytearray)):
        payload = payload.decode("utf-8", "replace")
    if isinstance(payload, str):
        payload = json.loads(payload)
    hsps = []
    for hit in payload.get("hits", []) or []:
        for h in hit.get("hit_hsps", []) or []:
            q_from, q_to = int(h["hsp_query_from"]), int(h["hsp_query_to"])
            h_from, h_to = int(h["hsp_hit_from"]), int(h["hsp_hit_to"])
            hsps.append({
                "acc":       hit.get("hit_acc", ""),
                "desc":      hit.get("hit_desc", "") or "",
                "os":        hit.get("hit_os", "") or "",
                "q_from":    min(q_from, q_to),
                "q_to":      max(q_from, q_to),
                "h_from":    h_from,
                "h_to":      h_to,
                # EBI reports this as hsp_strand ("plus/plus", "plus/minus"); fall back
                # to coordinate order, which is unambiguous, when it is absent
                "h_strand":  _hit_strand(h, h_from, h_to),
                "qseq":      h.get("hsp_qseq", "") or "",
                "hseq":      h.get("hsp_hseq", "") or "",
                "identity":  float(h.get("hsp_identity", 0) or 0),
                "bits":      float(h.get("hsp_bit_score", 0) or 0),
                "evalue":    float(h.get("hsp_expect", 0) or 0),
                "align_len": int(h.get("hsp_align_len", 0) or 0)
                             or len(h.get("hsp_qseq", "") or ""),
            })
    return hsps


# ── Classification: self / transcript / duplicate ─────────────────────────────

def _looks_like_transcript(desc, keywords):
    """True when a hit description names an mRNA/cDNA/transcript record."""
    low = desc.lower()
    return any(k in low for k in keywords)


def _query_blocks(hsp):
    """Ungapped query intervals covered by an HSP, as 0-based half-open [start, end)."""
    blocks = []
    q = hsp["q_from"] - 1          # to 0-based
    run_start = None
    for qc in hsp["qseq"]:
        if qc == "-":
            # gap in the query: close any open run, subject coordinate advances only
            if run_start is not None:
                blocks.append((run_start, q))
                run_start = None
            continue
        if run_start is None:
            run_start = q
        q += 1
    if run_start is not None:
        blocks.append((run_start, q))
    return blocks


def _overlap(a, b):
    """Length of the overlap between two half-open intervals."""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def _exon_coincidence(blocks, exons):
    """Fraction of a query footprint that falls inside annotated exons."""
    total = sum(e - s for s, e in blocks)
    if not total or not exons:
        return 0.0
    inside = sum(_overlap(b, e) for b in blocks for e in exons)
    return inside / total


def _skips_introns(blocks, exons):
    """True when a footprint covers >= 2 exons while leaving the introns between them bare.

    This is the signature of a transcript: an mRNA record aligns to exon after exon
    and never to the intervening intronic sequence. It deliberately requires two
    exons, because a hit confined to a *single* exon is evidence of nothing — a
    paralogous gene's coding exon looks exactly the same, which is precisely the
    gene-family case this screen exists to find.
    """
    if len(exons) < 2 or not blocks:
        return False
    ordered = sorted(exons)
    touched = [i for i, e in enumerate(ordered)
               if any(_overlap(b, e) > 0 for b in blocks)]
    if len(touched) < 2:
        return False
    # every intron lying between the first and last touched exon must be uncovered
    for i in range(touched[0], touched[-1]):
        intron = (ordered[i][1], ordered[i + 1][0])
        if intron[1] <= intron[0]:
            continue
        if any(_overlap(b, intron) > 0 for b in blocks):
            return False
    return True


def classify_hsps(hsps, query_len, exons=None, cfg=None):
    """Tag each HSP 'self', 'transcript' or 'duplicate'; returns the same list.

    Classification is per *accession*, not per HSP: BLAST splits an mRNA hit into
    one HSP per exon, so no single HSP can show that a record skips introns.
    """
    cfg = cfg or load_config()
    sh, tr = cfg["self_hit"], cfg["transcripts"]
    exons = sorted(exons or [])

    by_acc = {}
    for h in hsps:
        by_acc.setdefault(h["acc"], []).append(h)

    for acc, group in by_acc.items():
        blocks = [b for h in group for b in _query_blocks(h)]
        covered = sum(e - s for s, e in _merge(blocks))
        coverage = 100.0 * covered / max(1, query_len)
        best_identity = max(h["identity"] for h in group)
        desc = next((h["desc"] for h in group if h["desc"]), "")

        klass = "duplicate"
        # The query IS a slice of this genome, so its own locus returns a
        # near-full-length near-identical hit. Recognised without needing the
        # region's provenance, since a user may paste their own FASTA.
        if best_identity >= sh["self_identity"] and coverage >= sh["self_coverage"]:
            klass = "self"
        elif tr["exclude_transcripts"]:
            keyword_hit = tr["use_keywords"] and _looks_like_transcript(desc, tr["keywords"])
            # Structural evidence must be positive proof of intron skipping AND
            # near-identity: our own transcript is ~100% identical to our exons,
            # whereas a paralog is not.
            structural_hit = (
                _exon_coincidence(blocks, exons) >= tr["exon_coincidence"]
                and _skips_introns(blocks, exons)
                and best_identity >= tr.get("transcript_identity", 99.0)
            )
            if keyword_hit or structural_hit:
                klass = "transcript"
        for h in group:
            h["klass"] = klass
    return hsps


def _merge(intervals):
    """Merge overlapping half-open intervals."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def apply_post_filter(hsps, cfg=None):
    """Drop low-scoring 'duplicate' HSPs; returns (kept, n_dropped).

    Screen A runs BLAST at a permissive E-value so that short perfectly-identical
    segments survive, then removes noise here where the thresholds are ours.
    """
    cfg = cfg or load_config()
    pf = cfg["post_filter"]
    kept, dropped = [], 0
    for h in hsps:
        if h.get("klass") != "duplicate":
            kept.append(h)
            continue
        if (h["identity"] >= pf["min_identity"]
                and h["align_len"] >= pf["min_align_len"]
                and h["bits"] >= pf["min_bitscore"]):
            kept.append(h)
        else:
            dropped += 1
    return kept, dropped


def duplicates(hsps):
    """The HSPs that represent a real duplicated segment elsewhere in the genome."""
    return [h for h in hsps if h.get("klass") == "duplicate"]


def self_spans(hsps, pad=0):
    """Subject coordinates of the query's own locus, per accession.

    Screen A identifies these as a side effect of recognising self-hits, and Screen B
    needs them: a spacer trivially matches its own on-target site, and ENA carries
    several independent submissions of the same genome (six for C. elegans), so
    without this every guide reports roughly one spurious hit per assembly.
    """
    spans = {}
    for h in hsps:
        if h.get("klass") != "self":
            continue
        lo, hi = min(h["h_from"], h["h_to"]), max(h["h_from"], h["h_to"])
        spans.setdefault(h["acc"], []).append((lo - pad, hi + pad))
    return {acc: _merge(v) for acc, v in spans.items()}


def in_self_span(acc, span, spans):
    """True when a subject span lies inside the query's own locus on that accession."""
    if not span or not spans:
        return False
    return any(span[0] >= lo and span[1] <= hi for lo, hi in spans.get(acc, []))


def collapse_loci(sites):
    """Collapse sites that are the same locus reported by different submissions.

    ENA holds several independent assemblies per species (six for C. elegans), so a
    single off-target is reported once per assembly — measured on snt-1, one real
    site appeared as 5 hits. Two sites with the same matched subject sequence and
    PAM are treated as one locus, and the contributing accessions are kept on the
    representative so nothing is silently discarded.

    Two genuinely distinct loci with byte-identical sequence would merge; that is
    the accepted cost, and such a case is reported anyway by the region screen's
    identical-segment warning.
    """
    by_sig = {}
    for s in sites:
        sig = (s.get("subject_seq", ""), s.get("pam"), s.get("mismatches"), s.get("gaps"))
        if sig in by_sig:
            by_sig[sig]["accessions"].append(s["acc"])
            continue
        rep = dict(s)
        rep["accessions"] = [s["acc"]]
        by_sig[sig] = rep
    return list(by_sig.values())


def own_locus_accessions(hsps):
    """Accessions that are records OF the query's own locus, not other loci.

    Transcript records of our own gene must be excluded wholesale rather than by
    coordinate: an mRNA is our own exonic sequence, so every spacer in an exon
    matches it. Measured on snt-1, one mRNA record (L15302.1) accounted for 87 of
    110 apparent spacer off-targets.
    """
    return {h["acc"] for h in hsps if h.get("klass") in ("self", "transcript")}


# ── Alignment window scoring ──────────────────────────────────────────────────

def _column_index(hsp):
    """Map region coordinate -> alignment column for every ungapped query base."""
    index = {}
    q = hsp["q_from"] - 1
    for col, qc in enumerate(hsp["qseq"]):
        if qc == "-":
            continue
        index[q] = col
        q += 1
    return index


def score_window(hsp, start, end, three_prime_side, cfg=None):
    """Score one oligo footprint against an HSP alignment.

    start/end are 0-based half-open region coordinates. three_prime_side names
    which end of that interval is the oligo's 3' terminus ('right' = end-1,
    'left' = start) — passed explicitly because forward primers, reverse primers
    and minus-strand guides each place it differently.

    Returns None when the HSP does not fully cover the oligo, else a dict with
    mismatch/gap counts, the 3'-end verdict, and whether it is a perfect match.
    """
    cfg = cfg or load_config()
    sc = cfg["scoring"]
    index = _column_index(hsp)
    cols = [index.get(p) for p in range(start, end)]
    # Partial coverage cannot be scored honestly — a missing column is not a match
    if any(c is None for c in cols):
        return None

    qs, hs = hsp["qseq"], hsp["hseq"]
    mismatches = gaps = 0
    per_base = []          # True where the aligned subject base matches, 5'->3' in region order
    for c in cols:
        qc, hc = qs[c].upper(), hs[c].upper()
        if hc == "-":
            gaps += 1
            per_base.append(False)
        elif qc != hc:
            mismatches += 1
            per_base.append(False)
        else:
            per_base.append(True)

    # Orient 5'->3' along the oligo before taking its 3'-terminal window
    oriented = per_base if three_prime_side == "right" else per_base[::-1]
    tp = int(sc["three_prime_len"])
    three_prime_ok = all(oriented[-tp:]) if len(oriented) >= tp else all(oriented)

    total = mismatches + gaps
    return {
        "acc":             hsp["acc"],
        "desc":            hsp["desc"],
        "mismatches":      mismatches,
        "gaps":            gaps,
        "three_prime_ok":  three_prime_ok,
        "perfect":         total == 0,
        "flagged":         three_prime_ok and total <= int(sc["max_mismatch"]),
        "subject_span":    _subject_span(hsp, start, end),
        "h_strand":        hsp["h_strand"],
        # the matched subject bases; used to collapse the same locus reported from
        # several independent genome submissions of the same species
        "subject_seq":     "".join(hs[c] for c in cols).upper(),
    }


def _subject_span(hsp, start, end):
    """Subject coordinates aligned to a region interval, as (low, high) 1-based."""
    index = _column_index(hsp)
    cols = [index.get(p) for p in range(start, end) if index.get(p) is not None]
    if not cols:
        return None
    # Walk the subject once, counting non-gap characters up to each column
    hs = hsp["hseq"]
    step = 1 if hsp["h_strand"] == "+" else -1
    pos = hsp["h_from"]
    col_to_subject = {}
    for col, hc in enumerate(hs):
        if hc == "-":
            continue
        col_to_subject[col] = pos
        pos += step
    mapped = [col_to_subject[c] for c in cols if c in col_to_subject]
    if not mapped:
        return None
    return (min(mapped), max(mapped))


# ── Guides ────────────────────────────────────────────────────────────────────

def guide_footprint(row, guide_length=20, pam_len=3):
    """Region coords of a guide's spacer and PAM, plus which end is its 3' terminus.

    Works from the columns the reagents TSV actually stores (guide_strand,
    pam_fwd_start): find_guides' guide_fwd_start/end are not persisted but are
    fully determined by these.
    """
    pam_start = int(row["pam_fwd_start"])
    if str(row["guide_strand"]) == "+":
        spacer = (pam_start - guide_length, pam_start)
        three_prime_side = "right"          # spacer 3' end abuts the PAM on its right
    else:
        spacer = (pam_start + pam_len, pam_start + pam_len + guide_length)
        three_prime_side = "left"           # PAM lies to the left in forward coordinates
    return spacer, (pam_start, pam_start + pam_len), three_prime_side


def offtarget_pam(hsp, pam_span, strand, pam="NGG"):
    """Read the off-target's PAM out of the alignment; None when not covered.

    hseq is already oriented to pair with qseq, so the subject bases opposite our
    PAM are the off-target's PAM in the same relative sense. For a minus-strand
    guide the forward-coordinate bases must be reverse-complemented to read the
    PAM on the guide strand.
    """
    index = _column_index(hsp)
    cols = [index.get(p) for p in range(*pam_span)]
    if any(c is None for c in cols):
        return None
    bases = "".join(hsp["hseq"][c] for c in cols).upper()
    if "-" in bases:
        return None
    if strand == "-":
        bases = reverse_complement(bases)
    return bases


def pam_matches(bases, pam="NGG"):
    """True when PAM bases satisfy the IUPAC PAM pattern."""
    return bases is not None and re.fullmatch(iupac_to_regex(pam), bases) is not None


def screen_guide(hsps, row, guide_length=20, pam="NGG", cfg=None):
    """Screen one guide row against duplicated segments (Screen A).

    Returns {'sites': [...], 'n_total', 'n_identical', 'n_pam_unverified'}. A site
    is reported only when the spacer's 3' end is intact and mismatches are within
    threshold; whether the off-target also carries a PAM is reported per site
    rather than assumed.
    """
    cfg = cfg or load_config()
    spacer_span, pam_span, tp_side = guide_footprint(row, guide_length, len(pam))
    strand = str(row["guide_strand"])
    sites = []
    for h in duplicates(hsps):
        hit = score_window(h, spacer_span[0], spacer_span[1], tp_side, cfg)
        if not hit or not hit["flagged"]:
            continue
        bases = offtarget_pam(h, pam_span, strand, pam)
        hit["pam"] = bases
        hit["pam_ok"] = pam_matches(bases, pam)
        hit["pam_unverified"] = bases is None
        hit["screen"] = "region"
        # A second Cas9 site needs its own PAM; without one the near-match cannot cut,
        # so only PAM-bearing or PAM-unknown hits are reported as sites
        if hit["pam_ok"] or hit["pam_unverified"]:
            sites.append(hit)
    sites = collapse_loci(sites)
    return {
        "sites":            sites,
        "n_total":          len(sites),
        "n_identical":      sum(1 for s in sites if s["perfect"] and s["pam_ok"]),
        "n_pam_unverified": sum(1 for s in sites if s["pam_unverified"]),
    }


# ── Screen B: concatenated spacer query ───────────────────────────────────────

def build_spacer_query(spacers, separator_len=25):
    """Join spacer+PAM blocks with N-runs into one query; returns (seq, block_len)."""
    if not spacers:
        return "", 0
    block_len = max(len(s) for s in spacers)
    sep = "N" * int(separator_len)
    # Pad short blocks so every block starts at a predictable offset
    padded = [s.upper().ljust(block_len, "N") for s in spacers]
    return sep.join(padded), block_len


def block_of(query_pos, block_len, separator_len, n_blocks):
    """Which concatenated block a 1-based query coordinate falls in, or None.

    Returns (index, offset_within_block). A coordinate inside a separator, or an
    HSP straddling one, yields None — blastn will not gap-extend across an N-run,
    so this should not occur and is treated as unmappable rather than guessed.
    """
    stride = block_len + int(separator_len)
    p = int(query_pos) - 1
    idx, off = divmod(p, stride)
    if idx >= n_blocks or off >= block_len:
        return None
    return idx, off


def screen_spacer_hits(hsps, spacers, block_len, separator_len, pam="NGG", cfg=None,
                       exclude_spans=None, exclude_accessions=None):
    """Map Screen B HSPs back to spacers; returns {spacer_index: [sites]}.

    The query here is spacer+PAM, so the PAM columns are part of the alignment and
    an off-target's PAM is verified the same way as in Screen A.

    exclude_spans (from self_spans()) removes each guide's own on-target site. It is
    not optional in practice: without it every guide reports a perfect PAM-bearing
    "off-target" for every copy of the source genome in the database.
    """
    cfg = cfg or load_config()
    n = len(spacers)
    by_spacer = {}
    n_self = 0
    for h in hsps:
        start = block_of(h["q_from"], block_len, separator_len, n)
        end = block_of(h["q_to"], block_len, separator_len, n)
        # Must fall wholly inside one block; anything else is unmappable
        if start is None or end is None or start[0] != end[0]:
            continue
        idx = start[0]
        spacer_len = len(spacers[idx]) - len(pam)
        # Re-anchor the HSP onto block-local coordinates so score_window can be reused
        local = dict(h)
        local["q_from"] = start[1] + 1
        local["q_to"] = end[1] + 1
        local["klass"] = "duplicate"
        hit = score_window(local, 0, spacer_len, "right", cfg)
        if not hit or not hit["flagged"]:
            continue
        bases = offtarget_pam(local, (spacer_len, spacer_len + len(pam)), "+", pam)
        hit["pam"] = bases
        hit["pam_ok"] = pam_matches(bases, pam)
        hit["pam_unverified"] = bases is None
        hit["screen"] = "spacer"
        # drop the guide's own on-target site, in whichever assembly it was matched,
        # and any record that is itself our locus (our gene's own mRNA entries)
        if (h["acc"] in (exclude_accessions or set())
                or in_self_span(h["acc"], hit["subject_span"], exclude_spans)):
            n_self += 1
            continue
        if hit["pam_ok"] or hit["pam_unverified"]:
            by_spacer.setdefault(idx, []).append(hit)
    # one entry per distinct locus, not per database record
    return {k: collapse_loci(v) for k, v in by_spacer.items()}, n_self


# ── Primers: predicted spurious amplicons ─────────────────────────────────────

def predict_amplicons(hsps, fwd_span, rev_span, cfg=None):
    """Predict spurious products where BOTH primers bind the same duplicated segment.

    A lone primer in a duplicate is harmless; a product needs both primers in the
    same duplicate, correctly oriented, within max_product. Returns a list of
    {acc, desc, product_size, fwd, rev}.
    """
    cfg = cfg or load_config()
    max_product = int(cfg["scoring"]["max_product"])
    out = []
    for h in duplicates(hsps):
        # Forward primer's 3' end points right; the reverse primer's points left
        f = score_window(h, fwd_span[0], fwd_span[1], "right", cfg)
        r = score_window(h, rev_span[0], rev_span[1], "left", cfg)
        if not f or not r or not f["flagged"] or not r["flagged"]:
            continue
        if not f["subject_span"] or not r["subject_span"]:
            continue
        lo = min(f["subject_span"][0], r["subject_span"][0])
        hi = max(f["subject_span"][1], r["subject_span"][1])
        size = hi - lo + 1
        if size > max_product:
            continue
        out.append({
            "acc":          h["acc"],
            "desc":         h["desc"],
            "product_size": size,
            "fwd":          f,
            "rev":          r,
            "perfect":      f["perfect"] and r["perfect"],
        })
    return out


# ── Region-level summary ──────────────────────────────────────────────────────

def identical_segments(hsps, cfg=None):
    """Non-self HSPs that are 100% identical over a long stretch (segmental duplication).

    Every oligo inside such a window is suspect regardless of its own score, so
    this is surfaced per site rather than per oligo.
    """
    cfg = cfg or load_config()
    min_len = int(cfg["identical"]["identical_min_len"])
    return [h for h in duplicates(hsps)
            if h["identity"] >= 99.999 and h["align_len"] >= min_len]


def summarise(n_total, n_identical, n_pam_unverified):
    """Compact human-readable detail string for the reagents TSV / UI."""
    if not n_total:
        return "no sites found"
    parts = []
    if n_identical:
        parts.append("{} identical".format(n_identical))
    near = n_total - n_identical - n_pam_unverified
    if near > 0:
        parts.append("{} near".format(near))
    if n_pam_unverified:
        parts.append("{} PAM unverified".format(n_pam_unverified))
    return "{} site{} ({})".format(n_total, "" if n_total == 1 else "s", ", ".join(parts))
