"""
plasmid_assembly.py

SapTrap donor-plasmid reagent design.  Stateless, no Shiny imports.
Dual-use: imported by reagents_server.py and runnable from the CLI.

Produces:
  - Order-sheet rows: annealed-oligo pairs (arm ≤ OLIGO_MAX nt) or
    PCR amplification primers with SapI sites (arm > OLIGO_MAX nt).
  - Optional linear GenBank records (one per fragment) formatted for ApE's
    Golden Gate assembler.

SapTrap connector overhangs verified from:
  Schwartz & Jorgensen 2016, Genetics 202:1277 (PMID 26837755), Fig 2A
  Dickinson et al. 2018, microPublication Biology (PMID 32550377), Fig 1A

Matt Rich, 2025
"""

import io
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from crispr_util import reverse_complement, CODONS, SYN_CODONS
from reagent_sequences import (
    calc_tm,
    _primer3_anchored_single,
    _primer3_single,
    _tm_growth_bind_regions,
)

# ── constants ──────────────────────────────────────────────────────────────────

SAPI_SITE = "GCTCTTC"   # recognition sequence (5'→3')
SAPI_RC   = "GAAGAGC"   # reverse complement of recognition (for bottom-strand site)
SAPI_SPR  = "A"         # 1-nt spacer between recognition site and overhang in primers

# Maximum total oligo length for ordering as an annealed pair vs. PCR amplicon.
# IDT / Twist routinely synthesize up to ~90 nt single-stranded oligos;
# 87 nt leaves room for the 3-nt overhang tail.
OLIGO_MAX = 87

# SapI connector overhangs (3-nt, top strand 5'→3' at each junction).
# Both Schwartz-2016 and Dickinson-SEC use this identical set for the
# site-specific fragments (sgRNA, 5'HA, 3'HA).
# Format: (left_overhang, right_overhang) reading 5'→3' along the donor.
SAPTRAP_OVERHANGS = {
    "sgrna": ("TTG", "GTT"),   # spacer → U6/scaffold junctions
    "ha5":   ("TGG", "GCG"),   # vector → 5'HA; 5'HA → C-tag connector
    "ha3":   ("ACG", "GTA"),   # N-tag connector → 3'HA; 3'HA → vector
}

# Per-fragment colors for ApE visualization (ApEinfo hex format)
FRAG_COLORS = {
    "sgrna": "#4a90d9",   # blue
    "ha5":   "#2ecc71",   # green
    "ha3":   "#e67e22",   # orange
}
FRAG_LABELS = {"sgrna": "sgRNA", "ha5": "5' HA", "ha3": "3' HA"}

LINKER_SEQ = "nnn"   # separates fragments in the combined GenBank record


# ── low-level helpers ──────────────────────────────────────────────────────────

def _rc(seq):
    """Reverse complement (uppercase)."""
    return reverse_complement(seq.upper())


def find_internal_sapi(seq):
    """Return positions (0-indexed) of internal SapI sites (both strands).

    Reports hits for GCTCTTC (+ strand) and GAAGAGC (- strand = RC recognition).
    Positions are start indices of the recognition sequence.
    """
    seq_up = seq.upper()
    rc_site = SAPI_RC
    hits = []
    # forward strand
    i = 0
    while True:
        p = seq_up.find(SAPI_SITE, i)
        if p == -1:
            break
        hits.append((p, "+"))
        i = p + 1
    # reverse strand
    i = 0
    while True:
        p = seq_up.find(rc_site, i)
        if p == -1:
            break
        hits.append((p, "-"))
        i = p + 1
    return hits


# ── SapI site domestication ────────────────────────────────────────────────────

# GC-matched swap: change G↔C or A↔T (preserves GC count)
_GC_SWAP = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A'}


def _arm_codon_lookup(arm, is_left_arm):
    """Map each exonic (uppercase) position in arm to its codon.

    Uses the same case-encoding convention as the reagents TSV:
      uppercase = exonic/coding, lowercase = intronic/non-coding.

    Reading direction follows the arm orientation relative to the insertion site:
      is_left_arm=True  → codons grouped from the RIGHT end (insertion-proximal)
      is_left_arm=False → codons grouped from the LEFT end (insertion-proximal)

    Returns dict: {pos: (sorted_codon_positions, codon_str_5prime_to_3prime)}
    Incomplete triplets at the distal end are omitted.
    """
    ups = [i for i in range(len(arm)) if arm[i].isupper()]
    n = len(ups)
    lookup = {}

    if is_left_arm:
        # Group from the right: ups[-3:] → codon 0 (most proximal), etc.
        for codon_num in range(n // 3):
            end_idx   = n - 3 * codon_num
            start_idx = end_idx - 3
            cod_pos   = ups[start_idx:end_idx]   # left-to-right (5'→3')
            cod_str   = ''.join(arm[p].upper() for p in cod_pos)
            for p in cod_pos:
                lookup[p] = (cod_pos, cod_str)
    else:
        # Group from the left: ups[:3] → codon 0 (most proximal), etc.
        for codon_num in range(n // 3):
            start_idx = 3 * codon_num
            end_idx   = start_idx + 3
            cod_pos   = ups[start_idx:end_idx]
            cod_str   = ''.join(arm[p].upper() for p in cod_pos)
            for p in cod_pos:
                lookup[p] = (cod_pos, cod_str)

    return lookup


def domesticate_sapi(arm_seq, is_left_arm=True):
    """Mutate all internal SapI recognition sites in arm_seq.

    Applies the same three-tier ladder used for PAM disruption:
      1. Synonymous codon change (coding/uppercase position) that disrupts the site.
      2. GC-matched single-base swap at a non-coding (lowercase) position.
      3. GC-matched single-base swap at a coding (uppercase) position (non-synonymous).
      Fallback: any single-base change that disrupts the site.

    arm_seq uses uppercase for exonic bases, lowercase for intronic — same as the
    left_arm/right_arm columns of the .reagents.tsv.

    is_left_arm: True for the 5′ HA (insertion at the right end of the arm),
                 False for the 3′ HA (insertion at the left end).

    Returns (new_arm_seq, changes) where changes is a list of
    (position, original_char, new_char, method_str).
    """
    arm = list(arm_seq)
    changes = []

    for _guard in range(20):   # at most 20 iterations (one per hit)
        hits = find_internal_sapi(''.join(arm))
        if not hits:
            break
        hit_pos, _ = hits[0]
        site_pos = list(range(hit_pos, hit_pos + 7))

        codon_lu = _arm_codon_lookup(''.join(arm), is_left_arm)

        def _disrupted(test_arm):
            """Check that the 7-mer at hit_pos is no longer a SapI site."""
            s = ''.join(test_arm[hit_pos:hit_pos + 7]).upper()
            return s != SAPI_SITE and s != SAPI_RC

        applied = False

        # ── 1. Synonymous codon change ────────────────────────────────────────
        seen = {}
        for p in site_pos:
            if p in codon_lu:
                cod_pos, cod_str = codon_lu[p]
                key = tuple(cod_pos)
                if key not in seen:
                    seen[key] = (cod_pos, cod_str)

        for cod_pos, orig_codon in seen.values():
            aa = CODONS.get(orig_codon)
            if not aa or aa == '_':
                continue
            for alt in SYN_CODONS.get(aa, []):
                if alt == orig_codon:
                    continue
                # at least one differing base must be inside the SapI site
                if not any(cod_pos[k] in site_pos for k in range(3)
                           if orig_codon[k] != alt[k]):
                    continue
                test = arm[:]
                for k in range(3):
                    if orig_codon[k] != alt[k]:
                        test[cod_pos[k]] = alt[k]   # uppercase preserved
                if _disrupted(test):
                    for k in range(3):
                        if orig_codon[k] != alt[k]:
                            changes.append((cod_pos[k], arm[cod_pos[k]], alt[k], 'syn'))
                            arm[cod_pos[k]] = alt[k]
                    applied = True
                    break
            if applied:
                break
        if applied:
            continue

        # ── 2. GC-matched change at intronic (lowercase) position ─────────────
        for p in site_pos:
            if arm[p].islower():
                new_b = _GC_SWAP[arm[p].upper()].lower()
                if new_b != arm[p]:
                    test = arm[:]
                    test[p] = new_b
                    if _disrupted(test):
                        changes.append((p, arm[p], new_b, 'noncoding'))
                        arm[p] = new_b
                        applied = True
                        break
        if applied:
            continue

        # ── 3. GC-matched change at exonic (uppercase) position ───────────────
        for p in site_pos:
            if arm[p].isupper():
                new_b = _GC_SWAP[arm[p].upper()]
                if new_b != arm[p]:
                    test = arm[:]
                    test[p] = new_b
                    if _disrupted(test):
                        changes.append((p, arm[p], new_b, 'mut'))
                        arm[p] = new_b
                        applied = True
                        break
        if applied:
            continue

        # ── Fallback: any change that works ───────────────────────────────────
        for p in site_pos:
            for b in 'ACGT':
                new_b = b if arm[p].isupper() else b.lower()
                if new_b == arm[p]:
                    continue
                test = arm[:]
                test[p] = new_b
                if _disrupted(test):
                    changes.append((p, arm[p], new_b, 'forced'))
                    arm[p] = new_b
                    applied = True
                    break
            if applied:
                break
        if not applied:
            break   # no change found — shouldn't happen with valid DNA

    return ''.join(arm), changes


# ── annealed oligo pair ────────────────────────────────────────────────────────

def annealed_oligo_pair(seq, left_oh, right_oh):
    """Build a top/bottom annealed oligo pair with 3-nt 5′ overhang tails.

    After annealing the two oligos, the duplex presents:
      - 5′ overhang on top strand (left end): left_oh
      - 5′ overhang on bottom strand (right end): RC(right_oh)
    which are the two SapTrap connector overhangs for this fragment.

    top    = left_oh + seq                   (5'→3')
    bottom = RC(right_oh) + RC(seq)          (5'→3')
    """
    top = left_oh + seq
    bot = _rc(right_oh) + _rc(seq)
    return top, bot


def sgrna_oligos(spacer):
    """Annealed oligo pair for cloning a 20-nt spacer into a SapTrap sgRNA vector.

    Follows Schwartz & Jorgensen 2016 convention:
      top = TTG + spacer   (5'→3', presents 5'-TTG overhang)
      bot = AAC + RC(spacer)  (5'→3', presents 5'-AAC overhang at right end)

    Note: design_guides_across_region.py uses the reversed convention (aac/ttg)
    — this is the published-standard orientation.
    """
    left_oh, right_oh = SAPTRAP_OVERHANGS["sgrna"]   # TTG, GTT
    return annealed_oligo_pair(spacer, left_oh, right_oh)


# ── PCR amplification primers ──────────────────────────────────────────────────

def sapi_amplification_primers(arm_seq, left_oh, right_oh, tm_target=60.0,
                               nominal_length=None, tolerance=0.2,
                               pinned_side='rev'):
    """Design SapI-tailed PCR primers to amplify arm_seq from WT genomic template.

    arm_seq should be the MUTATED arm (from the reagents TSV) so that PAM-
    disrupting mutations land in the synthesized primer rather than requiring
    them to be present in the genomic template. When tolerance > 0, arm_seq is
    expected to be widened beyond nominal_length on the far/outer side only
    (see reagent_sequences.truncate_arms_with_tolerance) — the insertion-
    adjacent boundary must never move regardless of widening.

    pinned_side: 'rev' when the LAST base of arm_seq is insertion-adjacent
    (must be exact — left/5'HA arm case, far/outer end is position 0); 'fwd'
    when the FIRST base is insertion-adjacent (right/3'HA arm case, far/outer
    end is the last base). The pinned end is designed exactly at that fixed
    position via primer3; the opposite (flexible) end is searched within
    ±tolerance of nominal_length, letting primer3 pick a better-quality
    primer (Tm, hairpins, GC clamp) than a fixed boundary would allow.

    Primer structure (tail always present regardless of exact/flexible end):
      fwd: 5'-GCTCTTCA-{left_oh}-{binding_region}-3'
      rev: 5'-GCTCTTCA-{RC(right_oh)}-{RC(binding_region)}-3'

    After SapI digestion the amplicon presents:
      left end : 5'-{left_oh} overhang (top strand)
      right end: 5'-{RC(right_oh)} overhang (bottom strand) — mates with {right_oh}

    Falls back to the old hand-rolled Tm-growth loop, run at exactly
    nominal_length (not the widened arm_seq), if primer3 finds no valid
    design for either end.

    Returns ((fwd_suffix, fwd_seq), (rev_suffix, rev_seq), meta), where
    meta = {"effective_arm_length": int, "used_fallback": bool}.
    """
    arm = arm_seq.upper()
    tail = SAPI_SITE + SAPI_SPR  # GCTCTTCA
    n = len(arm)
    if nominal_length is None:
        nominal_length = n

    if pinned_side == 'rev':
        pinned_pos, pinned_pick_left = n - 1, False
        flex_pick_left = True
        lo = max(0, n - int(nominal_length * (1 + tolerance)))
        hi = max(lo, n - int(round(nominal_length * (1 - tolerance))))
    else:
        pinned_pos, pinned_pick_left = 0, True
        flex_pick_left = False
        hi = min(n - 1, int(nominal_length * (1 + tolerance)) - 1)
        lo = min(hi, int(round(nominal_length * (1 - tolerance))) - 1)
        lo = max(lo, 0)

    pinned_result = _primer3_anchored_single(arm, pinned_pos, pinned_pick_left, tm_target)

    # A too-narrow search window (e.g. tolerance=0, or a short nominal_length)
    # can't fit even a minimal primer — skip straight to fallback rather than
    # letting primer3 raise on an unsatisfiable SEQUENCE_INCLUDED_REGION.
    _MIN_WINDOW = 30
    region_len = hi - lo + 1
    if region_len < _MIN_WINDOW:
        flex_seq = None
    else:
        flex_seq = _primer3_single(
            arm, region_start=lo, region_len=region_len,
            excluded_region=[pinned_pos, 0], pick_left=flex_pick_left,
            primer_opt_tm=tm_target,
        )

    used_fallback = False
    effective_arm_length = nominal_length

    if pinned_result and flex_seq:
        # primer3's PRIMER_LEFT/RIGHT_0_SEQUENCE are already final, correctly
        # oriented oligo sequences (RIGHT is already the reverse-complement
        # binding to the top strand) — no further RC needed here.
        pinned_seq = pinned_result[0]
        if pinned_side == 'rev':
            fwd_bind, rev_bind = flex_seq, pinned_seq
            start = arm.find(fwd_bind)
            effective_arm_length = n - start if start >= 0 else nominal_length
        else:
            fwd_bind, rev_bind = pinned_seq, flex_seq
            binding_template = reverse_complement(rev_bind)   # back to genomic (top-strand) substring for lookup
            end = arm.find(binding_template)
            effective_arm_length = (end + len(rev_bind)) if end >= 0 else nominal_length
    else:
        used_fallback = True
        exact_arm = arm[-nominal_length:] if pinned_side == 'rev' else arm[:nominal_length]
        fwd_bind, rev_bind = _tm_growth_bind_regions(exact_arm, tm_target)
        effective_arm_length = nominal_length

    fwd_seq = tail + left_oh        + fwd_bind
    rev_seq = tail + _rc(right_oh)  + rev_bind

    meta = {"effective_arm_length": effective_arm_length, "used_fallback": used_fallback}
    return ("_F", fwd_seq), ("_R", rev_seq), meta


# ── fragment dispatcher ────────────────────────────────────────────────────────

def build_fragment(seq, left_oh, right_oh, tm_target=60.0,
                   nominal_length=None, tolerance=0.2, pinned_side='rev'):
    """Dispatch to annealed-oligo pair or PCR primers based on sequence length.

    nominal_length/tolerance/pinned_side are forwarded to
    sapi_amplification_primers (PCR-primer branch only; the oligo-pair branch
    is always exact, no search/fallback concept applies there).

    Returns a dict with keys:
      kind   : "oligo_pair" or "primers"
      top    : top oligo / fwd primer sequence (5'→3')
      bottom : bottom oligo / rev primer sequence (5'→3') — field name "bottom"
      fwd    : alias for top (primers only)  — field name "fwd"
      rev    : alias for bottom (primers only)
      meta   : {"effective_arm_length", "used_fallback"} — only present for "primers"
    All sequences are returned regardless of which is the primary label.
    """
    if len(seq) <= OLIGO_MAX:
        top, bot = annealed_oligo_pair(seq, left_oh, right_oh)
        return {"kind": "oligo_pair", "top": top, "bottom": bot,
                "fwd": top, "rev": bot}
    else:
        (_, fwd), (_, rev), meta = sapi_amplification_primers(
            seq, left_oh, right_oh, tm_target,
            nominal_length=nominal_length, tolerance=tolerance, pinned_side=pinned_side,
        )
        return {"kind": "primers", "top": fwd, "bottom": rev,
                "fwd": fwd, "rev": rev, "meta": meta}


# ── GenBank record builder ─────────────────────────────────────────────────────

def sapi_flanked_record(seq, left_oh, right_oh, name, description=""):
    """Build a linear SeqRecord with flanking SapI sites for ApE Golden Gate simulation.

    Sequence structure (top strand, 5'→3'):
      GCTCTTCA + {left_oh} + {seq} + {RC(right_oh)} + A + GAAGAGC

    After SapI digestion this produces a fragment with:
      left  5′ overhang: left_oh
      right 5′ overhang: RC(right_oh)  (mates with right_oh on adjacent fragment)

    Features annotate the SapI sites, overhangs, and insert body.
    Works for any fragment length; short fragments (sgRNA, short HA) are
    represented in full-duplex form so ApE can simulate the digest.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    rc_right_oh = _rc(right_oh)

    # full double-strandable sequence
    full = SAPI_SITE + SAPI_SPR + left_oh + seq.upper() + rc_right_oh + SAPI_SPR + SAPI_RC
    record = SeqRecord(
        Seq(full),
        id=name,
        name=name[:16],           # GenBank NAME field ≤16 chars
        description=description or name,
        annotations={"molecule_type": "DNA"},
    )

    # Positions (0-indexed, half-open)
    p0 = 0
    p_loh_start  = len(SAPI_SITE) + len(SAPI_SPR)   # 8
    p_insert_start = p_loh_start + len(left_oh)       # 11
    p_insert_end   = p_insert_start + len(seq)
    p_roh_end     = p_insert_end + len(rc_right_oh)   # = p_insert_end + 3
    p_end         = len(full)

    def _feat(start, end, label, strand=1, ftype="misc_feature"):
        return SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type=ftype,
            qualifiers={"label": [label]},
        )

    record.features = [
        _feat(p0, p_loh_start, "SapI", strand=1),
        _feat(p_loh_start, p_insert_start, "left overhang ({})".format(left_oh), strand=1),
        _feat(p_insert_start, p_insert_end, name, strand=1),
        _feat(p_insert_end, p_roh_end, "right overhang ({})".format(right_oh), strand=1),
        _feat(p_roh_end + len(SAPI_SPR), p_end, "SapI", strand=-1),
    ]
    return record


# ── combined GenBank record ────────────────────────────────────────────────────

def combine_saptrap_records(records_by_type, site_name):
    """Concatenate sgRNA, 5'HA, and 3'HA records into one annotated GenBank record.

    Fragments are joined by LINKER_SEQ ("nnn") spacers.  Each fragment gets a
    colored spanning feature (ApEinfo_fwdcolor / ApEinfo_revcolor) so ApE
    displays each section in a distinct color; sub-features from each individual
    record are re-indexed and included beneath it.

    records_by_type : dict {"sgrna": SeqRecord, "ha5": SeqRecord, "ha3": SeqRecord}
    site_name       : string used as the record id / GenBank LOCUS name
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    order = ["sgrna", "ha5", "ha3"]
    combined_seq = ""
    offsets = {}

    for i, key in enumerate(order):
        offsets[key] = len(combined_seq)
        combined_seq += str(records_by_type[key].seq)
        if i < len(order) - 1:
            combined_seq += LINKER_SEQ

    record = SeqRecord(
        Seq(combined_seq),
        id=site_name,
        name=site_name[:16],
        description="{} SapTrap fragments: sgRNA + 5'HA + 3'HA".format(site_name),
        annotations={"molecule_type": "DNA"},
    )

    features = []
    for i, key in enumerate(order):
        off  = offsets[key]
        rec  = records_by_type[key]
        color = FRAG_COLORS[key]
        flen  = len(rec.seq)

        # Spanning colored label feature for this fragment
        features.append(SeqFeature(
            FeatureLocation(off, off + flen, strand=1),
            type="misc_feature",
            qualifiers={
                "label":              [FRAG_LABELS[key]],
                "ApEinfo_fwdcolor":   [color],
                "ApEinfo_revcolor":   [color],
            },
        ))

        # Re-index sub-features from the individual record
        for f in rec.features:
            features.append(SeqFeature(
                FeatureLocation(
                    off + f.location.start,
                    off + f.location.end,
                    strand=f.location.strand,
                ),
                type=f.type,
                qualifiers=dict(f.qualifiers),
            ))

        # Linker feature between fragments
        if i < len(order) - 1:
            ls = off + flen
            features.append(SeqFeature(
                FeatureLocation(ls, ls + len(LINKER_SEQ), strand=1),
                type="misc_feature",
                qualifiers={
                    "label":            ["linker"],
                    "ApEinfo_fwdcolor": ["#aaaaaa"],
                    "ApEinfo_revcolor": ["#aaaaaa"],
                },
            ))

    record.features = features
    return record


# ── per-site orchestration ─────────────────────────────────────────────────────

def build_saptrap_site(row, left_arm, right_arm, tm_target=60.0, genbank=False,
                       arm_length=None, tolerance=0.2):
    """Design all SapTrap reagents for one selected guide / tag site.

    row       : a row from the .reagents.tsv (pandas Series or dict-like)
    left_arm  : mutated left arm — exactly arm_length for short (oligo-pair)
                arms, or widened via truncate_arms_with_tolerance for long
                (PCR-primer) arms; see reagent_sequences.truncate_arms_with_tolerance
    right_arm : mutated right arm, same convention as left_arm
    tm_target : Tm target (°C) for PCR primer binding region
    genbank   : if True, also build SeqRecords for ApE export
    arm_length: the nominal (requested) arm length — used as the search-window
                center for the PCR-primer path; defaults to len(left_arm)/
                len(right_arm) if not given (i.e. no widening happened)
    tolerance : fraction of arm_length the far/outer primer end may search
                within, when the PCR-primer (long-arm) path is used

    Returns (oligo_rows, records) where:
      oligo_rows : list of dicts {name, sequence, kind, notes}
      records    : list of SeqRecords (empty if genbank=False)
    """
    spacer   = str(row["spacer"]).upper()
    amino    = str(row["amino_acid"])
    ridx     = int(row["residue_index"])
    site     = "{}{}".format(amino, ridx)
    recut    = str(row.get("recut_block_method", "") or "")
    has_mut  = recut in ("syn_1", "mut_1")

    oligo_rows = []
    records    = []
    frag_recs  = {}   # keyed "sgrna"/"ha5"/"ha3"; combined into one record if genbank

    # ── sgRNA oligo pair ──────────────────────────────────────────────────────
    sgrna_top, sgrna_bot = sgrna_oligos(spacer)
    oligo_rows += [
        {"name": "{}_sgRNA_top".format(site), "sequence": sgrna_top,
         "kind": "oligo_top", "notes": "TTG overhang; anneal with _bot"},
        {"name": "{}_sgRNA_bot".format(site), "sequence": sgrna_bot,
         "kind": "oligo_bot", "notes": "AAC overhang; anneal with _top"},
    ]
    if genbank:
        frag_recs["sgrna"] = sapi_flanked_record(
            spacer, *SAPTRAP_OVERHANGS["sgrna"], "{}_sgRNA".format(site))

    # ── 5′ homology arm ───────────────────────────────────────────────────────
    left_oh5, right_oh5 = SAPTRAP_OVERHANGS["ha5"]

    # domesticate internal SapI sites for short (oligo) arms before designing oligos
    left_arm_use = left_arm
    if len(left_arm) <= OLIGO_MAX and find_internal_sapi(left_arm):
        left_arm_use, dom5_changes = domesticate_sapi(left_arm, is_left_arm=True)
    else:
        dom5_changes = []

    left_nominal = arm_length if arm_length is not None else len(left_arm)
    frag5 = build_fragment(left_arm_use, left_oh5, right_oh5, tm_target,
                           nominal_length=left_nominal, tolerance=tolerance,
                           pinned_side='rev')

    if frag5["kind"] == "oligo_pair":
        dom5_note = ("SapI site domesticated: " +
                     ", ".join("pos{} {}→{} ({})".format(p, o, n, m)
                               for p, o, n, m in dom5_changes)
                     if dom5_changes else "")
        oligo_rows += [
            {"name": "{}_5HA_top".format(site), "sequence": frag5["top"],
             "kind": "oligo_top",
             "notes": "TGG overhang; anneal with _bot" +
                      ("; " + dom5_note if dom5_note else "")},
            {"name": "{}_5HA_bot".format(site), "sequence": frag5["bottom"],
             "kind": "oligo_bot", "notes": "CGC overhang; anneal with _top"},
        ]
    else:
        mut_note = ("PAM-disrupting mutation in rev primer — verify primer covers it"
                    if has_mut else "")
        meta5 = frag5["meta"]
        meta5_note = "effective_arm_length={}{}".format(
            meta5["effective_arm_length"],
            "; fallback=True" if meta5["used_fallback"] else "",
        )
        oligo_rows += [
            {"name": "{}_5HA_F".format(site), "sequence": frag5["fwd"],
             "kind": "primer_fwd", "notes": "SapI-tailed fwd; TGG overhang; " + meta5_note},
            {"name": "{}_5HA_R".format(site), "sequence": frag5["rev"],
             "kind": "primer_rev", "notes": "SapI-tailed rev; GCG overhang" +
                                            ("; " + mut_note if mut_note else "") +
                                            "; " + meta5_note},
        ]

    if genbank:
        frag_recs["ha5"] = sapi_flanked_record(
            left_arm_use.upper(), left_oh5, right_oh5, "{}_5HA".format(site))

    # ── 3′ homology arm ───────────────────────────────────────────────────────
    left_oh3, right_oh3 = SAPTRAP_OVERHANGS["ha3"]

    right_arm_use = right_arm
    if len(right_arm) <= OLIGO_MAX and find_internal_sapi(right_arm):
        right_arm_use, dom3_changes = domesticate_sapi(right_arm, is_left_arm=False)
    else:
        dom3_changes = []

    right_nominal = arm_length if arm_length is not None else len(right_arm)
    frag3 = build_fragment(right_arm_use, left_oh3, right_oh3, tm_target,
                           nominal_length=right_nominal, tolerance=tolerance,
                           pinned_side='fwd')

    if frag3["kind"] == "oligo_pair":
        dom3_note = ("SapI site domesticated: " +
                     ", ".join("pos{} {}→{} ({})".format(p, o, n, m)
                               for p, o, n, m in dom3_changes)
                     if dom3_changes else "")
        oligo_rows += [
            {"name": "{}_3HA_top".format(site), "sequence": frag3["top"],
             "kind": "oligo_top",
             "notes": "ACG overhang; anneal with _bot" +
                      ("; " + dom3_note if dom3_note else "")},
            {"name": "{}_3HA_bot".format(site), "sequence": frag3["bottom"],
             "kind": "oligo_bot", "notes": "CGT overhang; anneal with _top"},
        ]
    else:
        mut_note = ("PAM-disrupting mutation in fwd primer — verify primer covers it"
                    if has_mut else "")
        meta3 = frag3["meta"]
        meta3_note = "effective_arm_length={}{}".format(
            meta3["effective_arm_length"],
            "; fallback=True" if meta3["used_fallback"] else "",
        )
        oligo_rows += [
            {"name": "{}_3HA_F".format(site), "sequence": frag3["fwd"],
             "kind": "primer_fwd", "notes": "SapI-tailed fwd; ACG overhang" +
                                            ("; " + mut_note if mut_note else "") +
                                            "; " + meta3_note},
            {"name": "{}_3HA_R".format(site), "sequence": frag3["rev"],
             "kind": "primer_rev", "notes": "SapI-tailed rev; GTA overhang; " + meta3_note},
        ]

    if genbank:
        frag_recs["ha3"] = sapi_flanked_record(
            right_arm_use.upper(), left_oh3, right_oh3, "{}_3HA".format(site))

    if genbank:
        records.append(combine_saptrap_records(frag_recs, site))

    return oligo_rows, records


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from argparse import ArgumentParser
    import pandas as pd
    from Bio import SeqIO

    parser = ArgumentParser(
        description="Design SapTrap donor reagents from a design_tag_reagents TSV."
    )
    parser.add_argument("--reagents_tsv", required=True, help="Path to .reagents.tsv")
    parser.add_argument("--arm_length",   type=int, default=500,
                        help="Homology arm length in bp (default: 500)")
    parser.add_argument("--tm",           type=float, default=60.0,
                        help="Primer binding Tm target °C for PCR arms (default: 60)")
    parser.add_argument("--tolerance",    type=float, default=0.2,
                        help="Fraction of arm_length the far/outer PCR primer end may "
                             "search within, for arms > OLIGO_MAX (default: 0.2)")
    parser.add_argument("--genbank",      action="store_true",
                        help="Write a GenBank file per site for ApE Golden Gate simulation")
    parser.add_argument("--output",       required=True,
                        help="Output path prefix; writes {prefix}_saptrap-oligos.tsv"
                             " and (with --genbank) {prefix}_{site}_saptrap.gb")
    args = parser.parse_args()

    from reagent_sequences import truncate_arms, truncate_arms_with_tolerance

    df = pd.read_csv(args.reagents_tsv, sep="\t")

    all_rows = []
    for _, row in df.iterrows():
        try:
            if args.arm_length > OLIGO_MAX:
                left, right = truncate_arms_with_tolerance(
                    str(row["left_arm"]), str(row["right_arm"]),
                    args.arm_length, args.tolerance,
                )
            else:
                left, right = truncate_arms(
                    str(row["left_arm"]), str(row["right_arm"]), args.arm_length
                )
        except ValueError as e:
            print("WARNING: {} — skipping row".format(e), file=sys.stderr)
            continue

        # check for internal SapI sites in the arms/spacer
        spacer = str(row["spacer"]).upper()
        for seq_name, seq_val in [("spacer", spacer), ("5'HA", left), ("3'HA", right)]:
            hits = find_internal_sapi(seq_val)
            if hits:
                print("WARNING: internal SapI site in {} for site {}{} — "
                      "this fragment must be domesticated.".format(
                          seq_name, row["amino_acid"], int(row["residue_index"])),
                      file=sys.stderr)

        oligo_rows, records = build_saptrap_site(
            row, left, right, tm_target=args.tm, genbank=args.genbank,
            arm_length=args.arm_length, tolerance=args.tolerance,
        )
        for orow in oligo_rows:
            if "fallback=True" in orow.get("notes", ""):
                print("WARNING: {} — primer3 design failed; fell back to Tm-growth "
                      "primers at the exact {} bp arm".format(
                          orow["name"], args.arm_length), file=sys.stderr)
        all_rows.extend(oligo_rows)

        if args.genbank and records:
            # one combined .gb per site (sgRNA + 5'HA + 3'HA in a single record)
            rec = records[0]
            gb_path = "{}_{}_saptrap.gb".format(args.output, rec.id)
            with open(gb_path, "w") as fh:
                SeqIO.write(rec, fh, "genbank")
            print("Wrote {}".format(gb_path), file=sys.stderr)

    out_df = pd.DataFrame(all_rows)
    tsv_path = "{}_saptrap-oligos.tsv".format(args.output)
    out_df.to_csv(tsv_path, sep="\t", index=False)
    print("Wrote {} rows to {}".format(len(out_df), tsv_path), file=sys.stderr)
