"""
crispr_util.py

Shared utilities for CRISPR guide design and PAM-disruption donor construction.
Extracted so parse_genewise.py, find_guides.py, design_tag_reagents.py, and
design_guides_across_region.py all share one copy.

Matt Rich, 2025
"""

import re

# ── IUPAC / sequence helpers ──────────────────────────────────────────────────

IUPAC_DICT = {
    'A': 'A',    'C': 'C',    'G': 'G',    'T': 'T',    'U': 'U',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]','D': '[AGT]',
    'H': '[ACT]','V': '[ACG]','N': '[ACGT]',
}

IUPAC_COMPLEMENT = {
    'A':'T','T':'A','C':'G','G':'C','N':'N','U':'A',
    'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
    'B':'V','V':'B','D':'H','H':'D',
}


def iupac_to_regex(pam):
    """Convert an IUPAC PAM string to a regex pattern string."""
    return ''.join(IUPAC_DICT.get(c.upper(), '') for c in pam)


def reverse_complement(seq):
    """Reverse complement a DNA sequence (ACGT + N only)."""
    dna_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(dna_dict[x.upper()] for x in seq[::-1])


def reverse_complement_iupac(pam):
    """Reverse complement of an IUPAC PAM string (preserves IUPAC codes)."""
    return ''.join(IUPAC_COMPLEMENT.get(c.upper(), c) for c in pam[::-1])


# ── Codon tables ──────────────────────────────────────────────────────────────

CODONS = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_',
    'TGC':'C','TGT':'C','TGA':'_','TGG':'W',
}

# Reverse map: amino acid → list of synonymous codons (excluding stops)
SYN_CODONS = {}
for _codon, _aa in CODONS.items():
    if _aa == '_':
        continue
    SYN_CODONS.setdefault(_aa, []).append(_codon)


def gc_count(seq):
    return sum(1 for b in seq.upper() if b in 'GC')


# ── Guide finding ─────────────────────────────────────────────────────────────

def find_guides(seq, pam='NGG', guide_length=20, cut_offset=3):
    """
    Find all guide RNA spacers + cut sites in *seq* on both strands.

    PAM is 3′ of the guide spacer on the non-template (protospacer) strand.
    SpCas9 cuts cut_offset bp upstream of the PAM (default 3).

    All coordinates are 0-based on the *forward* (input) sequence.

    Returns a list of dicts:
      strand         '+' or '-'
      spacer         guide sequence (5′→3′ on the guide strand)
      pam_seq        PAM sequence (on the guide strand)
      pam_fwd_start  start of the PAM region in fwd sequence coords
      guide_fwd_start, guide_fwd_end  fwd coords of the spacer
      cut_pos        fwd coord of the DSB (first base of right-side fragment)
      three_prime_C  True if 3′ base of spacer is C (quality flag)
    """
    pam_regex = iupac_to_regex(pam)
    pam_len = len(pam)
    L = len(seq)
    results = []

    # --- Forward strand (+): spacer immediately 5′ of PAM ---
    for m in re.finditer('(?=({}))'.format(pam_regex), seq.upper()):
        pam_start = m.start()
        guide_start = pam_start - guide_length
        if guide_start < 0:
            continue
        spacer = seq[guide_start:pam_start].upper()
        results.append({
            'strand': '+',
            'spacer': spacer,
            'pam_seq': seq[pam_start:pam_start + pam_len].upper(),
            'pam_fwd_start': pam_start,
            'guide_fwd_start': guide_start,
            'guide_fwd_end': pam_start,
            # Cut is cut_offset bp upstream (5′) of PAM on + strand
            'cut_pos': pam_start - cut_offset,
            'three_prime_C': spacer[-1] == 'C' if spacer else False,
        })

    # --- Reverse strand (−): search on RC sequence, convert coords ---
    seq_rc = reverse_complement(seq)
    for m in re.finditer('(?=({}))'.format(pam_regex), seq_rc.upper()):
        rc_pam_start = m.start()
        rc_guide_start = rc_pam_start - guide_length
        if rc_guide_start < 0:
            continue
        spacer_rc = seq_rc[rc_guide_start:rc_pam_start].upper()

        # PAM in fwd coords: RC positions [rc_pam_start, rc_pam_start+pam_len)
        # correspond to fwd positions [L-rc_pam_start-pam_len, L-rc_pam_start)
        fwd_pam_start = L - rc_pam_start - pam_len
        # Guide in fwd coords: immediately 3′ of fwd_pam on the plus strand
        fwd_guide_start = fwd_pam_start + pam_len
        fwd_guide_end = fwd_guide_start + guide_length

        # Cut: cut_offset bp 3′ of fwd_pam end (i.e., into the fwd guide region)
        # (mirrors the '+' formula reflected across the PAM)
        cut_pos = fwd_pam_start + pam_len + cut_offset

        results.append({
            'strand': '-',
            'spacer': spacer_rc,
            'pam_seq': seq_rc[rc_pam_start:rc_pam_start + pam_len].upper(),
            'pam_fwd_start': fwd_pam_start,
            'guide_fwd_start': fwd_guide_start,
            'guide_fwd_end': fwd_guide_end,
            'cut_pos': cut_pos,
            'three_prime_C': spacer_rc[-1] == 'C' if spacer_rc else False,
        })

    return results


def _pam_disrupted(seq, pam, pam_fwd_start, pam_len, strand):
    """Return True if the PAM at pam_fwd_start is no longer intact in seq."""
    fwd_seq = seq[pam_fwd_start:pam_fwd_start + pam_len].upper()
    check_seq = reverse_complement(fwd_seq) if strand == '-' else fwd_seq
    return not bool(re.fullmatch(iupac_to_regex(pam), check_seq))


# ── Frame lookup ──────────────────────────────────────────────────────────────

def build_frame_lookup(cds_df, dna):
    """
    Build a dict mapping every coding genomic position to codon information.

    Keys: 0-based genomic position (int)
    Values: dict with keys:
      codon_start  – 0-based start of the codon in the genomic sequence
      pos_in_codon – position within the codon (0, 1, or 2)
      codon        – 3-character codon string (may span intron if split)
      aa           – single-letter amino acid (None if split codon)
      split        – True if codon spans an intron boundary

    Positions in introns are absent from the dict.
    """
    lookup = {}
    for _, exon in cds_df.iterrows():
        start = int(exon['start'])   # 0-indexed
        stop  = int(exon['stop'])    # 0-indexed; last base of last full codon
        frame = int(exon['frame'])
        frame_prime = 3 if frame == 0 else frame  # see parse_genewise.py for explanation

        for x in range(start + frame_prime, stop + 3, 3):
            codon_start = x - 3
            # split if codon crosses either exon boundary
            is_split = (codon_start < start) or (x - 1 > stop)

            if not is_split:
                codon = dna[codon_start:x].upper()
                aa    = CODONS.get(codon, 'X')
                for k, pos in enumerate(range(codon_start, x)):
                    lookup[pos] = {
                        'codon_start': codon_start,
                        'pos_in_codon': k,
                        'codon': codon,
                        'aa': aa,
                        'split': False,
                    }
            else:
                # Mark exon-side bases of the split codon; codon/aa unreliable
                for pos in range(start, x):
                    k = pos - codon_start  # position within the (split) codon
                    lookup[pos] = {
                        'codon_start': codon_start,
                        'pos_in_codon': k,
                        'codon': None,
                        'aa': None,
                        'split': True,
                    }
    return lookup


# ── PAM disruption ────────────────────────────────────────────────────────────

def disrupt_pam(seq, pam, pam_fwd_start, strand, frame_lookup,
                seed_len=15, min_seed_mm=2):
    """
    Return a mutated copy of *seq* in which the guide can no longer be re-cut
    after HDR, along with a description and the method used.  Every mutation is
    protein-preserving (synonymous or intronic); if re-cutting cannot be blocked
    without changing the protein the guide is rejected (returns None) so the
    caller can discard it.

    Ladder (applied in order, returns on first success):
      1. 'syn_1'    — one synonymous codon change that disrupts the PAM.
      2. 'mut_1'    — one silent (intronic) single-base change within the PAM.
      3. 'syn_seed' — synonymous codon change(s) placing >= min_seed_mm mismatches
                      within the seed_len nt of the protospacer proximal to the PAM
                      (blocks Cas9 re-binding without touching the PAM).

    Non-synonymous PAM changes are never made: if none of the above succeed
    (only possible when the seed has too little synonymous freedom, e.g. a run of
    Trp/Met codons), the function returns None and the guide should be discarded.

    Parameters
    ----------
    seq           : str  full genomic sequence (local coords, forward strand)
    pam           : str  IUPAC PAM string (e.g. 'NGG')
    pam_fwd_start : int  0-based start of the PAM in fwd seq coords
    strand        : str  '+' or '-' (guide strand)
    frame_lookup  : dict output of build_frame_lookup()
    seed_len      : int  protospacer seed length (bp proximal to PAM)
    min_seed_mm   : int  synonymous mismatches required in the seed for 'syn_seed'

    Returns
    -------
    (mutated_seq, description, method) on success, or None if impossible.
    method is one of 'syn_1', 'mut_1', 'syn_seed'.
    """
    pam_len = len(pam)
    pam_positions = list(range(pam_fwd_start, pam_fwd_start + pam_len))

    # Helper: apply a set of base substitutions {pos: new_base} to seq
    def _apply(changes):
        s = list(seq)
        for pos, base in changes.items():
            s[pos] = base
        return ''.join(s)

    # Collect which codon(s) overlap the PAM and are safely coding (not split)
    codons_in_pam = {}   # codon_start → codon info dict
    for pos in pam_positions:
        info = frame_lookup.get(pos)
        if info and not info['split']:
            cs = info['codon_start']
            if cs not in codons_in_pam:
                codons_in_pam[cs] = info

    # ── Step 1: single synonymous codon change disrupting the PAM ─────────────
    for cs, info in codons_in_pam.items():
        orig_codon = info['codon']
        aa = info['aa']
        for alt in SYN_CODONS.get(aa, []):
            if alt == orig_codon:
                continue
            # substitute only the codon bases that differ
            changes = {cs + k: alt[k] for k in range(3)
                       if orig_codon[k] != alt[k]}
            if not changes:
                continue
            new_seq = _apply(changes)
            if _pam_disrupted(new_seq, pam, pam_fwd_start, pam_len, strand):
                desc = '{}>{} at genomic pos {}'.format(orig_codon, alt, cs)
                return (new_seq, desc, 'syn_1')

    # ── Step 2: single silent (intronic) base change within the PAM ───────────
    # Coding PAM positions are skipped here — any coding PAM base change would be
    # non-synonymous; those are handled protein-preservingly via the seed instead.
    for pos in pam_positions:
        if pos in frame_lookup:
            continue   # coding → skip (non-synonymous)
        orig = seq[pos].upper()
        for base in 'ACGT':
            if base == orig:
                continue
            new_seq = _apply({pos: base})
            if _pam_disrupted(new_seq, pam, pam_fwd_start, pam_len, strand):
                desc = 'pos {} {}>{} (intronic)'.format(pos, orig, base)
                return (new_seq, desc, 'mut_1')

    # ── Step 3: synonymous mismatches in the protospacer seed ─────────────────
    return _disrupt_seed(seq, pam_fwd_start, pam_len, strand,
                         frame_lookup, seed_len, min_seed_mm)


def _disrupt_seed(seq, pam_fwd_start, pam_len, strand, frame_lookup,
                  seed_len, min_mm):
    """Introduce >= min_mm synonymous mismatches within the seed_len bp of the
    protospacer proximal to the PAM.  Returns (mutated_seq, desc, 'syn_seed') or
    None if that many synonymous mismatches can't be placed (e.g. Trp/Met run).

    Seed positions are enumerated nearest-to-PAM first so mismatches are placed
    as close to the PAM as possible (most disruptive to Cas9 re-binding).
    """
    pam_end = pam_fwd_start + pam_len
    # protospacer seed positions (fwd coords), nearest the PAM first
    if strand == '+':
        seed_positions = list(range(pam_fwd_start - 1, pam_fwd_start - 1 - seed_len, -1))
    else:
        seed_positions = list(range(pam_end, pam_end + seed_len))
    seed_positions = [p for p in seed_positions if 0 <= p < len(seq)]
    seed_set = set(seed_positions)

    # coding, non-split codons overlapping the seed, ordered by proximity to PAM
    codons = []
    seen = set()
    for p in seed_positions:
        info = frame_lookup.get(p)
        if info and not info['split'] and info['codon_start'] not in seen:
            seen.add(info['codon_start'])
            codons.append(info)

    changes = {}          # pos → new base
    mm_positions = set()   # seed positions actually changed
    descs = []
    for info in codons:
        if len(mm_positions) >= min_mm:
            break
        cs = info['codon_start']
        orig_codon, aa = info['codon'], info['aa']
        # pick the synonymous alt adding the most NEW seed mismatches (fewest total edits on ties)
        best = None   # (n_new_seed_mm, -n_total_changes, alt, codon_changes)
        for alt in SYN_CODONS.get(aa, []):
            if alt == orig_codon:
                continue
            cchanges = {cs + k: alt[k] for k in range(3) if orig_codon[k] != alt[k]}
            new_mm = [p for p in cchanges if p in seed_set and p not in mm_positions]
            if not new_mm:
                continue
            cand = (len(new_mm), -len(cchanges), alt, cchanges)
            if best is None or cand[:2] > best[:2]:
                best = cand
        if best is None:
            continue
        _, _, alt, cchanges = best
        changes.update(cchanges)
        mm_positions.update(p for p in cchanges if p in seed_set)
        descs.append('{}>{} at genomic pos {}'.format(orig_codon, alt, cs))

    if len(mm_positions) < min_mm:
        return None
    new_seq = list(seq)
    for pos, base in changes.items():
        new_seq[pos] = base
    desc = 'seed syn ({} mismatches in {} nt proximal to PAM): {}'.format(
        len(mm_positions), seed_len, '; '.join(descs))
    return (''.join(new_seq), desc, 'syn_seed')
