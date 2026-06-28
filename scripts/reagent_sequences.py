"""
reagent_sequences.py

Pure sequence-assembly logic for CRISPR editing reagents: truncate precomputed
homology arms, build ssODN / PCR primers, calculate Tm, and render ASCII diagrams.

All functions are stateless — no Shiny imports.  Designed for dual use:
in-process calls from reagents_server.py and standalone CLI.

Matt Rich, 2025
"""

import csv
import os
import sys
from html import escape as _esc

from Bio.SeqUtils import MeltingTemp as _mt

sys.path.insert(0, os.path.dirname(__file__))
from crispr_util import reverse_complement, CODONS


def _rc_case(seq):
    """Reverse complement preserving uppercase/lowercase case encoding."""
    _comp = {'A':'T','C':'G','G':'C','T':'A','N':'N',
             'a':'t','c':'g','g':'c','t':'a','n':'n'}
    return ''.join(_comp.get(c, c) for c in seq[::-1])


# ── Tag library ───────────────────────────────────────────────────────────────

def load_tags(tags_path):
    """Load tables/tags.tsv → OrderedDict {name: dna_sequence}."""
    tags = {}
    with open(tags_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            name = row.get('name', '').strip()
            seq  = row.get('dna_sequence', '').strip()
            if name and seq:
                tags[name] = seq
    return tags


# ── Arm handling ──────────────────────────────────────────────────────────────

def truncate_arms(left_arm, right_arm, arm_length):
    """Truncate precomputed arms to arm_length bp from the insertion junction.

    Takes from the 3' end of left_arm and 5' end of right_arm so that the bases
    closest to the insertion site (including PAM-disrupting mutations) are retained.
    Raises ValueError when arm_length exceeds the stored arm length.
    """
    max_len = min(len(left_arm), len(right_arm))
    if arm_length > len(left_arm) or arm_length > len(right_arm):
        raise ValueError(
            'Requested arm_length {} exceeds stored arm length {} / {} — '
            'rerun the reagents pipeline with a larger arm_length.'.format(
                arm_length, len(left_arm), len(right_arm))
        )
    return left_arm[-arm_length:], right_arm[:arm_length]


# ── Tm calculation ────────────────────────────────────────────────────────────

def calc_tm(seq):
    """Nearest-neighbour Tm (°C) for a DNA sequence; returns 0.0 if too short."""
    seq = seq.upper().replace('U', 'T')
    if len(seq) < 8:
        return 0.0
    try:
        return float(_mt.Tm_NN(seq))
    except Exception:
        return 0.0


# ── ssODN assembly ────────────────────────────────────────────────────────────

def build_ssodn(left_arm, right_arm, insert_seq, cut_pos, insert_pos, consider_strand):
    """Build a single-stranded oligo donor from truncated arms and an insert.

    Strand selection (when consider_strand=True):
      cut 5' of insert (cut_pos < insert_pos) → reverse-complement the ssODN
      cut at or 3' of insert                  → keep forward strand

    Returns (sequence, total_length, strand_used '+'/'-').
    """
    oligo = left_arm + insert_seq + right_arm
    if consider_strand and cut_pos < insert_pos:
        oligo = _rc_case(oligo)   # preserves exonic/intronic case encoding
        strand = '-'
    else:
        strand = '+'
    return oligo, len(oligo), strand


# ── PCR primer design ─────────────────────────────────────────────────────────

def design_pcr_primers(left_arm, right_arm, template, tm_target, phos,
                       cut_pos, insert_pos):
    """Design PCR primers that add homology arms to a tag template.

    Forward primer = left_arm + first K bases of template where Tm(template[:K]) >= tm_target.
    Reverse primer = RC(right_arm) + RC(last K bases of template).
    Tm is calculated on the template-binding region only (not the arm overhang).

    Phosphorylation (for lambda-exo ssDNA): same strand logic as ssODN.
      cut 5' of insert → phosphorylate forward primer (to produce antisense ssDNA)
      cut at/3' of insert → phosphorylate reverse primer (to produce sense ssDNA)
    /5Phos/ IDT modifier is prepended to the sequence when phos=True.

    Returns ((fwd_suffix, fwd_seq), (rev_suffix, rev_seq)).
    fwd/rev_suffix are '_F' / '_R' (already include phos tag in seq, not suffix).
    """
    template = template.upper()

    # Grow template-binding region until Tm >= target (or exhaust template)
    k_fwd = 12
    while k_fwd < len(template) and calc_tm(template[:k_fwd]) < tm_target:
        k_fwd += 1

    k_rev = 12
    while k_rev < len(template) and calc_tm(template[-k_rev:]) < tm_target:
        k_rev += 1

    fwd_bind = template[:k_fwd]
    rev_bind  = reverse_complement(template[-k_rev:])

    fwd_seq = left_arm + fwd_bind
    rev_seq = _rc_case(right_arm) + rev_bind   # preserve exonic/intronic case in arm overhang

    # Phosphorylate the primer that produces the desired repair strand
    if phos:
        if cut_pos < insert_pos:
            # want antisense ssDNA after lambda exo → phos fwd primer
            fwd_seq = '/5Phos/' + fwd_seq
        else:
            # want sense ssDNA → phos rev primer
            rev_seq = '/5Phos/' + rev_seq

    return ('_F', fwd_seq), ('_R', rev_seq)


# ── Translation helper ────────────────────────────────────────────────────────

def _build_aa_line(seq_raw, local_ins):
    """Return a string of same length as seq_raw with AAs centered in coding triplets.

    Uppercase = exonic (coding); lowercase = intronic (skipped).
    Codons are read outward from local_ins (the insert/junction point).
    """
    ann = [' '] * len(seq_raw)

    # Right side: codons read forward from local_ins
    i = local_ins
    while i < len(seq_raw):
        cod_pos = []
        j = i
        while j < len(seq_raw) and len(cod_pos) < 3:
            ch = seq_raw[j]
            if ch.isupper():
                cod_pos.append(j)
                j += 1
            elif ch.islower():
                j += 1   # skip intronic
            else:
                break
        if len(cod_pos) == 3:
            codon = ''.join(seq_raw[k].upper() for k in cod_pos)
            ann[cod_pos[1]] = CODONS.get(codon, 'X')
            i = j
        else:
            break

    # Left side: codons read backward from local_ins
    i = local_ins - 1
    while i >= 0:
        cod_pos = []
        j = i
        while j >= 0 and len(cod_pos) < 3:
            ch = seq_raw[j]
            if ch.isupper():
                cod_pos.insert(0, j)   # prepend to keep left-to-right order
                j -= 1
            elif ch.islower():
                j -= 1
            else:
                break
        if len(cod_pos) == 3:
            codon = ''.join(seq_raw[k].upper() for k in cod_pos)
            ann[cod_pos[1]] = CODONS.get(codon, 'X')
            i = j
        else:
            break

    return ''.join(ann)


# ── ASCII diagram ─────────────────────────────────────────────────────────────

def ascii_diagram(guide_row, left_arm, right_arm, left_wt=None, right_wt=None):
    """Return an HTML string for embedding in a <pre> block.

    Rows (top to bottom):
      1. WT genomic sequence (PAM restored, blue; extends into wt context)
      2. Translation (single-letter AA centered in each coding triplet)
      3. Repair template (arm with PAM mutations; mutated PAM red, context dimmed)
      4. Annotation (^ insert site green+bold, | cut site)
      5. Guide (>>>> PAM for + strand; PAM <<<< for - strand)

    left_wt / right_wt: stored arms trimmed to arm_length+10 for WT display context.
    If omitted, same as left_arm / right_arm (no extra context shown).
    """
    insert_pos    = int(guide_row['insert_pos'])
    cut_pos       = int(guide_row['cut_pos'])
    pam_fwd_start = int(guide_row['pam_fwd_start'])
    pam_seq       = str(guide_row['pam_seq'])
    pam_len       = len(pam_seq)
    spacer_len    = len(str(guide_row['spacer']))
    strand        = str(guide_row['guide_strand'])
    recut_method  = str(guide_row.get('recut_block_method', 'none') or 'none')
    mutation_desc = str(guide_row.get('mutation_desc', '') or '')
    # a PAM mutation exists when method is syn_1 or mut_1 (not insertion or none)
    mutated = recut_method in ('syn_1', 'mut_1')
    dist    = abs(cut_pos - insert_pos)

    # Display window is based on the WT arms (wider than truncated repair arms)
    l_wt = left_wt if left_wt is not None else left_arm
    r_wt = right_wt if right_wt is not None else right_arm
    MAX  = 100   # hard cap per side
    l_show = l_wt[-MAX:]
    r_show = r_wt[:MAX]

    win_len   = len(l_show) + len(r_show)
    local_ins = len(l_show)

    gds             = insert_pos - local_ins
    local_cut       = cut_pos       - gds
    local_pam_start = pam_fwd_start - gds

    def _in_pam(i):
        return local_pam_start <= i < local_pam_start + pam_len

    # WT fwd-strand bases at PAM positions (left-to-right in genomic coords).
    # pam_seq is 5'→3' on the guide strand.
    # For + strand: guide strand = fwd strand, so pam_seq IS the fwd bases.
    # For - strand: guide strand = rev strand, so fwd bases = RC(pam_seq).
    wt_pam_fwd = pam_seq.upper() if strand == '+' else reverse_complement(pam_seq).upper()

    # Arm offset within the display window
    left_ctx  = local_ins - len(left_arm)    # chars of WT context to the left of the arm
    right_ctx = len(r_show) - len(right_arm) # chars of WT context to the right of the arm

    pfx = '   '   # 3-char left prefix (space or '...')

    # ── Row 1: WT sequence ────────────────────────────────────────────────────
    wt_raw = l_show + r_show   # mutated arm; we restore PAM positions below
    wt_parts = [pfx]
    for i, ch in enumerate(wt_raw):
        if _in_pam(i):
            k = i - local_pam_start
            base = wt_pam_fwd[k] if k < len(wt_pam_fwd) else ch.upper()
            wt_parts.append('<span style="color:#4a90d9">{}</span>'.format(_esc(base)))
        else:
            wt_parts.append(_esc(ch))
    wt_line = ''.join(wt_parts)

    # ── Row 2: translation (AAs centered in coding triplets) ─────────────────
    aa_raw  = _build_aa_line(wt_raw, local_ins)
    aa_line = pfx + _esc(aa_raw)

    # ── Row 3: repair template ────────────────────────────────────────────────
    # Context bases (outside the arm) shown in dim gray; arm bases shown normally
    # with mutated PAM positions in red, unmutated PAM in blue.
    def _arm_ch(i):
        """Character from the repair arm at display position i."""
        if local_ins > i >= left_ctx:
            return left_arm[i - left_ctx]
        if local_ins <= i < local_ins + len(right_arm):
            return right_arm[i - local_ins]
        return wt_raw[i]   # context: same as WT

    rep_parts = [pfx]
    for i, ch in enumerate(wt_raw):
        in_ctx = i < left_ctx or i >= local_ins + len(right_arm)
        arm_ch = _arm_ch(i)
        if in_ctx:
            rep_parts.append('<span style="color:#bbb">{}</span>'.format(_esc(arm_ch)))
        elif _in_pam(i) and mutated:
            rep_parts.append('<span style="color:#c0392b">{}</span>'.format(_esc(arm_ch)))
        elif _in_pam(i):
            rep_parts.append('<span style="color:#4a90d9">{}</span>'.format(_esc(arm_ch)))
        else:
            rep_parts.append(_esc(arm_ch))
    rep_line = ''.join(rep_parts)

    # ── Row 4: annotation (^ insert, | cut) ──────────────────────────────────
    # Both markers are shifted left by 0.5ch so they appear between two bases
    # rather than under the right-hand one (monospace cells are whole-width only).
    _HALF = 'position:relative;left:-0.5ch'
    ann_parts = [pfx]
    for i in range(win_len):
        if i == local_ins and 0 <= local_ins < win_len:
            ann_parts.append(
                '<span style="color:#27ae60;font-weight:bold;{}">' \
                '^</span>'.format(_HALF)
            )
        elif i == local_cut and 0 <= local_cut < win_len:
            ann_parts.append('<span style="{}">|</span>'.format(_HALF))
        else:
            ann_parts.append(' ')
    ann_line = ''.join(ann_parts)

    # ── Row 5: guide line ─────────────────────────────────────────────────────
    guide_slots = [''] * win_len
    if strand == '+':
        # + strand: >>>>PAM  (PAM label = pam_seq left-to-right in fwd coords)
        pam_label = pam_seq.upper()
        for k in range(spacer_len):
            p = local_pam_start - spacer_len + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#888">&gt;</span>'
        for k, ch in enumerate(pam_label):
            p = local_pam_start + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#4a90d9">{}</span>'.format(_esc(ch))
    else:
        # - strand: PAM<<<<  (PAM label reversed so letters match fwd-strand positions)
        # pam_seq is 5'→3' on the minus strand (right-to-left in fwd coords).
        # Reversing gives the same letters left-to-right in fwd display order.
        pam_label = pam_seq[::-1].lower()
        for k, ch in enumerate(pam_label):
            p = local_pam_start + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#4a90d9">{}</span>'.format(_esc(ch))
        for k in range(spacer_len):
            p = local_pam_start + pam_len + k
            if 0 <= p < win_len:
                guide_slots[p] = '<span style="color:#888">&lt;</span>'

    guide_line = pfx + ''.join(s if s else ' ' for s in guide_slots)

    # ── Footer ────────────────────────────────────────────────────────────────
    if mutation_desc and recut_method not in ('none', ''):
        recut_str = '{}: {}'.format(_esc(recut_method), _esc(mutation_desc))
    elif recut_method == 'none':
        recut_str = 'no PAM disruption'
    else:
        recut_str = _esc(recut_method)

    footer = 'cut&#8596;insert: {} bp  PAM: {}  strand: {}  recut: {}'.format(
        dist, _esc(pam_seq.upper()), strand, recut_str
    )

    return '\n'.join([wt_line, aa_line, rep_line, ann_line, guide_line, footer])


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == '__main__':
    from argparse import ArgumentParser
    import pandas as pd

    parser = ArgumentParser(
        description='Assemble editing reagents from a design_tag_reagents TSV.'
    )
    parser.add_argument('--reagents_tsv',  required=True, help='Path to .reagents.tsv')
    parser.add_argument('--repair_type',   default='ssodn',
                        choices=['ha', 'ssodn', 'amplicon'],
                        help='Repair strategy (default: ssodn)')
    parser.add_argument('--arm_length',    type=int, default=60,
                        help='Homology arm length in bp (default: 60)')
    parser.add_argument('--insert_seq',    default='', help='Tag/insert DNA sequence')
    parser.add_argument('--consider_strand', action='store_true', default=True,
                        help='Apply repair-strand logic for ssODN (default: on)')
    parser.add_argument('--pcr_template',  default='', help='Template sequence for amplicon')
    parser.add_argument('--pcr_tm',        type=float, default=60.0,
                        help='Primer Tm target (°C) for amplicon (default: 60)')
    parser.add_argument('--pcr_phos',      action='store_true',
                        help='Add 5-phosphate to one PCR primer for lambda-exo')
    parser.add_argument('--output',        required=True, help='Output TSV path')
    args = parser.parse_args()

    df   = pd.read_csv(args.reagents_tsv, sep='\t')
    rows = []
    for _, row in df.iterrows():
        try:
            left, right = truncate_arms(row['left_arm'], row['right_arm'], args.arm_length)
        except ValueError as e:
            print('WARNING: {}'.format(e), file=sys.stderr)
            continue

        site  = '{}{}'.format(row['amino_acid'], row['residue_index'])
        crrna = str(row['spacer'])
        cut   = int(row['cut_pos'])
        ins   = int(row['insert_pos'])

        if args.repair_type == 'ssodn':
            seq, total, strand = build_ssodn(
                left, right, args.insert_seq.upper(), cut, ins, args.consider_strand
            )
            rows.append({'name': '{}_ssODN'.format(site), 'sequence': seq,
                         'total_len': total, 'strand': strand})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

        elif args.repair_type == 'amplicon' and args.pcr_template:
            (_, fwd), (_, rev) = design_pcr_primers(
                left, right, args.pcr_template.upper(), args.pcr_tm,
                args.pcr_phos, cut, ins
            )
            rows.append({'name': '{}_F'.format(site), 'sequence': fwd})
            rows.append({'name': '{}_R'.format(site), 'sequence': rev})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

        else:  # ha only
            rows.append({'name': '{}_5HA'.format(site), 'sequence': left})
            rows.append({'name': '{}_3HA'.format(site), 'sequence': right})
            rows.append({'name': '{}_crRNA'.format(site), 'sequence': crrna})

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep='\t', index=False)
    print('Wrote {} rows to {}'.format(len(out_df), args.output), file=sys.stderr)
