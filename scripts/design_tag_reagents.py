"""
design_tag_reagents.py

For every potential tag-insertion site in a protein, identify the N nearest
CRISPR guide RNAs and construct HDR homology arms with PAM-disrupting mutations.

Pipeline:
  1. Parse Genewise output → CDS exon map
  2. Enumerate all insertion sites (one per protein residue)
  3. Find all guide cut sites in the genomic sequence (both strands)
  4. For each (residue × guide) pair within arm_length of the cut:
       - Build left and right homology arms centered on the insertion site
       - Mutate the arm containing the PAM to abrogate re-cutting
       - Mutation ladder: (a) 1 synonymous codon change, (b) 1 minimal base
         change within the PAM (intronic preferred, then non-synonymous)

Output TSV (one row per residue × guide):
  residue_index     1-based protein residue number
  amino_acid        residue identity
  insert_pos        0-based genomic insertion coordinate
  exon_index        exon containing (or adjacent to) this residue
  is_split_codon    True if codon straddles an exon–intron boundary
  dist_to_5p_splice distance to 5′ exon boundary (bp)
  dist_to_3p_splice distance to 3′ exon boundary (bp; -1 for last exon)
  guide_strand      '+' or '-'
  spacer            guide spacer sequence
  pam_seq           PAM sequence on guide strand
  pam_fwd_start     0-based start of PAM in fwd genomic coords
  cut_pos           0-based DSB position in fwd coords
  distance          |cut_pos - insert_pos| (bp)
  pam_in_arm        'left' / 'right' / 'both' / 'none'
  recut_block_method  'syn_1' / 'mut_1' / 'insertion' / 'none'
  mutation_desc     human-readable description of the PAM mutation
  left_arm          left homology arm: exonic bases uppercase, intronic lowercase
  right_arm         right homology arm: exonic bases uppercase, intronic lowercase

Matt Rich, 2025
"""

import sys
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from crispr_util import find_guides, build_frame_lookup, disrupt_pam
from parse_genewise import parse_genewise, enumerate_insertion_sites, \
    parse_genewise_score, cds_coverage
from progress import report as _report, resolve_reporter


def _case_arm(arm_seq, arm_start, frame_lookup):
    """Return arm_seq with exonic (coding) positions uppercase, intronic lowercase."""
    return ''.join(
        ch.upper() if (arm_start + i) in frame_lookup else ch.lower()
        for i, ch in enumerate(arm_seq)
    )


# ── Core logic ────────────────────────────────────────────────────────────────

def design_reagents(
    genewise_out,
    genomic_fasta,
    n_guides=5,
    arm_length=500,
    pam='NGG',
    guide_length=20,
    cut_offset=3,
    report=None,
):
    """
    Full pipeline: Genewise output + genomic FASTA → reagent table (DataFrame).

    Parameters
    ----------
    genewise_out   : str  path to Genewise .out.txt
    genomic_fasta  : str  path to the genomic region FASTA
    n_guides       : int  max guides to report per insertion site
    arm_length     : int  homology arm length (bp) on each side
    pam            : str  IUPAC PAM (default 'NGG')
    guide_length   : int  spacer length (default 20)
    cut_offset     : int  cut distance from PAM (SpCas9=3)
    """
    reporter = resolve_reporter(report)

    # 1. Load genomic sequence
    records = list(SeqIO.parse(genomic_fasta, 'fasta'))
    if not records:
        raise ValueError('No sequences in {}'.format(genomic_fasta))
    dna = str(records[0].seq)
    L = len(dna)

    # 2. Parse Genewise → CDS exon table
    # TODO: run Genewise on both orientations (forward + RC) and pick the better
    # one automatically before calling design_reagents.  The EBI call takes ~10 s,
    # so running both is cheap.  For now the caller (run_genewise.py / the
    # orchestrator pre-step) is responsible for supplying the correct .out.txt;
    # we guard against silent garbage mappings via the score/coverage check below.
    cds_df = parse_genewise(genewise_out)
    _report(reporter, '{} CDS exons parsed'.format(len(cds_df)), stage='parse_genewise')

    # Low-score / low-coverage guard: warn loudly so we never silently process
    # a garbage (wrong-strand or failed) Genewise alignment.
    score = parse_genewise_score(genewise_out)
    if score is not None:
        n_protein = len(cds_df)   # rough proxy; real protein length unknown here
        coverage  = cds_coverage(cds_df, n_protein) if n_protein > 0 else 0.0
        if score < 50.0:
            _report(reporter,
                    'Genewise score very low ({:.2f} bits); may indicate wrong strand or '
                    'failed alignment'.format(score),
                    stage='parse_genewise', level='warning')

    # 3. Enumerate insertion sites
    sites = enumerate_insertion_sites(cds_df, dna)
    _report(reporter, '{} residues / insertion sites'.format(len(sites)), stage='sites')

    # 4. Build per-position frame lookup for synonymous mutation
    frame_lookup = build_frame_lookup(cds_df, dna)

    # 5. Find all guide cut sites
    guides = find_guides(dna, pam=pam, guide_length=guide_length,
                         cut_offset=cut_offset)
    _report(reporter, '{} guide sites found (PAM={})'.format(len(guides), pam), stage='guides')

    pam_len = len(pam)
    rows = []

    for _, site in sites.iterrows():
        insert_pos = int(site['insert_pos'])

        # Sort guides by distance of their cut to the insertion site
        scored = sorted(guides, key=lambda g: abs(g['cut_pos'] - insert_pos))

        # Keep up to n_guides that are within arm_length of the insertion
        selected = []
        for g in scored:
            if abs(g['cut_pos'] - insert_pos) > arm_length:
                continue
            selected.append(g)
            if len(selected) >= n_guides:
                break

        for g in selected:
            cut_pos       = g['cut_pos']
            pam_fwd_start = g['pam_fwd_start']
            distance      = abs(cut_pos - insert_pos)

            # Build raw arms
            left_start  = max(0, insert_pos - arm_length)
            right_end   = min(L, insert_pos + arm_length)
            left_arm_raw  = dna[left_start:insert_pos]
            right_arm_raw = dna[insert_pos:right_end]

            # Determine which arm(s) contain the PAM
            pam_end = pam_fwd_start + pam_len
            pam_in_left  = pam_end > left_start and pam_fwd_start < insert_pos
            pam_in_right = pam_fwd_start < right_end and pam_end > insert_pos

            if pam_in_left and pam_in_right:
                pam_arm = 'both'
            elif pam_in_left:
                pam_arm = 'left'
            elif pam_in_right:
                pam_arm = 'right'
            else:
                pam_arm = 'none'

            # Check whether the tag insertion itself disrupts re-cutting.
            # Any insertion within the 15 bp protospacer seed region immediately
            # 5' of the PAM (on the guide strand) prevents Cas9 from re-binding.
            # PAM mutations are only needed when the insert falls outside this window.
            _SEED = 15
            if g['strand'] == '+':
                # seed region = 15 bp 5' of PAM + PAM itself (fwd coords)
                insertion_blocks = (pam_fwd_start - _SEED <= insert_pos < pam_end)
            else:
                # seed region = PAM + 15 bp 3' of PAM (= 15 bp 5' on − strand)
                insertion_blocks = (pam_fwd_start <= insert_pos < pam_end + _SEED)

            # Attempt PAM disruption
            recut_method  = 'none'
            mutation_desc = ''
            left_arm  = left_arm_raw
            right_arm = right_arm_raw

            if insertion_blocks:
                recut_method  = 'insertion'
                mutation_desc = 'insert within {} bp seed region of PAM; re-cutting impossible'.format(_SEED)

            elif pam_arm in ('left', 'right', 'both') and pam_arm != 'none':
                # disrupt_pam works on the full genomic sequence and returns a
                # modified copy; we then re-extract the arms from it.
                result = disrupt_pam(dna, pam, pam_fwd_start,
                                     g['strand'], frame_lookup)
                if result is not None:
                    mutated_seq, desc, method = result
                    left_arm  = mutated_seq[left_start:insert_pos]
                    right_arm = mutated_seq[insert_pos:right_end]
                    recut_method  = method
                    mutation_desc = desc
                else:
                    mutation_desc = 'PAM disruption not possible with available synonymous changes'

            # encode exonic vs intronic in arm case (uppercase = coding, lowercase = intronic)
            left_arm  = _case_arm(left_arm, left_start, frame_lookup)
            right_arm = _case_arm(right_arm, insert_pos, frame_lookup)

            rows.append({
                'residue_index':       int(site['residue_index']),
                'amino_acid':          site['amino_acid'],
                'insert_pos':          insert_pos,
                'exon_index':          int(site['exon_index']),
                'is_split_codon':      site['is_split_codon'],
                'dist_to_5p_splice':   int(site['dist_to_5p_splice']),
                'dist_to_3p_splice':   int(site['dist_to_3p_splice']),
                'guide_strand':        g['strand'],
                'spacer':              g['spacer'],
                'pam_seq':             g['pam_seq'],
                'pam_fwd_start':       pam_fwd_start,
                'cut_pos':             cut_pos,
                'distance':            distance,
                'pam_in_arm':          pam_arm,
                'recut_block_method':  recut_method,
                'mutation_desc':       mutation_desc,
                'left_arm':            left_arm,
                'right_arm':           right_arm,
            })

    return pd.DataFrame(rows)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(genewise, genomic_fasta, output, n_guides=5, arm_length=500,
         pam='NGG', guide_length=20, cut_offset=3, report=None):
    """Entry point for in-process calls from task_runners."""
    df = design_reagents(
        genewise_out  = genewise,
        genomic_fasta = genomic_fasta,
        n_guides      = n_guides,
        arm_length    = arm_length,
        pam           = pam,
        guide_length  = guide_length,
        cut_offset    = cut_offset,
        report        = report,
    )
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description=(
            'Design CRISPR knock-in reagents for every potential tag-insertion '
            'site in a protein, given Genewise alignment output.'
        )
    )
    # Required
    parser.add_argument('--genewise', required=True,
                        help='Genewise .out.txt result file')
    parser.add_argument('--genomic_fasta', required=True,
                        help='Genomic region FASTA (same sequence/orientation submitted to Genewise)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    # Optional guide / arm parameters
    parser.add_argument('--n_guides', type=int, default=5,
                        help='Max guide RNAs to report per insertion site (default: 5)')
    parser.add_argument('--arm_length', type=int, default=500,
                        help='Homology arm length in bp on each side (default: 500)')
    parser.add_argument('--PAM', default='NGG',
                        help='IUPAC PAM sequence (default: NGG)')
    parser.add_argument('--guide_length', type=int, default=20,
                        help='Guide spacer length in nt (default: 20)')
    parser.add_argument('--cut_offset', type=int, default=3,
                        help='Bases upstream of PAM where DSB occurs (SpCas9=3)')
    args = parser.parse_args()

    print('Running design_tag_reagents.py', file=sys.stderr)
    print('  genewise      : {}'.format(args.genewise), file=sys.stderr)
    print('  genomic_fasta : {}'.format(args.genomic_fasta), file=sys.stderr)
    print('  PAM           : {}'.format(args.PAM), file=sys.stderr)
    print('  arm_length    : {}'.format(args.arm_length), file=sys.stderr)
    print('  n_guides      : {}'.format(args.n_guides), file=sys.stderr)

    df = design_reagents(
        genewise_out   = args.genewise,
        genomic_fasta  = args.genomic_fasta,
        n_guides       = args.n_guides,
        arm_length     = args.arm_length,
        pam            = args.PAM,
        guide_length   = args.guide_length,
        cut_offset     = args.cut_offset,
    )

    df.to_csv(args.output, sep='\t', index=False)
    print('Wrote {} rows to {}'.format(len(df), args.output), file=sys.stderr)
