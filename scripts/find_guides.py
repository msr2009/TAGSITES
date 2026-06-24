"""
find_guides.py

Find all potential CRISPR guide RNA cut sites in a genomic sequence (both strands).

This is a thin CLI around crispr_util.find_guides().  It does NOT perform
off-target (BLAST) filtering — the input is a local-region FASTA, not a whole
genome, so off-target analysis requires a separate whole-genome BLAST database.

Output TSV columns:
  strand          '+' or '-'
  spacer          guide spacer (5′→3′ on guide strand)
  pam_seq         PAM sequence on the guide strand
  pam_fwd_start   0-based start of the PAM in the fwd genomic sequence
  guide_fwd_start 0-based start of the spacer in fwd coords
  guide_fwd_end   0-based end (exclusive) of the spacer in fwd coords
  cut_pos         0-based fwd coord of the DSB (first base of right fragment)
  three_prime_C   True if the 3′ base of the spacer is C (lower-activity flag)

Matt Rich, 2025
"""

import sys
from Bio import SeqIO

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from crispr_util import find_guides


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Enumerate CRISPR guide cut sites in a genomic region (both strands).'
    )
    parser.add_argument('--genomic_fasta', '--input_file', required=True, dest='FASTA',
                        help='FASTA file of the genomic region')
    parser.add_argument('--PAM', default='NGG', dest='PAM',
                        help='IUPAC PAM sequence (default: NGG)')
    parser.add_argument('--guide_length', type=int, default=20, dest='GUIDE_LENGTH',
                        help='Guide spacer length in nt (default: 20)')
    parser.add_argument('--cut_offset', type=int, default=3, dest='CUT_OFFSET',
                        help='Bases upstream of PAM where cut occurs (SpCas9=3, default: 3)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    args = parser.parse_args()

    records = list(SeqIO.parse(args.FASTA, 'fasta'))
    if not records:
        sys.exit('ERROR: no sequences found in {}'.format(args.FASTA))
    seq = str(records[0].seq)

    guides = find_guides(seq,
                         pam=args.PAM,
                         guide_length=args.GUIDE_LENGTH,
                         cut_offset=args.CUT_OFFSET)

    print('Found {} potential guide sites (PAM={}, guide_length={})'.format(
        len(guides), args.PAM, args.GUIDE_LENGTH), file=sys.stderr)

    header = ['strand', 'spacer', 'pam_seq', 'pam_fwd_start',
              'guide_fwd_start', 'guide_fwd_end', 'cut_pos', 'three_prime_C']

    with open(args.output, 'w') as fh:
        print('\t'.join(header), file=fh)
        for g in guides:
            print('\t'.join(str(g[col]) for col in header), file=fh)

    print('Wrote guide table to {}'.format(args.output), file=sys.stderr)
