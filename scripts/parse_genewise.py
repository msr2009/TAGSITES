"""
parse_genewise.py

Parse Genewise alignment output and enumerate all possible tag insertion sites
(codon boundaries in the genomic sequence, one per protein residue).

Based on math prototyped in FORTESTING_DELETEWHENREADY/parse_genewise.ipynb.

Key insight (from notebook cell 5):
  The Genewise GFF 'frame' column is 1-indexed in the file; we convert to 0-indexed.
  'frame' here means: the number of bases in the CURRENT exon that complete the
  split codon from the previous exon (i.e., the GFF3 'phase').
  frame=0  → full codon starts at exon start; remap to frame_prime=3 for math
  frame=1  → 1 base in this exon finishes the previous split codon
  frame=2  → 2 bases in this exon finish the previous split codon

  tag_sites per exon = range(exon_start + frame_prime, exon_stop + 3, 3)
  Each value x is the genomic insert position (base after the 3′ end of a codon).
  Codon sequence = dna[x-3 : x].

Output TSV columns:
  residue_index   1-based protein residue number
  amino_acid      single-letter code (X if codon spans intron = split codon)
  insert_pos      0-based genomic position where tag inserts (between x-1 and x)
  exon_index      0-based exon number
  is_split_codon  whether this codon straddles an intron boundary
  dist_to_5p_splice  bases from insert_pos to 5′ end of current exon
  dist_to_3p_splice  bases from insert_pos to 3′ end of current exon last codon
                     (set to -1 for the last exon = no downstream intron)

Matt Rich, 2025
"""

import re
import sys
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, __file__.rsplit('/', 1)[0])  # allow import from scripts/
from crispr_util import CODONS


# ── Score extraction ──────────────────────────────────────────────────────────

def parse_genewise_score(out_txt):
    """
    Extract the top-level alignment score from a Genewise .out.txt file.

    Returns the float score from the line:
        Score NNN.NN bits over entire alignment

    Returns None if the line is not found.
    NOTE: the EBI REST API does not write this header line; use
    parse_genewise_gff_score() instead for reliable score extraction.
    """
    with open(out_txt) as fh:
        for line in fh:
            m = re.match(r'^Score\s+([\d.]+)\s+bits', line)
            if m:
                return float(m.group(1))
    return None


def parse_genewise_gff_score(out_txt):
    """
    Extract the alignment score from the GFF 'match' row score column.

    The EBI REST Genewise API does not write a 'Score NNN bits' header line,
    but the score is present in the GFF section as the score of the 'match'
    feature (column 6).  This function works for both EBI output and the
    local standalone format; parse_genewise_score() only works for standalone.

    Returns float score, or None if no match row is found.
    """
    with open(out_txt) as fh:
        content = fh.read()
    for section in content.split('//'):
        if '\tGeneWise\t' not in section:
            continue
        for line in section.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 9 and parts[2] == 'match':
                try:
                    return float(parts[5])
                except ValueError:
                    pass
    return None


def cds_coverage(cds_df, protein_length):
    """
    Estimate what fraction of the protein is covered by the Genewise CDS exons.

    Computes sum-of-CDS-nucleotides / 3 / protein_length.  Values near 1.0
    indicate a good alignment; values near 0 indicate a garbage / wrong-strand
    result.  protein_length = 0 returns 0.0.
    """
    if protein_length <= 0 or cds_df.empty:
        return 0.0
    total_nt = sum(int(row['stop']) - int(row['start']) + 1
                   for _, row in cds_df.iterrows())
    return (total_nt / 3) / protein_length


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_genewise(out_txt):
    """
    Parse a Genewise .out.txt file and return a DataFrame of CDS exons.

    Columns: name, source, type, start (0-indexed), stop (0-indexed),
             score, strand, frame (as int), note
    Only rows with type == 'cds' are returned, sorted by start.

    The GFF block is identified by the presence of '\\tGeneWise\\t' — this is
    more robust than the hardcoded index-7 approach in the prototype notebook.
    """
    with open(out_txt) as fh:
        content = fh.read()

    # Find the //…// section that contains GeneWise GFF rows
    gff_section = None
    for section in content.split('//'):
        if '\tGeneWise\t' in section:
            gff_section = section
            break

    if gff_section is None:
        raise ValueError(
            'No GFF section (containing \\tGeneWise\\t) found in: {}'.format(out_txt)
        )

    rows = []
    for line in gff_section.strip().split('\n'):
        line = line.strip()
        if not line or '\tGeneWise\t' not in line:
            continue
        parts = line.split('\t')
        if len(parts) < 9:
            continue
        rows.append({
            'name':   parts[0],
            'source': parts[1],
            'type':   parts[2],
            'start':  int(parts[3]) - 1,   # 1-indexed → 0-indexed
            'stop':   int(parts[4]) - 1,   # 1-indexed → 0-indexed
            'score':  float(parts[5]) if parts[5] not in ('.', '') else 0.0,
            'strand': parts[6],
            'frame':  int(parts[7]) if parts[7] not in ('.', '') else 0,
            'note':   parts[8],
        })

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError('GFF section found but no rows parsed from: {}'.format(out_txt))

    cds = df[df['type'] == 'cds'].copy().sort_values('start').reset_index(drop=True)
    return cds


# ── Insertion-site enumeration ─────────────────────────────────────────────────

def enumerate_insertion_sites(cds_df, dna):
    """
    Map every protein residue to a genomic tag-insertion position.

    Parameters
    ----------
    cds_df : pd.DataFrame  output of parse_genewise()
    dna    : str           genomic sequence (forward strand, 0-indexed)

    Returns
    -------
    pd.DataFrame with one row per residue, columns as described in the module
    docstring.  Also returns a second item: the CDS DataFrame sorted by exon
    start (same as input, provided for convenience).
    """
    dna_upper = dna.upper()
    n_exons = len(cds_df)
    rows = []
    residue_idx = 0

    for exon_i, exon in cds_df.iterrows():
        start = int(exon['start'])
        stop  = int(exon['stop'])
        frame = int(exon['frame'])
        # frame_prime = offset to the END of the first codon in this exon
        # (= 3 when frame=0, meaning the exon starts a clean codon)
        frame_prime = 3 if frame == 0 else frame

        is_last_exon = (exon_i == n_exons - 1)

        for x in range(start + frame_prime, stop + 3, 3):
            codon_start = x - 3

            # A codon is split if it crosses either the 5′ or 3′ exon boundary:
            #   5′ split: codon starts before the exon (codon_start < exon start)
            #     → "completing split": the previous exon left a partial codon; this
            #       is the valid exonic insertion site for that split codon boundary.
            #   3′ split: last codon base (x-1) is past the last exon base (stop)
            #     → "originating split": the codon extends into the intron; skip it.
            #       The next exon's completing split is the valid insertion site.
            if x - 1 > stop:
                continue   # intronic insertion; next exon handles this boundary

            residue_idx += 1
            is_split = codon_start < start

            if is_split:
                aa = 'X'
            else:
                codon = dna_upper[codon_start:x]
                aa = CODONS.get(codon, 'X')

            # Distances to splice sites
            # dist_to_5p_splice: bases from exon start to this insert position
            #   (small → close to 5′ splice site; =3 for first clean codon)
            dist_5p = x - start

            # dist_to_3p_splice: bases from insert_pos to end of exon last codon.
            #   Value of 0 means tag inserts exactly at the 3′ exon boundary.
            #   Negative means the insert is PAST the boundary (3′ split codon).
            #   99999 is the sentinel for the last exon (no downstream intron).
            if is_last_exon:
                dist_3p = 99999
            else:
                dist_3p = (stop + 1) - x   # 0 when x = stop+1; negative for 3′ splits

            rows.append({
                'residue_index':    residue_idx,
                'amino_acid':       aa,
                'insert_pos':       x,        # 0-based; insert between x-1 and x
                'exon_index':       int(exon_i),
                'is_split_codon':   is_split,
                'dist_to_5p_splice': dist_5p,
                'dist_to_3p_splice': dist_3p,
            })

    return pd.DataFrame(rows)


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description=(
            'Parse Genewise output and enumerate tag insertion sites '
            '(one per protein residue).'
        )
    )
    parser.add_argument('--genewise', required=True,
                        help='Genewise .out.txt result file')
    parser.add_argument('--fasta', required=True,
                        help='Genomic region FASTA (same sequence submitted to Genewise)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    parser.add_argument('--verify', action='store_true',
                        help='Print codon-verification table to stderr for spot-checking')
    args = parser.parse_args()

    # Load genomic sequence
    records = list(SeqIO.parse(args.fasta, 'fasta'))
    if not records:
        sys.exit('ERROR: no sequences found in {}'.format(args.fasta))
    dna = str(records[0].seq)

    # Parse Genewise
    cds_df = parse_genewise(args.genewise)
    print('Parsed {} CDS exons from {}'.format(len(cds_df), args.genewise),
          file=sys.stderr)
    print(cds_df.to_string(), file=sys.stderr)

    # Enumerate insertion sites
    sites = enumerate_insertion_sites(cds_df, dna)
    sites.to_csv(args.output, sep='\t', index=False)
    print('Wrote {} insertion sites to {}'.format(len(sites), args.output),
          file=sys.stderr)

    if args.verify:
        print('\nCodon verification (first 10 residues):', file=sys.stderr)
        for _, row in sites.head(10).iterrows():
            x = row['insert_pos']
            codon = dna[x - 3:x].upper()
            print('  residue {:4d}  insert_pos {:6d}  codon {}  aa {}{}'.format(
                row['residue_index'], x, codon, row['amino_acid'],
                '  [SPLIT]' if row['is_split_codon'] else ''),
                file=sys.stderr)
