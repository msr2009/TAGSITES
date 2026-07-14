"""
junction_to_residue.py

Converts a short genomic/CDS junction sequence (e.g. a CRISPR knock-in cut site
or homology-arm sequence quoted in a paper/CGC page) into a residue position,
by locating it in a gene's coding sequence (CDS) and reporting which codon it
falls in.

Parameters:
    - cds_fasta (str): path to a FASTA file containing the gene's CDS (coding
      sequence only, e.g. NCBI's fasta_cds_na output for a RefSeq mRNA)
    - junction_seq (str): the short junction/homology sequence to locate

Returns:
    - prints the CDS nucleotide offset, residue number, and codon-frame status
      (whether the junction falls exactly on a codon boundary, i.e. a clean
      in-frame insertion site, or mid-codon)

Matt Rich, 7/2026
"""

from Bio import SeqIO
from Bio.Seq import Seq


def load_cds(cds_fasta):
    """Read the first record of a CDS FASTA file as an uppercase string."""
    rec = next(SeqIO.parse(cds_fasta, "fasta"))
    return str(rec.seq).upper()


def locate_junction(cds, junction_seq):
    """Find junction_seq in cds on either strand; returns (nt_index, strand) or (None, None)."""
    junction_seq = junction_seq.upper()
    rc = str(Seq(junction_seq).reverse_complement())
    for strand, seq in (("+", junction_seq), ("-", rc)):
        i = cds.find(seq)
        if i != -1:
            return i, strand
    return None, None


def nt_to_residue(nt_index):
    """Convert a 0-based CDS nucleotide offset to (codon_index, frame_offset).

    frame_offset == 0 means the junction starts exactly at a codon boundary
    (a clean in-frame insertion site); nonzero means it falls mid-codon.
    """
    return nt_index // 3, nt_index % 3


def main(cds_fasta, junction_seq):
    cds = load_cds(cds_fasta)
    protein = str(Seq(cds).translate(to_stop=True))

    nt_index, strand = locate_junction(cds, junction_seq)
    if nt_index is None:
        print(f"NOT FOUND: '{junction_seq}' does not occur in this CDS on either strand.")
        return None

    codon_index, frame_offset = nt_to_residue(nt_index)
    residue = codon_index + 1

    print(f"Found on {strand} strand at CDS nt {nt_index} (0-based)")
    print(f"Protein context: ...{protein[max(0, codon_index - 8):codon_index]} | "
          f"{protein[codon_index:codon_index + 8]}...")
    if frame_offset == 0:
        print(f"Clean in-frame insertion site: between residue {residue - 1} and residue {residue}")
    else:
        print(f"Falls mid-codon (offset {frame_offset}): within residue {residue}")
    return {"nt_index": nt_index, "strand": strand, "residue": residue, "frame_offset": frame_offset}


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--cds", action="store", type=str, dest="CDS_FASTA",
        help="path to gene CDS FASTA (e.g. NCBI fasta_cds_na output)", required=True)
    parser.add_argument("--junction", action="store", type=str, dest="JUNCTION_SEQ",
        help="junction/homology sequence to locate", required=True)

    args = parser.parse_args()

    main(args.CDS_FASTA, args.JUNCTION_SEQ)
