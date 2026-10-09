"""
genbank_input.py

Use a user-supplied GenBank file (genomic region + annotated CDS or exon features) as the gene
model in place of Genewise. Writes the same three files genewise_remote.main() produces
({outprefix}.genewise.out.txt, .genewise_genomic.fa, .genewise_orientation.txt) so reagent
design is unchanged, after checking the model's translation against the input protein with the
same rule as a Genewise model (cds_check.check_model: length difference fails).

A GenBank file with no CDS/exon features is only a genomic sequence: it is written as
{outprefix}.genomic.fa and the caller falls through to Genewise.

CLI:
  python genbank_input.py --genbank region.gb --protein_fasta prot.fa --outprefix data/X/X_genewise
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cds_check
from parse_genewise import GROUND_TRUTH_SCORE, write_genewise_gff
from progress import report as _report, resolve_reporter
from site_selection_util import save_fasta

GENBANK_EXTENSIONS = (".gb", ".gbk", ".genbank", ".gbff")
MODEL_FEATURES = ("CDS", "exon")


def is_genbank(path):
    """True if the path has a GenBank file extension."""
    return str(path).lower().endswith(GENBANK_EXTENSIONS)


def read_record(gb_path):
    """Return the first record of a GenBank file."""
    from Bio import SeqIO

    return next(SeqIO.parse(gb_path, "genbank"))


def has_gene_model(gb_path):
    """True if the first record has a CDS or exon feature."""
    return any(f.type in MODEL_FEATURES for f in read_record(gb_path).features)


def genbank_to_fasta(gb_path, fasta_path):
    """Write the record's sequence as a FASTA (for a GenBank with no gene model); return path."""
    rec = read_record(gb_path)
    save_fasta(rec.id, str(rec.seq).upper(), fasta_path)
    return fasta_path


def _feature_positions(feature, trim_first):
    """Genomic positions (0-based) of a feature's bases in transcript order, and its strand."""
    strand = feature.location.strand or 1
    positions = []
    for part in sorted(feature.location.parts, key=lambda p: int(p.start)):
        positions.extend(range(int(part.start), int(part.end)))
    if strand == -1:
        positions.reverse()  # transcript order runs from the high coordinate down
    return positions[trim_first:], strand


def _candidates(rec):
    """Yield (label, positions, strand): every CDS feature, else all exons as one model."""
    cds = [f for f in rec.features if f.type == "CDS"]
    for i, f in enumerate(cds, start=1):
        # codon_start (1-3) marks the first complete codon of a CDS that begins mid-codon
        trim = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1
        positions, strand = _feature_positions(f, trim)
        label = f.qualifiers.get("gene", f.qualifiers.get("locus_tag", [str(i)]))[0]
        yield f"CDS {label}", positions, strand
    exons = [f for f in rec.features if f.type == "exon"]
    if exons and not cds:
        strand = exons[0].location.strand or 1
        positions = sorted({p for f in exons for p in range(int(f.location.start),
                                                            int(f.location.end))})
        if strand == -1:
            positions.reverse()
        yield "exon features", positions, strand


def _trim_to_query(positions, odna, query_seq):
    """Trim an exon model (which may include UTR) to the query's reading frame, if it is found."""
    transcript = "".join(odna[p] for p in positions)
    from Bio.Seq import Seq

    for frame in range(3):
        usable = (len(transcript) - frame) // 3 * 3
        protein = str(Seq(transcript[frame : frame + usable]).translate())
        hit = protein.find(query_seq.rstrip("*"))
        if hit >= 0:
            start = frame + 3 * hit
            return positions[start : start + 3 * len(query_seq.rstrip("*"))]
    return positions


def _spans_df(positions):
    """Collapse transcript-ordered positions into exon rows (0-based start/stop, GFF phase)."""
    rows = []
    run_start = prev = positions[0]
    for p in positions[1:]:
        if p == prev + 1:
            prev = p
        else:
            rows.append([run_start, prev])
            run_start = prev = p
    rows.append([run_start, prev])
    df = pd.DataFrame(rows, columns=["start", "stop"])
    # phase = bases at the start of this exon that finish the previous exon's split codon
    before = (df["stop"] - df["start"] + 1).cumsum().shift(fill_value=0)
    df["frame"] = (3 - before % 3) % 3
    return df


def _rank(check):
    """Sort key for candidate models: exact match, then fewest substitutions, then length gap."""
    if check["status"] == "ok":
        return (0, 0)
    if check["status"] == "warn":
        return (1, len(check["mismatches"]))
    return (2, abs(check["model_len"] - check["query_len"]))


def genbank_to_genewise(gb_path, protein_fasta, outprefix, report=None):
    """Convert a GenBank gene model to the Genewise output files; raise ValueError if it fails."""
    reporter = resolve_reporter(report)
    query_seq = cds_check.read_protein(protein_fasta)
    rec = read_record(gb_path)
    dna = str(rec.seq).upper()
    flipped = str(rec.seq.reverse_complement()).upper()  # tolerates IUPAC codes

    best = None
    for label, positions, strand in _candidates(rec):
        # downstream code assumes the coding strand is "+": a minus-strand model is mapped
        # onto the reverse-complemented sequence
        odna = dna if strand == 1 else flipped
        if strand == -1:
            positions = [len(dna) - 1 - p for p in positions]
        if query_seq:
            positions = _trim_to_query(positions, odna, query_seq) if label == "exon features" \
                else positions
        if not positions:
            continue
        df = _spans_df(positions)
        check = cds_check.compare_to_query(cds_check.translate_model(df, odna), query_seq)
        _report(reporter, f"GenBank {label}: {check['status']} "
                          f"({check['model_len']} vs {check['query_len']} residues)",
                stage="genbank_input")
        if best is None or _rank(check) < _rank(best[0]):
            best = (check, df, odna, label)

    if best is None:
        raise ValueError("No CDS or exon features found in {}. {}".format(
            os.path.basename(gb_path), cds_check.MODEL_MISMATCH_HELP))
    check, df, odna, label = best
    if check["status"] == "fail":
        raise ValueError("GenBank gene model ({}): {} {}".format(
            label, check["message"], cds_check.MODEL_MISMATCH_HELP))
    if check["status"] == "warn":
        _report(reporter, check["message"], stage="genbank_input", level="warning")

    out_txt = f"{outprefix}.genewise.out.txt"
    genomic_fa = f"{outprefix}.genewise_genomic.fa"
    save_fasta(rec.id, odna, genomic_fa)
    write_genewise_gff(out_txt, df, rec.id, f"genbank_input.py ({label})")
    with open(f"{outprefix}.genewise_orientation.txt", "w") as fh:
        fh.write("+\n")  # the written sequence is already on the coding strand
    _report(reporter, f"GenBank gene model ({label}, {len(df)} exons) → {out_txt}",
            stage="genbank_input")
    return {
        "orientation": "+",
        "out_txt": out_txt,
        "genomic_fa": genomic_fa,
        "winner_score": GROUND_TRUTH_SCORE,
        "warning": check["message"] or None,
    }


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Use a GenBank gene model in place of Genewise.")
    parser.add_argument("--genbank", required=True, help="GenBank file: genomic region + CDS/exons")
    parser.add_argument("--protein_fasta", required=True, help="Input protein FASTA (or PDB)")
    parser.add_argument("--outprefix", required=True, help="Prefix of the Genewise-style outputs")
    args, _ = parser.parse_known_args()

    if has_gene_model(args.genbank):
        genbank_to_genewise(args.genbank, args.protein_fasta, args.outprefix)
    else:
        # sequence only: caller runs Genewise on this FASTA
        fa = genbank_to_fasta(args.genbank, f"{args.outprefix}.genomic.fa")
        print(f"No CDS/exon features in {args.genbank}; sequence written to {fa}", file=sys.stderr)
