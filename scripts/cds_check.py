"""
cds_check.py

Validate a gene model (Genewise output or a user GenBank) against the query protein.

Residue numbers in the reagent table come from translating the model's CDS, so any indel
between model and query silently shifts every downstream residue (DBL-1: two unspliced short
introns added 42 residues and Q239 was reported as F239). The model is globally aligned to the
query and the alignment is parsed: any gap (model-only or query-only residues) is a hard
failure, even when the lengths happen to agree. Substitutions alone are a warning that lists the
positions (the reagents are still built on the right codons) up to MAX_SUBSTITUTION_FRACTION of
the query length, and a failure above it.
"""

import json
import os

MODEL_MISMATCH_HELP = (
    "Confirm the genomic region contains the gene isoform of the input protein, or upload a "
    "GenBank (.gb) file of the genomic region with the CDS (or exon) features annotated "
    "to replace the Genewise gene model."
)
MAX_SUBSTITUTION_FRACTION = 0.02  # substitutions above this share of the query length fail
MAX_LISTED = 25  # mismatches / indels spelled out in a message before truncating


def read_protein(path):
    """Return the first protein sequence of a FASTA (or PDB) file as a string, '' if none."""
    if not path or not os.path.exists(path):
        return ""
    if path.endswith(".pdb"):
        from site_selection_util import get_sequence

        return str(get_sequence(path))
    from Bio import SeqIO

    recs = list(SeqIO.parse(path, "fasta"))
    return str(recs[0].seq) if recs else ""


def translate_model(cds_df, dna):
    """Translate the joined CDS spans (0-based inclusive start/stop) to a protein string."""
    from Bio.Seq import Seq

    cds = "".join(dna[int(r.start) : int(r.stop) + 1] for r in cds_df.itertuples())
    cds = cds[: len(cds) // 3 * 3]  # drop a trailing partial codon
    protein = str(Seq(cds).translate())
    # terminal stop is not a residue; an internal stop is selenocysteine (U)
    return protein.rstrip("*").replace("*", "U")


def align_to_query(query_seq, model_seq):
    """Globally align the model to the query; return aligned rows (query, model) as strings."""
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "global"  # whole-protein comparison: a local alignment would hide bad termini
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    # BLOSUM62 has no selenocysteine (U); U -> C keeps every position unchanged
    alignment = aligner.align(query_seq.replace("U", "C"), model_seq.replace("U", "C"))[0]
    return alignment[0], alignment[1]


def parse_alignment(aligned_query, aligned_model, query_seq, model_seq):
    """Walk the alignment columns into substitutions, model-only insertions, query-only deletions."""
    mismatches, insertions, deletions = [], [], []
    qpos = mpos = 0  # residues consumed so far (1-based after increment)
    for q, m in zip(aligned_query, aligned_model):
        if q != "-" and m != "-":
            qpos += 1
            mpos += 1
            # compare the original letters, not the U->C alignment surrogates
            if query_seq[qpos - 1] != model_seq[mpos - 1]:
                mismatches.append([qpos, query_seq[qpos - 1], model_seq[mpos - 1]])
        elif m != "-":
            # model residue with no query partner: extra sequence (e.g. an unspliced intron)
            mpos += 1
            if insertions and insertions[-1]["after"] == qpos:
                insertions[-1]["length"] += 1
            else:
                insertions.append({"after": qpos, "length": 1})
        else:
            # query residue the model lacks (e.g. a missed exon)
            qpos += 1
            if deletions and deletions[-1]["stop"] == qpos - 1:
                deletions[-1]["stop"] = qpos
            else:
                deletions.append({"start": qpos, "stop": qpos})
    return mismatches, insertions, deletions


def compare_to_query(model_seq, query_seq):
    """Align model to query and judge the parsed result: 'ok', 'warn' (substitutions only) or 'fail'."""
    query_seq = query_seq.rstrip("*")
    aligned_query, aligned_model = align_to_query(query_seq, model_seq)
    mismatches, insertions, deletions = parse_alignment(
        aligned_query, aligned_model, query_seq, model_seq
    )
    result = {
        "status": "ok",
        "query_len": len(query_seq),
        "model_len": len(model_seq),
        "mismatches": mismatches,
        "insertions": insertions,
        "deletions": deletions,
        "message": "",
    }
    if insertions or deletions:
        # any gap shifts the numbering of every downstream residue, even when lengths agree
        result["status"] = "fail"
        parts = [
            "model has {} extra residue{} after query residue {}".format(
                i["length"], "s" if i["length"] != 1 else "", i["after"]
            )
            for i in insertions[:MAX_LISTED]
        ] + [
            "model lacks query residues {}".format(
                d["start"] if d["start"] == d["stop"] else "{}-{}".format(d["start"], d["stop"])
            )
            for d in deletions[:MAX_LISTED]
        ]
        result["message"] = "Gene model translates to {} residues, input protein has {}: {}.".format(
            len(model_seq), len(query_seq), "; ".join(parts)
        )
    elif mismatches:
        frac = len(mismatches) / max(len(query_seq), 1)
        detail = "{} residue{} ({:.1%} of {}): {}".format(
            len(mismatches), "s" if len(mismatches) != 1 else "", frac, len(query_seq),
            format_mismatches(mismatches),
        )
        if frac > MAX_SUBSTITUTION_FRACTION:
            # too divergent to be the same protein (wrong gene, isoform or strain)
            result["status"] = "fail"
            result["message"] = "Gene model differs from the input protein at {}; the limit is " \
                "{:.0%}.".format(detail, MAX_SUBSTITUTION_FRACTION)
        else:
            result["status"] = "warn"
            result["message"] = "Gene model differs from the input protein at {}.".format(detail)
    return result


def format_mismatches(mismatches):
    """Format mismatches as 'Q239>F, ...' (query residue, position, model residue)."""
    shown = ", ".join("{}{}>{}".format(q, pos, m) for pos, q, m in mismatches[:MAX_LISTED])
    extra = len(mismatches) - MAX_LISTED
    return shown + (", ... (+{} more)".format(extra) if extra > 0 else "")


def check_model(cds_df, dna, query_seq):
    """Check a gene model against the query; raise ValueError on 'fail', else return result."""
    result = compare_to_query(translate_model(cds_df, dna), query_seq)
    if result["status"] == "fail":
        raise ValueError("{} {}".format(result["message"], MODEL_MISMATCH_HELP))
    return result


def write_sidecar(path, result):
    """Write a check result as JSON for the reagents tab."""
    with open(path, "w") as fh:
        json.dump(result, fh)


def read_sidecar(path):
    """Return a check result written by write_sidecar, or None if absent or unreadable."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None
