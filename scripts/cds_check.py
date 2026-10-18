"""
cds_check.py

Validate a gene model (Genewise output or a user GenBank) against the query protein.

Residue numbers in the reagent table come from translating the model's CDS, so any indel
between model and query silently shifts every downstream residue (DBL-1: two unspliced short
introns added 42 residues and Q239 was reported as F239). The model is globally aligned to the
query and the alignment is parsed. Substitutions and gaps (model-only or query-only residues) are
warnings that list the positions; each is judged against its own limit (MAX_SUBSTITUTION_FRACTION,
MAX_GAP_FRACTION of the query length) and either one over its limit fails, so a model with a gap
cannot hide a divergent stretch (unc-5: a 1-residue gap plus 23 substitutions in 947 residues).
With gaps the result carries a residue_map (model position -> query position, None for
model-only residues) so reagent numbering stays in input-protein coordinates; query residues
the model lacks have no genomic codon and get no reagents.
"""

import json
import os

MODEL_MISMATCH_HELP = (
    "Confirm the genomic region contains the gene isoform of the input protein, or upload a "
    "GenBank (.gb) file of the genomic region with the CDS (or exon) features annotated "
    "to replace the Genewise gene model."
)
MAX_SUBSTITUTION_FRACTION = 0.02  # substitutions above this share of the query length fail
MAX_GAP_FRACTION = 0.02  # gap residues (model-only + query-only) above this share fail
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


def build_residue_map(aligned_query, aligned_model):
    """Map each model residue (1-based) to its query residue number, None if model-only."""
    residue_map = {}
    qpos = mpos = 0
    for q, m in zip(aligned_query, aligned_model):
        # advance each sequence only on its own non-gap columns
        if q != "-":
            qpos += 1
        if m != "-":
            mpos += 1
            residue_map[mpos] = qpos if q != "-" else None
    return residue_map


def compare_to_query(model_seq, query_seq):
    """Align model to query and judge the parsed result: 'ok', 'warn' or 'fail'."""
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
        "residue_map": None,
        "message": "",
    }
    qlen = max(len(query_seq), 1)
    gap_residues = sum(i["length"] for i in insertions) + sum(
        d["stop"] - d["start"] + 1 for d in deletions
    )
    # gaps and substitutions are judged independently: either one over its limit fails
    too_gappy = gap_residues / qlen > MAX_GAP_FRACTION
    sub_frac = len(mismatches) / qlen
    too_divergent = sub_frac > MAX_SUBSTITUTION_FRACTION
    messages = []
    if insertions or deletions:
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
        msg = "Gene model translates to {} residues, input protein has {}: {}.".format(
            len(model_seq), len(query_seq), "; ".join(parts)
        )
        if too_gappy:
            msg += " The gap limit is {:.0%} of the protein.".format(MAX_GAP_FRACTION)
        messages.append(msg)
    if mismatches:
        detail = "{} residue{} ({:.1%} of {}): {}".format(
            len(mismatches), "s" if len(mismatches) != 1 else "", sub_frac, len(query_seq),
            format_mismatches(mismatches),
        )
        msg = "Gene model differs from the input protein at {}.".format(detail)
        if too_divergent:
            # too divergent to be the same protein (wrong gene, isoform or strain)
            msg = msg[:-1] + "; the substitution limit is {:.0%}.".format(MAX_SUBSTITUTION_FRACTION)
        messages.append(msg)

    if too_gappy or too_divergent:
        result["status"] = "fail"
    elif insertions or deletions or mismatches:
        result["status"] = "warn"
        if insertions or deletions:
            # a gap shifts the model's numbering downstream, so renumber through the alignment;
            # JSON keys are strings and only residues whose number differs need an entry
            result["residue_map"] = {
                str(m): q
                for m, q in build_residue_map(aligned_query, aligned_model).items()
                if q != m
            }
            messages.append(
                "Residue numbers follow the input protein; query residues the model lacks "
                "get no reagents."
            )
    result["message"] = " ".join(messages)
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
