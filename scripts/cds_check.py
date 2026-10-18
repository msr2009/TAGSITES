"""
cds_check.py

Validate a gene model (UCSC, NCBI RefSeq, GFF3, Genewise or a user GenBank) against the query
protein.

Residue numbers in the reagent table come from translating the model's CDS, so any indel
between model and query silently shifts every downstream residue (DBL-1: two unspliced short
introns added 42 residues and Q239 was reported as F239). The model is globally aligned to the
query and the alignment is parsed into substitutions and gaps (model-only or query-only residues).

A model fails when its differences look like a wrong exon rather than genome variation:
  - two substitutions closer than CLUSTER_WINDOW residues (a different alternative exon: Sh
    P08510-1 on UCSC's alignment had E422D N423S A432G T449R, 0.6% of the protein)
  - any one gap run longer than MAX_GAP_RUN residues (an unspliced intron or a missed exon)
  - a substitution within GAP_SUB_WINDOW residues of a gap (a shifted exon boundary)
  - substitutions or gap residues above MAX_SUBSTITUTION_FRACTION / MAX_GAP_FRACTION
A sweep of 150 UniProt isoforms (38 genes, 6 species) found every wrong-exon model clustered,
while real genome/UniProt differences were isolated residues (Act5C A261S + I331V, 70 apart),
so a fraction limit alone passes wrong exons. Isolated substitutions and short gaps warn.
With gaps the result carries a residue_map (model position -> query position, None for
model-only residues) so reagent numbering stays in input-protein coordinates; query residues
the model lacks have no genomic codon and get no reagents. describe() renders any result as a
located, sequence-level report.
"""

import json
import os

MODEL_MISMATCH_HELP = (
    "To supply the gene model yourself, upload a GenBank file (.gb/.gbk/.genbank/.gbff) of the "
    "genomic region in the 'Genomic region' field, annotated with a CDS feature for this "
    "isoform (a join over its exons) — or exon features for this isoform only. It replaces the "
    "automatic model and is checked the same way."
)
MAX_SUBSTITUTION_FRACTION = 0.02  # substitutions above this share of the query length fail
MAX_GAP_FRACTION = 0.02  # gap residues (model-only + query-only) above this share fail
CLUSTER_WINDOW = 15  # two substitutions closer than this (in query residues) fail
MAX_GAP_RUN = 3  # a single insertion or deletion run longer than this fails
GAP_SUB_WINDOW = 10  # a substitution within this many residues of a gap fails
MAX_LISTED = 25  # mismatches / indels spelled out in a message before truncating
CONTEXT = 5  # residues of context either side of a region in describe()
WRAP = 60  # alignment columns per line in describe()
LOGGED_CANDIDATES = 3  # candidate models summarised in the log per source


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


def _plural(n, word):
    """'1 residue' / '2 residues'."""
    return "{} {}{}".format(n, word, "" if n == 1 else "s")


def _gap_residues(insertions, deletions):
    """Total model-only plus query-only residues."""
    return sum(i["length"] for i in insertions) + sum(d["stop"] - d["start"] + 1 for d in deletions)


def _failure_regions(mismatches, insertions, deletions):
    """Stretches whose differences look like a wrong exon: [{start, stop, text}] in query residues."""
    regions = []
    # substitutions closer than CLUSTER_WINDOW chain into one cluster
    groups = []
    for m in mismatches:
        if groups and m[0] - groups[-1][-1][0] < CLUSTER_WINDOW:
            groups[-1].append(m)
        else:
            groups.append([m])
    for g in groups:
        if len(g) < 2:
            continue
        lo, hi = g[0][0], g[-1][0]
        regions.append(
            {
                "start": lo,
                "stop": hi,
                "text": (
                    "Residues {}–{}: {} within {} residues ({}). Two or more substitutions this close "
                    "together suggest the model uses a different exon.".format(
                        lo, hi, _plural(len(g), "substitution"), hi - lo + 1, format_mismatches(g)
                    )
                ),
            }
        )
    # gap runs: (first query residue touched, last, description)
    gaps = [
        (
            i["after"],
            i["after"] + 1,
            i["length"],
            "model has {} after input residue {}".format(
                _plural(i["length"], "extra residue"), i["after"]
            ),
        )
        for i in insertions
    ]
    gaps += [
        (
            d["start"],
            d["stop"],
            d["stop"] - d["start"] + 1,
            "model lacks input residue {}".format(d["start"])
            if d["start"] == d["stop"]
            else "model lacks input residues {}–{} ({} residues)".format(
                d["start"], d["stop"], d["stop"] - d["start"] + 1
            ),
        )
        for d in deletions
    ]
    for lo, hi, n, text in gaps:
        if n > MAX_GAP_RUN:
            regions.append(
                {
                    "start": lo,
                    "stop": hi,
                    "text": (
                        f"Residues {lo}–{hi}: {text}; the limit is {MAX_GAP_RUN} per gap. A gap this long suggests an "
                        "unspliced intron or a missed exon."
                    ),
                }
            )
            continue
        # a short gap is acceptable unless a substitution sits beside it
        near = [m for m in mismatches if lo - GAP_SUB_WINDOW <= m[0] <= hi + GAP_SUB_WINDOW]
        if near:
            a, b = min(lo, near[0][0]), max(hi, near[-1][0])
            regions.append(
                {
                    "start": a,
                    "stop": b,
                    "text": (
                        "Residues {}–{}: {} next to {} ({}). A gap beside substitutions suggests a "
                        "shifted exon boundary.".format(
                            a, b, text, _plural(len(near), "substitution"), format_mismatches(near)
                        )
                    ),
                }
            )
    # overlapping regions describe one problem: merge them
    regions.sort(key=lambda r: r["start"])
    merged = []
    for r in regions:
        if merged and r["start"] <= merged[-1]["stop"] + CONTEXT:
            merged[-1]["stop"] = max(merged[-1]["stop"], r["stop"])
            merged[-1]["text"] += "\n" + r["text"]
        else:
            merged.append(dict(r))
    return merged


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
        "failures": [],
        "aligned_query": aligned_query,
        "aligned_model": aligned_model,
        "message": "",
    }
    qlen = max(len(query_seq), 1)
    failures = _failure_regions(mismatches, insertions, deletions)
    # the fraction limits still cap scattered differences
    gap_frac = _gap_residues(insertions, deletions) / qlen
    sub_frac = len(mismatches) / qlen
    if gap_frac > MAX_GAP_FRACTION:
        failures.append(
            {
                "start": None,
                "stop": None,
                "text": (
                    f"Gap residues are {gap_frac:.1%} of the protein; the limit is {MAX_GAP_FRACTION:.0%}."
                ),
            }
        )
    if sub_frac > MAX_SUBSTITUTION_FRACTION:
        # too divergent to be the same protein (wrong gene, isoform or strain)
        failures.append(
            {
                "start": None,
                "stop": None,
                "text": (
                    f"Substitutions are {sub_frac:.1%} of the protein; the limit is {MAX_SUBSTITUTION_FRACTION:.0%}."
                ),
            }
        )
    result["failures"] = failures

    if failures:
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
    result["message"] = describe(result)
    return result


def summary_line(check, label=""):
    """One-line header of a check: source, verdict, lengths and difference counts."""
    verb = {"ok": "matches", "warn": "differs from", "fail": "does not match"}[check["status"]]
    ins = sum(i["length"] for i in check["insertions"])
    dels = sum(d["stop"] - d["start"] + 1 for d in check["deletions"])
    # a run count only where there are gap residues
    runs = lambda n, k: _plural(n, k[0]) + (" in " + _plural(k[1], "gap") if n else "")
    return "Gene model{} {} the input protein ({} vs {} residues): {}, {}, {}.".format(
        " (" + label + ")" if label else "",
        verb,
        check["model_len"],
        check["query_len"],
        _plural(len(check["mismatches"]), "substitution"),
        runs(ins, ("inserted residue", len(check["insertions"]))),
        runs(dels, ("missing residue", len(check["deletions"]))),
    )


def _alignment_block(aligned_query, aligned_model, start, stop):
    """Aligned input/model rows over query residues start..stop (plus CONTEXT), with markers."""
    # each column's query position (the last one consumed, for gap columns)
    cols, qpos = [], 0
    for q in aligned_query:
        if q != "-":
            qpos += 1
        cols.append(qpos)
    lo, hi = start - CONTEXT, stop + CONTEXT
    idx = [i for i, p in enumerate(cols) if lo <= p <= hi]
    if not idx:
        return ""
    c0, c1 = idx[0], idx[-1]
    # residue numbers at the block start for each row
    q_before = sum(1 for c in aligned_query[:c0] if c != "-")
    m_before = sum(1 for c in aligned_model[:c0] if c != "-")
    # a long region shows its first and last line only
    starts = list(range(c0, c1 + 1, WRAP))
    if len(starts) > 2:
        starts = [starts[0], starts[-1]]
    lines = []
    prev = None  # start column of the previous line shown
    for w in starts:
        if prev is not None:
            skipped = w - (prev + WRAP)
            if skipped > 0:
                lines.append(f"  ... {skipped} alignment columns omitted ...")
                q_before += sum(1 for c in aligned_query[prev + WRAP : w] if c != "-")
                m_before += sum(1 for c in aligned_model[prev + WRAP : w] if c != "-")
        prev = w
        qseg = aligned_query[w : min(w + WRAP, c1 + 1)]
        mseg = aligned_model[w : min(w + WRAP, c1 + 1)]
        # * substitution, ^ gap in either row
        marks = "".join(
            "^" if "-" in (a, b) else "*" if a != b else " " for a, b in zip(qseg, mseg)
        )
        qn = sum(1 for c in qseg if c != "-")
        mn = sum(1 for c in mseg if c != "-")
        lines.append(f"  input {q_before + 1:>5} {qseg} {q_before + qn}")
        lines.append(f"  model {m_before + 1:>5} {mseg} {m_before + mn}")
        lines.append("  {:>11} {}".format("", marks.rstrip()))
        q_before += qn
        m_before += mn
    return "\n".join(lines)


def describe(check, label=""):
    """Located, sequence-level report of a check: header, failing regions, other differences."""
    if check["status"] == "ok":
        return summary_line(check, label)
    lines = [summary_line(check, label)]
    failures = check.get("failures") or []
    regions = [f for f in failures if f["start"] is not None]
    aq, am = check.get("aligned_query"), check.get("aligned_model")
    for r in regions[:MAX_LISTED]:
        lines.append(r["text"])
        if aq and am:
            lines.append(_alignment_block(aq, am, r["start"], r["stop"]))
    if regions and aq and am:
        lines.append("(* substitution, ^ gap)")
    if len(regions) > MAX_LISTED:
        lines.append(f"... and {len(regions) - MAX_LISTED} more failing region(s).")
    lines += [f["text"] for f in failures if f["start"] is None]
    # differences outside the failing regions
    inside = lambda p: any(r["start"] <= p <= r["stop"] for r in regions)
    subs = [m for m in check["mismatches"] if not inside(m[0])]
    gaps = [
        "{} extra after {}".format(i["length"], i["after"])
        for i in check["insertions"]
        if not inside(i["after"])
    ]
    gaps += [
        "missing {}".format(
            d["start"] if d["start"] == d["stop"] else "{}–{}".format(d["start"], d["stop"])
        )
        for d in check["deletions"]
        if not inside(d["start"])
    ]
    if subs or gaps:
        parts = []
        if subs:
            parts.append(
                "{} ({})".format(
                    _plural(len(subs), "isolated substitution"), format_mismatches(subs)
                )
            )
        if gaps:
            parts.append(
                "{} ({})".format(_plural(len(gaps), "short gap"), "; ".join(gaps[:MAX_LISTED]))
            )
        lines.append(
            ("Other differences: " if regions else "Differences: ") + ", ".join(parts) + "."
        )
    if check["status"] == "warn" and (check["insertions"] or check["deletions"]):
        lines.append(
            "Residue numbers follow the input protein; input residues the model lacks "
            "get no reagents."
        )
    return "\n".join(line for line in lines if line)


def format_mismatches(mismatches):
    """Format mismatches as 'Q239>F, ...' (query residue, position, model residue)."""
    shown = ", ".join(f"{q}{pos}>{m}" for pos, q, m in mismatches[:MAX_LISTED])
    extra = len(mismatches) - MAX_LISTED
    return shown + (f", ... (+{extra} more)" if extra > 0 else "")


def check_model(cds_df, dna, query_seq, label=""):
    """Check a gene model against the query; raise ValueError on 'fail', else return result."""
    result = compare_to_query(translate_model(cds_df, dna), query_seq)
    if label:
        result["message"] = describe(result, label)
    if result["status"] == "fail":
        raise ValueError("{}\n{}".format(result["message"], MODEL_MISMATCH_HELP))
    return result


def rank(check):
    """Sort key for candidate models: ok, then warn, then fail; fewest differences within each."""
    order = {"ok": 0, "warn": 1, "fail": 2}[check["status"]]
    return (
        order,
        len(check["mismatches"]) + _gap_residues(check["insertions"], check["deletions"]),
        abs(check["model_len"] - check["query_len"]),
    )


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
