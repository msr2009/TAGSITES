"""ncbi_gene_model.py — gene model candidates from every NCBI RefSeq transcript of the gene.

Second source for the reagents gene model (see gene_model.py), tried when no UCSC UniProt
alignment translates exactly. UCSC lifts a UniProt isoform through one RefSeq transcript and
sometimes picks the wrong one (Sh P08510-1: the mutually exclusive exon at chrX:17946598 instead
of 17945218); here every protein-coding RefSeq transcript of the gene, curated (NM_) and
predicted (XM_), is translated from the same genome and compared with the input protein.
In a sweep of 150 UniProt isoforms this recovered an exact model for all 21 isoforms whose UCSC
alignment used the wrong exon.

  UniProt entry                      accession -> NCBI GeneID (GeneID cross-reference)
  Datasets /gene/id/{id}             reference genomic accession (NC_...), span and assembly
  efetch nuccore fasta               one plus-strand sequence: the gene plus flank
  efetch nuccore ft                  the annotation's CDS features there (genomic intervals)
  Datasets /gene/id/{id}/product_report   protein -> transcript name and RefSeq Select flag

The CDS comes from the genome annotation's own intervals, not by mapping the mRNA's CDS range
through its exons: a RefSeq mRNA can differ from (or align only partly to) the genome (chicken
TP53 NM_205264.1: 1555 nt, 965 nt aligned), which would shift the frame.

Works for any organism NCBI annotates, including those without a UCSC assembly. Candidates use
the model dict of ucsc_gene_model.item_to_model and are written by gene_model.write_model;
region.json carries the NC_ accession as chrom (phastcons_track aliases it to the UCSC name and
refuses a track whose codons do not read back identically).

CLI writes a model only on an exact translation: exit 0 when written, 2 otherwise.

Matt Rich, 2026
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import cds_check
from crispr_util import reverse_complement
from parse_genewise import spans_from_positions
from progress import report as _report, resolve_reporter

STOP_CODONS = ("TAA", "TAG", "TGA")


def gene_id_for(accession):
    """NCBI GeneID from a UniProt entry's cross-references, or '' if it has none."""
    from uniprot_api import _gene_id_xref, fetch_entry

    return _gene_id_xref(fetch_entry(accession.split("-")[0])) or ""


def _cds_positions(intervals):
    """1-based genomic positions of a CDS feature's bases in transcript order, and its strand."""
    # minus-strand features list each interval high -> low
    strand = "-" if any(b > e for b, e in intervals) else "+"
    out = []
    for b, e in intervals:
        out.extend(range(b, e + 1) if strand == "+" else range(b, e - 1, -1))
    return out, strand


def _protein_id(qualifiers):
    """The RefSeq protein accession of a CDS feature ('ref|NP_728120.1|' -> 'NP_728120.1')."""
    ids = qualifiers.get("protein_id") or [""]
    ref = [i for i in ids if i.startswith("ref|")]
    return (ref or ids)[0].replace("ref|", "").strip("|")


def cds_to_model(intervals, codon_start, chrom, region, flank_bp):
    """One annotated CDS -> the model dict of ucsc_gene_model.item_to_model.

    region is (lo1, plus-strand sequence) covering the CDS; the model's dna is the CDS span
    plus flank_bp on the coding strand, with the stop codon left out of the CDS (the Genewise
    convention the reagent design expects).
    """
    rlo, rseq = region
    rhi = rlo + len(rseq) - 1
    gpos, strand = _cds_positions(intervals)
    gpos = gpos[codon_start - 1:]  # a 5'-partial CDS starts mid-codon
    start, stop = min(gpos), max(gpos)
    tlo, thi = max(rlo, start - flank_bp), min(rhi, stop + flank_bp)
    dna = rseq[tlo - rlo : thi - rlo + 1]
    if strand == "-":
        dna = reverse_complement(dna)
        idx = [thi - p for p in gpos]
    else:
        idx = [p - tlo for p in gpos]
    # drop a terminal stop codon, as UCSC and Genewise models exclude it
    if len(idx) >= 3 and "".join(dna[i] for i in idx[-3:]) in STOP_CODONS:
        idx = idx[:-3]
    low, high = start - tlo, thi - stop  # genomic-low / genomic-high flank actually obtained
    return {
        "cds_df": spans_from_positions(idx), "dna": dna,
        "chrom": chrom, "start": start, "stop": stop, "strand": strand,
        "flank_5p": low if strand == "+" else high,
        "flank_3p": high if strand == "+" else low,
    }


def _preference(c):
    """Tie-break among equally good candidates: RefSeq Select/MANE, then NM_, then longest span."""
    m = c["model"]
    return (not c["select"], c["name"].startswith("XM_"), -(m["stop"] - m["start"]), c["name"])


def evaluate(protein_fa, taxid, accession="", flank_bp=2000, report=None, deadline_s=None):
    """Every annotated RefSeq CDS of the protein's gene as a candidate model, best first.

    Returns (candidates, gene_id); candidates follow ucsc_gene_model.evaluate's dict shape.
    Raises LookupError(reason) when the gene or its coding transcripts cannot be found.
    """
    import ensembl_rest
    import http_retry
    import ncbi_rest
    from uniprot_api import resolve_accession

    reporter = resolve_reporter(report)
    query = cds_check.read_protein(protein_fa).rstrip("*")
    if not query:
        raise LookupError(f"no protein sequence in {protein_fa}")
    acc, how = resolve_accession(protein_fa, accession, taxid)
    if not acc:
        raise LookupError(how)
    gene_id = gene_id_for(acc)
    if not gene_id:
        raise LookupError(f"UniProt {acc} has no NCBI GeneID cross-reference")

    cfg = ensembl_rest.load_config()
    budget = deadline_s or ensembl_rest.timeouts()["total_timeout_s"]
    deadline = http_retry.deadline_from(budget)
    locs = ncbi_rest.gene_report(gene_id, cfg, deadline=deadline)
    if not locs:
        raise LookupError(f"NCBI GeneID {gene_id} has no annotated genomic location")
    ref = locs[0]
    # the gene plus flank: every CDS of the gene and its homology-arm sequence
    seq, lo, hi = ncbi_rest.fetch_region(ref["accession"], ref["start"], ref["end"], "+",
                                         flank_bp, cfg, deadline=deadline)
    xref = f"GeneID:{gene_id}"
    cds = [f for f in ncbi_rest.feature_table(ref["accession"], lo, hi, cfg, deadline=deadline)
           if f["type"] == "CDS" and xref in f["qualifiers"].get("db_xref", [])]
    if not cds:
        raise LookupError(f"NCBI GeneID {gene_id} has no annotated CDS on {ref['accession']}")
    products = ncbi_rest.product_report(gene_id, cfg, deadline=deadline)
    _report(reporter, f"NCBI GeneID {gene_id} ({ref.get('symbol') or acc}): {len(cds)} RefSeq "
                      f"CDS on {ref['accession']} ({ref['assembly']})", stage="ncbi_model")

    candidates = []
    for f in cds:
        protein = _protein_id(f["qualifiers"])
        product = products.get(protein, {})
        name = product.get("transcript") or protein
        try:
            codon_start = int((f["qualifiers"].get("codon_start") or ["1"])[0])
            model = cds_to_model(f["intervals"], codon_start, ref["accession"], (lo, seq),
                                 flank_bp)
        except (ValueError, IndexError) as e:
            _report(reporter, f"{name}: {e}", stage="ncbi_model", level="warning")
            continue
        check = cds_check.compare_to_query(
            cds_check.translate_model(model["cds_df"], model["dna"]), query)
        candidates.append({
            "source": "ncbi", "name": name,
            "label": f"NCBI RefSeq {name} ({protein}, GeneID {gene_id})",
            "model": model, "check": check, "select": bool(product.get("select_category")),
            "fasta_source": "ncbi", "assembly": ref["assembly"],
            "gff_label": f"ncbi_gene_model.py ({ref['assembly']} {name} GeneID {gene_id})",
            "region": {"source": "ncbi", "genome": ref["assembly"], "accession": acc,
                       "accession_from": how, "gene_id": gene_id, "transcript_id": name,
                       "protein_id": protein, "equivalent_transcripts": []},
        })
    candidates.sort(key=lambda c: (cds_check.rank(c["check"]), _preference(c)))
    # transcripts with the same CDS translation are equivalent to the chosen one
    exact = [c["name"] for c in candidates if c["check"]["status"] == "ok"]
    if exact:
        candidates[0]["region"]["equivalent_transcripts"] = exact
    # the closest few are enough to show what was compared
    for c in candidates[:cds_check.LOGGED_CANDIDATES]:
        _report(reporter, cds_check.summary_line(c["check"], c["label"]), stage="ncbi_model")
    if len(candidates) > cds_check.LOGGED_CANDIDATES:
        _report(reporter, f"... {len(candidates) - cds_check.LOGGED_CANDIDATES} further "
                          "candidate(s) differ more", stage="ncbi_model")
    return candidates, gene_id


def build(protein_fa, taxid, outprefix, accession="", flank_bp=2000, report=None,
          deadline_s=None):
    """Write the exact NCBI RefSeq model; raise LookupError(reason) when no transcript is exact."""
    from gene_model import write_model

    candidates, gene_id = evaluate(protein_fa, taxid, accession, flank_bp, report=report,
                                   deadline_s=deadline_s)
    if not candidates or candidates[0]["check"]["status"] != "ok":
        raise LookupError(f"no RefSeq transcript of GeneID {gene_id} translates to the input "
                          "protein" + (" — closest: " + cds_check.summary_line(
                              candidates[0]["check"], candidates[0]["label"])
                              if candidates else ""))
    out_txt, genomic_fa, region_json = write_model(outprefix, candidates[0])
    return {"source": "ncbi", "out_txt": out_txt, "genomic_fa": genomic_fa,
            "region_json": region_json, "transcript": candidates[0]["name"], "gene_id": gene_id}


def main(protein_fa, taxid, outprefix, accession="", flank_bp=2000, report=None):
    """CLI entry: build the model; return 0 on success, 2 when NCBI has none."""
    try:
        build(protein_fa, taxid, outprefix, accession, flank_bp, report=report)
    except LookupError as e:
        print(f"NCBI gene model unavailable: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Gene model + flanked region from the NCBI RefSeq "
                                        "transcript that translates exactly to the protein.")
    parser.add_argument("--protein_fasta", "--input_file", dest="protein_fasta", required=True,
                        help="protein FASTA (or PDB) of the run")
    parser.add_argument("--taxid", default="", help="species taxid (used for the CRC64 lookup)")
    parser.add_argument("--accession", "--uniprot_accession", dest="accession", default="",
                        help="UniProt accession (default: FASTA header, then CRC64 lookup)")
    parser.add_argument("--outprefix", required=True,
                        help="output prefix, e.g. {working_dir}/{run_name}_genewise")
    parser.add_argument("--flank_bp", type=int, default=2000,
                        help="genomic flank kept either side of the transcript (default 2000)")
    args, _ = parser.parse_known_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    sys.exit(main(args.protein_fasta, args.taxid, args.outprefix, args.accession, args.flank_bp))
