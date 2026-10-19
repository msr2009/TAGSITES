"""gene_model.py — the one gene-model step for reagent design, shared by the app and the CLI.

Picks a single gene model for the run's protein and writes it as the Genewise-format files
({outprefix}.genewise.out.txt / .genewise_genomic.fa / .genewise_orientation.txt) that
design_tag_reagents reads. Sources are tried in this order:

  1. a user GenBank with CDS/exon features         genbank_input.genbank_to_genewise
  2. UCSC's UniProt alignments (every isoform)      ucsc_gene_model.evaluate
                                                   (REST; the local index in bulk mode)
  3. every NCBI RefSeq transcript of the gene       ncbi_gene_model.evaluate (remote only)
     or every GFF3 transcript of the gene           genewise_bulk.evaluate   (bulk mode)
  4. the closest candidate of 2-3 whose differences cds_check only warns about (isolated
     substitutions, short gaps: genome variation, not a wrong exon)
  5. fallback: EBI Genewise on the genomic region  run_genewise.main (remote mode)

An exact translation wins at steps 2-3 (searched in parallel remotely; select_model, also
previewed by the app at protein input). Only one source wins; models are never merged. Every
source is held to the same cds_check rules (clustered substitutions, long gaps or gaps beside
substitutions fail). A UCSC or NCBI model needs no genomic upload; an uploaded plain region is
used only by Genewise.

ensure_model() returns {"source", "out_txt", "genomic_fa", "region_json", "note"}, or
Genewise's {"ebi_status": ...} sentinel while its EBI jobs are still queued (app path).
The CLI (used by run_tag_sites_from_json.py) exits 0 when a model was written, 1 otherwise.

Matt Rich, 2026
"""

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from progress import report as _report, resolve_reporter

SOURCE_FILE = "{outprefix}.model_source.json"  # which source won, read by the Reagents tab


def _clear(outprefix):
    """Remove a previous run's model files, so a stale region.json never outlives its model."""
    for suffix in (".genewise.out.txt", ".genewise_genomic.fa", ".genewise_orientation.txt",
                   ".region.json", ".model_source.json"):
        try:
            os.remove(outprefix + suffix)
        except FileNotFoundError:
            pass


def _record(outprefix, result):
    """Write {outprefix}.model_source.json and return result unchanged."""
    with open(SOURCE_FILE.format(outprefix=outprefix), "w") as fh:
        json.dump({k: v for k, v in result.items() if k != "ebi_status"}, fh)
    return result


def _files(outprefix, source, note="", region_json=None):
    """Result dict for a model written under outprefix."""
    return {"source": source, "out_txt": outprefix + ".genewise.out.txt",
            "genomic_fa": outprefix + ".genewise_genomic.fa",
            "region_json": region_json, "note": note}


def read_source(outprefix):
    """The recorded model source for a run, or None (models predating this file)."""
    try:
        with open(SOURCE_FILE.format(outprefix=outprefix)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def model_fasta(cand):
    """A candidate's flanked genomic window as FASTA text, coding strand (the app's region)."""
    from fetch_genomic_sequence import _format_fasta

    model = cand["model"]
    if not cand.get("fasta_source"):
        return f">{cand['fasta_id']}\n{model['dna']}\n"
    plus = model["strand"] == "+"
    # flanked region bounds, 1-based inclusive, for the FASTA header
    lo = model["start"] - (model["flank_5p"] if plus else model["flank_3p"])
    hi = model["stop"] + (model["flank_3p"] if plus else model["flank_5p"])
    return _format_fasta(cand["fasta_source"], cand["assembly"], model["chrom"], lo, hi,
                         model["strand"], model["dna"])


def write_model(outprefix, cand):
    """Write a candidate model (see ucsc_gene_model.evaluate) as the Genewise-format files.

    Returns (out_txt, genomic_fa, region_json). The FASTA header is the auto-fetch form
    >source:assembly:chrom:lo:hi:strand (fasta_source set) or >{fasta_id} (GFF3); region.json
    carries the candidate's region keys plus the chromosome coordinates phastCons needs.
    """
    from parse_genewise import write_genewise_gff

    model = cand["model"]
    out_txt = f"{outprefix}.genewise.out.txt"
    genomic_fa = f"{outprefix}.genewise_genomic.fa"
    region_json = f"{outprefix}.region.json"
    with open(genomic_fa, "w") as fh:
        fh.write(model_fasta(cand))
    write_genewise_gff(out_txt, model["cds_df"], model["chrom"], cand["gff_label"])
    with open(f"{outprefix}.genewise_orientation.txt", "w") as fh:
        fh.write("+\n")  # the written sequence is already on the coding strand
    with open(region_json, "w") as fh:
        json.dump(dict(cand["region"], chrom=model["chrom"], start=model["start"],
                       stop=model["stop"], strand=model["strand"], flank_5p=model["flank_5p"],
                       flank_3p=model["flank_3p"]), fh)
    return out_txt, genomic_fa, region_json


SOURCE_ORDER = {"ucsc": 0, "ncbi": 1, "gff3": 2}  # tie-break between equally close candidates


def _closest(candidates):
    """The best candidate overall (exact, then fewest differences, then source order), or None."""
    import cds_check

    ranked = sorted(candidates, key=lambda c: (cds_check.rank(c["check"]),
                                               SOURCE_ORDER[c["source"]]))
    return ranked[0] if ranked else None


def _loci(candidates):
    """Distinct genomic loci (chrom, start, stop) of candidates, overlapping spans merged."""
    loci = []
    for c in sorted(candidates, key=lambda c: (c["model"]["chrom"], c["model"]["start"])):
        m = c["model"]
        if loci and loci[-1][0] == m["chrom"] and m["start"] <= loci[-1][2]:
            loci[-1][2] = max(loci[-1][2], m["stop"])
        else:
            loci.append([m["chrom"], m["start"], m["stop"]])
    return loci


def _use(outprefix, cand, note, reporter):
    """Write a winning candidate, log it, and record it as the run's model."""
    import cds_check

    _clear(outprefix)
    _, _, region_json = write_model(outprefix, cand)
    level = "warning" if cand["check"]["status"] == "warn" else "info"
    _report(reporter, cds_check.describe(cand["check"], cand["label"]), stage="gene_model",
            level=level)
    return _record(outprefix, _files(outprefix, cand["source"], note, region_json))


def _ucsc(protein_fa, taxid, accession, flank_bp, bulk, report, deadline_s):
    """UCSC candidates as (candidates, genome, reason or None); never raises."""
    import ucsc_gene_model

    try:
        ucsc, genome = ucsc_gene_model.evaluate(protein_fa, taxid, accession=accession,
                                                flank_bp=flank_bp,
                                                backend="local" if bulk else "remote",
                                                report=report, deadline_s=deadline_s)
    except Exception as e:
        return [], "", f"UCSC: {e}"
    exact = ucsc and ucsc[0]["check"]["status"] == "ok"
    return ucsc, genome, None if exact else f"UCSC {genome}: no alignment translates exactly"


def _ncbi(protein_fa, taxid, accession, flank_bp, report, deadline_s):
    """NCBI RefSeq candidates as (candidates, gene_id, reason or None); never raises."""
    import ncbi_gene_model

    try:
        ncbi, gene_id = ncbi_gene_model.evaluate(protein_fa, taxid, accession=accession,
                                                 flank_bp=flank_bp, report=report,
                                                 deadline_s=deadline_s)
    except Exception as e:
        return [], "", f"NCBI: {e}"
    exact = ncbi and ncbi[0]["check"]["status"] == "ok"
    return ncbi, gene_id, None if exact else (f"NCBI GeneID {gene_id}: no RefSeq transcript "
                                              "translates exactly")


def select_model(protein_fa, taxid, accession="", flank_bp=2000, bulk=False, report=None,
                 deadline_s=None):
    """Steps 2-4 of the module docstring: pick the best UCSC/NCBI (or GFF3) candidate.

    Returns {"winner", "note", "best", "reasons", "gene_id", "genome"}: winner is the chosen
    candidate (exact, or the closest that only warns) or None; best is the closest candidate
    even when it fails. Remotely, UCSC and NCBI are searched in parallel (~5 s instead of the
    sum, ~10 s for Sh) and the order is applied to the combined results. The app previews
    this at protein input; the reagents task repeats it.
    """
    import cds_check

    reporter = resolve_reporter(report)
    out = {"winner": None, "note": "", "best": None, "reasons": [], "gene_id": "", "genome": ""}

    # 2-3. remote: UCSC and NCBI concurrently. Bulk: UCSC-local first, since its alignments
    # locate the gene for the GFF3 search
    if bulk:
        ucsc, out["genome"], ucsc_reason = _ucsc(protein_fa, taxid, accession, flank_bp, True,
                                                 report, deadline_s)
        other, other_reason = _gff3(protein_fa, ucsc, report)
    else:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=2) as pool:
            ucsc_job = pool.submit(_ucsc, protein_fa, taxid, accession, flank_bp, False, report,
                                   deadline_s)
            ncbi_job = pool.submit(_ncbi, protein_fa, taxid, accession, flank_bp, report,
                                   deadline_s)
            ucsc, out["genome"], ucsc_reason = ucsc_job.result()
            other, out["gene_id"], other_reason = ncbi_job.result()
    candidates = ucsc + other
    out["reasons"] = [r for r in (ucsc_reason, other_reason) if r]
    # errors (not mere mismatches) are worth a log line
    for r in out["reasons"]:
        if r.startswith(("UCSC:", "NCBI:", "GFF3:")):
            _report(reporter, f"no {r.split(':')[0]} gene model ({r.split(': ', 1)[1]})",
                    stage="gene_model", level="warning")

    loci = _loci(ucsc)
    if len(loci) > 1:
        # UniProt also aligns to a paralog: only an exact UCSC model is trusted, and the
        # closest non-exact one must come from the gene-locked NCBI/GFF3 candidates
        _report(reporter, f"UCSC aligns {ucsc[0]['name'].split('-')[0]} to {len(loci)} loci ("
                          + ", ".join(f"{c}:{a}-{b}" for c, a, b in loci) + "); a "
                          "non-exact UCSC model will not be used", stage="gene_model",
                level="warning")

    # an exact translation wins, UCSC before NCBI/GFF3
    if ucsc and ucsc[0]["check"]["status"] == "ok":
        out["winner"] = out["best"] = ucsc[0]
        out["note"] = f"UCSC {out['genome']} (UniProt {ucsc[0]['name']})"
        return out
    if other and other[0]["check"]["status"] == "ok":
        out["winner"] = out["best"] = other[0]
        out["note"] = (f"local GFF3 transcript {other[0]['name']}" if bulk else
                       f"NCBI RefSeq {other[0]['name']} (GeneID {out['gene_id']})")
        return out

    # 4. the closest candidate, if its differences look like genome variation; UCSC alignments
    # at several loci are left out, since the closest of them may be a paralog
    trusted = [c for c in candidates if c["source"] != "ucsc" or len(loci) == 1]
    best = _closest(trusted) or _closest(candidates)
    out["best"] = best
    if best is not None and best["check"]["status"] == "warn" and best in trusted:
        out["winner"] = best
        out["note"] = (f"closest model: {best['label']} "
                       f"({cds_check.summary_line(best['check'])})")
    return out


def _gff3(protein_fa, ucsc, report):
    """Bulk mode's GFF3 candidates near the UCSC locus, as (candidates, reason or None)."""
    import genewise_bulk

    # the UCSC alignments locate the gene for the closest-transcript search, unless they sit
    # at more than one locus (a paralog)
    best = _closest(ucsc) if len(_loci(ucsc)) == 1 else None
    locus = (best["model"]["chrom"], best["model"]["start"], best["model"]["stop"],
             best["model"]["strand"]) if best else None
    try:
        gff = genewise_bulk.evaluate(protein_fa, locus=locus, report=report)
    except Exception as e:
        return [], f"GFF3: {e}"
    exact = gff and gff[0]["check"]["status"] == "ok"
    return gff, None if exact else "GFF3: no transcript translates exactly"


def ensure_model(protein_fa, genomic_fa, taxid, outprefix, email="", accession="",
                 flank_bp=2000, report=None, job_id_cb=None, resume_job_ids=None):
    """Write the run's gene model (see the module docstring); raise ValueError if none qualifies."""
    import cds_check
    import genbank_input
    import providers

    reporter = resolve_reporter(report)
    bulk = providers.region_from_annotation()

    # 1. user-curated exons always win (raises if they do not match the protein)
    if genomic_fa and genbank_input.is_genbank(genomic_fa) \
            and genbank_input.has_gene_model(genomic_fa):
        _clear(outprefix)
        genbank_input.genbank_to_genewise(genomic_fa, protein_fa, outprefix, report=report)
        return _record(outprefix, _files(outprefix, "genbank",
                                         note=f"GenBank {os.path.basename(genomic_fa)}"))

    # 2-4. UCSC, NCBI/GFF3, closest warn
    sel = select_model(protein_fa, taxid, accession=accession, flank_bp=flank_bp, bulk=bulk,
                       report=report)
    if sel["winner"] is not None:
        uploaded = "; the uploaded genomic region was not needed" if genomic_fa else ""
        return _use(outprefix, sel["winner"], sel["note"] + uploaded, reporter)
    best = sel["best"]
    closest = ("\nClosest candidate:\n" + cds_check.describe(best["check"], best["label"])) \
        if best else ""
    tried = "; ".join(sel["reasons"])

    # 5. fallback: Genewise on the supplied region (remote mode only)
    _clear(outprefix)
    if bulk:
        raise ValueError(f"No gene model qualifies ({tried}).{closest}\n"
                         + cds_check.MODEL_MISMATCH_HELP)
    if not genomic_fa:
        raise ValueError(f"No gene model qualifies ({tried}).{closest}\nFetch or upload the "
                         "genomic region so Genewise can build one. "
                         + cds_check.MODEL_MISMATCH_HELP)
    if best is not None:
        _report(reporter, "no UCSC or NCBI model qualifies; closest:\n"
                          + cds_check.describe(best["check"], best["label"]),
                stage="gene_model", level="warning")
    import run_genewise

    if genbank_input.is_genbank(genomic_fa):
        # sequence-only GenBank: Genewise needs a FASTA
        genomic_fa = genbank_input.genbank_to_fasta(genomic_fa, outprefix + ".genomic.fa")
    result = run_genewise.main(protein_fa, genomic_fa, email, outprefix, report=report,
                               job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)
    # only genewise's own {"ebi_status": ...} sentinel means "not done yet"
    if isinstance(result, dict) and "ebi_status" in result:
        return result
    return _record(outprefix, _files(outprefix, "genewise", f"EBI Genewise ({tried})"))


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Write the run's gene model: GenBank exons > UCSC "
                                        "UniProt alignment > NCBI RefSeq (or local GFF3) "
                                        "transcripts > Genewise.")
    parser.add_argument("--protein_fasta", "--input_file", dest="protein_fasta", required=True)
    parser.add_argument("--genomic_fasta", default="",
                        help="genomic region (FASTA or GenBank); needed only for the fallback")
    parser.add_argument("--taxid", default="")
    parser.add_argument("--accession", "--uniprot_accession", dest="accession", default="")
    parser.add_argument("--email", default="", help="required by EBI for the Genewise fallback")
    parser.add_argument("--outprefix", required=True,
                        help="output prefix, e.g. {working_dir}/{run_name}_genewise")
    parser.add_argument("--flank_bp", type=int, default=2000)
    args, _ = parser.parse_known_args()
    try:
        out = ensure_model(args.protein_fasta, args.genomic_fasta, args.taxid, args.outprefix,
                           email=args.email, accession=args.accession, flank_bp=args.flank_bp)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    print(json.dumps(out))
    sys.exit(0 if "ebi_status" not in out else 1)
