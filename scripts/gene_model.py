"""gene_model.py — the one gene-model step for reagent design, shared by the app and the CLI.

Picks a single gene model for the run's protein, in this order, and writes it as the
Genewise-format files ({outprefix}.genewise.out.txt / .genewise_genomic.fa /
.genewise_orientation.txt) that design_tag_reagents reads:

  1. a user GenBank with CDS/exon features         genbank_input.genbank_to_genewise
  2. UCSC's UniProt alignment for the accession    ucsc_gene_model.build
                                                   (REST; the local index in bulk mode)
  3. fallback: EBI Genewise on the genomic region  run_genewise.main (remote mode)
            or the local GFF3 exact match           genewise_bulk.main  (bulk mode)

Only one source wins; models are never merged. Every source is held to cds_check (exact
translation for UCSC; GenBank and Genewise keep their warn/fail limits). A UCSC model needs
no genomic upload; an uploaded plain region is used only when UCSC has no model.

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


def ensure_model(protein_fa, genomic_fa, taxid, outprefix, email="", accession="",
                 flank_bp=2000, report=None, job_id_cb=None, resume_job_ids=None):
    """Write the run's gene model (GenBank > UCSC > Genewise/GFF3); raise ValueError if none."""
    import genbank_input
    import providers
    import ucsc_gene_model

    reporter = resolve_reporter(report)
    bulk = providers.region_from_annotation()

    # 1. user-curated exons always win (raises if they do not match the protein)
    if genomic_fa and genbank_input.is_genbank(genomic_fa) \
            and genbank_input.has_gene_model(genomic_fa):
        _clear(outprefix)
        genbank_input.genbank_to_genewise(genomic_fa, protein_fa, outprefix, report=report)
        return _record(outprefix, _files(outprefix, "genbank",
                                         note=f"GenBank {os.path.basename(genomic_fa)}"))

    # 2. UCSC's UniProt alignment, exact translation only
    try:
        res = ucsc_gene_model.build(protein_fa, taxid, outprefix, accession=accession,
                                    flank_bp=flank_bp, backend="local" if bulk else "remote",
                                    report=report)
        note = f"UCSC {res['genome']} (UniProt {res['accession']})"
        if genomic_fa:
            note += "; the uploaded genomic region was not needed"
        return _record(outprefix, _files(outprefix, "ucsc", note, res["region_json"]))
    except Exception as e:
        reason = str(e)
        _report(reporter, f"no UCSC gene model ({reason}); falling back", stage="gene_model",
                level="warning")

    # 3. fallback: local GFF3 in bulk mode, EBI Genewise on the supplied region otherwise
    _clear(outprefix)
    if bulk:
        import genewise_bulk

        genewise_bulk.main(protein_fa, genomic_fa, email, outprefix, report=report)
        return _record(outprefix, _files(outprefix, "gff3", "local GFF3 transcript",
                                         outprefix + ".region.json"))
    if not genomic_fa:
        raise ValueError(f"No gene model: {reason}. Fetch or upload the genomic region so "
                         "Genewise can build one.")
    import run_genewise

    if genbank_input.is_genbank(genomic_fa):
        # sequence-only GenBank: Genewise needs a FASTA
        genomic_fa = genbank_input.genbank_to_fasta(genomic_fa, outprefix + ".genomic.fa")
    result = run_genewise.main(protein_fa, genomic_fa, email, outprefix, report=report,
                               job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)
    # only genewise's own {"ebi_status": ...} sentinel means "not done yet"
    if isinstance(result, dict) and "ebi_status" in result:
        return result
    return _record(outprefix, _files(outprefix, "genewise",
                                     f"EBI Genewise (UCSC: {reason})"))


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Write the run's gene model: GenBank exons > UCSC "
                                        "UniProt alignment > Genewise (or local GFF3).")
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
