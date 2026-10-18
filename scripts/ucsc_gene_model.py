"""ucsc_gene_model.py — gene model + flanked genomic region from UCSC's UniProt alignments.

First choice for the reagents gene model (see gene_model.py): once the run's protein has a
UniProt accession, UCSC's unipAliSwissprot / unipAliTrembl tracks (UniProt sequences aligned
to the genome, on the 122 native assemblies) give that isoform's CDS blocks directly.
Their blocks are exactly the CDS, stop codon excluded, the same convention as Genewise.

  remote   /search by accession -> getData/track over the hit -> getData/sequence of the
           flanked region (3 REST calls, a few seconds; keyless)
  local    the bulk-mode accession index and genome (ucsc_local.py), no network; also finds
           an alignment by the protein's CRC64 when no accession is known

The model is accepted only if it translates EXACTLY to the input protein; anything else is
a failure and the caller falls back (Genewise remotely, the GFF3 match in bulk mode). Writes
the same files as genbank_input.genbank_to_genewise (so design_tag_reagents, genbank_export
and the Reagents tab are unchanged) plus {outprefix}.region.json with the chromosome
coordinates the phastCons track needs:

  {outprefix}.genewise.out.txt          GFF-embedded exon table (parse_genewise reads it)
  {outprefix}.genewise_genomic.fa       flanked region on the coding strand,
                                        header >ucsc:{genome}:{chrom}:{lo}:{hi}:{strand}
  {outprefix}.genewise_orientation.txt  "+"
  {outprefix}.region.json               genome, track, accession, chrom, span, flanks

Exit status 0 when a model was written, 2 when UCSC has none (the reason is printed).

Matt Rich, 2026
"""

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import cds_check
from crispr_util import reverse_complement
from parse_genewise import spans_from_positions, write_genewise_gff
from progress import report as _report, resolve_reporter

MAX_CANDIDATES = 4  # alignments tried per accession (each costs one sequence fetch remotely)


def _ints(text):
    """'183,917,689,' -> [183, 917, 689]."""
    return [int(x) for x in str(text).strip(",").split(",") if x != ""]


def _rank(item, accession):
    """Sort key: the exact isoform, then the canonical (isMain) entry, Swiss-Prot first."""
    track = item.get("_track", "")
    exact = item.get("name") == accession
    main = bool(item.get("isMain")) and item.get("acc") == accession.split("-")[0]
    return (not exact, not main, track != "unipAliSwissprot", item.get("name", ""))


def item_to_model(item, flank_bp=2000, backend="remote", genome="", cfg=None, deadline=None):
    """One UniProt alignment -> {cds_df, dna, chrom, start, stop, strand, flank_5p, flank_3p}.

    dna is the flanked region on the coding strand; cds_df indexes into it (0-based
    inclusive, GFF phase). start/stop are the alignment's own 1-based span.
    """
    chrom, c0, c1, strand = item["chrom"], int(item["chromStart"]), int(item["chromEnd"]), \
        item.get("strand", "+")
    lo0 = max(0, c0 - flank_bp)
    hi = c1 + flank_bp
    if backend == "local":
        import genome_regions
        import ucsc_local

        lchrom = ucsc_local.local_chrom(chrom)
        hi = min(hi, len(genome_regions._load_genome()[lchrom]))
        region = genome_regions.extract_sequence(lchrom, lo0 + 1, hi, "+").upper()
    else:
        import offtarget_blat

        size = item.get("chromSize")
        if size:
            hi = min(hi, int(size))
        region = offtarget_blat.fetch_sequence(genome, chrom, lo0, hi, deadline=deadline).upper()
    if len(region) != hi - lo0:
        raise RuntimeError(f"got {len(region)} bp for {chrom}:{lo0}-{hi}, expected {hi - lo0}")

    # CDS bases in chromosome order, then in transcript order on the coding strand
    pos = [c0 + s + i for s, n in zip(_ints(item["chromStarts"]), _ints(item["blockSizes"]))
           for i in range(n)]
    pos = [p - lo0 for p in pos]
    n = len(region)
    if strand == "-":
        region = reverse_complement(region)
        pos = sorted(n - 1 - p for p in pos)
    low, high = c0 - lo0, hi - c1  # genomic-low / genomic-high flank actually obtained
    return {
        "cds_df": spans_from_positions(pos), "dna": region,
        "chrom": chrom, "start": c0 + 1, "stop": c1, "strand": strand,
        "flank_5p": low if strand == "+" else high,
        "flank_3p": high if strand == "+" else low,
    }


def extend_termini(model, check):
    """Model with missing terminal codons restored, or None if the gaps are not terminal.

    UniProt alignments sometimes stop one or two codons short of the protein's ends (e.g.
    the initial Met). When the only differences are such terminal deletions, the adjacent
    bases on the coding strand are taken as those codons; the caller re-checks the result
    for an exact translation, so a wrong guess is never accepted.
    """
    if check["mismatches"] or check["insertions"] or not check["deletions"]:
        return None
    qlen = check["query_len"]
    head = tail = 0
    for d in check["deletions"]:
        # only deletions touching residue 1 or the last residue are terminal
        if d["start"] == 1:
            head = d["stop"]
        elif d["stop"] == qlen:
            tail = d["stop"] - d["start"] + 1
        else:
            return None
    bases = [p for r in model["cds_df"].itertuples() for p in range(r.start, r.stop + 1)]
    first, last = bases[0] - 3 * head, bases[-1] + 3 * tail
    if first < 0 or last >= len(model["dna"]):
        return None
    bases = list(range(first, bases[0])) + bases + list(range(bases[-1] + 1, last + 1))
    return dict(model, cds_df=spans_from_positions(bases))


def _remote_items(genome, accession, cfg, deadline, reporter):
    """REST: alignment items for an accession (and its isoforms) on an assembly."""
    import ucsc_rest

    hits = ucsc_rest.search_accession(genome, accession, cfg, deadline=deadline)
    if not hits:
        return []
    names = {h[4] for h in hits}
    items = []
    # one getData/track call per (track, chrom) span covering that track's hits there
    spans = {}
    for track, chrom, start, end, _ in hits:
        lo, hi = spans.get((track, chrom), (start, end))
        spans[(track, chrom)] = (min(lo, start), max(hi, end))
    for (track, chrom), (start, end) in spans.items():
        ucsc_rest._pause(deadline)
        for it in ucsc_rest.track_items(genome, track, chrom, start - 1, end, cfg,
                                        deadline=deadline):
            if it.get("name") in names:
                items.append(dict(it, _track=track))
    _report(reporter, f"UCSC {genome}: {len(items)} UniProt alignment(s) for {accession}",
            stage="ucsc_model")
    return items


def _write(outprefix, model, genome, track, accession, how):
    """Write the Genewise-format files plus region.json; return their paths."""
    from fetch_genomic_sequence import _format_fasta

    out_txt = f"{outprefix}.genewise.out.txt"
    genomic_fa = f"{outprefix}.genewise_genomic.fa"
    region_json = f"{outprefix}.region.json"
    # flanked region bounds, 1-based inclusive, for the FASTA header
    lo = model["start"] - (model["flank_5p"] if model["strand"] == "+" else model["flank_3p"])
    hi = model["stop"] + (model["flank_3p"] if model["strand"] == "+" else model["flank_5p"])
    with open(genomic_fa, "w") as fh:
        fh.write(_format_fasta("ucsc", genome, model["chrom"], lo, hi, model["strand"],
                               model["dna"]))
    write_genewise_gff(out_txt, model["cds_df"], model["chrom"],
                       f"ucsc_gene_model.py ({genome} {track} {accession})")
    with open(f"{outprefix}.genewise_orientation.txt", "w") as fh:
        fh.write("+\n")  # the written sequence is already on the coding strand
    with open(region_json, "w") as fh:
        json.dump({"source": "ucsc", "genome": genome, "track": track,
                   "accession": accession, "accession_from": how,
                   "transcript_id": accession, "chrom": model["chrom"],
                   "start": model["start"], "stop": model["stop"], "strand": model["strand"],
                   "flank_5p": model["flank_5p"], "flank_3p": model["flank_3p"],
                   "equivalent_transcripts": []}, fh)
    return out_txt, genomic_fa, region_json


def build(protein_fa, taxid, outprefix, accession="", flank_bp=2000, backend="remote",
          report=None, deadline_s=None):
    """Write the UCSC gene model for the run's protein; raise LookupError(reason) if none.

    Returns {source, out_txt, genomic_fa, region_json, genome, accession, track}.
    """
    import http_retry
    from Bio.SeqUtils.CheckSum import crc64
    from uniprot_api import resolve_accession

    reporter = resolve_reporter(report)
    query = cds_check.read_protein(protein_fa).rstrip("*")
    if not query:
        raise LookupError(f"no protein sequence in {protein_fa}")

    if backend == "local":
        import ucsc_local

        genome = ucsc_local.settings()[0]
        if not ucsc_local.available():
            raise LookupError("local UCSC index not built (python scripts/reference_data.py "
                              "--only ucsc)")
        # bulk input usually carries its accession; the CRC64 lookup is local, not UniProt
        acc, how = resolve_accession(protein_fa, accession, taxid, use_checksum=False)
        items = ucsc_local.lookup(acc, crc64(query).replace("CRC-", ""))
        cfg = deadline = None
    else:
        import ensembl_rest
        import offtarget_blat
        import offtarget_screen

        genome = offtarget_blat.assembly_for_taxid(taxid, offtarget_screen.load_config()) \
            if taxid else ""
        if not genome:
            raise LookupError(f"no UCSC assembly configured for taxid {taxid or '(none)'}")
        acc, how = resolve_accession(protein_fa, accession, taxid)
        if not acc:
            raise LookupError(how)
        cfg = ensembl_rest.load_config()
        budget = deadline_s or ensembl_rest.timeouts()["total_timeout_s"]
        deadline = http_retry.deadline_from(budget)
        items = _remote_items(genome, acc, cfg, deadline, reporter)

    label = acc or "the input sequence"
    if not items:
        raise LookupError(f"UCSC {genome} has no UniProt alignment for {label}")

    reasons = []
    for item in sorted(items, key=lambda it: _rank(it, acc or ""))[:MAX_CANDIDATES]:
        name, track = item.get("name", ""), item.get("_track", "")
        try:
            if backend != "local":
                import ucsc_rest

                ucsc_rest._pause(deadline)
            model = item_to_model(item, flank_bp, backend, genome, cfg, deadline)
        except Exception as e:
            reasons.append(f"{name}: {e}")
            continue
        check = cds_check.compare_to_query(
            cds_check.translate_model(model["cds_df"], model["dna"]), query)
        extended = extend_termini(model, check) if check["status"] != "ok" else None
        if extended is not None:
            echeck = cds_check.compare_to_query(
                cds_check.translate_model(extended["cds_df"], extended["dna"]), query)
            if echeck["status"] == "ok":
                _report(reporter, f"{name}: alignment lacked terminal codon(s); restored from "
                                  "the adjacent genomic bases (exact translation)",
                        stage="ucsc_model")
                model, check = extended, echeck
        if check["status"] != "ok":
            reasons.append(f"{name} ({track}) does not translate to the input protein "
                           f"({check['model_len']} vs {check['query_len']} aa, "
                           f"{len(check['mismatches'])} substitution(s))")
            continue
        out_txt, genomic_fa, region_json = _write(outprefix, model, genome, track, name, how)
        _report(reporter, f"Gene model from UCSC {genome} {track} {name} "
                          f"({model['chrom']}:{model['start']}-{model['stop']}:{model['strand']}, "
                          f"{len(model['cds_df'])} exons; accession from {how})",
                stage="ucsc_model")
        return {"source": "ucsc", "out_txt": out_txt, "genomic_fa": genomic_fa,
                "region_json": region_json, "genome": genome, "accession": name,
                "track": track}
    raise LookupError(f"UCSC {genome}: no alignment of {label} translates to the input protein"
                      + (" — " + "; ".join(reasons) if reasons else ""))


def main(protein_fa, taxid, outprefix, accession="", flank_bp=2000, backend="remote",
         report=None):
    """CLI entry: build the model; return 0 on success, 2 when UCSC has none."""
    try:
        build(protein_fa, taxid, outprefix, accession, flank_bp, backend, report=report)
    except LookupError as e:
        print(f"UCSC gene model unavailable: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Gene model + flanked region from UCSC's UniProt "
                                        "alignments (exact translation required).")
    parser.add_argument("--protein_fasta", "--input_file", dest="protein_fasta", required=True,
                        help="protein FASTA (or PDB) of the run")
    parser.add_argument("--taxid", default="", help="species taxid (selects the UCSC assembly)")
    parser.add_argument("--accession", "--uniprot_accession", dest="accession", default="",
                        help="UniProt accession (default: FASTA header, then CRC64 lookup)")
    parser.add_argument("--outprefix", required=True,
                        help="output prefix, e.g. {working_dir}/{run_name}_genewise")
    parser.add_argument("--flank_bp", type=int, default=2000,
                        help="genomic flank kept either side of the CDS (default 2000)")
    parser.add_argument("--backend", choices=("remote", "local"), default="",
                        help="remote REST or the local bulk index (default: from "
                             "batch.config.json backends.genewise)")
    args, _ = parser.parse_known_args()
    backend = args.backend
    if not backend:
        import providers

        backend = "local" if providers.region_from_annotation() else "remote"
    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    sys.exit(main(args.protein_fasta, args.taxid, args.outprefix, args.accession,
                  args.flank_bp, backend))
