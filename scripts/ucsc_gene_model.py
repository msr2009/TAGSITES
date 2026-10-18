"""ucsc_gene_model.py — gene model + flanked genomic region from UCSC's UniProt alignments.

First choice for the reagents gene model (see gene_model.py): once the run's protein has a
UniProt accession, UCSC's unipAliSwissprot / unipAliTrembl tracks (UniProt sequences aligned
to the genome, on the 122 native assemblies) give that isoform's CDS blocks directly.
Their blocks are exactly the CDS, stop codon excluded, the same convention as Genewise.

  remote   /search by accession -> getData/track over the hit -> getData/sequence of the
           flanked region (3 REST calls, a few seconds; keyless)
  local    the bulk-mode accession index and genome (ucsc_local.py), no network; also finds
           an alignment by the protein's CRC64 when no accession is known

evaluate() translates EVERY alignment of the accession and its isoforms (UCSC sometimes lifts
an isoform through the wrong RefSeq transcript, i.e. the wrong alternative exon: Sh P08510-1)
and returns each as a candidate with its cds_check result; gene_model.py compares them with
the NCBI RefSeq candidates. build() (the CLI) writes a model only on an exact translation.
Writes the same files as genbank_input.genbank_to_genewise (so design_tag_reagents,
genbank_export and the Reagents tab are unchanged) plus {outprefix}.region.json with the
chromosome coordinates the phastCons track needs:

  {outprefix}.genewise.out.txt          GFF-embedded exon table (parse_genewise reads it)
  {outprefix}.genewise_genomic.fa       flanked region on the coding strand,
                                        header >ucsc:{genome}:{chrom}:{lo}:{hi}:{strand}
  {outprefix}.genewise_orientation.txt  "+"
  {outprefix}.region.json               genome, track, accession, chrom, span, flanks

Exit status 0 when a model was written, 2 when UCSC has none (the reason is printed).

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

def _ints(text):
    """'183,917,689,' -> [183, 917, 689]."""
    return [int(x) for x in str(text).strip(",").split(",") if x != ""]


def _rank(item, accession):
    """Sort key: the exact isoform, then the canonical (isMain) entry, Swiss-Prot first."""
    track = item.get("_track", "")
    exact = item.get("name") == accession
    main = bool(item.get("isMain")) and item.get("acc") == accession.split("-")[0]
    return (not exact, not main, track != "unipAliSwissprot", item.get("name", ""))


def item_to_model(item, flank_bp=2000, backend="remote", genome="", cfg=None, deadline=None,
                  region=None):
    """One UniProt alignment -> {cds_df, dna, chrom, start, stop, strand, flank_5p, flank_3p}.

    dna is the flanked region on the coding strand; cds_df indexes into it (0-based
    inclusive, GFF phase). start/stop are the alignment's own 1-based span. region is an
    optional pre-fetched (lo0, plus-strand sequence) covering the flanked span.
    """
    chrom, c0, c1, strand = item["chrom"], int(item["chromStart"]), int(item["chromEnd"]), \
        item.get("strand", "+")
    lo0 = max(0, c0 - flank_bp)
    hi = c1 + flank_bp
    if region is not None:
        rlo, rseq = region
        hi = min(hi, rlo + len(rseq))
        region = rseq[lo0 - rlo : hi - rlo]
    elif backend == "local":
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


def _fetch_regions(items, flank_bp, genome, deadline):
    """Remote: one plus-strand sequence per chromosome covering every item; {chrom: (lo0, seq)}."""
    import offtarget_blat
    import ucsc_rest

    spans = {}
    for it in items:
        lo, hi = spans.get(it["chrom"], (int(it["chromStart"]), int(it["chromEnd"])))
        spans[it["chrom"]] = (min(lo, int(it["chromStart"])), max(hi, int(it["chromEnd"])))
    out = {}
    for chrom, (lo, hi) in spans.items():
        size = next((int(it["chromSize"]) for it in items
                     if it["chrom"] == chrom and it.get("chromSize")), None)
        lo0, hi = max(0, lo - flank_bp), hi + flank_bp
        if size:
            hi = min(hi, size)
        ucsc_rest._pause(deadline)
        out[chrom] = (lo0, offtarget_blat.fetch_sequence(genome, chrom, lo0, hi,
                                                         deadline=deadline).upper())
    return out


def evaluate(protein_fa, taxid, accession="", flank_bp=2000, backend="remote", report=None,
             deadline_s=None):
    """Every UCSC UniProt alignment of the run's protein as a candidate model.

    Returns (candidates, genome); each candidate is {source, name, label, model, check,
    fasta_source, assembly, gff_label, region} (see gene_model.write_model), best first:
    exact translations, then the fewest differences, then _rank's isoform order.
    Raises LookupError(reason) when UCSC has no alignment.
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
    items = sorted(items, key=lambda it: _rank(it, acc or ""))
    regions = _fetch_regions(items, flank_bp, genome, deadline) if backend != "local" else {}

    candidates = []
    for order, item in enumerate(items):
        name, track = item.get("name", ""), item.get("_track", "")
        try:
            model = item_to_model(item, flank_bp, backend, genome, cfg, deadline,
                                  region=regions.get(item["chrom"]))
        except Exception as e:
            _report(reporter, f"{name}: {e}", stage="ucsc_model", level="warning")
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
        candidates.append({
            "source": "ucsc", "name": name, "label": f"UCSC {genome} {track} {name}",
            "model": model, "check": check, "order": order,
            "fasta_source": "ucsc", "assembly": genome,
            "gff_label": f"ucsc_gene_model.py ({genome} {track} {name})",
            "region": {"source": "ucsc", "genome": genome, "track": track, "accession": name,
                       "accession_from": how, "transcript_id": name,
                       "equivalent_transcripts": []},
        })
    candidates.sort(key=lambda c: (cds_check.rank(c["check"]), c["order"]))
    # the closest few are enough to show what was compared
    for c in candidates[:cds_check.LOGGED_CANDIDATES]:
        _report(reporter, cds_check.summary_line(c["check"], c["label"]), stage="ucsc_model")
    if len(candidates) > cds_check.LOGGED_CANDIDATES:
        _report(reporter, f"... {len(candidates) - cds_check.LOGGED_CANDIDATES} further "
                          "candidate(s) differ more", stage="ucsc_model")
    return candidates, genome


def build(protein_fa, taxid, outprefix, accession="", flank_bp=2000, backend="remote",
          report=None, deadline_s=None):
    """Write the UCSC gene model for the run's protein; raise LookupError(reason) if none is exact.

    Returns {source, out_txt, genomic_fa, region_json, genome, accession, track}.
    """
    from gene_model import write_model

    reporter = resolve_reporter(report)
    candidates, genome = evaluate(protein_fa, taxid, accession, flank_bp, backend,
                                  report=report, deadline_s=deadline_s)
    exact = [c for c in candidates if c["check"]["status"] == "ok"]
    if not exact:
        best = candidates[0]["check"] if candidates else None
        raise LookupError(f"UCSC {genome}: no alignment translates to the input protein"
                          + (" — closest: " + cds_check.summary_line(best, candidates[0]["label"])
                             if best else ""))
    c = exact[0]
    out_txt, genomic_fa, region_json = write_model(outprefix, c)
    m = c["model"]
    _report(reporter, f"Gene model from {c['label']} ({m['chrom']}:{m['start']}-{m['stop']}:"
                      f"{m['strand']}, {len(m['cds_df'])} exons; accession from "
                      f"{c['region']['accession_from']})", stage="ucsc_model")
    return {"source": "ucsc", "out_txt": out_txt, "genomic_fa": genomic_fa,
            "region_json": region_json, "genome": genome, "accession": c["name"],
            "track": c["region"]["track"]}


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
