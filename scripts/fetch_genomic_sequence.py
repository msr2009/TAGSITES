"""
fetch_genomic_sequence.py

Automate genomic-sequence extraction for a gene, given an organism taxid and a
gene symbol. Replaces manual genomic FASTA upload/paste.

Sources are tried in order (ensembl.config.json -> sources.order, default
ncbi, ucsc, ensembl) and the first that returns a sequence wins:
  ncbi     NCBI Datasets + efetch (ncbi_rest.py), any annotated NCBI genome;
           measured at 0.3-0.8 s per lookup, so it goes first
  ucsc     UCSC REST (ucsc_rest.py) for the assemblies in offtarget.config.json;
           4-8 s, since UCSC asks for ~1 request/s across its three calls
  ensembl  Ensembl REST (ensembl_rest.py), the original and now last resort,
           since it is often slow or down
One wall-clock deadline covers the whole chain; every source but the last is
also capped at sources.per_source_budget_s so a stalled first source cannot
starve the rest. A source that times out, errors, or does not know the
organism/gene hands over to the next.

Ensembl pipeline: symbol -> gene id (xref_symbol) -> coordinates (lookup_id) ->
flanked genomic FASTA (fetch_region_fasta). flank_bp of extra sequence is
fetched on each side of the gene body so downstream reagent design
(scripts/design_tag_reagents.py) has enough flanking sequence for homology
arms; the default (2000) covers the Design Reagents tab's max allowed
arm_length (2000bp).

Matt Rich, 2026
"""

import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent))
import ensembl_rest
import http_retry
import ncbi_rest
import offtarget_blat
import offtarget_screen
import ucsc_rest
from progress import report as _report, resolve_reporter

SOURCES = ("ncbi", "ucsc", "ensembl")


def _pick_most_complete_gene(candidates, deadline=None):
    """When a symbol resolves to >1 gene id, prefer the one with the largest genomic span.

    Genewise needs the FULL gene body to align the protein against; a
    truncated/partial duplicate annotation (e.g. on a patch/alt-haplotype
    contig with an incomplete gene model) has a shorter genomic span than the
    complete one — confirmed directly for zebrafish tp53 (11,579 bp primary
    vs. 8,702 bp duplicate) and gata1a (6,943 bp vs. 5,054 bp). This targets
    sequence completeness directly instead of guessing from Ensembl's
    assembly/region naming conventions, which vary by species and aren't a
    documented API contract (the previous approach here).

    Returns (gene_id, coords dict) for the winner, reusing the lookup_id
    call already made while comparing spans — the caller doesn't need to
    re-fetch it. coords is None if every lookup_id call failed, in which
    case the caller should fetch it fresh for the first candidate.

    At most max_gene_candidates ids are checked: this loop used to be unbounded,
    so a symbol matching N genes cost N serial round trips. Candidates are sorted
    before truncating so the choice stays reproducible. The cap can in principle
    miss a longer span beyond the cutoff, but symbols with more matches than that
    are LRG-style duplicates, which gene_candidates() has already filtered.
    """
    cap = ensembl_rest.load_config()["limits"]["max_gene_candidates"]
    scored = []
    for gid in sorted(candidates)[:cap]:
        try:
            # stop early rather than spend the remaining budget on ranking, and
            # fall back to the best span found so far
            http_retry.remaining(deadline)
            info = ensembl_rest.lookup_id(gid, deadline=deadline)
        except TimeoutError:
            break
        except Exception:
            continue
        span = info["end"] - info["start"] + 1
        scored.append((span, gid, info))
    if not scored:
        return candidates[0], None
    scored.sort(key=lambda t: t[0], reverse=True)
    _, gene_id, info = scored[0]
    return gene_id, info


def _fetch_ensembl(taxid, gene_symbol, flank_bp, reporter, deadline):
    """Ensembl REST: symbol -> gene id -> coordinates -> flanked FASTA; returns (fasta, meta)."""
    species = ensembl_rest.resolve_species_slug(taxid, deadline=deadline)
    if species is None:
        raise ValueError(f"Organism (taxid {taxid}) not found on Ensembl.")
    _report(reporter, f"Resolved taxid {taxid} -> {species}", stage="ensembl_species")

    xrefs = ensembl_rest.xref_symbol(species, gene_symbol, deadline=deadline)
    candidates = ensembl_rest.gene_candidates(xrefs)
    if not candidates:
        raise ValueError(f"No gene found for symbol '{gene_symbol}' in {species}.")

    ambiguous = len(candidates) > 1
    if ambiguous:
        gene_id, coords = _pick_most_complete_gene(candidates, deadline=deadline)
        if coords is None:
            coords = ensembl_rest.lookup_id(gene_id, deadline=deadline)
        _report(reporter, f"{gene_symbol} -> {gene_id} (chose 1 of {len(candidates)} "
                           f"matches by longest genomic span: {candidates})",
               stage="ensembl_xref")
    else:
        gene_id = candidates[0]
        coords = ensembl_rest.lookup_id(gene_id, deadline=deadline)
        _report(reporter, f"{gene_symbol} -> {gene_id}", stage="ensembl_xref")

    seq_region = coords["seq_region_name"]
    start      = coords["start"]
    end        = coords["end"]
    strand     = coords["strand"]
    assembly   = coords.get("assembly_name", "")
    _report(reporter, f"{gene_id}: {seq_region}:{start}-{end}:{strand} ({assembly})",
           stage="ensembl_lookup")

    fasta_text = ensembl_rest.fetch_region_fasta(
        species, seq_region, start, end, strand,
        expand_5prime=flank_bp, expand_3prime=flank_bp, deadline=deadline,
    )
    _report(reporter, f"fetched {len(fasta_text)} chars of genomic FASTA "
                       f"(+/- {flank_bp}bp flank)", stage="ensembl_fetch")

    meta = {
        "species":    species,
        "gene_id":    gene_id,
        "assembly":   assembly,
        "region":     f"{seq_region}:{start}-{end}:{strand}",
        "ambiguous":  ambiguous,
        "candidates": candidates if ambiguous else [],
        "source":     "ensembl",
    }
    return fasta_text, meta


def _strand_int(strand):
    """'+'/'-' -> 1/-1, Ensembl's convention, used in meta['region'] and FASTA headers."""
    return -1 if strand == "-" else 1


def _format_fasta(source, assembly, chrom, lo, hi, strand, seq, width=60):
    """FASTA text with an Ensembl-style header: >source:assembly:chrom:start:end:strand."""
    header = f">{source}:{assembly}:{chrom}:{lo}:{hi}:{_strand_int(strand)}"
    body = "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
    return f"{header}\n{body}\n"


def _fetch_ucsc(taxid, gene_symbol, flank_bp, reporter, deadline):
    """UCSC REST: symbol -> gene span + strand -> flanked sequence; returns (fasta, meta)."""
    cfg = ensembl_rest.load_config()
    genome = offtarget_blat.assembly_for_taxid(taxid, offtarget_screen.load_config())
    if not genome:
        raise ValueError(f"no UCSC assembly configured for taxid {taxid}")

    coords = ucsc_rest.find_gene(genome, gene_symbol, cfg, deadline=deadline)
    if coords is None:
        raise ValueError(f"No gene found for symbol '{gene_symbol}' in {genome}.")
    chrom, start, end, strand = coords["chrom"], coords["start"], coords["end"], coords["strand"]
    ambiguous = len(coords["candidates"]) > 1
    _report(reporter, f"{gene_symbol} -> {chrom}:{start}-{end}:{strand} ({genome}, "
                       f"{coords['track']})" + (f"; chose 1 of {len(coords['candidates'])} "
                       f"loci by longest span: {coords['candidates']}" if ambiguous else ""),
            stage="ucsc_lookup")

    seq, lo, hi = ucsc_rest.fetch_region(genome, chrom, start, end, strand, flank_bp,
                                         deadline=deadline)
    fasta_text = _format_fasta("ucsc", genome, chrom, lo, hi, strand, seq)
    _report(reporter, f"fetched {len(seq)} bp of genomic sequence (+/- {flank_bp}bp flank)",
            stage="ucsc_fetch")

    meta = {
        "species":    genome,
        "gene_id":    gene_symbol,
        "assembly":   genome,
        "region":     f"{chrom}:{start}-{end}:{_strand_int(strand)}",
        "ambiguous":  ambiguous,
        "candidates": coords["candidates"] if ambiguous else [],
        "source":     "ucsc",
    }
    return fasta_text, meta


def _fetch_ncbi(taxid, gene_symbol, flank_bp, reporter, deadline):
    """NCBI Datasets + efetch: symbol -> location -> flanked FASTA; returns (fasta, meta)."""
    cfg = ensembl_rest.load_config()
    coords = ncbi_rest.find_gene(taxid, gene_symbol, cfg, deadline=deadline)
    if coords is None:
        raise ValueError(f"No gene found for symbol '{gene_symbol}' in taxid {taxid} on NCBI.")
    acc, start, end, strand = coords["accession"], coords["start"], coords["end"], coords["strand"]
    ambiguous = len(coords["candidates"]) > 1
    _report(reporter, f"{gene_symbol} -> GeneID {coords['gene_id']} {acc}:{start}-{end}:{strand} "
                       f"({coords['assembly']})" + (f"; chose 1 of {len(coords['candidates'])} "
                       f"loci by longest span: {coords['candidates']}" if ambiguous else ""),
            stage="ncbi_lookup")

    seq, lo, hi = ncbi_rest.fetch_region(acc, start, end, strand, flank_bp, cfg,
                                         deadline=deadline)
    fasta_text = _format_fasta("ncbi", coords["assembly"], acc, lo, hi, strand, seq)
    _report(reporter, f"fetched {len(seq)} bp of genomic sequence (+/- {flank_bp}bp flank)",
            stage="ncbi_fetch")

    meta = {
        "species":    coords.get("taxname") or f"taxid {taxid}",
        "gene_id":    f"GeneID:{coords['gene_id']}",
        "assembly":   coords["assembly"],
        "region":     f"{acc}:{start}-{end}:{_strand_int(strand)}",
        "ambiguous":  ambiguous,
        "candidates": coords["candidates"] if ambiguous else [],
        "source":     "ncbi",
    }
    return fasta_text, meta


_FETCHERS = {"ncbi": _fetch_ncbi, "ucsc": _fetch_ucsc, "ensembl": _fetch_ensembl}


def _short_error(e):
    """One-line reason for a failed source, without urllib3's nested exception text."""
    if isinstance(e, requests.exceptions.HTTPError) and e.response is not None:
        return f"HTTP {e.response.status_code}"
    if isinstance(e, (requests.exceptions.Timeout, TimeoutError)):
        return "timed out"
    if isinstance(e, requests.exceptions.ConnectionError):
        return "connection failed"
    return str(e)


def source_settings():
    """(order, per_source_budget_s) from ensembl.config.json's sources block."""
    block = ensembl_rest.load_config().get("sources") or {}
    return list(block.get("order") or SOURCES), block.get("per_source_budget_s")


def fetch_genomic_sequence(taxid, gene_symbol, flank_bp=2000, report=None,
                           total_timeout=None, sources=None):
    """Fetch flanked genomic FASTA for a gene; returns (fasta_text, meta dict).

    Tries each source in turn (sources, else ensembl.config.json's sources.order)
    and returns the first success; meta['source'] names it and meta['fallbacks']
    lists what each earlier source reported. Raises ValueError when every source
    fails (the message lists each one) — callers should fall back to manual
    upload/paste. With a single source its own exception is re-raised unchanged.

    total_timeout is a wall-clock budget for the whole chain (default from
    ensembl.config.json); TimeoutError is raised once it is spent. Without it
    the per-request timeouts are per-socket only, so a slow-dribbling response
    never trips them — see issue #64.
    """
    reporter = resolve_reporter(report)
    if total_timeout is None:
        total_timeout = ensembl_rest.timeouts()["total_timeout_s"]
    deadline = http_retry.deadline_from(total_timeout)
    order, per_source = source_settings()
    if sources:
        order = list(sources)
    unknown = [s for s in order if s not in _FETCHERS]
    if unknown:
        raise ValueError(f"unknown genomic source(s): {unknown}; choose from {list(SOURCES)}")

    failures = []
    for i, source in enumerate(order):
        # every source but the last is capped, so a stalled one leaves time for the rest
        src_deadline = deadline
        if per_source and i < len(order) - 1:
            capped = http_retry.deadline_from(per_source)
            src_deadline = capped if deadline is None else min(deadline, capped)
        try:
            http_retry.remaining(deadline)
            fasta_text, meta = _FETCHERS[source](taxid, gene_symbol, flank_bp, reporter,
                                                 src_deadline)
        except Exception as e:
            if len(order) == 1:
                raise
            failures.append((source, e))
            _report(reporter, f"{source} failed ({_short_error(e)}); "
                               f"{'trying next source' if i < len(order) - 1 else 'no sources left'}",
                    stage="genomic_source", level="warning")
            # the overall budget is gone: later sources would only raise immediately
            if deadline is not None and http_retry.deadline_from(0) >= deadline:
                break
            continue
        meta["fallbacks"] = [f"{s}: {_short_error(e)}" for s, e in failures]
        return fasta_text, meta

    summary = "; ".join(f"{s}: {_short_error(e)}" for s, e in failures)
    # only a spent overall budget counts as a timeout; one slow source is just a failure
    if deadline is not None and http_retry.deadline_from(0) >= deadline:
        raise TimeoutError(f"timed out fetching genomic sequence for '{gene_symbol}' ({summary})")
    raise ValueError(f"Genomic sequence for '{gene_symbol}' (taxid {taxid}) not found by any "
                     f"source — {summary}")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Fetch flanked genomic FASTA for a gene via UCSC, NCBI or Ensembl REST."
    )
    parser.add_argument("--taxid", required=True,
                        help="NCBI taxonomy ID of the organism")
    parser.add_argument("--gene_symbol", required=True,
                        help="gene symbol to look up (e.g. unc-18)")
    parser.add_argument("--flank_bp", type=int, default=2000,
                        help="bp of flanking sequence to fetch on each side (default: 2000)")
    parser.add_argument("--output", required=True,
                        help="output FASTA file path")
    parser.add_argument("--timeout", type=float, default=None,
                        help="wall-clock budget in seconds for the whole fetch "
                             "(default: total_timeout_s from ensembl.config.json)")
    parser.add_argument("--source", default="auto", choices=("auto",) + SOURCES,
                        help="genomic source; auto tries sources.order from ensembl.config.json "
                             "(default: auto)")
    args, unknowns = parser.parse_known_args()

    fasta_text, meta = fetch_genomic_sequence(
        args.taxid, args.gene_symbol, args.flank_bp, total_timeout=args.timeout,
        sources=None if args.source == "auto" else [args.source],
    )

    print(f"source:   {meta['source']}", file=sys.stderr)
    for note in meta["fallbacks"]:
        print(f"skipped:  {note}", file=sys.stderr)
    print(f"species:  {meta['species']}", file=sys.stderr)
    print(f"gene_id:  {meta['gene_id']}", file=sys.stderr)
    print(f"region:   {meta['region']}", file=sys.stderr)
    print(f"assembly: {meta['assembly']}", file=sys.stderr)
    if meta["ambiguous"]:
        print(f"WARNING: {len(meta['candidates'])} gene matches for "
              f"'{args.gene_symbol}'; chose {meta['gene_id']} (longest genomic "
              f"span) — verify this is correct. All matches: {meta['candidates']}",
              file=sys.stderr)

    with open(args.output, "w") as f:
        f.write(fasta_text)
    print(f"Wrote genomic FASTA to {args.output}", file=sys.stderr)
