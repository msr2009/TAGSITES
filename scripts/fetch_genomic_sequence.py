"""
fetch_genomic_sequence.py

Automate genomic-sequence extraction for a gene, given an organism taxid and a
gene symbol, using the Ensembl REST API. Replaces manual genomic FASTA
upload/paste for any organism Ensembl hosts (species resolution is dynamic —
see ensembl_rest.resolve_species_slug()).

Pipeline: symbol -> gene id (xref_symbol) -> coordinates (lookup_id) ->
flanked genomic FASTA (fetch_region_fasta). flank_bp of extra sequence is
fetched on each side of the gene body so downstream reagent design
(scripts/design_tag_reagents.py) has enough flanking sequence for homology
arms; the default (2000) covers the Design Reagents tab's max allowed
arm_length (2000bp).

Matt Rich, 2026
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ensembl_rest
import http_retry
from progress import report as _report, resolve_reporter


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


def fetch_genomic_sequence(taxid, gene_symbol, flank_bp=2000, report=None,
                           total_timeout=None):
    """Fetch flanked genomic FASTA for a gene; returns (fasta_text, meta dict).

    Raises ValueError if the organism isn't on Ensembl or the gene symbol
    isn't found — callers should fall back to manual upload/paste.

    total_timeout is a wall-clock budget for the whole operation (default from
    ensembl.config.json); TimeoutError is raised once it is spent. Without it
    the per-request timeouts are per-socket only, so a slow-dribbling response
    never trips them — see issue #64.
    """
    reporter = resolve_reporter(report)
    if total_timeout is None:
        total_timeout = ensembl_rest.timeouts()["total_timeout_s"]
    deadline = http_retry.deadline_from(total_timeout)

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
    }
    return fasta_text, meta


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Fetch flanked genomic FASTA for a gene via the Ensembl REST API."
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
    args, unknowns = parser.parse_known_args()

    fasta_text, meta = fetch_genomic_sequence(args.taxid, args.gene_symbol, args.flank_bp,
                                              total_timeout=args.timeout)

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
