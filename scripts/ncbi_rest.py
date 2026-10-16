"""ncbi_rest.py — gene symbol -> flanked genomic FASTA via NCBI Datasets v2 + E-utilities efetch.

First source in fetch_genomic_sequence's fallback chain (the fastest), and the broadest: any
organism with an annotated NCBI genome, including E. coli and non-preset taxa.

  Datasets /gene/symbol/{symbol}/taxon/{taxid}   -> accession, range, orientation
  efetch db=nuccore seq_start/seq_stop/strand      -> flanked FASTA, already in the
                                                      gene's orientation (strand=2
                                                      returns the reverse complement)

Datasets rejects taxa above species (E. coli 562), so ncbi.taxid_aliases in
ensembl.config.json maps those to a reference strain. An optional NCBI_API_KEY
from the environment raises the E-utilities rate limit; it is never stored.

Matt Rich, 2026
"""

import os
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, str(Path(__file__).parent))
import http_retry

DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/{taxid}"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
KEY_ENV = "NCBI_API_KEY"


def _get(url, params, cfg, deadline, read_timeout=None):
    """GET an NCBI endpoint through the shared retry/deadline wrapper; return the Response."""
    t = cfg["timeouts"]
    ret = cfg["retries"]
    key = os.environ.get(KEY_ENV, "").strip()
    if key:
        params = dict(params, api_key=key)
    return http_retry.request_with_retries(
        "get", url, retries=ret["retries"], retry_wait=ret["retry_wait_s"], deadline=deadline,
        params=params, timeout=(t["connect_timeout_s"], read_timeout or t["read_timeout_s"]),
    )


def resolve_taxid(taxid, cfg):
    """Apply ncbi.taxid_aliases (e.g. E. coli 562 -> K-12 MG1655 511145)."""
    aliases = (cfg.get("ncbi") or {}).get("taxid_aliases") or {}
    return str(aliases.get(str(taxid), taxid))


def gene_locations(taxid, symbol, cfg, deadline=None):
    """Every genomic location of an exact symbol match on each gene's first annotation.

    The first annotation is the reference assembly (GRCh38 ahead of T2T, the
    newest zebrafish build first). Returns a list of dicts: gene_id, accession,
    chrom, start, end (1-based inclusive), strand ('+'/'-'), assembly, taxname.
    """
    url = DATASETS_URL.format(symbol=quote(symbol, safe=""), taxid=resolve_taxid(taxid, cfg))
    payload = _get(url, {}, cfg, deadline).json()
    out = []
    for report in payload.get("reports") or []:
        gene = report.get("gene") or {}
        # the endpoint can return synonym matches; keep only the symbol itself
        if str(gene.get("symbol", "")).lower() != symbol.lower():
            continue
        annotations = gene.get("annotations") or []
        if not annotations:
            continue
        ann = annotations[0]
        for loc in ann.get("genomic_locations") or []:
            rng = loc.get("genomic_range") or {}
            if not rng.get("begin") or not loc.get("genomic_accession_version"):
                continue
            out.append({
                "gene_id":   str(gene.get("gene_id", "")),
                "accession": loc["genomic_accession_version"],
                "chrom":     loc.get("sequence_name") or loc["genomic_accession_version"],
                "start":     int(rng["begin"]),
                "end":       int(rng["end"]),
                "strand":    "-" if rng.get("orientation") == "minus" else "+",
                "assembly":  ann.get("assembly_name", ""),
                "taxname":   gene.get("taxname", ""),
            })
    return out


def find_gene(taxid, symbol, cfg, deadline=None):
    """Resolve a symbol to its longest genomic location; returns a coords dict or None.

    coords is one gene_locations() entry plus candidates, every location as
    'gene_id@accession:start-end', longest first.
    """
    locs = gene_locations(taxid, symbol, cfg, deadline=deadline)
    if not locs:
        return None
    # longest genomic span wins, as in fetch_genomic_sequence._pick_most_complete_gene
    locs.sort(key=lambda d: (-(d["end"] - d["start"]), d["gene_id"], d["accession"]))
    best = dict(locs[0])
    best["candidates"] = [f"{d['gene_id']}@{d['accession']}:{d['start']}-{d['end']}"
                          for d in locs]
    return best


def fetch_region(accession, start, end, strand, flank, cfg, deadline=None):
    """Flanked region as a sequence string, in the gene's orientation; returns (seq, start, end).

    start/end are 1-based inclusive; the returned pair is the flanked region,
    clamped at the sequence start (efetch clamps the far end itself).
    """
    lo = max(1, start - flank)
    hi = end + flank
    resp = _get(EFETCH_URL, {
        "db": "nuccore", "id": accession, "seq_start": lo, "seq_stop": hi,
        "strand": 2 if strand == "-" else 1, "rettype": "fasta", "retmode": "text",
        "tool": "tagsites",
    }, cfg, deadline, read_timeout=cfg["timeouts"]["fasta_read_timeout_s"])
    lines = resp.text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"NCBI efetch returned no FASTA for {accession}:{lo}-{hi}")
    seq = "".join(lines[1:]).upper()
    if not seq:
        raise RuntimeError(f"NCBI efetch returned an empty sequence for {accession}:{lo}-{hi}")
    return seq, lo, lo + len(seq) - 1
