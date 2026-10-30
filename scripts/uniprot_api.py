"""uniprot_api.py — shared UniProt REST helpers, usable without the Shiny app.

Extracted from existing_AF_model.py so uniprot_features.py can reuse the same
checksum-lookup logic without importing AFDB-specific code. resolve_accession() is the one
place a run's protein is mapped to a UniProt accession (uniprot_features.py,
ucsc_gene_model.py, phastcons_track.py).

Matt Rich, 2026
"""

import os
import re

import requests
from Bio.SeqUtils.CheckSum import crc64

from http_retry import request_with_retries

ACCESSION_SIDECAR = "{run}.uniprot_accession.txt"  # written by structure_remote on a match


def checksum_lookup_response_to_accession(resp_json):
    """Extract the first matching UniProt accession from a checksum-search response, or None."""
    results = resp_json.get("results", [])
    if results:
        return results[0]["primaryAccession"]
    return None


def checksum_lookup(seq, taxid=None):
    """Query UniProt REST API for an exact sequence match by CRC64 checksum.

    Returns the first matching UniProt accession, or None.
    """
    checksum = crc64(seq).replace("CRC-", "")  # Biopython adds "CRC-" prefix; UniProt expects bare hex
    query = f"checksum:{checksum}"
    if str(taxid) not in ("", "1", "1.0", "None", None):
        query += f" AND taxonomy_id:{taxid}"
    # batch-scale traffic sees transient 5xx/429 here; retrying avoids a needless BLAST fallback
    resp = request_with_retries(
        "get", "https://rest.uniprot.org/uniprotkb/search",
        params={"query": query, "format": "json", "fields": "accession", "size": 1},
        timeout=15,
    )
    return checksum_lookup_response_to_accession(resp.json())


def parse_accession(text):
    """UniProt accession (isoform suffix kept, e.g. 'Q26261-2') at the start of text, or ''.

    Accepts a bare accession or a UniProt FASTA id ('sp|Q26261|UNC5_CAEEL').
    """
    from site_selection_util import uniprot_accession_regex

    text = (text or "").strip().split()[0] if (text or "").strip() else ""
    parts = text.split("|")
    # sp|ACC|NAME / tr|ACC|NAME: the accession is the middle field
    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
        text = parts[1]
    m = uniprot_accession_regex(text)
    if m is None:
        return ""
    iso = re.match(r"-\d+", text[m.end():])
    return m.group(0) + (iso.group(0) if iso else "")


def _fasta_header(path):
    """First FASTA header (without '>') of a file, or '' if unreadable."""
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    return line[1:].strip()
    except OSError:
        pass
    return ""


def resolve_accession(protein_fa, accession="", taxid="", run_dir=None, run_name=None,
                      use_checksum=True):
    """Best UniProt accession for a run's protein; returns (accession, how) or ('', reason).

    Order: explicit accession > the FASTA header > user_{name}.fa's header (the AFDB
    pre-search replaces {run}.fa with a run-named copy) > CRC64 checksum lookup >
    the {run}.uniprot_accession.txt sidecar structure_remote writes on a BLAST/CRC match.
    """
    acc = parse_accession(accession)
    if acc:
        return acc, "given"
    acc = parse_accession(_fasta_header(protein_fa)) if protein_fa else ""
    if acc:
        return acc, "FASTA header"
    folder = run_dir or (os.path.dirname(protein_fa) if protein_fa else "")
    stem = os.path.splitext(os.path.basename(protein_fa or ""))[0]
    stem = stem.removeprefix("user_")
    name = run_name or stem
    user_fa = os.path.join(folder, f"user_{stem}.fa")
    if os.path.exists(user_fa):
        acc = parse_accession(_fasta_header(user_fa))
        if acc:
            return acc, "original FASTA header"
    if use_checksum and protein_fa and os.path.exists(protein_fa):
        from cds_check import read_protein

        seq = read_protein(protein_fa)
        try:
            acc = checksum_lookup(seq.rstrip("*"), taxid) if seq else None
        except Exception:
            acc = None
        if acc:
            return acc, "CRC64 checksum"
    sidecar = os.path.join(folder, ACCESSION_SIDECAR.format(run=name)) if name else ""
    if sidecar and os.path.exists(sidecar):
        with open(sidecar) as fh:
            acc = parse_accession(fh.read())
        if acc:
            return acc, "AFDB search match"
    return "", "no UniProt accession found (header, checksum or BLAST match)"


def fetch_entry(accession):
    """Fetch the full UniProt entry (JSON) for an accession. Raises on HTTP error."""
    resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=15)
    resp.raise_for_status()
    return resp.json()


def fetch_fasta(accession):
    """Fetch a single accession's sequence as a plain string. Raises on HTTP error."""
    resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", timeout=15)
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    return "".join(lines[1:])


def _protein_name(entry):
    """Best-effort protein name from a UniProt entry: recommended name, else first
    submitted name (TrEMBL entries usually lack a recommendedName), else ''."""
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if rec:
        return rec
    submitted = desc.get("submittedName", [])
    if submitted:
        return submitted[0].get("fullName", {}).get("value", "")
    return ""


def _gene_id_xref(entry):
    """NCBI GeneID cross-reference from a UniProt entry, or None. GeneID is organism-
    agnostic and present for nearly every annotated gene, making it a reliable key for
    finding sibling entries regardless of species-specific databases like WormBase."""
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GeneID":
            return xref.get("id")
    return None


def fetch_computationally_mapped_isoforms(base_acc, gene_id):
    """Return UniProt's "computationally mapped potential isoform sequences" for a gene.

    These are Ensembl-predicted transcripts imported into TrEMBL as their own full
    entries (distinct primary accessions, not "-N" suffixes) because they were never
    curated into an ALTERNATIVE PRODUCTS comment. Found by searching for other entries
    sharing base_acc's GeneID cross-reference. Returns [] if gene_id is falsy or the
    search fails — callers should treat that as "none found", not an error.
    """
    if not gene_id:
        return []
    try:
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f"xref:GeneID-{gene_id}", "format": "json",
                    "fields": "accession,protein_name,sequence", "size": 50},
            timeout=15,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])
    except Exception:
        return []

    records = []
    for entry in results:
        acc = entry.get("primaryAccession")
        seq = entry.get("sequence", {}).get("value", "")
        if not acc or acc == base_acc or not seq:
            continue
        records.append({
            "accession":    acc,
            "isoform":      acc,
            "isoform_name": _protein_name(entry) or acc,
            "sequence":     seq,
            "is_canonical": False,
            "computationally_mapped": True,
        })
    return records


def fetch_isoform_sequences(accession):
    """Return per-isoform sequence records for a UniProt gene: curated ALTERNATIVE
    PRODUCTS isoforms plus computationally mapped potential isoforms (Ensembl-predicted
    TrEMBL entries sharing the same GeneID — see fetch_computationally_mapped_isoforms).

    accession may carry an isoform suffix (e.g. "P04637-2"); it is stripped since the
    isoform-suffixed entry endpoint returns no comments. Returns [] only on a fetch
    error for the base entry itself — callers should fall back to another
    isoform-detection strategy in that case. A single-isoform gene with no
    computationally mapped siblings returns a length-1 list (canonical only).
    """
    base_acc = accession.split("-")[0]
    try:
        entry = fetch_entry(base_acc)
    except Exception:
        return []

    canonical_seq = entry.get("sequence", {}).get("value", "")
    canonical_acc = base_acc
    canonical_name = _protein_name(entry) or base_acc

    records = []
    alt = next(
        (c for c in entry.get("comments", []) if c.get("commentType") == "ALTERNATIVE PRODUCTS"),
        None,
    )
    if alt and alt.get("isoforms"):
        for iso in alt["isoforms"]:
            iso_ids = iso.get("isoformIds", [])
            if not iso_ids:
                continue
            iso_acc  = iso_ids[0]  # e.g. "P04637-2"
            synonyms = [s["value"] for s in iso.get("synonyms", [])]
            iso_name = synonyms[0] if synonyms else iso.get("name", {}).get("value", "")
            status   = iso.get("isoformSequenceStatus", "")
            is_canonical = status == "Displayed"

            if is_canonical:
                seq = canonical_seq
                canonical_acc, canonical_name = iso_acc, iso_name or canonical_name
            else:
                try:
                    seq = fetch_fasta(iso_acc)
                except Exception:
                    continue

            records.append({
                "accession":    iso_acc,
                "isoform":      iso_acc.rsplit("-", 1)[-1] if "-" in iso_acc else "1",
                "isoform_name": iso_name,
                "sequence":     seq,
                "is_canonical": is_canonical,
            })

    if not any(r["is_canonical"] for r in records):
        records.insert(0, {
            "accession":    canonical_acc,
            "isoform":      "1",
            "isoform_name": canonical_name,
            "sequence":     canonical_seq,
            "is_canonical": True,
        })

    seen_seqs = {r["sequence"] for r in records}
    gene_id = _gene_id_xref(entry)
    for rec in fetch_computationally_mapped_isoforms(base_acc, gene_id):
        if rec["sequence"] in seen_seqs:
            continue  # already covered by a curated ALTERNATIVE PRODUCTS isoform
        seen_seqs.add(rec["sequence"])
        records.append(rec)

    return records
