"""uniprot_api.py — shared UniProt REST helpers, usable without the Shiny app.

Extracted from existing_AF_model.py so uniprot_features.py can reuse the same
checksum-lookup logic without importing AFDB-specific code.

Matt Rich, 2026
"""

import requests
from Bio.SeqUtils.CheckSum import crc64


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
    resp = requests.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={"query": query, "format": "json", "fields": "accession", "size": 1},
        timeout=15,
    )
    resp.raise_for_status()
    return checksum_lookup_response_to_accession(resp.json())


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
