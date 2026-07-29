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
