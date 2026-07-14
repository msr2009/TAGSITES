"""ensembl_rest.py — thin requests-based wrapper for the Ensembl REST API.

Used to resolve a gene symbol (in a given species) to genomic coordinates and
fetch flanked genomic FASTA sequence, without requiring the user to manually
upload a genomic file. No import-time network calls; safe to import anywhere.

Species resolution (taxid -> Ensembl species slug) is dynamic: FAST_PATH_SPECIES
is a cache seed for the app's known organisms (zero extra network calls), and
resolve_species_slug() falls back to Ensembl's /info/species discovery endpoint
for anything else, so any Ensembl-hosted organism works without a code change.
"""

import requests

BASE_URL = "https://rest.ensembl.org"

# division search order for dynamic species resolution — Vertebrates first
# since it's the default (no division param) call and covers most requests
DIVISIONS = [
    "EnsemblVertebrates", "EnsemblMetazoa", "EnsemblFungi",
    "EnsemblPlants", "EnsemblProtists", "EnsemblBacteria",
]

# cache seed, not an allowlist — covers the app's preset organisms with no
# network round-trip; add a new organism here directly, or let
# resolve_species_slug() find it dynamically via /info/species
FAST_PATH_SPECIES = {
    7227:   "drosophila_melanogaster",
    6239:   "caenorhabditis_elegans",
    559292: "saccharomyces_cerevisiae",
    7955:   "danio_rerio",
    10090:  "mus_musculus",
    10116:  "rattus_norvegicus",
    8364:   "xenopus_tropicalis",
    9606:   "homo_sapiens",
    562:    "escherichia_coli_str_k_12_substr_mg1655_gca_000005845",
}

# per-division species listings, populated lazily by list_species()
_division_cache = {}


def xref_symbol(species, symbol):
    """GET /xrefs/symbol/{species}/{symbol}; return list of {type, id} dicts."""
    resp = requests.get(
        f"{BASE_URL}/xrefs/symbol/{species}/{symbol}",
        params={"content-type": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json()


def lookup_id(ensembl_id):
    """GET /lookup/id/{id}; return dict with seq_region_name, start, end, strand, assembly_name."""
    resp = requests.get(
        f"{BASE_URL}/lookup/id/{ensembl_id}",
        params={"content-type": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json()


def fetch_region_fasta(species, seq_region, start, end, strand, expand_5prime=0, expand_3prime=0):
    """GET /sequence/region/{species}/{region}; return FASTA text."""
    region = f"{seq_region}:{start}-{end}:{strand}"
    resp = requests.get(
        f"{BASE_URL}/sequence/region/{species}/{region}",
        params={
            "expand_5prime": expand_5prime,
            "expand_3prime": expand_3prime,
            "content-type": "text/x-fasta",
        },
        timeout=60,
    )
    resp.raise_for_status()
    return resp.text


def list_species(division=None):
    """GET /info/species (optionally filtered by division); return list of species dicts."""
    params = {"content-type": "application/json"}
    if division:
        params["division"] = division
    resp = requests.get(f"{BASE_URL}/info/species", params=params, timeout=60)
    resp.raise_for_status()
    return resp.json().get("species", [])


def gene_candidates(xref_json):
    """Extract all type=='gene' ids from an xref_symbol response, LRG_* duplicates filtered out.

    LRG_ ids (Locus Reference Genomic) are a curated duplicate identifier for the
    same gene, not a distinct locus — observed in practice for several human
    genes (BRCA1, TP53, EGFR, ...). Filtering them out here resolves most
    real-world multi-match cases with no extra network call.
    """
    genes = [e["id"] for e in xref_json if e.get("type") == "gene" and e.get("id")]
    real = [g for g in genes if not g.startswith("LRG_")]
    return real or genes


def first_gene_id(xref_json):
    """Extract the best type=='gene' entry's id from an xref_symbol response, or None."""
    candidates = gene_candidates(xref_json)
    return candidates[0] if candidates else None


def resolve_species_slug(taxid):
    """Resolve a taxid to an Ensembl species slug, or None if not found on Ensembl.

    Checks FAST_PATH_SPECIES first (no network call); otherwise searches each
    division's /info/species listing in turn, caching each division's full
    listing so repeated lookups within a process don't re-fetch it.
    """
    taxid = int(taxid)
    if taxid in FAST_PATH_SPECIES:
        return FAST_PATH_SPECIES[taxid]

    taxid_str = str(taxid)
    for division in DIVISIONS:
        if division not in _division_cache:
            _division_cache[division] = list_species(division=division)
        matches = sorted(
            (sp["name"] for sp in _division_cache[division]
             if str(sp.get("taxon_id")) == taxid_str),
        )
        if matches:
            return matches[0]
    return None
