"""ensembl_rest.py — thin requests-based wrapper for the Ensembl REST API.

Used to resolve a gene symbol (in a given species) to genomic coordinates and
fetch flanked genomic FASTA sequence, without requiring the user to manually
upload a genomic file. No import-time network calls; safe to import anywhere.

Species resolution (taxid -> Ensembl species slug) is dynamic: FAST_PATH_SPECIES
is a cache seed for the app's known organisms (zero extra network calls), and
resolve_species_slug() falls back to Ensembl's /info/species discovery endpoint
for anything else, so any Ensembl-hosted organism works without a code change.

All request functions accept an optional deadline (from http_retry.deadline_from)
so the whole symbol -> FASTA operation is bounded in wall-clock time rather than
only per-socket; see ensembl.config.json and issue #64.
"""

import json
import sys
from pathlib import Path

import requests  # noqa: F401  — kept as the patch seam for tests

sys.path.insert(0, str(Path(__file__).parent))
import http_retry

CONFIG_PATH = Path(__file__).parent.parent / "ensembl.config.json"

_config_cache = None


def load_config(path=None):
    """Load ensembl.config.json (or a given path), cached after the first read."""
    global _config_cache
    if path is not None:
        with open(path) as f:
            return json.load(f)
    if _config_cache is None:
        with open(CONFIG_PATH) as f:
            _config_cache = json.load(f)
    return _config_cache


def timeouts():
    """The timeouts block from the config."""
    return load_config()["timeouts"]


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

# taxid -> slug, or None for "searched and not on Ensembl". The negative half is
# the point: without it, re-clicking Fetch for an unsupported organism re-walks
# every division each time. Both caches are module-level and so shared across
# sessions, which is intentional — this is immutable public reference data — and
# bounded by the division count and the number of distinct taxids seen.
_species_slug_cache = {}


def _get(path, params, read_timeout=None, deadline=None):
    """GET {BASE_URL}{path} through the shared retry/deadline wrapper.

    Single choke point so every Ensembl call gets the same bounded retries and
    the same clamping against the operation's remaining time budget.
    """
    cfg = timeouts()
    ret = load_config()["retries"]
    read = read_timeout if read_timeout is not None else cfg["read_timeout_s"]
    return http_retry.request_with_retries(
        "get", f"{BASE_URL}{path}",
        retries=ret["retries"], retry_wait=ret["retry_wait_s"], deadline=deadline,
        params=params, timeout=(cfg["connect_timeout_s"], read),
    )


def xref_symbol(species, symbol, deadline=None):
    """GET /xrefs/symbol/{species}/{symbol}; return list of {type, id} dicts."""
    resp = _get(f"/xrefs/symbol/{species}/{symbol}",
                {"content-type": "application/json"}, deadline=deadline)
    return resp.json()


def lookup_id(ensembl_id, deadline=None):
    """GET /lookup/id/{id}; return dict with seq_region_name, start, end, strand, assembly_name."""
    resp = _get(f"/lookup/id/{ensembl_id}",
                {"content-type": "application/json"}, deadline=deadline)
    return resp.json()


def fetch_region_fasta(species, seq_region, start, end, strand, expand_5prime=0,
                       expand_3prime=0, deadline=None):
    """GET /sequence/region/{species}/{region}; return FASTA text."""
    region = f"{seq_region}:{start}-{end}:{strand}"
    resp = _get(
        f"/sequence/region/{species}/{region}",
        {
            "expand_5prime": expand_5prime,
            "expand_3prime": expand_3prime,
            "content-type": "text/x-fasta",
        },
        read_timeout=timeouts()["fasta_read_timeout_s"], deadline=deadline,
    )
    return resp.text


def list_species(division=None, deadline=None):
    """GET /info/species (optionally filtered by division); return list of species dicts."""
    params = {"content-type": "application/json"}
    if division:
        params["division"] = division
    resp = _get("/info/species", params,
                read_timeout=timeouts()["fasta_read_timeout_s"], deadline=deadline)
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


def resolve_species_slug(taxid, deadline=None):
    """Resolve a taxid to an Ensembl species slug, or None if not found on Ensembl.

    Checks FAST_PATH_SPECIES first (no network call), then the slug cache, then
    searches each division's /info/species listing in turn, caching each
    division's full listing so repeated lookups within a process don't re-fetch
    it. The search stops after max_divisions and honors the deadline, so a stalled
    Ensembl can't turn this into a six-request wait.
    """
    taxid = int(taxid)
    if taxid in FAST_PATH_SPECIES:
        return FAST_PATH_SPECIES[taxid]
    if taxid in _species_slug_cache:
        return _species_slug_cache[taxid]

    taxid_str = str(taxid)
    for division in DIVISIONS[:load_config()["limits"]["max_divisions"]]:
        # raises TimeoutError rather than starting a request with no budget left
        http_retry.remaining(deadline)
        if division not in _division_cache:
            _division_cache[division] = list_species(division=division, deadline=deadline)
        matches = sorted(
            (sp["name"] for sp in _division_cache[division]
             if str(sp.get("taxon_id")) == taxid_str),
        )
        if matches:
            _species_slug_cache[taxid] = matches[0]
            return matches[0]
    _species_slug_cache[taxid] = None
    return None
