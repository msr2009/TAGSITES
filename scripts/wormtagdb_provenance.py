"""
wormtagdb_provenance.py

Fetches and caches source evidence (CGC strain pages, paper full text) for rows
in a wormtagdb internal-tags export, so a human/agent can annotate the actual
tag insertion position by reading the cached text. This script only does the
mechanical fetch/cache/tracking step; extracting a position from prose is a
separate manual annotation pass.

Parameters:
    - input_xlsx (str): path to wormtagdb internal-tags .xlsx export
    - output_dir (str): directory for cache files and the working CSV
    - limit (int or None): if set, take a pilot subset instead of all rows

Returns:
    - working CSV written to <output_dir>/wormtagdb_provenance.csv
    - raw CGC/paper text cached under <output_dir>/cache/

Matt Rich, 7/2026
"""

import os
import re
import csv

import pandas as pd
import requests

CGC_STRAIN_URL = "https://cgc.umn.edu/strain/{}"
EUROPEPMC_FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML"

TASK_COLUMNS = [
    "strain", "gene", "allele", "wormbase_id",
    "cgc_fetched", "cgc_cache_path",
    "pmcid", "pub_fetched", "pub_cache_path",
    "tier", "position_estimate", "confidence", "evidence_quote", "source_url", "notes",
]


def load_wormtagdb(xlsx_path):
    """Read the wormtagdb internal-tags sheet into a DataFrame."""
    df = pd.read_excel(xlsx_path, na_values=[], keep_default_na=False)
    return df


def parse_publication_ids(accession_str):
    """Pull PMID/PMCID/DOI out of an Accession cell like 'DOI:... PMCID:... PMID:...'."""
    if not isinstance(accession_str, str):
        return {"pmid": None, "pmcid": None, "doi": None}
    pmid = re.search(r"PMID:(\S+)", accession_str)
    pmcid = re.search(r"PMCID:(\S+)", accession_str)
    doi = re.search(r"DOI:(\S+)", accession_str)
    return {
        "pmid": pmid.group(1) if pmid else None,
        "pmcid": pmcid.group(1) if pmcid else None,
        "doi": doi.group(1) if doi else None,
    }


def fetch_cgc_page(strain_id, cache_dir):
    """GET the CGC strain page and cache it; returns (fetched, cache_path)."""
    if not strain_id or strain_id == "NA":
        return False, None
    cache_path = os.path.join(cache_dir, f"cgc_{strain_id}.html")
    if os.path.exists(cache_path):
        return True, cache_path
    resp = requests.get(CGC_STRAIN_URL.format(strain_id), timeout=30)
    if resp.status_code != 200:
        return False, None
    with open(cache_path, "w") as f:
        f.write(resp.text)
    return True, cache_path


def fetch_paper_fulltext(pmcid, cache_dir):
    """GET Europe PMC full-text XML for a PMCID and cache it; returns (fetched, cache_path)."""
    if not pmcid:
        return False, None
    cache_path = os.path.join(cache_dir, f"paper_{pmcid}.xml")
    if os.path.exists(cache_path):
        return True, cache_path
    resp = requests.get(EUROPEPMC_FULLTEXT_URL.format(pmcid), timeout=30)
    if resp.status_code != 200 or not resp.text.strip():
        return False, None
    with open(cache_path, "w") as f:
        f.write(resp.text)
    return True, cache_path


def select_pilot_rows(df, limit):
    """Pick a pilot subset mixing CGC-sourced and non-CGC (NA strain) rows."""
    if limit is None or limit >= len(df):
        return df
    cgc_rows = df[df["Strain"] != "NA"]
    other_rows = df[df["Strain"] == "NA"]
    n_cgc = min(len(cgc_rows), limit // 2)
    n_other = limit - n_cgc
    return pd.concat([cgc_rows.head(n_cgc), other_rows.head(n_other)])


def build_task_table(df, cache_dir):
    """Fetch CGC/paper evidence for each row and assemble the working task table."""
    rows_out = []
    for _, row in df.iterrows():
        strain = row.get("Strain")
        ids = parse_publication_ids(row.get("Accession"))

        cgc_fetched, cgc_path = fetch_cgc_page(strain, cache_dir)
        pub_fetched, pub_path = fetch_paper_fulltext(ids["pmcid"], cache_dir)

        rows_out.append({
            "strain": strain,
            "gene": row.get("Gene"),
            "allele": row.get("Allele"),
            "wormbase_id": row.get("WormBase ID"),
            "cgc_fetched": cgc_fetched,
            "cgc_cache_path": cgc_path,
            "pmcid": ids["pmcid"],
            "pub_fetched": pub_fetched,
            "pub_cache_path": pub_path,
            "tier": "",
            "position_estimate": "",
            "confidence": "",
            "evidence_quote": "",
            "source_url": "",
            "notes": "",
        })
    return rows_out


def main(input_xlsx, output_dir, limit=None):
    cache_dir = os.path.join(output_dir, "cache")
    os.makedirs(cache_dir, exist_ok=True)

    df = load_wormtagdb(input_xlsx)
    df = select_pilot_rows(df, limit)

    task_rows = build_task_table(df, cache_dir)

    out_csv = os.path.join(output_dir, "wormtagdb_provenance.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TASK_COLUMNS)
        writer.writeheader()
        writer.writerows(task_rows)

    return out_csv


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--input", action="store", type=str, dest="INPUT_XLSX",
        help="path to wormtagdb internal-tags .xlsx export", required=True)
    parser.add_argument("--output-dir", action="store", type=str, dest="OUTPUT_DIR",
        help="directory for cache files and the working CSV", required=True)
    parser.add_argument("--limit", action="store", type=int, dest="LIMIT",
        help="take a pilot subset of N rows instead of all rows", default=None)

    args = parser.parse_args()

    main(args.INPUT_XLSX, args.OUTPUT_DIR, args.LIMIT)
