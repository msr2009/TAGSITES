"""
domains_bulk.py

Local backend for domain annotation: looks up precomputed InterPro matches
in local_store.sqlite3 (built by scripts/local_store.py --build-domains from
scripts/reference_data.py's protein2ipr.filtered.tsv) instead of submitting an
EBI iprscan5 job. Falls back to raising FileNotFoundError with a clear message
if the local index hasn't been built — scripts/providers.py has no automatic
fallback to domains_remote itself; that's the batch driver's job (Phase C),
per batch.config.json's documented "config with fallback" model.

protein2ipr.dat has no human-readable source-database name column — only the
member database's own signature accession (e.g. "PF10325" for Pfam,
"PTHR23021" for PANTHER). _PREFIX_TO_SOURCE maps the well-known InterPro
member-database signature prefixes to the same source names
domains_remote.py's EBI TSV parsing already produces, so downstream code
(utils/results.py) doesn't need to know which backend produced a row. Only
"Pfam" rows are actually displayed by the app today (see
utils/results.py's _FEATURE_SOURCES / _DOMAIN_STRUCT_SOURCES) — everything
else is written through for completeness (domain-count summaries, future use)
but won't appear in the structure/domain panels either way.
"""

import sys
from pathlib import Path

from site_selection_util import read_fasta, uniprot_accession_regex

sys.path.insert(0, str(Path(__file__).parent))
from local_store import lookup_by_crc64, lookup_domains_by_accession, open_index
from progress import report as _report, resolve_reporter

# InterPro member-database signature-accession prefixes -> source name.
# Order matters: check longer/more-specific prefixes before shorter ones.
_PREFIX_TO_SOURCE = [
    ("PF",      "Pfam"),
    ("PTHR",    "PANTHER"),
    ("SSF",     "SUPERFAMILY"),
    ("PS",      "PROSITE"),
    ("SM",      "SMART"),
    ("PR",      "PRINTS"),
    ("PIRSF",   "PIRSF"),
    ("cd",      "CDD"),
    ("G3DSA:",  "Gene3D"),
    ("MF_",     "HAMAP"),
    ("TIGR",    "TIGRFAM"),
    ("NF",      "NCBIfam"),
    ("PD",      "ProDom"),
]


def _source_for(external_db_match_id):
    """Map an InterPro member-database signature accession to its source name;
    falls back to "InterPro" for prefixes not in the table above."""
    for prefix, source in _PREFIX_TO_SOURCE:
        if external_db_match_id.startswith(prefix):
            return source
    return "InterPro"


def _resolve_accession(fasta_in, conn):
    """Resolve fasta_in to a UniProt accession: the FASTA record id if it
    already looks like one, else a CRC64 checksum lookup against the local
    proteins index. Returns None if neither resolves.
    """
    name, seq = read_fasta(fasta_in)
    if uniprot_accession_regex(name):
        return name
    from Bio.SeqUtils.CheckSum import crc64
    checksum = crc64(str(seq)).replace("CRC-", "")
    matches = lookup_by_crc64(checksum, conn=conn)
    return matches[0]["accession"] if matches else None


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Look up precomputed InterPro domain matches for fasta_in's protein and
    write the same (source, start, stop, description) TSV domains_remote.py
    produces. email/job_id_cb/resume_job_ids are accepted but unused (no
    network job is submitted); kept so providers.resolve("domains") can call
    either backend with the identical call.
    """
    reporter = resolve_reporter(report)
    conn = open_index()
    try:
        accession = _resolve_accession(fasta_in, conn)
        if accession is None:
            _report(reporter,
                    "Could not resolve a UniProt accession for this sequence "
                    "in the local index — no domain annotation available "
                    "from the bulk backend.", stage="domains_bulk", level="warning")
            with open(outputfile, "w"):
                pass
            return

        rows = lookup_domains_by_accession(accession, conn=conn)
    finally:
        conn.close()

    domain_rows = [
        (_source_for(r["external_db_match_id"]), str(r["start"]), str(r["stop"]), r["description"])
        for r in rows
    ]
    with open(outputfile, "w") as f_out:
        for row in domain_rows:
            print("\t".join(row), file=f_out)

    descriptions = sorted({row[3] for row in domain_rows if row[3]})
    if domain_rows:
        shown = ", ".join(descriptions[:5])
        if len(descriptions) > 5:
            shown += ", ..."
        summary = f"Found {len(domain_rows)} domain hit(s) for {accession} (local index): {shown}"
    else:
        summary = f"Found 0 domain hits for {accession} (local index)"
    _report(reporter, summary, stage="done")
