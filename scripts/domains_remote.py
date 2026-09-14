"""
domains_remote.py

EBI InterProScan5 backend for domain annotation — moved out of call_interpro.py
verbatim as part of the local/remote backend split (see scripts/providers.py).
This is the default backend; a future domains_bulk.py will look up precomputed
InterPro matches for a bulk-downloaded proteome instead of submitting a queued job.

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py
"""

import os
import sys
from pathlib import Path

from site_selection_util import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, timed_poll_adapter


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Submit sequence to InterProScan5, parse TSV result, write output.

    job_id_cb(index, jid) and resume_job_ids follow the same convention as
    blast_orthologs.main(): this task makes one EBI submission, tagged index 0.
    """
    # imported lazily (not at module top) to avoid a load-time circular import —
    # call_interpro.py imports providers.py, which imports this module on demand
    from call_interpro import iprscan_tsv_to_domains

    reporter = resolve_reporter(report)
    name, seq = read_fasta(fasta_in)

    resume_id = (resume_job_ids or [None])[0]
    if resume_id:
        _report(reporter, "Checking previously-submitted InterProScan5 job…", stage="iprscan_submit")
        state, payload = ebi_rest.resume_job(ebi_rest.IPRSCAN5, resume_id, "tsv")
        if state == "pending":
            return {"ebi_status": "pending", "detail": payload}
        if state == "expired":
            return {"ebi_status": "expired", "detail": payload}
        tsv_bytes = payload
    else:
        params = {
            "email":    email,
            "stype":    "p",        # EBI iprscan5 uses 'p' for protein, not 'protein'
            "sequence": str(seq),   # Biopython Seq objects must be coerced to str
            "goterms":  "true",     # must be strings, not Python bools
            "pathways": "true",
        }

        _report(reporter, "Submitting InterProScan5 job…", stage="iprscan_submit")
        poll_cb = ebi_rest.combined_poll_cb(
            ebi_rest.indexed_job_id_cb(job_id_cb, 0),
            timed_poll_adapter(reporter, stage="iprscan_submit"),
        )
        job_id = ebi_rest.run_job(ebi_rest.IPRSCAN5, params, poll_cb=poll_cb)

        # fetch TSV result and save intermediate file (mirrors old naming: name.interpro.tsv.tsv)
        tsv_bytes = ebi_rest.fetch_result(ebi_rest.IPRSCAN5, job_id, "tsv")

    intermediate = os.path.join(workingdir, f"{name}.interpro.tsv.tsv")
    with open(intermediate, "wb") as f:
        f.write(tsv_bytes)

    domain_rows = iprscan_tsv_to_domains(tsv_bytes.decode())
    with open(outputfile, "w") as f_out:
        for row in domain_rows:
            print("\t".join(row), file=f_out)

    descriptions = sorted({row[3] for row in domain_rows if row[3]})
    if domain_rows:
        shown = ", ".join(descriptions[:5])
        if len(descriptions) > 5:
            shown += ", ..."
        summary = f"Found {len(domain_rows)} domain hit(s): {shown}"
    else:
        summary = "Found 0 domain hits"
    _report(reporter, summary, stage="done")
