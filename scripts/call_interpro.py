"""
call_interpro.py

CLI entry point for InterPro domain annotation. The actual lookup is delegated
to a backend resolved by scripts/providers.py — scripts/domains_remote.py (an
EBI iprscan5 job) by default, or a bulk-precomputed lookup when
batch.config.json sets backends.domains to a local mode. With no config file
present this always resolves to "remote", so the Shiny app's behavior here is
unchanged.

Parameters:
    - fasta_in (str): name of fasta file containing seq
    - email (str): EBI-registered email
    - workingdir (str): directory for intermediate files
    - clients_folder (str): (unused; retained for CLI compatibility)
    - outputfile (str): output filename (BED-like TSV)

Returns:
    - outputfile written with domain annotations

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py; backend-split 2026
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import resolve


def iprscan_tsv_to_domains(tsv_text):
    """Parse raw InterProScan5 TSV text into (source, start, stop, description) rows.

    Keeps only source(col3), start(col6), stop(col7), description(col5); rows with
    fewer than 8 tab-separated fields are malformed and skipped.
    """
    rows = []
    for line in tsv_text.splitlines():
        l = line.strip().split("\t")
        if len(l) < 8:
            continue
        rows.append((l[3], l[6], l[7], l[5]))
    return rows


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Run domain annotation via the configured backend (EBI InterProScan5 by
    default); same signature/return value as before the backend split.

    job_id_cb(index, jid) and resume_job_ids follow the same convention as
    blast_orthologs.main(): the remote backend makes one EBI submission, tagged
    index 0.
    """
    backend_main = resolve("domains")
    return backend_main(fasta_in, email, workingdir, clients_folder, outputfile,
                         report=report, job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", "--input_file", action="store", type=str, dest="FASTA_IN",
        help="name of fasta file containing seq.", required=True)
    parser.add_argument("--email", action="store", type=str, dest="EMAIL",
        help="email address, required by EBI job submission.", required=True)
    parser.add_argument("--dir", "--working_dir", action="store", type=str, dest="WORKINGDIR",
        help="working directory for output", required=True)
    parser.add_argument("--clients-folder", action="store", type=str, dest="CLIENTS_FOLDER",
        help="(unused; retained for CLI compatibility)", default="scripts/")
    parser.add_argument("--output", action="store", type=str, dest="OUTPUT",
        help="output file name")

    args, unknowns = parser.parse_known_args()

    main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.CLIENTS_FOLDER, args.OUTPUT)
