"""
call_interpro.py

A script to perform an InterPro search via the EBI iprscan5 REST API.

Parameters:
    - fasta_in (str): name of fasta file containing seq
    - email (str): EBI-registered email
    - workingdir (str): directory for intermediate files
    - clients_folder (str): (unused; retained for CLI compatibility)
    - outputfile (str): output filename (BED-like TSV)

Returns:
    - outputfile written with domain annotations

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py
"""

import os
import sys
from pathlib import Path

from site_selection_util import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, poll_adapter


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None):
    """Submit sequence to InterProScan5, parse TSV result, write output."""
    reporter = resolve_reporter(report)
    name, seq = read_fasta(fasta_in)

    params = {
        "email":    email,
        "stype":    "p",        # EBI iprscan5 uses 'p' for protein, not 'protein'
        "sequence": str(seq),   # Biopython Seq objects must be coerced to str
        "goterms":  "true",     # must be strings, not Python bools
        "pathways": "true",
    }

    _report(reporter, "Submitting InterProScan5 job…", stage="iprscan_submit")
    job_id = ebi_rest.run_job(ebi_rest.IPRSCAN5, params, poll_cb=poll_adapter(reporter))

    # fetch TSV result and save intermediate file (mirrors old naming: name.interpro.tsv.tsv)
    tsv_bytes = ebi_rest.fetch_result(ebi_rest.IPRSCAN5, job_id, "tsv")
    intermediate = os.path.join(workingdir, f"{name}.interpro.tsv.tsv")
    with open(intermediate, "wb") as f:
        f.write(tsv_bytes)

    # parse and reformat: keep source(col3), start(col6), stop(col7), description(col5)
    with open(outputfile, "w") as f_out:
        for line in open(intermediate):
            l = line.strip().split("\t")
            if len(l) < 8:
                continue
            print("\t".join([l[3], l[6], l[7], l[5]]), file=f_out)


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
