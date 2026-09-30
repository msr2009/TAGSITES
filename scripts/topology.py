"""
topology.py

CLI entry point for protein topology annotation (transmembrane helices,
signal peptides, inside/outside orientation). Delegates to a backend resolved
by scripts/providers.py — scripts/topology_deeptmhmm.py (local DeepTMHMM) by
default, since no remote topology backend exists yet. A future
topology_remote.py could wrap a Phobius-only EBI submission without any
caller here changing.

Not part of task_definitions.json's default_tasks — DeepTMHMM is far more
expensive per protein than the other default analyses (an ESM1b transformer
pass plus a 5-model CRF ensemble), so it is opt-in.

Parameters:
    - fasta_in (str): name of fasta file containing seq
    - email (str): unused by the local backend; kept for CLI/signature parity
    - workingdir (str): directory for intermediate files
    - clients_folder (str): (unused; retained for CLI compatibility)
    - outputfile (str): output filename (BED-like TSV)

Returns:
    - outputfile written with topology region annotations
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import resolve


def main(fasta_in, email, workingdir, clients_folder, outputfile, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Run topology annotation via the configured backend (local DeepTMHMM by
    default); same signature as every other task script.
    """
    backend_main = resolve("topology", default="deeptmhmm")
    return backend_main(fasta_in, email, workingdir, clients_folder, outputfile,
                         report=report, job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", "--input_file", action="store", type=str, dest="FASTA_IN",
        help="name of fasta file containing seq.", required=True)
    parser.add_argument("--email", action="store", type=str, dest="EMAIL",
        help="unused by the local backend; kept for CLI/signature parity", default="")
    parser.add_argument("--dir", "--working_dir", action="store", type=str, dest="WORKINGDIR",
        help="working directory for output", required=True)
    parser.add_argument("--clients-folder", action="store", type=str, dest="CLIENTS_FOLDER",
        help="(unused; retained for CLI compatibility)", default="scripts/")
    parser.add_argument("--output", action="store", type=str, dest="OUTPUT",
        help="output file name")

    args, unknowns = parser.parse_known_args()

    main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.CLIENTS_FOLDER, args.OUTPUT)
