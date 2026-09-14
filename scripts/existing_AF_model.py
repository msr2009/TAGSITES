"""
existing_AF_model.py

CLI entry point for AlphaFold structure lookup. The actual lookup is
delegated to a backend resolved by scripts/providers.py —
scripts/structure_remote.py (UniProt checksum → BLAST fallback → AFDB
download) by default, or a bulk-precomputed lookup when batch.config.json
sets backends.structure to a local mode. With no config file present this
always resolves to "remote", so the Shiny app's behavior here is unchanged.

Parameters:
    - fasta_in (str): path to FASTA file, OR a raw UniProt accession
    - email (str): EBI-registered email
    - workingdir (str): directory for output files
    - name (str): run name (prefix for output files)
    - taxid (str|int): taxonomy ID to constrain search
    - evalue (float): E-value threshold for BLAST hit (fallback only)
    - percentid (float): %ID threshold for BLAST hit (fallback only)
    - clients_folder (str): (unused; retained for CLI compatibility)
    - report: optional progress reporter (see scripts/progress.py)

Returns:
    - path to the AF2 FASTA (.fa), or 1 if not found

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py; backend-split 2026
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import resolve

# re-exported for backward compatibility: these pure parsers used to live here
# directly and are covered by tests/test_existing_af_model.py against this
# module's name, even though the logic that calls them now lives in
# structure_remote.py (the default "structure" backend).
from uniprot_api import checksum_lookup_response_to_accession
from structure_remote import blast_tsv_line_to_afdb_hit, is_afdb_not_found


def search_AFDB(fasta_in, email, workingdir, name, taxid, evalue, percentid,
                clients_folder, report=None):
    """Run structure lookup via the configured backend (UniProt/AFDB by
    default); same signature/return value as before the backend split.

    Returns the path to the downloaded AF2 FASTA, or 1 if not found.
    """
    backend_main = resolve("structure")
    return backend_main(fasta_in, email, workingdir, name, taxid, evalue, percentid,
                         clients_folder, report=report)


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", "--input_file", action="store", type=str, dest="FASTA_IN",
        help="path to FASTA file, or a UniProt accession", required=True)
    parser.add_argument("--email", action="store", type=str, dest="EMAIL",
        help="email address, required by EBI job submission.", required=True)
    parser.add_argument("--dir", "--working_dir", action="store", type=str, dest="WORKINGDIR",
        help="working directory for output", required=True)
    parser.add_argument("--name", "--run_name", action="store", type=str, dest="NAME",
        help="name for output", required=True)
    parser.add_argument("--taxid", action="store", type=str, dest="TAXID",
        help="UniProt taxid to limit BLAST search to", default="1")
    parser.add_argument("--evalue", action="store", type=float, dest="EVALUE",
        help="E-value threshold for BLAST hit (1e-100)", default=1e-100)
    parser.add_argument("--percent_id", action="store", type=float, dest="PERCENTID",
        help="Identity threshold for BLAST hit (99)", default=99)
    parser.add_argument("--clients_folder", action="store", type=str, dest="CLIENTS_FOLDER",
        help="(unused; retained for CLI compatibility)", default="./scripts/")

    args, unknowns = parser.parse_known_args()

    result = search_AFDB(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME,
                          args.TAXID, args.EVALUE, args.PERCENTID, args.CLIENTS_FOLDER)
    # search_AFDB returns 1 (not 0/None) when no AFDB model was found; propagate
    # that as the process exit code so callers (run_tag_sites_from_json.py) can
    # tell "no match" apart from "successfully found and downloaded a model".
    sys.exit(1 if result == 1 else 0)
