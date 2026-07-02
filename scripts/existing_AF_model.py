"""
existing_AF2_model.py

Check whether the AlphaFold database has a predicted structure for a protein.

1) Tries an exact CRC64 checksum lookup against UniProt (fast, no queue).
2) Falls back to BLASTing against UniProt if no exact match is found.
3) Downloads the PDB file from the AFDB for the matched accession.

Parameters:
    - fasta_in (str): path to FASTA file, OR a raw UniProt accession
    - email (str): EBI-registered email
    - workingdir (str): directory for output files
    - name (str): run name (prefix for output files)
    - taxid (str|int): taxonomy ID to constrain search
    - evalue (float): E-value threshold for BLAST hit (fallback only)
    - percentid (float): %ID threshold for BLAST hit (fallback only)
    - clients_folder (str): (unused; retained for CLI compatibility)
    - progress_cb: optional callback(job_id, status) for EBI poll status

Returns:
    - path to the AF2 FASTA (.fa), or 1 if not found

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py
"""

import sys
import requests
from pathlib import Path
from Bio.SeqUtils.CheckSum import crc64

from site_selection_util import get_sequence, uniprot_accession_regex, save_fasta

sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, timed_poll_adapter


def _uniprot_checksum_lookup(seq, taxid=None):
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
    results = resp.json().get("results", [])
    if results:
        return results[0]["primaryAccession"]
    return None


def search_AFDB(fasta_in, email, workingdir, name, taxid, evalue, percentid,
                clients_folder, report=None):
    """Checksum lookup → BLAST fallback → AFDB download → write PDB + FASTA.

    Returns the path to the downloaded AF2 FASTA, or 1 if not found.
    """
    reporter = resolve_reporter(report)
    match_accession = ""
    match_eval = 1e-200
    match_id = 100.0

    outfile_prefix = f"{workingdir}/{name}.AF"

    # if fasta_in looks like a UniProt accession, skip all searching
    if uniprot_accession_regex(fasta_in) is not None:
        match_accession = fasta_in
    else:
        seq = get_sequence(fasta_in)

        # fast path: exact sequence match via CRC64 checksum (no queue, milliseconds)
        _report(reporter, "Checking UniProt for exact sequence match…", stage="afdb_checksum")
        try:
            acc = _uniprot_checksum_lookup(seq, taxid)
        except Exception as e:
            _report(reporter, f"checksum lookup failed ({e}); falling back to BLAST",
                    stage="afdb_checksum", level="warning")
            acc = None

        if acc:
            _report(reporter, f"Exact UniProt match: {acc} — skipping BLAST", stage="afdb_checksum")
            match_accession = acc
        else:
            # fallback: BLAST against UniProt to find closest homolog
            _report(reporter, "No exact match; submitting NCBI BLAST job for AFDB lookup…",
                    stage="afdb_submit")
            params = {
                "email":      email,
                "program":    "blastp",
                "stype":      "protein",
                "sequence":   seq,
                "database":   "uniprotkb",
                "outformat":  "tsv",
                "alignments": 1,
                "scores":     1,
                "exp":        "1e-5",
            }
            if str(taxid) not in ("", "1", "1.0"):
                params["taxids"] = str(taxid)

            blast_job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, params,
                                            poll_cb=timed_poll_adapter(reporter, stage="afdb_blast"))
            tsv_bytes = ebi_rest.fetch_result(ebi_rest.NCBIBLAST, blast_job_id, "tsv")

            tsv_path = f"{outfile_prefix}.ncbiblast.tsv.tsv"
            with open(tsv_path, "wb") as f:
                f.write(tsv_bytes)

            lines = open(tsv_path).readlines()
            if len(lines) < 2:
                _report(reporter, "no BLAST hits found for AFDB lookup", stage="afdb_nohit")
                return 1

            l = lines[1].strip().split("\t")
            match_id        = float(l[7])
            match_eval      = float(l[9])
            match_accession = l[2]

    # try to download AlphaFold model for the matched accession
    if match_id >= percentid and match_eval <= evalue and match_accession != "":
        try:
            pdb_bytes = ebi_rest.dbfetch("afdb", match_accession, "pdb", "raw")
            pdb_text = pdb_bytes.decode("utf-8", errors="replace")
        except Exception as e:
            _report(reporter, f"dbfetch failed for {match_accession}: {e}",
                    stage="afdb_fetch", level="error")
            return 1

        # the EBI dbfetch returns an error message as plain text when not found
        if "ERROR" in pdb_text[:50]:
            _report(reporter, f"PDB not found in AFDB for {match_accession}",
                    stage="afdb_fetch", level="warning")
            return 1

        pdb_path = f"{outfile_prefix}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_text)
        _report(reporter, f"saved {match_accession} PDB to {pdb_path}", stage="afdb_save")

        fasta_path = f"{outfile_prefix}.fa"
        save_fasta(name, get_sequence(pdb_path), fasta_path)
        _report(reporter, f"saved FASTA from PDB to {fasta_path}", stage="afdb_save")
        return fasta_path
    else:
        _report(reporter, f"no BLAST hit better than E={evalue} and %ID={percentid}",
                stage="afdb_nohit", level="warning")
        return 1


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

    search_AFDB(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME,
                args.TAXID, args.EVALUE, args.PERCENTID, args.CLIENTS_FOLDER)
