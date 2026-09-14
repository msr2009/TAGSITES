"""
structure_remote.py

UniProt/AFDB backend for structure lookup — moved out of existing_AF_model.py
verbatim as part of the local/remote backend split (see scripts/providers.py).
This is the default backend; a future structure_bulk.py will look up a
bulk-downloaded AlphaFold DB proteome tarball instead of the per-protein
checksum-lookup / BLAST-fallback / dbfetch-download sequence below.

1) Tries an exact CRC64 checksum lookup against UniProt (fast, no queue).
2) Falls back to BLASTing against UniProt if no exact match is found.
3) Downloads the PDB file from the AFDB for the matched accession.

Matt Rich, 4/2024 / updated 2026 — EBI REST calls via ebi_rest.py; backend-split 2026
"""

import sys
from pathlib import Path

from site_selection_util import get_sequence, uniprot_accession_regex, save_fasta
from uniprot_api import checksum_lookup as _uniprot_checksum_lookup

sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, timed_poll_adapter


def blast_tsv_line_to_afdb_hit(line):
    """Parse one raw NCBIBLAST-TSV data line into an AFDB candidate hit dict."""
    l = line.strip().split("\t")
    return {"accession": l[2], "percent_id": float(l[7]), "evalue": float(l[9])}


def is_afdb_not_found(pdb_text):
    """EBI dbfetch returns a plain-text 'ERROR ...' message when an AFDB accession isn't found."""
    return "ERROR" in pdb_text[:50]


def main(fasta_in, email, workingdir, name, taxid, evalue, percentid,
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
                "alignments": 5,
                "scores":     5,
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

            hit = blast_tsv_line_to_afdb_hit(lines[1])
            match_id        = hit["percent_id"]
            match_eval      = hit["evalue"]
            match_accession = hit["accession"]

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
        if is_afdb_not_found(pdb_text):
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
