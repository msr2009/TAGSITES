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


def clean_sequence(seq):
    """Strip stop-codon asterisks and whitespace so checksum and BLAST see the bare protein."""
    return "".join(seq.replace("*", "").split())


def fetch_afdb_pdb(accession, reporter):
    """Download the AFDB PDB text for an accession, or None if unavailable."""
    try:
        pdb_bytes = ebi_rest.dbfetch("afdb", accession, "pdb", "raw")
        pdb_text = pdb_bytes.decode("utf-8", errors="replace")
    except Exception as e:
        _report(reporter, f"dbfetch failed for {accession}: {e}",
                stage="afdb_fetch", level="error")
        return None

    # the EBI dbfetch returns an error message as plain text when not found
    if is_afdb_not_found(pdb_text):
        _report(reporter, f"PDB not found in AFDB for {accession}",
                stage="afdb_fetch", level="warning")
        return None
    return pdb_text


def main(fasta_in, email, workingdir, name, taxid, evalue, percentid,
         clients_folder, report=None):
    """Checksum lookup → BLAST fallback → AFDB download → write PDB + FASTA.

    Returns the path to the downloaded AF2 FASTA, or 1 if not found.
    """
    reporter = resolve_reporter(report)
    candidates = []  # (accession, percent_id, evalue), best BLAST hit first

    outfile_prefix = f"{workingdir}/{name}.AF"

    # if fasta_in looks like a UniProt accession, skip all searching
    if uniprot_accession_regex(fasta_in) is not None:
        candidates.append((fasta_in, 100.0, 1e-200))
    else:
        seq = clean_sequence(get_sequence(fasta_in))

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
            candidates.append((acc, 100.0, 1e-200))
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

            # keep every hit: identical sequences tie, and the first may lack an AFDB model
            for line in lines[1:]:
                if line.strip():
                    hit = blast_tsv_line_to_afdb_hit(line)
                    candidates.append((hit["accession"], hit["percent_id"], hit["evalue"]))

    # try candidates in order; the first passing thresholds with an AFDB model wins
    for match_accession, match_id, match_eval in candidates:
        if match_accession == "" or match_id < percentid or match_eval > evalue:
            continue
        pdb_text = fetch_afdb_pdb(match_accession, reporter)
        if pdb_text is None:
            continue

        pdb_path = f"{outfile_prefix}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_text)
        _report(reporter, f"saved {match_accession} PDB to {pdb_path}", stage="afdb_save")

        fasta_path = f"{outfile_prefix}.fa"
        save_fasta(name, get_sequence(pdb_path), fasta_path)
        _report(reporter, f"saved FASTA from PDB to {fasta_path}", stage="afdb_save")
        return fasta_path

    _report(reporter, f"no hit with an AFDB model at E<={evalue} and %ID>={percentid}",
            stage="afdb_nohit", level="warning")
    return 1
