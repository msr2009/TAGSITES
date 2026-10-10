"""
structure_bulk.py

Local backend for structure lookup: reads a PDB directly out of the
AFDB proteome tarball extracted by reference_data.fetch_afdb() into
out_dir/afdb_pdbs/AF-<accession>-F1-model_v6.pdb.gz, instead of a UniProt
checksum lookup + BLAST fallback + EBI dbfetch download.

Accession resolution mirrors domains_bulk.py: fasta_in's FASTA header if it
already looks like a UniProt accession, else a raw accession string, else a
CRC64 checksum lookup against local_store.py's proteins table (the same key
uniprot_api.checksum_lookup() uses over the network).

AFDB only has predictions for canonical UniProt entries — a non-canonical
isoform's accession won't have its own AF-*.pdb.gz member. That's a known,
documented limitation (see the plan's Phase B), not something this backend
tries to work around; it returns 1 (not found) for those, exactly like
structure_remote.py's "no BLAST hit good enough" branch.

A model is used only when its sequence EQUALS the input sequence. task_runners.afdb_presearch
swaps the model's sequence in as the reference for every task, so a model built on an older
UniProt release would silently move modifications, scores and conservation onto the wrong
protein (28 of 28,626 proteome entries, three of them a different gene model, <45%
identical). A mismatch returns 1, the same "no model" outcome, and the input sequence stays.
Selenocysteine is compared as cysteine because the PDB carries it as CYS.
"""

import gzip
import os
import sys
from pathlib import Path

from site_selection_util import get_sequence, uniprot_accession_regex, save_fasta

sys.path.insert(0, str(Path(__file__).parent))
from local_store import lookup_by_crc64, _reference_dir
from progress import report as _report, resolve_reporter


def _resolve_accession(fasta_in):
    """Resolve fasta_in to a UniProt accession: fasta_in itself if it already
    looks like one (mirrors structure_remote.main()'s accession-input shortcut),
    else the FASTA record id if that looks like one, else a CRC64 lookup
    against the local proteins index. Returns None if nothing resolves.
    """
    if uniprot_accession_regex(fasta_in):
        return fasta_in

    from site_selection_util import read_fasta
    name, seq = read_fasta(fasta_in)
    if uniprot_accession_regex(name):
        return name

    from Bio.SeqUtils.CheckSum import crc64
    checksum = crc64(str(seq)).replace("CRC-", "")
    matches = lookup_by_crc64(checksum)
    return matches[0]["accession"] if matches else None


def main(fasta_in, email, workingdir, name, taxid, evalue, percentid,
         clients_folder, report=None):
    """Look up a local AFDB structure for fasta_in's protein and write the
    same {name}.AF.pdb / {name}.AF.fa pair structure_remote.py produces.
    email/taxid/evalue/percentid are accepted but unused (no BLAST fallback
    is run locally); kept so providers.resolve("structure") can call either
    backend identically. Returns the path to the AF2 FASTA, or 1 if not found.
    """
    reporter = resolve_reporter(report)
    accession = _resolve_accession(fasta_in)
    # a FASTA path carries the sequence the model has to match; a bare accession does not
    expected = None
    if not uniprot_accession_regex(fasta_in):
        from site_selection_util import read_fasta
        expected = str(read_fasta(fasta_in)[1])
    if accession is None:
        _report(reporter,
                "Could not resolve a UniProt accession for this sequence in "
                "the local index — no AFDB structure available from the bulk "
                "backend.", stage="afdb_bulk", level="warning")
        return 1

    ref_dir, _ = _reference_dir()
    pdb_gz_path = ref_dir / "afdb_pdbs" / f"AF-{accession}-F1-model_v6.pdb.gz"
    if not pdb_gz_path.exists():
        _report(reporter, f"no local AFDB structure for {accession}",
                stage="afdb_bulk", level="warning")
        return 1

    outfile_prefix = f"{workingdir}/{name}.AF"
    pdb_path = f"{outfile_prefix}.pdb"
    with gzip.open(pdb_gz_path, "rt") as f_in, open(pdb_path, "w") as f_out:
        f_out.write(f_in.read())
    _report(reporter, f"copied local AFDB structure for {accession} -> {pdb_path}",
            stage="afdb_save")

    model_seq = str(get_sequence(pdb_path))
    if expected is not None and model_seq.replace("U", "C") != expected.replace("U", "C"):
        _report(reporter, f"AFDB model for {accession} is a different sequence version "
                          f"({len(model_seq)} aa vs {len(expected)} aa); not used",
                stage="afdb_bulk", level="warning")
        os.remove(pdb_path)
        return 1

    fasta_path = f"{outfile_prefix}.fa"
    save_fasta(name, model_seq, fasta_path)
    _report(reporter, f"saved FASTA from PDB to {fasta_path}", stage="afdb_save")
    return fasta_path
