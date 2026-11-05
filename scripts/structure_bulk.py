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

When the tarball lacks an accession (it holds one model per gene, about 19.7k for C. elegans,
while AFDB also models the other UniProt entries), the model is fetched from AFDB through
structure_remote.fetch_afdb_pdb and cached next to the tarball's files, so the next lookup is
local. A failed fetch is the same "no model" outcome as before.

UniProt often keeps one sequence under two entries, and AFDB may model only the one the
reference proteome does not list (unc-25: Q7K7A3 in the proteome, the model under G5EDB7).
So after the protein's own accession, every accession with the identical sequence in
local_store's sequence_aliases index is tried, local files before AFDB fetches; the model's
accession is then recorded in {name}.uniprot_accession.txt.

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
from local_store import afdb_model_path, aliases_for_sequence, lookup_by_crc64
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


def _fetch_missing_model(accession, pdb_gz_path, reporter):
    """Fetch a model absent from the tarball and cache it as gzip; return True on success."""
    from structure_remote import fetch_afdb_pdb

    # AFDB models canonical entries only, so an isoform accession cannot be fetched
    if "-" in accession:
        return False
    pdb_text = fetch_afdb_pdb(accession, reporter)
    if pdb_text is None:
        return False
    # write through a temp file so an interrupted run never leaves a truncated cache entry
    tmp = Path(str(pdb_gz_path) + ".tmp")
    with gzip.open(tmp, "wt") as f_out:
        f_out.write(pdb_text)
    os.replace(tmp, pdb_gz_path)
    _report(reporter, f"fetched AFDB structure for {accession} (not in the local tarball)",
            stage="afdb_bulk")
    return True


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
    accessions = [accession] if accession else []
    # other entries with this exact sequence: AFDB often models only one of them
    if expected is not None:
        accessions += [a for a in _aliases(expected) if a not in accessions]
    if not accessions:
        _report(reporter,
                "Could not resolve a UniProt accession for this sequence in "
                "the local index — no AFDB structure available from the bulk "
                "backend.", stage="afdb_bulk", level="warning")
        return 1

    # local models first (no network), then AFDB fetches in the same accession order
    for acc, fetch in [(a, False) for a in accessions] + [(a, True) for a in accessions]:
        pdb_gz_path = afdb_model_path(acc)
        if fetch:
            # a file present now was already tried in the local pass
            if pdb_gz_path.exists() or not _fetch_missing_model(acc, pdb_gz_path, reporter):
                continue
        elif not pdb_gz_path.exists():
            continue
        fasta_path = _use_model(pdb_gz_path, acc, workingdir, name, expected, reporter)
        if fasta_path is None:
            continue
        if acc != accessions[0]:
            _report(reporter, f"AFDB has no model for {accessions[0]}; used {acc}, an entry "
                              "with the identical sequence", stage="afdb_bulk")
            # record the model's accession beside the run (as structure_remote does on a match)
            from structure_remote import _record_accession
            _record_accession(workingdir, name, acc)
        return fasta_path

    _report(reporter, f"no AFDB structure for {', '.join(accessions)}",
            stage="afdb_bulk", level="warning")
    return 1


def _aliases(seq):
    """Accessions in the local alias index with exactly this sequence ([] without the index)."""
    try:
        return aliases_for_sequence(seq)
    except Exception:
        return []


def _use_model(pdb_gz_path, accession, workingdir, name, expected, reporter):
    """Write {name}.AF.pdb / .fa from a model whose sequence equals expected; else None."""
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
        return None

    fasta_path = f"{outfile_prefix}.fa"
    save_fasta(name, model_seq, fasta_path)
    _report(reporter, f"saved FASTA from PDB to {fasta_path}", stage="afdb_save")
    return fasta_path
