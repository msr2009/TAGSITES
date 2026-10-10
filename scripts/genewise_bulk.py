"""
genewise_bulk.py

Local gene-model backend. No Genewise runs: the CDS exons are read from the GFF3 index
(scripts/genome_regions.py), the matching region is cut from the local genome FASTA, and the
transcript is chosen by translating every annotated CDS and matching the input protein exactly
(genome_regions.find_transcripts_by_sequence). No genomic region needs to be supplied.

Output files keep Genewise's names and format ({outprefix}.genewise.out.txt,
{outprefix}.genewise_genomic.fa, {outprefix}.genewise_orientation.txt) so the reagent code
downstream is unchanged. The alignment score column is a fixed sentinel
(parse_genewise.GROUND_TRUTH_SCORE), since annotation-derived exons are ground truth. Orientation
is always "+": the extracted DNA is already reverse-complemented to the coding strand.

Indexes built before the translation table existed fall back to the UniProt route: accession or
checksum in local_store -> WormBase locus name -> same-prefix GFF3 transcripts, each translated
and compared (genome_regions.resolve_transcript_for_accession). Raises LookupError when nothing
translates to the input protein.
"""

import json
import sys
from pathlib import Path

from site_selection_util import get_sequence, save_fasta, uniprot_accession_regex, read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from local_store import lookup_by_accession, lookup_by_crc64, open_index as open_proteins_index
from genome_regions import (
    find_transcripts_by_sequence, get_transcript_region, resolve_transcript_for_accession,
    open_index as open_genome_index,
)
from parse_genewise import write_genewise_gff
from progress import report as _report, resolve_reporter


def _resolve_accession(protein_fasta, proteins_conn):
    """Resolve protein_fasta to a UniProt accession: the string itself if it
    already looks like one, else its FASTA/PDB record id, else a CRC64
    checksum lookup against local_store. Returns None if nothing resolves.
    """
    if uniprot_accession_regex(protein_fasta):
        return protein_fasta

    if protein_fasta.endswith(".pdb"):
        seq = get_sequence(protein_fasta)
        name = None
    else:
        name, seq = read_fasta(protein_fasta)
        if uniprot_accession_regex(name):
            return name

    from Bio.SeqUtils.CheckSum import crc64
    checksum = crc64(str(seq)).replace("CRC-", "")
    matches = lookup_by_crc64(checksum, conn=proteins_conn)
    return matches[0]["accession"] if matches else None


def _flank_bp():
    """batch_run.genomic_flank_bp (default 2000): bases of genome kept around the gene.

    Homology arms are up to 1 kb either side of an insertion site, so a gene without flank
    would make every N- and C-terminal site fail the short-arm check.
    """
    from providers import _load_config

    return int(_load_config().get("batch_run", {}).get("genomic_flank_bp", 2000) or 0)


def _local_store_row(protein_fasta):
    """local_store row for the input (accession, record id or checksum); None if absent or unbuilt."""
    try:
        conn = open_proteins_index()
    except FileNotFoundError:
        return None
    try:
        accession = _resolve_accession(protein_fasta, conn)
        return lookup_by_accession(accession, conn=conn) if accession else None
    finally:
        conn.close()


def _input_sequence(protein_fasta, protein_row):
    """The protein sequence to match: the input file's own, else the local_store entry's."""
    if Path(protein_fasta).exists():
        return str(get_sequence(protein_fasta) if protein_fasta.endswith(".pdb")
                   else read_fasta(protein_fasta)[1])
    if protein_row:
        return protein_row["sequence"]
    raise LookupError(f"{protein_fasta!r} is neither a file nor a known UniProt accession")


def main(protein_fasta, genomic_fasta, email, outprefix, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Write Genewise-format files for the GFF3 transcript whose CDS translates to the protein.

    genomic_fasta/email/job_id_cb/resume_job_ids are accepted but unused (the region comes from
    the local genome, no network job is submitted) so providers.resolve("genewise") can call
    either backend identically. batch_run.genomic_flank_bp of genome is added around the gene
    (clipped at chromosome ends) and recorded, with the transcript coordinates and any
    equivalent transcripts, in {outprefix}.region.json. Raises LookupError if no transcript matches.
    """
    reporter = resolve_reporter(report)

    protein_row = _local_store_row(protein_fasta)
    seq = _input_sequence(protein_fasta, protein_row)
    label = protein_row["accession"] if protein_row else Path(protein_fasta).stem

    _report(reporter, f"Resolving GFF3 transcript for {label}…", stage="genewise_bulk")

    genome_conn = open_genome_index()
    try:
        # exact translation match against every annotated CDS
        matches = find_transcripts_by_sequence(seq, conn=genome_conn)
        # older indexes lack the translation table: use the UniProt locus name instead
        if not matches and protein_row and protein_row.get("wormbase_gene"):
            transcript_id = resolve_transcript_for_accession(
                protein_row["accession"], protein_row["wormbase_gene"], seq, conn=genome_conn,
            )
            matches = [transcript_id] if transcript_id else []
        if not matches:
            raise LookupError(
                f"no GFF3 transcript translates to {label}'s sequence; supply a GenBank gene "
                "model, or set backends.genewise to \"remote\" and upload a genomic region"
            )
        # transcripts with the same CDS (UTR/splice variants) are equivalent; take the longest
        transcript_id = matches[0]
        region = get_transcript_region(transcript_id, conn=genome_conn, flank=_flank_bp())
    finally:
        genome_conn.close()

    _report(reporter, f"Resolved to transcript {transcript_id} ({region['chrom']}:"
                       f"{region['start']}-{region['stop']}, strand {region['strand']})",
            stage="genewise_bulk")

    winner_fa = f"{outprefix}.genewise_genomic.fa"
    save_fasta(f"{transcript_id}_genomic", region["dna"], winner_fa)

    winner_out = f"{outprefix}.genewise.out.txt"
    write_genewise_gff(winner_out, region["cds_df"], region["chrom"],
                       "genome_regions.py (GFF3-derived)")

    with open(f"{outprefix}.region.json", "w") as fh:
        json.dump({"transcript_id": transcript_id, "chrom": region["chrom"],
                   "start": region["start"], "stop": region["stop"],
                   "strand": region["strand"], "flank_5p": region["flank_5p"],
                   "flank_3p": region["flank_3p"],
                   "equivalent_transcripts": matches}, fh)

    orient_file = f"{outprefix}.genewise_orientation.txt"
    with open(orient_file, "w") as fh:
        # always "+": the extracted DNA is already oriented to the coding
        # strand (see module docstring), unlike genewise_remote.py's RC
        # submission which reflects the *original* input region's orientation
        fh.write("+\n")

    _report(reporter, f"winner → {winner_out}  genomic → {winner_fa}  orientation → {orient_file}",
            stage="genewise_select")

    return {
        "orientation": "+",
        "out_txt": winner_out,
        "genomic_fa": winner_fa,
        "transcript_id": transcript_id,
        "warning": None,
    }
