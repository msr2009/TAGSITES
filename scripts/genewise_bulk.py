"""
genewise_bulk.py

Local backend for CDS-exon inference: reads the exon structure straight out
of the WormBase GFF3 annotation (via scripts/genome_regions.py) for any
isoform present in it, instead of submitting two EBI Genewise jobs (forward +
reverse-complement) to infer that structure genewise_remote.py's way.

Matching a protein to a specific GFF3 transcript is not a plain string match:
local_store.py's wormbase_gene (e.g. "C10C5.1g") is the isoform-lettered
locus name, but WormBase can have several numbered transcript versions under
that same name (splice/UTR variants with no further distinguishing hint) —
verified in a 300-accession sample, ~8% of name-prefix matches were
ambiguous and ~1% had none at all. genome_regions.resolve_transcript_for_accession()
resolves this by translating each same-prefix candidate's CDS and comparing
to the UniProt protein sequence already in local_store — 296/300 (98.7%)
resolved this way in that sample; the remainder (in that sample: two
non-nuclear-code mitochondrial genes plus two others) raise LookupError so
a batch driver can fall back to genewise_remote.py for them, per
batch.config.json's documented "config with fallback" model (the actual
per-task fallback dispatch is Phase C work, not yet implemented — this
module only ever does the local lookup or raises).

Writes the same {outprefix}.genewise.out.txt / {outprefix}.genewise_genomic.fa
/ {outprefix}.genewise_orientation.txt files genewise_remote.py's winning
orientation produces, in the same GFF-embedded-score format
parse_genewise.py's parser expects — with the alignment score column set to
a fixed sentinel (see parse_genewise.GROUND_TRUTH_SCORE) since there's no genewise bitscore
for annotation-derived exons, only exact ground truth. Orientation is always
"+" because the extracted genomic FASTA is already reverse-complemented to
the coding strand by genome_regions.extract_sequence() when needed, unlike
genewise_remote.py's RC submission which reflects the *original* input
region's orientation.
"""

import sys
from pathlib import Path

from site_selection_util import get_sequence, save_fasta, uniprot_accession_regex, read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from local_store import lookup_by_accession, lookup_by_crc64, open_index as open_proteins_index
from genome_regions import resolve_transcript_for_accession, open_index as open_genome_index
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


def main(protein_fasta, genomic_fasta, email, outprefix, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Look up the GFF3-derived CDS exon structure for protein_fasta's
    protein and write the same {outprefix}.genewise.out.txt /
    {outprefix}.genewise_genomic.fa / {outprefix}.genewise_orientation.txt
    files genewise_remote.main() produces. genomic_fasta/email/job_id_cb/
    resume_job_ids are accepted but unused (no genomic region needs to be
    supplied — it's read from the local genome FASTA via the resolved
    transcript's span — and no network job is submitted); kept so
    providers.resolve("genewise") can call either backend identically.

    Raises LookupError if the protein's accession can't be resolved, or no
    GFF3 transcript resolves to it — see this module's docstring for the
    intended caller behaviour (fall back to genewise_remote.py).
    """
    reporter = resolve_reporter(report)

    proteins_conn = open_proteins_index()
    try:
        accession = _resolve_accession(protein_fasta, proteins_conn)
        if accession is None:
            raise LookupError(
                f"could not resolve a UniProt accession for {protein_fasta!r} "
                "in the local index"
            )
        protein_row = lookup_by_accession(accession, conn=proteins_conn)
        if protein_row is None or not protein_row.get("wormbase_gene"):
            raise LookupError(
                f"accession {accession!r} has no WormBase cross-reference in the local index"
            )
    finally:
        proteins_conn.close()

    _report(reporter, f"Resolving GFF3 transcript for {accession}…", stage="genewise_bulk")

    genome_conn = open_genome_index()
    try:
        transcript_id = resolve_transcript_for_accession(
            accession, protein_row["wormbase_gene"], protein_row["sequence"], conn=genome_conn,
        )
        if transcript_id is None:
            raise LookupError(
                f"no GFF3 transcript under locus {protein_row['wormbase_gene']!r} "
                f"translates to accession {accession}'s sequence"
            )

        from genome_regions import get_transcript_region
        region = get_transcript_region(transcript_id, conn=genome_conn)
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
