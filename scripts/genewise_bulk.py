"""
genewise_bulk.py

Local gene-model backend, the bulk-mode FALLBACK: gene_model.ensure_model first tries UCSC's
UniProt alignment from the local index (ucsc_local.py / ucsc_gene_model.py, same rules as the
remote path), and calls this only when that has no exact model. No Genewise runs: the CDS
exons are read from the GFF3 index (scripts/genome_regions.py), the matching region is cut
from the local genome FASTA, and the transcript is chosen by translating every annotated CDS
and matching the input protein exactly (genome_regions.find_transcripts_by_sequence). With no
exact match, every transcript at the gene's locus is compared with cds_check and the closest
is used only if it passes the same cluster/gap rules as the UCSC and NCBI models (evaluate()).
No genomic region needs to be supplied.

Output files keep Genewise's names and format ({outprefix}.genewise.out.txt,
{outprefix}.genewise_genomic.fa, {outprefix}.genewise_orientation.txt) so the reagent code
downstream is unchanged. The alignment score column is a fixed sentinel
(parse_genewise.GROUND_TRUTH_SCORE), since annotation-derived exons are ground truth. Orientation
is always "+": the extracted DNA is already reverse-complemented to the coding strand.

Indexes built before the translation table existed fall back to the UniProt route: accession or
checksum in local_store -> WormBase locus name -> same-prefix GFF3 transcripts, each translated
and compared (genome_regions.resolve_transcript_for_accession). Raises LookupError when no
transcript qualifies.
"""

import sys
from pathlib import Path

from site_selection_util import get_sequence, read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from local_store import lookup_by_accession, lookup_by_crc64, open_index as open_proteins_index
from genome_regions import (
    _load_genome, find_transcripts_by_sequence, get_transcript_region,
    resolve_transcript_for_accession, open_index as open_genome_index,
)
from uniprot_api import parse_accession
from progress import report as _report, resolve_reporter


def _resolve_accession(protein_fasta, proteins_conn):
    """Resolve protein_fasta to a UniProt accession: the string itself if it
    already looks like one, else its FASTA/PDB record id, else a CRC64
    checksum lookup against local_store. Returns None if nothing resolves.
    """
    acc = parse_accession(protein_fasta)
    if acc:
        return acc

    if protein_fasta.endswith(".pdb"):
        seq = get_sequence(protein_fasta)
        name = None
    else:
        name, seq = read_fasta(protein_fasta)
        acc = parse_accession(name)
        if acc:
            return acc

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


def _candidate(transcript_id, region, check, equivalent):
    """A GFF3 transcript in the candidate shape of ucsc_gene_model.evaluate."""
    return {
        "source": "gff3", "name": transcript_id, "label": f"GFF3 transcript {transcript_id}",
        "model": region, "check": check, "fasta_source": None,
        "fasta_id": f"{transcript_id}_genomic", "gff_label": "genome_regions.py (GFF3-derived)",
        "region": {"source": "gff3", "transcript_id": transcript_id,
                   "equivalent_transcripts": equivalent},
    }


def evaluate(protein_fasta, locus=None, report=None):
    """GFF3 transcripts as candidate models for the protein, best first.

    An exact translation (CRC64 index, else the WormBase locus name) returns at once. Without
    one, every transcript at the gene's locus is translated and compared with cds_check, the
    same rules as the UCSC and NCBI candidates. locus is (chrom, start, stop, strand), 1-based,
    in UCSC or GFF3 chromosome names (from gene_model's UCSC candidates); the WormBase locus
    name's transcripts are added when local_store knows the protein. Raises LookupError when
    no transcript can be found.
    """
    import cds_check
    from genome_regions import transcripts_overlapping

    reporter = resolve_reporter(report)
    protein_row = _local_store_row(protein_fasta)
    seq = _input_sequence(protein_fasta, protein_row).rstrip("*")
    label = protein_row["accession"] if protein_row else Path(protein_fasta).stem
    _report(reporter, f"Resolving GFF3 transcript for {label}…", stage="genewise_bulk")

    flank = _flank_bp()
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
        if matches:
            # transcripts with the same CDS (UTR/splice variants) are equivalent; take the longest
            region = get_transcript_region(matches[0], conn=genome_conn, flank=flank)
            check = cds_check.compare_to_query(
                cds_check.translate_model(region["cds_df"], region["dna"]), seq)
            return [_candidate(matches[0], region, check, matches)]

        # no exact match: every transcript of the gene, judged like the other sources
        ids = []
        if locus:
            chrom, start, stop, strand = locus
            genome = _load_genome()
            if chrom not in genome:
                import ucsc_local

                chrom = ucsc_local.local_chrom(chrom)
            ids += transcripts_overlapping(chrom, start, stop, strand, conn=genome_conn)
        if protein_row and protein_row.get("wormbase_gene"):
            ids += [r[0] for r in genome_conn.execute(
                "SELECT transcript_id FROM transcript_spans WHERE transcript_id LIKE ?",
                (protein_row["wormbase_gene"] + ".%",)).fetchall()]
        ids = list(dict.fromkeys(ids))
        if not ids:
            raise LookupError(f"no GFF3 transcript translates to {label}'s sequence and no "
                              "gene locus is known to search for the closest one")
        candidates = []
        for transcript_id in ids:
            try:
                region = get_transcript_region(transcript_id, conn=genome_conn, flank=flank)
            except (ValueError, KeyError):
                continue
            check = cds_check.compare_to_query(
                cds_check.translate_model(region["cds_df"], region["dna"]), seq)
            candidates.append(_candidate(transcript_id, region, check, []))
    finally:
        genome_conn.close()
    candidates.sort(key=lambda c: (cds_check.rank(c["check"]),
                                   -(c["model"]["stop"] - c["model"]["start"])))
    # the closest few are enough to show what was compared
    for c in candidates[:cds_check.LOGGED_CANDIDATES]:
        _report(reporter, cds_check.summary_line(c["check"], c["label"]), stage="genewise_bulk")
    if len(candidates) > cds_check.LOGGED_CANDIDATES:
        _report(reporter, f"... {len(candidates) - cds_check.LOGGED_CANDIDATES} further "
                          "candidate(s) differ more", stage="genewise_bulk")
    return candidates


def main(protein_fasta, genomic_fasta, email, outprefix, report=None,
         job_id_cb=None, resume_job_ids=None, locus=None):
    """Write Genewise-format files for the best GFF3 transcript of the protein.

    genomic_fasta/email/job_id_cb/resume_job_ids are accepted but unused (the region comes from
    the local genome, no network job is submitted) so providers.resolve("genewise") can call
    either backend identically. batch_run.genomic_flank_bp of genome is added around the gene
    (clipped at chromosome ends) and recorded, with the transcript coordinates and any
    equivalent transcripts, in {outprefix}.region.json. An exact translation wins; otherwise
    the closest transcript is used when cds_check only warns. Raises LookupError, with the
    located differences of the closest transcript, when none qualifies.
    """
    import cds_check
    from gene_model import write_model

    reporter = resolve_reporter(report)
    candidates = evaluate(protein_fasta, locus=locus, report=report)
    best = candidates[0] if candidates else None
    if best is None or best["check"]["status"] == "fail":
        raise LookupError(
            "no GFF3 transcript qualifies as the gene model"
            + ("\n" + cds_check.describe(best["check"], best["label"]) if best else "")
            + "\n" + cds_check.MODEL_MISMATCH_HELP)
    if best["check"]["status"] == "warn":
        _report(reporter, cds_check.describe(best["check"], best["label"]),
                stage="genewise_bulk", level="warning")
    winner_out, winner_fa, _ = write_model(outprefix, best)
    region = best["model"]
    _report(reporter, f"Resolved to transcript {best['name']} ({region['chrom']}:"
                       f"{region['start']}-{region['stop']}, strand {region['strand']})",
            stage="genewise_bulk")
    _report(reporter, f"winner → {winner_out}  genomic → {winner_fa}", stage="genewise_select")
    return {
        "orientation": "+",
        "out_txt": winner_out,
        "genomic_fa": winner_fa,
        "transcript_id": best["name"],
        "warning": None,
    }
