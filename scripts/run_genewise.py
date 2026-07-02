"""
run_genewise.py

Run the EBI Genewise REST client on BOTH orientations of a genomic region FASTA
(forward and reverse-complement) and automatically select the better alignment.

Rationale
---------
Genewise only searches the forward strand of the submitted DNA sequence.  If the
gene of interest is encoded on the minus strand, the forward submission yields a
garbage low-score alignment (score ~1–2 bits, <5% coverage) while the
reverse-complement submission yields the correct result (~1000+ bits).

The score difference is large and unambiguous, so running both orientations
(~10 s each via the EBI REST API) and picking the winner is the cleanest approach.

TODO: integrate this call directly into the pipeline pre-step so that the user
only needs to supply a genomic region FASTA and the orchestrator handles both
orientations automatically, exactly as existing_AF_model.py is called for pLDDT.

Usage
-----
    python scripts/run_genewise.py \\
        --protein_fasta src-1.fa \\
        --genomic_fasta src-1_genomic.fa \\
        --email your@email.com \\
        --outprefix results/src-1

Outputs
-------
    <outprefix>.genewise.out.txt      Genewise output from the winning orientation
    <outprefix>.genewise_genomic.fa   Genomic FASTA in the winning orientation
    <outprefix>.rc.fa                 RC genomic FASTA (always written; for inspection)

The chosen orientation ('+' or 'rc') is printed to stdout and to
<outprefix>.genewise_orientation.txt.

Matt Rich, 2025
"""

import os
import re
import shutil
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Allow importing sibling scripts directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crispr_util import reverse_complement
from parse_genewise import parse_genewise_score, cds_coverage, parse_genewise
from progress import report as _report, resolve_reporter, timed_poll_adapter
import ebi_rest

# Score / coverage thresholds
LOW_SCORE_WARN   = 50.0   # bits – below this a warning is printed
LOW_COVER_WARN   = 0.50   # fraction of protein length covered by CDS


# ── EBI client wrapper ────────────────────────────────────────────────────────

def run_genewise_client(protein_fa, genomic_fa, email, outfile_prefix, report=None):
    """Submit protein+DNA FASTAs to EBI Genewise REST API; write .out.txt and return its path."""
    reporter = resolve_reporter(report)

    with open(protein_fa) as fh:
        asequence = fh.read()
    with open(genomic_fa) as fh:
        bsequence = fh.read()

    params = {
        "email":     email,
        "asequence": asequence,
        "bsequence": bsequence,
        "gff":       "true",  # embed GFF section in output; required by parse_genewise
    }

    _report(reporter, f"Submitting Genewise job for {os.path.basename(protein_fa)}", stage='genewise_run')
    job_id = ebi_rest.run_job(ebi_rest.GENEWISE, params,
                              poll_cb=timed_poll_adapter(reporter, stage='genewise_run'))

    result_bytes = ebi_rest.fetch_result(ebi_rest.GENEWISE, job_id, "out")
    expected = f"{outfile_prefix}.out.txt"
    with open(expected, "wb") as fh:
        fh.write(result_bytes)

    if not os.path.exists(expected):
        raise FileNotFoundError(f"Expected Genewise output not found: {expected}")
    return expected


# ── Orientation selection ─────────────────────────────────────────────────────

def select_orientation(fwd_out_txt, fwd_genomic_fa,
                       rc_out_txt,  rc_genomic_fa,
                       protein_length=None):
    """
    Compare two Genewise results (forward and reverse-complement orientations)
    and return information about the better one.

    This function is intentionally network-free and directly testable: it only
    reads existing .out.txt files and does not call the EBI client.

    Parameters
    ----------
    fwd_out_txt    : str  path to Genewise output for forward orientation
    fwd_genomic_fa : str  path to the forward genomic FASTA
    rc_out_txt     : str  path to Genewise output for RC orientation
    rc_genomic_fa  : str  path to the RC genomic FASTA
    protein_length : int or None  used for coverage calculation (optional)

    Returns
    -------
    dict with keys:
      orientation   '+' or 'rc'
      out_txt       path to the winning .out.txt
      genomic_fa    path to the winning genomic FASTA
      fwd_score     float score for forward run
      rc_score      float score for RC run
      winner_score  float score of the winner
      warning       str or None  populated when the winner looks suspiciously poor
    """
    fwd_score = parse_genewise_score(fwd_out_txt) or 0.0
    rc_score  = parse_genewise_score(rc_out_txt)  or 0.0

    if fwd_score >= rc_score:
        orientation = '+'
        out_txt     = fwd_out_txt
        genomic_fa  = fwd_genomic_fa
        winner_score = fwd_score
    else:
        orientation = 'rc'
        out_txt     = rc_out_txt
        genomic_fa  = rc_genomic_fa
        winner_score = rc_score

    # Coverage check (optional)
    warning = None
    try:
        cds_df = parse_genewise(out_txt)
        plen   = protein_length or max(1, len(cds_df))
        cover  = cds_coverage(cds_df, plen)
    except Exception:
        cover = None

    if winner_score < LOW_SCORE_WARN:
        warning = (
            'WARNING: best Genewise score is very low ({:.2f} bits, orientation={}).\n'
            'This may indicate a failed alignment or a wrong genomic region.\n'
            'Forward score: {:.2f}   RC score: {:.2f}'.format(
                winner_score, orientation, fwd_score, rc_score)
        )
    elif cover is not None and cover < LOW_COVER_WARN:
        warning = (
            'WARNING: CDS coverage of protein is low ({:.0%}, orientation={}).\n'
            'Expected ~100% for a correct alignment.\n'
            'Forward score: {:.2f}   RC score: {:.2f}'.format(
                cover, orientation, fwd_score, rc_score)
        )

    return {
        'orientation':  orientation,
        'out_txt':      out_txt,
        'genomic_fa':   genomic_fa,
        'fwd_score':    fwd_score,
        'rc_score':     rc_score,
        'winner_score': winner_score,
        'warning':      warning,
    }


# ── RC FASTA writer ───────────────────────────────────────────────────────────

def write_rc_fasta(genomic_fa, rc_fa_path):
    """
    Write the reverse-complement of the first record in *genomic_fa* to
    *rc_fa_path*.  Returns the path written.
    """
    records = list(SeqIO.parse(genomic_fa, 'fasta'))
    if not records:
        raise ValueError('No sequences in {}'.format(genomic_fa))
    rec = records[0]
    rc_seq = reverse_complement(str(rec.seq))
    rc_record = SeqRecord(
        Seq(rc_seq),
        id=rec.id + '_rc',
        description='reverse complement of ' + rec.description,
    )
    with open(rc_fa_path, 'w') as fh:
        SeqIO.write(rc_record, fh, 'fasta')
    return rc_fa_path


# ── Main entry point ──────────────────────────────────────────────────────────

def main(protein_fasta, genomic_fasta, email, outprefix, report=None):
    """
    Full both-strand Genewise run + orientation selection.

    Writes:
      <outprefix>.rc.fa                 RC of the genomic FASTA
      <outprefix>.fwd.out.txt           Genewise result for forward orientation
      <outprefix>.rc.out.txt            Genewise result for RC orientation
      <outprefix>.genewise.out.txt      Winner
      <outprefix>.genewise_genomic.fa   Genomic FASTA in winning orientation
      <outprefix>.genewise_orientation.txt  '+' or 'rc'
    """
    reporter = resolve_reporter(report)

    # Handle PDB input: extract protein FASTA if needed
    if protein_fasta.endswith('.pdb'):
        prot_fa = outprefix + '.protein.fa'
        _script_dir = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, _script_dir)
        from site_selection_util import get_sequence, save_fasta
        seq = get_sequence(protein_fasta)
        save_fasta(os.path.basename(outprefix), seq, prot_fa)
        protein_fasta = prot_fa
        _report(reporter, 'extracted protein FASTA from PDB → {}'.format(prot_fa),
                stage='genewise_prep')

    # Write RC genomic FASTA
    rc_fa = outprefix + '.rc.fa'
    write_rc_fasta(genomic_fasta, rc_fa)
    _report(reporter, 'wrote RC FASTA → {}'.format(rc_fa), stage='genewise_prep')

    # Run Genewise on both orientations
    _report(reporter, 'running Genewise (forward orientation)…', stage='genewise_fwd')
    fwd_out = run_genewise_client(protein_fasta, genomic_fasta,
                                  email, outprefix + '.fwd', report=reporter)

    _report(reporter, 'running Genewise (reverse-complement orientation)…', stage='genewise_rc')
    rc_out  = run_genewise_client(protein_fasta, rc_fa,
                                  email, outprefix + '.rc', report=reporter)

    # Select best orientation
    result = select_orientation(fwd_out, genomic_fasta, rc_out, rc_fa)

    _report(reporter,
            'forward score: {:.2f} bits  RC score: {:.2f} bits  selected: {} ({:.2f} bits)'.format(
                result['fwd_score'], result['rc_score'],
                result['orientation'], result['winner_score']),
            stage='genewise_select')

    if result['warning']:
        _report(reporter, result['warning'], stage='genewise_select', level='warning')

    # Write canonical output files
    winner_out = outprefix + '.genewise.out.txt'
    winner_fa  = outprefix + '.genewise_genomic.fa'
    shutil.copy(result['out_txt'], winner_out)
    shutil.copy(result['genomic_fa'], winner_fa)

    orient_file = outprefix + '.genewise_orientation.txt'
    with open(orient_file, 'w') as fh:
        fh.write(result['orientation'] + '\n')

    _report(reporter, 'winner → {}  genomic → {}  orientation → {}'.format(
        winner_out, winner_fa, orient_file), stage='genewise_select')

    return result


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description=(
            'Run Genewise on forward + RC orientations and select the '
            'biologically correct one.  Requires an EBI account e-mail.'
        )
    )
    parser.add_argument('--protein_fasta', required=True,
                        help='Protein sequence FASTA (or PDB; sequence is extracted)')
    parser.add_argument('--genomic_fasta', required=True,
                        help='Genomic region FASTA (orientation unknown)')
    parser.add_argument('--email', required=True,
                        help='E-mail address registered with EBI')
    parser.add_argument('--outprefix', required=True,
                        help='Output file prefix (directory must exist)')
    args = parser.parse_args()

    main(
        protein_fasta = args.protein_fasta,
        genomic_fasta = args.genomic_fasta,
        email         = args.email,
        outprefix     = args.outprefix,
    )
