"""
run_genewise.py

CLI entry point for Genewise-based CDS-exon inference. The actual work is
delegated to a backend resolved by scripts/providers.py —
scripts/genewise_remote.py (both-orientation EBI Genewise jobs) by default, or
a WormBase-GFF3-derived lookup when batch.config.json sets backends.genewise
to a local mode. With no config file present this always resolves to
"remote", so the Shiny app's behavior here is unchanged.

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

Matt Rich, 2025 / backend-split 2026
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import resolve

# re-exported for backward compatibility: these pure, network-free helpers used
# to live here directly and are covered by tests/test_run_genewise.py against
# this module's name, even though the logic that calls them now lives in
# genewise_remote.py (the default "genewise" backend).
from genewise_remote import select_orientation, write_rc_fasta


def main(protein_fasta, genomic_fasta, email, outprefix, report=None,
         job_id_cb=None, resume_job_ids=None):
    """Run Genewise-based CDS inference via the configured backend (both-
    orientation EBI Genewise by default); same signature/return value as
    before the backend split.
    """
    backend_main = resolve("genewise")
    return backend_main(protein_fasta, genomic_fasta, email, outprefix,
                         report=report, job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)


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
