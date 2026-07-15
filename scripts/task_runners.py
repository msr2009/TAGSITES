"""task_runners.py — in-process wrappers that map a run-JSON args dict to each
analysis script's main() function.

Used by the Shiny progress module to run tasks in-process (no subprocess)
instead of launching them as shell commands.  The CLI interface of every
underlying script is unchanged; only the internal call path differs.

Each runner signature: runner(args_dict, report=None, job_id_cb=None, resume_job_ids=None)
  args_dict      — the "args" sub-dict from the run JSON for one task
  report         — optional progress reporter, forwarded to the script's main()
  job_id_cb(index, jid)  — optional; tags a freshly-submitted EBI job ID by its
                           sequential position among the task's EBI submissions
  resume_job_ids — optional; the job ID list persisted from a previous attempt,
                   used to reattach to still-pending or already-finished EBI jobs
                   instead of resubmitting. Runners that don't call EBI ignore both.
"""

import os
import sys
from pathlib import Path

# add scripts/ to sys.path so sibling imports work when called from the app
_SCRIPTS = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS))


def _str(v, default=""):
    """Coerce a value to str; return default for empty/None."""
    if v is None or v == "":
        return default
    return str(v)


def _float(v, default=0.0):
    """Coerce to float with fallback."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _int(v, default=0):
    """Coerce to int with fallback."""
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def _bool(v, default=False):
    """Coerce to bool; handles strings like 'True'/'False'."""
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.lower() not in ("false", "0", "")
    return bool(v) if v is not None else default


# ── AFDB pre-step ────────────────────────────────────────────────────────────

def afdb_presearch(tasks, report=None):
    """Run AFDB lookup before tasks start and swap the reference .fa for all tasks.

    Finds the plddt task whose pdb arg is empty, runs the AFDB search, and if a
    model is found: saves the user's .fa aside, puts the AF-derived .fa in its
    place (so all tasks share the same reference sequence), and updates the plddt
    task's pdb arg in-place.  Other task args need no update because their
    input_file paths already point to the canonical .fa path.

    Returns True if the swap happened, False if no AFDB lookup was needed or the
    lookup returned no hit.
    """
    import existing_AF_model

    plddt_task = next(
        (t for t in tasks if t["analysis"] == "plddt" and not t["args"].get("pdb")),
        None,
    )
    if plddt_task is None:
        return False

    args = plddt_task["args"]
    wd   = _str(args.get("working_dir"))
    name = _str(args.get("run_name"))

    result = existing_AF_model.search_AFDB(
        fasta_in      = _str(args.get("input_file") or args.get("fasta")),
        email         = _str(args.get("email")),
        workingdir    = wd,
        name          = name,
        taxid         = _str(args.get("taxid") or args.get("species_taxid") or "1"),
        evalue        = _float(args.get("evalue"), 1e-10),
        percentid     = _float(args.get("percent_id"), 99),
        clients_folder= str(_SCRIPTS) + "/",
        report        = report,
    )

    if result == 1:
        return False  # no hit; plddt will fail on its own, others use user's .fa

    # swap .fa files so all tasks use the AF sequence as reference
    orig_fa = os.path.join(wd, name + ".fa")
    user_fa = os.path.join(wd, "user_" + name + ".fa")
    af_fa   = os.path.join(wd, name + ".AF.fa")
    af_pdb  = os.path.join(wd, name + ".AF.pdb")

    if os.path.exists(orig_fa):
        os.rename(orig_fa, user_fa)
    if os.path.exists(af_fa):
        os.rename(af_fa, orig_fa)

    # point the plddt task at the downloaded PDB
    plddt_task["args"]["pdb"] = af_pdb
    return True


# ── analysis runners ─────────────────────────────────────────────────────────

def run_blast(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run blast_orthologs.main() from a JSON args dict."""
    import blast_orthologs
    return blast_orthologs.main(
        fasta_in       = _str(args.get("input_file") or args.get("fasta")),
        email          = _str(args.get("email")),
        workingdir     = _str(args.get("working_dir")),
        name           = _str(args.get("run_name")),
        output         = _str(args.get("output")),
        n              = _int(args.get("max_hits"), 100),
        evalue         = _float(args.get("evalue"), 1e-10),
        db             = _str(args.get("db"), "uniprotkb"),
        length_percent = _float(args.get("length"), 0),
        align_full_seqs= _bool(args.get("align_full_seqs"), True),
        taxid          = _str(args.get("taxid"), "1"),
        taxid_file     = _str(args.get("taxid_file")) or None,
        clients_folder = str(_SCRIPTS) + "/",
        exclude_paralogs = _bool(args.get("exclude_paralogs"), False),
        report         = report,
        job_id_cb      = job_id_cb,
        resume_job_ids = resume_job_ids,
    )


def run_plddt(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run extract_from_pdb.main(). AFDB lookup is handled by afdb_presearch() before tasks start."""
    import extract_from_pdb

    pdb_path = _str(args.get("pdb"))
    output   = _str(args.get("output"))

    if not pdb_path:
        raise RuntimeError("No PDB path set for pLDDT — AFDB lookup may have failed or no PDB was provided.")

    extract_from_pdb.main(pdb_path, output)


def run_modifications(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run regex_sites.main() from a JSON args dict."""
    import regex_sites
    fasta_in  = _str(args.get("input_file") or args.get("fasta"))
    sites_file = _str(args.get("sites_file"))
    output    = _str(args.get("output"))
    with open(output, "w") as fout:
        regex_sites.main(fasta_in, sites_file, fout)


def run_domains(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run call_interpro.main() from a JSON args dict."""
    import call_interpro
    return call_interpro.main(
        fasta_in      = _str(args.get("input_file") or args.get("fasta")),
        email         = _str(args.get("email")),
        workingdir    = _str(args.get("working_dir")),
        clients_folder= str(_SCRIPTS) + "/",
        outputfile    = _str(args.get("output")),
        report        = report,
        job_id_cb     = job_id_cb,
        resume_job_ids = resume_job_ids,
    )


def run_scores(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run calculate_protein_scores.main() from a JSON args dict."""
    import calculate_protein_scores
    fasta_in   = _str(args.get("input_file") or args.get("fasta"))
    scores_file = _str(args.get("scores_file"))
    output     = _str(args.get("output"))
    window     = _int(args.get("window"), 1)
    with open(output, "w") as fout:
        calculate_protein_scores.main(fasta_in, fout, window, scores_file)


def run_reagents(args, report=None, job_id_cb=None, resume_job_ids=None):
    """Run design_tag_reagents.main() from a JSON args dict.

    Runs the Genewise pre-step if needed (empty 'genewise' arg).
    """
    import design_tag_reagents
    import run_genewise

    genewise_path = _str(args.get("genewise"))
    genomic_fa    = _str(args.get("genomic_fasta"))

    protein_fa = _str(args.get("input_file") or args.get("fasta"))

    if not genewise_path and genomic_fa:
        # run both-strand Genewise and pick the winner
        email       = _str(args.get("email"))
        working_dir = _str(args.get("working_dir"))
        run_name    = _str(args.get("run_name"))
        outprefix   = os.path.join(working_dir, run_name + "_genewise")
        genewise_result = run_genewise.main(protein_fa, genomic_fa, email, outprefix, report=report,
                                            job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)
        # only genewise's own {"ebi_status": "pending"|"expired"} sentinel means
        # "not done yet" — its success return (select_orientation's score dict)
        # is also a dict, so it must not be mistaken for that sentinel here.
        if isinstance(genewise_result, dict) and "ebi_status" in genewise_result:
            return genewise_result
        genewise_path = outprefix + ".genewise.out.txt"
        args = dict(args, genewise=genewise_path,
                    genomic_fasta=outprefix + ".genewise_genomic.fa")

    # Compute protein_length from the protein FASTA for coverage validation
    protein_length = None
    if protein_fa and os.path.exists(protein_fa):
        from Bio import SeqIO as _SeqIO
        recs = list(_SeqIO.parse(protein_fa, 'fasta'))
        if recs:
            protein_length = len(recs[0].seq)

    design_tag_reagents.main(
        genewise            = genewise_path,
        genomic_fasta       = _str(args.get("genomic_fasta")),
        output              = _str(args.get("output")),
        protein_length      = protein_length,
        n_guides            = _int(args.get("n_guides"), 5),
        arm_length          = _int(args.get("arm_length"), 1000),
        pam                 = _str(args.get("PAM"), "NGG"),
        guide_length        = _int(args.get("guide_length"), 20),
        cut_offset          = _int(args.get("cut_offset"), 3),
        insert_sequence     = _str(args.get("insert_sequence")),
        internal_threshold  = _int(args.get("internal_threshold"), 500),
        primer_opt_tm       = _float(args.get("primer_opt_tm"), 60.0),
        product_opt_size    = _int(args.get("product_opt_size"), 200),
        flank_min           = _int(args.get("flank_min"), 50),
        flank_max           = _int(args.get("flank_max"), 150),
        report              = report,
    )


# map analysis type → runner function (used by progress_server and task_runners)
TASK_RUNNERS = {
    "blast":         run_blast,
    "plddt":         run_plddt,
    "modifications": run_modifications,
    "domains":       run_domains,
    "scores":        run_scores,
    "reagents":      run_reagents,
}
