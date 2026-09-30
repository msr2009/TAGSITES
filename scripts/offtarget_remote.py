"""
offtarget_remote.py

EBI blastn backend for the off-target / primer-specificity screen.

Submits at most two jobs per run and writes a sidecar JSON that both the pipeline
and the Shiny UI read:

  Screen A  the whole genomic region as one query -> duplicated segments elsewhere
            in the same species (primer co-amplification, guides in repeats)
  Screen B  every spacer+PAM concatenated into one query, separated by N-runs
            -> scattered guide near-matches with no regional homology

Why one big query instead of one per spacer: EBI's ncbiblast takes a single query
per job, so per-spacer searching would mean dozens of jobs. A multi-kb region query
also seeds robustly at the nucleotide word size, whereas a bare 20 nt query is
marginal there.

The sidecar exists because the Shiny UI re-designs genotyping primers on demand
from the stored homology arms; caching the parsed HSPs lets the UI screen those
primers locally instead of submitting another BLAST.

Named *_remote.py to match the scripts/providers.py <analysis>_<mode> convention,
so a local blastn backend can be dropped in later. Note conservation_local.py
already synthesises EBI-shaped JSON from DIAMOND tabular output — the same adapter
trick would work here.

Matt Rich, 2025
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
import offtarget_screen as ots
from progress import report as _report, resolve_reporter, timed_poll_adapter

# EBI submission indices for the reagents task. Genewise already owns 0 and 1
# (see genewise_remote.py), so these must not collide.
JOB_INDEX_REGION = 2
JOB_INDEX_SPACER = 3

_RESULT_TYPE = "json"       # must match between resume_job() and fetch_result()


def database_for_taxid(taxid, cfg):
    """Pick the ENA division for a taxid; returns (database, is_fallback)."""
    dbs = cfg["databases"]
    key = str(taxid).strip()
    if key in dbs:
        return dbs[key], False
    return dbs.get("default", "em_all"), True


def _blast_params(email, sequence, database, taxid, evalue, wordsize, alignments):
    """Build the ncbiblast POST body for a nucleotide search."""
    params = {
        "email":      email,
        "program":    "blastn",
        "stype":      "dna",
        "sequence":   str(sequence),
        "database":   database,
        "outformat":  _RESULT_TYPE,
        # alignments and scores must agree or EBI returns an empty hit list
        "alignments": alignments,
        "scores":     alignments,
        "exp":        ebi_rest.fmt_exp(evalue),
        "wordsize":   wordsize,
        "task":       "blastn",     # traditional blastn; megablast is too strict for paralogs
    }
    # taxid 1 (root) is the project-wide "no scope" sentinel; see structure_remote.py
    if str(taxid).strip() not in ("", "1", "1.0", "None"):
        params["taxids"] = str(taxid).strip()
    return params


def _submit(sequence, email, database, taxid, evalue, wordsize, alignments,
            reporter, stage, job_id_cb, job_index, resume_id, wordsize_fallback=None):
    """Run one blastn job (or resume it); returns (state, payload)."""
    if resume_id:
        _report(reporter, "Checking previously-submitted blastn job…", stage=stage)
        return ebi_rest.resume_job(ebi_rest.NCBIBLAST, resume_id, _RESULT_TYPE)

    params = _blast_params(email, sequence, database, taxid, evalue, wordsize, alignments)
    poll_cb = ebi_rest.combined_poll_cb(
        ebi_rest.indexed_job_id_cb(job_id_cb, job_index),
        timed_poll_adapter(reporter, stage=stage),
    )
    try:
        job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, params, poll_cb=poll_cb)
    except Exception as e:
        # A word size below EBI's nucleotide default may be rejected outright; retry
        # once at the documented value rather than losing the whole screen
        if wordsize_fallback is None or int(wordsize_fallback) == int(wordsize):
            raise
        _report(reporter, "wordsize {} rejected ({}); retrying at {}".format(
            wordsize, e, wordsize_fallback), stage=stage, level="warning")
        params["wordsize"] = wordsize_fallback
        job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, params, poll_cb=poll_cb)
    return "finished", ebi_rest.fetch_result(ebi_rest.NCBIBLAST, job_id, _RESULT_TYPE)


def run_screens(region_seq, spacers, email, taxid, sidecar_path, exons=None,
                pam="NGG", cfg=None, report=None, job_id_cb=None, resume_job_ids=None,
                run_spacer_screen=True):
    """Run both screens and write the sidecar; returns the sidecar dict.

    Returns {"ebi_status": ...} instead when a job is still queued or has expired,
    matching the sentinel that progress_server.py uses to resume a task.
    """
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    resume = list(resume_job_ids or [])

    def _resume_at(i):
        return resume[i] if len(resume) > i and resume[i] else None

    database, is_fallback = database_for_taxid(taxid, cfg)
    if is_fallback:
        _report(reporter, "No ENA division mapped for taxid {}; using {} — transcript "
                          "divisions are in scope and hits will be noisier".format(taxid, database),
                stage="offtarget", level="warning")

    blast, pf, sq = cfg["blast"], cfg["post_filter"], cfg["spacer_query"]

    # ── Screen A: region homology ────────────────────────────────────────────
    _report(reporter, "Submitting region off-target blastn ({} bp vs {}, taxid {})".format(
        len(region_seq), database, taxid or "unscoped"), stage="offtarget_region")
    state, payload = _submit(
        region_seq, email, database, taxid, blast["evalue_region"],
        blast["wordsize_region"], blast["alignments"], reporter,
        "offtarget_region", job_id_cb, JOB_INDEX_REGION, _resume_at(JOB_INDEX_REGION))
    if state in ("pending", "expired"):
        return {"ebi_status": state, "detail": payload}

    region_hsps = ots.parse_blast_json(payload)
    ots.classify_hsps(region_hsps, len(region_seq), exons or [], cfg)
    region_hsps, n_dropped = ots.apply_post_filter(region_hsps, cfg)
    dups = ots.duplicates(region_hsps)
    identical = ots.identical_segments(region_hsps, cfg)
    n_self = sum(1 for h in region_hsps if h["klass"] == "self")
    n_tx = sum(1 for h in region_hsps if h["klass"] == "transcript")
    _report(reporter, "Region screen: {} duplicated segment(s) ({} identical), "
                      "{} self, {} transcript, {} below thresholds".format(
                          len(dups), len(identical), n_self, n_tx, n_dropped),
            stage="offtarget_region")
    if identical:
        _report(reporter, "{} identical duplicated segment(s) >= {} bp — every oligo in "
                          "those windows will act on both copies".format(
                              len(identical), cfg["identical"]["identical_min_len"]),
                stage="offtarget_region", level="warning")

    # ── Screen B: concatenated spacers ──────────────────────────────────────
    spacer_hits = {}
    spacer_list = list(spacers or [])[: int(sq["max_spacers"])]
    block_len = 0
    if run_spacer_screen and spacer_list:
        query, block_len = ots.build_spacer_query(spacer_list, sq["separator_len"])
        _report(reporter, "Submitting spacer off-target blastn ({} spacers, {} bp)".format(
            len(spacer_list), len(query)), stage="offtarget_spacer")
        state, payload = _submit(
            query, email, database, taxid, blast["evalue_spacer"],
            blast["wordsize_spacer"], blast["alignments"], reporter,
            "offtarget_spacer", job_id_cb, JOB_INDEX_SPACER, _resume_at(JOB_INDEX_SPACER),
            wordsize_fallback=blast["wordsize_fallback"])
        if state in ("pending", "expired"):
            return {"ebi_status": state, "detail": payload}
        raw = ots.parse_blast_json(payload)
        spacer_hits = ots.screen_spacer_hits(
            raw, spacer_list, block_len, sq["separator_len"], pam, cfg)
        _report(reporter, "Spacer screen: {} of {} spacers have a near-match".format(
            len(spacer_hits), len(spacer_list)), stage="offtarget_spacer")

    sidecar = {
        "_meta": {
            "taxid":         str(taxid),
            "database":      database,
            "database_is_fallback": is_fallback,
            "region_len":    len(region_seq),
            "pam":           pam,
            "blast":         blast,
            "post_filter":   pf,
            "n_below_threshold": n_dropped,
            "spacer_block_len":  block_len,
            "spacer_separator_len": sq["separator_len"],
        },
        "duplicates": dups,
        "identical":  identical,
        # kept for visibility but never counted; a paralog whose only ENA record is
        # an mRNA would land here, which is the known blind spot of this screen
        "excluded":   [h for h in region_hsps if h["klass"] in ("self", "transcript")],
        "spacers":    spacer_list,
        "spacer_hits": {str(k): v for k, v in spacer_hits.items()},
    }
    Path(sidecar_path).write_text(json.dumps(sidecar, indent=1))
    return sidecar


def load_sidecar(path):
    """Read a previously written sidecar, or None when absent/unreadable."""
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None
