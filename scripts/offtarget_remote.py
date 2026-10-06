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

# EBI validates `exp` against this exact list of strings — it is NOT free-form. Note
# ebi_rest.fmt_exp() must NOT be used here: it rewrites "1000" as "1e3", which EBI
# rejects ('Value for "exp" is not valid'). Verified against the live service.
EXP_ALLOWED = {
    "1e-200", "1e-100", "1e-50", "1e-10", "1e-5", "1e-4", "1e-3",
    "1e-2", "1e-1", "1.0", "10", "100", "1000",
}


def _exp_value(value):
    """Coerce a configured E-value to one of EBI's accepted literals."""
    s = str(value).strip()
    if s in EXP_ALLOWED:
        return s
    # fall back to the loosest threshold rather than failing the whole screen
    raise ValueError(
        "offtarget.config.json E-value {!r} is not one of EBI's accepted values: "
        "{}".format(s, ", ".join(sorted(EXP_ALLOWED))))


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
        "exp":        _exp_value(evalue),
        "wordsize":   wordsize,
        "task":       "blastn",     # traditional blastn; megablast is too strict for paralogs
    }
    # taxid 1 (root) is the project-wide "no scope" sentinel; see structure_remote.py
    if str(taxid).strip() not in ("", "1", "1.0", "None"):
        params["taxids"] = str(taxid).strip()
    return params


def _submit(sequence, email, database, taxid, evalue, wordsize, alignments,
            reporter, stage, job_id_cb, job_index, resume_id, wordsize_fallback=None,
            poll_max_interval=20):
    """Run one blastn job (or resume it); returns (state, payload)."""
    if resume_id:
        _report(reporter, "Checking previously-submitted blastn job…", stage=stage)
        return ebi_rest.resume_job(ebi_rest.NCBIBLAST, resume_id, _RESULT_TYPE)

    params = _blast_params(email, sequence, database, taxid, evalue, wordsize, alignments)
    poll_cb = ebi_rest.combined_poll_cb(
        ebi_rest.indexed_job_id_cb(job_id_cb, job_index),
        timed_poll_adapter(reporter, stage=stage),
    )
    # run_job backs off to a 60 s poll by default, which adds up to a minute of dead
    # waiting after EBI has already finished. These jobs run for minutes, so a tighter
    # cap costs a handful of extra status checks and removes most of that lag.
    try:
        job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, params, poll_cb=poll_cb,
                                  max_interval=poll_max_interval)
    except Exception as e:
        # A word size below EBI's nucleotide default may be rejected outright; retry
        # once at the documented value rather than losing the whole screen
        if wordsize_fallback is None or int(wordsize_fallback) == int(wordsize):
            raise
        _report(reporter, "wordsize {} rejected ({}); retrying at {}".format(
            wordsize, e, wordsize_fallback), stage=stage, level="warning")
        params["wordsize"] = wordsize_fallback
        job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, params, poll_cb=poll_cb,
                                  max_interval=poll_max_interval)
    return "finished", ebi_rest.fetch_result(ebi_rest.NCBIBLAST, job_id, _RESULT_TYPE)


def _resume_at(resume_job_ids, i):
    """The persisted EBI job id for submission index i, or None."""
    resume = list(resume_job_ids or [])
    return resume[i] if len(resume) > i and resume[i] else None


def run_region_screen(region_seq, email, taxid, exons=None, cfg=None, report=None,
                      job_id_cb=None, resume_job_ids=None):
    """Screen A: the whole region as one query, finding duplicated segments.

    Returns a dict of results, or {"ebi_status": ...} when the job is still queued
    or has expired — the sentinel progress_server.py uses to resume a task.
    """
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    database, is_fallback = database_for_taxid(taxid, cfg)
    if is_fallback:
        _report(reporter, "No ENA division mapped for taxid {}; using {} — transcript "
                          "divisions are in scope and hits will be noisier".format(taxid, database),
                stage="offtarget_region", level="warning")

    blast = cfg["blast"]
    _report(reporter, "Submitting region off-target blastn ({} bp vs {}, taxid {})".format(
        len(region_seq), database, taxid or "unscoped"), stage="offtarget_region")
    state, payload = _submit(
        region_seq, email, database, taxid, blast["evalue_region"],
        blast["wordsize_region"], blast["alignments"], reporter,
        "offtarget_region", job_id_cb, JOB_INDEX_REGION,
        _resume_at(resume_job_ids, JOB_INDEX_REGION),
        poll_max_interval=int(blast.get("poll_max_interval", 20)))
    if state in ("pending", "expired"):
        return {"ebi_status": state, "detail": payload}

    hsps = ots.parse_blast_json(payload)
    ots.classify_hsps(hsps, len(region_seq), exons or [], cfg)
    hsps, n_dropped = ots.apply_post_filter(hsps, cfg)
    dups = ots.duplicates(hsps)
    identical = ots.identical_segments(hsps, cfg)
    n_self = sum(1 for h in hsps if h["klass"] == "self")
    n_tx = sum(1 for h in hsps if h["klass"] == "transcript")
    _report(reporter, "Region screen: {} duplicated segment(s) ({} identical), "
                      "{} self, {} transcript, {} below thresholds".format(
                          len(dups), len(identical), n_self, n_tx, n_dropped),
            stage="offtarget_region")
    if identical:
        _report(reporter, "{} identical duplicated segment(s) >= {} bp — every oligo in "
                          "those windows will act on both copies".format(
                              len(identical), cfg["identical"]["identical_min_len"]),
                stage="offtarget_region", level="warning")

    return {
        "database":     database,
        "is_fallback":  is_fallback,
        "duplicates":   dups,
        "identical":    identical,
        # kept for visibility but never counted; a paralog whose only ENA record is
        # an mRNA would land here, which is the known blind spot of this screen
        "excluded":     [h for h in hsps if h["klass"] in ("self", "transcript")],
        # Screen B needs these to drop each guide's own on-target site
        "self_spans":   ots.self_spans(hsps),
        # records that ARE our locus (own-gene mRNAs, every assembly copy)
        "own_accessions": sorted(ots.own_locus_accessions(hsps)),
        "n_below_threshold": n_dropped,
    }


def run_spacer_screen(spacers, email, taxid, pam="NGG", cfg=None, report=None,
                      job_id_cb=None, resume_job_ids=None, self_spans=None,
                      own_accessions=None):
    """Screen B: all spacers in one query, finding scattered near-matches.

    Call this with only the guides that actually reach the output — a region holds
    hundreds of candidates but the user sees a handful per site, and query length
    drives both EBI queue time and hit volume.
    """
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    sq, blast = cfg["spacer_query"], cfg["blast"]
    spacer_list = list(spacers or [])[: int(sq["max_spacers"])]
    if not spacer_list:
        return {"spacer_hits": {}, "spacers": [], "block_len": 0}
    if len(spacers or []) > len(spacer_list):
        _report(reporter, "Screening the first {} of {} spacers (max_spacers)".format(
            len(spacer_list), len(spacers)), stage="offtarget_spacer", level="warning")

    database, _ = database_for_taxid(taxid, cfg)
    query, block_len = ots.build_spacer_query(spacer_list, sq["separator_len"])
    _report(reporter, "Submitting spacer off-target blastn ({} spacers, {} bp, "
                      "wordsize {})".format(len(spacer_list), len(query),
                                            blast["wordsize_spacer"]),
            stage="offtarget_spacer")
    state, payload = _submit(
        query, email, database, taxid, blast["evalue_spacer"],
        blast["wordsize_spacer"], blast.get("alignments_spacer", blast["alignments"]),
        reporter,
        "offtarget_spacer", job_id_cb, JOB_INDEX_SPACER,
        _resume_at(resume_job_ids, JOB_INDEX_SPACER),
        wordsize_fallback=blast["wordsize_fallback"],
        poll_max_interval=int(blast.get("poll_max_interval", 20)))
    if state in ("pending", "expired"):
        return {"ebi_status": state, "detail": payload}

    raw = ots.parse_blast_json(payload)
    hits, n_self = ots.screen_spacer_hits(
        raw, spacer_list, block_len, sq["separator_len"], pam, cfg,
        exclude_spans=self_spans, exclude_accessions=set(own_accessions or []))
    if not self_spans:
        _report(reporter, "No self-locus coordinates available, so each guide's own "
                          "on-target site is counted as a hit — counts are inflated",
                stage="offtarget_spacer", level="warning")
    _report(reporter, "Spacer screen: {} of {} spacers have an off-target near-match "
                      "({} on-target matches excluded)".format(
                          len(hits), len(spacer_list), n_self),
            stage="offtarget_spacer")

    # blastn rarely extends a 23 nt query over the PAM columns, so settle the
    # undecided sites by fetching their flanks instead of assuming either way
    if cfg.get("pam_check", {}).get("resolve_by_fetch", True):
        hits, pam_stats = resolve_pams(
            hits, pam, reporter, "offtarget_spacer",
            max_fetches=int(cfg.get("pam_check", {}).get("max_fetches", 200)))
        _report(reporter, "Spacer screen after PAM check: {} of {} spacers have a "
                          "PAM-bearing or unresolved off-target".format(
                              len(hits), len(spacer_list)), stage="offtarget_spacer")
    return {"spacer_hits": hits, "spacers": spacer_list, "block_len": block_len}


def resolve_pams(sites_by_key, pam="NGG", reporter=None, stage="offtarget_spacer",
                 max_fetches=200):
    """Fetch each unverified site's flank from ENA and settle whether it has a PAM.

    A near-match with no PAM cannot be cut, so once a PAM is known to be absent the
    site is dropped. Sites whose fetch fails stay unverified rather than being
    assumed either way. Returns (sites_by_key, stats).

    Cheap in practice: sites are already collapsed to distinct loci, so this is a
    handful of small ranged requests (12 for the snt-1 test region), cached by
    accession and coordinates.
    """
    cache = {}
    stats = {"fetched": 0, "confirmed": 0, "rejected": 0, "failed": 0}
    out = {}
    for key, sites in (sites_by_key or {}).items():
        kept = []
        for s in sites:
            if not s.get("pam_unverified"):
                kept.append(s)
                continue
            span = ots.pam_fetch_span(s, len(pam))
            if span is None or stats["fetched"] >= max_fetches:
                kept.append(s)
                continue
            start, end, needs_rc = span
            ck = (s["acc"], start, end)
            if ck not in cache:
                cache[ck] = ebi_rest.ena_subsequence(s["acc"], start, end)
                stats["fetched"] += 1
            bases = cache[ck]
            if needs_rc and bases:
                bases = ots.reverse_complement(bases)
            if not ots.apply_fetched_pam(s, bases, pam):
                stats["failed"] += 1
                kept.append(s)          # still unknown, so still reported
                continue
            if s["pam_ok"]:
                stats["confirmed"] += 1
                kept.append(s)
            else:
                stats["rejected"] += 1   # no PAM -> cannot cut -> not an off-target
        if kept:
            out[key] = kept
    if reporter and stats["fetched"]:
        _report(reporter, "PAM check: fetched {} flank(s) — {} confirmed, {} dropped for "
                          "no PAM, {} still unknown".format(
                              stats["fetched"], stats["confirmed"], stats["rejected"],
                              stats["failed"]),
                stage=stage)
    return out, stats


def write_sidecar(path, region, spacer, taxid, region_len, pam, cfg):
    """Write the parsed hits beside the reagents TSV for the UI to reuse."""
    sq = cfg["spacer_query"]
    sidecar = {
        "_meta": {
            "taxid":                 str(taxid),
            "database":              region.get("database", ""),
            "database_is_fallback":  region.get("is_fallback", False),
            "backend":               region.get("backend", ""),
            "region_len":            region_len,
            "pam":                   pam,
            "blast":                 cfg["blast"],
            "post_filter":           cfg["post_filter"],
            "n_below_threshold":     region.get("n_below_threshold", 0),
            "spacer_block_len":      spacer.get("block_len", 0),
            "spacer_separator_len":  sq["separator_len"],
        },
        "duplicates":  region.get("duplicates", []),
        "identical":   region.get("identical", []),
        "excluded":    region.get("excluded", []),
        # keyed by chromosome for a local screen; lets the UI mark the intended amplicon
        "self_spans":  region.get("self_spans", {}) if region.get("backend") == "local" else {},
        "spacers":     spacer.get("spacers", []),
        "spacer_hits": {str(k): v for k, v in (spacer.get("spacer_hits") or {}).items()},
    }
    Path(path).write_text(json.dumps(sidecar, indent=1))
    return sidecar


def load_sidecar(path):
    """Read a previously written sidecar, or None when absent/unreadable."""
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None
