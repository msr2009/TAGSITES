"""ebi_rest.py — thin requests-based wrapper for EBI REST job services.

Each service follows the same REST pattern:
  POST {base_url}/run             → jobId string
  GET  {base_url}/status/{jobId} → QUEUED | RUNNING | FINISHED | ERROR | FAILURE
  GET  {base_url}/result/{jobId}/{resultType} → raw bytes

Use run_job() for the common submit-poll-fetch workflow.
No import-time network calls; safe to import anywhere.
"""

import sys
import time
from pathlib import Path

# requests is no longer called directly here (http_retry owns that), but it stays
# imported as the patch seam the tests use: monkeypatching ebi_rest.requests.get
# patches the one shared module object http_retry also calls through.
import requests  # noqa: F401

sys.path.insert(0, str(Path(__file__).parent))
import ebi_health
import http_retry

# canonical base URLs for each EBI REST service
NCBIBLAST = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
IPRSCAN5  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
CLUSTALO  = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
GENEWISE  = "https://www.ebi.ac.uk/Tools/services/rest/genewise"

DBFETCH_BASE = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"
ENA_FASTA_BASE = "https://www.ebi.ac.uk/ena/browser/api/fasta"


# The retry machinery moved to http_retry.py so ensembl_rest.py can share it.
# Re-exported under the original names: every call site in this file and any
# external reference keeps working, and the default deadline=None means the
# behavior here is byte-for-byte what it was.
RETRYABLE_EXCEPTIONS = http_retry.RETRYABLE_EXCEPTIONS
RETRYABLE_STATUS = http_retry.RETRYABLE_STATUS
_retry_after_seconds = http_retry.retry_after_seconds
_request_with_retries = http_retry.request_with_retries


def submit(base_url, params):
    """POST params to {base_url}/run; return the jobId string."""
    try:
        resp = _request_with_retries("post", f"{base_url}/run", data=params, timeout=60)
    except requests.exceptions.HTTPError as exc:
        # surface EBI's own reason ("temporarily blocked ...", a bad parameter) instead
        # of a bare "400 Client Error"; a fair-use block also feeds the app's EBI banner
        body = getattr(exc.response, "text", "") or ""
        desc = ebi_health.error_description(body)
        if ebi_health.is_block_message(body):
            # EBI gives a placeholder email (example.org, ...) the same refusal: that is
            # the submitter's address, not a fair-use block on the server
            if ebi_health._usable_email(params.get("email")) is None:
                desc += (" (EBI sends this for placeholder addresses such as example.org;"
                         " enter a real email address)")
            else:
                ebi_health.note("blocked", desc, email=params.get("email"))
        if desc:
            raise requests.exceptions.HTTPError(
                f"{exc} — EBI says: {desc}", response=exc.response) from exc
        raise
    return resp.text.strip()


def get_status(base_url, job_id):
    """GET current status string for a submitted job."""
    resp = _request_with_retries("get", f"{base_url}/status/{job_id}", timeout=30)
    return resp.text.strip()


# EBI's status endpoint asks its own SLURM REST API about the job; when that internal API
# rate-limits EBI (HTTP 429, "Maximum 10 requests/second"), EBI wraps it in a plain HTTP 400
# with a SlurmException body. Seen 2026-11-14 on ncbiblast and iprscan5 alike, persisting for
# a lone manual request, so it is EBI-wide load rather than our polling rate. The job itself
# is unaffected, so this means "status unknown right now", not "bad request".
_BACKEND_BUSY_MARKERS = ("SlurmException", "rate_limit_exceeded")


def _backend_busy(exc):
    """True when an HTTPError is EBI's wrapped internal rate-limit error, not a real 400."""
    resp = getattr(exc, "response", None)
    if resp is None or resp.status_code != 400:
        return False
    try:
        body = resp.text or ""
    except Exception:
        return False
    return any(m in body for m in _BACKEND_BUSY_MARKERS)


def _status_or_none(base_url, job_id):
    """Job status string, or None while EBI's backend is rate-limited (try again later)."""
    try:
        return get_status(base_url, job_id)
    except requests.exceptions.HTTPError as exc:
        if _backend_busy(exc):
            ebi_health.note("status_busy")
            print(f"EBI status for {job_id} unavailable (EBI backend rate-limited); "
                  f"will poll again", file=sys.stderr)
            return None
        raise


def fetch_result(base_url, job_id, result_type):
    """GET one result type for a finished job; return raw bytes."""
    resp = _request_with_retries("get", f"{base_url}/result/{job_id}/{result_type}", timeout=120)
    return resp.content


def dbfetch(db, accession, fmt="fasta", style="raw"):
    """Fetch a record from EBI dbfetch; return raw bytes."""
    url = f"{DBFETCH_BASE}/{db}/{accession}/{fmt}/{style}"
    resp = _request_with_retries("get", url, timeout=60)
    return resp.content


def ena_subsequence(accession, start, end):
    """Fetch bases [start, end] (1-based inclusive) of an ENA record as a plain string.

    Uses the ENA Browser API's ?range= form. dbfetch itself has no subsequence
    support (an "ACC:start-end" id returns "No entries found"), and the Browser
    API's ?start=&end= parameters are silently ignored — they return the WHOLE
    record, which for a chromosome is tens of megabytes. Only ?range= actually
    slices, so do not "simplify" this to the other spelling.

    Returns "" when the range cannot be fetched, so callers can treat a failure as
    "unknown" rather than as an answer.
    """
    url = f"{ENA_FASTA_BASE}/{accession}?range={int(start)}-{int(end)}"
    try:
        resp = _request_with_retries("get", url, timeout=60)
    except Exception:
        return ""
    lines = resp.text.splitlines()
    return "".join(ln.strip() for ln in lines if ln and not ln.startswith(">")).upper()


def fmt_exp(value):
    """Format an E-value float as the string EBI expects (e.g. 1e-10, not 1e-10 with zero-padded exp).

    EBI BLAST rejects Python's default float formatting when the exponent is zero-padded
    (e.g. '1e-05'). This function strips the zero-padding from the exponent.
    """
    s = f"{float(value):e}"  # e.g. '1.000000e-05'
    # convert to compact form: strip mantissa trailing zeros, remove + in exponent
    mantissa, exp = s.split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    exp_sign = "-" if exp.startswith("-") else ""
    exp_digits = exp.lstrip("+-").lstrip("0") or "0"
    if mantissa == "1":
        return f"1e{exp_sign}{exp_digits}"
    return f"{mantissa}e{exp_sign}{exp_digits}"


def run_job(base_url, params, poll_cb=None, poll_interval=5, backoff=1.5, max_interval=20,
            max_wait=7200):
    """Submit a job, poll until FINISHED, return the jobId.

    poll_cb(job_id, status_str) is called after each status check when provided —
    use it to capture the jobId and surface intermediate status to callers.
    Raises RuntimeError if the job ends in ERROR or FAILURE, or if it is still
    QUEUED/RUNNING after max_wait seconds of polling (default 2h) — a wall-clock
    safety valve so a job stuck at the EBI end doesn't poll forever; pass
    max_wait=None to disable it and poll indefinitely as before.

    max_interval caps the exponential backoff. It was 60 s, which meant a job
    finishing at 356 s was not noticed until 401 s — pure dead time after EBI had
    already finished. Measured region blastn runs are 173 s (snt-1), 260 s
    (snap-29) and 764 s (col-103), and Genewise/InterProScan/Clustal Omega are
    comparable, so every task was paying up to a minute per job and the reagents
    task pays it four times over (two Genewise + two blastn). At 20 s the lag
    drops to ~17 s for the cost of a few extra cheap status checks.
    """
    job_id = submit(base_url, params)
    if poll_cb:
        poll_cb(job_id, "QUEUED")

    interval = poll_interval
    elapsed = 0.0
    current = "QUEUED"
    while True:
        time.sleep(interval)
        elapsed += interval
        # an unreadable status (EBI backend rate-limited) keeps the last known one;
        # max_wait still bounds how long that can go on
        current = _status_or_none(base_url, job_id) or current
        if poll_cb:
            poll_cb(job_id, current)
        if current == "FINISHED":
            return job_id
        if current in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"EBI job {job_id} ended with status: {current}")
        if max_wait is not None and elapsed >= max_wait:
            raise RuntimeError(
                f"EBI job {job_id} did not finish within {max_wait}s "
                f"(last status: {current})"
            )
        # exponential backoff up to max_interval
        interval = min(interval * backoff, max_interval)


def resume_job(base_url, job_id, result_type):
    """Check a previously-submitted job's status once (no polling loop) and act on it.

    Returns a (state, payload) tuple instead of raising, so callers can branch on
    plain data:
      ("finished", <result bytes>)      job is done; payload is the fetched result
      ("pending",  <raw status string>) still QUEUED/RUNNING; try again later
      ("expired",  <raw status string>) ERROR/FAILURE/NOT_FOUND; the job is gone
    Use this to reattach to a job submitted in an earlier session instead of
    resubmitting it from scratch.
    """
    status = _status_or_none(base_url, job_id)
    # EBI's backend is rate-limited: the job is still there, just unreadable for now
    if status is None:
        return "pending", "RUNNING"
    if status == "FINISHED":
        return "finished", fetch_result(base_url, job_id, result_type)
    if status in {"QUEUED", "RUNNING"}:
        return "pending", status
    return "expired", status


def indexed_job_id_cb(job_id_cb, index):
    """Wrap a job_id_cb(index, jid) callback into a poll_cb(job_id, status)-shaped one.

    Lets a script tag a freshly-submitted job ID with its sequential position
    (0, 1, 2, …) among the EBI calls a task makes, so callers pass it straight
    into run_job()'s poll_cb= without hand-rolling a closure.
    """
    if not job_id_cb:
        return None
    return lambda jid, status: job_id_cb(index, jid)


def combined_poll_cb(*callbacks):
    """Combine several poll_cb(job_id, status) callbacks into one that calls each in turn."""
    def _cb(job_id, status):
        for cb in callbacks:
            if cb:
                cb(job_id, status)
    return _cb
