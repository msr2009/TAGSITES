"""ebi_rest.py — thin requests-based wrapper for EBI REST job services.

Each service follows the same REST pattern:
  POST {base_url}/run             → jobId string
  GET  {base_url}/status/{jobId} → QUEUED | RUNNING | FINISHED | ERROR | FAILURE
  GET  {base_url}/result/{jobId}/{resultType} → raw bytes

Use run_job() for the common submit-poll-fetch workflow.
No import-time network calls; safe to import anywhere.
"""

import time
from datetime import datetime
from email.utils import parsedate_to_datetime
import requests

# canonical base URLs for each EBI REST service
NCBIBLAST = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
IPRSCAN5  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
CLUSTALO  = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
GENEWISE  = "https://www.ebi.ac.uk/Tools/services/rest/genewise"

DBFETCH_BASE = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"


RETRYABLE_EXCEPTIONS = (requests.exceptions.Timeout, requests.exceptions.ConnectionError)

# HTTP statuses worth retrying: 429 (rate limited) and the common transient 5xx codes.
# Batch-scale traffic hits these routinely; a plain 4xx (bad request, not found, ...)
# still raises immediately since retrying it can't help.
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


def _retry_after_seconds(resp):
    """Parse a response's Retry-After header (delta-seconds or HTTP-date); None if absent/unparsable."""
    value = resp.headers.get("Retry-After")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        try:
            dt = parsedate_to_datetime(value)
            return max(0.0, (dt - datetime.now(dt.tzinfo)).total_seconds())
        except Exception:
            return None


def _request_with_retries(method, url, retries=3, retry_wait=5, **kwargs):
    """Call requests.<method>(url, **kwargs), retrying on transient timeout/connection errors
    and on 429/5xx responses (honoring Retry-After when the server sends one).

    EBI's REST endpoints occasionally hang past the read timeout, and under batch-scale
    load return 429/503; retrying the same idempotent GET/POST a few times clears most
    of these transparently. A non-retryable status (e.g. a plain 4xx) still raises
    immediately via raise_for_status().
    """
    last_exc = None
    for attempt in range(retries + 1):
        try:
            resp = getattr(requests, method)(url, **kwargs)
            if resp.status_code in RETRYABLE_STATUS and attempt < retries:
                wait = _retry_after_seconds(resp)
                time.sleep(wait if wait is not None else retry_wait)
                continue
            resp.raise_for_status()
            return resp
        except RETRYABLE_EXCEPTIONS as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(retry_wait)
    raise last_exc


def submit(base_url, params):
    """POST params to {base_url}/run; return the jobId string."""
    resp = _request_with_retries("post", f"{base_url}/run", data=params, timeout=60)
    return resp.text.strip()


def get_status(base_url, job_id):
    """GET current status string for a submitted job."""
    resp = _request_with_retries("get", f"{base_url}/status/{job_id}", timeout=30)
    return resp.text.strip()


def fetch_result(base_url, job_id, result_type):
    """GET one result type for a finished job; return raw bytes."""
    resp = _request_with_retries("get", f"{base_url}/result/{job_id}/{result_type}", timeout=120)
    return resp.content


def dbfetch(db, accession, fmt="fasta", style="raw"):
    """Fetch a record from EBI dbfetch; return raw bytes."""
    url = f"{DBFETCH_BASE}/{db}/{accession}/{fmt}/{style}"
    resp = _request_with_retries("get", url, timeout=60)
    return resp.content


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


def run_job(base_url, params, poll_cb=None, poll_interval=5, backoff=1.5, max_interval=60,
            max_wait=7200):
    """Submit a job, poll until FINISHED, return the jobId.

    poll_cb(job_id, status_str) is called after each status check when provided —
    use it to capture the jobId and surface intermediate status to callers.
    Raises RuntimeError if the job ends in ERROR or FAILURE, or if it is still
    QUEUED/RUNNING after max_wait seconds of polling (default 2h) — a wall-clock
    safety valve so a job stuck at the EBI end doesn't poll forever; pass
    max_wait=None to disable it and poll indefinitely as before.
    """
    job_id = submit(base_url, params)
    if poll_cb:
        poll_cb(job_id, "QUEUED")

    interval = poll_interval
    elapsed = 0.0
    while True:
        time.sleep(interval)
        elapsed += interval
        current = get_status(base_url, job_id)
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
    status = get_status(base_url, job_id)
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
