"""ebi_rest.py — thin requests-based wrapper for EBI REST job services.

Each service follows the same REST pattern:
  POST {base_url}/run             → jobId string
  GET  {base_url}/status/{jobId} → QUEUED | RUNNING | FINISHED | ERROR | FAILURE
  GET  {base_url}/result/{jobId}/{resultType} → raw bytes

Use run_job() for the common submit-poll-fetch workflow.
No import-time network calls; safe to import anywhere.
"""

import time
import requests

# canonical base URLs for each EBI REST service
NCBIBLAST = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
IPRSCAN5  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
CLUSTALO  = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
GENEWISE  = "https://www.ebi.ac.uk/Tools/services/rest/genewise"

DBFETCH_BASE = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"


def submit(base_url, params):
    """POST params to {base_url}/run; return the jobId string."""
    resp = requests.post(f"{base_url}/run", data=params, timeout=60)
    resp.raise_for_status()
    return resp.text.strip()


def get_status(base_url, job_id):
    """GET current status string for a submitted job."""
    resp = requests.get(f"{base_url}/status/{job_id}", timeout=30)
    resp.raise_for_status()
    return resp.text.strip()


def fetch_result(base_url, job_id, result_type):
    """GET one result type for a finished job; return raw bytes."""
    resp = requests.get(f"{base_url}/result/{job_id}/{result_type}", timeout=120)
    resp.raise_for_status()
    return resp.content


def dbfetch(db, accession, fmt="fasta", style="raw"):
    """Fetch a record from EBI dbfetch; return raw bytes."""
    url = f"{DBFETCH_BASE}/{db}/{accession}/{fmt}/{style}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
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


def run_job(base_url, params, poll_cb=None, poll_interval=5, backoff=1.5, max_interval=60):
    """Submit a job, poll until FINISHED, return the jobId.

    poll_cb(job_id, status_str) is called after each status check when provided —
    use it to capture the jobId and surface intermediate status to callers.
    Raises RuntimeError if the job ends in ERROR or FAILURE.
    """
    job_id = submit(base_url, params)
    if poll_cb:
        poll_cb(job_id, "QUEUED")

    interval = poll_interval
    while True:
        time.sleep(interval)
        current = get_status(base_url, job_id)
        if poll_cb:
            poll_cb(job_id, current)
        if current == "FINISHED":
            return job_id
        if current in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"EBI job {job_id} ended with status: {current}")
        # exponential backoff up to max_interval
        interval = min(interval * backoff, max_interval)
