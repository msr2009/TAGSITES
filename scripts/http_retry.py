"""http_retry.py — shared bounded-retry wrapper around requests, with an
optional wall-clock deadline.

Moved out of ebi_rest.py so the Ensembl client can reuse it. The retry behavior
is unchanged from that original; the deadline support is new and inert unless a
caller passes one (deadline=None reproduces the previous behavior exactly).

Why a deadline and not just requests' timeout=: requests' timeout is a
per-socket connect/read timeout, not a budget for the operation. A server that
dribbles bytes resets the read timer on every chunk and can stream for minutes
without ever tripping timeout=60. Clamping each request's read timeout against
the time remaining turns a chain of per-socket timeouts into a real deadline,
because no single request can outlive the budget and remaining() refuses to
start one that has nothing left.

Matt Rich, 2026
"""

import time
from datetime import datetime
from email.utils import parsedate_to_datetime
import requests

RETRYABLE_EXCEPTIONS = (requests.exceptions.Timeout, requests.exceptions.ConnectionError)

# HTTP statuses worth retrying: 429 (rate limited) and the common transient 5xx codes.
# Batch-scale traffic hits these routinely; a plain 4xx (bad request, not found, ...)
# still raises immediately since retrying it can't help.
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


def deadline_from(total_timeout):
    """Absolute monotonic instant an operation must finish by; None means no cap."""
    if total_timeout is None:
        return None
    return time.monotonic() + float(total_timeout)


def remaining(deadline):
    """Seconds left before deadline, or None when uncapped; raises when spent."""
    if deadline is None:
        return None
    left = deadline - time.monotonic()
    if left <= 0:
        raise TimeoutError("timed out: the time budget for this operation is spent")
    return left


def clamp_timeout(timeout, deadline):
    """Shrink a requests timeout so no single socket wait outlives the deadline.

    Accepts either a scalar read timeout or a (connect, read) tuple, and returns
    the same shape. This is necessary but NOT sufficient on its own: requests'
    read timeout is per-read, so a server dribbling bytes resets it forever and
    never trips it however low it is set. _read_within_deadline is what actually
    bounds that case.
    """
    left = remaining(deadline)
    if left is None or timeout is None:
        return timeout
    if isinstance(timeout, tuple):
        connect, read = timeout
        return (min(connect, left), min(read, left))
    return min(timeout, left)


def _read_within_deadline(resp, deadline, chunk_size=65536):
    """Consume a streamed response body, raising TimeoutError if the budget runs out.

    This is the piece that makes the deadline real. requests' timeout= bounds
    each individual socket read, not the transfer, so a response arriving one
    byte every few seconds never times out at any setting — the exact shape of
    the hang in issue #64. Checking the wall clock between chunks catches it.

    The body is stitched back into the Response so callers keep using .json() /
    .text / .content as if it had never been streamed.
    """
    chunks = []
    for chunk in resp.iter_content(chunk_size):
        chunks.append(chunk)
        remaining(deadline)
    resp._content = b"".join(chunks)
    resp._content_consumed = True
    return resp


def retry_after_seconds(resp):
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


def request_with_retries(method, url, retries=3, retry_wait=5, deadline=None, **kwargs):
    """Call requests.<method>(url, **kwargs), retrying on transient timeout/connection errors
    and on 429/5xx responses (honoring Retry-After when the server sends one).

    EBI's REST endpoints occasionally hang past the read timeout, and under batch-scale
    load return 429/503; retrying the same idempotent GET/POST a few times clears most
    of these transparently. A non-retryable status (e.g. a plain 4xx) still raises
    immediately via raise_for_status().

    With a deadline, each attempt's timeout is clamped to the time remaining, the
    body is streamed so the transfer itself is bounded in wall-clock time, and a
    backoff sleep that would overshoot the budget is skipped in favour of raising —
    otherwise the retries themselves become the thing that blows the deadline.
    """
    last_exc = None
    for attempt in range(retries + 1):
        try:
            if "timeout" in kwargs:
                kwargs["timeout"] = clamp_timeout(kwargs["timeout"], deadline)
            if deadline is not None:
                kwargs["stream"] = True
            resp = getattr(requests, method)(url, **kwargs)
            if resp.status_code in RETRYABLE_STATUS and attempt < retries:
                wait = retry_after_seconds(resp)
                wait = wait if wait is not None else retry_wait
                resp.close()
                _sleep_within(wait, deadline)
                continue
            resp.raise_for_status()
            # status comes from the headers, so it is checked before spending any
            # of the budget on a body that may never finish arriving
            if deadline is not None:
                try:
                    _read_within_deadline(resp, deadline)
                finally:
                    resp.close()
            return resp
        except RETRYABLE_EXCEPTIONS as exc:
            last_exc = exc
            if attempt < retries:
                _sleep_within(retry_wait, deadline)
    raise last_exc


def _sleep_within(wait, deadline):
    """Sleep for wait seconds, raising instead if that would overshoot the deadline."""
    left = remaining(deadline)
    if left is not None and wait >= left:
        raise TimeoutError(
            "timed out: retry backoff would exceed the time budget for this operation")
    time.sleep(wait)
