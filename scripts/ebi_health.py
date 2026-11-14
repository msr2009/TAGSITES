"""ebi_health.py — quick check that EBI's Job Dispatcher is reachable and accepting jobs.

Two probes, both free of fair-use cost (neither creates a job), ~1 s in total:
  1. GET {service}/parameters for every service the app uses: reachability + latency.
  2. POST an EMPTY form to clustalo/run. EBI rejects it either way, but the reason
     tells the two cases apart: a validation error means submissions are accepted,
     "temporarily blocked from submitting" means EBI refuses this submitter.

The empty POST needs a REAL email (TAGSITES_EBI_PROBE_EMAIL, or --email on the CLI);
without one it is skipped and only the GET probes and real-job observations count.
EBI answers "temporarily blocked" to any reserved placeholder domain (example.org,
example.com, ...) whatever the IP, so a placeholder probe email reports a block that
does not exist. That false positive was shipped on 2026-11-14; such addresses are
now refused here (_usable_email).

check() caches its verdict process-wide, so every Shiny session shares one probe
per TTL. ebi_rest also calls note() when a real job hits trouble (a submission
refused, the status backend rate-limited), and those recent observations are
folded into the verdict.

The empty POST runs on a much slower clock than the GETs: an address under a
fair-use block should not keep knocking, even with requests that create no job,
in case EBI counts the attempts. It re-runs every SUBMIT_TTL_BLOCKED_S while
blocked, every SUBMIT_TTL_S otherwise (a real refused job still flags a new
block at once via note()), and on check(force=True) no more than once per
SUBMIT_MIN_GAP_S, so repeated "Check again" clicks cannot hammer EBI.

Levels: "ok", "slow" (reachable but sluggish), "degraded" (some services or the
status backend misbehaving), "blocked" (submissions refused), "down" (unreachable).

CLI:  python scripts/ebi_health.py   (exit 0 when ok/slow, 1 otherwise)
"""

import argparse
import concurrent.futures
import os
import re
import sys
import threading
import time

import requests

REST = "https://www.ebi.ac.uk/Tools/services/rest"
SERVICES = ("ncbiblast", "iprscan5", "clustalo", "genewise")
# reserved / placeholder domains EBI treats as blocked submitters (see docstring)
_RESERVED_EMAIL = re.compile(
    r"@(?:[^@]*\.)?(?:example\.(?:org|com|net)|invalid|test|localhost)$", re.IGNORECASE
)


def _usable_email(email):
    """The email if it can probe EBI meaningfully, else None (missing or placeholder)."""
    email = (email or "").strip()
    if "@" not in email or _RESERVED_EMAIL.search(email):
        return None
    return email


# submission-probe email; None disables the empty-POST test
PROBE_EMAIL = _usable_email(os.environ.get("TAGSITES_EBI_PROBE_EMAIL"))

TIMEOUT_S = 10  # per probe request
SLOW_S = 4.0  # a /parameters call slower than this counts as sluggish
TTL_OK_S = 300  # re-probe a healthy EBI every 5 min
TTL_BAD_S = 60  # re-probe sooner while something is wrong (GETs only)
SUBMIT_TTL_S = 1800  # empty-POST re-test while not blocked
SUBMIT_TTL_BLOCKED_S = 900  # empty-POST re-test while blocked
SUBMIT_MIN_GAP_S = 60  # floor between empty POSTs, even when forced
NOTE_WINDOW_S = 900  # job-time observations older than this are dropped

_LEVEL_RANK = {"ok": 0, "slow": 1, "degraded": 2, "blocked": 3, "down": 4}

_lock = threading.Lock()
_cache = {"at": 0.0, "result": None}  # GET-probe verdict, shared by every email
_notes = {}  # kind -> (time, detail)
_submit = {}  # email -> {"at", "accepting", "detail"}: last empty-POST outcome per email


def note(kind, detail="", email=None):
    """Record a problem a real EBI job ran into (kind: 'blocked' or 'status_busy')."""
    with _lock:
        _notes[kind] = (time.time(), detail)
        # a refused real job is a fresh submission test for its email: no need to probe
        if kind == "blocked" and _usable_email(email):
            _submit[_usable_email(email)] = {
                "at": time.time(),
                "accepting": False,
                "detail": detail,
            }
        # a fresh observation makes the cached verdict stale
        _cache["at"] = 0.0


def error_description(text):
    """The <description> text of an EBI XML error body, or '' when there is none."""
    m = re.search(r"<description>(.*?)</description>", text or "", re.DOTALL)
    return m.group(1).strip() if m else ""


def is_block_message(text):
    """True when an EBI error body is the fair-use submission block."""
    return "temporarily blocked" in (text or "").lower()


def _probe_parameters(service):
    """(service, seconds, error-or-None) for one GET /parameters."""
    t0 = time.time()
    try:
        r = requests.get(f"{REST}/{service}/parameters", timeout=TIMEOUT_S)
        dt = time.time() - t0
        if r.status_code != 200:
            return service, dt, f"HTTP {r.status_code}"
        return service, dt, None
    except requests.RequestException as exc:
        return service, time.time() - t0, type(exc).__name__


def _probe_submission(email):
    """(accepting, detail) from an empty POST to clustalo/run; accepting is None if unknown."""
    try:
        r = requests.post(f"{REST}/clustalo/run", data={"email": email}, timeout=TIMEOUT_S)
    except requests.RequestException as exc:
        return None, type(exc).__name__
    desc = error_description(r.text)
    if is_block_message(r.text):
        return False, desc
    # a 4xx validation complaint is the healthy answer to an empty form
    if r.status_code < 500:
        return True, desc
    return None, f"HTTP {r.status_code}"


def _submission_due(email, force):
    """True when the empty-POST test for this email should run now (see module docstring)."""
    with _lock:
        last = _submit.get(email)
    # never tested with this email: test once
    if last is None:
        return True
    since = time.time() - last["at"]
    if since < SUBMIT_MIN_GAP_S:
        return False
    if force:
        return True
    return since > (SUBMIT_TTL_BLOCKED_S if last["accepting"] is False else SUBMIT_TTL_S)


def probe():
    """Run the GET probes now (uncached); return the reachability verdict dict."""
    with concurrent.futures.ThreadPoolExecutor(len(SERVICES)) as pool:
        params = list(pool.map(_probe_parameters, SERVICES))

    failed = [(s, err) for s, _, err in params if err]
    slow = [(s, dt) for s, dt, err in params if not err and dt > SLOW_S]
    latency = max((dt for _, dt, err in params if not err), default=None)

    level, problems = "ok", []
    # every service unreachable: nothing else matters
    if len(failed) == len(SERVICES):
        level = "down"
        problems.append(
            "EBI's job services are not responding ("
            + ", ".join(f"{s}: {e}" for s, e in failed)
            + ")."
        )
    else:
        if failed:
            level = "degraded"
            problems.append("Not responding: " + ", ".join(f"{s} ({e})" for s, e in failed) + ".")
        if slow:
            level = max(level, "slow", key=_LEVEL_RANK.get)
            problems.append(
                "Slow to respond: " + ", ".join(f"{s} ({dt:.1f} s)" for s, dt in slow) + "."
            )
    return {"level": level, "problems": problems, "latency_s": latency, "checked_at": time.time()}


def _add_submission(result, email):
    """Fold the email's last empty-POST outcome into a verdict (returns a new dict)."""
    out = dict(result, problems=list(result["problems"]), accepting=None, email=email)
    if email is None:
        return out
    with _lock:
        last = dict(_submit.get(email) or {})
    out["accepting"] = last.get("accepting")
    if last.get("accepting") is False and out["level"] != "down":
        out["level"] = "blocked"
        out["problems"].insert(
            0,
            f"EBI is refusing new jobs for {email}: "
            + (last.get("detail") or "temporarily blocked"),
        )
    return out


def _merge_notes(result):
    """Fold recent job-time observations into a probe verdict (returns a new dict)."""
    out = dict(result, problems=list(result["problems"]))
    now = time.time()
    with _lock:
        notes = {k: v for k, v in _notes.items() if now - v[0] < NOTE_WINDOW_S}
    if "blocked" in notes and out["level"] not in ("blocked", "down"):
        # the probe can miss a block that applies to one service only
        mins = max(1, int((now - notes["blocked"][0]) / 60))
        out["level"] = "blocked"
        out["problems"].append(
            f"A job submission was refused {mins} min ago: "
            f"{notes['blocked'][1] or 'temporarily blocked'}"
        )
    if "status_busy" in notes and _LEVEL_RANK[out["level"]] < _LEVEL_RANK["degraded"]:
        mins = max(1, int((now - notes["status_busy"][0]) / 60))
        out["level"] = "degraded"
        out["problems"].append(
            f"EBI's job-status backend was overloaded {mins} min ago "
            "(jobs still run, but progress updates may stall)."
        )
    return out


def check(force=False, email=None):
    """Cached verdict: {'level', 'problems', 'latency_s', 'accepting', 'email', 'checked_at'}.

    email: the submitter to test (e.g. the app's Email field); falls back to
    TAGSITES_EBI_PROBE_EMAIL. Placeholder or missing -> the submission test is skipped
    and 'accepting' is None.
    """
    email = _usable_email(email) or PROBE_EMAIL
    with _lock:
        cached, at = _cache["result"], _cache["at"]
    ttl = TTL_OK_S if cached and cached["level"] == "ok" else TTL_BAD_S
    # re-run the GETs when forced, never probed, or the cached verdict has expired
    if force or cached is None or time.time() - at > ttl:
        cached = probe()
        with _lock:
            _cache["result"], _cache["at"] = cached, time.time()
    # the empty POST runs on its own, slower, per-email clock
    if email and _submission_due(email, force):
        accepting, detail = _probe_submission(email)
        with _lock:
            _submit[email] = {"at": time.time(), "accepting": accepting, "detail": detail}
    return _merge_notes(_add_submission(cached, email))


def summary(result):
    """One-line human description of a verdict."""
    head = {
        "ok": "EBI services are working normally.",
        "slow": "EBI services are responding slowly; analyses may take longer than usual.",
        "degraded": "EBI services are partly impaired; some analyses may fail or stall.",
        "blocked": "EBI is not accepting new jobs; EBI-backed "
        "analyses (BLAST, InterPro, Clustal Omega, Genewise) will fail.",
        "down": "EBI services are unreachable; EBI-backed analyses will fail.",
    }[result["level"]]
    return head


def main():
    """CLI: print the verdict and exit non-zero when EBI is not usable."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--email",
        help="real email for the submission probe "
        "(default: TAGSITES_EBI_PROBE_EMAIL; omitted -> probe skipped)",
    )
    args = ap.parse_args()
    # a placeholder would report a false block: refuse it instead of probing
    if args.email and _usable_email(args.email) is None:
        ap.error(f"{args.email!r} is a placeholder address; EBI answers it as blocked")
    t0 = time.time()
    res = check(force=True, email=args.email)
    print(f"[{res['level'].upper()}] {summary(res)}")
    for p in res["problems"]:
        print(f"  - {p}")
    if res["email"] is None:
        print("  submission test skipped: no real email (pass --email to run it)")
    else:
        verdict = {True: "accepted", False: "refused"}.get(res["accepting"], "unknown")
        print(f"  submissions for {res['email']}: {verdict}")
    lat = res["latency_s"]
    print(
        f"  probe took {time.time() - t0:.1f} s"
        + (f"; slowest service {lat:.2f} s" if lat is not None else "")
    )
    return 0 if res["level"] in ("ok", "slow") else 1


if __name__ == "__main__":
    sys.exit(main())
