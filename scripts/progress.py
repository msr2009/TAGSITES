"""progress.py — structured progress reporting for TAGSITES analysis scripts.

Convention: every progress line is  [stage] message
  level info    →  [stage] message
  level warning →  [stage] WARNING: message
  level error   →  [stage] ERROR: message

parse_line() is the inverse: split off the bracketed stage token.
Lines without a leading '[' are passed through verbatim (library noise, tracebacks).
"""

import sys


def _format_line(message, stage, level):
    """Build the '[stage] message' string from components."""
    prefix = f"[{stage}] " if stage else ""
    if level == "warning":
        return f"{prefix}WARNING: {message}"
    if level == "error":
        return f"{prefix}ERROR: {message}"
    return f"{prefix}{message}"


def report(reporter, message, stage="", level="info"):
    """Emit one progress event through the given reporter callable."""
    reporter(_format_line(message, stage, level), stage, level)


def make_stderr_reporter():
    """Return a reporter that writes formatted lines to stderr (standalone default)."""
    def _reporter(line, stage, level):
        print(line, file=sys.stderr, flush=True)
    return _reporter


def make_status_reporter(wd, rn, tid, log_lines):
    """Return a reporter that appends to log_lines and writes stage+log to the status file.

    Imports run_status lazily so this module stays usable without it on sys.path.
    """
    import run_status as _rs

    def _reporter(line, stage, level):
        log_lines.append(line)
        _rs.update_task(wd, rn, tid, stage=stage, log="\n".join(log_lines))
    return _reporter


def resolve_reporter(reporter):
    """Return reporter unchanged, or a stderr reporter when None."""
    return reporter if reporter is not None else make_stderr_reporter()


def poll_adapter(reporter):
    """Return a poll_cb(job_id, status_str) that routes EBI poll events through reporter."""
    def _cb(job_id, status_str):
        report(reporter, f"job {job_id}: {status_str}", stage="ebi_poll")
    return _cb


def parse_line(line):
    """Parse '[stage] message' → (stage, message, level); untagged lines → ('', line, 'info')."""
    line = line.rstrip("\n")
    if not line.startswith("["):
        return "", line, "info"
    close = line.find("]")
    if close < 0:
        return "", line, "info"
    stage = line[1:close]
    rest = line[close + 1:].lstrip(" ")
    if rest.startswith("WARNING: "):
        return stage, rest[9:], "warning"
    if rest.startswith("ERROR: "):
        return stage, rest[7:], "error"
    return stage, rest, "info"
