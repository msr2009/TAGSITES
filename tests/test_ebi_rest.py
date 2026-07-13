"""
test_ebi_rest.py — offline unit tests for scripts/ebi_rest.py

Covers: fmt_exp, indexed_job_id_cb, combined_poll_cb (pure), submit/get_status/
fetch_result/dbfetch (requests mocked), run_job (submit/get_status/time.sleep
mocked to drive success/failure/backoff paths), resume_job (three branches).
"""

from unittest.mock import Mock

import pytest

import ebi_rest


# ── fmt_exp ────────────────────────────────────────────────────────────────────

class TestFmtExp:

    def test_simple_negative_exponent(self):
        assert ebi_rest.fmt_exp(1e-10) == "1e-10"

    def test_no_zero_padding(self):
        assert ebi_rest.fmt_exp(1e-5) == "1e-5"

    def test_non_unit_mantissa(self):
        assert ebi_rest.fmt_exp(5e-50) == "5e-50"

    def test_multi_digit_mantissa(self):
        assert ebi_rest.fmt_exp(1.5e-100) == "1.5e-100"

    def test_string_input_coerced_to_float(self):
        assert ebi_rest.fmt_exp("1e-10") == "1e-10"

    def test_positive_exponent(self):
        assert ebi_rest.fmt_exp(1e5) == "1e5"


# ── indexed_job_id_cb ─────────────────────────────────────────────────────────

class TestIndexedJobIdCb:

    def test_returns_none_when_no_callback(self):
        assert ebi_rest.indexed_job_id_cb(None, 0) is None

    def test_wraps_callback_with_index(self):
        calls = []
        cb = ebi_rest.indexed_job_id_cb(lambda i, jid: calls.append((i, jid)), 2)
        cb("job-123", "RUNNING")
        assert calls == [(2, "job-123")]

    def test_status_argument_ignored(self):
        calls = []
        cb = ebi_rest.indexed_job_id_cb(lambda i, jid: calls.append((i, jid)), 0)
        cb("job-1", "QUEUED")
        cb("job-1", "FINISHED")
        assert calls == [(0, "job-1"), (0, "job-1")]


# ── combined_poll_cb ──────────────────────────────────────────────────────────

class TestCombinedPollCb:

    def test_calls_all_callbacks_in_order(self):
        calls = []
        cb = ebi_rest.combined_poll_cb(
            lambda jid, s: calls.append(("a", jid, s)),
            lambda jid, s: calls.append(("b", jid, s)),
        )
        cb("job-1", "RUNNING")
        assert calls == [("a", "job-1", "RUNNING"), ("b", "job-1", "RUNNING")]

    def test_none_callbacks_skipped(self):
        calls = []
        cb = ebi_rest.combined_poll_cb(None, lambda jid, s: calls.append((jid, s)), None)
        cb("job-1", "FINISHED")
        assert calls == [("job-1", "FINISHED")]

    def test_no_callbacks_is_noop(self):
        cb = ebi_rest.combined_poll_cb()
        cb("job-1", "RUNNING")  # must not raise


# ── submit / get_status / fetch_result / dbfetch ──────────────────────────────

def _fake_response(text=None, content=None, status_ok=True):
    resp = Mock()
    resp.text = text
    resp.content = content
    if status_ok:
        resp.raise_for_status = Mock()
    else:
        resp.raise_for_status = Mock(side_effect=ebi_rest.requests.HTTPError("bad status"))
    return resp


class TestSubmit:

    def test_posts_to_run_endpoint_and_strips_job_id(self, monkeypatch):
        captured = {}

        def fake_post(url, data=None, timeout=None):
            captured["url"] = url
            captured["data"] = data
            return _fake_response(text="  job-42  \n")

        monkeypatch.setattr(ebi_rest.requests, "post", fake_post)
        job_id = ebi_rest.submit("https://example.org/svc", {"email": "a@b.com"})
        assert job_id == "job-42"
        assert captured["url"] == "https://example.org/svc/run"
        assert captured["data"] == {"email": "a@b.com"}

    def test_raises_on_bad_status(self, monkeypatch):
        monkeypatch.setattr(ebi_rest.requests, "post",
                            lambda *a, **k: _fake_response(text="x", status_ok=False))
        with pytest.raises(Exception):
            ebi_rest.submit("https://example.org/svc", {})


class TestGetStatus:

    def test_gets_status_endpoint(self, monkeypatch):
        captured = {}

        def fake_get(url, timeout=None):
            captured["url"] = url
            return _fake_response(text="RUNNING\n")

        monkeypatch.setattr(ebi_rest.requests, "get", fake_get)
        status = ebi_rest.get_status("https://example.org/svc", "job-42")
        assert status == "RUNNING"
        assert captured["url"] == "https://example.org/svc/status/job-42"


class TestFetchResult:

    def test_gets_result_endpoint_and_returns_bytes(self, monkeypatch):
        captured = {}

        def fake_get(url, timeout=None):
            captured["url"] = url
            return _fake_response(content=b"raw-bytes")

        monkeypatch.setattr(ebi_rest.requests, "get", fake_get)
        result = ebi_rest.fetch_result("https://example.org/svc", "job-42", "tsv")
        assert result == b"raw-bytes"
        assert captured["url"] == "https://example.org/svc/result/job-42/tsv"


class TestDbfetch:

    def test_builds_expected_url(self, monkeypatch):
        captured = {}

        def fake_get(url, timeout=None):
            captured["url"] = url
            return _fake_response(content=b">seq\nMKV\n")

        monkeypatch.setattr(ebi_rest.requests, "get", fake_get)
        result = ebi_rest.dbfetch("uniprotkb", "P12345", "fasta", "raw")
        assert result == b">seq\nMKV\n"
        assert captured["url"] == f"{ebi_rest.DBFETCH_BASE}/uniprotkb/P12345/fasta/raw"


# ── run_job ────────────────────────────────────────────────────────────────────

class TestRunJob:

    def test_success_returns_job_id_and_calls_poll_cb(self, monkeypatch):
        monkeypatch.setattr(ebi_rest, "submit", lambda base, params: "job-1")
        monkeypatch.setattr(ebi_rest.time, "sleep", lambda s: None)
        statuses = iter(["RUNNING", "RUNNING", "FINISHED"])
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: next(statuses))

        seen = []
        job_id = ebi_rest.run_job("https://x", {}, poll_cb=lambda jid, s: seen.append(s))

        assert job_id == "job-1"
        assert seen == ["QUEUED", "RUNNING", "RUNNING", "FINISHED"]

    def test_raises_runtime_error_on_ebi_error_status(self, monkeypatch):
        monkeypatch.setattr(ebi_rest, "submit", lambda base, params: "job-2")
        monkeypatch.setattr(ebi_rest.time, "sleep", lambda s: None)
        statuses = iter(["RUNNING", "ERROR"])
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: next(statuses))

        with pytest.raises(RuntimeError, match="job-2.*ERROR"):
            ebi_rest.run_job("https://x", {})

    @pytest.mark.parametrize("terminal_status", ["ERROR", "FAILURE", "NOT_FOUND"])
    def test_raises_on_each_terminal_failure_status(self, monkeypatch, terminal_status):
        monkeypatch.setattr(ebi_rest, "submit", lambda base, params: "job-3")
        monkeypatch.setattr(ebi_rest.time, "sleep", lambda s: None)
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: terminal_status)

        with pytest.raises(RuntimeError):
            ebi_rest.run_job("https://x", {})

    def test_backoff_grows_and_caps_at_max_interval(self, monkeypatch):
        monkeypatch.setattr(ebi_rest, "submit", lambda base, params: "job-4")
        sleeps = []
        monkeypatch.setattr(ebi_rest.time, "sleep", lambda s: sleeps.append(s))
        statuses = iter(["RUNNING"] * 4 + ["FINISHED"])
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: next(statuses))

        ebi_rest.run_job("https://x", {}, poll_interval=5, backoff=2, max_interval=15)

        # 5 -> 10 -> 15 (capped, would be 20) -> 15 (capped) -> 15
        assert sleeps == [5, 10, 15, 15, 15]

    def test_no_poll_cb_does_not_raise(self, monkeypatch):
        monkeypatch.setattr(ebi_rest, "submit", lambda base, params: "job-5")
        monkeypatch.setattr(ebi_rest.time, "sleep", lambda s: None)
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: "FINISHED")

        assert ebi_rest.run_job("https://x", {}) == "job-5"


# ── resume_job ─────────────────────────────────────────────────────────────────

class TestResumeJob:

    def test_finished_fetches_and_returns_result(self, monkeypatch):
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: "FINISHED")
        monkeypatch.setattr(ebi_rest, "fetch_result", lambda base, jid, rt: b"payload")

        state, payload = ebi_rest.resume_job("https://x", "job-1", "tsv")
        assert state == "finished"
        assert payload == b"payload"

    @pytest.mark.parametrize("status", ["QUEUED", "RUNNING"])
    def test_pending_statuses_return_pending(self, monkeypatch, status):
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: status)

        state, payload = ebi_rest.resume_job("https://x", "job-1", "tsv")
        assert state == "pending"
        assert payload == status

    @pytest.mark.parametrize("status", ["ERROR", "FAILURE", "NOT_FOUND"])
    def test_terminal_failure_statuses_return_expired(self, monkeypatch, status):
        monkeypatch.setattr(ebi_rest, "get_status", lambda base, jid: status)

        state, payload = ebi_rest.resume_job("https://x", "job-1", "tsv")
        assert state == "expired"
        assert payload == status
