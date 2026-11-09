"""
test_uniprot_api_fetch.py — offline unit tests for uniprot_api.fetch_fastas

Covers: record parsing, de-duplication, 1000-per-call chunking, whole-batch 400 fallback,
and per-accession fallback for accessions missing from a batch response.
"""

from unittest.mock import Mock

import pytest
import requests

import uniprot_api


def _fasta(records):
    """Build a UniProt-style multi-record FASTA from {acc: seq}."""
    return "".join(f">sp|{a}|{a}_X desc OS=x\n{s[:5]}\n{s[5:]}\n" for a, s in records.items())


def _resp(text):
    """Mock response carrying text."""
    return Mock(text=text)


@pytest.fixture
def batch_calls(monkeypatch):
    """Patch the batch HTTP call to serve from SEQS; record each call's accession list."""
    state = {"calls": [], "seqs": {}}

    def fake(method, url, params=None, **kw):
        accs = params["accessions"].split(",")
        state["calls"].append(accs)
        return _resp(_fasta({a: state["seqs"][a] for a in sorted(accs) if a in state["seqs"]}))

    monkeypatch.setattr(uniprot_api, "request_with_retries", fake)
    return state


class TestParse:

    def test_multiline_sequence_joined(self):
        out = uniprot_api._parse_uniprot_fasta(_fasta({"P12345": "MKVLAAGGHH"}))
        assert out == {"P12345": "MKVLAAGGHH"}

    def test_isoform_accession_kept(self):
        out = uniprot_api._parse_uniprot_fasta(">sp|P04637-2|P53_HUMAN iso\nMEEP\n")
        assert out == {"P04637-2": "MEEP"}


class TestFetchFastas:

    def test_returns_sequences_keyed_by_accession(self, batch_calls):
        batch_calls["seqs"] = {"A1": "MKVLAAGGHH", "B2": "MSTNPKPQRK"}
        assert uniprot_api.fetch_fastas(["B2", "A1"]) == batch_calls["seqs"]

    def test_duplicates_requested_once(self, batch_calls):
        batch_calls["seqs"] = {"A1": "MKVLAAGGHH"}
        uniprot_api.fetch_fastas(["A1", "A1", "A1"])
        assert batch_calls["calls"] == [["A1"]]

    def test_chunks_at_batch_limit(self, batch_calls):
        accs = [f"A{i}" for i in range(2500)]
        batch_calls["seqs"] = {a: "MKVLAAGGHH" for a in accs}
        out = uniprot_api.fetch_fastas(accs)
        assert [len(c) for c in batch_calls["calls"]] == [1000, 1000, 500]
        assert len(out) == 2500

    def test_empty_list_makes_no_call(self, batch_calls):
        assert uniprot_api.fetch_fastas([]) == {}
        assert batch_calls["calls"] == []

    def test_rejected_batch_falls_back_per_accession(self, monkeypatch):
        def fake(method, url, params=None, **kw):
            raise requests.exceptions.HTTPError("400")

        singles = {"A1": "MKVL", "B2": "MSTN"}
        monkeypatch.setattr(uniprot_api, "request_with_retries", fake)
        monkeypatch.setattr(uniprot_api, "fetch_fasta",
                            lambda a: singles[a] if a in singles else (_ for _ in ()).throw(
                                ValueError("bad")))
        # the malformed id fails its single fetch and is dropped; the rest are recovered
        assert uniprot_api.fetch_fastas(["A1", "BAD", "B2"]) == singles

    def test_missing_from_batch_uses_single_fetch(self, batch_calls, monkeypatch):
        batch_calls["seqs"] = {"A1": "MKVLAAGGHH"}  # OLD1 (secondary id) not returned
        monkeypatch.setattr(uniprot_api, "fetch_fasta", lambda a: "MSTNPKPQRK")
        out = uniprot_api.fetch_fastas(["A1", "OLD1"])
        # result is keyed by the requested accession, not the primary id
        assert out == {"A1": "MKVLAAGGHH", "OLD1": "MSTNPKPQRK"}

    def test_unfetchable_accession_absent(self, batch_calls, monkeypatch):
        def boom(a):
            raise requests.exceptions.HTTPError("404")

        monkeypatch.setattr(uniprot_api, "fetch_fasta", boom)
        assert uniprot_api.fetch_fastas(["NOPE1"]) == {}
