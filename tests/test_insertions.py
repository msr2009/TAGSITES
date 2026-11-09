"""Tests for query-insertion detection (utils/insertions.py) and its scoring credit."""

import pandas as pd

from utils.insertions import (INSERTION_SOURCE, absent_counts, insertion_rows,
                              parse_description)
from utils.scoring import _handle_column_ramp, _insertion_positions

QUERY = "ACDEFGHIKLMN"
# 12 orthologs lack query columns 2-4 (D, E, F)
DEEP = "AC---GHIKLMN"
# fragments that start at column 6 never span the insertion
FRAGMENT = "------HIKLMN"


def _write(tmp_path, n_deep=12, n_frag=2, header=True):
    """Write a FASTA alignment and its .jsd; return (jsd_path, aln_path)."""
    seqs = [("Q", QUERY)]
    seqs += [(f"deep{i}", DEEP) for i in range(n_deep)]
    seqs += [(f"frag{i}", FRAGMENT) for i in range(n_frag)]
    aln = tmp_path / "x.aln"
    aln.write_text("".join(f">{n}\n{s}\n" for n, s in seqs))
    jsd = tmp_path / "x.jsd"
    lines = ["# reference sequence: Q\n", "# align_column_number\tscore\tamino acid\n"] if header else []
    lines += [f"{i}\t0.5\t{aa}\n" for i, aa in enumerate(QUERY)]
    jsd.write_text("".join(lines))
    return str(jsd), str(aln)


def test_absent_counts_deep_insertion(tmp_path):
    """Columns 2-4 are lacked by every spanning ortholog."""
    n_absent, n_span = absent_counts(*_write(tmp_path))
    assert list(n_absent[2:5]) == [12, 12, 12]
    assert list(n_absent[:2]) == [0, 0] and list(n_absent[5:]) == [0] * 7


def test_fragments_do_not_count_as_absent(tmp_path):
    """A fragment that ends or starts away from a column neither spans nor lacks it."""
    n_absent, n_span = absent_counts(*_write(tmp_path))
    assert n_span[0] == 12       # fragments start at column 6
    assert n_span[2] == 12       # not 14
    assert n_span[8] == 14       # fragments span here


def test_headerless_jsd_finds_query_by_length(tmp_path):
    """A database-export .jsd has no header; the query is the row matching its length."""
    with_header = absent_counts(*_write(tmp_path))
    without = absent_counts(*_write(tmp_path, header=False))
    assert list(with_header[0]) == list(without[0])


def test_missing_reference_returns_none(tmp_path):
    """An unknown reference name or an unreadable alignment is no data, not an error."""
    jsd, aln = _write(tmp_path)
    (tmp_path / "x.jsd").write_text("# reference sequence: nope\n0\t0.5\tA\n")
    assert absent_counts(jsd, aln) is None
    assert absent_counts(jsd, str(tmp_path / "missing.aln")) is None


def test_insertion_rows_runs_and_description(tmp_path):
    """NaN residues with absent orthologs collapse into one run; others are excluded."""
    jsd, aln = _write(tmp_path)
    rows = insertion_rows(jsd, aln, {1, 3, 4, 5})
    assert len(rows) == 1                      # residue 1 has no absent orthologs
    source, start, stop, desc = rows[0]
    assert (source, start, stop) == (INSERTION_SOURCE, 3, 5)
    assert parse_description(desc) == (12, 1.0)


def test_parse_description_rejects_other_text():
    """Only the writer's own format parses."""
    assert parse_description("Transmembrane") is None


def _ctx(desc):
    """Scoring context with five residues, NaN conservation at 2-4 and one insertion row."""
    aa = pd.DataFrame({"pos": [1, 2, 3, 4, 5], "BROAD_blast": [0.9, None, None, None, 0.5]})
    rng = pd.DataFrame([[INSERTION_SOURCE, 2, 4, desc]],
                       columns=["source", "start", "stop", "description"])
    return {"aa_df": aa, "range_df": rng, "seq_len": 5}


PARAMS = {"analysis": "blast", "low": 0.635, "high": 0.831,
          "insertion": {"min_absent": 10, "min_absent_fraction": 0.75}}
DEEP_DESC = "absent from 79 spanning orthologs (79% of spanning)"
SHALLOW_DESC = "absent from 2 spanning orthologs (100% of spanning)"


def test_deep_insertion_earns_full_point():
    """Deep insertions score 1.0; ordinary residues keep the ramp."""
    out = _handle_column_ramp([1, 2, 3, 4, 5], PARAMS, _ctx(DEEP_DESC))
    assert out[0] == 0.0 and out[4] == 1.0
    assert out[1:4] == [1.0, 1.0, 1.0]


def test_shallow_insertion_earns_nothing():
    """Two absent hits is noise: the residues stay missing (NaN)."""
    out = _handle_column_ramp([1, 2, 3, 4, 5], PARAMS, _ctx(SHALLOW_DESC))
    assert all(v != v for v in out[1:4])


def test_fraction_gate():
    """Enough absent orthologs still fail when they are a minority of the spanning ones."""
    desc = "absent from 12 spanning orthologs (40% of spanning)"
    assert _insertion_positions(_ctx(desc), PARAMS["insertion"]) == set()
    assert _insertion_positions(_ctx(DEEP_DESC), PARAMS["insertion"]) == {2, 3, 4}


def test_no_insertion_block_keeps_old_behaviour():
    """Without the config block missing conservation stays NaN."""
    params = {k: v for k, v in PARAMS.items() if k != "insertion"}
    out = _handle_column_ramp([1, 2, 3, 4, 5], params, _ctx(DEEP_DESC))
    assert all(v != v for v in out[1:4])
