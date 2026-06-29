"""
test_results_loader.py — offline unit tests for utils/results.py

Covers: _translate_phobius_desc, _merge_phobius_intervals,
        load_data_from_json (smoke test with temp output files),
        load_run_metadata, _guess_analysis_type, _hex_to_rgb,
        _isoform_class_key.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

# utils/ is at the repo root, not in scripts/ which conftest adds
sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from utils.results import (
    _translate_phobius_desc,
    _merge_phobius_intervals,
    load_data_from_json,
    load_run_metadata,
    _guess_analysis_type,
    _hex_to_rgb,
    _isoform_class_key,
)
from config import RESULTS_TYPE_DICT


# ── _translate_phobius_desc ───────────────────────────────────────────────────

class TestTranslatePhobiusDesc:

    def test_cytoplasm(self):
        assert _translate_phobius_desc("cytoplasm region") == "Cytoplasmic"

    def test_extracellular(self):
        assert _translate_phobius_desc("extracellular region") == "Extracellular"

    def test_transmembrane(self):
        assert _translate_phobius_desc("embedded in the membrane") == "Transmembrane"

    def test_signal_peptide(self):
        assert _translate_phobius_desc("signal peptide region") == "Signal peptide"

    def test_signal_sequence(self):
        assert _translate_phobius_desc("signal sequence") == "Signal peptide"

    def test_case_insensitive(self):
        assert _translate_phobius_desc("CYTOPLASM region") == "Cytoplasmic"

    def test_unknown_desc_returned_unchanged(self):
        assert _translate_phobius_desc("some unknown description") == "some unknown description"


# ── _merge_phobius_intervals ──────────────────────────────────────────────────

class TestMergePhobiusIntervals:

    def _df(self, rows):
        return pd.DataFrame(rows, columns=["source", "start", "stop", "description"])

    def test_single_interval_passes_through(self):
        df = self._df([("Phobius", 1, 10, "Cytoplasmic")])
        merged = _merge_phobius_intervals(df)
        assert len(merged) == 1
        assert merged.iloc[0]["start"] == 1
        assert merged.iloc[0]["stop"] == 10

    def test_overlapping_same_desc_merged(self):
        df = self._df([
            ("Phobius", 1, 10, "Cytoplasmic"),
            ("Phobius", 8, 15, "Cytoplasmic"),
        ])
        merged = _merge_phobius_intervals(df)
        assert len(merged) == 1
        assert merged.iloc[0]["stop"] == 15

    def test_contiguous_intervals_merged(self):
        df = self._df([
            ("Phobius", 1, 5, "Cytoplasmic"),
            ("Phobius", 6, 10, "Cytoplasmic"),
        ])
        merged = _merge_phobius_intervals(df)
        assert len(merged) == 1

    def test_non_overlapping_different_desc_kept_separate(self):
        df = self._df([
            ("Phobius", 1, 5, "Cytoplasmic"),
            ("Phobius", 20, 30, "Transmembrane"),
        ])
        merged = _merge_phobius_intervals(df)
        assert len(merged) == 2

    def test_source_column_set_to_phobius(self):
        df = self._df([("Phobius", 1, 5, "Cytoplasmic")])
        merged = _merge_phobius_intervals(df)
        assert (merged["source"] == "Phobius").all()

    def test_columns_correct(self):
        df = self._df([("Phobius", 1, 5, "Cytoplasmic")])
        merged = _merge_phobius_intervals(df)
        assert set(merged.columns) == {"source", "start", "stop", "description"}


# ── _guess_analysis_type ──────────────────────────────────────────────────────

class TestGuessAnalysisType:

    def test_blast_in_name(self):
        assert _guess_analysis_type("WORM_blast") == "blast"

    def test_plddt_in_name(self):
        assert _guess_analysis_type("HUM_plddt") == "plddt"

    def test_sasa_suffix_maps_to_plddt(self):
        assert _guess_analysis_type("HUM_plddt_sasa") == "plddt"

    def test_scores_in_name(self):
        assert _guess_analysis_type("ANY_scores") == "scores"

    def test_unknown_returns_unknown(self):
        assert _guess_analysis_type("RANDOM_XYZ") == "unknown"


# ── _hex_to_rgb ───────────────────────────────────────────────────────────────

class TestHexToRgb:

    def test_black(self):
        assert _hex_to_rgb("#000000") == (0, 0, 0)

    def test_white(self):
        assert _hex_to_rgb("#ffffff") == (255, 255, 255)

    def test_known_color(self):
        r, g, b = _hex_to_rgb("#2e7d32")
        assert r == 0x2e
        assert g == 0x7d
        assert b == 0x32

    def test_uppercase_accepted(self):
        assert _hex_to_rgb("#FF0000") == (255, 0, 0)


# ── _isoform_class_key ────────────────────────────────────────────────────────

class TestIsoformClassKey:

    def test_constitutive_prefix(self):
        assert _isoform_class_key("constitutive exon") == "constitutive"

    def test_constitutive_case_insensitive(self):
        assert _isoform_class_key("Constitutive region") == "constitutive"

    def test_unique_prefix(self):
        assert _isoform_class_key("unique to isoform A") == "unique"

    def test_unique_case_insensitive(self):
        assert _isoform_class_key("UNIQUE segment") == "unique"

    def test_anything_else_is_intermediate(self):
        assert _isoform_class_key("intermediate region") == "intermediate"

    def test_unknown_defaults_to_intermediate(self):
        assert _isoform_class_key("some unrecognised description") == "intermediate"


# ── load_run_metadata ─────────────────────────────────────────────────────────

class TestLoadRunMetadata:

    def _make_json(self, tmp_path, fasta_content=None, pdb=""):
        fa = tmp_path / "run1.fa"
        seq = "MSTGKKVL"
        if fasta_content:
            fa.write_text(fasta_content)
        else:
            fa.write_text(f">run1\n{seq}\n")
        j = {
            "global": {
                "email": "t@t.com", "run_name": "run1",
                "working_dir": str(tmp_path) + "/",
                "input_file": str(fa),
                "pdb": pdb, "scripts_folder": "./scripts/",
                "selected_sites": [],
            },
            "tasks": {}
        }
        return j, seq

    def test_query_seq_extracted(self, tmp_path):
        j, seq = self._make_json(tmp_path)
        meta = load_run_metadata(j)
        assert meta["query_seq"] == seq

    def test_seq_len_correct(self, tmp_path):
        j, seq = self._make_json(tmp_path)
        meta = load_run_metadata(j)
        assert meta["seq_len"] == len(seq)

    def test_pdb_path_empty_when_none(self, tmp_path):
        j, _ = self._make_json(tmp_path)
        meta = load_run_metadata(j)
        assert meta["pdb_path"] == ""

    def test_pdb_path_returned_when_exists(self, tmp_path):
        pdb = tmp_path / "run1.pdb"
        pdb.write_text("ATOM ...")
        j, _ = self._make_json(tmp_path, pdb=str(pdb))
        meta = load_run_metadata(j)
        assert meta["pdb_path"] == str(pdb)

    def test_accepts_dict(self, tmp_path):
        j, seq = self._make_json(tmp_path)
        meta = load_run_metadata(j)
        assert meta["query_seq"] == seq

    def test_accepts_json_string(self, tmp_path):
        j, seq = self._make_json(tmp_path)
        meta = load_run_metadata(json.dumps(j))
        assert meta["query_seq"] == seq

    def test_missing_fasta_returns_empty_seq(self, tmp_path):
        j = {
            "global": {
                "email": "t@t.com", "run_name": "run1",
                "working_dir": str(tmp_path) + "/",
                "input_file": str(tmp_path / "missing.fa"),
                "pdb": "", "scripts_folder": "./", "selected_sites": [],
            },
            "tasks": {}
        }
        meta = load_run_metadata(j)
        assert meta["query_seq"] == ""
        assert meta["seq_len"] == 0


# ── load_data_from_json — smoke test ──────────────────────────────────────────

class TestLoadDataFromJson:

    def _minimal_continuous_tsv(self, tmp_path, n=5):
        """Write a minimal continuous output TSV (pos, score per row)."""
        out = tmp_path / "run1.A_scores.tsv"
        lines = "\n".join(f"pos{i+1}\t{float(i)}" for i in range(n))
        out.write_text(lines + "\n")
        return str(out)

    def _minimal_range_tsv(self, tmp_path):
        """Write a minimal range output TSV (source, start, stop, description)."""
        out = tmp_path / "run1.B_mods.txt"
        out.write_text("Phospho\t10\t12\tphosphorylation\n")
        return str(out)

    def test_continuous_task_produces_aa_df(self, tmp_path):
        out = self._minimal_continuous_tsv(tmp_path)
        j = {
            "global": {},
            "tasks": {
                "A_scores": {"type": "scores", "args": {"output": out}}
            }
        }
        aa_df, range_df, alns = load_data_from_json(j, RESULTS_TYPE_DICT)
        assert "A_scores" in aa_df.columns
        assert len(aa_df) == 5

    def test_range_task_produces_range_df(self, tmp_path):
        out = self._minimal_range_tsv(tmp_path)
        j = {
            "global": {},
            "tasks": {
                "B_mods": {"type": "modifications", "args": {"output": out}}
            }
        }
        aa_df, range_df, alns = load_data_from_json(j, RESULTS_TYPE_DICT)
        assert len(range_df) == 1
        assert range_df.iloc[0]["description"] == "phosphorylation"

    def test_missing_output_skipped(self, tmp_path):
        j = {
            "global": {},
            "tasks": {
                "A_blast": {"type": "blast", "args": {"output": str(tmp_path / "missing.jsd")}}
            }
        }
        aa_df, range_df, alns = load_data_from_json(j, RESULTS_TYPE_DICT)
        assert len(aa_df) == 0

    def test_accepts_json_string(self, tmp_path):
        out = self._minimal_continuous_tsv(tmp_path)
        j = {"global": {}, "tasks": {"A_scores": {"type": "scores", "args": {"output": out}}}}
        aa_df, _, _ = load_data_from_json(json.dumps(j), RESULTS_TYPE_DICT)
        assert "A_scores" in aa_df.columns

    def test_accepts_dict(self, tmp_path):
        out = self._minimal_continuous_tsv(tmp_path)
        j = {"global": {}, "tasks": {"A_scores": {"type": "scores", "args": {"output": out}}}}
        aa_df, _, _ = load_data_from_json(j, RESULTS_TYPE_DICT)
        assert "A_scores" in aa_df.columns

    def test_empty_tasks_returns_empty_dfs(self):
        j = {"global": {}, "tasks": {}}
        aa_df, range_df, alns = load_data_from_json(j, RESULTS_TYPE_DICT)
        assert len(aa_df) == 0
        assert len(range_df) == 0
        assert alns == []
