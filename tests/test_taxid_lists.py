"""
test_taxid_lists.py — offline validation of taxid_lists/*.txt against the
format resolve_taxids() expects (see scripts/site_selection_util.py):
one taxid per line, optional trailing '#' comment, blank lines allowed.
"""

from pathlib import Path

import pytest

from site_selection_util import resolve_taxids

_TAXID_LISTS_DIR = Path(__file__).parent.parent / "taxid_lists"
_LIST_FILES = sorted(_TAXID_LISTS_DIR.glob("*.txt"))


def _taxids(path):
    """Parse a taxid list file the same way resolve_taxids() does."""
    with open(path) as f:
        return [ln.split("#")[0].strip() for ln in f if ln.split("#")[0].strip()]


class TestTaxidListsDirectory:

    def test_directory_exists(self):
        assert _TAXID_LISTS_DIR.is_dir()

    def test_at_least_one_list_present(self):
        assert len(_LIST_FILES) > 0


@pytest.mark.parametrize("path", _LIST_FILES, ids=lambda p: p.stem)
class TestEachTaxidList:

    def test_not_empty(self, path):
        assert len(_taxids(path)) > 0

    def test_all_entries_are_numeric(self, path):
        for t in _taxids(path):
            assert t.isdigit(), f"non-numeric taxid {t!r} in {path.name}"

    def test_no_duplicate_taxids(self, path):
        ids = _taxids(path)
        assert len(ids) == len(set(ids)), f"duplicate taxids in {path.name}"

    def test_readable_by_resolve_taxids(self, path):
        merged = resolve_taxids("", str(path))
        assert merged == ",".join(_taxids(path))


class TestBroadMetazoanPanel:

    PATH = _TAXID_LISTS_DIR / "broad_metazoan_panel.txt"

    def test_file_exists(self):
        assert self.PATH.exists()

    def test_expected_taxid_count(self):
        assert len(_taxids(self.PATH)) == 51

    def test_contains_human(self):
        assert "9606" in _taxids(self.PATH)

    def test_contains_c_elegans_relatives(self):
        ids = _taxids(self.PATH)
        # C. briggsae and C. remanei, used as close outgroups to C. elegans
        assert "6238" in ids
        assert "31234" in ids
