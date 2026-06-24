"""
conftest.py — pytest configuration for TAGSITES

Registers:
  - scripts/ on sys.path so tests can import pipeline modules directly
  - DATA fixture pointing at tests/data/
  - 'network' marker + --run-network CLI flag to gate EBI-dependent tests
  - TAGSITES_EMAIL env var / --email flag for integration tests
"""

import os
import sys
from pathlib import Path

import pytest

# Make pipeline modules importable without installing
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def DATA():
    """Path to the committed test fixture directory."""
    return REPO_ROOT / "tests" / "data"


@pytest.fixture(scope="session")
def email(request):
    """EBI e-mail for network tests: --email flag or TAGSITES_EMAIL env var."""
    return request.config.getoption("--email") or os.environ.get("TAGSITES_EMAIL", "")


# ── CLI options ───────────────────────────────────────────────────────────────

def pytest_addoption(parser):
    parser.addoption(
        "--run-network",
        action="store_true",
        default=False,
        help="Run network-dependent integration tests (requires EBI access).",
    )
    parser.addoption(
        "--email",
        action="store",
        default="",
        help="EBI e-mail address for network tests (or set TAGSITES_EMAIL env var).",
    )


# ── markers ───────────────────────────────────────────────────────────────────

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "network: mark test as requiring network access and an EBI e-mail.",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-network"):
        skip = pytest.mark.skip(reason="Pass --run-network to run network tests.")
        for item in items:
            if item.get_closest_marker("network"):
                item.add_marker(skip)
