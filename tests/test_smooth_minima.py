"""
test_smooth_minima.py — offline unit tests for scripts/pLDDT_minima.py

Bug fix tested:
  smooth_iterative(data, windowsize) previously called list(int), raising
  TypeError. Fixed to wrap int windowsize in [windowsize].

Also tests:
  split_minima  — recursively clusters nearby minima
  find_insertion_sites — finds local minima in smoothed data
"""

import math
import numpy as np
import pytest

from pLDDT_minima import smooth_iterative, split_minima, find_insertion_sites


# ── smooth_iterative ──────────────────────────────────────────────────────────

class TestSmoothIterative:

    def test_int_windowsize_no_longer_raises(self):
        """Regression: smooth_iterative([...], 3) must not raise TypeError."""
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = smooth_iterative(data, 3)
        assert result is not None

    def test_int_windowsize_correct_values(self):
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = smooth_iterative(data, 3)
        expected = np.array([2.0, 3.0, 4.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_list_windowsize_still_works(self):
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = smooth_iterative(data, [3])
        expected = np.array([2.0, 3.0, 4.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_output_length(self):
        data = list(range(10))
        result = smooth_iterative(data, 3)
        assert len(result) == len(data) - 3 + 1  # == 8

    def test_flat_signal_unchanged(self):
        """Smoothing a constant signal should return the same constant."""
        data = [5.0] * 10
        result = smooth_iterative(data, 3)
        np.testing.assert_array_almost_equal(result, [5.0] * 8)

    def test_iterative_two_passes(self):
        """Passing a list of two window sizes applies two rounds of smoothing."""
        data = [1.0, 2.0, 3.0, 2.0, 1.0, 2.0, 3.0, 2.0, 1.0]
        result_double = smooth_iterative(data, [3, 3])
        result_single = smooth_iterative(data, 3)
        # Two passes should produce a shorter, smoother array
        assert len(result_double) < len(result_single)


# ── split_minima ──────────────────────────────────────────────────────────────

class TestSplitMinima:

    def test_single_minimum_returned_as_is(self):
        result = split_minima([5], min_dist=10)
        assert len(result) == 1

    def test_well_separated_minima_preserved(self):
        """Two minima far apart should survive as separate groups."""
        result = split_minima([0, 100], min_dist=10)
        assert len(result) == 2

    def test_close_minima_merged(self):
        """Two minima too close together should be merged into one."""
        result = split_minima([10, 11, 12], min_dist=100)
        assert len(result) == 1


# ── find_insertion_sites ──────────────────────────────────────────────────────

class TestFindInsertionSites:

    def _signal_with_two_minima(self):
        """Create a signal with clear minima at positions ~20 and ~80."""
        signal = np.ones(100) * 80.0
        signal[20] = 10.0   # first minimum
        signal[80] = 15.0   # second minimum
        return signal

    def test_finds_minima(self):
        signal = self._signal_with_two_minima()
        sites = find_insertion_sites(signal, window_size=1, num_minima=5,
                                     min_dist=30, max_score=100)
        assert len(sites) > 0

    def test_respects_max_score(self):
        """Sites with smoothed score above max_score should be excluded."""
        signal = np.ones(100) * 80.0
        signal[20] = 70.0   # above threshold 60 → excluded
        sites = find_insertion_sites(signal, window_size=1, num_minima=5,
                                     min_dist=5, max_score=60)
        # The minimum at 70 is above max_score=60, so should not appear
        assert all(signal[int(s)] < 60 for s in sites)

    def test_returns_list(self):
        signal = self._signal_with_two_minima()
        sites = find_insertion_sites(signal, window_size=1, num_minima=5,
                                     min_dist=5)
        assert isinstance(sites, list)

    def test_int_window_size_accepted(self):
        """find_insertion_sites passes window_size to smooth_iterative;
        int windowsize must work after the bug fix."""
        signal = self._signal_with_two_minima()
        # Should not raise
        sites = find_insertion_sites(signal, window_size=3, num_minima=2,
                                     min_dist=10, max_score=100)
        assert isinstance(sites, list)
