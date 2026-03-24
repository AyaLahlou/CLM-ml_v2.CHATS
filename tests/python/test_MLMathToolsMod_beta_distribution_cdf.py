"""
Functional tests for MLMathToolsMod :: beta_distribution_cdf.

Returns the regularised incomplete beta function I_x(a,b), i.e. the CDF
of the beta distribution at x: F(x; a, b).
"""

import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_beta_distribution_cdf.exe'


class TestBetaCDF:
    """Tests for MLMathToolsMod :: beta_distribution_cdf(a, b, x) → F(x;a,b)."""

    def test_lower_bound(self):
        """F(0; a,b) = 0."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 3.0, 'x': 0.0})
        assert abs(out['beta_cdf']) < 1e-10

    def test_upper_bound(self):
        """F(1; a,b) = 1."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 3.0, 'x': 1.0})
        assert abs(out['beta_cdf'] - 1.0) < 1e-10

    def test_midpoint_symmetric(self):
        """For symmetric Beta(a,a), F(0.5) = 0.5."""
        out = run_fortran(EXE, {'a': 3.0, 'b': 3.0, 'x': 0.5})
        assert abs(out['beta_cdf'] - 0.5) < 1e-7

    def test_monotone_increasing(self):
        """CDF is non-decreasing."""
        xs = [0.1, 0.3, 0.5, 0.7, 0.9]
        vals = [run_fortran(EXE, {'a': 2.0, 'b': 2.0, 'x': x})['beta_cdf']
                for x in xs]
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i + 1], \
                f"CDF not monotone: F({xs[i]})={vals[i]} > F({xs[i+1]})={vals[i+1]}"

    def test_beta22_half(self):
        """B(2,2): F(0.5) = 0.5 (by symmetry), F(0.25) = 5/32 = 0.15625.
        B(2,2) CDF: F(x) = 3x^2 - 2x^3."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 2.0, 'x': 0.25})
        expected = 3 * 0.25**2 - 2 * 0.25**3   # = 0.15625
        assert abs(out['beta_cdf'] - expected) < 1e-7
