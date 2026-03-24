"""
Functional tests for MLMathToolsMod :: log_gamma_function.

Returns ln(Gamma(x)).
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_log_gamma_function.exe'


class TestLogGamma:
    """Tests for MLMathToolsMod :: log_gamma_function(x) → ln(Gamma(x))."""

    def test_at_one(self):
        """ln(Gamma(1)) = ln(1) = 0."""
        out = run_fortran(EXE, {'x': 1.0})
        assert abs(out['gammaln']) < 1e-10

    def test_at_two(self):
        """ln(Gamma(2)) = ln(1!) = 0."""
        out = run_fortran(EXE, {'x': 2.0})
        assert abs(out['gammaln']) < 1e-10

    def test_at_three(self):
        """ln(Gamma(3)) = ln(2!) = ln(2)."""
        out = run_fortran(EXE, {'x': 3.0})
        assert abs(out['gammaln'] - math.log(2.0)) < 1e-7

    def test_at_half(self):
        """ln(Gamma(0.5)) = 0.5*ln(pi)."""
        out = run_fortran(EXE, {'x': 0.5})
        assert abs(out['gammaln'] - 0.5 * math.log(math.pi)) < 1e-7

    def test_positive_for_x_less_than_one(self):
        """Gamma(x) > 1 (hence ln > 0) for 0 < x < 1."""
        out = run_fortran(EXE, {'x': 0.3})
        assert out['gammaln'] > 0.0

    def test_increases_for_large_x(self):
        """ln(Gamma) is monotonically increasing for x > 1."""
        g5  = run_fortran(EXE, {'x': 5.0})['gammaln']
        g10 = run_fortran(EXE, {'x': 10.0})['gammaln']
        assert g10 > g5
