"""
Functional tests for MLMathToolsMod :: beta_function.

Returns B(a, b) = Gamma(a)*Gamma(b) / Gamma(a+b).
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_beta_function.exe'


class TestBetaFunction:
    """Tests for MLMathToolsMod :: beta_function(a, b) → B(a,b)."""

    def test_unit_arguments(self):
        """B(1,1) = 1."""
        out = run_fortran(EXE, {'a': 1.0, 'b': 1.0})
        assert abs(out['beta'] - 1.0) < 1e-10

    def test_half_integer(self):
        """B(0.5, 0.5) = pi."""
        out = run_fortran(EXE, {'a': 0.5, 'b': 0.5})
        assert abs(out['beta'] - math.pi) < 1e-7

    def test_integer_arguments(self):
        """B(m,n) = (m-1)!(n-1)! / (m+n-1)! for integers.
        B(2,3) = 1!*2!/4! = 1*2/24 = 1/12."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 3.0})
        assert abs(out['beta'] - 1.0 / 12.0) < 1e-10

    def test_symmetry(self):
        """B(a,b) = B(b,a)."""
        ab = run_fortran(EXE, {'a': 2.0, 'b': 5.0})['beta']
        ba = run_fortran(EXE, {'a': 5.0, 'b': 2.0})['beta']
        assert abs(ab - ba) < 1e-12

    def test_positive(self):
        """B(a,b) > 0 for a,b > 0."""
        out = run_fortran(EXE, {'a': 3.7, 'b': 2.1})
        assert out['beta'] > 0.0
