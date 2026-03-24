"""
Functional tests for MLMathToolsMod :: quadratic.

Solves ax^2 + bx + c = 0 for its two roots r1, r2.
"""

import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_quadratic.exe'


class TestQuadratic:
    """Tests for MLMathToolsMod :: quadratic(a, b, c) → r1, r2."""

    def test_integer_roots(self):
        """x^2 - 3x + 2 = 0  →  roots 1, 2."""
        out = run_fortran(EXE, {'a': 1.0, 'b': -3.0, 'c': 2.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - 1.0) < 1e-10
        assert abs(roots[1] - 2.0) < 1e-10

    def test_repeated_root(self):
        """x^2 - 2x + 1 = 0  →  double root at 1."""
        out = run_fortran(EXE, {'a': 1.0, 'b': -2.0, 'c': 1.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - 1.0) < 1e-10
        assert abs(roots[1] - 1.0) < 1e-10

    def test_negative_roots(self):
        """x^2 + 5x + 6 = 0  →  roots -2, -3."""
        out = run_fortran(EXE, {'a': 1.0, 'b': 5.0, 'c': 6.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - (-3.0)) < 1e-10
        assert abs(roots[1] - (-2.0)) < 1e-10

    def test_vieta_product(self):
        """Product of roots = c/a (Vieta's formula)."""
        a, b, c = 2.0, -7.0, 3.0
        out = run_fortran(EXE, {'a': a, 'b': b, 'c': c})
        assert abs(out['r1'] * out['r2'] - c / a) < 1e-10

    def test_vieta_sum(self):
        """Sum of roots = -b/a (Vieta's formula)."""
        a, b, c = 3.0, -11.0, 6.0
        out = run_fortran(EXE, {'a': a, 'b': b, 'c': c})
        assert abs(out['r1'] + out['r2'] - (-b / a)) < 1e-10
