"""
Functional tests for MLMathToolsMod :: beta_distribution_pdf.

Returns the beta distribution PDF: f(x; a, b).
"""

import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_beta_distribution_pdf.exe'


class TestBetaPDF:
    """Tests for MLMathToolsMod :: beta_distribution_pdf(a, b, x) → f(x;a,b)."""

    def test_uniform_distribution(self):
        """Beta(1,1) is the uniform distribution: f(x)=1 for all x in [0,1]."""
        for x in [0.1, 0.5, 0.9]:
            out = run_fortran(EXE, {'a': 1.0, 'b': 1.0, 'x': x})
            assert abs(out['beta_pdf'] - 1.0) < 1e-10, \
                f"Beta(1,1) PDF at x={x}: expected 1.0, got {out['beta_pdf']}"

    def test_symmetric_at_midpoint(self):
        """For symmetric Beta(a,a), PDF at x=0.5 equals PDF at x=0.5 (trivially).
        More useful: f(x; a,b) = f(1-x; b,a) — symmetry property."""
        a, b, x = 2.0, 3.0, 0.3
        f_ab   = run_fortran(EXE, {'a': a,   'b': b,   'x': x    })['beta_pdf']
        f_ba   = run_fortran(EXE, {'a': b,   'b': a,   'x': 1-x  })['beta_pdf']
        assert abs(f_ab - f_ba) < 1e-10

    def test_nonnegative(self):
        """PDF is non-negative."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 5.0, 'x': 0.4})
        assert out['beta_pdf'] >= 0.0

    def test_mode_beta22(self):
        """Beta(2,2) has mode at x=0.5 where PDF = 1.5."""
        out = run_fortran(EXE, {'a': 2.0, 'b': 2.0, 'x': 0.5})
        assert abs(out['beta_pdf'] - 1.5) < 1e-10
