"""
Functional tests for MLMathToolsMod.

Subroutines tested:
  quadratic           - solve ax^2 + bx + c = 0
  tridiag             - solve a tridiagonal system F*u = r
  log_gamma_function  - ln(Gamma(x))
  beta_function       - B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)
  beta_distribution_pdf  - f(x; a,b) beta PDF
  beta_distribution_cdf  - F(x; a,b) beta CDF

Each section has:
  - Reference-value tests: known analytical results verified to tight tolerance
  - Property tests: mathematical invariants (symmetry, monotonicity, bounds)
"""

import math
import pytest
from utils import run_fortran, run_fortran_array


# ──────────────────────────────────────────────────────────────────────────────
# quadratic
# ──────────────────────────────────────────────────────────────────────────────

class TestQuadratic:
    """Tests for MLMathToolsMod :: quadratic(a, b, c) → r1, r2."""

    def test_integer_roots(self):
        """x^2 - 3x + 2 = 0  →  roots 1, 2."""
        out = run_fortran('test_quadratic.exe', {'a': 1.0, 'b': -3.0, 'c': 2.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - 1.0) < 1e-10
        assert abs(roots[1] - 2.0) < 1e-10

    def test_repeated_root(self):
        """x^2 - 2x + 1 = 0  →  double root at 1."""
        out = run_fortran('test_quadratic.exe', {'a': 1.0, 'b': -2.0, 'c': 1.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - 1.0) < 1e-10
        assert abs(roots[1] - 1.0) < 1e-10

    def test_negative_roots(self):
        """x^2 + 5x + 6 = 0  →  roots -2, -3."""
        out = run_fortran('test_quadratic.exe', {'a': 1.0, 'b': 5.0, 'c': 6.0})
        roots = sorted([out['r1'], out['r2']])
        assert abs(roots[0] - (-3.0)) < 1e-10
        assert abs(roots[1] - (-2.0)) < 1e-10

    def test_vieta_product(self):
        """Product of roots = c/a (Vieta's formula)."""
        a, b, c = 2.0, -7.0, 3.0
        out = run_fortran('test_quadratic.exe', {'a': a, 'b': b, 'c': c})
        assert abs(out['r1'] * out['r2'] - c / a) < 1e-10

    def test_vieta_sum(self):
        """Sum of roots = -b/a (Vieta's formula)."""
        a, b, c = 3.0, -11.0, 6.0
        out = run_fortran('test_quadratic.exe', {'a': a, 'b': b, 'c': c})
        assert abs(out['r1'] + out['r2'] - (-b / a)) < 1e-10


# ──────────────────────────────────────────────────────────────────────────────
# tridiag
# ──────────────────────────────────────────────────────────────────────────────

class TestTridiag:
    """Tests for MLMathToolsMod :: tridiag(a, b, c, r, u, n)."""

    def test_3x3_symmetric(self):
        """Solve the 3x3 tridiagonal system with known solution u=[1,1,1].

        Matrix:  diag=[2,2,2],  off-diag=[-1,-1]
        RHS:     r=[1,0,1]
        Solution: u=[1,1,1]  (verify: 2-1=1, -1+2-1=0, -1+2=1)
        """
        inputs = {
            'n': 3,
            'a': [0.0, -1.0, -1.0],   # a(1) not referenced
            'b': [2.0,  2.0,  2.0],
            'c': [-1.0, -1.0, 0.0],   # c(3) not referenced
            'r': [1.0,  0.0,  1.0],
        }
        u = run_fortran_array('test_tridiag.exe', inputs, 'u')
        for val in u:
            assert abs(val - 1.0) < 1e-12

    def test_diagonal_system(self):
        """Pure diagonal system (a=c=0): solution is r_i/b_i."""
        inputs = {
            'n': 4,
            'a': [0.0, 0.0, 0.0, 0.0],
            'b': [2.0, 3.0, 4.0, 5.0],
            'c': [0.0, 0.0, 0.0, 0.0],
            'r': [6.0, 9.0, 12.0, 15.0],
        }
        u = run_fortran_array('test_tridiag.exe', inputs, 'u')
        expected = [3.0, 3.0, 3.0, 3.0]
        for i, (val, exp) in enumerate(zip(u, expected)):
            assert abs(val - exp) < 1e-12, f"u[{i+1}]: got {val}, expected {exp}"

    def test_identity_rhs(self):
        """Identity matrix (b=1, a=c=0): solution equals rhs."""
        r_vals = [7.0, 3.0, 5.0]
        inputs = {
            'n': 3,
            'a': [0.0, 0.0, 0.0],
            'b': [1.0, 1.0, 1.0],
            'c': [0.0, 0.0, 0.0],
            'r': r_vals,
        }
        u = run_fortran_array('test_tridiag.exe', inputs, 'u')
        for val, exp in zip(u, r_vals):
            assert abs(val - exp) < 1e-12


# ──────────────────────────────────────────────────────────────────────────────
# log_gamma_function
# ──────────────────────────────────────────────────────────────────────────────

class TestLogGamma:
    """Tests for MLMathToolsMod :: log_gamma_function(x) → ln(Gamma(x))."""

    def test_at_one(self):
        """ln(Gamma(1)) = ln(1) = 0."""
        out = run_fortran('test_log_gamma.exe', {'x': 1.0})
        assert abs(out['gammaln']) < 1e-10

    def test_at_two(self):
        """ln(Gamma(2)) = ln(1!) = 0."""
        out = run_fortran('test_log_gamma.exe', {'x': 2.0})
        assert abs(out['gammaln']) < 1e-10

    def test_at_three(self):
        """ln(Gamma(3)) = ln(2!) = ln(2)."""
        out = run_fortran('test_log_gamma.exe', {'x': 3.0})
        assert abs(out['gammaln'] - math.log(2.0)) < 1e-7

    def test_at_half(self):
        """ln(Gamma(0.5)) = 0.5*ln(pi)."""
        out = run_fortran('test_log_gamma.exe', {'x': 0.5})
        assert abs(out['gammaln'] - 0.5 * math.log(math.pi)) < 1e-7

    def test_positive_for_x_less_than_one(self):
        """Gamma(x) > 1 (hence ln > 0) for 0 < x < 1."""
        out = run_fortran('test_log_gamma.exe', {'x': 0.3})
        assert out['gammaln'] > 0.0

    def test_increases_for_large_x(self):
        """ln(Gamma) is monotonically increasing for x > 1."""
        g5  = run_fortran('test_log_gamma.exe', {'x': 5.0})['gammaln']
        g10 = run_fortran('test_log_gamma.exe', {'x': 10.0})['gammaln']
        assert g10 > g5


# ──────────────────────────────────────────────────────────────────────────────
# beta_function
# ──────────────────────────────────────────────────────────────────────────────

class TestBetaFunction:
    """Tests for MLMathToolsMod :: beta_function(a, b) → B(a,b)."""

    def test_unit_arguments(self):
        """B(1,1) = 1."""
        out = run_fortran('test_beta_function.exe', {'a': 1.0, 'b': 1.0})
        assert abs(out['beta'] - 1.0) < 1e-10

    def test_half_integer(self):
        """B(0.5, 0.5) = pi."""
        out = run_fortran('test_beta_function.exe', {'a': 0.5, 'b': 0.5})
        assert abs(out['beta'] - math.pi) < 1e-7

    def test_integer_arguments(self):
        """B(m,n) = (m-1)!(n-1)! / (m+n-1)! for integers.
        B(2,3) = 1!*2!/4! = 1*2/24 = 1/12."""
        out = run_fortran('test_beta_function.exe', {'a': 2.0, 'b': 3.0})
        assert abs(out['beta'] - 1.0 / 12.0) < 1e-10

    def test_symmetry(self):
        """B(a,b) = B(b,a)."""
        ab = run_fortran('test_beta_function.exe', {'a': 2.0, 'b': 5.0})['beta']
        ba = run_fortran('test_beta_function.exe', {'a': 5.0, 'b': 2.0})['beta']
        assert abs(ab - ba) < 1e-12

    def test_positive(self):
        """B(a,b) > 0 for a,b > 0."""
        out = run_fortran('test_beta_function.exe', {'a': 3.7, 'b': 2.1})
        assert out['beta'] > 0.0


# ──────────────────────────────────────────────────────────────────────────────
# beta_distribution_pdf
# ──────────────────────────────────────────────────────────────────────────────

class TestBetaPDF:
    """Tests for MLMathToolsMod :: beta_distribution_pdf(a, b, x) → f(x;a,b)."""

    def test_uniform_distribution(self):
        """Beta(1,1) is the uniform distribution: f(x)=1 for all x in [0,1]."""
        for x in [0.1, 0.5, 0.9]:
            out = run_fortran('test_beta_pdf.exe', {'a': 1.0, 'b': 1.0, 'x': x})
            assert abs(out['beta_pdf'] - 1.0) < 1e-10, \
                f"Beta(1,1) PDF at x={x}: expected 1.0, got {out['beta_pdf']}"

    def test_symmetric_at_midpoint(self):
        """For symmetric Beta(a,a), PDF at x=0.5 equals PDF at x=0.5 (trivially).
        More useful: f(x; a,b) = f(1-x; b,a) — symmetry property."""
        a, b, x = 2.0, 3.0, 0.3
        f_ab   = run_fortran('test_beta_pdf.exe', {'a': a,   'b': b,   'x': x    })['beta_pdf']
        f_ba   = run_fortran('test_beta_pdf.exe', {'a': b,   'b': a,   'x': 1-x  })['beta_pdf']
        assert abs(f_ab - f_ba) < 1e-10

    def test_nonnegative(self):
        """PDF is non-negative."""
        out = run_fortran('test_beta_pdf.exe', {'a': 2.0, 'b': 5.0, 'x': 0.4})
        assert out['beta_pdf'] >= 0.0

    def test_mode_beta22(self):
        """Beta(2,2) has mode at x=0.5 where PDF = 1.5."""
        out = run_fortran('test_beta_pdf.exe', {'a': 2.0, 'b': 2.0, 'x': 0.5})
        assert abs(out['beta_pdf'] - 1.5) < 1e-10


# ──────────────────────────────────────────────────────────────────────────────
# beta_distribution_cdf
# ──────────────────────────────────────────────────────────────────────────────

class TestBetaCDF:
    """Tests for MLMathToolsMod :: beta_distribution_cdf(a, b, x) → F(x;a,b)."""

    def test_lower_bound(self):
        """F(0; a,b) = 0."""
        out = run_fortran('test_beta_cdf.exe', {'a': 2.0, 'b': 3.0, 'x': 0.0})
        assert abs(out['beta_cdf']) < 1e-10

    def test_upper_bound(self):
        """F(1; a,b) = 1."""
        out = run_fortran('test_beta_cdf.exe', {'a': 2.0, 'b': 3.0, 'x': 1.0})
        assert abs(out['beta_cdf'] - 1.0) < 1e-10

    def test_midpoint_symmetric(self):
        """For symmetric Beta(a,a), F(0.5) = 0.5."""
        out = run_fortran('test_beta_cdf.exe', {'a': 3.0, 'b': 3.0, 'x': 0.5})
        assert abs(out['beta_cdf'] - 0.5) < 1e-7

    def test_monotone_increasing(self):
        """CDF is non-decreasing."""
        xs = [0.1, 0.3, 0.5, 0.7, 0.9]
        vals = [run_fortran('test_beta_cdf.exe',
                            {'a': 2.0, 'b': 2.0, 'x': x})['beta_cdf']
                for x in xs]
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i + 1], \
                f"CDF not monotone: F({xs[i]})={vals[i]} > F({xs[i+1]})={vals[i+1]}"

    def test_beta22_half(self):
        """B(2,2): F(0.5) = 0.5 (by symmetry), F(0.25) = 5/32 = 0.15625.
        B(2,2) CDF: F(x) = 3x^2 - 2x^3."""
        out = run_fortran('test_beta_cdf.exe', {'a': 2.0, 'b': 2.0, 'x': 0.25})
        expected = 3 * 0.25**2 - 2 * 0.25**3   # = 0.15625
        assert abs(out['beta_cdf'] - expected) < 1e-7
