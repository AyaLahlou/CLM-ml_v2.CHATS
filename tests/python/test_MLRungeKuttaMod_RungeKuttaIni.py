"""
Functional tests for MLRungeKuttaMod :: RungeKuttaIni.

Initialises the Butcher tableau (a, b, c) for the selected Runge-Kutta method.
The method is a compile-time constant in MLclm_varctl:
  runge_kutta_type = 41  →  4th-order classical Kutta method, nrk = 4

4th-order Kutta Butcher tableau
--------------------------------
Weights b:
  b = [1/6, 2/6, 2/6, 1/6]             sum = 1  (consistency condition)

Stage times c:
  c = [0, 1/2, 1/2, 1]                 c(1)=0 by convention

Stage coefficients a (lower triangle only; diagonal and upper = spval ≈ 1e36):
  a(2,1) = 1/2
  a(3,1) = 0,   a(3,2) = 1/2
  a(4,1) = 0,   a(4,2) = 0,   a(4,3) = 1

Consistency condition (row-sum): c(i) = Σ a(i,j) for j < i, i ≥ 2.

Reference: Kutta (1901); Butcher (1996).
"""

import pytest
from utils import run_fortran

EXE  = 'test_MLRungeKuttaMod_RungeKuttaIni.exe'
NRK  = 4        # nrk = runge_kutta_type/10 = 41/10 = 4
SPVAL = 1.0e36  # sentinel for unused Butcher-tableau entries


def _out():
    """Run the executable (no inputs) and return the output dict."""
    return run_fortran(EXE, {})


# ──────────────────────────────────────────────────────────────────────────────
# Exact coefficient checks — 4th-order classical Kutta
# ──────────────────────────────────────────────────────────────────────────────

class TestKutta4thOrderCoefficients:

    def test_b_weights_exact(self):
        """b = [1/6, 2/6, 2/6, 1/6] for 4th-order Kutta."""
        out = _out()
        expected = [1/6, 2/6, 2/6, 1/6]
        for i, exp in enumerate(expected, start=1):
            assert abs(out[f'b_{i}'] - exp) < 1e-15, \
                f"b_{i}: got {out[f'b_{i}']}, expected {exp}"

    def test_c_times_exact(self):
        """c = [0, 1/2, 1/2, 1] for 4th-order Kutta."""
        out = _out()
        expected = [0.0, 0.5, 0.5, 1.0]
        for i, exp in enumerate(expected, start=1):
            assert abs(out[f'c_{i}'] - exp) < 1e-15, \
                f"c_{i}: got {out[f'c_{i}']}, expected {exp}"

    def test_a_lower_triangle_exact(self):
        """Lower-triangle a-coefficients for 4th-order Kutta."""
        out = _out()
        lower = {
            (2, 1): 0.5,
            (3, 1): 0.0,
            (3, 2): 0.5,
            (4, 1): 0.0,
            (4, 2): 0.0,
            (4, 3): 1.0,
        }
        for (i, j), exp in lower.items():
            key = f'a_{i}_{j}'
            assert abs(out[key] - exp) < 1e-15, \
                f"a[{i},{j}]: got {out[key]}, expected {exp}"

    def test_a_upper_triangle_is_spval(self):
        """Upper triangle (and diagonal) of a must be spval (unused sentinel)."""
        out = _out()
        for i in range(1, NRK + 1):
            for j in range(i, NRK + 1):   # diagonal and above
                key = f'a_{i}_{j}'
                assert out[key] > 1e35, \
                    f"a[{i},{j}] should be spval (~1e36), got {out[key]}"


# ──────────────────────────────────────────────────────────────────────────────
# Mathematical consistency conditions
# ──────────────────────────────────────────────────────────────────────────────

class TestButcherConsistency:

    def test_sum_b_equals_one(self):
        """Consistency: sum of weights b = 1."""
        out = _out()
        b_sum = sum(out[f'b_{i}'] for i in range(1, NRK + 1))
        assert abs(b_sum - 1.0) < 1e-14, f"sum(b) = {b_sum}, expected 1.0"

    def test_c1_equals_zero(self):
        """First stage time c(1) must be 0 (evaluated at t_n)."""
        out = _out()
        assert abs(out['c_1']) < 1e-15, f"c_1 = {out['c_1']}, expected 0"

    def test_c_equals_row_sum_of_a(self):
        """Row-sum consistency: c(i) = Σ a(i,j) for j = 1..i-1, i ≥ 2."""
        out = _out()
        for i in range(2, NRK + 1):
            row_sum = sum(
                out[f'a_{i}_{j}']
                for j in range(1, i)
                if out[f'a_{i}_{j}'] < 1e35   # skip spval
            )
            ci = out[f'c_{i}']
            assert abs(row_sum - ci) < 1e-14, \
                f"Row-sum a[{i},:] = {row_sum} ≠ c[{i}] = {ci}"

    def test_c_last_equals_one(self):
        """For 4th-order Kutta: c(nrk) = 1 (last stage at t_n + h)."""
        out = _out()
        assert abs(out[f'c_{NRK}'] - 1.0) < 1e-15, \
            f"c_{NRK} = {out[f'c_{NRK}']}, expected 1.0"
