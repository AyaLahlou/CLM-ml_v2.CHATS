"""
Functional tests for MLCanopyTurbulenceMod :: GetBeta.

Computes beta = u*/u(h), the ratio of friction velocity to wind speed at
the canopy top, as a function of atmospheric stability.
See Bonan et al. (2018), eqs. (A22)-(A24).

Stability is parameterised by LcL = Lc/obu (canopy length scale / Obukhov length):
  LcL = 0  → neutral
  LcL < 0  → unstable (enhanced mixing)
  LcL > 0  → stable   (suppressed mixing)

The defining equation that GetBeta satisfies:
  beta * phi_m(LcL * beta^2) = beta_neutral

where phi_m is the Monin-Obukhov momentum stability function:
  stable   (zeta >= 0):  phi_m(zeta) = 1 + 5*zeta
  unstable (zeta <  0):  phi_m(zeta) = (1 - 16*zeta)^(-1/4)
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyTurbulenceMod_GetBeta.exe'


def _phim(zeta):
    """Python reference for Monin-Obukhov phi_m."""
    if zeta >= 0.0:
        return 1.0 + 5.0 * zeta
    else:
        return (1.0 - 16.0 * zeta) ** (-0.25)


class TestGetBeta:
    """Tests for MLCanopyTurbulenceMod :: GetBeta."""

    # --- Reference value tests ---

    def test_neutral_beta_equals_beta_neutral(self):
        """At LcL=0 (neutral), beta = beta_neutral exactly.
        Proof: unstable branch gives beta = sqrt((-bb+sqrt(bb^2+4*beta_n^4))/2)
        with bb=0 → beta = sqrt(beta_n^2) = beta_n.
        """
        for bn in [0.2, 0.3, 0.35]:
            out = run_fortran(EXE, {'beta_neutral': bn, 'LcL': 0.0})
            assert abs(out['beta'] - bn) < 1e-12, \
                f"beta_neutral={bn}: beta={out['beta']}, expected {bn}"

    def test_defining_equation_satisfied(self):
        """beta * phi_m(LcL * beta^2) = beta_neutral for all stability regimes."""
        cases = [
            (0.3, 0.0),    # neutral
            (0.3, -1.0),   # unstable
            (0.3,  1.0),   # stable
            (0.25, -0.5),
            (0.30,  0.5),
            (0.35, -2.0),
        ]
        for bn, lcl in cases:
            out = run_fortran(EXE, {'beta_neutral': bn, 'LcL': lcl})
            beta = out['beta']
            zeta = lcl * beta ** 2
            lhs = beta * _phim(zeta)
            assert abs(lhs - bn) < 1e-6, \
                f"beta_neutral={bn}, LcL={lcl}: beta*phim={lhs}, expected {bn}"

    # --- Property tests ---

    def test_unstable_beta_gt_beta_neutral(self):
        """For LcL < 0 (unstable): beta > beta_neutral
        (enhanced turbulence increases u*/u(h)).
        """
        for lcl in [-0.5, -1.0, -2.0]:
            out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl})
            assert out['beta'] > 0.3, \
                f"LcL={lcl}: beta={out['beta']} should exceed beta_neutral=0.3"

    def test_stable_beta_lt_beta_neutral(self):
        """For LcL > 0 (stable): beta < beta_neutral
        (suppressed turbulence reduces u*/u(h)).
        """
        for lcl in [0.5, 1.0, 2.0]:
            out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl})
            assert out['beta'] < 0.3, \
                f"LcL={lcl}: beta={out['beta']} should be less than beta_neutral=0.3"

    def test_beta_positive(self):
        """beta must be positive for all stability regimes."""
        for lcl in [-2.0, -1.0, 0.0, 0.5, 1.0, 2.0]:
            out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl})
            assert out['beta'] > 0.0, f"LcL={lcl}: beta={out['beta']} should be > 0"

    def test_beta_monotone_decreasing_with_stability(self):
        """beta decreases monotonically as LcL increases (more stable → less mixing)."""
        lcls = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
        betas = [run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl})['beta']
                 for lcl in lcls]
        for i in range(len(betas) - 1):
            assert betas[i] > betas[i + 1], \
                f"beta not decreasing: beta({lcls[i]})={betas[i]} <= beta({lcls[i+1]})={betas[i+1]}"
