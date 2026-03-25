"""
Functional tests for MLCanopyTurbulenceMod :: GetPsiRSL.

Computes RSL-modified stability functions psi for momentum and scalars between
height za and the canopy top hc.  See Bonan et al. (2018), appendix A2.

Combined formula:
  psim = -psim_MO(za) + psim_MO(hc) + psihat_m(za) - psihat_m(hc) + vkc/beta
  psic = -psic_MO(za) + psic_MO(hc) + psihat_c(za) - psihat_c(hc)

where psim_MO / psic_MO are the Monin-Obukhov psi functions and psihat_* are
the RSL correction functions (from the look-up table).

Key exact property — za = hc
------------------------------
When za = hc both the Monin-Obukhov and RSL terms cancel:
  psim(za=hc) = vkc / beta   (vkc = 0.4, the von Kármán constant)
  psic(za=hc) = 0.0

This holds for any stability (any obu ≠ 0) and any beta/PrSc.

psim2 output:
  psim2 = psim_MO((hc - disp) / obu)
which can be independently computed via the Python reference below.
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyTurbulenceMod_GetPsiRSL.exe'
VKC = 0.4   # von Kármán constant (clm_varcon)


def _psim_mo(zeta):
    """Python reference for Monin-Obukhov psi_m(zeta)."""
    if zeta >= 0.0:
        return -5.0 * zeta
    else:
        x = (1.0 - 16.0 * zeta) ** 0.25
        return (2.0 * math.log((1.0 + x) / 2.0)
                + math.log((1.0 + x**2) / 2.0)
                - 2.0 * math.atan(x)
                + math.pi / 2.0)


# ──────────────────────────────────────────────────────────────────────────────
# Exact cancellation property:  za = hc
# ──────────────────────────────────────────────────────────────────────────────

class TestCancellationAtCanopyHeight:
    """When za = hc the Monin-Obukhov and RSL terms cancel exactly."""

    def _run_at_hc(self, obu, beta=0.3, PrSc=0.5):
        return run_fortran(EXE, {'za': 20.0, 'hc': 20.0, 'disp': 12.0,
                                 'obu': obu, 'beta': beta, 'PrSc': PrSc})

    def test_psim_equals_vkc_over_beta_unstable(self):
        """Unstable (obu=-100): psim = vkc/beta = 0.4/0.3."""
        beta = 0.3
        out = self._run_at_hc(obu=-100.0, beta=beta)
        assert abs(out['psim'] - VKC / beta) < 1e-10, \
            f"psim={out['psim']}, expected {VKC/beta}"

    def test_psim_equals_vkc_over_beta_stable(self):
        """Stable (obu=+100): psim = vkc/beta."""
        beta = 0.25
        out = self._run_at_hc(obu=100.0, beta=beta)
        assert abs(out['psim'] - VKC / beta) < 1e-10, \
            f"psim={out['psim']}, expected {VKC/beta}"

    def test_psic_equals_zero_unstable(self):
        """Unstable: psic = 0 when za = hc."""
        out = self._run_at_hc(obu=-100.0)
        assert abs(out['psic']) < 1e-10, f"psic={out['psic']}, expected 0"

    def test_psic_equals_zero_stable(self):
        """Stable: psic = 0 when za = hc."""
        out = self._run_at_hc(obu=100.0)
        assert abs(out['psic']) < 1e-10, f"psic={out['psic']}, expected 0"

    def test_psim_scales_with_beta(self):
        """psim = vkc/beta → inversely proportional to beta."""
        for beta in [0.20, 0.25, 0.30, 0.35]:
            out = self._run_at_hc(obu=-100.0, beta=beta)
            assert abs(out['psim'] - VKC / beta) < 1e-10, \
                f"beta={beta}: psim={out['psim']}, expected {VKC/beta}"


# ──────────────────────────────────────────────────────────────────────────────
# psim2 reference check
# ──────────────────────────────────────────────────────────────────────────────

class TestPsim2Reference:
    """psim2 = psim_MO((hc - disp) / obu), computable in Python."""

    def test_psim2_unstable(self):
        """psim2 matches Python Monin-Obukhov reference (unstable)."""
        hc, disp, obu = 20.0, 12.0, -100.0
        expected = _psim_mo((hc - disp) / obu)
        out = run_fortran(EXE, {'za': 30.0, 'hc': hc, 'disp': disp,
                                'obu': obu, 'beta': 0.3, 'PrSc': 0.5})
        assert abs(out['psim2'] - expected) < 1e-10, \
            f"psim2={out['psim2']}, expected {expected}"

    def test_psim2_stable(self):
        """psim2 matches Python Monin-Obukhov reference (stable)."""
        hc, disp, obu = 20.0, 12.0, 100.0
        expected = _psim_mo((hc - disp) / obu)
        out = run_fortran(EXE, {'za': 30.0, 'hc': hc, 'disp': disp,
                                'obu': obu, 'beta': 0.3, 'PrSc': 0.5})
        assert abs(out['psim2'] - expected) < 1e-10, \
            f"psim2={out['psim2']}, expected {expected}"


# ──────────────────────────────────────────────────────────────────────────────
# Property tests — above-canopy measurements (za > hc)
# ──────────────────────────────────────────────────────────────────────────────

class TestProperties:

    def test_psim_increases_with_height_unstable(self):
        """Unstable: psim increases as za increases above hc
        (more wind-profile correction needed at greater height)."""
        base = {'hc': 20.0, 'disp': 12.0, 'obu': -100.0,
                'beta': 0.3, 'PrSc': 0.5}
        psim_vals = [
            run_fortran(EXE, {'za': za, **base})['psim']
            for za in [20.0, 22.0, 25.0, 30.0, 40.0]
        ]
        for i in range(len(psim_vals) - 1):
            assert psim_vals[i] <= psim_vals[i + 1], \
                f"psim not increasing: psim({20+2*i})={psim_vals[i]}, " \
                f"psim({20+2*(i+1)})={psim_vals[i+1]}"

    def test_psim_hat2_finite(self):
        """psim_hat2 (RSL correction at hc) should be a finite, reasonable value."""
        out = run_fortran(EXE, {'za': 30.0, 'hc': 20.0, 'disp': 12.0,
                                'obu': -100.0, 'beta': 0.3, 'PrSc': 0.5})
        assert math.isfinite(out['psim_hat2']), "psim_hat2 is not finite"
        assert abs(out['psim_hat2']) < 10.0, \
            f"psim_hat2={out['psim_hat2']} seems unreasonably large"

    def test_psim_psic_finite_for_all_cases(self):
        """psim and psic are finite for typical atmospheric conditions."""
        cases = [
            {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu': -100.0,
             'beta': 0.3, 'PrSc': 0.5},
            {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu':  100.0,
             'beta': 0.3, 'PrSc': 0.5},
            {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu': 1000.0,
             'beta': 0.3, 'PrSc': 0.5},
        ]
        for inp in cases:
            out = run_fortran(EXE, inp)
            assert math.isfinite(out['psim']), f"psim not finite: {inp}"
            assert math.isfinite(out['psic']), f"psic not finite: {inp}"
