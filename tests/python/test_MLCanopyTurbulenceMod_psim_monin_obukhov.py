"""
Functional tests for MLCanopyTurbulenceMod :: psim_monin_obukhov and psic_monin_obukhov.

Monin-Obukhov psi integrated stability correction functions
(Bonan et al. 2018, eq. A12-A13):

  psi_m(zeta):
    stable   (zeta >= 0):  -5*zeta
    unstable (zeta <  0):  2*ln((1+x)/2) + ln((1+x^2)/2) - 2*atan(x) + pi/2
                           where x = (1-16*zeta)^(1/4)

  psi_c(zeta):
    stable   (zeta >= 0):  -5*zeta
    unstable (zeta <  0):  2*ln((1+x^2)/2)
                           where x = (1-16*zeta)^(1/4)

Reference values:
  zeta=0:  psi_m=0, psi_c=0  (neutral)
  zeta=1:  psi_m=-5, psi_c=-5 (stable, linear)
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyTurbulenceMod_psim_monin_obukhov.exe'


class TestPsiMO:
    """Tests for psim_monin_obukhov and psic_monin_obukhov."""

    # --- Reference value tests (neutral and stable — analytic formulas) ---

    def test_neutral(self):
        """At zeta=0 (neutral): psim=0, psic=0 (no stability correction)."""
        out = run_fortran(EXE, {'zeta': 0.0})
        assert abs(out['psim']) < 1e-12
        assert abs(out['psic']) < 1e-12

    def test_stable_zeta1(self):
        """At zeta=1 (stable): psim = psic = -5."""
        out = run_fortran(EXE, {'zeta': 1.0})
        assert abs(out['psim'] - (-5.0)) < 1e-12
        assert abs(out['psic'] - (-5.0)) < 1e-12

    def test_stable_zeta_half(self):
        """At zeta=0.5 (stable): psim = psic = -2.5."""
        out = run_fortran(EXE, {'zeta': 0.5})
        assert abs(out['psim'] - (-2.5)) < 1e-12
        assert abs(out['psic'] - (-2.5)) < 1e-12

    def test_unstable_psim_reference(self):
        """At zeta=-1: hand-calculated psim.
        x = 17^(1/4) = 2.03054...
        psim = 2*ln((1+x)/2) + ln((1+x^2)/2) - 2*atan(x) + pi/2
        """
        x = 17.0 ** 0.25
        expected = (2.0 * math.log((1.0 + x) / 2.0)
                    + math.log((1.0 + x * x) / 2.0)
                    - 2.0 * math.atan(x)
                    + math.pi / 2.0)
        out = run_fortran(EXE, {'zeta': -1.0})
        assert abs(out['psim'] - expected) < 1e-10, \
            f"psim(-1): got {out['psim']}, expected {expected}"

    def test_unstable_psic_reference(self):
        """At zeta=-1: psic = 2*ln((1+x^2)/2), x = 17^(1/4)."""
        x = 17.0 ** 0.25
        expected = 2.0 * math.log((1.0 + x * x) / 2.0)
        out = run_fortran(EXE, {'zeta': -1.0})
        assert abs(out['psic'] - expected) < 1e-10, \
            f"psic(-1): got {out['psic']}, expected {expected}"

    # --- Property tests ---

    def test_psi_negative_stable(self):
        """In stable conditions, both psi functions are negative."""
        for zeta in [0.1, 0.5, 1.0, 2.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert out['psim'] < 0.0, f"psim({zeta}) = {out['psim']} should be < 0"
            assert out['psic'] < 0.0, f"psic({zeta}) = {out['psic']} should be < 0"

    def test_psi_positive_unstable(self):
        """In unstable conditions, both psi functions are positive."""
        for zeta in [-0.1, -0.5, -1.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert out['psim'] > 0.0, f"psim({zeta}) = {out['psim']} should be > 0"
            assert out['psic'] > 0.0, f"psic({zeta}) = {out['psic']} should be > 0"

    def test_psim_equals_psic_stable(self):
        """For zeta >= 0, both use the same formula: psi = -5*zeta."""
        for zeta in [0.0, 0.3, 1.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert abs(out['psim'] - out['psic']) < 1e-12

    def test_psim_greater_than_psic_unstable(self):
        """For zeta < 0: psim > psic (greater correction for momentum)."""
        for zeta in [-0.5, -1.0, -2.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert out['psim'] > out['psic'], \
                f"zeta={zeta}: psim={out['psim']} should exceed psic={out['psic']}"
