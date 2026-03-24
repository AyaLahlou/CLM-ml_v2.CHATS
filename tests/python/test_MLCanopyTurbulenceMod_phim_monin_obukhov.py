"""
Functional tests for MLCanopyTurbulenceMod :: phim_monin_obukhov and phic_monin_obukhov.

Monin-Obukhov phi stability correction functions (Bonan et al. 2018, eq. A10-A11):

  phi_m(zeta):
    stable   (zeta >= 0):  1 + 5*zeta
    unstable (zeta <  0):  (1 - 16*zeta)^(-1/4)

  phi_c(zeta):
    stable   (zeta >= 0):  1 + 5*zeta
    unstable (zeta <  0):  (1 - 16*zeta)^(-1/2)

Reference values:
  zeta=0:  phi_m=1, phi_c=1  (neutral)
  zeta=1:  phi_m=6, phi_c=6  (stable, linear)
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyTurbulenceMod_phim_monin_obukhov.exe'


class TestPhiMO:
    """Tests for phim_monin_obukhov and phic_monin_obukhov."""

    # --- Reference value tests (neutral and stable — analytic formulas) ---

    def test_neutral(self):
        """At zeta=0 (neutral), phim=1 and phic=1."""
        out = run_fortran(EXE, {'zeta': 0.0})
        assert abs(out['phim'] - 1.0) < 1e-12
        assert abs(out['phic'] - 1.0) < 1e-12

    def test_stable_zeta1(self):
        """At zeta=1, phim = phic = 1 + 5*1 = 6."""
        out = run_fortran(EXE, {'zeta': 1.0})
        assert abs(out['phim'] - 6.0) < 1e-12
        assert abs(out['phic'] - 6.0) < 1e-12

    def test_stable_zeta_half(self):
        """At zeta=0.5, phim = phic = 1 + 5*0.5 = 3.5."""
        out = run_fortran(EXE, {'zeta': 0.5})
        assert abs(out['phim'] - 3.5) < 1e-12
        assert abs(out['phic'] - 3.5) < 1e-12

    def test_unstable_phim(self):
        """At zeta=-1: phim = (1+16)^(-1/4) = 17^(-0.25)."""
        out = run_fortran(EXE, {'zeta': -1.0})
        expected_phim = 17.0 ** (-0.25)
        assert abs(out['phim'] - expected_phim) < 1e-10, \
            f"phim(-1): got {out['phim']}, expected {expected_phim}"

    def test_unstable_phic(self):
        """At zeta=-1: phic = (1+16)^(-1/2) = 1/sqrt(17)."""
        out = run_fortran(EXE, {'zeta': -1.0})
        expected_phic = 1.0 / math.sqrt(17.0)
        assert abs(out['phic'] - expected_phic) < 1e-10, \
            f"phic(-1): got {out['phic']}, expected {expected_phic}"

    # --- Property tests ---

    def test_phim_greater_than_phic_unstable(self):
        """For zeta < 0: phim > phic (momentum less suppressed than scalars
        in unstable conditions: (1-16z)^(-1/4) > (1-16z)^(-1/2) for z<0)."""
        for zeta in [-0.5, -1.0, -2.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert out['phim'] > out['phic'], \
                f"zeta={zeta}: phim={out['phim']} should exceed phic={out['phic']}"

    def test_phim_equals_phic_stable(self):
        """For zeta >= 0, both phi functions use the same linear formula."""
        for zeta in [0.0, 0.5, 1.0, 2.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert abs(out['phim'] - out['phic']) < 1e-12, \
                f"zeta={zeta}: phim={out['phim']} should equal phic={out['phic']}"

    def test_phi_increases_stably(self):
        """phi increases monotonically with zeta in stable range."""
        zetas = [0.0, 0.2, 0.5, 1.0, 2.0]
        phim_vals = [run_fortran(EXE, {'zeta': z})['phim'] for z in zetas]
        for i in range(len(phim_vals) - 1):
            assert phim_vals[i] < phim_vals[i + 1], \
                f"phim not increasing: phi({zetas[i]})={phim_vals[i]} >= phi({zetas[i+1]})={phim_vals[i+1]}"

    def test_phi_greater_than_zero(self):
        """phi functions must be positive."""
        for zeta in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            out = run_fortran(EXE, {'zeta': zeta})
            assert out['phim'] > 0.0
            assert out['phic'] > 0.0
