"""
Functional tests for Monin-Obukhov stability functions in
MLCanopyTurbulenceMod.

Functions tested (made public for testing):
  phim_monin_obukhov(zeta)   →  phi_m(zeta)  momentum stability correction
  phic_monin_obukhov(zeta)   →  phi_c(zeta)  scalar  stability correction
  psim_monin_obukhov(zeta)   →  psi_m(zeta)  momentum integrated correction
  psic_monin_obukhov(zeta)   →  psi_c(zeta)  scalar  integrated correction

Formulae (Bonan et al. 2018, eq. A10-A13):

  phi_m(zeta):
    stable   (zeta >= 0):  1 + 5*zeta
    unstable (zeta <  0):  (1 - 16*zeta)^(-1/4)

  phi_c(zeta):
    stable   (zeta >= 0):  1 + 5*zeta
    unstable (zeta <  0):  (1 - 16*zeta)^(-1/2)

  psi_m(zeta):
    stable   (zeta >= 0):  -5*zeta
    unstable (zeta <  0):  2*ln((1+x)/2) + ln((1+x^2)/2) - 2*atan(x) + pi/2
                           where x = (1-16*zeta)^(1/4)

  psi_c(zeta):
    stable   (zeta >= 0):  -5*zeta
    unstable (zeta <  0):  2*ln((1+x^2)/2)
                           where x = (1-16*zeta)^(1/4)

Reference values:
  zeta=0:  phi_m=1, phi_c=1, psi_m=0, psi_c=0  (neutral)
  zeta=1:  phi_m=6, phi_c=6, psi_m=-5, psi_c=-5 (stable, linear)
"""

import math
import pytest
from utils import run_fortran


# ──────────────────────────────────────────────────────────────────────────────
# phim / phic
# ──────────────────────────────────────────────────────────────────────────────

class TestPhiMO:
    """Tests for phim_monin_obukhov and phic_monin_obukhov."""

    # --- Reference value tests (neutral and stable — analytic formulas) ---

    def test_neutral(self):
        """At zeta=0 (neutral), phim=1 and phic=1."""
        out = run_fortran('test_phim_mo.exe', {'zeta': 0.0})
        assert abs(out['phim'] - 1.0) < 1e-12
        assert abs(out['phic'] - 1.0) < 1e-12

    def test_stable_zeta1(self):
        """At zeta=1, phim = phic = 1 + 5*1 = 6."""
        out = run_fortran('test_phim_mo.exe', {'zeta': 1.0})
        assert abs(out['phim'] - 6.0) < 1e-12
        assert abs(out['phic'] - 6.0) < 1e-12

    def test_stable_zeta_half(self):
        """At zeta=0.5, phim = phic = 1 + 5*0.5 = 3.5."""
        out = run_fortran('test_phim_mo.exe', {'zeta': 0.5})
        assert abs(out['phim'] - 3.5) < 1e-12
        assert abs(out['phic'] - 3.5) < 1e-12

    def test_unstable_phim(self):
        """At zeta=-1: phim = (1+16)^(-1/4) = 17^(-0.25)."""
        out = run_fortran('test_phim_mo.exe', {'zeta': -1.0})
        expected_phim = 17.0 ** (-0.25)
        assert abs(out['phim'] - expected_phim) < 1e-10, \
            f"phim(-1): got {out['phim']}, expected {expected_phim}"

    def test_unstable_phic(self):
        """At zeta=-1: phic = (1+16)^(-1/2) = 1/sqrt(17)."""
        out = run_fortran('test_phim_mo.exe', {'zeta': -1.0})
        expected_phic = 1.0 / math.sqrt(17.0)
        assert abs(out['phic'] - expected_phic) < 1e-10, \
            f"phic(-1): got {out['phic']}, expected {expected_phic}"

    # --- Property tests ---

    def test_phim_greater_than_phic_unstable(self):
        """For zeta < 0: phim > phic (momentum less suppressed than scalars
        in unstable conditions: (1-16z)^(-1/4) > (1-16z)^(-1/2) for z<0)."""
        for zeta in [-0.5, -1.0, -2.0]:
            out = run_fortran('test_phim_mo.exe', {'zeta': zeta})
            assert out['phim'] > out['phic'], \
                f"zeta={zeta}: phim={out['phim']} should exceed phic={out['phic']}"

    def test_phim_equals_phic_stable(self):
        """For zeta >= 0, both phi functions use the same linear formula."""
        for zeta in [0.0, 0.5, 1.0, 2.0]:
            out = run_fortran('test_phim_mo.exe', {'zeta': zeta})
            assert abs(out['phim'] - out['phic']) < 1e-12, \
                f"zeta={zeta}: phim={out['phim']} should equal phic={out['phic']}"

    def test_phi_increases_stably(self):
        """phi increases monotonically with zeta in stable range."""
        zetas = [0.0, 0.2, 0.5, 1.0, 2.0]
        phim_vals = [run_fortran('test_phim_mo.exe', {'zeta': z})['phim'] for z in zetas]
        for i in range(len(phim_vals) - 1):
            assert phim_vals[i] < phim_vals[i + 1], \
                f"phim not increasing: phi({zetas[i]})={phim_vals[i]} >= phi({zetas[i+1]})={phim_vals[i+1]}"

    def test_phi_greater_than_zero(self):
        """phi functions must be positive."""
        for zeta in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            out = run_fortran('test_phim_mo.exe', {'zeta': zeta})
            assert out['phim'] > 0.0
            assert out['phic'] > 0.0


# ──────────────────────────────────────────────────────────────────────────────
# psim / psic
# ──────────────────────────────────────────────────────────────────────────────

class TestPsiMO:
    """Tests for psim_monin_obukhov and psic_monin_obukhov."""

    # --- Reference value tests (neutral and stable — analytic formulas) ---

    def test_neutral(self):
        """At zeta=0 (neutral): psim=0, psic=0 (no stability correction)."""
        out = run_fortran('test_psim_mo.exe', {'zeta': 0.0})
        assert abs(out['psim']) < 1e-12
        assert abs(out['psic']) < 1e-12

    def test_stable_zeta1(self):
        """At zeta=1 (stable): psim = psic = -5."""
        out = run_fortran('test_psim_mo.exe', {'zeta': 1.0})
        assert abs(out['psim'] - (-5.0)) < 1e-12
        assert abs(out['psic'] - (-5.0)) < 1e-12

    def test_stable_zeta_half(self):
        """At zeta=0.5 (stable): psim = psic = -2.5."""
        out = run_fortran('test_psim_mo.exe', {'zeta': 0.5})
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
        out = run_fortran('test_psim_mo.exe', {'zeta': -1.0})
        assert abs(out['psim'] - expected) < 1e-10, \
            f"psim(-1): got {out['psim']}, expected {expected}"

    def test_unstable_psic_reference(self):
        """At zeta=-1: psic = 2*ln((1+x^2)/2), x = 17^(1/4)."""
        x = 17.0 ** 0.25
        expected = 2.0 * math.log((1.0 + x * x) / 2.0)
        out = run_fortran('test_psim_mo.exe', {'zeta': -1.0})
        assert abs(out['psic'] - expected) < 1e-10, \
            f"psic(-1): got {out['psic']}, expected {expected}"

    # --- Property tests ---

    def test_psi_negative_stable(self):
        """In stable conditions, both psi functions are negative."""
        for zeta in [0.1, 0.5, 1.0, 2.0]:
            out = run_fortran('test_psim_mo.exe', {'zeta': zeta})
            assert out['psim'] < 0.0, f"psim({zeta}) = {out['psim']} should be < 0"
            assert out['psic'] < 0.0, f"psic({zeta}) = {out['psic']} should be < 0"

    def test_psi_positive_unstable(self):
        """In unstable conditions, both psi functions are positive."""
        for zeta in [-0.1, -0.5, -1.0]:
            out = run_fortran('test_psim_mo.exe', {'zeta': zeta})
            assert out['psim'] > 0.0, f"psim({zeta}) = {out['psim']} should be > 0"
            assert out['psic'] > 0.0, f"psic({zeta}) = {out['psic']} should be > 0"

    def test_psim_equals_psic_stable(self):
        """For zeta >= 0, both use the same formula: psi = -5*zeta."""
        for zeta in [0.0, 0.3, 1.0]:
            out = run_fortran('test_psim_mo.exe', {'zeta': zeta})
            assert abs(out['psim'] - out['psic']) < 1e-12

    def test_psim_greater_than_psic_unstable(self):
        """For zeta < 0: psim > psic (greater correction for momentum)."""
        for zeta in [-0.5, -1.0, -2.0]:
            out = run_fortran('test_psim_mo.exe', {'zeta': zeta})
            assert out['psim'] > out['psic'], \
                f"zeta={zeta}: psim={out['psim']} should exceed psic={out['psic']}"
