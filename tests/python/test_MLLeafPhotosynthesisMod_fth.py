"""
Functional tests for MLLeafPhotosynthesisMod :: fth.

High-temperature inhibition of photosynthesis:
  fth(tl, hd, se, c) = c / (1 + exp((-hd + se*tl) / (rgas*tl)))

Physical background
-------------------
  With c = fth25(hd, se): normalised so fth(298.15) = 1.0
  Decreases sharply at high temperature (inhibition)

Typical CLM values:
  rgas              = 8.31446 J/K/mol
  vcmaxhd_noacclim  = 150000 J/mol
  vcmaxse_noacclim  = 490    J/mol/K
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLLeafPhotosynthesisMod_fth.exe'

RGAS = 8.31446
TFRZ = 273.15
T25  = TFRZ + 25.0   # 298.15 K
HD   = 150000.0
SE   = 490.0


class TestFth:
    """Tests for MLLeafPhotosynthesisMod :: fth(tl, hd, se, c) → inhibition factor."""

    def _fth25_python(self, hd, se):
        return 1.0 + math.exp((-hd + se * T25) / (RGAS * T25))

    def _expected(self, tl, hd, se, c):
        return c / (1.0 + math.exp((-hd + se * tl) / (RGAS * tl)))

    # --- Reference value tests ---

    def test_at_25C_with_c_equals_fth25(self):
        """When c = fth25(hd, se), fth(298.15) = 1.0 (normalisation).
        fth(298.15) = fth25 / (1 + exp(...at 298.15)) = fth25 / fth25 = 1.
        """
        c = self._fth25_python(HD, SE)
        out = run_fortran(EXE, {'tl': T25, 'hd': HD, 'se': SE, 'c': c})
        assert abs(out['ans'] - 1.0) < 1e-10, \
            f"fth(298.15, {HD}, {SE}, fth25) = {out['ans']}, expected 1.0"

    def test_reference_at_35C(self):
        """Hand-computed fth at 35°C with c=1 (no normalisation)."""
        tl = T25 + 10.0  # 308.15 K
        c = 1.0
        expected = self._expected(tl, HD, SE, c)
        out = run_fortran(EXE, {'tl': tl, 'hd': HD, 'se': SE, 'c': c})
        assert abs(out['ans'] - expected) < 1e-10, \
            f"fth(308.15): got {out['ans']}, expected {expected}"

    # --- Property tests ---

    def test_decreases_at_high_temperature(self):
        """fth (normalised) decreases above ~40-45°C (high-T inhibition)."""
        c = self._fth25_python(HD, SE)
        fth_25  = run_fortran(EXE, {'tl': T25,    'hd': HD, 'se': SE, 'c': c})['ans']
        fth_40  = run_fortran(EXE, {'tl': 313.15, 'hd': HD, 'se': SE, 'c': c})['ans']
        fth_50  = run_fortran(EXE, {'tl': 323.15, 'hd': HD, 'se': SE, 'c': c})['ans']
        # At 25°C = 1; at higher T the denominator grows, so fth decreases
        assert fth_50 < fth_40, \
            f"fth not decreasing with T: fth(40)={fth_40}, fth(50)={fth_50}"
        assert fth_40 < fth_25 * 1.1, \
            "fth at 40°C should not exceed fth at 25°C by much"

    def test_c_scales_output_linearly(self):
        """fth is linear in c: fth(tl, hd, se, 2*c) = 2 * fth(tl, hd, se, c)."""
        tl = 310.0
        c0 = 1.0
        out1 = run_fortran(EXE, {'tl': tl, 'hd': HD, 'se': SE, 'c': c0})
        out2 = run_fortran(EXE, {'tl': tl, 'hd': HD, 'se': SE, 'c': 2.0 * c0})
        assert abs(out2['ans'] - 2.0 * out1['ans']) < 1e-10
