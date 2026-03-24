"""
Functional tests for MLLeafPhotosynthesisMod :: fth25.

Scaling factor for photosynthesis temperature inhibition at 25°C (298.15 K):
  fth25(hd, se) = 1 + exp( (-hd + se * 298.15) / (rgas * 298.15) )

Physical background
-------------------
  - Always > 1 for typical CLM parameter values
  - Used as the normalisation constant c in fth()

Typical CLM values:
  rgas              = 8.31446 J/K/mol
  vcmaxhd_noacclim  = 150000 J/mol
  vcmaxse_noacclim  = 490    J/mol/K
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLLeafPhotosynthesisMod_fth25.exe'

RGAS = 8.31446
TFRZ = 273.15
T25  = TFRZ + 25.0   # 298.15 K
HD   = 150000.0
SE   = 490.0


class TestFth25:
    """Tests for MLLeafPhotosynthesisMod :: fth25(hd, se) → scaling factor."""

    def _expected(self, hd, se):
        return 1.0 + math.exp((-hd + se * T25) / (RGAS * T25))

    # --- Reference value tests ---

    def test_vcmax_standard_params(self):
        """fth25 with standard vcmax parameters (hd=150000, se=490)."""
        out = run_fortran(EXE, {'hd': HD, 'se': SE})
        expected = self._expected(HD, SE)
        assert abs(out['ans'] - expected) < 1e-10, \
            f"fth25({HD},{SE}): got {out['ans']}, expected {expected}"

    # --- Property tests ---

    def test_always_greater_than_one(self):
        """fth25 = 1 + exp(...) >= 1 always."""
        for hd, se in [(150000.0, 490.0), (200000.0, 490.0), (150000.0, 550.0)]:
            out = run_fortran(EXE, {'hd': hd, 'se': se})
            assert out['ans'] > 1.0, f"fth25({hd},{se}) = {out['ans']} should be > 1"

    def test_increases_with_se(self):
        """Larger se → larger (-hd + se*T25) → larger exponent → larger fth25."""
        out_low  = run_fortran(EXE, {'hd': HD, 'se': 450.0})
        out_high = run_fortran(EXE, {'hd': HD, 'se': 530.0})
        assert out_high['ans'] > out_low['ans']

    def test_decreases_with_hd(self):
        """Larger hd → more negative (-hd + se*T25) → smaller exponent → smaller fth25."""
        out_low  = run_fortran(EXE, {'hd': 140000.0, 'se': SE})
        out_high = run_fortran(EXE, {'hd': 160000.0, 'se': SE})
        assert out_high['ans'] < out_low['ans']
