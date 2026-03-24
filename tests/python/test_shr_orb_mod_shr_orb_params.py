"""
Functional tests for shr_orb_mod :: shr_orb_params.

Computes orbital parameters for a given calendar year using secular
variations in orbital elements (Berger 1978 / CESM orbit module).

For modern epoch (~2000 AD):
  eccen  ≈ 0.0167  (eccentricity, very small → nearly circular)
  obliq  ≈ 23.44°  (obliquity)
  obliqr = obliq in radians ≈ 0.4091 rad
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_shr_orb_mod_shr_orb_params.exe'


class TestShrOrbParams:
    """Tests for shr_orb_mod :: shr_orb_params."""

    # --- Reference value tests for year 2000 ---

    def test_eccentricity_modern(self):
        """Eccentricity for 2000 AD should be ~0.0167 (well-known value)."""
        out = run_fortran(EXE, {'iyear_AD': 2000})
        assert 0.015 < out['eccen'] < 0.020, \
            f"eccen={out['eccen']}, expected ~0.0167"

    def test_obliquity_modern(self):
        """Obliquity for 2000 AD should be ~23.44° (~0.409 rad)."""
        out = run_fortran(EXE, {'iyear_AD': 2000})
        assert 23.0 < out['obliq'] < 24.0, \
            f"obliq={out['obliq']}°, expected ~23.44°"
        # Check radians consistency
        assert abs(out['obliqr'] - math.radians(out['obliq'])) < 1e-6

    def test_obliqr_mvelpp_relationship(self):
        """mvelpp = mvelp in radians + pi."""
        out = run_fortran(EXE, {'iyear_AD': 2000})
        mvelp_rad = math.radians(out['mvelp'])
        # mvelpp wraps: abs(mvelpp - (mvelp_rad + pi)) should be 0 or 2*pi
        diff = abs(out['mvelpp'] - (mvelp_rad + math.pi))
        diff = min(diff, abs(diff - 2.0 * math.pi))
        assert diff < 1e-6, f"mvelpp relationship: diff={diff}"

    # --- Property tests ---

    def test_eccentricity_range(self):
        """Eccentricity is always in [0, 1)."""
        for year in [1950, 2000, 2050, 2007]:
            out = run_fortran(EXE, {'iyear_AD': year})
            assert 0.0 <= out['eccen'] < 0.1, \
                f"eccen out of range: {out['eccen']} for year {year}"

    def test_obliquity_plausible(self):
        """Obliquity varies slowly; modern epoch: 22°-25°."""
        for year in [1950, 2000, 2050]:
            out = run_fortran(EXE, {'iyear_AD': year})
            assert 22.0 < out['obliq'] < 25.0, \
                f"obliq={out['obliq']}° for year {year}"

    def test_obliqr_positive(self):
        """Obliquity in radians must be positive."""
        out = run_fortran(EXE, {'iyear_AD': 2007})
        assert out['obliqr'] > 0.0
