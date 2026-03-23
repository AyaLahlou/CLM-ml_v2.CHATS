"""
Functional tests for shr_orb_mod.

Subroutines tested:
  shr_orb_params(iyear_AD) → eccen, obliq, mvelp, obliqr, lambm0, mvelpp
  shr_orb_decl(calday, ...) → delta (solar declination), eccf (distance factor)

Physical background
-------------------
shr_orb_params:
  Computes orbital parameters for a given year using secular variations
  in orbital elements (Berger 1978 / CESM orbit module).
  For modern epoch (~2000 AD):
    eccen  ≈ 0.0167  (eccentricity, very small → nearly circular)
    obliq  ≈ 23.44°  (obliquity)
    obliqr = obliq in radians ≈ 0.4091 rad

shr_orb_decl:
  Computes solar declination (delta, radians) for a given calendar day.
  Key physical facts:
    - delta ≈ 0 at vernal equinox (day ~80 in NH) and autumnal equinox (~266)
    - delta ≈ +obliqr at summer solstice (day ~172, max ~23.44°)
    - delta ≈ -obliqr at winter solstice (day ~355, min ~-23.44°)
    - eccf ≥ 0 always (never negative)
    - eccf ≈ 1 + 2*eccen*cos(true_anomaly); varies slightly around 1
"""

import math
import pytest
from utils import run_fortran


# ──────────────────────────────────────────────────────────────────────────────
# shr_orb_params
# ──────────────────────────────────────────────────────────────────────────────

class TestShrOrbParams:
    """Tests for shr_orb_mod :: shr_orb_params."""

    # --- Reference value tests for year 2000 ---

    def test_eccentricity_modern(self):
        """Eccentricity for 2000 AD should be ~0.0167 (well-known value)."""
        out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2000})
        assert 0.015 < out['eccen'] < 0.020, \
            f"eccen={out['eccen']}, expected ~0.0167"

    def test_obliquity_modern(self):
        """Obliquity for 2000 AD should be ~23.44° (~0.409 rad)."""
        out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2000})
        assert 23.0 < out['obliq'] < 24.0, \
            f"obliq={out['obliq']}°, expected ~23.44°"
        # Check radians consistency
        assert abs(out['obliqr'] - math.radians(out['obliq'])) < 1e-6

    def test_obliqr_mvelpp_relationship(self):
        """mvelpp = mvelp in radians + pi."""
        out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2000})
        mvelp_rad = math.radians(out['mvelp'])
        # mvelpp wraps: abs(mvelpp - (mvelp_rad + pi)) should be 0 or 2*pi
        diff = abs(out['mvelpp'] - (mvelp_rad + math.pi))
        diff = min(diff, abs(diff - 2.0 * math.pi))
        assert diff < 1e-6, f"mvelpp relationship: diff={diff}"

    # --- Property tests ---

    def test_eccentricity_range(self):
        """Eccentricity is always in [0, 1)."""
        for year in [1950, 2000, 2050, 2007]:
            out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': year})
            assert 0.0 <= out['eccen'] < 0.1, \
                f"eccen out of range: {out['eccen']} for year {year}"

    def test_obliquity_plausible(self):
        """Obliquity varies slowly; modern epoch: 22°-25°."""
        for year in [1950, 2000, 2050]:
            out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': year})
            assert 22.0 < out['obliq'] < 25.0, \
                f"obliq={out['obliq']}° for year {year}"

    def test_obliqr_positive(self):
        """Obliquity in radians must be positive."""
        out = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2007})
        assert out['obliqr'] > 0.0


# ──────────────────────────────────────────────────────────────────────────────
# shr_orb_decl
# ──────────────────────────────────────────────────────────────────────────────

class TestShrOrbDecl:
    """Tests for shr_orb_mod :: shr_orb_decl."""

    def _run(self, iyear_AD, calday):
        return run_fortran('test_shr_orb_decl.exe',
                           {'iyear_AD': iyear_AD, 'calday': float(calday)})

    # --- Reference value tests ---

    def test_vernal_equinox_near_zero_declination(self):
        """Near the vernal equinox (day ~80), delta should be close to 0."""
        out = self._run(2007, 80.0)
        # Within ~5° of zero (0.087 rad)
        assert abs(out['delta']) < 0.10, \
            f"delta at vernal equinox: {math.degrees(out['delta']):.2f}°"

    def test_summer_solstice_max_declination(self):
        """Near summer solstice (day ~172), delta ≈ +obliqr."""
        out = self._run(2007, 172.0)
        # Should be close to max obliquity (~23.44° = 0.409 rad)
        assert 0.36 < out['delta'] < 0.43, \
            f"delta at summer solstice: {math.degrees(out['delta']):.2f}°"

    def test_winter_solstice_min_declination(self):
        """Near winter solstice (day ~355), delta ≈ -obliqr."""
        out = self._run(2007, 355.0)
        assert -0.43 < out['delta'] < -0.36, \
            f"delta at winter solstice: {math.degrees(out['delta']):.2f}°"

    def test_autumnal_equinox_near_zero(self):
        """Near the autumnal equinox (day ~266), delta ≈ 0."""
        out = self._run(2007, 266.0)
        assert abs(out['delta']) < 0.10, \
            f"delta at autumnal equinox: {math.degrees(out['delta']):.2f}°"

    def test_eccf_positive(self):
        """Earth-sun distance factor eccf must always be positive."""
        for calday in [1.0, 80.0, 172.0, 266.0, 355.0]:
            out = self._run(2007, calday)
            assert out['eccf'] > 0.0, f"eccf={out['eccf']} at calday={calday}"

    def test_eccf_near_unity(self):
        """eccf varies near 1.0 (Earth's orbit is nearly circular, e~0.017)."""
        for calday in [1.0, 80.0, 172.0, 266.0]:
            out = self._run(2007, calday)
            # Max variation: approximately 1 ± 2*eccen ≈ 1 ± 0.034
            assert 0.93 < out['eccf'] < 1.07, \
                f"eccf={out['eccf']} at calday={calday} is out of expected range"

    # --- Property tests ---

    def test_delta_bounded_by_obliquity(self):
        """Declination must lie within [-obliqr, +obliqr]."""
        out_params = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2007})
        obliqr = out_params['obliqr']
        for calday in [1.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0]:
            out = self._run(2007, calday)
            assert -obliqr - 0.001 <= out['delta'] <= obliqr + 0.001, \
                f"delta={out['delta']} rad at day {calday} exceeds obliqr={obliqr}"

    def test_declination_varies_over_year(self):
        """Declination changes sign between summer and winter."""
        summer = self._run(2007, 172.0)['delta']
        winter = self._run(2007, 355.0)['delta']
        assert summer > 0.0, "Summer solstice declination should be positive"
        assert winter < 0.0, "Winter solstice declination should be negative"

    def test_eccen_consistent_across_drivers(self):
        """Eccentricity from shr_orb_decl driver matches shr_orb_params driver."""
        out_params = run_fortran('test_shr_orb_params.exe', {'iyear_AD': 2007})
        out_decl   = self._run(2007, 80.0)
        assert abs(out_params['eccen'] - out_decl['eccen']) < 1e-12
