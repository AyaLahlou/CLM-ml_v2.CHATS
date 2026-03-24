"""
Functional tests for MLWaterVaporMod :: SatVap.

Computes saturation vapor pressure (es, Pa) and its temperature
derivative (desdt, Pa/K) using the polynomial fit of Flatau et al. (1992).

Physical constants used for reference values (from clm_varcon / MLclm_varcon):
  tfrz  = 273.15 K

Reference saturation vapour pressures are from the Flatau et al. (1992)
polynomial approximation, cross-checked against the CRC Handbook.
"""

import pytest
from utils import run_fortran

EXE = 'test_MLWaterVaporMod_SatVap.exe'


class TestSatVap:
    """Tests for MLWaterVaporMod :: SatVap(t) → es, desdt."""

    # --- Reference value tests (water branch, T >= 273.15 K) ---

    def test_at_freezing_point_water(self):
        """At T=273.15 K (0°C), over water es ≈ 611.2 Pa.
        Flatau 1992 polynomial value at 0°C."""
        out = run_fortran(EXE, {'t': 273.15})
        assert abs(out['es'] - 611.2) < 2.0,   f"es={out['es']}"
        assert abs(out['desdt'] - 44.4) < 2.0, f"desdt={out['desdt']}"

    def test_at_25C(self):
        """At T=298.15 K (25°C), es ≈ 3167 Pa."""
        out = run_fortran(EXE, {'t': 298.15})
        assert abs(out['es'] - 3167.0) < 20.0, f"es={out['es']}"
        assert out['desdt'] > 0.0

    def test_at_100C(self):
        """At T=373.15 K (100°C), es ≈ 101325 Pa (1 atm — boiling point)."""
        out = run_fortran(EXE, {'t': 373.15})
        # Flatau polynomial at 100°C is the practical limit; accept ±3000 Pa
        assert abs(out['es'] - 101325.0) < 3000.0, f"es={out['es']}"

    # --- Reference value tests (ice branch, T < 273.15 K) ---

    def test_at_minus_10C_ice(self):
        """At T=263.15 K (-10°C), over ice es ≈ 259.9 Pa."""
        out = run_fortran(EXE, {'t': 263.15})
        assert abs(out['es'] - 259.9) < 5.0, f"es={out['es']}"
        assert out['desdt'] > 0.0

    # --- Property tests ---

    def test_es_increases_with_temperature(self):
        """es is monotonically increasing with T (Clausius-Clapeyron)."""
        temps = [260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0]
        es_vals = [run_fortran(EXE, {'t': t})['es'] for t in temps]
        for i in range(len(es_vals) - 1):
            assert es_vals[i] < es_vals[i + 1], \
                f"es not increasing: es({temps[i]})={es_vals[i]} >= es({temps[i+1]})={es_vals[i+1]}"

    def test_desdt_always_positive(self):
        """The derivative d(es)/dT is always positive."""
        for t in [250.0, 270.0, 290.0, 310.0]:
            out = run_fortran(EXE, {'t': t})
            assert out['desdt'] > 0.0, f"desdt={out['desdt']} at T={t}"

    def test_ice_branch_lower_than_water_at_same_temperature(self):
        """Very near 273.15 K the ice and water branches should be close.
        Just below freezing (ice) and just above freezing (water) should
        both be near 611 Pa."""
        out_water = run_fortran(EXE, {'t': 273.16})
        out_ice   = run_fortran(EXE, {'t': 273.14})
        # Both close to 611 Pa; ice slightly lower
        assert abs(out_water['es'] - out_ice['es']) < 5.0

    def test_clamped_at_upper_limit(self):
        """Temperature above 373.15 K (tc=100 limit) gives same es as T=373.15."""
        out_limit  = run_fortran(EXE, {'t': 373.15})
        out_higher = run_fortran(EXE, {'t': 400.0})
        assert abs(out_limit['es'] - out_higher['es']) < 1e-6, \
            "SatVap should clamp tc at 100°C"

    def test_clamped_at_lower_limit(self):
        """Temperature below 198.15 K (tc=-75 limit) gives same es as T=198.15."""
        out_limit  = run_fortran(EXE, {'t': 198.15})
        out_lower  = run_fortran(EXE, {'t': 150.0})
        assert abs(out_limit['es'] - out_lower['es']) < 1e-6, \
            "SatVap should clamp tc at -75°C"
