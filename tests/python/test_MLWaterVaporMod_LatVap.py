"""
Functional tests for MLWaterVaporMod :: LatVap.

Computes the latent heat of vaporisation or sublimation (J/mol).

Physical constants used for reference values (from clm_varcon / MLclm_varcon):
  tfrz  = 273.15 K
  hvap  = 2.5010e6 J/kg   (latent heat of evaporation)
  hsub  = 2.8347e6 J/kg   (latent heat of sublimation)
  mmh2o = 18.02e-3 kg/mol (molar mass of water)
"""

import pytest
from utils import run_fortran

EXE = 'test_MLWaterVaporMod_LatVap.exe'


class TestLatVap:
    """Tests for MLWaterVaporMod :: LatVap(t) → lambda (J/mol)."""

    # Physical constants from clm_varcon / MLclm_varcon
    _hvap  = 2.5010e6    # J/kg
    _hsub  = 2.8347e6    # J/kg
    _mmh2o = 18.02e-3    # kg/mol

    def test_above_freezing_reference(self):
        """For T > 273.15, lambda = hvap * mmh2o."""
        expected = self._hvap * self._mmh2o
        out = run_fortran(EXE, {'t': 300.0})
        assert abs(out['lambda'] - expected) < 1.0, \
            f"lambda={out['lambda']}, expected={expected}"

    def test_below_freezing_reference(self):
        """For T <= 273.15, lambda = hsub * mmh2o."""
        expected = self._hsub * self._mmh2o
        out = run_fortran(EXE, {'t': 260.0})
        assert abs(out['lambda'] - expected) < 1.0, \
            f"lambda={out['lambda']}, expected={expected}"

    def test_sublimation_greater_than_evaporation(self):
        """Latent heat of sublimation > latent heat of evaporation."""
        liq = run_fortran(EXE, {'t': 280.0})['lambda']
        ice = run_fortran(EXE, {'t': 260.0})['lambda']
        assert ice > liq, \
            f"Expected lambda_ice={ice} > lambda_liq={liq}"

    def test_step_at_freezing(self):
        """LatVap jumps discontinuously at the freezing point."""
        above = run_fortran(EXE, {'t': 273.16})['lambda']
        below = run_fortran(EXE, {'t': 273.14})['lambda']
        assert below > above, \
            "lambda should be higher below freezing (sublimation > evaporation)"

    def test_constant_above_freezing(self):
        """lambda is constant for all T > tfrz (not temperature-dependent)."""
        l1 = run_fortran(EXE, {'t': 280.0})['lambda']
        l2 = run_fortran(EXE, {'t': 320.0})['lambda']
        assert abs(l1 - l2) < 1e-6

    def test_constant_below_freezing(self):
        """lambda is constant for all T <= tfrz."""
        l1 = run_fortran(EXE, {'t': 250.0})['lambda']
        l2 = run_fortran(EXE, {'t': 200.0})['lambda']
        assert abs(l1 - l2) < 1e-6
