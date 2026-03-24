"""
Functional tests for MLLeafPhotosynthesisMod :: ft.

Arrhenius-type photosynthesis temperature response, normalised to 1.0 at 25°C:
  ft(tl, ha) = exp( ha / (rgas * 298.15) * (1 - 298.15 / tl) )

Physical background
-------------------
  - Normalised to 1.0 at tl = 298.15 K (25°C) for any ha
  - Monotonically increasing for ha > 0 (activation energy)

Typical CLM values used:
  rgas  = 8.31446 J/K/mol
  tfrz  = 273.15 K
  kcha  = 79430   J/mol (kc activation energy)
  vcmaxha_noacclim = 65330 J/mol
  jmaxha_noacclim  = 43540 J/mol
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLLeafPhotosynthesisMod_ft.exe'

RGAS = 8.31446   # J/K/mol
TFRZ = 273.15    # K
T25  = TFRZ + 25.0  # 298.15 K

VC_HA = 65330.0    # vcmax activation energy (J/mol)
JM_HA = 43540.0    # jmax activation energy (J/mol)
KC_HA = 79430.0    # kc activation energy (J/mol)


class TestFt:
    """Tests for MLLeafPhotosynthesisMod :: ft(tl, ha) → Arrhenius response."""

    # --- Reference value tests ---

    def test_at_25C_any_ha(self):
        """ft(298.15, ha) == 1.0 for any activation energy (normalisation identity)."""
        for ha in [VC_HA, JM_HA, KC_HA, 30000.0, 100000.0]:
            out = run_fortran(EXE, {'tl': T25, 'ha': ha})
            assert abs(out['ans'] - 1.0) < 1e-10, \
                f"ft(298.15, {ha}) = {out['ans']}, expected 1.0"

    def test_vcmax_at_35C(self):
        """ft(308.15, vcmaxha) against hand-calculation of Arrhenius formula.
        ft = exp(65330 / (8.31446*298.15) * (1 - 298.15/308.15))
        """
        tl = T25 + 10.0   # 35°C
        ha = VC_HA
        expected = math.exp(ha / (RGAS * T25) * (1.0 - T25 / tl))
        out = run_fortran(EXE, {'tl': tl, 'ha': ha})
        assert abs(out['ans'] - expected) < 1e-10, \
            f"ft(308.15, {ha}): got {out['ans']}, expected {expected}"

    def test_above_one_for_high_temp(self):
        """ft > 1 for T > 25°C when ha > 0 (Arrhenius response increases with T)."""
        out = run_fortran(EXE, {'tl': 310.0, 'ha': VC_HA})
        assert out['ans'] > 1.0

    def test_below_one_for_low_temp(self):
        """ft < 1 for T < 25°C when ha > 0."""
        out = run_fortran(EXE, {'tl': 285.0, 'ha': VC_HA})
        assert out['ans'] < 1.0

    # --- Property tests ---

    def test_monotone_increasing_with_temperature(self):
        """ft increases monotonically with temperature for ha > 0."""
        temps = [280.0, 290.0, 298.15, 308.0, 318.0]
        vals  = [run_fortran(EXE, {'tl': t, 'ha': VC_HA})['ans']
                 for t in temps]
        for i in range(len(vals) - 1):
            assert vals[i] < vals[i + 1], \
                f"ft not increasing: ft({temps[i]})={vals[i]} >= ft({temps[i+1]})={vals[i+1]}"

    def test_zero_ha_gives_one(self):
        """With ha=0 the Arrhenius exponent is 0, so ft=1 at all temperatures."""
        for tl in [280.0, 300.0, 320.0]:
            out = run_fortran(EXE, {'tl': tl, 'ha': 0.0})
            assert abs(out['ans'] - 1.0) < 1e-10
