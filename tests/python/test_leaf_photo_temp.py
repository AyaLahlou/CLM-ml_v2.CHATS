"""
Functional tests for the photosynthesis temperature response functions in
MLLeafPhotosynthesisMod.

Functions tested (made public for testing):
  ft(tl, ha)          → Arrhenius temperature response, =1 at 25°C
  fth(tl, hd, se, c)  → high-temperature inhibition
  fth25(hd, se)       → scaling factor: fth denominator evaluated at 25°C

Physical background
-------------------
ft(tl, ha) = exp( ha / (rgas * 298.15) * (1 - 298.15 / tl) )

  - Normalised to 1.0 at tl = 298.15 K (25°C) for any ha
  - Monotonically increasing for ha > 0 (activation energy)

fth25(hd, se) = 1 + exp( (-hd + se * 298.15) / (rgas * 298.15) )

  - Always > 1 for typical CLM parameter values

fth(tl, hd, se, c) = c / ( 1 + exp( (-hd + se*tl) / (rgas*tl) ) )

  - With c = fth25(hd, se): normalised so fth(298.15) = 1.0
  - Decreases sharply at high temperature (inhibition)

Typical CLM values used:
  rgas  = 8.31446 J/K/mol
  tfrz  = 273.15 K
  kcha  = 79430   J/mol (kc activation energy)
  vcmaxha_noacclim = 65330 J/mol
  jmaxha_noacclim  = 43540 J/mol
  vcmaxhd_noacclim = 150000 J/mol
  vcmaxse_noacclim = 490   J/mol/K
"""

import math
import pytest
from utils import run_fortran


# Physical constants (from MLclm_varcon / clm_varcon)
RGAS = 8.31446   # J/K/mol
TFRZ = 273.15    # K
T25  = TFRZ + 25.0  # 298.15 K

# Standard CLM parameter sets
VC_HA = 65330.0    # vcmax activation energy (J/mol)
JM_HA = 43540.0    # jmax activation energy (J/mol)
KC_HA = 79430.0    # kc activation energy (J/mol)
HD    = 150000.0   # deactivation energy (J/mol)
SE    = 490.0      # entropy term (J/mol/K)


# ──────────────────────────────────────────────────────────────────────────────
# ft
# ──────────────────────────────────────────────────────────────────────────────

class TestFt:
    """Tests for MLLeafPhotosynthesisMod :: ft(tl, ha) → Arrhenius response."""

    # --- Reference value tests ---

    def test_at_25C_any_ha(self):
        """ft(298.15, ha) == 1.0 for any activation energy (normalisation identity)."""
        for ha in [VC_HA, JM_HA, KC_HA, 30000.0, 100000.0]:
            out = run_fortran('test_ft.exe', {'tl': T25, 'ha': ha})
            assert abs(out['ans'] - 1.0) < 1e-10, \
                f"ft(298.15, {ha}) = {out['ans']}, expected 1.0"

    def test_vcmax_at_35C(self):
        """ft(308.15, vcmaxha) against hand-calculation of Arrhenius formula.
        ft = exp(65330 / (8.31446*298.15) * (1 - 298.15/308.15))
        """
        tl = T25 + 10.0   # 35°C
        ha = VC_HA
        expected = math.exp(ha / (RGAS * T25) * (1.0 - T25 / tl))
        out = run_fortran('test_ft.exe', {'tl': tl, 'ha': ha})
        assert abs(out['ans'] - expected) < 1e-10, \
            f"ft(308.15, {ha}): got {out['ans']}, expected {expected}"

    def test_above_one_for_high_temp(self):
        """ft > 1 for T > 25°C when ha > 0 (Arrhenius response increases with T)."""
        out = run_fortran('test_ft.exe', {'tl': 310.0, 'ha': VC_HA})
        assert out['ans'] > 1.0

    def test_below_one_for_low_temp(self):
        """ft < 1 for T < 25°C when ha > 0."""
        out = run_fortran('test_ft.exe', {'tl': 285.0, 'ha': VC_HA})
        assert out['ans'] < 1.0

    # --- Property tests ---

    def test_monotone_increasing_with_temperature(self):
        """ft increases monotonically with temperature for ha > 0."""
        temps = [280.0, 290.0, 298.15, 308.0, 318.0]
        vals  = [run_fortran('test_ft.exe', {'tl': t, 'ha': VC_HA})['ans']
                 for t in temps]
        for i in range(len(vals) - 1):
            assert vals[i] < vals[i + 1], \
                f"ft not increasing: ft({temps[i]})={vals[i]} >= ft({temps[i+1]})={vals[i+1]}"

    def test_zero_ha_gives_one(self):
        """With ha=0 the Arrhenius exponent is 0, so ft=1 at all temperatures."""
        for tl in [280.0, 300.0, 320.0]:
            out = run_fortran('test_ft.exe', {'tl': tl, 'ha': 0.0})
            assert abs(out['ans'] - 1.0) < 1e-10


# ──────────────────────────────────────────────────────────────────────────────
# fth25
# ──────────────────────────────────────────────────────────────────────────────

class TestFth25:
    """Tests for MLLeafPhotosynthesisMod :: fth25(hd, se) → scaling factor."""

    def _expected(self, hd, se):
        return 1.0 + math.exp((-hd + se * T25) / (RGAS * T25))

    # --- Reference value tests ---

    def test_vcmax_standard_params(self):
        """fth25 with standard vcmax parameters (hd=150000, se=490)."""
        out = run_fortran('test_fth25.exe', {'hd': HD, 'se': SE})
        expected = self._expected(HD, SE)
        assert abs(out['ans'] - expected) < 1e-10, \
            f"fth25({HD},{SE}): got {out['ans']}, expected {expected}"

    # --- Property tests ---

    def test_always_greater_than_one(self):
        """fth25 = 1 + exp(...) >= 1 always."""
        for hd, se in [(150000.0, 490.0), (200000.0, 490.0), (150000.0, 550.0)]:
            out = run_fortran('test_fth25.exe', {'hd': hd, 'se': se})
            assert out['ans'] > 1.0, f"fth25({hd},{se}) = {out['ans']} should be > 1"

    def test_increases_with_se(self):
        """Larger se → larger (-hd + se*T25) → larger exponent → larger fth25."""
        out_low  = run_fortran('test_fth25.exe', {'hd': HD, 'se': 450.0})
        out_high = run_fortran('test_fth25.exe', {'hd': HD, 'se': 530.0})
        assert out_high['ans'] > out_low['ans']

    def test_decreases_with_hd(self):
        """Larger hd → more negative (-hd + se*T25) → smaller exponent → smaller fth25."""
        out_low  = run_fortran('test_fth25.exe', {'hd': 140000.0, 'se': SE})
        out_high = run_fortran('test_fth25.exe', {'hd': 160000.0, 'se': SE})
        assert out_high['ans'] < out_low['ans']


# ──────────────────────────────────────────────────────────────────────────────
# fth
# ──────────────────────────────────────────────────────────────────────────────

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
        out = run_fortran('test_fth.exe', {'tl': T25, 'hd': HD, 'se': SE, 'c': c})
        assert abs(out['ans'] - 1.0) < 1e-10, \
            f"fth(298.15, {HD}, {SE}, fth25) = {out['ans']}, expected 1.0"

    def test_reference_at_35C(self):
        """Hand-computed fth at 35°C with c=1 (no normalisation)."""
        tl = T25 + 10.0  # 308.15 K
        c = 1.0
        expected = self._expected(tl, HD, SE, c)
        out = run_fortran('test_fth.exe', {'tl': tl, 'hd': HD, 'se': SE, 'c': c})
        assert abs(out['ans'] - expected) < 1e-10, \
            f"fth(308.15): got {out['ans']}, expected {expected}"

    # --- Property tests ---

    def test_decreases_at_high_temperature(self):
        """fth (normalised) decreases above ~40-45°C (high-T inhibition)."""
        c = self._fth25_python(HD, SE)
        fth_25  = run_fortran('test_fth.exe', {'tl': T25,  'hd': HD, 'se': SE, 'c': c})['ans']
        fth_40  = run_fortran('test_fth.exe', {'tl': 313.15, 'hd': HD, 'se': SE, 'c': c})['ans']
        fth_50  = run_fortran('test_fth.exe', {'tl': 323.15, 'hd': HD, 'se': SE, 'c': c})['ans']
        # At 25°C = 1; at higher T the denominator grows, so fth decreases
        assert fth_50 < fth_40, \
            f"fth not decreasing with T: fth(40)={fth_40}, fth(50)={fth_50}"
        assert fth_40 < fth_25 * 1.1, \
            "fth at 40°C should not exceed fth at 25°C by much"

    def test_c_scales_output_linearly(self):
        """fth is linear in c: fth(tl, hd, se, 2*c) = 2 * fth(tl, hd, se, c)."""
        tl = 310.0
        c0 = 1.0
        out1 = run_fortran('test_fth.exe', {'tl': tl, 'hd': HD, 'se': SE, 'c': c0})
        out2 = run_fortran('test_fth.exe', {'tl': tl, 'hd': HD, 'se': SE, 'c': 2.0 * c0})
        assert abs(out2['ans'] - 2.0 * out1['ans']) < 1e-10
