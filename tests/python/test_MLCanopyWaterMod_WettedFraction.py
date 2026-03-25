"""
Functional tests for MLCanopyWaterMod :: CalcWettedFraction.

Computes wetted (fwet) and dry (fdry) fractions for a single canopy layer.
Uses module constants from MLclm_varcon (default values):
  dewmx                      = 0.1   kg H2O / m2 leaf (max interception)
  fwet_exponent              = 0.67  power-law exponent
  maximum_leaf_wetted_fraction = 0.05  cap on fwet

Formulas:
  h2ocanmx = dewmx * dpai
  fwet = clamp((h2ocan / h2ocanmx)^fwet_exponent, 0, max_fwet)
  fdry = (1 - fwet) * dlai / dpai

Special cases:
  dpai = 0:        fwet = fdry = 0
  h2ocan = 0:      fwet = 0,  fdry = dlai/dpai
  h2ocan ≥ threshold:  fwet capped at max_fwet = 0.05

Cap threshold: h2ocan_thresh = dewmx * dpai * max_fwet^(1/fwet_exponent)
               = 0.1 * dpai * 0.05^(1/0.67) ≈ 0.1 * dpai * 0.00706
So fwet is capped for nearly any measurable water loading.
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyWaterMod_WettedFraction.exe'

# Module constants from MLclm_varcon
DEWMX     = 0.1
EXPONENT  = 0.67
MAX_FWET  = 0.05


def _fwet_ref(h2ocan, dpai):
    """Python reference for fwet (before cap)."""
    if dpai <= 0.0:
        return 0.0
    h2ocanmx = DEWMX * dpai
    raw = max(h2ocan / h2ocanmx, 0.0) ** EXPONENT
    return min(raw, MAX_FWET)


def _fdry_ref(h2ocan, dpai, dlai):
    """Python reference for fdry."""
    if dpai <= 0.0:
        return 0.0
    fwet = _fwet_ref(h2ocan, dpai)
    return (1.0 - fwet) * dlai / dpai


# ──────────────────────────────────────────────────────────────────────────────
# Reference value tests
# ──────────────────────────────────────────────────────────────────────────────

class TestReferenceValues:

    def test_dry_leaf_fwet_zero(self):
        """h2ocan = 0: fwet = 0, fdry = dlai/dpai."""
        dpai, dlai = 1.0, 0.8
        out = run_fortran(EXE, {'h2ocan': 0.0, 'dpai': dpai, 'dlai': dlai})
        assert abs(out['fwet']) < 1e-15
        assert abs(out['fdry'] - dlai / dpai) < 1e-12

    def test_dry_leaf_all_plant_area_is_leaf(self):
        """h2ocan=0, dpai=dlai: fdry = 1.0."""
        out = run_fortran(EXE, {'h2ocan': 0.0, 'dpai': 0.5, 'dlai': 0.5})
        assert abs(out['fwet']) < 1e-15
        assert abs(out['fdry'] - 1.0) < 1e-12

    def test_at_capacity_fwet_capped(self):
        """h2ocan = dewmx*dpai: raw fwet = 1.0 → capped at max_fwet = 0.05."""
        dpai, dlai = 1.0, 0.8
        out = run_fortran(EXE, {'h2ocan': DEWMX * dpai, 'dpai': dpai, 'dlai': dlai})
        assert abs(out['fwet'] - MAX_FWET) < 1e-12

    def test_above_capacity_still_capped(self):
        """h2ocan >> h2ocanmx: fwet still capped at max_fwet."""
        out = run_fortran(EXE, {'h2ocan': 10.0, 'dpai': 1.0, 'dlai': 0.8})
        assert abs(out['fwet'] - MAX_FWET) < 1e-12

    def test_small_h2ocan_below_cap(self):
        """Very small h2ocan: fwet follows power-law without hitting cap."""
        h2ocan, dpai, dlai = 1e-5, 1.0, 0.8
        expected_fwet = _fwet_ref(h2ocan, dpai)
        expected_fdry = _fdry_ref(h2ocan, dpai, dlai)
        assert expected_fwet < MAX_FWET, "Test setup: should be below cap"
        out = run_fortran(EXE, {'h2ocan': h2ocan, 'dpai': dpai, 'dlai': dlai})
        assert abs(out['fwet'] - expected_fwet) < 1e-12
        assert abs(out['fdry'] - expected_fdry) < 1e-12

    def test_zero_dpai_gives_zero_fractions(self):
        """dpai = 0 (no plant area): fwet = fdry = 0."""
        out = run_fortran(EXE, {'h2ocan': 0.0, 'dpai': 0.0, 'dlai': 0.0})
        assert abs(out['fwet']) < 1e-15
        assert abs(out['fdry']) < 1e-15


# ──────────────────────────────────────────────────────────────────────────────
# Property tests
# ──────────────────────────────────────────────────────────────────────────────

class TestProperties:

    def test_fwet_in_unit_interval(self):
        """fwet ∈ [0, 1] for all inputs."""
        cases = [
            {'h2ocan': 0.0,  'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 0.05, 'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 0.5,  'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 0.0,  'dpai': 0.0, 'dlai': 0.0},
        ]
        for inp in cases:
            out = run_fortran(EXE, inp)
            assert 0.0 <= out['fwet'] <= 1.0, f"fwet={out['fwet']} out of [0,1]: {inp}"

    def test_fwet_never_exceeds_max_fwet(self):
        """fwet never exceeds maximum_leaf_wetted_fraction = 0.05."""
        cases = [
            {'h2ocan': 0.001, 'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 0.1,   'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 1.0,   'dpai': 1.0, 'dlai': 0.8},
            {'h2ocan': 100.0, 'dpai': 1.0, 'dlai': 0.8},
        ]
        for inp in cases:
            out = run_fortran(EXE, inp)
            assert out['fwet'] <= MAX_FWET + 1e-12, \
                f"fwet={out['fwet']} exceeds max_fwet={MAX_FWET}: {inp}"

    def test_fwet_monotone_increasing_with_h2ocan(self):
        """fwet increases (weakly) as h2ocan increases."""
        h2ocans = [0.0, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1.0]
        fwets = [
            run_fortran(EXE, {'h2ocan': h, 'dpai': 1.0, 'dlai': 0.8})['fwet']
            for h in h2ocans
        ]
        for i in range(len(fwets) - 1):
            assert fwets[i] <= fwets[i + 1] + 1e-14, \
                f"fwet not increasing: fwet({h2ocans[i]})={fwets[i]}, " \
                f"fwet({h2ocans[i+1]})={fwets[i+1]}"

    def test_fdry_plus_wet_fraction_le_one(self):
        """fdry + fwet*dpai/dlai ≤ 1 (dry and wet add up consistently)."""
        dpai, dlai = 1.0, 0.8
        # fdry = (1-fwet)*dlai/dpai, so fwet + fdry*dpai/dlai = fwet + (1-fwet) = 1
        for h2ocan in [0.0, 1e-5, 0.1, 1.0]:
            out = run_fortran(EXE, {'h2ocan': h2ocan, 'dpai': dpai, 'dlai': dlai})
            check = out['fwet'] + out['fdry'] * dpai / dlai
            assert abs(check - 1.0) < 1e-12, \
                f"h2ocan={h2ocan}: fwet + fdry*dpai/dlai = {check} ≠ 1"

    def test_reference_formula_matches(self):
        """Python reference formula matches Fortran for several h2ocan values."""
        dpai, dlai = 2.0, 1.5
        for h2ocan in [0.0, 1e-6, 0.001]:
            out = run_fortran(EXE, {'h2ocan': h2ocan, 'dpai': dpai, 'dlai': dlai})
            assert abs(out['fwet'] - _fwet_ref(h2ocan, dpai)) < 1e-12
            assert abs(out['fdry'] - _fdry_ref(h2ocan, dpai, dlai)) < 1e-12
