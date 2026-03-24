"""
Functional tests for MLLeafPhotosynthesisMod :: RealizedRate.

Computes gross photosynthesis agross from three rate-limiting terms:
  ac  - Rubisco-limited rate
  aj  - RuBP regeneration-limited rate
  ap  - Product-limited (C3) or CO2-limited (C4) rate

Two colimitation modes (colim_type in MLclm_varctl):
  0: simple minimum   – agross = min(ac, aj) for C3
  1: smooth quadratic – uses curvature parameter theta = 0.98 (C3) or 0.80 (C4)

Physical background
-------------------
For colim_type=1 (co-limitation), the C3 smooth minimum is the smaller root of:
  theta * x^2 - (ac + aj) * x + ac * aj = 0   (theta = colim_c3a = 0.98)

Key properties:
  - agross(colim_type=1) ≤ min(ac, aj) for C3  (smoothing never exceeds minimum)
  - agross → min(ac, aj) when one rate ≪ the other
  - agross > 0 for positive inputs
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLLeafPhotosynthesisMod_RealizedRate.exe'

# CLM default co-limitation parameters (from MLclm_varcon)
COLIM_C3A = 0.98   # C3 curvature for Ac-Aj co-limitation
COLIM_C4A = 0.80   # C4 curvature for Ac-Aj co-limitation
COLIM_C4B = 0.95   # C4 curvature for Ai-Ap co-limitation


def _colim_smooth(theta, ac, aj):
    """Python reference for co-limited minimum (smaller quadratic root)."""
    b = -(ac + aj)
    c = ac * aj
    disc = b * b - 4.0 * theta * c
    return (-b - math.sqrt(disc)) / (2.0 * theta)


# ──────────────────────────────────────────────────────────────────────────────
# colim_type = 0  (simple minimum)
# ──────────────────────────────────────────────────────────────────────────────

class TestRealizedRateMinimum:
    """Tests for colim_type=0: agross = min(ac, aj) for C3."""

    def test_c3_ac_limits(self):
        """C3, ac < aj: agross = ac."""
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': 5.0, 'aj': 15.0,
                                'ap': 20.0, 'colim_type_in': 0})
        assert abs(out['agross'] - 5.0) < 1e-12

    def test_c3_aj_limits(self):
        """C3, aj < ac: agross = aj."""
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': 15.0, 'aj': 8.0,
                                'ap': 20.0, 'colim_type_in': 0})
        assert abs(out['agross'] - 8.0) < 1e-12

    def test_c3_equal_rates(self):
        """C3, ac = aj: agross = that value."""
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': 10.0, 'aj': 10.0,
                                'ap': 20.0, 'colim_type_in': 0})
        assert abs(out['agross'] - 10.0) < 1e-12

    def test_c4_ap_limits(self):
        """C4: agross = min(ac, aj, ap) when ap is smallest."""
        out = run_fortran(EXE, {'c3psn': 0.0, 'ac': 10.0, 'aj': 12.0,
                                'ap': 4.0, 'colim_type_in': 0})
        assert abs(out['agross'] - 4.0) < 1e-12

    def test_c4_ac_limits(self):
        """C4: agross = ac when ac is the minimum."""
        out = run_fortran(EXE, {'c3psn': 0.0, 'ac': 3.0, 'aj': 12.0,
                                'ap': 8.0, 'colim_type_in': 0})
        assert abs(out['agross'] - 3.0) < 1e-12


# ──────────────────────────────────────────────────────────────────────────────
# colim_type = 1  (smooth co-limitation)
# ──────────────────────────────────────────────────────────────────────────────

class TestRealizedRateColimited:
    """Tests for colim_type=1: smooth quadratic co-limitation."""

    def test_c3_equal_rates_reference(self):
        """C3, ac=aj=10: co-limited agross from quadratic formula.
        0.98*x^2 - 20*x + 100 = 0  →  smaller root ≈ 8.761
        """
        ac = aj = 10.0
        expected = _colim_smooth(COLIM_C3A, ac, aj)
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                'ap': 20.0, 'colim_type_in': 1})
        assert abs(out['agross'] - expected) < 1e-10, \
            f"agross={out['agross']}, expected={expected}"

    def test_c3_asymmetric_reference(self):
        """C3, ac=10, aj=5: co-limited value from quadratic.
        0.98*x^2 - 15*x + 50 = 0  →  smaller root ≈ 4.906
        """
        ac, aj = 10.0, 5.0
        expected = _colim_smooth(COLIM_C3A, ac, aj)
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                'ap': 20.0, 'colim_type_in': 1})
        assert abs(out['agross'] - expected) < 1e-10, \
            f"agross={out['agross']}, expected={expected}"

    # --- Property tests ---

    def test_c3_colim_leq_minimum(self):
        """C3: co-limited agross ≤ min(ac, aj) always (smooth min < true min)."""
        cases = [(10.0, 5.0), (5.0, 10.0), (10.0, 10.0), (2.0, 20.0)]
        for ac, aj in cases:
            out = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                    'ap': 20.0, 'colim_type_in': 1})
            assert out['agross'] <= min(ac, aj) + 1e-12, \
                f"ac={ac}, aj={aj}: agross={out['agross']} > min={min(ac,aj)}"

    def test_c3_approaches_min_when_dominated(self):
        """C3: when one rate << the other, co-limited ≈ true minimum."""
        ac, aj = 1.0, 1000.0
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                'ap': 9999.0, 'colim_type_in': 1})
        assert abs(out['agross'] - ac) / ac < 1e-4, \
            f"Expected agross ≈ {ac}, got {out['agross']}"

    def test_c3_positive_for_positive_inputs(self):
        """C3: agross > 0 when ac, aj > 0."""
        out = run_fortran(EXE, {'c3psn': 1.0, 'ac': 5.0, 'aj': 8.0,
                                'ap': 20.0, 'colim_type_in': 1})
        assert out['agross'] > 0.0

    def test_colim_vs_minimum_ordering(self):
        """co-limited rate < minimum rate for equal inputs (smooth min effect)."""
        ac = aj = 12.0
        out_min   = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                      'ap': 20.0, 'colim_type_in': 0})
        out_colim = run_fortran(EXE, {'c3psn': 1.0, 'ac': ac, 'aj': aj,
                                      'ap': 20.0, 'colim_type_in': 1})
        assert out_colim['agross'] < out_min['agross'], \
            "Co-limited rate should be strictly less than minimum when ac=aj"
