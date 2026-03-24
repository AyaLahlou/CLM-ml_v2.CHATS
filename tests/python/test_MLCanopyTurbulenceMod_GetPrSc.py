"""
Functional tests for MLCanopyTurbulenceMod :: GetPrSc.

Computes the effective Prandtl/Schmidt number (PrSc) at the canopy top
as a function of atmospheric stability.
See Bonan et al. (2018), eqs. (A25), (A34).

Base formula (for sparse_canopy_type = 0):
  PrSc = Pr0 + Pr1 * tanh(Pr2 * LcL)

where the default constants from MLclm_varcon are:
  Pr0 = 0.5  (neutral value)
  Pr1 = 0.3  (amplitude)
  Pr2 = 2.0  (length scale)

When sparse_canopy_type = 1, the formula is blended toward 1.0 for sparse canopies:
  PrSc = (1 - beta_neutral/beta_neutral_max) * 1 + (beta_neutral/beta_neutral_max) * PrSc_base

Physical interpretation:
  LcL = 0 (neutral):   PrSc = 0.5
  LcL → +∞ (stable):  PrSc → 0.5 + 0.3 = 0.8
  LcL → -∞ (unstable): PrSc → 0.5 - 0.3 = 0.2
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyTurbulenceMod_GetPrSc.exe'

# Default constants from MLclm_varcon
PR0 = 0.5
PR1 = 0.3
PR2 = 2.0
BETA_NEUTRAL_MAX = 0.35   # from MLclm_varcon


def _prsc_base(lcl):
    """Python reference for PrSc base formula."""
    return PR0 + PR1 * math.tanh(PR2 * lcl)


def _prsc_sparse(bn, lcl):
    """Python reference for PrSc with sparse canopy adjustment."""
    base = _prsc_base(lcl)
    w = bn / BETA_NEUTRAL_MAX
    return (1.0 - w) * 1.0 + w * base


class TestGetPrSc:
    """Tests for MLCanopyTurbulenceMod :: GetPrSc."""

    # --- Reference value tests (sparse_canopy_type = 0) ---

    def test_neutral_gives_pr0(self):
        """At LcL=0 (neutral), PrSc = Pr0 = 0.5 (sparse_canopy_type=0)."""
        out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': 0.0,
                                'sparse_canopy_type_in': 0})
        assert abs(out['PrSc'] - PR0) < 1e-12, \
            f"PrSc={out['PrSc']}, expected Pr0={PR0}"

    def test_stable_reference(self):
        """At LcL=1, PrSc = 0.5 + 0.3*tanh(2) (sparse_canopy_type=0)."""
        lcl = 1.0
        expected = _prsc_base(lcl)
        out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl,
                                'sparse_canopy_type_in': 0})
        assert abs(out['PrSc'] - expected) < 1e-12, \
            f"PrSc={out['PrSc']}, expected {expected}"

    def test_unstable_reference(self):
        """At LcL=-1, PrSc = 0.5 + 0.3*tanh(-2) (sparse_canopy_type=0)."""
        lcl = -1.0
        expected = _prsc_base(lcl)
        out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl,
                                'sparse_canopy_type_in': 0})
        assert abs(out['PrSc'] - expected) < 1e-12, \
            f"PrSc={out['PrSc']}, expected {expected}"

    # --- Property tests (sparse_canopy_type = 0) ---

    def test_increases_with_stability(self):
        """PrSc increases monotonically with LcL (more stable → larger Pr)."""
        lcls = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
        vals = [run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl,
                                  'sparse_canopy_type_in': 0})['PrSc']
                for lcl in lcls]
        for i in range(len(vals) - 1):
            assert vals[i] < vals[i + 1], \
                f"PrSc not increasing: PrSc({lcls[i]})={vals[i]} >= PrSc({lcls[i+1]})={vals[i+1]}"

    def test_bounded_below_by_pr0_minus_pr1(self):
        """PrSc > Pr0 - Pr1 = 0.2 for any LcL (asymptote of tanh)."""
        lower = PR0 - PR1
        for lcl in [-5.0, -2.0, -1.0, -0.5]:
            out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl,
                                    'sparse_canopy_type_in': 0})
            assert out['PrSc'] > lower - 1e-10, \
                f"LcL={lcl}: PrSc={out['PrSc']} below lower bound {lower}"

    def test_bounded_above_by_pr0_plus_pr1(self):
        """PrSc < Pr0 + Pr1 = 0.8 for any LcL (asymptote of tanh)."""
        upper = PR0 + PR1
        for lcl in [0.5, 1.0, 2.0, 5.0]:
            out = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': lcl,
                                    'sparse_canopy_type_in': 0})
            assert out['PrSc'] < upper + 1e-10, \
                f"LcL={lcl}: PrSc={out['PrSc']} above upper bound {upper}"

    # --- Sparse canopy adjustment (sparse_canopy_type = 1) ---

    def test_sparse_neutral_reference(self):
        """At LcL=0, sparse_canopy_type=1: PrSc blended with 1.0.
        PrSc = (1 - bn/bmax)*1 + (bn/bmax)*0.5
        For bn=0.3, bmax=0.35: PrSc = (1-6/7)*1 + (6/7)*0.5 ≈ 0.571
        """
        bn = 0.3
        expected = _prsc_sparse(bn, 0.0)
        out = run_fortran(EXE, {'beta_neutral': bn, 'LcL': 0.0,
                                'sparse_canopy_type_in': 1})
        assert abs(out['PrSc'] - expected) < 1e-12, \
            f"PrSc={out['PrSc']}, expected {expected}"

    def test_sparse_adjustment_gt_base(self):
        """sparse_canopy_type=1 gives larger PrSc than base formula at neutral."""
        out_base   = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': 0.0,
                                       'sparse_canopy_type_in': 0})
        out_sparse = run_fortran(EXE, {'beta_neutral': 0.3, 'LcL': 0.0,
                                       'sparse_canopy_type_in': 1})
        # base = 0.5; sparse blends toward 1.0, so sparse > base
        assert out_sparse['PrSc'] > out_base['PrSc'], \
            "Sparse canopy adjustment should increase PrSc at neutral"

    def test_sparse_full_canopy_equals_base(self):
        """When beta_neutral = beta_neutral_max, weight = 1: sparse adjustment
        reduces to the base formula."""
        bn = BETA_NEUTRAL_MAX  # 0.35
        lcl = 0.5
        expected = _prsc_base(lcl)
        out = run_fortran(EXE, {'beta_neutral': bn, 'LcL': lcl,
                                'sparse_canopy_type_in': 1})
        assert abs(out['PrSc'] - expected) < 1e-12, \
            f"PrSc={out['PrSc']}, expected base={expected}"
