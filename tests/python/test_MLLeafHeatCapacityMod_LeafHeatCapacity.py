"""
Functional tests for MLLeafHeatCapacityMod :: CalcLeafHeatCapacity.

Computes leaf heat capacity (J/K/m2 leaf) from specific leaf area sla (m2/gC).
See Bonan et al. (2018) Geosci. Model Dev. eq. (A29).

Formula (using module constants from MLclm_varcon and clm_varcon):
  lma          = (1/sla) * 0.001             [kg C/m2]  (sla in m2/gC)
  dry_weight   = lma / fcarbon               [kg DM/m2] (fcarbon = 0.5)
  fresh_weight = dry_weight / (1 - fwater)   [kg FM/m2] (fwater = 0.7)
  leaf_water   = fwater * fresh_weight       [kg H2O/m2]
  cpleaf       = cpbio * dry_weight + cpliq * leaf_water  [J/K/m2]

Module constant defaults:
  cpbio   = 4188 / 3  ≈ 1396.0  J/kg/K  (specific heat dry biomass)
  fcarbon = 0.5                          (C fraction of dry mass)
  fwater  = 0.7                          (water fraction of fresh mass)
  cpliq   = 4188.0    J/kg/K            (from clm_varcon, specific heat water)

Simplified closed form:
  cpleaf = (1/sla) * 0.001 / 0.5 * (cpbio + fwater/(1-fwater) * cpliq)
         = 0.002/sla * (1396 + (7/3)*4188)
         ≈ 22.336 / sla   [J/K/m2 leaf]

Physical interpretation:
  Higher sla (thinner leaves) → lower cpleaf (less thermal mass per m2).
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLLeafHeatCapacityMod_LeafHeatCapacity.exe'

# Module constants
CPBIO   = 4188.0 / 3.0   # J/kg/K  (MLclm_varcon :: cpbio)
FCARBON = 0.5             # kg C / kg DM  (MLclm_varcon :: fcarbon)
FWATER  = 0.7             # kg H2O / kg FM  (MLclm_varcon :: fwater)
CPLIQ   = 4188.0          # J/kg/K  (clm_varcon :: cpliq)


def _cpleaf_ref(sla):
    """Python reference for CalcLeafHeatCapacity."""
    lma          = (1.0 / sla) * 0.001          # kg C/m2
    dry_weight   = lma / FCARBON                # kg DM/m2
    fresh_weight = dry_weight / (1.0 - FWATER)  # kg FM/m2
    leaf_water   = FWATER * fresh_weight        # kg H2O/m2
    return CPBIO * dry_weight + CPLIQ * leaf_water


# ──────────────────────────────────────────────────────────────────────────────
# Reference value tests
# ──────────────────────────────────────────────────────────────────────────────

class TestReferenceValues:

    def test_typical_broadleaf_sla(self):
        """sla=0.04 m2/gC (broadleaf): cpleaf matches Python reference."""
        sla = 0.04
        out = run_fortran(EXE, {'sla': sla})
        expected = _cpleaf_ref(sla)
        assert abs(out['cpleaf'] - expected) < 1e-6, \
            f"cpleaf={out['cpleaf']}, expected={expected}"

    def test_needleleaf_sla(self):
        """sla=0.02 m2/gC (needleleaf, high LMA): cpleaf matches reference."""
        sla = 0.02
        out = run_fortran(EXE, {'sla': sla})
        expected = _cpleaf_ref(sla)
        assert abs(out['cpleaf'] - expected) < 1e-6, \
            f"cpleaf={out['cpleaf']}, expected={expected}"

    def test_grass_sla(self):
        """sla=0.06 m2/gC (grass, low LMA): cpleaf matches reference."""
        sla = 0.06
        out = run_fortran(EXE, {'sla': sla})
        expected = _cpleaf_ref(sla)
        assert abs(out['cpleaf'] - expected) < 1e-6, \
            f"cpleaf={out['cpleaf']}, expected={expected}"

    def test_known_approximate_value(self):
        """cpleaf ≈ 22.336 / sla for default constants."""
        # Derive: 0.002/sla * (cpbio + fwater/(1-fwater)*cpliq)
        coeff = 0.002 * (CPBIO + FWATER / (1.0 - FWATER) * CPLIQ)
        for sla in [0.02, 0.04, 0.06, 0.10]:
            out = run_fortran(EXE, {'sla': sla})
            expected = coeff / sla
            assert abs(out['cpleaf'] - expected) < 1e-6, \
                f"sla={sla}: cpleaf={out['cpleaf']}, expected ≈ {expected}"


# ──────────────────────────────────────────────────────────────────────────────
# Property tests
# ──────────────────────────────────────────────────────────────────────────────

class TestProperties:

    def test_cpleaf_positive(self):
        """cpleaf > 0 for any positive sla."""
        for sla in [0.01, 0.04, 0.10]:
            out = run_fortran(EXE, {'sla': sla})
            assert out['cpleaf'] > 0.0, f"sla={sla}: cpleaf={out['cpleaf']} not positive"

    def test_cpleaf_decreases_with_sla(self):
        """Higher sla (thinner leaves) → lower cpleaf (less thermal mass per m2)."""
        slas = [0.01, 0.02, 0.04, 0.06, 0.10]
        cpleafs = [run_fortran(EXE, {'sla': s})['cpleaf'] for s in slas]
        for i in range(len(cpleafs) - 1):
            assert cpleafs[i] > cpleafs[i + 1], \
                f"cpleaf not decreasing: cpleaf({slas[i]})={cpleafs[i]}, " \
                f"cpleaf({slas[i+1]})={cpleafs[i+1]}"

    def test_cpleaf_inversely_proportional_to_sla(self):
        """cpleaf ∝ 1/sla: doubling sla halves cpleaf."""
        sla1, sla2 = 0.02, 0.04
        out1 = run_fortran(EXE, {'sla': sla1})
        out2 = run_fortran(EXE, {'sla': sla2})
        ratio = out1['cpleaf'] / out2['cpleaf']
        expected_ratio = sla2 / sla1   # = 2.0
        assert abs(ratio - expected_ratio) < 1e-10, \
            f"cpleaf ratio={ratio}, expected {expected_ratio} (inverse-sla scaling)"

    def test_water_dominates_heat_capacity(self):
        """Leaf water contributes more heat capacity than dry biomass.
        fwater/(1-fwater) * cpliq / cpbio = (7/3 * 4188) / (4188/3) = 7
        So water contribution = 7× dry-biomass contribution.
        """
        sla = 0.04
        lma = (1.0 / sla) * 0.001
        dry_weight = lma / FCARBON
        fresh_weight = dry_weight / (1.0 - FWATER)
        leaf_water = FWATER * fresh_weight

        cp_dry   = CPBIO  * dry_weight
        cp_water = CPLIQ  * leaf_water

        assert cp_water > cp_dry, \
            f"Water heat capacity ({cp_water:.2f}) should exceed dry-mass ({cp_dry:.2f})"
        assert abs(cp_water / cp_dry - 7.0) < 1e-10, \
            f"Expected ratio cp_water/cp_dry = 7, got {cp_water/cp_dry}"
