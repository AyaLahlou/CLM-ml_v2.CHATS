"""
Functional tests for MLCanopyNitrogenProfileMod :: CalcNitrogenScale.

Computes nitrogen scaling factors nscale_sun and nscale_sha for a single
canopy layer, integrating the exponential leaf nitrogen profile over
sunlit and shaded portions.

Two leaf_optics_type cases:

  leaf_optics_type = 0  (Bonan et al. 2021):
    fn       = exp(-kn*pai_above) * (1 - exp(-kn*dpai)) / kn
    fn_sun   = clump_fac/(kn + kb*clump_fac) * exp(-kn*pai_above)
               * tbi * (1 - exp(-(kn + kb*clump_fac)*dpai))
    fn_sha   = fn - fn_sun
    nscale_sun = fn_sun / (fracsun * dpai)
    nscale_sha = fn_sha / ((1 - fracsun) * dpai)

  leaf_optics_type = 1  (thin-layer approximation):
    nscale_sun = nscale_sha = fn / dpai

Key consistency property:
  nscale_sun * fracsun + nscale_sha * (1 - fracsun) = fn / dpai
  (weighted mean of nscale equals the layer-integrated nitrogen factor)
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyNitrogenProfileMod_NitrogenScale.exe'


def _fn_ref(kn, pai_above, dpai):
    """Total nitrogen integral over the layer."""
    return math.exp(-kn * pai_above) * (1.0 - math.exp(-kn * dpai)) / kn


def _nscale_ref(kn, pai_above, dpai, kb, clump_fac, fracsun, tbi, lot):
    """Python reference for nscale_sun and nscale_sha."""
    fn = _fn_ref(kn, pai_above, dpai)
    if lot == 1:
        ns = fn / dpai
        return ns, ns
    # lot == 0
    fn_sun = (clump_fac / (kn + kb * clump_fac)
              * math.exp(-kn * pai_above) * tbi
              * (1.0 - math.exp(-(kn + kb * clump_fac) * dpai)))
    fn_sha = fn - fn_sun
    nscale_sun = fn_sun / (fracsun * dpai)
    nscale_sha = fn_sha / ((1.0 - fracsun) * dpai)
    return nscale_sun, nscale_sha


# ──────────────────────────────────────────────────────────────────────────────
# Reference value tests
# ──────────────────────────────────────────────────────────────────────────────

class TestReferenceValues:

    def test_type0_top_layer_matches_reference(self):
        """leaf_optics_type=0, canopy top (pai_above=0, tbi=1): match Python."""
        inp = dict(kn=0.3, pai_above=0.0, dpai=0.5, kb=0.5,
                   clump_fac=0.8, fracsun=0.4, tbi=1.0, leaf_optics_type_in=0)
        out = run_fortran(EXE, inp)
        exp_sun, exp_sha = _nscale_ref(
            inp['kn'], inp['pai_above'], inp['dpai'], inp['kb'],
            inp['clump_fac'], inp['fracsun'], inp['tbi'], 0)
        assert abs(out['nscale_sun'] - exp_sun) < 1e-10, \
            f"nscale_sun={out['nscale_sun']}, expected {exp_sun}"
        assert abs(out['nscale_sha'] - exp_sha) < 1e-10, \
            f"nscale_sha={out['nscale_sha']}, expected {exp_sha}"

    def test_type0_mid_canopy_matches_reference(self):
        """leaf_optics_type=0, mid-canopy (pai_above=1.0, tbi<1): match Python."""
        tbi = math.exp(-0.5 * 1.0)  # approximate transmittance at pai_above=1
        inp = dict(kn=0.3, pai_above=1.0, dpai=0.4, kb=0.5,
                   clump_fac=0.8, fracsun=0.3, tbi=tbi, leaf_optics_type_in=0)
        out = run_fortran(EXE, inp)
        exp_sun, exp_sha = _nscale_ref(
            inp['kn'], inp['pai_above'], inp['dpai'], inp['kb'],
            inp['clump_fac'], inp['fracsun'], inp['tbi'], 0)
        assert abs(out['nscale_sun'] - exp_sun) < 1e-10
        assert abs(out['nscale_sha'] - exp_sha) < 1e-10

    def test_type1_sun_equals_shade(self):
        """leaf_optics_type=1: nscale_sun == nscale_sha (thin-layer)."""
        inp = dict(kn=0.3, pai_above=0.0, dpai=0.5, kb=0.5,
                   clump_fac=0.8, fracsun=0.4, tbi=1.0, leaf_optics_type_in=1)
        out = run_fortran(EXE, inp)
        assert abs(out['nscale_sun'] - out['nscale_sha']) < 1e-15, \
            f"type=1: nscale_sun={out['nscale_sun']} ≠ nscale_sha={out['nscale_sha']}"

    def test_type1_equals_fn_over_dpai(self):
        """leaf_optics_type=1: nscale = fn/dpai."""
        kn, pai_above, dpai = 0.3, 0.5, 0.4
        inp = dict(kn=kn, pai_above=pai_above, dpai=dpai, kb=0.5,
                   clump_fac=0.8, fracsun=0.4, tbi=0.8, leaf_optics_type_in=1)
        out = run_fortran(EXE, inp)
        expected = _fn_ref(kn, pai_above, dpai) / dpai
        assert abs(out['nscale_sun'] - expected) < 1e-10, \
            f"nscale_sun={out['nscale_sun']}, expected fn/dpai={expected}"


# ──────────────────────────────────────────────────────────────────────────────
# Property tests
# ──────────────────────────────────────────────────────────────────────────────

class TestProperties:

    def _run(self, **kw):
        inp = dict(kn=0.3, pai_above=0.0, dpai=0.5, kb=0.5,
                   clump_fac=0.8, fracsun=0.4, tbi=1.0, leaf_optics_type_in=0)
        inp.update(kw)
        return run_fortran(EXE, inp)

    def test_nscale_sun_positive(self):
        """nscale_sun > 0 for typical inputs."""
        out = self._run()
        assert out['nscale_sun'] > 0.0, f"nscale_sun={out['nscale_sun']} not positive"

    def test_nscale_sha_positive(self):
        """nscale_sha > 0 for typical inputs."""
        out = self._run()
        assert out['nscale_sha'] > 0.0, f"nscale_sha={out['nscale_sha']} not positive"

    def test_weighted_mean_consistency(self):
        """Weighted mean of nscale equals fn/dpai (consistency condition).

        nscale_sun * fracsun + nscale_sha * (1-fracsun) = fn / dpai
        """
        for fracsun in [0.2, 0.4, 0.6]:
            kn, pai_above, dpai = 0.3, 0.0, 0.5
            inp = dict(kn=kn, pai_above=pai_above, dpai=dpai, kb=0.5,
                       clump_fac=0.8, fracsun=fracsun, tbi=1.0,
                       leaf_optics_type_in=0)
            out = run_fortran(EXE, inp)
            fn = _fn_ref(kn, pai_above, dpai)
            weighted = out['nscale_sun'] * fracsun + out['nscale_sha'] * (1.0 - fracsun)
            expected = fn / dpai
            assert abs(weighted - expected) < 1e-10, \
                f"fracsun={fracsun}: weighted nscale={weighted}, expected {expected}"

    def test_nscale_decreases_with_depth(self):
        """Deeper layers (larger pai_above) have smaller nscale (Beer's law)."""
        pai_aboves = [0.0, 0.5, 1.0, 2.0]
        nscales = [self._run(pai_above=p)['nscale_sun'] for p in pai_aboves]
        for i in range(len(nscales) - 1):
            assert nscales[i] > nscales[i + 1], \
                f"nscale_sun not decreasing with depth: pai_above={pai_aboves[i]}"

    def test_nscale_decreases_with_kn(self):
        """Larger kn → steeper nitrogen gradient → lower nscale at top layer."""
        # At pai_above=0: nscale = fn/dpai = (1-exp(-kn*dpai))/(kn*dpai)
        # This decreases as kn increases (for fixed dpai)
        kns = [0.1, 0.3, 0.5, 1.0]
        # Use type=1 for clean formula
        nscales = [self._run(kn=k, leaf_optics_type_in=1)['nscale_sun'] for k in kns]
        for i in range(len(nscales) - 1):
            assert nscales[i] > nscales[i + 1], \
                f"nscale not decreasing with kn: kn={kns[i]}"

    def test_type0_and_type1_agree_for_thin_layer(self):
        """For very thin dpai, type=0 and type=1 give similar nscale (weighted mean)."""
        # For very small dpai, fn ≈ exp(-kn*pai_above) * dpai
        # Both types should give approximately exp(-kn*pai_above)
        # They won't be identical, but their weighted means should both
        # equal fn/dpai
        kn, dpai = 0.3, 0.01
        fracsun = 0.4
        inp0 = dict(kn=kn, pai_above=0.0, dpai=dpai, kb=0.5,
                    clump_fac=0.8, fracsun=fracsun, tbi=1.0, leaf_optics_type_in=0)
        inp1 = dict(**inp0, leaf_optics_type_in=1)
        out0 = run_fortran(EXE, inp0)
        out1 = run_fortran(EXE, inp1)
        fn = _fn_ref(kn, 0.0, dpai)
        wm0 = out0['nscale_sun'] * fracsun + out0['nscale_sha'] * (1.0 - fracsun)
        wm1 = out1['nscale_sun'] * fracsun + out1['nscale_sha'] * (1.0 - fracsun)
        assert abs(wm0 - fn / dpai) < 1e-10
        assert abs(wm1 - fn / dpai) < 1e-10
