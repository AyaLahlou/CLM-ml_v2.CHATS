"""
Parent-level integration tests for MLCanopyNitrogenProfileMod :: CanopyNitrogenProfile.

This driver replicates the full multi-layer nitrogen-profile loop using the
public CalcNitrogenScale helper and validates the central consistency condition:

    sum_ic( vcmax25_profile(ic) * dpai(ic) ) = vcmax25top * (1 - exp(-kn*totalPAI)) / kn

This is exactly the assertion that CanopyNitrogenProfile makes internally (the
call to endrun when the sum disagrees).  It is a mathematical property of the
Beer's-law exponential nitrogen decay: the discrete sum over arbitrary layers
always reproduces the closed-form integral, regardless of the number of layers
or the distribution of dpai / fracsun / kb / tbi.

Single-layer tests additionally cross-check against CalcNitrogenScale directly.

Conventions (matching CanopyNitrogenProfile and the test driver):
  ic = 1       →  bottom layer (largest pai_above)
  ic = ncan    →  top layer    (pai_above = 0)
  totalPAI     = lai + sai
"""

import math
import pytest
from utils import run_fortran

EXE = 'test_MLCanopyNitrogenProfileMod_CanopyNitrogenProfile.exe'


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _fn(kn, pai_above, dpai):
    """Nitrogen integral for one layer."""
    return math.exp(-kn * pai_above) * (1.0 - math.exp(-kn * dpai)) / kn


def _analytical(vcmax25top, kn, total_pai):
    """Analytical canopy integral of vcmax25."""
    return vcmax25top * (1.0 - math.exp(-kn * total_pai)) / kn


def _run(ncan, dpai_list, fracsun_list, kb_list, tbi_list,
         clump_fac=1.0, kn=0.3, vcmax25top=60.0,
         lai=None, sai=0.0, lot=0):
    """
    Run the Fortran driver and return (out_dict, nscale_sun, nscale_sha, vcmax25).
    Arrays are 1-indexed (ic=1 bottom, ic=ncan top), matching Fortran convention.
    """
    assert len(dpai_list)    == ncan
    assert len(fracsun_list) == ncan
    assert len(kb_list)      == ncan
    assert len(tbi_list)     == ncan

    if lai is None:
        lai = sum(dpai_list)

    # Pad to ncan_max=10
    pad = [0.0] * (10 - ncan)

    inp = dict(
        ncan=ncan,
        dpai    = dpai_list    + pad,
        fracsun = fracsun_list + pad,
        kb      = kb_list      + pad,
        tbi     = tbi_list     + pad,
        clump_fac=clump_fac,
        kn=kn, vcmax25top=vcmax25top,
        lai=lai, sai=sai,
        leaf_optics_type_in=lot,
    )
    out = run_fortran(EXE, inp)

    nscale_sun = [out[f'nscale_sun_{ic:02d}']      for ic in range(1, ncan + 1)]
    nscale_sha = [out[f'nscale_sha_{ic:02d}']      for ic in range(1, ncan + 1)]
    vcmax25     = [out[f'vcmax25_profile_{ic:02d}'] for ic in range(1, ncan + 1)]

    return out, nscale_sun, nscale_sha, vcmax25


def _tbi_uniform(ncan, dpai, kb, clump_fac=1.0):
    """Compute tbi array for a uniform-dpai canopy (ic=ncan is top)."""
    # pai_above at each layer (accumulated from top)
    tbi  = [0.0] * ncan
    pai  = 0.0
    for ic in range(ncan, 0, -1):         # top → bottom
        tbi[ic - 1] = math.exp(-kb * clump_fac * pai)
        pai += dpai
    return tbi                            # index 0 = bottom, ncan-1 = top


# ─────────────────────────────────────────────────────────────────────────────
# Integration sum tests (the core parent-level property)
# ─────────────────────────────────────────────────────────────────────────────

class TestIntegrationSum:
    """Verify sum(vcmax25_profile * dpai) == analytical for all canopy configs."""

    ATOL = 1e-8   # generous tolerance: Fortran float64 rounding

    def _check_sum(self, out):
        assert abs(out['numerical'] - out['analytical']) < self.ATOL, (
            f"Integration error: numerical={out['numerical']:.6e}, "
            f"analytical={out['analytical']:.6e}, "
            f"diff={abs(out['numerical']-out['analytical']):.2e}"
        )

    def test_3layer_uniform_type0(self):
        """3-layer uniform canopy (leaf_optics_type=0) satisfies integration sum."""
        n = 3; dpai = 0.5; kb = 0.5; cf = 1.0
        tbi = _tbi_uniform(n, dpai, kb, cf)
        out, *_ = _run(n,
                       dpai_list   =[dpai]*n,
                       fracsun_list=[0.20, 0.35, 0.50],
                       kb_list     =[kb]*n,
                       tbi_list    =tbi,
                       clump_fac=cf, kn=0.3, vcmax25top=60.0)
        self._check_sum(out)

    def test_3layer_uniform_type1(self):
        """3-layer uniform canopy (leaf_optics_type=1) satisfies integration sum."""
        n = 3; dpai = 0.5
        tbi = [0.6065, 0.7788, 1.0]  # not actually used for type=1
        out, *_ = _run(n,
                       dpai_list   =[dpai]*n,
                       fracsun_list=[0.20, 0.35, 0.50],
                       kb_list     =[0.5]*n,
                       tbi_list    =tbi,
                       lot=1)
        self._check_sum(out)

    def test_5layer_nonuniform_dpai(self):
        """5-layer non-uniform dpai canopy satisfies integration sum."""
        dpai_list   = [0.2, 0.3, 0.4, 0.4, 0.2]  # bottom to top
        fracsun_list= [0.10, 0.20, 0.30, 0.40, 0.55]
        kb_list     = [0.45, 0.50, 0.50, 0.52, 0.48]
        cf = 0.85
        n = len(dpai_list)
        # Compute tbi (ic=5 is top)
        tbi = [0.0] * n
        pai = 0.0
        for ic in range(n, 0, -1):
            tbi[ic-1] = math.exp(-kb_list[ic-1] * cf * pai)
            pai += dpai_list[ic-1]
        out, *_ = _run(n, dpai_list, fracsun_list, kb_list, tbi,
                       clump_fac=cf, kn=0.25, vcmax25top=80.0,
                       lai=sum(dpai_list))
        self._check_sum(out)

    def test_1layer_type0(self):
        """Single-layer canopy satisfies integration sum."""
        out, *_ = _run(1,
                       dpai_list   =[1.5],
                       fracsun_list=[0.4],
                       kb_list     =[0.5],
                       tbi_list    =[1.0],
                       kn=0.3, vcmax25top=60.0, lai=1.5)
        self._check_sum(out)

    def test_large_kn(self):
        """High kn (steep gradient) still satisfies integration sum."""
        n = 5; dpai = 0.4
        tbi = _tbi_uniform(n, dpai, 0.5, 1.0)
        out, *_ = _run(n,
                       dpai_list   =[dpai]*n,
                       fracsun_list=[0.15, 0.25, 0.35, 0.45, 0.55],
                       kb_list     =[0.5]*n,
                       tbi_list    =tbi,
                       kn=1.5, vcmax25top=50.0)
        self._check_sum(out)

    def test_small_kn(self):
        """Small kn (shallow gradient) still satisfies integration sum."""
        n = 4; dpai = 0.5
        tbi = _tbi_uniform(n, dpai, 0.5, 1.0)
        out, *_ = _run(n,
                       dpai_list   =[dpai]*n,
                       fracsun_list=[0.2, 0.3, 0.4, 0.5],
                       kb_list     =[0.5]*n,
                       tbi_list    =tbi,
                       kn=0.05, vcmax25top=70.0)
        self._check_sum(out)

    def test_analytical_value_matches_python(self):
        """numerical and analytical both match Python formula."""
        kn = 0.3; vcmax25top = 60.0; total_pai = 1.5
        n = 3
        tbi = _tbi_uniform(n, 0.5, 0.5, 1.0)
        out, *_ = _run(n,
                       dpai_list   =[0.5]*n,
                       fracsun_list=[0.20, 0.35, 0.50],
                       kb_list     =[0.5]*n,
                       tbi_list    =tbi,
                       kn=kn, vcmax25top=vcmax25top, lai=total_pai)
        expected = _analytical(vcmax25top, kn, total_pai)
        assert abs(out['analytical'] - expected) < 1e-12, \
            f"analytical={out['analytical']}, expected Python={expected}"
        assert abs(out['numerical'] - expected) < self.ATOL, \
            f"numerical={out['numerical']}, expected Python={expected}"


# ─────────────────────────────────────────────────────────────────────────────
# Property tests
# ─────────────────────────────────────────────────────────────────────────────

class TestProperties:

    def test_type1_sun_equals_shade_all_layers(self):
        """leaf_optics_type=1: nscale_sun == nscale_sha for every layer."""
        n = 4; dpai = 0.4
        tbi = [0.4724, 0.5696, 0.6873, 1.0]  # approximate Beer's law values
        _, ns, nsh, _ = _run(n,
                             dpai_list   =[dpai]*n,
                             fracsun_list=[0.15, 0.25, 0.35, 0.5],
                             kb_list     =[0.5]*n,
                             tbi_list    =tbi,
                             lot=1)
        for ic, (s, sh) in enumerate(zip(ns, nsh), 1):
            assert abs(s - sh) < 1e-15, \
                f"type=1: layer {ic}: nscale_sun={s} ≠ nscale_sha={sh}"

    def test_nscale_decreases_bottom_to_top(self):
        """Canopy top (largest index) has highest nscale; bottom has lowest.

        Beer's law: exp(-kn*pai_above) is largest at the top (pai_above=0).
        """
        n = 4; dpai = 0.5; kn = 0.4
        tbi = _tbi_uniform(n, dpai, 0.5, 1.0)
        _, ns, _, _ = _run(n,
                           dpai_list   =[dpai]*n,
                           fracsun_list=[0.2, 0.3, 0.4, 0.5],
                           kb_list     =[0.5]*n,
                           tbi_list    =tbi,
                           kn=kn, lot=1)   # type=1 gives clean exponential
        # ic=4 is top (nscale highest), ic=1 is bottom (nscale lowest)
        for ic in range(n - 1):
            assert ns[ic] < ns[ic + 1], \
                f"nscale_sun[{ic+1}]={ns[ic]:.4f} should be < nscale_sun[{ic+2}]={ns[ic+1]:.4f}"

    def test_single_layer_matches_calc_nitrogen_scale(self):
        """1-layer driver output matches the single-layer CalcNitrogenScale driver."""
        from utils import run_fortran as rf
        kn, dpai, kn_pai, fracsun, cf, kb = 0.3, 1.0, 0.0, 0.4, 0.8, 0.5
        tbi = 1.0

        # Parent-level driver (1-layer)
        out_parent, ns_p, nsh_p, _ = _run(
            1,
            dpai_list   =[dpai],
            fracsun_list=[fracsun],
            kb_list     =[kb],
            tbi_list    =[tbi],
            clump_fac=cf, kn=kn, lai=dpai, lot=0,
        )

        # Single-layer CalcNitrogenScale driver
        out_scalar = rf('test_MLCanopyNitrogenProfileMod_NitrogenScale.exe', dict(
            kn=kn, pai_above=0.0, dpai=dpai,
            kb=kb, clump_fac=cf, fracsun=fracsun, tbi=tbi,
            leaf_optics_type_in=0,
        ))

        assert abs(ns_p[0]  - out_scalar['nscale_sun']) < 1e-13, \
            f"nscale_sun mismatch: parent={ns_p[0]}, scalar={out_scalar['nscale_sun']}"
        assert abs(nsh_p[0] - out_scalar['nscale_sha']) < 1e-13, \
            f"nscale_sha mismatch: parent={nsh_p[0]}, scalar={out_scalar['nscale_sha']}"

    def test_vcmax25_profile_positive(self):
        """All vcmax25_profile values are positive for non-zero dpai."""
        n = 3
        tbi = _tbi_uniform(n, 0.5, 0.5, 1.0)
        _, _, _, vcmax25 = _run(n,
                                dpai_list   =[0.5]*n,
                                fracsun_list=[0.2, 0.35, 0.5],
                                kb_list     =[0.5]*n,
                                tbi_list    =tbi)
        for ic, v in enumerate(vcmax25, 1):
            assert v > 0.0, f"vcmax25_profile[{ic}]={v} not positive"

    def test_integration_sum_is_less_than_vcmax25top_times_total_pai(self):
        """The integrated canopy vcmax25 is less than vcmax25top * totalPAI.

        The exponential profile concentrates nitrogen at the top, so the
        canopy average is less than the top-of-canopy value.
        """
        kn = 0.3; vcmax25top = 60.0; total_pai = 2.0
        n = 4; dpai = 0.5
        tbi = _tbi_uniform(n, dpai, 0.5, 1.0)
        out, *_ = _run(n,
                       dpai_list   =[dpai]*n,
                       fracsun_list=[0.15, 0.25, 0.35, 0.50],
                       kb_list     =[0.5]*n,
                       tbi_list    =tbi,
                       kn=kn, vcmax25top=vcmax25top, lai=total_pai)
        assert out['numerical'] < vcmax25top * total_pai, (
            f"Integrated vcmax25={out['numerical']:.4f} should be less than "
            f"vcmax25top*totalPAI={vcmax25top*total_pai:.4f}"
        )

    def test_more_layers_same_result(self):
        """Refining 3 layers into 6 equal sub-layers gives the same integrated vcmax25."""
        kn = 0.3; vcmax25top = 60.0; fracsun_val = 0.4; kb = 0.5; cf = 1.0

        # 3 layers, each dpai=0.5
        n3 = 3; dpai3 = 0.5
        tbi3 = _tbi_uniform(n3, dpai3, kb, cf)
        out3, *_ = _run(n3,
                        dpai_list   =[dpai3]*n3,
                        fracsun_list=[fracsun_val]*n3,
                        kb_list     =[kb]*n3,
                        tbi_list    =tbi3,
                        clump_fac=cf, kn=kn, vcmax25top=vcmax25top, lot=1)

        # 6 layers, each dpai=0.25
        n6 = 6; dpai6 = 0.25
        tbi6 = _tbi_uniform(n6, dpai6, kb, cf)
        out6, *_ = _run(n6,
                        dpai_list   =[dpai6]*n6,
                        fracsun_list=[fracsun_val]*n6,
                        kb_list     =[kb]*n6,
                        tbi_list    =tbi6,
                        clump_fac=cf, kn=kn, vcmax25top=vcmax25top,
                        lai=n6*dpai6, lot=1)

        assert abs(out3['numerical'] - out6['numerical']) < 1e-10, (
            f"3-layer sum={out3['numerical']:.8f}, "
            f"6-layer sum={out6['numerical']:.8f}, "
            f"diff={abs(out3['numerical']-out6['numerical']):.2e}"
        )
