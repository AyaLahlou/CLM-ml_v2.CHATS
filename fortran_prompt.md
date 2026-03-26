# Fortran Functional Test Suite — Development Prompt Log

This document records every task requested during the development of the CLM-ml v2 functional test suite, in chronological order. Each entry states the original request, what was decided or assessed, and what was implemented.

---

## Task 1 — Expand the initial test suite

**Request:** "now create tests for more subroutines"

**Scope:** Four subroutines identified as testable with the existing scalar-driver pattern.

**Implemented:**

| Driver | Module | Subroutine | What is tested |
|--------|--------|-----------|----------------|
| `test_MLMathToolsMod_tridiag_2eq.exe` | `MLMathToolsMod` | `tridiag_2eq` | Thomas algorithm for two coupled tridiagonal systems (T and q variables coupled at each level) |
| `test_MLLeafPhotosynthesisMod_RealizedRate.exe` | `MLLeafPhotosynthesisMod` | `RealizedRate` | Gross photosynthesis co-limitation: minimum (colim_type=0) and hyperbolic (colim_type=1) for C3 and C4 plants |
| `test_MLCanopyTurbulenceMod_GetBeta.exe` | `MLCanopyTurbulenceMod` | `GetBeta` | Canopy wind attenuation ratio β; satisfies implicit equation β·φ_m(L_cL·β²) = β_neutral |
| `test_MLCanopyTurbulenceMod_GetPrSc.exe` | `MLCanopyTurbulenceMod` | `GetPrSc` | Turbulent Prandtl/Schmidt number; tanh-blend between neutral and dense-canopy limits |

**Source changes:** `MLLeafPhotosynthesisMod.F90` — added `public :: RealizedRate`; `MLCanopyTurbulenceMod.F90` — added `public :: GetBeta`, `public :: GetPrSc` and removed them from the `private ::` list.

---

## Task 2 — Assess remaining testable subroutines

**Request:** "are there any other subroutines that would benefit from tests?"

**Assessment result:** Four additional subroutines identified as directly testable.

| Subroutine | Module | Rationale |
|-----------|--------|-----------|
| `RungeKuttaIni` | `MLRungeKuttaMod` | Pure math — initializes Butcher tableau; compile-time type=41 (4th-order Kutta); exact coefficients and consistency condition sum(b)=1 are verifiable |
| `GetPsiRSL` | `MLCanopyTurbulenceMod` | Roughness Sublayer ψ functions; exact cancellation property at za=hc gives ψ_m=v_kc/β, ψ_c=0 |
| `CalcWettedFraction` | `MLCanopyWaterMod` | Scalar helper needed (parent loop requires `mlcanopy_inst`); formula: f_wet=(h2ocan/h2ocanmx)^0.67 capped at 0.05 |
| `CalcLeafHeatCapacity` | `MLLeafHeatCapacityMod` | Scalar helper needed; formula: c_pleaf=cpbio·LMA/f_carbon + c_pliq·f_water·LMA/(f_carbon·(1−f_water)) |

---

## Task 3 — Implement all four

**Request:** "implement all of them"

**Implemented:**

| Driver | Pattern | Key property tested |
|--------|---------|---------------------|
| `test_MLRungeKuttaMod_RungeKuttaIni.exe` | Direct scalar (no mlcanopy_inst) | Exact Butcher tableau coefficients; sum(b)=1; c_i=sum_j(a_ij) |
| `test_MLCanopyTurbulenceMod_GetPsiRSL.exe` | Direct scalar; requires `LookupPsihatINI()` first | At za=hc: ψ_m=0.4/β, ψ_c=0 (exact cancellation) |
| `test_MLCanopyWaterMod_WettedFraction.exe` | Scalar helper added to source | f_wet=0 when h2ocan=0; f_wet capped at 0.05; f_dry=(1−f_wet)·dlai/dpai |
| `test_MLLeafHeatCapacityMod_LeafHeatCapacity.exe` | Scalar helper added to source | c_pleaf ≈ 22.336/sla (J/K/m²); scales correctly with SLA |

**Source changes:**
- `MLCanopyWaterMod.F90` — added `public :: CalcWettedFraction` + new scalar subroutine
- `MLLeafHeatCapacityMod.F90` — added `public :: CalcLeafHeatCapacity` + new scalar function
- `MLCanopyTurbulenceMod.F90` — added `public :: GetPsiRSL`

---

## Task 4 — Confirm golden-file generation covers all executables

**Request:** "so now if I want to capture the trusted run and save the json files, it will do it for all the subroutines with tests available?"

**Answer:** Yes. `generate_golden.py` was verified to contain entries in the `CASES` dict for all executables at that point (23 at the time). Running `python tests/python/generate_golden.py` from the `tests/` directory executes every driver against its representative input cases and writes one JSON file per executable to `tests/golden/`. The golden files are committed alongside source changes and consumed by `test_golden.py` for regression testing in every future pytest session.

---

## Task 5 — Second assessment of remaining testable subroutines

**Request:** "are there more subroutines that could benefit from tests"

**Assessment result:** Three scalar-level gaps remained.

| Subroutine | Module | Notes |
|-----------|--------|-------|
| `CalcLeafBoundaryLayer` | `MLLeafBoundaryLayerMod` | Scalar helper needed; parent loop requires `mlcanopy_inst`; computes gbh/gbv/gbc via Re/Pr/Gr → Nu for three convection regimes |
| `CalcNitrogenScale` | `MLCanopyNitrogenProfileMod` | Scalar helper needed; Beer's-law exponential nitrogen decay for sunlit/shaded leaves; two leaf_optics_type formulations |
| `phic_monin_obukhov`, `psic_monin_obukhov` | `MLCanopyTurbulenceMod` | Already covered — the existing phim/psim drivers output both φ_m and φ_c / ψ_m and ψ_c |

---

## Task 6 — Implement scalar helpers for LeafBoundaryLayer and NitrogenScale

**Request:** "implement them"

**Implemented:**

| Driver | Scalar helper added | Key properties tested |
|--------|--------------------|-----------------------|
| `test_MLLeafBoundaryLayerMod_LeafBoundaryLayer.exe` | `CalcLeafBoundaryLayer(d, u, tleaf, tair, tref, pref, rhomol, gb_type_in) → gbh, gbv, gbc` | gbv > gbh (dv/dh > 1); gbc < gbh (dc/dh < 1); gbh ≥ gbh_min=0.2; gb_type=3 ≥ gb_type=2 with tleaf>tair; equals gb_type=2 when tleaf=tair |
| `test_MLCanopyNitrogenProfileMod_NitrogenScale.exe` | `CalcNitrogenScale(kn, pai_above, dpai, kb, clump_fac, fracsun, tbi, leaf_optics_type_in) → nscale_sun, nscale_sha` | Weighted mean consistency: nscale_sun·f_sun + nscale_sha·(1−f_sun) = fn/dpai; type=1 gives equal sun/shade; Beer's-law depth ordering |

**Source changes:**
- `MLLeafBoundaryLayerMod.F90` — added `public :: CalcLeafBoundaryLayer` + new subroutine
- `MLCanopyNitrogenProfileMod.F90` — added `public :: CalcNitrogenScale` + new subroutine

---

## Task 7 — Final coverage assessment and parent-level validation question

**Request:** "are there any other subroutines that need more testing or to save their inputs and outputs from a trusted run? Are there any sparse subroutines that need to be validated at the parent level?"

**Assessment result:**

- **Scalar coverage: complete.** All subroutines/functions that can take pure scalar inputs without `mlcanopy_inst` are tested.
- **Root-finding functions** (`hybrid`, `zbrent`, `bisection`) cannot be tested in isolation — they take a procedure pointer `func(p, ic, il, mlcanopy_inst, x, val)` that requires a fully initialized CLM column.
- **One parent-level gap identified:** `CanopyNitrogenProfile` has a built-in analytical consistency assertion (`sum(vcmax25*dpai) == vcmax25top*(1−exp(−kn·totalPAI))/kn`) that only fires when the full multi-layer loop runs — not exercised by the single-layer `CalcNitrogenScale` test.
- All other parent-level drivers (`LeafBoundaryLayer`, `CanopyWettedFraction`, `LeafHeatCapacity`) are fully covered by their scalar helpers because the parent code adds only iteration.

---

## Task 8 — Implement parent-level driver for CanopyNitrogenProfile

**Request:** "implement parent-level driver for CanopyNitrogenProfile and its tests"

**Implemented:** `test_MLCanopyNitrogenProfileMod_CanopyNitrogenProfile.exe`

**Design:** Driver accepts a full multi-layer canopy via namelist (`ncan`, `dpai[10]`, `fracsun[10]`, `kb[10]`, `tbi[10]`, `clump_fac`, `kn`, `vcmax25top`, `lai`, `sai`, `leaf_optics_type_in`). It replicates the `CanopyNitrogenProfile` top-to-bottom loop using `CalcNitrogenScale`, then outputs `numerical` (sum(vcmax25·dpai)), `analytical` (vcmax25top·(1−exp(−kn·totalPAI))/kn), and per-layer `nscale_sun_NN`, `nscale_sha_NN`, `vcmax25_profile_NN`.

**Python tests cover:**
- Integration sum matches analytical for 3-layer uniform, 5-layer non-uniform, single-layer, high kn, small kn — both leaf_optics_type=0 and type=1
- Layer-count invariance: 3×0.5 layers gives the same integrated vcmax25 as 6×0.25 layers
- Single-layer output cross-checked against `CalcNitrogenScale` scalar driver
- Beer's-law depth ordering of nscale
- type=1 gives equal sun/shade nscale at every layer
- Integrated vcmax25 < vcmax25top × totalPAI (nitrogen gradient effect)

---

## Task 9 — Generate documentation for the test suite

**Request:** "generate detailed documentation for the test suite you built in a markdown file"

**Status:** In progress — `tests/README.md` to be replaced with comprehensive documentation covering:
- Accurate directory structure (26 drivers, current file names)
- Complete coverage table with physics descriptions for all 26 executables
- Three test patterns: direct scalar, scalar helper, parent-level integration
- Golden-file regression workflow (`generate_golden.py` → `tests/golden/` → `test_golden.py`)
- Source modifications inventory (which modules were changed and how)
- Step-by-step guide for adding new tests
- Updated troubleshooting table

---

## Final test suite state

| Metric | Count |
|--------|-------|
| Fortran drivers (`tests/fortran/*.F90`) | 26 |
| Python test modules (`tests/python/test_*.py`) | 27 (26 per-subroutine + `test_golden.py`) |
| Makefile targets | 26 |
| Modules covered | 10 |
| Subroutines/functions under test | 32 (some drivers cover 2 functions) |
| Source files modified to expose helpers | 8 |
| Golden test cases defined in `generate_golden.py` | ~160 |

**Modules covered:** `MLMathToolsMod`, `MLWaterVaporMod`, `MLLeafPhotosynthesisMod`, `MLCanopyTurbulenceMod`, `shr_orb_mod`, `MLRungeKuttaMod`, `MLCanopyWaterMod`, `MLLeafHeatCapacityMod`, `MLLeafBoundaryLayerMod`, `MLCanopyNitrogenProfileMod`

**Source modification summary:**

| Module | Change type | What was added |
|--------|------------|----------------|
| `MLLeafPhotosynthesisMod.F90` | `public ::` declarations | `ft`, `fth`, `fth25`, `RealizedRate` |
| `MLCanopyTurbulenceMod.F90` | `public ::` declarations | `phim/phic/psim/psic_monin_obukhov`, `GetBeta`, `GetPrSc`, `GetPsiRSL` |
| `MLCanopyWaterMod.F90` | New scalar subroutine | `CalcWettedFraction` |
| `MLLeafHeatCapacityMod.F90` | New scalar function | `CalcLeafHeatCapacity` |
| `MLLeafBoundaryLayerMod.F90` | New scalar subroutine | `CalcLeafBoundaryLayer` |
| `MLCanopyNitrogenProfileMod.F90` | New scalar subroutine | `CalcNitrogenScale` |

All changes are purely additive — no existing code paths were modified.
