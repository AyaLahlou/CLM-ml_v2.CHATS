#!/usr/bin/env python3
"""
Generate golden (trusted-run) data for all Fortran test executables.

Usage (from the tests/ directory):
    python python/generate_golden.py

Or from anywhere:
    python tests/python/generate_golden.py

For each executable a JSON file is written to tests/golden/<exe_name>.json.
Each JSON file records the inputs and outputs of a representative set of
test cases captured from a verified, trusted build of the code.

Re-run this script whenever the Fortran code is intentionally changed and
the new outputs should become the new trusted reference.  Commit the updated
JSON files alongside the code change so that future test_golden.py runs
compare against the correct baseline.
"""

import json
import os
import sys
from datetime import datetime, timezone

# Allow 'from utils import ...' regardless of where the script is invoked from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import run_fortran  # noqa: E402

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'golden')

# ---------------------------------------------------------------------------
# Representative input cases for each executable
# Each entry in the list is a plain dict of namelist inputs.
# ---------------------------------------------------------------------------

CASES = {
    # ── MLMathToolsMod ──────────────────────────────────────────────────────
    'test_MLMathToolsMod_quadratic.exe': [
        {'a':  1.0, 'b':  -3.0, 'c':  2.0},   # roots 1, 2
        {'a':  1.0, 'b':  -2.0, 'c':  1.0},   # double root 1
        {'a':  1.0, 'b':   5.0, 'c':  6.0},   # roots -3, -2
        {'a':  2.0, 'b':  -7.0, 'c':  3.0},   # Vieta product check
        {'a':  3.0, 'b': -11.0, 'c':  6.0},   # Vieta sum check
    ],
    'test_MLMathToolsMod_tridiag.exe': [
        # 3×3 symmetric, solution [1,1,1]
        {'n': 3,
         'a': [0.0,  -1.0, -1.0],
         'b': [2.0,   2.0,  2.0],
         'c': [-1.0, -1.0,  0.0],
         'r': [1.0,   0.0,  1.0]},
        # 4×4 diagonal, solution [3,3,3,3]
        {'n': 4,
         'a': [0.0,  0.0,  0.0,  0.0],
         'b': [2.0,  3.0,  4.0,  5.0],
         'c': [0.0,  0.0,  0.0,  0.0],
         'r': [6.0,  9.0, 12.0, 15.0]},
        # 3×3 identity, solution = rhs
        {'n': 3,
         'a': [0.0, 0.0, 0.0],
         'b': [1.0, 1.0, 1.0],
         'c': [0.0, 0.0, 0.0],
         'r': [7.0, 3.0, 5.0]},
    ],
    'test_MLMathToolsMod_log_gamma_function.exe': [
        {'x': 1.0},
        {'x': 2.0},
        {'x': 3.0},
        {'x': 0.5},
        {'x': 0.3},
        {'x': 5.0},
        {'x': 10.0},
    ],
    'test_MLMathToolsMod_beta_function.exe': [
        {'a': 1.0, 'b': 1.0},
        {'a': 0.5, 'b': 0.5},
        {'a': 2.0, 'b': 3.0},
        {'a': 2.0, 'b': 5.0},
        {'a': 5.0, 'b': 2.0},
        {'a': 3.7, 'b': 2.1},
    ],
    'test_MLMathToolsMod_beta_distribution_pdf.exe': [
        {'a': 1.0, 'b': 1.0, 'x': 0.1},
        {'a': 1.0, 'b': 1.0, 'x': 0.5},
        {'a': 1.0, 'b': 1.0, 'x': 0.9},
        {'a': 2.0, 'b': 3.0, 'x': 0.3},
        {'a': 3.0, 'b': 2.0, 'x': 0.7},
        {'a': 2.0, 'b': 5.0, 'x': 0.4},
        {'a': 2.0, 'b': 2.0, 'x': 0.5},
    ],
    'test_MLMathToolsMod_beta_distribution_cdf.exe': [
        {'a': 2.0, 'b': 3.0, 'x': 0.0},
        {'a': 2.0, 'b': 3.0, 'x': 1.0},
        {'a': 3.0, 'b': 3.0, 'x': 0.5},
        {'a': 2.0, 'b': 2.0, 'x': 0.1},
        {'a': 2.0, 'b': 2.0, 'x': 0.3},
        {'a': 2.0, 'b': 2.0, 'x': 0.5},
        {'a': 2.0, 'b': 2.0, 'x': 0.7},
        {'a': 2.0, 'b': 2.0, 'x': 0.9},
        {'a': 2.0, 'b': 2.0, 'x': 0.25},
    ],
    # ── MLWaterVaporMod ─────────────────────────────────────────────────────
    'test_MLWaterVaporMod_SatVap.exe': [
        {'t': 273.15},   # 0°C, freezing point (water branch)
        {'t': 298.15},   # 25°C
        {'t': 373.15},   # 100°C (upper clamp limit)
        {'t': 263.15},   # -10°C (ice branch)
        {'t': 198.15},   # -75°C (lower clamp limit)
        {'t': 150.0},    # below lower clamp – should equal 198.15 K result
        {'t': 400.0},    # above upper clamp – should equal 373.15 K result
        {'t': 280.0},
        {'t': 310.0},
    ],
    'test_MLWaterVaporMod_LatVap.exe': [
        {'t': 300.0},    # above freezing (evaporation)
        {'t': 260.0},    # below freezing (sublimation)
        {'t': 273.16},   # just above freezing
        {'t': 273.14},   # just below freezing
        {'t': 280.0},
        {'t': 320.0},
        {'t': 250.0},
        {'t': 200.0},
    ],
    # ── MLLeafPhotosynthesisMod ─────────────────────────────────────────────
    'test_MLLeafPhotosynthesisMod_ft.exe': [
        {'tl': 298.15, 'ha':  65330.0},   # 25°C, vcmax ha → ft = 1
        {'tl': 308.15, 'ha':  65330.0},   # 35°C, vcmax
        {'tl': 308.15, 'ha':  43540.0},   # 35°C, jmax
        {'tl': 308.15, 'ha':  79430.0},   # 35°C, kc
        {'tl': 285.0,  'ha':  65330.0},   # below 25°C
        {'tl': 310.0,  'ha':  65330.0},
        {'tl': 280.0,  'ha':  65330.0},
        {'tl': 290.0,  'ha':  65330.0},
        {'tl': 318.0,  'ha':  65330.0},
        {'tl': 300.0,  'ha':  0.0},       # ha=0 → ft = 1
    ],
    'test_MLLeafPhotosynthesisMod_fth25.exe': [
        {'hd': 150000.0, 'se': 490.0},
        {'hd': 200000.0, 'se': 490.0},
        {'hd': 150000.0, 'se': 550.0},
        {'hd': 140000.0, 'se': 490.0},
        {'hd': 160000.0, 'se': 490.0},
        {'hd': 150000.0, 'se': 450.0},
        {'hd': 150000.0, 'se': 530.0},
    ],
    'test_MLLeafPhotosynthesisMod_fth.exe': [
        {'tl': 298.15, 'hd': 150000.0, 'se': 490.0, 'c': 1.0},
        {'tl': 308.15, 'hd': 150000.0, 'se': 490.0, 'c': 1.0},
        {'tl': 313.15, 'hd': 150000.0, 'se': 490.0, 'c': 1.0},
        {'tl': 323.15, 'hd': 150000.0, 'se': 490.0, 'c': 1.0},
        {'tl': 310.0,  'hd': 150000.0, 'se': 490.0, 'c': 1.0},
        {'tl': 310.0,  'hd': 150000.0, 'se': 490.0, 'c': 2.0},
    ],
    # ── MLCanopyTurbulenceMod ───────────────────────────────────────────────
    'test_MLCanopyTurbulenceMod_phim_monin_obukhov.exe': [
        {'zeta':  0.0},
        {'zeta':  0.5},
        {'zeta':  1.0},
        {'zeta':  2.0},
        {'zeta': -0.5},
        {'zeta': -1.0},
        {'zeta': -2.0},
        {'zeta':  0.2},
    ],
    'test_MLCanopyTurbulenceMod_psim_monin_obukhov.exe': [
        {'zeta':  0.0},
        {'zeta':  0.5},
        {'zeta':  1.0},
        {'zeta':  2.0},
        {'zeta': -0.1},
        {'zeta': -0.5},
        {'zeta': -1.0},
        {'zeta': -2.0},
        {'zeta':  0.3},
    ],
    # ── shr_orb_mod ─────────────────────────────────────────────────────────
    'test_shr_orb_mod_shr_orb_params.exe': [
        {'iyear_AD': 2000},
        {'iyear_AD': 2007},
        {'iyear_AD': 1950},
        {'iyear_AD': 2050},
    ],
    'test_shr_orb_mod_shr_orb_decl.exe': [
        {'iyear_AD': 2007, 'calday':   1.0},   # Jan 1
        {'iyear_AD': 2007, 'calday':  80.0},   # vernal equinox
        {'iyear_AD': 2007, 'calday': 172.0},   # summer solstice
        {'iyear_AD': 2007, 'calday': 266.0},   # autumnal equinox
        {'iyear_AD': 2007, 'calday': 355.0},   # winter solstice
        {'iyear_AD': 2000, 'calday':  80.0},
    ],
    # ── MLMathToolsMod :: tridiag_2eq ───────────────────────────────────────
    'test_MLMathToolsMod_tridiag_2eq.exe': [
        # n=1, decoupled: T=d1/b11, q=d2/b22
        {'n': 1, 'b11': [2.0], 'b22': [3.0], 'd1': [4.0], 'd2': [6.0]},
        # n=1, coupled 2x2: T=1.6, q=1.8
        {'n': 1,
         'b11': [2.0], 'b12': [1.0], 'b21': [1.0], 'b22': [3.0],
         'd1': [5.0], 'd2': [7.0]},
        # n=2, decoupled tridiag: T=[1,1], q=[1,1]
        {'n': 2,
         'a1': [0.0, -1.0], 'b11': [2.0, 2.0], 'b12': [0.0, 0.0],
         'c1': [-1.0, 0.0], 'd1': [1.0, 1.0],
         'a2': [0.0, -1.0], 'b21': [0.0, 0.0], 'b22': [3.0, 3.0],
         'c2': [-1.0, 0.0], 'd2': [2.0, 2.0]},
        # n=2, coupled diagonal: T=[1.6,1.6], q=[1.8,1.8]
        {'n': 2,
         'a1': [0.0, 0.0], 'b11': [2.0, 2.0], 'b12': [1.0, 1.0],
         'c1': [0.0, 0.0], 'd1': [5.0, 5.0],
         'a2': [0.0, 0.0], 'b21': [1.0, 1.0], 'b22': [3.0, 3.0],
         'c2': [0.0, 0.0], 'd2': [7.0, 7.0]},
    ],
    # ── MLLeafPhotosynthesisMod :: RealizedRate ─────────────────────────────
    'test_MLLeafPhotosynthesisMod_RealizedRate.exe': [
        # colim_type=0, C3: minimum rate
        {'c3psn': 1.0, 'ac':  5.0, 'aj': 15.0, 'ap': 20.0, 'colim_type_in': 0},
        {'c3psn': 1.0, 'ac': 15.0, 'aj':  8.0, 'ap': 20.0, 'colim_type_in': 0},
        {'c3psn': 1.0, 'ac': 10.0, 'aj': 10.0, 'ap': 20.0, 'colim_type_in': 0},
        # colim_type=0, C4: minimum rate including ap
        {'c3psn': 0.0, 'ac': 10.0, 'aj': 12.0, 'ap':  4.0, 'colim_type_in': 0},
        # colim_type=1, C3: co-limited
        {'c3psn': 1.0, 'ac': 10.0, 'aj': 10.0, 'ap': 20.0, 'colim_type_in': 1},
        {'c3psn': 1.0, 'ac': 10.0, 'aj':  5.0, 'ap': 20.0, 'colim_type_in': 1},
        {'c3psn': 1.0, 'ac':  1.0, 'aj': 20.0, 'ap': 30.0, 'colim_type_in': 1},
        # colim_type=1, C4
        {'c3psn': 0.0, 'ac': 10.0, 'aj':  8.0, 'ap':  6.0, 'colim_type_in': 1},
    ],
    # ── MLCanopyTurbulenceMod :: GetBeta ────────────────────────────────────
    'test_MLCanopyTurbulenceMod_GetBeta.exe': [
        {'beta_neutral': 0.3,  'LcL':  0.0},   # neutral
        {'beta_neutral': 0.3,  'LcL': -1.0},   # unstable
        {'beta_neutral': 0.3,  'LcL':  1.0},   # stable
        {'beta_neutral': 0.3,  'LcL': -0.5},
        {'beta_neutral': 0.3,  'LcL':  0.5},
        {'beta_neutral': 0.25, 'LcL':  0.0},
        {'beta_neutral': 0.35, 'LcL': -2.0},
        {'beta_neutral': 0.3,  'LcL':  2.0},
    ],
    # ── MLCanopyTurbulenceMod :: GetPrSc ────────────────────────────────────
    'test_MLCanopyTurbulenceMod_GetPrSc.exe': [
        {'beta_neutral': 0.3,  'LcL':  0.0, 'sparse_canopy_type_in': 0},
        {'beta_neutral': 0.3,  'LcL':  1.0, 'sparse_canopy_type_in': 0},
        {'beta_neutral': 0.3,  'LcL': -1.0, 'sparse_canopy_type_in': 0},
        {'beta_neutral': 0.3,  'LcL':  0.5, 'sparse_canopy_type_in': 0},
        {'beta_neutral': 0.3,  'LcL': -0.5, 'sparse_canopy_type_in': 0},
        {'beta_neutral': 0.3,  'LcL':  0.0, 'sparse_canopy_type_in': 1},
        {'beta_neutral': 0.35, 'LcL':  0.5, 'sparse_canopy_type_in': 1},
    ],
    # ── MLRungeKuttaMod :: RungeKuttaIni ────────────────────────────────────
    # runge_kutta_type=41 (4th-order Kutta) is a compile-time parameter;
    # no runtime inputs are needed — one case with empty inputs suffices.
    'test_MLRungeKuttaMod_RungeKuttaIni.exe': [
        {},   # no inputs; outputs are the Butcher tableau for 4th-order Kutta
    ],
    # ── MLCanopyTurbulenceMod :: GetPsiRSL ──────────────────────────────────
    'test_MLCanopyTurbulenceMod_GetPsiRSL.exe': [
        # za = hc: exact cancellation → psim = vkc/beta = 0.4/0.3, psic = 0
        {'za': 20.0, 'hc': 20.0, 'disp': 12.0, 'obu': -100.0,
         'beta': 0.3, 'PrSc': 0.5},
        # za above canopy, unstable
        {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu': -100.0,
         'beta': 0.3, 'PrSc': 0.5},
        # za above canopy, stable
        {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu':  100.0,
         'beta': 0.3, 'PrSc': 0.5},
        # neutral (large |obu|)
        {'za': 30.0, 'hc': 20.0, 'disp': 12.0, 'obu': 1000.0,
         'beta': 0.3, 'PrSc': 0.5},
        # different beta and PrSc
        {'za': 25.0, 'hc': 20.0, 'disp': 12.0, 'obu': -50.0,
         'beta': 0.25, 'PrSc': 0.6},
    ],
    # ── MLCanopyWaterMod :: CalcWettedFraction ──────────────────────────────
    'test_MLCanopyWaterMod_WettedFraction.exe': [
        # dry leaf
        {'h2ocan': 0.0,    'dpai': 1.0, 'dlai': 0.8},
        # very small h2ocan (below cap threshold)
        {'h2ocan': 1.0e-5, 'dpai': 1.0, 'dlai': 0.8},
        # at capacity (h2ocan = dewmx*dpai = 0.1): fwet hits cap 0.05
        {'h2ocan': 0.1,    'dpai': 1.0, 'dlai': 0.8},
        # above capacity: fwet still capped at 0.05
        {'h2ocan': 1.0,    'dpai': 1.0, 'dlai': 0.8},
        # dlai = dpai (all leaves): fdry = 1 - fwet
        {'h2ocan': 0.0,    'dpai': 0.5, 'dlai': 0.5},
        # dpai = 0: fwet = fdry = 0
        {'h2ocan': 0.0,    'dpai': 0.0, 'dlai': 0.0},
    ],
    # ── MLLeafHeatCapacityMod :: CalcLeafHeatCapacity ───────────────────────
    'test_MLLeafHeatCapacityMod_LeafHeatCapacity.exe': [
        {'sla': 0.04},   # typical broadleaf (~558 J/K/m2)
        {'sla': 0.02},   # needle-leaf  (higher LMA → larger cpleaf)
        {'sla': 0.06},   # grass/forb   (lower LMA → smaller cpleaf)
        {'sla': 0.01},   # very low SLA
        {'sla': 0.10},   # very high SLA
    ],
}


def generate(exe_name: str, cases: list) -> dict:
    """Run one executable over all its cases and return the golden record."""
    records = []
    for inputs in cases:
        outputs = run_fortran(exe_name, inputs)
        records.append({'inputs': inputs, 'outputs': outputs})
    return {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'executable': exe_name,
        'cases': records,
    }


def main():
    os.makedirs(GOLDEN_DIR, exist_ok=True)

    failed = []
    for exe, cases in CASES.items():
        print(f'  {exe} ... ', end='', flush=True)
        try:
            data = generate(exe, cases)
            out_path = os.path.join(GOLDEN_DIR, exe.replace('.exe', '.json'))
            with open(out_path, 'w') as f:
                json.dump(data, f, indent=2)
                f.write('\n')
            print(f'OK  ({len(cases)} cases → {os.path.basename(out_path)})')
        except Exception as exc:
            print(f'FAILED: {exc}')
            failed.append(exe)

    print()
    if failed:
        print(f'ERROR: {len(failed)} executable(s) failed:')
        for name in failed:
            print(f'  {name}')
        sys.exit(1)
    else:
        print(f'Golden data written to {os.path.abspath(GOLDEN_DIR)}/')
        print('Commit the JSON files alongside your code changes.')


if __name__ == '__main__':
    main()
