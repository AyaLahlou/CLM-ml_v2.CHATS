"""
Golden (trusted-run) regression tests for all Fortran subroutines.

For each executable, tests/golden/<exe>.json holds the inputs and outputs
captured from a verified, trusted build of the code.  These tests re-run
every case with the same inputs and verify that the outputs match the golden
values within a tight tolerance (1e-14 relative, 1e-300 absolute – effectively
bit-for-bit agreement on doubles).

If an intentional code change alters numerical output, regenerate the golden
files with:

    cd tests/
    python python/generate_golden.py

then commit the updated JSON files alongside the code change.
"""

import json
import os

import pytest

# Allow running directly from tests/python/ or from tests/
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import run_fortran  # noqa: E402

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'golden')

# Tolerance for regression comparison.
# Chosen tight enough to catch any meaningful numerical drift while still
# allowing benign last-bit differences across identical compiler flags.
REL_TOL = 1e-14
ABS_TOL = 1e-300   # near zero guard


def _load_all_cases():
    """
    Collect (exe_name, case_index, inputs, expected_outputs) tuples from
    every golden JSON file found in tests/golden/.

    Returns an empty list (no failures) when no JSON files exist yet –
    the test module then has zero parametrised test items and pytest reports
    them as "no tests collected" rather than an error.
    """
    params = []
    if not os.path.isdir(GOLDEN_DIR):
        return params

    for fname in sorted(os.listdir(GOLDEN_DIR)):
        if not fname.endswith('.json'):
            continue
        path = os.path.join(GOLDEN_DIR, fname)
        with open(path) as f:
            data = json.load(f)
        exe = data['executable']
        for idx, case in enumerate(data['cases']):
            params.append(pytest.param(
                exe,
                case['inputs'],
                case['outputs'],
                id=f'{exe}[{idx}]',
            ))
    return params


_ALL_CASES = _load_all_cases()

# ---------------------------------------------------------------------------
# Skip the entire module gracefully if golden files have not been generated
# yet (first-time checkout before running generate_golden.py).
# ---------------------------------------------------------------------------
if not _ALL_CASES:
    pytest.skip(
        'No golden data found in tests/golden/.  '
        'Run  python tests/python/generate_golden.py  first.',
        allow_module_level=True,
    )


@pytest.mark.parametrize('exe, inputs, expected', _ALL_CASES)
def test_golden_regression(exe, inputs, expected):
    """
    Re-run one Fortran test case and compare every output against its
    golden (trusted-run) value.
    """
    actual = run_fortran(exe, inputs)

    for key, golden_val in expected.items():
        assert key in actual, (
            f'{exe}: output key "{key}" missing.\n'
            f'  inputs   = {inputs}\n'
            f'  got keys = {sorted(actual.keys())}'
        )
        got = actual[key]
        # Relative-then-absolute tolerance (mirrors math.isclose semantics)
        diff = abs(got - golden_val)
        tol  = REL_TOL * abs(golden_val) + ABS_TOL
        assert diff <= tol, (
            f'{exe}: output "{key}" differs from golden value.\n'
            f'  inputs  = {inputs}\n'
            f'  got     = {got!r}\n'
            f'  golden  = {golden_val!r}\n'
            f'  |diff|  = {diff:.3e}  (tol={tol:.3e})'
        )
