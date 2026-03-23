"""
pytest configuration for CLM-ml Fortran functional test suite.

Automatically builds all test executables once per pytest session by
running ``make all`` in the tests/ directory before any tests run.

Usage (from tests/python/):
    python -m pytest -v
    python -m pytest test_water_vapor.py -v
    python -m pytest -k "satvap" -v

If the build fails, the entire session fails with a clear error message
showing the compiler output.
"""

import os
import subprocess
import sys

import pytest


@pytest.fixture(scope="session", autouse=True)
def build_test_executables():
    """Build all Fortran test executables before the test session starts."""
    tests_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    bin_dir = os.path.join(tests_dir, 'bin')

    result = subprocess.run(
        ['make', 'all'],
        cwd=tests_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        pytest.fail(
            "Fortran test build failed.\n\n"
            "Common causes:\n"
            "  - Main model not built yet: run 'make' in offline_executable/ first\n"
            "  - NETCDF_PATH not set in your environment\n"
            "  - nvfortran not on PATH\n\n"
            f"--- make stdout ---\n{result.stdout}\n"
            f"--- make stderr ---\n{result.stderr}"
        )

    # Verify at least one executable was produced
    exes = [f for f in os.listdir(bin_dir) if f.endswith('.exe')]
    if not exes:
        pytest.fail(
            f"Make succeeded but no .exe files found in {bin_dir}.\n"
            f"--- make stdout ---\n{result.stdout}"
        )
