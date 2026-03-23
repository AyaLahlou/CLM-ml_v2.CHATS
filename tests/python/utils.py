"""
Shared utilities for CLM-ml Fortran functional tests.

Each Fortran test driver reads a Fortran namelist from stdin and writes
key=value pairs to stdout.  This module provides:

  run_fortran(exe_name, inputs)
      Build a namelist string, pipe it to the Fortran executable, and
      return a dict mapping output key names to float values.

  run_fortran_array(exe_name, inputs, keys)
      Like run_fortran, but also collects indexed array outputs of the
      form  u_1=..., u_2=..., u_3=...  and assembles them into a list.
"""

import os
import re
import subprocess

# Directory containing the compiled test executables
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')


def _build_namelist(inputs: dict) -> str:
    """Convert a dict to a Fortran namelist string piped to stdin."""
    lines = ['&inputs']
    for key, val in inputs.items():
        if isinstance(val, (list, tuple)):
            # Format array: key = v1, v2, v3
            vals = ', '.join(f'{v:.17g}' if isinstance(v, float) else str(v)
                             for v in val)
            lines.append(f'  {key} = {vals}')
        elif isinstance(val, float):
            lines.append(f'  {key} = {val:.17g}')
        else:
            lines.append(f'  {key} = {val}')
    lines.append('/')
    return '\n'.join(lines) + '\n'


def run_fortran(exe_name: str, inputs: dict) -> dict:
    """
    Run a Fortran test executable with the given inputs.

    Parameters
    ----------
    exe_name : str
        Executable file name, e.g. 'test_SatVap.exe'.  Must exist in
        tests/bin/ (built by ``make all`` in tests/).
    inputs : dict
        Namelist variable names → values.  Floats, ints, and lists are
        all supported.

    Returns
    -------
    dict
        Parsed output: key → float.  Only lines of the form
        ``KEY=VALUE`` (no whitespace in key) are collected; all other
        output lines (Fortran warnings written via iulog) are ignored.

    Raises
    ------
    RuntimeError
        If the executable exits with a non-zero return code.
    FileNotFoundError
        If the executable does not exist in tests/bin/.
    """
    exe_path = os.path.join(BIN_DIR, exe_name)
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(
            f"Executable not found: {exe_path}\n"
            f"Run 'make all' in the tests/ directory first."
        )

    namelist = _build_namelist(inputs)
    result = subprocess.run(
        [exe_path],
        input=namelist,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"{exe_name} exited with code {result.returncode}.\n"
            f"--- stdin ---\n{namelist}\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )

    outputs = {}
    # Match lines like:  r1= 1.2345678901234567E+00
    pattern = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)=\s*(.+)$')
    for line in result.stdout.splitlines():
        m = pattern.match(line.strip())
        if m:
            try:
                outputs[m.group(1)] = float(m.group(2))
            except ValueError:
                pass  # Non-numeric output line (e.g. Fortran warning text)

    return outputs


def run_fortran_array(exe_name: str, inputs: dict, base_key: str) -> list:
    """
    Run a Fortran test executable and collect indexed array output.

    Fortran array outputs are written as::

        u_1= 1.0000...
        u_2= 2.0000...
        u_3= 3.0000...

    This function assembles them into a Python list [u_1, u_2, u_3, ...].

    Parameters
    ----------
    exe_name : str
        Executable name.
    inputs : dict
        Namelist inputs (must include ``n`` to indicate array length).
    base_key : str
        The base array name (e.g. 'u' to collect u_1, u_2, ...).

    Returns
    -------
    list of float
        Array values in order.
    """
    outputs = run_fortran(exe_name, inputs)
    n = int(inputs.get('n', 0))
    result = []
    for i in range(1, n + 1):
        key = f'{base_key}_{i}'
        if key not in outputs:
            raise KeyError(
                f"Expected output key '{key}' not found in {exe_name} output.\n"
                f"Available keys: {list(outputs.keys())}"
            )
        result.append(outputs[key])
    return result
