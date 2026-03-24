"""
Functional tests for MLMathToolsMod :: tridiag.

Solves the tridiagonal system F*u = r where F is defined by vectors a, b, c.
"""

import pytest
from utils import run_fortran, run_fortran_array

EXE = 'test_MLMathToolsMod_tridiag.exe'


class TestTridiag:
    """Tests for MLMathToolsMod :: tridiag(a, b, c, r, u, n)."""

    def test_3x3_symmetric(self):
        """Solve the 3x3 tridiagonal system with known solution u=[1,1,1].

        Matrix:  diag=[2,2,2],  off-diag=[-1,-1]
        RHS:     r=[1,0,1]
        Solution: u=[1,1,1]  (verify: 2-1=1, -1+2-1=0, -1+2=1)
        """
        inputs = {
            'n': 3,
            'a': [0.0, -1.0, -1.0],   # a(1) not referenced
            'b': [2.0,  2.0,  2.0],
            'c': [-1.0, -1.0, 0.0],   # c(3) not referenced
            'r': [1.0,  0.0,  1.0],
        }
        u = run_fortran_array(EXE, inputs, 'u')
        for val in u:
            assert abs(val - 1.0) < 1e-12

    def test_diagonal_system(self):
        """Pure diagonal system (a=c=0): solution is r_i/b_i."""
        inputs = {
            'n': 4,
            'a': [0.0, 0.0, 0.0, 0.0],
            'b': [2.0, 3.0, 4.0, 5.0],
            'c': [0.0, 0.0, 0.0, 0.0],
            'r': [6.0, 9.0, 12.0, 15.0],
        }
        u = run_fortran_array(EXE, inputs, 'u')
        expected = [3.0, 3.0, 3.0, 3.0]
        for i, (val, exp) in enumerate(zip(u, expected)):
            assert abs(val - exp) < 1e-12, f"u[{i+1}]: got {val}, expected {exp}"

    def test_identity_rhs(self):
        """Identity matrix (b=1, a=c=0): solution equals rhs."""
        r_vals = [7.0, 3.0, 5.0]
        inputs = {
            'n': 3,
            'a': [0.0, 0.0, 0.0],
            'b': [1.0, 1.0, 1.0],
            'c': [0.0, 0.0, 0.0],
            'r': r_vals,
        }
        u = run_fortran_array(EXE, inputs, 'u')
        for val, exp in zip(u, r_vals):
            assert abs(val - exp) < 1e-12
