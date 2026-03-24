"""
Functional tests for MLMathToolsMod :: tridiag_2eq.

Solves two coupled tridiagonal systems simultaneously:

  a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
  a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)

where T and q are two unknowns at each layer (e.g. air temperature and
water-vapour mole fraction in a canopy model).

Note: a*(1) and c*(n) are not referenced.
"""

import pytest
from utils import run_fortran

EXE = 'test_MLMathToolsMod_tridiag_2eq.exe'


def _run(n, a1, b11, b12, c1, d1, a2, b21, b22, c2, d2):
    """Helper: build inputs dict and run the executable."""
    return run_fortran(EXE, {
        'n': n,
        'a1': a1, 'b11': b11, 'b12': b12, 'c1': c1, 'd1': d1,
        'a2': a2, 'b21': b21, 'b22': b22, 'c2': c2, 'd2': d2,
    })


def _t(out, n):
    return [out[f't_{i}'] for i in range(1, n + 1)]


def _q(out, n):
    return [out[f'q_{i}'] for i in range(1, n + 1)]


# ──────────────────────────────────────────────────────────────────────────────
# Single-layer (n=1) tests — reduces to a 2×2 linear system
# ──────────────────────────────────────────────────────────────────────────────

class TestSingleLayer:

    def test_decoupled_diagonal(self):
        """n=1, no coupling (b12=b21=0): T=d1/b11, q=d2/b22."""
        out = _run(1,
                   a1=[0.0], b11=[2.0], b12=[0.0], c1=[0.0], d1=[4.0],
                   a2=[0.0], b21=[0.0], b22=[3.0], c2=[0.0], d2=[6.0])
        assert abs(out['t_1'] - 2.0) < 1e-12
        assert abs(out['q_1'] - 2.0) < 1e-12

    def test_coupled_2x2_system(self):
        """n=1, coupled 2×2: [2 1; 1 3]*[T;q] = [5;7] → T=1.6, q=1.8.
        Solution: det=5, T=(3*5-1*7)/5=8/5, q=(2*7-1*5)/5=9/5."""
        out = _run(1,
                   a1=[0.0], b11=[2.0], b12=[1.0], c1=[0.0], d1=[5.0],
                   a2=[0.0], b21=[1.0], b22=[3.0], c2=[0.0], d2=[7.0])
        assert abs(out['t_1'] - 1.6) < 1e-12
        assert abs(out['q_1'] - 1.8) < 1e-12

    def test_identity_system(self):
        """n=1, identity matrix: [1 0; 0 1]*[T;q] = [d1;d2] → T=d1, q=d2."""
        out = _run(1,
                   a1=[0.0], b11=[1.0], b12=[0.0], c1=[0.0], d1=[7.0],
                   a2=[0.0], b21=[0.0], b22=[1.0], c2=[0.0], d2=[3.0])
        assert abs(out['t_1'] - 7.0) < 1e-12
        assert abs(out['q_1'] - 3.0) < 1e-12


# ──────────────────────────────────────────────────────────────────────────────
# Multi-layer tests
# ──────────────────────────────────────────────────────────────────────────────

class TestMultiLayer:

    def test_two_layer_decoupled_tridiag(self):
        """n=2, no cross-coupling: two independent tridiag systems.

        Eq 1 (T):  [2 -1; -1 2]*[T1;T2] = [1;1]  → T1=T2=1
        Eq 2 (q):  [3 -1; -1 3]*[q1;q2] = [2;2]  → q1=q2=1
        """
        out = _run(2,
                   a1=[0.0, -1.0], b11=[2.0, 2.0], b12=[0.0, 0.0],
                   c1=[-1.0, 0.0], d1=[1.0, 1.0],
                   a2=[0.0, -1.0], b21=[0.0, 0.0], b22=[3.0, 3.0],
                   c2=[-1.0, 0.0], d2=[2.0, 2.0])
        assert abs(out['t_1'] - 1.0) < 1e-12
        assert abs(out['t_2'] - 1.0) < 1e-12
        assert abs(out['q_1'] - 1.0) < 1e-12
        assert abs(out['q_2'] - 1.0) < 1e-12

    def test_two_layer_uniform_coupled_diagonal(self):
        """n=2, coupled but no tridiag off-diagonals (a=c=0):
        same 2×2 system at each layer independently.
        [2 1; 1 3]*[T;q] = [5;7] → T=1.6, q=1.8 at each layer.
        """
        out = _run(2,
                   a1=[0.0, 0.0], b11=[2.0, 2.0], b12=[1.0, 1.0],
                   c1=[0.0, 0.0], d1=[5.0, 5.0],
                   a2=[0.0, 0.0], b21=[1.0, 1.0], b22=[3.0, 3.0],
                   c2=[0.0, 0.0], d2=[7.0, 7.0])
        for i in [1, 2]:
            assert abs(out[f't_{i}'] - 1.6) < 1e-12, f"T layer {i}"
            assert abs(out[f'q_{i}'] - 1.8) < 1e-12, f"q layer {i}"

    def test_three_layer_decoupled_identity(self):
        """n=3, identity system (no coupling, no off-diag): T=d1, q=d2."""
        d1_vals = [2.0, 4.0, 6.0]
        d2_vals = [1.0, 3.0, 5.0]
        out = _run(3,
                   a1=[0.0]*3, b11=[1.0]*3, b12=[0.0]*3,
                   c1=[0.0]*3, d1=d1_vals,
                   a2=[0.0]*3, b21=[0.0]*3, b22=[1.0]*3,
                   c2=[0.0]*3, d2=d2_vals)
        for i, (exp_t, exp_q) in enumerate(zip(d1_vals, d2_vals), start=1):
            assert abs(out[f't_{i}'] - exp_t) < 1e-12
            assert abs(out[f'q_{i}'] - exp_q) < 1e-12

    def test_solution_satisfies_equations(self):
        """Verify the solution satisfies the original coupled system.

        3-layer system with both tridiag coupling and cross-coupling.
        Check: A1*T + B11*T + B12*q + C1*T = D1 (row-wise).
        """
        n = 3
        a1  = [0.0, -0.5, -0.5]
        b11 = [3.0,  3.0,  3.0]
        b12 = [0.5,  0.5,  0.5]
        c1  = [-0.5, -0.5, 0.0]
        d1  = [2.0,  1.0,  2.0]
        a2  = [0.0, -0.3, -0.3]
        b21 = [0.3,  0.3,  0.3]
        b22 = [2.0,  2.0,  2.0]
        c2  = [-0.3, -0.3, 0.0]
        d2  = [1.5,  0.8,  1.5]

        out = _run(n, a1, b11, b12, c1, d1, a2, b21, b22, c2, d2)
        T = _t(out, n)
        q = _q(out, n)

        for i in range(n):  # 0-based indexing for Python lists
            # Equation 1 residual at layer i+1
            t_prev = T[i - 1] if i > 0     else 0.0
            t_next = T[i + 1] if i < n - 1 else 0.0
            q_prev = q[i - 1] if i > 0     else 0.0
            q_next = q[i + 1] if i < n - 1 else 0.0

            res1 = (a1[i] * t_prev
                    + b11[i] * T[i]
                    + b12[i] * q[i]
                    + c1[i] * t_next
                    - d1[i])
            res2 = (a2[i] * q_prev
                    + b21[i] * T[i]
                    + b22[i] * q[i]
                    + c2[i] * q_next
                    - d2[i])

            assert abs(res1) < 1e-11, f"Layer {i+1} eq1 residual: {res1}"
            assert abs(res2) < 1e-11, f"Layer {i+1} eq2 residual: {res2}"
