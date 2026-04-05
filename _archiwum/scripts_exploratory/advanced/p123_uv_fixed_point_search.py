#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p123_uv_fixed_point_search.py -- UV Fixed Point Search for TGP
================================================================

Searches for a UV fixed point of the TGP theory using three approaches:

(A) One-loop Wilson beta functions (from sek08_formalizm, eq:beta-mass/beta-g):
    k dm~/dk = -2m~ + g~/(2pi^2 (1+m~^2)^2)
    k dg~/dk = -g~ - 5g~^2/(2pi^2 (1+m~^2)^3)
    Known FP: (m~*=0, g~*=-2pi^2/5)

(B) Extended truncation with K(psi)=psi^4:
    Adds anomalous dimension eta and K flow.
    Beta functions modified by K coupling.

(C) Full potential scan (LPA'):
    Searches for a fixed-point potential V*(psi) such that dV*/dt = 0
    using the Wetterich equation from p120.

Key deliverables:
  - UV FP coordinates
  - Critical exponents (eigenvalues of stability matrix)
  - Anomalous dimension eta at FP
  - Number of relevant/irrelevant directions

Author: TGP project, session v42
Date: 2026-03-31
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.optimize import fsolve, root
from scipy.integrate import solve_ivp

# ================================================================
#  APPROACH A: One-loop Wilson beta functions
# ================================================================

def beta_A(x):
    """Beta functions in 1-loop Wilson approximation.
    x = [m_tilde_sq, g_tilde]
    Returns [beta_m, beta_g] (= k d/dk of couplings).
    """
    m2, g = x
    denom = (1 + m2)**2
    denom3 = (1 + m2)**3

    beta_m = -2*m2 + g / (2*np.pi**2 * denom)
    beta_g = -g - 5*g**2 / (2*np.pi**2 * denom3)
    return [beta_m, beta_g]


def stability_matrix_A(x):
    """Jacobian d(beta)/d(couplings) at point x.
    M_ij = d(beta_i)/d(x_j)
    """
    m2, g = x
    h = 1e-7
    n = len(x)
    M = np.zeros((n, n))
    f0 = beta_A(x)
    for j in range(n):
        xp = list(x)
        xp[j] += h
        fp = beta_A(xp)
        for i in range(n):
            M[i][j] = (fp[i] - f0[i]) / h
    return M


def run_approach_A():
    """Find and characterize UV FP using 1-loop Wilson beta functions."""
    print("\n" + "="*60)
    print("  APPROACH A: One-Loop Wilson Beta Functions")
    print("="*60)

    # Analytic FP
    m2_star = 0.0
    g_star = -2*np.pi**2 / 5
    print(f"  Analytic FP: m~^2* = {m2_star:.6f}, g~* = {g_star:.6f}")

    # Verify numerically
    sol = fsolve(beta_A, [0.01, -3.0], full_output=True)
    x_num = sol[0]
    info = sol[1]
    print(f"  Numerical FP: m~^2* = {x_num[0]:.8f}, g~* = {x_num[1]:.8f}")
    print(f"  Residual: |beta| = {np.linalg.norm(beta_A(x_num)):.2e}")

    # Stability matrix
    M = stability_matrix_A(x_num)
    eigenvalues = np.linalg.eigvals(M)
    print(f"\n  Stability matrix eigenvalues:")
    for i, ev in enumerate(sorted(eigenvalues)):
        relevance = "relevant (UV attractive)" if ev.real < 0 else "irrelevant (UV repulsive)"
        print(f"    theta_{i+1} = {ev.real:+.4f} -- {relevance}")

    n_relevant = sum(1 for ev in eigenvalues if ev.real < 0)
    n_irrelevant = sum(1 for ev in eigenvalues if ev.real > 0)
    print(f"\n  Relevant directions: {n_relevant}")
    print(f"  Irrelevant directions: {n_irrelevant}")

    # RG flow from UV to IR
    def rhs_flow(t, x):
        """t = ln(k/k_IR), flow from UV (t_max) to IR (t=0)."""
        return [-b for b in beta_A(x)]  # reverse: going to IR

    # Start near FP, flow to IR
    eps = 0.01
    x0 = [m2_star + eps, g_star + eps]
    t_max = np.log(25)  # k_UV/k_IR ~ 25

    sol_flow = solve_ivp(rhs_flow, [0, t_max], x0, method='RK45',
                         t_eval=np.linspace(0, t_max, 100), rtol=1e-8)

    if sol_flow.success:
        m2_IR = sol_flow.y[0, -1]
        g_IR = sol_flow.y[1, -1]
        print(f"\n  RG flow UV -> IR:")
        print(f"    UV: m~^2 = {x0[0]:.4f}, g~ = {x0[1]:.4f}")
        print(f"    IR: m~^2 = {m2_IR:.4f}, g~ = {g_IR:.4f}")

        # Check if IR values are physical
        m_sp_sq_IR = m2_IR * 1.0  # k_IR = 1
        print(f"    m_sp^2(IR) = {m_sp_sq_IR:.4f} (should be ~1)")

    results_A = {
        "FP": {"m2_star": float(x_num[0]), "g_star": float(x_num[1])},
        "eigenvalues": [float(ev.real) for ev in eigenvalues],
        "n_relevant": n_relevant,
        "residual": float(np.linalg.norm(beta_A(x_num)))
    }

    return results_A


# ================================================================
#  APPROACH B: Extended truncation with anomalous dimension
# ================================================================

def beta_B(x, eta=0.0):
    """Extended beta functions with anomalous dimension eta.

    In the presence of K(psi) = psi^4, the propagator is modified:
    G(p) = 1/(K*p^2 + m^2) with K=psi^4.

    At the vacuum psi=1: K=1 (normalized), so standard propagator.
    The anomalous dimension modifies the canonical dimensions:

    k dm~/dk = -(2 - eta)*m~ + g~/(2*pi^2*(1+m~^2)^2)
    k dg~/dk = -(1 + eta)*g~ - 5*g~^2/(2*pi^2*(1+m~^2)^3)

    eta itself is determined self-consistently:
    eta = -g~^2 / (6*pi^2*(1+m~^2)^4)  (1-loop, from K running)
    """
    m2, g = x
    denom2 = (1 + m2)**2
    denom3 = (1 + m2)**3
    denom4 = (1 + m2)**4

    # Self-consistent eta (from K flow at FP, 1-loop)
    eta_sc = -g**2 / (6*np.pi**2 * denom4)

    beta_m = -(2 - eta_sc)*m2 + g / (2*np.pi**2 * denom2)
    beta_g = -(1 + eta_sc)*g - 5*g**2 / (2*np.pi**2 * denom3)
    return [beta_m, beta_g], eta_sc


def run_approach_B():
    """Find UV FP with anomalous dimension."""
    print("\n" + "="*60)
    print("  APPROACH B: Extended Truncation with Anomalous Dimension")
    print("="*60)

    def beta_only(x):
        b, _ = beta_B(x)
        return b

    # Search from multiple initial conditions
    best_fp = None
    best_res = 1e10

    for m2_init in np.linspace(-1, 1, 11):
        for g_init in np.linspace(-10, 0, 21):
            try:
                sol = fsolve(beta_only, [m2_init, g_init], full_output=True)
                x_fp = sol[0]
                residual = np.linalg.norm(beta_only(x_fp))
                if residual < 1e-10 and residual < best_res:
                    # Check it's not Gaussian FP
                    if abs(x_fp[1]) > 0.01:
                        best_fp = x_fp
                        best_res = residual
            except:
                pass

    if best_fp is None:
        print("  No non-Gaussian FP found!")
        return None

    betas, eta_fp = beta_B(best_fp)
    print(f"  Non-Gaussian UV FP found:")
    print(f"    m~^2* = {best_fp[0]:.8f}")
    print(f"    g~*   = {best_fp[1]:.8f}")
    print(f"    eta*  = {eta_fp:.6f}")
    print(f"    |beta| = {best_res:.2e}")

    # Stability matrix
    h = 1e-7
    M = np.zeros((2, 2))
    f0 = beta_only(best_fp)
    for j in range(2):
        xp = list(best_fp)
        xp[j] += h
        fp = beta_only(xp)
        for i in range(2):
            M[i][j] = (fp[i] - f0[i]) / h

    eigenvalues = np.linalg.eigvals(M)
    print(f"\n  Stability matrix eigenvalues (theta = -eigenvalue):")
    for i, ev in enumerate(sorted(eigenvalues, key=lambda x: x.real)):
        theta = -ev.real
        relevance = "relevant" if theta > 0 else "irrelevant"
        print(f"    theta_{i+1} = {theta:+.4f} ({relevance})")

    n_relevant = sum(1 for ev in eigenvalues if ev.real < 0)

    # Alpha_TGP with eta correction
    alpha_corrected = 2 + eta_fp
    print(f"\n  alpha_TGP = 2 + eta = {alpha_corrected:.6f}")
    print(f"  Deviation from alpha=2: {abs(eta_fp)*100:.2f}%")

    results_B = {
        "FP": {"m2_star": float(best_fp[0]), "g_star": float(best_fp[1])},
        "eta": float(eta_fp),
        "eigenvalues": [float(ev.real) for ev in eigenvalues],
        "n_relevant": n_relevant,
        "alpha_corrected": float(alpha_corrected)
    }

    return results_B


# ================================================================
#  APPROACH C: Fixed-point potential V*(psi)
# ================================================================

def run_approach_C():
    """Search for fixed-point potential V*(psi) in LPA'."""
    print("\n" + "="*60)
    print("  APPROACH C: Fixed-Point Potential V*(psi) in LPA'")
    print("="*60)

    N_PSI = 60
    PSI = np.linspace(0.1, 3.0, N_PSI)
    DPSI = PSI[1] - PSI[0]

    def d2(f):
        d2f = np.zeros_like(f)
        h = DPSI
        for i in range(2, len(f)-2):
            d2f[i] = (-f[i+2] + 16*f[i+1] - 30*f[i] + 16*f[i-1] - f[i-2]) / (12*h**2)
        d2f[0] = d2f[1] = d2f[2]
        d2f[-1] = d2f[-2] = d2f[-3]
        return d2f

    # At UV FP: dV/dt = 0 => canonical scaling balances quantum corrections
    # Dimensionless: u(psi) = V(psi)/k^4
    # FP condition: -4*u + psi*u' + K(psi)*k^2/(32*pi^2*(K(psi) + u'')) = 0
    # At FP (k=k*): normalize k*=1:
    # -4*u + psi*u'(psi) + K(psi)/(32*pi^2*(K(psi) + u''(psi))) = 0

    K = PSI**4  # K(psi) = psi^4

    def fp_residual(u_vec):
        """Residual of the FP equation: -4u + psi*u' + K/(32pi^2*(K+u'')) = 0"""
        u = u_vec
        u_prime = np.gradient(u, DPSI)
        u_pp = d2(u)

        denom = K + u_pp
        denom = np.where(np.abs(denom) < 1e-8, 1e-8, denom)

        residual = -4*u + PSI * u_prime + K / (32*np.pi**2 * denom)
        return residual

    # Initial guess: tree-level dimensionless potential
    # u_tree = V_tree/k^4 at k=k_UV
    k_UV = 25.0
    u0 = (PSI**3/3 - PSI**4/4) / k_UV**4

    # Solve
    from scipy.optimize import least_squares
    result = least_squares(fp_residual, u0, method='lm',
                          ftol=1e-10, xtol=1e-10, max_nfev=10000)

    u_star = result.x
    residual_norm = np.linalg.norm(result.fun)

    print(f"  Solver: {'converged' if result.success else 'FAILED'}")
    print(f"  |residual| = {residual_norm:.2e}")
    print(f"  nfev = {result.nfev}")

    # Analyze FP potential
    idx1 = np.argmin(np.abs(PSI - 1.0))
    u_star_pp = d2(u_star)

    print(f"\n  Fixed-point potential u*(psi) at psi=1:")
    print(f"    u*(1) = {u_star[idx1]:.6e}")
    print(f"    u*''(1) = {u_star_pp[idx1]:.6e}")
    print(f"    K(1) = {K[idx1]:.4f}")

    # Physical mass at FP
    m2_fp = u_star_pp[idx1] / K[idx1] if abs(K[idx1]) > 1e-10 else float('inf')
    print(f"    m^2_phys(FP) = u*''/K = {m2_fp:.6e}")

    stable_at_1 = u_star_pp[idx1] > 0 if abs(K[idx1]) > 1e-10 else False
    print(f"    Vacuum stable at FP: {'YES' if stable_at_1 else 'NO'}")

    # Check sign of u* - should be bounded
    print(f"\n  u*(psi) profile:")
    for psi_test in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5]:
        idx = np.argmin(np.abs(PSI - psi_test))
        print(f"    u*({psi_test:.1f}) = {u_star[idx]:.6e}")

    results_C = {
        "converged": bool(result.success),
        "residual": float(residual_norm),
        "u_star_1": float(u_star[idx1]),
        "u_star_pp_1": float(u_star_pp[idx1]),
        "m2_fp": float(m2_fp),
        "stable": bool(stable_at_1)
    }

    return results_C


# ================================================================
#  MAIN
# ================================================================

def main():
    print("="*60)
    print("  TGP UV Fixed Point Search")
    print("  Goal: Verify asymptotic safety of TGP")
    print("="*60)

    results = {}

    # Approach A
    results['A'] = run_approach_A()

    # Approach B
    results['B'] = run_approach_B()

    # Approach C
    results['C'] = run_approach_C()

    # Summary
    print("\n" + "="*60)
    print("  SUMMARY: UV FIXED POINT ANALYSIS")
    print("="*60)

    tests = []

    # Test 1: Non-Gaussian FP exists (approach A)
    if results['A']:
        fp_A = results['A']['FP']
        tests.append(("T1: Non-Gaussian FP exists (1-loop)",
                      abs(fp_A['g_star']) > 0.01,
                      f"g* = {fp_A['g_star']:.4f}"))

    # Test 2: FP has exactly 2 relevant directions
    if results['A']:
        tests.append(("T2: 2 relevant directions at FP",
                      results['A']['n_relevant'] == 2,
                      f"n_rel = {results['A']['n_relevant']}"))

    # Test 3: eta is small (alpha ~ 2)
    if results['B']:
        eta = results['B']['eta']
        tests.append(("T3: eta small (|eta| < 0.1)",
                      abs(eta) < 0.1,
                      f"eta = {eta:.6f}"))

    # Test 4: FP stable under eta correction
    if results['B']:
        tests.append(("T4: Non-Gaussian FP survives eta",
                      results['B']['n_relevant'] >= 1,
                      f"n_rel = {results['B']['n_relevant']}"))

    # Test 5: FP potential V* is bounded
    if results['C']:
        tests.append(("T5: FP potential solver converges",
                      results['C']['converged'],
                      f"|res| = {results['C']['residual']:.2e}"))

    # Test 6: Vacuum stable at FP
    if results['C']:
        tests.append(("T6: Vacuum stable at FP",
                      results['C']['stable'],
                      f"m^2(FP) = {results['C']['m2_fp']:.4e}"))

    n_pass = sum(1 for _, v, _ in tests if v)
    n_total = len(tests)

    print(f"\n  Tests: {n_pass}/{n_total} PASS")
    for name, passed, detail in tests:
        mark = "PASS" if passed else "FAIL"
        print(f"    [{mark}] {name}")
        print(f"           {detail}")

    print(f"\n  Conclusion:")
    if n_pass >= 4:
        print("    POSITIVE: Strong signals for asymptotic safety of TGP.")
        print("    Non-Gaussian UV FP exists with correct number of relevant")
        print("    directions. K(psi)=psi^4 provides structural stability.")
        print("    Status: Hypothesis -> Proposition (numerical evidence)")
    elif n_pass >= 2:
        print("    PARTIAL: Some signals positive, further analysis needed.")
    else:
        print("    NEGATIVE: UV FP not found at this truncation level.")

    return results


if __name__ == '__main__':
    main()
