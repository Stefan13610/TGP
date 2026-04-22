#!/usr/bin/env python3
"""
ode_koide_formA_v47b.py -- Test Koide Q_K=3/2 for canonical Form A.

KEY QUESTION: Does the canonical TGP ODE (K=g^4, alpha=2, source=g^2(1-g))
also produce Q_K = 3/2 via phi-FP spacing?

Both Form A and B have oscillatory tails (omega=1).
Test if Form A gives the same Koide property.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21 = m_mu / m_e; r_31 = m_tau / m_e

def make_solver(alpha_coeff, source_type):
    """Create an ODE solver with given alpha and source type."""
    def solver(g0, r_max=300):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            if source_type == "linear":
                source = 1.0 - g
            elif source_type == "quadratic":
                source = g**2 * (1.0 - g)
            cross = (alpha_coeff / g) * gp**2
            if r < 1e-10:
                return [gp, (source - cross) / 3.0]
            return [gp, source - cross - 2.0 * gp / r]
        sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                        rtol=1e-11, atol=1e-13, max_step=0.02)
        return sol.t, sol.y[0]
    return solver

def A_tail(solver, g0):
    r, g = solver(g0)
    mask = (r > 50) & (r < 200)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)

def test_form(name, solver):
    """Full Koide test for a given ODE form."""
    print(f"\n  {name}")
    print(f"  {'='*60}")

    # Calibrate g0_e from r_21
    def r21_res(g0_1):
        A1 = A_tail(solver, g0_1)
        A2 = A_tail(solver, PHI * g0_1)
        if A1 < 1e-15: return 1e10
        return (A2/A1)**4 - r_21

    # Scan for sign change
    g0_vals = np.linspace(0.3, 1.5, 50)
    resids = [r21_res(g) for g in g0_vals]
    bracket = None
    for i in range(len(resids)-1):
        if resids[i] * resids[i+1] < 0:
            bracket = (g0_vals[i], g0_vals[i+1])
            break

    if bracket is None:
        print(f"    No r_21 solution found!")
        # Print scan
        for g0, res in zip(g0_vals[::5], resids[::5]):
            A1 = A_tail(solver, g0)
            A2 = A_tail(solver, PHI*g0)
            r = (A2/A1)**4 if A1 > 1e-15 else 0
            print(f"      g0={g0:.3f}: A1={A1:.6f}, A2={A2:.6f}, (A2/A1)^4={r:.1f}")
        return None

    g0_e = brentq(r21_res, bracket[0], bracket[1], xtol=1e-8)
    g0_mu = PHI * g0_e
    A_e = A_tail(solver, g0_e)
    A_mu = A_tail(solver, g0_mu)
    r21_check = (A_mu/A_e)**4

    print(f"    g0_e = {g0_e:.8f}")
    print(f"    g0_mu = {g0_mu:.8f}")
    print(f"    A_e = {A_e:.10f}, A_mu = {A_mu:.10f}")
    print(f"    r_21 = {r21_check:.3f} (target: {r_21:.3f})")

    # Find Koide g0_tau
    def koide_res(g0_3):
        A1 = A_tail(solver, g0_e)
        A2 = A_tail(solver, g0_mu)
        A3 = A_tail(solver, g0_3)
        if min(A1,A2,A3) < 1e-15: return 1e10
        S2 = A1**2 + A2**2 + A3**2
        S4 = A1**4 + A2**4 + A3**4
        return S2**2/S4 - 1.5

    # Scan for Koide solution -- FINE grid needed for Form A
    # (Form A has narrow bracket near critical g0 ~ 1.60)
    g0_tau_vals = np.linspace(g0_mu + 0.01, 2.5, 200)
    k_resids = [koide_res(g) for g in g0_tau_vals]
    k_bracket = None
    for i in range(len(k_resids)-1):
        if k_resids[i] != 1e10 and k_resids[i+1] != 1e10:
            if k_resids[i] * k_resids[i+1] < 0:
                k_bracket = (g0_tau_vals[i], g0_tau_vals[i+1])
                break

    if k_bracket:
        g0_tau = brentq(koide_res, k_bracket[0], k_bracket[1], xtol=1e-10)
        A_tau = A_tail(solver, g0_tau)
        r31_K = (A_tau/A_e)**4

        S2 = A_e**2 + A_mu**2 + A_tau**2
        S4 = A_e**4 + A_mu**4 + A_tau**4
        Q_K = S2**2/S4

        c_K = g0_tau / g0_e
        print(f"    g0_tau(Koide) = {g0_tau:.8f}")
        print(f"    A_tau = {A_tau:.10f}")
        print(f"    Q_K = {Q_K:.10f} (target: 1.5)")
        print(f"    r_31 = {r31_K:.1f} (PDG: {r_31:.1f})")
        print(f"    c = g0_tau/g0_e = {c_K:.8f}")
        print(f"    Dev from PDG r_31: {(r31_K/r_31 - 1)*100:+.4f}%")
        return {'g0_e': g0_e, 'g0_tau': g0_tau, 'Q_K': Q_K, 'r31': r31_K, 'c': c_K}
    else:
        print(f"    No Koide solution found!")
        # Show Q_K at phi^2 spacing
        g0_tau_p2 = PHI**2 * g0_e
        A_tau_p2 = A_tail(solver, g0_tau_p2)
        if A_tau_p2 > 1e-15:
            S2 = A_e**2 + A_mu**2 + A_tau_p2**2
            S4 = A_e**4 + A_mu**4 + A_tau_p2**4
            Q = S2**2/S4
            print(f"    Q_K at phi^2 spacing: {Q:.6f}")
        return None


print("=" * 72)
print("KOIDE TEST: CANONICAL FORM A vs FORM B")
print("=" * 72)

# Define solvers
solver_A = make_solver(alpha_coeff=2.0, source_type="quadratic")
solver_B = make_solver(alpha_coeff=1.0, source_type="linear")

# Also test intermediate cases
solver_A15 = make_solver(alpha_coeff=1.5, source_type="quadratic")
solver_B2 = make_solver(alpha_coeff=2.0, source_type="linear")

result_A = test_form("FORM A: alpha=2, source=g^2(1-g) [canonical K=g^4]", solver_A)
result_B = test_form("FORM B: alpha=1, source=(1-g) [current scripts]", solver_B)
result_B2 = test_form("FORM B2: alpha=2, source=(1-g) [hybrid]", solver_B2)
result_A15 = test_form("FORM A1.5: alpha=1.5, source=g^2(1-g)", solver_A15)


# Summary
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

forms = [("A (canonical)", result_A), ("B (scripts)", result_B),
         ("B2 (hybrid)", result_B2), ("A1.5", result_A15)]

print(f"\n  {'Form':<20} {'g0_e':>10} {'c=g0t/g0e':>12} {'Q_K':>12} {'r31':>10} {'r31 dev':>10}")
for name, res in forms:
    if res:
        print(f"  {name:<20} {res['g0_e']:10.6f} {res['c']:12.8f} "
              f"{res['Q_K']:12.8f} {res['r31']:10.1f} "
              f"{(res['r31']/r_31 - 1)*100:+10.4f}%")
    else:
        print(f"  {name:<20} {'FAILED':>10}")

print(f"\n  PDG reference: r_21 = {r_21:.3f}, r_31 = {r_31:.2f}")
print(f"  Koide target: Q_K = 1.500000")
