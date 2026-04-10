#!/usr/bin/env python3
"""
ode_koide_formA_exact_v47b.py -- Find exact Koide solution for canonical Form A.

DISCOVERY: Form A has a NARROW bracket for Koide near g0_tau ~ 1.57,
just below the critical g0 ~ 1.60 where solitons collapse.
The original scan (linspace 1.2-3.0, 50 pts) missed it.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21 = m_mu / m_e; r_31 = m_tau / m_e


def solver_A(g0, r_max=300):
    """Canonical Form A: alpha=2, source=g^2(1-g)."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-11, atol=1e-13, max_step=0.02)
    return sol.t, sol.y[0]


def solver_B(g0, r_max=300):
    """Form B: alpha=1, source=(1-g)."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = 1.0 - g
        cross = (1.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-11, atol=1e-13, max_step=0.02)
    return sol.t, sol.y[0]


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


def find_koide(name, solver):
    """Full Koide analysis for given ODE form."""
    print("=" * 65)
    print("  %s" % name)
    print("=" * 65)

    # 1. Calibrate g0_e from r_21
    def r21_res(g0_1):
        A1 = A_tail(solver, g0_1)
        A2 = A_tail(solver, PHI * g0_1)
        if A1 < 1e-15:
            return 1e10
        return (A2 / A1)**4 - r_21

    g0_e = brentq(r21_res, 0.8, 0.95, xtol=1e-10)
    g0_mu = PHI * g0_e
    A_e = A_tail(solver, g0_e)
    A_mu = A_tail(solver, g0_mu)

    print("  g0_e  = %.10f" % g0_e)
    print("  g0_mu = %.10f  (phi * g0_e)" % g0_mu)
    print("  A_e   = %.12f" % A_e)
    print("  A_mu  = %.12f" % A_mu)
    print("  r_21  = %.6f (target %.6f)" % ((A_mu / A_e)**4, r_21))
    print()

    # 2. Find Koide g0_tau - FINE scan to catch narrow bracket
    def koide_res(g0_3):
        A3 = A_tail(solver, g0_3)
        if A3 < 1e-15:
            return 1e10
        S2 = A_e**2 + A_mu**2 + A3**2
        S4 = A_e**4 + A_mu**4 + A3**4
        return S2**2 / S4 - 1.5

    # Scan with 200 pts in [g0_mu+0.01, 2.5] to catch narrow brackets
    g0_lo = g0_mu + 0.01
    g0_hi = 2.5
    g0_scan = np.linspace(g0_lo, g0_hi, 200)
    resids = []
    for g0 in g0_scan:
        res = koide_res(g0)
        resids.append(res)

    bracket = None
    for i in range(len(resids) - 1):
        if resids[i] != 1e10 and resids[i + 1] != 1e10:
            if resids[i] * resids[i + 1] < 0:
                bracket = (g0_scan[i], g0_scan[i + 1])
                break

    if bracket is None:
        # Try also scanning below g0_mu (shouldn't be physical but check)
        print("  No Koide bracket in [%.3f, %.3f]" % (g0_lo, g0_hi))
        # Diagnostic: print Q_K at key points
        for g0_3 in np.linspace(g0_mu + 0.05, min(g0_hi, 1.6), 10):
            A3 = A_tail(solver, g0_3)
            if A3 > 1e-15:
                S2 = A_e**2 + A_mu**2 + A3**2
                S4 = A_e**4 + A_mu**4 + A3**4
                Q = S2**2 / S4
                print("    g0=%.4f A3=%.8f Q_K=%.6f" % (g0_3, A3, Q))
        return None

    g0_tau = brentq(koide_res, bracket[0], bracket[1], xtol=1e-12)
    A_tau = A_tail(solver, g0_tau)

    S2 = A_e**2 + A_mu**2 + A_tau**2
    S4 = A_e**4 + A_mu**4 + A_tau**4
    Q_K = S2**2 / S4
    r31_K = (A_tau / A_e)**4
    c_K = g0_tau / g0_e

    print("  KOIDE SOLUTION FOUND:")
    print("  g0_tau = %.10f" % g0_tau)
    print("  A_tau  = %.12f" % A_tau)
    print("  Q_K    = %.12f  (target 1.500000)" % Q_K)
    print("  r_31   = %.2f  (PDG: %.2f)" % (r31_K, r_31))
    print("  c = g0_tau/g0_e = %.10f" % c_K)
    print("  Dev from PDG r_31: %+.6f%%" % ((r31_K / r_31 - 1) * 100))

    return {
        "g0_e": g0_e, "g0_mu": g0_mu, "g0_tau": g0_tau,
        "A_e": A_e, "A_mu": A_mu, "A_tau": A_tau,
        "Q_K": Q_K, "r31": r31_K, "c": c_K,
    }


# ====================================================================
print()
print("KOIDE FROM CANONICAL TGP ODE: FORM A vs FORM B")
print("(Fine scan to capture narrow bracket near critical g0)")
print()

res_A = find_koide("FORM A: alpha=2, source=g^2(1-g) [CANONICAL TGP]", solver_A)
print()
res_B = find_koide("FORM B: alpha=1, source=(1-g) [used in earlier scripts]", solver_B)

# ====================================================================
print()
print("=" * 65)
print("  SUMMARY")
print("=" * 65)
print()

header = "  %-22s %10s %12s %14s %10s %10s"
row    = "  %-22s %10.6f %12.8f %14.10f %10.1f %+10.4f%%"
print(header % ("Form", "g0_e", "c=g0t/g0e", "Q_K", "r_31", "r31 dev"))
for label, res in [("A (canonical)", res_A), ("B (scripts)", res_B)]:
    if res:
        print(row % (label, res["g0_e"], res["c"], res["Q_K"],
                      res["r31"], (res["r31"] / r_31 - 1) * 100))
    else:
        print("  %-22s %10s" % (label, "FAILED"))

print()
print("  PDG: r_21 = %.3f, r_31 = %.2f" % (r_21, r_31))
print("  Koide target: Q_K = 1.500000")

if res_A:
    print()
    print("  CRITICAL g0 ANALYSIS (Form A):")
    print("  g0_tau = %.6f" % res_A["g0_tau"])
    print("  g0_crit approx 1.600 (soliton collapse threshold)")
    margin = 1.600 - res_A["g0_tau"]
    print("  Margin: %.4f (%.2f%% below critical)" % (margin, margin / 1.600 * 100))

    print()
    print("  CONCLUSION: Canonical Form A (alpha=2, K=g^4) DOES produce")
    print("  Koide Q_K = 3/2, but g0_tau sits just below the critical")
    print("  threshold where solitons cease to exist.")
    print("  This suggests a deep connection between the Koide condition")
    print("  and the soliton existence boundary in TGP.")
