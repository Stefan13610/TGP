#!/usr/bin/env python3
"""
koide_alpha_window_v47b.py -- For which alpha values does Koide Q_K=3/2 exist?

The collapse threshold g0_crit = (2a+4)/(2a+1) shrinks with alpha.
The muon g0_mu = phi*g0_e is fixed by r_21 calibration.
Tau requires g0_tau between g0_mu and g0_crit.

KEY QUESTION: Is alpha=2 (canonical TGP) special for Koide?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21 = m_mu / m_e


def make_solver(alpha, r_max=300):
    """ODE solver with given alpha, source=g^2(1-g)."""
    def solver(g0):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            source = g**2 * (1.0 - g)
            cross = (alpha / g) * gp**2
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


def find_koide_for_alpha(alpha):
    """Find Koide solution for given alpha. Returns dict or None."""
    solver = make_solver(alpha)
    g0_crit = (2 * alpha + 4) / (2 * alpha + 1)

    # Calibrate g0_e from r_21
    def r21_res(g0_1):
        A1 = A_tail(solver, g0_1)
        if g0_1 * PHI >= g0_crit:
            return 1e10
        A2 = A_tail(solver, PHI * g0_1)
        if A1 < 1e-15:
            return 1e10
        return (A2 / A1)**4 - r_21

    # Scan for bracket
    g0_scan = np.linspace(0.5, min(0.98, g0_crit / PHI - 0.01), 30)
    resids = [r21_res(g) for g in g0_scan]
    bracket = None
    for i in range(len(resids) - 1):
        if resids[i] * resids[i + 1] < 0 and abs(resids[i]) < 1e8:
            bracket = (g0_scan[i], g0_scan[i + 1])
            break

    if bracket is None:
        return None

    g0_e = brentq(r21_res, bracket[0], bracket[1], xtol=1e-8)
    g0_mu = PHI * g0_e
    A_e = A_tail(solver, g0_e)
    A_mu = A_tail(solver, g0_mu)

    if A_e < 1e-15 or A_mu < 1e-15:
        return None

    # Find Koide g0_tau in (g0_mu, g0_crit)
    def koide_res(g0_3):
        A3 = A_tail(solver, g0_3)
        if A3 < 1e-15:
            return 1e10
        S2 = A_e**2 + A_mu**2 + A3**2
        S4 = A_e**4 + A_mu**4 + A3**4
        return S2**2 / S4 - 1.5

    g0_lo = g0_mu + 0.005
    g0_hi = g0_crit - 0.001
    if g0_hi <= g0_lo:
        return {"g0_e": g0_e, "g0_mu": g0_mu, "g0_crit": g0_crit,
                "status": "NO_ROOM", "margin": g0_crit - g0_mu}

    # Scan (coarser for speed)
    g0_tau_scan = np.linspace(g0_lo, g0_hi, 60)
    k_resids = []
    for g0 in g0_tau_scan:
        k_resids.append(koide_res(g0))

    k_bracket = None
    for i in range(len(k_resids) - 1):
        if k_resids[i] != 1e10 and k_resids[i + 1] != 1e10:
            if k_resids[i] * k_resids[i + 1] < 0:
                k_bracket = (g0_tau_scan[i], g0_tau_scan[i + 1])
                break

    if k_bracket is None:
        # Report Q_K range
        valid_Q = [1.5 + r for r in k_resids if abs(r) < 100]
        Q_min = min(valid_Q) if valid_Q else 0
        Q_max = max(valid_Q) if valid_Q else 0
        return {"g0_e": g0_e, "g0_mu": g0_mu, "g0_crit": g0_crit,
                "status": "NO_KOIDE", "margin": g0_crit - g0_mu,
                "Q_range": (Q_min, Q_max)}

    g0_tau = brentq(koide_res, k_bracket[0], k_bracket[1], xtol=1e-10)
    A_tau = A_tail(solver, g0_tau)
    S2 = A_e**2 + A_mu**2 + A_tau**2
    S4 = A_e**4 + A_mu**4 + A_tau**4
    Q_K = S2**2 / S4
    r31 = (A_tau / A_e)**4

    return {
        "g0_e": g0_e, "g0_mu": g0_mu, "g0_tau": g0_tau,
        "g0_crit": g0_crit, "Q_K": Q_K, "r31": r31,
        "c": g0_tau / g0_e,
        "margin": g0_crit - g0_tau,
        "margin_pct": (g0_crit - g0_tau) / g0_crit * 100,
        "status": "KOIDE_FOUND",
    }


# ================================================================
print("=" * 70)
print("KOIDE EXISTENCE WINDOW IN alpha")
print("=" * 70)
print()
print("  %-6s %-8s %-8s %-8s %-8s %-10s %-8s %-10s" %
      ("alpha", "g0_crit", "g0_e", "g0_mu", "g0_tau", "margin", "Q_K", "status"))

alpha_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

results = []
for alpha in alpha_vals:
    res = find_koide_for_alpha(alpha)
    results.append((alpha, res))

    if res is None:
        print("  %-6.1f  CALIBRATION FAILED" % alpha)
    elif res["status"] == "KOIDE_FOUND":
        print("  %-6.1f %-8.4f %-8.4f %-8.4f %-8.4f %-10.4f %-8.6f FOUND" %
              (alpha, res["g0_crit"], res["g0_e"], res["g0_mu"],
               res["g0_tau"], res["margin"], res["Q_K"]))
    elif res["status"] == "NO_ROOM":
        print("  %-6.1f %-8.4f %-8.4f %-8.4f %-8s %-10.4f %-8s NO_ROOM" %
              (alpha, res["g0_crit"], res["g0_e"], res["g0_mu"],
               "---", res["margin"], "---"))
    elif res["status"] == "NO_KOIDE":
        Q_lo, Q_hi = res.get("Q_range", (0, 0))
        print("  %-6.1f %-8.4f %-8.4f %-8.4f %-8s %-10.4f %-8s Q:[%.2f,%.2f]" %
              (alpha, res["g0_crit"], res["g0_e"], res["g0_mu"],
               "---", res["margin"], "---", Q_lo, Q_hi))

# Summary
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

koide_alphas = [a for a, r in results
                if r and r["status"] == "KOIDE_FOUND"]
if koide_alphas:
    print("\n  Koide Q_K = 3/2 exists for alpha in [%.1f, %.1f]" %
          (min(koide_alphas), max(koide_alphas)))
    print("  Canonical TGP: alpha = 2.0")
    print()

    # Show margins
    print("  Margin (g0_crit - g0_tau) vs alpha:")
    for alpha, res in results:
        if res and res["status"] == "KOIDE_FOUND":
            print("    alpha=%.1f: margin=%.4f (%.2f%%)" %
                  (alpha, res["margin"], res["margin_pct"]))

    print()
    print("  The margin SHRINKS with increasing alpha.")
    print("  At some alpha_max, the tau soliton can no longer fit")
    print("  below the collapse threshold.")
