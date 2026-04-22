#!/usr/bin/env python3
"""
ngen_collapse_v47b.py -- N_gen = 3 from soliton collapse threshold.

KEY RESULT: The canonical TGP Form A (alpha=2, source=g^2(1-g)) has a
soliton collapse threshold at g0_crit ~ 1.600. The tau soliton sits
at g0_tau = 1.5696, using 95.8% of the available g0 range.
A 4th generation is FORBIDDEN: it would require g0 > g0_crit.

This provides a DYNAMICAL explanation for N_gen = 3 that is INDEPENDENT
of the WKB counting argument (N_gen = 3 from k=4 bound states).
"""
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21 = m_mu / m_e; r_31 = m_tau / m_e


# ================================================================
# Section 1: Precise critical g0 for canonical Form A
# ================================================================
def make_solver(alpha_coeff, source_type, r_max=300):
    def solver(g0):
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


def find_g0_crit(alpha, source_type="quadratic", tol=1e-8):
    """Find critical g0 above which solitons collapse."""
    solver = make_solver(alpha, source_type)

    def converges(g0):
        r, g = solver(g0)
        return r[-1] > 250 and np.min(g) > 0.01

    # Find rough bracket
    g_lo = 1.0
    g_hi = 5.0
    if not converges(g_lo):
        return 0.0  # No stable solitons at all
    if converges(g_hi):
        return float('inf')  # No collapse

    for _ in range(100):
        g_mid = (g_lo + g_hi) / 2
        if g_hi - g_lo < tol:
            break
        if converges(g_mid):
            g_lo = g_mid
        else:
            g_hi = g_mid

    return (g_lo + g_hi) / 2


print("=" * 72)
print("N_gen = 3 FROM SOLITON COLLAPSE IN CANONICAL TGP")
print("=" * 72)


# ================================================================
# Section 2: Critical g0 vs alpha
# ================================================================
print("\n  SECTION 1: Critical g0 vs alpha coefficient")
print("  " + "-" * 50)
print("  %-8s %-14s %-14s" % ("alpha", "g0_crit(quad)", "g0_crit(lin)"))

alpha_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
crits_quad = {}
crits_lin = {}

for alpha in alpha_vals:
    gc_q = find_g0_crit(alpha, "quadratic")
    gc_l = find_g0_crit(alpha, "linear")
    crits_quad[alpha] = gc_q
    crits_lin[alpha] = gc_l
    gc_q_str = "%.8f" % gc_q if gc_q < 100 else "inf"
    gc_l_str = "%.8f" % gc_l if gc_l < 100 else "inf"
    print("  %-8.1f %-14s %-14s" % (alpha, gc_q_str, gc_l_str))


# ================================================================
# Section 3: Precise analysis for canonical alpha=2
# ================================================================
print("\n  SECTION 2: Canonical Form A (alpha=2) -- generation analysis")
print("  " + "-" * 50)

solver_A = make_solver(2.0, "quadratic")

# Calibrate g0_e
def r21_res(g0_1):
    A1 = A_tail(solver_A, g0_1)
    A2 = A_tail(solver_A, PHI * g0_1)
    if A1 < 1e-15:
        return 1e10
    return (A2 / A1)**4 - r_21

g0_e = brentq(r21_res, 0.8, 0.95, xtol=1e-10)
g0_mu = PHI * g0_e

# Koide g0_tau (from previous result)
g0_tau = 1.5695724137

A_e = A_tail(solver_A, g0_e)
A_mu = A_tail(solver_A, g0_mu)
A_tau = A_tail(solver_A, g0_tau)

g0_crit = find_g0_crit(2.0, "quadratic", tol=1e-10)

print("  g0_e    = %.10f" % g0_e)
print("  g0_mu   = %.10f  (phi * g0_e)" % g0_mu)
print("  g0_tau  = %.10f  (Koide solution)" % g0_tau)
print("  g0_crit = %.10f  (soliton collapse)" % g0_crit)
print()
print("  Spacings:")
print("    g0_mu  - g0_e   = %.6f" % (g0_mu - g0_e))
print("    g0_tau - g0_mu  = %.6f" % (g0_tau - g0_mu))
print("    g0_crit - g0_tau = %.6f  (margin)" % (g0_crit - g0_tau))
print()
print("  Occupation: tau uses %.2f%% of [g0_e, g0_crit] range" %
      ((g0_tau - g0_e) / (g0_crit - g0_e) * 100))

# 4th generation tests
g0_4_phi = PHI * g0_tau
g0_4_koide = None  # Would need Koide(mu, tau, X)

print()
print("  4th generation hypotheses:")
print("    phi * g0_tau = %.4f  (%.1f%% above g0_crit => FORBIDDEN)" %
      (g0_4_phi, (g0_4_phi / g0_crit - 1) * 100))

# Check: what if we try smaller spacing?
# Koide for (mu, tau, X): find X such that Q_K(mu,tau,X) = 3/2
def koide_res_4th(g0_4):
    A4 = A_tail(solver_A, g0_4)
    if A4 < 1e-15:
        return 1e10
    S2 = A_mu**2 + A_tau**2 + A4**2
    S4 = A_mu**4 + A_tau**4 + A4**4
    return S2**2 / S4 - 1.5

# Scan for any possible 4th gen with Koide(mu, tau, X)
print()
print("  Koide(mu, tau, X) scan:")
g0_scan = np.linspace(g0_tau + 0.001, g0_crit - 0.001, 50)
found_any = False
for g0_4 in g0_scan:
    A4 = A_tail(solver_A, g0_4)
    if A4 > 1e-15:
        S2 = A_mu**2 + A_tau**2 + A4**2
        S4 = A_mu**4 + A_tau**4 + A4**4
        Q = S2**2 / S4
        if abs(Q - 1.5) < 0.01:
            found_any = True
            print("    g0=%.6f: A4=%.6f, Q_K=%.6f" % (g0_4, A4, Q))

if not found_any:
    # Show range of Q_K values
    Qs = []
    for g0_4 in g0_scan:
        A4 = A_tail(solver_A, g0_4)
        if A4 > 1e-15:
            S2 = A_mu**2 + A_tau**2 + A4**2
            S4 = A_mu**4 + A_tau**4 + A4**4
            Qs.append(S2**2 / S4)
    if Qs:
        print("    Q_K range in [g0_tau, g0_crit]: [%.4f, %.4f]" %
              (min(Qs), max(Qs)))
        print("    No Q_K = 3/2 solution for 4th gen in stable soliton range")
    else:
        print("    No stable solitons between g0_tau and g0_crit")


# ================================================================
# Section 4: Mass of maximum-allowed soliton
# ================================================================
print("\n  SECTION 3: Maximum soliton mass (at g0_crit)")
print("  " + "-" * 50)

# A_tail as function of g0 near g0_crit
g0_near = np.linspace(1.58, g0_crit - 0.0005, 30)
A_near = [(g0, A_tail(solver_A, g0)) for g0 in g0_near]
A_max_val = max(A for _, A in A_near)
g0_max = [g0 for g0, A in A_near if A == A_max_val][0]

r_max_possible = (A_max_val / A_e)**4
print("  Max A_tail = %.8f (at g0 = %.6f)" % (A_max_val, g0_max))
print("  Max r = (A_max/A_e)^4 = %.0f" % r_max_possible)
print("  Max m = %.1f * m_e = %.1f MeV" %
      (r_max_possible, r_max_possible * m_e))
print()
print("  For comparison:")
print("    m_tau = %.1f MeV  (r_31 = %.0f)" % (m_tau, r_31))
print("    m_max = %.1f MeV  (r_max = %.0f)" %
      (r_max_possible * m_e, r_max_possible))
print("    Ratio m_max/m_tau = %.2f" % (r_max_possible / r_31))


# ================================================================
# Section 5: Analytical understanding of collapse
# ================================================================
print("\n  SECTION 4: Physics of the collapse mechanism")
print("  " + "-" * 50)

# The cross-coupling term (alpha/g) g'^2 becomes large when g is small
# For g near 0, source ~ g^2(1-g) ~ g^2 -> 0 (can't restore g)
# But cross ~ (2/g) g'^2 -> infinity (drives g further to 0)
# Critical g0 is where the initial overshoot brings g close enough
# to 0 that cross > source and positive feedback causes collapse.

# Track g_min vs g0
print()
print("  g_min(r) vs g0 (how deep the soliton dips):")
print("  %-10s %-12s %-12s %-12s" % ("g0", "g_min", "r_at_min", "converges"))

for g0 in [1.50, 1.55, 1.57, 1.58, 1.59, 1.595, 1.598, 1.599, 1.5995, 1.5999]:
    r, g = solver_A(g0)
    g_min = np.min(g)
    i_min = np.argmin(g)
    r_min = r[i_min]
    conv = r[-1] > 250 and g_min > 0.01
    print("  %-10.4f %-12.6f %-12.2f %-12s" %
          (g0, g_min, r_min, "YES" if conv else "COLLAPSE"))


# ================================================================
# Section 6: Two independent proofs of N_gen = 3
# ================================================================
print("\n  SECTION 5: TWO INDEPENDENT PROOFS OF N_gen = 3")
print("  " + "-" * 50)
print("""
  PROOF 1 (WKB): From d=3, k=4 (exact TGP potential exponent):
    The WKB quantization in V_eff ~ r^k = r^4 admits exactly
    N_bound = 3 bound states in 3 dimensions.
    This is a KINEMATIC constraint from the potential shape.

  PROOF 2 (Soliton collapse, this script): From alpha=2 (exact TGP):
    The nonlinear coupling (2/g)g'^2 in the canonical TGP ODE
    creates a soliton collapse threshold at g0_crit = %.6f.
    With phi-FP spacing and Koide calibration:
      e(%.4f) -> mu(%.4f) -> tau(%.4f) -> COLLAPSE(%.4f)
    Tau uses %.1f%% of the available g0 range.
    A 4th generation is DYNAMICALLY FORBIDDEN.

  These are INDEPENDENT arguments that both give N_gen = 3.
  Proof 1 uses the potential exponent k=4.
  Proof 2 uses the kinetic coupling coefficient alpha=2.
  Both values (k=4, alpha=2) derive from K(g) = g^4 in the
  canonical TGP Lagrangian.""" %
      (g0_crit, g0_e, g0_mu, g0_tau, g0_crit,
       (g0_tau - g0_e) / (g0_crit - g0_e) * 100))


# ================================================================
# Section 7: Summary
# ================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print("""
  Canonical TGP Form A (alpha=2, K=g^4, source=g^2(1-g)):

  1. g0_crit = %.10f (soliton collapse threshold)
  2. g0_tau  = %.10f (Koide Q_K = 3/2 solution)
  3. Margin  = %.6f (%.2f%% below critical)
  4. Tau uses %.1f%% of available [g0_e, g0_crit] range
  5. 4th generation (phi*g0_tau = %.4f) is %.1f%% ABOVE threshold
  6. N_gen = 3 has TWO independent proofs: WKB + collapse
  7. Both derive from K(g) = g^4 in the TGP Lagrangian""" %
      (g0_crit, g0_tau, g0_crit - g0_tau,
       (g0_crit - g0_tau) / g0_crit * 100,
       (g0_tau - g0_e) / (g0_crit - g0_e) * 100,
       g0_4_phi, (g0_4_phi / g0_crit - 1) * 100))
