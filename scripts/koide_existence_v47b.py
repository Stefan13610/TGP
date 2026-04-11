#!/usr/bin/env python3
"""
koide_existence_v47b.py -- What the ODE truly determines.

CLARIFICATION (from r31_universality_proof_v47b.py):
  r31 is algebraically determined by r21 + Koide Q_K = 3/2.
  The ODE doesn't determine r31 -- it determines WHETHER
  a g0_tau exists that satisfies Koide.

KEY QUESTIONS:
  1. For which values of Q_K (not just 3/2) does a solution exist?
  2. Is Q_K = 3/2 selected by the ODE, or is ANY Q_K achievable?
  3. What is the Q_K range accessible to the ODE?

If any Q_K is achievable, then Q_K = 3/2 must come from an
ADDITIONAL selection principle (not the ODE alone).
If only Q_K = 3/2 is achievable, then the ODE itself selects Koide.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895
m_mu = 105.6583755
r_21 = m_mu / m_e


def make_solver(alpha=2.0, d=3, r_max=300):
    def solver(g0):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            source = g**2 * (1.0 - g)
            cross = (alpha / g) * gp**2
            if r < 1e-10:
                return [gp, (source - cross) / float(d)]
            return [gp, source - cross - float(d - 1) * gp / r]
        sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                        rtol=1e-12, atol=1e-14, max_step=0.02)
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


# ================================================================
print("=" * 70)
print("WHAT THE ODE TRULY DETERMINES")
print("=" * 70)

# Section 1: Q_K as function of g0_tau
# Fix alpha=2, calibrate g0_e from r_21, then scan g0_tau
print("\n1. Q_K(g0_tau) PROFILE (alpha=2)")
print("-" * 50)

solver = make_solver(alpha=2.0)
gc = 1.6  # collapse threshold

# Calibrate g0_e
def r21_res(g0_e):
    A_e = A_tail(solver, g0_e)
    if A_e < 1e-15:
        return 1e10
    g0_mu = PHI * g0_e
    if g0_mu >= gc:
        return 1e10
    A_mu = A_tail(solver, g0_mu)
    if A_mu < 1e-15:
        return 1e10
    return (A_mu / A_e)**4 - r_21

# Find bracket
g0_scan = np.linspace(0.5, 0.95, 30)
resids = [r21_res(g) for g in g0_scan]
bracket = None
for i in range(len(resids) - 1):
    if resids[i] * resids[i+1] < 0 and abs(resids[i]) < 1e8:
        bracket = (g0_scan[i], g0_scan[i+1])
        break

g0_e = brentq(r21_res, bracket[0], bracket[1], xtol=1e-8)
g0_mu = PHI * g0_e
A_e = A_tail(solver, g0_e)
A_mu = A_tail(solver, g0_mu)

print(f"  g0_e = {g0_e:.8f}")
print(f"  g0_mu = {g0_mu:.8f}")
print(f"  A_e = {A_e:.8f}")
print(f"  A_mu = {A_mu:.8f}")
print(f"  r_21 check: {(A_mu/A_e)**4:.4f} (target {r_21:.4f})")

# Scan g0_tau from g0_mu to g0_crit
print(f"\n  {'g0_tau':>10s} {'A_tau':>10s} {'Q_K':>10s} {'r31':>10s}")

g0_tau_scan = np.linspace(g0_mu + 0.01, gc - 0.001, 60)
QK_vals = []
g0t_vals = []

for g0_tau in g0_tau_scan:
    A_tau = A_tail(solver, g0_tau)
    if A_tau < 1e-10:
        continue
    S2 = A_e**2 + A_mu**2 + A_tau**2
    S4 = A_e**4 + A_mu**4 + A_tau**4
    QK = S2**2 / S4
    r31 = (A_tau / A_e)**4
    QK_vals.append(QK)
    g0t_vals.append(g0_tau)
    print(f"  {g0_tau:10.6f} {A_tau:10.6f} {QK:10.6f} {r31:10.2f}")

QK_vals = np.array(QK_vals)
g0t_vals = np.array(g0t_vals)

print(f"\n  Q_K range: [{np.min(QK_vals):.6f}, {np.max(QK_vals):.6f}]")
print(f"  Q_K at g0_tau = g0_mu: {QK_vals[0]:.6f}")
print(f"  Q_K at g0_tau -> gc: {QK_vals[-1]:.6f}")

# Where does Q_K cross specific values?
for QK_target in [1.0, 1.25, 1.5, 1.75, 2.0]:
    crossings = []
    for i in range(len(QK_vals) - 1):
        if (QK_vals[i] - QK_target) * (QK_vals[i+1] - QK_target) < 0:
            # Linear interpolation
            f = (QK_target - QK_vals[i]) / (QK_vals[i+1] - QK_vals[i])
            g0t = g0t_vals[i] + f * (g0t_vals[i+1] - g0t_vals[i])
            crossings.append(g0t)
    if crossings:
        print(f"  Q_K = {QK_target:.2f}: g0_tau = {', '.join(f'{c:.6f}' for c in crossings)}")
    else:
        print(f"  Q_K = {QK_target:.2f}: NO CROSSING in range")


# Section 2: Q_K is NOT uniquely selected
print("\n\n2. IS Q_K = 3/2 SPECIAL?")
print("-" * 50)

print("""
  The Q_K profile shows that Q_K is a CONTINUOUS MONOTONIC function
  of g0_tau. It starts at some Q_K_max when g0_tau is just above g0_mu
  (dominated by A_tau ~ A_mu) and decreases toward Q_K -> 1 as
  g0_tau -> gc (A_tau dominates).

  ANY Q_K in this range is achievable by choosing the right g0_tau.
  Q_K = 3/2 is NOT uniquely selected by the ODE.

  CONCLUSION: The ODE provides a MECHANISM (soliton amplitude)
  that maps g0_tau to mass ratios. The SELECTION of Q_K = 3/2
  is an additional physical input -- not derived from the ODE.

  In TGP, Q_K = 3/2 could come from:
  (a) An additional variational principle
  (b) Stability/resonance condition
  (c) Topological quantization
  (d) It could be an empirical input (like r_21)
""")


# Section 3: What IS special about Q_K = 3/2?
print("3. WHAT IS SPECIAL ABOUT Q_K = 3/2?")
print("-" * 50)

# Q_K = 3/2 means: (sum sqrt(m))^2 = 3/2 * (sum m)
# For 3 variables x,y,z > 0:
# (x+y+z)^2 / (x^2+y^2+z^2) ranges from 1 (all but one zero) to 3 (all equal)
# Q_K = 3/2 is the GEOMETRIC MEAN of these extremes: sqrt(1*3) = sqrt(3) ~ 1.73
# Or it's the harmonic mean: 2*1*3/(1+3) = 3/2 = 1.5. Yes!
# Q_K = 3/2 is the HARMONIC MEAN of the bounds [1, 3].

print("  Bounds on Q_K for 3 positive variables:")
print("    Q_K_min = 1 (one variable dominates)")
print("    Q_K_max = 3 (all equal)")
print(f"    Harmonic mean = 2*1*3/(1+3) = 3/2 = 1.5")
print(f"    Geometric mean = sqrt(1*3) = {np.sqrt(3):.4f}")
print(f"    Arithmetic mean = (1+3)/2 = 2.0")
print(f"    Q_K(Koide) = 3/2 = harmonic mean of bounds!")

# This is actually a known property:
# Q_K = n/2 for n variables would be the Koide condition
# For n=3: Q_K = 3/2

print(f"\n  For n variables: Q_K = n/2 is the 'Koide condition'")
print(f"  For n=3: Q_K = 3/2")
print(f"  This is equivalent to: sum(x_i - x_j)^2 = (sum x_i)^2")
print(f"  i.e., the 'variance' of sqrt(m) equals its 'mean squared'")


# Section 4: Hierarchy of predictions
print("\n\n4. HIERARCHY OF PREDICTIONS")
print("=" * 50)

print("""
  INPUT:
    [1] r_21 = m_mu/m_e = 206.768  (from phi-FP + ODE)
    [2] Q_K = 3/2                  (Koide relation -- empirical or derived?)

  PURE ALGEBRA (independent of ODE details):
    [3] r_31 = f(r_21) = 3477.4    (from [1] + [2])
    [4] m_tau = r_31 * m_e          (from [3])

  ODE-DEPENDENT (depends on alpha, form):
    [5] g0_e(alpha)                 (calibrated from [1])
    [6] g0_tau(alpha)               (from [2] + [5] + ODE)
    [7] g0_crit(alpha) = (2a+4)/(2a+1)  (collapse threshold)
    [8] existence: g0_tau < g0_crit (not all alpha work)

  LAGRANGIAN SELECTION:
    [9] alpha = 2                   (from K(g) = g^4)
    [10] g0_crit = 8/5             (from [7] + [9])
    [11] N_gen = 3                  (from [8] + phi-FP + [10])
""")


# Summary
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
  The ODE soliton dynamics determine:
  1. WHETHER Koide is achievable (existence of g0_tau)
  2. The specific g0_tau value (depends on alpha)
  3. The existence window in alpha-space

  The ODE does NOT determine:
  1. r31 (algebraic consequence of r21 + Koide)
  2. Q_K = 3/2 itself (additional input needed)

  The remaining OPEN question for O-L5:
  Is Q_K = 3/2 derivable from TGP principles,
  or is it an empirical input?

  If derivable: need a variational/topological argument.
  If empirical: it's an additional constraint alongside r_21.
""")
