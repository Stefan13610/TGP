#!/usr/bin/env python3
"""
koide_observer_critical_v47b.py -- PATH 18b: CRITICAL RE-ANALYSIS
                                   of the "observer effect" claim.

PROBLEM IDENTIFIED:
  In koide_observer_v47b.py, g0_tau was FOUND BY ENFORCING Q_K(A^2)=3/2.
  Therefore Q_K(A^4)=3/2 is trivially the same condition restated.
  The "observation chain" analysis may be tautological.

CRITICAL TESTS:
  Test A: TAUTOLOGY CHECK
    If Q_K(x_i) = 3/2 for some x_i, what is Q_K(x_i^2)?
    Is Q_K(x^2)=3/2 a GENERIC consequence of Q_K(x)=3/2?
    Answer: NO. Q_K(x^2) != Q_K(x) in general.
    But: Q_K as defined in the code is (sum sqrt(x))^2 / (3*sum(x)).
    If x_i = A_i^2, then Q_K(x_i) = (sum A_i)^2/(3*sum A_i^2).
    The script sets THIS to 3/2.
    Then Q_K(A^4) = Q_K(x_i^2 where x_i=A_i^2)
                  = (sum sqrt(A_i^4))^2 / (3*sum A_i^4)
                  = (sum A_i^2)^2 / (3*sum A_i^4)
    But original Koide on sqrt(m) is: (sum sqrt(m))^2/(3*sum m)
    With m=A^4: (sum A^2)^2/(3*sum A^4) -- THIS IS THE STANDARD KOIDE.
    And the script enforced (sum A)^2/(3*sum A^2) = 3/2 -- DIFFERENT!

    Wait -- let me re-check what was enforced vs what was found.

  Test B: NON-TRIVIAL QUESTION
    If the ODE map g0->A is given, and we pick g0 triplets with phi-FP
    spacing, does Q_K(mass) generically approach 3/2?
    Or is it specific to one particular g0_e?
    THIS would be a genuine "observer effect".

  Test C: WHAT DOES THE NONLINEAR MAP ACTUALLY DO?
    For RANDOM triplets (not phi-FP spaced), does the ODE map
    g0->A->A^2->A^4 push Q_K toward any particular value?
    If squaring is "special," random A should cluster near Q_K=3/2
    after squaring.

  Test D: COLLAPSE PROXIMITY EFFECT
    G(g0) = A^2/(g0-1)^2 varies by 3.2x.
    Is this variation SUFFICIENT and NECESSARY for Koide?
    What if G were slightly different?
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0


def solver_A(g0, r_max=400):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0]


def A_tail(g0, r_min=50, r_max=300):
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


def QK(values):
    """Standard Koide: Q_K = (sum sqrt(x))^2 / (3 * sum(x))"""
    v = np.array(values, dtype=float)
    v = np.abs(v)
    v = v[v > 0]
    if len(v) < 3:
        return 0
    sq = np.sqrt(v)
    return np.sum(sq)**2 / (3.0 * np.sum(v))


# ================================================================
print("=" * 70)
print("PATH 18b: CRITICAL RE-ANALYSIS OF OBSERVER EFFECT")
print("=" * 70)

# ================================================================
print("\n" + "=" * 70)
print("TEST A: TAUTOLOGY CHECK -- WHAT WAS ACTUALLY ENFORCED?")
print("=" * 70)

# Reproduce the original calculation
g0_e = 0.86770494
g0_mu = PHI * g0_e

A_e = A_tail(g0_e)
A_mu = A_tail(g0_mu)

print(f"\n  Given: g0_e = {g0_e:.8f}, g0_mu = phi*g0_e = {g0_mu:.8f}")
print(f"  A_e = {A_e:.8f}, A_mu = {A_mu:.8f}")

# The original script enforced: (sum A_i)^2 / (sum A_i^2) = 3/2
# Wait -- let me re-read the original code carefully.
# Line 121-122: xe, xm, xt = A_e**2, A_mu**2, At**2
#               return (xe+xm+xt)**2/(xe**2+xm**2+xt**2) - 1.5
# So x_i = A_i^2, and the residual is (sum x)^2/(sum x^2) - 1.5
# Note: this is NOT the standard Q_K!
# Standard Q_K = (sum sqrt(x))^2/(3*sum(x))
# The residual is: (sum x)^2/(sum x^2) which equals 3*Q_K(x^2)
# So it enforced: 3*Q_K(x_i^2) = 3/2 where x_i = A_i^2
# i.e., Q_K(A_i^4) = 1/2 ???
# No wait: Q_K(m) = (sum sqrt(m))^2/(3*sum m)
# With m = A^4: Q_K(m) = (sum A^2)^2/(3*sum A^4)
# The residual is (sum A^2)^2/(sum A^4) = 1.5
# So: (sum A^2)^2/(sum A^4) = 3/2
# And Q_K(m) = (sum A^2)^2/(3*sum A^4) = 1.5/3 = 1/2 ???
#
# NO -- standard Koide Q is 2/3, not 3/2!
# Different normalization conventions. Let me check.
#
# Convention 1 (original Koide): Q = (sum sqrt(m))^2 / (sum m) = 2
#   Koide's original: this equals 2 (not 2/3 or 3/2)
#   Actually: sum(sqrt(m))^2 / sum(m) = 2 means Q/3 = 2/3
#
# Convention used in code QK_func:
#   QK = (sum sqrt(v))^2 / sum(v)  -- NO factor of 3!
#   Wait, line 101-102 of original: return np.sum(sq)**2 / np.sum(v)
#   NO factor of 3!
#
# But the residual (line 122): (xe+xm+xt)^2/(xe^2+xm^2+xt^2) - 1.5
#   This is NOT QK_func applied to anything standard.
#   This is: (sum x)^2 / (sum x^2) where x = A^2
#   = (sum A^2)^2 / (sum A^4) = 1.5
#   Note: (sum A^2)^2 / (sum A^4) = 3 * [(sum A^2)^2 / (3*sum A^4)]
#   The term in brackets IS the standard Q_K of the masses (with /3).
#   But wait, the QK_func in the code doesn't have /3.
#   Let me just be explicit.

# What was enforced:
print("\n  What the original script ENFORCED (line 122):")
print("  (A_e^2 + A_mu^2 + A_tau^2)^2 / (A_e^4 + A_mu^4 + A_tau^4) = 1.5")
print("  This is equivalent to: sum(sqrt(m))^2 / sum(m) = 1.5")
print("  i.e., the STANDARD Koide relation.")

# Find g0_tau using the original method
def find_g0_tau_original():
    def resid(g0t):
        At = A_tail(g0t)
        xe, xm, xt = A_e**2, A_mu**2, At**2
        return (xe+xm+xt)**2/(xe**2+xm**2+xt**2) - 1.5
    return brentq(resid, 1.50, 1.59)

g0_tau = find_g0_tau_original()
A_tau = A_tail(g0_tau)

# Now check what QK_func returns for different levels
# QK_func computes: (sum sqrt(v))^2 / sum(v)  [no factor of 3]
A_vals = [A_e, A_mu, A_tau]
A2_vals = [a**2 for a in A_vals]
A4_vals = [a**4 for a in A_vals]

# The code's QK_func on A^4:
# (sum sqrt(A^4))^2 / sum(A^4) = (sum A^2)^2 / sum(A^4) = 1.5 BY CONSTRUCTION
qk_A4 = QK(A4_vals)
print(f"\n  QK_func(A^4) = {qk_A4:.6f}  (= 3/2 by construction? {abs(qk_A4-0.5)<0.001})")
print(f"  Note: QK_func has NO /3 factor.")
print(f"  QK_func(A^4) = (sum A^2)^2 / (sum A^4)")

# Now the residual enforced (sum A^2)^2/(sum A^4) = 1.5
# But QK_func(A^4) = (sum sqrt(A^4))^2/sum(A^4) = (sum A^2)^2/sum(A^4)
# So QK_func(A^4) should = 1.5 by construction!
ratio = sum(A2_vals)**2 / sum(A4_vals)
print(f"  Direct: (sum A^2)^2/(sum A^4) = {ratio:.6f}")
print(f"\n  VERDICT: QK_func(A^4) = 1.5 IS the enforced condition. TAUTOLOGICAL.")

# What about QK_func on other levels?
qk_g0 = QK([g0_e, g0_mu, g0_tau])
qk_dev = QK([abs(g0_e-1), abs(g0_mu-1), abs(g0_tau-1)])
qk_A = QK(A_vals)
qk_A2 = QK(A2_vals)
print(f"\n  QK_func values at different levels:")
print(f"  g0:    {qk_g0:.6f}")
print(f"  |g0-1|:{qk_dev:.6f}")
print(f"  A:     {qk_A:.6f}")
print(f"  A^2:   {qk_A2:.6f}")
print(f"  A^4:   {qk_A4:.6f}  <-- enforced = 1.5")
print(f"\n  The output table in PATH 18 RELABELED these as:")
print(f"  A^4 line said 'Q_K = 1.500000 actual masses <-- KOIDE!'")
print(f"  But this was the DEFINITION, not a finding!")

# ================================================================
print("\n" + "=" * 70)
print("TEST A2: IS THERE ANYTHING NON-TRIVIAL IN THE Q_K CHAIN?")
print("=" * 70)

# The chain g0 -> A -> A^2 -> A^4 with Q_K at each level IS non-trivial
# in the sense that it shows how Q_K transforms under the nonlinear map.
# But the Q_K values at g0, A, A^2 levels are CONSEQUENCES of the
# enforced Q_K(A^4)=3/2 combined with the specific nonlinear map.
# They don't provide new information -- they're just the same constraint
# viewed in different coordinates.

# THE ONE POTENTIALLY NON-TRIVIAL THING:
# If we DON'T enforce Koide, but instead use phi-FP for ALL THREE:
#   g0_e, g0_mu = phi*g0_e, g0_tau = phi*g0_mu = phi^2*g0_e
# Does the ODE map produce something close to Q_K = 3/2?

print("\n  Non-enforced test: g0 = g0_e * phi^{0,1,2}")
g0_1 = g0_e
g0_2 = PHI * g0_1
g0_3 = PHI * g0_2
A_1 = A_tail(g0_1)
A_2 = A_tail(g0_2)
A_3 = A_tail(g0_3)
print(f"  g0: {g0_1:.6f}, {g0_2:.6f}, {g0_3:.6f}")
print(f"  A:  {A_1:.8f}, {A_2:.8f}, {A_3:.8f}")
m_1, m_2, m_3 = A_1**4, A_2**4, A_3**4
print(f"  m:  {m_1:.8e}, {m_2:.8e}, {m_3:.8e}")
qk_phi_fp = QK([m_1, m_2, m_3])
print(f"  QK(masses) for pure phi-FP spacing = {qk_phi_fp:.6f}")
print(f"  Deviation from 3/2: {qk_phi_fp - 1.5:+.6f} ({(qk_phi_fp-1.5)/1.5*100:+.2f}%)")
print(f"\n  If this were close to 3/2, the ODE map + phi-FP would")
print(f"  'naturally produce' Koide. If not, Koide requires g0_tau != phi^2*g0_e.")

# Also check: what is the actual g0_tau vs phi^2*g0_e?
g0_tau_phiFP = PHI**2 * g0_e
print(f"\n  g0_tau (Koide-enforced) = {g0_tau:.8f}")
print(f"  g0_tau (phi^2*g0_e)     = {g0_tau_phiFP:.8f}")
print(f"  Difference: {g0_tau - g0_tau_phiFP:+.8f} ({(g0_tau/g0_tau_phiFP - 1)*100:+.4f}%)")


# ================================================================
print("\n" + "=" * 70)
print("TEST B: DOES ODE MAP PUSH Q_K TOWARD 3/2?")
print("=" * 70)

# Scan over g0_e values with phi-FP spacing
# For each g0_e, compute g0_mu = phi*g0_e, find A values, compute Q_K(mass)
print("\n  Scanning g0_e with phi-FP spacing (g0_mu = phi*g0_e, g0_tau = phi^2*g0_e):")
print(f"  {'g0_e':>10s} {'g0_mu':>10s} {'g0_tau':>10s} {'QK(mass)':>10s} {'QK-3/2':>10s}")

for g0e_test in np.arange(0.50, 0.96, 0.05):
    g0m = PHI * g0e_test
    g0t = PHI * g0m
    if g0t >= G0_CRIT:
        continue
    try:
        Ae = A_tail(g0e_test)
        Am = A_tail(g0m)
        At = A_tail(g0t)
        if Ae > 0 and Am > 0 and At > 0:
            qk = QK([Ae**4, Am**4, At**4])
            print(f"  {g0e_test:10.4f} {g0m:10.4f} {g0t:10.4f} {qk:10.6f} {qk-1.5:+10.6f}")
    except:
        pass

# ================================================================
print("\n" + "=" * 70)
print("TEST C: DOES SQUARING GENERICALLY PUSH Q_K TOWARD 3/2?")
print("=" * 70)

# Take random triplets of A values, compute Q_K(A^n) for n=1,2,4
# If squaring is "special," Q_K(A^2) or Q_K(A^4) should cluster near 3/2
np.random.seed(42)
N = 10000
qk_A_rand = np.zeros(N)
qk_A2_rand = np.zeros(N)
qk_A4_rand = np.zeros(N)

for i in range(N):
    # Random positive values with large hierarchy (like solitons)
    A_rand = np.sort(np.random.exponential(1.0, 3))
    qk_A_rand[i] = QK(A_rand)
    qk_A2_rand[i] = QK(A_rand**2)
    qk_A4_rand[i] = QK(A_rand**4)

print(f"\n  {N} random triplets from Exp(1):")
print(f"  {'Level':>10s} {'mean(QK)':>10s} {'std(QK)':>10s} {'%near3/2':>10s}")
for name, arr in [("A", qk_A_rand), ("A^2", qk_A2_rand), ("A^4", qk_A4_rand)]:
    near = np.mean(np.abs(arr - 1.5) < 0.05) * 100
    print(f"  {name:>10s} {np.mean(arr):10.4f} {np.std(arr):10.4f} {near:10.1f}%")

# Same but with power-law distributed A (large hierarchy)
qk_A_pl = np.zeros(N)
qk_A2_pl = np.zeros(N)
qk_A4_pl = np.zeros(N)

for i in range(N):
    # Power-law: strong hierarchy like leptons
    logA = np.sort(np.random.uniform(-3, 1, 3))
    A_rand = 10.0**logA
    qk_A_pl[i] = QK(A_rand)
    qk_A2_pl[i] = QK(A_rand**2)
    qk_A4_pl[i] = QK(A_rand**4)

print(f"\n  {N} random triplets from 10^Uniform(-3,1) (strong hierarchy):")
print(f"  {'Level':>10s} {'mean(QK)':>10s} {'std(QK)':>10s} {'%near3/2':>10s}")
for name, arr in [("A", qk_A_pl), ("A^2", qk_A2_pl), ("A^4", qk_A4_pl)]:
    near = np.mean(np.abs(arr - 1.5) < 0.05) * 100
    print(f"  {name:>10s} {np.mean(arr):10.4f} {np.std(arr):10.4f} {near:10.1f}%")

print("\n  If squaring pushes Q_K toward 3/2, we'd see Q_K(A^2) or Q_K(A^4)")
print("  closer to 3/2 than Q_K(A). If not, squaring has no special role.")


# ================================================================
print("\n" + "=" * 70)
print("TEST D: WHAT THE NONLINEAR MAP g0->A ACTUALLY DOES TO Q_K")
print("=" * 70)

# The real question: for the ODE map, does Q_K change in a specific way?
# Compare Q_K at g0 level vs A level vs mass level
# for the phi-FP triplets

print("\n  For phi-FP spacing with different g0_e:")
print(f"  {'g0_e':>8s} {'QK(g0)':>10s} {'QK(A)':>10s} {'QK(A^2)':>10s} {'QK(A^4)':>10s} {'QK(mass)':>10s}")

for g0e_test in np.arange(0.50, 0.95, 0.05):
    g0m = PHI * g0e_test
    g0t = PHI * g0m
    if g0t >= G0_CRIT:
        continue
    try:
        Ae = A_tail(g0e_test)
        Am = A_tail(g0m)
        At = A_tail(g0t)
        if Ae > 0 and Am > 0 and At > 0:
            qk_g = QK([g0e_test, g0m, g0t])
            qk_a = QK([Ae, Am, At])
            qk_a2 = QK([Ae**2, Am**2, At**2])
            qk_a4 = QK([Ae**4, Am**4, At**4])
            qk_m = QK([Ae**4, Am**4, At**4])
            print(f"  {g0e_test:8.4f} {qk_g:10.6f} {qk_a:10.6f} {qk_a2:10.6f} {qk_a4:10.6f} {qk_m:10.6f}")
    except:
        pass


# ================================================================
print("\n" + "=" * 70)
print("TEST E: THE REAL NON-TRIVIAL QUESTION")
print("=" * 70)
print("""
  The GENUINE observer effect question would be:

  Given the TGP ODE with its specific nonlinearity,
  does the map g0 -> A_tail have a special property such that
  phi-FP spaced g0 values GENERICALLY produce Q_K ~ 3/2 for masses?

  From Test B, we know Q_K(mass) for pure phi-FP spacing.
  If Q_K(mass) ~ 3/2 naturally, then Koide IS emergent.
  If Q_K(mass) != 3/2, then g0_tau must be selected by the mechanism.

  The actual g0_tau deviates from phi^2*g0_e. The selection
  mechanism is structural: d=3 -> 2 tail components -> chi2(2)
  -> CV=1 -> Q_K=3/2 (see rem:T-cv1-origin, thm:T-QK-CV).
""")

# What is the Q_K for the exact phi-FP triplet (no Koide enforcement)?
g0_1_exact = g0_e
g0_2_exact = PHI * g0_e
g0_3_exact = PHI**2 * g0_e

A1 = A_tail(g0_1_exact)
A2 = A_tail(g0_2_exact)
A3 = A_tail(g0_3_exact)

qk_natural = QK([A1**4, A2**4, A3**4])
print(f"  NATURAL Q_K (pure phi-FP, no Koide enforcement):")
print(f"  g0 = {g0_1_exact:.6f}, {g0_2_exact:.6f}, {g0_3_exact:.6f}")
print(f"  A  = {A1:.8f}, {A2:.8f}, {A3:.8f}")
print(f"  Q_K(mass) = {qk_natural:.6f}")
print(f"  Deviation from 3/2: {(qk_natural-1.5)/1.5*100:+.4f}%")

if abs(qk_natural - 1.5) < 0.05:
    print(f"\n  RESULT: Q_K ~ 3/2 NATURALLY from phi-FP + ODE!")
    print(f"  This WOULD support the observer effect interpretation.")
else:
    print(f"\n  RESULT: Q_K != 3/2 for pure phi-FP.")
    print(f"  Koide requires g0_tau != phi^2*g0_e.")
    print(f"  The 'observer effect' does NOT produce Koide generically.")


# ================================================================
print("\n" + "=" * 70)
print("TEST F: NONLINEAR ENHANCEMENT -- PREDICTIVE OR DESCRIPTIVE?")
print("=" * 70)

# G(g0) = A^2/(g0-1)^2
# The claim was G varies by 3.2x. But what does this mean?
# Let's compute G for many g0 values to see if it has universal features

print("\n  G(g0) = A_tail(g0)^2 / (g0-1)^2 across the ODE:")
print(f"  {'g0':>8s} {'A':>12s} {'(g0-1)^2':>12s} {'G':>12s} {'dG/dg0':>12s}")

g0_range = np.arange(0.50, 1.58, 0.05)
G_vals = []
g0_list = []
for g0 in g0_range:
    if abs(g0 - 1.0) < 0.02:
        continue
    try:
        A = A_tail(g0)
        if A > 0:
            G = A**2 / (g0-1)**2
            G_vals.append(G)
            g0_list.append(g0)
            print(f"  {g0:8.4f} {A:12.8f} {(g0-1)**2:12.8f} {G:12.6f}")
    except:
        pass

# What G values would produce Q_K = 3/2?
# If A_i^2 = G_i * (g0_i - 1)^2, then sqrt(m_i) = G_i * (g0_i-1)^2
# Q_K depends on ratios of G values AND g0 deviations
print(f"\n  G at lepton soliton positions:")
G_e = A_e**2 / (g0_e - 1)**2
G_mu = A_mu**2 / (g0_mu - 1)**2
G_tau = A_tau**2 / (g0_tau - 1)**2
print(f"  G_e   = {G_e:.6f}  (g0_e   = {g0_e:.6f})")
print(f"  G_mu  = {G_mu:.6f}  (g0_mu  = {g0_mu:.6f})")
print(f"  G_tau = {G_tau:.6f}  (g0_tau = {g0_tau:.6f})")
print(f"  G_tau/G_e = {G_tau/G_e:.4f}")

# If G were constant (linear regime), what Q_K?
G_const = G_e  # use electron value
dev_e = abs(g0_e - 1)
dev_mu = abs(g0_mu - 1)
dev_tau = abs(g0_tau - 1)
m_linear = [G_const**2 * dev_e**4, G_const**2 * dev_mu**4, G_const**2 * dev_tau**4]
qk_linear = QK(m_linear)
print(f"\n  If G = G_e = const: Q_K(mass) = {qk_linear:.6f}")
print(f"  (This is just Q_K of |g0-1|^4 values)")

# Use actual G values
m_actual = [G_e**2 * dev_e**4, G_mu**2 * dev_mu**4, G_tau**2 * dev_tau**4]
qk_actual = QK(m_actual)
print(f"  With actual G values: Q_K(mass) = {qk_actual:.6f}")
print(f"  (Should be 3/2 = {1.5:.6f})")


# ================================================================
print("\n" + "=" * 70)
print("FINAL CRITICAL ASSESSMENT")
print("=" * 70)
print("""
  PATH 18 TAUTOLOGY DIAGNOSIS:

  1. TAUTOLOGICAL (no new content):
     - Q_K(A^4) = 3/2 was ENFORCED by choosing g0_tau (line 118-123)
     - The "observation chain" showing Q_K varies from 2.955 to 1.500
       is just the same constraint viewed in different coordinates
     - Saying "Koide is at the A^2 level" just means sqrt(m) = A^2,
       which is the DEFINITION of the mass-amplitude relation

  2. PARTIALLY NON-TRIVIAL (descriptive, not predictive):
     - The G(g0) enhancement curve IS a real property of the ODE
     - It tells us HOW the nonlinearity redistributes values
     - But it doesn't PREDICT Q_K = 3/2

  3. THE GENUINE QUESTION (addressed in Test E):
     - Does pure phi-FP spacing + ODE produce Q_K ~ 3/2?
     - If yes: Koide is emergent (real observer effect)
     - If no: Koide is an additional input (g0_tau tuned)
""")
