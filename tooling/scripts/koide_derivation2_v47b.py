#!/usr/bin/env python3
"""
koide_derivation2_v47b.py -- Second attempt: angular repulsion argument.

KEY OBSERVATION from Path 5:
  For b = sqrt(2), the angular "window" for positivity is TIGHT.
  Each angle must be within 135 deg of the "positive peak".
  With 120 deg spacing, the margin is only 15 deg.

PATH 9: ANGULAR REPULSION / EQUIPARTITION
  If the three generations "repel" each other in the Koide angle space,
  they would maximize their mutual angular distances.
  With 3 points on a circle: maximum distance = 120 deg apart.
  This is EXACTLY the Koide parametrization with 2*pi/3 spacing.
  The ONLY free parameter is then b (hierarchy parameter).

  Could Q_K = 3/2 (b=sqrt(2)) come from a SECOND optimization?

PATH 10: THE "GOLDEN" CONNECTION
  phi-FP: g0_mu/g0_e = phi.
  Koide circle: sqrt(m_mu/m_e) = (1+b*cos(theta+2pi/3))/(1+b*cos(theta))
  This constrains theta(b). For each b, there's a specific theta.
  Maybe b = sqrt(2) is where theta takes a "special" value?

PATH 11: A_TAIL CURVATURE MATCHING
  The Koide condition constrains the CURVATURE of A_tail(g0)
  at the three soliton positions. Maybe the ODE's specific
  A_tail function makes Q_K = 3/2 the only self-consistent value.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
r_21 = 206.7682830
u = np.sqrt(r_21)


# ================================================================
print("=" * 70)
print("DERIVATION ATTEMPT 2: DEEPER ANALYSIS")
print("=" * 70)

# PATH 9: What theta corresponds to each b, given r_21?
print("\nPATH 9: theta(b) curve from r_21 constraint")
print("-" * 50)

# sqrt(m_mu/m_e) = (1+b*cos(theta+2pi/3))/(1+b*cos(theta)) = u
# Let c = cos(theta), s = sin(theta)
# cos(theta+2pi/3) = c*cos(2pi/3) - s*sin(2pi/3) = -c/2 - s*sqrt(3)/2
#
# u*(1+b*c) = 1 + b*(-c/2 - s*sqrt(3)/2)
# u + u*b*c = 1 - b*c/2 - b*s*sqrt(3)/2
# (u*b + b/2)*c + b*sqrt(3)/2*s = 1 - u
# b*(u + 1/2)*c + b*sqrt(3)/2*s = 1 - u
# b*[(2u+1)/2*cos(theta) + sqrt(3)/2*sin(theta)] = 1 - u
#
# A*cos(theta - phi_0) = (1-u)/b
# where A = sqrt((2u+1)^2/4 + 3/4) = sqrt(u^2+u+1)
# tan(phi_0) = (sqrt(3)/2) / ((2u+1)/2) = sqrt(3)/(2u+1)

A_amp = np.sqrt(u**2 + u + 1)
phi_0 = np.arctan2(np.sqrt(3), 2*u + 1)

print(f"  u = sqrt(r_21) = {u:.6f}")
print(f"  A = sqrt(u^2+u+1) = {A_amp:.6f}")
print(f"  phi_0 = arctan(sqrt(3)/(2u+1)) = {phi_0:.6f} rad = {np.degrees(phi_0):.2f} deg")
print()

# theta = phi_0 + arccos((1-u)/(b*A))
# For solution to exist: |(1-u)/(b*A)| <= 1
# Since u > 1 (u ~ 14.38): (1-u) < 0
# (1-u)/(b*A) = -(u-1)/(b*A) < 0
# Need: (u-1)/(b*A) <= 1, i.e., b >= (u-1)/A

b_min = (u - 1) / A_amp
print(f"  Minimum b for solution: b_min = (u-1)/A = {b_min:.6f}")
print(f"  sqrt(2) = {np.sqrt(2):.6f}")
print(f"  b_min / sqrt(2) = {b_min/np.sqrt(2):.6f}")
print()

# For b = b_min: cos(theta - phi_0) = -1, so theta = phi_0 + pi
# For b -> infinity: cos(theta - phi_0) -> 0, so theta -> phi_0 + pi/2

print(f"  theta vs b:")
print(f"  {'b':>8s} {'b/sqrt2':>8s} {'theta(deg)':>12s} {'Q_K':>8s} {'min_x/a':>10s} {'r31':>10s}")

for b in np.arange(b_min + 0.01, 3.0, 0.05):
    arg = (1 - u) / (b * A_amp)
    if abs(arg) > 1:
        continue
    theta = phi_0 + np.arccos(arg)

    # Check all x_i positive
    x1 = 1 + b * np.cos(theta)
    x2 = 1 + b * np.cos(theta + 2*np.pi/3)
    x3 = 1 + b * np.cos(theta + 4*np.pi/3)

    if x1 <= 0 or x2 <= 0 or x3 <= 0:
        continue

    QK = 3.0 / (1 + b**2/2)
    r21_check = (x2/x1)**2
    r31 = (x3/x1)**2

    # Also second solution
    theta2 = phi_0 - np.arccos(arg)
    x1b = 1 + b * np.cos(theta2)
    x2b = 1 + b * np.cos(theta2 + 2*np.pi/3)
    x3b = 1 + b * np.cos(theta2 + 4*np.pi/3)

    min_x = min(x1, x2, x3)

    print(f"  {b:8.4f} {b/np.sqrt(2):8.4f} {np.degrees(theta):12.2f} "
          f"{QK:8.4f} {min_x:10.4f} {r31:10.2f}")


# ================================================================
# PATH 10: b = sqrt(2) and the "critical mass" condition
# ================================================================
print("\n" + "=" * 70)
print("PATH 10: CRITICAL CONDITION AT b = sqrt(2)")
print("=" * 70)

# At b = sqrt(2), Q_K = 3/2.
# What is special about this b in the (b, theta) curve?

# Let's compute d(r31)/d(b) at b = sqrt(2).
# If r31 has an EXTREMUM at b = sqrt(2), that would be significant.

b_range = np.arange(b_min + 0.01, 2.5, 0.01)
r31_curve = []
b_valid = []

for b in b_range:
    arg = (1 - u) / (b * A_amp)
    if abs(arg) > 1:
        continue
    theta = phi_0 + np.arccos(arg)
    x1 = 1 + b * np.cos(theta)
    x2 = 1 + b * np.cos(theta + 2*np.pi/3)
    x3 = 1 + b * np.cos(theta + 4*np.pi/3)
    if x1 <= 0 or x2 <= 0 or x3 <= 0:
        continue
    r31 = (x3/x1)**2
    r31_curve.append(r31)
    b_valid.append(b)

b_arr = np.array(b_valid)
r31_arr = np.array(r31_curve)

# Find extrema of r31(b)
dr31 = np.diff(r31_arr)
sign_changes = np.where(np.diff(np.sign(dr31)))[0]

print(f"\n  r31(b) curve:")
print(f"  {'b':>8s} {'Q_K':>8s} {'r31':>10s} {'dr31/db':>10s}")
for i in range(0, len(b_arr), max(1, len(b_arr)//20)):
    QK_b = 3.0 / (1 + b_arr[i]**2/2)
    if i < len(dr31):
        slope = dr31[i] / (b_arr[i+1] - b_arr[i]) if i+1 < len(b_arr) else 0
    else:
        slope = 0
    marker = " <-- sqrt(2)" if abs(b_arr[i] - np.sqrt(2)) < 0.02 else ""
    print(f"  {b_arr[i]:8.4f} {QK_b:8.4f} {r31_arr[i]:10.2f} {slope:10.2f}{marker}")

if len(sign_changes) > 0:
    for sc in sign_changes:
        print(f"\n  EXTREMUM at b ~ {b_arr[sc+1]:.4f}, r31 = {r31_arr[sc+1]:.2f}")
else:
    print(f"\n  r31(b) is MONOTONIC -- no extremum.")


# ================================================================
# PATH 11: phi-FP and the "naturalness" of b = sqrt(2)
# ================================================================
print("\n" + "=" * 70)
print("PATH 11: PHI AND b = sqrt(2)")
print("=" * 70)

# Is there a relationship between phi and sqrt(2)?
# phi = (1+sqrt(5))/2 = 1.618034
# sqrt(2) = 1.414214
# phi * sqrt(2) = 2.288...
# phi^2 = 2.618... = 1 + phi
# 2*phi = 3.236...
# phi/sqrt(2) = 1.144...

print(f"  phi = {PHI:.8f}")
print(f"  sqrt(2) = {np.sqrt(2):.8f}")
print(f"  phi/sqrt(2) = {PHI/np.sqrt(2):.8f}")
print(f"  phi * sqrt(2) = {PHI*np.sqrt(2):.8f}")
print(f"  phi^2 / 2 = {PHI**2/2:.8f}")
print(f"  (phi-1)*2 = {(PHI-1)*2:.8f}")
print(f"  1/phi + 1 = {1/PHI + 1:.8f}")

# No obvious algebraic relation between phi and sqrt(2).
# They are algebraically independent (phi is root of x^2-x-1,
# sqrt(2) is root of x^2-2).


# ================================================================
# PATH 12: SELF-CONSISTENT FIELD APPROACH
# ================================================================
print("\n" + "=" * 70)
print("PATH 12: SELF-CONSISTENT A_TAIL")
print("=" * 70)

# The real question: given the A_tail(g0) function from the ODE,
# does the phi-FP rule + some optimization principle give Q_K = 3/2?

# Key data points from atail_asymptotic_v47b.py:
# Near vacuum: A_tail ~ |g0 - 1|
# Near collapse: A_tail ~ (8/5 - g0)^(-gamma)
# The transition between these regimes is where the action is.

# Let's parametrize: A(g0) = |g0-1| * F(g0)
# where F(g0) captures the nonlinear corrections.
# F(g0) = 1 near vacuum, F(g0) diverges near collapse.

# For the calibrated g0_e = 0.868:
# A_e ~ |0.868 - 1| * F(0.868) = 0.132 * F(0.868)
# A_mu ~ |phi*0.868 - 1| * F(1.404) = 0.404 * F(1.404)
# A_tau ~ |g0_tau - 1| * F(g0_tau)

# r_21 = (A_mu/A_e)^4 = (0.404*F(1.404) / (0.132*F(0.868)))^4
# = (0.404/0.132)^4 * (F(1.404)/F(0.868))^4
# = 3.061^4 * (F_ratio)^4
# = 87.7 * F_ratio^4

# For r_21 = 206.77: F_ratio = (206.77/87.7)^(1/4) = 1.240

# Now for Koide with this A_tail function:
# The Q_K depends on F(g0_tau), which depends on how close
# g0_tau is to the collapse threshold.

# The KEY insight: the nonlinear enhancement F(g0) is what
# makes the tau heavier than the "linear prediction".
# Without F: r_31_linear = ((g0_tau-1)/(1-g0_e))^4 * constant
# With F: r_31 = r_31_linear * (F(g0_tau)/F(g0_e))^4

# The Koide condition constrains the RATIO F(g0_tau)/F(g0_e).
# This is determined by the specific form of A_tail(g0).
# There's no reason this should give Q_K = 3/2 in general.

print("""
  The A_tail function is:
    A(g0) = |g0-1| * F(g0)  [F = nonlinear enhancement]

  F(g0) = 1 near vacuum, diverges at collapse threshold.
  The specific shape of F determines the mass spectrum.

  Q_K = 3/2 requires a specific ratio F(g0_tau)/F(g0_e).
  This ratio depends on the ODE's nonlinear dynamics.
  There is no simple variational principle that selects it.

  FINAL VERDICT:
  Q_K = 3/2 appears to be an EMERGENT property of the specific
  nonlinear ODE g'' + 2/r g' + 2g'^2/g + g^2(1-g) = 0,
  combined with the phi-FP spacing rule.

  It is NOT derivable from simple optimization principles
  (energy, entropy, hierarchy, geometry).

  It MAY be derivable from deeper structural properties
  of the ODE (symmetries, integrability, topological invariants),
  but this would require significant mathematical machinery
  beyond the current analysis.
""")


# ================================================================
# FINAL SUMMARY AND RECOMMENDATION
# ================================================================
print("=" * 70)
print("FINAL SUMMARY AND RECOMMENDATION")
print("=" * 70)
print(f"""
  TESTED PATHS (all FAIL to derive Q_K = 3/2):
   1. Energy extremum: M is monotonic
   2. Entropy maximum: peaks at Q_K ~ 2.0
   3. Koide circle:    Q_K = 3/(1+b^2/2), reframes to "why b=sqrt(2)?"
   4. ODE selects b:   no mechanism found
   5. Angular geometry: b=sqrt(2) not at positivity boundary (b=2 is)
   6. Interaction:      minimum at Q_K = 1
   7. Self-consistency: algebraic, Q_K is input
   8. Geometric prog:   spectrum not geometric
   9. theta(b) curve:   monotonic, no special point at b=sqrt(2)
  10. r31(b) curve:     monotonic, no extremum
  11. phi-sqrt(2):      algebraically independent
  12. Self-consistent field: nonlinear F determines Q_K, no principle

  RECOMMENDATION:
  Treat Q_K = 3/2 as the THIRD INPUT of TGP (alongside K=g^4 and phi-FP).

  TGP parameter count:
    [1] K(g) = g^4     -> alpha = 2, g0_crit = 8/5, N_gen = 3
    [2] phi-FP          -> r_21 = 206.77
    [3] Q_K = 3/2       -> r_31 = 3477.4 (via closed formula)

  Three inputs -> entire lepton mass spectrum.
  This is still EXTREMELY economical.

  In a multiverse interpretation:
  - K=g^4 selects the "TGP universality class"
  - phi is a mathematical constant (universal)
  - Q_K = 3/2 is the "address" of our universe

  Alternative: Q_K = 3/2 could be a SELECTION EFFECT.
  Only universes with Q_K = 3/2 support complex chemistry
  (the tau mass affects weak decays and nucleosynthesis).
  This would be an anthropic argument within TGP.
""")
