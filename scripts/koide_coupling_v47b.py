#!/usr/bin/env python3
"""
koide_coupling_v47b.py -- PATH 14: WHY self/cross = 4?

From koide_field_balance_v47b.py we found:
  Q_K = 3/2 <=> sum(A_i^4) = 4 * sum_{i<j}(A_i^2 * A_j^2)

The number 4 appears. In TGP, alpha = 2 and:
  - K(g) = g^4            (alpha+2 = 4)
  - g^(alpha+1) = g^3     (canonical variable)
  - m ~ A^4               (mass-amplitude relation is FOURTH power)

Is the "4" in self/cross = 4 related to the FOURTH-power mass relation?

KEY REFORMULATION:
  Let m_i = c * A_i^4 (the TGP mass relation).
  Then A_i^2 = (m_i/c)^(1/2) = sqrt(m_i) * c^(-1/2).
  And A_i^4 = m_i / c.

  sum(A_i^4) = (1/c) * sum(m_i)
  sum(A_i^2 * A_j^2) = (1/c) * sum(sqrt(m_i)*sqrt(m_j))

  Wait: A_i^2 * A_j^2 = sqrt(m_i/c) * sqrt(m_j/c) = sqrt(m_i*m_j)/c

  self/cross = sum(m_i) / sum(sqrt(m_i*m_j))
             = sum(m_i) / sum_{i<j}(sqrt(m_i*m_j))

  For Q_K = 3/2:
  (sum sqrt(m))^2 = (3/2) * sum(m)
  sum(m) + 2*sum(sqrt(m_i*m_j)) = (3/2)*sum(m)
  2*sum(sqrt(m_i*m_j)) = (1/2)*sum(m)
  sum(sqrt(m_i*m_j)) = sum(m)/4

  So self/cross = sum(m) / (sum(m)/4) = 4. TAUTOLOGY CONFIRMED.

  But the PHYSICAL CONTENT is in the mass-amplitude relation m ~ A^4.
  If m ~ A^n instead, we'd get a DIFFERENT condition.

  For general power law m ~ A^n:
  A_i = (m_i/c)^(1/n)
  sum(A_i^4) = sum((m_i/c)^(4/n))
  sum(A_i^2*A_j^2) = sum((m_i*m_j)^(2/n))/c^(4/n)

  For n=4: A_i^4 = m_i/c, A_i^2*A_j^2 = sqrt(m_i*m_j)/c
  self/cross = sum(m) / sum(sqrt(m_i*m_j)) = 4 (for Koide)

  For n=2: A_i^4 = m_i^2/c^2, A_i^2*A_j^2 = m_i*m_j/c^2
  self/cross = sum(m^2) / sum(m_i*m_j) -- DIFFERENT condition

  So the "4" is ENTANGLED with the fourth-power mass relation.
  The TGP-specific prediction is: "self = 4*cross" in A-space,
  which via m ~ A^4 becomes exactly Koide Q_K = 3/2.

  NEW QUESTION: Is there a PHYSICAL reason why self = 4*cross
  in the soliton tail amplitude space?

  IDEA: In the ODE g'' + 2/r g' + (alpha/g)g'^2 + g^2(1-g) = 0,
  the nonlinear coupling term (alpha/g)g'^2 has coefficient alpha=2.
  The tail amplitude A satisfies matching conditions at the core-tail
  boundary. Could the coupling between different soliton tails
  (via the shared vacuum) impose self = 4*cross?

  DEEPER IDEA: The three solitons share the SAME vacuum g=1.
  Their tails overlap in this vacuum. The condition for the vacuum
  to be self-consistent under the perturbations from all three solitons
  might be what gives self = 4*cross.

  TEST: For different values of alpha, does the ODE still give Q_K = 3/2?
  We already know r_31 is universal wrt alpha (F_alpha_canonical).
  But is Q_K also universal? YES -- because Q_K comes from the mass ratios
  which come from A^4 ratios which are calibrated to r_21 and r_31.

  So the real question: WHY does the phi-FP + r_21 + ODE combination
  produce amplitudes satisfying self = 4*cross?
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2


def solver_general(g0, alpha, r_max=400):
    """Solve canonical Form A ODE for general alpha."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0]


def A_tail_general(g0, alpha, r_min=50, r_max=300):
    r, g = solver_general(g0, alpha, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


def calibrate_g0e(alpha, r21_target=206.7683, g0e_lo=0.5, g0e_hi=0.99):
    """Find g0_e such that (A_mu/A_e)^4 = r21, with g0_mu = phi*g0_e."""
    def residual(g0e):
        g0m = PHI * g0e
        gc = 2*(alpha+2) / (2*(alpha+2) - 3)
        if g0m >= gc:
            return 1e10
        Ae = A_tail_general(g0e, alpha)
        Am = A_tail_general(g0m, alpha)
        if Ae < 1e-15:
            return 1e10
        return (Am/Ae)**4 - r21_target
    try:
        return brentq(residual, g0e_lo, g0e_hi, xtol=1e-10)
    except:
        return None


# ================================================================
print("=" * 70)
print("PATH 14: WHY SELF/CROSS = 4?")
print("=" * 70)

# ================================================================
# Section 1: Verify self/cross = 4 is equivalent to Q_K = 3/2
# ================================================================
print("\nSECTION 1: ALGEBRAIC EQUIVALENCE")
print("-" * 50)

print("""
  For variables x_i = A_i^2 (proportional to sqrt(m_i)):

  Q_K = (sum x_i)^2 / sum(x_i^2)

  Let S = sum(x_i), P = sum(x_i^2), C = sum_{i<j}(x_i*x_j)
  Then S^2 = P + 2C, so Q_K = (P + 2C)/P = 1 + 2C/P

  Q_K = 3/2 => 2C/P = 1/2 => P = 4C => self = 4 * cross. QED.

  Now: via m_i = k * A_i^4 (TGP mass relation):
    x_i = A_i^2 = (m_i/k)^(1/2) = sqrt(m_i/k)
    P = sum(m_i/k), C = sum(sqrt(m_i*m_j)/k)
    self/cross = sum(m_i) / sum(sqrt(m_i*m_j))

  This ratio depends ONLY on mass ratios, not on k.
  So self/cross = 4 is a UNIVERSAL CONSTRAINT on mass ratios
  that is INDEPENDENT of the amplitude-mass power law!

  Wait -- that can't be right. Let me recheck for general power n.
""")

# For general mass-amplitude law m = k * A^n:
# x_i = A_i^2 = (m_i/k)^(2/n)
# P = sum(x_i^2) = sum((m_i/k)^(4/n))
# C = sum_{i<j}(x_i*x_j) = sum_{i<j}((m_i*m_j)^(2/n))/k^(4/n)
#
# Actually the Q_K formula uses sqrt(m_i), not A_i^2.
# sqrt(m_i) = sqrt(k) * A_i^(n/2)
# Q_K = (sum sqrt(m_i))^2 / sum(m_i)
# This is fixed by the actual masses, regardless of n.
# Q_K = 3/2 is a statement about MASSES, period.
#
# But the "self = 4*cross in A-space" is ALSO about masses,
# because A_i^2 = sqrt(m_i/k), so A_i^4 = m_i/k.
# self/cross in A-space = sum(m_i) / sum(sqrt(m_i*m_j)) = 4.
#
# This is the SAME condition regardless of n, because the
# A^2 vs sqrt(m) relation is always: A^2 proportional to m^(1/2)
# when m = k * A^4 (TGP). But for m = k * A^n:
# A^2 = (m/k)^(2/n), so A^4 = (m/k)^(4/n)
# self = sum(A^4) = sum(m^(4/n))/k^(4/n)
# cross = sum(A_i^2*A_j^2) = sum((m_i*m_j)^(2/n))/k^(4/n)
#
# For n=4: self = sum(m)/k, cross = sum(sqrt(m_i*m_j))/k
#   self/cross = sum(m)/sum(sqrt(m_i*m_j)) = 4 (Koide)
#
# For n=2: self = sum(m^2)/k^2, cross = sum(m_i*m_j)/k^2
#   self/cross = sum(m^2)/sum(m_i*m_j)
#   If Koide holds: sum(m) = (2/3)(sum sqrt(m))^2
#   sum(m^2)/sum(m_i*m_j) is NOT automatically 4.

print("  For n=4 (TGP): self/cross in A-space = sum(m)/sum(sqrt(m_i*m_j))")
print("  Koide => this ratio = 4.")
print()
print("  For n=2 (hypothetical): self/cross = sum(m^2)/sum(m_i*m_j)")

# Compute for actual masses
m_e = 0.511  # MeV
m_mu = 105.658
m_tau = 1776.86

# n=4 (TGP)
self_n4 = m_e + m_mu + m_tau
cross_n4 = np.sqrt(m_e*m_mu) + np.sqrt(m_e*m_tau) + np.sqrt(m_mu*m_tau)
print(f"\n  n=4: self/cross = {self_n4:.2f}/{cross_n4:.2f} = {self_n4/cross_n4:.6f}")

# n=2
self_n2 = m_e**2 + m_mu**2 + m_tau**2
cross_n2 = m_e*m_mu + m_e*m_tau + m_mu*m_tau
print(f"  n=2: self/cross = {self_n2:.2f}/{cross_n2:.2f} = {self_n2/cross_n2:.6f}")

# n=6
self_n6 = m_e**(4/6) + m_mu**(4/6) + m_tau**(4/6)
cross_n6 = (m_e*m_mu)**(2/6) + (m_e*m_tau)**(2/6) + (m_mu*m_tau)**(2/6)
print(f"  n=6: self/cross = {self_n6:.4f}/{cross_n6:.4f} = {self_n6/cross_n6:.6f}")

# n=1
self_n1 = m_e**4 + m_mu**4 + m_tau**4
cross_n1 = (m_e*m_mu)**2 + (m_e*m_tau)**2 + (m_mu*m_tau)**2
print(f"  n=1: self/cross = {self_n1:.2e}/{cross_n1:.2e} = {self_n1/cross_n1:.6f}")


# ================================================================
# Section 2: PROTON ANALOGY - three components with "fractional" properties
# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: PROTON ANALOGY -- FRACTIONAL DECOMPOSITION")
print("=" * 70)

print("""
  Proton: charge = +2/3 + 2/3 - 1/3 = +1
  The "whole" (proton charge = 1) decomposes into fractional parts.

  Lepton family: what quantity is "whole"?
  Candidates:
    (a) sum(sqrt(m_i)) = 3*a  (from Koide parametrization)
    (b) sum(m_i) = 3*a^2*(1+b^2/2)  (from Koide)
    (c) Q_K = 3/2 (the "charge" of the lepton family)

  For Koide with b=sqrt(2):
    sum(sqrt(m)) = 3a
    sum(m) = 3a^2 * 2 = 6a^2
    Q_K = 9a^2 / 6a^2 = 3/2

  The "fractional" sqrt(m_i) values:
""")

sqm = np.array([np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)])
a_koide = np.sum(sqm) / 3
print(f"  a = sum(sqrt(m)) / 3 = {a_koide:.6f}")
print(f"  sqrt(m_e)/a   = {sqm[0]/a_koide:.6f}  (fractional part)")
print(f"  sqrt(m_mu)/a  = {sqm[1]/a_koide:.6f}  (fractional part)")
print(f"  sqrt(m_tau)/a = {sqm[2]/a_koide:.6f}  (fractional part)")
print(f"  Sum = {np.sum(sqm/a_koide):.6f}  (= 3, by construction)")

# In Koide parametrization: sqrt(m_i)/a = 1 + b*cos(theta_i)
# where theta_i = theta + 2pi(i-1)/3
# So the "fractional" values are 1 + sqrt(2)*cos(theta_i)
b_koide = np.sqrt(2)
theta_K = np.arctan2(np.sqrt(3)*(sqm[1]/a_koide - sqm[0]/a_koide),
                     (sqm[0]/a_koide + sqm[1]/a_koide - 2))

# Compute from actual parametrization
print(f"\n  Koide decomposition (b = sqrt(2)):")
for i, name in enumerate(['e', 'mu', 'tau']):
    frac = sqm[i]/a_koide
    delta_i = frac - 1  # deviation from "democracy"
    print(f"    {name}: sqrt(m)/a = {frac:.6f} = 1 + {delta_i:+.6f}")

print(f"\n  Sum of deviations: {np.sum(sqm/a_koide - 1):.10f} (= 0 by construction)")
print(f"  Sum of deviations^2: {np.sum((sqm/a_koide - 1)**2):.6f}")
print(f"  2/3 * sum(dev^2) = b^2 = {2/3*np.sum((sqm/a_koide - 1)**2):.6f}")
print(f"  sqrt(2)^2 = {2:.6f}")

# KEY: the deviations sum to 0 (like quark charges?)
# AND their "energy" (sum dev^2) equals (3/2)*b^2 = 3.
# So Q_K = 3/2 means "the deviations have unit variance" (CV=1)


# ================================================================
# Section 3: THE VACUUM AS MEDIATOR -- field theory perspective
# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: VACUUM AS MEDIATOR")
print("=" * 70)

print("""
  Physical picture (user's hypothesis refined):

  1. The TGP vacuum (g=1) is a SHARED MEDIUM for all solitons.
  2. Each soliton deforms the vacuum with amplitude A_i.
  3. The vacuum's response is governed by the linearized ODE
     (spherical Bessel equation around g=1).
  4. The three solitons must be COMPATIBLE with the same vacuum.

  Condition: the total "stress" on the vacuum from all three
  solitons must satisfy a consistency relation.

  The soliton tail h_i(r) ~ A_i * sin(r - delta_i) / r
  satisfies the linearized equation. But the FULL nonlinear equation
  has terms O(h^2) that couple different modes.

  At second order: the nonlinear coupling generates terms
  proportional to A_i * A_j. These must be self-consistent
  with the vacuum equation.

  The condition for self-consistency at O(A^2):
    sum of second-order corrections must vanish
    (or: source for h^(2) must have no secular growth)

  This is like the SOLVABILITY CONDITION in perturbation theory!
""")

# Let me compute the nonlinear coupling strength
# h'' + 2/r h' - h = N(h)  where N(h) are nonlinear terms
# N(h) = 2h'^2/(1+h) - 2h^2 - h^3
# At second order: N^(2)(h) = 2h'^2 - 2h^2 (dropping h^3, using 1/(1+h)~1-h)
#
# If h = sum_i h_i = sum_i A_i sin(r-d_i)/r, then
# h'^2 = (sum A_i [cos(r-d_i)/r - sin(r-d_i)/r^2])^2
# The leading term at large r:
# h'^2 ~ (sum A_i cos(r-d_i)/r)^2 = sum_i,j A_iA_j cos(r-d_i)cos(r-d_j)/r^2
# h^2 ~ (sum A_i sin(r-d_i)/r)^2 = sum_i,j A_iA_j sin(r-d_i)sin(r-d_j)/r^2
#
# N^(2) = 2h'^2 - 2h^2
# = 2*sum_{i,j} A_iA_j [cos(r-d_i)cos(r-d_j) - sin(r-d_i)sin(r-d_j)] / r^2
# = 2*sum_{i,j} A_iA_j cos(2r - d_i - d_j) / r^2  [WRONG: this uses cos-cos - sin-sin = cos(sum)]
# Wait: cos(a)cos(b) - sin(a)sin(b) = cos(a+b)
# So N^(2) = 2*sum_{i,j} A_iA_j cos((r-d_i)+(r-d_j)) / r^2
# = 2*sum_{i,j} A_iA_j cos(2r - d_i - d_j) / r^2
#
# The source N^(2) has frequency 2 (double frequency).
# The free equation h'' + 2/r h' - h = 0 has characteristic exponents
# related to frequency 1 (for the oscillatory case).
# At frequency 2, the particular solution exists without secular growth
# IF the frequency 2 is not a resonance of the operator.
#
# Actually, the linearized operator is h'' + 2/r h' - h = 0
# which gives exponential (not oscillatory) tails.
# The oscillatory behavior comes from the canonical variable f = g^3.
# Let me redo this in the f-variable.

print("  In canonical variable f = g^(alpha+1) = g^3:")
print("  Linearized: eta'' + 2/r eta' + eta = 0  (oscillatory!)")
print("  Solutions: eta ~ sin(r)/r, cos(r)/r")
print("  Frequency omega = 1.")
print()
print("  At second order, nonlinear source has frequency 2.")
print("  The Green's function at frequency 2 is non-singular")
print("  (2 is not a resonance of the spherical Bessel operator).")
print("  So there's NO solvability condition at O(A^2)!")
print()
print("  The condition would arise at a HIGHER order,")
print("  or from a DIFFERENT mechanism (e.g., energy conservation).")


# ================================================================
# Section 4: ENERGY PARTITION TEST
# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: ENERGY PARTITION -- does the ODE predict Q_K?")
print("=" * 70)

# Test: for different alpha, compute the ACTUAL Q_K from the ODE
# when g0_tau is determined by phi-FP (g0_tau = phi^2 * g0_e)
# versus when g0_tau gives Q_K = 3/2.

# The key question: does any natural rule for g0_tau
# (not the Koide rule) give Q_K = 3/2?

print(f"\n  Testing natural rules for g0_tau:")
print(f"  g0_e  = 0.86770, g0_mu = phi*g0_e = {PHI*0.86770:.5f}")
print()

g0_e = 0.86770494
g0_mu = PHI * g0_e
A_e = A_tail_general(g0_e, 2)
A_mu = A_tail_general(g0_mu, 2)

rules = {
    'phi^2 * g0_e': PHI**2 * g0_e,
    'phi * g0_mu': PHI * g0_mu,
    '2 * g0_e': 2 * g0_e,
    'g0_mu + (g0_mu - g0_e)': g0_mu + (g0_mu - g0_e),
    'sqrt(2) * g0_mu / g0_e * g0_e': np.sqrt(2) * g0_mu,
    'g0_crit - (g0_crit - g0_mu)^2/(g0_crit - g0_e)': None,
}

# Compute the "geometric" rule: g0_crit minus something
gc = 8/5
rules['g0_crit - (gc-g0_mu)^2/(gc-g0_e)'] = gc - (gc - g0_mu)**2 / (gc - g0_e)

# Energy equipartition: E_tau = E_mu  or  E_tau = E_e + E_mu
# First compute soliton energies properly
def soliton_field_energy(g0, alpha=2):
    r, g = solver_general(g0, alpha, r_max=300)
    mask = (r > 0.001) & (r < 250)
    rf = r[mask]
    gf = g[mask]
    gp = np.gradient(gf, rf)
    # (g-1)^2 as proxy for field energy density
    integrand = (gf - 1.0)**2 * rf**2
    return trapz(integrand, rf)

E_field_e = soliton_field_energy(g0_e)
E_field_mu = soliton_field_energy(g0_mu)
print(f"  Field energy: E_e = {E_field_e:.4f}, E_mu = {E_field_mu:.4f}")

# What g0_tau gives E_tau = E_e + E_mu?  (energy conservation)
print(f"\n  Searching: E_tau = E_e + E_mu = {E_field_e + E_field_mu:.4f}")
def energy_balance(g0t):
    return soliton_field_energy(g0t) - (E_field_e + E_field_mu)

# Scan
print(f"\n  {'rule':>45s} {'g0_tau':>8s} {'Q_K':>8s} {'r31':>10s}")
for name, g0t in rules.items():
    if g0t is None or g0t >= gc or g0t <= g0_mu:
        continue
    At = A_tail_general(g0t, 2)
    if At < 1e-10:
        continue
    r31 = (At/A_e)**4
    xe, xm, xt = A_e**2, A_mu**2, At**2
    QK = (xe+xm+xt)**2 / (xe**2+xm**2+xt**2)
    marker = " <-- KOIDE!" if abs(QK - 1.5) < 0.02 else ""
    print(f"  {name:>45s} {g0t:8.5f} {QK:8.4f} {r31:10.2f}{marker}")

# Energy-balance rule
for g0t in np.arange(1.42, 1.59, 0.002):
    Ef = soliton_field_energy(g0t)
    if abs(Ef - (E_field_e + E_field_mu)) < 1.0:
        At = A_tail_general(g0t, 2)
        r31 = (At/A_e)**4
        xe, xm, xt = A_e**2, A_mu**2, At**2
        QK = (xe+xm+xt)**2 / (xe**2+xm**2+xt**2)
        print(f"  {'E_tau = E_e + E_mu':>45s} {g0t:8.5f} {QK:8.4f} {r31:10.2f}"
              + (" <-- KOIDE!" if abs(QK - 1.5) < 0.02 else ""))


# ================================================================
# Section 5: DEMOCRATIC DEFICIT (CV = 1) from field theory
# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: DEMOCRATIC DEFICIT / CV = 1")
print("=" * 70)

print("""
  Q_K = 3/2 <=> CV(sqrt(m)) = 1 <=> b = sqrt(2)
  <=> "variance of sqrt(m) = mean^2"

  In field theory terms:
    <A^2> = (A_e^2 + A_mu^2 + A_tau^2) / 3  (mean of A^2)
    Var(A^2) = <(A^2)^2> - <A^2>^2
             = (A_e^4 + A_mu^4 + A_tau^4)/3 - ((A_e^2+A_mu^2+A_tau^2)/3)^2

    CV(A^2) = sqrt(Var) / <A^2>

    Q_K = 3/2 <=> CV(A^2) = 1 <=> Var = <A^2>^2
    <=> <(A^2)^2> = 2 * <A^2>^2  (i.e., kurtosis condition)

  This is a STATISTICAL property of the amplitude distribution.
  It says: the three amplitudes are distributed such that
  their coefficient of variation is EXACTLY 1.

  Physical meaning: the fluctuations in soliton "strength"
  (measured by A^2) are as large as the mean strength itself.
  This is the MAXIMUM DISORDER compatible with all amplitudes > 0
  under the constraint of 3-fold structure... NO, it's not the maximum.
  CV can be larger.

  But CV = 1 IS special: it's the value where the distribution
  is "scale-free" in a certain sense. The standard deviation
  equals the mean.
""")

# Compute the actual CV
A_vals = np.array([A_e, A_mu, A_tau_val := A_tail_general(g0_tau_koide := 1.56962743, 2)])
x_vals = A_vals**2
mean_x = np.mean(x_vals)
var_x = np.var(x_vals)
cv_x = np.sqrt(var_x) / mean_x

print(f"  Amplitude values: A_e={A_vals[0]:.6f}, A_mu={A_vals[1]:.6f}, A_tau={A_vals[2]:.6f}")
print(f"  x = A^2 values:  {x_vals[0]:.6f}, {x_vals[1]:.6f}, {x_vals[2]:.6f}")
print(f"  Mean(x) = {mean_x:.6f}")
print(f"  Var(x)  = {var_x:.6f}")
print(f"  CV(x)   = {cv_x:.6f}  (Q_K = 3/2 requires CV = 1.0)")
print(f"  CV^2    = {cv_x**2:.6f}")


# ================================================================
# Section 6: THE PROTON-LEPTON PARALLEL (quantitative)
# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: QUANTITATIVE PROTON-LEPTON PARALLEL")
print("=" * 70)

print("""
  Proton:  3 quarks with charges q_i = (+2/3, +2/3, -1/3)
    sum(q_i) = +1 (integer "neutrality")
    sum(q_i^2) = 4/9 + 4/9 + 1/9 = 1
    Q_charge = sum(q_i)^2 / sum(q_i^2) = 1/1 = 1

  Lepton family: 3 generations with "charges" x_i = sqrt(m_i/m_e)
    sum(x_i) = S
    sum(x_i^2) = P (= sum(m_i/m_e))
    Q_K = S^2/P = 3/2

  In Koide parametrization with b=sqrt(2):
    x_i = a(1 + sqrt(2)*cos(theta_i))
    where theta_i are 120 deg apart.

  The "fractional charges" are:
    f_i = x_i / sum(x_j) = (1 + sqrt(2)*cos(theta_i)) / 3
""")

# Compute lepton "fractional charges"
sqm_norm = np.array([np.sqrt(m_e/m_e), np.sqrt(m_mu/m_e), np.sqrt(m_tau/m_e)])
S = np.sum(sqm_norm)
f_i = sqm_norm / S

print(f"  Lepton 'fractional charges':  f_e = {f_i[0]:.6f}, f_mu = {f_i[1]:.6f}, f_tau = {f_i[2]:.6f}")
print(f"  Sum = {np.sum(f_i):.6f}")
print(f"  sum(f_i^2) = {np.sum(f_i**2):.6f}")
print(f"  Q = (sum f_i)^2 / (3*sum(f_i^2)) = {1/(3*np.sum(f_i**2)):.6f}")

# Quark analogy
q_proton = np.array([2/3, 2/3, -1/3])
print(f"\n  Quark charges: q = {q_proton}")
print(f"  Sum = {np.sum(q_proton):.6f}")
print(f"  sum(q^2) = {np.sum(q_proton**2):.6f}")
print(f"  Q = (sum q)^2 / (3*sum(q^2)) = {np.sum(q_proton)**2/(3*np.sum(q_proton**2)):.6f}")

# So Q_proton = 1/3 * (1)^2 / (9/9) = 1/3.
# Q_lepton = 1/3 * (1)^2 / (sum f^2).
# Q_K = 3/2 gives sum(f^2) = 2/9.

# Intersting: the proton has Q = 1/(3*sum(q^2)) = 1/(3*1) = 1/3
# The lepton family has Q = 1/(3*sum(f^2)) = 1/(3*2/9) = 9/6 = 3/2 ... wait
# Q_K = (sum sqrt(m))^2 / sum(m) = S^2/P
# f_i = sqrt(m_i)/S, so S = 1 in normalized units
# sum(f_i^2) = sum(m_i)/S^2 = P/S^2 = 1/Q_K
# So 1/(3*sum(f_i^2)) = Q_K/3 = 1/2.

print(f"\n  Interesting parallel:")
print(f"    Q_proton_charge = (sum q)^2 / (3*sum(q^2)) = 1/3")
print(f"    Q_lepton_mass   = Q_K/3 = 1/2")
print(f"    Both are simple fractions!")
print(f"    Proton: 1/3 (from 3 colors)")
print(f"    Lepton: 1/2 (from ... what?)")


# ================================================================
# FINAL SYNTHESIS
# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS")
print("=" * 70)

print(f"""
  The user's proton analogy is STRUCTURALLY apt:

  1. PROTON: 3 quarks, fractional charges (2/3, 2/3, -1/3), sum = 1
     Enforced by: SU(3) color confinement + charge quantization.

  2. LEPTON FAMILY: 3 generations, "fractional" sqrt(m) values,
     Q_K = sum(sqrt(m))^2 / sum(m) = 3/2.
     Enforced by: ??? (the open question)

  WHAT WE LEARNED:
  - Complex neutrality: FAILS (phases don't cancel)
  - Field charge neutrality: FAILS (all Q_i < 0)
  - Self = 4*cross: TRUE but tautological (= Q_K = 3/2)
  - Energy partition: no natural rule gives Q_K = 3/2
  - CV = 1: the cleanest reformulation (democratic deficit)

  THE DEEPEST INSIGHT from the proton analogy:
  In the proton, the NUMBER 3 (quarks) and charge quantization
  are enforced by the GAUGE GROUP SU(3).
  In TGP, the NUMBER 3 (generations) is enforced by the ODE dynamics
  (soliton collapse threshold). Could there be an ANALOGOUS
  symmetry principle that fixes Q_K = 3/2?

  POSSIBLE DIRECTION:
  The soliton ODE with alpha=2 has a specific symmetry structure.
  The canonical variable f=g^3 transforms it to Lane-Emden type.
  Lane-Emden equations have scaling symmetries (Lie group structure).
  Could a discrete subgroup of this scaling symmetry
  constrain the mass ratios to satisfy Q_K = 3/2?

  This would be the "SU(3) of TGP" -- not a gauge group,
  but a symmetry of the soliton equation that quantizes
  the mass spectrum.
""")
