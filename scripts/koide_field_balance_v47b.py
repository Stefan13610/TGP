#!/usr/bin/env python3
"""
koide_field_balance_v47b.py -- PATH 13: FIELD BALANCE / SOLITON NEUTRALITY

HYPOTHESIS (user's insight):
  The creation of lighter solitons DISTURBS the TGP field in a way
  that REQUIRES the existence of heavier solitons for balance.
  The three generations are like quarks in a proton:
    proton charge: +2/3 + 2/3 - 1/3 = +1  (three inseparable parts)
    lepton family: e + mu + tau = "neutral" (three inseparable generations)

  This "field balance" or "soliton neutrality" condition
  could be what fixes Q_K = 3/2.

WHAT TO TEST:
  Each soliton has an oscillatory tail: h(r) ~ A_i * sin(r - delta_i) / r
  The amplitude A_i and phase delta_i depend on g0_i.

  Possible neutrality conditions:
  (a) Sum of amplitudes: sum(A_i) = 0 (impossible, all positive)
  (b) Sum of signed amplitudes: sum(s_i * A_i) = 0 (with signs)
  (c) Complex neutrality: sum(A_i * exp(i*delta_i)) = 0
  (d) Energy balance: sum(A_i^2) = const or sum(A_i^2 * f(delta_i)) = 0
  (e) Quadratic sum rule: sum(A_i^2) = (sum(A_i))^2 / 3  [<=> Q_K = 3/2!]
  (f) Tail integral: integral of g-1 over all space sums to zero
  (g) Virial-like: sum of "field charge" Q_i proportional to A_i vanishes

  KEY INSIGHT:
  sqrt(m_i) ~ A_i^2 (from A_tail^4 = mass ratio).
  Actually m_i/m_e = (A_i/A_e)^4, so sqrt(m_i/m_e) = (A_i/A_e)^2.
  Let x_i = A_i^2. Then m_i ~ x_i^2.
  Q_K = (sum sqrt(m_i))^2 / sum(m_i) = (sum x_i)^2 / sum(x_i^2).

  Q_K = 3/2 means:
    (sum x_i)^2 / sum(x_i^2) = 3/2
    2*(x1+x2+x3)^2 = 3*(x1^2+x2^2+x3^2)
    2*S^2 = 3*P  where S = sum(x_i), P = sum(x_i^2)
    Expanding: x1^2+x2^2+x3^2 = 2*(x1*x2+x1*x3+x2*x3)/1
    Actually: 2(x1+x2+x3)^2 = 3(x1^2+x2^2+x3^2)
    => 2(x1^2+x2^2+x3^2+2x1x2+2x1x3+2x2x3) = 3(x1^2+x2^2+x3^2)
    => x1^2+x2^2+x3^2 = 4(x1x2+x1x3+x2x3)

  So Q_K = 3/2 <=> sum(x_i^2) = 4*sum(x_i*x_j) for i<j
  In terms of amplitudes: "self-energy = 4 * cross-energy"

  This looks like a COUPLING CONDITION between the three solitons!
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0


def solver_A(g0, r_max=400):
    """Solve canonical Form A ODE (alpha=2)."""
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


def extract_tail(g0, r_min=50, r_max=300):
    """Extract A_tail and delta (amplitude and phase) from soliton tail."""
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0, 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf  # r * (g-1) ~ A*sin(r - delta)
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    A = np.sqrt(bc[0]**2 + bc[1]**2)
    delta = np.arctan2(-bc[0], bc[1])  # phase: h ~ A*sin(r - delta)
    return A, delta


def field_integral(g0, r_max=300):
    """Compute integral of (g-1) * r^2 dr  (total 'field charge')."""
    r, g = solver_A(g0, r_max=r_max+50)
    mask = r < r_max
    rf = r[mask]
    gf = g[mask]
    integrand = (gf - 1.0) * rf**2
    return trapz(integrand, rf)


def field_integral_sq(g0, r_max=300):
    """Compute integral of (g-1)^2 * r^2 dr  (field energy proxy)."""
    r, g = solver_A(g0, r_max=r_max+50)
    mask = r < r_max
    rf = r[mask]
    gf = g[mask]
    integrand = (gf - 1.0)**2 * rf**2
    return trapz(integrand, rf)


# ================================================================
print("=" * 70)
print("PATH 13: FIELD BALANCE / SOLITON NEUTRALITY")
print("=" * 70)

# Calibration
g0_e = 0.86770494
g0_mu = PHI * g0_e

print(f"\n  g0_e  = {g0_e:.8f}")
print(f"  g0_mu = {g0_mu:.8f}")

# Find g0_tau from Koide Q_K = 3/2
# We know from previous: g0_tau ~ 1.5696
# Let's refine it
A_e, delta_e = extract_tail(g0_e)
A_mu, delta_mu = extract_tail(g0_mu)

print(f"\n  A_e    = {A_e:.8f},  delta_e  = {delta_e:.4f} rad = {np.degrees(delta_e):.2f} deg")
print(f"  A_mu   = {A_mu:.8f},  delta_mu = {delta_mu:.4f} rad = {np.degrees(delta_mu):.2f} deg")

r_21 = (A_mu/A_e)**4
print(f"  r_21   = {r_21:.4f}")


# ================================================================
# Section 1: Extract tails for a range of g0_tau candidates
# ================================================================
print("\n" + "=" * 70)
print("SECTION 1: TAIL PROPERTIES vs g0_tau")
print("=" * 70)

g0_tau_range = np.arange(1.42, 1.60, 0.005)
print(f"\n  {'g0_tau':>8s} {'A_tau':>10s} {'delta_tau':>10s} {'Q_K':>8s} "
      f"{'sum_A':>10s} {'sum_A2':>10s} {'cross':>10s} {'self/cross':>10s}")

for g0t in g0_tau_range:
    if g0t >= G0_CRIT:
        continue
    A_t, delta_t = extract_tail(g0t)
    if A_t < 1e-10:
        continue

    # x_i = A_i^2 (proportional to sqrt(m_i))
    x_e = A_e**2
    x_mu = A_mu**2
    x_tau = A_t**2

    S = x_e + x_mu + x_tau
    P = x_e**2 + x_mu**2 + x_tau**2
    cross = x_e*x_mu + x_e*x_tau + x_mu*x_tau
    QK = S**2 / P  # = (sum sqrt(m))^2 / sum(m) since m ~ x^2

    # Q_K = 3/2 <=> P = 4*cross  (self = 4*cross)
    ratio = P / cross if cross > 0 else 0

    marker = " <--" if abs(QK - 1.5) < 0.01 else ""
    print(f"  {g0t:8.4f} {A_t:10.6f} {delta_t:10.4f} {QK:8.4f} "
          f"{S:10.6f} {P:10.6f} {cross:10.6f} {ratio:10.4f}{marker}")


# ================================================================
# Section 2: Complex neutrality test
# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: COMPLEX NEUTRALITY  sum(A_i * exp(i*delta_i))")
print("=" * 70)

print(f"\n  Phases:")
print(f"    delta_e  = {delta_e:.6f} rad = {np.degrees(delta_e):.2f} deg")
print(f"    delta_mu = {delta_mu:.6f} rad = {np.degrees(delta_mu):.2f} deg")

# Find g0_tau for Koide
def QK_func(g0t):
    At, _ = extract_tail(g0t)
    xe, xm, xt = A_e**2, A_mu**2, At**2
    S = xe + xm + xt
    P = xe**2 + xm**2 + xt**2
    return S**2/P - 1.5

# Bracket search
g0_tau_koide = brentq(QK_func, 1.50, 1.59)
A_tau, delta_tau = extract_tail(g0_tau_koide)

print(f"    delta_tau = {delta_tau:.6f} rad = {np.degrees(delta_tau):.2f} deg")
print(f"    g0_tau(Koide) = {g0_tau_koide:.8f}")

# Complex sum
z_e = A_e * np.exp(1j * delta_e)
z_mu = A_mu * np.exp(1j * delta_mu)
z_tau = A_tau * np.exp(1j * delta_tau)

z_sum = z_e + z_mu + z_tau
print(f"\n  Complex sum: z_e + z_mu + z_tau")
print(f"    = {z_sum.real:.6f} + {z_sum.imag:.6f}i")
print(f"    |z_sum| = {abs(z_sum):.6f}")
print(f"    |z_sum| / |z_tau| = {abs(z_sum)/abs(z_tau):.6f}")

# Is it zero? If so, that's complex neutrality!
# Also check: what g0_tau would make it zero?
print(f"\n  Is complex neutrality achievable?")
print(f"  For z_sum = 0: need z_tau = -(z_e + z_mu)")
z_target = -(z_e + z_mu)
print(f"  Target: A_tau = {abs(z_target):.6f}, delta_tau = {np.degrees(np.angle(z_target)):.2f} deg")
print(f"  Actual: A_tau = {A_tau:.6f}, delta_tau = {np.degrees(delta_tau):.2f} deg")


# ================================================================
# Section 3: Field charge (integral of g-1)
# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: FIELD CHARGE Q_i = integral (g-1)*r^2 dr")
print("=" * 70)

Q_e = field_integral(g0_e)
Q_mu = field_integral(g0_mu)

print(f"  Q_e  = {Q_e:.6f}")
print(f"  Q_mu = {Q_mu:.6f}")

# Scan g0_tau
print(f"\n  {'g0_tau':>8s} {'Q_tau':>12s} {'Q_sum':>12s} {'Q_K':>8s}")
for g0t in np.arange(1.45, 1.59, 0.005):
    Q_t = field_integral(g0t)
    At, _ = extract_tail(g0t)
    xe, xm, xt = A_e**2, A_mu**2, At**2
    S = xe + xm + xt
    P = xe**2 + xm**2 + xt**2
    QK = S**2/P
    Q_sum = Q_e + Q_mu + Q_t
    marker = " <-- Koide" if abs(QK - 1.5) < 0.01 else ""
    print(f"  {g0t:8.4f} {Q_t:12.4f} {Q_sum:12.4f} {QK:8.4f}{marker}")

Q_tau = field_integral(g0_tau_koide)
print(f"\n  At Koide point:")
print(f"  Q_e + Q_mu + Q_tau = {Q_e + Q_mu + Q_tau:.6f}")
print(f"  Q_e   = {Q_e:.6f}")
print(f"  Q_mu  = {Q_mu:.6f}")
print(f"  Q_tau = {Q_tau:.6f}")

# Check ratios
print(f"  Q_mu/Q_e = {Q_mu/Q_e:.6f}")
print(f"  Q_tau/Q_e = {Q_tau/Q_e:.6f}")
print(f"  Q_tau/Q_mu = {Q_tau/Q_mu:.6f}")


# ================================================================
# Section 4: Field energy (integral of (g-1)^2)
# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: FIELD ENERGY E_i = integral (g-1)^2 * r^2 dr")
print("=" * 70)

E_e = field_integral_sq(g0_e)
E_mu = field_integral_sq(g0_mu)

print(f"  E_e  = {E_e:.6f}")
print(f"  E_mu = {E_mu:.6f}")

print(f"\n  {'g0_tau':>8s} {'E_tau':>12s} {'E_sum':>12s} {'E_rat':>10s} {'Q_K':>8s}")
for g0t in np.arange(1.45, 1.59, 0.005):
    E_t = field_integral_sq(g0t)
    At, _ = extract_tail(g0t)
    xe, xm, xt = A_e**2, A_mu**2, At**2
    S = xe + xm + xt
    P = xe**2 + xm**2 + xt**2
    QK = S**2/P
    E_sum = E_e + E_mu + E_t
    E_cross = E_e*E_mu + E_e*E_t + E_mu*E_t
    E_rat = E_sum**2 / (E_e**2 + E_mu**2 + E_t**2) if E_t > 0 else 0
    marker = " <-- Koide" if abs(QK - 1.5) < 0.01 else ""
    print(f"  {g0t:8.4f} {E_t:12.4f} {E_sum:12.4f} {E_rat:10.4f} {QK:8.4f}{marker}")

E_tau = field_integral_sq(g0_tau_koide)
print(f"\n  At Koide point:")
print(f"  E_e + E_mu + E_tau = {E_e + E_mu + E_tau:.6f}")
print(f"  E ratio = (sum E)^2 / sum(E^2) = {(E_e+E_mu+E_tau)**2/(E_e**2+E_mu**2+E_tau**2):.6f}")


# ================================================================
# Section 5: THE KEY TEST - "self = 4 * cross" in amplitude space
# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: SELF-ENERGY = 4 * CROSS-ENERGY CONDITION")
print("=" * 70)

print("""
  Q_K = 3/2 is EQUIVALENT to:
    sum(A_i^4) = 4 * sum_{i<j}(A_i^2 * A_j^2)

  i.e., "self-interaction = 4 * cross-interaction"
  in the space of tail amplitudes.

  Physical interpretation:
  Each soliton's tail creates a "field disturbance" with amplitude A_i.
  The "self-energy" of these disturbances is sum(A_i^4).
  The "cross-energy" (interference between pairs) is sum(A_i^2 * A_j^2).

  Q_K = 3/2 requires their ratio to be EXACTLY 4.

  WHY 4? Is this related to the structure of the TGP ODE?
""")

x_e, x_mu, x_tau = A_e**2, A_mu**2, A_tau**2
self_E = x_e**2 + x_mu**2 + x_tau**2
cross_E = x_e*x_mu + x_e*x_tau + x_mu*x_tau

print(f"  A_e^2   = {x_e:.8f}")
print(f"  A_mu^2  = {x_mu:.8f}")
print(f"  A_tau^2 = {x_tau:.8f}")
print(f"  self   = sum(A_i^4)         = {self_E:.8f}")
print(f"  cross  = sum(A_i^2*A_j^2)   = {cross_E:.8f}")
print(f"  self/cross = {self_E/cross_E:.6f}  (Koide requires: 4.0)")
print(f"  Q_K = (sum x)^2/sum(x^2) = {(x_e+x_mu+x_tau)**2/self_E:.6f}")


# ================================================================
# Section 6: Deeper - phase differences and tail structure
# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: PHASE STRUCTURE OF THREE SOLITONS")
print("=" * 70)

print(f"\n  Soliton tail: (g-1)*r ~ A_i * sin(r - delta_i)")
print(f"  Phase differences (mod 2pi):")
d12 = (delta_mu - delta_e) % (2*np.pi)
d13 = (delta_tau - delta_e) % (2*np.pi)
d23 = (delta_tau - delta_mu) % (2*np.pi)
print(f"    delta_mu - delta_e  = {d12:.4f} rad = {np.degrees(d12):.2f} deg")
print(f"    delta_tau - delta_e = {d13:.4f} rad = {np.degrees(d13):.2f} deg")
print(f"    delta_tau - delta_mu = {d23:.4f} rad = {np.degrees(d23):.2f} deg")
print(f"    Sum of phase diffs  = {d12+d23:.4f} rad = {np.degrees(d12+d23):.2f} deg")

# Are the phases related to 2pi/3 spacing?
print(f"\n  2*pi/3 = {2*np.pi/3:.4f} rad = 120.00 deg")
print(f"  pi/3   = {np.pi/3:.4f} rad = 60.00 deg")

# Check Koide circle: in Koide parametrization, the angles ARE 2pi/3 apart
# Is there a connection between the ODE phases and Koide angles?
# Koide: sqrt(m_i) = a(1 + b*cos(theta_K + 2pi(i-1)/3))
# ODE:   (g_i-1)*r ~ A_i*sin(r - delta_i)
# The Koide angles theta_K are in "mass space"
# The ODE phases delta_i are in "field space"
# Are they related?

theta_K = 2.316617  # from previous computation (for b=sqrt(2))
print(f"\n  Koide parametrization angle: theta_K = {theta_K:.4f} rad = {np.degrees(theta_K):.2f} deg")
print(f"  Koide angles: {np.degrees(theta_K):.2f}, {np.degrees(theta_K+2*np.pi/3):.2f}, {np.degrees(theta_K+4*np.pi/3):.2f}")
print(f"  ODE phases:   {np.degrees(delta_e):.2f}, {np.degrees(delta_mu):.2f}, {np.degrees(delta_tau):.2f}")


# ================================================================
# Section 7: TAIL OVERLAP INTEGRAL — the interference test
# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: TAIL OVERLAP INTEGRAL (interference)")
print("=" * 70)

# If solitons are at different positions r_i in real space,
# their tails overlap. The overlap integral between solitons i and j
# depends on A_i, A_j, delta_i, delta_j, and separation.
#
# But in TGP, the three generations are different solutions to the
# SAME ODE -- they don't coexist in the same physical space.
# Instead, they represent different "vacua" or field configurations.
#
# The "balance" is more abstract: it's in the PARAMETER space (g0).
# Each g0 defines a soliton with specific A and delta.
# The phi-FP rule constrains g0_mu = phi * g0_e.
# What constrains g0_tau?

# Let's compute the "tail wavefunction" overlap
# <psi_i | psi_j> where psi_i(r) = A_i * sin(r - delta_i) / r
# This integral over r from R to infinity:
# integral A_i*A_j * sin(r-d_i)*sin(r-d_j) / r^2 dr
# = (A_i*A_j/2) * integral [cos(d_i-d_j) - cos(2r-d_i-d_j)] / r^2 dr
# For large R, the oscillating part averages out:
# ~ (A_i*A_j/2) * cos(d_i-d_j) / R

print(f"\n  Asymptotic overlap integrals (proportional to):")
print(f"    <e|e>   ~ A_e^2 / 2 = {A_e**2/2:.8f}")
print(f"    <mu|mu> ~ A_mu^2/2  = {A_mu**2/2:.8f}")
print(f"    <tau|tau> ~ A_tau^2/2 = {A_tau**2/2:.8f}")
print(f"\n    <e|mu>  ~ A_e*A_mu/2 * cos(d12) = {A_e*A_mu/2*np.cos(d12):.8f}")
print(f"    <e|tau> ~ A_e*A_tau/2 * cos(d13) = {A_e*A_tau/2*np.cos(d13):.8f}")
print(f"    <mu|tau> ~ A_mu*A_tau/2* cos(d23) = {A_mu*A_tau/2*np.cos(d23):.8f}")

diag_sum = A_e**2 + A_mu**2 + A_tau**2
off_sum = A_e*A_mu*np.cos(d12) + A_e*A_tau*np.cos(d13) + A_mu*A_tau*np.cos(d23)
print(f"\n    Sum diagonal  = {diag_sum/2:.8f}")
print(f"    Sum off-diag  = {off_sum/2:.8f}")
print(f"    Total overlap = {(diag_sum + off_sum)/2:.8f}")
print(f"    diag/off ratio = {diag_sum/off_sum:.6f}")

# What if total overlap = 0?  diag + off = 0 means cos-weighted cross terms
# cancel the self terms.  That would be a real "neutrality" condition.


# ================================================================
# Section 8: "FIELD DEBT" INTERPRETATION
# ================================================================
print("\n" + "=" * 70)
print("SECTION 8: FIELD DEBT INTERPRETATION")
print("=" * 70)

print("""
  User's hypothesis: lighter solitons create a "field debt"
  that heavier solitons must repay.

  In TGP, the soliton modifies g from vacuum (g=1):
    - For g0 < 1 (electron): field is SUPPRESSED below vacuum
    - For g0 > 1 (muon, tau): field is ENHANCED above vacuum

  The "debt" picture: electron pulls g below 1 (creates deficit),
  muon and tau push g above 1 (create surplus).

  Balance: total field displacement = 0?
  Or: total field "action" = 0?
""")

# Sign of field displacement
print(f"  g0_e  = {g0_e:.4f} {'< 1 (deficit)' if g0_e < 1 else '> 1 (surplus)'}")
print(f"  g0_mu = {g0_mu:.4f} {'< 1 (deficit)' if g0_mu < 1 else '> 1 (surplus)'}")
print(f"  g0_tau= {g0_tau_koide:.4f} {'< 1 (deficit)' if g0_tau_koide < 1 else '> 1 (surplus)'}")

print(f"\n  Field charges (integral of (g-1)*r^2 dr):")
print(f"  Q_e   = {Q_e:.6f} {'(negative: deficit)' if Q_e < 0 else '(positive: surplus)'}")
print(f"  Q_mu  = {Q_mu:.6f} {'(negative: deficit)' if Q_mu < 0 else '(positive: surplus)'}")
print(f"  Q_tau = {Q_tau:.6f} {'(negative: deficit)' if Q_tau < 0 else '(positive: surplus)'}")

# What g0_tau would give Q_e + Q_mu + Q_tau = 0?
print(f"\n  For neutrality: Q_e + Q_mu + Q_tau = 0")
print(f"  Need Q_tau = -{Q_e + Q_mu:.6f}")

# Search for g0_tau that gives field neutrality
print(f"\n  Searching for g0_tau with sum(Q_i) = 0:")
q_target = -(Q_e + Q_mu)
print(f"  Target Q_tau = {q_target:.6f}")

for g0t in np.arange(1.42, 1.595, 0.005):
    Qt = field_integral(g0t)
    At, _ = extract_tail(g0t)
    xe, xm, xt = A_e**2, A_mu**2, At**2
    S = xe + xm + xt
    P = xe**2 + xm**2 + xt**2
    QK = S**2/P
    Q_sum = Q_e + Q_mu + Qt
    marker = ""
    if abs(QK - 1.5) < 0.01:
        marker = " <-- Koide"
    if abs(Q_sum) < abs(q_target)*0.1:
        marker += " <-- ~neutral"
    print(f"    g0_tau={g0t:.4f}: Q_tau={Qt:.4f}, sum={Q_sum:.4f}, Q_K={QK:.4f}{marker}")


# ================================================================
# Section 9: Energy budget -- total soliton energy
# ================================================================
print("\n" + "=" * 70)
print("SECTION 9: TOTAL SOLITON ENERGY vs Q_K")
print("=" * 70)

def soliton_energy(g0, r_max=300):
    """Compute the full soliton energy functional."""
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > 0.001) & (r < r_max)
    rf = r[mask]
    gf = g[mask]

    # Gradient of g
    gp = np.gradient(gf, rf)

    # Kinetic: g'^2 * r^2
    T = gp**2 * rf**2

    # Potential from g^2(1-g): integral of g^3/3 - g^4/4
    # V(g) = g^3/3 - g^4/4  (antiderivative of g^2(1-g))
    V = (gf**3/3.0 - gf**4/4.0) * rf**2

    T_int = trapz(T, rf)
    V_int = trapz(V, rf)
    return T_int, V_int, T_int + V_int

T_e, V_e, E_e_tot = soliton_energy(g0_e)
T_mu, V_mu, E_mu_tot = soliton_energy(g0_mu)

print(f"  Soliton energies (kinetic, potential, total):")
print(f"  e:   T={T_e:.4f}, V={V_e:.4f}, E={E_e_tot:.4f}")
print(f"  mu:  T={T_mu:.4f}, V={V_mu:.4f}, E={E_mu_tot:.4f}")

print(f"\n  {'g0_tau':>8s} {'E_tot':>10s} {'E_sum':>10s} {'E_ratio':>10s} {'Q_K':>8s}")
for g0t in np.arange(1.45, 1.595, 0.005):
    Tt, Vt, Et = soliton_energy(g0t)
    At, _ = extract_tail(g0t)
    xe, xm, xt = A_e**2, A_mu**2, At**2
    S = xe + xm + xt
    P = xe**2 + xm**2 + xt**2
    QK = S**2/P
    E_sum = E_e_tot + E_mu_tot + Et
    E_rat = Et / E_mu_tot if E_mu_tot != 0 else 0
    marker = " <-- Koide" if abs(QK - 1.5) < 0.01 else ""
    print(f"  {g0t:8.4f} {Et:10.4f} {E_sum:10.4f} {E_rat:10.4f} {QK:8.4f}{marker}")


# ================================================================
# FINAL ANALYSIS
# ================================================================
print("\n" + "=" * 70)
print("FINAL ANALYSIS: FIELD BALANCE HYPOTHESIS")
print("=" * 70)

print(f"""
  RESULTS:

  1. COMPLEX NEUTRALITY (sum A_i * exp(i*delta_i) = 0):
     |z_sum|/|z_tau| = {abs(z_sum)/abs(z_tau):.4f}
     {'CLOSE TO ZERO!' if abs(z_sum)/abs(z_tau) < 0.1 else 'NOT zero.'}

  2. FIELD CHARGE NEUTRALITY (sum Q_i = 0):
     Q_e + Q_mu + Q_tau = {Q_e + Q_mu + Q_tau:.4f}
     {'CLOSE TO ZERO!' if abs(Q_e + Q_mu + Q_tau) < abs(Q_tau)*0.1 else 'NOT zero at Koide point.'}

  3. SELF/CROSS RATIO:
     sum(A_i^4) / sum(A_i^2*A_j^2) = {self_E/cross_E:.4f}
     (Q_K = 3/2 requires exactly 4.0)
     {'CONFIRMED!' if abs(self_E/cross_E - 4.0) < 0.01 else 'CHECK VALUE'}

  4. PHASE STRUCTURE:
     Phase differences: {np.degrees(d12):.1f}, {np.degrees(d13):.1f}, {np.degrees(d23):.1f} deg
     (Koide circle has 120 deg spacing)

  KEY INSIGHT:
  Q_K = 3/2 means "self-interaction = 4 * cross-interaction" in
  the space of soliton tail amplitudes. This is a specific
  COUPLING CONDITION between the three generations.

  The question is: does the TGP ODE dynamics enforce this ratio?
""")
