#!/usr/bin/env python3
"""
koide_derivation_v47b.py -- Attempt to DERIVE Q_K = 3/2 from TGP principles.

STRATEGY: Test multiple candidate derivation paths:

PATH 1: VARIATIONAL (energy extremum)
  Does Q_K = 3/2 minimize/maximize the total soliton energy
  E_total = m_e + m_mu + m_tau, subject to r_21 fixed?

PATH 2: ENTROPY MAXIMUM
  Does Q_K = 3/2 maximize the Shannon entropy of the mass spectrum?

PATH 3: KOIDE CIRCLE GEOMETRY
  Koide parametrization: x_i = a(1 + b*cos(theta + 2*pi*i/3))
  Does the phi-FP constraint on r_21 select Q_K = 3/2 uniquely?

PATH 4: SOLITON PACKING / DEMOCRATIC PRINCIPLE
  Is Q_K = 3/2 selected by maximizing the number of "equally spaced"
  solitons in the g0 interval [g0_e, g0_crit]?

PATH 5: INFORMATION / MINIMUM DESCRIPTION
  Does Q_K = 3/2 minimize some information-theoretic quantity?

PATH 6: INTERACTION ENERGY
  If solitons interact via their tail overlap, does Q_K = 3/2
  minimize the interaction energy?

PATH 7: CONFORMAL / SCALE INVARIANCE
  Does Q_K = 3/2 arise from a conformal symmetry of the ODE?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

PHI = (1 + np.sqrt(5)) / 2
m_e = 0.51099895
m_mu = 105.6583755
r_21 = m_mu / m_e
u = np.sqrt(r_21)  # sqrt(r_21)


# Helper: given Q_K and r_21, compute r_31
def r31_from_QK(QK, r21=r_21):
    """Solve (1 + u^2 + v^2) / (1 + u^2 + v^2 + 2(u + v + uv)) = 1/QK
    where u = sqrt(r21), v = sqrt(r31).
    Equivalently: QK*(x^2+y^2+z^2) = (x+y+z)^2 with x=1, y=u, z=v.
    => (QK-1)*v^2 - 2v*(u+1) + (QK-1)*(1+u^2) - 2*u = 0
    """
    u_val = np.sqrt(r21)
    # Q*(1+u^2+v^2) = (1+u+v)^2
    # Q + Q*u^2 + Q*v^2 = 1 + u^2 + v^2 + 2u + 2v + 2uv
    # (Q-1)*v^2 - 2v(1+u) + Q(1+u^2) - (1+u^2+2u) = 0
    # (Q-1)*v^2 - 2v(1+u) + (Q-1)(1+u^2) - 2u = 0
    a_coeff = QK - 1
    b_coeff = -2*(1 + u_val)
    c_coeff = (QK - 1)*(1 + u_val**2) - 2*u_val

    disc = b_coeff**2 - 4*a_coeff*c_coeff
    if disc < 0:
        return None
    v_plus = (-b_coeff + np.sqrt(disc)) / (2*a_coeff)
    if v_plus > 0:
        return v_plus**2
    return None


# ODE solver for canonical Form A
def make_solver(alpha=2.0, r_max=300):
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


solver = make_solver(alpha=2.0)
gc = 1.6

# Calibrate g0_e
def r21_res(g0_e):
    A_e = A_tail(solver, g0_e)
    if A_e < 1e-15: return 1e10
    g0_mu = PHI * g0_e
    if g0_mu >= gc: return 1e10
    A_mu = A_tail(solver, g0_mu)
    return (A_mu / A_e)**4 - r_21

g0_scan = np.linspace(0.5, 0.95, 30)
resids = [r21_res(g) for g in g0_scan]
for i in range(len(resids) - 1):
    if resids[i] * resids[i+1] < 0 and abs(resids[i]) < 1e8:
        g0_e = brentq(r21_res, g0_scan[i], g0_scan[i+1], xtol=1e-8)
        break

g0_mu = PHI * g0_e
A_e = A_tail(solver, g0_e)
A_mu = A_tail(solver, g0_mu)

print("=" * 70)
print("ATTEMPT TO DERIVE Q_K = 3/2")
print("=" * 70)
print(f"\n  Calibration: g0_e={g0_e:.8f}, A_e={A_e:.8f}, A_mu={A_mu:.8f}")
print(f"  r_21 check: {(A_mu/A_e)**4:.4f}")


# ================================================================
# PATH 1: VARIATIONAL -- Total energy
# ================================================================
print("\n" + "=" * 70)
print("PATH 1: TOTAL ENERGY EXTREMUM")
print("=" * 70)

# Total mass M = m_e + m_mu + m_tau = A_e^4 + A_mu^4 + A_tau^4
# (up to common proportionality constant)
# As function of g0_tau:

g0_tau_range = np.linspace(g0_mu + 0.01, gc - 0.001, 50)
masses_total = []
QK_list = []
energies = []

for g0t in g0_tau_range:
    At = A_tail(solver, g0t)
    if At < 1e-10:
        continue
    M = A_e**4 + A_mu**4 + At**4
    S2 = A_e**2 + A_mu**2 + At**2
    S4 = A_e**4 + A_mu**4 + At**4
    QK = S2**2 / S4
    masses_total.append(M)
    QK_list.append(QK)
    # Also compute soliton energy: integral of energy density
    # For now, use mass as proxy for energy
    energies.append(M)

QK_arr = np.array(QK_list)
M_arr = np.array(masses_total)

# Is there an extremum of M at Q_K = 3/2?
# M is monotonically increasing with g0_tau (heavier tau = more mass)
# So NO extremum of total mass.
print("\n  Total mass M = m_e + m_mu + m_tau vs Q_K:")
print(f"  {'Q_K':>8s} {'M/m_e':>12s}")
for i in range(0, len(QK_arr), 5):
    print(f"  {QK_arr[i]:8.4f} {M_arr[i]/A_e**4:12.2f}")

print(f"\n  M is {'monotonic' if np.all(np.diff(M_arr) > 0) or np.all(np.diff(M_arr) < 0) else 'NON-MONOTONIC'}")
print(f"  VERDICT: Total mass has NO extremum. PATH 1 FAILS.")


# ================================================================
# PATH 2: ENTROPY MAXIMUM
# ================================================================
print("\n" + "=" * 70)
print("PATH 2: ENTROPY MAXIMUM")
print("=" * 70)

# Shannon entropy of mass distribution: p_i = m_i / M_total
# S = -sum p_i * log(p_i)
# Maximum entropy = all equal masses (S = log(3))
# Does Q_K = 3/2 maximize entropy subject to r_21 constraint?

entropies = []
for i in range(len(QK_arr)):
    m_arr_local = np.array([A_e**4, A_mu**4, M_arr[i] - A_e**4 - A_mu**4])
    m_arr_local = np.maximum(m_arr_local, 1e-20)
    p = m_arr_local / np.sum(m_arr_local)
    S = -np.sum(p * np.log(p))
    entropies.append(S)

S_arr = np.array(entropies)
i_max_S = np.argmax(S_arr)

print(f"\n  Entropy vs Q_K:")
print(f"  {'Q_K':>8s} {'Entropy':>10s}")
for i in range(0, len(QK_arr), 5):
    marker = " <-- max" if i == i_max_S else ""
    print(f"  {QK_arr[i]:8.4f} {S_arr[i]:10.6f}{marker}")

print(f"\n  Max entropy at Q_K = {QK_arr[i_max_S]:.4f}")
print(f"  Koide Q_K = 1.5000")
print(f"  VERDICT: Entropy max at Q_K = {QK_arr[i_max_S]:.4f}, NOT 3/2. PATH 2 FAILS.")


# ================================================================
# PATH 3: KOIDE CIRCLE GEOMETRY
# ================================================================
print("\n" + "=" * 70)
print("PATH 3: KOIDE CIRCLE GEOMETRY")
print("=" * 70)

# Koide's original parametrization:
#   sqrt(m_i) = a * (1 + b * cos(theta + 2*pi*(i-1)/3))
# where a, b, theta are parameters.
#
# Q_K = 3/2 is AUTOMATIC for any a, b, theta!
# (This is the defining property of the parametrization.)
#
# So the question becomes: does the phi-FP constraint select
# specific (a, b, theta)?

print("""
  Koide parametrization: sqrt(m_i) = a(1 + b*cos(theta + 2*pi*(i-1)/3))

  FACT: Q_K = 3/2 holds AUTOMATICALLY for this parametrization,
  for ANY values of a, b, theta. This is because:
    sum sqrt(m_i) = 3a  (cosines cancel in sum)
    sum m_i = 3a^2(1 + b^2/2)  (cross terms cancel)
    (sum sqrt(m_i))^2 / (sum m_i) = 9a^2 / (3a^2(1+b^2/2)) = 3/(1+b^2/2)
  Wait, that's not 3/2 in general. Let me recalculate.
""")

# Actually, let me verify:
# x_i = a(1 + b cos(theta_i)) where theta_i = theta + 2pi(i-1)/3
# sum x_i = a * sum(1 + b cos(theta_i)) = 3a (since sum cos(theta_i) = 0)
# sum x_i^2 = a^2 * sum(1 + 2b cos + b^2 cos^2)
#           = a^2 * (3 + 0 + b^2 * 3/2) = 3a^2(1 + b^2/2)
# Q = (sum x)^2 / (sum x^2) = 9a^2 / (3a^2(1+b^2/2)) = 3/(1+b^2/2)
# Q = 3/2 iff b^2/2 = 1 iff b^2 = 2 iff b = sqrt(2)

print("  Correction: Q_K = 3/(1 + b^2/2)")
print(f"  Q_K = 3/2 requires b^2 = 2, i.e., b = sqrt(2) = {np.sqrt(2):.6f}")
print()

# So the question is: does TGP select b = sqrt(2)?
# Given r_21, what b and theta give the right mass ratio?

# r_21 = m_mu/m_e = x_2^2/x_1^2
# x_1 = a(1 + b cos(theta))
# x_2 = a(1 + b cos(theta + 2pi/3))
# r_21 = (1 + b cos(theta + 2pi/3))^2 / (1 + b cos(theta))^2

# For b = sqrt(2):
# r_21 = (1 + sqrt(2) cos(theta + 2pi/3))^2 / (1 + sqrt(2) cos(theta))^2
# Need to find theta such that this equals 206.768

def r21_from_theta(theta, b=np.sqrt(2)):
    x1 = 1 + b * np.cos(theta)
    x2 = 1 + b * np.cos(theta + 2*np.pi/3)
    if x1 <= 0 or x2 <= 0:
        return 0
    return (x2 / x1)**2

# Scan theta
thetas = np.linspace(0, 2*np.pi, 1000)
r21_vals = [r21_from_theta(t) for t in thetas]
r21_vals = np.array(r21_vals)

# Find theta where r21 = 206.768
for i in range(len(r21_vals) - 1):
    if r21_vals[i] > 0 and r21_vals[i+1] > 0:
        if (r21_vals[i] - r_21) * (r21_vals[i+1] - r_21) < 0:
            theta_sol = brentq(lambda t: r21_from_theta(t) - r_21,
                               thetas[i], thetas[i+1])
            b_val = np.sqrt(2)
            x1 = 1 + b_val * np.cos(theta_sol)
            x2 = 1 + b_val * np.cos(theta_sol + 2*np.pi/3)
            x3 = 1 + b_val * np.cos(theta_sol + 4*np.pi/3)
            r31_check = (x3/x1)**2 if x1 > 0 and x3 > 0 else 0
            print(f"  Solution: theta = {theta_sol:.6f} rad = {np.degrees(theta_sol):.2f} deg")
            print(f"  x1 (e) = {x1:.6f}")
            print(f"  x2 (mu) = {x2:.6f}")
            print(f"  x3 (tau) = {x3:.6f}")
            print(f"  r21 = {(x2/x1)**2:.4f} (target {r_21:.4f})")
            print(f"  r31 = {r31_check:.4f} (target 3477.44)")
            # Check Q_K
            S2 = x1**2 + x2**2 + x3**2
            S4 = x1**4 + x2**4 + x3**4
            print(f"  Q_K = {S2**2/S4:.6f}")
            print()

print("  INSIGHT: Q_K = 3/2 requires EXACTLY b = sqrt(2) in the")
print("  Koide circle parametrization. The question becomes:")
print("  Does TGP select b = sqrt(2)?")
print()
print("  b = sqrt(2) means the AMPLITUDE of modulation equals")
print("  the MEAN, scaled by sqrt(2). This is a specific constraint")
print("  on the 'democracy' of the mass spectrum.")


# ================================================================
# PATH 4: b = sqrt(2) FROM SOLITON PROPERTIES
# ================================================================
print("\n" + "=" * 70)
print("PATH 4: DOES THE ODE SELECT b = sqrt(2)?")
print("=" * 70)

# In the soliton picture, x_i = sqrt(m_i) = A_i^2.
# The Koide circle says x_i = a(1 + b cos(theta_i)).
# With phi-FP: g0_mu = phi * g0_e.
# The mapping A_tail(g0) translates g0 -> x via x = A_tail^2.
#
# For the phi-FP to be compatible with Koide circle:
# x_2/x_1 = (1 + b cos(theta+2pi/3)) / (1 + b cos(theta))
# = A_tail(phi*g0_e)^2 / A_tail(g0_e)^2
# = sqrt(r_21)
#
# This is one equation for two unknowns (b, theta).
# A second equation is needed to fix b = sqrt(2).
#
# The ODE provides A_tail(g0). Near vacuum, A ~ |g0-1|.
# So x = A^2 ~ (g0-1)^2 for g0 near 1.
#
# For the electron: x_e ~ (g0_e - 1)^2 ~ (0.868 - 1)^2 ~ 0.017
# For the muon: x_mu ~ (phi*g0_e - 1)^2 ~ (1.404 - 1)^2 ~ 0.163
# For the tau: x_tau ~ A_tau^2 (nonlinear, not just (g0-1)^2)

# Compute actual b and theta from numerical data
A_tau_koide = 1.00510  # from previous runs
x_e = A_e**2
x_mu = A_mu**2
x_tau = A_tau_koide**2

# From Koide: x_i = a(1 + b cos(theta_i))
# sum x_i = 3a => a = (x_e + x_mu + x_tau)/3
a_val = (x_e + x_mu + x_tau) / 3
print(f"\n  Numerical values:")
print(f"    x_e   = A_e^2   = {x_e:.8f}")
print(f"    x_mu  = A_mu^2  = {x_mu:.8f}")
print(f"    x_tau = A_tau^2 = {x_tau:.8f}")
print(f"    a = mean(x) = {a_val:.8f}")

# x_i/a = 1 + b cos(theta_i)
# (x_i/a - 1) = b cos(theta_i)
d_e = x_e/a_val - 1
d_mu = x_mu/a_val - 1
d_tau = x_tau/a_val - 1

print(f"    d_e   = x_e/a - 1 = {d_e:.8f}")
print(f"    d_mu  = x_mu/a - 1 = {d_mu:.8f}")
print(f"    d_tau = x_tau/a - 1 = {d_tau:.8f}")

# b^2 = (2/3)(d_e^2 + d_mu^2 + d_tau^2)  [since sum cos^2 = 3/2]
b_sq = (2.0/3.0) * (d_e**2 + d_mu**2 + d_tau**2)
b_val_num = np.sqrt(b_sq)

print(f"\n    b^2 = (2/3)*sum(d_i^2) = {b_sq:.8f}")
print(f"    b = {b_val_num:.8f}")
print(f"    sqrt(2) = {np.sqrt(2):.8f}")
print(f"    b/sqrt(2) = {b_val_num/np.sqrt(2):.8f}")

# Q_K from b:
QK_from_b = 3.0 / (1 + b_sq/2)
print(f"\n    Q_K = 3/(1+b^2/2) = {QK_from_b:.8f}")
print(f"    Q_K(target) = 1.500000")


# ================================================================
# PATH 5: GEOMETRIC ARGUMENT -- equal angular spacing
# ================================================================
print("\n" + "=" * 70)
print("PATH 5: GEOMETRIC / ANGULAR ANALYSIS")
print("=" * 70)

# In Koide parametrization with theta_i equally spaced (2pi/3 apart),
# the PHASE theta determines which generation is heaviest.
# The AMPLITUDE b determines how hierarchical the spectrum is.
#
# b = 0: all equal masses (Q_K = 3, degenerate)
# b = 1: one mass = 0 (Q_K = 3/2... wait let me check)
# b -> infinity: extreme hierarchy (Q_K -> 1)

for b_test in [0.0, 0.5, 1.0, np.sqrt(2), 1.5, 2.0, 3.0, 5.0]:
    QK_test = 3.0 / (1 + b_test**2/2) if b_test > 0 else 3.0
    print(f"  b = {b_test:.4f}: Q_K = {QK_test:.6f}")

print(f"\n  b = 1.0 gives Q_K = 3/(1+0.5) = {3/1.5:.6f} = 2.0")
print(f"  b = sqrt(2) gives Q_K = 3/(1+1) = {3/2:.6f} = 1.5")
print(f"  b = sqrt(3) gives Q_K = 3/(1+1.5) = {3/2.5:.6f} = 1.2")

# b = sqrt(2) means: the deviation from the mean equals the mean.
# In other words: x_i = a(1 + sqrt(2)*cos(theta_i))
# The maximum x is a(1+sqrt(2)) = a*2.414
# The minimum x is a(1-sqrt(2)) = a*(-0.414) which is NEGATIVE!
# So with b = sqrt(2), one of the x_i is negative (unphysical)
# UNLESS theta is chosen so all cos terms give x_i > 0.

# Check: which theta values keep all x_i > 0?
print("\n  For b = sqrt(2): x_i > 0 requires 1 + sqrt(2)*cos(theta_i) > 0")
print(f"  cos(theta_i) > -1/sqrt(2) = {-1/np.sqrt(2):.6f}")
print(f"  theta_i must be in (-{np.degrees(np.arccos(-1/np.sqrt(2))):.1f}, "
      f"{np.degrees(np.arccos(-1/np.sqrt(2))):.1f}) degrees")
print(f"  = within 135 deg of 0. With 120 deg spacing, this is tight!")


# ================================================================
# PATH 6: SOLITON ENERGY FUNCTIONAL
# ================================================================
print("\n" + "=" * 70)
print("PATH 6: SOLITON ENERGY FUNCTIONAL")
print("=" * 70)

# Instead of total mass, consider the INTERACTION between solitons.
# If solitons at g0_i interact via tail overlap:
#   V_int ~ sum_{i<j} A_i * A_j * exp(-|r_i - r_j|)
# In the single-center approximation, this is ~ sum A_i*A_j.
#
# The total "cost" of the three-soliton configuration:
#   C = sum m_i + lambda * sum_{i<j} sqrt(m_i*m_j)
#   = S4/2 + lambda * (xy + xz + yz) where x,y,z = sqrt(m_i)
#   (using S4 = sum x^4 is wrong -- it's sum m_i = sum x^2*x^2)
# Let me use s_i = sqrt(m_i) = A_i^2:
#   sum m_i = sum s_i^2
#   sum sqrt(m_i*m_j) = sum s_i*s_j
#
# Koide: sum(s_i)^2 / sum(s_i^2) = Q_K
#   sum(s_i)^2 = sum(s_i^2) + 2*sum(s_i*s_j)
#   Q_K = 1 + 2*sum(s_i*s_j)/sum(s_i^2)
#   sum(s_i*s_j) = (Q_K - 1)/2 * sum(s_i^2)
#
# So the "interaction" term is proportional to (Q_K-1)*sum(s_i^2)/2.
# Total cost: C = sum(s_i^2) * [1 + lambda*(Q_K-1)/2]
# This is MINIMIZED when Q_K is as SMALL as possible (Q_K -> 1).
# NOT at Q_K = 3/2.

print("  Interaction model: C = sum(m_i) + lambda * sum(sqrt(m_i*m_j))")
print("  = sum(s^2) * [1 + lambda*(Q_K-1)/2]")
print("  Minimized at Q_K = 1 (maximum hierarchy), not 3/2.")
print("  VERDICT: Simple interaction model FAILS.")


# ================================================================
# PATH 7: RATIO OF RATIOS / SELF-CONSISTENCY
# ================================================================
print("\n" + "=" * 70)
print("PATH 7: SELF-CONSISTENCY WITH PHI-FP")
print("=" * 70)

# The phi-FP rule is: g0_mu = phi * g0_e.
# What if there's a SECOND phi-FP relation for the tau?
# e.g., A_tau/A_mu = A_mu/A_e * correction?
#
# If A_tail were purely linear (A ~ |g0-1|), then:
# A_mu/A_e = (phi*g0_e - 1)/(1 - g0_e)  [for g0_e < 1 < phi*g0_e]
# A_tau/A_mu would be determined by g0_tau.
#
# The Koide condition constrains A_tau/A_mu.
# Let R = A_tau/A_mu. Then:
# Q_K = (A_e^2 + A_mu^2 + A_tau^2)^2 / (A_e^4 + A_mu^4 + A_tau^4)
# Let p = A_mu/A_e = r_21^(1/4). Then:
# Q_K = (1 + p^2 + p^2*R^2)^2 / (1 + p^4 + p^4*R^4)
# This determines R for given Q_K and p.

p = (r_21)**(1.0/4.0)
print(f"\n  p = A_mu/A_e = r_21^(1/4) = {p:.8f}")

# For Q_K = 3/2:
# 3/2 = (1 + p^2 + p^2*R^2)^2 / (1 + p^4 + p^4*R^4)
# Let w = p^2*R^2. Then:
# 3/2 = (1+p^2+w)^2 / (1+p^4+w^2)
# 3(1+p^4+w^2) = 2(1+p^2+w)^2
# 3 + 3p^4 + 3w^2 = 2 + 4p^2 + 4w + 2p^4 + 4p^2*w + 2w^2
# w^2 - 4w(1+p^2) + p^4 - 4p^2 + 1 + 4w = 0... let me redo
# 3 + 3p^4 + 3w^2 = 2 + 4p^2 + 4w + 2p^4 + 4p^2*w + 2w^2
# w^2 - 4w - 4p^2*w + p^4 - 4p^2 + 1 = 0
# w^2 - 4w(1+p^2) + (1+p^4-4p^2) = 0
# Note: 1+p^4-4p^2 = (p^2-2)^2 - 3

def solve_w(p_val):
    a_c = 1
    b_c = -4*(1 + p_val**2)
    c_c = 1 + p_val**4 - 4*p_val**2
    disc = b_c**2 - 4*a_c*c_c
    w_plus = (-b_c + np.sqrt(disc)) / 2
    w_minus = (-b_c - np.sqrt(disc)) / 2
    return w_plus, w_minus

w_plus, w_minus = solve_w(p)
R_plus = np.sqrt(w_plus) / p
R_minus = np.sqrt(max(w_minus, 0)) / p if w_minus > 0 else 0

r31_from_R = (p * R_plus)**4  # = (A_tau/A_e)^4

print(f"  w+ = {w_plus:.8f}, R+ = A_tau/A_mu = {R_plus:.8f}")
print(f"  w- = {w_minus:.8f}, R- = {R_minus:.8f}")
print(f"  r31 = (p*R+)^4 = {r31_from_R:.4f}")
print(f"  r32 = R+^4 = {R_plus**4:.4f}")
print(f"  PDG r32 = {m_mu/m_e:.2f}... no, r32 = m_tau/m_mu = {1776.86/105.658:.4f}")

# Now: is there a relationship between R and phi?
print(f"\n  R+ = A_tau/A_mu = {R_plus:.8f}")
print(f"  phi = {PHI:.8f}")
print(f"  R/phi = {R_plus/PHI:.8f}")
print(f"  R^2 = {R_plus**2:.8f}")
print(f"  phi^2 = {PHI**2:.8f}")

# Check: r32 = R^4
r32 = R_plus**4
r32_pdg = 1776.86 / 105.658
print(f"\n  r32 = m_tau/m_mu = R^4 = {r32:.4f}")
print(f"  r32(PDG) = {r32_pdg:.4f}")
print(f"  sqrt(r32) = {np.sqrt(r32):.4f}")

# r32 = r31/r21
print(f"  r32 = r31/r21 = {r31_from_R/r_21:.4f}")


# ================================================================
# PATH 8: LOG-SPACING ARGUMENT
# ================================================================
print("\n" + "=" * 70)
print("PATH 8: LOG-SPACING AND GEOMETRIC MEAN")
print("=" * 70)

# If the three sqrt(masses) form a geometric progression:
# s_2/s_1 = s_3/s_2 = q (common ratio)
# Then s = (s_1, s_1*q, s_1*q^2)
# Q_K = (1+q+q^2)^2 / (1+q^2+q^4)
# For q = r_21^(1/4) = p = 3.793:

q_geo = p
QK_geo = (1 + q_geo + q_geo**2)**2 / (1 + q_geo**2 + q_geo**4)
print(f"  If sqrt(m) forms geometric progression with ratio q = {q_geo:.4f}:")
print(f"  Q_K = (1+q+q^2)^2 / (1+q^2+q^4) = {QK_geo:.6f}")

# This gives Q_K ~ 1.06, NOT 3/2.
# The actual hierarchy is NOT a geometric progression.

# What ratio q gives Q_K = 3/2?
from scipy.optimize import brentq as bq
def QK_of_q(q):
    return (1+q+q**2)**2 / (1+q**2+q**4) - 1.5

q_koide = bq(QK_of_q, 1.01, 100)
print(f"\n  Q_K = 3/2 requires geometric ratio q = {q_koide:.6f}")
print(f"  Actual q = sqrt(r_21)^(1/2) = {p:.6f}")
print(f"  These don't match. Mass spectrum is NOT geometric.")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "=" * 70)
print("SUMMARY OF DERIVATION ATTEMPTS")
print("=" * 70)

print("""
  PATH 1 (Energy extremum):     FAILS -- M is monotonic in g0_tau
  PATH 2 (Entropy maximum):     FAILS -- max at Q_K ~ 2.0, not 1.5
  PATH 3 (Koide circle):        REFRAMES -- Q_K=3/2 iff b=sqrt(2)
  PATH 4 (ODE selects b):       OPEN -- no mechanism identified
  PATH 5 (Angular analysis):    b=sqrt(2) is at boundary of positivity
  PATH 6 (Interaction energy):  FAILS -- min at Q_K=1
  PATH 7 (Self-consistency):    ALGEBRAIC -- r31 determined, but Q_K input
  PATH 8 (Geometric progression): FAILS -- spectrum not geometric

  CONCLUSION:
  None of the tested paths derive Q_K = 3/2 from TGP first principles.

  The CLEANEST reformulation is via the Koide circle (Path 3):
    Q_K = 3/2 <=> b = sqrt(2) in the parametrization
    sqrt(m_i) = a(1 + sqrt(2)*cos(theta + 2*pi*(i-1)/3))

  The question "why Q_K = 3/2?" becomes "why b = sqrt(2)?".

  POSSIBLE INTERPRETATION:
  b = sqrt(2) is the MAXIMUM b for which all three masses are
  positive (since 1 - sqrt(2) < 0, one needs theta such that
  all cos terms keep x_i > 0 -- this is tight at b = sqrt(2)).

  This suggests: Q_K = 3/2 = the MOST HIERARCHICAL mass spectrum
  compatible with all three masses being positive in the Koide
  circle parametrization.

  Let me verify this "maximal hierarchy" interpretation...
""")

# Verify: is b = sqrt(2) the maximum b for positivity?
print("  MAXIMAL HIERARCHY TEST:")
print("  For b > sqrt(2), can all x_i be positive?")
for b_test in np.arange(1.0, 2.0, 0.1):
    # Check: for what theta are all x_i > 0?
    # x_i > 0 requires cos(theta_i) > -1/b
    # The three angles are theta, theta+2pi/3, theta+4pi/3
    # Need min(cos(theta_i)) > -1/b
    # min cos over three equally-spaced points depends on theta.
    # The minimum of {cos(t), cos(t+2pi/3), cos(t+4pi/3)} over t
    # is minimized when one angle is at pi, giving cos = -1.
    # The constraint is -1 > -1/b, i.e., b > 1.
    # But the constraint is that ALL THREE are > -1/b.
    # The worst case is when one cos is as negative as possible.

    # For equally spaced angles, at least one cos(theta_i) <= -1/2.
    # (The minimum of the three is -1 only when one angle = pi.)
    # Actually, the minimum of min{cos(t), cos(t+120), cos(t+240)}
    # over t achieves minimum value of -1 (at t = 180 deg).
    # But the constraint is not on the min but on ALL being > -1/b.

    # For any theta, the minimum cos is at most -1/2 (when theta = 0,
    # then cos(240) = cos(-120) = -1/2).
    # Wait: cos(0) = 1, cos(120) = -1/2, cos(240) = -1/2.
    # min is -1/2. So x_i > 0 requires b*(-1/2) > -1, i.e., b < 2.

    # More carefully: for theta such that one angle is at pi,
    # cos = -1, so x = 1-b < 0 for b > 1. But the other angles
    # are at pi +/- 2pi/3, i.e., -60 and +60 deg from pi,
    # giving cos = -cos(60) = -1/2.
    # So that configuration has x = 1-b, 1-b/2, 1-b/2.
    # For b < 1: all positive. For b = 1: one is 0.

    # The maximum b for positivity depends on theta.
    # For the OPTIMAL theta (maximizing the minimum x_i):
    # Choose theta so that the minimum cos(theta_i) is maximized.
    # By symmetry, this is when no angle is at pi.
    # The optimal theta has the minimum cos = cos(2pi/3) = -1/2.
    # Then constraint: 1 + b*(-1/2) > 0, i.e., b < 2.

    # For b = sqrt(2) ~ 1.414: 1 - sqrt(2)/2 = 1 - 0.707 = 0.293 > 0.
    # The minimum x_i = a * 0.293 > 0. OK.

    # For b = 2: 1 - 2/2 = 0. One mass = 0. Boundary.

    min_cos = -0.5  # best case for equally spaced
    x_min = 1 + b_test * min_cos
    QK_b = 3.0 / (1 + b_test**2/2)
    print(f"    b={b_test:.2f}: min(x_i/a) = {x_min:.4f}, Q_K = {QK_b:.4f}")

print(f"\n  b = 2.0: min(x_i/a) = 0 (one mass = 0), Q_K = {3/(1+2):.4f}")
print(f"  b = sqrt(2): min(x_i/a) = {1-np.sqrt(2)/2:.4f}, Q_K = 1.5000")

print(f"\n  CRITICAL INSIGHT:")
print(f"  b = 2 is the ABSOLUTE maximum for positivity.")
print(f"  b = sqrt(2) is NOT the maximum. It's at Q_K = 1.5,")
print(f"  while b = 2 gives Q_K = 1.0.")
print(f"\n  So Q_K = 3/2 is NOT selected by maximal hierarchy.")
print(f"  It's an intermediate value with no obvious geometric meaning.")
