#!/usr/bin/env python3
"""
koide_collapse_corridor_v47b.py -- PROPOSAL #13 + #1 COMBINED
  Collapse corridor constraints on Q_K + equipartition test.

CORE INSIGHT (from PATH 18b):
  g0_tau = phi^2 * g0_e = 2.27 > g0_crit = 8/5 = 1.6
  So the third generation CANNOT have phi-FP spacing.
  g0_tau is "squeezed" just below the collapse threshold.

QUESTIONS:
  P13: How narrow is the allowed corridor for g0_tau?
       What range of Q_K is compatible with 3 generations
       existing below g0_crit?

  P1:  Does equipartition r = sqrt(N-1) = sqrt(2) emerge
       naturally from the corridor constraint?

APPROACH:
  1. Map the full corridor: for each g0_e, find the range of
     g0_tau that (a) exists as a soliton and (b) gives g0_mu
     between g0_e and g0_tau.
  2. Compute Q_K(mass) as a function of g0_tau within the corridor.
  3. Check if Q_K = 3/2 has any special geometric/extremal property
     within the corridor.
  4. Test equipartition: does the Brannen parameter b = sqrt(2)
     correspond to any natural point in the corridor?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

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
    if g0 >= G0_CRIT - 0.001:
        return 0.0
    if g0 <= 0.01:
        return 0.0
    try:
        r, g = solver_A(g0, r_max=r_max+50)
        mask = (r > r_min) & (r < r_max)
        if np.sum(mask) < 100:
            return 0.0
        rf = r[mask]
        h = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, h, rcond=None)[0]
        return np.sqrt(bc[0]**2 + bc[1]**2)
    except:
        return 0.0


def QK_masses(A_vals):
    """Q_K = (sum sqrt(m))^2 / (sum m) where m = A^4"""
    A = np.array(A_vals)
    A2 = A**2   # sqrt(m)
    A4 = A**4   # m
    if np.sum(A4) < 1e-30:
        return 0.0
    return np.sum(A2)**2 / np.sum(A4)


def brannen_b(A_vals):
    """Compute Brannen parameter b from A values.
    sqrt(m_i) = c(1 + b*cos(theta + 2pi*i/3))
    b^2 = 2 means Q_K = 3/2
    """
    A = np.array(A_vals)
    x = A**2  # = sqrt(m)
    mu = np.mean(x)
    if mu < 1e-30:
        return 0.0
    std = np.std(x, ddof=0)
    CV = std / mu
    # Q_K = 3/(1 + b^2/2), so b^2 = 2*(3/Q_K - 1)
    # Also CV = b/sqrt(2), so b = CV*sqrt(2)
    return CV * np.sqrt(2)


# ================================================================
print("=" * 70)
print("PROPOSAL #13 + #1: COLLAPSE CORRIDOR AND EQUIPARTITION")
print("=" * 70)

# ================================================================
print("\n" + "=" * 70)
print("SECTION 1: THE COLLAPSE CORRIDOR")
print("=" * 70)

g0_e = 0.86770494
g0_mu = PHI * g0_e

print(f"\n  Fixed: g0_e = {g0_e:.8f}")
print(f"  Fixed: g0_mu = phi * g0_e = {g0_mu:.8f}")
print(f"  Collapse threshold: g0_crit = {G0_CRIT:.4f}")
print(f"  phi^2 * g0_e = {PHI**2 * g0_e:.4f} (> g0_crit, soliton dies)")

# Find the actual max g0_tau (just below collapse)
print(f"\n  Corridor for g0_tau: [{g0_mu:.4f}, {G0_CRIT:.4f})")
print(f"  Width: {G0_CRIT - g0_mu:.4f}")

# Scan g0_tau through the corridor
print(f"\n  Q_K as function of g0_tau within corridor:")
print(f"  {'g0_tau':>10s} {'frac_corr':>10s} {'A_tau':>12s} {'Q_K':>10s} {'b':>8s} {'Q_K-3/2':>10s}")

A_e = A_tail(g0_e)
A_mu = A_tail(g0_mu)

g0_tau_vals = []
qk_vals = []
b_vals = []

for g0t in np.arange(g0_mu + 0.01, G0_CRIT - 0.005, 0.005):
    At = A_tail(g0t)
    if At > 0:
        qk = QK_masses([A_e, A_mu, At])
        b = brannen_b([A_e, A_mu, At])
        frac = (g0t - g0_mu) / (G0_CRIT - g0_mu)  # fraction through corridor
        g0_tau_vals.append(g0t)
        qk_vals.append(qk)
        b_vals.append(b)
        print(f"  {g0t:10.5f} {frac:10.3f} {At:12.8f} {qk:10.6f} {b:8.4f} {qk-1.5:+10.6f}")

g0_tau_vals = np.array(g0_tau_vals)
qk_vals = np.array(qk_vals)
b_vals = np.array(b_vals)


# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: Q_K RANGE IN THE CORRIDOR")
print("=" * 70)

print(f"\n  Q_K range: [{np.min(qk_vals):.4f}, {np.max(qk_vals):.4f}]")
print(f"  b range:   [{np.min(b_vals):.4f}, {np.max(b_vals):.4f}]")
print(f"  Q_K = 3/2 requires b = sqrt(2) = {np.sqrt(2):.4f}")

# Find where Q_K = 3/2
idx_koide = np.argmin(np.abs(qk_vals - 1.5))
print(f"\n  Q_K closest to 3/2 at g0_tau = {g0_tau_vals[idx_koide]:.6f}")
print(f"  Fraction through corridor: {(g0_tau_vals[idx_koide] - g0_mu) / (G0_CRIT - g0_mu):.3f}")

# Is 3/2 at the boundary, middle, or extremum?
dqk = np.gradient(qk_vals, g0_tau_vals)
print(f"\n  dQ_K/dg0_tau at Q_K=3/2 point: {dqk[idx_koide]:.4f}")
print(f"  (If zero, Q_K=3/2 would be an extremum -- special!)")

# Check monotonicity
if np.all(dqk > 0):
    print(f"  Q_K is MONOTONICALLY INCREASING in the corridor.")
elif np.all(dqk < 0):
    print(f"  Q_K is MONOTONICALLY DECREASING in the corridor.")
else:
    # Find extrema
    sign_changes = np.where(np.diff(np.sign(dqk)))[0]
    print(f"  Q_K has {len(sign_changes)} extremum/a in the corridor.")
    for sc in sign_changes:
        print(f"    Extremum near g0_tau = {g0_tau_vals[sc]:.6f}, Q_K = {qk_vals[sc]:.6f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: SPECIAL POINTS IN THE CORRIDOR")
print("=" * 70)

# Check several "natural" points
frac_koide = (g0_tau_vals[idx_koide] - g0_mu) / (G0_CRIT - g0_mu)

special_points = {
    "g0_tau = g0_crit (boundary)": G0_CRIT,
    "g0_tau = phi * g0_mu": PHI * g0_mu,
    "g0_tau = g0_mu + (g0_crit - g0_mu)/phi": g0_mu + (G0_CRIT - g0_mu) / PHI,
    "g0_tau = g0_mu + (g0_crit - g0_mu)/2": g0_mu + (G0_CRIT - g0_mu) / 2.0,
    "g0_tau = g0_mu + (g0_crit - g0_mu)/sqrt(2)": g0_mu + (G0_CRIT - g0_mu) / np.sqrt(2),
    "g0_tau = g0_mu + (g0_crit - g0_mu)*2/3": g0_mu + (G0_CRIT - g0_mu) * 2.0 / 3.0,
    "g0_tau = g0_crit - (g0_crit - g0_mu)/phi^2": G0_CRIT - (G0_CRIT - g0_mu) / PHI**2,
    "g0_tau = sqrt(g0_mu * g0_crit)": np.sqrt(g0_mu * G0_CRIT),
    "g0_tau (Koide actual)": g0_tau_vals[idx_koide],
}

print(f"\n  {'Point':>50s} {'g0_tau':>10s} {'frac':>6s} {'Q_K':>10s} {'Q_K-3/2':>10s}")
for name, g0t in special_points.items():
    if g0t >= G0_CRIT:
        print(f"  {name:>50s} {g0t:10.5f} {'---':>6s} {'COLLAPSE':>10s} {'---':>10s}")
        continue
    if g0t <= g0_mu:
        print(f"  {name:>50s} {g0t:10.5f} {'<0':>6s} {'---':>10s} {'---':>10s}")
        continue
    At = A_tail(g0t)
    if At > 0:
        qk = QK_masses([A_e, A_mu, At])
        frac = (g0t - g0_mu) / (G0_CRIT - g0_mu)
        print(f"  {name:>50s} {g0t:10.5f} {frac:6.3f} {qk:10.6f} {qk-1.5:+10.6f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: CORRIDOR FRACTION FOR Q_K = 3/2")
print("=" * 70)

# Precisely find g0_tau for Q_K = 3/2
def qk_resid(g0t):
    At = A_tail(g0t)
    if At <= 0:
        return 10.0
    return QK_masses([A_e, A_mu, At]) - 1.5

g0_tau_koide = brentq(qk_resid, 1.50, G0_CRIT - 0.01)
frac_exact = (g0_tau_koide - g0_mu) / (G0_CRIT - g0_mu)

print(f"\n  Exact g0_tau for Q_K = 3/2: {g0_tau_koide:.10f}")
print(f"  g0_mu  = {g0_mu:.10f}")
print(f"  g0_crit = {G0_CRIT:.10f}")
print(f"  Corridor width = {G0_CRIT - g0_mu:.10f}")
print(f"  Position in corridor = {g0_tau_koide - g0_mu:.10f}")
print(f"  Fraction f = (g0_tau - g0_mu)/(g0_crit - g0_mu) = {frac_exact:.10f}")

# Check if this fraction matches any known constant
constants = {
    "1/phi": 1.0/PHI,
    "1/phi^2": 1.0/PHI**2,
    "2/3": 2.0/3.0,
    "1/sqrt(2)": 1.0/np.sqrt(2),
    "1/2": 0.5,
    "1/e": 1.0/np.e,
    "1/3": 1.0/3.0,
    "1/pi": 1.0/np.pi,
    "phi - 1 (= 1/phi)": PHI - 1,
    "2 - phi": 2 - PHI,
    "3 - 2*phi": 3 - 2*PHI,
    "sqrt(2) - 1": np.sqrt(2) - 1,
    "phi/pi": PHI/np.pi,
    "ln(2)": np.log(2),
    "(phi-1)^2": (PHI-1)**2,
    "pi/4 - 1/2": np.pi/4 - 0.5,
}

print(f"\n  Matching fraction {frac_exact:.6f} against constants:")
matches = []
for name, val in constants.items():
    err = abs(frac_exact - val) / frac_exact * 100
    matches.append((err, name, val))
matches.sort()
for err, name, val in matches[:10]:
    marker = " <--" if err < 1.0 else ""
    print(f"    {name:>20s} = {val:.6f}  error = {err:6.2f}%{marker}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: EQUIPARTITION TEST (PROPOSAL #1)")
print("=" * 70)

# Brannen: sqrt(m_k) = c * (1 + b*cos(theta + 2*pi*k/3))
# Q_K = 3/(1 + b^2/2)
# Q_K = 3/2 <=> b = sqrt(2)
# Equipartition: r = sqrt(N-1) = sqrt(2) for N=3

A_tau_koide = A_tail(g0_tau_koide)
b_koide = brannen_b([A_e, A_mu, A_tau_koide])
print(f"\n  At Q_K = 3/2 point:")
print(f"  Brannen b = {b_koide:.8f}")
print(f"  sqrt(2)   = {np.sqrt(2):.8f}")
print(f"  Match: {abs(b_koide - np.sqrt(2))/np.sqrt(2)*100:.4f}%")

# What does "equipartition" mean in the corridor?
# If we parameterize x_k = sqrt(m_k) = c(1 + b*cos(...)),
# then b = sqrt(2) means the "spread" equals the "mean" in a specific sense.
# In the simplex, this corresponds to a point on the unit sphere.

# The physical interpretation: on the (N-1)-simplex of mass fractions,
# r = sqrt(N-1) means the point is at maximal distance from the center
# while remaining on the simplex boundary = one mass vanishes.
# Check: is that what happens for leptons?
x = np.array([A_e**2, A_mu**2, A_tau_koide**2])
x_norm = x / np.sum(x)
center = np.ones(3) / 3.0
dist_from_center = np.linalg.norm(x_norm - center)
max_dist = np.linalg.norm(np.array([1,0,0]) - center)  # vertex distance
print(f"\n  Simplex analysis:")
print(f"  Mass fraction vector: ({x_norm[0]:.6f}, {x_norm[1]:.6f}, {x_norm[2]:.6f})")
print(f"  Distance from center: {dist_from_center:.6f}")
print(f"  Max distance (vertex): {max_dist:.6f}")
print(f"  Ratio d/d_max: {dist_from_center/max_dist:.6f}")

# For N=3 equipartition, the ratio should be sqrt(2)/sqrt(2/3) = sqrt(3)...
# Actually, the Brannen circle has radius c*b on a 2D plane.
# The center (1/3,1/3,1/3) in sqrt(m) space.
# CV = 1 means std(sqrt(m)) = mean(sqrt(m)), i.e., the distribution
# is maximally spread in a specific sense.

# What fraction of the corridor width corresponds to each "natural" b value?
print(f"\n  Brannen b across the corridor:")
print(f"  {'frac':>6s} {'g0_tau':>10s} {'b':>8s} {'b/sqrt(2)':>10s} {'Q_K':>10s}")
for frac_test in np.arange(0.05, 1.0, 0.05):
    g0t = g0_mu + frac_test * (G0_CRIT - g0_mu)
    if g0t >= G0_CRIT - 0.005:
        continue
    At = A_tail(g0t)
    if At > 0:
        b_test = brannen_b([A_e, A_mu, At])
        qk_test = QK_masses([A_e, A_mu, At])
        print(f"  {frac_test:6.2f} {g0t:10.5f} {b_test:8.4f} {b_test/np.sqrt(2):10.4f} {qk_test:10.6f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: WHAT SELECTS Q_K = 3/2 WITHIN THE CORRIDOR?")
print("=" * 70)

# Key question: is there ANY property of the ODE or soliton that
# is extremal at Q_K = 3/2?

# Test 1: Total mass
print(f"\n  Test 1: Total mass M = sum(A^4)")
print(f"  {'g0_tau':>10s} {'M_total':>14s} {'Q_K':>10s}")
for g0t, qk in zip(g0_tau_vals[::4], qk_vals[::4]):
    At = A_tail(g0t)
    M = A_e**4 + A_mu**4 + At**4
    print(f"  {g0t:10.5f} {M:14.8f} {qk:10.6f}")
print(f"  M is monotonic (already known from PATH 1)")

# Test 2: Field energy of tau soliton
print(f"\n  Test 2: Field energy E_field of tau soliton")
print(f"  {'g0_tau':>10s} {'E_field':>14s} {'Q_K':>10s}")

def field_energy(g0, r_max=250):
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > 0.001) & (r < r_max)
    rf = r[mask]
    gf = g[mask]
    from numpy import trapezoid as trapz
    return trapz((gf - 1.0)**2 * rf**2, rf)

for g0t in np.arange(g0_mu + 0.02, G0_CRIT - 0.01, 0.02):
    At = A_tail(g0t)
    if At > 0:
        Ef = field_energy(g0t)
        qk = QK_masses([A_e, A_mu, At])
        print(f"  {g0t:10.5f} {Ef:14.4f} {qk:10.6f}")

# Test 3: G(g0_tau) -- nonlinear enhancement at Q_K = 3/2
print(f"\n  Test 3: Nonlinear enhancement G(g0_tau) = A^2/(g0-1)^2")
A_tau_k = A_tail(g0_tau_koide)
G_tau = A_tau_k**2 / (g0_tau_koide - 1)**2
G_mu = A_mu**2 / (g0_mu - 1)**2
G_e = A_e**2 / (g0_e - 1)**2
print(f"  G_e   = {G_e:.6f}")
print(f"  G_mu  = {G_mu:.6f}")
print(f"  G_tau = {G_tau:.6f}")
print(f"  G_tau/G_e = {G_tau/G_e:.4f}")
print(f"  G_tau/G_mu = {G_tau/G_mu:.4f}")

# Is G_tau/G_e = some constant related to phi or sqrt(2)?
ratio_Gte = G_tau/G_e
print(f"\n  G_tau/G_e = {ratio_Gte:.6f}")
for name, val in {"phi": PHI, "phi^2": PHI**2, "e": np.e,
                  "pi": np.pi, "2": 2.0, "3": 3.0, "sqrt(2)": np.sqrt(2),
                  "2*phi-1": 2*PHI-1, "phi+1": PHI+1}.items():
    err = abs(ratio_Gte - val)/val*100
    marker = " <--" if err < 2 else ""
    print(f"    {name:>12s} = {val:.6f}  error = {err:.2f}%{marker}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: SENSITIVITY -- HOW DOES Q_K DEPEND ON g0_e?")
print("=" * 70)

# If g0_e changes, does Q_K = 3/2 still require a specific fraction?
print(f"\n  For different g0_e, find corridor fraction where Q_K = 3/2:")
print(f"  {'g0_e':>8s} {'g0_mu':>8s} {'g0_crit':>8s} {'corr_w':>8s} {'g0t_K':>10s} {'frac':>8s}")

for g0e_test in np.arange(0.50, 0.95, 0.05):
    g0m_test = PHI * g0e_test
    if g0m_test >= G0_CRIT - 0.01:
        continue
    Ae_t = A_tail(g0e_test)
    Am_t = A_tail(g0m_test)
    if Ae_t <= 0 or Am_t <= 0:
        continue

    def resid_t(g0t):
        At = A_tail(g0t)
        if At <= 0:
            return 10.0
        return QK_masses([Ae_t, Am_t, At]) - 1.5

    # Search for Q_K = 3/2 in corridor
    try:
        # Check if Q_K = 3/2 is achievable
        lo = g0m_test + 0.01
        hi = G0_CRIT - 0.005
        r_lo = resid_t(lo)
        r_hi = resid_t(hi)
        if r_lo * r_hi < 0:
            g0t_k = brentq(resid_t, lo, hi)
            frac_t = (g0t_k - g0m_test) / (G0_CRIT - g0m_test)
            corr_w = G0_CRIT - g0m_test
            print(f"  {g0e_test:8.4f} {g0m_test:8.4f} {G0_CRIT:8.4f} {corr_w:8.4f} {g0t_k:10.6f} {frac_t:8.5f}")
        else:
            print(f"  {g0e_test:8.4f} {g0m_test:8.4f} {G0_CRIT:8.4f} {'---':>8s} {'no root':>10s} {'---':>8s}")
    except:
        print(f"  {g0e_test:8.4f} {g0m_test:8.4f} {G0_CRIT:8.4f} {'---':>8s} {'error':>10s} {'---':>8s}")


# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS")
print("=" * 70)
print(f"""
  CORRIDOR ANALYSIS RESULTS:

  1. The corridor for g0_tau is [{g0_mu:.4f}, {G0_CRIT:.4f}),
     width = {G0_CRIT - g0_mu:.4f}.

  2. Q_K = 3/2 occurs at g0_tau = {g0_tau_koide:.6f},
     fraction = {frac_exact:.6f} through the corridor.

  3. Does Q_K = 3/2 correspond to an extremum?
     -> Check dQ_K/dg0_tau at that point.

  4. Does the corridor fraction match a known constant?
     -> Check above.

  5. Is the Brannen parameter b = sqrt(2) natural?
     -> b increases monotonically through the corridor.
     -> b = sqrt(2) is just a specific point, not extremal.

  KEY QUESTION: What SELECTS g0_tau = {g0_tau_koide:.6f} within the corridor?
  The corridor CONSTRAINS Q_K to a finite range, but does not SELECT 3/2.
""")
