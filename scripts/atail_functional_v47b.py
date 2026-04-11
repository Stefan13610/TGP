#!/usr/bin/env python3
"""
atail_functional_v47b.py -- Study A_tail(g0) functional form for Form A.

Goal: understand WHY the specific functional form of A_tail produces
Koide Q_K = 3/2 with phi-FP spacing.

Key structure:
- A_tail has minimum near g0 = 1 (vacuum)
- Grows as g0 -> 0 (linear regime)
- DIVERGES as g0 -> 8/5 (collapse threshold)
- Electron sits below vacuum, muon and tau above
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0  # = 1.6 exactly


def solver_A(g0, r_max=300):
    """Canonical Form A ODE solver."""
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


def A_tail(g0):
    r, g = solver_A(g0)
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
# Section 1: A_tail on fine grid
# ================================================================
print("=" * 65)
print("A_tail(g0) FUNCTIONAL FORM -- CANONICAL FORM A")
print("=" * 65)

g0_grid = np.concatenate([
    np.linspace(0.2, 0.95, 16),
    np.linspace(0.96, 1.04, 9),
    np.linspace(1.05, 1.55, 20),
    np.linspace(1.56, 1.598, 15),
])
g0_grid = np.sort(np.unique(g0_grid))

print("\n  g0 vs A_tail:")
print("  %-10s %-14s %-10s %-10s" % ("g0", "A_tail", "|g0-1|", "8/5-g0"))
results = []
for g0 in g0_grid:
    A = A_tail(g0)
    results.append((g0, A))
    d_vac = abs(g0 - 1.0)
    d_crit = G0_CRIT - g0
    print("  %-10.5f %-14.8f %-10.4f %-10.4f" % (g0, A, d_vac, d_crit))


# ================================================================
# Section 2: Asymptotic near vacuum (g0 -> 1)
# ================================================================
print("\n  ASYMPTOTIC: g0 -> 1 (vacuum)")
print("  " + "-" * 50)

for g0, A in results:
    d = abs(g0 - 1.0)
    if 0.005 < d < 0.2 and A > 1e-10:
        ratio = A / d
        print("  g0=%.4f  delta=%.4f  A/delta=%.6f" % (g0, d, ratio))

# Near vacuum: A_tail is approximately proportional to |g0-1|
# but with corrections


# ================================================================
# Section 3: Asymptotic near collapse (g0 -> 8/5)
# ================================================================
print("\n  ASYMPTOTIC: g0 -> 8/5 (collapse)")
print("  " + "-" * 50)

deltas = []
log_As = []
for g0, A in results:
    d = G0_CRIT - g0
    if 0 < d < 0.05 and A > 1e-10:
        deltas.append(d)
        log_As.append(np.log(A))
        print("  g0=%.5f  delta=%.5f  A=%.6f" % (g0, d, A))

if len(deltas) > 3:
    deltas = np.array(deltas)
    log_As = np.array(log_As)
    log_ds = np.log(deltas)
    p = np.polyfit(log_ds, log_As, 1)
    gamma = -p[0]
    print("\n  Power law: A_tail = C * (8/5 - g0)^(-%.4f)" % gamma)
    print("  (exponent gamma = %.4f)" % gamma)


# ================================================================
# Section 4: A_tail in canonical variable f = g^3
# ================================================================
print("\n  CANONICAL VARIABLE f = g^3")
print("  " + "-" * 50)

# In terms of f: g = f^(1/3), g0 = f0^(1/3)
# A_tail for g corresponds to A_tail_f for f (with factor 3)
# Since g = 1 + h, f = (1+h)^3 = 1 + 3h + ..., so A_tail_f = 3*A_tail_g

g0_e = 0.8676968610
g0_mu = PHI * g0_e
g0_tau = 1.5695724137

f0_e = g0_e**3
f0_mu = g0_mu**3
f0_tau = g0_tau**3
f0_crit = G0_CRIT**3  # = 512/125 = 4.096

print("  f0_e   = g0_e^3   = %.8f" % f0_e)
print("  f0_mu  = g0_mu^3  = %.8f" % f0_mu)
print("  f0_tau = g0_tau^3 = %.8f" % f0_tau)
print("  f0_crit = (8/5)^3 = %.8f = 512/125" % f0_crit)
print()
print("  Ratios in f-space:")
print("    f0_mu/f0_e   = %.8f" % (f0_mu / f0_e))
print("    f0_tau/f0_e  = %.8f" % (f0_tau / f0_e))
print("    f0_crit/f0_e = %.8f" % (f0_crit / f0_e))
print("    phi^3        = %.8f" % (PHI**3))
print("    f0_mu/f0_e vs phi^3: %.6f" % (f0_mu / f0_e / PHI**3))


# ================================================================
# Section 5: Ratio A(phi*g0)/A(g0) as function of g0
# ================================================================
print("\n  RATIO R(g0) = A(phi*g0) / A(g0)")
print("  " + "-" * 50)

# This ratio determines r_21 via R^4 = r_21
# If R were constant, all g0_e values would give the same r_21
# The specific g0_e is where R^4 = 206.768

g0_test = np.linspace(0.4, 0.95, 20)
print("  %-8s %-12s %-12s %-12s" % ("g0", "A(g0)", "A(phi*g0)", "R^4"))
for g0 in g0_test:
    A1 = A_tail(g0)
    g0_2 = PHI * g0
    if g0_2 < G0_CRIT:
        A2 = A_tail(g0_2)
        if A1 > 1e-10:
            R = A2 / A1
            print("  %-8.4f %-12.8f %-12.8f %-12.2f" % (g0, A1, A2, R**4))


# ================================================================
# Section 6: Key insight - separation of scales
# ================================================================
print("\n  SCALE SEPARATION ANALYSIS")
print("  " + "-" * 50)

A_e = A_tail(g0_e)
A_mu = A_tail(g0_mu)
A_tau = A_tail(g0_tau)

# In log space
log_Ae = np.log(A_e)
log_Amu = np.log(A_mu)
log_Atau = np.log(A_tau)

print("  ln(A_e)   = %.6f" % log_Ae)
print("  ln(A_mu)  = %.6f" % log_Amu)
print("  ln(A_tau) = %.6f" % log_Atau)
print()
print("  Spacings in log-A:")
print("    ln(A_mu) - ln(A_e)  = %.6f" % (log_Amu - log_Ae))
print("    ln(A_tau) - ln(A_mu) = %.6f" % (log_Atau - log_Amu))
print("    Ratio of spacings: %.4f" % (
    (log_Atau - log_Amu) / (log_Amu - log_Ae)))
print()

# Check if spacings relate to phi
print("  Comparison with phi:")
print("    ln(phi) = %.6f" % np.log(PHI))
print("    Spacing ratio vs 1/phi = %.6f" % (
    (log_Atau - log_Amu) / (log_Amu - log_Ae) / (1/PHI)))

# x_i = A_i^2 for Koide
x_e = A_e**2
x_mu = A_mu**2
x_tau = A_tau**2

print()
print("  Koide variables x = A^2 (proportional to sqrt(mass)):")
print("    x_e   = %.10f" % x_e)
print("    x_mu  = %.10f" % x_mu)
print("    x_tau = %.10f" % x_tau)
print("    Q_K = (sum x)^2 / (sum x^2) = %.10f" % (
    (x_e + x_mu + x_tau)**2 / (x_e**2 + x_mu**2 + x_tau**2)))


# ================================================================
# Section 7: Summary
# ================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print("""
  A_tail(g0) for canonical Form A has structure:
  - Minimum near g0 = 1 (vacuum): A_tail -> 0 linearly
  - Below vacuum (g0 < 1): electron regime
  - Above vacuum (g0 > 1): muon, tau regime
  - Diverges at g0 -> 8/5: A_tail ~ (8/5 - g0)^(-gamma)

  The tau soliton sits in the DIVERGENCE regime,
  drawing its large mass from proximity to collapse.

  Koide Q_K = 3/2 constrains A_tau relative to A_e, A_mu.
  The specific value 3/2 arises from the nonlinear ODE
  dynamics -- not from any simple approximation.
""")
