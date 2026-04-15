#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_barrier_structural.py -- Structural analysis of g0_crit(d)
"""

import numpy as np
from scipy.optimize import curve_fit
from math import gamma, sqrt, pi, log, exp, e

PHI = (1+sqrt(5))/2

# High-precision values:
g0c = {
    1: 4/3,
    2: 1.732438098578,
    3: 2.206189384431,
    4: 2.764548578180,
    5: 3.418676157100,
    6: 4.181058454400,
    7: 5.065648932900,
    8: 6.088025137600,
}

ds = sorted(g0c.keys())

print("=" * 70)
print("  STRUCTURAL ANALYSIS OF g0_crit(d)")
print("=" * 70)

# 1. Consecutive ratios
print("\n1. CONSECUTIVE RATIOS")
print("  If g0_crit(d+1) = R(d) * g0_crit(d), what is R(d)?")
print(f"  {'d':>3s}  {'g0_crit':>14s}  {'R(d)':>14s}")
for i in range(len(ds)-1):
    d1, d2 = ds[i], ds[i+1]
    r = g0c[d2]/g0c[d1]
    print(f"  {d1:3d}  {g0c[d1]:14.10f}  {r:14.10f}")
print(f"  {ds[-1]:3d}  {g0c[ds[-1]]:14.10f}")

# 2. Power law
print("\n2. POWER LAW: g0c = a * d^b")
log_d = np.log(np.array(ds, dtype=float))
log_g = np.log(np.array([g0c[d] for d in ds]))
coeffs = np.polyfit(log_d, log_g, 1)
b_fit, a_log = coeffs
a_fit = exp(a_log)
print(f"  a = {a_fit:.6f}, b = {b_fit:.6f}")
for d in ds:
    pred = a_fit * d**b_fit
    print(f"    d={d}: pred={pred:.6f}, actual={g0c[d]:.6f}, err={abs(pred-g0c[d])/g0c[d]*100:.2f}%")

# 3. g0c/d
print("\n3. g0_crit/d:")
for d in ds:
    print(f"    d={d}: g0c/d = {g0c[d]/d:.8f}")

# 4. Gamma function
print("\n4. GAMMA FUNCTION CONNECTION")
for d in ds:
    val_gamma = gamma(d/2.0 + 1) / gamma(d/2.0)
    ratio = g0c[d] / val_gamma
    print(f"  d={d}: G({d}/2+1)/G({d}/2) = {val_gamma:.6f}, g0c/val = {ratio:.6f}")

# 5. Solid angle
print("\n5. SOLID ANGLE S_d = 2*pi^(d/2) / Gamma(d/2)")
for d in ds:
    S_d = 2 * pi**(d/2.0) / gamma(d/2.0)
    V_d = pi**(d/2.0) / gamma(d/2.0 + 1) if d > 0 else 1
    print(f"  d={d}: S_d={S_d:.4f}, V_d(unit ball)={V_d:.4f}, g0c*S_d={g0c[d]*S_d:.4f}, g0c/V_d={g0c[d]/V_d:.4f}")

# 6. Deviations from (4/pi)*sqrt(d)
print("\n6. DEVIATIONS FROM (4/pi)*sqrt(d)")
for d in ds:
    pred = (4/pi) * sqrt(d)
    delta = g0c[d] - pred
    delta_rel = delta / pred
    print(f"  d={d}: g0c={g0c[d]:.6f}, pred={pred:.6f}, delta={delta:+.6f} ({delta_rel:+.4f})")

# 7. Exact constraint
print("\n7. CONSTRAINT: g0_crit(1) = 4/3 exactly")
print(f"  (4/pi)*sqrt(1) = {4/pi:.6f} != 4/3 = {4/3:.6f}")
print(f"  So g_bar = (4/pi)*sqrt(d) fails at d=1 by {abs(4/pi - 4/3)/(4/3)*100:.2f}%")

# 8. Fit g0c = C * (d + delta)^alpha
print("\n8. FIT: g0_crit(d) = C * (d + delta)^alpha")
def model_power(d, C, delta, alpha):
    return C * (d + delta)**alpha

ds_arr = np.array(ds, dtype=float)
g0c_arr = np.array([g0c[d] for d in ds])
popt, _ = curve_fit(model_power, ds_arr, g0c_arr, p0=[1.0, 0.1, 0.7], maxfev=10000)
C, delta, alpha = popt
print(f"  C = {C:.8f}, delta = {delta:.8f}, alpha = {alpha:.8f}")
for d in ds:
    pred = model_power(d, C, delta, alpha)
    print(f"  d={d}: pred={pred:.8f}, actual={g0c[d]:.8f}, err={abs(pred-g0c[d])/g0c[d]*100:.6f}%")

print(f"\n  delta candidates:")
for name, val in [("1/pi", 1/pi), ("1/3", 1/3), ("1/4", 1/4), ("1/e", 1/e),
                   ("pi/10", pi/10), ("1/pi^2", 1/pi**2), ("0", 0),
                   ("1/(2*pi)", 1/(2*pi)), ("phi-1", PHI-1), ("1/phi", 1/PHI)]:
    print(f"    {name:10s} = {val:.8f}  (diff {abs(delta-val):.4e})")

print(f"\n  alpha candidates:")
for name, val in [("1/2", 0.5), ("3/4", 0.75), ("2/3", 2/3),
                   ("1/sqrt(2)", 1/sqrt(2)), ("ln(2)", log(2)),
                   ("pi/4", pi/4), ("1/phi", 1/PHI)]:
    print(f"    {name:10s} = {val:.8f}  (diff {abs(alpha-val):.4e})")

# 9. Test: g0c(d) = (4/3) * prod_{k=1}^{d-1} R_k
print("\n9. PRODUCT STRUCTURE")
print("  g0c(d) = g0c(1) * R_1 * R_2 * ... * R_{d-1}")
print(f"  {'d':>3s}  {'R_{d-1}':>14s}  {'1+2/(3d-1)':>14s}  {'diff%':>8s}")
for i in range(len(ds)-1):
    d1, d2 = ds[i], ds[i+1]
    r = g0c[d2]/g0c[d1]
    # Is R related to d?
    # Try R = 1 + 2/(3d-1), R = 1 + 1/d, R = (d+1)/d, etc
    r_try = 1 + 2/(3*d1 - 1)
    diff = abs(r - r_try)/r * 100
    print(f"  {d1:3d}  {r:14.10f}  {r_try:14.10f}  {diff:8.4f}")

# 10. The key ratio g0_crit(3)/g0_crit(2) ~ 4/pi
print("\n10. THE KEY RATIO g0_crit(3)/g0_crit(2)")
g2, g3 = g0c[2], g0c[3]
r32 = g3/g2
print(f"  g0c(3)/g0c(2) = {r32:.12f}")
print(f"  4/pi           = {4/pi:.12f}")
print(f"  diff           = {r32 - 4/pi:.6e} ({abs(r32-4/pi)/(4/pi)*100:.6f}%)")

# What if: g0c(d) = (4/pi)^(d-2) * sqrt(3) * correction?
print("\n11. GEOMETRIC PROGRESSION WITH BASE 4/pi")
print("  If g0c(d) = sqrt(3) * (4/pi)^(d-2):")
for d in ds:
    pred = sqrt(3) * (4/pi)**(d-2)
    err = abs(pred - g0c[d])/g0c[d]*100
    print(f"  d={d}: pred={pred:.8f}, actual={g0c[d]:.8f}, err={err:.4f}%")

# 12. What about g0c(d) = (4/pi)^(d-1) * (pi/3)?
# g0c(1) = pi/3 * 1 = pi/3 = 1.047 != 4/3
# No, that doesn't work either.
# g0c(1) = 4/3, g0c(d) = (4/3) * product
# The product ratio is ~1.27-1.30, decreasing
#
# What if we use the Q formulation directly?
# Q_d = pi*g0c(d)/4
# Q_1 = pi/3 (exact)
# Q_2 = 1.360654 (close to pi*sqrt(3)/4 = 1.360350)
# Q_3 = 1.732737 (close to sqrt(3) = 1.732051)
#
# Note: Q_3/Q_2 = 1.273459 = 4/pi to 0.017%
# And Q_2/Q_1 = 1.29933 = 3*sqrt(3)/4 * (4/pi) / 1 hmm

print("\n12. Q VALUES (Q_d = pi*g0c/4)")
for d in ds:
    Q = pi * g0c[d] / 4
    print(f"  d={d}: Q_{d} = {Q:.10f}")

print("\n  Q ratios:")
for i in range(len(ds)-1):
    d1, d2 = ds[i], ds[i+1]
    Q1 = pi*g0c[d1]/4
    Q2 = pi*g0c[d2]/4
    r = Q2/Q1
    print(f"  Q_{d2}/Q_{d1} = {r:.10f}")

# 13. Summary
print("\n" + "="*70)
print("  SUMMARY")
print("="*70)
print(f"""
  g0_crit(d) for substrate ODE in d dimensions:

  d=1: {g0c[1]:.12f}  = 4/3 EXACTLY (conservation law)
  d=2: {g0c[2]:.12f}  ~ sqrt(3)   (0.022%)
  d=3: {g0c[3]:.12f}  ~ 4*sqrt(3)/pi (0.040%)

  Ratios:
    g0c(3)/g0c(2) = {g0c[3]/g0c[2]:.10f} ~ 4/pi = {4/pi:.10f} (0.017%)
    g0c(2)/g0c(1) = {g0c[2]/g0c[1]:.10f}

  Best fit: g0c(d) = {C:.6f} * (d + {delta:.6f})^{alpha:.6f}

  The (4/pi)*Q_d framework:
    Q_1 = pi/3  (EXACT)
    Q_2 = {pi*g0c[2]/4:.8f}  ~ pi*sqrt(3)/4 = {pi*sqrt(3)/4:.8f} (0.022%)
    Q_3 = {pi*g0c[3]/4:.8f}  ~ sqrt(3)      = {sqrt(3):.8f} (0.040%)

  INTERPRETATION:
    The pattern g_bar = (4/pi)*Q_d with Q_d ~ sqrt(d)-like
    is SUGGESTIVE but not EXACT (deviations 0.02-0.04%).

    The deviations are REAL (bracket precision 1e-11).

    However, the closeness to clean algebraic values
    (especially sqrt(3) and 4/pi ratio) suggests there MAY be
    an exact formula with small corrections from the nonlinear
    (1/g)g'^2 term in the ODE.
""")
