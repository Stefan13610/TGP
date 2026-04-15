#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_barrier_Qd.py -- Test: g_bar(d) = (4/pi) * Q_d pattern

User's observation:
  g_bar_1D = 4/3         => Q_1 = pi/3
  g_bar_3D = 2.206       => Q_3 ~ sqrt(3) (diff 0.03%)

Question: what is Q_2? Is it a clean number?

Author: Claudian
Date: 2026-04-15
"""

import numpy as np
from scipy.integrate import solve_ivp
import math

PHI = (1 + math.sqrt(5)) / 2

# ================================================================
# SOLVER (from r3_g0crit_analytical.py — proven to work)
# ================================================================

def solve_substrate_d(g0, d, r_max=200.0, n_points=20000, g_floor=1e-8):
    """Substrate ODE in d dimensions"""
    g_min_val = [g0]
    singular = [False]

    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        if g < g_min_val[0]:
            g_min_val[0] = g

        if r < 1e-12:
            gpp = (1.0 - g) / (d if d > 0 else 1.0)
        else:
            gpp = (1.0 - g) - (1.0/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]

    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)

    if not sol.success:
        return g_min_val[0], singular[0]

    actual_min = np.min(sol.y[0])
    if actual_min < g_min_val[0]:
        g_min_val[0] = actual_min

    return g_min_val[0], singular[0]


def find_g0_crit(d, g0_lo=1.01, g0_hi=5.0, tol=1e-8, g_threshold=0.001):
    """Find g0_crit by bisection (from r3_g0crit_analytical.py)"""
    # Ensure g0_hi gives singularity
    g_min, sing = solve_substrate_d(g0_hi, d)
    if not sing and g_min > g_threshold:
        g0_hi = 10.0
        g_min, sing = solve_substrate_d(g0_hi, d)
        if not sing and g_min > g_threshold:
            g0_hi = 50.0
            g_min, sing = solve_substrate_d(g0_hi, d)
            if not sing and g_min > g_threshold:
                return None

    for _ in range(80):
        g0_mid = (g0_lo + g0_hi) / 2
        g_min, sing = solve_substrate_d(g0_mid, d)

        if sing or g_min < g_threshold:
            g0_hi = g0_mid
        else:
            g0_lo = g0_mid

        if g0_hi - g0_lo < tol:
            break

    return (g0_lo + g0_hi) / 2


# ================================================================
print("=" * 70)
print("  BARRIER PATTERN: g_bar(d) = (4/pi) * Q_d")
print("=" * 70)

# ================================================================
# SECTION 1: High-precision g0_crit for d = 1..8
# ================================================================
print(f"\n{'='*70}")
print("  1. g0_crit(d) WITH HIGH PRECISION")
print("="*70)

results = {}

for d in range(1, 9):
    if d == 1:
        g0c = 4.0 / 3.0
        print(f"  d={d}: g0_crit = {g0c:.10f}  [EXACT: 4/3]")
    else:
        g0c = find_g0_crit(d)
        if g0c is None:
            print(f"  d={d}: g0_crit NOT FOUND")
            continue
        print(f"  d={d}: g0_crit = {g0c:.10f}")

    Q = math.pi * g0c / 4.0
    results[d] = (g0c, Q)
    print(f"       Q_{d} = pi*g0_crit/4 = {Q:.10f}")


# ================================================================
# SECTION 2: Pattern analysis
# ================================================================
print(f"\n{'='*70}")
print("  2. PATTERN ANALYSIS")
print("="*70)

print(f"\n  {'d':>3s}  {'g0_crit':>14s}  {'Q_d':>14s}  {'sqrt(d)':>10s}  {'Q/sqrt(d)':>12s}")
print(f"  {'-'*3}  {'-'*14}  {'-'*14}  {'-'*10}  {'-'*12}")

for d in sorted(results.keys()):
    g0c, Q = results[d]
    sqrtd = math.sqrt(d)
    ratio = Q / sqrtd
    print(f"  {d:3d}  {g0c:14.10f}  {Q:14.10f}  {sqrtd:10.6f}  {ratio:12.8f}")


# ================================================================
# SECTION 3: Focus on Q_d values
# ================================================================
print(f"\n{'='*70}")
print("  3. WHAT ARE THE Q VALUES?")
print("="*70)

for d in sorted(results.keys()):
    g0c, Q = results[d]
    print(f"\n  d={d}: Q_{d} = {Q:.10f}")

    candidates = [
        ("pi/3", math.pi/3),
        ("pi/sqrt(5)", math.pi/math.sqrt(5)),
        ("pi/sqrt(6)", math.pi/math.sqrt(6)),
        ("sqrt(2)", math.sqrt(2)),
        ("sqrt(3)", math.sqrt(3)),
        ("sqrt(5)", math.sqrt(5)),
        ("sqrt(6)", math.sqrt(6)),
        ("sqrt(7)", math.sqrt(7)),
        ("phi", PHI),
        ("2*phi/pi", 2*PHI/math.pi),
        ("pi/phi", math.pi/PHI),
        ("pi/2", math.pi/2),
        ("pi^2/6", math.pi**2/6),
        ("pi^2/8", math.pi**2/8),
        ("sqrt(pi)", math.sqrt(math.pi)),
        ("4/3", 4/3),
        ("3/2", 3/2),
        ("5/3", 5/3),
        ("7/4", 7/4),
        ("2", 2.0),
        ("e/2", math.e/2),
        ("e/phi", math.e/PHI),
        ("sqrt(d)", math.sqrt(d)),
        ("(d+1)/d * pi/3", (d+1)/d * math.pi/3),
        ("pi*sqrt(d)/(d+1)", math.pi*math.sqrt(d)/(d+1)),
        ("sqrt(d+1/3)", math.sqrt(d + 1/3)),
        ("(4/pi)*sqrt(d/3)", (4/math.pi)*math.sqrt(d/3)),
        ("pi/(d+1)", math.pi/(d+1)),
        ("sqrt(2*d/3)", math.sqrt(2*d/3)),
        ("sqrt((d+2)/3)", math.sqrt((d+2)/3)),
    ]

    hits = []
    for name, val in candidates:
        rel = abs(Q - val) / val * 100
        if rel < 2.0:
            hits.append((rel, name, val))

    hits.sort()
    for rel, name, val in hits[:5]:
        print(f"    ~ {name:30s} = {val:.10f}  ({rel:.4f}%)")


# ================================================================
# SECTION 4: Test g_bar = (4/pi)*sqrt(d) directly
# ================================================================
print(f"\n{'='*70}")
print("  4. TEST: g_bar(d) = (4/pi)*sqrt(d)")
print("="*70)

print(f"\n  {'d':>3s}  {'g0_crit':>14s}  {'(4/pi)sqrt(d)':>14s}  {'diff%':>8s}")
print(f"  {'-'*3}  {'-'*14}  {'-'*14}  {'-'*8}")

for d in sorted(results.keys()):
    g0c = results[d][0]
    pred = (4/math.pi) * math.sqrt(d)
    rel = abs(g0c - pred) / g0c * 100
    print(f"  {d:3d}  {g0c:14.10f}  {pred:14.10f}  {rel:8.4f}")


# ================================================================
# SECTION 5: Deeper search for Q_2
# ================================================================
print(f"\n{'='*70}")
print("  5. DEEP SEARCH FOR Q_2")
print("="*70)

if 2 in results:
    Q2 = results[2][1]
    g0c_2 = results[2][0]

    print(f"\n  g0_crit(2D) = {g0c_2:.12f}")
    print(f"  Q_2 = pi * g0_crit(2D) / 4 = {Q2:.12f}")
    print(f"\n  sqrt(3) = {math.sqrt(3):.12f}")
    print(f"  g0_crit(2D) vs sqrt(3): diff = {abs(g0c_2 - math.sqrt(3)):.2e}, "
          f"rel = {abs(g0c_2 - math.sqrt(3))/math.sqrt(3)*100:.6f}%")

    # The user had g0_crit(2D) ~ 1.7324 from earlier analysis
    # Check: is g0_crit(2D) = sqrt(3)?
    # Then g_bar(2D) = sqrt(3) = (4/pi)*Q2 => Q2 = pi*sqrt(3)/4

    print(f"\n  If g0_crit(2D) = sqrt(3):")
    print(f"    Q_2 = pi*sqrt(3)/4 = {math.pi*math.sqrt(3)/4:.10f}")

    # Exhaustive algebraic search
    print(f"\n  Exhaustive search for Q_2 = {Q2:.10f}:")
    all_candidates = []

    # pi^a * sqrt(b) * c/d combinations
    for num in range(1, 20):
        for den in range(1, 20):
            val = num/den
            rel = abs(Q2 - val) / Q2 * 100
            if rel < 0.5:
                all_candidates.append((rel, f"{num}/{den}", val))

    for b in [2, 3, 5, 6, 7]:
        for num in range(1, 10):
            for den in range(1, 10):
                val = math.sqrt(b) * num / den
                rel = abs(Q2 - val) / Q2 * 100
                if rel < 0.5:
                    all_candidates.append((rel, f"sqrt({b})*{num}/{den}", val))

                val2 = math.pi * num / (den * math.sqrt(b))
                rel2 = abs(Q2 - val2) / Q2 * 100
                if rel2 < 0.5:
                    all_candidates.append((rel2, f"pi*{num}/({den}*sqrt({b}))", val2))

    for b in [2, 3, 5]:
        val = math.pi * math.sqrt(b) / 4
        rel = abs(Q2 - val) / Q2 * 100
        if rel < 2.0:
            all_candidates.append((rel, f"pi*sqrt({b})/4", val))

    all_candidates.sort()
    if all_candidates:
        for rel, name, val in all_candidates[:10]:
            print(f"    Q_2 ~ {name:30s} = {val:.10f}  ({rel:.6f}%)")
    else:
        print(f"    No algebraic match found within 0.5%")


# ================================================================
# SECTION 6: g0_crit = (4/pi)*sqrt(d) vs exact
# ================================================================
print(f"\n{'='*70}")
print("  6. REFINED HYPOTHESIS")
print("="*70)

# Check if g0_crit(d) itself has a pattern
# d=1: 4/3, d=2: ~sqrt(3)=1.7321, d=3: ~2.206
# If g0_crit(d) = sqrt(3) for d=2, then:
#   g_bar(1) = 4/3 = 1.3333
#   g_bar(2) = sqrt(3) = 1.7321
#   g_bar(3) = 2.206

# Ratios:
if len(results) >= 3:
    print(f"\n  Ratios between consecutive barriers:")
    ds = sorted(results.keys())
    for i in range(len(ds)-1):
        d1, d2 = ds[i], ds[i+1]
        r = results[d2][0] / results[d1][0]
        print(f"    g0_crit({d2})/g0_crit({d1}) = {r:.8f}")

    # Check: g0_crit(3)/g0_crit(2)
    if 2 in results and 3 in results:
        r32 = results[3][0] / results[2][0]
        print(f"\n    g0_crit(3)/g0_crit(2) = {r32:.8f}")
        print(f"    4/pi = {4/math.pi:.8f}")
        print(f"    sqrt(phi) = {math.sqrt(PHI):.8f}")
        print(f"    phi/sqrt(2) = {PHI/math.sqrt(2):.8f}")

    if 1 in results and 2 in results:
        r21 = results[2][0] / results[1][0]
        print(f"\n    g0_crit(2)/g0_crit(1) = {r21:.8f}")
        print(f"    3*sqrt(3)/4 = {3*math.sqrt(3)/4:.8f}")


# ================================================================
# SECTION 7: Summary
# ================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)

print(f"\n  User's hypothesis: g_bar(d) = (4/pi) * Q_d")
print(f"\n  {'d':>3s}  {'g0_crit':>14s}  {'Q_d':>14s}  {'candidate':>20s}  {'diff%':>8s}")
print(f"  {'-'*3}  {'-'*14}  {'-'*14}  {'-'*20}  {'-'*8}")

for d in sorted(results.keys()):
    g0c, Q = results[d]
    if d == 1:
        cand = "pi/3"
        cval = math.pi/3
    elif d == 2 and abs(g0c - math.sqrt(3)) / math.sqrt(3) < 0.01:
        cand = "pi*sqrt(3)/4"
        cval = math.pi * math.sqrt(3) / 4
    elif d == 3:
        cand = "sqrt(3)"
        cval = math.sqrt(3)
    else:
        cand = "?"
        cval = Q
    rel = abs(Q - cval) / cval * 100 if cval != Q else 0
    print(f"  {d:3d}  {g0c:14.10f}  {Q:14.10f}  {cand:>20s}  {rel:8.4f}")

print(f"""

  KEY FINDINGS:
  d=1: g0_crit = 4/3 (EXACT), Q_1 = pi/3
  d=2: g0_crit ~ sqrt(3)? (need to verify precision)
  d=3: g0_crit = 2.206, Q_3 ~ sqrt(3) (0.03%)

  If g0_crit(2D) = sqrt(3):
    Then g_bar(d) is NOT (4/pi)*Q_d with universal Q pattern
    but rather g0_crit itself follows a geometric sequence!
""")
