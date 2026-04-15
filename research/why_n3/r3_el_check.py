#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_el_check.py — Check the E-L equations for general K(g)

KEY QUESTION: Does g₀_crit depend on the choice of K(g)?

The substrate ODE: g'' + (1/g)g'² + (2/r)g' = 1 - g

IS THIS the E-L equation of L = K(g)·[g'²/2 + (1-g)²/2] with K=g²?

Let me derive carefully.
"""

import numpy as np
from scipy.integrate import solve_ivp
import math

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941

# ================================================================
# E-L DERIVATION
# ================================================================
print("=" * 70)
print("  E-L DERIVATION FOR L = K(g)·[g'²/2 + (1-g)²/2]·r²")
print("=" * 70)

print("""
  L = K(g)·[g'²/2 + (1-g)²/2]·r²

  ∂L/∂g' = K·g'·r²
  d/dr(∂L/∂g') = K'·g'²·r² + K·g''·r² + 2K·g'·r

  ∂L/∂g = K'·[g'²/2 + (1-g)²/2]·r² - K·(1-g)·r²

  E-L: ∂L/∂g - d/dr(∂L/∂g') = 0

  K'g'²/2 + K'(1-g)²/2 - K(1-g) - K'g'² - Kg'' - 2Kg'/r = 0

  Kg'' + K'g'²/2 + 2Kg'/r = K'(1-g)²/2 - K(1-g)

  Divide by K:
  g'' + [K'/(2K)]g'² + (2/r)g' = [K'/(2K)](1-g)² - (1-g)

  For K = g^{2α}: K'/K = 2α/g, K'/(2K) = α/g

  g'' + (α/g)g'² + (2/r)g' = (α/g)(1-g)² - (1-g)
                             = (1-g)·[α(1-g)/g - 1]
                             = (1-g)·[α - αg - g]/g
                             = (1-g)·[α - (α+1)g]/g

  CHECK α=1:
  RHS = (1-g)(1 - 2g)/g = (1 - 3g + 2g²)/g
  This is NOT the substrate RHS = (1-g)!

  So substrate ODE ≠ E-L from K=g² action.
  The substrate ODE comes from a DIFFERENT action.
""")

# What action gives the substrate ODE?
print("What action gives g'' + (1/g)g'² + (2/r)g' = 1-g ?")
print("""
  Rewrite: g·g'' + g'² + (2g/r)g' = g(1-g)
  Note: g·g'' + g'² = d/dr(g·g') = (gg')'
  ... no: (gg')' = g'² + gg''  ✓

  So: (gg')' + (2g/r)g' = g(1-g)
  And (gg')' + 2gg'/r = g(1-g)
  d/dr(r²·g·g') = r²·g(1-g)  ... no:
  d/dr(r²·g·g') = r²(g'² + gg'') + 2r·g·g'
                = r²[(gg')' + 2gg'/r] - r²g'² + r²g'²  ... same thing

  The equation can come from:
  L_sub = g·g'²/2 + V_sub(g)  where V_sub'' = ...

  Actually, let me try L_sub = (1/2)g(g')² + g²/2 - g³/3
  Then ∂L/∂g' = g·g'
  d/dr(∂L/∂g') = g'² + g·g''

  ∂L/∂g = g'²/2 + g - g²

  E-L: g'²/2 + g - g² - g'² - g'' = 0 ... wrong sign/structure

  Let me use proper Lagrangian with r² factor:
  L = [g(g')²/2 + U(g)]·r²

  ∂L/∂g' = g·g'·r²
  d/dr(...) = g'²·r² + g·g''·r² + 2g·g'·r

  ∂L/∂g = [g'²/2 + U'(g)]·r²

  E-L: g'²/2 + U' - g'² - gg'' - 2gg'/r = 0
  gg'' + g'²/2 + 2gg'/r = U'(g)

  Divide by g: g'' + g'²/(2g) + 2g'/r = U'(g)/g

  Compare with substrate: g'' + g'²/g + 2g'/r = 1-g

  The g'² coefficient is 1/(2g) from E-L but 1/g from substrate.
  These are DIFFERENT unless we modify the kinetic term.

  For L = [f(g)·g'²/2 + U(g)]·r²:
  E-L: f·g'' + f'g'²/2 + 2fg'/r + f'g'²/2 = U'
  Wait: ∂L/∂g' = f·g'·r², d/dr = f'g'²r² + fg''r² + 2fg'r
        ∂L/∂g = f'g'²/2·r² + U'r²
  E-L: f'g'²/2 + U' - f'g'² - fg'' - 2fg'/r = 0
  → fg'' + f'g'²/2 + 2fg'/r = U'

  Divide by f: g'' + f'/(2f)·g'² + 2g'/r = U'/f

  For substrate: f'/(2f) = 1/g → f'/f = 2/g → f = g²
  And U'/f = U'/g² = 1-g → U' = g²(1-g) = g² - g³ → U = g³/3 - g⁴/4

  So substrate ODE comes from:
  L_substrate = [g²·g'²/2 + g³/3 - g⁴/4]·r²

  NOT from g²·[g'²/2 + (1-g)²/2]!
""")

# Verify
print("=" * 70)
print("  VERIFICATION: Substrate = E-L of g²g'²/2 + g³/3 - g⁴/4")
print("=" * 70)

print("""
  L = g²·g'²/2 + g³/3 - g⁴/4  (times r²)

  f(g) = g², U(g) = g³/3 - g⁴/4
  f'/2f = 2g/(2g²) = 1/g  ✓
  U'/f = (g² - g³)/g² = 1 - g  ✓

  Substrate ODE: g'' + (1/g)g'² + (2/r)g' = 1-g  ✓✓✓

  So substrate potential is V_sub(g) = g³/3 - g⁴/4
  NOT V = g²(1-g)²/2 = g²/2 - g³ + g⁴/2!

  The canonical K=g⁴ case:
  L_can = g⁴·g'²/2 + U_can(g)

  If we want g'' + (2/g)g'² + (2/r)g' = RHS:
    f'/(2f) = 4g³/(2g⁴) = 2/g  ✓
    RHS = U_can'/g⁴

  But what IS U_can(g)?
""")

# The key question: what potential V gives the "right" ODE?
# The substrate has V_sub = g³/3 - g⁴/4
# V_sub(1) = 1/3 - 1/4 = 1/12
# V_sub'(g) = g² - g³ = g²(1-g) → zero at g=0, g=1
# V_sub(0) = 0

# For the metric singularity: g→0
# V_sub(0) = 0, V_sub(1) = 1/12
# The "energy" in the 1D problem: E = g²g'²/2 + g³/3 - g⁴/4
# At g=g₀, g'=0: E₀ = g₀³/3 - g₀⁴/4
# At g_min, g'=0: E_min = g_min³/3 - g_min⁴/4
# Conservation: E₀ = E_min (in 1D)
# g₀³/3 - g₀⁴/4 = g_min³/3 - g_min⁴/4
# This is F(g₀) = F(g_min) where F(x) = x³/3 - x⁴/4

# Wait, but earlier we had F(x) = 2x³/3 - x⁴/2
# The factor of 2 comes from the conservation law derivation.
# Let me re-derive...

# 1D ODE: g'' + (1/g)g'² = 1 - g
# Multiply by g²g':
# g²g'g'' + g(g')³ = g²(1-g)g'
# d/dr[g²g'²/2] = g²(1-g)g' = d/dr[g³/3 - g⁴/4]
# So E = g²g'²/2 - g³/3 + g⁴/4 = const

# At g₀: E₀ = -g₀³/3 + g₀⁴/4
# At g_min: E_min = -g_min³/3 + g_min⁴/4
# → g_min³/3 - g_min⁴/4 = g₀³/3 - g₀⁴/4
# → F(g_min) = F(g₀) where F(x) = x³/3 - x⁴/4
# Note: 2F(x) = 2x³/3 - x⁴/2 which matches the earlier formula! ✓

print("=" * 70)
print("  GENERAL ODE WITH CORRECT LAGRANGIAN")
print("=" * 70)

print("""
  The CORRECT Lagrangian for general α is:
  L = f(g)·g'²/2 + U(g)  with f(g) = g^{2α}

  ODE: g'' + α/g·g'² + (2/r)g' = U'(g)/g^{2α}

  For substrate (α=1): U' = g²(1-g), U = g³/3 - g⁴/4
  For α=2: we need to choose U₂(g) such that the ODE is physical.

  Key insight: the PHYSICS determines U(g), not just K(g)!

  If U = g³/3 - g⁴/4 for ALL α:
    ODE: g'' + α/g·g'² + (2/r)g' = (g² - g³)/g^{2α}
       = g^{2-2α} - g^{3-2α}

    α=1: g^0 - g^1 = 1 - g  ✓
    α=2: g^{-2} - g^{-1} = (1-g)/g²
    α=3: g^{-4} - g^{-3} = (1-g)/g⁴ · 1/g^... hmm

  WAIT! For α=2 with U = g³/3 - g⁴/4:
    RHS = (g² - g³)/g⁴ = 1/g² - 1/g = (1-g)/g²

  So the canonical ODE I had: g'' + (2/g)g'² + (2/r)g' = (1-g)/g²
  is WRONG! I had (3/g) instead of (2/g).

  The correct ODE with f=g⁴, U=g³/3-g⁴/4:
    g'' + (2/g)g'² + (2/r)g' = (1-g)/g²

  And with f=g², U=g³/3-g⁴/4:
    g'' + (1/g)g'² + (2/r)g' = 1-g

  THESE TWO SHARE THE SAME POTENTIAL U!
  The only difference is the kinetic coupling α/g.
""")

# Let me compute g₀_crit for the CORRECT canonical ODE
# g'' + (2/g)g'² + (2/r)g' = (1-g)/g²

print("=" * 70)
print("  g₀_crit FOR CORRECT E-L EQUATIONS (same U, different α)")
print("=" * 70)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}: PASS  {detail}")
    else:
        FAIL += 1
        print(f"  ✗ {name}: FAIL  {detail}")
    return condition

def solve_alpha(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    """
    ODE: g'' + (α/g)g'² + ((d-1)/r)g' = (1-g)·g^{2-2α}

    This comes from L = g^{2α}·g'²/2 + g³/3 - g⁴/4
    (same potential U for all α, different kinetic coupling)
    """
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
            # g'' + (d-1)g'' = RHS(g₀) = (1-g₀)·g₀^{2-2α}
            gpp = (1-g) * g**(2 - 2*alpha) / d
        else:
            # RHS = (g² - g³) / g^{2α} = g^{2-2α} - g^{3-2α} = (1-g)g^{2-2α}
            rhs_val = (1-g) * g**(2 - 2*alpha)
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1)/r) * gp
        return [gp, gpp]

    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)

    if sol.success:
        actual_min = np.min(sol.y[0])
        if actual_min < g_min_val[0]:
            g_min_val[0] = actual_min

    return g_min_val[0], singular[0], sol


def find_g0_crit_alpha(alpha, d=3, g0_lo=1.01, g0_hi=5.0, tol=1e-7):
    """Find g₀_crit for given α."""
    g_threshold = 0.005

    g_min, sing, sol = solve_alpha(g0_hi, alpha, d)
    is_bad = sing or g_min < g_threshold or not sol.success

    if not is_bad:
        g0_hi = 10.0
        g_min, sing, sol = solve_alpha(g0_hi, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success
        if not is_bad:
            return None

    for _ in range(70):
        g0_mid = (g0_lo + g0_hi) / 2
        g_min, sing, sol = solve_alpha(g0_mid, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success

        if is_bad:
            g0_hi = g0_mid
        else:
            g0_lo = g0_mid

        if g0_hi - g0_lo < tol:
            break

    return (g0_lo + g0_hi) / 2


# Compute g₀_crit for various α
print()
alphas = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
g0_crits_alpha = {}

for alpha in alphas:
    g0c = find_g0_crit_alpha(alpha, d=3)
    if g0c is not None:
        g0_crits_alpha[alpha] = g0c
        print(f"  α = {alpha:.1f}: g₀_crit(3D) = {g0c:.6f}")
    else:
        print(f"  α = {alpha:.1f}: NO SINGULARITY FOUND")

# 1D conservation law: does g₀_crit(1D) change with α?
print()
print("  1D case (no damping):")
for alpha in [1.0, 2.0, 3.0]:
    g0c_1d = find_g0_crit_alpha(alpha, d=1)
    if g0c_1d is not None:
        print(f"  α = {alpha:.1f}: g₀_crit(1D) = {g0c_1d:.6f}")

# Analysis: 1D conservation law for general α
print("""
  1D conservation law for general α:
  ODE: g'' + (α/g)g'² = (1-g)g^{2-2α}

  Let q = g'^2. Then dq/dg + (2α/g)q = 2(1-g)g^{2-2α}
  Integrating factor: g^{2α}
  d/dg[g^{2α}·q] = 2(1-g)·g^{2-2α}·g^{2α} = 2(1-g)g² = 2g² - 2g³
  g^{2α}·q = 2g³/3 - g⁴/2 + C
  q = g^{-2α}·[2g³/3 - g⁴/2 + C]

  At g=g₀: q=0 → C = g₀⁴/2 - 2g₀³/3
  At g_min: q=0 → same condition as before!

  F(g_min) = F(g₀) where F(x) = 2x³/3 - x⁴/2

  g_min = 0 → F(0) = 0 = F(g₀) → g₀ = 4/3

  CONCLUSION: g₀_crit(1D) = 4/3 FOR ALL α!

  The α dependence appears ONLY through the damping term!
""")

check("T1: g₀_crit(1D) independent of α",
      all(abs(find_g0_crit_alpha(a, d=1) - 4/3) < 0.01 for a in [1, 2, 3]),
      "conservation law: F(g₀) = 0 is α-independent")


# Key test: does g₀_crit(3D) depend on α?
print()
print("  KEY QUESTION: Does g₀_crit(3D) depend on α?")
if len(g0_crits_alpha) >= 2:
    vals = list(g0_crits_alpha.values())
    spread = max(vals) - min(vals)
    print(f"  Spread of g₀_crit across α: {spread:.6f}")

    if spread < 0.01:
        print("  → g₀_crit is INDEPENDENT of α! (spread < 0.01)")
        print("  → The barrier is a property of the POTENTIAL, not the kinetic term!")
    else:
        print(f"  → g₀_crit VARIES with α (spread = {spread:.4f})")
        print("  → The kinetic coupling affects the barrier")

# Generation counting for each α
print()
print("  GENERATION COUNTING for each α:")
for alpha, g0c in sorted(g0_crits_alpha.items()):
    n_max = math.log(g0c / G0_E) / math.log(PHI) if g0c > G0_E else 0
    N = max(0, math.floor(n_max) + 1)
    print(f"  α={alpha:.1f}: g₀_crit={g0c:.4f}, n_max={n_max:.3f}, N_gen={N}")


# Mass integral for α=2 (correct ODE)
print()
print("=" * 70)
print("  MASS FUNCTION FOR α=2 (CORRECT E-L)")
print("=" * 70)

if 2.0 in g0_crits_alpha:
    g0c_a2 = g0_crits_alpha[2.0]
    masses_a2 = []
    g0_vals_a2 = []

    for g0 in np.concatenate([
        np.arange(0.5, 1.0, 0.1),
        np.arange(1.0, min(g0c_a2 - 0.02, 3.0), 0.1),
        [g0c_a2 - 0.05, g0c_a2 - 0.02, g0c_a2 - 0.01]
    ]):
        if abs(g0 - 1.0) < 0.01:
            continue

        g_min, sing, sol = solve_alpha(g0, 2.0, d=3)
        if sing or not sol.success:
            continue

        r = sol.t
        g = sol.y[0]
        gp = sol.y[1]

        # Energy density for α=2: ε = g⁴g'²/2 + g³/3 - g⁴/4
        # Wait: total energy = kinetic + potential
        # E = g^{2α}·g'²/2 + U(g) where U = g³/3 - g⁴/4
        # But U can be < 0 for g < 1. Need to subtract vacuum energy.
        # U(1) = 1/3 - 1/4 = 1/12
        # ε = g⁴g'²/2 + (g³/3 - g⁴/4) - 1/12

        eps = g**4 * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
        integrand = eps * r**2
        mass = 4 * np.pi * np.trapezoid(integrand, r)

        if mass > 0:  # Only positive masses
            g0_vals_a2.append(g0)
            masses_a2.append(mass)

    print(f"\n  g₀_crit(α=2, 3D) = {g0c_a2:.6f}")
    print(f"\n  {'g₀':>8s}  {'m(g₀)':>12s}")
    for g0, m in zip(g0_vals_a2, masses_a2):
        note = ""
        if abs(g0 - G0_E) < 0.02:
            note = "~ electron"
        elif abs(g0 - PHI * G0_E) < 0.02:
            note = "~ muon"
        elif abs(g0 - 1.729) < 0.05:
            note = "~ tau(Koide)"
        print(f"  {g0:8.4f}  {m:12.4f}  {note}")

    if len(masses_a2) > 3:
        from scipy.interpolate import interp1d
        g0_arr = np.array(g0_vals_a2)
        m_arr = np.array(masses_a2)

        if g0_arr[0] < 0.87 and g0_arr[-1] > 1.41:
            m_interp = interp1d(g0_arr, m_arr, kind='cubic', fill_value='extrapolate')
            m_e = float(m_interp(0.869))
            m_mu = float(m_interp(1.407))

            if m_e > 0 and m_mu > 0:
                ratio = m_mu / m_e
                A_e = abs(0.869 - 1)
                A_mu = abs(1.407 - 1)
                alpha_eff = math.log(ratio) / (2*math.log(A_mu/A_e))

                print(f"\n  Mass ratio m_μ/m_e = {ratio:.2f}  (exp: 206.77)")
                print(f"  Effective α from integral: {alpha_eff:.4f}")

                # With this α_eff, check N=3
                A_tau_exp = A_e * (3477)**(1/(2*alpha_eff)) if alpha_eff > 0 else 999
                g0_tau = 1 + A_tau_exp
                A_4 = A_e * (58477)**(1/(2*alpha_eff)) if alpha_eff > 0 else 999
                g0_4 = 1 + A_4

                print(f"\n  With α_eff={alpha_eff:.3f}:")
                print(f"    g₀^τ = {g0_tau:.4f}  vs  g₀_crit = {g0c_a2:.4f}  "
                      f"{'✓' if g0_tau < g0c_a2 else '✗'}")
                print(f"    g₀^(4) = {g0_4:.4f}  vs  g₀_crit = {g0c_a2:.4f}  "
                      f"{'✓' if g0_4 > g0c_a2 else '✗'}")

                check("T2: α=2 correct E-L: τ below barrier",
                      g0_tau < g0c_a2,
                      f"g₀^τ = {g0_tau:.4f}")

                check("T3: α=2 correct E-L: 4th above barrier",
                      g0_4 > g0c_a2,
                      f"g₀^(4) = {g0_4:.4f}")


# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'=' * 70}")
print("  SUMMARY")
print("=" * 70)

print(f"""
  KEY FINDING: SUBSTRATE ODE COMES FROM SPECIFIC LAGRANGIAN

  L = g^{{2α}}·g'²/2 + g³/3 - g⁴/4  (times r²)

  The potential U(g) = g³/3 - g⁴/4 is FIXED (same for all α).
  Only the kinetic coupling g^{{2α}} changes with α.

  1D THEOREM: g₀_crit(1D) = 4/3 FOR ALL α
  (Conservation law is α-independent!)

  3D: g₀_crit depends on α through the damping term.
""")

if g0_crits_alpha:
    for alpha, g0c in sorted(g0_crits_alpha.items()):
        print(f"  α={alpha:.1f}: g₀_crit(3D) = {g0c:.6f}")

print(f"\n{'=' * 70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS + FAIL}")
print(f"{'=' * 70}")
