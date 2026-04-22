#!/usr/bin/env python3
"""
gcrit_energy_proof_v47b.py -- Energy dissipation approach to g0_crit for general d.

STRATEGY:
In canonical variable f = g^(alpha+1), the ODE becomes:
  f'' + (d-1)/r * f' = S(f)
where S(f) = (alpha+1) * f^((alpha+2)/(alpha+1)) * (1 - f^(1/(alpha+1)))

Define potential V(f) such that S(f) = -dV/df:
  V(f) = -(alpha+1)^2/(alpha+2) * f^((alpha+2)/(alpha+1))
         +(alpha+1)^2/(2*alpha+3) * f^((2*alpha+3)/(alpha+1))
  [Since S = (a+1)*f^p*(1-f^q) with p=(a+2)/(a+1), q=1/(a+1)]

Energy: E = f'^2/2 + V(f)
  dE/dr = -(d-1)/r * f'^2   (friction)

For d=1: no friction => E conserved => V(f0) = V(0) = 0 => gives g0_crit.
For d>1: V(f0) = V(0) + integral of friction = integral of (d-1)/r * f'^2 dr

QUESTION: Does the friction integral have a clean form that gives the formula?
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp

# ================================================================
# Section 1: Define the canonical ODE in f-variable
# ================================================================

def make_f_solver(alpha, d, r_max=500):
    """Solve the ODE in canonical variable f = g^(alpha+1)."""
    a1 = alpha + 1  # shorthand

    def solver(f0):
        def rhs(r, y):
            f, fp = y
            f = max(f, 1e-15)
            # S(f) = (a+1) * f^((a+2)/(a+1)) * (1 - f^(1/(a+1)))
            f_pa = f ** ((alpha + 2) / a1)
            f_q = f ** (1.0 / a1)
            source = a1 * f_pa * (1.0 - f_q)
            if d == 1:
                return [fp, source]
            if r < 1e-10:
                return [fp, source / float(d)]
            return [fp, source - float(d - 1) * fp / r]

        sol = solve_ivp(rhs, (1e-6, r_max), [f0, 0.0],
                        rtol=1e-12, atol=1e-14, max_step=0.01,
                        dense_output=True)
        return sol.t, sol.y[0], sol.y[1]
    return solver


def V_potential(f, alpha):
    """Potential V(f) such that S(f) = -dV/df."""
    a1 = alpha + 1
    # Integrate S(f) = (a+1)*f^((a+2)/(a+1))*(1 - f^(1/(a+1)))
    # S = (a+1)*f^p - (a+1)*f^(p+q) where p=(a+2)/(a+1), q=1/(a+1)
    # p+q = (a+3)/(a+1)
    # V = -integral of S df = -(a+1)*f^(p+1)/(p+1) + (a+1)*f^(p+q+1)/(p+q+1)
    # p+1 = (2a+3)/(a+1), p+q+1 = (2a+4)/(a+1)

    exp1 = (2*alpha + 3) / a1
    exp2 = (2*alpha + 4) / a1

    V = -(a1**2 / (2*alpha + 3)) * f**exp1 + (a1**2 / (2*alpha + 4)) * f**exp2
    return V


def g0_crit_formula(alpha, d):
    return (2*alpha + 4) / (2*alpha + 4 - d)


# ================================================================
# Section 2: Compute friction integral for critical solutions
# ================================================================
print("=" * 70)
print("ENERGY DISSIPATION ANALYSIS FOR g0_crit")
print("=" * 70)

print("\nSection 1: Verify V(f) definition")
print("-" * 50)

# Check V(0) = 0
for alpha in [1.0, 2.0, 3.0]:
    V0 = V_potential(1e-15, alpha)
    print(f"  alpha={alpha:.1f}: V(0) = {V0:.2e}")

# Check -dV/df = S(f) numerically
print("\n  Verify -dV/df = S(f):")
for alpha in [1.0, 2.0, 3.0]:
    a1 = alpha + 1
    for f_test in [0.5, 1.0, 2.0, 3.0]:
        df = 1e-8
        dVdf = (V_potential(f_test + df, alpha) - V_potential(f_test - df, alpha)) / (2*df)
        S = a1 * f_test**((alpha+2)/a1) * (1 - f_test**(1/a1))
        print(f"  alpha={alpha:.0f}, f={f_test:.1f}: S={S:.8f}, -dV/df={-dVdf:.8f}, diff={abs(S+dVdf):.2e}")


# ================================================================
# Section 3: For each (alpha, d), find critical solution and compute
#            friction integral vs V(f0)
# ================================================================
print("\n\nSection 2: Friction integral for critical solutions")
print("=" * 70)

def find_critical_and_friction(alpha, d, tol=1e-8):
    """Find g0_crit numerically and compute friction integral."""
    a1 = alpha + 1
    gc = g0_crit_formula(alpha, d)
    fc = gc ** a1  # f0_crit

    solver = make_f_solver(alpha, d, r_max=500)

    # Bisect to find critical f0
    def converges(f0):
        r, f, fp = solver(f0)
        return r[-1] > 400 and np.min(f) > 0.01

    f_lo = fc * 0.8
    f_hi = fc * 1.2

    # Ensure brackets
    if not converges(f_lo):
        f_lo = 0.5
    while converges(f_hi) and f_hi < fc * 3:
        f_hi *= 1.5

    for _ in range(80):
        f_mid = (f_lo + f_hi) / 2
        if f_hi - f_lo < tol * fc:
            break
        if converges(f_mid):
            f_lo = f_mid
        else:
            f_hi = f_mid

    f0_num = (f_lo + f_hi) / 2
    g0_num = f0_num ** (1.0 / a1)

    # Now compute energy quantities for the critical solution
    # Use f0 slightly below critical so it oscillates (doesn't collapse)
    f0_use = f_lo  # just below critical
    r, f, fp = solver(f0_use)

    # V at f0
    V_f0 = V_potential(f0_use, alpha)
    V_0 = 0.0  # V(0) = 0

    # Friction integral: integral of (d-1)/r * fp^2 dr
    # Use trapezoidal rule
    if d == 1:
        friction = 0.0
    else:
        # Avoid r=0 singularity
        mask = r > 0.01
        r_m = r[mask]
        fp_m = fp[mask]
        integrand = (d - 1) / r_m * fp_m**2
        friction = trapz(integrand, r_m)

    # Energy at start: E0 = fp(0)^2/2 + V(f0) = V(f0) since fp(0)=0
    E0 = V_f0

    # Energy at r->inf: E_inf = V(0) = 0 (for oscillating solution around f=1)
    # Actually for f->1 (vacuum), V(1) = ?
    V_1 = V_potential(1.0, alpha)

    # For critical solution: f starts at f0, goes toward 0 or oscillates
    # The critical solution is the boundary -- it approaches 0 asymptotically
    # So E_inf -> V(0) = 0

    # Energy balance: E0 = E_inf + friction_dissipated
    # V(f0) = 0 + friction (if critical sol goes to f=0)
    # But V(f0) < 0 for our potential, and friction > 0
    # So: V(f0) + friction_dissipated = 0  ???
    # Actually: dE/dr = -(d-1)/r * fp^2 <= 0
    # E decreases. E0 = V(f0). E_inf = 0.
    # So V(f0) = friction_dissipated (friction is the energy lost)
    # friction_dissipated = -integral of (d-1)/r * fp^2 dr
    # V(f0) = -friction_integral
    # friction_integral = -V(f0)

    return {
        'alpha': alpha, 'd': d,
        'g0_crit_formula': gc,
        'g0_crit_num': g0_num,
        'f0_crit': f0_num,
        'V_f0': V_f0,
        'V_1': V_1,
        'friction': friction,
        'V_f0_plus_friction': V_f0 + friction,  # should be ~0 if E_inf=0
        'ratio_fric_V': friction / abs(V_f0) if abs(V_f0) > 1e-15 else 0,
    }


# Test cases
test_cases = [
    (1.0, 1), (2.0, 1), (3.0, 1),  # d=1: friction=0
    (1.0, 2), (2.0, 2), (3.0, 2),
    (1.0, 3), (1.5, 3), (2.0, 3), (2.5, 3), (3.0, 3),
    (2.0, 4), (2.0, 5),
]

print(f"\n  {'alpha':>5s} {'d':>3s} {'g0c_form':>10s} {'g0c_num':>10s} {'V(f0)':>12s} {'friction':>12s} {'V+fric':>12s} {'fric/|V|':>10s}")

results = []
for alpha, d in test_cases:
    res = find_critical_and_friction(alpha, d)
    results.append(res)
    print(f"  {alpha:5.1f} {d:3d} {res['g0_crit_formula']:10.6f} {res['g0_crit_num']:10.6f} "
          f"{res['V_f0']:12.6f} {res['friction']:12.6f} {res['V_f0_plus_friction']:12.6f} {res['ratio_fric_V']:10.4f}")


# ================================================================
# Section 4: Analyze the V(f0) formula at g0_crit
# ================================================================
print("\n\nSection 3: V(f0_crit) analytical formula")
print("=" * 70)

print("\n  V(f0) = -(a+1)^2/(2a+3) * f0^((2a+3)/(a+1)) + (a+1)^2/(2a+4) * f0^((2a+4)/(a+1))")
print("\n  At g0_crit = (2a+4)/(2a+4-d), f0_crit = g0_crit^(a+1):")

for alpha in [1.0, 2.0, 3.0]:
    a1 = alpha + 1
    print(f"\n  alpha = {alpha:.0f}:")
    for d in [1, 2, 3, 4]:
        if d >= 2*alpha + 4:
            continue
        gc = (2*alpha + 4) / (2*alpha + 4 - d)
        fc = gc ** a1
        V = V_potential(fc, alpha)

        # For d=1: V(f0_crit) should be 0
        # For d>1: V(f0_crit) should be negative (energy lost to friction)
        print(f"    d={d}: g0c={gc:.6f}, f0c={fc:.6f}, V(f0c)={V:.10f}")


# ================================================================
# Section 5: Key question -- does V(f0_crit) = 0 for all d?
# ================================================================
print("\n\nSection 4: Does V(f0_crit) = 0 for general d?")
print("=" * 70)

print("\n  V(f0) = 0 condition (the d=1 proof condition):")
print("  -(a+1)^2/(2a+3) * f0^((2a+3)/(a+1)) + (a+1)^2/(2a+4) * f0^((2a+4)/(a+1)) = 0")
print("  => f0^(1/(a+1)) = (2a+4)/(2a+3)")
print("  => g0 = (2a+4)/(2a+3)  -- this is the d=1 formula!")
print()
print("  For d>1, V(f0_crit) != 0. The friction compensates.")
print("  The question is: does the friction integral have a closed form?")

# Let's check what ratio friction/|V(f0)| is
print("\n  Friction ratio analysis:")
for res in results:
    if res['d'] > 1 and abs(res['V_f0']) > 1e-10:
        d = res['d']
        alpha = res['alpha']
        # The d=1 critical g0 is (2a+4)/(2a+3)
        # The general d critical g0 is (2a+4)/(2a+4-d)
        # V(f0_crit) for g0=(2a+4)/(2a+4-d) is nonzero
        gc = res['g0_crit_formula']
        gc_d1 = (2*alpha + 4) / (2*alpha + 3)

        print(f"  alpha={alpha:.1f}, d={d}: fric/|V|={res['ratio_fric_V']:.6f}, "
              f"g0c/g0c_d1={gc/gc_d1:.6f}")


# ================================================================
# Section 6: Alternative -- work directly in g-variable
# ================================================================
print("\n\nSection 5: Direct g-variable energy analysis")
print("=" * 70)

print("""
  Alternative: work in original g-variable.
  ODE: g'' + (d-1)/r*g' + (alpha/g)*g'^2 = g^2(1-g)

  Multiply by g^(2*alpha)*g':
  d/dr[g^(2a)*g'^2/2] + (d-1)/r * g^(2a)*g'^2 = g^(2a+2)*(1-g)*g'

  Define: H = g^(2a)*g'^2/2 (kinetic)
          W(g) = integral of g^(2a+2)*(1-g) dg
               = g^(2a+3)/(2a+3) - g^(2a+4)/(2a+4)

  Energy: dH/dr = dW/dg * g' - (d-1)/r * g^(2a)*g'^2
          d/dr[H - W(g)] = -(d-1)/r * g^(2a)*g'^2

  Define E_g = H - W(g) = g^(2a)*g'^2/2 - W(g)
  At r=0: E_g(0) = -W(g0) = -(g0^(2a+3)/(2a+3) - g0^(2a+4)/(2a+4))

  Critical solution: g -> 0 as r -> inf (for collapse boundary)
  E_g(inf) = 0 - W(0) = 0

  So: -W(g0) + integral of (d-1)/r * g^(2a)*g'^2 dr = 0
  i.e.: W(g0) = integral of (d-1)/r * g^(2a)*g'^2 dr
""")

# Compute this in g-variable
def find_critical_g_energy(alpha, d, r_max=500, tol=1e-8):
    """Find critical g0 and compute energy balance in g-variable."""
    def solve_g(g0):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            source = g**2 * (1.0 - g)
            cross = (alpha / g) * gp**2
            if d == 1:
                return [gp, source - cross]
            if r < 1e-10:
                return [gp, (source - cross) / float(d)]
            return [gp, source - cross - float(d - 1) * gp / r]
        sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                        rtol=1e-12, atol=1e-14, max_step=0.01)
        return sol.t, sol.y[0], sol.y[1]

    gc = g0_crit_formula(alpha, d)

    def converges(g0):
        r, g, gp = solve_g(g0)
        return r[-1] > 400 and np.min(g) > 0.01

    g_lo = max(1.001, gc - 0.3)
    g_hi = gc + 0.3
    if not converges(g_lo):
        g_lo = 1.001
    while converges(g_hi):
        g_hi += 0.5

    for _ in range(80):
        g_mid = (g_lo + g_hi) / 2
        if g_hi - g_lo < tol:
            break
        if converges(g_mid):
            g_lo = g_mid
        else:
            g_hi = g_mid

    g0_num = g_lo  # just below critical
    r, g, gp = solve_g(g0_num)

    # W(g0)
    W_g0 = g0_num**(2*alpha+3)/(2*alpha+3) - g0_num**(2*alpha+4)/(2*alpha+4)

    # Friction integral in g-variable
    if d > 1:
        mask = r > 0.01
        r_m = r[mask]
        g_m = g[mask]
        gp_m = gp[mask]
        integrand = (d-1) / r_m * g_m**(2*alpha) * gp_m**2
        friction_g = trapz(integrand, r_m)
    else:
        friction_g = 0.0

    return {
        'g0': g0_num, 'gc_form': gc,
        'W_g0': W_g0,
        'friction_g': friction_g,
        'balance': W_g0 - friction_g,  # should be ~0
        'ratio': friction_g / W_g0 if W_g0 > 1e-15 else 0,
    }


print(f"\n  {'alpha':>5s} {'d':>3s} {'g0c':>10s} {'W(g0)':>14s} {'friction':>14s} {'W-fric':>14s} {'fric/W':>10s}")

for alpha, d in test_cases:
    res = find_critical_g_energy(alpha, d)
    print(f"  {alpha:5.1f} {d:3d} {res['gc_form']:10.6f} {res['W_g0']:14.10f} {res['friction_g']:14.10f} "
          f"{res['balance']:14.2e} {res['ratio']:10.6f}")


# ================================================================
# Section 7: Analyze the friction ratio
# ================================================================
print("\n\nSection 6: Friction ratio = fric/W(g0) patterns")
print("=" * 70)

# For d=1: ratio = 0 (no friction)
# For d=2: ratio = ?
# For d=3: ratio = ?
# Pattern?

print("\n  Fixed alpha=2, varying d:")
for d in [1, 2, 3, 4, 5, 6, 7]:
    if d >= 8:
        continue
    res = find_critical_g_energy(2.0, d)
    gc = res['gc_form']
    W = res['W_g0']
    F = res['friction_g']

    # W(g0_crit) for d=1 is 0 by definition
    # For general d, W(g0_crit) = g0c^(2a+3)/(2a+3) - g0c^(2a+4)/(2a+4)
    # = g0c^(2a+3) * [1/(2a+3) - g0c/(2a+4)]
    # At g0c=(2a+4)/(2a+4-d):
    # = gc^(2a+3) * [1/(2a+3) - (2a+4)/((2a+4-d)*(2a+4))]
    # = gc^(2a+3) * [(2a+4-d)*(2a+4) - (2a+3)*(2a+4)] / [(2a+3)*(2a+4)*(2a+4-d)]
    # Wait, let me recalculate...
    # = gc^(2a+3) * [1/(2a+3) - gc/(2a+4)]
    # gc = N/(N-d) where N = 2a+4
    # = (N/(N-d))^(N-1) * [1/(N-1) - N/((N-d)*N)]
    # = (N/(N-d))^(N-1) * [1/(N-1) - 1/(N-d)]
    # = (N/(N-d))^(N-1) * [(N-d-N+1)/((N-1)*(N-d))]
    # = (N/(N-d))^(N-1) * [(1-d)/((N-1)*(N-d))]
    # = (N/(N-d))^(N-1) * (1-d) / ((N-1)*(N-d))

    # For d=1: W = gc^(N-1) * 0 / (...) = 0. Correct!
    # For d>1: W < 0 since (1-d) < 0. But W should be positive...
    # Let me recheck W = g0^(2a+3)/(2a+3) - g0^(2a+4)/(2a+4)
    # For g0 > 1 (which g0_crit is): g0^(2a+4) > g0^(2a+3), and 2a+4 > 2a+3
    # W = g0^(2a+3)[1/(2a+3) - g0/(2a+4)]
    # For g0 < (2a+4)/(2a+3) = gc(d=1): the term is positive
    # For g0 > gc(d=1): W < 0

    # For d>1: gc(d) > gc(d=1), so W(gc(d)) < 0
    # This means the critical solution for d>1 has W < 0
    # Energy balance: E0 = -W(g0) (positive), friction drains it to 0

    print(f"    d={d}: gc={gc:.6f}, W(gc)={W:.10f}, fric={F:.10f}, balance={W-F:.2e}")


# ================================================================
# Section 8: W(g0_crit) analytical formula
# ================================================================
print("\n\nSection 7: W(g0_crit) analytical formula")
print("=" * 70)

print("\n  W(g0) = g0^(N-1)/(N-1) - g0^N/N   where N = 2*alpha + 4")
print("  At g0 = N/(N-d):")
print("  W = (N/(N-d))^(N-1)/(N-1) - (N/(N-d))^N/N")
print("  W = (N/(N-d))^(N-1) * [1/(N-1) - N/((N-d)*N)]")
print("  W = (N/(N-d))^(N-1) * [1/(N-1) - 1/(N-d)]")
print("  W = (N/(N-d))^(N-1) * (1-d)/((N-1)*(N-d))")

print("\n  Verification:")
for alpha in [1.0, 2.0, 3.0]:
    N = 2*alpha + 4
    for d in [1, 2, 3]:
        if d >= N:
            continue
        gc = N / (N - d)
        W_num = gc**(N-1)/(N-1) - gc**N/N
        W_formula = (N/(N-d))**(N-1) * (1-d) / ((N-1)*(N-d))
        print(f"  alpha={alpha:.0f}, d={d}: W_num={W_num:.10f}, W_formula={W_formula:.10f}, diff={abs(W_num-W_formula):.2e}")


# ================================================================
# Section 9: Check if friction = W analytically
# ================================================================
print("\n\nSection 8: CRITICAL TEST -- Does friction integral = W(g0_crit)?")
print("=" * 70)

print("\n  If the energy balance holds exactly:")
print("  friction = W(g0_crit) = (N/(N-d))^(N-1) * (1-d)/((N-1)*(N-d))")
print("  For d=1: friction = 0, W = 0. Trivially true.")
print("  For d>1: friction > 0, W < 0... wait, that means -W(g0) = friction")
print()
print("  Correction: E0 = -W(g0) + 0 = -W(g0) (since g'(0)=0, H(0)=0)")
print("  E_inf = 0. dE/dr <= 0. So E0 = total_friction_loss.")
print("  -W(g0) = integral (d-1)/r * g^(2a)*g'^2 dr")
print()

print("  Recheck with corrected signs:")
for alpha, d in [(2.0, 1), (2.0, 2), (2.0, 3), (2.0, 4), (2.0, 5),
                  (1.0, 2), (1.0, 3), (3.0, 2), (3.0, 3)]:
    if d >= 2*alpha + 4:
        continue
    res = find_critical_g_energy(alpha, d)
    neg_W = -res['W_g0']
    fric = res['friction_g']
    print(f"  alpha={alpha:.0f}, d={d}: -W(g0c)={neg_W:14.10f}, friction={fric:14.10f}, "
          f"ratio={fric/neg_W:.8f}" if neg_W > 1e-15 else
          f"  alpha={alpha:.0f}, d={d}: -W(g0c)={neg_W:14.10f}, friction={fric:14.10f}")


print("\n\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
  Energy balance for critical soliton:
    -W(g0_crit) = integral_0^inf (d-1)/r * g^(2a)*g'^2 dr

  For d=1: Both sides = 0. This IS the proof (W=0 => g0_crit formula).
  For d>1: Both sides positive. The friction integral is a functional
           of the solution trajectory g(r), not just g0.

  The formula g0_crit = (2a+4)/(2a+4-d) means:
    W(g0_crit) = (N/(N-d))^(N-1) * (1-d)/((N-1)*(N-d))

  This value must EXACTLY equal the friction integral.
  If we can show this independently, we have the proof.
""")
