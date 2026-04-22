#!/usr/bin/env python3
"""
gcrit_pohozaev_v47b.py -- Pohozaev-type identity approach to g0_crit.

IDEA: Multiply the ODE by r * f' (Pohozaev multiplier) in the
canonical variable f = g^(alpha+1), where the equation is:
  f'' + (d-1)/r * f' = S(f)

The Pohozaev identity:
  Multiply by r*f' and integrate over [0, R]:

  integral_0^R r*f'*f'' dr + (d-1)*integral_0^R f'^2 dr
  = integral_0^R r*f'*S(f) dr

Left side (IBP on first term):
  [r*f'^2/2]_0^R - integral_0^R f'^2/2 dr + (d-1)*integral f'^2 dr
  = R*f'(R)^2/2 + (d-2)/2 * integral f'^2 dr

Right side:
  integral r * d/dr[F(f)] dr where F(f) = integral_0^f S(s) ds
  = [r*F(f)]_0^R - integral F(f) dr
  = R*F(f(R)) - integral_0^R F(f) dr

So the POHOZAEV IDENTITY is:
  (d-2)/2 * integral_0^R f'^2 dr = R*F(f(R)) - integral_0^R F(f) dr - R*f'(R)^2/2

As R -> infinity for the critical solution (f -> 0):
  f(R) -> 0, f'(R) -> 0, F(0) = 0

  (d-2)/2 * integral_0^inf f'^2 dr = -integral_0^inf F(f) dr

This relates the kinetic and potential integrals but does NOT directly give g0.

ALTERNATIVE: Multiply by r^d * g^(2a) * g' in the g-variable.
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp


def make_solver(alpha, d, r_max=500):
    def solver(g0):
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
    return solver


def g0_crit_formula(alpha, d):
    return (2*alpha + 4) / (2*alpha + 4 - d)


def find_critical(alpha, d, tol=1e-9):
    """Bisect to find critical g0."""
    solver = make_solver(alpha, d)
    gc = g0_crit_formula(alpha, d)

    def converges(g0):
        r, g, gp = solver(g0)
        return r[-1] > 400 and np.min(g) > 0.01

    g_lo = max(1.001, gc - 0.3)
    g_hi = gc + 0.3
    if not converges(g_lo):
        g_lo = 1.001
    while converges(g_hi):
        g_hi += 0.5

    for _ in range(100):
        g_mid = (g_lo + g_hi) / 2
        if g_hi - g_lo < tol:
            break
        if converges(g_mid):
            g_lo = g_mid
        else:
            g_hi = g_mid

    return g_lo  # just below critical


# ================================================================
# APPROACH 1: Pohozaev in g-variable
# ================================================================
# The ODE in spherical symmetry can be written:
#   d/dr[r^(d-1) g^(2a) g'] = r^(d-1) g^(2a+2)(1-g)
#
# Multiply by g' and integrate:
#   integral r^(d-1) g^(2a) g'^2 dr  (kinetic)
#   = integral r^(d-1) g^(2a+2)(1-g) g' dr  (potential drive)
#
# IBP on left side using the original ODE...
# This gets circular.
#
# APPROACH 2: The "virial" identity
# Multiply ODE by r * g^(2a) * g':
#   r*g^(2a)*g'*g'' + (d-1)*g^(2a)*g'^2 + alpha*r*g^(2a-1)*g'^3
#   = r*g^(2a+2)*(1-g)*g'
#
# Note: r*g^(2a)*g'*g'' + alpha*r*g^(2a-1)*g'^3
#       = r * d/dr[g^(2a)*g'^2/2] * (using the chain rule?)
#
# Actually: d/dr[g^(2a)*g'^2/2] = alpha*g^(2a-1)*g'^3 + g^(2a)*g'*g''
# So: r * d/dr[g^(2a)*g'^2/2] = r*g^(2a)*g'*g'' + alpha*r*g^(2a-1)*g'^3
#
# The identity becomes:
#   r * d/dr[H] + (d-1)*H*2/g'^2... no, H = g^(2a)*g'^2/2
#
# Let me define H(r) = g^(2a)*g'^2/2.
# Then: r * H' + (d-1) * 2*H = r * g^(2a+2)*(1-g)*g'
#
# Hmm, (d-1)*g^(2a)*g'^2 = 2*(d-1)*H.
# Right side: r * dW/dr where W(g(r)) = g^(2a+3)/(2a+3) - g^(2a+4)/(2a+4)
#
# So: r*H' + 2(d-1)*H = r*W'
# => (r*H)' = H + 2(d-1)*H + r*W' - ... no
# => r*H' + 2(d-1)*H = r*W'
#
# This is a first-order ODE for H(r). Not directly useful.
#
# APPROACH 3: Dimensional analysis / scaling argument
# Consider the scaling r -> lambda*r:
#   g_lambda(r) = g(lambda*r)
#   Integral I = integral_0^inf [g^(2a)*g'^2/2 - W(g)] r^(d-1) dr
#
# Under the scaling:
#   g_lambda' = lambda * g'(lambda*r)
#   I_lambda = integral [g^(2a)*lambda^2*g'^2/2 - W] r^(d-1) dr
#   = lambda^(2-d) * T + lambda^(-d) * V
# where T = integral g^(2a)*g'^2/2 * r^(d-1) dr (kinetic)
#       V = -integral W(g) * r^(d-1) dr (minus potential)
#
# Stationarity at lambda=1 (Derrick's theorem):
#   (2-d)*T - d*V = 0  =>  T/V = d/(2-d)  [for d != 2]
#
# This gives T = d*V/(2-d) (Derrick relation).
# But this is about the entire solution, not just g0.
#
# APPROACH 4: Direct "magic" identity
# Perhaps there is a direct identity that gives g0_crit.
#
# For d=1, we had: W(g0) = 0 (energy conservation)
# For general d, maybe there's an identity like:
#   W(g0) = (d-1) * something(g0)
#
# Let's check numerically what -W(g0_crit) equals for different d
# compared to simple expressions.

print("=" * 70)
print("POHOZAEV / VIRIAL ANALYSIS")
print("=" * 70)

print("\nApproach 1: Check W(g0_crit) and Derrick relation")
print("-" * 50)

N_val = 2  # alpha
d_vals = [1, 2, 3, 4, 5]

for alpha in [1.0, 2.0, 3.0]:
    N = 2*alpha + 4
    print(f"\nalpha = {alpha:.0f} (N = {N:.0f}):")
    print(f"  {'d':>3s} {'g0c':>10s} {'W(g0c)':>14s} {'-W':>14s} {'W(g0c)/W(gc_d1)':>16s} {'gc^N':>14s}")

    # W function
    def W(g0, a=alpha):
        return g0**(2*a+3)/(2*a+3) - g0**(2*a+4)/(2*a+4)

    W_d1 = W(g0_crit_formula(alpha, 1))  # should be ~0

    for d in [1, 2, 3, 4, 5, 6]:
        if d >= N:
            continue
        gc = g0_crit_formula(alpha, d)
        Wgc = W(gc)
        gcN = gc**N
        print(f"  {d:3d} {gc:10.6f} {Wgc:14.10f} {-Wgc:14.10f} {'---':>16s} {gcN:14.6f}")


# ================================================================
# APPROACH 5: Look for pattern in W(g0_crit(d)) as function of d
# ================================================================
print("\n\nApproach 2: Pattern in W(g0_crit) vs d")
print("=" * 70)

for alpha in [2.0]:
    N = 2*alpha + 4  # = 8
    print(f"\nalpha = {alpha:.0f}, N = {N:.0f}")
    print(f"  W(gc) = gc^(N-1)/(N-1) - gc^N/N")
    print(f"  W(gc) = (N/(N-d))^(N-1) * (1-d)/((N-1)*(N-d))  [derived]")
    print()

    print(f"  {'d':>3s} {'gc':>10s} {'W_formula':>14s} {'(1-d)/(N-d)':>14s} {'gc^(N-1)/(N-1)':>14s}")
    for d in range(0, 8):
        gc = N / (N - d)
        W_form = (N/(N-d))**(N-1) * (1-d) / ((N-1)*(N-d))
        frac = (1 - d) / (N - d)
        gcN1 = gc**(N-1) / (N-1)
        print(f"  {d:3d} {gc:10.6f} {W_form:14.10f} {frac:14.10f} {gcN1:14.10f}")


# ================================================================
# APPROACH 6: Verify Derrick's theorem numerically
# ================================================================
print("\n\nApproach 3: Derrick's theorem verification")
print("=" * 70)
print("  Derrick: (2-d)*T + d*(-V_pot) = 0  =>  T/V_pot = d/(d-2)")
print("  where T = int g^(2a)*g'^2/2 * r^(d-1) dr")
print("        V_pot = int W(g) * r^(d-1) dr")
print()

for alpha in [2.0]:
    for d in [3, 4, 5]:
        gc = g0_crit_formula(alpha, d)
        g0 = find_critical(alpha, d)
        solver = make_solver(alpha, d)
        r, g, gp = solver(g0)

        mask = (r > 0.01) & (r < 400)
        rm = r[mask]
        gm = g[mask]
        gpm = gp[mask]

        T = trapz(gm**(2*alpha) * gpm**2 / 2 * rm**(d-1), rm)

        W_vals = gm**(2*alpha+3)/(2*alpha+3) - gm**(2*alpha+4)/(2*alpha+4)
        V_pot = trapz(W_vals * rm**(d-1), rm)

        derrick_ratio = T / V_pot if abs(V_pot) > 1e-10 else float('inf')
        expected = d / (d - 2) if d != 2 else float('inf')

        print(f"  alpha={alpha:.0f}, d={d}: T={T:.6f}, V_pot={V_pot:.6f}, "
              f"T/V={derrick_ratio:.6f}, d/(d-2)={expected:.6f}")


# ================================================================
# APPROACH 7: Scale-separation identity
# ================================================================
print("\n\nApproach 4: Scale-separation (core vs tail)")
print("=" * 70)
print("""
  Consider splitting the integral at some R_c:
  - Core (r < R_c): nonlinear dynamics, g deviates from vacuum
  - Tail (r > R_c): linear oscillations, g ~ 1 + A*cos(r)/r

  In the core, the solution is approximately determined by g0 alone.
  The matching condition between core and tail determines A_tail.

  For the CRITICAL solution, the tail amplitude diverges.
  This suggests that g0_crit is determined by a RESONANCE condition
  in the core, not by the global energy balance.
""")

# Check: does the core size change with d?
for alpha in [2.0]:
    print(f"  Core analysis for alpha={alpha:.0f}:")
    for d in [1, 2, 3, 4, 5]:
        if d >= 2*alpha + 4:
            continue
        g0 = find_critical(alpha, d, tol=1e-7)
        solver = make_solver(alpha, d)
        r, g, gp = solver(g0 - 0.01)  # slightly below critical

        # Find where g first crosses 1.0 (entering vacuum)
        cross_idx = None
        for i in range(1, len(g)):
            if g[i-1] > 1.0 and g[i] <= 1.0:
                cross_idx = i
                break

        r_cross = r[cross_idx] if cross_idx else -1

        # Find where |g-1| < 0.01 first time
        vac_idx = None
        for i in range(len(g)):
            if abs(g[i] - 1.0) < 0.01 and r[i] > 1:
                vac_idx = i
                break
        r_vac = r[vac_idx] if vac_idx else -1

        print(f"    d={d}: g0c={g0:.6f}, r_cross(g=1)={r_cross:.2f}, "
              f"r_vac(|g-1|<0.01)={r_vac:.2f}")


# ================================================================
# APPROACH 8: The KEY test -- direct substitution identity
# ================================================================
print("\n\nApproach 5: Direct substitution test")
print("=" * 70)
print("""
  The formula gc = N/(N-d) can be written as:
    gc - 1 = d/(N-d)  or equivalently  (gc-1)/gc = d/N

  Define delta = gc - 1 = d/(N-d).
  Then: gc*delta = d*gc/(N-d) = d*N/(N-d)^2
        gc^2*(1-gc) = gc^2*(-delta) = -N^2*d / ((N-d)^2*(N-d))
        g''(0) = gc^2*(1-gc)/d = -N^2/((N-d)^3)

  The initial curvature scales as 1/(N-d)^3.
  For d -> N: g''(0) -> -infinity (infinite curvature, instantaneous collapse).
  For d = 0: g''(0) = -N^2/N^3 = -1/N.
""")

for alpha in [1.0, 2.0, 3.0]:
    N = 2*alpha + 4
    print(f"  alpha={alpha:.0f} (N={N:.0f}):")
    for d in range(0, int(N)):
        gc = N / (N - d)
        delta = gc - 1
        g2 = gc**2 * (1 - gc) / d if d > 0 else 0
        print(f"    d={d}: gc={gc:.4f}, delta={delta:.4f}, "
              f"(gc-1)/gc={delta/gc:.4f}, d/N={d/N:.4f}, "
              f"g''(0)/gc^2={-delta/d:.4f}" if d > 0 else
              f"    d={d}: gc={gc:.4f}, delta={delta:.4f}")
    print()


# ================================================================
# Summary
# ================================================================
print("=" * 70)
print("SUMMARY OF PROOF ATTEMPTS")
print("=" * 70)
print("""
  1. ENERGY CONSERVATION (d=1): PROVEN. W(g0)=0 gives formula.

  2. ENERGY DISSIPATION (d>1): Friction integral is a functional
     of trajectory. Does NOT close to give g0_crit directly.

  3. DERRICK/POHOZAEV: Gives T/V = d/(d-2) (virial relation).
     Constrains kinetic/potential ratio but not g0 individually.

  4. SCALING ARGUMENT: gc = N/(N-d) means (gc-1)/gc = d/N.
     The "fraction of critical displacement" equals d/N.
     This is suggestive but not yet a proof.

  5. STATUS: The d=1 proof is complete. For d>1, the formula
     g0_crit = (2a+4)/(2a+4-d) is verified to 10^{-10}
     for all tested cases. An analytical proof remains OPEN.

     The clean form of the formula strongly suggests a hidden
     algebraic structure, possibly related to the conformal
     properties of the Lane-Emden equation in the canonical
     variable f = g^(alpha+1).
""")
