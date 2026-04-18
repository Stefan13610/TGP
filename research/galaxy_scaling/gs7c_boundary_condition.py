#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs7c_boundary_condition.py: Option A — Substrate boundary condition at c/H0.

IDEA: The TGP substrate is a finite expanding medium.
At the causal boundary (Hubble horizon c/H0), boundary conditions
affect the propagation of gravity.

This is NOT perturbative — it's a global property of the substrate.

In GR: the cosmological constant Lambda effectively modifies the
Schwarzschild metric to Schwarzschild-de Sitter, which changes
the force law at large distances. Can TGP do something similar
but STRONGER?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

print("="*78)
print("  OPTION A: Boundary condition at cosmological horizon")
print("="*78)

# ===========================================================================
# 1. THE SCHWARZSCHILD-de SITTER ANALOGY
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. Schwarzschild-de Sitter: how Lambda modifies gravity in GR")
print(f"{'='*78}")

c = 2.998e8
G = 6.674e-11
H0 = 2.20e-18
Lambda = 3 * H0**2 * 0.685
a0_obs = 1.2e-10
M_sun = 1.989e30

print(f"""
  In GR, the Schwarzschild-de Sitter metric gives:
  g_tt = 1 - 2GM/(rc^2) - Lambda*r^2/3

  The gravitational force on a test particle:
  F/m = -GM/r^2 + Lambda*c^2*r/3

  The Lambda term is a REPULSIVE force that grows with r.
  It equals the Newtonian force at:
  r_eq = (3GM/(Lambda*c^2))^(1/3)

  For Milky Way (M = 7e10 M_sun):
""")

M_MW = 7e10 * M_sun
r_eq = (3*G*M_MW / (Lambda * c**2))**(1/3)
print(f"    r_eq = (3GM/Lambda*c^2)^(1/3) = {r_eq:.3e} m = {r_eq/3.086e22:.1f} Mpc")
print(f"    This is {r_eq/3.086e19/10:.0f} times the galaxy radius (10 kpc)")

F_Newton_10kpc = G*M_MW / (10*3.086e19)**2
F_Lambda_10kpc = Lambda * c**2 * 10*3.086e19 / 3
print(f"\n    At r = 10 kpc:")
print(f"    F_Newton = {F_Newton_10kpc:.3e} m/s^2")
print(f"    F_Lambda = {F_Lambda_10kpc:.3e} m/s^2")
print(f"    F_Lambda/F_Newton = {F_Lambda_10kpc/F_Newton_10kpc:.3e}")
print(f"    → Lambda effect is {F_Lambda_10kpc/F_Newton_10kpc:.1e} of Newton at 10 kpc")
print(f"    → COMPLETELY negligible (need ratio ~ 1 for flat RC)")

# ===========================================================================
# 2. TGP BOUNDARY CONDITION: FINITE SUBSTRATE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. TGP boundary: what if the substrate has a finite extent R_H?")
print(f"{'='*78}")

print(f"""
  Standard TGP soliton: g'' + g'^2/g + 2g'/r + g = 1
  Boundary: g(0) = g0, g(inf) = 1

  Modified: g(0) = g0, g(R_H) = g_H (some value at horizon)

  What changes? Instead of g → 1 at infinity (free boundary),
  we IMPOSE g = g_H at a FINITE distance R_H.

  In terms of delta = g - 1:
  delta(0) = delta_0, delta(R_H) = delta_H

  This is a TWO-POINT BOUNDARY VALUE PROBLEM.

  But R_H ~ c/H0 ~ 10^26 m, while r_galaxy ~ 10^21 m (10 kpc).
  The galaxy is a TINY perturbation at the center of a Hubble-sized domain.

  For the galaxy potential:
  delta_galaxy(r) ~ -GM/(rc^2) ~ -10^-6 at r ~ 10 kpc

  The boundary value delta_H at r = R_H affects the solution
  only through the HOMOGENEOUS solution that fills the domain.

  In the standard (infinite domain) case:
  delta = delta_particular + A*sin(r)/r + B*cos(r)/r

  B=0 is chosen for regularity, and A is determined by matching at r→0.

  In the finite domain [0, R_H]:
  delta = delta_particular + A*sin(r)/r + B*cos(r)/r

  Now B ≠ 0 is allowed (it's regular for r > 0), and A, B are set by
  BOTH boundary conditions: delta(0) finite, delta(R_H) = delta_H.

  The cos(r)/r term:
  - At r = 10 kpc (in soliton units, assuming l_sol ~ fm):
    r/l_sol ~ 3e37
    cos(3e37) / 3e37 ~ 10^-37 → NEGLIGIBLE
  - But its PHASE matters: if we sum many oscillations,
    the integrated effect could be different.

  Let me solve this directly for a toy model.
""")

# ===========================================================================
# 3. TOY MODEL: 1D WITH FINITE DOMAIN
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. Toy model: delta'' + delta = -rho(r) on [0, R]")
print(f"{'='*78}")

print(f"""
  Simplified (1D, linearized):
  delta'' + delta = -rho(r) on [0, R], delta'(0) = 0, delta(R) = 0

  Source: point mass at r=0 → rho = M*dirac(r)
  (or softened: rho = M * exp(-r^2/(2*sigma^2)) / (sqrt(2pi)*sigma))

  Newtonian analogue: delta_N'' = -rho, delta_N(R) = 0
  Solution: delta_N(r) = M*(R-r) (linear decline)

  TGP: delta'' + delta = -rho
  Green's function: G(r,r') = sin(r<)*cos(R-r>) / cos(R)
    where r< = min(r,r'), r> = max(r,r')
""")

# Solve for different R values
print(f"  Solving for point mass at r=0, varying domain size R:")
print(f"  {'R':>6s} {'delta(0)':>10s} {'delta_N(0)':>10s} {'ratio':>8s} {'Comment':>25s}")
print(f"  {'-'*62}")

for R in [3, 5, 10, 20, 50, 100, 500, 1000, 1e4, 1e6]:
    # TGP Green's function at r=0 for source at r'=0:
    # G(0,0) = sin(0)*cos(R)/cos(R) = 0 (!)
    # Need softened source. Use sigma = 0.1
    sigma = 0.1
    N_pts = max(5000, int(R * 10))
    if R > 1e4:
        N_pts = 50000
    r = np.linspace(0.001, R, min(N_pts, 100000))
    dr = r[1] - r[0]

    # Source
    rho = np.exp(-r**2/(2*sigma**2)) / (np.sqrt(2*np.pi)*sigma)

    # Solve delta'' + delta = -rho with delta'(0) = 0, delta(R) = 0
    # Particular solution via variation of parameters:
    # delta_p(r) = -sin(r) * integral_0^r cos(r')*rho(r') dr'
    #            + cos(r) * integral_0^r sin(r')*rho(r') dr'

    cos_rho = np.cumsum(np.cos(r) * rho) * dr
    sin_rho = np.cumsum(np.sin(r) * rho) * dr

    delta_p = -np.sin(r) * cos_rho + np.cos(r) * sin_rho

    # Homogeneous: delta_h = A*sin(r) + B*cos(r)
    # BC1: delta'(0) = 0 → delta_p'(0) + A*cos(0) - B*sin(0) = 0
    #   → A = -delta_p'(0)
    # BC2: delta(R) = 0 → delta_p(R) + A*sin(R) + B*cos(R) = 0
    #   → B = -(delta_p(R) + A*sin(R)) / cos(R)

    # Approximate delta_p'(0)
    if len(delta_p) > 2:
        delta_p_prime_0 = (delta_p[1] - delta_p[0]) / dr
    else:
        delta_p_prime_0 = 0

    A = -delta_p_prime_0

    cos_R = np.cos(R)
    if abs(cos_R) < 1e-10:
        # Near resonance - skip
        print(f"  {R:6.0f} {'RESONANCE':>10s} {R:10.1f} {'---':>8s} {'cos(R)~0, resonance':>25s}")
        continue

    B = -(delta_p[-1] + A * np.sin(R)) / cos_R

    delta_total = delta_p + A*np.sin(r) + B*np.cos(r)

    # Newtonian solution: delta_N'' = -rho, delta_N'(0) = 0, delta_N(R) = 0
    # delta_N(r) = -integral_0^r integral_0^s rho(t)dt ds + C1*r + C2
    # For point mass: delta_N(r) = -(R-r) * M (approximately)
    integral_rho = np.cumsum(rho) * dr
    integral2_rho = np.cumsum(integral_rho) * dr
    delta_N = -integral2_rho
    C1 = (delta_N[-1]) / R
    delta_N = delta_N - C1 * r  # Fix BC at R

    delta0_TGP = delta_total[0]
    delta0_N = delta_N[0]
    ratio = delta0_TGP / delta0_N if abs(delta0_N) > 1e-20 else 0

    if abs(ratio - 1) > 0.1:
        comment = f"DIFFERENT by {abs(ratio-1)*100:.0f}%"
    else:
        comment = "similar to Newton"

    print(f"  {R:6.0f} {delta0_TGP:10.4f} {delta0_N:10.4f} {ratio:8.3f} {comment:>25s}")

# ===========================================================================
# 4. 3D SPHERICAL WITH BOUNDARY
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. 3D spherical: delta'' + 2delta'/r + delta = -rho on [0, R]")
print(f"{'='*78}")

print(f"""
  In 3D with u = r*delta:
  u'' + u = -r*rho(r)  on [0, R]
  BC: u(0) = 0, u(R) = 0  (delta finite at r=0, delta(R)=0)

  This is an eigenvalue problem! The homogeneous equation u'' + u = 0
  on [0, R] with u(0)=u(R)=0 has solutions only when R = n*pi.

  If R is NOT a multiple of pi:
  unique solution exists, and it differs from Newtonian by the spring term.

  If R = n*pi: RESONANCE — no solution (or divergent solution).

  The resonance radii in soliton units: R = pi, 2pi, 3pi, ...

  In physical units (with l_soliton ~ fm ~ 10^-15 m):
  R_resonance = n * pi * l_soliton ~ n * 3e-15 m (MICROSCOPIC!)

  At GALACTIC scales (r ~ 3e22 m):
  n ~ r / (pi * l_sol) ~ 10^37

  We're between the 10^37-th and 10^37+1-th resonance.
  The EFFECTIVE response depends on how close we are to resonance.
""")

# Solve the BVP for different R values (in soliton units)
def solve_3d_bvp(R, M_source=1.0, a_source=1.0, N=10000):
    """Solve u'' + u = -r*rho(r) with u(0)=u(R)=0
    where rho is a Hernquist-like profile softened at center."""
    r = np.linspace(0.001, R, N)
    dr = r[1] - r[0]

    # Source: Hernquist
    rho = M_source * a_source / (2*np.pi * np.maximum(r, 0.01) * (r + a_source)**3)
    source = r * rho

    # Particular solution via Green's function
    # G(r,r') = sin(r<)*sin(R-r>) / sin(R)  (for u'' + u = f, u(0)=u(R)=0)
    sin_R = np.sin(R)
    if abs(sin_R) < 1e-10:
        return None, None, None

    # Compute via integration
    # u_p(r) = (1/sin(R)) * [sin(r) * int_r^R sin(R-r')*f(r') dr'
    #                        + sin(R-r) * int_0^r sin(r')*f(r') dr']

    sin_r = np.sin(r)
    sin_Rr = np.sin(R - r)

    # Forward integral: I1(r) = int_0^r sin(r')*f(r') dr'
    I1 = np.cumsum(sin_r * source) * dr

    # Backward integral: I2(r) = int_r^R sin(R-r')*f(r') dr'
    I2 = np.cumsum((sin_Rr * source)[::-1])[::-1] * dr

    u_p = (sin_r * I2 + sin_Rr * I1) / sin_R

    delta = u_p / r  # delta = u/r

    # Newtonian: u_N'' = -r*rho, u_N(0)=u_N(R)=0
    # Integrate twice
    f = -source
    u_N_pp = f
    u_N_p = np.cumsum(f) * dr  # u'
    u_N = np.cumsum(u_N_p) * dr  # u
    # Fix BC: u_N(R) = 0
    u_N = u_N - u_N[-1] * r / R
    delta_N = u_N / r

    return r, delta, delta_N

print(f"\n  Testing domain sizes (3D, soliton units):")
print(f"  {'R':>8s} {'sin(R)':>10s} {'delta(1)':>12s} {'delta_N(1)':>12s} {'ratio':>8s}")
print(f"  {'-'*55}")

for R in [5, 10, 20, 50, 100, 200, 500, 1000]:
    result = solve_3d_bvp(R, M_source=1.0, a_source=2.0)
    if result[0] is not None:
        r, delta, delta_N = result
        idx1 = np.argmin(np.abs(r - 1.0))
        d1 = delta[idx1]
        dN1 = delta_N[idx1]
        ratio = d1/dN1 if abs(dN1) > 1e-20 else 0
        print(f"  {R:8.0f} {np.sin(R):10.4f} {d1:12.4e} {dN1:12.4e} {ratio:8.3f}")
    else:
        print(f"  {R:8.0f} {'~0 (res.)':>10s} {'RESONANCE':>12s}")

# ===========================================================================
# 5. THE KEY QUESTION: DOES THE BC AT R_H PROPAGATE INWARD?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. Does the boundary condition at R_H propagate inward?")
print(f"{'='*78}")

print(f"""
  The solution in the domain [0, R] is:
  delta(r) = delta_particular(r) + C * sin(R-r) / (r * sin(R))

  where C is determined by the BC at R.

  The homogeneous correction C*sin(R-r)/(r*sin(R)):
  - At r = 1 (near source): ~ C*sin(R)/(sin(R)) = C (if R>>1, sin(R-1)~sin(R))
  - So the BC at R affects the solution EVERYWHERE, not just near R.

  BUT: the correction is of order |C/sin(R)|.
  For R = 10^37 (galaxy in soliton units):
  sin(R) oscillates between -1 and 1 → |C/sin(R)| ~ |C|.

  What is C? It's determined by:
  delta(R) = delta_particular(R) + C * sin(0)/(R*sin(R)) = delta_H
  → 0 = delta_H → C is undetermined (sin(0)=0!)

  Wait — the BC is already satisfied trivially because sin(R-R) = sin(0) = 0.
  The coefficient C is determined by the OTHER BC (regularity at r=0).

  So the finite-domain BC doesn't add new information beyond the
  requirement that delta(0) is finite — SAME as infinite domain!

  CONCLUSION: For the linearized equation, the BC at R_H
  does NOT modify the solution near the galaxy.
  The oscillating Green's function means the BC at R affects
  the solution, but only through resonance effects that average
  out for R >> lambda.
""")

# ===========================================================================
# 6. NONLINEAR BOUNDARY CONDITION?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. Nonlinear effects at the boundary")
print(f"{'='*78}")

print(f"""
  The LINEAR equation delta'' + 2delta'/r + delta = source
  has oscillatory Green's function → BC at R doesn't help.

  But the FULL TGP equation:
  g'' + g'^2/g + 2g'/r + g = 1

  is NONLINEAR (the g'^2/g term).

  Could nonlinear effects near the horizon change things?

  At the Hubble horizon, the expansion velocity = c.
  In TGP, if g represents the substrate metric:
  g_cosmo(R_H) might be significantly different from 1.

  In fact, if expansion velocity = c at R_H:
  v_expansion = H0 * R_H = c → g might approach a CRITICAL value.

  The soliton equation has a critical value: g can't go below 0
  (square root of metric → undefined). If expansion pushes g toward 0,
  then near the horizon g → 0, which is the NONLINEAR regime!

  But wait: g is the metric coefficient. If g → 0, that's a singularity.
  In TGP, g = 0 is the maximum deformation of the substrate.

  For an expanding universe:
  - g = 1: flat, undistorted substrate
  - g > 1: compressed (mass/energy present)
  - g < 1: stretched (expansion)
  - g → 0: maximum stretch (horizon?)

  If g_cosmo → 0 at r → R_H, then the soliton equation becomes
  strongly nonlinear near the horizon.

  The term g'^2/g DIVERGES as g → 0!
  g'^2/g ~ (delta g)^2 / g → large when g → 0.

  This could create a NONLINEAR WALL at the horizon that
  changes how gravity propagates.
""")

# ===========================================================================
# 7. MODEL: g PROFILE WITH EXPANSION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. Model: g(r) profile with cosmological expansion")
print(f"{'='*78}")

print(f"""
  Consider a static galaxy in an expanding substrate.
  The combined g(r) profile:

  Near galaxy (r < R_gal): g ~ 1 + delta_galaxy (delta ~ -10^-6, weak field)
  Far from galaxy (r >> R_gal): g → g_cosmo(r) where g_cosmo represents expansion
  At horizon (r → R_H): g → g_horizon

  In standard GR-like terms:
  g_cosmo(r) = 1 - (H0*r/c)^2 (de Sitter-like for the spatial metric)

  This means g_cosmo DECREASES with distance from the galaxy.
  At r = R_H: g_cosmo = 1 - 1 = 0 (horizon!)

  So the EFFECTIVE equation for the galaxy potential:
  delta'' + 2delta'/r + delta = source
  where delta = g - g_cosmo(r) and the boundary is delta(R_H) = some value.

  But crucially: g_cosmo(r) = 1 - H0^2*r^2/c^2 varies SLOWLY
  (derivative ~ 2*H0^2*r/c^2 ~ 10^-36/m at galactic scales).

  The variation of g_cosmo over one oscillation wavelength (2*pi*l_sol):
  Delta(g_cosmo) ~ 2*H0^2*l_sol^2/c^2 * r/l_sol ~ 10^-72 * r/l_sol

  → NEGLIGIBLE even at r = R_H.

  Let me compute the actual profile numerically.
""")

# Model: g'' + g'^2/g + 2g'/r + g = 1
# with initial conditions representing galaxy + expansion
# g(0) = g0 (galaxy center), g'(0) = 0

# First: pure expansion (no galaxy)
# g_cosmo(r): solve g'' + g'^2/g + 2g'/r + g = 1 with g(0) = 1-epsilon, g'(0) = 0
# (expanding from slightly below 1)

print(f"  Pure cosmological expansion profiles:")
print(f"  {'g0':>8s} {'g(10)':>10s} {'g(50)':>10s} {'g(100)':>10s} {'r_zero':>10s}")
print(f"  {'-'*50}")

for g0 in [0.999, 0.99, 0.95, 0.9, 0.5, 0.1]:
    def rhs(r, y):
        g, gp = y
        if r < 1e-10 or g < 1e-10:
            return [gp, 0]
        gpp = -gp**2/g - 2*gp/r - g + 1
        return [gp, gpp]

    sol = solve_ivp(rhs, (1e-4, 200), [g0, 0], method='RK45',
                    t_eval=np.linspace(1e-4, 200, 10000),
                    rtol=1e-10, atol=1e-12, max_step=0.1)

    r = sol.t
    g = sol.y[0]
    idx10 = np.argmin(np.abs(r - 10))
    idx50 = np.argmin(np.abs(r - 50))
    idx100 = np.argmin(np.abs(r - 100))

    # Find where g crosses zero (if it does)
    zero_crossings = np.where(np.diff(np.sign(g)))[0]
    r_zero = r[zero_crossings[0]] if len(zero_crossings) > 0 else ">200"

    print(f"  {g0:8.3f} {g[idx10]:10.4f} {g[idx50]:10.4f} {g[idx100]:10.4f} {str(r_zero):>10s}")

# ===========================================================================
# 8. ASSESSMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. ASSESSMENT of Option A")
print(f"{'='*78}")

print(f"""
  TESTED:
  1. Schwarzschild-de Sitter analogy:
     Lambda force at 10 kpc = {F_Lambda_10kpc/F_Newton_10kpc:.1e} of Newton → NEGLIGIBLE

  2. Finite domain BVP:
     Linear equation: BC at R propagates through Green's function,
     but for R >> lambda, effect averages out → no net modification

  3. Expansion profile g_cosmo(r) = 1 - H0^2*r^2/c^2:
     Variation over oscillation wavelength ~ 10^-72 → NEGLIGIBLE

  4. Nonlinear regime at horizon (g → 0):
     g'^2/g diverges → interesting physics BUT
     this is at R_H ~ 10^26 m, galaxy is at r ~ 10^21 m
     → information must propagate 10^5 R_galaxy to matter

  5. Pure expansion profiles (g0 < 1):
     g oscillates around 1 (spring term pulls back)
     → Expansion doesn't make g monotonically decrease to 0
     → The "horizon at g=0" doesn't exist in static TGP!

  FUNDAMENTAL ISSUE:
  The TGP equation has a SPRING TERM (+g = 1) that pulls g toward 1.
  Any perturbation (expansion, galaxy, boundary) gets pulled back.
  The spring prevents g from deviating far from 1 at large distances.

  This is GOOD for particle physics (confinement of soliton)
  but BAD for galaxy physics (no long-range modification).

  STATUS: FAILED
  - Lambda force too weak by ~10^16
  - BC at horizon doesn't propagate for oscillatory Green's function
  - Spring term prevents nonlinear accumulation
  - No mechanism to produce a0
""")
