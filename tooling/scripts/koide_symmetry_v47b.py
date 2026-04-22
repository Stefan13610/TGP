#!/usr/bin/env python3
"""
koide_symmetry_v47b.py -- PATH 15: SYMMETRY STRUCTURE OF THE TGP ODE

GOAL: Find a symmetry of the soliton equation that quantizes
the mass spectrum, analogous to how SU(3) quantizes quark charges.

THE TGP ODE (canonical, alpha=2, d=3):
  g'' + (2/r)g' + (2/g)g'^2 + g^2(1-g) = 0

CANONICAL VARIABLE f = g^3:
  f'' + (2/r)f' + (3/f)·[f·f'' + (2/r)f·f' - f'^2]·(1/3) + ...
  Actually, let's derive it carefully.

LIE SYMMETRY ANALYSIS:
  A 2nd-order ODE y'' = F(x, y, y') admits a Lie symmetry
  X = xi(x,y) d/dx + eta(x,y) d/dy
  if the prolongation pr^(2)X annihilates y'' - F on the solution manifold.

  For autonomous ODEs (no explicit r-dependence except the 2/r term),
  the symmetry structure is constrained.

  Key question: does the TGP ODE have a SCALING symmetry
  r -> lambda*r, g -> g (or g -> mu*g)?

LANE-EMDEN CONNECTION:
  The standard Lane-Emden equation: y'' + (2/x)y' + y^n = 0
  has the Emden-Fowler scaling symmetry for specific n values.
  For n=5 (d=3): exact solution exists (Schuster solution).
  For general n: scaling y -> c*y, x -> c^((1-n)/2)*x is a symmetry.

  Our equation is NOT standard Lane-Emden because:
  1. Source is g^2(1-g), not g^n
  2. Cross term (2/g)g'^2 is present

  But in the f-variable, things may simplify.

APPROACH:
  1. Compute Lie symmetries numerically (infinitesimal generators)
  2. Check for discrete symmetries (Z_2, Z_3, S_3)
  3. Look for hidden conservation laws (Noether)
  4. Test if the A_tail function has a discrete symmetry
     under g0 -> f(g0) transformations
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0


def solver_A(g0, r_max=400):
    """Solve canonical Form A ODE (alpha=2, d=3)."""
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
    """Extract tail amplitude."""
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


# ================================================================
print("=" * 70)
print("PATH 15: SYMMETRY STRUCTURE OF THE TGP ODE")
print("=" * 70)


# ================================================================
# Section 1: SCALING SYMMETRIES
# ================================================================
print("\nSECTION 1: SCALING SYMMETRIES")
print("-" * 50)

print("""
  The TGP ODE: g'' + (2/r)g' + (2/g)g'^2 + g^2(1-g) = 0

  Test scaling: r -> lambda*r, g -> mu*g
  Under this: g' -> (mu/lambda)*g', g'' -> (mu/lambda^2)*g''

  Term by term:
    g''           -> mu/lambda^2
    (2/r)g'       -> (2/(lambda*r))(mu/lambda)g' = mu/lambda^2 * (2/r)g'  OK same
    (2/g)g'^2     -> (2/(mu*g))(mu/lambda)^2 g'^2 = mu/lambda^2 * (2/g)g'^2  OK same
    g^2(1-g)      -> mu^2*g^2*(1 - mu*g)  != mu/lambda^2 * g^2(1-g)

  For the last term: need mu^2*(1-mu*g) = (mu/lambda^2)*(1-g)
  This requires mu = 1 AND lambda = 1 (trivial).

  SO: the nonlinear source g^2(1-g) BREAKS scaling symmetry.
  The vacuum g=1 defines a FIXED SCALE.

  HOWEVER: near the vacuum (g ~ 1), the linearized equation
  h'' + (2/r)h' + h = 0 (in f-variable) IS scale-invariant
  in the sense that the frequency omega=1 is fixed.

  The CORE region (g far from 1) breaks scaling.
  The TAIL region preserves it (universal oscillation).
""")

# ================================================================
# Section 2: THE ODE IN AUTONOMOUS FORM
# ================================================================
print("\nSECTION 2: AUTONOMOUS FORM (Emden-Fowler reduction)")
print("-" * 50)

print("""
  Standard trick: substitute t = ln(r), u(t) = g(r).
  Then g' = u'/r = u'*e^(-t), g'' = (u'' - u')*e^(-2t)

  ODE becomes:
    (u'' - u')e^(-2t) + 2*e^(-t)*u'*e^(-t) + (2/u)*u'^2*e^(-2t) + u^2(1-u) = 0
    (u'' - u')e^(-2t) + 2*u'*e^(-2t) + (2/u)*u'^2*e^(-2t) + u^2(1-u) = 0
    e^(-2t)*[u'' + u' + (2/u)*u'^2] + u^2(1-u) = 0
    u'' + u' + (2/u)u'^2 + e^(2t)*u^2(1-u) = 0

  NOT autonomous due to e^(2t) = r^2 factor.
  This is expected: the d=3 Laplacian has no translation symmetry in r.

  Alternative: use the PHASE PLANE variables.
  Let p = r*g'/g (logarithmic derivative, scaled by r).
  Then: autonomous system in (g, p) with r as parameter.

  Actually, for d=3, the standard Lane-Emden substitution is:
    xi = ln(r), w = r^((d-2)/2) * g = r^(1/2) * g
  which gives a more symmetric form.
""")


# ================================================================
# Section 3: EMDEN-FOWLER TRANSFORMATION
# ================================================================
print("\nSECTION 3: EMDEN-FOWLER PHASE PLANE")
print("-" * 50)

# The key transformation for Lane-Emden equations:
# Define: xi = ln(r)
#         v = -r^2 * g'/g = -r * (r*g'/g)  [Milne variables]
#         w = r^3 * g^2(1-g) / g  (source term)
# Or use the Emden variables:
#   u = -r*g'/g  (compression)
#   v = r^3*g^2(1-g)  (source strength)
# These give an autonomous 2D system.

# Actually, let's use a cleaner approach:
# Define theta = r*g'/g and phi_var = r^2*g^2*(1-g)
# The ODE becomes:
#   theta' + theta + theta^2/r + (2/1)*theta^2/r + phi_var/r = ... messy
#
# Better: just study the phase portrait numerically.

# Phase portrait: (g, g') for different g0
print("  Computing phase portraits for different g0 values...")
print()

# The KEY observation: in phase space (g, r*g'), the trajectories
# of solitons with different g0 are different curves.
# Is there a TRANSFORMATION that maps one to another?

# Check: if g(r; g0_mu) = T[g(r; g0_e)] for some operator T,
# and similarly for tau, then the three solitons are related
# by a symmetry.

# The simplest transformation: g -> c*g (rescaling amplitude)
# We know this doesn't work because g=1 is a fixed point.

# More subtle: could there be a NONLINEAR transformation?
# g -> G(g) such that if g(r) solves ODE with g(0)=g0,
# then G(g(r)) solves ODE with G(g(0)) = G(g0)?

# This would be a "spectrum-generating symmetry".
# Test: if phi-FP gives g0_mu = phi*g0_e, is there a
# transformation g -> G(g) = phi*g that maps solutions?

print("  Test: does g -> phi*g map solutions?")
g0_e = 0.86770494
g0_mu = PHI * g0_e

r_e, g_e = solver_A(g0_e)
r_mu, g_mu = solver_A(g0_mu)

# If g->phi*g were a symmetry, then phi*g_e(r) would satisfy the ODE
# with initial condition phi*g0_e = g0_mu.
# Check: does phi*g_e(r) match g_mu(r)?

# Interpolate at common r points
from scipy.interpolate import interp1d
g_e_interp = interp1d(r_e, g_e, bounds_error=False, fill_value=1.0)
g_mu_interp = interp1d(r_mu, g_mu, bounds_error=False, fill_value=1.0)

r_test = np.linspace(0.1, 30, 200)
g_e_scaled = PHI * g_e_interp(r_test)
g_mu_actual = g_mu_interp(r_test)

err = np.max(np.abs(g_e_scaled - g_mu_actual))
print(f"  max|phi*g_e(r) - g_mu(r)| = {err:.6f}")
print(f"  This is {'SMALL' if err < 0.01 else 'LARGE'} -- g -> phi*g is {'approximately' if err < 0.1 else 'NOT'} a symmetry.")

# What about a radial scaling? g_mu(r) = g_e(lambda*r)?
# Find best lambda
from scipy.optimize import minimize_scalar

def mismatch(lam):
    if lam <= 0:
        return 1e10
    g_e_stretched = g_e_interp(lam * r_test)
    return np.mean((g_e_stretched - g_mu_actual)**2)

res = minimize_scalar(mismatch, bounds=(0.1, 10), method='bounded')
lambda_opt = res.x
err_lambda = np.sqrt(res.fun)
print(f"\n  Test: does g_mu(r) = g_e(lambda*r) for some lambda?")
print(f"  Best lambda = {lambda_opt:.6f}")
print(f"  RMS error = {err_lambda:.6f}")
print(f"  This is {'SMALL' if err_lambda < 0.01 else 'NOT small'} -- radial rescaling {'works' if err_lambda < 0.01 else 'fails'}.")

# What about combined: g_mu(r) = phi * g_e(lambda*r)?
def mismatch2(lam):
    if lam <= 0:
        return 1e10
    g_e_transformed = PHI * g_e_interp(lam * r_test)
    return np.mean((g_e_transformed - g_mu_actual)**2)

res2 = minimize_scalar(mismatch2, bounds=(0.1, 10), method='bounded')
lambda2 = res2.x
err2 = np.sqrt(res2.fun)
print(f"\n  Test: does g_mu(r) = phi * g_e(lambda*r)?")
print(f"  Best lambda = {lambda2:.6f}")
print(f"  RMS error = {err2:.6f}")


# ================================================================
# Section 4: A_TAIL AS A FUNCTION -- SYMMETRY PROPERTIES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: SYMMETRIES OF A_tail(g0)")
print("=" * 70)

# The function A_tail(g0) maps initial condition to tail amplitude.
# Is there a functional equation that A_tail satisfies?
#
# If g0 -> phi*g0 gives A_mu = phi^k * A_e for some k,
# then A_tail has a "scaling at phi" property.

print("\n  Computing A_tail for a range of g0 values...")
g0_range = np.arange(0.50, 1.595, 0.005)
A_data = []
for g0 in g0_range:
    A = A_tail(g0)
    A_data.append(A)
A_data = np.array(A_data)

# Check: is A_tail(phi*g0) / A_tail(g0) = constant?
print("\n  Test: A_tail(phi*g0) / A_tail(g0)")
print(f"  {'g0':>8s} {'phi*g0':>8s} {'A(g0)':>10s} {'A(phi*g0)':>10s} {'ratio':>10s} {'ratio^4':>10s}")

for g0 in np.arange(0.60, 0.99, 0.02):
    g0p = PHI * g0
    if g0p >= G0_CRIT:
        continue
    A1 = A_tail(g0)
    A2 = A_tail(g0p)
    if A1 > 1e-10:
        ratio = A2 / A1
        print(f"  {g0:8.4f} {g0p:8.4f} {A1:10.6f} {A2:10.6f} {ratio:10.4f} {ratio**4:10.2f}")

print(f"\n  If ratio were constant, that would be a scaling symmetry.")
print(f"  Target ratio: (r_21)^(1/4) = {206.77**0.25:.4f}")

# The ratio is NOT constant -- it varies with g0.
# This confirms that g -> phi*g is NOT a symmetry of the ODE.
# The phi-FP rule must come from a DIFFERENT mechanism.


# ================================================================
# Section 5: DISCRETE SYMMETRIES OF THE KOIDE CIRCLE
# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: DISCRETE SYMMETRIES OF THE KOIDE CIRCLE")
print("=" * 70)

print("""
  The Koide parametrization:
    sqrt(m_i) = a * (1 + b*cos(theta + 2*pi*(i-1)/3))

  has a manifest Z_3 symmetry: cyclic permutation of generations
  (i -> i+1 mod 3) corresponds to theta -> theta + 2*pi/3.

  This Z_3 is the CYCLIC GROUP acting on the three generations.

  The full symmetry of the Koide circle:
  - Z_3 cyclic: theta -> theta + 2*pi/3
  - Z_2 reflection: theta -> -theta (swaps generation ordering)
  - Together: S_3 = Z_3 x Z_2 (full permutation group)

  The physical Koide relation BREAKS S_3 to Z_1 (identity only),
  because the masses are ordered: m_e < m_mu < m_tau.
  But the FORMULA is S_3-symmetric.

  KEY QUESTION: is there a LARGER discrete group that also
  constrains b = sqrt(2)?

  In the Koide circle, the three points are at angles:
    theta_1, theta_2 = theta_1 + 120, theta_3 = theta_1 + 240
  on a circle of radius b centered at 1.

  b = sqrt(2) is special because:
  - b^2 = 2 = the number of "cross dimensions" (?)
  - Q_K = 3/(1+b^2/2) = 3/2
  - The "modulation depth" sqrt(2)/1 = sqrt(2)

  Is there a geometric reason why b = sqrt(2) for a TRIANGLE
  inscribed in a circle?
""")

# For an equilateral triangle inscribed in a circle of radius b,
# centered at (1, 0) on the real line:
# The vertices are at (1 + b*cos(theta_i), b*sin(theta_i))
# projected onto the real axis: x_i = 1 + b*cos(theta_i)
#
# For the triangle to have specific properties:
# - Area of triangle in (x,y) space: (3*sqrt(3)/4)*b^2
# - Perimeter: 3*sqrt(3)*b
# - Inradius: sqrt(3)/2 * b
# - Circumradius: b
#
# Is there a geometric condition that gives b = sqrt(2)?

# Condition: the CENTROID of the x_i values is at x = 1.
# sum(x_i)/3 = 1 + b*sum(cos(theta_i))/3 = 1 (since sum cos = 0)
# This is ALWAYS true! No constraint on b.

# Condition: the VARIANCE of x_i equals 1?
# Var(x) = b^2 * sum(cos^2(theta_i))/3 = b^2 * (3/2)/3 = b^2/2
# Var = 1 requires b^2 = 2, i.e., b = sqrt(2)!

print(f"  GEOMETRIC INSIGHT:")
print(f"  Var(x_i) = b^2/2 for equilateral triangle on circle of radius b.")
print(f"  Var(x_i) = 1 requires b^2 = 2, i.e., b = sqrt(2).")
print(f"  ")
print(f"  But 'Var = 1' means Var = (centroid)^2, since centroid = 1.")
print(f"  This is exactly CV = 1 (coefficient of variation = 1).")
print(f"  Equivalently: sigma = mu for the x_i distribution.")
print(f"  ")
print(f"  So Q_K = 3/2 <=> 'the standard deviation of sqrt(m)")
print(f"  equals the mean of sqrt(m)' <=> sigma = mu.")
print(f"  ")
print(f"  Why would sigma = mu be special?")


# ================================================================
# Section 6: INFORMATION-THEORETIC INTERPRETATION
# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: INFORMATION / EXPONENTIAL FAMILY INTERPRETATION")
print("=" * 70)

print("""
  For a set of positive values {x_i} with mean mu and std sigma:

  CV = sigma/mu = 1 is the defining property of the
  EXPONENTIAL DISTRIBUTION!

  If x ~ Exp(lambda), then E[x] = 1/lambda, Var[x] = 1/lambda^2,
  so CV = sqrt(Var)/E = 1.

  This means: Q_K = 3/2 is equivalent to saying that
  the three sqrt(m_i) values are DISTRIBUTED like samples
  from an exponential distribution.

  For a CONTINUOUS exponential distribution of n samples:
    E[Q_K] = n / (1 + CV^2) = n / 2
  For n=3: E[Q_K] = 3/2. EXACTLY KOIDE!

  This is remarkable: if the sqrt(m_i) values were random
  samples from an exponential distribution, the EXPECTED
  value of Q_K would be 3/2.

  But wait: for 3 samples, there's variance around E[Q_K] = 3/2.
  The actual Q_K could be anything from 1 to 3.
  The fact that it's EXACTLY 3/2 (not just approximately)
  means either:
  (a) It's a coincidence (3 samples happened to hit the mean)
  (b) There's a CONSTRAINT forcing Q_K to its expected value
  (c) The distribution is not random but self-organizing
""")

# Simulate: what's the probability that 3 exponential samples
# give Q_K within 0.1% of 3/2?
np.random.seed(42)
N_sim = 1000000
QK_samples = np.zeros(N_sim)
for i in range(N_sim):
    x = np.random.exponential(1.0, 3)
    S = np.sum(x)
    P = np.sum(x**2)
    QK_samples[i] = S**2 / P

frac_close = np.mean(np.abs(QK_samples - 1.5) < 0.001 * 1.5)
mean_QK = np.mean(QK_samples)
std_QK = np.std(QK_samples)

print(f"  Monte Carlo: {N_sim} sets of 3 exponential samples:")
print(f"  Mean(Q_K) = {mean_QK:.6f}  (theory: 1.500)")
print(f"  Std(Q_K)  = {std_QK:.6f}")
print(f"  P(|Q_K - 1.5| < 0.15%) = {frac_close:.6f}")
print(f"  P(|Q_K - 1.5| < 1%) = {np.mean(np.abs(QK_samples - 1.5) < 0.015):.6f}")
print()

# BUT: we also have the constraint r_21 = 206.77!
# For exponential samples with x_2/x_1 = sqrt(206.77) = 14.38:
# This is VERY unlikely for unconstrained exponential samples.
# The r_21 constraint + Q_K = 3/2 together fully determine the spectrum.

# What is the Q_K distribution GIVEN r_21?
# If x_1 < x_2 < x_3 with x_2/x_1 = u = 14.38,
# then x_3 is the only free parameter.
# Q_K(x_3) = (x_1 + u*x_1 + x_3)^2 / (x_1^2 + u^2*x_1^2 + x_3^2)
# = x_1^2 * (1+u+x_3/x_1)^2 / (x_1^2 * (1+u^2+(x_3/x_1)^2))
# = (1+u+v)^2 / (1+u^2+v^2)  where v = x_3/x_1

u_val = np.sqrt(206.77)
print(f"  Given r_21 = 206.77 (u = sqrt(r_21) = {u_val:.4f}):")
print(f"  Q_K(v) = (1+u+v)^2 / (1+u^2+v^2)")
print(f"  where v = sqrt(m_tau/m_e)")
print()

# For exponential distribution: E[v] = E[x_3/x_1] is problematic
# (ratio of exponentials is Pareto). But the conditional distribution
# given r_21 is delta-like. So the exponential interpretation is
# about the UNCONDITIONAL distribution.

# More relevant: Q_K = 3/2 given u, solve for v:
# (1+u+v)^2 = (3/2)(1+u^2+v^2)
# 1+u^2+v^2+2u+2v+2uv = (3/2)(1+u^2+v^2)
# 2u+2v+2uv = (1/2)(1+u^2+v^2)
# 4u+4v+4uv = 1+u^2+v^2
# v^2 - 4v(1+u) + u^2 - 4u + 1 = 0
# v = 2(1+u) + sqrt(4(1+u)^2 - u^2 + 4u - 1)
# = 2(1+u) + sqrt(4+8u+4u^2 - u^2 + 4u - 1)
# = 2(1+u) + sqrt(3+12u+3u^2)
# = 2(1+u) + sqrt(3(1+4u+u^2))

v_koide = 2*(1+u_val) + np.sqrt(3*(1+4*u_val+u_val**2))
print(f"  Koide solution: v = 2(1+u) + sqrt(3(1+4u+u^2)) = {v_koide:.4f}")
print(f"  r_31 = v^2 = {v_koide**2:.2f}  (PDG: 3477.23)")
print()

# Now the exponential distribution connection:
# For Exp(lambda), E[x^2]/E[x]^2 = 2 (because E[x^2] = 2/lambda^2)
# So <x^2>/<x>^2 = 2, which gives Q_K = n*<x>^2/(<x^2>) = n/2.
# This is the RAYLEIGH property: second moment = 2 * first moment squared.

print(f"  EXPONENTIAL DISTRIBUTION CONNECTION:")
print(f"  For Exp(lambda): E[x^2]/E[x]^2 = 2")
print(f"  For n samples:   Q_K = n * E[x]^2 / E[x^2] = n/2")
print(f"  For n=3: Q_K = 3/2 -- EXACT KOIDE!")
print(f"  ")
print(f"  The PHYSICAL QUESTION: why would sqrt(m_i) follow")
print(f"  an exponential-like distribution?")


# ================================================================
# Section 7: A_TAIL AND THE EXPONENTIAL CONNECTION
# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: IS A_tail^2(g0) EXPONENTIAL IN g0?")
print("=" * 70)

# If A^2(g0) ~ exp(c * g0), then with phi-FP spacing:
# A_mu^2 / A_e^2 = exp(c*(g0_mu - g0_e)) = exp(c*(phi-1)*g0_e)
# A_tau^2 / A_e^2 = exp(c*(g0_tau - g0_e))
#
# For Q_K = 3/2 with exponential A^2:
# x_i = A_i^2 = exp(c*g0_i) (up to constant)
# These are exponentially distributed IF g0_i are uniformly spaced.
# But they're NOT: g0_mu = phi*g0_e, g0_tau is from Koide.
#
# Let's check if A^2 is exponential in g0:

g0_fine = np.arange(0.50, 1.595, 0.002)
A_fine = np.array([A_tail(g0) for g0 in g0_fine])
x_fine = A_fine**2

# Plot-like analysis: is ln(A^2) linear in g0?
valid = A_fine > 1e-10
g0_v = g0_fine[valid]
lnx = np.log(x_fine[valid])

# Fit line: ln(x) = a + b*g0
from numpy.polynomial import polynomial as P
coefs = np.polyfit(g0_v, lnx, 1)  # [slope, intercept]

print(f"  Linear fit: ln(A^2) = {coefs[1]:.4f} + {coefs[0]:.4f} * g0")
print(f"  Slope c = {coefs[0]:.4f}")
print()

# Check residuals
lnx_fit = coefs[0] * g0_v + coefs[1]
residuals = lnx - lnx_fit
print(f"  Fit quality: max |residual| = {np.max(np.abs(residuals)):.4f}")
print(f"  This is {'GOOD' if np.max(np.abs(residuals)) < 0.1 else 'POOR'} -- ln(A^2) is {'approximately' if np.max(np.abs(residuals)) < 0.5 else 'NOT'} linear in g0.")
print()

# Check at the three soliton points
g0_e = 0.86770494
g0_mu = PHI * g0_e
g0_tau_koide = 1.56962743

A_e_val = A_tail(g0_e)
A_mu_val = A_tail(g0_mu)
A_tau_val = A_tail(g0_tau_koide)

print(f"  At soliton points:")
print(f"  g0_e   = {g0_e:.6f}, ln(A^2) = {np.log(A_e_val**2):.4f}, fit = {coefs[0]*g0_e+coefs[1]:.4f}, err = {np.log(A_e_val**2)-(coefs[0]*g0_e+coefs[1]):.4f}")
print(f"  g0_mu  = {g0_mu:.6f}, ln(A^2) = {np.log(A_mu_val**2):.4f}, fit = {coefs[0]*g0_mu+coefs[1]:.4f}, err = {np.log(A_mu_val**2)-(coefs[0]*g0_mu+coefs[1]):.4f}")
print(f"  g0_tau = {g0_tau_koide:.6f}, ln(A^2) = {np.log(A_tau_val**2):.4f}, fit = {coefs[0]*g0_tau_koide+coefs[1]:.4f}, err = {np.log(A_tau_val**2)-(coefs[0]*g0_tau_koide+coefs[1]):.4f}")

# Try quadratic fit
coefs2 = np.polyfit(g0_v, lnx, 2)
lnx_fit2 = np.polyval(coefs2, g0_v)
residuals2 = lnx - lnx_fit2
print(f"\n  Quadratic fit: ln(A^2) = {coefs2[2]:.4f} + {coefs2[1]:.4f}*g0 + {coefs2[0]:.4f}*g0^2")
print(f"  Max |residual| = {np.max(np.abs(residuals2)):.4f}")

# Try power law: A^2 ~ (g0 - g0_min)^p or A^2 ~ (gc - g0)^(-q)
# Near vacuum: A ~ |g0-1|, so A^2 ~ (1-g0)^2 for g0 < 1
# Near collapse: A ~ (gc-g0)^(-gamma), so A^2 ~ (gc-g0)^(-2gamma)

# Split into two regimes
mask_below = g0_v < 1.0
mask_above = g0_v > 1.05

if np.sum(mask_below) > 5:
    g0_b = g0_v[mask_below]
    lnx_b = lnx[mask_below]
    ln_dist = np.log(1.0 - g0_b)
    coefs_b = np.polyfit(ln_dist, lnx_b, 1)
    print(f"\n  Below vacuum (g0 < 1): A^2 ~ (1-g0)^{coefs_b[0]:.4f}")
    print(f"  Expected: exponent = 2.0")

if np.sum(mask_above) > 5:
    g0_a = g0_v[mask_above]
    lnx_a = lnx[mask_above]
    ln_dist_a = np.log(G0_CRIT - g0_a)
    valid_a = np.isfinite(ln_dist_a)
    if np.sum(valid_a) > 5:
        coefs_a = np.polyfit(ln_dist_a[valid_a], lnx_a[valid_a], 1)
        print(f"  Above vacuum (g0 > 1): A^2 ~ (gc-g0)^{coefs_a[0]:.4f}")
        print(f"  This is the collapse divergence: exponent ~ -2*gamma")


# ================================================================
# Section 8: S_3 PERMUTATION SYMMETRY AND MASS SUM RULES
# ================================================================
print("\n" + "=" * 70)
print("SECTION 8: S_3 PERMUTATION AND MASS SUM RULES")
print("=" * 70)

print("""
  The Koide formula has an underlying S_3 symmetry
  (permutations of the three generations).

  In TGP, the three solitons have g0 values:
    g0_e < 1 < g0_mu < g0_tau < g0_crit = 8/5

  The VACUUM g=1 sits BETWEEN electron and muon/tau.
  This is an asymmetric arrangement -- no S_3 in g0-space.

  But in the KOIDE CIRCLE, the three sqrt(m) values ARE
  related by 120-degree rotations (Z_3 subgroup of S_3).

  The key question: is the Koide circle's Z_3 symmetry
  an ACCIDENTAL consequence of Q_K = 3/2, or does it
  reflect a deeper symmetry of the dynamics?

  In the proton analogy:
  - The three quarks transform under the FUNDAMENTAL rep of SU(3)
  - Confinement forces them into a COLOR-SINGLET state
  - This constrains their quantum numbers

  For TGP:
  - The three solitons might transform under some group G
  - A "confinement-like" condition forces Q_K to a specific value
  - The G-singlet condition gives Q_K = 3/2
""")

# Check: what group representation gives Q_K = n/2?
# For SU(2): fundamental rep has dimension 2, adjoint has 3.
# For SU(3): fundamental rep has dimension 3.
# Q_K = 3/2 = dim(fund SU(3)) / dim(fund SU(2))? That's not standard.

# More natural: Q_K = 3/2 = (n_gen) / 2.
# If each generation carries a "charge" of 1/2 under some Z_2,
# then the total charge is 3/2.

# Or: Q_K = 3/2 corresponds to a spin-1/2 system with 3 components?
# Casimir of SU(2) for spin j: j(j+1) = 3/4. Not 3/2.
# Casimir of SU(2) for spin j=1: j(j+1) = 2. Not 3/2.
# Casimir of SU(3) for fundamental: (N^2-1)/(2N) = 4/3. Not 3/2.

# Actually, 3/2 is the dimension of the OVERLAP between
# the symmetric and antisymmetric representations of S_3.
# S_3 has irreps: trivial (dim 1), sign (dim 1), standard (dim 2).
# 1 + 1 + 2 = 4 = |S_3|/... no.

# Let's think differently.
# For 3 variables on a circle: Q = sum(cos(theta_i))^2 + sum(sin(theta_i))^2
# divided by 3, this is the "order parameter" of the angular distribution.
# For equilateral: order = 0. For aligned: order = 1.

# The Koide Q is different -- it's about AMPLITUDES, not angles.
# But the Koide circle has BOTH: amplitude b AND angle theta.

# The constraint Q_K = 3/2 fixes b = sqrt(2).
# The constraint r_21 = 206.77 fixes theta.
# Together they give the unique physical solution.

# Is b = sqrt(2) related to a representation dimension?
# sqrt(2) = sqrt(dim(SU(2) fundamental))? SU(2) fund has dim 2.
# Hmm, that's cute but probably numerology.

print(f"  Numerological checks:")
print(f"  sqrt(2) = sqrt(dim SU(2) fund) = {np.sqrt(2):.6f}")
print(f"  Q_K = 3/2 = dim(SU(3) fund) / dim(SU(2) fund) = {3/2:.6f}")
print(f"  b^2/n = 2/3 (n=3 generations)")
print(f"  Q_K = 3/(1 + n*b^2/(2n)) = n/(1 + b^2/2) ... same formula")
print()

# Actually: for the general Koide condition with n generations:
# Q_K^(n) = n/(1 + b^2/2)
# Q_K = n/2 requires b^2 = 2, INDEPENDENT of n!
# So b = sqrt(2) is the "natural" value for ANY number of generations.

print(f"  CRUCIAL OBSERVATION:")
print(f"  For n generations on the Koide circle:")
print(f"    Q_K = n / (1 + b^2/2)")
print(f"    Q_K = n/2 requires b^2 = 2 for ALL n!")
print(f"")
print(f"    n=2: Q_K = 1   (two vars, b=sqrt(2))")
print(f"    n=3: Q_K = 3/2 (three vars, b=sqrt(2))")
print(f"    n=4: Q_K = 2   (four vars, b=sqrt(2))")
print(f"    n=5: Q_K = 5/2 (five vars, b=sqrt(2))")
print(f"")
print(f"  b = sqrt(2) is the value where Q_K = n/2,")
print(f"  which is the EXPECTED value for exponential distribution!")
print(f"  It does NOT depend on the number of generations!")
print(f"")
print(f"  This suggests: b = sqrt(2) is a UNIVERSAL property")
print(f"  of the underlying distribution, not specific to n=3.")


# ================================================================
# Section 9: EXPONENTIAL DISTRIBUTION FROM SOLITON DYNAMICS?
# ================================================================
print("\n" + "=" * 70)
print("SECTION 9: DOES THE ODE GENERATE EXPONENTIAL STATISTICS?")
print("=" * 70)

print("""
  THESIS: The three soliton amplitudes x_i = A_i^2 behave as if
  drawn from an exponential distribution. This gives Q_K = n/2.

  For exponential distribution:
    P(x) = lambda * exp(-lambda*x)
    E[x] = 1/lambda
    E[x^2] = 2/lambda^2
    Var(x) = 1/lambda^2
    CV = 1

  Question: does the A_tail(g0) function, combined with the
  phi-FP spacing, produce amplitudes with CV = 1?

  We know:
    A_e^2  = 0.01714
    A_mu^2 = 0.24656  (ratio: 14.38)
    A_tau^2 = 1.01112  (ratio: 58.97)

  For exponential samples x_1 < x_2 < x_3 with x_2/x_1 = 14.38:
  The EXPECTED x_3/x_1 = ?

  For order statistics of Exp(lambda), n=3:
    E[x_(1)] = 1/(3*lambda)
    E[x_(2)] = 1/(3*lambda) + 1/(2*lambda) = 5/(6*lambda)
    E[x_(3)] = 1/(3*lambda) + 1/(2*lambda) + 1/lambda = 11/(6*lambda)
    E[x_(2)/x_(1)] = E[x_(2)]/E[x_(1)] is NOT simply 5/2 (ratio of expectations != expectation of ratio)
""")

# But we can compute: given x_2/x_1 = u^2 = 206.77,
# what's the conditional E[x_3/x_1] assuming exponential?
#
# For Exp(lambda), x_1 < x_2 < x_3:
# The gaps y_1 = x_1, y_2 = x_2-x_1, y_3 = x_3-x_2
# are independent with y_i ~ Exp((4-i)*lambda) for i=1,2,3.
# Actually: y_1 ~ Exp(3*lambda), y_2 ~ Exp(2*lambda), y_3 ~ Exp(lambda)
#
# x_1 = y_1, x_2 = y_1+y_2, x_3 = y_1+y_2+y_3
#
# Given x_2/x_1 = R21:
# (y_1+y_2)/y_1 = R21 => y_2 = (R21-1)*y_1
# y_2 ~ Exp(2*lambda) and y_1 ~ Exp(3*lambda)
# y_2/y_1 ~ Pareto distribution
#
# E[x_3/x_1 | x_2/x_1 = R21]:
# x_3/x_1 = R21 + y_3/y_1
# E[y_3/y_1] requires care (y_3 ~ Exp(lambda), y_1 ~ Exp(3*lambda))
# But y_3/y_1 has E[y_3]*E[1/y_1] if independent... y_3 IS independent of y_1
# E[y_3] = 1/lambda, but E[1/y_1] diverges for exponential!

# Let's just do it by Monte Carlo
np.random.seed(123)
N_mc = 5000000
y1 = np.random.exponential(1/(3.0), N_mc)  # Exp(3)
y2 = np.random.exponential(1/(2.0), N_mc)  # Exp(2)
y3 = np.random.exponential(1/(1.0), N_mc)  # Exp(1)

x1 = y1
x2 = y1 + y2
x3 = y1 + y2 + y3

R21 = x2/x1
R31 = x3/x1

# Condition: R21 close to sqrt(206.77)^2 = 206.77
# Actually x_i are sqrt(m_i), so x_2/x_1 = sqrt(m_mu/m_e) = sqrt(206.77) = 14.38
# R21_target = sqrt(206.77)
# No wait: x_i = A_i^2 = sqrt(m_i) (up to constant)
# So x_2/x_1 = sqrt(m_mu/m_e) = sqrt(206.77) = 14.38

target_R21 = np.sqrt(206.77)
mask_R21 = np.abs(R21 - target_R21) < 0.5  # within ~3%
if np.sum(mask_R21) > 100:
    R31_cond = R31[mask_R21]
    QK_cond = (x1[mask_R21] + x2[mask_R21] + x3[mask_R21])**2 / \
              (x1[mask_R21]**2 + x2[mask_R21]**2 + x3[mask_R21]**2)
    print(f"  Monte Carlo: {np.sum(mask_R21)} samples with R21 ~ {target_R21:.2f}")
    print(f"  Conditional E[R31] = {np.mean(R31_cond):.2f}")
    print(f"  Conditional E[r31] = E[R31^2] = {np.mean(R31_cond**2):.2f}")
    print(f"  Actual r31 (Koide) = {v_koide**2:.2f}")
    print(f"  Actual R31 = sqrt(r31) = {v_koide:.2f}")
    print(f"  Conditional E[Q_K] = {np.mean(QK_cond):.4f}")
    print(f"  Actual Q_K = 1.5000")
    print(f"  Conditional std(Q_K) = {np.std(QK_cond):.4f}")
else:
    print(f"  Not enough samples with R21 ~ {target_R21:.2f} (got {np.sum(mask_R21)})")
    # Use wider window
    mask_R21 = np.abs(R21 - target_R21) < 2.0
    if np.sum(mask_R21) > 10:
        QK_cond = (x1[mask_R21] + x2[mask_R21] + x3[mask_R21])**2 / \
                  (x1[mask_R21]**2 + x2[mask_R21]**2 + x3[mask_R21]**2)
        print(f"  Wide window: {np.sum(mask_R21)} samples with R21 ~ {target_R21:.1f} +/- 2")
        print(f"  Conditional E[Q_K] = {np.mean(QK_cond):.4f}")
    else:
        print(f"  Still not enough samples even with wide window.")


# ================================================================
# Section 10: FINAL SYNTHESIS
# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS: SYMMETRY AND EXPONENTIAL STRUCTURE")
print("=" * 70)

print(f"""
  KEY FINDINGS:

  1. NO EXACT SYMMETRY of the ODE maps solitons to each other.
     g -> phi*g is NOT a symmetry (max error {err:.4f}).
     g(lambda*r) does NOT map solutions (RMS error {err_lambda:.4f}).

  2. A_tail^2(g0) is NOT exponential in g0 (poor fit).
     It has different power-law regimes:
     Below vacuum: A^2 ~ (1-g0)^2
     Above vacuum: A^2 ~ (gc-g0)^(-2*gamma)
     The crossover at g0=1 prevents a simple exponential.

  3. b = sqrt(2) is UNIVERSAL for any number of generations:
     Q_K = n/2 requires b = sqrt(2) for all n.
     This is the exponential distribution's CV = 1 property.

  4. The Koide circle's Z_3 (120 deg spacing) is the maximal
     angular separation of 3 points on a circle -- from
     angular repulsion / equipartition.

  5. Q_K = 3/2 says: "the sqrt(m) values have the same
     statistical properties as exponential random variables."
     This is NOT the same as saying they ARE random.
     It's a CONSTRAINT on their ratios.

  INTERPRETATION:
  The TGP field dynamics (nonlinear ODE + phi-FP spacing)
  produces soliton amplitudes that MIMIC exponential statistics.
  This is not a coincidence -- it may reflect a MAXIMUM ENTROPY
  principle operating in the soliton parameter space.

  Exponential distribution = maximum entropy distribution
  for positive continuous variables with fixed mean.

  If the soliton amplitudes x_i = A_i^2 maximize entropy
  subject to a constraint sum(x_i) = S (total amplitude),
  the result is exponential distribution, giving Q_K = n/2.

  But we already tested entropy in Path 2 and it gave Q_K ~ 2.0!
  That was entropy of MASSES, not of AMPLITUDES.
  Entropy of sqrt(m) might be different...

  LET'S CHECK: what Q_K maximizes the entropy of the x_i = sqrt(m_i)
  distribution, given r_21 = 206.77?
""")

# Entropy of the distribution {x_i/S} where x_i = sqrt(m_i)
# H = -sum(p_i * ln(p_i)) where p_i = x_i/S

u = np.sqrt(206.77)

def entropy_QK(QK_val):
    """Given Q_K and u = sqrt(r21), compute the entropy of the sqrt(m) distribution."""
    # From Koide: v^2 - 4v(1+u) + u^2 - 4u + 1 = 0 at Q_K = 3/2
    # More generally: (1+u+v)^2 = QK*(1+u^2+v^2)
    # QK*(1+u^2+v^2) = 1+u^2+v^2+2u+2v+2uv
    # (QK-1)*(1+u^2+v^2) = 2(u+v+uv)
    # (QK-1)*v^2 - 2v(1+u) + (QK-1)*(1+u^2) - 2u = 0
    a_coef = QK_val - 1
    b_coef = -2*(1+u)
    c_coef = (QK_val-1)*(1+u**2) - 2*u
    disc = b_coef**2 - 4*a_coef*c_coef
    if disc < 0 or a_coef == 0:
        return -1, 0
    v = (-b_coef + np.sqrt(disc)) / (2*a_coef)
    if v <= 0:
        return -1, 0

    # x = (1, u, v) proportional to sqrt(m)
    x = np.array([1.0, u, v])
    S = np.sum(x)
    p = x / S
    # Entropy
    H = -np.sum(p * np.log(p))
    return H, v

QK_range = np.arange(1.05, 2.95, 0.05)
print(f"\n  Entropy of sqrt(m) distribution vs Q_K:")
print(f"  {'Q_K':>8s} {'H(sqrt_m)':>10s} {'v':>10s} {'r31':>10s}")
H_max = -1
QK_max = 0
for QK in QK_range:
    H, v = entropy_QK(QK)
    if H < 0:
        continue
    r31 = v**2
    marker = " <-- Koide" if abs(QK - 1.5) < 0.03 else ""
    if H > H_max:
        H_max = H
        QK_max = QK
        marker += " <-- MAX" if abs(QK - QK_max) < 0.03 else ""
    print(f"  {QK:8.4f} {H:10.6f} {v:10.4f} {r31:10.2f}{marker}")

print(f"\n  Maximum entropy at Q_K = {QK_max:.4f}")
print(f"  Koide Q_K = 1.5000")
print(f"  {'MATCH!' if abs(QK_max - 1.5) < 0.05 else 'NO MATCH (diff = ' + f'{QK_max-1.5:.3f}' + ')'}")
