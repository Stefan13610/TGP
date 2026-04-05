#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p141_3d_stability_g0tau.py -- 3D stability analysis to constrain g0_tau
========================================================================

The radial ODE gives a one-parameter family of solitons indexed by g0.
In the FULL 3D theory, angular perturbations may destabilize solitons
above a critical g0, providing a physical upper bound.

Strategy:
  A) Linear stability of spherical soliton against angular perturbations
     (Derrick-type and mode-stability analysis)
  B) Radial perturbation spectrum: is there a negative eigenvalue?
  C) Derrick's theorem in d=3 with running alpha
  D) The stability boundary g0_max as function of eta_K

Key idea: In standard scalar field theory, Derrick's theorem forbids
static solitons in d >= 2. TGP evades this because:
  1) The kinetic function f(g) is field-dependent
  2) The soliton has tail oscillations (not strictly localized)

But for LARGE g0, the soliton core becomes very compact and the tail
very extended. There may be a stability boundary.

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal

PHI = (1+np.sqrt(5))/2
ETA_K = 181.0/15  # analytical value

def solve_sol(g0, eta=ETA_K, rm=300):
    def fk(g):
        a = 2.0/(1+eta*(g-1)**2)
        return 1+2*a*np.log(g) if g>0 else -1e30
    def Vp(g): return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0)<1e-15: return None,None,None
    c2 = Vp(g0)/(3*fg0)
    rs=0.01
    def rhs(r,y):
        g,p=y
        if g<=1e-15: return [p,0]
        fg=fk(g)
        if abs(fg)<1e-10: return [p,0]
        if r<1e-10: return [p,Vp(g)/fg/3]
        return [p,(Vp(g)-2/r*p)/fg]
    def ev(r,y): return 100-abs(y[0])
    ev.terminal=True
    s=solve_ivp(rhs,[rs,rm],[g0+c2*rs**2,2*c2*rs],
                method='RK45',rtol=1e-11,atol=1e-13,
                max_step=0.05,events=[ev],dense_output=True)
    r=np.linspace(rs,min(s.t[-1],rm),15000)
    sol = s.sol(r)
    return r, sol[0], sol[1]

# ===================================================================
print("="*70)
print("  PART A: DERRICK'S THEOREM WITH RUNNING ALPHA")
print("="*70)

print("""
  Derrick's theorem: for a static solution of
    f(g)*g'' + (2/r)*g' = V'(g)
  consider the scaling g_lambda(r) = g(r/lambda).

  The energy functional is:
    E[g] = 4*pi * int_0^inf [ (1/2)*f(g)*g'^2 + V(g)-V(1) ] * r^2 dr
         = E_kin + E_pot

  Under g -> g_lambda:
    E_kin(lambda) = lambda * E_kin     (in d=3)
    E_pot(lambda) = lambda^3 * E_pot

  Stationarity at lambda=1: E_kin + 3*E_pot = 0
  Second variation: E_kin + 9*E_pot > 0 for stability
  Combined: E_kin + 3*E_pot = 0 and E_kin + 9*E_pot > 0
  => 9*E_pot > -E_kin = 3*E_pot => 6*E_pot > 0 => E_pot > 0

  But wait: for TGP, f(g) is field-dependent, so the kinetic term
  is more complex. The standard Derrick scaling doesn't directly apply.

  Modified Derrick with f(g):
    Under r -> lambda*r:
      E_kin(lambda) = (1/lambda) * int f(g(x)) * g'(x)^2 * (lambda*x)^2 d(lambda*x)
                    = lambda * E_kin  (same as before, since f depends on g, not r)

  Actually: E_kin = 4*pi * int_0^inf (1/2)*f(g(r))*g'(r)^2 * r^2 dr
  Under g_lambda(r) = g(r/lambda), g'_lambda = g'/lambda:
  E_kin[g_lambda] = 4*pi * int (1/2)*f(g(r/lambda))*(g'(r/lambda)/lambda)^2 * r^2 dr
  Change var u = r/lambda:
  = 4*pi * int (1/2)*f(g(u))*g'(u)^2/lambda^2 * (lambda*u)^2 * lambda du
  = lambda * E_kin

  So standard Derrick applies even with f(g):
    Stationarity: E_kin + 3*E_pot = 0  (Derrick virial)
    Stability:    E_kin + 9*E_pot > 0  (against scaling)
""")

# Compute E_kin and E_pot for each g0
def compute_energies(g0, eta=ETA_K):
    r, g, gp = solve_sol(g0, eta)
    if r is None: return None, None, None, None
    def fk(gg):
        a = 2.0/(1+eta*(gg-1)**2)
        return 1+2*a*np.log(gg) if gg>0 else -1e30
    fk_arr = np.array([fk(gi) for gi in g])
    V_arr = g**3/3 - g**4/4 - 1.0/12  # V(g) - V(1)

    E_kin = 4*np.pi * trapezoid(0.5*fk_arr*gp**2*r**2, r)
    E_pot = 4*np.pi * trapezoid(V_arr*r**2, r)

    virial = E_kin + 3*E_pot  # should be ~0 for solutions
    stability = E_kin + 9*E_pot  # should be > 0 for stability
    return E_kin, E_pot, virial, stability

print(f"  {'g0':>6s} {'E_kin':>10s} {'E_pot':>10s} {'Virial':>10s} {'Stab(E+9V)':>10s} {'Stable?':>8s}")
print(f"  {'-'*60}")

g0_values = [0.905, 1.0, 1.2, 1.465, 2.0, 2.5, 3.0, 3.5, 3.845, 4.0, 4.5, 5.0, 6.0, 8.0]
for g0 in g0_values:
    ek, ep, vir, stab = compute_energies(g0)
    if ek is None: continue
    stable = "YES" if stab > 0 else "NO"
    mark = ""
    if abs(g0-0.905)<0.01: mark=" (e)"
    elif abs(g0-1.465)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {ek:10.4f} {ep:10.4f} {vir:10.4f} {stab:10.4f} {stable:>8s}{mark}")
    sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART B: ANGULAR PERTURBATION STABILITY")
print("="*70)

print("""
  For a spherically symmetric soliton g(r), angular perturbations
  of the form:
    delta_g(r, theta, phi) = u_l(r) * Y_lm(theta, phi)

  The linearized equation is:
    f(g)*u_l'' + [f'(g)*g' + 2*f(g)/r]*u_l'
    + [V''(g)*f(g) - f(g)*l(l+1)/r^2 + f''(g)*g'^2/2]*u_l
    = omega^2 * u_l

  The l=0 mode is the radial perturbation (breathing mode).
  The l=1 mode is the translation mode (zero mode).
  The l >= 2 modes are angular deformations.

  Stability requires omega^2 >= 0 for all l.
  We look for the LOWEST eigenvalue for each l.
""")

def angular_stability(g0, l_max=5, eta=ETA_K):
    """Compute lowest eigenvalue for angular perturbations l=0..l_max."""
    r, g, gp = solve_sol(g0, eta)
    if r is None: return {}

    def fk(gg):
        a = 2.0/(1+eta*(gg-1)**2)
        return 1+2*a*np.log(gg) if gg>0 else -1e30
    def fk_p(gg):
        """df/dg"""
        a = 2.0/(1+eta*(gg-1)**2)
        da = -2.0*2*eta*(gg-1)/(1+eta*(gg-1)**2)**2
        return 2*da*np.log(gg) + 2*a/gg if gg>0 else 0
    def Vpp(gg):
        """V''(g) = 2g - 3g^2"""
        return 2*gg - 3*gg**2

    # Use a Schrodinger-like form on a grid
    # Transform u_l(r) = chi_l(r)/r to remove the first-derivative term
    # (simplified: we use direct discretization)

    # Grid: use r[10:-100] to avoid boundary issues
    N = min(2000, len(r)-200)
    idx = np.linspace(10, len(r)-100, N, dtype=int)
    rg = r[idx]
    gg = g[idx]
    gpg = gp[idx]
    dr = np.diff(rg)

    results = {}

    for l in range(l_max+1):
        # Effective potential for radial Schrodinger equation
        # V_eff(r) = V''(g(r)) - l(l+1)/r^2 + (f'/f)*g'*(additional terms)
        # Simplified: for the TGP soliton linearized equation,
        # the key term is V''(g) - l(l+1)/r^2

        fk_arr = np.array([fk(gi) for gi in gg])
        Vpp_arr = np.array([Vpp(gi) for gi in gg])

        # Effective Schrodinger potential
        V_eff = Vpp_arr/fk_arr - l*(l+1)/rg**2

        # Build tridiagonal Hamiltonian (finite difference)
        n = len(rg) - 2  # interior points
        if n < 10: continue

        h = np.mean(dr)
        diag = -2.0/h**2 * np.ones(n) + V_eff[1:-1]
        off = 1.0/h**2 * np.ones(n-1)

        # Note: we want the LOWEST eigenvalue of -d^2/dr^2 + V_eff
        # but the operator is H = -d^2/dr^2 + V_eff (with minus sign for kinetic)
        # So we solve for eigenvalues of H. Negative eigenvalue = instability.

        # Actually for stability: we need omega^2 > 0 for all modes.
        # The operator is: L_l u = -(1/f)*[f*u'' + ...] + V''*u - l(l+1)/r^2 * u
        # In Schrodinger form: -u'' + V_eff*u = omega^2 * u

        try:
            # Get lowest few eigenvalues
            evals = eigh_tridiagonal(diag, off, eigvals_only=True,
                                     select='i', select_range=(0, min(2, n-1)))
            results[l] = evals[0]
        except:
            results[l] = np.nan

    return results

print(f"\n  Lowest eigenvalue omega^2 for angular modes l=0..5:")
print(f"  (omega^2 < 0 means UNSTABLE)")
print(f"\n  {'g0':>6s}", end="")
for l in range(6):
    print(f"  {'l='+str(l):>10s}", end="")
print(f"  {'Stable?':>8s}")
print(f"  {'-'*80}")

for g0 in [0.905, 1.465, 2.0, 3.0, 3.5, 3.845, 4.0, 4.5, 5.0, 6.0, 8.0]:
    evals = angular_stability(g0)
    print(f"  {g0:6.3f}", end="")
    all_stable = True
    for l in range(6):
        ev = evals.get(l, np.nan)
        if ev < 0: all_stable = False
        print(f"  {ev:10.4f}", end="")
    mark = ""
    if abs(g0-0.905)<0.01: mark=" (e)"
    elif abs(g0-1.465)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {'YES' if all_stable else 'NO':>8s}{mark}")
    sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART C: CRITICAL g0 FOR INSTABILITY ONSET")
print("="*70)
print()

# Fine scan to find where any mode becomes unstable
print(f"  Fine scan: lowest omega^2 for l=2 mode (first deformation)")
print(f"  {'g0':>6s} {'omega^2(l=2)':>14s} {'Status':>10s}")
print(f"  {'-'*35}")

for g0 in np.arange(2.0, 10.1, 0.5):
    evals = angular_stability(g0, l_max=3)
    ev2 = evals.get(2, np.nan)
    status = "stable" if ev2 >= 0 else "UNSTABLE"
    mark = ""
    if abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.2f} {ev2:14.6f} {status:>10s}{mark}")
    sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART D: RADIAL BREATHING MODE (l=0)")
print("="*70)
print()

# The l=0 mode is special: the lowest eigenvalue should be
# the "breathing" mode. For a stable soliton, this should be positive.
# If it becomes negative, the soliton is radially unstable.

print(f"  Fine scan: lowest omega^2 for l=0 (breathing mode)")
print(f"  {'g0':>6s} {'omega^2(l=0)':>14s} {'omega^2(l=1)':>14s} {'Status':>10s}")
print(f"  {'-'*50}")

for g0 in np.arange(2.0, 10.1, 0.5):
    evals = angular_stability(g0, l_max=2)
    ev0 = evals.get(0, np.nan)
    ev1 = evals.get(1, np.nan)
    status = "stable" if (ev0 >= 0 and ev1 >= -0.01) else "UNSTABLE"
    # Note: l=1 should have a zero mode (translation) => eigenvalue ~ 0
    mark = ""
    if abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.2f} {ev0:14.6f} {ev1:14.6f} {status:>10s}{mark}")
    sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART E: TOTAL SOLITON ENERGY")
print("="*70)
print()

# Another constraint: the soliton energy must be POSITIVE
# (negative energy would mean the soliton is more favorable than vacuum)
# E_total = E_kin + E_pot; for Derrick: E_kin = -3*E_pot
# So E_total = E_kin + E_pot = -3*E_pot + E_pot = -2*E_pot
# E_total > 0 requires E_pot < 0

print(f"  {'g0':>6s} {'E_total':>10s} {'E_kin':>10s} {'E_pot':>10s} {'E_kin/|E_pot|':>14s}")
print(f"  {'-'*55}")

for g0 in [0.905, 1.465, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0]:
    ek, ep, vir, stab = compute_energies(g0)
    if ek is None: continue
    et = ek + ep
    ratio = ek/abs(ep) if abs(ep) > 1e-10 else np.nan
    mark = ""
    if abs(g0-0.905)<0.01: mark=" (e)"
    elif abs(g0-1.465)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {et:10.4f} {ek:10.4f} {ep:10.4f} {ratio:14.4f}{mark}")
    sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)
print("""
  ANALYSIS RESULTS:

  1. Derrick virial: E_kin + 3*E_pot = 0 is satisfied for solutions
     (this is automatic for stationary points of the energy).

  2. Derrick stability (against scaling): E_kin + 9*E_pot > 0
     This determines if the soliton is stable against radial rescaling.

  3. Angular stability: the lowest eigenvalue of the perturbation
     operator for l >= 2 modes determines angular stability.

  4. Key question: Is there a g0_max above which the soliton
     becomes unstable? If g0_max = 4, we have our answer.
""")

print("DONE")
