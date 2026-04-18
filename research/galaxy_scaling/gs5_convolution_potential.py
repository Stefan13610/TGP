#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs5_convolution_potential.py: Averaged N-soliton potential via convolution.

THE SELF-TERM PROBLEM (from gs4):
  |∇g|²/g includes self-terms |∇δᵢ|² which are localized at each soliton.
  These are already part of the measured soliton mass → can't count them
  as "extra dark matter."

  The CORRECT approach: compute the N-soliton potential DIRECTLY and
  extract the rotation curve. The soliton profile g(r) already includes
  the nonlinear g'²/g term — it's the FULL TGP solution.

CLEANER METHOD:
  For a spherically symmetric galaxy with soliton density n(r'):
    δ_total(r) = ∫ n(r') × δ_soliton(|r - r'|) d³r'

  This is a 3D CONVOLUTION of n(r) with the soliton profile δ(r).
  Due to spherical symmetry, reduces to 1D integral.

  The key: the soliton tail δ ~ A·sin(r)/r oscillates.
  When convolved with a smooth density, the oscillations partially cancel.
  The RESIDUAL after cancellation determines the rotation curve.

  Compare:
  A) Newtonian soliton (δ ~ 1/r, no oscillation): gives standard 1/r potential
  B) TGP soliton (δ ~ sin(r)/r): gives MODIFIED potential from convolution

  If the TGP convolution gives a potential that falls slower than 1/r,
  the rotation curve will be flatter than Keplerian.

PLAN:
  1. Compute single soliton profile (from ODE)
  2. Spherical convolution with galaxy density n(r)
  3. Compare TGP vs Newtonian potential and rotation curve
  4. Multiple galaxy profiles (Hernquist, isothermal, exponential)
  5. Scale analysis: when does TGP deviate from Newton?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, quad, simpson
from scipy.interpolate import interp1d
import time

print("="*78)
print("  CONVOLUTION POTENTIAL: N-soliton averaged field")
print("="*78)

# ==========================================================================
# STEP 1: Soliton profile
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 1: Single soliton profile (TGP vs Newtonian)")
print(f"{'='*78}")

def tgp_rhs(r, y):
    g, gp = y
    if g <= 1e-10:
        return [gp, 0]
    if r < 1e-10:
        gpp = (1 - g - gp**2/g) / 3
    else:
        gpp = -gp**2/g - 2*gp/r - g + 1
    return [gp, gpp]

g0 = 1.1  # Weak soliton (galaxy particle in weak field)
r_max_sol = 300
sol = solve_ivp(tgp_rhs, [1e-6, r_max_sol], [g0, 0],
                method='RK45', max_step=0.02, rtol=1e-11, atol=1e-13,
                dense_output=True)

r_s = sol.t
g_s = sol.y[0]
delta_s = g_s - 1
gp_s = sol.y[1]

# Create fine interpolation
delta_func = interp1d(r_s, delta_s, kind='cubic', fill_value=0, bounds_error=False)
gp_func = interp1d(r_s, gp_s, kind='cubic', fill_value=0, bounds_error=False)

# Also create "Newtonian equivalent": same core mass, but 1/r tail (no oscillation)
# Newtonian enclosed mass = -r²·g' evaluated at large r
# For the soliton, M_enc(r) oscillates. The "asymptotic Newtonian mass" is the
# time-average (phase-average) of M_enc over one oscillation period
mask_far = r_s > 20
M_enc_far = -r_s[mask_far]**2 * gp_s[mask_far]
M_newton_avg = np.mean(M_enc_far)  # average over oscillations

# Newtonian profile: δ_N(r) = M_newton_avg / r (1/r potential)
delta_newton = lambda r: M_newton_avg / np.maximum(r, 0.01)

# Check the tail
print(f"  g₀ = {g0}, δ₀ = {g0-1:.4f}")
print(f"  Asymptotic Newtonian mass: M_avg = {M_newton_avg:.4f}")
print(f"\n  {'r':>6s} {'delta_TGP':>10s} {'delta_Newton':>12s} {'ratio':>8s} {'r*delta_TGP':>11s}")
print(f"  {'---':>6s} {'---':>10s} {'---':>12s} {'---':>8s} {'---':>11s}")
for r in [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 80, 100, 150, 200]:
    dt = delta_func(r)
    dn = delta_newton(r)
    rat = dt/dn if abs(dn) > 1e-15 else 0
    print(f"  {r:6.0f} {dt:10.5f} {dn:12.5f} {rat:8.3f} {r*dt:11.5f}")

# ==========================================================================
# STEP 2: Spherical convolution integral
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 2: Spherical convolution — math")
print(f"{'='*78}")

print("""
  For a spherically symmetric density n(r') and soliton profile delta(s):

    delta_total(r) = integral n(r') * <delta(|r-r'|)>_angle * 4pi r'^2 dr'

  Angular average reduces to radial integral.
  For Newtonian delta=M/s: shell theorem gives M_enc/r.
  For TGP delta=A*sin(s)/s: angular average gives 2A*sin(r)*sin(r')/(r*r').

  REMARKABLE: oscillation is PRESERVED after angular averaging!
  Total field ~ (sin(r)/r) * integral(n(r')*sin(r')*r' dr')

  Real soliton has core at small s. Computing numerically...
""")

# ==========================================================================
# STEP 3: Numerical convolution
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 3: Numerical convolution for various galaxy profiles")
print(f"{'='*78}")

def angular_avg_delta(r, rprime, delta_func_local, n_quad=200):
    """
    Compute the angular average of delta(|r - r'|) for given r and r'.
    <delta> = (1/(2*r*r')) * integral_{|r-r'|}^{r+r'} delta(s) * s ds
    """
    if r < 1e-6 or rprime < 1e-6:
        return delta_func_local(abs(r - rprime))

    s_min = abs(r - rprime)
    s_max = r + rprime

    s_pts = np.linspace(s_min + 1e-8, s_max, n_quad)
    integrand = delta_func_local(s_pts) * s_pts
    integral = np.trapezoid(integrand, s_pts)

    return integral / (2 * r * rprime)


def compute_convolved_potential(r_eval, n_func, r_max_galaxy, delta_func_local,
                                 n_rprime=200):
    """
    Compute delta_total(r) = integral n(r') * <delta(|r-r'|)>_angle * 4*pi*r'^2 dr'

    n_func: number density as function of r'
    r_max_galaxy: maximum radius of galaxy
    delta_func_local: single soliton profile
    """
    delta_total = np.zeros_like(r_eval)

    rprime_pts = np.linspace(0.1, r_max_galaxy, n_rprime)
    dr = rprime_pts[1] - rprime_pts[0]

    for i, r in enumerate(r_eval):
        integrand = np.zeros(n_rprime)
        for j, rp in enumerate(rprime_pts):
            avg_delta = angular_avg_delta(r, rp, delta_func_local)
            integrand[j] = n_func(rp) * avg_delta * 4 * np.pi * rp**2
        delta_total[i] = np.trapezoid(integrand, rprime_pts)

    return delta_total


# Galaxy profiles
def hernquist_density(r, N_total, a):
    """Hernquist profile: n(r) = N/(2*pi) * a / (r * (r+a)^3)"""
    return N_total / (2 * np.pi) * a / (r * (r + a)**3)

def isothermal_density(r, N_total, r_max):
    """Isothermal: n(r) = n0 / r^2, truncated at r_max"""
    # N_total = 4*pi*n0*r_max → n0 = N_total/(4*pi*r_max)
    n0 = N_total / (4 * np.pi * r_max)
    return np.where(r < r_max, n0 / np.maximum(r, 0.1)**2, 0)

def uniform_density(r, N_total, R):
    """Uniform sphere: n = const for r < R"""
    n0 = N_total / (4/3 * np.pi * R**3)
    return np.where(r < R, n0, 0)


# Parameters
N_total = 1000  # number of solitons
R_galaxy = 30   # soliton units

# Evaluation points
r_eval = np.linspace(1, 4 * R_galaxy, 60)

print(f"  Computing convolutions (N={N_total}, R_gal={R_galaxy})...")
print(f"  This takes some time (numerical integration)...\n")

# Compute for TGP and Newton, for Hernquist profile
n_hern = lambda r: hernquist_density(r, N_total, R_galaxy / 3)  # a = R/3

t0 = time.time()
delta_TGP = compute_convolved_potential(r_eval, n_hern, 3*R_galaxy, delta_func, n_rprime=100)
t1 = time.time()
delta_N = compute_convolved_potential(r_eval, n_hern, 3*R_galaxy, delta_newton, n_rprime=100)
t2 = time.time()

print(f"  TGP convolution: {t1-t0:.1f}s")
print(f"  Newton convolution: {t2-t1:.1f}s")

# Rotation curves: v^2 = -r * d(delta)/dr (since delta = 2*Phi/c^2, v^2 = r*a = -r*Phi')
# Numerical derivative
ddelta_dr_TGP = np.gradient(delta_TGP, r_eval)
ddelta_dr_N = np.gradient(delta_N, r_eval)

v2_TGP = -r_eval * ddelta_dr_TGP / 2  # factor 1/2 from delta = 2*Phi/c^2
v2_N = -r_eval * ddelta_dr_N / 2

v_TGP = np.sqrt(np.maximum(v2_TGP, 0))
v_N = np.sqrt(np.maximum(v2_N, 0))

# ==========================================================================
# STEP 4: Results — TGP vs Newton
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 4: Rotation curves — TGP vs Newton (Hernquist galaxy)")
print(f"{'='*78}")

# Normalize to peak velocity
v_peak_N = np.max(v_N[r_eval > 3])
v_peak_TGP = np.max(v_TGP[r_eval > 3])

print(f"\n  Peak v_Newton = {v_peak_N:.4f} at r = {r_eval[np.argmax(v_N)]:.1f}")
print(f"  Peak v_TGP    = {v_peak_TGP:.4f} at r = {r_eval[np.argmax(v_TGP)]:.1f}")

print(f"\n  {'r':>6s} {'delta_TGP':>10s} {'delta_N':>10s} {'v_TGP':>8s} {'v_Newton':>8s} {'v_T/v_N':>7s} {'v_TGP/peak':>10s} {'v_N/peak':>8s}")
print(f"  {'---'*30}")

for i in range(0, len(r_eval), max(1, len(r_eval)//25)):
    r = r_eval[i]
    dt = delta_TGP[i]
    dn = delta_N[i]
    vt = v_TGP[i]
    vn = v_N[i]
    rat = vt/vn if vn > 0.001 else 0
    vt_norm = vt/v_peak_TGP if v_peak_TGP > 0 else 0
    vn_norm = vn/v_peak_N if v_peak_N > 0 else 0
    print(f"  {r:6.1f} {dt:10.5f} {dn:10.5f} {vt:8.4f} {vn:8.4f} {rat:7.3f} {vt_norm:10.3f} {vn_norm:8.3f}")

# ==========================================================================
# STEP 5: Flatness analysis
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 5: Flatness analysis")
print(f"{'='*78}")

# Power law fit to v(r) in outer region (r > R_galaxy)
mask_outer = (r_eval > R_galaxy) & (v_TGP > 0) & (v_N > 0)
if np.sum(mask_outer) > 5:
    log_r_out = np.log(r_eval[mask_outer])
    log_v_TGP_out = np.log(v_TGP[mask_outer])
    log_v_N_out = np.log(v_N[mask_outer])

    # Remove any NaN/inf
    valid_T = np.isfinite(log_v_TGP_out)
    valid_N = np.isfinite(log_v_N_out)

    if np.sum(valid_T) > 3:
        alpha_TGP = np.polyfit(log_r_out[valid_T], log_v_TGP_out[valid_T], 1)[0]
    else:
        alpha_TGP = float('nan')

    if np.sum(valid_N) > 3:
        alpha_N = np.polyfit(log_r_out[valid_N], log_v_N_out[valid_N], 1)[0]
    else:
        alpha_N = float('nan')

    print(f"  Power law fit v ∝ r^alpha for r > {R_galaxy} (outer region):")
    print(f"  TGP:    alpha = {alpha_TGP:.3f}")
    print(f"  Newton: alpha = {alpha_N:.3f}")
    print(f"  Keplerian: alpha = -0.500")
    print(f"  Flat:      alpha =  0.000")
    print(f"\n  TGP is {'FLATTER' if alpha_TGP > alpha_N else 'STEEPER'} than Newton by {alpha_TGP - alpha_N:+.3f}")

    if alpha_TGP > -0.2:
        print(f"  → NEARLY FLAT rotation curve! TGP mechanism works!")
    elif alpha_TGP > -0.4:
        print(f"  → SIGNIFICANTLY flatter than Keplerian")
    else:
        print(f"  → Still close to Keplerian")
else:
    print(f"  Not enough valid data points for fit")
    alpha_TGP = float('nan')
    alpha_N = float('nan')

# ==========================================================================
# STEP 6: Different galaxy profiles
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 6: Different galaxy density profiles")
print(f"{'='*78}")

profiles = [
    ("Hernquist",   lambda r: hernquist_density(r, N_total, R_galaxy/3)),
    ("Isothermal",  lambda r: isothermal_density(r, N_total, R_galaxy)),
    ("Uniform",     lambda r: uniform_density(r, N_total, R_galaxy)),
]

for prof_name, n_func in profiles:
    print(f"\n  Profile: {prof_name}")

    delta_T = compute_convolved_potential(r_eval, n_func, 3*R_galaxy, delta_func, n_rprime=80)
    delta_Nw = compute_convolved_potential(r_eval, n_func, 3*R_galaxy, delta_newton, n_rprime=80)

    dd_T = np.gradient(delta_T, r_eval)
    dd_N = np.gradient(delta_Nw, r_eval)

    v2_T = -r_eval * dd_T / 2
    v2_Nw = -r_eval * dd_N / 2

    vT = np.sqrt(np.maximum(v2_T, 0))
    vNw = np.sqrt(np.maximum(v2_Nw, 0))

    # Power law in outer region
    mask_o = (r_eval > R_galaxy) & (vT > 0) & (vNw > 0)
    if np.sum(mask_o) > 5:
        lr = np.log(r_eval[mask_o])
        lv_T = np.log(vT[mask_o])
        lv_N = np.log(vNw[mask_o])

        valid_T = np.isfinite(lv_T)
        valid_N = np.isfinite(lv_N)

        aT = np.polyfit(lr[valid_T], lv_T[valid_T], 1)[0] if np.sum(valid_T)>3 else float('nan')
        aN = np.polyfit(lr[valid_N], lv_N[valid_N], 1)[0] if np.sum(valid_N)>3 else float('nan')

        vp_T = np.max(vT[r_eval > 3])
        vp_N = np.max(vNw[r_eval > 3])
        v_ratio_2R = vT[np.argmin(np.abs(r_eval-2*R_galaxy))] / vp_T if vp_T > 0 else 0
        v_ratio_2R_N = vNw[np.argmin(np.abs(r_eval-2*R_galaxy))] / vp_N if vp_N > 0 else 0

        print(f"    alpha_TGP = {aT:.3f}, alpha_Newton = {aN:.3f}, Delta = {aT-aN:+.3f}")
        print(f"    v(2R)/v_peak: TGP = {v_ratio_2R:.3f}, Newton = {v_ratio_2R_N:.3f}")
    else:
        print(f"    Insufficient data for power law fit")

# ==========================================================================
# STEP 7: Scale dependence — vary R_galaxy at fixed N
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 7: Scale dependence — how does flattening depend on R_galaxy?")
print(f"{'='*78}")

N_fix = 500
R_values = [5, 10, 20, 40, 60]

print(f"\n  Fixed N = {N_fix}")
print(f"  {'R_gal':>6s} {'alpha_TGP':>10s} {'alpha_N':>10s} {'Delta':>8s} {'R/lambda':>8s}")
print(f"  {'---'*20}")

for R_g in R_values:
    r_ev = np.linspace(1, 3*R_g, 40)
    n_h = lambda r, Rg=R_g: hernquist_density(r, N_fix, Rg/3)

    dT = compute_convolved_potential(r_ev, n_h, 2*R_g, delta_func, n_rprime=60)
    dN = compute_convolved_potential(r_ev, n_h, 2*R_g, delta_newton, n_rprime=60)

    ddT = np.gradient(dT, r_ev)
    ddN = np.gradient(dN, r_ev)

    v2T = -r_ev * ddT / 2
    v2N = -r_ev * ddN / 2

    vT = np.sqrt(np.maximum(v2T, 0))
    vN = np.sqrt(np.maximum(v2N, 0))

    mask_o = (r_ev > R_g) & (vT > 0) & (vN > 0)
    if np.sum(mask_o) > 5:
        lr = np.log(r_ev[mask_o])
        aT = np.polyfit(lr[np.isfinite(np.log(vT[mask_o]))],
                        np.log(vT[mask_o])[np.isfinite(np.log(vT[mask_o]))], 1)[0]
        aN = np.polyfit(lr[np.isfinite(np.log(vN[mask_o]))],
                        np.log(vN[mask_o])[np.isfinite(np.log(vN[mask_o]))], 1)[0]
        delta_alpha = aT - aN
        R_over_lam = R_g / (2*np.pi)
        print(f"  {R_g:6d} {aT:10.3f} {aN:10.3f} {delta_alpha:+8.3f} {R_over_lam:8.2f}")
    else:
        print(f"  {R_g:6d}     ---        ---      ---   {R_g/(2*np.pi):8.2f}")

# ==========================================================================
# STEP 8: The key comparison — potential envelope
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 8: Potential envelope — does TGP fall slower than 1/r?")
print(f"{'='*78}")

# For the main Hernquist galaxy (N=1000, R=30):
# Compare |delta_TGP(r)| vs |delta_Newton(r)| at large r

print(f"\n  Hernquist galaxy, N={N_total}, R_gal={R_galaxy}")
print(f"  Potential at large r (normalized to r=R_gal):")
print(f"\n  {'r':>6s} {'|dT|':>10s} {'|dN|':>10s} {'|dT/dN|':>8s} {'r*|dT|':>10s} {'r*|dN|':>10s} {'r*ratio':>8s}")
print(f"  {'---'*25}")

for i in range(0, len(r_eval), max(1, len(r_eval)//20)):
    r = r_eval[i]
    dt = abs(delta_TGP[i])
    dn = abs(delta_N[i])
    rat = dt/dn if dn > 1e-10 else 0
    print(f"  {r:6.1f} {dt:10.5f} {dn:10.5f} {rat:8.3f} {r*dt:10.5f} {r*dn:10.5f} {rat:8.3f}")

# Check: does r*delta_TGP stay constant (1/r potential) or grow (slower than 1/r)?
r_far = r_eval[r_eval > R_galaxy]
if len(r_far) > 3:
    rd_T = r_far * np.abs(delta_TGP[r_eval > R_galaxy])
    rd_N = r_far * np.abs(delta_N[r_eval > R_galaxy])

    # Fit r*|delta| ∝ r^β → |delta| ∝ r^(β-1)
    # β = 0: 1/r potential (Newtonian)
    # β > 0: slower than 1/r → flatter rotation curve
    valid_T = rd_T > 0
    valid_N = rd_N > 0

    if np.sum(valid_T) > 3:
        beta_T = np.polyfit(np.log(r_far[valid_T]), np.log(rd_T[valid_T]), 1)[0]
    else:
        beta_T = 0

    if np.sum(valid_N) > 3:
        beta_N = np.polyfit(np.log(r_far[valid_N]), np.log(rd_N[valid_N]), 1)[0]
    else:
        beta_N = 0

    print(f"\n  Potential envelope fit: |delta| ∝ r^(beta-1)")
    print(f"  TGP:    beta = {beta_T:.3f} → |delta| ∝ r^{{{beta_T-1:.3f}}}")
    print(f"  Newton: beta = {beta_N:.3f} → |delta| ∝ r^{{{beta_N-1:.3f}}}")
    print(f"  Newtonian: beta = 0 → |delta| ∝ r^{{-1}}")

    if beta_T > beta_N + 0.05:
        print(f"\n  → TGP potential decays SLOWER than Newton by r^{{{beta_T-beta_N:.3f}}}")
        print(f"  → This produces FLATTER rotation curves ✓")
    elif beta_T > beta_N - 0.05:
        print(f"\n  → TGP and Newton have SIMILAR decay rates")
    else:
        print(f"\n  → TGP potential decays FASTER than Newton (unexpected)")

# ==========================================================================
# STEP 9: SUMMARY
# ==========================================================================
print(f"\n{'='*78}")
print(f"  SUMMARY")
print(f"{'='*78}")

print(f"""
  METHOD: Compute the averaged potential of a galaxy modeled as
  N solitons with TGP profile delta(r) ~ sin(r)/r.
  Compare with Newtonian (delta ~ 1/r).

  The convolution of the soliton profile with a smooth galaxy density
  gives the COLLECTIVE potential at each radius.

  KEY FINDINGS:

  1. POTENTIAL PROFILE:
     TGP: |delta| ∝ r^{{{beta_T-1:.3f}}} (outer region)
     Newton: |delta| ∝ r^{{{beta_N-1:.3f}}} (outer region)
     Difference in decay: {beta_T - beta_N:+.3f}

  2. ROTATION CURVE:
     TGP alpha = {alpha_TGP:.3f}
     Newton alpha = {alpha_N:.3f}
     (Keplerian = -0.500, flat = 0.000)

  3. PHYSICAL MECHANISM:
     The oscillating soliton tail sin(r)/r has a DIFFERENT convolution
     behavior than the monotonic 1/r Newtonian potential.

     For 1/r: convolution with smooth density → shell theorem → M_enc/r
     For sin(r)/r: convolution preserves oscillatory structure
     → After angular averaging: ~ sin(r)·integral(n(r')·sin(r')·r' dr')

     The sin(r)·sin(r') factor creates COHERENT constructive interference
     at scales r ≈ n·pi (soliton tail resonances).

  4. SCALING:
     The TGP effect depends on R_galaxy/lambda_tail.
     When R >> lambda (= 2*pi ≈ 6.3 soliton units):
     many soliton oscillation periods fit inside → averaging smooths
     but systematic ENVELOPE remains.
""")
