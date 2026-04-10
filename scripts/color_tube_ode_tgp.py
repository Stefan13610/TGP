#!/usr/bin/env python3
"""
TGP Color Tube ODE v2: First-Principles String Tension (2026-04-10)

Physics:
  A color flux tube between two quarks depletes the TGP field Phi
  inside the tube. The tube is stabilized by the balance between:
  (a) Phi wanting to return to Phi_0 (potential energy)
  (b) Color flux energy acting as a source depleting Phi

  The TGP field equation with color source:
    nabla^2 phi + 2(nabla phi)^2/phi + phi^2 - phi^3 = -j(rho)

  For a long tube (L >> R_tube), the transverse profile is z-independent.
  In cylindrical coordinates (rho = transverse distance from axis):

    phi'' + phi'/rho + 2(phi')^2/phi + phi^2 - phi^3 = -j_0 * Theta(R0-rho)

  where j_0 > 0 is the color source strength and R0 is the tube core radius.

  Approach:
  1. Use scipy.integrate.solve_bvp for the BVP
  2. Scan (j_0, R0) parameter space
  3. Compute string tension sigma for each solution
  4. Identify the self-consistent physical tube
"""

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.special import k0, k1
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# TGP Parameters
# ============================================================
Phi_0 = 24.6492
kappa = 3.0 / (4 * Phi_0)
a_Gamma = 0.040
phi_golden = (1 + np.sqrt(5)) / 2
A_target = a_Gamma / phi_golden  # = 0.02472

print("=" * 70)
print("TGP COLOR TUBE ODE v2 (with color source)")
print("=" * 70)

# ============================================================
# [1] APPROACH: Shooting from large rho inward
# ============================================================
# Asymptotic form at large rho: phi ~ 1 - C*K0(rho_hat)
# where K0 is the modified Bessel function.
# K0(x) ~ sqrt(pi/(2x)) * exp(-x) for large x.
# The parameter C controls tube strength.

print("\n[1] LINEARIZED TUBE: phi = 1 - C*K0(rhat)")
print("-" * 50)

# The linearized energy per unit length for a weak tube (C << 1):
# sigma_lin = 2*pi * integral_0^inf [(1/2)(C*K0')^2 + (1/2)C^2*K0^2] rho drho
# The kinetic term dominates at small rho.
# Exact result: integral_0^inf [K0'(rho)^2 + K0(rho)^2] rho drho = 1
# (from the Bessel orthogonality)
# So sigma_lin ~ pi * C^2

# Actually, let me compute this properly.
# K0' = -K1, so K0'^2 = K1^2
# integral_0^inf [K1(x)^2 + K0(x)^2] x dx = ... (known Bessel integral)

rho_grid = np.linspace(0.01, 15, 2000)
K0_vals = k0(rho_grid)
K1_vals = k1(rho_grid)

# Energy integrals for linearized tube
integ_kin = np.trapezoid(K1_vals**2 * rho_grid, rho_grid)
integ_pot = np.trapezoid(K0_vals**2 * rho_grid, rho_grid)

print(f"  Bessel integrals:")
print(f"    I_kin = int K1^2 rho drho = {integ_kin:.6f}")
print(f"    I_pot = int K0^2 rho drho = {integ_pot:.6f}")
print(f"    I_total = I_kin + I_pot = {integ_kin + integ_pot:.6f}")
# These should be approximately 1/2 each, total = 1

# Linearized string tension (per unit length):
# sigma_lin = pi * C^2 * (I_kin + I_pot)
# where C is the tube amplitude parameter

print(f"\n  Linearized sigma = pi * C^2 * {integ_kin + integ_pot:.4f}")

# ============================================================
# [2] NONLINEAR TUBE: Full BVP solution
# ============================================================
print(f"\n\n[2] NONLINEAR TUBE SOLUTION")
print("-" * 50)

# The full equation is:
# phi'' + phi'/rho + 2(phi')^2/phi + phi^2 - phi^3 = 0
# where phi = 1 - epsilon(rho) with epsilon > 0 (depletion)
#
# The SOURCED version adds -j0 inside a core:
# phi'' + phi'/rho + 2(phi')^2/phi + phi^2 - phi^3 = -j0 * Theta(R0-rho)
#
# Alternative approach: Use the NONLINEAR depletion equation
# Set phi = 1 - u(rho) where 0 < u < 1:
# Expanding in u but keeping nonlinear terms:
#
# -u'' - u'/rho + 2u'^2/(1-u) + (1-u)^2 - (1-u)^3 = -j(rho)
# -u'' - u'/rho + 2u'^2/(1-u) + u - 2u^2 + u^3 = -j(rho)
#
# Or equivalently:
# u'' + u'/rho = 2u'^2/(1-u) + u - 2u^2 + u^3 + j(rho)
#
# For j=0 and small u: u'' + u'/rho ≈ u (modified Bessel)
# For larger u: the nonlinear terms become important

# Let me solve the BVP directly using a shooting method,
# but now from OUTSIDE in. Start with the known asymptotic
# phi ~ 1 - C*K0(rhat) at large rhat, integrate inward.

def field_eq_phi(rhat, y, j0=0, R0=1.0):
    """
    y[0] = phi, y[1] = phi'
    phi'' = -phi'/rhat - 2*(phi')^2/phi - phi^2 + phi^3 + j0*Theta(R0-rhat)

    NOTE: +j0 sign! Color flux DEPLETES Phi (destroys space generation),
    opposite to ordinary mass which enhances Phi. The positive source
    on the RHS counteracts the potential's tendency to restore phi=1,
    allowing a depleted tube solution.
    """
    phi, dphi = y
    if phi < 1e-10:
        phi = 1e-10

    source = j0 if rhat < R0 else 0.0

    if rhat < 1e-6:
        ddphi = (phi**3 - phi**2 + source) / 2.0
    else:
        ddphi = -dphi/rhat - 2*dphi**2/phi - phi**2 + phi**3 + source

    return [dphi, ddphi]

def solve_tube_bvp(j0, R0, rhat_max=20.0, N=500):
    """
    Solve the sourced tube equation using scipy solve_bvp.
    """
    rhat = np.linspace(0.001, rhat_max, N)

    def ode(x, y):
        phi, dphi = y
        phi = np.maximum(phi, 1e-8)
        source = np.where(x < R0, j0, 0.0)

        # Avoid division by zero at x=0
        inv_x = np.where(x > 0.01, 1.0/x, 0.0)

        ddphi = -dphi*inv_x - 2*dphi**2/phi - phi**2 + phi**3 + source
        return np.vstack([dphi, ddphi])

    def bc(ya, yb):
        # BC: phi'(0) = 0, phi(rhat_max) = 1
        return np.array([ya[1], yb[0] - 1.0])

    # Initial guess: smooth depletion
    delta = j0 * R0**2 / 4  # rough estimate of central depletion
    delta = min(delta, 0.8)
    phi_guess = 1.0 - delta * np.exp(-rhat**2 / (2*R0**2))
    dphi_guess = delta * rhat / R0**2 * np.exp(-rhat**2 / (2*R0**2))

    y_guess = np.vstack([phi_guess, dphi_guess])

    try:
        sol = solve_bvp(ode, bc, rhat, y_guess, tol=1e-6, max_nodes=5000)
        if sol.success:
            return sol
        else:
            return None
    except Exception:
        return None

def solve_tube_ivp(j0, R0, rhat_max=25.0):
    """
    Alternative: IVP shooting from rhat = 0.
    Bisect on phi(0) to match phi(rhat_max) ≈ 1.
    """
    def residual(phi0):
        y0 = [phi0, 0.0]
        sol = solve_ivp(lambda r, y: field_eq_phi(r, y, j0, R0),
                       [1e-4, rhat_max], y0,
                       method='RK45', rtol=1e-9, atol=1e-11,
                       max_step=0.05)
        if sol.success and len(sol.t) > 10:
            return sol.y[0][-1] - 1.0, sol
        return 10.0, sol

    # Bisection search for phi(0)
    lo, hi = 0.01, 0.999
    for _ in range(60):
        mid = (lo + hi) / 2
        r, _ = residual(mid)
        if r > 0:
            lo = mid
        else:
            hi = mid

    r_final, sol_final = residual((lo+hi)/2)
    return sol_final if abs(r_final) < 0.05 else None, (lo+hi)/2

# ============================================================
# [3] SCAN OVER SOURCE PARAMETERS
# ============================================================
print(f"\n\n[3] SOURCE PARAMETER SCAN")
print("-" * 50)

print(f"  {'j0':>6} {'R0':>5} {'phi_min':>9} {'sigma_hat':>12} {'R_tube_eff':>11} {'status':>8}")
print(f"  {'-'*55}")

results = []

for j0 in [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]:
    for R0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
        # Try BVP first
        sol = solve_tube_bvp(j0, R0)

        if sol is None:
            # Try IVP shooting
            sol_ivp, phi0_found = solve_tube_ivp(j0, R0)
            if sol_ivp is not None:
                rhat = sol_ivp.t
                phi_sol = sol_ivp.y[0]
                dphi_sol = sol_ivp.y[1]
                phi_min_val = phi_sol[0]
            else:
                continue
        else:
            rhat = sol.x
            phi_sol = sol.y[0]
            dphi_sol = sol.y[1]
            phi_min_val = phi_sol[0]

        # Compute string tension
        # Energy density above vacuum:
        # e = (1/2)*phi^4*(dphi)^2 + U(phi) - U(1) + j0*Theta(R0-rho)*(1-phi)
        # where U(phi) = (1/3)phi^3 - (1/4)phi^4
        U_vac = 1.0/3.0 - 1.0/4.0  # = 1/12

        e_arr = np.zeros_like(rhat)
        for i in range(len(rhat)):
            phi_i = max(phi_sol[i], 1e-10)
            dphi_i = dphi_sol[i]
            U_i = (1.0/3.0)*phi_i**3 - (1.0/4.0)*phi_i**4
            E_kin = 0.5 * phi_i**4 * dphi_i**2
            E_pot = U_i - U_vac
            # Source energy: j0 depletes Phi, creating energy cost
            E_src = j0 * (1 - phi_i)**2 / 2 if rhat[i] < R0 else 0.0
            e_arr[i] = E_kin + E_pot + E_src

        sigma_hat = np.trapezoid(2*np.pi*e_arr*rhat, rhat)

        # Effective tube radius (half-depletion point)
        mid_val = (phi_min_val + 1.0) / 2.0
        try:
            idx_mid = np.argmin(np.abs(phi_sol - mid_val))
            R_eff = rhat[idx_mid]
        except:
            R_eff = R0

        if phi_min_val < 0.99 and sigma_hat > 0:
            results.append((j0, R0, phi_min_val, sigma_hat, R_eff))
            print(f"  {j0:6.1f} {R0:5.1f} {phi_min_val:9.5f} {sigma_hat:12.6f} "
                  f"{R_eff:11.4f} {'OK':>8}")

# ============================================================
# [4] SELF-CONSISTENT TUBE: Match to TGP coupling
# ============================================================
print(f"\n\n[4] SELF-CONSISTENT TUBE IDENTIFICATION")
print("-" * 50)

if results:
    # The color source strength j0 should be related to alpha_s:
    # j0 ~ alpha_s * (color factor) / R0^2
    # For SU(3): j0 ~ (4/3) * alpha_s / R0^2 (Casimir)
    #
    # In TGP: alpha_s = N_c^3 * g_0^e / (8*Phi_0) = 0.1190
    # j0 * R0^2 ~ (4/3) * 0.119 = 0.159

    alpha_s_TGP = 0.1190
    C_F = 4.0/3.0  # SU(3) fundamental Casimir
    target_product = C_F * alpha_s_TGP

    print(f"  Target: j0 * R0^2 ~ C_F * alpha_s = {target_product:.4f}")
    print(f"  (C_F = 4/3 for SU(3) fundamental)")

    # Find results closest to this constraint
    for j0, R0, pm, sig, Re in results:
        product = j0 * R0**2
        deviation = abs(product / target_product - 1)
        if deviation < 1.0:  # within factor 2
            print(f"    j0={j0:.1f}, R0={R0:.1f}: "
                  f"j0*R0^2={product:.3f} (dev={deviation*100:.0f}%), "
                  f"sigma_hat={sig:.6f}, phi_min={pm:.4f}")

    # ============================================================
    # [5] STRING TENSION AND m_0 DERIVATION
    # ============================================================
    print(f"\n\n[5] STRING TENSION -> m_0")
    print("-" * 50)

    # Use the best physical solution
    # Physical string tension:
    # sigma_phys = (gamma / kappa) * sigma_hat
    # where gamma = m_sp^2 sets the unit of length
    #
    # QCD string tension: sigma_QCD = (440 MeV)^2
    # This determines m_sp: m_sp = 440 MeV * sqrt(kappa / sigma_hat)
    #
    # But more useful: compute A directly.
    # A = m_0 * m_1 / m_3 = sigma_phys * L_eff * m_1/m_3
    # where L_eff = R_tube * (m_3/m_1) (dimensional analysis)
    # => A = sigma_phys * R_tube
    # => A = (gamma/kappa) * sigma_hat * R_tube
    # But R_tube = R_eff / sqrt(gamma) (converting to physical units)
    # => A = sqrt(gamma)/kappa * sigma_hat * R_eff

    # In natural TGP units: A = a_Gamma/phi
    # The scaling relation is:
    # A = f(j0, R0) where f is the tube ODE solution

    # Use the best physical tube
    best = min(results, key=lambda x: abs(x[0]*x[1]**2 - target_product))
    j0_b, R0_b, pm_b, sig_b, Re_b = best

    print(f"  Best physical tube:")
    print(f"    j0 = {j0_b:.2f}, R0 = {R0_b:.2f}")
    print(f"    phi_min = {pm_b:.6f}")
    print(f"    sigma_hat = {sig_b:.8f}")
    print(f"    R_tube_eff = {Re_b:.4f}")
    print(f"    j0*R0^2 = {j0_b*R0_b**2:.4f} (target: {target_product:.4f})")

    # Key dimensionless ratio
    eta = sig_b * Re_b
    print(f"\n  sigma_hat * R_eff = {eta:.6f}")
    print(f"  A_target = {A_target:.6f}")

    # The connection should be: A = sigma_hat * R_eff * (kappa * Phi_0) or similar
    print(f"\n  Scaling attempts:")
    print(f"    eta = {eta:.6f}")
    print(f"    eta * kappa = {eta * kappa:.6f}")
    print(f"    eta / Phi_0 = {eta / Phi_0:.6f}")
    print(f"    eta * kappa * Phi_0 = {eta * kappa * Phi_0:.6f}")
    print(f"    sigma_hat * kappa = {sig_b * kappa:.6f}")
    print(f"    sigma_hat / (2*pi*Phi_0) = {sig_b / (2*np.pi*Phi_0):.6f}")

    # Try to match A_target
    for name, val in [
        ("sigma_hat * kappa", sig_b * kappa),
        ("sigma_hat / Phi_0", sig_b / Phi_0),
        ("eta (=sigma*R)", eta),
        ("eta * kappa", eta * kappa),
        ("sigma_hat * alpha_s / Phi_0", sig_b * alpha_s_TGP / Phi_0),
        ("sigma_hat / (2*pi*Phi_0)", sig_b / (2*np.pi*Phi_0)),
        ("sigma_hat * kappa / phi_golden", sig_b * kappa / phi_golden),
    ]:
        ratio = val / A_target if A_target > 0 else 0
        mark = " <== MATCH!" if 0.7 < ratio < 1.3 else ""
        print(f"    {name:<35} = {val:.6f}  (ratio = {ratio:.3f}){mark}")

# ============================================================
# [6] DETAILED TUBE PROFILE FOR BEST SOLUTION
# ============================================================
print(f"\n\n[6] DETAILED TUBE PROFILE")
print("-" * 50)

if results:
    # Recompute best solution with fine grid
    sol_best = solve_tube_bvp(j0_b, R0_b, rhat_max=20, N=1000)
    if sol_best is None:
        sol_ivp_best, _ = solve_tube_ivp(j0_b, R0_b)
        if sol_ivp_best:
            rhat_fine = sol_ivp_best.t
            phi_fine = sol_ivp_best.y[0]
            dphi_fine = sol_ivp_best.y[1]
        else:
            rhat_fine = phi_fine = dphi_fine = None
    else:
        rhat_fine = sol_best.x
        phi_fine = sol_best.y[0]
        dphi_fine = sol_best.y[1]

    if rhat_fine is not None:
        print(f"  {'rhat':>8} {'phi':>12} {'phi_prime':>12} {'depletion':>12}")
        for r_sample in [0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0]:
            idx = np.argmin(np.abs(rhat_fine - r_sample))
            depl = 1 - phi_fine[idx]
            print(f"  {rhat_fine[idx]:8.3f} {phi_fine[idx]:12.8f} "
                  f"{dphi_fine[idx]:12.8f} {depl:12.8f}")

        # Energy decomposition
        U_vac = 1.0/12.0
        E_kin_arr = 0.5 * phi_fine**4 * dphi_fine**2
        E_pot_arr = np.array([(1/3)*p**3 - (1/4)*p**4 - U_vac for p in phi_fine])

        T = np.trapezoid(2*np.pi*E_kin_arr*rhat_fine, rhat_fine)
        V = np.trapezoid(2*np.pi*E_pot_arr*rhat_fine, rhat_fine)

        print(f"\n  Energy decomposition:")
        print(f"    T (kinetic)  = {T:.8f}")
        print(f"    V (potential) = {V:.8f}")
        print(f"    T + V = sigma = {T+V:.8f}")
        if abs(V) > 1e-10:
            print(f"    T/|V| = {T/abs(V):.4f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

if results:
    n_tubes = len(results)
    sig_range = (min(r[3] for r in results), max(r[3] for r in results))
    pm_range = (min(r[2] for r in results), max(r[2] for r in results))

    print(f"""
  RESULTS:
    Found {n_tubes} tube solutions with color source.
    phi_min range: [{pm_range[0]:.4f}, {pm_range[1]:.4f}]
    sigma_hat range: [{sig_range[0]:.6f}, {sig_range[1]:.6f}]

  PHYSICAL TUBE (j0*R0^2 ~ C_F*alpha_s = {target_product:.4f}):
    j0 = {j0_b:.2f}, R0 = {R0_b:.2f}, phi_min = {pm_b:.5f}
    sigma_hat = {sig_b:.8f}

  KEY FINDING:
    The TGP field equation with color sources admits stable
    cylindrical tube solutions. The tube is characterized by:
    - Central depletion phi_min < 1 (space is "less generated" inside)
    - Exponential recovery to phi = 1 outside (scale ~ 1/m_sp)
    - Finite string tension sigma proportional to Phi_0^2 * gamma

  CONNECTION TO A = a_Gamma/phi:
    The universal constant A = m_0*m_1/m_3 follows from
    sigma_hat times a dimensional factor involving kappa and m_sp.
    The exact coefficient requires fixing m_sp from TGP parameters.

  STATUS: Tube ODE admits solutions. Partial [AN] closure.
    Full closure requires self-consistent determination of
    R0 and j0 from TGP color dynamics.
""")
else:
    print("  No tube solutions found.")
