#!/usr/bin/env python3
"""
tgp_strong_field_solver.py — Full strong-field TGP soliton profile
===================================================================

Solves the static, spherically symmetric TGP vacuum field equation:

    φ'' + (2/r)φ' + 2(φ')²/φ + β φ² - γ φ³ = 0

with boundary conditions:
    φ(r → ∞) = 1   (vacuum)
    φ(r → 0) → ∞   (black hole / frozen interior)

The metric is:
    ds² = -(c₀²/φ) dt² + φ δᵢⱼ dxⁱ dxʲ

Key outputs:
    - Full φ(r) profile from center to infinity
    - Photon sphere radius: minimum of impact parameter b(r)
    - Shadow size comparison with GR (Schwarzschild)

Uses the substitution u = φ³ which linearizes the kinetic operator:
    D_kin[φ] = ∇²φ + 2(∇φ)²/φ = (1/3φ²) ∇²(φ³) = (1/3φ²) ∇²u

So the ODE becomes:
    u'' + (2/r) u' + 3φ² (β φ² - γ φ³) = 0
    u'' + (2/r) u' + 3β u^{4/3} - 3γ u^{5/3} = 0   (where φ = u^{1/3})

For the BH solution we shoot from r_min with large φ₀ inward boundary.

Author: Mateusz Serafin (with Claude)
Date: April 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import sys

# ============================================================
# 1. Field equation in original variable φ
# ============================================================

def field_rhs(r, y, beta, gamma_coeff):
    """
    RHS of the ODE system:
        y[0] = φ(r)
        y[1] = φ'(r)

    φ'' = -(2/r)φ' - 2(φ')²/φ - β φ² + γ φ³
    """
    phi, dphi = y
    if phi < 1e-30:
        phi = 1e-30  # safety floor

    ddphi = -(2.0/r) * dphi - 2.0 * dphi**2 / phi \
            - beta * phi**2 + gamma_coeff * phi**3

    return [dphi, ddphi]


def field_rhs_u(r, y, beta, gamma_coeff):
    """
    RHS using u = φ³ substitution (linearized kinetic term):
        y[0] = u(r)
        y[1] = u'(r)

    u'' = -(2/r)u' - 3β u^{4/3} + 3γ u^{5/3}
    """
    u, du = y
    if u < 1e-30:
        u = 1e-30

    ddu = -(2.0/r) * du - 3.0 * beta * u**(4.0/3.0) \
          + 3.0 * gamma_coeff * u**(5.0/3.0)

    return [du, ddu]


# ============================================================
# 2. Shooting from large r inward (BH profile)
# ============================================================

def shoot_bh_profile(M_param, beta=1.0, gamma_coeff=1.0,
                     r_max=500.0, r_min=0.01, n_points=10000):
    """
    Shoot from r_max inward with Newtonian asymptotic:
        φ(r) ≈ 1 + 2U/c₀² = 1 - r_s/r   (weak field, r_s = 2GM/c₀²)

    Parameters
    ----------
    M_param : float
        Dimensionless mass parameter r_s (Schwarzschild radius in code units)
    beta, gamma_coeff : float
        Self-interaction parameters (β = γ for vacuum condition)
    r_max : float
        Outer boundary (should be >> r_s)
    r_min : float
        Inner boundary (smallest radius to integrate to)
    n_points : int
        Number of output points

    Returns
    -------
    r_arr : array, φ_arr : array, dφ_arr : array
    """
    r_s = M_param

    # Asymptotic boundary condition at r_max (weak field)
    # φ = exp(2U/c²) ≈ exp(-r_s/r) for exponential metric
    phi_inf = np.exp(-r_s / r_max)
    dphi_inf = (r_s / r_max**2) * np.exp(-r_s / r_max)

    y0 = [phi_inf, dphi_inf]

    r_eval = np.linspace(r_max, r_min, n_points)

    sol = solve_ivp(
        lambda r, y: field_rhs(r, y, beta, gamma_coeff),
        [r_max, r_min],
        y0,
        method='RK45',
        t_eval=r_eval,
        rtol=1e-12,
        atol=1e-14,
        max_step=0.1
    )

    if not sol.success:
        print(f"  Integration warning: {sol.message}")

    return sol.t, sol.y[0], sol.y[1]


def shoot_bh_profile_exact(M_param, beta=1.0, gamma_coeff=1.0,
                           r_max=500.0, r_min=0.01, n_points=10000):
    """
    Same as above but using the EXACT exponential boundary condition
    at r_max, and the full nonlinear ODE.

    The exponential metric gives φ = exp(-r_s/r) exactly in the
    linearized limit. For the full nonlinear solution, we use this
    as the starting point and let the ODE evolve.

    Key insight: the exponential metric IS the exact weak-field solution.
    The question is whether the full nonlinear ODE deviates from it
    in the strong field (r ~ r_s).
    """
    r_s = M_param

    # --- Exact exponential ansatz as outer BC ---
    # φ_exp(r) = exp(-r_s/r)
    # φ'_exp(r) = (r_s/r²) exp(-r_s/r)
    # φ''_exp(r) = (-2r_s/r³ + r_s²/r⁴) exp(-r_s/r)
    #
    # Check: does φ_exp satisfy the full ODE?
    # Substituting into φ'' + (2/r)φ' + 2(φ')²/φ + βφ² - γφ³ = 0:
    #
    # Term 1: φ'' = (-2r_s/r³ + r_s²/r⁴) exp(-r_s/r)
    # Term 2: (2/r)φ' = (2r_s/r³) exp(-r_s/r)
    # Term 3: 2(φ')²/φ = 2(r_s²/r⁴) exp(-r_s/r)
    # Sum of kinetic: (r_s²/r⁴ + 2r_s²/r⁴) exp(-r_s/r) = 3r_s²/r⁴ exp(-r_s/r)
    # Term 4: βφ² - γφ³ = β exp(-2r_s/r) - γ exp(-3r_s/r)
    #
    # For the ODE to be satisfied: 3r_s²/r⁴ exp(-r_s/r) + β exp(-2r_s/r) - γ exp(-3r_s/r) = 0
    # This is NOT zero in general!
    #
    # So the exponential metric is NOT the exact solution of the full nonlinear ODE.
    # The deviation is strongest where r ~ r_s.

    phi_start = np.exp(-r_s / r_max)
    dphi_start = (r_s / r_max**2) * np.exp(-r_s / r_max)

    y0 = [phi_start, dphi_start]

    r_eval = np.linspace(r_max, r_min, n_points)

    # Event: stop if φ goes negative or blows up
    def phi_zero(r, y):
        return y[0] - 1e-10
    phi_zero.terminal = True
    phi_zero.direction = -1

    def phi_huge(r, y):
        return y[0] - 1e6
    phi_huge.terminal = True
    phi_huge.direction = 1

    sol = solve_ivp(
        lambda r, y: field_rhs(r, y, beta, gamma_coeff),
        [r_max, r_min],
        y0,
        method='Radau',  # stiff-capable
        t_eval=r_eval,
        rtol=1e-12,
        atol=1e-15,
        max_step=0.05,
        events=[phi_zero, phi_huge]
    )

    return sol.t, sol.y[0], sol.y[1], sol


# ============================================================
# 3. Photon sphere and shadow from φ(r)
# ============================================================

def compute_photon_sphere(r_arr, phi_arr):
    """
    In TGP isotropic coordinates, the metric is:
        ds² = -(c₀²/φ) dt² + φ δᵢⱼ dxⁱdxʲ

    For null geodesics, the impact parameter is:
        b(r) = r √(g_rr / |g_tt|) = r √(φ · φ) = r · φ(r)

    Wait — let me be more careful.

    g_tt = -c₀²/φ,  g_rr = φ  (isotropic)

    For a photon orbit in the equatorial plane:
        b = r √(g_φφ / |g_tt|)

    In isotropic coords: g_φφ = φ r² (since g_ij = φ δ_ij in spherical → g_φφ = φ r²)

    So: b(r) = r √(φ r² / (c₀²/φ)) = ... wait, let me redo.

    The effective potential approach for null geodesics:

    ds² = -A(r) dt² + B(r) dr² + C(r) dΩ²

    In our case (isotropic coords):
        A(r) = c₀²/φ
        B(r) = φ
        C(r) = φ r²

    Photon sphere: d/dr[C(r)/A(r)] = 0

    C/A = φ r² / (c₀²/φ) = φ² r² / c₀²

    d/dr[φ² r²] = 2φ φ' r² + 2φ² r = 2φ r (φ' r + φ) = 0

    Since φ > 0 and r > 0:
        φ'(r_ph) r_ph + φ(r_ph) = 0
        → φ'(r_ph) = -φ(r_ph) / r_ph

    Impact parameter at photon sphere:
        b_crit = √(C/A)|_{r_ph} = φ(r_ph) r_ph / c₀
    """
    # Find where φ'·r + φ = 0, i.e., d/dr[r·φ] = 0
    # Numerically: compute h(r) = r · φ(r), find its minimum

    h = r_arr * phi_arr

    # Find local minimum of h (photon sphere)
    idx_min = np.argmin(h[1:]) + 1  # skip r=0

    if idx_min <= 1 or idx_min >= len(h) - 1:
        return None, None, None

    # Refine with parabolic interpolation
    r_ph = r_arr[idx_min]
    b_crit = h[idx_min]  # = r_ph * φ(r_ph), up to c₀ factor
    phi_ph = phi_arr[idx_min]

    return r_ph, b_crit, phi_ph


def shadow_ratio_vs_schwarzschild(r_ph, b_crit, r_s):
    """
    Compare TGP shadow with Schwarzschild.

    Schwarzschild in Schwarzschild coords:
        R_ph = 3GM/c² = (3/2) r_s
        b_crit_GR = 3√3 GM/c² = (3√3/2) r_s

    Shadow angular size ratio:
        θ_TGP / θ_GR = b_crit_TGP / b_crit_GR
    """
    b_crit_gr = (3.0 * np.sqrt(3.0) / 2.0) * r_s
    ratio = b_crit / b_crit_gr
    return ratio


# ============================================================
# 4. Analytical checks
# ============================================================

def check_exponential_residual(r_arr, r_s, beta, gamma_coeff):
    """
    Check: how well does φ_exp(r) = exp(-r_s/r) satisfy the full ODE?

    Residual = φ'' + (2/r)φ' + 2(φ')²/φ + β φ² - γ φ³

    For exact solution, residual = 0.
    """
    phi = np.exp(-r_s / r_arr)
    dphi = (r_s / r_arr**2) * phi
    ddphi = (-2.0 * r_s / r_arr**3 + r_s**2 / r_arr**4) * phi

    kinetic = ddphi + (2.0 / r_arr) * dphi + 2.0 * dphi**2 / phi
    potential = beta * phi**2 - gamma_coeff * phi**3
    residual = kinetic + potential

    return residual, phi


def exponential_photon_sphere(r_s):
    """
    Analytical photon sphere for the exponential metric φ = exp(-r_s/r).

    Condition: φ'·r + φ = 0
    → (r_s/r²)·exp(-r_s/r)·r + exp(-r_s/r) = 0
    → (r_s/r + 1) exp(-r_s/r) = 0
    → r_s/r + 1 = 0  (impossible for r > 0!)

    Wait — this means d/dr[r·φ] = d/dr[r·exp(-r_s/r)]
    = exp(-r_s/r) + r·(r_s/r²)·exp(-r_s/r)
    = exp(-r_s/r)·(1 + r_s/r)

    This is ALWAYS POSITIVE for r > 0!

    So r·φ = r·exp(-r_s/r) is monotonically increasing — NO MINIMUM.

    That means the exponential metric in isotropic coordinates
    has NO photon sphere at all!

    This is consistent with the fact that the Yilmaz metric
    (exponential isotropic) is known to lack a photon sphere
    in isotropic coordinates. The "photon sphere" in the BH shadow
    paper was found by transforming to curvature (Schwarzschild-like)
    coordinates.

    Let me redo the analysis in curvature coordinates.
    """
    # In curvature coordinates R = r·√φ = r·exp(-r_s/(2r)):
    # The metric becomes (approximately):
    #   ds² = -exp(-r_s/r(R)) dt² + ... dR² + R² dΩ²
    # The photon sphere condition is d/dR[R²/|g_tt(R)|] = 0
    #
    # In isotropic coords the photon sphere condition is different:
    #   d/dr [ C(r)/A(r) ] = 0 where C = φ r², A = c₀²/φ
    #   d/dr [ φ²r² ] = 0
    #   2φφ'r² + 2φ²r = 0
    #   φ'r + φ = 0
    #   φ'/φ = -1/r
    #   For φ = exp(-r_s/r): φ'/φ = r_s/r²
    #   So: r_s/r² = -1/r → r = -r_s (no solution for r > 0)
    #
    # CORRECT: in isotropic coordinates, the exponential metric
    # has NO photon sphere because g_tt never reaches zero.
    # The shadow analysis must be done in curvature coordinates.

    print("Exponential metric in isotropic coords: NO photon sphere")
    print("(g_tt = -1/φ = -exp(r_s/r) never reaches zero)")
    print("Must transform to curvature coordinates for shadow analysis.")
    print()

    # Curvature coordinate transformation: R = r·√φ = r·exp(-r_s/(2r))
    # Photon sphere in curvature coords: R_ph such that
    #   d/dR [R² · g_tt_curv(R)] = 0
    #
    # From the BH shadow paper: R_ph = r_s/2, b_crit = (r_s/2)·e
    # Shadow ratio = e/(3√3) ≈ 0.523

    R_ph_curv = r_s / 2.0
    b_crit_curv = (r_s / 2.0) * np.e
    b_crit_gr = (3.0 * np.sqrt(3.0) / 2.0) * r_s
    ratio = b_crit_curv / b_crit_gr

    print(f"Curvature-coordinate analytical results (exponential metric):")
    print(f"  R_ph = r_s/2 = {R_ph_curv:.4f}")
    print(f"  b_crit = (r_s/2)·e = {b_crit_curv:.4f}")
    print(f"  b_crit_GR = (3√3/2)·r_s = {b_crit_gr:.4f}")
    print(f"  Shadow ratio = {ratio:.4f} ({(1-ratio)*100:.1f}% smaller)")

    return R_ph_curv, b_crit_curv, ratio


# ============================================================
# 5. Transform to curvature coordinates and find photon sphere
# ============================================================

def isotropic_to_curvature(r_arr, phi_arr, dphi_arr):
    """
    Transform from isotropic to curvature (areal) coordinates.

    Isotropic metric: ds² = -(c₀²/φ) dt² + φ (dr² + r²dΩ²)

    Curvature radius: R = r·√φ  (so that g_ΩΩ = R²)

    The metric in curvature coordinates:
        ds² = -A(R) dt² + B(R) dR² + R² dΩ²

    where A(R) = c₀²/φ(r(R)) and B(R) is determined by the
    coordinate transformation.

    Photon sphere condition: d/dR [R²/A(R)] = 0
    → d/dR [R² φ(r(R))] = 0
    """
    # R = r √φ
    R_arr = r_arr * np.sqrt(phi_arr)

    # A(R) = c₀²/φ = 1/φ  (setting c₀=1)
    A_arr = 1.0 / phi_arr

    # Photon sphere: d/dR [R²/A] = d/dR [R² φ] = 0
    # In terms of isotropic r: R² φ = r² φ · φ = r² φ²
    # So we need: d/dr [r² φ²] · dr/dR = 0
    # Since dr/dR ≠ 0 (monotonic transform), we need:
    # d/dr [r² φ²] = 0
    # → 2r φ² + 2r² φ φ' = 0
    # → φ + r φ' = 0  (same condition!)
    #
    # Hmm, same condition. But we can also work directly with R.

    # h(r) = R² · φ = r² φ²
    h = r_arr**2 * phi_arr**2

    # For GR comparison, h should have a minimum at the photon sphere.
    # Let's compute numerically.

    # Actually the photon sphere condition in curvature coords is
    # d/dR [R²/A] = d/dR [R² φ] = 0
    # R² φ = r² φ · φ = r² φ²
    # So minimum of r² φ² (in r-space, noting R is monotonic in r)

    # Check: for exponential φ = exp(-r_s/r):
    # r² φ² = r² exp(-2r_s/r)
    # d/dr = 2r exp(-2r_s/r) + r²·(2r_s/r²)·exp(-2r_s/r)
    #       = exp(-2r_s/r)·(2r + 2r_s)
    # This is always > 0! No minimum again!
    #
    # But wait — in the BH shadow paper we found a photon sphere.
    # The issue is: in Schwarzschild coords the photon sphere exists
    # because g_tt has a different radial dependence.
    #
    # The problem is that in the FULL nonlinear solution,
    # φ(r) may deviate from exp(-r_s/r) and create a photon sphere
    # that the exponential approximation misses.

    return R_arr, A_arr, h


def find_photon_sphere_curvature(r_arr, phi_arr):
    """
    Find photon sphere using the curvature-coordinate condition.

    In the Schwarzschild-like form ds² = -f(R)dt² + f(R)⁻¹dR² + R²dΩ²,
    the photon sphere is at f'(R_ph)·R_ph = 2·f(R_ph).

    But TGP's curvature-coord metric is NOT Schwarzschild form.
    Starting from isotropic:
        ds² = -(1/φ) dt² + φ(dr² + r²dΩ²)

    Transform R = r√φ, then:
        dR = √φ(1 + rφ'/(2φ)) dr
        dr = dR / [√φ(1 + rφ'/(2φ))]

    ds² = -(1/φ) dt² + φ·dR²/[φ(1+rφ'/(2φ))²] + R²dΩ²
        = -(1/φ) dt² + dR²/(1+rφ'/(2φ))² + R²dΩ²

    So: A(R) = 1/φ,  B(R) = 1/(1+rφ'/(2φ))²

    Photon sphere in general static spherical metric with R²dΩ²:
        condition: A'(R)/A(R) = 2/R

    In isotropic r:
        A = 1/φ,  R = r√φ
        dA/dr = -φ'/φ²
        dR/dr = √φ + rφ'/(2√φ) = √φ(1 + rφ'/(2φ))

        A'(R) = (dA/dr)/(dR/dr) = [-φ'/φ²] / [√φ(1 + rφ'/(2φ))]

        A'/A = [-φ'/φ²] / [√φ(1 + rφ'/(2φ))] · φ
             = -φ'/(φ√φ(1 + rφ'/(2φ)))

        2/R = 2/(r√φ)

    Condition: -φ'/(φ√φ(1 + rφ'/(2φ))) = 2/(r√φ)

    Simplify: -φ'/(φ(1 + rφ'/(2φ))) = 2/r

    → -rφ' = 2φ(1 + rφ'/(2φ))
    → -rφ' = 2φ + rφ'
    → -2rφ' = 2φ
    → φ'(r_ph) = -φ(r_ph)/r_ph

    SAME CONDITION as before! This is coordinate-invariant as it should be.

    The issue is that for φ_exp = exp(-r_s/r), we get:
        φ'/φ = r_s/r²  and we need  φ'/φ = -1/r
        → r_s/r = -1 (impossible)

    CONCLUSION: The pure exponential metric (linearized TGP solution)
    has NO photon sphere at all. The "48% smaller shadow" from the
    BH paper was derived by a coordinate-change trick that implicitly
    assumed a different metric form.

    THE REAL QUESTION: Does the FULL NONLINEAR solution φ(r) deviate
    from exp(-r_s/r) enough to create a photon sphere?

    Let's find out numerically.
    """
    # Condition: φ' + φ/r = 0, i.e., d/dr[r·φ] = 0
    rp = r_arr * phi_arr

    # Look for minimum
    drp = np.diff(rp)
    sign_changes = np.where(drp[:-1] * drp[1:] < 0)[0]

    results = []
    for idx in sign_changes:
        # Check if this is a minimum (drp goes from negative to positive)
        if drp[idx] < 0 and drp[idx+1] > 0:
            r_ph = r_arr[idx+1]
            phi_ph = phi_arr[idx+1]
            b_crit = r_ph * phi_ph  # = r·φ at photon sphere
            R_ph = r_ph * np.sqrt(phi_ph)  # curvature coordinate
            results.append({
                'r_ph_iso': r_ph,
                'R_ph_curv': R_ph,
                'phi_ph': phi_ph,
                'b_crit': b_crit,
            })

    return results


# ============================================================
# 6. Main analysis
# ============================================================

def analyze_strong_field(r_s=1.0, beta=1.0, gamma_coeff=1.0,
                         r_max=200.0, r_min=0.005, n_points=50000):
    """
    Full strong-field analysis:
    1. Check if exponential is exact (compute residual)
    2. Solve full nonlinear ODE
    3. Find photon sphere (if any)
    4. Compare with GR
    """
    print("=" * 70)
    print("TGP Strong-Field Analysis")
    print(f"  r_s = {r_s},  beta = {beta},  gamma = {gamma_coeff}")
    print(f"  Integration: r ∈ [{r_min}, {r_max}],  N = {n_points}")
    print("=" * 70)

    # --- Step 1: Analytical exponential check ---
    print("\n--- Step 1: Exponential metric residual ---")
    r_test = np.linspace(r_min, r_max, 10000)
    residual, phi_exp = check_exponential_residual(r_test, r_s, beta, gamma_coeff)

    # Relative residual (normalized by |potential term|)
    pot_scale = np.abs(beta * phi_exp**2 - gamma_coeff * phi_exp**3)
    pot_scale = np.maximum(pot_scale, 1e-30)
    rel_residual = np.abs(residual) / pot_scale

    print(f"  Max |residual| = {np.max(np.abs(residual)):.6e}")
    print(f"  Max |rel residual| = {np.max(rel_residual):.6e}")

    # Where is residual largest?
    idx_max = np.argmax(np.abs(residual))
    r_max_res = r_test[idx_max]
    print(f"  Largest residual at r = {r_max_res:.4f} (r/r_s = {r_max_res/r_s:.4f})")
    print(f"  φ_exp(r_max_res) = {phi_exp[idx_max]:.6f}")

    # --- Step 1b: Analytical photon sphere of exponential ---
    print("\n--- Step 1b: Exponential photon sphere analysis ---")
    exponential_photon_sphere(r_s)

    # --- Step 2: Solve full nonlinear ODE ---
    print("\n--- Step 2: Full nonlinear ODE integration ---")
    r_arr, phi_arr, dphi_arr, sol = shoot_bh_profile_exact(
        r_s, beta, gamma_coeff, r_max, r_min, n_points
    )

    print(f"  Integration success: {sol.success}")
    print(f"  Points: {len(r_arr)}")
    print(f"  r range: [{r_arr[-1]:.4f}, {r_arr[0]:.4f}]")
    print(f"  φ range: [{phi_arr.min():.6f}, {phi_arr.max():.6f}]")

    # Compare with exponential
    phi_exp_full = np.exp(-r_s / r_arr)
    delta = phi_arr - phi_exp_full
    rel_delta = delta / phi_exp_full

    print(f"\n  Deviation from exponential metric:")
    print(f"    max |φ - φ_exp| = {np.max(np.abs(delta)):.6e}")
    print(f"    max |Δφ/φ_exp|  = {np.max(np.abs(rel_delta)):.6e}")

    # Where is deviation largest?
    idx_max_dev = np.argmax(np.abs(rel_delta))
    print(f"    Largest deviation at r = {r_arr[idx_max_dev]:.4f} " +
          f"(r/r_s = {r_arr[idx_max_dev]/r_s:.4f})")
    print(f"    φ_ODE = {phi_arr[idx_max_dev]:.6f}, " +
          f"φ_exp = {phi_exp_full[idx_max_dev]:.6f}")

    # --- Step 3: Find photon sphere ---
    print("\n--- Step 3: Photon sphere search ---")

    # Check d/dr[r·φ] for the numerical solution
    rp = r_arr * phi_arr
    drp = np.gradient(rp, r_arr)

    # Is there a sign change in d/dr[r·φ]?
    n_sign_changes = np.sum(np.diff(np.sign(drp)) != 0)
    print(f"  d/dr[r·φ] sign changes: {n_sign_changes}")

    ph_results = find_photon_sphere_curvature(r_arr, phi_arr)

    if ph_results:
        print(f"  PHOTON SPHERE FOUND: {len(ph_results)} solution(s)")
        for i, res in enumerate(ph_results):
            b_gr = (3.0 * np.sqrt(3.0) / 2.0) * r_s
            ratio = res['b_crit'] / b_gr
            print(f"    [{i}] r_ph(iso) = {res['r_ph_iso']:.4f}, "
                  f"R_ph(curv) = {res['R_ph_curv']:.4f}, "
                  f"φ(r_ph) = {res['phi_ph']:.4f}")
            print(f"        b_crit = {res['b_crit']:.4f}, "
                  f"b_GR = {b_gr:.4f}, "
                  f"ratio = {ratio:.4f} ({(1-ratio)*100:.1f}% {'smaller' if ratio < 1 else 'larger'})")
    else:
        print("  NO PHOTON SPHERE FOUND in numerical solution")
        print("  (same as exponential analytical result)")
        print()
        print("  IMPLICATION: TGP exponential metric has NO shadow")
        print("  at all — photons are not captured. This is because")
        print("  g_tt = -1/φ never reaches zero (no horizon).")
        print()
        print("  This is QUALITATIVELY different from GR and means")
        print("  the BH shadow paper's 48% result was based on a")
        print("  coordinate transformation that assumed Schwarzschild-like")
        print("  behavior which the actual TGP metric does not have.")

    # --- Step 4: Check for circular orbits more carefully ---
    print("\n--- Step 4: Effective potential analysis ---")

    # For null geodesics with angular momentum L:
    # (dr/dλ)² = E² / B(r) · [1/A(r) - L²/(E²·C(r))]
    # = (1/φ) · [φ - L²/(E²·φ·r²)]  ... need to be careful
    #
    # Actually in isotropic coords:
    # ds² = -(1/φ)dt² + φ(dr² + r²dΩ²)
    #
    # Conserved: E = (1/φ)·dt/dλ, L = φ·r²·dθ/dλ
    # Null: 0 = -(1/φ)(dt/dλ)² + φ[(dr/dλ)² + r²(dθ/dλ)²]
    # → 0 = -E²φ + φ(dr/dλ)² + L²/φr²... wait
    #
    # Let me redo: g_μν (dx^μ/dλ)(dx^ν/dλ) = 0
    # -(1/φ)(dt/dλ)² + φ(dr/dλ)² + φ r²(dθ/dλ)² = 0
    #
    # E = -(1/φ)(dt/dλ) → dt/dλ = -Eφ
    # L = φ r²(dθ/dλ) → dθ/dλ = L/(φr²)
    #
    # -(1/φ)(Eφ)² + φ(dr/dλ)² + φr²·L²/(φr²)² = 0
    # -E²φ + φ(dr/dλ)² + L²/(φr²) = 0
    # (dr/dλ)² = E²  - L²/(φ²r²)
    #
    # Effective potential: V_eff = L²/(φ²r²)
    # Circular orbit: V_eff = E², dV_eff/dr = 0
    #
    # dV_eff/dr = L²·d/dr[1/(φ²r²)] = L²·[-2φ'/(φ³r²) - 2/(φ²r³)]
    #           = -2L²/(φ²r³)·[rφ'/φ + 1]
    #           = 0 when rφ'/φ + 1 = 0, i.e., φ' = -φ/r
    #
    # Same condition. Confirmed.

    # V_eff(r) = 1/(φ²r²) (dropping L²)
    V_eff = 1.0 / (phi_arr**2 * r_arr**2)

    # For a photon sphere, V_eff must have a LOCAL MAXIMUM
    dV = np.gradient(V_eff, r_arr)

    # Sign changes in dV (look for max: dV goes + to -)
    n_max = 0
    for i in range(len(dV)-1):
        if dV[i] > 0 and dV[i+1] < 0:
            n_max += 1
            print(f"  V_eff local MAX at r ≈ {r_arr[i]:.4f}, V_eff = {V_eff[i]:.6f}")

    if n_max == 0:
        print("  V_eff has NO local maximum → no unstable circular photon orbit")
        print("  V_eff is monotonically decreasing (r→0: V→∞, r→∞: V→0)")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Exponential residual (max): {np.max(np.abs(residual)):.4e}")
    print(f"  ODE deviation from exp (max rel): {np.max(np.abs(rel_delta)):.4e}")

    if ph_results:
        print(f"  Photon sphere: YES")
        ratio = ph_results[0]['b_crit'] / ((3*np.sqrt(3)/2)*r_s)
        print(f"  Shadow ratio vs GR: {ratio:.4f}")
    else:
        print(f"  Photon sphere: NO")
        print(f"  Shadow: TGP predicts NO classical BH shadow")
        print(f"          (fundamentally different from GR)")

    return {
        'r': r_arr, 'phi': phi_arr, 'dphi': dphi_arr,
        'phi_exp': phi_exp_full, 'delta': delta,
        'residual': residual, 'r_test': r_test,
        'photon_spheres': ph_results,
    }


# ============================================================
# 7. Parameter scan: does ANY β,γ produce a photon sphere?
# ============================================================

def scan_beta_gamma(r_s=1.0, r_max=200.0, r_min=0.01, n_points=20000):
    """
    Scan over β, γ parameter space to check if any combination
    produces a photon sphere.
    """
    print("\n" + "=" * 70)
    print("Parameter scan: β, γ → photon sphere?")
    print("=" * 70)

    betas = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    gammas = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    for beta in betas:
        for gamma_coeff in gammas:
            try:
                r_arr, phi_arr, dphi_arr, sol = shoot_bh_profile_exact(
                    r_s, beta, gamma_coeff, r_max, r_min, n_points
                )
                if not sol.success or len(r_arr) < 100:
                    continue

                ph = find_photon_sphere_curvature(r_arr, phi_arr)
                status = f"PHOTON SPHERE (n={len(ph)})" if ph else "no"

                # Also check deviation from exponential
                phi_exp = np.exp(-r_s / r_arr)
                max_dev = np.max(np.abs((phi_arr - phi_exp)/phi_exp))

                print(f"  β={beta:5.1f}, γ={gamma_coeff:5.1f}: "
                      f"max|Δφ/φ|={max_dev:.2e}, photon sphere: {status}")

                if ph:
                    for res in ph:
                        b_gr = (3*np.sqrt(3)/2)*r_s
                        print(f"    → r_ph={res['r_ph_iso']:.3f}, "
                              f"b/b_GR={res['b_crit']/b_gr:.3f}")
            except Exception as e:
                print(f"  β={beta:5.1f}, γ={gamma_coeff:5.1f}: FAILED ({e})")


# ============================================================
if __name__ == "__main__":
    quick = "--quick" in sys.argv

    if quick:
        results = analyze_strong_field(r_s=1.0, r_max=100.0, r_min=0.02,
                                       n_points=10000)
    else:
        results = analyze_strong_field(r_s=1.0, r_max=500.0, r_min=0.005,
                                       n_points=50000)
        scan_beta_gamma()

    print("\nDone.")
