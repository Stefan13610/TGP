"""
soliton_interaction.py -- Two-soliton interaction from field overlap
=====================================================================

PHYSICAL MOTIVATION:
--------------------
The n-body force law in TGP (pairwise V_2 + three-body V_3) is POSTULATED
in Path B as arising from Yukawa sources C_i at each particle.

Path C (P5, yukawa_from_defect.py) DERIVES C_eff from the classical defect
profile via projection onto the Yukawa Green's function.

This module COMPLETES the loop by computing the ACTUAL overlap between
two defect profiles and comparing with the Yukawa prediction.

LEADING-ORDER INTERACTION:
--------------------------
Two defects in linear superposition: g(x) = 1 - delta_1(|x|) - delta_2(|x-d|)

Expanding the TGP energy to second order in deltas:
  E = E_1 + E_2 + V_cross

where the cross-term (interaction) at leading order is:

  V_cross(d) = int [grad(delta_1) . grad(delta_2) + V''(1)*delta_1*delta_2] d^3x
             = int [-delta_1 * nabla^2(delta_2) + V''(1)*delta_1*delta_2] d^3x
             = -int delta_1 * [nabla^2 + |V''(1)|] * delta_2 d^3x

For standard TGP:  V''(1) = 2*beta - 3*gamma = -1 (for beta=gamma=1)
  The linearized equation is (nabla^2 + 1)*delta = 0 at large r.
  So V_cross ~ 0 at leading order (defects decouple far-field)!
  Nonzero V_cross comes from nonlinear (core overlap) terms.

For EFT (Yukawa propagator): V''_eff = +m_sp^2
  V_cross_EFT(d) = -int delta_1 * [nabla^2 - m_sp^2] * delta_2 d^3x
  With delta_2 ~ C*exp(-m*r)/r: (nabla^2 - m^2)*delta_2 = -4*pi*C*delta^3(x)
  => V_cross_EFT = 4*pi*C_1*C_2*delta_1(d) ~ exp(-m*d)/d  (YUKAWA!)

This shows WHY classical field theory gives oscillatory (decoupled)
interaction while the EFT mass gap produces Yukawa.

COMPUTED QUANTITIES:
-------------------
1. Overlap integral: S(d) = int delta_1(|x|) * delta_2(|x-d*z|) d^3x
   = 2*pi * int_0^inf int_{-inf}^{inf} delta_1(sqrt(rho^2+z^2))
     * delta_2(sqrt(rho^2+(z-d)^2)) * rho drho dz

2. Gradient cross: G(d) = int grad(delta_1) . grad(delta_2) d^3x

3. Classical interaction: V_cl(d) = G(d) - S(d)   [since V''(1) = -1]

4. Yukawa prediction: V_Y(d) = -4*pi * C_eff^2 * exp(-m_sp*d) / d

5. Far-field value: delta(d) = 1 - g(d), directly from defect profile
"""

import numpy as np
from .yukawa_from_defect import (
    solve_defect_standard, V_tgp, dV_tgp,
    compute_C_eff_projection,
)
from .tgp_field import screening_mass


def generate_defect_profile(g0, beta=1.0, gamma=1.0, r_max=60.0,
                            n_eval=5000):
    """
    Generate a single defect profile (standard TGP, no stabilization).

    Returns dict with r, g, gp arrays and metadata (C_eff_proj, m_sp, E_defect).
    """
    return solve_defect_standard(g0, beta, gamma, r_max=r_max,
                                 kinetic="full", n_eval=n_eval)


def _interp_delta(r_grid, g_grid, r_query):
    """
    Interpolate delta(r) = 1 - g(r) at arbitrary radii.

    For r > r_max: delta = 0 (vacuum).
    For r < r_min: delta = 1 - g(r_min) (flat core).
    """
    r_min = r_grid[0]
    r_max = r_grid[-1]
    delta_grid = 1.0 - g_grid

    result = np.zeros_like(r_query)

    mask_in = (r_query >= r_min) & (r_query <= r_max)
    if np.any(mask_in):
        result[mask_in] = np.interp(r_query[mask_in], r_grid, delta_grid)

    mask_lo = r_query < r_min
    if np.any(mask_lo):
        result[mask_lo] = delta_grid[0]

    return result


def _interp_gp(r_grid, gp_grid, r_query):
    """Interpolate g'(r) at arbitrary radii. Zero outside grid."""
    r_min = r_grid[0]
    r_max = r_grid[-1]

    result = np.zeros_like(r_query)
    mask_in = (r_query >= r_min) & (r_query <= r_max)
    if np.any(mask_in):
        result[mask_in] = np.interp(r_query[mask_in], r_grid, gp_grid)
    return result


def overlap_integral(defect, d, n_rho=300, n_z=600, rho_max=None, z_max=None):
    """
    Compute the overlap integral between two identical defects at separation d.

    S(d) = 2*pi * int_0^rho_max int_{-z_max}^{z_max}
           delta_1(sqrt(rho^2+z^2)) * delta_2(sqrt(rho^2+(z-d)^2)) * rho drho dz

    where defect 1 is at origin, defect 2 is at z = d.

    Returns S(d) (scalar).
    """
    r_grid = defect['r']
    g_grid = defect['g']
    r_max_def = r_grid[-1]

    if rho_max is None:
        rho_max = min(r_max_def * 0.7, 25.0)
    if z_max is None:
        z_max = d + r_max_def * 0.7

    rho_arr = np.linspace(0, rho_max, n_rho)
    z_arr = np.linspace(-r_max_def * 0.5, d + r_max_def * 0.5, n_z)
    RHO, Z = np.meshgrid(rho_arr, z_arr, indexing='ij')

    R1 = np.sqrt(RHO**2 + Z**2)
    R2 = np.sqrt(RHO**2 + (Z - d)**2)

    d1 = _interp_delta(r_grid, g_grid, R1.ravel()).reshape(R1.shape)
    d2 = _interp_delta(r_grid, g_grid, R2.ravel()).reshape(R2.shape)

    integrand = d1 * d2 * RHO  # * rho for cylindrical

    drho = rho_arr[1] - rho_arr[0]
    dz = z_arr[1] - z_arr[0]

    S = 2.0 * np.pi * np.sum(integrand) * drho * dz
    return float(S)


def gradient_cross_integral(defect, d, n_rho=300, n_z=600,
                            rho_max=None, z_max=None):
    """
    Compute gradient cross-term: G(d) = int grad(delta_1) . grad(delta_2) d^3x.

    Since grad(delta) = -g'(r) * r_hat, we have:
    grad(delta_1) . grad(delta_2) = g'_1 * g'_2 * (r_hat_1 . r_hat_2)

    In cylindrical:
      r_hat_1 . r_hat_2 = (rho^2 + z*(z-d)) / (r_1 * r_2)

    Returns G(d) (scalar).
    """
    r_grid = defect['r']
    g_grid = defect['g']
    gp_grid = defect['gp']
    r_max_def = r_grid[-1]

    if rho_max is None:
        rho_max = min(r_max_def * 0.7, 25.0)
    if z_max is None:
        z_max = d + r_max_def * 0.7

    rho_arr = np.linspace(0, rho_max, n_rho)
    z_arr = np.linspace(-r_max_def * 0.5, d + r_max_def * 0.5, n_z)
    RHO, Z = np.meshgrid(rho_arr, z_arr, indexing='ij')

    R1 = np.sqrt(RHO**2 + Z**2)
    R2 = np.sqrt(RHO**2 + (Z - d)**2)

    # grad(delta_i) = -g'(r_i) * r_hat_i
    # So grad(delta_1) . grad(delta_2) = g'_1 * g'_2 * cos(angle)
    gp1 = _interp_gp(r_grid, gp_grid, R1.ravel()).reshape(R1.shape)
    gp2 = _interp_gp(r_grid, gp_grid, R2.ravel()).reshape(R2.shape)

    R1_safe = np.maximum(R1, 1e-10)
    R2_safe = np.maximum(R2, 1e-10)

    # cos(angle) = r_hat_1 . r_hat_2
    cos_angle = (RHO**2 + Z * (Z - d)) / (R1_safe * R2_safe)

    integrand = gp1 * gp2 * cos_angle * RHO

    drho = rho_arr[1] - rho_arr[0]
    dz = z_arr[1] - z_arr[0]

    G = 2.0 * np.pi * np.sum(integrand) * drho * dz
    return float(G)


def classical_interaction(defect, d, beta=1.0, gamma=1.0,
                          n_rho=300, n_z=600):
    """
    Compute leading-order classical interaction energy:

    V_cl(d) = G(d) + V''(1)*S(d) = G(d) - (3*gamma - 2*beta)*S(d)

    For beta=gamma=1: V_cl = G(d) - S(d)

    Note: V''(1) = 2*beta - 3*gamma (negative for standard TGP).

    Returns dict with V_cl, G, S, and diagnostic info.
    """
    V_pp = 2.0 * beta - 3.0 * gamma  # V''(1) at vacuum

    S = overlap_integral(defect, d, n_rho, n_z)
    G = gradient_cross_integral(defect, d, n_rho, n_z)

    V_cl = G + V_pp * S  # = G - |V''(1)| * S for standard TGP

    return {
        'd': d,
        'V_classical': V_cl,
        'overlap_S': S,
        'gradient_G': G,
        'V_pp': V_pp,
    }


def far_field_value(defect, d):
    """
    Return delta(d) = 1 - g(d) directly from defect profile.

    This is the simplest way to see interaction:
    V ~ C_1 * delta_2(d) for well-separated solitons.
    """
    r_grid = defect['r']
    g_grid = defect['g']
    delta_grid = 1.0 - g_grid

    if d > r_grid[-1]:
        return 0.0
    return float(np.interp(d, r_grid, delta_grid))


def interaction_scan(defect, d_values, beta=1.0, gamma=1.0,
                     n_rho=300, n_z=600):
    """
    Scan classical interaction V_cl(d) and compare with Yukawa prediction.

    Returns dict with arrays for each separation.
    """
    d_values = np.asarray(d_values, dtype=float)
    m_sp = screening_mass(beta, gamma)
    C_eff = defect['C_eff_proj']

    results = {
        'd': d_values,
        'V_classical': np.zeros(len(d_values)),
        'overlap_S': np.zeros(len(d_values)),
        'gradient_G': np.zeros(len(d_values)),
        'V_yukawa': np.zeros(len(d_values)),
        'delta_at_d': np.zeros(len(d_values)),
        'C_eff': C_eff,
        'm_sp': m_sp,
    }

    for i, d in enumerate(d_values):
        res = classical_interaction(defect, d, beta, gamma, n_rho, n_z)
        results['V_classical'][i] = res['V_classical']
        results['overlap_S'][i] = res['overlap_S']
        results['gradient_G'][i] = res['gradient_G']

        # Yukawa prediction
        results['V_yukawa'][i] = -4.0 * np.pi * C_eff**2 * np.exp(-m_sp * d) / d

        # Direct far-field value
        results['delta_at_d'][i] = far_field_value(defect, d)

    return results


def yukawa_interaction(d, C_eff, m_sp):
    """Yukawa potential: V_Y(d) = -4*pi * C^2 * exp(-m*d) / d."""
    return -4.0 * np.pi * C_eff**2 * np.exp(-m_sp * d) / d


def yukawa_force(d, C_eff, m_sp):
    """Yukawa force: F_Y(d) = -dV_Y/dd = -4*pi*C^2*exp(-m*d)*(1/d^2 + m/d)."""
    return -4.0 * np.pi * C_eff**2 * np.exp(-m_sp * d) * (1.0 / d**2 + m_sp / d)
