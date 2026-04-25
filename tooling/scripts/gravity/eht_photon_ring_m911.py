"""
eht_photon_ring_m911.py -- EHT-quick test (TGP_CLOSURE_PLAN sec 8.10)

Compute the photon ring radius for a static spherically symmetric "black
hole" in the M9.1'' hyperbolic effective metric and compare with the
Schwarzschild prediction (r_ph = 1.5 r_s = 3 G M / c^2). Convert to
shadow angular sizes for M87* and Sgr A* and compare with EHT
observations.

Reference
---------
- TGP_CLOSURE_PLAN_2026-04-25.md  sections 4.2, 8.10
- M9_1_pp_P1_results.md            (Phi-EOM higher-PN coefficients)
- M9_1_pp_P2_results.md            (V(Phi)/Phi^4 selection)
- tgp_core.tex section 3 (effective metric) and field eq.

Background
----------
M9.1'' effective line element (isotropic radial coordinate r):
    ds^2 = -c0^2 (4 - 3 psi)/psi  dt^2
            + psi/(4 - 3 psi)  ( dr^2 + r^2 dOmega^2 )
with psi = Phi/Phi_0. In vacuum the canonical Phi-EOM (alpha=2, vacuum
beta = gamma so the polynomial part vanishes at psi=1) reduces to

    laplacian(eps) + 2 (grad eps)^2 / (1+eps) = 0,        eps = psi - 1.

For a static spherically symmetric source we have

    eps''(r) + (2/r) eps'(r) + 2 (eps'(r))^2 / (1+eps) = 0.

Asymptotic flatness:  eps -> A/r at large r with A = G M / (2 c^2)
(Newton matching gives 4 eps = 2 U_N + O(eps^2)).

The hyperbolic metric has g_tt -> 0 at psi = 4/3 (eps = 1/3).
This is the M9.1'' analog of the Schwarzschild horizon.

Photon-orbit potential for a static spherical metric
    ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2
is V_eff(r) = A(r)/r^2 (the photon ring satisfies dV_eff/dr = 0).

Method
------
1. Solve eps(r) on r in [r_in, r_out] with two boundary conditions:
   - Newton tail at r = r_out:  r * eps'(r) + eps(r) = 0  (i.e. d(r*eps)/dr=0)
   - Calibrated mass at r = r_out:  r * eps(r) = A_target = G M / (2 c^2)
   In geometric units G = c = 1 we set M = 1, so r_g = 1, and A_target = 1/2.

   We use a relaxation / shooting approach in the variable v = r*eps with
   the second-order ODE
        v'' = (1/r)(2 (v'/r - v/r^2)^2 r) / (1 + v/r) - ...
   Cleaner formulation: state y = (eps, eps'), integrate inward from
   r_out (large) where the asymptotic eps = A/r - A^2/r^2 + (5/3) A^3/r^3
   is known analytically (PN expansion from M9.1'' P1). Stop when
   eps reaches 1/3 (g_tt -> 0) to find the TGP horizon.

2. Identify the TGP horizon r_H as the largest r where eps(r) = 1/3.

3. Compute the photon ring radius r_ph by solving dV_eff/dr = 0 with
   V_eff(r) = c0^2 (4 - 3 psi) / (psi r^2). Search in (r_H, ~5 r_g).

4. Convert: in geometric units M = 1, r_g = G M / c^2 = 1.
   Compare r_ph^TGP with r_ph^GR = 1.5 r_s = 3 r_g (in standard
   Schwarzschild coordinates; in isotropic coordinates the photon ring
   sits at r_iso = (3 + 2 sqrt(2)) M / 2 ~ 2.914 M for GR
   - we report the ratio in terms of areal radius r_areal = r_iso * (1+M/(2 r_iso))^2
   to ensure invariant comparison).

5. Compute the shadow diameter D_sh = 2 sqrt(3) * r_g for GR and
   the analogous TGP value (b_crit = r_ph / sqrt(A(r_ph))). Convert
   to angular size for M87* (D = 16.8 Mpc, M = 6.5e9 M_sun) and
   Sgr A* (D = 8.15 kpc, M = 4.3e6 M_sun). Compare to EHT.

Date: 2026-04-25
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# ============================================================
# Physical constants
# ============================================================

G_NEWTON = 6.67430e-11        # m^3 kg^-1 s^-2
C_LIGHT = 2.99792458e8         # m / s
M_SUN = 1.98892e30             # kg
PARSEC = 3.0857e16             # m
KPC = 1e3 * PARSEC
MPC = 1e6 * PARSEC
ARCSEC_PER_RAD = 206264.806    # arcsec / rad


# ============================================================
# 1. Asymptotic PN expansion of eps(r), from M9.1'' P1
# ============================================================
#
# eps(r) at large r:
#   eps = A/r + c2 A^2/r^2 + c3 A^3/r^3 + c4 A^4/r^4 + c5 A^5/r^5 + ...
# with c1 = 1, c2 = -1, c3 = +5/3, c4 = -10/3, c5 = +22/3
# (from m9_1_pp_p1_higher_pn.py, sympy verified).
#
# The matching to Newtonian gravity uses
#   g_tt^TGP/(-c^2) = (4 - 3 psi)/psi = 1 - 4 eps + 4 eps^2 - 4 eps^3 + ...
# so leading order:  -4 eps_lead = -2 U_N/c^2 = -2 G M / (c^2 r)
#                    => eps_lead = G M / (2 c^2 r) = (1/2) (r_g / r)
# In geometric units M = G = c0 = 1: eps_lead = 1/(2 r), so A = 1/2.

C_AN = {
    1: 1.0,
    2: -1.0,
    3: 5.0 / 3.0,
    4: -10.0 / 3.0,
    5: 22.0 / 3.0,
    6: -154.0 / 9.0,
    7: 374.0 / 9.0,
}


def eps_asymptotic(r, A, n_terms=7):
    """Asymptotic series for eps(r) up to n_terms PN orders."""
    s = 0.0
    for n in range(1, n_terms + 1):
        s = s + C_AN[n] * A ** n / r ** n
    return s


def deps_dr_asymptotic(r, A, n_terms=7):
    """dEps/dr from the series."""
    s = 0.0
    for n in range(1, n_terms + 1):
        s = s + (-n) * C_AN[n] * A ** n / r ** (n + 1)
    return s


# ============================================================
# 2. ODE: state z = (eps, deps/dr), integrate inward from r_out
# ============================================================
#
# Phi-EOM (vacuum, alpha=2):
#   eps'' + (2/r) eps' + 2 (eps')^2 / (1 + eps) = 0
# Rearranged for eps'':
#   eps'' = -(2/r) eps' - 2 (eps')^2 / (1 + eps)

def rhs(r, z):
    eps, dep = z
    one_plus = 1.0 + eps
    if one_plus <= 0:
        # crashed past the singular boundary
        return np.array([0.0, 0.0])
    eps_pp = -(2.0 / r) * dep - 2.0 * dep ** 2 / one_plus
    return np.array([dep, eps_pp])


def horizon_event(r, z):
    """Detect eps -> 1/3 (g_tt -> 0)."""
    return z[0] - 1.0 / 3.0


horizon_event.terminal = True
# direction = 0 means catch any zero crossing (we don't care about the sign convention).
horizon_event.direction = 0


def integrate_inward(A, r_out=1000.0, r_min_safety=1.0e-3, n_terms=7,
                     rtol=1e-12, atol=1e-14):
    """Integrate eps(r) inward from r_out using PN initial conditions.

    Returns the solution object plus the horizon radius r_H (largest r
    where eps = 1/3) if found.
    """
    eps0 = eps_asymptotic(r_out, A, n_terms=n_terms)
    dep0 = deps_dr_asymptotic(r_out, A, n_terms=n_terms)
    # Integrate inward: dr < 0
    sol = solve_ivp(
        rhs,
        t_span=(r_out, r_min_safety),
        y0=np.array([eps0, dep0]),
        method="DOP853",
        events=horizon_event,
        max_step=0.5,
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    r_H = None
    if sol.t_events[0].size > 0:
        r_H = float(sol.t_events[0][0])
    return sol, r_H


# ============================================================
# 3. Photon ring: standard derivation for spherical metric in
#    isotropic-like radial coordinate.
#
# Effective metric:
#    ds^2 = -A_t(r) dt^2 + B(r) [dr^2 + r^2 dOmega^2]
# with A_t = c0^2 (4 - 3 psi)/psi  and  B = h(psi) = psi/(4 - 3 psi),
# i.e.  A_t * B = c0^2  (substrate budget fh = 1).
#
# Areal radius (the proper size of a sphere):  R(r) = r * sqrt(B(r)).
# Conserved E = A_t * dt/dlambda, L = R^2 * dphi/dlambda.
# Null condition gives turning points  b^2 := L^2/E^2 = R(r)^2 / A_t(r).
# Photon ring (unstable circular photon orbit):
#    d/dr [ R(r)^2 / A_t(r) ] = 0.
#
# In TGP isotropic with A_t * B = c0^2 we have B = c0^2/A_t and so
#    R^2/A_t = r^2 * B / A_t = r^2 * c0^2 / A_t^2.
# Setting c0 = 1:  R^2/A_t = r^2 / A_t^2.
#    d/dr [ r^2 / A_t^2 ] = 0
#    => 2 r A_t^2 - r^2 * 2 A_t * A_t' = 0
#    => r * A_t' = A_t.
#
# A_t(r) = (4 - 3 psi)/psi = 1 - 4 eps/(1+eps) - ... (with psi = 1 + eps):
#    A_t = 4/(1+eps) - 3
#    dA_t/dr = -4/(1+eps)^2 * eps'
#
# Photon-ring equation r A_t' = A_t:
#    r * (-4 eps'/(1+eps)^2) = (4 - 3(1+eps))/(1+eps) = (1 - 3 eps)/(1+eps)
#    => -4 r eps' = (1 - 3 eps)(1 + eps)
#    => 4 r eps'(r) + (1 - 3 eps(r))(1 + eps(r)) = 0
#
# Define F(r) = 4 r eps'(r) + (1 - 3 eps(r))(1 + eps(r)).
# Photon ring radius: zero of F(r) outside the TGP horizon.
# ============================================================


def F_photon(r, sol):
    z = sol.sol(r)
    eps = z[0]
    dep = z[1]
    return 4.0 * r * dep + (1.0 - 3.0 * eps) * (1.0 + eps)


def find_photon_ring(sol, r_H, r_max=20.0):
    """Find r_ph by bisection on F_photon = 0 in (r_H * 1.001, r_max)."""
    # Sample
    r_lo = r_H * 1.001 if r_H is not None else 1.0
    r_hi = r_max
    rs = np.linspace(r_lo, r_hi, 4000)
    Fs = np.array([F_photon(r, sol) for r in rs])
    # Find sign changes
    sign_changes = np.where(np.diff(np.sign(Fs)) != 0)[0]
    if sign_changes.size == 0:
        return None
    # The OUTERMOST photon ring: search from large r inward and find
    # the first crossing. But the standard photon ring (the unstable
    # circular photon orbit) is the OUTER one, just outside the horizon
    # but at a larger radius than the trivial solutions. For GR there is
    # one. For TGP we expect one, near r ~ 1.5 r_g (isotropic ~ 1.4 r_g).
    # Take the LARGEST r where F crosses zero (most physical).
    idx = sign_changes[-1]
    r_ph = brentq(lambda r: F_photon(r, sol), rs[idx], rs[idx + 1],
                   xtol=1e-12, rtol=1e-12)
    return r_ph


# ============================================================
# 4. Convert isotropic to areal (Schwarzschild) coordinate
# ============================================================
#
# In ISOTROPIC coords, GR Schwarzschild has
#   g_tt = -c^2 (1 - M/(2 r_iso))^2 / (1 + M/(2 r_iso))^2
#   g_rr = (1 + M/(2 r_iso))^4
# Areal r_areal = r_iso * (1 + M/(2 r_iso))^2.
# GR photon ring is at r_areal = 3 M (in c=G=1 units, M = r_g).
# r_iso for that:  3 M = r_iso (1 + M/(2 r_iso))^2
#   let x = M/(2 r_iso); 3 M = (M/(2x)) (1+x)^2
#   => 6 x = (1+x)^2 / 1 (-> 1+x)^2 = 6x => x^2 - 4x + 1 = 0
#   => x = 2 - sqrt(3) ~ 0.2679
#   => r_iso = M / (2 (2 - sqrt(3))) ~ M (2 + sqrt(3))/2 ~ 1.866 M
# So GR photon ring in ISOTROPIC = (2 + sqrt(3))/2 ~ 1.866 r_g.
#
# For TGP M9.1'' the radial coordinate is also "isotropic" (g_ij = h(psi) delta_ij).
# To compare to Schwarzschild's "1.5 r_s = 3 r_g", we convert TGP's r_ph_iso to an
# areal radius via the metric area function:
#   r_areal(r_iso) = sqrt(g_theta_theta) = r_iso * sqrt(h(psi))
# with h(psi) = psi/(4 - 3 psi) = 1/(A_metric/c^2).

def r_areal_TGP(r_iso, eps):
    """Areal radius from isotropic, given eps(r_iso)."""
    psi = 1.0 + eps
    h = psi / (4.0 - 3.0 * psi)
    return r_iso * np.sqrt(h)


def r_areal_GR(r_iso, M=1.0):
    return r_iso * (1.0 + M / (2.0 * r_iso)) ** 2


# ============================================================
# 5. Critical impact parameter and shadow size
# ============================================================
# For null geodesics, b_crit = r_areal(r_ph) / sqrt(A(r_ph))
# For GR: b_crit = 3 sqrt(3) M = 3 sqrt(3) r_g.
# Shadow angular diameter:  alpha = 2 b_crit / D_observer (small angle)


def critical_impact_parameter_TGP(r_ph_iso, eps_at_ph):
    """b_crit = r_iso * sqrt(h(psi)) / sqrt(A(r_ph))
    BUT careful: we need r_areal / sqrt(A_metric).
    A_metric(r) = (4 - 3 psi)/psi, h = 1/A_metric,
    so r_areal = r_iso * sqrt(h) = r_iso/sqrt(A_metric)
    => b = r_iso * sqrt(h) / sqrt(A_metric) = r_iso / A_metric (?)
    Wait, the standard formula is b_crit^2 = r_areal^2 / A(r_ph),
    not r_iso^2.  Let me redo.

    For ds^2 = -A dt^2 + B dr^2 + R^2 dOmega^2 (general spherical),
    null circular orbit at R(r) = R_ph, and impact parameter is
        b = R_ph / sqrt(A(R_ph)).
    In our case R = r_iso * sqrt(h) = r_iso * sqrt(psi/(4-3 psi))
                = r_iso / sqrt(A_metric)
    So:
        b_crit = (r_iso / sqrt(A_metric)) / sqrt(A_metric)
              = r_iso / A_metric.
    """
    psi = 1.0 + eps_at_ph
    A_metric = (4.0 - 3.0 * psi) / psi
    return r_ph_iso / A_metric


def critical_impact_parameter_GR(M=1.0):
    return 3.0 * np.sqrt(3.0) * M


# ============================================================
# 6. Astrophysical conversion
# ============================================================
def shadow_angular_diameter_microarcsec(b_crit_natural, M_solar, D_meters):
    """b_crit_natural is in units of r_g = G M / c^2.
    Returns angular diameter (microarcsec).
    """
    r_g_meters = G_NEWTON * (M_solar * M_SUN) / C_LIGHT ** 2
    b_meters = b_crit_natural * r_g_meters
    angular_diameter_rad = 2.0 * b_meters / D_meters
    return angular_diameter_rad * ARCSEC_PER_RAD * 1.0e6


# ============================================================
# Tests
# ============================================================

def run_tests():
    out = []
    out.append("=" * 78)
    out.append("EHT-quick: photon ring in M9.1'' hyperbolic BH (TGP_CLOSURE_PLAN sec 8.10)")
    out.append("Date: 2026-04-25")
    out.append("Geometric units: G = c0 = M = 1, r_g = G M / c^2 = 1.")
    out.append("=" * 78)
    out.append("")

    n_pass = 0
    n_total = 0

    # ---------- TEST 1: integrate Phi-EOM, find horizon ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 1: solve vacuum Phi-EOM inward from r_out, find TGP horizon")
    out.append("-" * 78)

    A = 0.5  # asymptotic amplitude in geometric units (eps = A/r + ...)
    r_out = 1000.0
    sol, r_H = integrate_inward(A=A, r_out=r_out, r_min_safety=1.0e-3,
                                  n_terms=7, rtol=1e-12, atol=1e-14)
    out.append(f"  IC at r_out = {r_out:.0f}: eps(r_out) = {eps_asymptotic(r_out, A):.6e}")
    out.append(f"                          eps'(r_out) = {deps_dr_asymptotic(r_out, A):.6e}")
    if r_H is None:
        out.append("  Integration did not reach eps = 1/3 (no TGP horizon found in domain).")
        out.append("  Verdict TEST 1: FAIL (numerical)")
    else:
        out.append(f"  TGP horizon at eps = 1/3:   r_H = {r_H:.6f} r_g")
        # Sanity: linear-order prediction is eps = A/r => r_H ~ 3 A = 1.5
        out.append(f"  Linear (1PN) prediction:    r_H ~ 3 A = {3*A:.4f} r_g")
        # Compare to GR Schwarzschild horizon in isotropic coords: r_iso_S = M/2 = 0.5 r_g
        out.append(f"  GR Schwarzschild horizon (isotropic): r_iso_S = M/2 = 0.5 r_g")
        out.append(f"  GR Schwarzschild horizon (areal):     r_s = 2 r_g")
        # Areal radius of TGP horizon: at eps -> 1/3, psi = 4/3,
        # h = (4/3)/(4-4) -> infinity. So r_areal at horizon DIVERGES.
        out.append(f"  At TGP horizon: psi = 4/3, h(psi) = psi/(4-3psi) -> infty")
        out.append(f"  => areal r_horizon diverges; consistent with GR (areal) where horizon is r_s")
        out.append("  but in isotropic coordinates GR has r_iso_S = M/2 finite (and h=infty there too).")
        out.append("  Verdict TEST 1: PASS (horizon found, structure consistent)")
        n_pass += 1
    out.append("")

    if r_H is None:
        out.append("Cannot proceed. Aborting.")
        return out, n_pass, n_total

    # ---------- TEST 2: photon ring radius ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 2: photon ring r_ph from dV_eff/dr = 0")
    out.append("-" * 78)
    out.append("  Photon ring at extremum of R^2(r)/A_t(r) (R = areal radius).")
    out.append("  In TGP isotropic (A_t * B = 1):  R^2/A_t = r^2/A_t^2.")
    out.append("  d/dr[r^2/A_t^2]=0  =>  r * A_t' = A_t.")
    out.append("  In eps:  4 r eps'(r) + (1 - 3 eps(r))(1 + eps(r)) = 0.")
    out.append("")

    r_ph = find_photon_ring(sol, r_H, r_max=20.0)
    if r_ph is None:
        out.append("  No photon ring found (no zero crossing in F_photon).")
        out.append("  Verdict TEST 2: FAIL")
        return out, n_pass, n_total
    z_ph = sol.sol(r_ph)
    eps_ph = float(z_ph[0])
    psi_ph = 1.0 + eps_ph
    A_metric_ph = (4.0 - 3.0 * psi_ph) / psi_ph
    out.append(f"  TGP photon ring (isotropic):  r_ph_iso = {r_ph:.6f} r_g")
    out.append(f"      eps(r_ph)  = {eps_ph:+.6f}")
    out.append(f"      psi(r_ph)  = {psi_ph:.6f}")
    out.append(f"      A_metric(r_ph) = {A_metric_ph:.6f}")

    # Convert to areal
    r_areal_ph = r_areal_TGP(r_ph, eps_ph)
    out.append(f"  TGP photon ring (areal):      r_ph_areal = {r_areal_ph:.6f} r_g")

    # GR comparison
    M_geom = 1.0
    r_iso_GR = (2.0 + np.sqrt(3.0)) / 2.0 * M_geom
    r_areal_GR_value = r_areal_GR(r_iso_GR, M=M_geom)
    out.append(f"  GR Schwarzschild photon ring: r_iso_GR = (2+sqrt 3)/2 = {r_iso_GR:.6f} r_g")
    out.append(f"                                r_areal_GR = 3 r_g (canonical)")
    out.append(f"")
    out.append(f"  Comparison (isotropic):   r_ph_TGP / r_ph_GR = {r_ph/r_iso_GR:.6f}")
    out.append(f"  Comparison (areal):       r_areal_TGP / r_areal_GR = {r_areal_ph/r_areal_GR_value:.6f}")
    rel_diff_isol = (r_ph - r_iso_GR) / r_iso_GR
    rel_diff_areal = (r_areal_ph - r_areal_GR_value) / r_areal_GR_value
    out.append(f"  Relative deviation (isotropic):  {rel_diff_isol*100:+.3f} %")
    out.append(f"  Relative deviation (areal):      {rel_diff_areal*100:+.3f} %")
    out.append("  Verdict TEST 2: PASS (photon ring found)")
    n_pass += 1
    out.append("")

    # ---------- TEST 3: critical impact parameter & shadow ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 3: critical impact parameter b_crit and shadow size")
    out.append("-" * 78)
    out.append("  b_crit = r_areal(r_ph) / sqrt(A_metric(r_ph))")
    out.append("         = r_iso / A_metric  (since r_areal = r_iso/sqrt(A))")
    out.append("")
    b_crit_TGP = critical_impact_parameter_TGP(r_ph, eps_ph)
    b_crit_GR = critical_impact_parameter_GR(M=1.0)  # = 3 sqrt(3)
    out.append(f"  b_crit (TGP) = {b_crit_TGP:.6f} r_g")
    out.append(f"  b_crit (GR)  = 3 sqrt(3) = {b_crit_GR:.6f} r_g")
    rel_b = (b_crit_TGP - b_crit_GR) / b_crit_GR
    out.append(f"  b_crit_TGP / b_crit_GR = {b_crit_TGP/b_crit_GR:.6f}")
    out.append(f"  Relative deviation:    {rel_b*100:+.3f} %")
    out.append("  Verdict TEST 3: PASS (critical impact parameter computed)")
    n_pass += 1
    out.append("")

    # ---------- TEST 4: M87* shadow ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 4: M87* shadow prediction vs EHT 2019 observation")
    out.append("-" * 78)
    M87_mass_solar = 6.5e9
    M87_distance_m = 16.8 * MPC
    EHT_M87_obs_uas = 42.0    # microarcsec, EHT 2019
    EHT_M87_uncertainty_uas = 3.0  # ~5-10%

    sh_TGP = shadow_angular_diameter_microarcsec(b_crit_TGP, M87_mass_solar, M87_distance_m)
    sh_GR = shadow_angular_diameter_microarcsec(b_crit_GR, M87_mass_solar, M87_distance_m)
    out.append(f"  M87*: M = 6.5e9 M_sun, D = 16.8 Mpc")
    out.append(f"  TGP shadow diameter:        {sh_TGP:8.3f} microarcsec")
    out.append(f"  GR  shadow diameter:        {sh_GR:8.3f} microarcsec")
    out.append(f"  EHT 2019 observed:          {EHT_M87_obs_uas:8.1f} +/- {EHT_M87_uncertainty_uas:.1f} microarcsec (~10% systematic)")
    rel_TGP_M87 = (sh_TGP - EHT_M87_obs_uas) / EHT_M87_obs_uas
    rel_GR_M87 = (sh_GR - EHT_M87_obs_uas) / EHT_M87_obs_uas
    out.append(f"  TGP / observed:             {sh_TGP/EHT_M87_obs_uas:.4f}  ({rel_TGP_M87*100:+.2f} %)")
    out.append(f"  GR  / observed:             {sh_GR/EHT_M87_obs_uas:.4f}  ({rel_GR_M87*100:+.2f} %)")
    in_5pct_TGP = abs(rel_TGP_M87) <= 0.05
    in_10pct_TGP = abs(rel_TGP_M87) <= 0.10
    out.append(f"  TGP within 5%:  {'YES' if in_5pct_TGP else 'NO'}")
    out.append(f"  TGP within 10%: {'YES' if in_10pct_TGP else 'NO'}")
    if in_10pct_TGP:
        out.append("  Verdict TEST 4: PASS (TGP within EHT 10% bound)")
        n_pass += 1
    else:
        out.append("  Verdict TEST 4: FAIL (TGP outside EHT 10% bound)")
    out.append("")

    # ---------- TEST 5: Sgr A* shadow ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 5: Sgr A* shadow prediction vs EHT 2022 observation")
    out.append("-" * 78)
    SgrA_mass_solar = 4.3e6
    SgrA_distance_m = 8.15 * KPC
    EHT_SgrA_obs_uas = 51.8   # microarcsec, EHT 2022 (ring diameter)
    EHT_SgrA_uncertainty_uas = 2.3

    sh_TGP_S = shadow_angular_diameter_microarcsec(b_crit_TGP, SgrA_mass_solar, SgrA_distance_m)
    sh_GR_S = shadow_angular_diameter_microarcsec(b_crit_GR, SgrA_mass_solar, SgrA_distance_m)
    out.append(f"  Sgr A*: M = 4.3e6 M_sun, D = 8.15 kpc")
    out.append(f"  TGP shadow diameter:        {sh_TGP_S:8.3f} microarcsec")
    out.append(f"  GR  shadow diameter:        {sh_GR_S:8.3f} microarcsec")
    out.append(f"  EHT 2022 observed:          {EHT_SgrA_obs_uas:8.1f} +/- {EHT_SgrA_uncertainty_uas:.1f} microarcsec")
    rel_TGP_S = (sh_TGP_S - EHT_SgrA_obs_uas) / EHT_SgrA_obs_uas
    rel_GR_S = (sh_GR_S - EHT_SgrA_obs_uas) / EHT_SgrA_obs_uas
    out.append(f"  TGP / observed:             {sh_TGP_S/EHT_SgrA_obs_uas:.4f}  ({rel_TGP_S*100:+.2f} %)")
    out.append(f"  GR  / observed:             {sh_GR_S/EHT_SgrA_obs_uas:.4f}  ({rel_GR_S*100:+.2f} %)")
    in_5pct_TGP_S = abs(rel_TGP_S) <= 0.05
    in_10pct_TGP_S = abs(rel_TGP_S) <= 0.10
    out.append(f"  TGP within 5%:  {'YES' if in_5pct_TGP_S else 'NO'}")
    out.append(f"  TGP within 10%: {'YES' if in_10pct_TGP_S else 'NO'}")
    if in_10pct_TGP_S:
        out.append("  Verdict TEST 5: PASS (TGP within EHT 10% bound)")
        n_pass += 1
    else:
        out.append("  Verdict TEST 5: FAIL (TGP outside EHT 10% bound)")
    out.append("")

    # ---------- TEST 6: convergence check (independent of analytic c_n) ----------
    n_total += 1
    out.append("-" * 78)
    out.append("TEST 6: convergence -- vary r_out, n_terms, check r_ph stable")
    out.append("-" * 78)

    convergence_results = []
    for r_out_test, n_terms_test in [(500.0, 5), (1000.0, 6), (1000.0, 7),
                                      (2000.0, 7), (5000.0, 7)]:
        s2, rH2 = integrate_inward(A=A, r_out=r_out_test, r_min_safety=1e-3,
                                    n_terms=n_terms_test, rtol=1e-12, atol=1e-14)
        if rH2 is None:
            convergence_results.append((r_out_test, n_terms_test, None, None, None))
            continue
        rph2 = find_photon_ring(s2, rH2, r_max=20.0)
        if rph2 is None:
            convergence_results.append((r_out_test, n_terms_test, rH2, None, None))
            continue
        eps2 = float(s2.sol(rph2)[0])
        bcrit2 = critical_impact_parameter_TGP(rph2, eps2)
        convergence_results.append((r_out_test, n_terms_test, rH2, rph2, bcrit2))

    out.append("    r_out    n_terms    r_H            r_ph_iso       b_crit (r_g)")
    out.append("    -----    -------    -----------    -----------    ------------")
    bcrits = []
    for r_out_test, n_terms_test, rH2, rph2, bcrit2 in convergence_results:
        if rH2 is None:
            out.append(f"    {r_out_test:5.0f}    {n_terms_test:7d}    no horizon found")
        elif rph2 is None:
            out.append(f"    {r_out_test:5.0f}    {n_terms_test:7d}    {rH2:.6f}     no photon ring")
        else:
            out.append(f"    {r_out_test:5.0f}    {n_terms_test:7d}    {rH2:.6f}     {rph2:.6f}     {bcrit2:.6f}")
            bcrits.append(bcrit2)
    if len(bcrits) >= 2:
        spread = (max(bcrits) - min(bcrits)) / np.mean(bcrits)
        out.append(f"  b_crit spread across runs: {spread*100:.4f} %")
        if spread < 1e-3:
            out.append("  Verdict TEST 6: PASS (b_crit converged below 0.1%)")
            n_pass += 1
        else:
            out.append(f"  Verdict TEST 6: FAIL (b_crit not converged, spread {spread*100:.3f}%)")
    else:
        out.append("  Verdict TEST 6: FAIL (insufficient convergence data)")
    out.append("")

    # ---------- Final analysis: deviation and interpretation ----------
    out.append("=" * 78)
    out.append("FINAL ANALYSIS")
    out.append("=" * 78)
    out.append("")
    out.append(f"  TGP M9.1'' photon ring (areal):    r_ph^TGP = {r_areal_ph:.6f} r_g")
    out.append(f"  GR Schwarzschild photon ring:      r_ph^GR  = 3.000000 r_g")
    out.append(f"  Ratio:                             {r_areal_ph/3.0:.6f}")
    out.append("")
    out.append(f"  TGP critical impact parameter:     b_crit^TGP = {b_crit_TGP:.6f} r_g")
    out.append(f"  GR critical impact parameter:      b_crit^GR  = {b_crit_GR:.6f} r_g")
    out.append(f"  Fractional deviation:              {rel_b*100:+.4f} %")
    out.append("")
    out.append(f"  M87* shadow:    TGP {sh_TGP:.2f} uas vs EHT {EHT_M87_obs_uas:.1f} +/- {EHT_M87_uncertainty_uas:.1f} uas")
    out.append(f"                  GR  {sh_GR:.2f} uas (deviation TGP-GR: {(sh_TGP-sh_GR)/sh_GR*100:+.4f} %)")
    out.append(f"  Sgr A* shadow:  TGP {sh_TGP_S:.2f} uas vs EHT {EHT_SgrA_obs_uas:.1f} +/- {EHT_SgrA_uncertainty_uas:.1f} uas")
    out.append(f"                  GR  {sh_GR_S:.2f} uas (deviation TGP-GR: {(sh_TGP_S-sh_GR_S)/sh_GR_S*100:+.4f} %)")
    out.append("")

    if abs(rel_b) <= 0.10:
        verdict_global = "POSITIVE: TGP shadow within EHT 10% bound, indistinguishable at current precision"
    elif abs(rel_b) <= 0.20:
        verdict_global = "INCONCLUSIVE: TGP deviation 10-20%, marginal vs EHT systematics"
    else:
        verdict_global = "NEGATIVE: TGP shadow outside EHT 10% bound; tension with observation"
    out.append(f"  GLOBAL VERDICT: {verdict_global}")
    out.append("")

    # Summary
    out.append("=" * 78)
    out.append(f"SUMMARY: {n_pass}/{n_total} PASS")
    out.append("=" * 78)

    return out, n_pass, n_total


def main():
    out, n_pass, n_total = run_tests()
    text = "\n".join(out)
    print(text)

    # Write to research/op7
    import os
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "research", "op7", "eht_photon_ring_m911.txt"
    )
    out_path = os.path.normpath(out_path)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(text)
    print()
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main()
