#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qnm_spectrum.py -- TGP Quasi-Normal Mode Spectrum (O17)
========================================================
Numerically computes the QNM spectrum for TGP (Teoria Generowanej Przestrzeni)
around compact objects and compares with standard GR Schwarzschild QNMs.

Theory reference: dodatekC_ringdown.tex, sek08

Ringdown equation (Theorem thm:ringdown, eq:ringdown):

    d^2 Phi_w / dr*^2 + [omega^2 - V_eff(r)] Phi_w = 0

Effective potential (eq:V-eff):

    V_eff(r) = (c0^2 Phi0 / f(r)) * [
        l(l+1)/r^2 + (9/4) f''/f + (59/16)(f'/f)^2
        + 4 f'/(r f) - 2 beta chi + 3 gamma chi^2
    ]

    where chi(r) = f(r)/Phi0

Tortoise coordinate (eq:tortoise):

    dr*/dr = sqrt(chi(r)) / c0 = sqrt(f(r)/Phi0) / c0

Background profile (eq:bg-static):

    f'' + (2/r)f' + alpha (f')^2/f + beta f^2/Phi0 - gamma f^3/Phi0^2 = 0

    with f'(0)=0, f(r->inf) = Phi0, alpha=2, beta=gamma (vacuum)

Boundary conditions for QNMs (Definition def:QNM):

    r* -> +inf:     outgoing, Phi_w ~ exp(+i k_inf r*),  k_inf = sqrt(omega^2 - c0^2 gamma)
    r* -> r*_min:   ingoing into frozen zone (c -> 0)

Key TGP features vs GR (Remark rem:vs-RW):
    1. Mass gap:  omega_min = c0 sqrt(gamma), no low-frequency modes
    2. Rising barrier: V_eff ~ 3 gamma c0^2 chi  for f >> Phi0
    3. Scalar breathing modes (l=0 allowed, unlike GR tensor l>=2)

Parameters from theory (sek04, sek05):
    beta = gamma  (vacuum condition)
    q = 8 pi G0 / c0^2
    gamma ~ 10^{-52} m^{-2}  (cosmological mass gap)
    Phi0 ~ 25

GR comparison:
    omega_GR ~ (1/r_s) [l + 1/2 - i(n+1/2)/sqrt(27)]   (Schutz-Will WKB)
    Regge-Wheeler potential: V_RW = (1-r_s/r)[l(l+1)/r^2 + (1-s^2)r_s/r^3]

Output:
    scripts/plots/qnm_spectrum.png  -- V_eff comparison, complex spectrum, deviations
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import os
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===================================================================
# PHYSICAL CONSTANTS (SI)
# ===================================================================
G_SI = 6.67430e-11       # m^3 kg^-1 s^-2
c_SI = 2.99792458e8      # m/s
M_sun = 1.98892e30       # kg

# ===================================================================
# TGP PARAMETERS (sek04_stale.tex, sek05_ciemna_energia.tex)
# ===================================================================
PHI0 = 25.0              # Background field (dimensionless)
C0 = c_SI                # Asymptotic speed of light
GAMMA_PHYS = 1.0e-52     # m^{-2}  (cosmological mass gap)
BETA_PHYS = GAMMA_PHYS   # vacuum condition beta = gamma
ALPHA = 2.0              # Fixed by theory (sek02_pole.tex)
Q_COUPLING = 8 * np.pi * G_SI / C0**2  # scalar charge (sek08)

# BH configurations for analysis
BH_MASSES_MSUN = [1.0, 3.0, 5.0, 10.0, 30.0, 50.0, 100.0]

BH_CONFIGS = [
    {"M_Msun": 1.0,  "label": r"$1\,M_\odot$"},
    {"M_Msun": 10.0, "label": r"$10\,M_\odot$"},
]

ELL = 2  # Dominant GW multipole

# Dimensionless coupling values for numerical study
# Physical gamma_hat ~ 10^{-45} is negligible; we also study enhanced values
GAMMA_HAT_VALUES = [0.001, 0.01, 0.1, 0.5]

# Source strengths (compactness parameter)
SOURCE_STRENGTHS = [0.5, 2.0]


# ===================================================================
# HELPERS
# ===================================================================

def schwarzschild_radius(M_Msun):
    """Schwarzschild radius r_s = 2GM/c^2 in meters."""
    return 2 * G_SI * M_Msun * M_sun / c_SI**2


def dimensionless_gamma(M_Msun):
    """gamma_hat = gamma * r_s^2 -- dimensionless TGP coupling at BH scale."""
    r_s = schwarzschild_radius(M_Msun)
    return GAMMA_PHYS * r_s**2


# ===================================================================
# 1. GR SCHWARZSCHILD QNM FREQUENCIES (tabulated + WKB)
# ===================================================================

def qnm_gr_schwarzschild(M_Msun, ell=2, n=0):
    """
    GR Schwarzschild QNM frequencies.

    Tabulated values from Berti, Cardoso, Starinets (2009).
    Dimensionless: omega * GM/c^3.

    WKB fallback (Schutz & Will 1985):
        omega * M_geo ~ (l + 1/2)/2 - i(n + 1/2)/(2 sqrt(27))

    Physical: omega_phys = omega_dimless * c^3 / (G M)
    """
    r_s = schwarzschild_radius(M_Msun)

    # Tabulated: omega_dimless * (GM/c^3)
    # Tensor (s=2):
    gr_tensor = {
        (2, 0): complex(0.37367, -0.08896),
        (2, 1): complex(0.34671, -0.27392),
        (2, 2): complex(0.30105, -0.47828),
        (3, 0): complex(0.59944, -0.09270),
        (3, 1): complex(0.58264, -0.28130),
    }
    # Scalar (s=0) -- relevant for TGP breathing mode:
    gr_scalar = {
        (2, 0): complex(0.48364, -0.09676),
        (2, 1): complex(0.46383, -0.29561),
        (0, 0): complex(0.11049, -0.10049),
        (1, 0): complex(0.29293, -0.09766),
    }

    # WKB fallback
    omega_wkb = complex((ell + 0.5) / 2.0, -(n + 0.5) / (2.0 * np.sqrt(27)))

    key = (ell, n)
    omega_tensor_dl = gr_tensor.get(key, omega_wkb)
    omega_scalar_dl = gr_scalar.get(key, omega_wkb)

    # Convert to physical units
    conv = c_SI**3 / (G_SI * M_Msun * M_sun)
    omega_tensor = omega_tensor_dl * conv
    omega_scalar = omega_scalar_dl * conv

    return {
        'omega_tensor': omega_tensor,
        'omega_scalar': omega_scalar,
        'omega_tensor_dimless': omega_tensor_dl,
        'omega_scalar_dimless': omega_scalar_dl,
        'f_tensor_Hz': omega_tensor.real / (2 * np.pi),
        'f_scalar_Hz': omega_scalar.real / (2 * np.pi),
        'r_s': r_s,
    }


# ===================================================================
# 2. GR REGGE-WHEELER POTENTIAL
# ===================================================================

def regge_wheeler_potential(r_over_rs, ell=2, spin=0):
    """
    Regge-Wheeler potential (units of 1/r_s^2):
        V_RW = (1 - 1/x) [l(l+1)/x^2 + (1-s^2)/x^3]
    where x = r/r_s.
    """
    x = r_over_rs
    return (1 - 1.0 / x) * (ell * (ell + 1) / x**2 + (1 - spin**2) / x**3)


def tortoise_schwarzschild(r_over_rs):
    """r*/r_s = r/r_s + ln(r/r_s - 1)"""
    x = r_over_rs
    return x + np.log(np.abs(x - 1.0))


# ===================================================================
# 3. TGP BACKGROUND AND EFFECTIVE POTENTIAL
# ===================================================================

def solve_tgp_background(gamma_hat, source_strength, sigma=0.5,
                         x_min=0.02, x_max=80.0, n_mesh=1200):
    """
    Solve dimensionless TGP background (eq:bg-static):
        f'' + (2/x)f' + alpha(f')^2/f + beta_hat f^2 - gamma_hat f^3 = source
    with beta_hat=gamma_hat, f'(0)=0, f(x_max)=1.
    Returns (x, f, fp, fpp) or None.
    """
    beta_hat = gamma_hat

    def source_gaussian(x, C, sig):
        return C * np.exp(-x**2 / (2 * sig**2)) / (2 * np.pi * sig**2)**1.5

    def bg_ode(x, y):
        f = np.maximum(y[0], 1e-15)
        fp = y[1]
        x_s = np.maximum(x, 1e-10)
        src = -source_gaussian(x, source_strength, sigma)
        fpp = (-(2.0 / x_s) * fp - ALPHA * fp**2 / f
               - beta_hat * f**2 + gamma_hat * f**3 + src)
        return np.vstack([fp, fpp])

    def bg_bc(ya, yb):
        return np.array([ya[1], yb[0] - 1.0])

    x_mesh = np.linspace(x_min, x_max, n_mesh)
    peak = source_strength / (4 * np.pi * sigma)
    f_g = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
    fp_g = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
    y_guess = np.vstack([f_g, fp_g])

    for tol, nodes in [(1e-6, 60000), (1e-4, 120000), (1e-3, 200000)]:
        sol = solve_bvp(fun=bg_ode, bc=bg_bc, x=x_mesh, y=y_guess,
                        tol=tol, max_nodes=nodes, verbose=0)
        if sol.success:
            break
        x_mesh = np.linspace(x_min, x_max, 2 * len(x_mesh))
        f_g = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
        fp_g = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
        y_guess = np.vstack([f_g, fp_g])

    if not sol.success:
        return None

    x = sol.x
    f = sol.y[0]
    fp = sol.y[1]
    rhs = bg_ode(x, sol.y)
    fpp = rhs[1]
    return x, f, fp, fpp


def compute_V_eff_tgp(x, f, fp, fpp, ell, gamma_hat):
    """
    TGP effective potential (eq:V-eff, dimensionless units 1/r_s^2):
        V_eff = (1/f) * [l(l+1)/r^2 + (9/4)f''/f + (59/16)(f'/f)^2
                         + 4f'/(rf) - 2 beta_hat chi + 3 gamma_hat chi^2]
    with chi = f (since Phi0=1 in dimensionless units), beta_hat = gamma_hat.
    """
    beta_hat = gamma_hat
    f_s = np.maximum(f, 1e-15)
    x_s = np.maximum(x, 1e-10)
    chi = f_s
    bracket = (ell * (ell + 1) / x_s**2
               + (9.0 / 4.0) * fpp / f_s
               + (59.0 / 16.0) * (fp / f_s)**2
               + 4.0 * fp / (x_s * f_s)
               - 2.0 * beta_hat * chi
               + 3.0 * gamma_hat * chi**2)
    return bracket / f_s


def compute_tortoise_tgp(x, f):
    """TGP tortoise: dr*/dr = sqrt(f) (dimensionless, eq:tortoise)."""
    integrand = np.sqrt(np.maximum(f, 1e-15))
    r_star = np.zeros_like(x)
    for i in range(1, len(x)):
        r_star[i] = (r_star[i - 1]
                     + 0.5 * (integrand[i] + integrand[i - 1]) * (x[i] - x[i - 1]))
    return r_star


def build_tgp_potential(x, f, fp, fpp, ell, gamma_hat):
    """Build V_eff(r*) spline. Returns (spline, rs_clean, V_inf) or (None,None,None)."""
    V = compute_V_eff_tgp(x, f, fp, fpp, ell, gamma_hat)
    r_star = compute_tortoise_tgp(x, f)

    mask = np.isfinite(V) & np.isfinite(r_star)
    rs = r_star[mask]
    Vc = V[mask]

    dr = np.diff(rs)
    mono = np.concatenate([[True], dr > 1e-15])
    rs, Vc = rs[mono], Vc[mono]

    skip = min(5, len(rs) // 10)
    rs, Vc = rs[skip:], Vc[skip:]

    if len(rs) < 30:
        return None, None, None

    cs = CubicSpline(rs, Vc, extrapolate=True)
    # V_inf = gamma_hat in weak field (Corollary cor:weak-ringdown)
    V_inf = gamma_hat
    return cs, rs, V_inf


# ===================================================================
# 4. SHOOTING METHOD FOR QNM
# ===================================================================

def shoot_qnm(omega, V_spline, rs_min, rs_max, V_inf):
    """
    Integrate Phi'' + [omega^2 - V(r*)] Phi = 0 from rs_min to rs_max.
    BC: Phi(rs_min)=0, Phi'(rs_min)=1  (evanescent at rising barrier).
    Returns (A_in, A_out).  QNM: A_in = 0.  (Definition def:QNM)
    """
    omega2 = omega**2

    def ode(rs, y):
        V = float(V_spline(rs))
        return [y[1], -(omega2 - V) * y[0]]

    y0 = [0.0 + 0.0j, 1.0 + 0.0j]
    try:
        sol = solve_ivp(ode, [rs_min, rs_max], y0, method='RK45',
                        rtol=1e-9, atol=1e-11,
                        max_step=(rs_max - rs_min) / 300)
        if not sol.success:
            return None, None
    except Exception:
        return None, None

    Phi_end = sol.y[0, -1]
    dPhi_end = sol.y[1, -1]

    k_inf = np.sqrt(omega2 - V_inf + 0j)
    if k_inf.real < 0:
        k_inf = -k_inf
    if abs(k_inf) < 1e-15:
        return None, None

    exp_p = np.exp(1j * k_inf * rs_max)
    exp_m = np.exp(-1j * k_inf * rs_max)
    A_out = (Phi_end * 1j * k_inf + dPhi_end) / (2 * 1j * k_inf * exp_p)
    A_in = (Phi_end * 1j * k_inf - dPhi_end) / (2 * 1j * k_inf * exp_m)
    return A_in, A_out


def find_qnm_newton(omega_guess, V_spline, rs_min, rs_max, V_inf,
                     max_iter=120, tol=1e-8, d_omega=1e-5):
    """Newton-Raphson for A_in(omega)=0 in complex omega plane."""
    omega = omega_guess
    for it in range(max_iter):
        A_in, A_out = shoot_qnm(omega, V_spline, rs_min, rs_max, V_inf)
        if A_in is None:
            return None, it, "integration failed"
        F = A_in / A_out if (A_out is not None and abs(A_out) > 1e-30) else A_in
        if abs(F) < tol:
            return omega, it, "converged"
        A_in_d, A_out_d = shoot_qnm(omega + d_omega, V_spline, rs_min, rs_max, V_inf)
        if A_in_d is None:
            return None, it, "Jacobian failed"
        F_d = A_in_d / A_out_d if (A_out_d is not None and abs(A_out_d) > 1e-30) else A_in_d
        dF = (F_d - F) / d_omega
        if abs(dF) < 1e-30:
            return None, it, "zero derivative"
        delta = F / dF
        step = min(1.0, 0.5 / max(abs(delta) / (abs(omega) + 0.1), 1.0))
        omega_new = omega - step * delta
        if omega_new.imag > 0:
            omega_new = omega_new.real - abs(omega_new.imag) * 1j
        omega = omega_new
    return omega, max_iter, "max iterations"


def wkb_qnm_tgp(V_spline, rs_arr, V_inf, n=0):
    """
    WKB estimate for TGP QNMs.

    The TGP potential (eq:V-eff) is fundamentally different from GR:
    - It rises monotonically toward the compact object (Corollary cor:strong-ringdown)
    - At infinity it floors at V_inf = gamma_hat (mass gap, Corollary cor:weak-ringdown)

    For a potential with a local maximum (barrier), standard WKB gives:
        omega^2 = V_peak - i(n+1/2) sqrt(-2 V''_peak)

    If the potential is monotonically decreasing (typical for TGP strong field),
    we look for a shoulder/inflection region where V transitions from steep
    to flat.  The QNM frequency is set by the curvature at this transition.

    Returns complex omega.
    """
    rs_test = np.linspace(rs_arr[0], rs_arr[-1], 3000)
    V_test = V_spline(rs_test)
    finite = np.isfinite(V_test)
    rs_test, V_test = rs_test[finite], V_test[finite]

    if len(rs_test) < 10:
        return None

    # Check for local maximum (GR-like barrier)
    dV = np.diff(V_test)
    sign_changes = np.where(dV[:-1] > 0)[0]  # rising followed by falling?
    has_peak = False
    for idx in sign_changes:
        if idx + 1 < len(dV) and dV[idx + 1] < 0:
            has_peak = True
            break

    if has_peak:
        # Standard WKB at peak
        idx_max = np.argmax(V_test)
        rs_peak = rs_test[idx_max]
        V_peak = V_test[idx_max]
        V_dd = float(V_spline.derivative(2)(rs_peak))
        if V_dd >= 0:
            dr = (rs_arr[-1] - rs_arr[0]) / 5000
            V_dd = float((V_spline(rs_peak + dr) + V_spline(rs_peak - dr)
                          - 2 * V_peak) / dr**2)
        if V_dd < 0 and V_peak > 0:
            omega_sq = V_peak - 1j * (n + 0.5) * np.sqrt(-2 * V_dd)
            omega = np.sqrt(omega_sq + 0j)
            if omega.real < 0:
                omega = -omega
            if omega.imag > 0:
                omega = omega.conjugate()
            return omega

    # Monotonically decreasing TGP potential: QNMs are "box modes"
    # trapped between the rising inner barrier and the mass-gap floor.
    # For a potential well with floor V_inf and steep inner wall,
    # approximate as particle-in-box with soft walls:
    #
    # omega_n^2 ~ V_inf + (n + 1/2)^2 * pi^2 / L_eff^2
    #
    # where L_eff is the effective width where V(r*) < some threshold V_th.
    # Damping: omega_I ~ -sqrt(V_inf) / Q_eff

    # Find effective cavity width at V_th = 2 * V_inf
    V_threshold = max(3.0 * V_inf, V_inf + 0.1)
    inside = V_test < V_threshold
    if inside.any():
        rs_inside = rs_test[inside]
        L_eff = rs_inside[-1] - rs_inside[0]
    else:
        L_eff = rs_test[-1] - rs_test[0]

    L_eff = max(L_eff, 1.0)

    omega_R_sq = V_inf + ((n + 0.5) * np.pi / L_eff)**2
    omega_R = np.sqrt(max(omega_R_sq, V_inf + 0.01))

    # Damping estimate: leakage through the mass-gap floor
    # omega_I ~ -1/(2 L_eff) for lowest modes (rough estimate)
    omega_I = -(n + 0.5) / (2.0 * L_eff)

    return complex(omega_R, omega_I)


def find_qnm_modes(V_spline, rs_arr, V_inf, ell, label, n_modes=2):
    """
    Find QNMs: WKB initial guess + Newton refinement + grid search.
    Returns list of (n, omega, converged_bool).
    """
    rs_min = rs_arr[0]
    rs_max = rs_arr[-1]

    # Truncate outer domain where V is flat
    test_rs = np.linspace(rs_arr[0], rs_arr[-1], 600)
    test_V = V_spline(test_rs)
    for k in range(len(test_rs) - 1, 0, -1):
        if abs(test_V[k] - test_V[-1]) > 0.01:
            rs_max = min(test_rs[k] + 15.0, rs_arr[-1])
            break

    # Handle l=0 V<0 region
    if ell == 0:
        test_V2 = V_spline(np.linspace(rs_min, rs_max, 800))
        test_rs2 = np.linspace(rs_min, rs_max, 800)
        for k in range(len(test_rs2) - 1):
            if test_V2[k] < 0 and test_V2[k + 1] >= 0:
                rs_min = test_rs2[k + 1]

    results = []
    for n in range(n_modes):
        # WKB initial guess
        omega_wkb = wkb_qnm_tgp(V_spline, rs_arr, V_inf, n=n)
        if omega_wkb is None:
            continue

        # Newton refinement
        omega_qnm, n_it, status = find_qnm_newton(
            omega_wkb, V_spline, rs_min, rs_max, V_inf,
            max_iter=100, tol=1e-9)

        if omega_qnm is not None and status == "converged":
            results.append((n, omega_qnm, True))
        else:
            results.append((n, omega_wkb, False))

    # Grid search as backup
    omega_gap = np.sqrt(max(V_inf, 1e-10))
    V_test = V_spline(np.linspace(rs_min, rs_max, 400))
    V_max_est = np.max(V_test[np.isfinite(V_test)])
    for n in range(n_modes):
        already = any(r[0] == n and r[2] for r in results)
        if already:
            continue
        for wr_try in np.linspace(omega_gap * 0.9, min(np.sqrt(V_max_est), 5.0), 6):
            found = False
            for wi_try in [-0.02, -0.05, -0.1, -0.3]:
                omega_guess = complex(wr_try, wi_try * (n + 1))
                omega_qnm, _, status = find_qnm_newton(
                    omega_guess, V_spline, rs_min, rs_max, V_inf,
                    max_iter=80, tol=1e-8)
                if omega_qnm is not None and status == "converged" and omega_qnm.real > 0.01:
                    # Check not duplicate
                    is_dup = any(abs(omega_qnm.real - r[1].real) / max(r[1].real, 0.01) < 0.02
                                for r in results if r[2])
                    if not is_dup:
                        results.append((n, omega_qnm, True))
                        found = True
                        break
            if found:
                break

    return results


# ===================================================================
# 5. COMPUTE FULL QNM ANALYSIS FOR A BH MASS
# ===================================================================

def compute_tgp_qnm_for_bh(M_Msun, ell=2):
    """
    Compute TGP QNM spectrum for a given BH mass.

    Physical gamma_hat ~ 10^{-45} gives negligible deviations.
    We also compute with enhanced gamma_hat values to map
    the sensitivity/detectability curve.

    Returns dict with GR and TGP results.
    """
    r_s = schwarzschild_radius(M_Msun)
    gamma_hat_phys = dimensionless_gamma(M_Msun)

    # GR reference
    gr_results = {}
    for n in range(2):
        gr_results[n] = qnm_gr_schwarzschild(M_Msun, ell=ell, n=n)

    # TGP with various gamma_hat values
    tgp_results = {}
    for gamma_hat in GAMMA_HAT_VALUES:
        for C in SOURCE_STRENGTHS:
            key = (gamma_hat, C)
            bg = solve_tgp_background(gamma_hat, C)
            if bg is None:
                continue
            x, f, fp, fpp = bg
            V_sp, rs_arr, V_inf = build_tgp_potential(x, f, fp, fpp, ell, gamma_hat)
            if V_sp is None:
                continue
            modes = find_qnm_modes(V_sp, rs_arr, V_inf, ell,
                                   f"M={M_Msun}, gh={gamma_hat:.1e}, C={C}")
            tgp_results[key] = {
                'modes': modes, 'V_spline': V_sp, 'rs_arr': rs_arr,
                'V_inf': V_inf, 'x': x, 'f': f, 'f_max': f.max(),
                'gamma_hat': gamma_hat, 'source': C,
            }

    return {
        'M_Msun': M_Msun, 'r_s': r_s,
        'gamma_hat_phys': gamma_hat_phys,
        'gr': gr_results, 'tgp': tgp_results,
    }


# ===================================================================
# 6. ANALYTIC DEVIATION ESTIMATE
# ===================================================================

def analytic_deviation(M_Msun, ell=2, n=0):
    """
    Analytic estimate of TGP QNM deviation from GR.

    From V_eff (eq:V-eff), the TGP correction in weak field (chi~1):
        delta_V ~ gamma_hat  (since 3 gamma chi^2 - 2 beta chi = gamma for chi=1)

    WKB: delta_omega/omega ~ delta_V_peak / (2 V_peak)
    For l=2 Schwarzschild: V_peak ~ 0.15/r_s^2
    => delta_omega/omega ~ gamma_hat / 0.3

    For strong-field amplification (Corollary cor:strong-ringdown):
        delta_V ~ 3 gamma_hat chi_max^2
        => delta_omega/omega ~ gamma_hat * chi_max^2 / 0.1

    Returns dict with estimates.
    """
    gamma_hat = dimensionless_gamma(M_Msun)
    r_s = schwarzschild_radius(M_Msun)
    delta_frac = gamma_hat / 0.3
    # chi_max needed for 1% deviation
    chi_max_1pct = np.sqrt(0.01 / max(gamma_hat, 1e-100))
    return {
        'gamma_hat': gamma_hat, 'r_s_m': r_s,
        'delta_omega_frac': delta_frac,
        'chi_max_1pct': chi_max_1pct,
    }


# ===================================================================
# 7. COMPUTE DEVIATION TABLE
# ===================================================================

def compute_deviations(all_results, ell=2):
    """
    Compute delta_omega/omega_GR for each TGP configuration vs GR.
    """
    table = []
    for res in all_results:
        M = res['M_Msun']
        for (gamma_hat, C), tgp in res['tgp'].items():
            for n_mode, omega_tgp, converged in tgp['modes']:
                gr = qnm_gr_schwarzschild(M, ell=ell, n=n_mode)
                # GR omega in r_s units: omega_phys * r_s/c = omega_dimless * 2
                omega_gr_rs = 2 * gr['omega_scalar_dimless']
                delta_re = (abs(omega_tgp.real - omega_gr_rs.real)
                            / max(abs(omega_gr_rs.real), 1e-10))
                delta_im = (abs(omega_tgp.imag - omega_gr_rs.imag)
                            / max(abs(omega_gr_rs.imag), 1e-10))
                table.append({
                    'M_Msun': M, 'gamma_hat': gamma_hat, 'source': C,
                    'n': n_mode, 'omega_tgp': omega_tgp,
                    'omega_gr_rs': omega_gr_rs, 'converged': converged,
                    'delta_re': delta_re, 'delta_im': delta_im,
                    'f_max': tgp['f_max'],
                })
    return table


# ===================================================================
# 8. PLOTTING
# ===================================================================

def create_plots(all_results, deviation_table, save_dir):
    """
    Generate 5-panel QNM analysis figure:
    (a) V_eff(r*): TGP vs GR Regge-Wheeler
    (b) QNM spectrum in complex omega plane
    (c) delta_f/f_GR vs BH mass (analytic)
    (d) TGP potential structure detail
    (e) Predictions summary text

    References: dodatekC_ringdown.tex, Remark rem:vs-RW
    """
    fig = plt.figure(figsize=(20, 20))
    fig.suptitle("TGP Quasi-Normal Mode Spectrum vs GR (O17)\n"
                 "Ref: dodatekC_ringdown.tex, Theorem thm:ringdown",
                 fontsize=16, y=0.995)

    colors_tgp = ['#2171b5', '#cb181d', '#31a354', '#e6550d']
    colors_gh = ['#74a9cf', '#2171b5', '#cb181d', '#e6550d']

    # ---------------------------------------------------------------
    # (a) V_eff(r*): TGP vs GR Regge-Wheeler
    # ---------------------------------------------------------------
    ax_a = fig.add_subplot(3, 2, 1)

    # GR Regge-Wheeler for l=2
    r_gr = np.linspace(1.02, 25, 2000)
    rs_gr = tortoise_schwarzschild(r_gr)
    for spin, sty, lbl in [(0, '--', 'GR scalar (s=0)'), (2, ':', 'GR tensor (s=2)')]:
        V_rw = regge_wheeler_potential(r_gr, ell=2, spin=spin)
        ax_a.plot(rs_gr, V_rw, 'k', ls=sty, lw=1.8, alpha=0.7, label=lbl)

    # TGP potentials
    ci = 0
    for res in all_results[:1]:
        for (gamma_hat, C), tgp in res['tgp'].items():
            if C != 2.0:
                continue
            rs = tgp['rs_arr']
            V = tgp['V_spline'](rs)
            mask = np.isfinite(V) & (V < 2.0)  # clip extreme values for readability
            lbl = f"TGP $\\hat{{\\gamma}}$={gamma_hat}"
            ax_a.plot(rs[mask], V[mask], color=colors_gh[ci % len(colors_gh)],
                      lw=2.0, label=lbl)
            ci += 1

    ax_a.axhline(0, color='gray', ls='-', lw=0.5, alpha=0.3)
    ax_a.set_xlabel(r'$r_* / r_s$', fontsize=12)
    ax_a.set_ylabel(r'$V_{\mathrm{eff}} \cdot r_s^2$', fontsize=12)
    ax_a.set_title('(a) Effective potential: TGP vs GR (l=2)', fontsize=13)
    ax_a.legend(fontsize=8, loc='upper right')
    ax_a.grid(True, ls=':', alpha=0.4)
    ax_a.set_xlim(-5, 40)
    ax_a.set_ylim(-0.05, 1.2)
    ax_a.annotate("TGP: rising barrier\n(no horizon, $V \\to \\infty$)",
                  xy=(1, 0.9), fontsize=8, style='italic', color='#2171b5')
    ax_a.annotate("GR: $V \\to 0$ at horizon\nand at $r \\to \\infty$",
                  xy=(15, 0.01), fontsize=8, style='italic', color='black')

    # ---------------------------------------------------------------
    # (b) QNM spectrum in complex omega plane
    # ---------------------------------------------------------------
    ax_b = fig.add_subplot(3, 2, 2)

    # GR reference points (both tensor and scalar, l=2)
    for n in range(2):
        gr = qnm_gr_schwarzschild(1.0, ell=2, n=n)
        w_t = gr['omega_tensor_dimless']
        w_s = gr['omega_scalar_dimless']
        mk = 'o' if n == 0 else 's'
        ax_b.scatter(w_t.real, w_t.imag, marker=mk, c='black', s=140,
                     edgecolors='black', linewidths=1.2, zorder=10,
                     label=f'GR tensor n={n}' if n == 0 else f'GR tensor n={n}')
        ax_b.scatter(w_s.real, w_s.imag, marker=mk, c='gray', s=100,
                     edgecolors='black', linewidths=0.8, zorder=9,
                     label=f'GR scalar n={n}')

    # TGP QNMs (dimensionless, converted to GM/c^3 units for comparison)
    ci = 0
    for res in all_results:
        for (gamma_hat, C), tgp in res['tgp'].items():
            if C != 2.0:
                continue
            for n_mode, omega, conv in tgp['modes']:
                omega_dl = omega / 2.0  # r_s units -> GM/c^3 units
                mk = 'D' if n_mode == 0 else '^'
                lbl = (f"TGP $\\hat{{\\gamma}}$={gamma_hat}, "
                       f"M={res['M_Msun']:.0f}$M_\\odot$, n={n_mode}")
                fc = colors_gh[ci % len(colors_gh)]
                ec = 'k' if conv else 'red'
                ax_b.scatter(omega_dl.real, omega_dl.imag, marker=mk,
                             c=fc, s=90, edgecolors=ec, linewidths=0.8,
                             zorder=8, label=lbl, alpha=0.8 if conv else 0.5)
            ci += 1

    ax_b.set_xlabel(r'$\omega_R \cdot GM/c^3$', fontsize=12)
    ax_b.set_ylabel(r'$\omega_I \cdot GM/c^3$', fontsize=12)
    ax_b.set_title(r'(b) QNM spectrum: complex $\omega$ plane (l=2)', fontsize=13)
    ax_b.legend(fontsize=6.5, loc='lower right', ncol=2)
    ax_b.grid(True, ls=':', alpha=0.4)

    # ---------------------------------------------------------------
    # (c) Deviation delta_f/f_GR vs BH mass (analytic)
    # ---------------------------------------------------------------
    ax_c = fig.add_subplot(3, 2, 3)

    masses = np.logspace(-0.3, 2.5, 200)
    for gi, gamma_val in enumerate([GAMMA_PHYS, 1e-40, 1e-30, 1e-20]):
        devs = []
        for M in masses:
            r_s = schwarzschild_radius(M)
            g_hat = gamma_val * r_s**2
            devs.append(g_hat / 0.3)
        ax_c.loglog(masses, devs, color=colors_gh[gi], lw=2.0,
                    label=f'$\\gamma = {gamma_val:.0e}$ m$^{{-2}}$')

    # Detector thresholds
    ax_c.axhline(0.01, color='purple', ls='--', lw=1.0, alpha=0.5)
    ax_c.annotate('LIGO O4 ~1%', xy=(200, 0.012), fontsize=9, color='purple', alpha=0.7)
    ax_c.axhline(1e-4, color='green', ls='--', lw=1.0, alpha=0.5)
    ax_c.annotate('LISA ~0.01%', xy=(200, 1.2e-4), fontsize=9, color='green', alpha=0.7)
    ax_c.axhline(1e-6, color='teal', ls='--', lw=1.0, alpha=0.5)
    ax_c.annotate('3G detectors ~ppm', xy=(200, 1.2e-6), fontsize=9, color='teal', alpha=0.7)

    # Numerical points
    for entry in deviation_table:
        if entry['converged'] and entry['n'] == 0:
            ax_c.scatter(entry['M_Msun'], max(entry['delta_re'], 1e-100),
                         marker='o', c='red', s=50, edgecolors='k',
                         linewidths=0.5, zorder=10)

    ax_c.set_xlabel(r'$M / M_\odot$', fontsize=12)
    ax_c.set_ylabel(r'$|\Delta\omega_R| / \omega_{R,\mathrm{GR}}$', fontsize=12)
    ax_c.set_title(r'(c) QNM deviation $\Delta f / f_{\mathrm{GR}}$ vs mass', fontsize=13)
    ax_c.legend(fontsize=9, loc='upper left')
    ax_c.grid(True, ls=':', alpha=0.4, which='both')
    ax_c.set_xlim(0.5, 300)
    ax_c.set_ylim(1e-100, 1e2)

    # ---------------------------------------------------------------
    # (d) TGP potential structure (multiple gamma_hat)
    # ---------------------------------------------------------------
    ax_d = fig.add_subplot(3, 2, 4)

    # Show TGP potential for one BH mass, multiple gamma_hat
    if all_results:
        res0 = all_results[0]
        ci = 0
        for (gamma_hat, C), tgp in res0['tgp'].items():
            if C != 2.0:
                continue
            rs = tgp['rs_arr']
            V = tgp['V_spline'](rs)
            mask = np.isfinite(V)
            # Normalize to show shape
            V_norm = V[mask]
            rs_norm = rs[mask]
            # Clip extreme values
            clip = V_norm < 5.0
            ax_d.semilogy(rs_norm[clip], np.maximum(V_norm[clip], 1e-6),
                          color=colors_gh[ci % len(colors_gh)], lw=2.0,
                          label=f'$\\hat{{\\gamma}}={gamma_hat}$, C={C}')
            # Mass gap floor
            ax_d.axhline(gamma_hat, color=colors_gh[ci % len(colors_gh)],
                         ls=':', lw=0.8, alpha=0.5)
            ci += 1

    ax_d.set_xlabel(r'$r_* / r_s$', fontsize=12)
    ax_d.set_ylabel(r'$V_{\mathrm{eff}}(r_*)$ [log scale]', fontsize=12)
    ax_d.set_title(f'(d) TGP potential structure (M={res0["M_Msun"]:.0f}$M_\\odot$, l=2)',
                   fontsize=13)
    ax_d.legend(fontsize=8, loc='upper right')
    ax_d.grid(True, ls=':', alpha=0.4, which='both')
    ax_d.annotate("Rising barrier\n(frozen zone, $c \\to 0$)",
                  xy=(2, 1.0), fontsize=9, style='italic', color='#cb181d')
    ax_d.annotate("Mass gap floor\n$V \\to c_0^2\\gamma$",
                  xy=(40, 0.002), fontsize=9, style='italic', color='#2171b5')

    # ---------------------------------------------------------------
    # (e) GR frequency table
    # ---------------------------------------------------------------
    ax_e = fig.add_subplot(3, 2, 5)
    ax_e.axis('off')

    lines = []
    lines.append("GR SCHWARZSCHILD QNM FREQUENCIES (l=2)")
    lines.append("=" * 55)
    lines.append(f"{'M/Msun':>8}  {'n':>3}  {'f_tensor [Hz]':>14}  "
                 f"{'f_scalar [Hz]':>14}  {'tau [ms]':>10}")
    lines.append("-" * 55)
    for M in [1.0, 5.0, 10.0, 30.0, 100.0]:
        for n in range(2):
            gr = qnm_gr_schwarzschild(M, ell=2, n=n)
            tau = 1000 / abs(gr['omega_tensor'].imag)
            lines.append(f"{M:>8.1f}  {n:>3}  {gr['f_tensor_Hz']:>14.1f}  "
                         f"{gr['f_scalar_Hz']:>14.1f}  {tau:>10.4f}")
    lines.append("")
    lines.append("TGP ANALYTIC DEVIATION ESTIMATES")
    lines.append("=" * 55)
    lines.append(f"{'M/Msun':>8}  {'gamma_hat':>14}  {'delta_w/w':>14}  "
                 f"{'chi(1%)':>12}")
    lines.append("-" * 55)
    for M in BH_MASSES_MSUN:
        est = analytic_deviation(M)
        lines.append(f"{M:>8.1f}  {est['gamma_hat']:>14.2e}  "
                     f"{est['delta_omega_frac']:>14.2e}  "
                     f"{est['chi_max_1pct']:>12.2e}")
    lines.append("")
    lines.append("Mass gap: f_min = c0*sqrt(gamma)/(2pi)")
    lines.append(f"        = {C0 * np.sqrt(GAMMA_PHYS) / (2*np.pi):.4e} Hz")

    text = "\n".join(lines)
    ax_e.text(0.02, 0.98, text, transform=ax_e.transAxes, fontsize=8.5,
              fontfamily='monospace', va='top',
              bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    # ---------------------------------------------------------------
    # (f) Key predictions summary
    # ---------------------------------------------------------------
    ax_f = fig.add_subplot(3, 2, 6)
    ax_f.axis('off')

    pred = []
    pred.append("TGP QNM KEY PREDICTIONS (O17)")
    pred.append("=" * 55)
    pred.append("")
    pred.append("1. MASS GAP (Remark rem:vs-RW, point 1):")
    pred.append(f"   omega_min = c0*sqrt(gamma) = {C0*np.sqrt(GAMMA_PHYS):.2e} rad/s")
    pred.append(f"   f_min = {C0*np.sqrt(GAMMA_PHYS)/(2*np.pi):.2e} Hz")
    pred.append("   No propagating modes below f_min.")
    pred.append("   GR: NO such gap (V_RW -> 0).")
    pred.append("")
    pred.append("2. RISING BARRIER (Cor. cor:strong-ringdown):")
    pred.append("   V_eff ~ 3*gamma*c0^2*chi  for chi >> 1")
    pred.append("   Potential GROWS toward compact object.")
    pred.append("   GR: V_RW -> 0 at horizon.")
    pred.append("")
    pred.append("3. BREATHING MODE (Remark rem:vs-RW, point 4):")
    pred.append("   TGP ringdown = scalar breathing mode.")
    pred.append("   l=0 allowed (GR tensor needs l>=2).")
    pred.append("")
    pred.append("4. PHYSICAL DEVIATION:")
    pred.append("   gamma ~ 10^{-52} m^{-2} =>")
    pred.append("   delta_f/f ~ 10^{-45} (stellar BH)")
    pred.append("   => UNDETECTABLE with current GW detectors.")
    pred.append("   Near-horizon chi amplification needed")
    pred.append("   for observable signatures.")
    pred.append("")
    pred.append("5. DETECTABILITY THRESHOLDS:")
    pred.append("   LIGO O4:    delta_f/f ~ 1%")
    pred.append("   LISA:       delta_f/f ~ 0.01%")
    pred.append("   3G (ET/CE): delta_f/f ~ 1 ppm")
    pred.append("   => Need gamma_eff > 10^{-8} m^{-2}")
    pred.append("      (via nonlinear chi_max amplification)")

    text2 = "\n".join(pred)
    ax_f.text(0.02, 0.98, text2, transform=ax_f.transAxes, fontsize=8.5,
              fontfamily='monospace', va='top',
              bbox=dict(boxstyle='round', facecolor='#e8f4fd', alpha=0.9))

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    path = os.path.join(save_dir, "qnm_spectrum.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved plot: {path}")
    plt.close(fig)
    return path


# ===================================================================
# 9. MAIN
# ===================================================================

def main():
    print("=" * 72)
    print("  TGP QUASI-NORMAL MODE SPECTRUM (O17)")
    print("  Theory of Generated Space -- Ringdown Analysis")
    print("  Reference: dodatekC_ringdown.tex, sek08")
    print("=" * 72)

    print(f"\n  PHYSICAL PARAMETERS:")
    print(f"    Phi0          = {PHI0}")
    print(f"    c0            = {C0:.6e} m/s")
    print(f"    gamma         = {GAMMA_PHYS:.2e} m^{{-2}}")
    print(f"    beta = gamma  (vacuum condition)")
    print(f"    alpha         = {ALPHA}")
    print(f"    q = 8piG/c^2  = {Q_COUPLING:.6e} m/kg")
    print(f"    Mass gap: omega_min = c0*sqrt(gamma) = {C0 * np.sqrt(GAMMA_PHYS):.4e} rad/s")
    print(f"              f_min     = {C0 * np.sqrt(GAMMA_PHYS) / (2*np.pi):.4e} Hz")

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # STEP 1: GR reference
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  STEP 1: GR Schwarzschild QNM reference frequencies (l=2)")
    print(f"{'='*72}")

    print(f"\n  {'M/M_sun':>8}  {'r_s [m]':>12}  {'Mode':>6}  {'f_tensor [Hz]':>14}  "
          f"{'f_scalar [Hz]':>14}  {'tau_tensor [ms]':>16}")
    print(f"  " + "-" * 80)

    for M in [1.0, 5.0, 10.0, 30.0, 100.0]:
        for n in range(2):
            gr = qnm_gr_schwarzschild(M, ell=2, n=n)
            tau = 1000 / abs(gr['omega_tensor'].imag)
            print(f"  {M:>8.1f}  {gr['r_s']:>12.1f}  (2,{n})  "
                  f"{gr['f_tensor_Hz']:>14.1f}  {gr['f_scalar_Hz']:>14.1f}  "
                  f"{tau:>16.4f}")

    # ------------------------------------------------------------------
    # STEP 2: TGP backgrounds, potentials, QNMs
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  STEP 2: TGP background profiles and QNM computation")
    print(f"{'='*72}")

    all_results = []
    for cfg in BH_CONFIGS:
        M = cfg['M_Msun']
        label = cfg['label']
        print(f"\n  --- {label} (M = {M} M_sun) ---")
        r_s = schwarzschild_radius(M)
        gamma_hat_phys = dimensionless_gamma(M)
        print(f"    r_s = {r_s:.2f} m")
        print(f"    gamma_hat (physical) = {gamma_hat_phys:.2e}")

        res = compute_tgp_qnm_for_bh(M, ell=ELL)
        all_results.append(res)

        for (gamma_hat, C), tgp in res['tgp'].items():
            print(f"\n    gamma_hat={gamma_hat}, C={C}:")
            print(f"      f_max = {tgp['f_max']:.4f}")
            for n_mode, omega, conv in tgp['modes']:
                Q = -omega.real / (2 * omega.imag) if omega.imag != 0 else float('inf')
                tag = "CONVERGED" if conv else "WKB est."
                print(f"      n={n_mode}: omega*r_s = {omega.real:.6f} "
                      f"{omega.imag:+.6f}i  (Q={Q:.2f}) [{tag}]")

    # ------------------------------------------------------------------
    # STEP 3: Deviations from GR
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  STEP 3: QNM deviations from GR")
    print(f"{'='*72}")

    deviation_table = compute_deviations(all_results, ell=ELL)

    print(f"\n  {'M':>5}  {'gamma_hat':>10}  {'C':>4}  {'n':>2}  "
          f"{'d_wR/wR':>12}  {'d_wI/wI':>12}  {'f_max':>7}  {'Conv':>5}")
    print(f"  " + "-" * 65)

    for e in deviation_table:
        tag = "Y" if e['converged'] else "N"
        print(f"  {e['M_Msun']:>5.0f}  {e['gamma_hat']:>10.3f}  "
              f"{e['source']:>4.1f}  {e['n']:>2}  "
              f"{e['delta_re']:>12.4e}  {e['delta_im']:>12.4e}  "
              f"{e['f_max']:>7.4f}  {tag:>5}")

    # ------------------------------------------------------------------
    # STEP 4: Analytic estimates
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  STEP 4: Analytic deviation estimates (physical gamma)")
    print(f"{'='*72}")

    print(f"\n  {'M/M_sun':>8}  {'r_s [m]':>12}  {'gamma_hat':>14}  "
          f"{'delta_w/w':>14}  {'chi_max(1%)':>14}")
    print(f"  " + "-" * 70)
    for M in BH_MASSES_MSUN:
        est = analytic_deviation(M)
        print(f"  {M:>8.1f}  {est['r_s_m']:>12.1f}  {est['gamma_hat']:>14.2e}  "
              f"{est['delta_omega_frac']:>14.2e}  {est['chi_max_1pct']:>14.2e}")

    # ------------------------------------------------------------------
    # STEP 5: Plots
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  STEP 5: Generating plots")
    print(f"{'='*72}")

    plot_path = create_plots(all_results, deviation_table, save_dir)

    # ------------------------------------------------------------------
    # STEP 6: Predictions table
    # ------------------------------------------------------------------
    print(f"\n{'='*72}")
    print("  PREDICTIONS TABLE: TGP vs GR QNM FREQUENCIES")
    print(f"{'='*72}")

    print(f"\n  +----------+-------+----------------+----------------+----------------+")
    print(f"  | M/M_sun  | Mode  | f_tensor [Hz]  | f_scalar [Hz]  | delta_f/f(TGP) |")
    print(f"  +----------+-------+----------------+----------------+----------------+")
    for M in [1.0, 10.0, 30.0, 100.0]:
        for n in range(2):
            gr = qnm_gr_schwarzschild(M, ell=2, n=n)
            est = analytic_deviation(M, ell=2, n=n)
            print(f"  | {M:>7.1f}  | (2,{n}) | {gr['f_tensor_Hz']:>13.1f}  | "
                  f"{gr['f_scalar_Hz']:>13.1f}  | {est['delta_omega_frac']:>13.2e}  |")
    print(f"  +----------+-------+----------------+----------------+----------------+")

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------
    print(f"""
{'='*72}
  SUMMARY: TGP QNM SPECTRUM (O17) -- KEY RESULTS
{'='*72}

  1. MASS GAP (Remark rem:vs-RW, Corollary cor:weak-ringdown):
     omega_min = c0 * sqrt(gamma) = {C0 * np.sqrt(GAMMA_PHYS):.4e} rad/s
     f_min = {C0 * np.sqrt(GAMMA_PHYS) / (2*np.pi):.4e} Hz
     No propagating scalar modes below f_min.
     GR: no mass gap (V_RW -> 0 at spatial infinity).

  2. RISING BARRIER (Corollary cor:strong-ringdown):
     V_eff ~ 3 gamma c0^2 chi  for chi = f/Phi0 >> 1
     Potential GROWS into compact object -- traps perturbations.
     Opposite to GR (V_RW -> 0 at horizon).

  3. BREATHING MODE (Remark rem:vs-RW, point 4):
     TGP ringdown eq. (eq:ringdown) is scalar -- breathing polarization.
     l=0 modes allowed (GR tensor requires l >= 2).

  4. PHYSICAL DEVIATION:
     gamma ~ 10^{{-52}} m^{{-2}} => gamma_hat ~ 10^{{-45}} (stellar BH)
     delta_f/f_GR ~ 10^{{-45}} : FAR below detector sensitivity.

  5. SPECIFIC GR FREQUENCIES (l=2, n=0, tensor):""")

    for M in [1.0, 10.0, 30.0]:
        gr = qnm_gr_schwarzschild(M, ell=2, n=0)
        tau = 1000 / abs(gr['omega_tensor'].imag)
        print(f"     M = {M:>5.1f} M_sun:  f = {gr['f_tensor_Hz']:>10.1f} Hz,  "
              f"tau = {tau:>8.4f} ms")

    print(f"""
  6. DETECTABILITY:
     For observable deviations (>1%), need near-horizon chi_max ~ 10^21
     through nonlinear TGP field amplification.

  Plot saved to: {plot_path}

{'='*72}
  DONE
{'='*72}
""")


if __name__ == "__main__":
    main()
