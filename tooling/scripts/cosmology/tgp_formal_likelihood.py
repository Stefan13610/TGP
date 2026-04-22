"""
tgp_formal_likelihood.py -- Formal cosmological likelihood: TGP vs LCDM
=========================================================================
H6 deliverable: chi-squared comparison with public data.

Data used:
  1. DESI DR1 BAO  (arXiv:2404.03002, Table 1)
     - D_M/r_d and D_H/r_d at 5 effective redshifts + D_V/r_d at z=0.30
     - Published correlation coefficients between D_M and D_H at each z
  2. Planck 2018 compressed CMB distance priors  (arXiv:1807.06209, Table 2)
     - Shift parameter R, acoustic scale l_A, baryon density omega_b
  3. Cosmic chronometer H(z) compilation  (Moresco+ 2022)

Models:
  - TGP:   1 free parameter  (Phi0), with beta=gamma=Phi0*H0^2/c0^2
  - LCDM:  2 free parameters (Omega_m, H0)

Outputs:
  - chi2 comparison table  (TGP vs LCDM)
  - Delta AIC / Delta BIC
  - Best-fit parameter posteriors
  - H(z), D_M(z), D_H(z) comparison plots

Usage:
  python tgp_formal_likelihood.py              # grid scan + best-fit
  python tgp_formal_likelihood.py --mcmc       # full MCMC
  python tgp_formal_likelihood.py --quick      # fast grid only (for testing)
"""

import os
import sys
import warnings
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import interp1d

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore', category=RuntimeWarning)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# Physical constants (SI)
# ============================================================================
c_SI      = 2.99792458e8       # m/s
G_SI      = 6.67430e-11        # m^3 kg^-1 s^-2
Mpc_m     = 3.08567758e22      # m per Mpc
km_per_Mpc = Mpc_m / 1e3      # conversion factor: H [s^-1] -> [km/s/Mpc]

# Planck 2018 fiducial values (Table 2, arXiv:1807.06209)
H0_PLANCK   = 67.36            # km/s/Mpc
OMm_PLANCK  = 0.3153
OMb_PLANCK  = 0.0493
OMr_PLANCK  = 9.14e-5          # radiation (photons + 3 massless nu)
Neff        = 3.044
T_CMB       = 2.7255           # K
z_star      = 1089.92          # decoupling redshift (Planck 2018)
r_d_PLANCK  = 147.09           # Mpc (sound horizon at drag, Planck 2018)

# ============================================================================
# DESI DR1 BAO data  (arXiv:2404.03002, Table 1)
# ============================================================================
# Each entry: (z_eff, observable_type, value, sigma)
# For z=0.30 (BGS): only D_V/r_d available
# For others: D_M/r_d and D_H/r_d with correlation rho

DESI_BAO = {
    # z_eff: {'DV': (val, sig)} or {'DM': (val, sig), 'DH': (val, sig), 'rho': corr}
    0.295: {'DV': (7.93, 0.15)},
    0.510: {'DM': (13.62, 0.25), 'DH': (20.98, 0.61), 'rho': -0.44},
    0.706: {'DM': (16.85, 0.32), 'DH': (20.08, 0.61), 'rho': -0.42},
    0.930: {'DM': (21.71, 0.28), 'DH': (17.88, 0.35), 'rho': -0.39},
    1.317: {'DM': (27.79, 0.69), 'DH': (13.82, 0.42), 'rho': -0.45},
    2.330: {'DM': (39.71, 0.94), 'DH': (8.52, 0.17), 'rho': -0.48},
}

# ============================================================================
# Planck 2018 compressed CMB distance priors
# (Planck Collaboration VI 2018, Table 2; Chen+ 2018 arXiv:1807.06209)
# ============================================================================
# R = sqrt(Omega_m * H0^2) * D_M(z*) / c
# l_A = pi * D_M(z*) / r_s(z*)
PLANCK_CMB = {
    'R':     (1.7502, 0.0046),
    'l_A':   (301.471, 0.090),
    'omega_b': (0.02236, 0.00015),
    # Correlation matrix (R, l_A, omega_b):
    'corr': np.array([
        [ 1.0000,  0.4600, -0.6600],
        [ 0.4600,  1.0000, -0.3300],
        [-0.6600, -0.3300,  1.0000],
    ]),
}

# Build Planck CMB covariance matrix
_planck_sig = np.array([PLANCK_CMB['R'][1], PLANCK_CMB['l_A'][1], PLANCK_CMB['omega_b'][1]])
PLANCK_COV = np.outer(_planck_sig, _planck_sig) * PLANCK_CMB['corr']
PLANCK_COV_INV = np.linalg.inv(PLANCK_COV)

# ============================================================================
# Cosmic chronometer H(z) data  (Moresco+ 2022 compilation)
# ============================================================================
DATA_HZ = np.array([
    [0.07,69.0,19.6],[0.10,69.0,12.0],[0.12,68.6,26.2],[0.17,83.0,8.0],
    [0.20,72.9,29.6],[0.27,77.0,14.0],[0.28,88.8,36.6],[0.35,82.7,8.4],
    [0.38,83.8,3.7],[0.40,95.0,17.0],[0.42,87.1,11.2],[0.44,82.6,7.8],
    [0.48,87.9,2.6],[0.51,90.4,1.9],[0.57,96.8,3.4],[0.59,98.5,3.7],
    [0.60,87.9,6.1],[0.61,97.3,2.1],[0.68,92.4,12.4],[0.73,97.3,7.0],
    [0.78,105.0,12.0],[0.80,113.1,15.1],[0.875,125.0,17.0],[0.88,90.0,40.0],
    [0.90,117.0,23.0],[1.04,154.0,20.0],[1.30,168.0,17.0],[1.36,160.0,33.6],
    [1.43,177.0,18.0],[1.53,140.0,14.0],[1.75,202.0,40.0],[2.33,224.0,8.0],
    [2.34,222.0,7.0],
])


# ############################################################################
#                          LCDM MODEL
# ############################################################################

def _H_lcdm_si(a, Om, H0_kms):
    """LCDM H(a) in s^-1."""
    H0_si = H0_kms * 1e3 / Mpc_m
    ODE = 1.0 - Om - OMr_PLANCK
    return H0_si * np.sqrt(OMr_PLANCK / a**4 + Om / a**3 + ODE)


def _H_lcdm_kms(z, Om, H0_kms):
    """LCDM H(z) in km/s/Mpc."""
    a = 1.0 / (1.0 + z)
    ODE = 1.0 - Om - OMr_PLANCK
    return H0_kms * np.sqrt(OMr_PLANCK * (1+z)**4 + Om * (1+z)**3 + ODE)


def _comoving_dist_lcdm(z, Om, H0_kms):
    """Comoving distance D_M(z) in Mpc for flat LCDM. Uses fast Gauss quadrature."""
    if z <= 0:
        return 0.0
    # 64-point Gauss-Legendre on [0, z]
    nodes, weights = np.polynomial.legendre.leggauss(64)
    zp = 0.5 * z * (nodes + 1.0)
    w = 0.5 * z * weights
    H_inv = 1.0 / _H_lcdm_kms(zp, Om, H0_kms)
    return (c_SI / 1e3) * np.sum(w * H_inv)


def _sound_horizon_lcdm(Om, H0_kms, Ob=None):
    """Sound horizon r_d at drag epoch.
    Uses improved fitting formula (Aubourg+ 2015, arXiv:1411.1074)
    calibrated against Planck 2018.  Accuracy: ~0.3% for standard cosmologies.
    Returns r_d in Mpc."""
    if Ob is None:
        # Keep omega_b = Ob*h^2 fixed at Planck value 0.02236
        h = H0_kms / 100.0
        Ob = 0.02236 / h**2
    h = H0_kms / 100.0
    omega_m = Om * h**2
    omega_b = Ob * h**2
    # Aubourg+ 2015 fitting formula (calibrated to Planck)
    r_d = 147.49 * (omega_b / 0.02230)**(-0.130) * (omega_m / 0.1426)**(-0.255)
    return r_d


def lcdm_predictions(z_arr, Om, H0_kms):
    """Compute LCDM predictions: H(z), D_M(z), D_H(z), D_V(z), all in Mpc or km/s/Mpc."""
    z_arr = np.asarray(z_arr, dtype=float)
    H_kms = _H_lcdm_kms(z_arr, Om, H0_kms)
    D_M = np.array([_comoving_dist_lcdm(z, Om, H0_kms) for z in z_arr])
    D_H = (c_SI / 1e3) / H_kms
    D_V = np.where(z_arr > 0, (D_M**2 * (c_SI/1e3) * z_arr / H_kms)**(1.0/3.0), 0.0)
    return H_kms, D_M, D_H, D_V


def _sound_horizon_star(Om, H0_kms):
    """Sound horizon at recombination z* (NOT at drag z_d).
    For CMB acoustic scale l_A = pi*D_M(z*)/r_s(z*).
    Uses scaling: r_s(z*)/r_d ~ 0.9819 (Planck 2018: 144.43/147.09)."""
    r_d = _sound_horizon_lcdm(Om, H0_kms)
    return r_d * 0.9819


def lcdm_cmb_observables(Om, H0_kms):
    """Compute CMB shift parameter R and acoustic scale l_A for LCDM."""
    D_M_star = _comoving_dist_lcdm(z_star, Om, H0_kms)
    R = np.sqrt(Om) * H0_kms / (c_SI / 1e3) * D_M_star
    r_s = _sound_horizon_star(Om, H0_kms)
    l_A = np.pi * D_M_star / r_s
    omega_b = 0.02236  # fixed
    return R, l_A, omega_b


# ############################################################################
#                          TGP MODEL
# ############################################################################
#
# Dimensionless Friedmann equation for TGP:
#
#   E^2(a) = H^2/H0^2 = Omega_m/a^3 + Omega_r/a^4 + Omega_psi(a)
#
# where Omega_psi(a) = Omega_DE0 * rho_psi(a)/rho_psi(1).
#
# The field psi satisfies the TGP Klein-Gordon equation:
#   psi'' + (3 + E'/E) psi' + (1/E^2) dV_hat/dpsi = 0
# with ' = d/d(ln a), and V_hat = V/V0 is the normalized potential.
#
# TGP potential shape (for beta = gamma):
#   f(psi) = 4*psi^3 - 3*psi^4     (normalized: f(1) = 1)
#   f'(psi) = 12*psi^2 - 12*psi^3 = 12*psi^2*(1 - psi)
#
# The field mass parameter mu^2 = V''(1)/H0^2 controls the field
# oscillation timescale.  For TGP: V''(1) = V0*(24 - 36) = -12*V0.
# So mu^2 = -12 * V0 / H0^2.
#
# Free parameters: (Phi0, Omega_m, H0).
#   - Phi0 controls the dimensionless mass: mu^2 = Phi0 * (some factor)
#   - In practice, mu ~ O(H0) means the field evolves on Hubble timescale
#
# IMPORTANT: Omega_DE0 = 1 - Omega_m - Omega_r (closure condition).
# TGP does NOT predict Omega_DE; it predicts w_DE(z) via the field dynamics.
# ############################################################################

def _f_pot(psi):
    """Normalized TGP potential shape: f(psi) = 4*psi^3 - 3*psi^4.
    f(1) = 1. This is U_eff(psi)/U_eff(1) for beta=gamma."""
    return 4.0 * psi**3 - 3.0 * psi**4


def _w_tgp(psi):
    """TGP driving function (dimensionless).
    From the exact TGP field equation: W_exact = c0^2 * gamma * [(7/3)*psi^2 - 2*psi^3]
    In dimensionless form: w(psi) = (7/3)*psi^2 - 2*psi^3.
    Note: w(1) = 7/3 - 2 = 1/3 != 0  (non-zero driving at vacuum).
    This is NOT dU_eff/dpsi -- it comes from the full nonlinear TGP substrate dynamics."""
    return (7.0 / 3.0) * psi**2 - 2.0 * psi**3


def _E2_tgp(a, psi, Om, ODE):
    """Dimensionless Friedmann: E^2 = H^2/H0^2.
    Dark energy provided by TGP field potential: Omega_DE(a) = ODE * f(psi(a))."""
    E2 = Om / a**3 + OMr_PLANCK / a**4 + ODE * _f_pot(psi)
    return max(E2, 1e-60)


def _ode_tgp_dimless(lna, state, Om, ODE, Phi0):
    """TGP exact field equation in dimensionless form.
    state = [psi, chi],  chi = d(psi)/d(ln a).

    From tgp_cosmo.py, the original ODE is:
       chi' = W(psi)/H^2 - 3*chi - 2*chi^2/psi - (H'/H)*chi
    where W = c0^2*gamma*[(7/3)psi^2 - 2psi^3].

    In dimensionless units (dividing by H0^2):
       W/H^2 = (c0^2*gamma/H0^2) * w(psi) / E^2 = Phi0 * w(psi) / E^2

    because gamma = Phi0*H0^2/c0^2  =>  c0^2*gamma/H0^2 = Phi0.

    So Phi0 directly controls the field evolution rate.
    """
    a = np.exp(lna)
    psi, chi = state
    psi = max(psi, 1e-10)

    E2 = _E2_tgp(a, psi, Om, ODE)

    # d(ln E^2)/d(ln a) via finite diff
    eps = 1e-5
    E2_p = _E2_tgp(a * (1 + eps), psi, Om, ODE)
    E2_m = _E2_tgp(a * (1 - eps), psi, Om, ODE)
    dlnE2 = (np.log(max(E2_p, 1e-60)) - np.log(max(E2_m, 1e-60))) / (2 * eps)
    Edot_over_E = 0.5 * dlnE2

    # TGP exact field equation (dimensionless):
    chi_prime = Phi0 * _w_tgp(psi) / E2 - 3.0 * chi - 2.0 * chi**2 / psi - Edot_over_E * chi

    return [chi, chi_prime]


def solve_tgp(Phi0, Om=OMm_PLANCK, H0_kms=H0_PLANCK, n_pts=3000):
    """Integrate TGP field equation in dimensionless form.
    Returns (a, psi, E, chi) where E = H/H0."""
    ODE = 1.0 - Om - OMr_PLANCK

    # Phi0 is the dimensionless coupling that controls field evolution rate.
    # From the exact TGP equation: W/H0^2 = Phi0 * w(psi) / E^2.
    # Phi0 >> 1: field evolves rapidly, psi departs from 1, w_DE departs from -1
    # Phi0 << 1: field frozen at psi=1, w_DE = -1 (Lambda-like limit)
    # Phi0 ~ O(1): intermediate quintessence-like behavior

    lna_span = (np.log(1.0 / (1.0 + 3000)), 0.0)
    lna_eval = np.linspace(lna_span[0], lna_span[1], n_pts)

    sol = solve_ivp(
        lambda lna, y: _ode_tgp_dimless(lna, y, Om, ODE, Phi0),
        lna_span, [1.0, 0.0],  # psi(early) = 1, chi(early) = 0
        t_eval=lna_eval, method="RK45",
        rtol=1e-10, atol=1e-12, max_step=0.01,
    )
    if sol.status != 0:
        return None

    a = np.exp(sol.t)
    psi = sol.y[0]
    chi = sol.y[1]
    E = np.array([np.sqrt(_E2_tgp(a[i], psi[i], Om, ODE))
                   for i in range(len(a))])
    return a, psi, E, chi


def tgp_w_de(a, psi, chi, mu2, ODE):
    """Compute w_DE(a) from field dynamics.
    w_DE = (KE - PE)/(KE + PE) where:
      PE = ODE * f(psi)
      KE ~ ODE * mu2 * chi^2 / (some normalization)
    For slow-roll (chi small): w_DE ~ -1.
    For exact: use the continuity equation approach."""
    w = np.full_like(a, -1.0)
    for i in range(1, len(a)):
        f_val = _f_pot(psi[i])
        if abs(f_val) > 1e-30:
            # From continuity: rho_DE'/rho_DE = -3(1+w_DE)
            # rho_DE proportional to f(psi), so:
            # d(ln f)/d(ln a) = -3(1+w_DE)
            # w_DE = -1 - (1/3) d(ln f)/d(ln a)
            f_prev = _f_pot(psi[i-1])
            if f_prev > 1e-30:
                dlna = np.log(a[i]) - np.log(a[i-1])
                if abs(dlna) > 1e-15:
                    dlnf = np.log(f_val / f_prev)
                    w[i] = -1.0 - dlnf / (3.0 * dlna)
    # Smooth
    w[0] = w[1]
    return w


def tgp_predictions(z_arr, Phi0, Om=OMm_PLANCK, H0_kms=H0_PLANCK):
    """Compute TGP: H(z) [km/s/Mpc], D_M(z) [Mpc], D_H(z) [Mpc], D_V(z) [Mpc]."""
    result = solve_tgp(Phi0, Om, H0_kms)
    if result is None:
        nan = np.full_like(z_arr, np.nan, dtype=float)
        return nan, nan, nan, nan

    a_sol, psi_sol, E_sol, chi_sol = result
    z_sol = 1.0 / a_sol - 1.0
    H_kms_sol = E_sol * H0_kms  # E = H/H0, so H = E*H0

    # Interpolate H(z)
    valid = np.isfinite(H_kms_sol) & (a_sol > 0)
    if valid.sum() < 10:
        nan = np.full_like(z_arr, np.nan, dtype=float)
        return nan, nan, nan, nan

    H_interp = interp1d(z_sol[valid], H_kms_sol[valid], kind='linear',
                        fill_value='extrapolate', bounds_error=False)
    H_kms = H_interp(z_arr)

    # D_M via cumulative trapezoidal integration on fine z grid
    z_fine = np.linspace(0, max(z_arr) * 1.01, 500)
    H_fine = H_interp(z_fine)
    integrand_fine = (c_SI / 1e3) / H_fine
    DM_cumul = np.concatenate([[0], np.cumsum(0.5 * (integrand_fine[:-1] + integrand_fine[1:]) * np.diff(z_fine))])
    DM_interp = interp1d(z_fine, DM_cumul, kind='linear', fill_value='extrapolate')
    D_M = DM_interp(z_arr)

    D_H = (c_SI / 1e3) / H_kms
    D_V = np.zeros_like(z_arr, dtype=float)
    for i in range(len(z_arr)):
        if z_arr[i] > 0 and H_kms[i] > 0 and D_M[i] > 0:
            D_V[i] = (D_M[i]**2 * (c_SI/1e3) * z_arr[i] / H_kms[i])**(1.0/3.0)
        else:
            D_V[i] = np.nan
    return H_kms, D_M, D_H, D_V


def tgp_sound_horizon(Phi0, Om=OMm_PLANCK, H0_kms=H0_PLANCK):
    """Sound horizon in TGP.
    Pre-recombination: field psi ~ 1 (frozen by Hubble friction at z >> 1).
    So r_d is the same as LCDM to O(delta_psi) ~ O(10^-5) at z_drag."""
    return _sound_horizon_lcdm(Om, H0_kms)


def tgp_cmb_observables(Phi0, Om=OMm_PLANCK, H0_kms=H0_PLANCK):
    """CMB shift parameter R and acoustic scale l_A for TGP."""
    result = solve_tgp(Phi0, Om, H0_kms, n_pts=3000)
    if result is None:
        return np.nan, np.nan, np.nan

    a_sol, psi_sol, E_sol, chi_sol = result
    z_sol = 1.0 / a_sol - 1.0
    H_kms_sol = E_sol * H0_kms

    valid = np.isfinite(H_kms_sol) & (a_sol > 0)
    H_interp = interp1d(z_sol[valid], H_kms_sol[valid], kind='linear',
                        fill_value='extrapolate', bounds_error=False)

    # D_M to z_star via 128-point Gauss-Legendre
    nodes, weights = np.polynomial.legendre.leggauss(128)
    zp = 0.5 * z_star * (nodes + 1.0)
    w = 0.5 * z_star * weights
    H_vals = np.array([float(H_interp(z)) for z in zp])
    H_vals = np.maximum(H_vals, 1e-10)
    D_M_star = (c_SI / 1e3) * np.sum(w / H_vals)

    R = np.sqrt(Om) * H0_kms / (c_SI / 1e3) * D_M_star
    r_s = _sound_horizon_star(Om, H0_kms)  # r_s at z*, same as LCDM pre-recombination
    l_A = np.pi * D_M_star / r_s
    omega_b = 0.02236  # fixed
    return R, l_A, omega_b


# ############################################################################
#                       CHI-SQUARED COMPUTATION
# ############################################################################

def chi2_bao(z_pred, DM_pred, DH_pred, DV_pred, r_d):
    """Chi-squared from DESI DR1 BAO data."""
    chi2 = 0.0
    # Build interpolators
    DM_interp = interp1d(z_pred, DM_pred, kind='linear', fill_value='extrapolate')
    DH_interp = interp1d(z_pred, DH_pred, kind='linear', fill_value='extrapolate')
    DV_interp = interp1d(z_pred, DV_pred, kind='linear', fill_value='extrapolate')

    for z_eff, data in DESI_BAO.items():
        if 'DV' in data and 'DM' not in data:
            # Single D_V measurement
            DV_mod = float(DV_interp(z_eff)) / r_d
            DV_obs, DV_sig = data['DV']
            chi2 += ((DV_mod - DV_obs) / DV_sig)**2
        elif 'DM' in data and 'DH' in data:
            # Correlated D_M/r_d and D_H/r_d
            DM_mod = float(DM_interp(z_eff)) / r_d
            DH_mod = float(DH_interp(z_eff)) / r_d
            DM_obs, DM_sig = data['DM']
            DH_obs, DH_sig = data['DH']
            rho = data['rho']
            # 2x2 covariance
            cov = np.array([
                [DM_sig**2, rho * DM_sig * DH_sig],
                [rho * DM_sig * DH_sig, DH_sig**2],
            ])
            delta = np.array([DM_mod - DM_obs, DH_mod - DH_obs])
            chi2 += float(delta @ np.linalg.solve(cov, delta))
    return chi2


def chi2_hz(z_pred, H_pred):
    """Chi-squared from CC H(z) data."""
    H_interp = interp1d(z_pred, H_pred, kind='linear', fill_value='extrapolate')
    H_mod = np.array([float(H_interp(z)) for z in DATA_HZ[:, 0]])
    return float(np.sum(((DATA_HZ[:, 1] - H_mod) / DATA_HZ[:, 2])**2))


def chi2_cmb(R_mod, lA_mod, ob_mod):
    """Chi-squared from Planck compressed CMB priors."""
    obs = np.array([PLANCK_CMB['R'][0], PLANCK_CMB['l_A'][0], PLANCK_CMB['omega_b'][0]])
    mod = np.array([R_mod, lA_mod, ob_mod])
    delta = mod - obs
    return float(delta @ PLANCK_COV_INV @ delta)


# ############################################################################
#                       TOTAL LIKELIHOOD
# ############################################################################

def total_chi2_lcdm(Om, H0_kms, verbose=False):
    """Total chi-squared for LCDM."""
    z_eval = np.linspace(0.01, 3.0, 100)
    H_kms, D_M, D_H, D_V = lcdm_predictions(z_eval, Om, H0_kms)
    r_d = _sound_horizon_lcdm(Om, H0_kms)

    c2_bao = chi2_bao(z_eval, D_M, D_H, D_V, r_d)
    c2_hz  = chi2_hz(z_eval, H_kms)
    R, l_A, ob = lcdm_cmb_observables(Om, H0_kms)
    c2_cmb = chi2_cmb(R, l_A, ob)
    total = c2_bao + c2_hz + c2_cmb

    if verbose:
        print(f"    LCDM(Om={Om:.4f}, H0={H0_kms:.2f}): "
              f"chi2_BAO={c2_bao:.2f}, chi2_Hz={c2_hz:.2f}, chi2_CMB={c2_cmb:.2f}, "
              f"TOTAL={total:.2f}  (r_d={r_d:.2f} Mpc)")
    return total, c2_bao, c2_hz, c2_cmb


def total_chi2_tgp(Phi0, Om=OMm_PLANCK, H0_kms=H0_PLANCK, verbose=False):
    """Total chi-squared for TGP."""
    z_eval = np.linspace(0.01, 3.0, 100)
    H_kms, D_M, D_H, D_V = tgp_predictions(z_eval, Phi0, Om, H0_kms)
    if np.any(np.isnan(H_kms)):
        return 1e10, 1e10, 1e10, 1e10
    r_d = tgp_sound_horizon(Phi0, Om, H0_kms)

    c2_bao = chi2_bao(z_eval, D_M, D_H, D_V, r_d)
    c2_hz  = chi2_hz(z_eval, H_kms)

    R, l_A, ob = tgp_cmb_observables(Phi0, Om, H0_kms)
    if np.isnan(R):
        return 1e10, 1e10, 1e10, 1e10
    c2_cmb = chi2_cmb(R, l_A, ob)
    total = c2_bao + c2_hz + c2_cmb

    if verbose:
        ODE = 1.0 - Om - OMr_PLANCK
        print(f"    TGP(Phi0={Phi0:.2f}, Om={Om:.4f}, H0={H0_kms:.2f}, ODE={ODE:.4f}): "
              f"chi2_BAO={c2_bao:.2f}, chi2_Hz={c2_hz:.2f}, chi2_CMB={c2_cmb:.2f}, "
              f"TOTAL={total:.2f}  (r_d={r_d:.2f} Mpc)")
    return total, c2_bao, c2_hz, c2_cmb


# ############################################################################
#                       GRID SCAN & OPTIMIZATION
# ############################################################################

def optimize_lcdm(verbose=True):
    """Find best-fit LCDM parameters (Om, H0)."""
    if verbose:
        print("\n  Optimizing LCDM ...")

    def neg_chi2(params):
        Om, H0 = params
        if not (0.1 < Om < 0.6 and 50 < H0 < 90):
            return 1e10
        try:
            total, _, _, _ = total_chi2_lcdm(Om, H0)
            return total
        except Exception:
            return 1e10

    # Multi-start Nelder-Mead (fast) from several initial points
    best_result = None
    for Om0, H0_0 in [(0.30, 67.5), (0.28, 70.0), (0.32, 68.0), (0.31, 72.0)]:
        res = minimize(neg_chi2, [Om0, H0_0], method='Nelder-Mead',
                       options={'xatol': 1e-5, 'fatol': 0.01, 'maxiter': 200})
        if best_result is None or res.fun < best_result.fun:
            best_result = res

    Om_best, H0_best = best_result.x
    chi2_best = best_result.fun

    if verbose:
        total_chi2_lcdm(Om_best, H0_best, verbose=True)
    return Om_best, H0_best, chi2_best


def optimize_tgp(verbose=True):
    """Find best-fit TGP parameter Phi0 (keeping Om, H0 at Planck values for now)."""
    if verbose:
        print("\n  Optimizing TGP (Phi0 scan) ...")

    # Coarse scan first
    Phi0_vals = np.logspace(0.3, 2.5, 40)  # 2 to ~300
    chi2_vals = np.zeros(len(Phi0_vals))
    for i, P in enumerate(Phi0_vals):
        try:
            chi2_vals[i], _, _, _ = total_chi2_tgp(P)
        except Exception:
            chi2_vals[i] = 1e10

    idx_best = np.argmin(chi2_vals)
    Phi0_init = Phi0_vals[idx_best]

    if verbose:
        print(f"    Coarse best: Phi0={Phi0_init:.2f}, chi2={chi2_vals[idx_best]:.2f}")

    # Fine optimization using bounded method (more robust)
    def neg_chi2_phi(P):
        if P < 1 or P > 500:
            return 1e10
        try:
            total, _, _, _ = total_chi2_tgp(P)
            return total
        except Exception:
            return 1e10

    result = minimize_scalar(neg_chi2_phi, bounds=(max(1, Phi0_init*0.3),
                                                    min(500, Phi0_init*3)),
                             method='bounded')
    Phi0_best = result.x
    chi2_best = result.fun

    if verbose:
        total_chi2_tgp(Phi0_best, verbose=True)
    return Phi0_best, chi2_best


def optimize_tgp_full(verbose=True):
    """Find best-fit TGP parameters (Phi0, Om, H0)."""
    if verbose:
        print("\n  Optimizing TGP (Phi0, Om, H0) ...")

    def neg_chi2(params):
        Phi0, Om, H0 = params
        if not (0.001 < Phi0 < 500 and 0.1 < Om < 0.6 and 50 < H0 < 90):
            return 1e10
        try:
            total, _, _, _ = total_chi2_tgp(Phi0, Om, H0)
            return total
        except Exception:
            return 1e10

    # Multi-start Nelder-Mead from several initial points
    best_result = None
    starts = [
        (0.001, 0.30, 67.5),  # Lambda-like (frozen field)
        (0.01, 0.31, 68.0),
        (0.1, 0.30, 70.0),
        (1.0, 0.29, 68.0),
        (5.0, 0.31, 67.5),
        (10.0, 0.30, 68.0),
        (50.0, 0.28, 70.0),
    ]
    for P0, Om0, H0_0 in starts:
        if verbose:
            print(f"    Trying Phi0={P0}, Om={Om0}, H0={H0_0} ...", end="", flush=True)
        res = minimize(neg_chi2, [P0, Om0, H0_0], method='Nelder-Mead',
                       options={'xatol': 1e-5, 'fatol': 0.01, 'maxiter': 300})
        if verbose:
            print(f" chi2={res.fun:.2f}")
        if best_result is None or res.fun < best_result.fun:
            best_result = res

    Phi0_best, Om_best, H0_best = best_result.x
    chi2_best = best_result.fun

    if verbose:
        total_chi2_tgp(Phi0_best, Om_best, H0_best, verbose=True)
    return Phi0_best, Om_best, H0_best, chi2_best


# ############################################################################
#                       INFORMATION CRITERIA
# ############################################################################

def compute_ic(chi2_lcdm, k_lcdm, chi2_tgp, k_tgp, N_data):
    """Compute AIC and BIC for both models."""
    AIC_lcdm = chi2_lcdm + 2 * k_lcdm
    AIC_tgp  = chi2_tgp  + 2 * k_tgp
    BIC_lcdm = chi2_lcdm + k_lcdm * np.log(N_data)
    BIC_tgp  = chi2_tgp  + k_tgp  * np.log(N_data)

    return {
        'AIC_lcdm': AIC_lcdm, 'AIC_tgp': AIC_tgp,
        'BIC_lcdm': BIC_lcdm, 'BIC_tgp': BIC_tgp,
        'dAIC': AIC_tgp - AIC_lcdm,
        'dBIC': BIC_tgp - BIC_lcdm,
    }


# ############################################################################
#                       PLOTTING
# ############################################################################

def plot_comparison(lcdm_params, tgp_params_best, save_dir):
    """Publication-quality comparison plots."""
    os.makedirs(save_dir, exist_ok=True)
    Om_l, H0_l = lcdm_params
    Phi0_t, Om_t, H0_t = tgp_params_best

    z_plot = np.linspace(0.01, 2.8, 300)

    H_l, DM_l, DH_l, DV_l = lcdm_predictions(z_plot, Om_l, H0_l)
    H_t, DM_t, DH_t, DV_t = tgp_predictions(z_plot, Phi0_t, Om_t, H0_t)
    rd_l = _sound_horizon_lcdm(Om_l, H0_l)
    rd_t = tgp_sound_horizon(Phi0_t, Om_t, H0_t)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # --- Panel 1: H(z) ---
    ax = axes[0, 0]
    ax.plot(z_plot, H_l, 'k-', lw=2, label=rf'$\Lambda$CDM ($\Omega_m={Om_l:.3f}, H_0={H0_l:.1f}$)')
    ax.plot(z_plot, H_t, 'C0-', lw=2, label=rf'TGP ($\Phi_0={Phi0_t:.1f}$)')
    ax.errorbar(DATA_HZ[:,0], DATA_HZ[:,1], yerr=DATA_HZ[:,2],
                fmt='ro', ms=3, elinewidth=0.8, capsize=2, alpha=0.6, label='CC data')
    ax.set_xlabel(r'$z$'); ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.set_title('Hubble parameter'); ax.legend(fontsize=8); ax.grid(True, ls=':', alpha=0.3)

    # --- Panel 2: D_M/r_d ---
    ax = axes[0, 1]
    ax.plot(z_plot, DM_l/rd_l, 'k-', lw=2, label=r'$\Lambda$CDM')
    ax.plot(z_plot, DM_t/rd_t, 'C0-', lw=2, label='TGP')
    for z_eff, data in DESI_BAO.items():
        if 'DM' in data:
            ax.errorbar(z_eff, data['DM'][0], yerr=data['DM'][1],
                        fmt='rs', ms=6, elinewidth=1.2, capsize=3, zorder=5)
    ax.set_xlabel(r'$z$'); ax.set_ylabel(r'$D_M/r_d$')
    ax.set_title('Comoving distance / sound horizon')
    ax.legend(fontsize=9); ax.grid(True, ls=':', alpha=0.3)

    # --- Panel 3: D_H/r_d ---
    ax = axes[1, 0]
    ax.plot(z_plot, DH_l/rd_l, 'k-', lw=2, label=r'$\Lambda$CDM')
    ax.plot(z_plot, DH_t/rd_t, 'C0-', lw=2, label='TGP')
    for z_eff, data in DESI_BAO.items():
        if 'DH' in data:
            ax.errorbar(z_eff, data['DH'][0], yerr=data['DH'][1],
                        fmt='rs', ms=6, elinewidth=1.2, capsize=3, zorder=5)
    ax.set_xlabel(r'$z$'); ax.set_ylabel(r'$D_H/r_d$')
    ax.set_title('Hubble distance / sound horizon')
    ax.legend(fontsize=9); ax.grid(True, ls=':', alpha=0.3)

    # --- Panel 4: D_V/r_d ---
    ax = axes[1, 1]
    ax.plot(z_plot, DV_l/rd_l, 'k-', lw=2, label=r'$\Lambda$CDM')
    ax.plot(z_plot, DV_t/rd_t, 'C0-', lw=2, label='TGP')
    for z_eff, data in DESI_BAO.items():
        if 'DV' in data:
            ax.errorbar(z_eff, data['DV'][0], yerr=data['DV'][1],
                        fmt='rs', ms=6, elinewidth=1.2, capsize=3, zorder=5)
    ax.set_xlabel(r'$z$'); ax.set_ylabel(r'$D_V/r_d$')
    ax.set_title('Volume-averaged distance / sound horizon')
    ax.legend(fontsize=9); ax.grid(True, ls=':', alpha=0.3)

    fig.suptitle('TGP vs $\\Lambda$CDM: BAO + H(z) comparison', fontsize=15, y=1.01)
    fig.tight_layout()
    p = os.path.join(save_dir, 'tgp_lcdm_bao_comparison.png')
    fig.savefig(p, dpi=180, bbox_inches='tight')
    print(f"  Saved {p}")
    plt.close(fig)


def plot_chi2_profile(save_dir):
    """Chi-squared profile as function of Phi0."""
    os.makedirs(save_dir, exist_ok=True)
    Phi0_arr = np.logspace(0.3, 2.3, 50)
    chi2_arr = np.zeros(len(Phi0_arr))
    chi2_bao_arr = np.zeros(len(Phi0_arr))
    chi2_hz_arr = np.zeros(len(Phi0_arr))
    chi2_cmb_arr = np.zeros(len(Phi0_arr))

    print("\n  Computing chi2 profile ...")
    for i, P in enumerate(Phi0_arr):
        try:
            chi2_arr[i], chi2_bao_arr[i], chi2_hz_arr[i], chi2_cmb_arr[i] = total_chi2_tgp(P)
        except Exception:
            chi2_arr[i] = np.nan

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Total chi2
    valid = np.isfinite(chi2_arr) & (chi2_arr < 1e9)
    ax1.plot(Phi0_arr[valid], chi2_arr[valid], 'C0-', lw=2)
    ax1.set_xscale('log')
    ax1.set_xlabel(r'$\Phi_0$', fontsize=13)
    ax1.set_ylabel(r'$\chi^2_{\rm total}$', fontsize=13)
    ax1.set_title(r'TGP $\chi^2$ profile (BAO + Hz + CMB)')
    ax1.grid(True, ls=':', alpha=0.3)
    idx_min = np.nanargmin(chi2_arr)
    ax1.axvline(Phi0_arr[idx_min], color='r', ls='--', lw=0.8,
                label=rf'best $\Phi_0 = {Phi0_arr[idx_min]:.1f}$')
    ax1.legend()

    # Component breakdown
    ax2.plot(Phi0_arr[valid], chi2_bao_arr[valid], 'C1-', lw=1.5, label=r'$\chi^2_{\rm BAO}$')
    ax2.plot(Phi0_arr[valid], chi2_hz_arr[valid], 'C2-', lw=1.5, label=r'$\chi^2_{H(z)}$')
    ax2.plot(Phi0_arr[valid], chi2_cmb_arr[valid], 'C3-', lw=1.5, label=r'$\chi^2_{\rm CMB}$')
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$\Phi_0$', fontsize=13)
    ax2.set_ylabel(r'$\chi^2$', fontsize=13)
    ax2.set_title(r'$\chi^2$ component breakdown')
    ax2.legend(); ax2.grid(True, ls=':', alpha=0.3)

    fig.tight_layout()
    p = os.path.join(save_dir, 'tgp_chi2_profile.png')
    fig.savefig(p, dpi=180, bbox_inches='tight')
    print(f"  Saved {p}")
    plt.close(fig)


def plot_residuals(lcdm_params, tgp_params_best, save_dir):
    """Residual plot: (model - data)/sigma for each BAO point."""
    os.makedirs(save_dir, exist_ok=True)
    Om_l, H0_l = lcdm_params
    Phi0_t, Om_t, H0_t = tgp_params_best

    z_eval = np.linspace(0.01, 3.0, 500)
    _, DM_l, DH_l, DV_l = lcdm_predictions(z_eval, Om_l, H0_l)
    _, DM_t, DH_t, DV_t = tgp_predictions(z_eval, Phi0_t, Om_t, H0_t)
    rd_l = _sound_horizon_lcdm(Om_l, H0_l)
    rd_t = tgp_sound_horizon(Phi0_t, Om_t, H0_t)

    DM_l_interp = interp1d(z_eval, DM_l/rd_l, fill_value='extrapolate')
    DH_l_interp = interp1d(z_eval, DH_l/rd_l, fill_value='extrapolate')
    DV_l_interp = interp1d(z_eval, DV_l/rd_l, fill_value='extrapolate')
    DM_t_interp = interp1d(z_eval, DM_t/rd_t, fill_value='extrapolate')
    DH_t_interp = interp1d(z_eval, DH_t/rd_t, fill_value='extrapolate')
    DV_t_interp = interp1d(z_eval, DV_t/rd_t, fill_value='extrapolate')

    labels = []
    resid_l = []
    resid_t = []

    for z_eff, data in sorted(DESI_BAO.items()):
        if 'DV' in data and 'DM' not in data:
            obs, sig = data['DV']
            rl = (float(DV_l_interp(z_eff)) - obs) / sig
            rt = (float(DV_t_interp(z_eff)) - obs) / sig
            labels.append(f'$D_V/r_d$ z={z_eff:.2f}')
            resid_l.append(rl); resid_t.append(rt)
        if 'DM' in data:
            obs, sig = data['DM']
            rl = (float(DM_l_interp(z_eff)) - obs) / sig
            rt = (float(DM_t_interp(z_eff)) - obs) / sig
            labels.append(f'$D_M/r_d$ z={z_eff:.2f}')
            resid_l.append(rl); resid_t.append(rt)
        if 'DH' in data:
            obs, sig = data['DH']
            rl = (float(DH_l_interp(z_eff)) - obs) / sig
            rt = (float(DH_t_interp(z_eff)) - obs) / sig
            labels.append(f'$D_H/r_d$ z={z_eff:.2f}')
            resid_l.append(rl); resid_t.append(rt)

    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(labels))
    ax.barh(y - 0.2, resid_l, height=0.35, color='gray', alpha=0.7, label=r'$\Lambda$CDM')
    ax.barh(y + 0.2, resid_t, height=0.35, color='C0', alpha=0.7, label='TGP')
    ax.axvline(0, color='k', lw=0.8)
    for ns in [-2, -1, 1, 2]:
        ax.axvline(ns, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r'Residual $({\rm model} - {\rm data})/\sigma$', fontsize=12)
    ax.set_title('BAO residuals: TGP vs $\\Lambda$CDM')
    ax.legend(fontsize=10); ax.grid(True, axis='x', ls=':', alpha=0.3)
    fig.tight_layout()
    p = os.path.join(save_dir, 'tgp_bao_residuals.png')
    fig.savefig(p, dpi=180, bbox_inches='tight')
    print(f"  Saved {p}")
    plt.close(fig)


# ############################################################################
#                       MCMC
# ############################################################################

def run_mcmc_formal(save_dir, nwalkers=24, nsteps=800, burn_in=300):
    """Full MCMC for TGP: fit (Phi0, Om, H0) to BAO + Hz + CMB."""
    os.makedirs(save_dir, exist_ok=True)

    try:
        import emcee
    except ImportError:
        print("  [WARN] emcee not installed. Skipping MCMC.")
        return None

    ndim = 3
    labels = [r'$\Phi_0$', r'$\Omega_m$', r'$H_0$']

    def log_prior(theta):
        Phi0, Om, H0 = theta
        if 1 < Phi0 < 500 and 0.1 < Om < 0.6 and 55 < H0 < 85:
            return 0.0
        return -np.inf

    def log_likelihood(theta):
        Phi0, Om, H0 = theta
        try:
            total, _, _, _ = total_chi2_tgp(Phi0, Om, H0)
            if not np.isfinite(total) or total > 1e9:
                return -np.inf
            return -0.5 * total
        except Exception:
            return -np.inf

    def log_posterior(theta):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        ll = log_likelihood(theta)
        return lp + ll

    # Initialize around best-fit
    p0_center = np.array([20.0, 0.31, 67.5])
    p0 = p0_center + 0.01 * p0_center * np.random.randn(nwalkers, ndim)
    p0[:, 0] = np.clip(p0[:, 0], 2, 490)
    p0[:, 1] = np.clip(p0[:, 1], 0.11, 0.59)
    p0[:, 2] = np.clip(p0[:, 2], 56, 84)

    print(f"\n  Running MCMC: {nwalkers} walkers x {nsteps} steps, 3 params ...")
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)
    sampler.run_mcmc(p0, nsteps, progress=True)

    flat = sampler.get_chain(discard=burn_in, flat=True)
    print(f"  Chain shape: {flat.shape}")

    for i, lab in enumerate(labels):
        q16, q50, q84 = np.percentile(flat[:, i], [16, 50, 84])
        print(f"  {lab}: {q50:.3f}  (+{q84-q50:.3f} / -{q50-q16:.3f})")

    # Trace plot
    chain = sampler.get_chain()
    fig, axes = plt.subplots(ndim, 1, figsize=(12, 8), sharex=True)
    for i in range(ndim):
        axes[i].plot(chain[:, :, i], alpha=0.2, lw=0.5)
        axes[i].set_ylabel(labels[i])
        axes[i].axvline(burn_in, color='r', ls='--', lw=0.8)
    axes[-1].set_xlabel('Step')
    fig.suptitle('TGP MCMC trace', fontsize=14)
    fig.tight_layout()
    p = os.path.join(save_dir, 'tgp_formal_mcmc_trace.png')
    fig.savefig(p, dpi=150); print(f"  Saved {p}"); plt.close(fig)

    # Corner plot
    try:
        import corner
        fig = corner.corner(flat, labels=labels,
                            quantiles=[0.16, 0.5, 0.84], show_titles=True)
        p = os.path.join(save_dir, 'tgp_formal_mcmc_corner.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"  Saved {p}"); plt.close(fig)
    except ImportError:
        print("  [WARN] corner not installed, skipping corner plot")

    return sampler, flat


# ############################################################################
#                          MAIN
# ############################################################################

def print_summary_table(chi2_l, k_l, chi2_t, k_t, N_data,
                        params_l, params_t):
    """Print publication-ready summary table."""
    ic = compute_ic(chi2_l, k_l, chi2_t, k_t, N_data)

    print("\n" + "=" * 72)
    print("  FORMAL LIKELIHOOD COMPARISON: TGP vs LCDM")
    print("=" * 72)
    print(f"  {'':20s}  {'LCDM':>15s}  {'TGP':>15s}")
    print("  " + "-" * 55)
    print(f"  {'Free parameters':20s}  {k_l:>15d}  {k_t:>15d}")
    print(f"  {'chi2_total':20s}  {chi2_l:>15.2f}  {chi2_t:>15.2f}")
    print(f"  {'chi2/N_data':20s}  {chi2_l/N_data:>15.3f}  {chi2_t/N_data:>15.3f}")
    print(f"  {'AIC':20s}  {ic['AIC_lcdm']:>15.2f}  {ic['AIC_tgp']:>15.2f}")
    print(f"  {'BIC':20s}  {ic['BIC_lcdm']:>15.2f}  {ic['BIC_tgp']:>15.2f}")
    print("  " + "-" * 55)
    print(f"  {'Delta chi2 (TGP-L)':20s}  {chi2_t - chi2_l:>15.2f}")
    print(f"  {'Delta AIC  (TGP-L)':20s}  {ic['dAIC']:>15.2f}")
    print(f"  {'Delta BIC  (TGP-L)':20s}  {ic['dBIC']:>15.2f}")
    print("  " + "-" * 55)

    # Interpretation
    dAIC = ic['dAIC']
    if dAIC < -10:
        interp = "Strong preference for TGP"
    elif dAIC < -6:
        interp = "Moderate preference for TGP"
    elif dAIC < -2:
        interp = "Weak preference for TGP"
    elif dAIC < 2:
        interp = "No significant difference"
    elif dAIC < 6:
        interp = "Weak preference for LCDM"
    elif dAIC < 10:
        interp = "Moderate preference for LCDM"
    else:
        interp = "Strong preference for LCDM"
    print(f"  {'Interpretation':20s}  {interp}")
    print("=" * 72)

    Om_l, H0_l = params_l
    Phi0_t, Om_t, H0_t = params_t
    ODE_t = 1.0 - Om_t - OMr_PLANCK
    print(f"\n  Best-fit LCDM:  Omega_m = {Om_l:.4f},  H0 = {H0_l:.2f} km/s/Mpc")
    print(f"  Best-fit TGP:   Phi0 = {Phi0_t:.2f} (mu^2),  Omega_m = {Om_t:.4f},  "
          f"H0 = {H0_t:.2f} km/s/Mpc")
    print(f"                  Omega_DE = {ODE_t:.4f}")
    print(f"  N_data = {N_data} (BAO: 11, Hz: {len(DATA_HZ)}, CMB: 3)")
    return ic


def main():
    quick = '--quick' in sys.argv
    do_mcmc = '--mcmc' in sys.argv

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            'plots')
    os.makedirs(save_dir, exist_ok=True)

    print("=" * 72)
    print("  TGP Formal Cosmological Likelihood Pipeline")
    print("  Data: DESI DR1 BAO + Planck 2018 CMB + CC H(z)")
    print("=" * 72)

    # Count data points:
    # DESI: 1 (DV at z=0.30) + 5*2 (DM+DH at 5 redshifts) = 11
    # CC H(z): 33
    # CMB: 3 (R, l_A, omega_b)
    N_data = 11 + len(DATA_HZ) + 3  # = 47

    # ---- Step 1: Optimize LCDM ----
    Om_l, H0_l, chi2_l_total = optimize_lcdm(verbose=True)
    chi2_l_total, c2l_bao, c2l_hz, c2l_cmb = total_chi2_lcdm(Om_l, H0_l, verbose=True)
    k_lcdm = 2  # (Omega_m, H0)

    if quick:
        # Quick TGP scan at fixed Om, H0
        Phi0_t, chi2_t_total = optimize_tgp(verbose=True)
        chi2_t_total, c2t_bao, c2t_hz, c2t_cmb = total_chi2_tgp(Phi0_t, verbose=True)
        tgp_best = (Phi0_t, OMm_PLANCK, H0_PLANCK)
        k_tgp = 1
    else:
        # Full TGP optimization (Phi0, Om, H0)
        Phi0_t, Om_t, H0_t, chi2_t_total = optimize_tgp_full(verbose=True)
        chi2_t_total, c2t_bao, c2t_hz, c2t_cmb = total_chi2_tgp(
            Phi0_t, Om_t, H0_t, verbose=True)
        tgp_best = (Phi0_t, Om_t, H0_t)
        k_tgp = 3  # (Phi0, Om, H0)

    # ---- Summary ----
    ic = print_summary_table(chi2_l_total, k_lcdm, chi2_t_total, k_tgp,
                             N_data, (Om_l, H0_l), tgp_best)

    # ---- Plots ----
    print("\n  Generating comparison plots ...")
    plot_comparison((Om_l, H0_l), tgp_best, save_dir)
    plot_residuals((Om_l, H0_l), tgp_best, save_dir)
    plot_chi2_profile(save_dir)

    # ---- MCMC ----
    if do_mcmc:
        run_mcmc_formal(save_dir)

    # ---- Print component breakdown ----
    print(f"\n  Chi2 breakdown:")
    print(f"    {'Component':15s}  {'LCDM':>10s}  {'TGP':>10s}  {'Delta':>10s}")
    print(f"    {'-'*50}")
    print(f"    {'BAO':15s}  {c2l_bao:>10.2f}  {c2t_bao:>10.2f}  {c2t_bao-c2l_bao:>+10.2f}")
    print(f"    {'H(z)':15s}  {c2l_hz:>10.2f}  {c2t_hz:>10.2f}  {c2t_hz-c2l_hz:>+10.2f}")
    print(f"    {'CMB':15s}  {c2l_cmb:>10.2f}  {c2t_cmb:>10.2f}  {c2t_cmb-c2l_cmb:>+10.2f}")
    print(f"    {'TOTAL':15s}  {chi2_l_total:>10.2f}  {chi2_t_total:>10.2f}  "
          f"{chi2_t_total-chi2_l_total:>+10.2f}")

    print("\n  Done.")
    return ic


if __name__ == "__main__":
    main()
