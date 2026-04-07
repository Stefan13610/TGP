"""
tgp_cosmo.py  --  Theory of Generated Space: consolidated cosmological script
===============================================================================
Merges w_de_exact.py (physical SI units, derived parameters) with
cosmological_evolution.py (growth rate, G_eff, MCMC).

Key design:
  - All TGP parameters derived from a single physical parameter Phi0
    via the vacuum condition beta = gamma = Phi0 H0^2 / c0^2.
  - tau0 = sqrt(Phi0) / (c0 sqrt(gamma))  (NOT free).
  - Physical SI units throughout.
  - EXACT nonlinear field equation (not linearized).
  - w_DE from field energy density (kinetic + potential), not from H^2 residual.
  - Growth rate f*sigma8 via G_eff(k,a) with m_sp^2 = gamma (for beta=gamma).
  - MCMC fits Phi0 (single physical parameter) to H(z) + f*sigma8(z) data.

Modes:
  python tgp_cosmo.py              -- scan Phi0 values, produce all plots
  python tgp_cosmo.py --mcmc       -- MCMC fit of Phi0 to observational data

Outputs (saved to scripts/plots/):
  tgp_w_de.png         -- w_DE(z) for different Phi0
  tgp_w0_wa.png        -- w0-wa plane with observational constraints
  tgp_psi.png          -- psi(z) field evolution
  tgp_growth.png       -- f*sigma8(z) and G_eff(k)
  tgp_mcmc_trace.png   -- MCMC trace plots
  tgp_mcmc_corner.png  -- MCMC posteriors
  tgp_mcmc_bestfit.png -- best-fit model vs data
"""

import os
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# Physical constants (SI)
# ============================================================================
c0    = 3.0e8            # m/s
G0    = 6.674e-11        # m^3 kg^-1 s^-2
H0_SI = 2.27e-18         # s^-1  (~ 70 km/s/Mpc)
H0_km = 70.0             # km/s/Mpc  (for data comparison)
kappa = 8.0 * np.pi * G0 / c0**4

# Cosmological density fractions today
Omega_m0 = 0.3
Omega_r0 = 9.0e-5

# ============================================================================
# TGP parameters from Phi0  (vacuum condition beta = gamma)
# ============================================================================
PHI0_VALUES = [25, 50, 80, 100, 115, 130, 150]


def tgp_params(Phi0):
    """Derive (beta, gamma, tau0) from Phi0 under beta = gamma."""
    gamma = Phi0 * H0_SI**2 / c0**2
    beta  = gamma
    tau0  = np.sqrt(Phi0) / (c0 * np.sqrt(gamma))
    return beta, gamma, tau0


# ============================================================================
# Potentials
# ============================================================================
def W_exact(psi, beta, gamma):
    """RHS of exact FRW field equation: c0^2 (gamma*psi - beta).

    Derived from correct action S[g] = int[1/2 g^4(nabla g)^2 + P(g)]
    with P(g) = (beta/7)g^7 - (gamma/8)g^8 and sqrt(-g_eff) = c0*psi.
    The FRW equation is: psi_tt + 3H psi_t + 3 psi_t^2/psi = c0^2(gamma*psi - beta).
    NOTE: kappa does NOT appear here -- it cancels because 1/kappa
    multiplies the entire gravitational sector in the action."""
    return c0**2 * (gamma * psi - beta)


def U_eff(psi, beta, gamma):
    """Action potential: P(psi) = (beta/7)psi^7 - (gamma/8)psi^8.

    Derived from P'(g) = K(g) * (gamma g^3 - beta g^2) with K = g^4,
    giving P(g) = (beta/7)g^7 - (gamma/8)g^8."""
    return (beta / 7.0) * psi**7 - (gamma / 8.0) * psi**8


# ============================================================================
# Friedmann equation  H^2(a, psi, dpsi_dt)
# ============================================================================
def rho_crit0():
    return 3.0 * H0_SI**2 / (8.0 * np.pi * G0)


def hubble_squared(a, psi, dpsi_dt, beta, gamma):
    """H^2 from Friedmann eq with matter, radiation, and field energy.
    NOTE: kappa cancels from the field energy (1/kappa on entire gravitational sector)."""
    rc0 = rho_crit0()
    rho_m = Omega_m0 * rc0 / a**3
    rho_r = Omega_r0 * rc0 / a**4
    kinetic = 0.5 * (dpsi_dt / c0)**2
    potential = U_eff(psi, beta, gamma)
    H2 = (8.0 * np.pi * G0 / 3.0) * (rho_m + rho_r + kinetic + potential)
    return max(H2, 1e-60)


def hubble_lcdm(a):
    """LCDM Hubble parameter H(a) in s^-1."""
    Omega_DE0 = 1.0 - Omega_m0 - Omega_r0
    return H0_SI * np.sqrt(Omega_r0 / a**4 + Omega_m0 / a**3 + Omega_DE0)


# ============================================================================
# ODE: exact field equation in ln(a)
# ============================================================================
def ode_exact(lna, state, beta, gamma):
    """State = [psi, chi] where chi = dpsi/d(ln a). Returns [chi, chi']."""
    a = np.exp(lna)
    psi, chi = state
    psi = max(psi, 1e-10)

    H2_est = hubble_squared(a, psi, 0.0, beta, gamma)
    H_est  = np.sqrt(H2_est)
    dpsi_dt_est = H_est * chi
    H2 = hubble_squared(a, psi, dpsi_dt_est, beta, gamma)

    eps = 1e-5
    H2_p = hubble_squared(a * (1 + eps), psi, dpsi_dt_est, beta, gamma)
    H2_m = hubble_squared(a * (1 - eps), psi, dpsi_dt_est, beta, gamma)
    dlnH2 = (np.log(max(H2_p, 1e-60)) - np.log(max(H2_m, 1e-60))) / (2.0 * eps)
    Hdot_over_H = 0.5 * dlnH2

    W = W_exact(psi, beta, gamma)
    # FRW field eq: psi_tt + 3H psi_t + 3 psi_t^2/psi = W(psi)
    # In chi = dpsi/dlna: chi' = W/H^2 - 3chi - 3chi^2/psi - (Hdot/H)chi
    chi_prime = W / H2 - 3.0 * chi - 3.0 * chi**2 / psi - Hdot_over_H * chi
    return [chi, chi_prime]


# ============================================================================
# Integration driver
# ============================================================================
def solve_field(Phi0, lna_span=None, n_pts=4000):
    """Integrate exact field equation. Returns (a, psi, chi, H, dpsi_dt)."""
    beta, gamma, _ = tgp_params(Phi0)
    if lna_span is None:
        lna_span = (np.log(1.0 / (1.0 + 1e4)), 0.0)

    lna_eval = np.linspace(lna_span[0], lna_span[1], n_pts)
    sol = solve_ivp(
        lambda lna, y: ode_exact(lna, y, beta, gamma),
        lna_span, [1.0, 0.0],
        t_eval=lna_eval, method="RK45",
        rtol=1e-10, atol=1e-12, max_step=0.01,
    )
    if sol.status != 0:
        print(f"  [warn] Integration Phi0={Phi0}: {sol.message}")

    a_arr   = np.exp(sol.t)
    psi_arr = sol.y[0]
    chi_arr = sol.y[1]
    H_arr      = np.zeros_like(a_arr)
    dpsi_dt_arr = np.zeros_like(a_arr)
    for i in range(len(a_arr)):
        H2 = hubble_squared(a_arr[i], psi_arr[i], 0.0, beta, gamma)
        H_arr[i] = np.sqrt(H2)
        dpsi_dt_arr[i] = H_arr[i] * chi_arr[i]
    return a_arr, psi_arr, chi_arr, H_arr, dpsi_dt_arr


# ============================================================================
# w_DE from field energy (kinetic + potential)
# ============================================================================
def compute_w_de(a_arr, psi_arr, dpsi_dt_arr, beta, gamma):
    """w_DE = p_DE / rho_DE from field kinetic + potential energy."""
    w = np.full_like(a_arr, -1.0)
    for i in range(len(a_arr)):
        kin = 0.5 * (dpsi_dt_arr[i] / c0)**2
        pot = U_eff(psi_arr[i], beta, gamma)
        rho = kin + pot
        if abs(rho) > 1e-60:
            w[i] = (kin - pot) / rho
    return w


# ============================================================================
# CPL fit: w(a) = w0 + wa(1-a)
# ============================================================================
def fit_w0_wa(a_arr, w_arr):
    mask = (a_arr > 0.3) & (a_arr < 1.0) & np.isfinite(w_arr) & (np.abs(w_arr) < 5)
    if mask.sum() < 5:
        return -1.0, 0.0
    A = np.column_stack([np.ones(mask.sum()), 1.0 - a_arr[mask]])
    res, _, _, _ = np.linalg.lstsq(A, w_arr[mask], rcond=None)
    return res[0], res[1]


# ============================================================================
# Scale-dependent G_eff(k, a)
# ============================================================================
def G_eff_ratio(k_H0, a, gamma):
    """
    G_eff(k,a)/G0 for beta=gamma.
    k_H0 is wavenumber in units of H0 (dimensionless).
    m_sp^2 = gamma (for beta=gamma), alpha_eff = q/(4pi) with q = 8piG0/c0^2.
    k_phys = k_H0 * H0_SI / a  (physical wavenumber).
    """
    q = 8.0 * np.pi * G0 / c0**2
    alpha_eff = q / (4.0 * np.pi)
    m_sp = np.sqrt(max(gamma, 0.0))
    # ratio = (a * m_sp / k_phys)^2, but m_sp is in "code" units matching k_H0
    # We need consistent units: m_sp has units from gamma = Phi0*H0^2/c0^2
    # In dimensionless form with H0=1: m_sp_dimless = m_sp / H0_SI
    # Actually m_sp^2 = gamma which has units of [s^-2 * m^-2 * ...] -- let's
    # keep it simple: use the ratio in H0 units
    m_sp_H0 = np.sqrt(gamma) / H0_SI  # in units of H0
    ratio = (a * m_sp_H0 / k_H0)**2
    return 1.0 + 2.0 * alpha_eff**2 / (1.0 + ratio)


# ============================================================================
# Growth rate f*sigma8
# ============================================================================
def compute_growth(a_arr, H_arr, gamma, k_ref=0.1):
    """
    Solve growth ODE with G_eff, return (a_g, fsigma8, f_growth).
    sigma8 normalised to 0.811 at a=1.
    """
    lna = np.log(a_arr)
    H_interp = interp1d(lna, H_arr, kind="cubic", fill_value="extrapolate")

    def growth_ode(lna_val, y):
        a = np.exp(lna_val)
        H = float(H_interp(lna_val))
        if H < 1e-30:
            H = 1e-30
        dl = 1e-5
        Hp = float(H_interp(lna_val + dl))
        Hm = float(H_interp(lna_val - dl))
        dlnH = (np.log(max(Hp, 1e-60)) - np.log(max(Hm, 1e-60))) / (2 * dl)
        Omega_m_a = Omega_m0 * H0_SI**2 / (a**3 * H**2)
        Geff = G_eff_ratio(k_ref, a, gamma)
        delta, dp = y
        dpp = -(2.0 + dlnH) * dp + 1.5 * Omega_m_a * Geff * delta
        return [dp, dpp]

    idx0 = 10
    sol = solve_ivp(growth_ode, (lna[idx0], lna[-1]), [1.0, 1.0],
                    t_eval=lna[idx0:], method="RK45", rtol=1e-8, atol=1e-10)
    a_g = np.exp(sol.t)
    delta_g = sol.y[0]
    dp_g = sol.y[1]
    f_growth = dp_g / delta_g
    norm = 0.811 / np.interp(1.0, a_g, delta_g)
    sigma8_a = delta_g * norm
    return a_g, f_growth * sigma8_a, f_growth


# ============================================================================
# Observational data (from cosmological_evolution.py)
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

DATA_FSIG8 = np.array([
    [0.02,0.428,0.0465],[0.10,0.370,0.130],[0.15,0.490,0.145],
    [0.17,0.510,0.060],[0.18,0.360,0.090],[0.38,0.440,0.060],
    [0.25,0.3512,0.0583],[0.37,0.4602,0.0378],[0.32,0.384,0.095],
    [0.44,0.413,0.080],[0.51,0.455,0.039],[0.57,0.441,0.043],
    [0.60,0.390,0.063],[0.61,0.457,0.052],[0.73,0.437,0.072],
    [0.78,0.380,0.040],[0.80,0.470,0.080],[0.85,0.315,0.095],
    [0.978,0.379,0.176],[1.05,0.280,0.080],[1.40,0.482,0.116],
    [1.48,0.462,0.045],[1.52,0.420,0.076],
])


# ============================================================================
# MCMC model predictions
# ============================================================================
def model_Hz(z_arr, Phi0):
    """Compute H(z) in km/s/Mpc from TGP for a given Phi0."""
    try:
        a, psi, chi, H, dpsi_dt = solve_field(Phi0, n_pts=1500)
        z_mod = 1.0 / a - 1.0
        # Convert H from s^-1 to km/s/Mpc:  H [s^-1] * (Mpc_in_m / 1000)
        Mpc_m = 3.0857e22
        H_km = H * Mpc_m / 1e3
        valid = np.isfinite(H_km) & (a > 0)
        interp = interp1d(z_mod[valid], H_km[valid], kind="linear",
                          fill_value="extrapolate", bounds_error=False)
        return interp(z_arr)
    except Exception:
        return np.full_like(z_arr, np.nan)


def model_fsigma8(z_arr, Phi0):
    """Compute f*sigma8(z) from TGP for a given Phi0."""
    try:
        _, gamma, _ = tgp_params(Phi0)
        a, psi, chi, H, dpsi_dt = solve_field(Phi0, n_pts=1500)
        a_g, fsig8, _ = compute_growth(a, H, gamma)
        z_g = 1.0 / a_g - 1.0
        valid = np.isfinite(fsig8) & (a_g > 0)
        if valid.sum() < 5:
            return np.full_like(z_arr, np.nan)
        interp = interp1d(z_g[valid], fsig8[valid], kind="linear",
                          fill_value="extrapolate", bounds_error=False)
        return interp(z_arr)
    except Exception:
        return np.full_like(z_arr, np.nan)


# ============================================================================
# MCMC
# ============================================================================
def log_prior_mcmc(theta):
    """Flat prior on Phi0 in [1, 200]."""
    (Phi0,) = theta
    if 1.0 < Phi0 < 200.0:
        return 0.0
    return -np.inf


def log_likelihood_mcmc(theta):
    (Phi0,) = theta
    chi2 = 0.0
    # H(z)
    Hz_mod = model_Hz(DATA_HZ[:, 0], Phi0)
    if np.any(np.isnan(Hz_mod)):
        return -np.inf
    chi2 += np.sum(((DATA_HZ[:, 1] - Hz_mod) / DATA_HZ[:, 2])**2)
    # f*sigma8
    fs8_mod = model_fsigma8(DATA_FSIG8[:, 0], Phi0)
    if np.any(np.isnan(fs8_mod)):
        return -np.inf
    chi2 += np.sum(((DATA_FSIG8[:, 1] - fs8_mod) / DATA_FSIG8[:, 2])**2)
    return -0.5 * chi2


def log_posterior_mcmc(theta):
    lp = log_prior_mcmc(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood_mcmc(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def run_mcmc(nwalkers=16, nsteps=1000, burn_in=300, save_dir=None):
    """Run MCMC fitting Phi0 to H(z) + f*sigma8(z) data."""
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    import emcee

    os.makedirs(save_dir, exist_ok=True)
    ndim = 1
    # Initial positions scattered around Phi0 ~ 20
    p0 = 20.0 + 2.0 * np.random.randn(nwalkers, ndim)
    p0 = np.clip(p0, 2.0, 190.0)

    print(f"\n  Running MCMC: {nwalkers} walkers, {nsteps} steps, fitting Phi0 ...")
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_mcmc)
    sampler.run_mcmc(p0, nsteps, progress=True)

    flat = sampler.get_chain(discard=burn_in, flat=True)
    print(f"  Chain shape after burn-in: {flat.shape}")

    q16, q50, q84 = np.percentile(flat[:, 0], [16, 50, 84])
    print(f"\n  Phi0 = {q50:.2f}  (+{q84-q50:.2f} / -{q50-q16:.2f})")
    Phi0_best = q50
    beta, gamma, tau0 = tgp_params(Phi0_best)
    print(f"  gamma = beta = {gamma:.4e}")
    print(f"  tau0 = {tau0:.4e} s")

    # --- Trace plot ---
    chain_full = sampler.get_chain()
    fig_t, ax_t = plt.subplots(figsize=(10, 4))
    ax_t.plot(chain_full[:, :, 0], alpha=0.2, color="C0", lw=0.5)
    ax_t.axvline(burn_in, color="r", ls="--", lw=0.8, label="burn-in")
    ax_t.set_xlabel("Step"); ax_t.set_ylabel(r"$\Phi_0$")
    ax_t.set_title("MCMC trace"); ax_t.legend(); ax_t.grid(True, ls=":", alpha=0.3)
    fig_t.tight_layout()
    p = os.path.join(save_dir, "tgp_mcmc_trace.png")
    fig_t.savefig(p, dpi=150); print(f"  Saved {p}"); plt.close(fig_t)

    # --- Corner / histogram ---
    try:
        import corner
        fig_c = corner.corner(flat, labels=[r"$\Phi_0$"],
                              quantiles=[0.16, 0.5, 0.84], show_titles=True)
        p = os.path.join(save_dir, "tgp_mcmc_corner.png")
        fig_c.savefig(p, dpi=150, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig_c)
    except ImportError:
        fig_h, ax_h = plt.subplots(figsize=(6, 4))
        ax_h.hist(flat[:, 0], bins=40, density=True, color="C0", alpha=0.7)
        ax_h.axvline(q50, color="r", label=f"median={q50:.1f}")
        ax_h.set_xlabel(r"$\Phi_0$"); ax_h.set_title("Posterior"); ax_h.legend()
        fig_h.tight_layout()
        p = os.path.join(save_dir, "tgp_mcmc_corner.png")
        fig_h.savefig(p, dpi=150); print(f"  Saved {p}"); plt.close(fig_h)

    # --- Best-fit vs data ---
    fig_bf, (ax_h, ax_f) = plt.subplots(1, 2, figsize=(14, 6))
    z_plot = np.linspace(0.01, 2.5, 100)

    n_draw = min(40, len(flat))
    for idx in np.random.choice(len(flat), n_draw, replace=False):
        Ph = flat[idx, 0]
        Hz = model_Hz(z_plot, Ph)
        if np.all(np.isfinite(Hz)):
            ax_h.plot(z_plot, Hz, color="C0", alpha=0.05, lw=0.5)
        fs = model_fsigma8(z_plot, Ph)
        if np.all(np.isfinite(fs)):
            ax_f.plot(z_plot, fs, color="C0", alpha=0.05, lw=0.5)

    Hz_best = model_Hz(z_plot, Phi0_best)
    fs_best = model_fsigma8(z_plot, Phi0_best)
    ax_h.plot(z_plot, Hz_best, "C1-", lw=2, label=f"TGP best ($\\Phi_0={Phi0_best:.1f}$)")
    ax_f.plot(z_plot, fs_best, "C1-", lw=2, label=f"TGP best ($\\Phi_0={Phi0_best:.1f}$)")

    # LCDM reference
    a_ref = 1.0 / (1.0 + z_plot)
    Mpc_m = 3.0857e22
    ax_h.plot(z_plot, hubble_lcdm(a_ref) * Mpc_m / 1e3, "k--", lw=1.2, label=r"$\Lambda$CDM")

    ax_h.errorbar(DATA_HZ[:,0], DATA_HZ[:,1], yerr=DATA_HZ[:,2],
                  fmt="ro", ms=3, elinewidth=0.8, capsize=2, label="Data")
    ax_f.errorbar(DATA_FSIG8[:,0], DATA_FSIG8[:,1], yerr=DATA_FSIG8[:,2],
                  fmt="ro", ms=3, elinewidth=0.8, capsize=2, label="Data")

    for ax, yl, tl in [(ax_h, r"$H(z)$ [km/s/Mpc]", "Hubble parameter"),
                        (ax_f, r"$f\sigma_8(z)$", "Growth rate")]:
        ax.set_xlabel(r"$z$"); ax.set_ylabel(yl); ax.set_title(tl)
        ax.legend(fontsize=9); ax.grid(True, ls=":", alpha=0.4)

    fig_bf.suptitle("TGP MCMC fit to observational data", fontsize=14, y=1.01)
    fig_bf.tight_layout()
    p = os.path.join(save_dir, "tgp_mcmc_bestfit.png")
    fig_bf.savefig(p, dpi=180, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig_bf)

    # Derived quantities
    a_bf, psi_bf, _, H_bf, dpsi_bf = solve_field(Phi0_best)
    w_bf = compute_w_de(a_bf, psi_bf, dpsi_bf, beta, gamma)
    w0, wa = fit_w0_wa(a_bf, w_bf)
    print(f"\n  Best-fit derived: w0 = {w0:.4f}, wa = {wa:.4f}, psi(today) = {psi_bf[-1]:.6f}")
    return sampler, flat, Phi0_best


# ============================================================================
# Plotting: scan mode
# ============================================================================
def plot_w_de(results, save_dir):
    """Plot w_DE(z) for different Phi0."""
    os.makedirs(save_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results)))
    for (Phi0, a, psi, w_de), col in zip(results, colors):
        z = 1.0 / a - 1.0
        m = (z > 0) & (z < 5) & np.isfinite(w_de)
        ax.plot(z[m], w_de[m], color=col, lw=1.8, label=rf"$\Phi_0 = {Phi0}$")
    ax.axhline(-1, color="k", ls="--", lw=1, label=r"$\Lambda$CDM")
    z_b = np.linspace(0, 3, 200)
    ax.fill_between(z_b, -1 - 0.05*(1+0.3*z_b), -1 + 0.05*(1+0.3*z_b),
                    color="orange", alpha=0.15, label="DESI (schematic)")
    ax.fill_between(z_b, -1 - 0.03*(1+0.2*z_b), -1 + 0.03*(1+0.2*z_b),
                    color="blue", alpha=0.10, label="Euclid (schematic)")
    ax.set_xlabel(r"$z$", fontsize=13); ax.set_ylabel(r"$w_{\rm DE}(z)$", fontsize=13)
    ax.set_title("TGP dark energy equation of state (exact equation)", fontsize=14)
    ax.set_xlim(0, 3); ax.set_ylim(-1.3, -0.7)
    ax.legend(fontsize=9, loc="lower left"); ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    p = os.path.join(save_dir, "tgp_w_de.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig)


def plot_w0_wa(results, save_dir):
    """Plot w0-wa plane."""
    os.makedirs(save_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results)))
    for (Phi0, a, psi, w_de), col in zip(results, colors):
        w0, wa = fit_w0_wa(a, w_de)
        ax.plot(w0, wa, "o", ms=10, color=col,
                label=rf"$\Phi_0={Phi0}$ ($w_0={w0:.3f}, w_a={wa:.3f}$)")
    ax.plot(-1, 0, "k*", ms=16, zorder=6, label=r"$\Lambda$CDM")
    th = np.linspace(0, 2*np.pi, 200)
    for ns, ls in [(1, "-"), (2, "--")]:
        ax.plot(-1 + ns*0.08*np.cos(th), ns*0.35*np.sin(th), color="gray", ls=ls, lw=0.8)
    ax.fill(-1 + 2*0.08*np.cos(th), 2*0.35*np.sin(th),
            color="gray", alpha=0.06, label=r"Planck+BAO $2\sigma$ (schematic)")
    ax.plot(-0.75 + 0.12*np.cos(th), -0.8 + 0.5*np.sin(th), "r-", lw=1)
    ax.fill(-0.75 + 0.12*np.cos(th), -0.8 + 0.5*np.sin(th),
            color="red", alpha=0.08, label=r"DESI 2024 hint (schematic)")
    ax.set_xlabel(r"$w_0$", fontsize=14); ax.set_ylabel(r"$w_a$", fontsize=14)
    ax.set_title(r"TGP predictions in the $w_0$--$w_a$ plane", fontsize=14)
    ax.set_xlim(-1.5, -0.4); ax.set_ylim(-2, 1)
    ax.legend(fontsize=8, loc="upper left"); ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    p = os.path.join(save_dir, "tgp_w0_wa.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig)


def plot_psi(results, save_dir):
    """Plot psi(z)."""
    os.makedirs(save_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results)))
    for (Phi0, a, psi, _), col in zip(results, colors):
        z = 1.0 / a - 1.0
        lab = rf"$\Phi_0 = {Phi0}$"
        m1 = (z > 0) & (z < 10)
        ax1.plot(z[m1], psi[m1], color=col, lw=1.8, label=lab)
        m2 = z > 0
        ax2.plot(z[m2], psi[m2], color=col, lw=1.8, label=lab)
    for ax in (ax1, ax2):
        ax.axhline(1, color="k", ls="--", lw=0.8)
        ax.set_xlabel(r"$z$", fontsize=13); ax.set_ylabel(r"$\psi(z)$", fontsize=13)
        ax.legend(fontsize=9); ax.grid(True, ls=":", alpha=0.4)
    ax1.set_xlim(0, 10); ax1.set_title("Field evolution (linear)", fontsize=14)
    ax2.set_xscale("log"); ax2.set_title("Field evolution (log)", fontsize=14)
    fig.tight_layout()
    p = os.path.join(save_dir, "tgp_psi.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig)


def plot_growth(results_growth, gamma_first, save_dir):
    """Plot f*sigma8(z) and G_eff(k)."""
    os.makedirs(save_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results_growth)))

    for (Phi0, a_g, fsig8), col in zip(results_growth, colors):
        z_g = 1.0 / a_g - 1.0
        m = (z_g > 0) & (z_g < 2)
        ax1.plot(z_g[m], fsig8[m], color=col, lw=1.5, label=rf"$\Phi_0={Phi0}$")

    # LCDM reference growth
    a_ref = np.linspace(0.01, 1.0, 2000)
    H_ref = np.array([hubble_lcdm(aa) for aa in a_ref])
    try:
        a_g0, fs0, _ = compute_growth(a_ref, H_ref, 0.0)
        z0 = 1.0 / a_g0 - 1.0
        m0 = (z0 > 0) & (z0 < 2)
        ax1.plot(z0[m0], fs0[m0], "k--", lw=1.2, label=r"$\Lambda$CDM")
    except Exception:
        pass

    # Data points
    ax1.errorbar(DATA_FSIG8[:,0], DATA_FSIG8[:,1], yerr=DATA_FSIG8[:,2],
                 fmt="ro", ms=4, elinewidth=0.8, capsize=2, label="Data", zorder=5)
    ax1.set_xlabel(r"$z$", fontsize=13); ax1.set_ylabel(r"$f\sigma_8(z)$", fontsize=13)
    ax1.set_title("Growth rate", fontsize=14)
    ax1.legend(fontsize=8); ax1.grid(True, ls=":", alpha=0.4)

    # G_eff(k) for first Phi0 at several epochs
    k_arr = np.logspace(-4, 0, 300)
    for lab, a_ep, ls in [("a=1.0", 1.0, "-"), ("a=0.5", 0.5, "--"), ("a=0.2", 0.2, ":")]:
        Geff = np.array([G_eff_ratio(k, a_ep, gamma_first) for k in k_arr])
        ax2.plot(k_arr, Geff, ls=ls, lw=1.5, label=lab)
    ax2.set_xscale("log"); ax2.axhline(1, color="k", ls="-", lw=0.5)
    ax2.set_xlabel(r"$k$ [$H_0$]", fontsize=13)
    ax2.set_ylabel(r"$G_{\rm eff}/G_0$", fontsize=13)
    ax2.set_title(r"Scale-dependent $G_{\rm eff}$", fontsize=14)
    ax2.legend(fontsize=9); ax2.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    p = os.path.join(save_dir, "tgp_growth.png")
    fig.savefig(p, dpi=180, bbox_inches="tight"); print(f"  Saved {p}"); plt.close(fig)


# ============================================================================
# Main: scan mode
# ============================================================================
def main_scan():
    print("=" * 65)
    print("TGP consolidated cosmological script -- scan mode")
    print("=" * 65)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    results_plot = []   # (Phi0, a, psi, w_de)
    results_growth = [] # (Phi0, a_g, fsig8)

    for Phi0 in PHI0_VALUES:
        beta, gamma, tau0 = tgp_params(Phi0)
        print(f"\n  Phi0 = {Phi0}:  gamma = {gamma:.4e},  tau0 = {tau0:.4e} s")

        a, psi, chi, H, dpsi_dt = solve_field(Phi0)
        w_de = compute_w_de(a, psi, dpsi_dt, beta, gamma)
        w0, wa = fit_w0_wa(a, w_de)
        print(f"    psi(today) = {psi[-1]:.8f},  w0 = {w0:.4f},  wa = {wa:.4f}")
        results_plot.append((Phi0, a, psi, w_de))

        # Growth rate
        try:
            a_g, fsig8, _ = compute_growth(a, H, gamma)
            results_growth.append((Phi0, a_g, fsig8))
        except Exception as e:
            print(f"    [warn] growth failed: {e}")

    print("\n  Generating plots ...")
    plot_w_de(results_plot, save_dir)
    plot_w0_wa(results_plot, save_dir)
    plot_psi(results_plot, save_dir)
    if results_growth:
        _, gamma0, _ = tgp_params(PHI0_VALUES[0])
        plot_growth(results_growth, gamma0, save_dir)

    # Summary table
    print("\n  " + "=" * 65)
    print(f"  {'Phi0':>6s}  {'gamma':>12s}  {'psi(0)':>10s}  {'w0':>8s}  {'wa':>8s}")
    print("  " + "-" * 65)
    for (Phi0, a, psi, w_de) in results_plot:
        w0, wa = fit_w0_wa(a, w_de)
        _, gam, _ = tgp_params(Phi0)
        print(f"  {Phi0:6d}  {gam:12.4e}  {psi[-1]:10.6f}  {w0:8.4f}  {wa:8.4f}")
    print("  " + "=" * 65)
    print("\nDone.")


def main_mcmc():
    print("=" * 65)
    print("TGP consolidated cosmological script -- MCMC mode")
    print("=" * 65)
    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    run_mcmc(nwalkers=16, nsteps=500, burn_in=150, save_dir=save_dir)
    print("\nMCMC Done.")


if __name__ == "__main__":
    if "--mcmc" in sys.argv:
        main_mcmc()
    else:
        main_scan()
