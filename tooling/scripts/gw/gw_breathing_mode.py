"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
gw_breathing_mode.py  --  Theory of Generated Space (TGP)
=========================================================
Compute and visualize gravitational wave properties in TGP theory.

Physics
-------
In TGP, the conformal metric g_mu_nu = Phi * eta_mu_nu generates scalar
(breathing mode) gravitational waves from perturbations delta_Phi of the
background field.  The perturbation obeys a massive Klein-Gordon equation:

    Dispersion relation:  omega^2 = c0^2 (k^2/a^2 + m_sp^2)

where the spin-0 mass parameter is:

    m_sp^2 = 3 gamma - 2 beta          [CORRECTED: no 1/Phi0 factor]

Key derived quantities:
    Phase velocity:   v_phi = c0 sqrt(1 + (a m_sp / k)^2)
    Group velocity:   v_g   = c0 / sqrt(1 + (a m_sp / k)^2)
    Low-freq cutoff:  f_cut = c0 m_sp / (2 pi)

Polarization structure (sek08, Theorem thm:no-tensor + Hypothesis hyp:disformal):
    - The CONFORMAL metric g_mu_nu = Phi * eta_mu_nu CANNOT generate tensor
      modes h_+, h_x.  This is proven (Theorem thm:no-tensor): a conformally
      flat metric has vanishing Weyl tensor, so no transverse-traceless
      perturbations exist.
    - The breathing mode (scalar, transverse trace) is the PRIMARY/fundamental
      GW polarization in TGP.
    - HOWEVER, the DISFORMAL extension:
          g_mu_nu = A(Phi) eta_mu_nu + B(Phi) partial_mu Phi partial_nu Phi / M*^4
      CAN generate tensor modes h_+, h_x from source anisotropy
      (Hypothesis hyp:disformal).  The disformal term breaks conformal
      flatness and introduces a preferred direction (the gradient of Phi),
      which sources transverse-traceless perturbations when the source
      itself is anisotropic (e.g., binary inspiral).
    - The ratio h_breathing / h_tensor depends on the disformal coupling B/M*^4.

Cosmological vs. nuclear scales:
    - With nuclear-scale r_0 ~ 1e-15 m, the hatted parameters give enormous
      m_sp and unphysically high f_cut.  This is an artifact of using
      dimensionless beta_hat = beta * r_0^2 at nuclear scale.
    - For COSMOLOGICAL GW predictions, the natural scale is:
          gamma_cosmo = Phi_0 * H_0^2 / c_0^2
      giving m_sp ~ H_0/c_0 ~ 1e-26 m^-1 and f_cut ~ 1e-18 Hz.
    - At LIGO frequencies (100 Hz), the deviation is ~ 1e-44, trivially
      satisfying GW170817 constraints.

Honest assessment:
    - With natural (cosmological) parameters, the TGP breathing mode is
      INDISTINGUISHABLE from GR at LIGO/Virgo/KAGRA frequencies.
    - The breathing mode is detectable only if an additional mechanism
      enhances m_sp beyond the cosmological scale.
    - The disformal coupling is NEEDED to reproduce the observed tensor
      polarization pattern in GW events.

Observational constraint (GW170817):
    |c_GW - c0| / c0 < 1e-15
    This bounds the combination (a m_sp / k)^2 at LIGO frequencies.

Outputs  (saved to scripts/plots/):
    gw_dispersion.png           -- omega(k) vs GR massless relation
    gw_velocities.png           -- phase and group velocity vs frequency
    gw_cutoff_vs_gamma.png      -- f_cut as function of gamma
    gw_gw170817.png             -- GW170817 exclusion region in parameter space
    gw_disformal_tensor.png     -- disformal tensor mode contribution
    gw_cosmo_assessment.png     -- cosmological-scale analysis
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════════════════
# Physical constants (SI)
# ═══════════════════════════════════════════════════════════════════════════
c0 = 3.0e8          # m/s
G0 = 6.674e-11      # m^3 kg^-1 s^-2
r0 = 1.0e-15        # m  (nuclear scale for dimensionless hat parameters)
Phi0 = 115.0         # background field (dimensionless; ~115 from TGP cosmology)
H0 = 2.2e-18        # Hubble constant [1/s] (~67.4 km/s/Mpc)

# ═══════════════════════════════════════════════════════════════════════════
# TGP parameters (dimensionless, hatted: beta_hat = beta * r0^2, etc.)
# ═══════════════════════════════════════════════════════════════════════════
BETA_HAT_DEFAULT  = 0.02
GAMMA_HAT_DEFAULT = 0.03


# ═══════════════════════════════════════════════════════════════════════════
# Core physics functions
# ═══════════════════════════════════════════════════════════════════════════
def m_sp_squared(beta_hat, gamma_hat):
    """
    Spin-0 mass parameter squared (in units of 1/r0^2).

    CORRECTED formula:  m_sp^2 = 3 gamma - 2 beta
    (NOT divided by Phi0 -- that was an earlier error in the manuscript.)

    Physical units: m_sp_phys^2 = (3 gamma_hat - 2 beta_hat) / r0^2
    """
    return 3.0 * gamma_hat - 2.0 * beta_hat


def m_sp_physical(beta_hat, gamma_hat):
    """m_sp in SI units [1/m], i.e. inverse Compton wavelength."""
    ms2 = m_sp_squared(beta_hat, gamma_hat)
    if ms2 <= 0:
        return 0.0
    return np.sqrt(ms2) / r0


def omega_dispersion(k, beta_hat, gamma_hat, a=1.0):
    """
    Dispersion relation:  omega(k) = c0 sqrt(k^2/a^2 + m_sp_phys^2)

    Parameters
    ----------
    k : array-like, wavenumber [1/m]
    a : float, scale factor (default 1 for today)

    Returns
    -------
    omega : angular frequency [rad/s]
    """
    m = m_sp_physical(beta_hat, gamma_hat)
    return c0 * np.sqrt((k / a)**2 + m**2)


def phase_velocity(k, beta_hat, gamma_hat, a=1.0):
    """v_phi = omega / k = c0 sqrt(1 + (a m_sp / k)^2)"""
    m = m_sp_physical(beta_hat, gamma_hat)
    return c0 * np.sqrt(1.0 + (a * m / k)**2)


def group_velocity(k, beta_hat, gamma_hat, a=1.0):
    """v_g = d omega / d k = c0 / sqrt(1 + (a m_sp / k)^2)"""
    m = m_sp_physical(beta_hat, gamma_hat)
    return c0 / np.sqrt(1.0 + (a * m / k)**2)


def f_cutoff(beta_hat, gamma_hat):
    """Low-frequency cutoff: f_cut = c0 m_sp / (2 pi)  [Hz]"""
    m = m_sp_physical(beta_hat, gamma_hat)
    return c0 * m / (2.0 * np.pi)


# ═══════════════════════════════════════════════════════════════════════════
# Cosmological-scale mass parameter
# ═══════════════════════════════════════════════════════════════════════════
def m_sp_cosmological():
    """
    Cosmological spin-0 mass parameter.

    With natural (cosmological) parameters:
        gamma_cosmo = Phi_0 * H_0^2 / c_0^2

    This gives m_sp ~ H_0 / c_0 ~ 1e-26 m^-1.

    Returns
    -------
    m_sp_cosmo : float, [1/m]
    gamma_cosmo : float, dimensionless
    f_cut_cosmo : float, [Hz]
    """
    gamma_cosmo = Phi0 * H0**2 / c0**2
    # For cosmological case, beta_cosmo ~ same order, take beta ~ gamma/2
    # so m_sp^2 = 3*gamma - 2*beta ~ 2*gamma (order of magnitude)
    m_sp_cosmo = np.sqrt(3.0 * gamma_cosmo) * (H0 / c0)
    # More directly: m_sp ~ sqrt(Phi_0) * H_0 / c_0
    m_sp_cosmo = np.sqrt(Phi0) * H0 / c0
    f_cut_cosmo = c0 * m_sp_cosmo / (2.0 * np.pi)
    return m_sp_cosmo, gamma_cosmo, f_cut_cosmo


# ═══════════════════════════════════════════════════════════════════════════
# Disformal tensor mode contribution
# ═══════════════════════════════════════════════════════════════════════════
def disformal_tensor_amplitude(h_breathing, B_over_Mstar4, source_anisotropy=1.0):
    """
    Estimate tensor mode amplitude from disformal coupling.

    The disformal metric extension:
        g_mu_nu = A(Phi) eta_mu_nu + B(Phi) (partial_mu Phi)(partial_nu Phi) / M*^4

    generates tensor modes h_+, h_x when the source is anisotropic.
    The tensor amplitude scales as:

        h_tensor ~ B / M*^4 * (grad Phi)^2 * Q_anisotropy * h_breathing

    where Q_anisotropy characterizes the source quadrupole anisotropy
    (Q ~ 1 for a binary inspiral, Q -> 0 for isotropic sources).

    Parameters
    ----------
    h_breathing : float or array, breathing mode strain amplitude
    B_over_Mstar4 : float, disformal coupling B/M*^4 [dimensionless in natural units]
    source_anisotropy : float, Q_anisotropy in [0, 1]

    Returns
    -------
    h_tensor : float or array, tensor mode amplitude (h_+ or h_x)
    """
    # The tensor mode is sourced by the disformal term acting on the
    # anisotropic part of the scalar field gradient.
    # Parametric estimate: h_tensor ~ B/M*^4 * source_anisotropy * h_breathing
    h_tensor = B_over_Mstar4 * source_anisotropy * h_breathing
    return h_tensor


def breathing_to_tensor_ratio(B_over_Mstar4, source_anisotropy=1.0):
    """
    Ratio h_breathing / h_tensor as function of disformal coupling.

    For GR-like behavior (tensor-dominated), need B/M*^4 >> 1.
    For pure TGP (breathing-dominated), B/M*^4 << 1.

    Returns
    -------
    ratio : float, h_breathing / h_tensor
    """
    if B_over_Mstar4 * source_anisotropy == 0:
        return np.inf
    return 1.0 / (B_over_Mstar4 * source_anisotropy)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1: Dispersion relation omega(k) vs GR
# ═══════════════════════════════════════════════════════════════════════════
def plot_dispersion(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Wavenumber range (in units of m_sp for clarity)
    m_ref = m_sp_physical(BETA_HAT_DEFAULT, GAMMA_HAT_DEFAULT)
    k_over_m = np.logspace(-1, 3, 500)
    k = k_over_m * m_ref

    # GR: massless  omega = c0 k
    omega_GR = c0 * k

    # TGP dispersion for several parameter choices
    param_sets = [
        (0.02, 0.03, r"$\hat\beta=0.02,\;\hat\gamma=0.03$"),
        (0.01, 0.02, r"$\hat\beta=0.01,\;\hat\gamma=0.02$"),
        (0.05, 0.10, r"$\hat\beta=0.05,\;\hat\gamma=0.10$"),
        (0.10, 0.10, r"$\hat\beta=0.10,\;\hat\gamma=0.10$"),
    ]
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(param_sets)))

    # GR reference
    ax.plot(k_over_m, omega_GR / (c0 * m_ref), "k--", lw=2,
            label=r"GR: $\omega = c_0 k$ (massless)")

    for (bh, gh, label), col in zip(param_sets, colors):
        ms2 = m_sp_squared(bh, gh)
        if ms2 <= 0:
            continue
        m_phys = np.sqrt(ms2) / r0
        k_loc = k_over_m * m_phys
        omega = omega_dispersion(k_loc, bh, gh)
        ax.plot(k_over_m, omega / (c0 * m_phys), color=col, lw=1.8, label=label)

    # Mark the mass gap
    ax.axhline(1.0, color="gray", ls=":", lw=0.8, alpha=0.6)
    ax.annotate(r"$\omega_{\rm min} = c_0 m_{\rm sp}$ (mass gap)",
                xy=(0.15, 1.02), fontsize=10, color="gray", style="italic")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$k / m_{\rm sp}$", fontsize=13)
    ax.set_ylabel(r"$\omega / (c_0 \, m_{\rm sp})$", fontsize=13)
    ax.set_title("TGP gravitational wave dispersion relation", fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, ls=":", alpha=0.4)
    ax.set_xlim(k_over_m[0], k_over_m[-1])

    fig.tight_layout()
    path = os.path.join(save_dir, "gw_dispersion.png")
    fig.savefig(path, dpi=180)
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2: Phase and group velocity vs frequency
# ═══════════════════════════════════════════════════════════════════════════
def plot_velocities(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    param_sets = [
        (0.02, 0.03, r"$\hat\beta=0.02,\;\hat\gamma=0.03$"),
        (0.01, 0.02, r"$\hat\beta=0.01,\;\hat\gamma=0.02$"),
        (0.05, 0.10, r"$\hat\beta=0.05,\;\hat\gamma=0.10$"),
    ]
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(param_sets)))

    for (bh, gh, label), col in zip(param_sets, colors):
        ms2 = m_sp_squared(bh, gh)
        if ms2 <= 0:
            continue
        f_cut = f_cutoff(bh, gh)

        # Frequency range: from just above cutoff to 10^4 * cutoff
        f_arr = np.logspace(np.log10(f_cut * 1.01), np.log10(f_cut * 1e4), 500)
        k_arr = 2.0 * np.pi * f_arr / c0  # approximate k from f for high-k
        # More precisely, k from dispersion: k = (1/c0) sqrt(omega^2 - c0^2 m_sp^2)
        m_phys = m_sp_physical(bh, gh)
        omega_arr = 2.0 * np.pi * f_arr
        k_arr = np.sqrt(omega_arr**2 - (c0 * m_phys)**2) / c0

        v_ph = phase_velocity(k_arr, bh, gh)
        v_gr = group_velocity(k_arr, bh, gh)

        ax1.plot(f_arr / f_cut, v_ph / c0, color=col, lw=1.8, label=label)
        ax2.plot(f_arr / f_cut, v_gr / c0, color=col, lw=1.8, label=label)

        # Mark cutoff
        ax1.axvline(1.0, color=col, ls=":", lw=0.6, alpha=0.5)
        ax2.axvline(1.0, color=col, ls=":", lw=0.6, alpha=0.5)

    # GR reference
    ax1.axhline(1.0, color="k", ls="--", lw=1.2, label=r"GR: $v_\varphi = c_0$")
    ax2.axhline(1.0, color="k", ls="--", lw=1.2, label=r"GR: $v_g = c_0$")

    ax1.set_xscale("log")
    ax1.set_xlabel(r"$f / f_{\rm cut}$", fontsize=13)
    ax1.set_ylabel(r"$v_\varphi / c_0$", fontsize=13)
    ax1.set_title("Phase velocity (superluminal)", fontsize=14)
    ax1.legend(fontsize=9)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.set_ylim(0.95, 3.0)

    ax2.set_xscale("log")
    ax2.set_xlabel(r"$f / f_{\rm cut}$", fontsize=13)
    ax2.set_ylabel(r"$v_g / c_0$", fontsize=13)
    ax2.set_title("Group velocity (subluminal, causal)", fontsize=14)
    ax2.legend(fontsize=9)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.set_ylim(0.0, 1.05)

    fig.suptitle("TGP gravitational wave velocities (breathing mode)", fontsize=15, y=1.02)
    fig.tight_layout()
    path = os.path.join(save_dir, "gw_velocities.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3: f_cut as function of gamma
# ═══════════════════════════════════════════════════════════════════════════
def plot_cutoff_vs_gamma(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    beta_vals = [0.00, 0.01, 0.05, 0.10]
    colors = plt.cm.cool(np.linspace(0.1, 0.9, len(beta_vals)))

    gamma_arr = np.linspace(0.001, 0.5, 500)

    for bh, col in zip(beta_vals, colors):
        f_arr = []
        g_valid = []
        for gh in gamma_arr:
            ms2 = m_sp_squared(bh, gh)
            if ms2 > 0:
                m_phys = np.sqrt(ms2) / r0
                f_arr.append(c0 * m_phys / (2.0 * np.pi))
                g_valid.append(gh)
        if len(f_arr) > 0:
            ax.plot(g_valid, f_arr, color=col, lw=1.8,
                    label=rf"$\hat\beta = {bh}$")

    # Reference frequency bands
    ax.axhspan(10, 1e4, color="orange", alpha=0.08, label="LIGO band (10-10$^4$ Hz)")
    ax.axhspan(1e-4, 1e-1, color="green", alpha=0.08, label="LISA band (0.1-100 mHz)")
    ax.axhspan(1e-9, 1e-7, color="blue", alpha=0.08, label="PTA band (1-100 nHz)")

    ax.set_xlabel(r"$\hat\gamma$", fontsize=13)
    ax.set_ylabel(r"$f_{\rm cut}$ [Hz]", fontsize=13)
    ax.set_yscale("log")
    ax.set_title(r"Low-frequency cutoff $f_{\rm cut} = c_0 m_{\rm sp} / (2\pi)$",
                 fontsize=14)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    path = os.path.join(save_dir, "gw_cutoff_vs_gamma.png")
    fig.savefig(path, dpi=180)
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 4: GW170817 constraint region
# ═══════════════════════════════════════════════════════════════════════════
def plot_gw170817_constraint(save_dir=None):
    """
    GW170817 measured |c_GW - c0|/c0 < 1e-15.  In TGP the group velocity is

        v_g / c0 = 1 / sqrt(1 + (a m_sp / k)^2)

    so  |v_g - c0|/c0 ~ (1/2)(a m_sp / k)^2  for small mass.

    At LIGO frequency f ~ 100 Hz:  k = 2 pi f / c0
    Constraint:  (a m_sp / k)^2 < 2e-15

    This bounds m_sp_phys, and therefore (3 gamma_hat - 2 beta_hat) / r0^2.
    """
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 7))

    # GW170817 parameters
    f_LIGO = 100.0       # Hz (characteristic frequency)
    k_LIGO = 2.0 * np.pi * f_LIGO / c0   # 1/m
    delta_v_bound = 1e-15  # |v_g - c0| / c0

    # Constraint: (m_sp_phys / k_LIGO)^2 < 2 * delta_v_bound  (a=1 today)
    m_sp_max = k_LIGO * np.sqrt(2.0 * delta_v_bound)

    # In terms of hat parameters: m_sp_phys^2 = (3 gamma_hat - 2 beta_hat) / r0^2
    # So:  3 gamma_hat - 2 beta_hat < m_sp_max^2 * r0^2
    constraint_rhs = m_sp_max**2 * r0**2

    print(f"  GW170817 constraint:")
    print(f"    f_LIGO            = {f_LIGO} Hz")
    print(f"    k_LIGO            = {k_LIGO:.6e} m^-1")
    print(f"    m_sp_max          = {m_sp_max:.6e} m^-1")
    print(f"    m_sp_max * r0     = {m_sp_max * r0:.6e}")
    print(f"    3 gamma - 2 beta  < {constraint_rhs:.6e}  (in hat units)")

    # Parameter space scan
    beta_arr = np.linspace(0, 0.5, 400)
    gamma_arr = np.linspace(0, 0.5, 400)
    B, G = np.meshgrid(beta_arr, gamma_arr)

    # m_sp^2 = 3 gamma - 2 beta (must be positive for propagating modes)
    MS2 = 3.0 * G - 2.0 * B
    M_PHYS_SQ = MS2 / r0**2

    # Group velocity deviation: |v_g - c0|/c0 ~ (1/2)(m_sp_phys / k_LIGO)^2
    # for m_sp << k
    delta_v = np.where(MS2 > 0, 0.5 * M_PHYS_SQ / k_LIGO**2, 0.0)

    # Plot: log10 of the velocity deviation
    delta_v_plot = np.where(delta_v > 0, np.log10(delta_v), -20)

    im = ax.pcolormesh(beta_arr, gamma_arr, delta_v_plot,
                       cmap="RdYlGn_r", shading="auto",
                       vmin=-20, vmax=-10)
    cbar = fig.colorbar(im, ax=ax, label=r"$\log_{10}\,|\Delta v_g / c_0|$")

    # GW170817 exclusion boundary: 3 gamma - 2 beta = constraint_rhs
    # gamma = (constraint_rhs + 2 beta) / 3
    beta_line = np.linspace(0, 0.5, 200)
    gamma_line = (constraint_rhs + 2.0 * beta_line) / 3.0
    ax.plot(beta_line, gamma_line, "w-", lw=2.5,
            label=rf"GW170817 bound ($|\Delta v|/c_0 < 10^{{-15}}$)")

    # Stability boundary: m_sp^2 = 0 => gamma = 2 beta / 3
    gamma_stab = 2.0 * beta_line / 3.0
    ax.plot(beta_line, gamma_stab, "k--", lw=1.5,
            label=r"Stability: $m_{\rm sp}^2 = 0$")

    # Shade excluded region (above the GW170817 line AND above stability)
    # "Excluded" means the deviation is TOO LARGE, i.e. above the constraint line
    ax.fill_between(beta_line, gamma_line, 0.5,
                    where=(gamma_line < 0.5),
                    color="red", alpha=0.15, label="Excluded by GW170817")

    # Shade tachyonic region (below stability line)
    ax.fill_between(beta_line, 0, gamma_stab,
                    color="gray", alpha=0.2, label=r"Tachyonic ($m_{\rm sp}^2 < 0$)")

    ax.set_xlabel(r"$\hat\beta$", fontsize=13)
    ax.set_ylabel(r"$\hat\gamma$", fontsize=13)
    ax.set_title("GW170817 constraint on TGP parameter space\n"
                 "(nuclear-scale hatted parameters)", fontsize=14)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.5)

    fig.tight_layout()
    path = os.path.join(save_dir, "gw_gw170817.png")
    fig.savefig(path, dpi=180)
    print(f"  Saved {path}")
    plt.close(fig)

    return constraint_rhs


# ═══════════════════════════════════════════════════════════════════════════
# Plot 5: Disformal tensor mode contribution
# ═══════════════════════════════════════════════════════════════════════════
def plot_disformal_tensor(save_dir=None):
    """
    Plot the ratio h_breathing / h_tensor as function of disformal coupling B/M*^4.

    Key physics (sek08, Hypothesis hyp:disformal):
    - The conformal metric g = Phi * eta gives ONLY breathing mode (proven).
    - The disformal extension g = A(Phi)eta + B(Phi)(dPhi)(dPhi)/M*^4
      generates tensor modes h_+, h_x from anisotropic sources.
    - For B/M*^4 >> 1: tensor modes dominate (GR-like polarization pattern).
    - For B/M*^4 << 1: breathing mode dominates (pure TGP prediction).
    - For B/M*^4 ~ 1: mixed polarization (testable with detector networks).
    """
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # --- Left panel: ratio h_breathing / h_tensor vs B/M*^4 ---
    B_values = np.logspace(-4, 4, 500)
    anisotropy_values = [1.0, 0.5, 0.1, 0.01]
    colors = plt.cm.magma(np.linspace(0.2, 0.8, len(anisotropy_values)))

    for Q, col in zip(anisotropy_values, colors):
        ratio = 1.0 / (B_values * Q)
        ax1.plot(B_values, ratio, color=col, lw=1.8,
                 label=rf"$Q_{{\rm aniso}} = {Q}$")

    ax1.axhline(1.0, color="gray", ls="--", lw=1.0, alpha=0.7)
    ax1.annotate("Equal amplitudes", xy=(1e-3, 1.2), fontsize=10,
                 color="gray", style="italic")

    # Shade regions
    ax1.axhspan(1, 1e8, color="blue", alpha=0.05)
    ax1.annotate("Breathing dominates\n(pure TGP)", xy=(1e-3, 1e3),
                 fontsize=10, color="blue", alpha=0.7, ha="center")
    ax1.axhspan(1e-8, 1, color="red", alpha=0.05)
    ax1.annotate("Tensor dominates\n(GR-like)", xy=(1e2, 1e-3),
                 fontsize=10, color="red", alpha=0.7, ha="center")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$B / M_*^4$ (disformal coupling)", fontsize=13)
    ax1.set_ylabel(r"$h_{\rm breathing} / h_{\rm tensor}$", fontsize=13)
    ax1.set_title("Polarization ratio vs disformal coupling", fontsize=14)
    ax1.legend(fontsize=9)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.set_ylim(1e-5, 1e5)

    # --- Right panel: effective polarization content vs B/M*^4 ---
    Q = 1.0  # maximally anisotropic source (binary)
    h_b = np.ones_like(B_values)  # normalized breathing = 1
    h_t = disformal_tensor_amplitude(h_b, B_values, Q)

    # Fractional power in each mode
    total_sq = h_b**2 + h_t**2
    frac_breathing = h_b**2 / total_sq
    frac_tensor = h_t**2 / total_sq

    ax2.fill_between(B_values, 0, frac_tensor, color="red", alpha=0.3,
                     label=r"Tensor ($h_+, h_\times$)")
    ax2.fill_between(B_values, frac_tensor, 1.0, color="blue", alpha=0.3,
                     label="Breathing (scalar)")
    ax2.plot(B_values, frac_breathing, "b-", lw=1.5)
    ax2.plot(B_values, frac_tensor, "r-", lw=1.5)

    ax2.axvline(1.0, color="gray", ls="--", lw=1.0)
    ax2.annotate(r"$B/M_*^4 = 1$", xy=(1.2, 0.5), fontsize=10, color="gray")

    ax2.set_xscale("log")
    ax2.set_xlabel(r"$B / M_*^4$ (disformal coupling)", fontsize=13)
    ax2.set_ylabel("Fractional GW power", fontsize=13)
    ax2.set_title(r"Polarization content ($Q_{\rm aniso} = 1$, binary source)",
                  fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.set_ylim(0, 1)

    fig.suptitle("TGP: Disformal coupling and tensor mode generation\n"
                 r"($g_{\mu\nu} = A(\Phi)\eta_{\mu\nu} + "
                 r"B(\Phi)\,\partial_\mu\Phi\,\partial_\nu\Phi / M_*^4$)",
                 fontsize=14, y=1.05)
    fig.tight_layout()
    path = os.path.join(save_dir, "gw_disformal_tensor.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 6: Cosmological-scale analysis and honest assessment
# ═══════════════════════════════════════════════════════════════════════════
def plot_cosmo_assessment(save_dir=None):
    """
    Compare nuclear-scale vs cosmological-scale predictions.

    Key result: with natural cosmological parameters, the TGP breathing mode
    is indistinguishable from GR at all current detector frequencies.
    """
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    m_sp_cosmo, gamma_cosmo, f_cut_cosmo = m_sp_cosmological()

    # --- Left panel: velocity deviation at different detector frequencies ---
    # For cosmological m_sp
    f_detectors = {
        "PTA (nHz)": np.logspace(-9, -7, 100),
        "LISA (mHz)": np.logspace(-4, -1, 100),
        "LIGO (Hz)": np.logspace(1, 4, 100),
        "Einstein Tel. (Hz)": np.logspace(0, 4, 100),
    }
    detector_colors = ["blue", "green", "orange", "red"]

    for (name, f_arr), col in zip(f_detectors.items(), detector_colors):
        k_arr = 2.0 * np.pi * f_arr / c0
        delta_v = 0.5 * (m_sp_cosmo / k_arr)**2
        ax1.plot(f_arr, delta_v, color=col, lw=1.8, label=name)

    # GW170817 bound
    ax1.axhline(1e-15, color="k", ls="--", lw=1.5, label=r"GW170817: $|\Delta v|/c_0 < 10^{-15}$")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Frequency [Hz]", fontsize=13)
    ax1.set_ylabel(r"$|\Delta v_g| / c_0$", fontsize=13)
    ax1.set_title(f"Velocity deviation with cosmological $m_{{\\rm sp}}$\n"
                  f"($m_{{\\rm sp}} = {m_sp_cosmo:.2e}$ m$^{{-1}}$, "
                  f"$f_{{\\rm cut}} = {f_cut_cosmo:.2e}$ Hz)",
                  fontsize=12)
    ax1.legend(fontsize=9, loc="upper right")
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.set_ylim(1e-55, 1e-5)

    # --- Right panel: comparison of nuclear vs cosmological predictions ---
    f_range = np.logspace(-20, 5, 1000)

    # Nuclear scale
    m_nuc = m_sp_physical(BETA_HAT_DEFAULT, GAMMA_HAT_DEFAULT)
    k_nuc = 2.0 * np.pi * f_range / c0
    delta_v_nuc = 0.5 * (m_nuc / k_nuc)**2

    # Cosmological scale
    delta_v_cosmo = 0.5 * (m_sp_cosmo / k_nuc)**2

    ax2.plot(f_range, delta_v_nuc, "r-", lw=1.8,
             label=f"Nuclear scale ($m_{{\\rm sp}} = {m_nuc:.1e}$ m$^{{-1}}$)")
    ax2.plot(f_range, delta_v_cosmo, "b-", lw=1.8,
             label=f"Cosmological ($m_{{\\rm sp}} = {m_sp_cosmo:.1e}$ m$^{{-1}}$)")

    ax2.axhline(1e-15, color="k", ls="--", lw=1.5, alpha=0.7,
                label="GW170817 bound")

    # Shade detector bands
    ax2.axvspan(10, 1e4, color="orange", alpha=0.08, label="LIGO")
    ax2.axvspan(1e-4, 1e-1, color="green", alpha=0.08, label="LISA")
    ax2.axvspan(1e-9, 1e-7, color="blue", alpha=0.08, label="PTA")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Frequency [Hz]", fontsize=13)
    ax2.set_ylabel(r"$|\Delta v_g| / c_0$", fontsize=13)
    ax2.set_title("Nuclear vs cosmological scale predictions", fontsize=14)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.set_ylim(1e-55, 1e5)
    ax2.set_xlim(1e-20, 1e5)

    fig.suptitle("TGP honest assessment: cosmological parameters predict\n"
                 "breathing mode indistinguishable from GR at detector frequencies",
                 fontsize=13, y=1.05)
    fig.tight_layout()
    path = os.path.join(save_dir, "gw_cosmo_assessment.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)

    return m_sp_cosmo, gamma_cosmo, f_cut_cosmo


# ═══════════════════════════════════════════════════════════════════════════
# Summary of observational constraints
# ═══════════════════════════════════════════════════════════════════════════
def print_summary(constraint_rhs):
    print("\n" + "=" * 72)
    print("  TGP Gravitational Wave Properties -- Summary")
    print("=" * 72)

    print(f"\n  Physical constants:")
    print(f"    c0   = {c0:.3e} m/s")
    print(f"    G0   = {G0:.3e} m^3/(kg s^2)")
    print(f"    r0   = {r0:.3e} m  (nuclear scale)")
    print(f"    H0   = {H0:.3e} s^-1  (Hubble constant)")
    print(f"    Phi0 = {Phi0}  (background field)")

    print(f"\n  Key formulae (CORRECTED):")
    print(f"    m_sp^2 = 3 gamma - 2 beta       [NOT divided by Phi0]")
    print(f"    omega^2 = c0^2 (k^2/a^2 + m_sp^2)")
    print(f"    v_phi = c0 sqrt(1 + (a m_sp/k)^2)    [superluminal]")
    print(f"    v_g   = c0 / sqrt(1 + (a m_sp/k)^2)  [subluminal, causal]")
    print(f"    f_cut = c0 m_sp / (2 pi)")

    print(f"\n  Polarization (sek08, Thm thm:no-tensor + Hyp hyp:disformal):")
    print(f"    - Conformal metric g = Phi*eta: breathing mode ONLY (proven)")
    print(f"      (Weyl tensor vanishes for conformally flat metrics)")
    print(f"    - Breathing mode is the PRIMARY/fundamental GW polarization")
    print(f"    - Disformal extension g = A*eta + B*(dPhi)(dPhi)/M*^4:")
    print(f"      CAN generate tensor modes h_+, h_x from source anisotropy")
    print(f"    - Ratio: h_breathing/h_tensor ~ 1/(B/M*^4 * Q_aniso)")
    print(f"    - For observed GW events (tensor-dominated), need B/M*^4 >> 1")

    print(f"\n  GW170817 constraint (nuclear-scale parameters):")
    print(f"    |c_GW - c0| / c0 < 1e-15")
    print(f"    => 3 gamma_hat - 2 beta_hat < {constraint_rhs:.3e}")
    print(f"    For beta = gamma:  gamma_hat < {constraint_rhs:.3e}")
    print(f"    => m_sp < {np.sqrt(constraint_rhs) / r0:.3e} m^-1")
    print(f"    => Compton wavelength lambda_sp > "
          f"{2 * np.pi * r0 / np.sqrt(constraint_rhs):.3e} m")

    # Cosmological analysis
    m_sp_cosmo, gamma_cosmo, f_cut_cosmo = m_sp_cosmological()
    print(f"\n  Cosmological-scale analysis:")
    print(f"    gamma_cosmo = Phi0 * H0^2 / c0^2 = {gamma_cosmo:.3e}")
    print(f"    m_sp_cosmo  = sqrt(Phi0) * H0/c0 = {m_sp_cosmo:.3e} m^-1")
    print(f"    f_cut_cosmo = {f_cut_cosmo:.3e} Hz")
    k_LIGO = 2.0 * np.pi * 100.0 / c0
    delta_v_cosmo = 0.5 * (m_sp_cosmo / k_LIGO)**2
    print(f"    delta_v at 100 Hz (LIGO)          = {delta_v_cosmo:.3e}")
    print(f"    (trivially satisfies GW170817: {delta_v_cosmo:.1e} << 1e-15)")

    # Example: default parameters (nuclear scale)
    bh, gh = BETA_HAT_DEFAULT, GAMMA_HAT_DEFAULT
    ms2 = m_sp_squared(bh, gh)
    print(f"\n  Example (nuclear-scale): beta_hat={bh}, gamma_hat={gh}")
    print(f"    m_sp^2 (hat)     = {ms2:.6f}")
    print(f"    m_sp (physical)  = {m_sp_physical(bh, gh):.6e} m^-1")
    print(f"    f_cut            = {f_cutoff(bh, gh):.6e} Hz")
    if ms2 > 0:
        lam_sp = 2 * np.pi / m_sp_physical(bh, gh)
        print(f"    Compton lambda   = {lam_sp:.6e} m")

    print(f"\n" + "-" * 72)
    print(f"  HONEST ASSESSMENT")
    print(f"-" * 72)
    print(f"    1. With natural (cosmological) parameters, the TGP breathing")
    print(f"       mode is INDISTINGUISHABLE from GR at LIGO/Virgo/KAGRA")
    print(f"       frequencies. The velocity deviation is ~ 1e-44, far below")
    print(f"       any conceivable measurement precision.")
    print(f"    2. The breathing mode becomes detectable ONLY if there is an")
    print(f"       additional mechanism that enhances m_sp beyond the")
    print(f"       cosmological Hubble scale (m_sp >> H_0/c_0 ~ 1e-26 m^-1).")
    print(f"    3. The disformal coupling (B/M*^4 term) is REQUIRED to")
    print(f"       reproduce the observed tensor polarization pattern in")
    print(f"       GW events like GW150914 and GW170817.")
    print(f"    4. TGP's unique prediction: an additional breathing mode")
    print(f"       component on top of tensor modes. Detectable with 3+")
    print(f"       detector networks via polarization decomposition.")
    print(f"    5. The nuclear-scale analysis (r_0 ~ 1e-15 m) gives extreme")
    print(f"       values because hatted parameters are rescaled by r_0^2.")
    print(f"       This is NOT the physically relevant regime for GW astronomy.")
    print("=" * 72)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 72)
    print("TGP gravitational wave breathing mode analysis")
    print("  (with disformal tensor modes and cosmological assessment)")
    print("=" * 72)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')

    print("\n  [1/6] Dispersion relation ...")
    plot_dispersion(save_dir=save_dir)

    print("\n  [2/6] Phase and group velocities ...")
    plot_velocities(save_dir=save_dir)

    print("\n  [3/6] Cutoff frequency vs gamma ...")
    plot_cutoff_vs_gamma(save_dir=save_dir)

    print("\n  [4/6] GW170817 constraint region ...")
    constraint_rhs = plot_gw170817_constraint(save_dir=save_dir)

    print("\n  [5/6] Disformal tensor mode analysis ...")
    plot_disformal_tensor(save_dir=save_dir)

    print("\n  [6/6] Cosmological-scale assessment ...")
    plot_cosmo_assessment(save_dir=save_dir)

    print_summary(constraint_rhs)

    print("\nDone.")


if __name__ == "__main__":
    main()
