"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
w_de_redshift.py  --  Theory of Generated Space (TGP)
=====================================================
Compute w_DE(z) from structure formation history in TGP.

Physics summary
---------------
TGP gives Lambda_eff = gamma / 56, with Phi_0 ~ 115 matching
Lambda_obs.  The dark energy equation of state w_DE deviates from -1
through three mechanisms:

  (1) Spatial variation:  Phi = Phi_0 (1 + delta_Phi), where delta_Phi
      tracks matter perturbations via the Newtonian potential Psi_N.
      This gives Lambda_eff(x) = gamma Phi(x) / 56.

  (2) Time evolution:  as structures grow, <delta_Phi^2> ~ D^2(z) <Psi^2>
      increases, so the volume-averaged Lambda_eff evolves.

  (3) Backreaction:  the nonlinear coupling beta Phi^2 - gamma Phi^3
      modifies the effective dark energy density and pressure.

The result is w_DE(z) = -1 + epsilon(z), where epsilon ~ 5.4e-9 from
lambda_eff_quantitative.py.  This is HONEST: the deviation is far below
any foreseeable observational threshold (current: ~0.03; future: ~0.01).

Outputs
-------
  scripts/plots/w_de_redshift.png   -- 4-panel diagnostic plot

Refs: lambda_eff_quantitative.py, growth_factor_tgp.py, w_de_exact.py
"""

import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================================
# Physical constants (SI)
# =====================================================================
c0 = 2.998e8           # m/s
G0 = 6.674e-11         # m^3 kg^-1 s^-2
H0_km = 67.4           # km/s/Mpc
Mpc_m = 3.0857e22      # m/Mpc
H0 = H0_km * 1e3 / Mpc_m   # s^-1  ~2.184e-18

Omega_m = 0.315
Omega_Lambda = 0.685
Lambda_obs = 1.11e-52   # m^-2  (observed)

# TGP parameters
Phi0 = 168.0 * Omega_Lambda   # = 115.1 (from P(1) = gamma/56)
gamma_tgp = Phi0 * H0**2 / c0**2   # m^-2

# Structure correction amplitude from lambda_eff_quantitative.py
# <delta_Phi^2 / Phi_0^2> dominated by cosmic-web Bardeen potential
sigma_Psi = 3e-5        # RMS Bardeen potential on large scales
psi2_0 = sigma_Psi**2   # = 9e-10
epsilon_0 = 6.0 * psi2_0   # delta_Lambda / Lambda at z=0 = 5.4e-9


# =====================================================================
# (1) LCDM growth factor D(z) via ODE
# =====================================================================
def E_squared(z):
    """E^2(z) = H^2(z) / H_0^2 for flat LCDM."""
    return Omega_m * (1 + z)**3 + Omega_Lambda


def growth_ode(lna, y):
    """
    Growth equation for delta_m in terms of ln(a):
        D'' + [2 + d ln H / d ln a] D' = (3/2) Omega_m(a) D
    State: y = [D, dD/d(ln a)]
    """
    a = np.exp(lna)
    z = 1.0 / a - 1.0
    D, Dp = y

    E2 = E_squared(z)
    Omega_m_a = Omega_m * (1 + z)**3 / E2

    # d ln H / d ln a = (1/2) d ln E^2 / d ln a
    # d E^2 / d ln a = d E^2 / da * a = -3 Omega_m (1+z)^3 * (-1/a^2) * a
    #                = -3 Omega_m / a^3 * a = -3 Omega_m (1+z)^3
    # Wait, let's be careful:
    # E^2 = Omega_m a^{-3} + Omega_Lambda
    # d E^2 / d ln a = -3 Omega_m a^{-3}
    dlnE2_dlna = -3.0 * Omega_m * (1 + z)**3 / E2
    dlnH_dlna = 0.5 * dlnE2_dlna

    Dpp = 1.5 * Omega_m_a * D - (2.0 + dlnH_dlna) * Dp
    return [Dp, Dpp]


def compute_growth_factor(z_arr):
    """
    Compute the linear growth factor D(z)/D(0) for flat LCDM.
    Returns D normalized to D(0)=1.
    """
    # Integrate from high z (matter domination: D ~ a) to z=0
    a_ini = 1e-4
    lna_span = (np.log(a_ini), 0.0)
    y0 = [a_ini, a_ini]  # D = a, dD/d(ln a) = a in matter era

    sol = solve_ivp(growth_ode, lna_span, y0, method='RK45',
                    rtol=1e-10, atol=1e-13,
                    dense_output=True)

    if not sol.success:
        raise RuntimeError(f"Growth factor ODE failed: {sol.message}")

    # Evaluate at requested redshifts
    a_eval = 1.0 / (1.0 + z_arr)
    lna_eval = np.log(a_eval)

    D_arr = np.zeros_like(z_arr)
    for i, lna in enumerate(lna_eval):
        state = sol.sol(lna)
        D_arr[i] = state[0]

    # Normalize to D(z=0) = 1
    D0 = sol.sol(0.0)[0]
    D_arr /= D0

    return D_arr


def compute_growth_rate(z_arr, D_arr):
    """
    Compute f(z) = d ln D / d ln a = -(1+z) dD/dz / D.
    Uses numerical differentiation.
    """
    # Use the approximation f ~ Omega_m(z)^{0.55} (Linder 2005)
    # This is sufficiently accurate for our purposes
    E2 = E_squared(z_arr)
    Omega_m_z = Omega_m * (1 + z_arr)**3 / E2
    f_arr = Omega_m_z**0.55
    return f_arr


# =====================================================================
# (2) Phi perturbation sourced by matter growth
# =====================================================================
def compute_delta_phi_variance(z_arr, D_arr):
    """
    In TGP, delta_Phi / Phi_0 ~ Psi_N (Newtonian potential).
    The large-scale Bardeen potential is approximately constant in
    matter domination and decays during Lambda domination:
        Psi(z) ~ Psi_0 * T(z)
    where T(z) accounts for the ISW effect.

    For volume-averaged <delta_Phi^2>:
        <delta_Phi^2>(z) / Phi_0^2 ~ <Psi^2>(z)

    On linear scales, <Psi^2> ~ const during matter domination
    (Psi = const for growing mode), but the structure formation
    increases delta_rho, which couples to delta_Phi through:
        delta_Phi / Phi_0 = Psi_N = -(3/2)(H_0/k)^2 Omega_m (1+z) delta_m D(z)

    The variance contribution from structure is:
        <delta_Phi^2> / Phi_0^2 ~ sigma_Psi^2 * [D(z)/D(0)]^2 * g(z)^2

    where g(z) is the Bardeen potential growth suppression:
        g(z) = D(z) * (1+z) * H_0^2 / H^2(z) * (5 Omega_m / 2)
        normalized so g(0) = 1 approximately.

    For simplicity and honesty, we use the dominant scaling:
        <delta_Phi^2>(z) / Phi_0^2 = psi2_0 * [D(z)]^2
    This captures the main effect: delta_Phi grows with structure.
    """
    psi2_z = psi2_0 * D_arr**2
    return psi2_z


# =====================================================================
# (3) Effective dark energy from Phi evolution
# =====================================================================
def compute_w_de(z_arr, D_arr, f_arr, psi2_z):
    """
    Compute w_DE(z) from the TGP Phi field.

    The effective cosmological constant in TGP is:
        Lambda_eff(z) = (gamma/56) [1 + 6 <delta_Phi^2>/Phi_0^2]
                      = Lambda_0 [1 + epsilon(z)]

    where epsilon(z) = 6 psi2_z = epsilon_0 * D^2(z).

    For a time-varying Lambda, the effective equation of state is:
        w_DE(z) = -1 - (1/3) d ln rho_DE / d ln(1+z)

    With rho_DE ~ Lambda_eff(z):
        d ln Lambda_eff / d ln(1+z) = d ln(1 + epsilon D^2) / d ln(1+z)
                                     ~ epsilon * d(D^2) / d ln(1+z)   (for epsilon << 1)
                                     = epsilon * 2 D * dD/d ln(1+z)
                                     = -epsilon * 2 D * dD/d ln a
                                     = -2 epsilon_0 D^2 * f(z)

    where f(z) = d ln D / d ln a is the growth rate.

    Therefore:
        w_DE(z) = -1 - (1/3) * (-2 epsilon_0 D^2(z) f(z))
                = -1 + (2/3) epsilon_0 D^2(z) f(z)

    At z=0: w_0 = -1 + (2/3) epsilon_0 f(0)
           with f(0) = Omega_m^{0.55} ~ 0.547
           w_0 - (-1) ~ (2/3) * 5.4e-9 * 0.547 ~ 1.97e-9
    """
    epsilon_z = epsilon_0 * D_arr**2

    # w_DE(z) = -1 + (2/3) * epsilon_0 * D^2(z) * f(z)
    deviation = (2.0 / 3.0) * epsilon_0 * D_arr**2 * f_arr
    w_de = -1.0 + deviation

    return w_de, deviation


# =====================================================================
# (4) CPL parametrization fit: w(z) = w_0 + w_a * z/(1+z)
# =====================================================================
def fit_cpl(z_arr, w_arr):
    """
    Fit w(z) = w_0 + w_a * z/(1+z) to the computed w_DE(z).
    Returns (w_0, w_a).
    """
    mask = (z_arr > 0.01) & (z_arr < 2.5) & np.isfinite(w_arr)
    z_fit = z_arr[mask]
    w_fit = w_arr[mask]

    def cpl_model(z, w0, wa):
        return w0 + wa * z / (1.0 + z)

    popt, pcov = curve_fit(cpl_model, z_fit, w_fit, p0=[-1.0, 0.0])
    return popt[0], popt[1]


# =====================================================================
# (5) delta_Lambda/Lambda evolution
# =====================================================================
def compute_delta_lambda(z_arr, D_arr):
    """
    delta_Lambda/Lambda(z) = 6 * <delta_Phi^2>(z) / Phi_0^2
                            = epsilon_0 * D^2(z)
    """
    return epsilon_0 * D_arr**2


# =====================================================================
# (6) 4-panel plot
# =====================================================================
def make_plot(z_arr, w_de, epsilon_z, w0_cpl, wa_cpl, delta_lambda_z,
              D_arr, f_arr, save_path):
    """Generate 4-panel diagnostic plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # ── Panel (a): w_DE(z) vs z ──────────────────────────────────────
    ax = axes[0, 0]

    ax.plot(z_arr, w_de, 'b-', lw=2.5, label=r'TGP: $w_{\rm DE}(z)$')
    ax.axhline(-1.0, color='k', ls='--', lw=1.2,
               label=r'$\Lambda$CDM: $w = -1$')

    # Observational bands (schematic)
    z_band = np.linspace(0, 3, 200)
    sigma_current = 0.05 * np.ones_like(z_band)  # current ~5%
    sigma_future = 0.01 * np.ones_like(z_band)    # future ~1%
    ax.fill_between(z_band, -1 - sigma_current, -1 + sigma_current,
                    color='orange', alpha=0.15,
                    label=r'Current bounds ($\pm 0.05$)')
    ax.fill_between(z_band, -1 - sigma_future, -1 + sigma_future,
                    color='blue', alpha=0.10,
                    label=r'Future reach ($\pm 0.01$)')

    ax.set_xlabel(r'Redshift $z$', fontsize=12)
    ax.set_ylabel(r'$w_{\rm DE}(z)$', fontsize=12)
    ax.set_title(r'(a) Dark energy equation of state $w_{\rm DE}(z)$',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(0, 3)
    ax.set_ylim(-1.08, -0.92)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, ls=':', alpha=0.4)

    # Annotate: TGP curve is indistinguishable from w=-1
    ax.text(0.03, 0.05,
            r'TGP prediction: $w_{\rm DE} = -1 + \mathcal{O}(10^{-9})$'
            '\nIndistinguishable from $\\Lambda$CDM\nat this scale',
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9),
            verticalalignment='bottom')

    # ── Panel (b): epsilon(z) = w + 1 on log scale ──────────────────
    ax = axes[0, 1]

    ax.semilogy(z_arr[z_arr > 0.001], epsilon_z[z_arr > 0.001],
                'b-', lw=2.5, label=r'$\varepsilon(z) = w_{\rm DE}(z) + 1$')

    # Detection thresholds
    ax.axhline(0.05, color='orange', ls='--', lw=1.5,
               label=r'Current: DESI+Planck+SNe ($\sim 0.05$)')
    ax.axhline(0.01, color='green', ls='--', lw=1.5,
               label=r'Future: Euclid+LSST ($\sim 0.01$)')
    ax.axhline(1e-3, color='purple', ls=':', lw=1.2,
               label=r'Optimistic future ($\sim 10^{-3}$)')

    ax.set_xlabel(r'Redshift $z$', fontsize=12)
    ax.set_ylabel(r'$\varepsilon(z) = w_{\rm DE} + 1$', fontsize=12)
    ax.set_title(r'(b) Deviation from $w = -1$ (log scale)',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(0, 3)
    ax.set_ylim(1e-11, 1)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, ls=':', alpha=0.4)

    # Gap annotation
    gap_orders = int(np.log10(0.01 / epsilon_z[0]))
    ax.annotate(
        f'Gap: ~{gap_orders} orders\nof magnitude',
        xy=(0.5, epsilon_z[len(z_arr)//6]),
        xytext=(1.5, 1e-5),
        fontsize=10, color='red',
        arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
        bbox=dict(boxstyle='round', fc='mistyrose', alpha=0.9))

    # ── Panel (c): CPL plane (w_0, w_a) with observational contours ──
    ax = axes[1, 0]

    # TGP prediction point
    ax.plot(w0_cpl, wa_cpl, 'r*', ms=20, zorder=10,
            label=rf'TGP ($w_0={w0_cpl:.2e},\;w_a={wa_cpl:.2e}$)')

    # Lambda CDM point
    ax.plot(-1.0, 0.0, 'ko', ms=10, zorder=10,
            label=r'$\Lambda$CDM ($w_0=-1, w_a=0$)')

    # Observational constraint ellipses (schematic)
    theta = np.linspace(0, 2 * np.pi, 300)

    # Planck 2018 + BAO (approximate 1-sigma, 2-sigma)
    for nsig, ls_ell, alpha_f in [(1, '-', 0.15), (2, '--', 0.08)]:
        ew0 = nsig * 0.08
        ewa = nsig * 0.36
        # Tilted ellipse (correlation ~ -0.5)
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        angle = -0.4  # rotation angle
        x_ell = -1.0 + ew0 * (np.cos(angle) * cos_t - np.sin(angle) * sin_t)
        y_ell = 0.0 + ewa * (np.sin(angle) * cos_t + np.cos(angle) * sin_t)
        ax.plot(x_ell, y_ell, color='gray', ls=ls_ell, lw=1.0)
    ax.fill(x_ell, y_ell, color='gray', alpha=0.08,
            label=r'Planck+BAO $2\sigma$ (schematic)')

    # DESI 2024 hint region (w0 ~ -0.7, wa ~ -1.0)
    ew0_d = 0.12
    ewa_d = 0.45
    x_desi = -0.75 + ew0_d * np.cos(theta)
    y_desi = -0.8 + ewa_d * np.sin(theta)
    ax.fill(x_desi, y_desi, color='red', alpha=0.08,
            label=r'DESI 2024 hint $1\sigma$ (schematic)')
    ax.plot(x_desi, y_desi, 'r-', lw=0.8)

    # Phantom divide and axes
    ax.axvline(-1.0, color='gray', ls=':', lw=0.5)
    ax.axhline(0.0, color='gray', ls=':', lw=0.5)

    ax.set_xlabel(r'$w_0$', fontsize=12)
    ax.set_ylabel(r'$w_a$', fontsize=12)
    ax.set_title(r'(c) CPL plane: $w(z) = w_0 + w_a\,z/(1+z)$',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(-1.5, -0.4)
    ax.set_ylim(-2.0, 1.0)
    ax.legend(fontsize=7, loc='upper left')
    ax.grid(True, ls=':', alpha=0.4)

    # Inset: zoom near Lambda CDM
    inset = ax.inset_axes([0.55, 0.55, 0.42, 0.42])
    inset.plot(w0_cpl, wa_cpl, 'r*', ms=15, zorder=10)
    inset.plot(-1.0, 0.0, 'ko', ms=8, zorder=10)

    zoom_range = 5e-8
    inset.set_xlim(-1.0 - zoom_range, -1.0 + zoom_range)
    inset.set_ylim(-zoom_range, zoom_range)
    inset.set_title('Zoom near $\\Lambda$CDM', fontsize=7)
    inset.ticklabel_format(style='sci', scilimits=(-2, 2))
    inset.tick_params(labelsize=6)
    inset.grid(True, ls=':', alpha=0.4)

    # ── Panel (d): delta_Lambda/Lambda(z) evolution ──────────────────
    ax = axes[1, 1]

    ax.semilogy(z_arr, delta_lambda_z, 'b-', lw=2.5,
                label=r'$\delta\Lambda/\Lambda = \varepsilon_0\,D^2(z)$')
    ax.axhline(epsilon_0, color='r', ls='--', lw=1.5,
               label=rf'$\varepsilon_0 = {epsilon_0:.1e}$ (today)')

    # Also plot D^2(z) scaled
    ax.semilogy(z_arr, D_arr**2 * epsilon_0 * 0.999, 'g:', lw=1.0,
                alpha=0.0)  # invisible, just for range

    ax.set_xlabel(r'Redshift $z$', fontsize=12)
    ax.set_ylabel(r'$\delta\Lambda / \Lambda$', fontsize=12)
    ax.set_title(r'(d) Fractional $\Lambda$ variation from structure',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(0, 3)
    ax.legend(fontsize=9)
    ax.grid(True, ls=':', alpha=0.4)

    # Growth factor on twin axis
    ax2 = ax.twinx()
    ax2.plot(z_arr, D_arr, 'g--', lw=1.5, alpha=0.6)
    ax2.set_ylabel(r'$D(z)/D(0)$', fontsize=10, color='green')
    ax2.tick_params(axis='y', labelcolor='green')

    ax.text(0.55, 0.35,
            r'$\delta\Lambda/\Lambda \sim 5\times10^{-9}$'
            '\nCompletely negligible\n'
            r'$\Lambda$ is effectively constant',
            transform=ax.transAxes, fontsize=9, ha='center',
            bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

    # ── Global ───────────────────────────────────────────────────────
    fig.suptitle(
        r'TGP: Dark energy equation of state $w_{\rm DE}(z)$ from structure formation',
        fontsize=15, y=1.01)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=180, bbox_inches='tight')
    print(f"\n  Plot saved to {save_path}")
    plt.close(fig)


# =====================================================================
# Main computation
# =====================================================================
def main():
    print()
    print("#" * 70)
    print("#  TGP: w_DE(z) from structure formation history")
    print("#" * 70)

    # ── Cosmological parameters ──────────────────────────────────────
    print("\n" + "=" * 70)
    print("  PARAMETERS")
    print("=" * 70)
    print(f"  H_0       = {H0:.4e} s^-1  ({H0_km} km/s/Mpc)")
    print(f"  Omega_m   = {Omega_m}")
    print(f"  Omega_L   = {Omega_Lambda}")
    print(f"  Lambda_obs= {Lambda_obs:.3e} m^-2")
    print(f"  Phi_0     = 168 * Omega_Lambda = {Phi0:.4f}")
    print(f"  gamma     = Phi_0 * H_0^2 / c_0^2 = {gamma_tgp:.4e} m^-2")
    print(f"  sigma_Psi = {sigma_Psi:.1e} (Bardeen potential RMS)")
    print(f"  epsilon_0 = 6 * sigma_Psi^2 = {epsilon_0:.3e}")

    # ── Compute growth factor D(z) ──────────────────────────────────
    print("\n" + "=" * 70)
    print("  GROWTH FACTOR D(z)")
    print("=" * 70)

    z_arr = np.linspace(0, 3, 2000)
    D_arr = compute_growth_factor(z_arr)
    f_arr = compute_growth_rate(z_arr, D_arr)

    print(f"\n  {'z':>6s}  {'D(z)/D(0)':>10s}  {'f(z)':>8s}")
    print("  " + "-" * 30)
    for z_val in [0, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]:
        idx = np.argmin(np.abs(z_arr - z_val))
        print(f"  {z_val:6.1f}  {D_arr[idx]:10.4f}  {f_arr[idx]:8.4f}")

    # ── Compute delta_Phi variance ──────────────────────────────────
    print("\n" + "=" * 70)
    print("  PHI PERTURBATION VARIANCE")
    print("=" * 70)

    psi2_z = compute_delta_phi_variance(z_arr, D_arr)

    print(f"\n  <delta_Phi^2/Phi_0^2>(z=0) = {psi2_z[0]:.3e}")
    print(f"  <delta_Phi^2/Phi_0^2>(z=1) = {psi2_z[np.argmin(np.abs(z_arr-1))]:.3e}")
    print(f"  <delta_Phi^2/Phi_0^2>(z=3) = {psi2_z[-1]:.3e}")

    # ── Compute w_DE(z) ─────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  DARK ENERGY EQUATION OF STATE")
    print("=" * 70)

    w_de, epsilon_z = compute_w_de(z_arr, D_arr, f_arr, psi2_z)

    print(f"\n  w_DE(z) = -1 + epsilon(z)")
    print(f"  epsilon(z) = (2/3) * epsilon_0 * D^2(z) * f(z)")
    print(f"\n  {'z':>6s}  {'w_DE':>15s}  {'epsilon = w+1':>15s}")
    print("  " + "-" * 40)
    for z_val in [0, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]:
        idx = np.argmin(np.abs(z_arr - z_val))
        print(f"  {z_val:6.1f}  {w_de[idx]:15.10f}  {epsilon_z[idx]:15.3e}")

    # ── CPL fit ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  CPL PARAMETRIZATION: w(z) = w_0 + w_a * z/(1+z)")
    print("=" * 70)

    w0_cpl, wa_cpl = fit_cpl(z_arr, w_de)

    print(f"\n  w_0 = {w0_cpl:.10f}")
    print(f"  w_a = {wa_cpl:.10e}")
    print(f"  w_0 - (-1) = {w0_cpl + 1:.3e}")
    print(f"  w_a - 0    = {wa_cpl:.3e}")

    # ── Comparison with observational bounds ────────────────────────
    print("\n" + "=" * 70)
    print("  COMPARISON WITH OBSERVATIONAL BOUNDS")
    print("=" * 70)

    sigma_w0_planck = 0.08      # Planck 2018 + BAO
    sigma_w0_desi = 0.05        # DESI DR1 2024
    sigma_w0_euclid = 0.01      # Euclid forecast
    sigma_w0_stage4 = 0.005     # Stage IV forecast

    dev_w0 = abs(w0_cpl + 1)

    print(f"\n  TGP prediction: |w_0 + 1| = {dev_w0:.3e}")
    print(f"\n  {'Experiment':30s}  {'sigma(w_0)':>12s}  {'Detectable?':>12s}")
    print("  " + "-" * 58)
    for name, sigma in [("Planck 2018 + BAO", sigma_w0_planck),
                        ("DESI DR1 (2024)", sigma_w0_desi),
                        ("Euclid (forecast)", sigma_w0_euclid),
                        ("Stage IV (forecast)", sigma_w0_stage4),
                        ("Ideal cosmic variance", 1e-3)]:
        detectable = "NO" if dev_w0 < sigma else "YES"
        ratio = dev_w0 / sigma
        print(f"  {name:30s}  {sigma:12.1e}  {detectable:>12s}  "
              f"(ratio = {ratio:.1e})")

    # ── delta_Lambda/Lambda ────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  DELTA_LAMBDA / LAMBDA EVOLUTION")
    print("=" * 70)

    delta_lambda_z = compute_delta_lambda(z_arr, D_arr)

    print(f"\n  delta_Lambda/Lambda(z=0) = {delta_lambda_z[0]:.3e}")
    print(f"  delta_Lambda/Lambda(z=1) = {delta_lambda_z[np.argmin(np.abs(z_arr-1))]:.3e}")
    print(f"  delta_Lambda/Lambda(z=3) = {delta_lambda_z[-1]:.3e}")

    # ── Summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  SUMMARY AND HONESTY STATEMENT")
    print("=" * 70)
    print(f"""
  TGP predicts w_DE(z) = -1 + epsilon(z), where:

    epsilon(z) = (2/3) * 6 * <Psi^2> * D^2(z) * f(z)

  Key results:
    w_0 = {w0_cpl:.10f}
    w_a = {wa_cpl:.3e}
    |w_0 + 1| = {dev_w0:.3e}

  HONEST ASSESSMENT:
  ------------------
  The TGP dark energy equation of state is INDISTINGUISHABLE from
  w = -1 (pure cosmological constant) by any current or foreseeable
  observation.

  The deviation |w_0 + 1| ~ {dev_w0:.0e} is:
    - {dev_w0/sigma_w0_desi:.0e}x below DESI sensitivity
    - {dev_w0/sigma_w0_euclid:.0e}x below Euclid forecast
    - {dev_w0/1e-3:.0e}x below cosmic variance limit

  This is NOT a failure of TGP -- it is a PREDICTION:
  TGP predicts that dark energy behaves as a true cosmological
  constant to extraordinary precision.  Any detected deviation
  from w = -1 at the percent level would FALSIFY TGP (or require
  new physics beyond the minimal beta=gamma sector).

  The structure formation correction delta_Lambda/Lambda ~ 5e-9
  is similarly negligible.  Lambda_eff is effectively constant
  across cosmic history.
""")

    # ── Generate plot ───────────────────────────────────────────────
    save_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots', "w_de_redshift.png")
    make_plot(z_arr, w_de, epsilon_z, w0_cpl, wa_cpl, delta_lambda_z,
              D_arr, f_arr, save_path)

    print("\nAll computations completed successfully.")


if __name__ == "__main__":
    main()
