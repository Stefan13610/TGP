"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
lambda_eff_estimation.py  --  Theory of Generated Space (TGP)
==============================================================
Estimate the effective cosmological constant Lambda_eff from TGP
parameters and compare with observed value.

Physics
-------
In TGP, the effective cosmological constant has two sources:

  1. Residual action potential at psi = 1:
       P(1) = beta/7 - gamma/8
     For beta = gamma:  P(1) = gamma/56

     NOTE: This uses the CORRECT action potential P(g) = (beta/7)g^7 - (gamma/8)g^8
     derived from P'(g) = K(g)*(gamma g^3 - beta g^2) with K=g^4.
     The old (incorrect) value was U(1) = gamma/12 from the field-equation potential.

  2. Backreaction from structure formation:
       Lambda_back = gamma * <(delta_Phi / Phi_0)^2>
                   = 4 * gamma * sigma_Phi^2
     where sigma_Phi^2 ~ 10^{-10} -- 10^{-9} is the variance
     of the gravitational potential.

The dominant contribution is (1), giving:
    Lambda_eff ~ gamma/56

The naturalness condition tau_0 ~ 1/H_0 implies:
    gamma ~ Phi_0 * H_0^2 / c_0^2

Therefore:
    Lambda_eff ~ Phi_0 * H_0^2 / (56 * c_0^2)
    Lambda_obs = 3 * Omega_Lambda * H_0^2 / c_0^2 ~ 2.1 * H_0^2 / c_0^2

The ratio Lambda_obs / Lambda_eff ~ 117 / Phi_0, so Phi_0 ~ 100-120
gives the correct order of magnitude WITHOUT fine-tuning.

Outputs (saved to tooling/scripts/plots/):
    lambda_eff_vs_Phi0.png        -- Lambda_eff / Lambda_obs vs Phi_0
    lambda_eff_components.png     -- Potential vs backreaction contributions
    lambda_eff_parameter_space.png -- Contour plot in (gamma, Phi_0) plane
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════════════════
# Physical constants (SI)
# ═══════════════════════════════════════════════════════════════════════════
c0 = 2.998e8          # m/s
G0 = 6.674e-11        # m^3 kg^-1 s^-2
H0_si = 2.27e-18      # s^-1  (H0 ~ 70 km/s/Mpc)
Omega_Lambda = 0.7     # observed dark energy fraction

# Derived
Lambda_obs = 3 * Omega_Lambda * H0_si**2 / c0**2  # m^-2
print(f"Lambda_obs = {Lambda_obs:.4e} m^-2")


# ═══════════════════════════════════════════════════════════════════════════
# TGP estimation functions
# ═══════════════════════════════════════════════════════════════════════════
def gamma_natural(Phi0):
    """Natural value of gamma from tau_0 ~ 1/H_0 condition.
    gamma = Phi_0 * H_0^2 / c_0^2
    """
    return Phi0 * H0_si**2 / c0**2


def Lambda_eff_potential(Phi0, beta_ratio=1.0):
    """Lambda_eff from residual action potential P(psi=1).

    P(g) = (beta/7)g^7 - (gamma/8)g^8  (correct action potential)
    P(1) = beta/7 - gamma/8
    For beta = beta_ratio * gamma:
        P(1) = gamma * (beta_ratio/7 - 1/8)

    Parameters
    ----------
    Phi0 : float or array, background field value
    beta_ratio : float, beta/gamma ratio (1.0 for vacuum condition)
    """
    gamma = gamma_natural(Phi0)
    return gamma * (beta_ratio / 7.0 - 1.0 / 8.0)


def Lambda_eff_backreaction(Phi0, sigma_Phi=3e-5):
    """Lambda_eff from structure formation backreaction.

    Lambda_back = 4 * gamma * sigma_Phi^2

    Parameters
    ----------
    Phi0 : float or array
    sigma_Phi : float, RMS of Newtonian potential / c^2
    """
    gamma = gamma_natural(Phi0)
    return 4.0 * gamma * sigma_Phi**2


def Lambda_eff_total(Phi0, beta_ratio=1.0, sigma_Phi=3e-5):
    """Total Lambda_eff = potential + backreaction."""
    return Lambda_eff_potential(Phi0, beta_ratio) + Lambda_eff_backreaction(Phi0, sigma_Phi)


def m_sp_from_gamma(gamma_val):
    """Scalar field mass from gamma (for beta=gamma).
    m_sp = sqrt(gamma) [m^-1]
    """
    return np.sqrt(gamma_val)


def compton_wavelength(gamma_val):
    """Compton wavelength in meters."""
    m = m_sp_from_gamma(gamma_val)
    return 2.0 * np.pi / m if m > 0 else np.inf


def hubble_radius():
    """Hubble radius c_0 / H_0 [m]."""
    return c0 / H0_si


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1: Lambda_eff / Lambda_obs vs Phi_0
# ═══════════════════════════════════════════════════════════════════════════
def plot_lambda_vs_Phi0(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 6))

    Phi0_arr = np.linspace(1, 250, 500)

    # Total Lambda_eff
    L_pot = Lambda_eff_potential(Phi0_arr)
    L_back = Lambda_eff_backreaction(Phi0_arr)
    L_tot = L_pot + L_back

    ratio = L_tot / Lambda_obs

    ax.plot(Phi0_arr, ratio, "b-", lw=2.5, label=r"$\Lambda_{\rm eff}/\Lambda_{\rm obs}$")
    ax.axhline(1.0, color="r", ls="--", lw=1.5, label=r"$\Lambda_{\rm eff} = \Lambda_{\rm obs}$")

    # Mark the crossing point (P(1) = gamma/56, so Phi0_match ~ 168*Omega_L)
    idx_cross = np.argmin(np.abs(ratio - 1.0))
    Phi0_match = Phi0_arr[idx_cross] if idx_cross > 0 else 115.0
    ax.axvline(Phi0_match, color="green", ls=":", lw=1.2, alpha=0.7)
    ax.annotate(rf"$\Phi_0 \approx {Phi0_match:.0f}$",
                xy=(Phi0_match, 1.0), xytext=(min(Phi0_match + 15, 85), 2.5),
                fontsize=11, color="green",
                arrowprops=dict(arrowstyle="->", color="green"))

    # Shade "natural" region
    ax.axvspan(50, 200, color="green", alpha=0.08, label=r"$\Phi_0 = \mathcal{O}(100)$ (natural)")

    ax.set_xlabel(r"$\Phi_0$ (background field value)", fontsize=13)
    ax.set_ylabel(r"$\Lambda_{\rm eff} / \Lambda_{\rm obs}$", fontsize=13)
    ax.set_title(r"TGP prediction: $\Lambda_{\rm eff}$ vs observed $\Lambda$", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, ls=":", alpha=0.4)
    ax.set_xlim(1, 250)
    ax.set_ylim(0, 4)

    fig.tight_layout()
    path = os.path.join(save_dir, "lambda_eff_vs_Phi0.png")
    fig.savefig(path, dpi=180)
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2: Potential vs backreaction contributions
# ═══════════════════════════════════════════════════════════════════════════
def plot_components(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    Phi0_arr = np.logspace(0, 2, 200)

    # Component magnitudes
    L_pot = Lambda_eff_potential(Phi0_arr)
    sigma_vals = [1e-5, 3e-5, 1e-4, 3e-4]

    ax1.plot(Phi0_arr, L_pot / Lambda_obs, "b-", lw=2.5,
             label=r"Potential: $P(1) = \gamma/56$")

    colors = plt.cm.Oranges(np.linspace(0.3, 0.9, len(sigma_vals)))
    for sigma, col in zip(sigma_vals, colors):
        L_back = Lambda_eff_backreaction(Phi0_arr, sigma)
        ax1.plot(Phi0_arr, L_back / Lambda_obs, color=col, ls="--", lw=1.5,
                 label=rf"Backreaction: $\sigma_\Phi = {sigma:.0e}$")

    ax1.axhline(1.0, color="r", ls=":", lw=1, alpha=0.7)
    ax1.set_xlabel(r"$\Phi_0$", fontsize=13)
    ax1.set_ylabel(r"contribution / $\Lambda_{\rm obs}$", fontsize=13)
    ax1.set_title("Components of $\\Lambda_{\\rm eff}$", fontsize=14)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend(fontsize=9)
    ax1.grid(True, ls=":", alpha=0.4)

    # Right panel: ratio potential / backreaction
    for sigma, col in zip(sigma_vals, colors):
        ratio = (1.0/56.0) / (4.0 * sigma**2) * np.ones_like(Phi0_arr)
        ax2.axhline(ratio[0], color=col, ls="-", lw=1.8,
                     label=rf"$\sigma_\Phi = {sigma:.0e}$: ratio = {ratio[0]:.1e}")

    ax2.set_xlabel(r"$\Phi_0$", fontsize=13)
    ax2.set_ylabel(r"$P(1) / \Lambda_{\rm back}$", fontsize=13)
    ax2.set_title("Potential / Backreaction ratio", fontsize=14)
    ax2.set_yscale("log")
    ax2.legend(fontsize=10)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.set_xlim(1, 100)

    fig.suptitle("TGP: Cosmological constant components", fontsize=15, y=1.02)
    fig.tight_layout()
    path = os.path.join(save_dir, "lambda_eff_components.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3: Parameter space contour plot
# ═══════════════════════════════════════════════════════════════════════════
def plot_parameter_space(save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 7))

    # gamma in units of H_0^2/c_0^2
    gamma_unit = H0_si**2 / c0**2  # base unit

    # Scan over gamma/gamma_unit and Phi_0
    g_factors = np.logspace(-1, 2, 300)
    Phi0_arr = np.linspace(1, 100, 300)
    G, P = np.meshgrid(g_factors, Phi0_arr)

    gamma_grid = G * gamma_unit
    # Lambda_eff = gamma / 56 (correct action potential P(1) = gamma/56)
    L_eff = gamma_grid / 56.0
    ratio = np.log10(L_eff / Lambda_obs)

    im = ax.pcolormesh(g_factors, Phi0_arr, ratio,
                       cmap="RdBu_r", shading="auto",
                       vmin=-2, vmax=2)
    cbar = fig.colorbar(im, ax=ax,
                        label=r"$\log_{10}(\Lambda_{\rm eff}/\Lambda_{\rm obs})$")

    # Contour: Lambda_eff = Lambda_obs
    cs = ax.contour(g_factors, Phi0_arr, ratio,
                    levels=[0], colors=["lime"], linewidths=2.5)
    ax.clabel(cs, fmt=r"$\Lambda_{\rm eff} = \Lambda_{\rm obs}$", fontsize=10)

    # Natural region: gamma = Phi_0 * gamma_unit (the tau_0 ~ 1/H_0 locus)
    # At this locus: gamma/gamma_unit = Phi_0
    # So g_factor = Phi_0
    ax.plot(Phi0_arr, Phi0_arr, "k--", lw=2,
            label=r"$\gamma = \Phi_0 H_0^2/c_0^2$ ($\tau_0 \sim 1/H_0$)")

    ax.set_xlabel(r"$\gamma / (H_0^2/c_0^2)$", fontsize=13)
    ax.set_ylabel(r"$\Phi_0$", fontsize=13)
    ax.set_xscale("log")
    ax.set_title(r"TGP: $\Lambda_{\rm eff}$ in parameter space ($\beta = \gamma$)",
                 fontsize=14)
    ax.legend(fontsize=10, loc="upper left")

    fig.tight_layout()
    path = os.path.join(save_dir, "lambda_eff_parameter_space.png")
    fig.savefig(path, dpi=180)
    print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════
def print_summary():
    print("\n" + "=" * 65)
    print("  TGP Effective Cosmological Constant -- Summary")
    print("=" * 65)

    print(f"\n  Physical constants:")
    print(f"    c_0             = {c0:.3e} m/s")
    print(f"    G_0             = {G0:.3e} m^3/(kg s^2)")
    print(f"    H_0             = {H0_si:.3e} s^-1  (~70 km/s/Mpc)")
    print(f"    Lambda_obs      = {Lambda_obs:.4e} m^-2")
    print(f"    R_Hubble        = {hubble_radius():.4e} m")

    print(f"\n  TGP formulae:")
    print(f"    gamma_natural   = Phi_0 * H_0^2 / c_0^2")
    print(f"    Lambda_eff      = gamma / 56     (for beta = gamma, correct action)")
    print(f"    Lambda_eff      = Phi_0 * H_0^2 / (56 c_0^2)")

    print(f"\n  Comparison with observations:")
    for Phi0 in [10, 25, 50, 80, 100, 115, 150]:
        gamma = gamma_natural(Phi0)
        L_eff = gamma / 56.0
        L_back = Lambda_eff_backreaction(Phi0)
        m_sp = m_sp_from_gamma(gamma)
        lam_c = compton_wavelength(gamma)
        R_H = hubble_radius()

        print(f"\n    Phi_0 = {Phi0:6.1f}:")
        print(f"      gamma           = {gamma:.4e} m^-2")
        print(f"      Lambda_eff(pot) = {L_eff:.4e} m^-2")
        print(f"      Lambda_eff(bk)  = {L_back:.4e} m^-2  (backreaction)")
        print(f"      L_eff/L_obs     = {L_eff/Lambda_obs:.4f}")
        print(f"      m_sp            = {m_sp:.4e} m^-1")
        print(f"      lambda_Compton  = {lam_c:.4e} m")
        print(f"      lambda_C / R_H  = {lam_c/R_H:.2f}")

    print(f"\n  Key insight:")
    print(f"    Phi_0 ~ 115 gives Lambda_eff ~ Lambda_obs")
    print(f"    This is O(100), NOT 10^122 -- no fine-tuning!")
    print(f"    Compton wavelength ~ Hubble radius (consistency)")

    print(f"\n  Cosmological constant problem resolution:")
    print(f"    Standard: Lambda requires 10^-122 fine-tuning in Planck units")
    print(f"    TGP:      Lambda ~ gamma ~ H_0^2/c_0^2 naturally")
    print(f"              The mass of the space-generating field is ~H_0/c_0")
    print(f"              Its range is ~Hubble radius (self-consistent)")
    print("=" * 65)


# ═══════════════════════════════════════════════════════════════════════════
# w_DE estimation: equation of state from slow-roll
# ═══════════════════════════════════════════════════════════════════════════
def estimate_w_DE(Phi0=115.0):
    """
    Estimate the dark energy equation of state w_DE.

    With the correct unified action S[g] = int[1/2 g^4 (nabla g)^2 + P(g)],
    the cosmological field equation is:

        psi'' + 3H psi' + 3(psi')^2/psi = c0^2 (gamma*psi - beta)

    At the vacuum psi = 1 (with beta = gamma):
        W(1) = c0^2 (gamma - beta) = 0   (true equilibrium!)
        W'(1) = c0^2 * gamma              (restoring force)

    Since W(1) = 0, the field sits at a stable minimum.
    Hubble friction (3H psi') damps any perturbation exponentially.
    The screening mass m^2 = gamma >> H^2, so psi is frozen at 1.

    Result: w_DE = -1 EXACTLY (pure cosmological constant behavior).
    This is confirmed by numerical ODE integration (see w_de_redshift.py).
    """
    gamma = gamma_natural(Phi0)

    # Correct cosmological potential: W(psi) = c0^2 * (gamma*psi - beta)
    # At vacuum psi = 1 (beta = gamma):
    W1 = 0.0          # gamma - beta = 0: true equilibrium
    Wp1 = gamma        # restoring force coefficient (in c0^2 units)

    # Screening mass vs Hubble rate
    # m^2 = gamma ~ H_0^2 / c_0^2, but m_phys = sqrt(gamma) * c_0 >> H_0
    # The field responds on timescale ~ 1/m_phys << 1/H_0
    # Any deviation from psi=1 is damped within a Hubble time

    # Slow-roll parameter: epsilon = 0 because W(1) = 0
    epsilon = 0.0

    print(f"\n  Dark energy EOS estimate (Phi_0 = {Phi0}):")
    print(f"    Correct action: S[g] = int[1/2 g^4 (nabla g)^2 + P(g)]")
    print(f"    Cosmological potential: W(psi) = c0^2 (gamma*psi - beta)")
    print(f"    W(1) = c0^2 (gamma - beta) = {W1:.4e}  (exact equilibrium)")
    print(f"    W'(1) = c0^2 * gamma = {Wp1:.4e}  (restoring force)")
    print(f"    -> psi = 1 is a STABLE vacuum: field frozen by Hubble friction")
    print(f"    -> w_DE = -1 exactly (cosmological constant)")
    print(f"    -> No slow-roll: field is AT the minimum, not rolling toward it")
    print(f"    -> Confirmed by numerical ODE (w_de_redshift.py)")

    return epsilon


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 65)
    print("TGP effective cosmological constant estimation")
    print("=" * 65)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')

    print("\n  [1/4] Lambda_eff vs Phi_0 ...")
    plot_lambda_vs_Phi0(save_dir=save_dir)

    print("\n  [2/4] Component analysis ...")
    plot_components(save_dir=save_dir)

    print("\n  [3/4] Parameter space ...")
    plot_parameter_space(save_dir=save_dir)

    print("\n  [4/4] Summary and w_DE estimate ...")
    print_summary()
    estimate_w_DE()

    print("\nDone.")


if __name__ == "__main__":
    main()
