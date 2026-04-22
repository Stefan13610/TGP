"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
effective_potential.py  --  Theory of Generated Space (TGP)
==========================================================
Compute the effective interaction energy V_eff(d) between two equal
point sources separated by distance d.

Physics
-------
The TGP energy functional (eq:energy-corrected, f(Phi) ≈ 1 for weak fields):

    E[Phi] = integral { (1/2) f(Phi) (grad Phi)^2
                       + (beta/3 Phi0) Phi^3
                       - (gamma/4 Phi0^2) Phi^4 } d^3x

Interaction energy V_eff(d) = E_total(d) - 2 * E_single, computed
analytically from overlap of single-source Yukawa profiles.

THREE REGIMES AT beta = gamma (vacuum condition)
-------------------------------------------------
The three contributions to V_eff(d) are:

    E_lin   = -4*pi*C^2 / d         (gradient, attractive)
    E_beta  = +8*pi*beta*C^2 / d^2  (quadratic, repulsive)
    E_gamma = -24*pi*gamma*C^3 / d^3 (cubic, confining)

For beta = gamma: the repulsive and confining terms share the SAME
coupling constant.  Three regimes arise NOT from separate couplings,
but from different powers of C:
  - E_beta ~ C^2 / d^2  (dominates at intermediate d when C < 1)
  - E_gamma ~ C^3 / d^3  (dominates at short d, suppressed by extra C)

Three regimes exist when beta > 9*C/2 (Proposition prop:trzy-rezimy-beta-gamma,
sek03).  For elementary particles C << 1, so this is trivially satisfied
for any beta > 0.  For macroscopic objects C >> 1, only regime I (gravity)
survives — consistent with observation.

Source sign: TGP convention (mass generates excess space, -qΦ₀ρ on RHS).

Outputs  (saved to scripts/plots/):
    effective_potential.png    -- V_eff(d) at beta=gamma (vacuum), three regimes
    regime_map.png             -- parameter space scan (beta=gamma, C)
"""

import os
import numpy as np
from scipy.integrate import simpson
from scipy.special import erfc
import matplotlib.pyplot as plt

# ── TGP parameters ─────────────────────────────────────────────────────────
# Gradient coupling constant.  At leading order in the weak-field expansion,
# alpha does not appear in V_eff(d) because f(Phi) = 1 + 2*alpha*ln(Phi/Phi0)
# reduces to f ≈ 1 for the linearized (Yukawa) profiles used here.  The
# constant is retained for use in higher-order or numerical computations.
ALPHA = 2.0

# Separation distances (dimensionless units)
D_MIN, D_MAX, N_D = 0.5, 30.0, 120


# ── Linearized Yukawa profile ──────────────────────────────────────────────
def yukawa_profile(r, C, m_eff):
    """
    Single-source linearized profile: delta(r) = C * exp(-m_eff * r) / r
    where delta = f - 1 = (Phi - Phi0) / Phi0.
    """
    r_safe = np.maximum(r, 1e-10)
    return C * np.exp(-m_eff * r_safe) / r_safe


def yukawa_gradient(r, C, m_eff):
    """d(delta)/dr for Yukawa profile."""
    r_safe = np.maximum(r, 1e-10)
    return -C * np.exp(-m_eff * r_safe) * (m_eff + 1.0 / r_safe) / r_safe


# ── Analytical V_eff for two Yukawa sources ───────────────────────────────
def veff_yukawa_analytic(d, C, m_eff, bg):
    """
    Analytical interaction energy for two identical Yukawa sources
    at distance d.  In the linearized regime:

    V_lin(d) = -C^2 * (2*pi / d) * exp(-m_eff * d) * (1 + 1/(m_eff*d))

    For the TGP sign convention (mass → excess space), the gradient
    interaction is attractive (negative), while the cubic self-interaction
    produces repulsion (positive) at intermediate scales.

    This function computes the numerical 2D integral for better accuracy.
    """
    # Cylindrical grid (exploit axial symmetry around z-axis)
    N_rho, N_z = 150, 300
    rho_max = max(30.0, 3 * d)
    z_max = max(40.0, 2 * d)
    rho = np.linspace(0.01, rho_max, N_rho)
    z = np.linspace(-z_max, z_max, N_z)
    RHO, Z = np.meshgrid(rho, z, indexing="ij")

    # Distances from each source (at z = ±d/2)
    r1 = np.sqrt(RHO**2 + (Z - d / 2)**2)
    r2 = np.sqrt(RHO**2 + (Z + d / 2)**2)

    # Profiles
    d1 = yukawa_profile(r1, C, m_eff)
    d2 = yukawa_profile(r2, C, m_eff)

    # Gradients (cylindrical components)
    # For source 1: grad_rho = (rho/r1) * d(delta)/dr1, grad_z = ((z-d/2)/r1) * d(delta)/dr1
    ddr1 = yukawa_gradient(r1, C, m_eff)
    ddr2 = yukawa_gradient(r2, C, m_eff)

    r1_safe = np.maximum(r1, 1e-10)
    r2_safe = np.maximum(r2, 1e-10)

    grad1_rho = RHO / r1_safe * ddr1
    grad1_z = (Z - d / 2) / r1_safe * ddr1
    grad2_rho = RHO / r2_safe * ddr2
    grad2_z = (Z + d / 2) / r2_safe * ddr2

    # E_lin: gradient cross-term  = integral { grad(d1) . grad(d2) } 2*pi*rho
    integrand_lin = (grad1_rho * grad2_rho + grad1_z * grad2_z) * 2 * np.pi * RHO

    # E_cubic: cross-terms from beta * (f^2 - f^3) expanded around f = 1 + d1 + d2
    # f = 1 + d1 + d2: f^2 - f^3 = -(d1+d2) - 3(d1+d2)^2/2 - ...
    # Cross-term: -3 * d1 * d2 at leading order, times bg
    integrand_cubic = -3.0 * bg * d1 * d2 * 2 * np.pi * RHO

    def integrate_2d(integrand):
        inner = simpson(integrand, x=z, axis=1)
        return simpson(inner, x=rho)

    E_lin = integrate_2d(integrand_lin)
    E_cubic = integrate_2d(integrand_cubic)
    V_total = E_lin + E_cubic

    return V_total, E_lin, E_cubic


# ── Simple analytical formula (valid for d >> 1/m_eff) ─────────────────────
def veff_three_regime(d, C, beta=None, gamma=None, vacuum=True):
    """
    Leading-order analytic TGP interaction potential from the energy functional
    decomposition (Theorem thm:three-regimes, eq. eq:Eint-decomp, sek08).

    The three contributions arise from overlap integrals of single-source
    profiles delta_Phi(r) = C/r (Proposition prop:trzy-rezimy-beta-gamma,
    sek08):

        E_lin(d)   = -4*pi*C^2 / d          (gradient, attractive, longest range)
        E_beta(d)  = +8*pi*beta*C^2 / d^2   (quadratic, repulsive, medium range)
        E_gamma(d) = -24*pi*gamma*C^3 / d^3 (cubic, confining, shortest range)

    Parameters
    ----------
    d : float
        Separation between the two sources.
    C : float
        Source strength (dimensionless charge).
    beta, gamma : float or None
        Self-interaction couplings.  When *vacuum=True* (default) the vacuum
        condition beta = gamma is enforced: pass a single value via *beta*
        (or *gamma*) and the other is set equal to it.  When *vacuum=False*
        both must be supplied explicitly (non-vacuum / cosmological density
        offset regime).
    vacuum : bool
        If True (default), enforce beta = gamma (vacuum condition).

    Notes
    -----
    The effective screening mass in the vacuum case is
        m_sp = sqrt(3*gamma - 2*beta) = sqrt(gamma)   [for beta = gamma].
    """
    # ── resolve vacuum vs. non-vacuum parameters ──
    if vacuum:
        # Vacuum condition: beta = gamma
        val = beta if beta is not None else (gamma if gamma is not None else 1.0)
        beta = val
        gamma = val
    else:
        if beta is None or gamma is None:
            raise ValueError(
                "Both beta and gamma must be supplied in non-vacuum mode "
                "(vacuum=False)."
            )
    d_safe = max(d, 0.1)
    # Gradient overlap (attractive)
    V_grad = -4 * np.pi * C**2 / d_safe
    # Beta term: overlap of delta1 * delta2 (repulsive, ∝ 1/d^2)
    V_beta = 8 * np.pi * beta * C**2 / d_safe**2
    # Gamma term: overlap of delta1 * delta2 * (delta1 + delta2) (confining, ∝ 1/d^3)
    V_gamma = -24 * np.pi * gamma * C**3 / d_safe**3

    return V_grad + V_beta + V_gamma, V_grad, V_beta, V_gamma


# ── Plotting ────────────────────────────────────────────────────────────────
def plot_effective_potential(save_dir=None):
    """
    Main figure: V_eff(d) demonstrating three regimes.

    KEY PHYSICS (Proposition prop:trzy-rezimy-beta-gamma, sek03):
    Three regimes exist ALREADY at beta = gamma (vacuum condition),
    provided  beta > 9*C/2.  For elementary particles (C << 1) this
    condition is trivially satisfied for any beta > 0.

    The three contributions are:
      E_lin   = -4*pi*C^2 / d           (gradient, attractive, 1/d)
      E_beta  = +8*pi*beta*C^2 / d^2    (repulsive, 1/d^2)
      E_gamma = -24*pi*gamma*C^3 / d^3  (confining, 1/d^3)

    For beta = gamma, E_beta and E_gamma are BOTH controlled by the
    SAME coupling constant.  This is NOT a problem: the 1/d^2 term
    has a C^2 prefactor while the 1/d^3 term has C^3, so for C << 1
    the beta-term (repulsion) dominates at intermediate d, and the
    gamma-term (confinement) dominates at short d.  The mechanism is:
    different powers of C, NOT different coupling constants.
    """
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    # ── Figure 1: Vacuum case beta = gamma (THE physical scenario) ──
    # Three regimes exist when beta > 9C/2  (trivially true for C << 1).
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: V_eff for several C values at fixed beta = gamma
    beta_vac = 5.0  # beta = gamma (vacuum)
    C_values = [
        (0.3, r"$C = 0.3,\;\beta/C = 16.7$  (3 regimes)"),
        (0.5, r"$C = 0.5,\;\beta/C = 10$  (3 regimes)"),
        (1.0, r"$C = 1.0,\;\beta/C = 5$  (3 regimes)"),
        (1.2, r"$C = 1.2,\;\beta/C = 4.2$  (1 regime, subcritical)"),
    ]
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(C_values)))

    d_arr = np.linspace(0.3, 30.0, 200)

    for (C, label), col in zip(C_values, colors):
        ratio = beta_vac / C
        print(f"  Vacuum case: beta=gamma={beta_vac}, C={C}, beta/C={ratio:.1f} ...")

        V_arr = np.array([veff_three_regime(d, C, beta=beta_vac)[0] for d in d_arr])

        # Normalize for visibility
        V_scale = np.max(np.abs(V_arr[np.isfinite(V_arr)])) or 1.0
        ax1.plot(d_arr, V_arr / V_scale, color=col, lw=2, label=label)

        # Mark zero crossings
        for i in range(len(V_arr) - 1):
            if V_arr[i] * V_arr[i + 1] < 0:
                d_cross = d_arr[i] - V_arr[i] * (d_arr[i+1] - d_arr[i]) / (V_arr[i+1] - V_arr[i])
                ax1.axvline(d_cross, color=col, ls=":", lw=0.8, alpha=0.5)

    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_xlabel(r"Separation $d / r_0$", fontsize=13)
    ax1.set_ylabel(r"$V_{\rm eff}(d) / |V|_{\max}$", fontsize=13)
    ax1.set_title(
        rf"Three regimes at $\beta = \gamma = {beta_vac}$ (vacuum)"
        "\n" + r"Critical: $\beta > 9C/2$ (Prop. 3, sek03)",
        fontsize=12)
    ax1.legend(fontsize=8)
    ax1.grid(True, ls=":", alpha=0.4)
    ax1.set_ylim(-1.5, 1.5)

    # Right panel: Energy decomposition for vacuum case
    C0 = 0.5
    V_arr = np.zeros(len(d_arr))
    Vg_arr = np.zeros(len(d_arr))
    Vb_arr = np.zeros(len(d_arr))
    Vc_arr = np.zeros(len(d_arr))
    for i, d in enumerate(d_arr):
        V, Vg, Vb, Vc = veff_three_regime(d, C0, beta=beta_vac)  # vacuum=True default
        V_arr[i] = V
        Vg_arr[i] = Vg
        Vb_arr[i] = Vb
        Vc_arr[i] = Vc

    V_scale = np.max(np.abs(V_arr[np.isfinite(V_arr)])) or 1.0
    ax2.plot(d_arr, V_arr / V_scale, "k-", lw=2.2,
             label=r"$V_{\rm eff}$ (total)")
    ax2.plot(d_arr, Vg_arr / V_scale, "b--", lw=1.5,
             label=r"$-4\pi C^2/d$ (gradient, attractive)")
    ax2.plot(d_arr, Vb_arr / V_scale, "r-.", lw=1.5,
             label=r"$+8\pi\beta C^2/d^2$ (repulsive)")
    ax2.plot(d_arr, Vc_arr / V_scale, "m:", lw=1.5,
             label=r"$-24\pi\beta C^3/d^3$ (confining)")
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set_xlabel(r"Separation $d / r_0$", fontsize=13)
    ax2.set_ylabel(r"Energy / $|V|_{\max}$", fontsize=13)
    ax2.set_title(
        rf"Decomposition: $\beta = \gamma = {beta_vac}$, $C = {C0}$"
        "\n" + r"Same coupling $\beta$; separation by powers of $C$",
        fontsize=12)
    ax2.legend(fontsize=8)
    ax2.grid(True, ls=":", alpha=0.4)
    ax2.set_ylim(-1.5, 1.5)

    # Annotate regimes
    ax2.annotate("I: gravity", xy=(18, -0.15), fontsize=11, color="blue",
                 ha="center", style="italic")
    ax2.annotate("II: repulsion", xy=(5, 0.3), fontsize=11, color="red",
                 ha="center", style="italic")
    ax2.annotate("III: confinement", xy=(1.5, -0.8), fontsize=11, color="purple",
                 ha="center", style="italic")

    fig.suptitle(
        r"TGP three-regime potential: $\beta = \gamma$ (vacuum condition)"
        "\n" + r"Three regimes from ONE coupling constant "
        r"($\sim C^2/d$ vs $\sim C^2/d^2$ vs $\sim C^3/d^3$)",
        fontsize=13, y=1.06)
    fig.tight_layout()
    path = os.path.join(save_dir, "effective_potential.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)


def plot_regime_map(save_dir=None):
    """
    Scan beta-gamma parameter space and show where V_eff changes sign
    (boundary between attractive and repulsive regimes).
    """
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)
    d_arr = np.linspace(D_MIN, D_MAX, 80)

    betas = np.linspace(0.005, 0.2, 25)
    C_vals = np.linspace(0.5, 5.0, 25)
    d_cross_map = np.full((len(betas), len(C_vals)), np.nan)
    V_min_map = np.full((len(betas), len(C_vals)), np.nan)

    for i, bg in enumerate(betas):
        # Screening mass in the vacuum case (beta = gamma):
        #   m_sp = sqrt(3*gamma - 2*beta) = sqrt(gamma)  when beta = gamma.
        m_sp = np.sqrt(bg)
        for j, C in enumerate(C_vals):
            # Vacuum condition: beta = gamma = bg
            V_arr = np.array([veff_three_regime(d, C, beta=bg)[0] for d in d_arr])
            V_min_map[i, j] = np.min(V_arr)

            # Zero crossing
            for k in range(len(V_arr) - 1):
                if V_arr[k] * V_arr[k + 1] < 0:
                    d_cross = d_arr[k] - V_arr[k] * (d_arr[k+1] - d_arr[k]) / (V_arr[k+1] - V_arr[k])
                    d_cross_map[i, j] = d_cross
                    break

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    extent = [C_vals[0], C_vals[-1], betas[0], betas[-1]]

    im1 = ax1.imshow(d_cross_map, origin="lower", aspect="auto", extent=extent,
                     cmap="YlOrRd")
    ax1.set_xlabel(r"Source strength $C$", fontsize=13)
    ax1.set_ylabel(r"$\hat\beta = \hat\gamma$", fontsize=13)
    ax1.set_title(r"Zero-crossing distance $d_{\rm cross}/r_0$", fontsize=13)
    fig.colorbar(im1, ax=ax1)

    im2 = ax2.imshow(V_min_map, origin="lower", aspect="auto", extent=extent,
                     cmap="RdBu_r", vmin=-0.5, vmax=0.5)
    ax2.set_xlabel(r"Source strength $C$", fontsize=13)
    ax2.set_ylabel(r"$\hat\beta = \hat\gamma$", fontsize=13)
    ax2.set_title(r"Minimum $V_{\rm eff}$ (well depth)", fontsize=13)
    fig.colorbar(im2, ax=ax2)

    # Add critical line beta = 9C/2 to both panels
    C_crit_line = np.linspace(C_vals[0], C_vals[-1], 100)
    beta_crit_line = 4.5 * C_crit_line
    for ax in (ax1, ax2):
        mask = beta_crit_line <= betas[-1]
        ax.plot(C_crit_line[mask], beta_crit_line[mask], "w--", lw=2.0,
                label=r"$\beta = 9C/2$ (critical)")
        ax.legend(fontsize=8, loc="upper left")

    fig.suptitle(r"TGP potential landscape ($\beta = \gamma$, vacuum)"
                 "\n" + r"Three regimes exist above the critical line $\beta > 9C/2$",
                 fontsize=13, y=1.05)
    fig.tight_layout()
    path = os.path.join(save_dir, "regime_map.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)


# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("TGP effective potential computation")
    print("=" * 60)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')

    plot_effective_potential(save_dir=save_dir)
    print()
    print("  Computing parameter space map ...")
    plot_regime_map(save_dir=save_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
