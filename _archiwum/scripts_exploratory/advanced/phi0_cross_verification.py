"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
phi0_cross_verification.py  --  Theory of Generated Space (TGP)
================================================================
Cross-verify the TGP parameter Phi_0 from multiple independent
observational constraints.

Phi_0 is the cosmological background value of the scalar field Phi
(space-density field).  It enters:
  - Dynamic constants: c(Phi)=c0*(Phi/Phi0)^{-1/2}, G(Phi)=G0*(Phi/Phi0)^{-1}
  - Dark energy:  Lambda_eff = gamma/12,  gamma ~ Phi0*H0^2/c0^2
  - Relaxation:   tau0 = sqrt(Phi0)/(c0*sqrt(gamma))
  - Scale-dep G:  G_eff(k,a) = G0*(1 + 2*alpha_eff^2/(1+(a*m_eff/k)^2))

Six constraints are evaluated:
  1. Dark energy  (Lambda_obs -> gamma -> Phi0  via  tau0 ~ 1/H0)
  2. BBN          (tau0*H_BBN >> 1  => lower bound on Phi0)
  3. CMB          (field frozen at recombination)
  4. GW speed     (|c_GW - c0|/c0 < 1e-15)
  5. PPN / Cassini (gamma_PPN - 1 < 2.3e-5)
  6. Galaxy clusters (effective G modification from X-ray / lensing)

Output:
  TGP/TGP_v1/plots/phi0_cross_verification.png
  Console summary of the consistent Phi_0 range.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import scipy.constants as const

# =====================================================================
# Physical constants from scipy.constants (SI)
# =====================================================================
c0    = const.c                    # 2.998e8  m/s
G0    = const.G                    # 6.674e-11  m^3 kg^-1 s^-2
hbar0 = const.hbar                 # 1.0546e-34 J s
k_B   = const.k                    # Boltzmann constant

# Cosmological parameters
H0_km  = 70.0                      # km/s/Mpc
Mpc_m  = 3.0857e22                 # metres per Mpc
H0_SI  = H0_km * 1e3 / Mpc_m      # s^-1  (~2.27e-18)

Omega_Lambda = 0.685               # Planck 2018
Omega_m      = 0.315
Omega_r      = 9.1e-5              # radiation today

rho_crit = 3.0 * H0_SI**2 / (8.0 * np.pi * G0)       # kg/m^3
rho_Lambda_obs = Omega_Lambda * rho_crit               # kg/m^3

# Observed Lambda in geometric units  [m^-2]
Lambda_obs = 3.0 * Omega_Lambda * H0_SI**2 / c0**2     # ~1.1e-52

# Planck length (invariant in TGP)
l_P = np.sqrt(hbar0 * G0 / c0**3)

# BBN / CMB reference epochs
T_BBN   = 1.0e6 * const.eV / k_B          # ~1 MeV  in K
z_BBN   = 3.0e9
H_BBN   = H0_SI * np.sqrt(Omega_r) * (1 + z_BBN)**2   # rad-dominated H

z_rec   = 1089.0                                        # recombination
H_rec   = H0_SI * np.sqrt(Omega_m * (1 + z_rec)**3
                           + Omega_r * (1 + z_rec)**4)

# Solar-system numbers
r_S_sun = 2.0 * G0 * 1.989e30 / c0**2    # Schwarzschild radius of Sun
r_AU    = const.au                          # 1 AU in metres


# =====================================================================
# TGP helper: gamma and tau0
# =====================================================================
def gamma_of(Phi0):
    """gamma = Phi0 * H0^2 / c0^2   [m^-2]."""
    return Phi0 * H0_SI**2 / c0**2


def tau0_of(Phi0):
    """tau0 = sqrt(Phi0) / (c0 * sqrt(gamma))   [s]."""
    g = gamma_of(Phi0)
    return np.sqrt(Phi0) / (c0 * np.sqrt(g))


# =====================================================================
# Constraint 1 : Dark energy   Lambda_eff = gamma/12 = Lambda_obs
# =====================================================================
def constraint_dark_energy(Phi0_arr):
    """
    Lambda_eff = gamma/12 = Phi0*H0^2 / (12*c0^2).
    Setting Lambda_eff = Lambda_obs  =>  Phi0_DE.

    Also expressible through energy density:
        rho_DE = Lambda_eff * c0^2 / (8*pi*G0)
               = Phi0 * H0^2 / (96*pi*G0)
    """
    Lambda_eff = gamma_of(Phi0_arr) / 12.0
    Phi0_central = Lambda_obs * 12.0 * c0**2 / H0_SI**2

    # 1-sigma from Planck  Omega_Lambda = 0.685 +/- 0.007
    Phi0_lo = (0.685 - 0.007) / 0.685 * Phi0_central
    Phi0_hi = (0.685 + 0.007) / 0.685 * Phi0_central

    return dict(
        Lambda_eff=Lambda_eff,
        Phi0_central=Phi0_central,
        Phi0_lo=Phi0_lo,
        Phi0_hi=Phi0_hi,
        label="Dark energy",
    )


# =====================================================================
# Constraint 2 : BBN   tau0 * H_BBN >> 1  (field frozen)
# =====================================================================
def constraint_bbn(Phi0_arr):
    """
    Require tau0 * H_BBN > N_freeze  (we take N_freeze = 10).
    tau0 = sqrt(Phi0)/(c0*sqrt(gamma)) = c0 / (H0*sqrt(Phi0))   [using gamma def]

    Actually:  tau0 = sqrt(Phi0)/(c0*sqrt(gamma))
      gamma = Phi0*H0^2/c0^2
      sqrt(gamma) = sqrt(Phi0)*H0/c0
      tau0 = sqrt(Phi0)/(c0 * sqrt(Phi0)*H0/c0) = 1/H0

    So tau0 = 1/H0 identically for the natural gamma!
    Then tau0*H_BBN = H_BBN/H0 >> 1, which is automatically satisfied
    since H_BBN ~ 10^{18} * H0.

    However, for gamma != Phi0*H0^2/c0^2 (i.e. gamma as free param),
    the condition tau0*H_BBN >> 1 gives a LOWER bound on Phi0.

    More physically: the BBN constraint |Delta G/G| < 0.13 combined
    with the field equation gives phi(z_BBN) ~ phi_0/(tau0*H_BBN)^2.
    Requiring this < 0.13:
        Phi0 > (0.13)^{-1} * (H0/H_BBN)^2  ...  extremely weak.

    For an interesting lower bound, consider the case where gamma
    is fixed by Lambda_obs (gamma = 12*Lambda_obs) independently:
        tau0 = sqrt(Phi0) / (c0 * sqrt(12*Lambda_obs))
        tau0 * H_BBN > 10
        Phi0 > (10 * c0 * sqrt(12*Lambda_obs) / H_BBN)^2
    """
    gamma_fixed = 12.0 * Lambda_obs        # from DE constraint
    # tau0(Phi0) = sqrt(Phi0) / (c0 * sqrt(gamma_fixed))
    tau0_arr = np.sqrt(Phi0_arr) / (c0 * np.sqrt(gamma_fixed))
    ratio = tau0_arr * H_BBN

    N_freeze = 10.0
    # Lower bound:  sqrt(Phi0_min) = N_freeze * c0 * sqrt(gamma_fixed) / H_BBN
    Phi0_min = (N_freeze * c0 * np.sqrt(gamma_fixed) / H_BBN)**2

    return dict(
        tau0_H_BBN=ratio,
        Phi0_min=Phi0_min,
        N_freeze=N_freeze,
        label="BBN (frozen field)",
    )


# =====================================================================
# Constraint 3 : CMB   field frozen at recombination
# =====================================================================
def constraint_cmb(Phi0_arr):
    """
    Same logic as BBN but at z_rec = 1089.
    Require tau0 * H_rec > N_freeze.
    With gamma from DE:
        Phi0 > (N_freeze * c0 * sqrt(12*Lambda_obs) / H_rec)^2

    Also: the ISW effect constrains time-variation of the potential
    Phi_dot/Phi < epsilon * H0  =>  additional bound.
    CMB anisotropy constrains  delta(G)/G < ~0.05  at recombination.
    """
    gamma_fixed = 12.0 * Lambda_obs
    tau0_arr = np.sqrt(Phi0_arr) / (c0 * np.sqrt(gamma_fixed))
    ratio = tau0_arr * H_rec

    N_freeze = 10.0
    Phi0_min_freeze = (N_freeze * c0 * np.sqrt(gamma_fixed) / H_rec)**2

    # ISW / delta-G constraint:  |1 - psi(z_rec)| < 0.05
    # psi deviation ~ 1/(tau0*H_rec)^2  for frozen field
    # |delta_psi| = 1/(tau0*H_rec)^2 < 0.05
    # tau0*H_rec > 1/sqrt(0.05) ~ 4.47
    # already weaker than N_freeze = 10 above

    # Combined: take the stronger one (N_freeze = 10)
    Phi0_min = Phi0_min_freeze

    return dict(
        tau0_H_rec=ratio,
        Phi0_min=Phi0_min,
        label="CMB (frozen at rec.)",
    )


# =====================================================================
# Constraint 4 : GW speed   |c_GW/c - 1| < 1e-15
# =====================================================================
def constraint_gw_speed(Phi0_arr):
    """
    GW170817 + GRB 170817A:  |c_GW/c - 1| < ~5e-16.

    In TGP (disformal metric):
        c_GW^2/c^2 = A(Phi) / (A + B * Phi_dot^2)

    For the natural (frozen-field) regime  Phi_dot ~ 0:
        c_GW = c  exactly.

    For a slowly evolving field:
        |c_GW/c - 1| ~ (1/2) * B/A * (Phi_dot/c)^2

    With  Phi_dot ~ delta_psi * Phi0 * H0  and  B/A ~ 1/Phi0^2:
        |c_GW/c - 1| ~ (delta_psi * H0)^2 / 2

    The field deviation delta_psi at late times:
        delta_psi ~ 1 / (tau0 * H0)^2

    For natural gamma:  tau0 = 1/H0  =>  delta_psi ~ 1  (too large?).
    But actually the dimensionless field equation gives
        delta_psi ~ gamma / (12 * H0^2/c0^2) = Phi0/12  ... need care.

    More careful:  B(Phi) in TGP has dimension [length^4],
        B ~ l_P^2 / Phi0^2  (from the action),  A ~ Phi/Phi0 ~ 1.
    Then:
        |delta c/c| = (1/2)*(l_P^2/Phi0^2)*(Phi_dot/c)^2
                     = (1/2)*(l_P * H0 * delta_psi)^2
    With l_P*H0 ~ 10^{-61},  this is  ~ 10^{-122} * delta_psi^2.

    So the GW bound gives:
        delta_psi^2 < 2 * 5e-16 / (l_P * H0)^2
    => Essentially NO constraint on Phi0  (trivially satisfied).

    For visualization, we parametrize as:
        |delta c/c|(Phi0) ~ (l_P * H0)^2 / (2 * Phi0)
    giving a formal (very weak) lower bound.
    """
    gw_bound = 5e-16

    # |delta c/c| ~ (l_P * H0)^2 / (2 * Phi0)   [crude parametrization]
    delta_c = (l_P * H0_SI)**2 / (2.0 * Phi0_arr)

    # formal lower bound (extremely small)
    Phi0_min = (l_P * H0_SI)**2 / (2.0 * gw_bound)

    return dict(
        delta_c_over_c=delta_c,
        gw_bound=gw_bound,
        Phi0_min=Phi0_min,
        label="GW speed (GW170817)",
    )


# =====================================================================
# Constraint 5 : PPN / Cassini   |gamma_PPN - 1| < 2.3e-5
# =====================================================================
def constraint_ppn(Phi0_arr):
    """
    PPN deviation in TGP from the gradient of Phi around the Sun:
        delta_Phi/Phi0 ~ G*M_sun / (c0^2 * r * Phi0)  = r_S / (2*r*Phi0)

    The PPN parameter deviation:
        |gamma_PPN - 1| ~ (r_S / (2*r_AU*Phi0))^2

    Cassini bound:  2.3e-5.

        Phi0 > r_S / (2*r_AU*sqrt(2.3e-5))
    """
    cassini = 2.3e-5
    dev = (r_S_sun / (2.0 * r_AU * Phi0_arr))**2

    Phi0_min = r_S_sun / (2.0 * r_AU * np.sqrt(cassini))

    return dict(
        gamma_ppn_dev=dev,
        cassini_bound=cassini,
        Phi0_min=Phi0_min,
        label="PPN / Cassini",
    )


# =====================================================================
# Constraint 6 : Galaxy clusters  (X-ray / lensing G modification)
# =====================================================================
def constraint_clusters(Phi0_arr):
    """
    Galaxy-cluster masses from X-ray gas profiles vs gravitational
    lensing constrain  |G_eff/G0 - 1| < ~0.1  on scales k ~ 1/Mpc.

    In TGP:
        G_eff(k) = G0 * (1 + 2*alpha_eff^2 / (1 + (m_sp/k)^2))

    where  m_sp = sqrt(gamma)  [spatial mass, m^{-1}]
    and    alpha_eff ~ 1/(2*Phi0)  [coupling].

    At cluster scales  k ~ 2*pi / (1 Mpc)  ~ 2e-22 m^{-1}:
        m_sp = sqrt(gamma) ~ sqrt(Phi0) * H0/c0  ~ sqrt(Phi0) * 7.6e-27 m^{-1}

    For Phi0 ~ 10-50:  m_sp ~ 2-5e-26 m^{-1}  <<  k_cluster.
    So  (m_sp/k)^2 ~ (3e-26/2e-22)^2 ~ 2e-8  <<  1.

    Then:
        G_eff/G0 - 1 ~ 2*alpha_eff^2 = 2/(4*Phi0^2) = 1/(2*Phi0^2)

    Constraint:  1/(2*Phi0^2) < 0.1  =>  Phi0 > sqrt(5) ~ 2.24.

    More conservatively with 5% precision:
        Phi0 > sqrt(10) ~ 3.16.
    """
    cluster_bound = 0.10          # |G_eff/G - 1| < 10%
    k_cluster = 2.0 * np.pi / (Mpc_m)   # ~2e-22 m^-1

    gamma_arr = gamma_of(Phi0_arr)
    m_sp = np.sqrt(gamma_arr)
    alpha_eff = 1.0 / (2.0 * Phi0_arr)

    # Full expression
    delta_G = 2.0 * alpha_eff**2 / (1.0 + (m_sp / k_cluster)**2)

    Phi0_min = 1.0 / np.sqrt(2.0 * cluster_bound)   # ~ 2.24

    return dict(
        delta_G=delta_G,
        cluster_bound=cluster_bound,
        Phi0_min=Phi0_min,
        label="Galaxy clusters",
    )


# =====================================================================
# Collect all constraints and find consistent range
# =====================================================================
def compute_all(Phi0_arr):
    c1 = constraint_dark_energy(Phi0_arr)
    c2 = constraint_bbn(Phi0_arr)
    c3 = constraint_cmb(Phi0_arr)
    c4 = constraint_gw_speed(Phi0_arr)
    c5 = constraint_ppn(Phi0_arr)
    c6 = constraint_clusters(Phi0_arr)
    return c1, c2, c3, c4, c5, c6


def find_consistent_range(c1, c2, c3, c4, c5, c6):
    """Return (Phi0_lo, Phi0_hi) where all constraints overlap."""
    # Lower bounds
    lower_bounds = {
        "BBN":      c2["Phi0_min"],
        "CMB":      c3["Phi0_min"],
        "GW speed": c4["Phi0_min"],
        "PPN":      c5["Phi0_min"],
        "Clusters": c6["Phi0_min"],
    }
    # Dark-energy band
    de_lo = c1["Phi0_lo"]
    de_hi = c1["Phi0_hi"]

    overall_lo = max(max(lower_bounds.values()), de_lo)
    overall_hi = de_hi  # DE is the tightest upper

    return overall_lo, overall_hi, lower_bounds


# =====================================================================
# Plotting
# =====================================================================
def make_plot(Phi0_arr, save_path):
    c1, c2, c3, c4, c5, c6 = compute_all(Phi0_arr)

    overall_lo, overall_hi, lower_bounds = find_consistent_range(
        c1, c2, c3, c4, c5, c6
    )

    fig, axes = plt.subplots(3, 2, figsize=(16, 17))

    # ----- Panel (0,0): Dark energy -----
    ax = axes[0, 0]
    ratio = c1["Lambda_eff"] / Lambda_obs
    ax.plot(Phi0_arr, ratio, "b-", lw=2.2,
            label=r"$\Lambda_{\rm eff}/\Lambda_{\rm obs}$")
    ax.axhline(1.0, color="red", ls="--", lw=1.5,
               label=r"$\Lambda_{\rm eff} = \Lambda_{\rm obs}$")
    ax.axvline(c1["Phi0_central"], color="green", ls=":", lw=2,
               label=rf"$\Phi_0 = {c1['Phi0_central']:.1f}$")
    ax.axvspan(c1["Phi0_lo"], c1["Phi0_hi"],
               color="green", alpha=0.18,
               label=rf"Planck $1\sigma$: [{c1['Phi0_lo']:.1f}, {c1['Phi0_hi']:.1f}]")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$\Lambda_{\rm eff}/\Lambda_{\rm obs}$", fontsize=12)
    ax.set_title("1. Dark energy", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_xlim(1, 80)
    ax.set_ylim(0, 5)
    ax.grid(True, ls=":", alpha=0.4)

    # ----- Panel (0,1): BBN -----
    ax = axes[0, 1]
    ax.semilogy(Phi0_arr, c2["tau0_H_BBN"], "b-", lw=2.2,
                label=r"$\tau_0 \cdot H_{\rm BBN}$")
    ax.axhline(c2["N_freeze"], color="red", ls="--", lw=1.5,
               label=rf"Freeze threshold $= {c2['N_freeze']:.0f}$")
    ax.axvline(c2["Phi0_min"], color="orange", ls=":", lw=2,
               label=rf"$\Phi_0^{{\min}} = {c2['Phi0_min']:.2e}$")
    ax.axvspan(c2["Phi0_min"], Phi0_arr[-1],
               color="green", alpha=0.08, label="Allowed")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$\tau_0 \cdot H_{\rm BBN}$", fontsize=12)
    ax.set_title("2. BBN (frozen field)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(1, 80)
    ax.grid(True, ls=":", alpha=0.4)

    # ----- Panel (1,0): CMB -----
    ax = axes[1, 0]
    ax.semilogy(Phi0_arr, c3["tau0_H_rec"], "b-", lw=2.2,
                label=r"$\tau_0 \cdot H_{\rm rec}$")
    ax.axhline(10.0, color="red", ls="--", lw=1.5,
               label="Freeze threshold = 10")
    ax.axvline(c3["Phi0_min"], color="orange", ls=":", lw=2,
               label=rf"$\Phi_0^{{\min}} = {c3['Phi0_min']:.2e}$")
    ax.axvspan(c3["Phi0_min"], Phi0_arr[-1],
               color="green", alpha=0.08, label="Allowed")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$\tau_0 \cdot H_{\rm rec}$", fontsize=12)
    ax.set_title("3. CMB (frozen at recombination)", fontsize=13,
                 fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(1, 80)
    ax.grid(True, ls=":", alpha=0.4)

    # ----- Panel (1,1): GW speed -----
    ax = axes[1, 1]
    ax.semilogy(Phi0_arr, c4["delta_c_over_c"], "b-", lw=2.2,
                label=r"$|\delta c_{\rm GW}/c|$ (TGP)")
    ax.axhline(c4["gw_bound"], color="red", ls="--", lw=1.5,
               label=rf"GW170817 bound: {c4['gw_bound']:.0e}")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$|c_{\rm GW}/c - 1|$", fontsize=12)
    ax.set_title("4. GW speed (GW170817)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(1, 80)
    ax.grid(True, ls=":", alpha=0.4)
    ax.text(0.55, 0.25,
            "Trivially satisfied:\n"
            r"$|\delta c/c| \sim (l_P H_0)^2/\Phi_0$"
            "\n"
            r"$\sim 10^{-122}/\Phi_0$",
            transform=ax.transAxes, fontsize=9, ha="center",
            bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

    # ----- Panel (2,0): PPN / Cassini -----
    ax = axes[2, 0]
    ax.loglog(Phi0_arr, c5["gamma_ppn_dev"], "b-", lw=2.2,
              label=r"$|\gamma_{\rm PPN} - 1|$ (TGP)")
    ax.axhline(c5["cassini_bound"], color="red", ls="--", lw=1.5,
               label=f"Cassini bound: {c5['cassini_bound']:.1e}")
    ax.axvline(c5["Phi0_min"], color="orange", ls=":", lw=2,
               label=rf"$\Phi_0^{{\min}} = {c5['Phi0_min']:.2e}$")
    ax.axvspan(c5["Phi0_min"], Phi0_arr[-1],
               color="green", alpha=0.08, label="Allowed")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$|\gamma_{\rm PPN} - 1|$", fontsize=12)
    ax.set_title("5. PPN / Cassini", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(1, 80)
    ax.grid(True, ls=":", alpha=0.4)

    # ----- Panel (2,1): Galaxy clusters -----
    ax = axes[2, 1]
    ax.loglog(Phi0_arr, c6["delta_G"], "b-", lw=2.2,
              label=r"$|G_{\rm eff}/G_0 - 1|$ (TGP)")
    ax.axhline(c6["cluster_bound"], color="red", ls="--", lw=1.5,
               label=f"Cluster bound: {c6['cluster_bound']:.0%}")
    ax.axvline(c6["Phi0_min"], color="orange", ls=":", lw=2,
               label=rf"$\Phi_0^{{\min}} = {c6['Phi0_min']:.2f}$")
    ax.axvspan(c6["Phi0_min"], Phi0_arr[-1],
               color="green", alpha=0.08, label="Allowed")
    ax.set_xlabel(r"$\Phi_0$", fontsize=12)
    ax.set_ylabel(r"$|G_{\rm eff}/G_0 - 1|$", fontsize=12)
    ax.set_title("6. Galaxy clusters (X-ray/lensing)", fontsize=13,
                 fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(1, 80)
    ax.grid(True, ls=":", alpha=0.4)

    # ----- Super-title -----
    fig.suptitle(
        r"TGP: Cross-verification of $\Phi_0$ from 6 independent constraints"
        "\n"
        rf"Consistent range: $\Phi_0 \in [{overall_lo:.1f},\; {overall_hi:.1f}]$",
        fontsize=15, y=1.01,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=180, bbox_inches="tight")
    print(f"  Plot saved to {save_path}")
    plt.close(fig)


# =====================================================================
# Console summary
# =====================================================================
def print_summary(Phi0_arr):
    c1, c2, c3, c4, c5, c6 = compute_all(Phi0_arr)
    overall_lo, overall_hi, lower_bounds = find_consistent_range(
        c1, c2, c3, c4, c5, c6
    )

    sep = "=" * 68
    print(sep)
    print("  TGP: Cross-verification of Phi_0")
    print(sep)

    print(f"\n  Physical constants (scipy.constants):")
    print(f"    c0    = {c0:.6e} m/s")
    print(f"    G0    = {G0:.4e} m^3 kg^-1 s^-2")
    print(f"    hbar0 = {hbar0:.6e} J s")
    print(f"    H0    = {H0_SI:.4e} s^-1  ({H0_km:.0f} km/s/Mpc)")
    print(f"    l_P   = {l_P:.4e} m")

    print(f"\n  Observed:")
    print(f"    Lambda_obs   = {Lambda_obs:.4e} m^-2")
    print(f"    rho_Lambda   = {rho_Lambda_obs:.4e} kg/m^3")
    print(f"    Omega_Lambda = {Omega_Lambda}")

    print(f"\n  --- Constraint 1: Dark energy ---")
    print(f"    Lambda_eff = gamma/12 = Phi0*H0^2/(12*c0^2)")
    print(f"    Matching Lambda_obs:")
    print(f"      Phi0 = {c1['Phi0_central']:.2f}")
    print(f"      1-sigma band: [{c1['Phi0_lo']:.2f}, {c1['Phi0_hi']:.2f}]")

    print(f"\n  --- Constraint 2: BBN (frozen field) ---")
    print(f"    tau0*H_BBN > {c2['N_freeze']:.0f}")
    print(f"    Lower bound: Phi0 > {c2['Phi0_min']:.4e}")
    print(f"    (trivially satisfied for Phi0 ~ O(10))")

    print(f"\n  --- Constraint 3: CMB (frozen at recombination) ---")
    print(f"    tau0*H_rec > 10")
    print(f"    Lower bound: Phi0 > {c3['Phi0_min']:.4e}")
    print(f"    (trivially satisfied for Phi0 ~ O(10))")

    print(f"\n  --- Constraint 4: GW speed (GW170817) ---")
    print(f"    |delta c/c| < {c4['gw_bound']:.0e}")
    print(f"    Formal lower bound: Phi0 > {c4['Phi0_min']:.4e}")
    print(f"    (|delta c/c| ~ (l_P*H0)^2/Phi0 ~ 10^{{-122}}/Phi0)")
    print(f"    Trivially satisfied for any Phi0.")

    print(f"\n  --- Constraint 5: PPN / Cassini ---")
    print(f"    |gamma_PPN - 1| < {c5['cassini_bound']:.1e}")
    print(f"    Lower bound: Phi0 > {c5['Phi0_min']:.4e}")
    idx25 = np.argmin(np.abs(Phi0_arr - c1["Phi0_central"]))
    print(f"    At Phi0={c1['Phi0_central']:.1f}: "
          f"|gamma_PPN-1| = {c5['gamma_ppn_dev'][idx25]:.2e}")

    print(f"\n  --- Constraint 6: Galaxy clusters ---")
    print(f"    |G_eff/G0 - 1| < {c6['cluster_bound']:.0%}")
    print(f"    Lower bound: Phi0 > {c6['Phi0_min']:.2f}")
    print(f"    At Phi0={c1['Phi0_central']:.1f}: "
          f"|delta G| = {c6['delta_G'][idx25]:.2e}")

    print(f"\n  {'=' * 50}")
    print(f"  CONSISTENT Phi_0 RANGE")
    print(f"  {'=' * 50}")
    print(f"    Lower bounds:")
    for name, val in sorted(lower_bounds.items(), key=lambda x: -x[1]):
        tag = " <-- binding" if val >= overall_lo * 0.99 else ""
        print(f"      {name:12s}: Phi0 > {val:.4e}{tag}")
    print(f"    Dark-energy band: [{c1['Phi0_lo']:.2f}, {c1['Phi0_hi']:.2f}]")
    print(f"\n    => Overall consistent range:")
    print(f"       Phi_0 in [{overall_lo:.2f}, {overall_hi:.2f}]")
    print(f"       Central value: Phi_0 = {c1['Phi0_central']:.2f}")
    print(f"\n  Key result: Phi_0 ~ O(10) is fixed by ONE observable")
    print(f"  (dark energy density) and is compatible with all other")
    print(f"  constraints.  No fine-tuning needed (no hierarchy problem).")
    print(sep)


# =====================================================================
# Main
# =====================================================================
def main():
    Phi0_arr = np.linspace(1.0, 80.0, 2000)

    # Output path: TGP/TGP_v1/plots/phi0_cross_verification.png
    save_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots', "phi0_cross_verification.png")

    print_summary(Phi0_arr)
    make_plot(Phi0_arr, save_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
