"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
three_regimes_quantitative.py  --  Theory of Generated Space (TGP), O18
========================================================================
Quantitative predictions for the three interaction regimes in TGP.

Maps the dimensionless TGP potential to physical (SI) units and computes
crossover radii, well depths, forces, equilibrium points, and sensitivity
to the beta/gamma ratio.

Physics  (sek03, sek08)
-----------------------
The effective potential between two masses M1, M2 at separation d is:

    V_eff(d) = E_lin(d) + E_beta(d) + E_gamma(d)

where (eq:Eint-decomp, sek08):

    E_lin(d)   = -(qPhi0)^2 M1 M2 / (4 pi d)                [attractive, Newton-like, 1/d]
    E_beta(d)  = 2 beta (qPhi0)^2 M1 M2 / ((4pi)^2 Phi0 d)
                 * ln(d / r0)                                 [repulsive, logarithmic]
    E_gamma(d) = gamma (qPhi0)^2 M1 M2 / ((4pi)^3 Phi0^2 d^3)
                 * K_gamma                                    [short-range well, 1/d^3]

In the overlap-integral (power-law) approximation used in the companion
scripts effective_potential.py and two_body_three_regimes.py, the leading
terms reduce to:

    E_lin   = -4 pi C^2 / d        (sek08, gradient overlap)
    E_beta  = +8 pi beta_hat C^2 / d^2   (cubic self-interaction)
    E_gamma = -24 pi gamma_hat C^3 / d^3  (quartic self-interaction)

with C = alpha_eff M / (4 pi r0), beta_hat = beta r0^2, gamma_hat = gamma r0^2.

Three regimes exist when beta_hat > 9 C / 2  (Proposition prop:trzy-rezimy-beta-gamma).

Parameters
----------
alpha = 2                            (fixed, gradient nonlinearity)
beta = gamma                         (vacuum condition, sek03)
q = 8 pi G0 / c0^2                  (source coupling, [m/kg])
Phi0 ~ 25                           (from Lambda_eff matching, sek08)
gamma ~ Phi0 H0^2 / c0^2            (from tau0 ~ 1/H0, sek08)
Lambda_eff = gamma / 12             (residual vacuum energy, sek08)

Outputs (saved to tooling/scripts/plots/):
    three_regimes_quantitative.png    -- 4-panel figure: V_eff, force, sensitivity
    three_regimes_force.png           -- standalone force plot with equilibria

Author: TGP consistency analysis
"""

import os
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# 1. Fundamental physical constants (SI)
#    (sek08: these define q, and through H0 they fix gamma)
# ============================================================================

G0       = 6.67430e-11      # m^3 kg^-1 s^-2   Newton constant
c0       = 2.99792458e8     # m s^-1            speed of light
hbar0    = 1.054571817e-34  # J s               reduced Planck constant
H0_si    = 2.268e-18        # s^-1              Hubble parameter (~70 km/s/Mpc)
Omega_Lambda = 0.6889       # Planck 2018       dark energy fraction

# Observed cosmological constant  (m^-2)
Lambda_obs = 3.0 * Omega_Lambda * H0_si**2 / c0**2

# Hubble radius  (m)
R_Hubble = c0 / H0_si

# ============================================================================
# 2. TGP parameters from fundamental constants  (sek03, sek08)
# ============================================================================

ALPHA = 2.0            # gradient nonlinearity exponent (fixed by the action)

# Source coupling  q = 8 pi G0 / c0^2   [m / kg]   (sek03, eq:field-eq-rhs)
q_coupling = 8.0 * np.pi * G0 / c0**2

# Background field value.  Lambda_eff = gamma/12 = Lambda_obs
# gives  gamma = 12 Lambda_obs.  Combined with
# gamma = Phi0 H0^2 / c0^2  =>  Phi0 = 12 Lambda_obs c0^2 / H0^2
#                                      = 12 * 3 Omega_Lambda = 36 Omega_Lambda ~ 24.8
PHI0 = 12.0 * Lambda_obs * c0**2 / H0_si**2
# (Should be ~ 25)

# gamma and beta from the naturalness condition  (sek08)
# gamma = Phi0 * H0^2 / c0^2   [m^-2]
GAMMA_PHYS = PHI0 * H0_si**2 / c0**2    # m^-2
BETA_PHYS  = GAMMA_PHYS                  # vacuum condition beta = gamma

# Effective cosmological constant from TGP
Lambda_eff_tgp = GAMMA_PHYS / 12.0

# Characteristic length scale r0: Compton wavelength of the spatial scalar
# m_sp = sqrt(gamma)  =>  r0 ~ 1/m_sp = 1/sqrt(gamma)   [m]
m_sp = np.sqrt(GAMMA_PHYS)       # spatial mass [m^-1]
r0   = 1.0 / m_sp                # characteristic length [m]

# Dimensionless couplings
BETA_HAT  = BETA_PHYS * r0**2     # = 1.0 by construction (beta = gamma, r0 = 1/sqrt(gamma))
GAMMA_HAT = GAMMA_PHYS * r0**2    # = 1.0

# Known physical scales for comparison
r_Bohr    = 5.292e-11     # m   Bohr radius
r_proton  = 8.75e-16      # m   proton charge radius
r_Planck  = np.sqrt(hbar0 * G0 / c0**3)   # m   Planck length

# ============================================================================
# 3. Dimensionless source strength C for various objects
#    C = alpha_eff * M / (4 pi r0)  where alpha_eff = q * Phi0
#    (sek08, eq:C-definition)
# ============================================================================

alpha_eff = q_coupling * PHI0  # effective coupling [m / kg] * [1] = [m / kg]


def source_strength_C(M):
    """
    Dimensionless source strength for a mass M.
    C = alpha_eff * M / (4 pi r0)     (sek08)
    """
    return alpha_eff * M / (4.0 * np.pi * r0)


# Particle masses  (kg)
m_electron = 9.1094e-31
m_proton   = 1.6726e-27
m_sun      = 1.989e30
m_earth    = 5.972e24

# ============================================================================
# 4. Analytical interaction potential (dimensionless units)
#    (sek08, Theorem thm:three-regimes)
# ============================================================================

def V_eff_dimless(d, C, beta_hat=BETA_HAT, gamma_hat=GAMMA_HAT):
    """
    Leading-order analytical interaction energy in dimensionless units
    (d measured in units of r0).

    V_eff(d) = E_lin + E_beta + E_gamma

    E_lin   = -4 pi C^2 / d                (gradient overlap, attractive)
    E_beta  = +8 pi beta_hat C^2 / d^2     (cubic self-interaction, repulsive)
    E_gamma = -24 pi gamma_hat C^3 / d^3   (quartic self-interaction, confining)

    (sek08, eq:Eint-decomp; Proposition prop:trzy-rezimy-beta-gamma, sek03)

    Returns: (V_total, E_lin, E_beta, E_gamma)
    """
    d_safe = np.maximum(np.asarray(d, dtype=float), 1e-30)
    E_lin   = -4.0 * np.pi * C**2 / d_safe
    E_beta  =  8.0 * np.pi * beta_hat * C**2 / d_safe**2
    E_gamma = -24.0 * np.pi * gamma_hat * C**3 / d_safe**3
    return E_lin + E_beta + E_gamma, E_lin, E_beta, E_gamma


def force_dimless(d, C, beta_hat=BETA_HAT, gamma_hat=GAMMA_HAT):
    """
    Force F = -dV_eff/dd in dimensionless units.

    F(d) = -4 pi C^2 / d^2  +  16 pi beta_hat C^2 / d^3
           - 72 pi gamma_hat C^3 / d^4

    (sek08, derivative of eq:Eint-decomp)
    """
    d_safe = np.maximum(np.asarray(d, dtype=float), 1e-30)
    return (-4.0 * np.pi * C**2 / d_safe**2
            + 16.0 * np.pi * beta_hat * C**2 / d_safe**3
            - 72.0 * np.pi * gamma_hat * C**3 / d_safe**4)


# ============================================================================
# 5. Crossover radii: d1 (Newton -> repulsion) and d2 (repulsion -> well)
#    (sek03, Proposition prop:trzy-rezimy-beta-gamma)
# ============================================================================

def find_crossover_radii(C, beta_hat=BETA_HAT, gamma_hat=GAMMA_HAT):
    """
    Find the two force-zero crossover radii from F(d) = 0.

    The force zeros satisfy (sek03):
        d^2 - 4 beta_hat d + 18 gamma_hat C = 0      [for beta_hat = gamma_hat]

    Solutions:
        d_2 = 2 beta_hat - sqrt(4 beta_hat^2 - 18 gamma_hat C)    (inner, well -> repulsion)
        d_1 = 2 beta_hat + sqrt(4 beta_hat^2 - 18 gamma_hat C)    (outer, repulsion -> gravity)

    Three regimes exist iff discriminant >= 0:
        4 beta_hat^2 - 18 gamma_hat C >= 0
        => beta_hat >= (9/2) * gamma_hat * C / beta_hat  =>  beta_hat^2 >= (9/2) gamma_hat C

    For beta_hat = gamma_hat this reduces to:  beta_hat > 9 C / 2.

    Returns
    -------
    (d_2, d_1) in dimensionless units (units of r0), or None if only 1 regime.
    d_2 < d_1: d_2 is the inner transition (short-range well boundary),
               d_1 is the outer transition (Newton gravity boundary).
    """
    disc = 4.0 * beta_hat**2 - 18.0 * gamma_hat * C
    if disc < 0:
        return None

    # For very small C (elementary particles), direct subtraction
    # 2*beta_hat - sqrt(4*beta_hat^2 - 18*gamma_hat*C) loses precision.
    # Use Taylor expansion: d_inner ~ 9*gamma_hat*C / (2*beta_hat)
    # when 18*gamma_hat*C << 4*beta_hat^2.
    ratio = 18.0 * gamma_hat * C / (4.0 * beta_hat**2)
    if ratio < 1e-10:
        # Taylor regime: avoid catastrophic cancellation
        d_inner = 9.0 * gamma_hat * C / (2.0 * beta_hat)
        d_outer = 4.0 * beta_hat - d_inner   # ~ 4*beta_hat for small C
    else:
        sqrt_disc = np.sqrt(disc)
        d_inner = 2.0 * beta_hat - sqrt_disc
        d_outer = 2.0 * beta_hat + sqrt_disc

    if d_inner <= 0:
        return None
    return d_inner, d_outer


def crossover_to_physical(d_dimless):
    """Convert dimensionless distance to physical distance [m]."""
    return d_dimless * r0


# ============================================================================
# 6. Well depth and equilibrium analysis
# ============================================================================

def well_depth(C, beta_hat=BETA_HAT, gamma_hat=GAMMA_HAT):
    """
    Compute the depth of the confining well at d_2 (inner equilibrium).

    Returns (d_min, V_min) in dimensionless units, or None if no well exists.
    """
    radii = find_crossover_radii(C, beta_hat, gamma_hat)
    if radii is None:
        return None
    d_inner, d_outer = radii

    # The local minimum of V_eff is between 0 and d_inner.
    # Find it by looking for F(d) = 0 with F' > 0 (stable equilibrium).
    # Actually, d_inner IS a force zero. Check sign of V_eff there and
    # find the actual minimum of V_eff between some small d and d_inner.
    d_scan = np.linspace(0.01, d_inner, 5000)
    V_scan = V_eff_dimless(d_scan, C, beta_hat, gamma_hat)[0]
    idx_min = np.argmin(V_scan)
    return d_scan[idx_min], V_scan[idx_min]


def find_equilibria(C, beta_hat=BETA_HAT, gamma_hat=GAMMA_HAT):
    """
    Find equilibrium points where F(d) = 0 and classify their stability.

    Returns list of (d_eq, stability) where stability is 'stable' or 'unstable'.
    """
    radii = find_crossover_radii(C, beta_hat, gamma_hat)
    if radii is None:
        return []

    d_inner, d_outer = radii
    equilibria = []

    # Check stability: dF/dd at the equilibrium point
    # Stable if dF/dd < 0 (restoring force); unstable if dF/dd > 0
    for d_eq in [d_inner, d_outer]:
        # Central finite difference for dF/dd
        dx = 1e-6
        f_plus  = force_dimless(d_eq + dx, C, beta_hat, gamma_hat)
        f_minus = force_dimless(d_eq - dx, C, beta_hat, gamma_hat)
        dF = (f_plus - f_minus) / (2.0 * dx)
        stability = 'stable' if dF < 0 else 'unstable'
        equilibria.append((d_eq, stability))

    return equilibria


# ============================================================================
# 7. Sensitivity analysis: how d1, d2 change with beta/gamma ratio
#    (sek03: vacuum condition is beta = gamma; explore deviations)
# ============================================================================

def sensitivity_analysis(C, n_ratios=200):
    """
    Scan the beta/gamma ratio from 0.5 to 2.0 and track crossover radii.

    At beta/gamma = 1 (vacuum condition), the canonical result holds.
    Deviations correspond to non-vacuum / cosmological density offset.

    Returns
    -------
    ratios : array of beta/gamma values
    d_inner_arr : array of d_2 values (NaN where no 3 regimes)
    d_outer_arr : array of d_1 values (NaN where no 3 regimes)
    """
    ratios = np.linspace(0.5, 2.0, n_ratios)
    d_inner_arr = np.full(n_ratios, np.nan)
    d_outer_arr = np.full(n_ratios, np.nan)

    for i, ratio in enumerate(ratios):
        bh = ratio * GAMMA_HAT
        gh = GAMMA_HAT
        result = find_crossover_radii(C, bh, gh)
        if result is not None:
            d_inner_arr[i], d_outer_arr[i] = result

    return ratios, d_inner_arr, d_outer_arr


# ============================================================================
# 8. Plotting
# ============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR   = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(PLOT_DIR, exist_ok=True)


def format_si(value, unit="m"):
    """Format a value in SI with appropriate prefix."""
    if value == 0:
        return f"0 {unit}"
    exp = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exp
    return f"{mantissa:.3f}e{exp:+d} {unit}"


def make_main_figure():
    """
    4-panel production figure:

    (a) V_eff(d) for an elementary particle showing all 3 regimes
    (b) Force F(d) with equilibrium points annotated
    (c) V_eff(d) for several C values (particle -> macroscopic)
    (d) Sensitivity of crossover radii to beta/gamma ratio
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 13))
    fig.suptitle(
        "TGP (O18): Quantitative Three-Regime Predictions\n"
        r"$\alpha = 2$, $\beta = \gamma$ (vacuum), "
        rf"$\Phi_0 = {PHI0:.2f}$, "
        rf"$\gamma = {GAMMA_PHYS:.4e}$ m$^{{-2}}$, "
        rf"$r_0 = {r0:.4e}$ m",
        fontsize=13, fontweight='bold', y=0.99
    )

    # --- Use a representative elementary particle: proton ---
    C_proton = source_strength_C(m_proton)

    # ------------------------------------------------------------------
    # Panel (a): V_eff(d) decomposition for proton, 3 regimes
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    radii_p = find_crossover_radii(C_proton)

    if radii_p is not None:
        d_inner, d_outer = radii_p
        d_max_plot = min(4.0 * d_outer, 1e4)
    else:
        d_inner, d_outer = 1.0, 5.0
        d_max_plot = 30.0

    d_arr = np.linspace(0.01, d_max_plot, 2000)
    V_tot, E_lin, E_beta, E_gamma = V_eff_dimless(d_arr, C_proton)

    # Normalise for visibility
    finite = np.isfinite(V_tot) & (np.abs(V_tot) < 1e30)
    Vsc = np.max(np.abs(V_tot[finite])) if np.any(finite) else 1.0
    Vsc = max(Vsc, 1e-50)

    ax.plot(d_arr, V_tot / Vsc, 'k-', lw=2.2, label=r"$V_{\rm eff}$ (total)")
    ax.plot(d_arr, E_lin / Vsc, 'b--', lw=1.3,
            label=r"$E_{\rm lin} \propto -1/d$ (gravity)")
    ax.plot(d_arr, E_beta / Vsc, 'r-.', lw=1.3,
            label=r"$E_\beta \propto +1/d^2$ (repulsion)")
    ax.plot(d_arr, E_gamma / Vsc, 'm:', lw=1.5,
            label=r"$E_\gamma \propto -1/d^3$ (confinement)")
    ax.axhline(0, color='k', lw=0.5)

    if radii_p is not None:
        ax.axvline(d_inner, color='green', ls=':', lw=1.2, alpha=0.7,
                   label=rf"$d_2 = {d_inner:.4f}\,r_0$ (well boundary)")
        ax.axvline(d_outer, color='orange', ls=':', lw=1.2, alpha=0.7,
                   label=rf"$d_1 = {d_outer:.4f}\,r_0$ (gravity onset)")
        # Shade regimes
        ax.axvspan(0, d_inner, alpha=0.04, color='purple')
        ax.axvspan(d_inner, d_outer, alpha=0.04, color='red')
        ax.axvspan(d_outer, d_max_plot, alpha=0.04, color='blue')
        # Labels
        x_well = d_inner * 0.4 if d_inner > 0.1 else 0.05
        ax.text(x_well, -0.6, "III\nwell", fontsize=9, ha='center',
                color='purple', style='italic')
        ax.text((d_inner + d_outer) / 2, 0.5, "II\nrepulsion", fontsize=9,
                ha='center', color='red', style='italic')
        ax.text(d_outer * 1.3, -0.3, "I\ngravity", fontsize=9,
                ha='center', color='blue', style='italic')

    ax.set_xlabel(r"Separation $d / r_0$", fontsize=12)
    ax.set_ylabel(r"Energy / $|V|_{\max}$", fontsize=12)
    ax.set_title(
        rf"(a) Potential decomposition: proton ($C = {C_proton:.4e}$)"
        "\n" + r"$\hat\beta > 9C/2$ satisfied $\Rightarrow$ 3 regimes",
        fontsize=11)
    ax.legend(fontsize=7, loc='best')
    ax.grid(True, ls=':', alpha=0.3)
    ax.set_ylim(-1.5, 1.5)

    # ------------------------------------------------------------------
    # Panel (b): Force F(d) and equilibrium points
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    F_arr = force_dimless(d_arr, C_proton)
    F_finite = np.isfinite(F_arr) & (np.abs(F_arr) < 1e30)
    Fsc = np.max(np.abs(F_arr[F_finite])) if np.any(F_finite) else 1.0
    Fsc = max(Fsc, 1e-50)

    ax.plot(d_arr, F_arr / Fsc, 'k-', lw=2, label=r"$F = -dV_{\rm eff}/dd$")
    ax.axhline(0, color='gray', lw=0.8)

    equilibria = find_equilibria(C_proton)
    for d_eq, stab in equilibria:
        marker = 'o' if stab == 'stable' else 's'
        color = 'green' if stab == 'stable' else 'red'
        d_phys = crossover_to_physical(d_eq)
        ax.axvline(d_eq, color=color, ls=':', lw=1, alpha=0.6)
        ax.plot(d_eq, 0, marker, color=color, ms=10, zorder=5,
                label=f"{stab}: $d = {d_eq:.4f}\\,r_0$ = {format_si(d_phys)}")

    ax.set_xlabel(r"Separation $d / r_0$", fontsize=12)
    ax.set_ylabel(r"$F / |F|_{\max}$", fontsize=12)
    ax.set_title(
        r"(b) Force $F(d)$ and equilibrium points" "\n"
        rf"proton, $C = {C_proton:.4e}$",
        fontsize=11)
    ax.legend(fontsize=7, loc='best')
    ax.grid(True, ls=':', alpha=0.3)
    ax.set_ylim(-1.5, 1.5)

    # ------------------------------------------------------------------
    # Panel (c): V_eff for several C values (log scale in C)
    # ------------------------------------------------------------------
    ax = axes[1, 0]

    C_values = [
        (source_strength_C(m_electron), r"electron", '#1f77b4'),
        (source_strength_C(m_proton),   r"proton",   '#ff7f0e'),
        (0.01,                          r"$C=0.01$", '#2ca02c'),
        (0.1,                           r"$C=0.1$",  '#d62728'),
        (1.0,                           r"$C=1.0$",  '#9467bd'),
        (10.0,                          r"$C=10$",   '#8c564b'),
    ]

    for C_val, name, color in C_values:
        radii_c = find_crossover_radii(C_val)
        satisfied = BETA_HAT > 4.5 * C_val

        if radii_c is not None:
            d_plot_max = min(4.0 * radii_c[1], 1e5)
        else:
            d_plot_max = max(30.0, 20.0 * C_val)
        d_c = np.linspace(0.01, d_plot_max, 1000)
        V_c = V_eff_dimless(d_c, C_val)[0]

        V_fin = V_c[np.isfinite(V_c) & (np.abs(V_c) < 1e30)]
        if len(V_fin) == 0:
            continue
        V_sc_c = np.max(np.abs(V_fin))
        V_sc_c = max(V_sc_c, 1e-50)

        n_reg = 3 if satisfied else 1
        ax.plot(d_c, V_c / V_sc_c, color=color, lw=1.5,
                label=f"{name}: {n_reg} regime(s)")

    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel(r"Separation $d / r_0$", fontsize=12)
    ax.set_ylabel(r"$V_{\rm eff} / |V|_{\max}$ (each normalised)", fontsize=12)
    ax.set_title(
        "(c) Potential across mass scales\n"
        r"3 regimes for $C < C_{\rm crit}$, 1 regime for $C > C_{\rm crit}$",
        fontsize=11)
    ax.legend(fontsize=7, loc='best')
    ax.grid(True, ls=':', alpha=0.3)
    ax.set_ylim(-1.5, 1.5)

    # ------------------------------------------------------------------
    # Panel (d): Sensitivity analysis -- crossover radii vs beta/gamma
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    C_sens = 0.01   # representative small C (elementary particle regime)
    ratios, d_in, d_out = sensitivity_analysis(C_sens, n_ratios=300)

    has_data = np.isfinite(d_in)
    if np.any(has_data):
        ax.plot(ratios[has_data], d_in[has_data], 'r-', lw=2,
                label=r"$d_2$ (well boundary)")
        ax.plot(ratios[has_data], d_out[has_data], 'b-', lw=2,
                label=r"$d_1$ (gravity onset)")
        ax.fill_between(ratios[has_data], d_in[has_data], d_out[has_data],
                        alpha=0.1, color='red', label='Regime II (repulsion)')

    ax.axvline(1.0, color='green', ls='--', lw=1.5,
               label=r"$\beta/\gamma = 1$ (vacuum)")

    # Also show sensitivity for C = 0.1
    _, d_in2, d_out2 = sensitivity_analysis(0.1, n_ratios=300)
    has2 = np.isfinite(d_in2)
    if np.any(has2):
        ax.plot(ratios[has2], d_in2[has2], 'r--', lw=1.2, alpha=0.6)
        ax.plot(ratios[has2], d_out2[has2], 'b--', lw=1.2, alpha=0.6,
                label=r"$C=0.1$ (dashed)")

    ax.set_xlabel(r"$\beta / \gamma$ ratio", fontsize=12)
    ax.set_ylabel(r"Crossover radius $d / r_0$", fontsize=12)
    ax.set_title(
        r"(d) Sensitivity: crossover radii vs $\beta/\gamma$" "\n"
        rf"Solid: $C = {C_sens}$, Dashed: $C = 0.1$",
        fontsize=11)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, ls=':', alpha=0.3)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    outpath = os.path.join(PLOT_DIR, "three_regimes_quantitative.png")
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {outpath}")
    return outpath


def make_force_figure():
    """
    Standalone force plot with equilibrium markers for representative C values.
    Uses dimensionless C values (not real particle masses) so that the
    three-regime structure is clearly visible on the plot.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    # Use C values that show clear three-regime structure
    C_force_values = [
        (0.001, "$C = 0.001$", '#1f77b4'),
        (0.01,  "$C = 0.01$",  '#ff7f0e'),
        (0.05,  "$C = 0.05$",  '#2ca02c'),
        (0.1,   "$C = 0.1$",   '#d62728'),
    ]

    has_any = False
    for C_val, name, color in C_force_values:
        radii = find_crossover_radii(C_val)
        if radii is None:
            continue
        has_any = True
        d_inner, d_outer = radii
        d_arr = np.linspace(0.001, 4.0 * d_outer, 3000)
        F = force_dimless(d_arr, C_val)

        F_fin = F[np.isfinite(F) & (np.abs(F) < 1e30)]
        Fsc = np.max(np.abs(F_fin)) if len(F_fin) > 0 else 1.0

        ax.plot(d_arr, F / Fsc, color=color, lw=2,
                label=f"{name}")

        # Mark equilibria
        for d_eq, stab in find_equilibria(C_val):
            marker = 'o' if stab == 'stable' else 's'
            mcolor = 'green' if stab == 'stable' else 'red'
            ax.plot(d_eq, 0, marker, color=mcolor, ms=8, zorder=5)

    ax.axhline(0, color='gray', lw=0.8)
    ax.set_xlabel(r"Separation $d / r_0$", fontsize=13)
    ax.set_ylabel(r"$F / |F|_{\max}$", fontsize=13)
    ax.set_title(
        "TGP: Force for representative source strengths\n"
        r"$F = -dV_{\rm eff}/dd$, equilibria: o = stable, s = unstable",
        fontsize=13)
    if has_any:
        ax.legend(fontsize=10)
    ax.grid(True, ls=':', alpha=0.4)
    ax.set_ylim(-1.5, 1.5)

    fig.tight_layout()
    outpath = os.path.join(PLOT_DIR, "three_regimes_force.png")
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {outpath}")
    return outpath


# ============================================================================
# 9. Print quantitative predictions table
# ============================================================================

def print_predictions_table():
    """
    Print a comprehensive table of quantitative predictions from TGP O18.
    """
    sep = "=" * 80
    print(f"\n{sep}")
    print("  TGP (O18) QUANTITATIVE PREDICTIONS: THREE INTERACTION REGIMES")
    print(sep)

    # --- Fundamental parameters ---
    print("\n  FUNDAMENTAL CONSTANTS (SI):")
    print(f"    G0             = {G0:.5e} m^3 kg^-1 s^-2")
    print(f"    c0             = {c0:.8e} m/s")
    print(f"    hbar           = {hbar0:.6e} J s")
    print(f"    H0             = {H0_si:.4e} s^-1  (~70 km/s/Mpc)")
    print(f"    Omega_Lambda   = {Omega_Lambda:.4f}")
    print(f"    Lambda_obs     = {Lambda_obs:.6e} m^-2")

    print("\n  TGP PARAMETERS:")
    print(f"    alpha          = {ALPHA}")
    print(f"    Phi0           = {PHI0:.4f}  (from Lambda_eff = Lambda_obs)")
    print(f"    q = 8piG/c^2   = {q_coupling:.6e} m/kg")
    print(f"    gamma          = {GAMMA_PHYS:.6e} m^-2")
    print(f"    beta           = {BETA_PHYS:.6e} m^-2  (= gamma, vacuum condition)")
    print(f"    Lambda_eff     = gamma/12 = {Lambda_eff_tgp:.6e} m^-2")
    print(f"    Lambda_eff/Lambda_obs = {Lambda_eff_tgp / Lambda_obs:.6f}")
    print(f"    m_sp = sqrt(gamma)     = {m_sp:.6e} m^-1")
    print(f"    r0 = 1/m_sp            = {r0:.6e} m")
    print(f"    r0 / R_Hubble          = {r0 / R_Hubble:.4f}")
    print(f"    beta_hat = beta*r0^2   = {BETA_HAT:.6f}")
    print(f"    gamma_hat = gamma*r0^2 = {GAMMA_HAT:.6f}")

    # --- Source strengths ---
    print(f"\n  SOURCE STRENGTHS  C = alpha_eff * M / (4 pi r0):")
    objects = [
        ("electron",  m_electron),
        ("proton",    m_proton),
        ("Earth",     m_earth),
        ("Sun",       m_sun),
    ]
    for name, mass in objects:
        C = source_strength_C(mass)
        satisfied = BETA_HAT > 4.5 * C
        ratio = BETA_HAT / (4.5 * C) if C > 0 else np.inf
        print(f"    {name:12s}: M = {mass:.4e} kg, C = {C:.6e}, "
              f"beta_hat/(9C/2) = {ratio:.4e}  "
              f"{'-> 3 regimes' if satisfied else '-> 1 regime (gravity only)'}")

    # --- Crossover radii ---
    print(f"\n  CROSSOVER RADII (d_2: well<->repulsion, d_1: repulsion<->gravity):")
    print(f"  {'Object':>12s}  {'C':>14s}  {'d_2 [r0]':>14s}  {'d_2 [m]':>16s}  "
          f"{'d_1 [r0]':>14s}  {'d_1 [m]':>16s}  {'Well depth':>14s}")
    print(f"  {'-'*106}")

    for name, mass in objects:
        C = source_strength_C(mass)
        radii = find_crossover_radii(C)
        if radii is not None:
            d_in, d_out = radii
            d_in_m = crossover_to_physical(d_in)
            d_out_m = crossover_to_physical(d_out)
            wd = well_depth(C)
            wd_str = f"{wd[1]:.6e}" if wd is not None else "---"
            print(f"  {name:>12s}  {C:>14.6e}  {d_in:>14.6e}  {d_in_m:>16.6e}  "
                  f"{d_out:>14.6e}  {d_out_m:>16.6e}  {wd_str:>14s}")
        else:
            print(f"  {name:>12s}  {C:>14.6e}  {'---':>14s}  {'---':>16s}  "
                  f"{'---':>14s}  {'---':>16s}  {'---':>14s}")

    # --- Comparison with known scales ---
    print(f"\n  COMPARISON WITH KNOWN PHYSICAL SCALES:")
    print(f"    Planck length   l_P        = {r_Planck:.6e} m")
    print(f"    Proton radius   r_p        = {r_proton:.6e} m")
    print(f"    Bohr radius     a_0        = {r_Bohr:.6e} m")
    print(f"    TGP length      r_0        = {r0:.6e} m")
    print(f"    Hubble radius   R_H        = {R_Hubble:.6e} m")

    C_proton = source_strength_C(m_proton)
    radii_p = find_crossover_radii(C_proton)
    if radii_p is not None:
        d2_m = crossover_to_physical(radii_p[0])
        d1_m = crossover_to_physical(radii_p[1])
        print(f"\n    Proton d_2 (well) = {d2_m:.6e} m")
        print(f"      d_2 / l_Planck  = {d2_m / r_Planck:.4e}")
        print(f"      d_2 / r_proton  = {d2_m / r_proton:.4e}")
        print(f"    Proton d_1 (grav) = {d1_m:.6e} m")
        print(f"      d_1 / a_Bohr    = {d1_m / r_Bohr:.4e}")
        print(f"      d_1 / R_Hubble  = {d1_m / R_Hubble:.4e}")

    # --- Critical source strength ---
    C_crit = 2.0 * BETA_HAT / 9.0
    M_crit = C_crit * 4.0 * np.pi * r0 / alpha_eff
    print(f"\n  CRITICAL SOURCE STRENGTH:")
    print(f"    C_crit = 2 beta_hat / 9 = {C_crit:.6e}")
    print(f"    M_crit = C_crit * 4 pi r0 / alpha_eff = {M_crit:.6e} kg")
    print(f"    M_crit / m_proton = {M_crit / m_proton:.4e}")
    print(f"    M_crit / m_sun    = {M_crit / m_sun:.4e}")
    print(f"    Objects with M < M_crit exhibit three regimes.")
    print(f"    Objects with M > M_crit show only Newtonian gravity.")

    # --- Equilibria ---
    print(f"\n  EQUILIBRIUM POINTS (proton):")
    eqs = find_equilibria(C_proton)
    for d_eq, stab in eqs:
        d_phys = crossover_to_physical(d_eq)
        print(f"    d = {d_eq:.6e} r0 = {d_phys:.6e} m  ({stab})")

    # --- Sensitivity highlights ---
    print(f"\n  SENSITIVITY TO beta/gamma RATIO (C = 0.01):")
    C_sens = 0.01
    for bg_ratio in [0.8, 0.9, 1.0, 1.1, 1.2, 1.5]:
        bh = bg_ratio * GAMMA_HAT
        radii_s = find_crossover_radii(C_sens, bh, GAMMA_HAT)
        if radii_s is not None:
            d_in_s, d_out_s = radii_s
            d_in_phys = crossover_to_physical(d_in_s)
            d_out_phys = crossover_to_physical(d_out_s)
            print(f"    beta/gamma = {bg_ratio:.2f}:  "
                  f"d_2 = {d_in_s:.6f} r0 ({d_in_phys:.4e} m),  "
                  f"d_1 = {d_out_s:.6f} r0 ({d_out_phys:.4e} m),  "
                  f"d_1/d_2 = {d_out_s/d_in_s:.2f}")
        else:
            print(f"    beta/gamma = {bg_ratio:.2f}:  no three regimes")

    print(f"\n{sep}")
    print("  END OF PREDICTIONS TABLE")
    print(sep)


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    print_predictions_table()
    outpath1 = make_main_figure()
    outpath2 = make_force_figure()
    print(f"\nAll done.")
    print(f"  Main plot:  {outpath1}")
    print(f"  Force plot: {outpath2}")
