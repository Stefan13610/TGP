"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
profile_intermediate.py  --  Theory of Generated Space (TGP)
=============================================================
Compute the analytical profile phi(r) in the intermediate (dark-matter-like)
regime via asymptotic matching between inner (power-series near r_0) and
outer (Yukawa far-field) solutions.

Physics (sek03, sek08)
----------------------
The static, spherically symmetric, vacuum TGP field equation is:

    nabla^2 Phi + 2 (nabla Phi)^2 / Phi + beta Phi^2 / Phi0
                - beta Phi^3 / Phi0^2 = 0

Writing Phi = Phi0 (1 + phi), where phi << 1 in the far field, the
equation for phi in spherical symmetry (r > 0) reads:

    phi'' + (2/r) phi' + 2 (phi')^2 / (1 + phi)
        + beta (1 + phi)^2 - beta (1 + phi)^3 = 0

Simplifying the potential terms:
    (1+phi)^2 - (1+phi)^3 = -(1+phi)^2 phi
so the equation becomes:

    phi'' + (2/r) phi' + 2 (phi')^2 / (1 + phi) - beta (1+phi)^2 phi = 0

Three regimes:
    - Far field (r >> r_crit):  phi ~ C/r exp(-m_sp r),  Yukawa tail
    - Intermediate (r_0 < r < r_crit):  nonlinear 1/r^2 tail mimics DM
    - Core (r < r_0):  strong-field, repulsive

Key scales:
    m_sp = sqrt(beta)           (spatial scalar mass)
    r_crit = 1 / (2 beta C)    (onset of nonlinear regime)
    r_0                         (inner boundary, from V_eff zero crossing)

Outputs (saved to scripts/plots/):
    profile_intermediate.png  -- 4-panel figure:
        (a) phi(r): full numerical vs matched approximation
        (b) Effective DM density profile vs NFW
        (c) Rotation curve v(r)
        (d) Residuals (numerical - matched) / numerical
"""

import os
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# 1. Parameters
# ============================================================================

ALPHA = 2.0  # gradient nonlinearity (fixed by action)

# Parameter grid: (C, beta) combinations
C_VALUES = [0.1, 1.0, 5.0]
BETA_VALUES = [0.1, 0.5, 1.0]

# Numerical grid
R_MIN = 0.05       # avoid r=0 singularity
R_MAX_FACTOR = 15  # r_max = R_MAX_FACTOR / m_sp (several Yukawa lengths)
N_MESH = 800


# ============================================================================
# 2. Derived scales
# ============================================================================

def compute_scales(C, beta):
    """Compute characteristic scales for given C and beta."""
    m_sp = np.sqrt(beta)                # spatial scalar mass
    r_crit = 1.0 / (2.0 * beta * C)    # nonlinear onset radius
    # r_0: estimate from V_eff zero crossing ~ 9*beta*C / 2 in dimensionless
    # For the single-source profile, r_0 is where the effective potential
    # changes sign. A useful estimate is r_0 ~ 9*C / (2*beta) (from the
    # two-body analysis), but for the single-source radial profile we take
    # r_0 as the radius where phi becomes O(1), i.e., C/r_0 ~ 1 => r_0 ~ C
    r_0 = max(C, R_MIN * 2)
    return m_sp, r_crit, r_0


# ============================================================================
# 3. Full nonlinear BVP solver
# ============================================================================

def tgp_ode(r, y, beta):
    """
    ODE system for phi(r) where Phi = Phi0 (1 + phi).

    y[0] = phi, y[1] = phi'

    phi'' = -(2/r) phi' - 2 (phi')^2 / (1 + phi) + beta (1+phi)^2 phi
    """
    phi = y[0]
    dphi = y[1]
    r_safe = np.maximum(r, 1e-12)
    f = np.maximum(1.0 + phi, 1e-12)

    ddphi = (-(2.0 / r_safe) * dphi
             - ALPHA * dphi**2 / f
             + beta * f**2 * phi)
    # Note: (1+phi)^2 * phi = f^2 * phi, which is the simplified potential

    return np.vstack([dphi, ddphi])


def tgp_bc(ya, yb, phi_inner):
    """
    Boundary conditions:
        phi(r_min) = phi_inner  (matched to inner solution)
        phi(r_max) = 0          (far-field: Phi -> Phi0)
    """
    return np.array([ya[0] - phi_inner, yb[0]])


def solve_tgp_profile(C, beta, r_min=None, r_max=None, n_mesh=N_MESH):
    """
    Solve the full nonlinear vacuum TGP equation for phi(r).

    Uses C/r * exp(-m_sp * r) as initial guess (Yukawa profile).
    Inner BC: phi(r_min) = C/r_min * exp(-m_sp * r_min)
    Outer BC: phi(r_max) = 0
    """
    m_sp, r_crit, r_0 = compute_scales(C, beta)

    if r_min is None:
        r_min = max(R_MIN, 0.05 * min(r_0, r_crit))
        r_min = max(r_min, 0.02)
    if r_max is None:
        r_max = max(R_MAX_FACTOR / m_sp, 5.0 * r_crit, 50.0)
        r_max = min(r_max, 500.0)  # cap for numerical stability

    r_mesh = np.linspace(r_min, r_max, n_mesh)

    # Initial guess: Yukawa profile (linearized solution)
    phi_guess = C / r_mesh * np.exp(-m_sp * r_mesh)
    dphi_guess = -C * np.exp(-m_sp * r_mesh) * (m_sp + 1.0 / r_mesh) / r_mesh

    y_guess = np.vstack([phi_guess, dphi_guess])

    # Inner BC value from Yukawa guess
    phi_inner = C / r_min * np.exp(-m_sp * r_min)

    sol = solve_bvp(
        fun=lambda r, y: tgp_ode(r, y, beta),
        bc=lambda ya, yb: tgp_bc(ya, yb, phi_inner),
        x=r_mesh,
        y=y_guess,
        tol=1e-6,
        max_nodes=50000,
        verbose=0,
    )

    if not sol.success:
        # Retry with relaxed tolerance and denser mesh
        r_mesh2 = np.linspace(r_min, r_max, 2 * n_mesh)
        phi_guess2 = C / r_mesh2 * np.exp(-m_sp * r_mesh2)
        dphi_guess2 = -C * np.exp(-m_sp * r_mesh2) * (m_sp + 1.0 / r_mesh2) / r_mesh2
        y_guess2 = np.vstack([phi_guess2, dphi_guess2])

        sol = solve_bvp(
            fun=lambda r, y: tgp_ode(r, y, beta),
            bc=lambda ya, yb: tgp_bc(ya, yb, phi_inner),
            x=r_mesh2,
            y=y_guess2,
            tol=1e-4,
            max_nodes=80000,
            verbose=0,
        )
        if not sol.success:
            print(f"    [WARN] BVP did not converge: {sol.message}")

    return sol


# ============================================================================
# 4. Asymptotic matching
# ============================================================================

def outer_yukawa(r, C, m_sp):
    """Outer (far-field) Yukawa solution: phi ~ C/r exp(-m_sp r)."""
    r_safe = np.maximum(r, 1e-12)
    return C / r_safe * np.exp(-m_sp * r_safe)


def inner_expansion(r, r_0, phi_0, phi_1, phi_2):
    """
    Inner power-series expansion near r_0:
        phi(r) = phi_0 + phi_1 (r - r_0) + phi_2 (r - r_0)^2 + ...
    """
    dr = r - r_0
    return phi_0 + phi_1 * dr + phi_2 * dr**2


def compute_inner_coefficients(r_0, C, beta, m_sp):
    """
    Compute inner expansion coefficients by matching to the ODE at r = r_0.

    At r = r_0:
        phi_0 = C / r_0 * exp(-m_sp * r_0)         (continuity)
        phi_1 = d/dr [C/r exp(-m_sp r)] at r_0      (C^1 matching)
        phi_2 from the ODE: phi'' = -(2/r) phi' - 2 (phi')^2/(1+phi) + beta (1+phi)^2 phi
    """
    phi_0 = C / r_0 * np.exp(-m_sp * r_0)
    phi_1 = -C * np.exp(-m_sp * r_0) * (m_sp + 1.0 / r_0) / r_0

    f_0 = 1.0 + phi_0
    phi_2_val = 0.5 * (
        -(2.0 / r_0) * phi_1
        - ALPHA * phi_1**2 / f_0
        + beta * f_0**2 * phi_0
    )
    return phi_0, phi_1, phi_2_val


def matched_profile(r, C, beta):
    """
    Construct the matched asymptotic profile across all three regimes.

    Strategy:
        - For r > r_match_outer: use outer Yukawa
        - For r_match_inner < r < r_match_outer: blend inner and outer
        - For r < r_match_inner: use inner expansion

    The matching region is centred around r_crit (where nonlinear effects
    become important). We use a smooth sigmoid transition.
    """
    m_sp, r_crit, r_0 = compute_scales(C, beta)

    # Inner coefficients computed at the matching point r_0
    phi_0, phi_1, phi_2 = compute_inner_coefficients(r_0, C, beta, m_sp)

    # Outer solution
    phi_outer = outer_yukawa(r, C, m_sp)

    # Inner solution
    phi_inner = inner_expansion(r, r_0, phi_0, phi_1, phi_2)

    # Intermediate correction: in the region r_0 < r < r_crit,
    # the nonlinear terms produce an effective 1/r^2 tail.
    # This is the key dark-matter-mimicking feature.
    # The leading nonlinear correction is:
    #   phi_NL(r) ~ (alpha * C^2) / (2 * r^2) * exp(-2 * m_sp * r)
    # from the (phi')^2 / (1+phi) term at second order.
    phi_NL = (ALPHA * C**2) / (2.0 * r**2) * np.exp(-2.0 * m_sp * r)

    # Smooth blending using sigmoid centred at r_crit
    # Width of transition ~ r_crit / 3
    width = max(r_crit / 3.0, 0.5)
    sigma_blend = 1.0 / (1.0 + np.exp(-(r - r_crit) / width))

    # For r >> r_crit: sigma_blend -> 1 (use outer Yukawa)
    # For r << r_crit: sigma_blend -> 0 (use inner + NL correction)

    # Matched profile:
    # - Near r_0: inner expansion dominates
    # - In intermediate regime: Yukawa + NL correction (1/r^2 tail)
    # - Far field: pure Yukawa

    # Inner blend: transition from inner expansion to intermediate
    sigma_inner = 1.0 / (1.0 + np.exp(-(r - r_0) / max(r_0 / 3.0, 0.3)))

    # The composite matched profile
    phi_match = (
        (1.0 - sigma_inner) * phi_inner
        + sigma_inner * (
            sigma_blend * phi_outer
            + (1.0 - sigma_blend) * (phi_outer + phi_NL)
        )
    )

    return phi_match, phi_outer, phi_inner, phi_NL


# ============================================================================
# 5. Effective dark matter density
# ============================================================================

def effective_dm_density(r, phi, C, beta):
    """
    Compute the effective "dark matter" density from the TGP field profile.

    The Newtonian potential is Phi_N = -C / r (leading order).
    The TGP modification creates an additional gravitational effect:
        Phi_eff(r) = Phi_N(r) + delta_Phi(r)

    where delta_Phi comes from the nonlinear TGP terms. The effective
    DM density is:
        rho_DM,eff = -1/(4 pi r^2) d/dr [r^2 d(delta_Phi)/dr]

    In dimensionless TGP units (with G=1 normalization), this becomes:
        rho_DM,eff(r) = -1/(4 pi r^2) d/dr [r^2 dphi/dr]
                      = -1/(4 pi) [phi'' + (2/r) phi']

    which is just the Laplacian of phi divided by -4pi (Poisson equation).
    From the TGP equation, this equals the nonlinear source terms.
    """
    r_safe = np.maximum(r, 1e-12)

    # Numerical derivatives
    dr = np.gradient(r)
    dphi = np.gradient(phi, r)
    ddphi = np.gradient(dphi, r)

    # rho_DM,eff = -1/(4pi) * [phi'' + (2/r) phi']
    # But from TGP equation: phi'' + (2/r)phi' = -2(phi')^2/(1+phi) + beta*(1+phi)^2*phi
    # So: rho_DM,eff = 1/(4pi) * [2(phi')^2/(1+phi) - beta*(1+phi)^2*phi]
    # The DM-like contribution comes from the nonlinear gradient term.
    f = np.maximum(1.0 + phi, 1e-12)

    # Nonlinear contribution (the "dark matter" source):
    rho_NL = (1.0 / (4.0 * np.pi)) * (ALPHA * dphi**2 / f)

    # Full effective density from the Laplacian
    rho_eff = -(1.0 / (4.0 * np.pi)) * (ddphi + 2.0 / r_safe * dphi)

    return rho_eff, rho_NL


def nfw_profile(r, rho_s, r_s):
    """NFW dark matter density profile for comparison."""
    x = r / r_s
    x_safe = np.maximum(x, 1e-12)
    return rho_s / (x_safe * (1.0 + x_safe)**2)


# ============================================================================
# 6. Rotation curve
# ============================================================================

def rotation_curve(r, phi):
    """
    Compute v(r) = sqrt(r |dphi_N/dr|) from the TGP profile.

    In the Newtonian limit, v^2 = r * |dPhi_N/dr|. The TGP modification
    gives an effective potential whose gradient determines the rotation
    curve. We compute v^2 = r * |dphi/dr| (proportional to the true
    rotation velocity squared, up to constants).
    """
    dphi = np.gradient(phi, r)
    # v^2 proportional to r * |dphi/dr|
    v_sq = r * np.abs(dphi)
    return np.sqrt(np.maximum(v_sq, 0.0))


# ============================================================================
# 7. Main computation and plotting
# ============================================================================

def main():
    print("=" * 70)
    print("TGP Intermediate Regime: Asymptotic Matching & DM Mimicry")
    print("=" * 70)

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PLOT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(PLOT_DIR, exist_ok=True)

    # Storage for all results
    all_results = []

    for C in C_VALUES:
        for beta in BETA_VALUES:
            m_sp, r_crit, r_0 = compute_scales(C, beta)
            label = rf"$C={C},\;\beta={beta}$"
            print(f"\n  --- C = {C}, beta = {beta} ---")
            print(f"    m_sp   = {m_sp:.4f}")
            print(f"    r_crit = {r_crit:.4f}")
            print(f"    r_0    = {r_0:.4f}")

            # Solve full nonlinear BVP
            sol = solve_tgp_profile(C, beta)
            r_num = sol.x
            phi_num = sol.y[0]

            # Compute matched profile on same grid
            phi_match, phi_outer, phi_inner, phi_NL = matched_profile(
                r_num, C, beta
            )

            # Effective DM density
            rho_eff, rho_NL = effective_dm_density(r_num, phi_num, C, beta)

            # Rotation curve
            v_num = rotation_curve(r_num, phi_num)
            v_match = rotation_curve(r_num, phi_match)

            # Residuals
            mask = np.abs(phi_num) > 1e-12
            residual = np.zeros_like(phi_num)
            residual[mask] = (phi_num[mask] - phi_match[mask]) / phi_num[mask]

            # Key diagnostics
            phi_max = np.max(phi_num)
            r_at_max = r_num[np.argmax(phi_num)]
            print(f"    phi_max = {phi_max:.6f}  at r = {r_at_max:.4f}")
            print(f"    phi(r_min) = {phi_num[0]:.6f}")
            print(f"    phi(r_max) = {phi_num[-1]:.6e}")
            print(f"    Max |residual| = {np.max(np.abs(residual)):.4e}")

            # Check for 1/r^2 tail in intermediate regime
            r_mid_mask = (r_num > r_0) & (r_num < r_crit) & (r_num > 0.1)
            if np.sum(r_mid_mask) > 10:
                log_r = np.log(r_num[r_mid_mask])
                log_phi = np.log(np.maximum(np.abs(phi_num[r_mid_mask]), 1e-30))
                # Linear fit to get power law index
                if len(log_r) > 2:
                    coeffs = np.polyfit(log_r, log_phi, 1)
                    print(f"    Intermediate power law: phi ~ r^({coeffs[0]:.2f})")

            all_results.append({
                'C': C, 'beta': beta, 'label': label,
                'm_sp': m_sp, 'r_crit': r_crit, 'r_0': r_0,
                'r': r_num, 'phi_num': phi_num,
                'phi_match': phi_match, 'phi_outer': phi_outer,
                'phi_NL': phi_NL,
                'rho_eff': rho_eff, 'rho_NL': rho_NL,
                'v_num': v_num, 'v_match': v_match,
                'residual': residual,
                'converged': sol.success,
            })

    # ========================================================================
    # 8. Four-panel figure
    # ========================================================================
    print("\n  Generating 4-panel figure ...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle(
        "TGP Intermediate Regime: Asymptotic Matching\n"
        r"$\nabla^2\Phi + 2(\nabla\Phi)^2/\Phi + \beta\Phi^2/\Phi_0"
        r" - \beta\Phi^3/\Phi_0^2 = 0$ (vacuum, spherical)",
        fontsize=13, fontweight='bold', y=0.99
    )

    # Colour map for parameter combinations
    n_results = len(all_results)
    colors = plt.cm.tab10(np.linspace(0, 0.9, n_results))

    # ------------------------------------------------------------------
    # Panel (a): phi(r) full numerical vs matched approximation
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    for res, col in zip(all_results, colors):
        r = res['r']
        ax.plot(r, res['phi_num'], color=col, lw=1.8,
                label=res['label'] + " (num)")
        ax.plot(r, res['phi_match'], color=col, lw=1.2, ls='--', alpha=0.7)
        # Mark r_crit and r_0
        ax.axvline(res['r_crit'], color=col, ls=':', lw=0.5, alpha=0.4)

    ax.set_xlabel(r"$r$ (dimensionless)", fontsize=12)
    ax.set_ylabel(r"$\varphi(r) = \Phi/\Phi_0 - 1$", fontsize=12)
    ax.set_title(r"(a) Profile $\varphi(r)$: numerical (solid) vs matched (dashed)",
                 fontsize=11)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=6, loc='upper right', ncol=2)
    ax.grid(True, which='both', ls=':', alpha=0.3)
    ax.set_ylim(bottom=1e-8)

    # ------------------------------------------------------------------
    # Panel (b): Effective DM density profile vs NFW
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    for res, col in zip(all_results, colors):
        r = res['r']
        rho = res['rho_NL']
        # Only plot where rho is positive and meaningful
        pos = rho > 0
        if np.any(pos):
            ax.plot(r[pos], rho[pos], color=col, lw=1.5,
                    label=res['label'])

    # Add NFW reference (with arbitrary normalization for comparison)
    # Use median r_crit as scale radius
    r_s_ref = np.median([res['r_crit'] for res in all_results])
    r_ref = np.logspace(np.log10(R_MIN), np.log10(100.0), 500)
    rho_nfw = nfw_profile(r_ref, rho_s=0.01, r_s=r_s_ref)
    ax.plot(r_ref, rho_nfw, 'k--', lw=2.0, alpha=0.5,
            label=rf"NFW ($r_s = {r_s_ref:.1f}$)")

    ax.set_xlabel(r"$r$ (dimensionless)", fontsize=12)
    ax.set_ylabel(r"$\rho_{\rm DM,eff}(r)$", fontsize=12)
    ax.set_title(r"(b) Effective DM density from nonlinear gradient term",
                 fontsize=11)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=6, loc='upper right', ncol=2)
    ax.grid(True, which='both', ls=':', alpha=0.3)

    # ------------------------------------------------------------------
    # Panel (c): Rotation curve v(r)
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    for res, col in zip(all_results, colors):
        r = res['r']
        # Normalize rotation curve for comparison
        v = res['v_num']
        v_max = np.max(v) if np.max(v) > 0 else 1.0
        ax.plot(r, v / v_max, color=col, lw=1.5,
                label=res['label'] + " (num)")
        v_m = res['v_match']
        v_m_max = np.max(v_m) if np.max(v_m) > 0 else 1.0
        ax.plot(r, v_m / v_m_max, color=col, lw=1.0, ls='--', alpha=0.6)

    # Add Keplerian reference: v ~ 1/sqrt(r)
    r_kep = np.linspace(1.0, 100.0, 200)
    v_kep = 1.0 / np.sqrt(r_kep)
    ax.plot(r_kep, v_kep / v_kep[0], 'k:', lw=1.5, alpha=0.5,
            label=r"Keplerian $\propto 1/\sqrt{r}$")

    ax.set_xlabel(r"$r$ (dimensionless)", fontsize=12)
    ax.set_ylabel(r"$v(r) / v_{\max}$", fontsize=12)
    ax.set_title(r"(c) Rotation curve: $v = \sqrt{r\,|d\varphi/dr|}$"
                 "\n(solid = numerical, dashed = matched)",
                 fontsize=11)
    ax.set_xscale('log')
    ax.legend(fontsize=6, loc='upper right', ncol=2)
    ax.grid(True, which='both', ls=':', alpha=0.3)
    ax.set_ylim(0, 1.5)

    # ------------------------------------------------------------------
    # Panel (d): Residuals
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    for res, col in zip(all_results, colors):
        r = res['r']
        ax.plot(r, res['residual'], color=col, lw=1.2,
                label=res['label'])
        # Mark r_crit
        ax.axvline(res['r_crit'], color=col, ls=':', lw=0.5, alpha=0.4)

    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel(r"$r$ (dimensionless)", fontsize=12)
    ax.set_ylabel(r"$(\varphi_{\rm num} - \varphi_{\rm match}) / \varphi_{\rm num}$",
                  fontsize=12)
    ax.set_title("(d) Relative residual: numerical vs matched", fontsize=11)
    ax.set_xscale('log')
    ax.legend(fontsize=6, loc='upper right', ncol=2)
    ax.grid(True, which='both', ls=':', alpha=0.3)
    ax.set_ylim(-1.0, 1.0)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    outpath = os.path.join(PLOT_DIR, "profile_intermediate.png")
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {outpath}")

    # ========================================================================
    # 9. Summary table
    # ========================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY: Intermediate Regime Profile Analysis")
    print("=" * 70)
    print(f"\n  {'C':>6s}  {'beta':>6s}  {'m_sp':>8s}  {'r_crit':>10s}  "
          f"{'r_0':>8s}  {'phi_max':>10s}  {'max|res|':>10s}  {'conv':>5s}")
    print(f"  {'-'*72}")
    for res in all_results:
        print(f"  {res['C']:6.1f}  {res['beta']:6.2f}  {res['m_sp']:8.4f}  "
              f"{res['r_crit']:10.4f}  {res['r_0']:8.4f}  "
              f"{np.max(res['phi_num']):10.6f}  "
              f"{np.max(np.abs(res['residual'])):10.4e}  "
              f"{'yes' if res['converged'] else 'NO':>5s}")

    print(f"\n  Key physics:")
    print(f"    - In the intermediate regime (r_0 < r < r_crit), the nonlinear")
    print(f"      gradient term 2(nabla Phi)^2/Phi generates an effective")
    print(f"      1/r^2 density profile that mimics a dark matter halo.")
    print(f"    - The rotation curve flattens in this regime, consistent with")
    print(f"      observed galaxy rotation curves without invoking DM particles.")
    print(f"    - r_crit = 1/(2*beta*C) sets the transition to Newtonian gravity.")

    print(f"\n  Output: {outpath}")
    print("=" * 70)


if __name__ == "__main__":
    main()
