"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
two_body_three_regimes.py  --  Theory of Generated Space (TGP)
===============================================================
Numerical verification that three regimes exist for beta = gamma > 9C/2.

This script validates Proposition prop:trzy-rezimy-beta-gamma:
for the vacuum-consistent case beta = gamma, the effective potential
V_eff(d) has three regimes (attraction -> repulsion -> confinement)
if and only if beta > 9C/2, where C is the dimensionless source strength.

Method:
  Analytical V_eff(d) from linearized overlap integrals:
    V_eff(d) = -4piC^2/d + 8pibetaC^2/d^2 - 24pibetaC^3/d^3

  Force F(d) = -dV/dd zeros give regime transitions.

Outputs (saved to tooling/scripts/plots/):
    three_regimes_beta_eq_gamma.png  -- V_eff and F for several C/beta ratios
    regime_phase_diagram.png         -- phase diagram in (C, beta) space
    critical_mass_scaling.png        -- M_crit vs beta

Author: TGP consistency analysis
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -- Output directory ------------------------------------------------------
PLOT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(PLOT_DIR, exist_ok=True)

# -- Analytical functions --------------------------------------------------

def V_eff(d, C, beta):
    """
    Effective potential for beta = gamma (vacuum condition).
    V_eff(d) = -4piC^2/d + 8pibetaC^2/d^2 - 24pibetaC^3/d^3
    """
    return -4*np.pi*C**2/d + 8*np.pi*beta*C**2/d**2 - 24*np.pi*beta*C**3/d**3


def force(d, C, beta):
    """
    Force F = -dV/dd.
    F(d) = -4piC^2/d^2 + 16pibetaC^2/d^3 - 72pibetaC^3/d^4
    """
    return -4*np.pi*C**2/d**2 + 16*np.pi*beta*C**2/d**3 - 72*np.pi*beta*C**3/d**4


def force_zeros(C, beta):
    """
    Zeros of the force: d^2 - 4betad + 18betaC = 0.
    Returns (d_rep, d_well) or None if no real roots.
    """
    discriminant = 16*beta**2 - 72*beta*C
    if discriminant < 0:
        return None
    sqrt_disc = np.sqrt(discriminant)
    d_rep = 2*beta - sqrt_disc/2    # smaller root (repulsion onset)
    d_well = 2*beta + sqrt_disc/2   # larger root (well onset)
    # Correcting: from d^2 - 4betad + 18betaC = 0
    # d = (4beta +- sqrt(16beta^2 - 72betaC)) / 2 = 2beta +- sqrt(4beta^2 - 18betaC)
    discriminant2 = 4*beta**2 - 18*beta*C
    if discriminant2 < 0:
        return None
    sqrt_disc2 = np.sqrt(discriminant2)
    d_rep = 2*beta - sqrt_disc2
    d_well = 2*beta + sqrt_disc2
    if d_rep <= 0:
        return None
    return d_rep, d_well


def critical_ratio():
    """The critical value: beta/C = 9/2."""
    return 9.0 / 2.0


# -- Plot 1: V_eff and force for several C/beta ratios -----------------------

def plot_three_regimes():
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        r"TGP: Three regimes for $\beta = \gamma$ (Prop. 3.X)",
        fontsize=14, fontweight='bold'
    )

    beta = 5.0
    # Cases: beta/C = 3 (subcritical), 4.5 (critical), 6, 10, 20 (supercritical)
    ratios = [3.0, 4.5, 6.0, 10.0, 20.0, 50.0]
    labels = [
        r"$\beta/C = 3$ (subcritical)",
        r"$\beta/C = 4.5$ (critical)",
        r"$\beta/C = 6$ (supercritical)",
        r"$\beta/C = 10$",
        r"$\beta/C = 20$",
        r"$\beta/C = 50$",
    ]

    for idx, (ratio, label) in enumerate(zip(ratios, labels)):
        ax = axes.flat[idx]
        C = beta / ratio
        d = np.linspace(0.3*C, max(8*beta, 30*C), 2000)

        V = V_eff(d, C, beta)
        F = force(d, C, beta)

        # Normalize for plotting
        V_scale = np.max(np.abs(V[np.isfinite(V)])) if np.any(np.isfinite(V)) else 1.0
        F_scale = np.max(np.abs(F[np.isfinite(F)])) if np.any(np.isfinite(F)) else 1.0

        ax.plot(d, V/V_scale, 'b-', linewidth=1.5, label=r'$V_{\rm eff}/|V|_{\max}$')
        ax.plot(d, F/F_scale, 'r--', linewidth=1.2, label=r'$F/|F|_{\max}$')
        ax.axhline(0, color='gray', linewidth=0.5)

        # Mark regime transitions
        zeros = force_zeros(C, beta)
        if zeros is not None:
            d_rep, d_well = zeros
            ax.axvline(d_rep, color='green', linestyle=':', alpha=0.7,
                       label=rf'$d_{{\rm rep}} = {d_rep:.2f}$')
            ax.axvline(d_well, color='orange', linestyle=':', alpha=0.7,
                       label=rf'$d_{{\rm well}} = {d_well:.2f}$')
            # Shade regimes
            ax.axvspan(0, d_rep, alpha=0.05, color='purple')   # III: confinement
            ax.axvspan(d_rep, d_well, alpha=0.05, color='red')  # II: repulsion
            ax.axvspan(d_well, d[-1], alpha=0.05, color='blue') # I: gravity

        is_super = ratio > critical_ratio()
        status = "OK 3 regimes" if is_super else ("critical" if abs(ratio - 4.5) < 0.01 else "FAIL 1 regime")
        ax.set_title(f"{label}\n{status}", fontsize=10)
        ax.set_xlabel(r'$d$')
        ax.legend(fontsize=7, loc='best')
        ax.set_ylim(-1.5, 1.5)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "three_regimes_beta_eq_gamma.png"), dpi=150)
    plt.close()
    print("Saved: three_regimes_beta_eq_gamma.png")


# -- Plot 2: Phase diagram in (C, beta) space --------------------------------

def plot_phase_diagram():
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    C_arr = np.linspace(0.01, 5.0, 500)
    beta_arr = np.linspace(0.01, 30.0, 500)
    CC, BB = np.meshgrid(C_arr, beta_arr)

    # Number of regimes: 3 if beta > 9C/2, else 1
    n_regimes = np.where(BB > 4.5 * CC, 3, 1)

    ax.contourf(CC, BB, n_regimes, levels=[0.5, 2.0, 3.5],
                colors=['#ffcccc', '#ccffcc'], alpha=0.7)
    ax.plot(C_arr, 4.5 * C_arr, 'k-', linewidth=2,
            label=r'$\beta = \frac{9}{2}C$ (critical)')

    # Mark physical regions
    ax.annotate("3 regimes\n(particles)", xy=(0.5, 15), fontsize=12,
                ha='center', color='darkgreen', fontweight='bold')
    ax.annotate("1 regime\n(gravity only)", xy=(3.0, 5), fontsize=12,
                ha='center', color='darkred', fontweight='bold')

    # Typical particle region
    ax.axvspan(0, 0.1, alpha=0.1, color='blue', label='Typical particle $C$')

    ax.set_xlabel(r'Dimensionless source strength $C = \alpha_{\rm eff} M / (4\pi r_0)$',
                  fontsize=11)
    ax.set_ylabel(r'Self-coupling $\beta = \gamma$', fontsize=11)
    ax.set_title(r'TGP Regime Phase Diagram ($\beta = \gamma$, vacuum condition)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "regime_phase_diagram.png"), dpi=150)
    plt.close()
    print("Saved: regime_phase_diagram.png")


# -- Plot 3: Transition scales vs C ---------------------------------------

def plot_transition_scales():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    beta = 5.0
    C_arr = np.linspace(0.01, 2*beta/9 - 0.01, 500)  # subcritical to just below critical

    d_reps = []
    d_wells = []
    for C in C_arr:
        zeros = force_zeros(C, beta)
        if zeros:
            d_reps.append(zeros[0])
            d_wells.append(zeros[1])
        else:
            d_reps.append(np.nan)
            d_wells.append(np.nan)

    d_reps = np.array(d_reps)
    d_wells = np.array(d_wells)

    ax1.plot(C_arr, d_reps, 'r-', linewidth=2, label=r'$d_{\rm rep}$ (repulsion onset)')
    ax1.plot(C_arr, d_wells, 'b-', linewidth=2, label=r'$d_{\rm well}$ (well onset)')
    ax1.fill_between(C_arr, d_reps, d_wells, alpha=0.1, color='red',
                     label='Regime II (repulsion)')
    C_crit = 2*beta/9
    ax1.axvline(C_crit, color='gray', linestyle='--',
                label=rf'$C_{{\rm crit}} = 2\beta/9 = {C_crit:.2f}$')
    ax1.set_xlabel(r'Source strength $C$', fontsize=11)
    ax1.set_ylabel(r'Transition scale $d$', fontsize=11)
    ax1.set_title(rf'Regime transitions ($\beta = \gamma = {beta}$)', fontsize=12)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Plot separation ratio
    ratio = d_wells / d_reps
    ax2.plot(C_arr/C_crit, ratio, 'k-', linewidth=2)
    ax2.set_xlabel(r'$C / C_{\rm crit}$', fontsize=11)
    ax2.set_ylabel(r'$d_{\rm well} / d_{\rm rep}$', fontsize=11)
    ax2.set_title('Scale separation between regimes', fontsize=12)
    ax2.axhline(1, color='gray', linestyle='--', alpha=0.5, label='Merger point')
    ax2.set_xlim(0, 1.05)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "critical_mass_scaling.png"), dpi=150)
    plt.close()
    print("Saved: critical_mass_scaling.png")


# -- Validation: check analytical condition numerically --------------------

def validate_condition():
    """
    Numerical validation: count force zero crossings for beta = gamma
    across a range of C/beta ratios and verify the critical ratio is 9/2.
    """
    beta = 10.0
    n_test = 200
    C_ratios = np.linspace(0.01, 1.0, n_test)  # C/beta from 0.01 to 1

    print("\n=== Numerical validation of beta > 9C/2 condition ===")
    print(f"beta = gamma = {beta}")
    print(f"Critical ratio C/beta = 2/9 ~ {2/9:.6f}")
    print("-" * 50)

    transition_found = None
    for C_over_beta in C_ratios:
        C = C_over_beta * beta
        d = np.linspace(0.01*C, 50*beta, 100000)
        F = force(d, C, beta)

        # Count sign changes
        sign_changes = np.sum(np.diff(np.sign(F)) != 0)
        n_regimes = 1 + sign_changes  # 1 regime + sign changes

        if n_regimes == 3 and transition_found is None:
            # Still 3 regimes
            pass
        elif n_regimes < 3 and transition_found is None:
            transition_found = C_over_beta
            # Check previous
            C_prev = (C_over_beta - 1.0/n_test) * beta
            zeros_prev = force_zeros(C_prev, beta)
            print(f"Transition from 3->1 regimes at C/beta ~ {C_over_beta:.6f}")
            print(f"Analytical prediction: C/beta = 2/9 ~ {2/9:.6f}")
            print(f"Relative error: {abs(C_over_beta - 2/9)/(2/9)*100:.2f}%")

    if transition_found is None:
        print("All tested ratios show 3 regimes (need larger C/beta range)")

    # Specific test cases
    print("\n--- Specific tests ---")
    test_cases = [
        (0.1, 5.0, True),   # beta/C = 50 >> 4.5: should have 3 regimes
        (1.0, 5.0, True),   # beta/C = 5 > 4.5: should have 3 regimes
        (1.2, 5.0, False),  # beta/C = 4.17 < 4.5: should have 1 regime
        (2.0, 5.0, False),  # beta/C = 2.5 < 4.5: should have 1 regime
        (0.01, 1.0, True),  # beta/C = 100 >> 4.5: fundamental particle
    ]

    all_pass = True
    for C, beta_test, expected_three in test_cases:
        zeros = force_zeros(C, beta_test)
        has_three = zeros is not None
        status = "OK" if has_three == expected_three else "FAIL FAIL"
        if has_three != expected_three:
            all_pass = False
        ratio = beta_test / C
        print(f"  C={C:.2f}, beta={beta_test:.1f}, beta/C={ratio:.1f}: "
              f"{'3 regimes' if has_three else '1 regime'} {status}")

    print(f"\nAll tests passed: {all_pass}")
    return all_pass


# -- Main ------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("TGP: Three-regime verification for beta = gamma")
    print("=" * 60)

    ok = validate_condition()

    plot_three_regimes()
    plot_phase_diagram()
    plot_transition_scales()

    print("\nDone. Plots saved to tooling/scripts/plots/")
    if not ok:
        print("WARNING: Some validation tests FAILED!")
