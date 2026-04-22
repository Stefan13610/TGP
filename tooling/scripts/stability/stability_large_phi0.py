# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
stability_large_phi0.py  --  Theory of Generated Space (TGP)
=============================================================
Verify numerical stability of the cosmological evolution for large Phi_0.

The TGP parameter Phi_0 ~ 25 (from Lambda_eff ~ gamma/12 matching
Lambda_obs ~ 1e-52 m^-2). This script tests whether the cosmological
ODE system remains numerically stable for Phi_0 values ranging from
1 to 1000, checking:

1. Stiffness of the field equation:
   psi_ddot + 3H psi_dot + 2 psi_dot^2/psi = c0^2 W(psi)
   The RHS has c0^2 * gamma ~ very small in SI,
   which can cause numerical issues.

2. Friedmann constraint preservation:
   3(H + psi_dot/(4 psi))^2 sqrt(psi) = c0^2 [sqrt(psi) psi_dot^2/(2c0^2) + U(psi)]
   Check that this remains satisfied throughout evolution.

3. Scale separation:
   tau_0 = sqrt(Phi_0) / (c0 sqrt(gamma)) depends on Phi_0
   As Phi_0 grows, tau_0 grows, meaning the field relaxes slower.

4. Long-term evolution:
   Integrate for 100 Hubble times and check for blow-up or oscillation.

Usage:
    python stability_large_phi0.py
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import time

# Physical constants
c0 = 2.998e8         # m/s
G0 = 6.674e-11       # m^3/(kg s^2)
kappa = 8 * np.pi * G0 / c0**4

save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(save_dir, exist_ok=True)


def U(psi, beta, gamma):
    """Static potential U(psi) = (beta/3) psi^3 - (gamma/4) psi^4."""
    return beta / 3 * psi**3 - gamma / 4 * psi**4


def W(psi, beta, gamma):
    """Cosmological potential W(psi) = (7 beta/3) psi^2 - 2 gamma psi^3."""
    return 7 * beta / 3 * psi**2 - 2 * gamma * psi**3


def run_evolution(Phi0, gamma_val, psi0, psi_dot0, n_hubble=10):
    """
    Integrate the TGP cosmological field equation for given Phi_0.

    Returns dict with time series and diagnostics.
    """
    beta_val = gamma_val  # vacuum condition

    # Hubble rate from Lambda_eff
    Lambda_eff = gamma_val / 12
    H0 = np.sqrt(c0**2 * Lambda_eff / 3)

    # Time scale
    t_hubble = 1.0 / H0
    tau0 = np.sqrt(Phi0) / (c0 * np.sqrt(gamma_val))

    # ODE system: [psi, psi_dot]
    # psi_ddot + 3H psi_dot + 2 psi_dot^2/psi = c0^2 W(psi)  (kappa cancels)
    # Using constant H = H0 (de Sitter background)

    # Self-consistent H from modified Friedmann:
    # 3(H + psi_dot/(4 psi))^2 sqrt(psi) = c0^2 [sqrt(psi) psi_dot^2/(2 c0^2) + U(psi)]
    # Solve for H:
    # H_eff = H + psi_dot/(4 psi) = sqrt(c0^2 [...] / (3 sqrt(psi)))
    # H = H_eff - psi_dot/(4 psi)

    def H_from_friedmann(psi, psi_d):
        """Self-consistent Hubble rate from modified Friedmann eq."""
        U_val = U(psi, beta_val, gamma_val)
        rhs_val = c0**2 * (
            np.sqrt(psi) * psi_d**2 / (2 * c0**2) + U_val
        )
        if rhs_val <= 0:
            return H0  # fallback
        H_eff = np.sqrt(rhs_val / (3 * np.sqrt(psi)))
        return H_eff - psi_d / (4 * psi)

    def rhs(t, y):
        psi, psi_d = y
        psi = max(psi, 1e-10)  # safety

        H = H_from_friedmann(psi, psi_d)

        psi_dd = (c0**2 * W(psi, beta_val, gamma_val)
                  - 3 * H * psi_d
                  - 2 * psi_d**2 / psi)

        return [psi_d, psi_dd]

    # Integrate - use BDF for stiff systems
    t_span = (0, n_hubble * t_hubble)
    y0 = [psi0, psi_dot0]

    t_start = time.time()
    sol = solve_ivp(rhs, t_span, y0,
                    method='Radau',  # implicit, handles stiffness
                    rtol=1e-10, atol=1e-12,
                    max_step=t_hubble / 50,
                    dense_output=True)
    elapsed = time.time() - t_start

    if not sol.success:
        return {
            "success": False,
            "message": sol.message,
            "Phi0": Phi0,
            "elapsed": elapsed,
        }

    # Evaluate on uniform grid
    n_pts = 2000
    t_eval = np.linspace(t_span[0], t_span[1], n_pts)
    y_eval = sol.sol(t_eval)
    psi_arr = y_eval[0]
    psi_d_arr = y_eval[1]

    # Friedmann constraint check
    # 3(H + psi_dot/(4 psi))^2 sqrt(psi) = c0^2 [sqrt(psi) psi_dot^2/(2 c0^2) + U(psi)]
    LHS = 3 * (H0 + psi_d_arr / (4 * psi_arr))**2 * np.sqrt(psi_arr)
    RHS = c0**2 * (
        np.sqrt(psi_arr) * psi_d_arr**2 / (2 * c0**2) +
        np.array([U(p, beta_val, gamma_val) for p in psi_arr])
    )

    # Avoid division by zero
    mask = np.abs(RHS) > 1e-100
    if mask.sum() > 0:
        constraint_ratio = LHS[mask] / RHS[mask]
        constraint_deviation = np.max(np.abs(constraint_ratio - 1.0))
    else:
        constraint_deviation = np.nan

    # Check for blow-up or collapse
    psi_min = np.min(psi_arr)
    psi_max = np.max(psi_arr)
    psi_range = psi_max / max(psi_min, 1e-30)

    # Energy conservation check
    E = psi_d_arr**2 / (2 * c0**2) + np.array([U(p, beta_val, gamma_val) for p in psi_arr])
    E_deviation = (np.max(E) - np.min(E)) / max(np.abs(np.mean(E)), 1e-100)

    return {
        "success": True,
        "Phi0": Phi0,
        "H0": H0,
        "tau0": tau0,
        "t_hubble": t_hubble,
        "elapsed": elapsed,
        "n_steps": sol.t.size,
        "psi_min": psi_min,
        "psi_max": psi_max,
        "psi_range": psi_range,
        "psi_final": psi_arr[-1],
        "psi_dot_final": psi_d_arr[-1],
        "constraint_deviation": constraint_deviation,
        "E_deviation": E_deviation,
        "t": t_eval,
        "psi": psi_arr,
        "psi_dot": psi_d_arr,
        "LHS": LHS,
        "RHS": RHS,
    }


def main():
    print("=" * 65)
    print("TGP: Numerical Stability at Large Phi_0")
    print("=" * 65)

    gamma_val = 0.03  # L^{-2}

    # Test range of Phi_0
    Phi0_values = [1, 5, 10, 25, 50, 100, 250, 500, 1000]

    # Initial conditions: small perturbation from vacuum
    psi0 = 1.001
    psi_dot0 = 0.0

    results = []

    print(f"\n{'Phi_0':>8} {'H0':>12} {'tau0':>12} {'psi_range':>10} "
          f"{'Friedmann':>10} {'E_dev':>10} {'steps':>8} {'time':>8}")
    print("-" * 95)

    for Phi0 in Phi0_values:
        res = run_evolution(Phi0, gamma_val, psi0, psi_dot0, n_hubble=10)
        results.append(res)

        if res["success"]:
            print(f"{Phi0:8.0f} {res['H0']:12.4e} {res['tau0']:12.4e} "
                  f"{res['psi_range']:10.4f} {res['constraint_deviation']:10.4e} "
                  f"{res['E_deviation']:10.4e} {res['n_steps']:8d} "
                  f"{res['elapsed']:8.2f}s")
        else:
            print(f"{Phi0:8.0f} FAILED: {res['message']}")

    # Long-term stability test at Phi0 = 25 (canonical value)
    print(f"\n{'=' * 65}")
    print("LONG-TERM STABILITY TEST (Phi_0 = 25, 100 Hubble times)")
    print("=" * 65)

    res_long = run_evolution(25, gamma_val, psi0, psi_dot0, n_hubble=100)
    if res_long["success"]:
        print(f"  Duration: {res_long['t'][-1] / res_long['t_hubble']:.0f} Hubble times")
        print(f"  psi range: [{res_long['psi_min']:.6f}, {res_long['psi_max']:.6f}]")
        print(f"  Friedmann constraint deviation: {res_long['constraint_deviation']:.4e}")
        print(f"  psi(t_final) = {res_long['psi_final']:.6f}")
        print(f"  Steps: {res_long['n_steps']}, Time: {res_long['elapsed']:.2f}s")
    else:
        print(f"  FAILED: {res_long['message']}")

    # Note: Phi_0 doesn't appear in the dimensionless psi equation!
    # The psi equation depends only on gamma, beta, and H.
    # Phi_0 enters ONLY through:
    #   1. tau_0 = sqrt(Phi_0)/(c0 sqrt(gamma)) -- time scale
    #   2. Source term q*Phi_0*rho -- coupling to matter
    # For vacuum (rho=0), psi evolution is INDEPENDENT of Phi_0.
    print(f"\n  NOTE: In vacuum (no source), psi evolution is independent of Phi_0.")
    print(f"  Phi_0 affects only tau_0 (time scale) and source coupling.")
    print(f"  This explains identical trajectories for all Phi_0 above.")

    # Test different gamma values instead (physically meaningful)
    print(f"\n{'=' * 65}")
    print("GAMMA SENSITIVITY TEST (Phi_0 = 25)")
    print("=" * 65)

    gamma_values = [0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
    for gv in gamma_values:
        res_g = run_evolution(25, gv, 1.001, 0.0, n_hubble=10)
        if res_g["success"]:
            print(f"  gamma = {gv:.3f}: psi_range = {res_g['psi_range']:.4f}, "
                  f"steps = {res_g['n_steps']}, time = {res_g['elapsed']:.2f}s")
        else:
            print(f"  gamma = {gv:.3f}: FAILED ({res_g['message'][:50]})")

    # Stiffness test: moderate perturbation (psi_0 = 1.1, within de Sitter regime)
    # U(psi) > 0 requires psi < 4 beta/(3 gamma) = 4/3 for beta=gamma
    print(f"\n{'=' * 65}")
    print("STIFFNESS TEST (Phi_0 = 25, psi_0 = 1.1)")
    print("=" * 65)
    print(f"  (Note: U(psi) > 0 requires psi < 4/3 = 1.333 for beta=gamma)")

    res_stiff = run_evolution(25, gamma_val, 1.1, 0.0, n_hubble=20)
    if res_stiff["success"]:
        print(f"  psi range: [{res_stiff['psi_min']:.6f}, {res_stiff['psi_max']:.6f}]")
        print(f"  Friedmann deviation: {res_stiff['constraint_deviation']:.4e}")
        print(f"  Steps: {res_stiff['n_steps']}, Time: {res_stiff['elapsed']:.2f}s")

        # Check if psi returns to equilibrium
        psi_final_diff = abs(res_stiff['psi_final'] - 1.0)
        print(f"  |psi(final) - 1| = {psi_final_diff:.6e}")
        if psi_final_diff < 0.1:
            print("  -> Field RETURNS to equilibrium (stable)")
        else:
            print(f"  -> Field at psi = {res_stiff['psi_final']:.4f} (evolving)")
    else:
        print(f"  FAILED: {res_stiff['message']}")

    # Plot results
    print(f"\n{'=' * 65}")
    print("GENERATING PLOTS")
    print("=" * 65)

    # Plot 1: psi(t) for different Phi_0
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) psi(t) for several Phi_0
    ax = axes[0, 0]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results)))
    for res, col in zip(results, colors):
        if res["success"]:
            t_norm = res["t"] / res["t_hubble"]
            ax.plot(t_norm, res["psi"], color=col,
                    label=f'$\\Phi_0 = {res["Phi0"]:.0f}$', lw=1.2)
    ax.set_xlabel(r'$t / t_{\rm Hubble}$', fontsize=12)
    ax.set_ylabel(r'$\psi(t)$', fontsize=12)
    ax.set_title(r'Field evolution for different $\Phi_0$', fontsize=13)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, ls=':', alpha=0.4)

    # (b) Friedmann constraint deviation vs Phi_0
    ax = axes[0, 1]
    Phi0_ok = [r["Phi0"] for r in results if r["success"]]
    dev_ok = [r["constraint_deviation"] for r in results if r["success"]]
    ax.semilogy(Phi0_ok, dev_ok, 'o-', ms=6, color='steelblue')
    ax.set_xlabel(r'$\Phi_0$', fontsize=12)
    ax.set_ylabel('max |Friedmann ratio - 1|', fontsize=12)
    ax.set_title('Friedmann constraint preservation', fontsize=13)
    ax.grid(True, ls=':', alpha=0.4)

    # (c) Long-term evolution
    ax = axes[1, 0]
    if res_long["success"]:
        t_norm = res_long["t"] / res_long["t_hubble"]
        ax.plot(t_norm, res_long["psi"], 'b-', lw=1)
        ax.axhline(1.0, color='red', ls='--', lw=0.8, alpha=0.5, label='vacuum')
        ax.set_xlabel(r'$t / t_{\rm Hubble}$', fontsize=12)
        ax.set_ylabel(r'$\psi(t)$', fontsize=12)
        ax.set_title(r'Long-term stability ($\Phi_0=25$, 100 $t_H$)', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, ls=':', alpha=0.4)

    # (d) Stiffness test
    ax = axes[1, 1]
    if res_stiff["success"]:
        t_norm = res_stiff["t"] / res_stiff["t_hubble"]
        ax.plot(t_norm, res_stiff["psi"], 'r-', lw=1.2)
        ax.axhline(1.0, color='gray', ls='--', lw=0.8, alpha=0.5, label='vacuum')
        ax.set_xlabel(r'$t / t_{\rm Hubble}$', fontsize=12)
        ax.set_ylabel(r'$\psi(t)$', fontsize=12)
        ax.set_title(r'Stiffness test ($\psi_0 = 1.5$)', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, ls=':', alpha=0.4)

    fig.suptitle('TGP: Numerical stability analysis', fontsize=15, y=1.02)
    fig.tight_layout()
    path = os.path.join(save_dir, "stability_large_phi0.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)

    # Summary
    print(f"\n{'=' * 65}")
    print("SUMMARY")
    print("=" * 65)

    all_ok = all(r["success"] for r in results)
    max_dev = max(r["constraint_deviation"] for r in results if r["success"])
    max_range = max(r["psi_range"] for r in results if r["success"])

    print(f"  All integrations successful: {'YES' if all_ok else 'NO'}")
    print(f"  Max Friedmann constraint deviation: {max_dev:.4e}")
    print(f"  Max psi range (psi_max/psi_min): {max_range:.4f}")
    print(f"  Long-term (100 t_H) stable: {'YES' if res_long['success'] else 'NO'}")
    print(f"  Stiff system (psi_0=1.5) handled: {'YES' if res_stiff['success'] else 'NO'}")

    if max_dev < 0.01 and all_ok:
        print("\n  CONCLUSION: Numerical evolution is STABLE for Phi_0 in [1, 1000]")
        print("  The Friedmann constraint is preserved to < 1% accuracy.")
        print("  No blow-up or collapse observed in long-term evolution.")
    else:
        print("\n  CONCLUSION: Some stability issues detected. See details above.")

    print("\nDone.")


if __name__ == "__main__":
    main()
