#!/usr/bin/env python3
"""
ex200 -- Lyapunov ratio TGP/Newton vs beta WITH 3-body forces (yukawa_feynman)
================================================================================

Like ex199 but uses the yukawa_feynman backend (V2 + V3) instead of
pairwise-only (V2). The 3-body irreducible force from Phi^4 nonlinearity
creates an additional repulsive barrier that should lower beta_crit.

Also tests convergence with increasing t_final at key beta values.
"""

import sys, os, time
import numpy as np

repo = r"C:\Users\Mateusz\Documents\ObsydnianMain\TGP\TGP_v1"
sys.path.insert(0, repo)

from nbody.lyapunov import (
    largest_lyapunov_exponent_benettin_leapfrog,
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_tgp_pairwise_softened,
    acceleration_jacobian_yukawa_feynman_analytic,
    random_velocities_for_excess_energy,
    scale_velocities_match_energy,
    total_mechanical_energy,
    pythagorean_three_body_burrau,
)
from nbody.dynamics_v2 import (
    forces_newton,
    potential_newton,
    potential_tgp,
    analytical_forces_tgp_pairwise,
)
from nbody.dynamics_backends import build_tgp_integration_pair

quick = "--quick" in sys.argv
SOFT = 1e-6
G = 4.0 * np.pi
JAC_EPS = 1e-5

if quick:
    N_QUAD = 10
    T_FINALS_SCAN = [1.5]
    T_FINALS_CONV = [1.0, 1.5, 2.0]
    DT = 0.05
    RENORM = 7
    BETAS = [0.05, 0.08, 0.20, 0.50, 1.0]
    CONV_BETAS = [0.08, 0.50]
else:
    N_QUAD = 14
    T_FINALS_SCAN = [2.0]
    T_FINALS_CONV = [1.0, 2.0, 3.0, 4.0, 6.0]
    DT = 0.04
    RENORM = 10
    BETAS = [0.02, 0.05, 0.08, 0.12, 0.20, 0.35, 0.50, 0.80, 1.0, 1.5]
    CONV_BETAS = [0.08, 0.35, 0.50, 1.0]


def run_one_beta(beta, t_final, backend="yukawa_feynman"):
    """Run Newton vs TGP at given beta and t_final."""
    gamma = beta
    pos0, vel0, C = pythagorean_three_body_burrau()

    # Newton
    def pot_N(p, c):
        return potential_newton(p, c, G=G, softening=SOFT)
    def acc_N(p, c):
        return forces_newton(p, c, G=G, softening=SOFT) / c[:, None]
    def jac_N(p, c):
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=SOFT)

    H_N = total_mechanical_energy(pos0, vel0, C, pot_N)

    # TGP
    acc_T, pot_T = build_tgp_integration_pair(
        backend,
        beta=beta, gamma=gamma, softening=SOFT,
        include_3body=(backend != "pairwise"),
        n_quad_feynman=N_QUAD,
    )

    # Match energy
    vel_T = scale_velocities_match_energy(pos0, vel0, C, pot_T, target_energy=H_N)
    if float(np.sum(vel_T**2)) < 1e-28:
        vel_T = random_velocities_for_excess_energy(
            pos0, C, pot_T, target_energy=H_N,
            rng=np.random.default_rng(200),
        )
    H_T = total_mechanical_energy(pos0, vel_T, C, pot_T)

    # Jacobians — wrap to match (pos, C) -> ndarray signature
    if backend == "pairwise":
        def jac_T(p, c):
            return acceleration_jacobian_tgp_pairwise_softened(
                p, c, beta=beta, gamma=gamma, softening=SOFT)
    else:
        # yukawa_feynman: analytic Jacobian with keyword args wrapped
        def jac_T(p, c):
            return acceleration_jacobian_yukawa_feynman_analytic(
                p, c, beta=beta, gamma=gamma, softening=SOFT,
                n_quad_feynman=N_QUAD)

    # Newton lambda_max
    try:
        lam_N, _ = largest_lyapunov_exponent_benettin_leapfrog(
            pos0, vel0, C, acc_N, t_final=t_final, dt=DT,
            renorm_every=RENORM, jac_eps=JAC_EPS,
            position_jacobian_fn=jac_N,
            rng=np.random.default_rng(200))
    except Exception as e:
        lam_N = np.nan

    # TGP lambda_max
    try:
        lam_T, _ = largest_lyapunov_exponent_benettin_leapfrog(
            pos0, vel_T, C, acc_T, t_final=t_final, dt=DT,
            renorm_every=RENORM, jac_eps=JAC_EPS,
            position_jacobian_fn=jac_T,
            rng=np.random.default_rng(200))
    except Exception as e:
        lam_T = np.nan

    ratio = lam_T / lam_N if (np.isfinite(lam_N) and lam_N > 0.001) else np.nan

    return {
        'beta': beta, 't_final': t_final, 'backend': backend,
        'lam_N': lam_N, 'lam_T': lam_T, 'ratio': ratio,
        'H_N': H_N, 'H_T': H_T,
        'energy_ok': abs(H_T - H_N) < 1e-3 * max(1.0, abs(H_N)),
    }


def main():
    print("=" * 78)
    print("ex200 -- Lyapunov beta scan WITH 3-body forces (yukawa_feynman)")
    print("=" * 78)
    print(f"  G_Newton = 4*pi, n_quad = {N_QUAD}, dt = {DT}")

    # ── Part 1: Beta scan with V2+V3 ──
    print(f"\n--- Part 1: Beta scan (yukawa_feynman, V2+V3) t={T_FINALS_SCAN[0]} ---")
    print(f"  {'beta':>6s} {'lam_N':>8s} {'lam_T':>8s} {'ratio':>7s} {'H_ok':>5s} {'verdict':>14s} {'time':>6s}")
    print(f"  {'-'*58}")

    scan_results = []
    for beta in BETAS:
        t0 = time.time()
        r = run_one_beta(beta, T_FINALS_SCAN[0], "yukawa_feynman")
        elapsed = time.time() - t0

        verdict = "N/A"
        if np.isfinite(r['ratio']):
            verdict = "SUPPRESS" if r['ratio'] < 0.9 else ("ENHANCE" if r['ratio'] > 1.1 else "COMPARABLE")
        h_ok = "OK" if r['energy_ok'] else "FAIL"

        print(f"  {r['beta']:6.3f} {r['lam_N']:8.4f} {r['lam_T']:8.4f} "
              f"{r['ratio']:7.3f} {h_ok:>5s} {verdict:>14s} {elapsed:5.1f}s")
        scan_results.append(r)

    # ── Part 2: Same scan with V2 only (for comparison) ──
    print(f"\n--- Part 2: Beta scan (pairwise V2 only) t={T_FINALS_SCAN[0]} ---")
    print(f"  {'beta':>6s} {'lam_N':>8s} {'lam_T':>8s} {'ratio':>7s} {'verdict':>14s}")
    print(f"  {'-'*48}")

    pair_results = []
    for beta in BETAS:
        t0 = time.time()
        r = run_one_beta(beta, T_FINALS_SCAN[0], "pairwise")
        elapsed = time.time() - t0

        verdict = "N/A"
        if np.isfinite(r['ratio']):
            verdict = "SUPPRESS" if r['ratio'] < 0.9 else ("ENHANCE" if r['ratio'] > 1.1 else "COMPARABLE")

        print(f"  {r['beta']:6.3f} {r['lam_N']:8.4f} {r['lam_T']:8.4f} "
              f"{r['ratio']:7.3f} {verdict:>14s}")
        pair_results.append(r)

    # ── Part 3: Convergence with t_final ──
    print(f"\n--- Part 3: Convergence (ratio vs t_final) at key betas ---")
    print(f"  {'beta':>6s} {'t_final':>7s} {'lam_N':>8s} {'lam_T(V3)':>10s} {'ratio(V3)':>10s} "
          f"{'lam_T(V2)':>10s} {'ratio(V2)':>10s}")
    print(f"  {'-'*70}")

    conv_results = []
    for beta in CONV_BETAS:
        for t in T_FINALS_CONV:
            r_v3 = run_one_beta(beta, t, "yukawa_feynman")
            r_v2 = run_one_beta(beta, t, "pairwise")

            print(f"  {beta:6.3f} {t:7.1f} {r_v3['lam_N']:8.4f} "
                  f"{r_v3['lam_T']:10.4f} {r_v3['ratio']:10.3f} "
                  f"{r_v2['lam_T']:10.4f} {r_v2['ratio']:10.3f}")
            conv_results.append({'beta': beta, 't': t,
                                'ratio_v3': r_v3['ratio'], 'ratio_v2': r_v2['ratio']})

    # ── Summary ──
    print(f"\n{'='*78}")
    print("SUMMARY: V2+V3 vs V2-only beta scan")
    print(f"{'='*78}")
    print(f"\n  {'beta':>6s} {'ratio(V3)':>10s} {'ratio(V2)':>10s} {'V3 effect':>12s}")
    print(f"  {'-'*42}")

    for s, p in zip(scan_results, pair_results):
        if np.isfinite(s['ratio']) and np.isfinite(p['ratio']):
            effect = "V3 helps" if s['ratio'] < p['ratio'] else "V3 hurts"
            print(f"  {s['beta']:6.3f} {s['ratio']:10.3f} {p['ratio']:10.3f} {effect:>12s}")

    # Count
    v3_suppress = sum(1 for r in scan_results if np.isfinite(r['ratio']) and r['ratio'] < 0.9)
    v2_suppress = sum(1 for r in pair_results if np.isfinite(r['ratio']) and r['ratio'] < 0.9)
    n_total = len(BETAS)

    print(f"\n  V2+V3 suppresses: {v3_suppress}/{n_total}")
    print(f"  V2-only suppresses: {v2_suppress}/{n_total}")

    if v3_suppress > v2_suppress:
        print(f"  -> V3 (3-body force) HELPS suppress chaos (as predicted by P1)")
    elif v3_suppress == v2_suppress:
        print(f"  -> V3 has NEGLIGIBLE impact on suppression count")
    else:
        print(f"  -> V3 REDUCES suppression (unexpected)")

    # Find beta_crit for each
    for label, results in [("V2+V3", scan_results), ("V2-only", pair_results)]:
        ratios_valid = [(r['beta'], r['ratio']) for r in results if np.isfinite(r['ratio'])]
        transition = None
        for i in range(len(ratios_valid)-1):
            b1, r1 = ratios_valid[i]
            b2, r2 = ratios_valid[i+1]
            if r1 > 1.0 and r2 < 1.0:
                # Linear interpolation for crossing
                beta_cross = b1 + (1.0 - r1) * (b2 - b1) / (r2 - r1)
                transition = beta_cross
                break
        if transition:
            print(f"  {label}: beta_crit ~ {transition:.3f}")


if __name__ == "__main__":
    main()
