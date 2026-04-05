#!/usr/bin/env python3
"""
ex199 -- Lyapunov ratio TGP/Newton vs beta (parameter scan)
=============================================================

Scans beta = gamma from small (weak barrier) to large (strong barrier)
to find if there exists a regime where TGP suppresses chaos (P1 thesis).

Uses: Burrau IC, pairwise V2, G = 4*pi (matched coupling), matched H.
"""

import sys, os, time
import numpy as np

repo = r"C:\Users\Mateusz\Documents\ObsydnianMain\TGP\TGP_v1"
sys.path.insert(0, repo)

from nbody.lyapunov import (
    largest_lyapunov_exponent_benettin_leapfrog,
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_tgp_pairwise_softened,
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
from nbody.pairwise import force_zeros_2body

quick = "--quick" in sys.argv
SOFT = 1e-6
G = 4.0 * np.pi  # Match TGP gradient coupling
JAC_EPS = 1e-5
T_FINAL = 2.0 if not quick else 1.5
DT = 0.04 if not quick else 0.05
RENORM = 10 if not quick else 7

BETAS = [0.02, 0.05, 0.08, 0.12, 0.20, 0.35, 0.50, 0.80, 1.0, 1.5, 2.0]
if quick:
    BETAS = [0.05, 0.08, 0.20, 0.50, 1.0, 2.0]


def run_at_beta(beta):
    gamma = beta
    pos0, vel0, C = pythagorean_three_body_burrau()

    # Newton (constant G=4pi, independent of beta)
    def pot_N(p, c):
        return potential_newton(p, c, G=G, softening=SOFT)
    def acc_N(p, c):
        return forces_newton(p, c, G=G, softening=SOFT) / c[:, None]
    def jac_N(p, c):
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=SOFT)

    H_N = total_mechanical_energy(pos0, vel0, C, pot_N)

    # TGP pairwise at this beta
    def pot_T(p, c):
        return potential_tgp(p, c, beta, gamma, SOFT)
    def acc_T(p, c):
        return analytical_forces_tgp_pairwise(p, c, beta, gamma, SOFT) / c[:, None]
    def jac_T(p, c):
        return acceleration_jacobian_tgp_pairwise_softened(p, c, beta=beta, gamma=gamma, softening=SOFT)

    # Match energy
    vel_T = scale_velocities_match_energy(pos0, vel0, C, pot_T, target_energy=H_N)
    if float(np.sum(vel_T**2)) < 1e-28:
        vel_T = random_velocities_for_excess_energy(
            pos0, C, pot_T, target_energy=H_N, rng=np.random.default_rng(199)
        )
    H_T = total_mechanical_energy(pos0, vel_T, C, pot_T)

    # Repulsive barrier distance (for equal mass approximation using min C)
    C_min = min(C)
    zeros = force_zeros_2body(C_min, beta, gamma)
    d_rep = zeros[0] if zeros and len(zeros) >= 1 else np.nan

    # Lyapunov
    rng = np.random.default_rng(199)
    try:
        lam_N, _ = largest_lyapunov_exponent_benettin_leapfrog(
            pos0, vel0, C, acc_N, t_final=T_FINAL, dt=DT,
            renorm_every=RENORM, jac_eps=JAC_EPS,
            position_jacobian_fn=jac_N, rng=np.random.default_rng(199))
    except Exception:
        lam_N = np.nan

    try:
        lam_T, _ = largest_lyapunov_exponent_benettin_leapfrog(
            pos0, vel_T, C, acc_T, t_final=T_FINAL, dt=DT,
            renorm_every=RENORM, jac_eps=JAC_EPS,
            position_jacobian_fn=jac_T, rng=np.random.default_rng(199))
    except Exception:
        lam_T = np.nan

    ratio = lam_T / lam_N if (np.isfinite(lam_N) and lam_N > 0.001) else np.nan

    return {
        'beta': beta, 'lam_N': lam_N, 'lam_T': lam_T,
        'ratio': ratio, 'd_rep': d_rep,
        'H_N': H_N, 'H_T': H_T,
        'dH': abs(H_T - H_N),
    }


def main():
    print("=" * 78)
    print("ex199 -- Lyapunov ratio TGP/Newton vs beta")
    print("=" * 78)
    print(f"  G_Newton = 4*pi = {G:.4f} (matched to TGP gradient coupling)")
    print(f"  IC: Burrau (3,4,5), t_final = {T_FINAL}, dt = {DT}")
    print()

    print(f"  {'beta':>6s} {'lam_N':>8s} {'lam_TGP':>8s} {'ratio':>7s} "
          f"{'d_rep':>7s} {'dH':>10s} {'verdict':>14s}")
    print(f"  {'-'*66}")

    results = []
    for beta in BETAS:
        t0 = time.time()
        r = run_at_beta(beta)
        elapsed = time.time() - t0

        if np.isfinite(r['ratio']):
            if r['ratio'] < 0.9:
                verdict = "SUPPRESS"
            elif r['ratio'] > 1.1:
                verdict = "ENHANCE"
            else:
                verdict = "COMPARABLE"
        else:
            verdict = "N/A"

        print(f"  {r['beta']:6.3f} {r['lam_N']:8.4f} {r['lam_T']:8.4f} "
              f"{r['ratio']:7.3f} {r['d_rep']:7.4f} {r['dH']:10.2e} "
              f"{verdict:>14s}  [{elapsed:.1f}s]")
        results.append(r)

    # Summary
    print(f"\n{'='*78}")
    print("ANALYSIS")
    print(f"{'='*78}")

    suppress = [r for r in results if np.isfinite(r['ratio']) and r['ratio'] < 0.9]
    enhance = [r for r in results if np.isfinite(r['ratio']) and r['ratio'] > 1.1]
    comparable = [r for r in results if np.isfinite(r['ratio']) and 0.9 <= r['ratio'] <= 1.1]

    print(f"  Suppress: {len(suppress)} / {len(results)}")
    print(f"  Enhance:  {len(enhance)} / {len(results)}")
    print(f"  Comparable: {len(comparable)} / {len(results)}")

    if suppress:
        betas_s = [r['beta'] for r in suppress]
        print(f"  Suppression at beta = {betas_s}")
        print(f"  -> P1 thesis may hold for LARGE beta (strong barrier)")

    if enhance:
        betas_e = [r['beta'] for r in enhance]
        print(f"  Enhancement at beta = {betas_e}")

    # Trend
    valid = [(r['beta'], r['ratio']) for r in results if np.isfinite(r['ratio'])]
    if len(valid) >= 3:
        betas_v = [v[0] for v in valid]
        ratios_v = [v[1] for v in valid]
        print(f"\n  Trend: ratio at min beta = {ratios_v[0]:.3f}, max beta = {ratios_v[-1]:.3f}")
        if ratios_v[-1] < ratios_v[0]:
            print(f"  -> DECREASING trend: stronger barrier reduces chaos enhancement")
            if ratios_v[-1] < 1.0:
                print(f"  -> At high enough beta, TGP SUPPRESSES chaos!")
        else:
            print(f"  -> No clear decreasing trend")


if __name__ == "__main__":
    main()
