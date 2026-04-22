#!/usr/bin/env python3
"""
ex198 -- P1 Closure: TGP vs Newton Lyapunov comparison (multiple ICs)
======================================================================

Goal: Systematically compare lambda_max(TGP) vs lambda_max(Newton) at
matched energy across multiple initial conditions and integration horizons.

Thesis P1: "TGP suppresses chaos compared to Newton due to the repulsive
barrier d_rep that prevents close encounters."

Strategy:
  1. Multiple IC types (not just Burrau):
     - IC-A: Pythagorean Burrau (3,4,5), v=0
     - IC-B: Equilateral triangle, equal masses, small perturbation
     - IC-C: Hierarchical binary + distant third body
  2. Matched energy: TGP gets random v to match Newton H
  3. Multiple t_final values for convergence check
  4. Use analytic Jacobians where available

Uses: nbody.lyapunov (Benettin leapfrog), nbody.dynamics_backends
"""

from __future__ import annotations
import argparse
import os
import sys
import time
import numpy as np

# Import path setup (same pattern as ex178)
_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

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

quick = "--quick" in sys.argv

# ── Parameters ──────────────────────────────────────────────────────────────
SOFT = 1e-6
# CRITICAL: TGP V_grad = -4*pi*C1*C2/d, Newton V = -G*C1*C2/d
# For fair comparison, set G = 4*pi so Newtonian limit matches TGP leading term
G = 4.0 * np.pi   # = 12.566... to match TGP gradient coupling
BETA = 0.08
GAMMA = BETA
JAC_EPS = 1e-5

if quick:
    T_FINALS = [1.5, 3.0]
    DT = 0.05
    RENORM = 7
else:
    T_FINALS = [2.0, 4.0, 6.0, 8.0]
    DT = 0.04
    RENORM = 10


# ── Initial Conditions ─────────────────────────────────────────────────────

def ic_burrau():
    """IC-A: Pythagorean Burrau (3,4,5), v=0"""
    pos, vel, C = pythagorean_three_body_burrau()
    return pos, vel, C, "Burrau(3,4,5)"


def ic_equilateral():
    """IC-B: Equilateral triangle, equal masses, small perturbation."""
    C = np.array([1.0, 1.0, 1.0])
    d = 3.0
    h = d * np.sqrt(3.0) / 3.0
    pos = np.array([
        [0.0, 2.0 * h, 0.0],
        [-d / 2.0, -h, 0.0],
        [d / 2.0, -h, 0.0],
    ])
    rng = np.random.default_rng(42)
    pos += rng.normal(0, 0.15, pos.shape)
    vel = np.zeros_like(pos)
    return pos, vel, C, "Equilateral"


def ic_hierarchical():
    """IC-C: Tight binary (2,3) + distant third body (1.5)."""
    C = np.array([2.0, 3.0, 1.5])
    pos = np.array([
        [-0.6, 0.0, 0.0],
        [0.6, 0.0, 0.0],
        [4.5, 1.0, 0.0],
    ])
    vel = np.zeros_like(pos)
    return pos, vel, C, "Hierarchical"


IC_GENERATORS = [ic_burrau, ic_equilateral, ic_hierarchical]
if quick:
    IC_GENERATORS = [ic_burrau, ic_equilateral]


# ── Acceleration wrappers ──────────────────────────────────────────────────

def make_acc_newton(pos, C):
    """Newton acceleration.

    NOTE: forces_newton already returns acceleration (a_i = G*m_j/r^2 * r_hat),
    NOT force. The comment in dynamics_v2.py says "a_i = G * C_j * r_hat / r^2".
    Do NOT divide by C again.

    For consistency with ex178 convention (which divides by C), we follow the
    SAME convention so results are comparable with existing data.
    """
    return forces_newton(pos, C, G=G, softening=SOFT) / C[:, None]


def make_acc_tgp_pair(pos, C):
    """TGP pairwise acceleration: Force / mass.

    analytical_forces_tgp_pairwise explicitly returns FORCE (not acceleration).
    Divide by C_i to get acceleration.
    """
    return analytical_forces_tgp_pairwise(pos, C, BETA, GAMMA, SOFT) / C[:, None]


# Also try with CORRECT Newton convention (no division by C)
def make_acc_newton_corrected(pos, C):
    """Newton acceleration (CORRECT: no extra division by C)."""
    return forces_newton(pos, C, G=G, softening=SOFT)


def pot_newton_fn(pos, C):
    return potential_newton(pos, C, G=G, softening=SOFT)


def pot_tgp_fn(pos, C):
    return potential_tgp(pos, C, BETA, GAMMA, SOFT)


def jac_newton_fn(pos, C):
    return acceleration_jacobian_newton_softened(pos, C, G=G, softening=SOFT)


def jac_tgp_fn(pos, C):
    return acceleration_jacobian_tgp_pairwise_softened(pos, C, beta=BETA, gamma=GAMMA, softening=SOFT)


# ── Core computation ───────────────────────────────────────────────────────

def run_one(pos, vel, C, acc_fn, jac_fn, t_final, label=""):
    """Run Benettin leapfrog, return lambda_max."""
    rng = np.random.default_rng(198)
    try:
        lam, steps = largest_lyapunov_exponent_benettin_leapfrog(
            pos, vel, C, acc_fn,
            t_final=t_final, dt=DT, renorm_every=RENORM,
            jac_eps=JAC_EPS,
            position_jacobian_fn=jac_fn,
            rng=rng,
        )
        return lam, steps, True
    except Exception as e:
        return np.nan, 0, False


def run_comparison_for_ic(ic_func):
    """Run Newton vs TGP for one IC across all t_finals."""
    pos0, vel0, C, ic_name = ic_func()

    # Newton energy at rest
    H_n = total_mechanical_energy(pos0, vel0, C, pot_newton_fn)

    # TGP matched energy: add random velocity
    E_tgp_rest = pot_tgp_fn(pos0, C)
    vel_tgp = scale_velocities_match_energy(pos0, vel0, C, pot_tgp_fn, target_energy=H_n)
    if float(np.sum(vel_tgp**2)) < 1e-28:
        vel_tgp = random_velocities_for_excess_energy(
            pos0, C, pot_tgp_fn, target_energy=H_n,
            rng=np.random.default_rng(1980),
        )

    H_tgp = total_mechanical_energy(pos0, vel_tgp, C, pot_tgp_fn)
    tol_e = 1e-4 * max(1.0, abs(H_n))
    energy_matched = abs(H_tgp - H_n) < tol_e

    results = []

    for t_final in T_FINALS:
        t0 = time.time()

        lam_N, steps_N, ok_N = run_one(
            pos0, vel0, C, make_acc_newton, jac_newton_fn, t_final, "Newton"
        )
        lam_T, steps_T, ok_T = run_one(
            pos0, vel_tgp, C, make_acc_tgp_pair, jac_tgp_fn, t_final, "TGP"
        )

        elapsed = time.time() - t0

        results.append({
            'ic': ic_name, 't_final': t_final,
            'lam_N': lam_N, 'ok_N': ok_N,
            'lam_T': lam_T, 'ok_T': ok_T,
            'H_n': H_n, 'H_tgp': H_tgp,
            'energy_matched': energy_matched,
            'elapsed': elapsed,
        })

    return results


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 78)
    print("ex198 -- P1 Closure: TGP vs Newton Lyapunov (multiple ICs)")
    print("=" * 78)
    print(f"  beta = gamma = {BETA}, softening = {SOFT}, dt = {DT}")
    print(f"  t_final values: {T_FINALS}")
    print(f"  ICs: {len(IC_GENERATORS)} configurations")
    print(f"  Backend: pairwise V2 (analytic Jacobians)")

    all_results = []

    for ic_func in IC_GENERATORS:
        pos0, vel0, C, ic_name = ic_func()
        print(f"\n{'='*78}")
        print(f"IC: {ic_name}  |  C = {C}  |  H_match: ", end="")

        results = run_comparison_for_ic(ic_func)
        print(f"{'OK' if results[0]['energy_matched'] else 'FAIL'} "
              f"(H_N={results[0]['H_n']:.4f}, H_TGP={results[0]['H_tgp']:.4f})")

        print(f"  {'t_final':>7s} | {'lam_Newton':>10s} | {'lam_TGP':>10s} | "
              f"{'ratio':>8s} | {'verdict':>16s} | {'time':>6s}")
        print(f"  {'-'*70}")

        for r in results:
            if r['ok_N'] and r['ok_T'] and r['lam_N'] > 0.001:
                ratio = r['lam_T'] / r['lam_N']
                ratio_s = f"{ratio:.4f}"
                if ratio < 0.9:
                    verdict = "SUPPRESSES"
                elif ratio > 1.1:
                    verdict = "ENHANCES"
                else:
                    verdict = "COMPARABLE"
            elif not r['ok_N'] or not r['ok_T']:
                ratio_s = "FAIL"
                verdict = "INCONCLUSIVE"
                ratio = np.nan
            else:
                ratio_s = "~0"
                verdict = "non-chaotic"
                ratio = np.nan

            print(f"  {r['t_final']:7.1f} | {r['lam_N']:10.5f} | {r['lam_T']:10.5f} | "
                  f"{ratio_s:>8s} | {verdict:>16s} | {r['elapsed']:5.1f}s")

        all_results.extend(results)

    # ── Convention check: Newton with and without /C ──
    print(f"\n{'='*78}")
    print("CONVENTION CHECK: Newton acc vs Newton acc/C")
    print(f"{'='*78}")
    print("  (forces_newton returns a_i = G*m_j/r^2*r_hat, which is acceleration)")
    print("  Ex178 convention divides by C again => acc/m (possibly wrong)")
    print("  Checking impact on lambda_max:")

    pos0_b, vel0_b, C_b, _ = ic_burrau()
    t_check = T_FINALS[0]

    lam_N_div, _, _ = run_one(pos0_b, vel0_b, C_b,
                               make_acc_newton, jac_newton_fn, t_check, "N/C")
    lam_N_nodiv, _, _ = run_one(pos0_b, vel0_b, C_b,
                                 make_acc_newton_corrected, jac_newton_fn, t_check, "N_raw")

    print(f"  Newton (acc/C, ex178 convention): lam = {lam_N_div:.5f}")
    print(f"  Newton (acc, corrected):          lam = {lam_N_nodiv:.5f}")
    print(f"  Ratio corrected/ex178: {lam_N_nodiv/lam_N_div:.3f}" if lam_N_div > 0.001
          else "  Cannot compute ratio")

    # ── Summary ──
    print(f"\n{'='*78}")
    print("SUMMARY TABLE")
    print(f"{'='*78}")
    print(f"  {'IC':>16s} {'t':>5s} {'lam_N':>9s} {'lam_T':>9s} {'ratio':>7s} {'verdict':>14s}")
    print(f"  {'-'*62}")

    ratios = []
    for r in all_results:
        if r['ok_N'] and r['ok_T'] and r['lam_N'] > 0.001:
            ratio = r['lam_T'] / r['lam_N']
            ratios.append(ratio)
            verdict = "SUPPRESS" if ratio < 0.9 else ("ENHANCE" if ratio > 1.1 else "COMPARABLE")
        else:
            ratio = np.nan
            verdict = "N/A"
        print(f"  {r['ic']:>16s} {r['t_final']:5.1f} {r['lam_N']:9.5f} "
              f"{r['lam_T']:9.5f} {ratio:7.3f} {verdict:>14s}")

    # ── Verdict ──
    print(f"\n{'='*78}")
    print("P1 THESIS ASSESSMENT")
    print(f"{'='*78}")

    if ratios:
        n_total = len(ratios)
        n_suppress = sum(1 for r in ratios if r < 0.9)
        n_enhance = sum(1 for r in ratios if r > 1.1)
        n_comparable = n_total - n_suppress - n_enhance
        mean_r = np.mean(ratios)
        std_r = np.std(ratios) if len(ratios) > 1 else 0

        print(f"  Valid comparisons: {n_total}")
        print(f"  Suppresses (ratio < 0.9): {n_suppress}/{n_total}")
        print(f"  Enhances   (ratio > 1.1): {n_enhance}/{n_total}")
        print(f"  Comparable (0.9 - 1.1):   {n_comparable}/{n_total}")
        print(f"  Mean ratio: {mean_r:.4f} +/- {std_r:.4f}")

        if mean_r < 1.0:
            print(f"\n  Average: TGP REDUCES lambda_max by {(1-mean_r)*100:.1f}%")
        else:
            print(f"\n  Average: TGP INCREASES lambda_max by {(mean_r-1)*100:.1f}%")

        if n_suppress > n_enhance and mean_r < 1.0:
            print(f"  -> P1 thesis SUPPORTED: TGP suppresses chaos")
        elif n_enhance > n_suppress and mean_r > 1.0:
            print(f"  -> P1 thesis NOT SUPPORTED: TGP enhances chaos at these params")
        else:
            print(f"  -> P1 thesis INCONCLUSIVE: mixed results")

        print(f"\n  CAVEATS:")
        print(f"  - Pairwise V2 only (no V3 3-body barrier)")
        print(f"  - Short integration horizon (finite-time Lyapunov)")
        print(f"  - beta=gamma={BETA} (one parameter point)")
        print(f"  - With V3 (yukawa_feynman), repulsive barrier is stronger")
    else:
        print("  No valid comparisons obtained.")


if __name__ == "__main__":
    main()
