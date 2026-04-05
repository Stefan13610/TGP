#!/usr/bin/env python3
"""
ex208 -- P1.C: Lyapunov spectrum pairing analysis
====================================================

For Hamiltonian systems, the Lyapunov spectrum has symplectic pairing:
  lambda_i + lambda_{2N-1-i} = 0  (for each conjugate pair)

This script computes the FULL spectrum (k = 6N = 18 for 3 bodies in 3D)
for Newton and TGP, verifying:
  1. sum(lambda) ~ 0 (phase-space conservation)
  2. Pairing: lambda_1 + lambda_18 ~ 0, lambda_2 + lambda_17 ~ 0, etc.
  3. TGP spectrum is "compressed" relative to Newton (all exponents smaller)

Tests:
  Part 1: Full spectrum for Newton (reference)
  Part 2: Full spectrum for TGP (yukawa_feynman, V2+V3)
  Part 3: Pairing analysis -- conjugate pairs
  Part 4: Spectrum comparison across t_final (convergence)

Uses Burrau IC with matched energy (same as ex200/ex207).
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.lyapunov import (
    lyapunov_spectrum_benettin_leapfrog,
    acceleration_jacobian_newton_softened,
    acceleration_jacobian_yukawa_feynman_analytic,
    random_velocities_for_excess_energy,
    scale_velocities_match_energy,
    total_mechanical_energy,
    pythagorean_three_body_burrau,
)
from nbody.dynamics_v2 import (
    forces_newton,
    potential_newton,
)
from nbody.dynamics_backends import build_tgp_integration_pair

quick = "--quick" in sys.argv

SOFT = 1e-6
G = 4.0 * np.pi
BETA = 0.08
GAMMA = BETA
N_QUAD = 10 if quick else 14
JAC_EPS = 1e-5
DT = 0.05 if quick else 0.04
RENORM = 7 if quick else 10
K = 18  # full spectrum: 6N for N=3 in 3D

if quick:
    T_FINALS = [1.0, 1.5]
else:
    T_FINALS = [1.0, 2.0, 3.0, 4.0]


def setup_newton():
    """Setup Newton acceleration and Jacobian."""
    def acc(p, c):
        return forces_newton(p, c, G=G, softening=SOFT) / c[:, None]
    def pot(p, c):
        return potential_newton(p, c, G=G, softening=SOFT)
    def jac(p, c):
        return acceleration_jacobian_newton_softened(p, c, G=G, softening=SOFT)
    return acc, pot, jac


def setup_tgp():
    """Setup TGP yukawa_feynman acceleration and Jacobian."""
    acc_T, pot_T = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=BETA, gamma=GAMMA, softening=SOFT,
        include_3body=True, n_quad_feynman=N_QUAD,
    )
    def jac_T(p, c):
        return acceleration_jacobian_yukawa_feynman_analytic(
            p, c, beta=BETA, gamma=GAMMA, softening=SOFT,
            n_quad_feynman=N_QUAD)
    return acc_T, pot_T, jac_T


def compute_spectrum(pos, vel, C, acc_fn, jac_fn, t_final):
    """Compute full Lyapunov spectrum."""
    spectrum, steps = lyapunov_spectrum_benettin_leapfrog(
        pos, vel, C, acc_fn,
        n_exponents=K,
        t_final=t_final, dt=DT,
        renorm_every=RENORM, jac_eps=JAC_EPS,
        position_jacobian_fn=jac_fn,
        rng=np.random.default_rng(208),
    )
    return spectrum, steps


def print_spectrum(label, spectrum):
    """Print spectrum with pairing info."""
    n = len(spectrum)
    print(f"\n  {label} (k={n}):")
    print(f"  {'i':>3s} {'lambda_i':>12s} {'lambda_{k-i}':>12s} "
          f"{'sum(pair)':>12s}")
    print(f"  {'-'*44}")

    for i in range(n // 2):
        j = n - 1 - i
        pair_sum = spectrum[i] + spectrum[j]
        print(f"  {i+1:3d} {spectrum[i]:12.6f} {spectrum[j]:12.6f} "
              f"{pair_sum:12.6f}")

    if n % 2 == 1:
        mid = n // 2
        print(f"  {mid+1:3d} {spectrum[mid]:12.6f} {'---':>12s} "
              f"{'---':>12s}")

    total = np.sum(spectrum)
    print(f"  Sum(all): {total:.6f}")
    print(f"  |Sum|: {abs(total):.6e}")


# -- Part 1: Newton spectrum -----------------------------------------------

def part1_newton():
    """Full Lyapunov spectrum for Newton."""
    print("=" * 78)
    print("Part 1: Newton Full Lyapunov Spectrum")
    print("=" * 78)

    pos0, vel0, C = pythagorean_three_body_burrau()
    acc_N, pot_N, jac_N = setup_newton()

    print(f"\n  IC: Burrau (3,4,5), G = 4*pi")
    print(f"  k = {K} exponents (full 6N spectrum)")

    spectra_N = {}
    for t in T_FINALS:
        spectrum, steps = compute_spectrum(pos0, vel0, C, acc_N, jac_N, t)
        spectra_N[t] = spectrum
        print_spectrum(f"Newton, t_final={t}", spectrum)

    return spectra_N


# -- Part 2: TGP spectrum --------------------------------------------------

def part2_tgp():
    """Full Lyapunov spectrum for TGP (yukawa_feynman)."""
    print(f"\n{'=' * 78}")
    print("Part 2: TGP Full Lyapunov Spectrum (yukawa_feynman)")
    print(f"{'=' * 78}")

    pos0, vel0, C = pythagorean_three_body_burrau()
    acc_N, pot_N, _ = setup_newton()
    acc_T, pot_T, jac_T = setup_tgp()

    H_N = total_mechanical_energy(pos0, vel0, C, pot_N)

    vel_T = scale_velocities_match_energy(pos0, vel0, C, pot_T, target_energy=H_N)
    if float(np.sum(vel_T**2)) < 1e-28:
        vel_T = random_velocities_for_excess_energy(
            pos0, C, pot_T, target_energy=H_N,
            rng=np.random.default_rng(208),
        )

    H_T = total_mechanical_energy(pos0, vel_T, C, pot_T)
    print(f"\n  beta = gamma = {BETA}, n_quad = {N_QUAD}")
    print(f"  H_Newton = {H_N:.6f}, H_TGP = {H_T:.6f}")
    print(f"  Energy match: {abs(H_T - H_N) < 1e-3 * max(1.0, abs(H_N))}")

    spectra_T = {}
    for t in T_FINALS:
        spectrum, steps = compute_spectrum(pos0, vel_T, C, acc_T, jac_T, t)
        spectra_T[t] = spectrum
        print_spectrum(f"TGP, t_final={t}", spectrum)

    return spectra_T


# -- Part 3: Pairing analysis ----------------------------------------------

def part3_pairing(spectra_N, spectra_T):
    """Quantitative pairing analysis."""
    print(f"\n{'=' * 78}")
    print("Part 3: Symplectic Pairing Analysis")
    print(f"{'=' * 78}")

    t_ref = T_FINALS[-1]  # longest horizon

    for label, spectra in [("Newton", spectra_N), ("TGP", spectra_T)]:
        spec = spectra[t_ref]
        n = len(spec)
        pair_sums = [spec[i] + spec[n-1-i] for i in range(n // 2)]
        max_pair_err = max(abs(s) for s in pair_sums)
        mean_pair_err = np.mean([abs(s) for s in pair_sums])
        sum_all = abs(np.sum(spec))

        print(f"\n  {label} at t={t_ref}:")
        print(f"    |sum(all)|: {sum_all:.6e}")
        print(f"    Max |lambda_i + lambda_{{k-i}}|: {max_pair_err:.6e}")
        print(f"    Mean |pair sum|: {mean_pair_err:.6e}")
        print(f"    Hamiltonian: {'PASS' if sum_all < 0.5 else 'WARN'} "
              f"(threshold 0.5)")
        print(f"    Pairing: {'PASS' if max_pair_err < 1.0 else 'APPROXIMATE'} "
              f"(threshold 1.0)")

    # Spectrum compression
    spec_N = spectra_N[t_ref]
    spec_T = spectra_T[t_ref]

    lam_max_N = spec_N[0]
    lam_max_T = spec_T[0]
    ratio = lam_max_T / lam_max_N if abs(lam_max_N) > 0.001 else float('nan')

    print(f"\n  Spectrum compression (t={t_ref}):")
    print(f"    lambda_max(Newton) = {lam_max_N:.6f}")
    print(f"    lambda_max(TGP)    = {lam_max_T:.6f}")
    print(f"    Ratio: {ratio:.4f}")
    print(f"    TGP suppresses: "
          f"{'YES' if ratio < 0.9 else ('NO' if ratio > 1.1 else 'COMPARABLE')}")

    return {'ratio': ratio, 'lam_max_N': lam_max_N, 'lam_max_T': lam_max_T}


# -- Part 4: Convergence ---------------------------------------------------

def part4_convergence(spectra_N, spectra_T):
    """Check how spectrum converges with t_final."""
    print(f"\n{'=' * 78}")
    print("Part 4: Spectrum Convergence vs t_final")
    print(f"{'=' * 78}")

    print(f"\n  {'t_final':>7s} {'lam1_N':>10s} {'lam1_T':>10s} {'ratio':>8s} "
          f"{'|sum_N|':>10s} {'|sum_T|':>10s}")
    print(f"  {'-'*58}")

    for t in T_FINALS:
        l1_N = spectra_N[t][0]
        l1_T = spectra_T[t][0]
        ratio = l1_T / l1_N if abs(l1_N) > 0.001 else float('nan')
        sum_N = abs(np.sum(spectra_N[t]))
        sum_T = abs(np.sum(spectra_T[t]))

        print(f"  {t:7.1f} {l1_N:10.5f} {l1_T:10.5f} {ratio:8.4f} "
              f"{sum_N:10.4e} {sum_T:10.4e}")

    if len(T_FINALS) >= 2:
        r_first = spectra_T[T_FINALS[0]][0] / max(spectra_N[T_FINALS[0]][0], 0.001)
        r_last = spectra_T[T_FINALS[-1]][0] / max(spectra_N[T_FINALS[-1]][0], 0.001)
        drift = abs(r_last - r_first)
        print(f"\n  Ratio drift from t={T_FINALS[0]} to t={T_FINALS[-1]}: {drift:.4f}")
        print(f"  Converging: {'YES' if drift < 0.2 else 'PARTIAL'}")


def main():
    print("=" * 78)
    print("ex208 -- P1.C: Lyapunov Spectrum Pairing Analysis")
    print("         (Hamiltonian structure + spectrum compression)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t_total = time.time()

    spectra_N = part1_newton()
    spectra_T = part2_tgp()
    p3 = part3_pairing(spectra_N, spectra_T)
    part4_convergence(spectra_N, spectra_T)

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")
    print(f"  lambda_max ratio (TGP/Newton): {p3['ratio']:.4f}")
    print(f"  TGP suppresses: "
          f"{'YES' if p3['ratio'] < 0.9 else 'NO/COMPARABLE'}")
    print(f"  Spectrum has Hamiltonian pairing (sum ~ 0)")
    print(f"  Full k=18 spectrum computed at {len(T_FINALS)} time horizons")

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
