#!/usr/bin/env python3
"""
ex202 -- P4: Multipole I_Y vs Feynman 2D (accuracy + speed benchmark)
======================================================================

Verifies the multipole (Legendre/Gegenbauer) expansion of the triple
Yukawa overlap integral I_Y against the exact Feynman parametrization
from three_body_force_exact.py.

Tests:
  Part 1: I_Y value comparison on a grid of geometries
  Part 2: Convergence with L_max
  Part 3: Force comparison (dI/dd) — multipole FD vs Feynman analytic
  Part 4: Speed benchmark (multipole vs Feynman)
  Part 5: 3-body forces on Burrau IC — multipole vs exact
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.multipole_triple_overlap import (
    yukawa_overlap_multipole,
    yukawa_overlap_multipole_with_derivatives,
    three_body_forces_multipole,
    three_body_potential_multipole,
)
from nbody.three_body_force_exact import (
    yukawa_overlap_exact,
    _dI_dd_components,
    three_body_forces_exact,
)
from nbody.tgp_field import default_beta_gamma, screening_mass

quick = "--quick" in sys.argv

BETA = 0.08
GAMMA = BETA
_, GAMMA_C = default_beta_gamma(BETA, GAMMA)
M = screening_mass(BETA, GAMMA_C)

if quick:
    GEOMETRIES = [
        ("equilateral_small", 1.0, 1.0, 1.0),
        ("equilateral_large", 5.0, 5.0, 5.0),
        ("isosceles", 2.0, 3.0, 3.0),
        ("scalene", 1.5, 2.5, 3.0),
    ]
    L_MAXES = [3, 5, 10, 15]
    N_QUAD_REF = 30
    N_RAD = 40
    N_BENCH = 5
else:
    GEOMETRIES = [
        ("equilateral_tiny", 0.3, 0.3, 0.3),
        ("equilateral_small", 1.0, 1.0, 1.0),
        ("equilateral_med", 3.0, 3.0, 3.0),
        ("equilateral_large", 5.0, 5.0, 5.0),
        ("isosceles_narrow", 1.0, 5.0, 5.0),
        ("isosceles_wide", 4.0, 3.0, 3.0),
        ("scalene", 1.5, 2.5, 3.0),
        ("nearly_collinear", 1.0, 2.0, 2.99),
        ("compact", 0.5, 0.6, 0.7),
    ]
    L_MAXES = [1, 2, 3, 5, 8, 10, 15, 20]
    N_QUAD_REF = 50
    N_RAD = 50
    N_BENCH = 20


def part1_value_comparison():
    """Compare I_Y values: multipole vs Feynman."""
    print("=" * 78)
    print("Part 1: I_Y Value Comparison (multipole vs Feynman)")
    print("=" * 78)
    print(f"  m = {M:.4f}, beta = {BETA}, gamma = {GAMMA}")
    print(f"  Feynman: n_quad = {N_QUAD_REF}")
    print(f"  Multipole: L_max = 10, n_rad = {N_RAD}")

    print(f"\n  {'geometry':>20s} {'d12':>5s} {'d13':>5s} {'d23':>5s} "
          f"{'I_feynman':>12s} {'I_multipole':>12s} {'rel_err':>10s} {'status':>8s}")
    print(f"  {'-'*88}")

    results = []
    for name, d12, d13, d23 in GEOMETRIES:
        I_feyn = yukawa_overlap_exact(d12, d13, d23, M, n_quad=N_QUAD_REF)
        I_mult = yukawa_overlap_multipole(d12, d13, d23, M, L_max=10, n_rad=N_RAD)

        if abs(I_feyn) > 1e-30:
            rel_err = abs(I_mult - I_feyn) / abs(I_feyn)
        else:
            rel_err = abs(I_mult - I_feyn)

        status = "PASS" if rel_err < 0.05 else ("WARN" if rel_err < 0.20 else "FAIL")

        print(f"  {name:>20s} {d12:5.2f} {d13:5.2f} {d23:5.2f} "
              f"{I_feyn:12.6e} {I_mult:12.6e} {rel_err:10.2e} {status:>8s}")

        results.append({
            'name': name, 'd12': d12, 'd13': d13, 'd23': d23,
            'I_feyn': I_feyn, 'I_mult': I_mult, 'rel_err': rel_err,
        })

    return results


def part2_convergence():
    """Convergence of multipole with L_max."""
    print(f"\n{'=' * 78}")
    print("Part 2: Convergence with L_max")
    print(f"{'=' * 78}")

    # Pick a representative geometry
    test_cases = [
        ("equilateral d=1", 1.0, 1.0, 1.0),
        ("equilateral d=3", 3.0, 3.0, 3.0),
        ("scalene", 1.5, 2.5, 3.0),
    ]

    results = []
    for name, d12, d13, d23 in test_cases:
        I_ref = yukawa_overlap_exact(d12, d13, d23, M, n_quad=N_QUAD_REF)
        print(f"\n  {name}: d = ({d12},{d13},{d23}), m*d_min = {M*min(d12,d13,d23):.2f}")
        print(f"  I_ref (Feynman) = {I_ref:.8e}")
        print(f"  {'L_max':>6s} {'I_multipole':>14s} {'rel_err':>12s} {'status':>8s}")
        print(f"  {'-'*44}")

        for L in L_MAXES:
            I_m = yukawa_overlap_multipole(d12, d13, d23, M, L_max=L, n_rad=N_RAD)
            re = abs(I_m - I_ref) / max(abs(I_ref), 1e-30)
            status = "PASS" if re < 0.01 else ("OK" if re < 0.05 else "LOW")
            print(f"  {L:6d} {I_m:14.8e} {re:12.4e} {status:>8s}")
            results.append({'name': name, 'L_max': L, 'rel_err': re})

    return results


def part3_derivative_comparison():
    """Compare derivatives dI/dd: multipole FD vs Feynman analytic."""
    print(f"\n{'=' * 78}")
    print("Part 3: Derivative Comparison (dI/dd)")
    print(f"{'=' * 78}")

    print(f"\n  {'geometry':>20s} {'deriv':>6s} {'Feynman':>12s} {'Multipole':>12s} "
          f"{'rel_err':>10s} {'status':>8s}")
    print(f"  {'-'*76}")

    results = []
    test_geoms = GEOMETRIES[:4] if quick else GEOMETRIES[:6]

    for name, d12, d13, d23 in test_geoms:
        # Feynman analytic derivatives
        dI12_f, dI13_f, dI23_f = _dI_dd_components(d12, d13, d23, M, n_quad=N_QUAD_REF)

        # Multipole FD derivatives
        _, dI12_m, dI13_m, dI23_m = yukawa_overlap_multipole_with_derivatives(
            d12, d13, d23, M, L_max=10, n_rad=N_RAD,
        )

        for label, vf, vm in [('dd12', dI12_f, dI12_m),
                               ('dd13', dI13_f, dI13_m),
                               ('dd23', dI23_f, dI23_m)]:
            if abs(vf) > 1e-30:
                re = abs(vm - vf) / abs(vf)
            else:
                re = abs(vm - vf)
            status = "PASS" if re < 0.05 else ("WARN" if re < 0.20 else "FAIL")
            print(f"  {name:>20s} {label:>6s} {vf:12.6e} {vm:12.6e} {re:10.2e} {status:>8s}")
            results.append({'name': name, 'deriv': label, 'rel_err': re})

    return results


def part4_speed_benchmark():
    """Speed comparison: multipole vs Feynman."""
    print(f"\n{'=' * 78}")
    print("Part 4: Speed Benchmark")
    print(f"{'=' * 78}")

    d12, d13, d23 = 2.0, 2.5, 3.0

    # Feynman
    t0 = time.time()
    for _ in range(N_BENCH):
        yukawa_overlap_exact(d12, d13, d23, M, n_quad=N_QUAD_REF)
    t_feyn = (time.time() - t0) / N_BENCH

    print(f"\n  Feynman 2D (n_quad={N_QUAD_REF}): {t_feyn*1000:.2f} ms/eval")

    for L in [5, 10, 15]:
        t0 = time.time()
        for _ in range(N_BENCH):
            yukawa_overlap_multipole(d12, d13, d23, M, L_max=L, n_rad=N_RAD)
        t_mult = (time.time() - t0) / N_BENCH
        speedup = t_feyn / t_mult if t_mult > 0 else float('inf')
        print(f"  Multipole (L_max={L:2d}, n_rad={N_RAD}): {t_mult*1000:.2f} ms/eval  "
              f"(speedup: {speedup:.1f}x)")

    # With derivatives
    t0 = time.time()
    for _ in range(N_BENCH):
        _dI_dd_components(d12, d13, d23, M, n_quad=N_QUAD_REF)
    t_feyn_d = (time.time() - t0) / N_BENCH

    t0 = time.time()
    for _ in range(N_BENCH):
        yukawa_overlap_multipole_with_derivatives(d12, d13, d23, M, L_max=10, n_rad=N_RAD)
    t_mult_d = (time.time() - t0) / N_BENCH

    print(f"\n  With derivatives:")
    print(f"    Feynman analytic: {t_feyn_d*1000:.2f} ms")
    print(f"    Multipole FD:    {t_mult_d*1000:.2f} ms  "
          f"(ratio: {t_mult_d/t_feyn_d:.1f}x)")


def part5_force_comparison():
    """Compare full 3-body forces on Burrau IC."""
    print(f"\n{'=' * 78}")
    print("Part 5: 3-Body Forces on Burrau IC")
    print(f"{'=' * 78}")

    # Burrau IC
    pos = np.array([
        [1.0, 3.0, 0.0],
        [-2.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ])
    C = np.array([3.0, 4.0, 5.0])

    n_quad = N_QUAD_REF
    L_max = 10

    # Exact forces (Feynman)
    t0 = time.time()
    F_exact = three_body_forces_exact(pos, C, BETA, GAMMA, n_quad=n_quad)
    t_exact = time.time() - t0

    # Multipole forces
    t0 = time.time()
    F_mult = three_body_forces_multipole(pos, C, BETA, GAMMA, L_max=L_max, n_rad=N_RAD)
    t_mult = time.time() - t0

    print(f"\n  Forces (Feynman, n_quad={n_quad}):")
    for i in range(3):
        print(f"    body {i+1}: [{F_exact[i,0]:+10.6f}, {F_exact[i,1]:+10.6f}, {F_exact[i,2]:+10.6f}]")

    print(f"\n  Forces (Multipole, L_max={L_max}):")
    for i in range(3):
        print(f"    body {i+1}: [{F_mult[i,0]:+10.6f}, {F_mult[i,1]:+10.6f}, {F_mult[i,2]:+10.6f}]")

    # Relative error
    F_max = np.max(np.abs(F_exact))
    rel_err = np.max(np.abs(F_mult - F_exact)) / max(F_max, 1e-30)
    print(f"\n  Max relative error: {rel_err:.4e}")
    print(f"  Time: Feynman={t_exact*1000:.1f}ms, Multipole={t_mult*1000:.1f}ms")

    # Potential energy comparison
    V_exact = 0.0
    from itertools import combinations
    for i, j, k in combinations(range(3), 3):
        rij = pos[j] - pos[i]
        rik = pos[k] - pos[i]
        rjk = pos[k] - pos[j]
        dij = np.linalg.norm(rij)
        dik = np.linalg.norm(rik)
        djk = np.linalg.norm(rjk)
        I_Y = yukawa_overlap_exact(dij, dik, djk, M, n_quad=n_quad)
        V_exact += (2.0 * BETA - 6.0 * GAMMA_C) * C[i] * C[j] * C[k] * I_Y

    V_mult = three_body_potential_multipole(pos, C, BETA, GAMMA, L_max=L_max, n_rad=N_RAD)

    print(f"\n  V3 (Feynman):   {V_exact:.8e}")
    print(f"  V3 (Multipole): {V_mult:.8e}")
    print(f"  Relative error: {abs(V_mult - V_exact)/max(abs(V_exact), 1e-30):.4e}")

    # Newton's 3rd law check: sum of forces should be zero
    F_sum_exact = np.sum(F_exact, axis=0)
    F_sum_mult = np.sum(F_mult, axis=0)
    print(f"\n  Sum(F) exact:    [{F_sum_exact[0]:+.2e}, {F_sum_exact[1]:+.2e}, {F_sum_exact[2]:+.2e}]")
    print(f"  Sum(F) multipole: [{F_sum_mult[0]:+.2e}, {F_sum_mult[1]:+.2e}, {F_sum_mult[2]:+.2e}]")

    status = "PASS" if rel_err < 0.05 else "FAIL"
    print(f"\n  Force comparison: {status}")

    return {'rel_err': rel_err, 'status': status}


def main():
    print("=" * 78)
    print("ex202 -- P4: Multipole I_Y vs Feynman 2D (Triple Overlap)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")
    print(f"  m = sqrt(3*gamma - 2*beta) = {M:.4f}")
    print(f"  beta = {BETA}, gamma = {GAMMA}")

    t_total = time.time()

    p1 = part1_value_comparison()
    p2 = part2_convergence()
    p3 = part3_derivative_comparison()
    part4_speed_benchmark()
    p5 = part5_force_comparison()

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")

    n_pass_1 = sum(1 for r in p1 if r['rel_err'] < 0.05)
    print(f"  Part 1 (values): {n_pass_1}/{len(p1)} within 5%")

    # L_max needed for 1% accuracy
    for name in set(r['name'] for r in p2):
        entries = [r for r in p2 if r['name'] == name and r['rel_err'] < 0.01]
        if entries:
            L_min = min(r['L_max'] for r in entries)
            print(f"  Part 2 ({name}): L_max >= {L_min} for 1% accuracy")

    n_pass_3 = sum(1 for r in p3 if r['rel_err'] < 0.10)
    print(f"  Part 3 (derivatives): {n_pass_3}/{len(p3)} within 10%")
    print(f"  Part 5 (forces): {p5['status']} (rel_err = {p5['rel_err']:.2e})")

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
