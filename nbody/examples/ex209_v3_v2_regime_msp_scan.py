#!/usr/bin/env python3
"""
ex209 -- V3/V2 ratio across m_sp regimes
==========================================

Validates when 3-body forces (V3) are significant vs negligible relative
to pairwise (V2) as a function of the dimensionless parameter t = m_sp * d.

Physics:
  - m_sp = sqrt(3*gamma - 2*beta); for beta=gamma: m_sp = sqrt(beta)
  - V2 ~ pairwise Yukawa + gradient + confining terms
  - V3 = (2*beta - 6*gamma)*C1*C2*C3*I_Y(d;m_sp) (irreducible 3-body overlap)
  - t << 1: Coulomb regime (physical TGP), V3/V2 saturates (constant)
  - t >> 1: Yukawa regime (unit-code TGP), V3/V2 exponentially suppressed

Tests:
  Part 1: V3/V2 vs t = m_sp*d for equilateral triangle (systematic scan)
  Part 2: V3/V2 vs beta at fixed d (connects to Lyapunov ex200/ex207)
  Part 3: Geometry dependence (equilateral vs isosceles vs linear)
  Part 4: C-dependence of V3/V2 (linear in C for equilateral)
  Part 5: CSV export of Part 1 + Part 2

Uses exact Feynman 2D integral (three_body_force_exact), not saddle-point.
"""

import sys
import os
import time
import csv
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.three_body_force_exact import (
    three_body_energy_exact,
    yukawa_overlap_exact,
    triangle_shape_coordinates,
    yukawa_overlap_geometry_rate,
)
from nbody.dynamics_v2 import potential_tgp
from nbody.tgp_field import screening_mass

quick = "--quick" in sys.argv

SOFT = 1e-6

# ── Part 1: V3/V2 vs t = m_sp * d ─────────────────────────────────────────

def part1_v3v2_vs_msp_d():
    """Scan V3/V2 ratio vs dimensionless parameter t = m_sp*d."""
    print("=" * 78)
    print("Part 1: V3/V2 vs t = m_sp*d (equilateral triangle)")
    print("=" * 78)

    C = 0.20
    C_arr = np.array([C, C, C])

    if quick:
        t_vals = np.array([0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
        n_quad = 12
    else:
        t_vals = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0,
                           1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])
        n_quad = 18

    results = []

    print(f"\n  C = {C}, equilateral d12=d13=d23=d")
    print(f"  n_quad = {n_quad}")
    print(f"\n  {'t=m*d':>8s} {'beta':>8s} {'d':>8s} "
          f"{'|V2|':>12s} {'|V3|':>12s} {'|V3/V2|':>12s} {'V3 sign':>8s}")
    print(f"  {'-'*72}")

    for t in t_vals:
        # Fix m_sp=1 and vary d=t (equivalent to fixing d=1 and varying m_sp=t)
        # With m_sp=1: beta=gamma=1 (since m_sp = sqrt(beta) for beta=gamma)
        m_sp = 1.0
        d = float(t)
        beta = m_sp ** 2  # = 1.0
        gamma = beta

        # Equilateral positions
        pos = np.array([
            [0.0, 0.0, 0.0],
            [d,   0.0, 0.0],
            [d/2, d*np.sqrt(3)/2, 0.0],
        ])

        V2 = potential_tgp(pos, C_arr, beta, gamma, SOFT)
        V3 = three_body_energy_exact(d, d, d, C, C, C,
                                      beta=beta, gamma=gamma, n_quad=n_quad)

        ratio = abs(V3 / V2) if abs(V2) > 1e-30 else float('nan')
        sign = "+" if V3 > 0 else "-"

        results.append({
            't': t, 'beta': beta, 'd': d,
            'V2': V2, 'V3': V3, 'ratio': ratio,
        })

        print(f"  {t:8.3f} {beta:8.4f} {d:8.4f} "
              f"{abs(V2):12.6e} {abs(V3):12.6e} {ratio:12.6e} {sign:>8s}")

    # Analysis
    if len(results) >= 2:
        r_small = results[0]['ratio']
        r_large = results[-1]['ratio']
        suppression = r_large / r_small if r_small > 1e-30 else float('nan')
        print(f"\n  Suppression factor (t={t_vals[-1]} vs t={t_vals[0]}): "
              f"{suppression:.4e}")

        # Find where V3/V2 crosses 1%
        for r in results:
            if r['ratio'] < 0.01:
                print(f"  V3/V2 < 1% at t = {r['t']:.2f}")
                break
        else:
            print(f"  V3/V2 > 1% for all t tested")

    return results


# ── Part 2: V3/V2 vs beta at fixed d ─────────────────────────────────────

def part2_v3v2_vs_beta():
    """Scan V3/V2 vs beta at fixed geometry (Burrau-like separations)."""
    print(f"\n{'=' * 78}")
    print("Part 2: V3/V2 vs beta (equilateral, d=3.0, C=0.20)")
    print(f"{'=' * 78}")

    C = 0.20
    C_arr = np.array([C, C, C])
    d = 3.0

    if quick:
        betas = [0.01, 0.05, 0.10, 0.20, 0.50, 1.0]
        n_quad = 12
    else:
        betas = [0.005, 0.01, 0.025, 0.05, 0.08, 0.10, 0.15,
                 0.20, 0.35, 0.50, 0.75, 1.0, 2.0, 5.0]
        n_quad = 18

    pos = np.array([
        [0.0, 0.0, 0.0],
        [d,   0.0, 0.0],
        [d/2, d*np.sqrt(3)/2, 0.0],
    ])

    results = []

    print(f"\n  d = {d}, C = {C}")
    print(f"\n  {'beta':>8s} {'m_sp':>8s} {'t=m*d':>8s} "
          f"{'|V2|':>12s} {'|V3|':>12s} {'|V3/V2|':>12s}")
    print(f"  {'-'*62}")

    for beta in betas:
        gamma = beta
        m_sp = screening_mass(beta, gamma)
        t = m_sp * d

        V2 = potential_tgp(pos, C_arr, beta, gamma, SOFT)
        V3 = three_body_energy_exact(d, d, d, C, C, C,
                                      beta=beta, gamma=gamma, n_quad=n_quad)

        ratio = abs(V3 / V2) if abs(V2) > 1e-30 else float('nan')

        results.append({
            'beta': beta, 'm_sp': m_sp, 't': t,
            'V2': V2, 'V3': V3, 'ratio': ratio,
        })

        print(f"  {beta:8.4f} {m_sp:8.4f} {t:8.3f} "
              f"{abs(V2):12.6e} {abs(V3):12.6e} {ratio:12.6e}")

    # Analysis: identify regime transition
    if len(results) >= 3:
        ratios = [r['ratio'] for r in results]
        t_vals = [r['t'] for r in results]

        # Find t where V3/V2 drops below 1%
        for i, r in enumerate(results):
            if r['ratio'] < 0.01 and r['t'] > 0.5:
                print(f"\n  V3/V2 < 1% at beta={r['beta']:.3f} "
                      f"(m_sp*d={r['t']:.2f})")
                break

        # Max ratio
        i_max = int(np.argmax(ratios))
        print(f"  Max |V3/V2| = {ratios[i_max]:.4e} at beta={results[i_max]['beta']}")

    return results


# ── Part 3: Geometry dependence ───────────────────────────────────────────

def part3_geometry():
    """Compare V3/V2 for different triangle geometries."""
    print(f"\n{'=' * 78}")
    print("Part 3: V3/V2 geometry dependence")
    print(f"{'=' * 78}")

    C = 0.20
    beta = 0.08
    gamma = beta
    m_sp = screening_mass(beta, gamma)
    n_quad = 12 if quick else 18

    # Different geometries at similar scale
    d_base = 3.0
    geometries = {
        "equilateral": (d_base, d_base, d_base),
        "isosceles_2:1": (d_base, d_base, d_base/2),
        "isosceles_1:2": (d_base, d_base, 2*d_base),
        "linear_1:1": (d_base, d_base, 2*d_base),  # collinear
        "compact": (d_base/2, d_base/2, d_base/2),
        "extended": (2*d_base, 2*d_base, 2*d_base),
    }

    print(f"\n  beta = gamma = {beta}, m_sp = {m_sp:.4f}, n_quad = {n_quad}")
    print(f"\n  {'Geometry':>20s} {'q1':>6s} {'q2':>6s} {'lambda':>8s} {'t_max':>8s} "
          f"{'|V3/V2|':>12s}")
    print(f"  {'-'*74}")

    results = {}
    for name, (d12, d13, d23) in geometries.items():
        # Check triangle inequality
        if (d12 + d13 <= d23 or d12 + d23 <= d13 or d13 + d23 <= d12):
            # Degenerate/collinear — V3 should still be computable
            pass

        shape = triangle_shape_coordinates(d12, d13, d23, m=m_sp)
        rate = yukawa_overlap_geometry_rate(d12, d13, d23, m=m_sp)

        # Build positions (in-plane)
        # Body 1 at origin, body 2 at (d12,0,0), body 3 determined by d13,d23
        cos_angle = (d12**2 + d13**2 - d23**2) / (2*d12*d13 + 1e-30)
        cos_angle = np.clip(cos_angle, -1, 1)
        sin_angle = np.sqrt(max(0, 1 - cos_angle**2))
        pos = np.array([
            [0.0, 0.0, 0.0],
            [d12, 0.0, 0.0],
            [d13*cos_angle, d13*sin_angle, 0.0],
        ])
        C_arr = np.array([C, C, C])

        V2 = potential_tgp(pos, C_arr, beta, gamma, SOFT)
        V3 = three_body_energy_exact(d12, d13, d23, C, C, C,
                                      beta=beta, gamma=gamma, n_quad=n_quad)

        ratio = abs(V3 / V2) if abs(V2) > 1e-30 else float('nan')
        results[name] = {
            'ratio': ratio,
            't_max': shape['t'],
            'q1': shape['q1'],
            'q2': shape['q2'],
            'lambda': rate['lambda'],
            'suppression_exponent': rate['suppression_exponent'],
            'V2': V2,
            'V3': V3,
        }

        print(f"  {name:>20s} {shape['q1']:6.3f} {shape['q2']:6.3f} "
              f"{rate['lambda']:8.3f} {shape['t']:8.3f} {ratio:12.6e}")

    # Verify: compact geometry has larger V3/V2 (closer bodies = stronger overlap)
    if 'compact' in results and 'extended' in results:
        r_compact = results['compact']['ratio']
        r_extended = results['extended']['ratio']
        print(f"\n  compact/extended ratio: {r_compact/r_extended:.2f}x "
              f"({'OK: compact > extended' if r_compact > r_extended else 'UNEXPECTED'})")

    return results


# ── Part 4: C-dependence ─────────────────────────────────────────────────

def part4_C_dependence():
    """V3/V2 scales linearly with C for equal-mass equilateral."""
    print(f"\n{'=' * 78}")
    print("Part 4: V3/V2 vs C (equilateral, beta=0.08, d=3.0)")
    print(f"{'=' * 78}")

    beta = 0.08
    gamma = beta
    d = 3.0
    n_quad = 12 if quick else 18

    if quick:
        C_vals = [0.05, 0.10, 0.20, 0.40]
    else:
        C_vals = [0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0]

    pos = np.array([
        [0.0, 0.0, 0.0],
        [d,   0.0, 0.0],
        [d/2, d*np.sqrt(3)/2, 0.0],
    ])

    print(f"\n  {'C':>8s} {'|V2|':>12s} {'|V3|':>12s} "
          f"{'|V3/V2|':>12s} {'V3/(C*ref)':>12s}")
    print(f"  {'-'*56}")

    results = []
    ref_ratio_per_C = None

    for C in C_vals:
        C_arr = np.array([C, C, C])
        V2 = potential_tgp(pos, C_arr, beta, gamma, SOFT)
        V3 = three_body_energy_exact(d, d, d, C, C, C,
                                      beta=beta, gamma=gamma, n_quad=n_quad)

        ratio = abs(V3 / V2) if abs(V2) > 1e-30 else float('nan')

        # For equilateral equal-mass: V3 ~ C^3, V2 ~ C^2, so V3/V2 ~ C
        ratio_per_C = ratio / C if C > 0 else float('nan')
        if ref_ratio_per_C is None:
            ref_ratio_per_C = ratio_per_C

        norm = ratio_per_C / ref_ratio_per_C if ref_ratio_per_C > 1e-30 else float('nan')

        results.append({'C': C, 'ratio': ratio, 'ratio_per_C': ratio_per_C})
        print(f"  {C:8.4f} {abs(V2):12.6e} {abs(V3):12.6e} "
              f"{ratio:12.6e} {norm:12.6f}")

    # Check linearity
    if len(results) >= 2:
        norms = [r['ratio_per_C'] / ref_ratio_per_C for r in results]
        max_dev = max(abs(n - 1.0) for n in norms)
        print(f"\n  V3/V2 ~ C linearity: max deviation = {max_dev:.4e} "
              f"({'PASS: linear' if max_dev < 0.05 else 'FAIL: nonlinear'})")

    return results


# ── Part 5: CSV export ───────────────────────────────────────────────────

def part5_csv(p1_results, p2_results):
    """Export Part 1 and Part 2 data."""
    print(f"\n{'=' * 78}")
    print("Part 5: CSV export")
    print(f"{'=' * 78}")

    outdir = os.path.join(os.path.dirname(__file__), "_outputs")
    os.makedirs(outdir, exist_ok=True)

    # Part 1 CSV
    f1 = os.path.join(outdir, "ex209_v3v2_vs_msp_d.csv")
    with open(f1, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_msp_d", "beta", "d", "V2", "V3", "abs_V3_over_V2"])
        for r in p1_results:
            w.writerow([f"{r['t']:.4f}", f"{r['beta']:.4f}",
                        f"{r['d']:.4f}", f"{r['V2']:.10e}",
                        f"{r['V3']:.10e}", f"{r['ratio']:.10e}"])
    print(f"  Written: {f1}")

    # Part 2 CSV
    f2 = os.path.join(outdir, "ex209_v3v2_vs_beta.csv")
    with open(f2, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["beta", "m_sp", "t_msp_d", "V2", "V3", "abs_V3_over_V2"])
        for r in p2_results:
            w.writerow([f"{r['beta']:.4f}", f"{r['m_sp']:.4f}",
                        f"{r['t']:.4f}", f"{r['V2']:.10e}",
                        f"{r['V3']:.10e}", f"{r['ratio']:.10e}"])
    print(f"  Written: {f2}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    print("=" * 78)
    print("ex209 -- V3/V2 ratio across m_sp regimes")
    print("         (when are 3-body forces significant?)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t0 = time.time()

    p1 = part1_v3v2_vs_msp_d()
    p2 = part2_v3v2_vs_beta()
    p3 = part3_geometry()
    p4 = part4_C_dependence()
    part5_csv(p1, p2)

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")

    # Key findings
    if p1:
        small_t = [r for r in p1 if r['t'] <= 0.5]
        large_t = [r for r in p1 if r['t'] >= 3.0]
        if small_t:
            avg_small = np.mean([r['ratio'] for r in small_t])
            print(f"  t << 1 (Coulomb): <|V3/V2|> = {avg_small:.4e}")
        if large_t:
            avg_large = np.mean([r['ratio'] for r in large_t])
            print(f"  t >> 1 (Yukawa):  <|V3/V2|> = {avg_large:.4e}")

    if p2:
        # Find beta range used in ex200/ex207
        lyap_betas = [r for r in p2 if 0.02 <= r['beta'] <= 0.35]
        if lyap_betas:
            min_r = min(r['ratio'] for r in lyap_betas)
            max_r = max(r['ratio'] for r in lyap_betas)
            print(f"  Lyapunov range (beta=0.02-0.35): |V3/V2| = "
                  f"{min_r:.4e} -- {max_r:.4e}")
            print(f"    -> V3 is {'SIGNIFICANT' if max_r > 0.01 else 'NEGLIGIBLE'} "
                  f"in Lyapunov scan regime")

    print(f"\n  Regime classification:")
    print(f"    t < 0.5:  Coulomb regime — V3/V2 approx constant, V3 perturbative")
    print(f"    t ~ 1:    Transition — V3 still measurable")
    print(f"    t > 3:    Yukawa screening — V3 exponentially suppressed")

    print(f"\n  Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
