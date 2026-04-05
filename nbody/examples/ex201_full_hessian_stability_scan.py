#!/usr/bin/env python3
"""
ex201 -- P2: Full Hessian Stability Mapping (TGP vs Newton)
============================================================

Systematic normal-mode analysis at TGP equilibrium configurations.
Compares all eigenvalues of the mass-weighted Hessian for Newton and TGP
across multiple geometries and a (beta, C) parameter grid.

Key questions:
  Q1: Which modes does TGP stabilize vs Newton? (breathing, bending, ...)
  Q2: At what beta does bifurcation occur (stability classification changes)?
  Q3: Does V3 (3-body) affect stability classification?
  Q4: How do results depend on C (source strength)?

Configurations tested:
  - Equilateral triangle (3-body, d_rep and d_well)
  - Collinear (Euler-type, 3-body)
  - Regular N-gon (4, 5, 6 bodies)

Uses: stability.normal_mode_analysis, stability.stability_comparison,
      stability.stability_bifurcation_scan
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.stability import (
    normal_mode_analysis,
    stability_comparison,
    stability_bifurcation_scan,
)
from nbody.equilibria import (
    equilateral_pairwise_equilibrium,
    collinear_pairwise_equilibrium_numerical,
    ngon_pairwise_equilibrium_numerical,
)
from nbody.configurations import (
    equilateral_triangle,
    collinear_equal,
    regular_ngon,
)
from nbody.dynamics_v2 import (
    potential_tgp,
    potential_newton,
    full_potential_tgp,
)

quick = "--quick" in sys.argv
SOFT = 1e-6
DX = 1e-5
G_NEWTON = 4.0 * np.pi  # Match TGP leading term

if quick:
    BETAS = [0.5, 1.0, 2.0, 4.0]
    CS = [0.05, 0.10, 0.20]
    NGONS = [5]
else:
    BETAS = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
    CS = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
    NGONS = [5, 7, 8]


# ── Helper: build equilibrium config for equilateral triangle ──────────

def equilateral_config(C, beta, gamma, which='d_rep'):
    """Return (positions, C_values) at pairwise equilibrium, or None."""
    eq = equilateral_pairwise_equilibrium(C, beta, gamma)
    if not eq['exists']:
        return None
    d = eq[which]
    if d is None or d <= 0:
        return None
    pos, Cv, _ = equilateral_triangle(d, C)
    return pos, Cv


def ngon_config(n, C, beta, gamma, which_eq=0):
    """Return (positions, C_values) at pairwise N-gon equilibrium, or None."""
    eq = ngon_pairwise_equilibrium_numerical(n, C, beta, gamma)
    if eq['n_equilibria'] == 0:
        return None
    idx = min(which_eq, eq['n_equilibria'] - 1)
    R_eq = eq['equilibria'][idx]['R']
    pos, Cv, _ = regular_ngon(n, R_eq, C)
    return pos, Cv


# ── Part 1: Detailed normal-mode analysis at a reference point ─────────

def part1_detailed_analysis():
    """Detailed mode-by-mode comparison at one (beta, C) point."""
    beta = 2.0
    C = 0.10
    gamma = beta

    print("=" * 78)
    print("Part 1: Detailed Normal-Mode Analysis")
    print("=" * 78)
    print(f"  beta = gamma = {beta}, C = {C}, G_Newton = 4*pi")

    eq = equilateral_pairwise_equilibrium(C, beta, gamma)
    if not eq['exists']:
        print(f"  No equilibrium at beta={beta}, C={C}")
        return {}

    results = {}

    for which, label in [('d_rep', 'Repulsive'), ('d_well', 'Well')]:
        d = eq[which]
        if d is None or d <= 0:
            continue

        print(f"\n  --- {label} equilibrium: d = {d:.4f} ---")
        pos, Cv, _ = equilateral_triangle(d, C)

        comp = stability_comparison(pos, Cv, beta, gamma,
                                     softening=SOFT, dx=DX)

        for theory, nma in [('Newton', comp['newton']), ('TGP', comp['tgp'])]:
            print(f"\n  {theory}:")
            print(f"    Classification: {nma['classification']}")
            print(f"    Zero modes: {nma['n_zero_modes']}")
            print(f"    Physical modes: {nma['n_stable']} stable + {nma['n_unstable']} unstable")
            print(f"    omega^2 = {nma['omega2_physical']}")
            if len(nma['omega_physical']) > 0:
                print(f"    omega   = {nma['omega_physical']}")

            # Mode characters
            print(f"    Mode details:")
            for k, mc in enumerate(nma['mode_characters']):
                stab_s = "STABLE" if mc['stable'] else "UNSTABLE"
                print(f"      [{k}] omega^2 = {mc['omega2']:.6f}, type = {mc['type']}, "
                      f"radial_frac = {mc['radial_fraction']:.3f}, {stab_s}")

        print(f"\n  Comparison:")
        print(f"    Stability differs: {comp['stability_differs']}")
        print(f"    Modes stabilized by TGP: {comp['n_modes_stabilized']}")
        print(f"    Modes destabilized by TGP: {comp['n_modes_destabilized']}")
        print(f"    omega^2 ratio (TGP/Newton): {comp['omega2_ratio']}")

        results[which] = comp

    return results


# ── Part 2: Beta scan at fixed C (equilateral triangle) ───────────────

def part2_beta_scan():
    """Scan beta for equilateral triangle at fixed C."""
    C = 0.10

    print(f"\n{'=' * 78}")
    print(f"Part 2: Beta Scan -- Equilateral Triangle (C = {C})")
    print(f"{'=' * 78}")

    print(f"\n  {'beta':>6s} {'d_rep':>7s} {'d_well':>7s} "
          f"{'N_class':>8s} {'T_class':>8s} "
          f"{'N_stab':>6s} {'T_stab':>6s} "
          f"{'stab+':>5s} {'stab-':>5s} {'time':>6s}")
    print(f"  {'-'*72}")

    results = []
    for beta in BETAS:
        gamma = beta
        eq = equilateral_pairwise_equilibrium(C, beta, gamma)
        if not eq['exists']:
            print(f"  {beta:6.2f}   --- no equilibrium (beta/C = {beta/C:.1f} < 4.5) ---")
            results.append({'beta': beta, 'has_eq': False})
            continue

        for which in ['d_rep', 'd_well']:
            d = eq[which]
            if d is None or d <= 0:
                continue

            t0 = time.time()
            pos, Cv, _ = equilateral_triangle(d, C)
            comp = stability_comparison(pos, Cv, beta, gamma,
                                         softening=SOFT, dx=DX)
            elapsed = time.time() - t0

            n_nma = comp['newton']
            t_nma = comp['tgp']

            d_rep_s = f"{eq['d_rep']:.4f}" if which == 'd_rep' else ""
            d_well_s = f"{eq['d_well']:.4f}" if which == 'd_well' else ""

            print(f"  {beta:6.2f} {d_rep_s:>7s} {d_well_s:>7s} "
                  f"{n_nma['classification']:>8s} {t_nma['classification']:>8s} "
                  f"{n_nma['n_stable']:>6d} {t_nma['n_stable']:>6d} "
                  f"{comp['n_modes_stabilized']:>5d} {comp['n_modes_destabilized']:>5d} "
                  f"{elapsed:5.1f}s")

            results.append({
                'beta': beta, 'which': which, 'd': d,
                'has_eq': True,
                'newton_class': n_nma['classification'],
                'tgp_class': t_nma['classification'],
                'newton_n_stable': n_nma['n_stable'],
                'tgp_n_stable': t_nma['n_stable'],
                'n_stabilized': comp['n_modes_stabilized'],
                'n_destabilized': comp['n_modes_destabilized'],
                'tgp_omega2': t_nma['omega2_physical'],
                'newton_omega2': n_nma['omega2_physical'],
            })

    return results


# ── Part 3: (beta, C) bifurcation scan ────────────────────────────────

def part3_bifurcation_scan():
    """Scan (beta, C) grid for equilateral triangle stability at d_well."""

    print(f"\n{'=' * 78}")
    print(f"Part 3: Bifurcation Map -- Equilateral Triangle d_well (beta, C) grid")
    print(f"{'=' * 78}")
    print(f"  (d_rep is always marginal/flat; d_well is the physically relevant one)")

    def config_d_well(C, beta, gamma):
        return equilateral_config(C, beta, gamma, 'd_well')

    t0 = time.time()
    results_well = stability_bifurcation_scan(
        config_d_well, BETAS, CS,
        softening=SOFT, dx=DX,
    )
    elapsed = time.time() - t0

    # Display d_well TGP classification
    header = f"  {'beta\\C':>8s}" + "".join(f" {C:>7.3f}" for C in CS)

    print(f"\n  d_well equilibrium -- TGP classification (time: {elapsed:.1f}s)")
    print(header)
    print(f"  {'-' * (8 + 8 * len(CS))}")

    idx = 0
    for beta in BETAS:
        row = f"  {beta:8.2f}"
        for C in CS:
            r = results_well[idx]; idx += 1
            if not r.get('has_equilibrium', False):
                row += f" {'---':>7s}"
            elif r.get('error', False):
                row += f" {'ERR':>7s}"
            else:
                cls = r['tgp_classification'][:3].upper()
                row += f" {cls:>7s}"
        print(row)

    # d_well Newton classification
    print(f"\n  d_well equilibrium -- Newton classification")
    print(header)
    print(f"  {'-' * (8 + 8 * len(CS))}")

    idx = 0
    for beta in BETAS:
        row = f"  {beta:8.2f}"
        for C in CS:
            r = results_well[idx]; idx += 1
            if not r.get('has_equilibrium', False):
                row += f" {'---':>7s}"
            elif r.get('error', False):
                row += f" {'ERR':>7s}"
            else:
                cls = r['newton_classification'][:3].upper()
                row += f" {cls:>7s}"
        print(row)

    # TGP stable mode count
    print(f"\n  d_well -- TGP stable mode count")
    print(header)
    print(f"  {'-' * (8 + 8 * len(CS))}")

    idx = 0
    for beta in BETAS:
        row = f"  {beta:8.2f}"
        for C in CS:
            r = results_well[idx]; idx += 1
            if not r.get('has_equilibrium', False):
                row += f" {'---':>7s}"
            elif r.get('error', False):
                row += f" {'ERR':>7s}"
            else:
                ns = r.get('tgp_n_stable', 0)
                row += f" {ns:>7d}"
        print(row)

    # Differences
    print(f"\n  Where TGP and Newton DIFFER (d_well):")
    differ_count = 0
    for r in results_well:
        if r.get('has_equilibrium') and not r.get('error') and r.get('stability_differs'):
            differ_count += 1
            print(f"    beta={r['beta']:.2f}, C={r['C']:.3f}: "
                  f"Newton={r['newton_classification']}, TGP={r['tgp_classification']}, "
                  f"modes stabilized={r['n_modes_stabilized']}, "
                  f"destabilized={r['n_modes_destabilized']}")

    if differ_count == 0:
        print(f"    (none -- classifications match everywhere)")

    return results_well


# ── Part 4: N-gon stability ───────────────────────────────────────────

def part4_ngon_stability():
    """Normal-mode analysis for regular N-gons at equilibrium."""
    C = 0.10
    beta = 2.0
    gamma = beta

    print(f"\n{'=' * 78}")
    print(f"Part 4: Regular N-gon Stability (C = {C}, beta = {beta})")
    print(f"{'=' * 78}")

    results = []
    for n in NGONS:
        eq = ngon_pairwise_equilibrium_numerical(n, C, beta, gamma)
        print(f"\n  --- {n}-gon: {eq['n_equilibria']} equilibria found ---")

        if eq['n_equilibria'] == 0:
            results.append({'n': n, 'has_eq': False})
            continue

        for eq_idx, eq_data in enumerate(eq['equilibria']):
            R_eq = eq_data['R']
            side = eq_data['side']
            t0 = time.time()

            pos, Cv, _ = regular_ngon(n, R_eq, C)
            comp = stability_comparison(pos, Cv, beta, gamma,
                                         softening=SOFT, dx=DX)
            elapsed = time.time() - t0

            n_nma = comp['newton']
            t_nma = comp['tgp']

            print(f"\n  Equilibrium {eq_idx}: R = {R_eq:.4f}, side = {side:.4f}")
            print(f"    Newton: {n_nma['classification']} "
                  f"({n_nma['n_stable']} stable, {n_nma['n_unstable']} unstable)")
            print(f"    TGP:    {t_nma['classification']} "
                  f"({t_nma['n_stable']} stable, {t_nma['n_unstable']} unstable)")
            print(f"    Modes stabilized by TGP: {comp['n_modes_stabilized']}")

            # Show mode characters for TGP
            print(f"    TGP modes:")
            for k, mc in enumerate(t_nma['mode_characters']):
                stab_s = "+" if mc['stable'] else "-"
                print(f"      [{k}] omega^2={mc['omega2']:+10.4f} "
                      f"{mc['type']:>11s} (radial={mc['radial_fraction']:.2f}) {stab_s}")

            results.append({
                'n': n, 'eq_idx': eq_idx, 'R': R_eq, 'side': side,
                'has_eq': True,
                'newton_class': n_nma['classification'],
                'tgp_class': t_nma['classification'],
                'tgp_modes': t_nma['mode_characters'],
                'n_stabilized': comp['n_modes_stabilized'],
                'elapsed': elapsed,
            })

    return results


# ── Part 5: Effect of V3 (3-body) on stability ────────────────────────

def part5_v3_effect():
    """Compare stability with and without V3 at select points."""
    C = 0.10
    beta = 2.0
    gamma = beta

    print(f"\n{'=' * 78}")
    print(f"Part 5: V3 (3-body) Effect on Stability (C = {C}, beta = {beta})")
    print(f"{'=' * 78}")

    eq = equilateral_pairwise_equilibrium(C, beta, gamma)
    if not eq['exists']:
        print("  No equilibrium")
        return {}

    results = {}
    for which in ['d_rep', 'd_well']:
        d = eq[which]
        if d is None or d <= 0:
            continue

        pos, Cv, _ = equilateral_triangle(d, C)

        print(f"\n  --- {which} (d = {d:.4f}) ---")

        # V2 only
        comp_v2 = stability_comparison(pos, Cv, beta, gamma,
                                        softening=SOFT, include_3body=False, dx=DX)
        # V2 + V3
        comp_v3 = stability_comparison(pos, Cv, beta, gamma,
                                        softening=SOFT, include_3body=True, dx=DX)

        t_v2 = comp_v2['tgp']
        t_v3 = comp_v3['tgp']

        print(f"    TGP (V2 only):  {t_v2['classification']} "
              f"({t_v2['n_stable']} stable, {t_v2['n_unstable']} unstable)")
        print(f"    TGP (V2 + V3):  {t_v3['classification']} "
              f"({t_v3['n_stable']} stable, {t_v3['n_unstable']} unstable)")

        print(f"    omega^2 (V2):  {t_v2['omega2_physical']}")
        print(f"    omega^2 (V3):  {t_v3['omega2_physical']}")

        if len(t_v2['omega2_physical']) > 0 and len(t_v3['omega2_physical']) > 0:
            n_common = min(len(t_v2['omega2_physical']), len(t_v3['omega2_physical']))
            shifts = t_v3['omega2_physical'][:n_common] - t_v2['omega2_physical'][:n_common]
            print(f"    delta omega^2 (V3 - V2): {shifts}")
            v3_stabilizes = np.sum(shifts > 0)
            v3_destabilizes = np.sum(shifts < 0)
            print(f"    V3 shifts {v3_stabilizes} modes toward stability, "
                  f"{v3_destabilizes} away")

        results[which] = {
            'v2_class': t_v2['classification'],
            'v3_class': t_v3['classification'],
            'v2_omega2': t_v2['omega2_physical'],
            'v3_omega2': t_v3['omega2_physical'],
        }

    return results


# ── Summary ────────────────────────────────────────────────────────────

def print_summary(p1, p2, p3, p4, p5):
    print(f"\n{'=' * 78}")
    print("SUMMARY: P2 Full Hessian Stability Analysis")
    print(f"{'=' * 78}")

    # From Part 2: count where TGP differs from Newton
    if p2:
        n_pts = len([r for r in p2 if r.get('has_eq')])
        n_stab = sum(1 for r in p2 if r.get('has_eq') and r.get('n_stabilized', 0) > 0)
        n_dstab = sum(1 for r in p2 if r.get('has_eq') and r.get('n_destabilized', 0) > 0)
        print(f"\n  Beta scan (equilateral, {n_pts} equilibrium points):")
        print(f"    TGP stabilizes modes at {n_stab} points")
        print(f"    TGP destabilizes modes at {n_dstab} points")

    # From Part 3: bifurcation count (d_well)
    if p3:
        n_valid = sum(1 for r in p3 if r.get('has_equilibrium') and not r.get('error'))
        n_differ = sum(1 for r in p3 if r.get('stability_differs'))
        n_tgp_stable = sum(1 for r in p3 if r.get('tgp_classification') == 'stable')
        n_newt_stable = sum(1 for r in p3 if r.get('newton_classification') == 'stable')
        n_newt_marginal = sum(1 for r in p3 if r.get('newton_classification') == 'marginal')
        print(f"\n  Bifurcation scan -- d_well ({n_valid} valid (beta,C) points):")
        print(f"    TGP stable: {n_tgp_stable}/{n_valid}")
        print(f"    Newton stable: {n_newt_stable}/{n_valid}")
        print(f"    Newton marginal (Earnshaw): {n_newt_marginal}/{n_valid}")
        print(f"    Classification differs: {n_differ}/{n_valid}")

    # From Part 4: N-gon
    if p4:
        for r in p4:
            if r.get('has_eq'):
                print(f"\n  {r['n']}-gon (R={r['R']:.4f}): "
                      f"Newton={r['newton_class']}, TGP={r['tgp_class']}, "
                      f"modes_stabilized={r['n_stabilized']}")

    # From Part 5: V3 effect
    if p5:
        print(f"\n  V3 (3-body) effect on stability:")
        for which, data in p5.items():
            print(f"    {which}: V2-only={data['v2_class']}, V2+V3={data['v3_class']}")

    print(f"\n  KEY FINDING:")
    # Determine overall conclusion
    if p3:
        if n_tgp_stable > n_newt_stable:
            print(f"    TGP produces MORE stable equilibria than Newton (P2 thesis SUPPORTED)")
        elif n_tgp_stable == n_newt_stable:
            print(f"    TGP and Newton have EQUAL stability count (P2 thesis INCONCLUSIVE)")
        else:
            print(f"    TGP produces FEWER stable equilibria than Newton (unexpected)")
    else:
        print(f"    Insufficient data for conclusion")


# ── Main ──────────────────────────────────────────────────────────────

def main():
    print("=" * 78)
    print("ex201 -- P2: Full Hessian Stability Mapping (TGP vs Newton)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")
    print(f"  Betas: {BETAS}")
    print(f"  C values: {CS}")
    print(f"  N-gons: {NGONS}")
    print(f"  G_Newton = 4*pi = {G_NEWTON:.4f}")

    t_total = time.time()

    p1 = part1_detailed_analysis()
    p2 = part2_beta_scan()
    p3 = part3_bifurcation_scan()
    p4 = part4_ngon_stability()
    p5 = part5_v3_effect()

    print_summary(p1, p2, p3, p4, p5)

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
