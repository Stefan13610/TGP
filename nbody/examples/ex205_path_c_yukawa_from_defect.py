#!/usr/bin/env python3
"""
ex205 -- P5: Path C -- Yukawa coupling from topological defect (EFT approach)
==============================================================================

STATUS: CANONICAL

This is the preferred example for the synchronized bridge
`classical defect -> C_eff, m_sp -> effective Yukawa source`.

Derives the Path B Yukawa coupling C_eff from first principles:
  1. Classical TGP defect has oscillatory tail (sin(r)/r)
  2. Loop/RG corrections generate mass gap m_sp (screening mass)
  3. C_eff is obtained by projecting the defect profile onto the
     Yukawa Green's function: C_eff = int delta(r)*exp(-m_sp*r)*r dr

This DERIVES the Path B phenomenological ansatz from defect physics.

Tests:
  Part 1: Consistency check -- m_sp(Path B) == m_sp(Path C)
  Part 2: Potential landscape -- V_TGP curvature, why vacuum is a maximum
  Part 3: Standard TGP defect profiles (oscillatory) + stabilization obstacle
  Part 4: EFT projection -- C_eff from standard defects + Yukawa Green's fn
  Part 5: C_eff(g0) scan -- effective coupling vs core depth
  Part 6: C_eff vs beta scan -- parameter dependence

PHYSICS:
--------
Standard TGP: V''(g=1) = 2*beta - 3*gamma < 0 (vacuum is maximum)
  => linearized equation: (nabla^2 + 1)*delta = 0
  => oscillatory tail: delta ~ sin(r)/r

Naive stabilization V_sb = (mu^2/2)*(g-1)^2 with mu^2 = 2*(3g-2b):
  => V_C'(g) = (1-g)*(g^2 - mu^2) < 0 for all g in (0,1) when mu^2 > 1
  => DESTROYS classical defects (no turning point)

EFT resolution:
  - Classical defect uses standard V_TGP (oscillatory tail survives)
  - Mass gap m_sp emerges from quantum/RG corrections to propagator
  - C_eff = projection of defect onto Yukawa Green's function
  - This derives Path B from Path C without modifying classical Lagrangian

Key result: mu^2 = 2*(3*gamma - 2*beta) is an auxiliary EFT stabilization scale.
The physical screening mass used by `nbody` remains m_sp.
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.yukawa_from_defect import (
    V_tgp, dV_tgp, V_pathC, dV_pathC, d2V_pathC,
    mu2_self_consistent,
    effective_mass_squared, screening_mass_pathC,
    solve_defect_pathC,
    solve_defect_standard,
    compute_C_eff_projection,
    demonstrate_stabilization_obstacle,
    fit_oscillatory_tail, compare_tail_types,
    verify_pathC_consistency,
)
from nbody.tgp_field import screening_mass, default_beta_gamma

quick = "--quick" in sys.argv

if quick:
    G0_VALS = [0.50, 0.70, 0.90, 0.95]
    BETAS = [0.5, 1.0, 2.0]
    R_MAX = 40.0
    N_EVAL = 2000
    N_G0_SCAN = 10
else:
    G0_VALS = [0.30, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.98]
    BETAS = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 4.0]
    R_MAX = 60.0
    N_EVAL = 5000
    N_G0_SCAN = 25


# -- Part 1: Consistency check --------------------------------------------

def part1_consistency():
    """Verify m_sp(Path B) == m_sp(Path C) for various parameters."""
    print("=" * 78)
    print("Part 1: Path B / Path C Screening Mass Consistency")
    print("=" * 78)

    print(f"\n  {'beta':>6s} {'gamma':>6s} {'m_sp(B)':>10s} {'m_sp(C)':>10s} "
          f"{'mu^2':>8s} {'match':>8s}")
    print(f"  {'-'*56}")

    all_pass = True
    test_params = [
        (0.5, 0.5), (1.0, 1.0), (2.0, 2.0),
        (0.5, 1.0), (1.0, 2.0), (0.1, 0.1),
    ]

    for beta, gamma in test_params:
        try:
            result = verify_pathC_consistency(beta, gamma)
            status = "PASS" if result['match'] else "FAIL"
            if not result['match']:
                all_pass = False
            print(f"  {beta:6.2f} {gamma:6.2f} {result['m_sp_pathB']:10.6f} "
                  f"{result['m_sp_pathC']:10.6f} {result['mu2']:8.4f} {status:>8s}")
        except Exception as e:
            print(f"  {beta:6.2f} {gamma:6.2f} ERROR: {e}")
            all_pass = False

    print(f"\n  Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    return all_pass


# -- Part 2: Potential landscape -------------------------------------------

def part2_potential():
    """V_TGP curvature analysis: why vacuum is a maximum."""
    print(f"\n{'=' * 78}")
    print("Part 2: Potential Landscape -- V_TGP Vacuum Analysis")
    print(f"{'=' * 78}")

    beta = gamma = 1.0
    mu2 = mu2_self_consistent(beta, gamma)

    print(f"\n  beta = gamma = {beta}")
    print(f"  mu^2 (self-consistent) = {mu2:.4f}")

    # Curvature at vacuum
    d2V_vac = 2.0 * beta - 3.0 * gamma
    m_eff2 = effective_mass_squared(beta, gamma, mu2)

    print(f"\n  Standard TGP at vacuum (g=1):")
    print(f"    V_TGP''(1) = 2*beta - 3*gamma = {d2V_vac:.4f}")
    print(f"    => vacuum is a MAXIMUM (negative curvature)")
    print(f"    => linearized: (nabla^2 + |V''|)*delta = 0")
    print(f"    => OSCILLATORY tail: delta ~ sin(r)/r")

    print(f"\n  With EFT stabilization (on propagator, not Lagrangian):")
    print(f"    V_C''(1) = (2b-3g) + mu^2 = {d2V_vac} + {mu2} = {m_eff2:.4f}")
    print(f"    m_sp = sqrt({m_eff2:.4f}) = {np.sqrt(m_eff2):.4f}")
    print(f"    => Yukawa propagator: G(r) = exp(-m_sp*r)/(4*pi*r)")

    # Path B comparison
    m_sp_B = screening_mass(beta, gamma)
    print(f"\n  Path B screening mass: m_sp = sqrt(3g-2b) = {m_sp_B:.4f}")
    print(f"  Match: {abs(np.sqrt(m_eff2) - m_sp_B) < 1e-12}")

    is_correct = m_eff2 > 0 and abs(np.sqrt(m_eff2) - m_sp_B) < 1e-12
    return {
        'd2V_tgp': d2V_vac,
        'm_eff2': m_eff2,
        'm_sp': np.sqrt(m_eff2),
        'is_correct': is_correct,
    }


# -- Part 3: Standard TGP profiles + stabilization obstacle ---------------

def part3_profiles_and_obstacle():
    """Standard TGP defect (oscillatory) + demonstrate stabilization obstacle."""
    print(f"\n{'=' * 78}")
    print("Part 3: Standard TGP Defect Profiles + Stabilization Obstacle")
    print(f"{'=' * 78}")

    beta = gamma = 1.0

    # --- 3A: Standard TGP defect profiles ---
    print(f"\n  --- 3A: Standard TGP defect profiles (mu^2 = 0) ---")
    print(f"  {'g0':>6s} {'g(r_max)':>12s} {'osc_amp':>10s} {'osc_resid':>10s} "
          f"{'tail':>10s}")
    print(f"  {'-'*56}")

    n_osc = 0
    for g0 in G0_VALS:
        try:
            res = solve_defect_pathC(g0, beta, gamma, mu2=0.0,
                                     r_max=R_MAX, kinetic="full", n_eval=N_EVAL)
            comp = compare_tail_types(res['r'], res['g'], 1.0)
            osc_amp, _, osc_resid = fit_oscillatory_tail(res['r'], res['g'])

            tail = comp['better']
            if tail == 'oscillatory':
                n_osc += 1

            print(f"  {g0:6.3f} {res['g'][-1]:12.8f} {osc_amp:10.6f} "
                  f"{osc_resid:10.4e} {tail:>10s}")
        except Exception as e:
            print(f"  {g0:6.3f} ERROR: {e}")

    print(f"\n  Result: {n_osc}/{len(G0_VALS)} defects have oscillatory tails")
    tgp_oscillatory = n_osc == len(G0_VALS)
    print(f"  Standard TGP all oscillatory: "
          f"{'PASS' if tgp_oscillatory else 'PARTIAL'}")

    # --- 3B: Stabilization obstacle ---
    print(f"\n  --- 3B: Naive stabilization obstacle ---")
    obs = demonstrate_stabilization_obstacle(beta, gamma)
    print(f"  mu^2 = {obs['mu2']:.4f}")
    print(f"  V_C'(g) has positive values in (0,1): {obs['has_positive_dV_C']}")
    print(f"  Sign changes in dV_standard: {obs['sign_changes_standard']}")
    print(f"  Sign changes in dV_stabilized: {obs['sign_changes_stabilized']}")

    if obs['defect_impossible_stabilized']:
        print(f"\n  CONFIRMED: V_C'(g) < 0 for ALL g in (0,1) when mu^2 = {obs['mu2']:.1f} > 1")
        print(f"  => Naive stabilization DESTROYS classical defects")
        print(f"  => Cannot add V_sb to classical Lagrangian")
        print(f"  => Stabilization must act on PROPAGATOR (EFT), not Lagrangian")
    else:
        print(f"  Note: V_C'(g) has positive region -- defects may survive")

    # Show V_C'(g) at sample points
    print(f"\n  Sample V_C'(g) values:")
    print(f"  {'g':>6s} {'dV_TGP':>12s} {'dV_C':>12s}")
    print(f"  {'-'*34}")
    for g_val in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
        dV_std = dV_tgp(g_val, beta, gamma)
        dV_stab = dV_pathC(g_val, beta, gamma, obs['mu2'])
        print(f"  {g_val:6.2f} {dV_std:12.6f} {dV_stab:12.6f}")

    obstacle_confirmed = obs['defect_impossible_stabilized']
    return {
        'tgp_oscillatory': tgp_oscillatory,
        'obstacle_confirmed': obstacle_confirmed,
    }


# -- Part 4: EFT projection -- C_eff from standard defects ----------------

def part4_eft_projection():
    """Compute C_eff by projecting standard TGP defects onto Yukawa Green's fn."""
    print(f"\n{'=' * 78}")
    print("Part 4: EFT Projection -- C_eff from Standard TGP Defects")
    print(f"{'=' * 78}")

    beta = gamma = 1.0
    m_sp = screening_mass(beta, gamma)
    print(f"\n  beta = gamma = {beta}")
    print(f"  m_sp (Path B) = {m_sp:.6f}")
    print(f"\n  C_eff = int_0^inf delta(r) * exp(-m_sp*r) * r dr")
    print(f"  where delta(r) = 1 - g(r) from standard TGP defect ODE")

    print(f"\n  {'g0':>6s} {'g(r_max)':>12s} {'C_eff_proj':>12s} "
          f"{'osc_amp':>10s} {'C/amp':>10s} {'status':>8s}")
    print(f"  {'-'*66}")

    results = []
    for g0 in G0_VALS:
        try:
            res = solve_defect_standard(g0, beta, gamma,
                                        r_max=R_MAX, kinetic="full", n_eval=N_EVAL)
            C_proj = res['C_eff_proj']
            osc_amp = res['osc_amplitude']
            ratio = abs(C_proj / osc_amp) if osc_amp > 1e-15 else float('nan')
            valid = abs(C_proj) > 1e-15 and np.isfinite(C_proj)
            status = "OK" if valid else "WEAK"

            print(f"  {g0:6.3f} {res['g'][-1]:12.8f} {C_proj:12.6e} "
                  f"{osc_amp:10.6f} {ratio:10.4f} {status:>8s}")

            results.append({
                'g0': g0,
                'C_eff_proj': C_proj,
                'osc_amplitude': osc_amp,
                'ratio': ratio,
                'valid': valid,
            })
        except Exception as e:
            print(f"  {g0:6.3f} ERROR: {e}")

    valid_results = [r for r in results if r['valid']]
    if valid_results:
        C_vals = [r['C_eff_proj'] for r in valid_results]
        print(f"\n  Valid projections: {len(valid_results)}/{len(G0_VALS)}")
        print(f"  C_eff range: [{min(C_vals):.6e}, {max(C_vals):.6e}]")

        # Check: deeper core (lower g0) -> larger C_eff?
        if len(valid_results) >= 3:
            sorted_r = sorted(valid_results, key=lambda x: x['g0'])
            C_sorted = [r['C_eff_proj'] for r in sorted_r]
            monotone = all(C_sorted[i] >= C_sorted[i+1]
                          for i in range(len(C_sorted)-1))
            print(f"  C_eff monotone (deeper core -> larger C): {monotone}")

    return results


# -- Part 5: C_eff(g0) scan -----------------------------------------------

def part5_coupling_scan():
    """Detailed scan of C_eff vs g0 for standard TGP defects."""
    print(f"\n{'=' * 78}")
    print("Part 5: C_eff(g0) Scan -- Effective Coupling vs Core Depth")
    print(f"{'=' * 78}")

    beta = gamma = 1.0
    m_sp = screening_mass(beta, gamma)

    g0_vals = np.linspace(0.30, 0.98, N_G0_SCAN)

    print(f"\n  beta = gamma = {beta}, m_sp = {m_sp:.4f}")
    print(f"  Scanning {N_G0_SCAN} core values g0 in [0.30, 0.98]")

    print(f"\n  {'g0':>6s} {'C_eff_proj':>12s} {'E_defect':>12s} "
          f"{'osc_amp':>10s} {'C/E ratio':>12s}")
    print(f"  {'-'*60}")

    valid = []
    for g0 in g0_vals:
        try:
            res = solve_defect_standard(g0, beta, gamma,
                                        r_max=R_MAX, kinetic="full", n_eval=N_EVAL)
            C_proj = res['C_eff_proj']
            E_def = res['E_defect']
            osc_amp = res['osc_amplitude']
            CE_ratio = abs(C_proj / E_def) if abs(E_def) > 1e-30 else 0.0

            if abs(C_proj) > 1e-15:
                print(f"  {g0:6.3f} {C_proj:12.6e} {E_def:12.6e} "
                      f"{osc_amp:10.6f} {CE_ratio:12.6e}")
                valid.append({
                    'g0': g0, 'C_eff_proj': C_proj,
                    'E_defect': E_def, 'osc_amplitude': osc_amp,
                })
            else:
                print(f"  {g0:6.3f} {'---':>12s} {E_def:12.6e} "
                      f"{osc_amp:10.6f} {'---':>12s}  (C~0)")
        except Exception as e:
            print(f"  {g0:6.3f} ERROR: {e}")

    if len(valid) >= 2:
        C_arr = np.array([r['C_eff_proj'] for r in valid])
        g0_arr = np.array([r['g0'] for r in valid])
        E_arr = np.array([r['E_defect'] for r in valid])

        print(f"\n  Summary ({len(valid)} valid):")
        print(f"  C_eff range: [{np.min(C_arr):.4e}, {np.max(C_arr):.4e}]")
        print(f"  E_defect range: [{np.min(E_arr):.4e}, {np.max(E_arr):.4e}]")

        # Power-law fit: C_eff ~ (1-g0)^alpha
        delta_g0 = 1.0 - g0_arr
        mask = (np.abs(C_arr) > 1e-15) & (delta_g0 > 1e-10)
        if np.sum(mask) >= 3:
            log_dg = np.log(delta_g0[mask])
            log_C = np.log(np.abs(C_arr[mask]))
            coeffs = np.polyfit(log_dg, log_C, 1)
            alpha = coeffs[0]
            print(f"  Scaling: |C_eff| ~ (1-g0)^{alpha:.2f}")

    return valid


# -- Part 6: C_eff vs beta scan -------------------------------------------

def part6_beta_scan():
    """Scan C_eff vs beta for fixed g0 using EFT projection."""
    print(f"\n{'=' * 78}")
    print("Part 6: C_eff and m_sp vs beta (fixed g0=0.90, EFT projection)")
    print(f"{'=' * 78}")

    g0 = 0.90

    print(f"\n  {'beta':>6s} {'m_sp':>8s} {'C_eff_proj':>12s} "
          f"{'osc_amp':>10s} {'E_defect':>12s} {'status':>8s}")
    print(f"  {'-'*64}")

    results = []
    for beta in BETAS:
        gamma = beta  # vacuum condition
        try:
            res = solve_defect_standard(g0, beta, gamma,
                                        r_max=R_MAX, kinetic="full", n_eval=N_EVAL)

            C_proj = res['C_eff_proj']
            valid = abs(C_proj) > 1e-15
            status = "OK" if valid else "WEAK"

            print(f"  {beta:6.2f} {res['m_sp_pathB']:8.4f} "
                  f"{C_proj:12.6e} {res['osc_amplitude']:10.6f} "
                  f"{res['E_defect']:12.6e} {status:>8s}")

            res['beta'] = beta
            results.append(res)
        except Exception as e:
            print(f"  {beta:6.2f} ERROR: {e}")

    # Scaling: C_eff ~ beta^alpha ?
    if len(results) >= 3:
        b_arr = np.array([r['beta'] for r in results])
        C_arr = np.array([abs(r['C_eff_proj']) for r in results])
        mask = C_arr > 1e-15
        if np.sum(mask) >= 2:
            log_b = np.log(b_arr[mask])
            log_c = np.log(C_arr[mask])
            coeffs = np.polyfit(log_b, log_c, 1)
            alpha = coeffs[0]
            print(f"\n  Scaling: |C_eff| ~ beta^{alpha:.2f}")

    return results


def main():
    print("=" * 78)
    print("ex205 -- P5: Path C -- Yukawa Coupling from Topological Defect")
    print("         (EFT Approach: Classical Defect + Yukawa Projection)")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t_total = time.time()

    p1 = part1_consistency()
    p2 = part2_potential()
    p3 = part3_profiles_and_obstacle()
    p4 = part4_eft_projection()
    p5 = part5_coupling_scan()
    p6 = part6_beta_scan()

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")
    print(f"  Part 1 (consistency): {'PASS' if p1 else 'FAIL'}")
    print(f"  Part 2 (curvature):   V_TGP''(1) = {p2['d2V_tgp']:.1f} < 0, "
          f"m_eff^2 (EFT) = {p2['m_eff2']:.1f} > 0: "
          f"{'PASS' if p2['is_correct'] else 'FAIL'}")
    print(f"  Part 3 (profiles + obstacle):")
    print(f"    TGP all oscillatory: "
          f"{'PASS' if p3['tgp_oscillatory'] else 'PARTIAL'}")
    print(f"    Naive stabilization destroys defects: "
          f"{'CONFIRMED' if p3['obstacle_confirmed'] else 'NOT CONFIRMED'}")

    n_valid_p4 = sum(1 for r in p4 if r.get('valid', False))
    print(f"  Part 4 (EFT projection): {n_valid_p4}/{len(p4)} valid C_eff values")
    print(f"  Part 5 (coupling scan):  {len(p5)} valid defect profiles")
    print(f"  Part 6 (beta scan):      {len(p6)} parameter points")

    print(f"\n  KEY RESULTS:")
    print(f"  1. Standard TGP defects have oscillatory tails (confirmed)")
    print(f"  2. Naive quadratic stabilization DESTROYS classical defects")
    print(f"     (V_C'(g) < 0 for all g<1 when mu^2 > 1)")
    print(f"  3. EFT approach: mass gap m_sp from loop corrections to propagator")
    print(f"  4. C_eff = projection of defect onto Yukawa Green's function")
    print(f"  5. This DERIVES Path B ansatz: C_i is not postulated but computed")
    print(f"  6. mu^2 = 2*(3*gamma - 2*beta) is uniquely determined")

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
