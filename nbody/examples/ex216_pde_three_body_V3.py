#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex216 -- Phase 2: Three-body force from 3D PDE (Layer 5 verification)
======================================================================

Extract the irreducible 3-body FORCE by comparing field gradients,
avoiding the energy-functional mismatch that plagues energy differences.

Method:
  1. Solve the 3D PDE with 3 sources (equilateral triangle, side d)
  2. Solve the 3D PDE with each pair (0,1) and (0,2)
  3. Force on particle 0: F = 12*pi*C * grad(u) at r_0
  4. 3-body force: F3 = F_total(123) - F_pair(01) - F_pair(02)
  5. Compare with exact Feynman 3-body force

Why this works:
  - The field u(x) is uniquely determined by the PDE (from the ODE action)
  - grad(u) at a source position is action-independent
  - Self-field cancels in the subtraction (same Gaussian at same position)
  - No energy functional is evaluated -> no action mismatch

Pass criteria:
  - Sign agreement (radial component)
  - Magnitude ratio within factor of 5 (limited by grid resolution)
  - All PDE solves converge

Depends on: tgp_pde_solver, three_body_force_exact, tgp_field
"""
import sys, os

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np
import time
from nbody.tgp_pde_solver import (
    TGPGrid, solve_field, extract_3body_force,
)
from nbody.three_body_force_exact import three_body_force_triplet_exact
from nbody.tgp_field import screening_mass


def equilateral_triangle(d: float) -> np.ndarray:
    """Vertices of equilateral triangle with side d, centered at origin."""
    R = d / np.sqrt(3.0)
    angles = np.array([np.pi / 2, np.pi / 2 + 2 * np.pi / 3,
                       np.pi / 2 + 4 * np.pi / 3])
    pos = np.zeros((3, 3))
    pos[:, 0] = R * np.cos(angles)
    pos[:, 1] = R * np.sin(angles)
    return pos


def main() -> None:
    print("=" * 70)
    print("ex216: Phase 2 -- Three-body FORCE from 3D PDE (Layer 5)")
    print("=" * 70)

    # -- Parameters --------------------------------------------------------
    C = 0.3
    beta = 1.0
    gamma = beta
    m_sp = screening_mass(beta, gamma)

    N = 64
    L = 20.0
    grid = TGPGrid(Nx=N, Lx=L)
    sigma = 1.5 * grid.dx

    print(f"  C = {C}, beta = {beta}, gamma = {gamma}, m_sp = {m_sp:.4f}")
    print(f"  Grid: {N}^3, L = {L}, dx = {grid.dx:.4f}")
    print()

    picard_kw = dict(max_iter=300, tol=1e-7, alpha=0.5)

    # -- Scan over separations ---------------------------------------------
    # Focus on d ~ 3-4 where F3 ~ exp(-2md) is well above grid noise.
    # At d >= 5 the 3-body signal drops below systematic errors.
    d_values = [3.0, 4.0, 5.0]
    d_max_safe = L / 2 - 3.0 / m_sp
    d_values = [d for d in d_values if d / np.sqrt(3.0) < d_max_safe]

    C_values = np.array([C, C, C])

    print(f"  {'d':>5s}  {'F3r_PDE':>12s}  {'F3r_Feyn':>12s}  "
          f"{'ratio':>8s}  {'|F3|/|F2|':>10s}  {'time':>6s}")
    print(f"  {'-'*5}  {'-'*12}  {'-'*12}  "
          f"{'-'*8}  {'-'*10}  {'-'*6}")

    results = []

    for d in d_values:
        positions = equilateral_triangle(d)

        # Radial unit vector: particle 0 toward center (= -pos[0]/|pos[0]|)
        r0 = positions[0]
        R0 = np.linalg.norm(r0)
        if R0 > 1e-10:
            r_hat = -r0 / R0   # inward (toward center)
        else:
            r_hat = np.array([0.0, -1.0, 0.0])

        t0 = time.time()

        # Extract 3-body force
        res = extract_3body_force(grid, positions, C_values, beta, gamma,
                                  sigma=sigma, picard_kw=picard_kw,
                                  verbose=False)
        dt = time.time() - t0

        F3_pde = res['F_3body']
        F3r_pde = np.dot(F3_pde, r_hat)  # radial component (positive = inward)

        # Total 2-body force magnitude for scale
        F2_total = res['F_pair_01'] + res['F_pair_02']
        # The pair forces include self-field, so F_total - F2 = F3
        # For scale, use |F_pair_01 + F_pair_02| as the 2-body force magnitude
        F2_mag = np.linalg.norm(F2_total)

        # Feynman reference
        F1_feyn, F2_feyn, F3_feyn = three_body_force_triplet_exact(
            positions[0], positions[1], positions[2],
            C, C, C, beta=beta, gamma=gamma, n_quad=60,
        )
        F3r_feyn = np.dot(F1_feyn, r_hat)

        if abs(F3r_feyn) > 1e-30:
            ratio = F3r_pde / F3r_feyn
        else:
            ratio = float('inf')

        if F2_mag > 1e-30:
            f3_over_f2 = np.linalg.norm(F3_pde) / F2_mag
        else:
            f3_over_f2 = 0.0

        conv = "OK" if res['all_converged'] else "FAIL"
        print(f"  {d:5.1f}  {F3r_pde:12.4e}  {F3r_feyn:12.4e}  "
              f"{ratio:8.3f}  {f3_over_f2:10.4f}  {dt:5.1f}s {conv}")

        results.append({
            'd': d,
            'F3r_pde': F3r_pde,
            'F3r_feyn': F3r_feyn,
            'ratio': ratio,
            'f3_over_f2': f3_over_f2,
            'converged': res['all_converged'],
            'F3_pde_vec': F3_pde,
            'F1_feyn_vec': F1_feyn,
        })

    # -- Check exponential decay -------------------------------------------
    print()
    if len(results) >= 2:
        d1, d2 = results[0]['d'], results[-1]['d']
        f1 = abs(results[0]['F3r_pde'])
        f2 = abs(results[-1]['F3r_pde'])
        if f1 > 1e-30 and f2 > 1e-30:
            m_eff = -np.log(f2 / f1) / (d2 - d1)
            # Feynman V3 ~ exp(-2md) for equilateral triangle
            m_expected = 2.0 * m_sp
            print(f"  F3 decay rate: m_eff = {m_eff:.3f} "
                  f"(expected ~ 2*m_sp = {m_expected:.3f} for 3-body)")
            dm_err = abs(m_eff - m_expected) / m_expected
            print(f"  Decay rate error: {dm_err:.1%}")

        # Same for Feynman
        ff1 = abs(results[0]['F3r_feyn'])
        ff2 = abs(results[-1]['F3r_feyn'])
        if ff1 > 1e-30 and ff2 > 1e-30:
            m_feyn = -np.log(ff2 / ff1) / (d2 - d1)
            print(f"  Feynman decay rate: {m_feyn:.3f}")

    # -- Summary -----------------------------------------------------------
    print()
    print("  ANALYSIS:")

    if not results:
        print("  No data points computed!")
        raise SystemExit(1)

    all_converged = all(r['converged'] for r in results)

    # Only check d-values where Feynman force is above noise threshold
    # F3 ~ exp(-2*m*d), grid noise ~ dx^2 * C^2 ~ 1e-3
    noise_floor = 1e-4
    strong_results = [r for r in results if abs(r['F3r_feyn']) > noise_floor]

    if not strong_results:
        print("  No d-values with sufficient signal above noise floor!")
        raise SystemExit(1)

    print(f"    d-values with F3_Feyn > {noise_floor:.0e}: "
          f"{[r['d'] for r in strong_results]}")

    # Sign agreement (both should be negative = attractive = inward)
    signs_agree = all(
        np.sign(r['F3r_pde']) == np.sign(r['F3r_feyn'])
        for r in strong_results
    )

    # Magnitude ratio (only for strong-signal d-values)
    ratios = [r['ratio'] for r in strong_results]
    if ratios:
        ok_ratio = all(0.2 < abs(r) < 5.0 for r in ratios)
        mean_ratio = np.mean(ratios)
    else:
        ok_ratio = False
        mean_ratio = float('inf')

    print(f"    Sign agreement: {'YES' if signs_agree else 'NO'}")
    print(f"    Ratios F3_PDE/F3_Feyn: [{', '.join(f'{r:.2f}' for r in ratios)}]")
    print(f"    Mean ratio: {mean_ratio:.3f}")
    print(f"    Ratios in [0.2, 5.0]: {'YES' if ok_ratio else 'NO'}")
    print(f"    All converged: {'YES' if all_converged else 'NO'}")
    print()

    # Key result
    print("  KEY RESULT:")
    if signs_agree and ok_ratio:
        print("    The 3D PDE solver independently reproduces the")
        print("    irreducible 3-body FORCE from the Feynman integral.")
        print("    Sign and magnitude agree -- Layer 5 VALIDATED.")
    elif signs_agree:
        print("    Signs agree; magnitude needs refinement (increase N).")
    else:
        print("    Sign or magnitude disagreement detected.")

    # Direction check
    print()
    print("  DIRECTION CHECK (3-body force vector at d={:.1f}):".format(
        results[0]['d']))
    F_pde = results[0]['F3_pde_vec']
    F_feyn = results[0]['F1_feyn_vec']
    cos_angle = (np.dot(F_pde, F_feyn)
                 / (np.linalg.norm(F_pde) * np.linalg.norm(F_feyn) + 1e-30))
    print(f"    F3_PDE  = [{F_pde[0]:.4e}, {F_pde[1]:.4e}, {F_pde[2]:.4e}]")
    print(f"    F3_Feyn = [{F_feyn[0]:.4e}, {F_feyn[1]:.4e}, {F_feyn[2]:.4e}]")
    print(f"    cos(angle) = {cos_angle:.6f}")

    # -- PASS/FAIL ---------------------------------------------------------
    ok = signs_agree and ok_ratio and all_converged

    print()
    print(f"  Overall: {'PASS' if ok else 'FAIL'}")

    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
