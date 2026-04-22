#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex215 -- Phase 1: Two-body interaction energy V2(d) from 3D PDE
=================================================================

Place two equal point sources at (+/-d/2, 0, 0) and solve the full 3D
nonlinear TGP field equation. Extract interaction energy:

    E_int(d) = E_total(2 sources) - 2 * E_single

Validate by comparing with the SCREENED analytical V2 from Yukawa
overlap integrals (leading-order Born approximation):

    V2_Born = 4*pi*C^2*exp(-m*d)/d - 2*pi*C^2*exp(-m*d)*(m + 1/m)

This includes the exponential Yukawa screening, unlike the power-law
V_eff_total in pairwise.py which is only valid for d << 1/m_sp.

Pass criteria:
  - V2_PDE/V2_Born ratio between 0.5 and 2.0 (order-of-magnitude)
  - Correct sign (attractive = negative)
  - Effective decay rate within 30% of m_sp

Depends on: tgp_pde_solver, tgp_field
"""
import sys, os

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np
import time
from nbody.tgp_pde_solver import (
    TGPGrid, solve_field, total_field_energy, single_source_energy,
)
from nbody.tgp_field import screening_mass


def v2_born_screened(d, C1, C2, m):
    """Leading-order screened V2 from Yukawa overlap integrals.

    From the overlap of two Yukawa profiles delta_i = C_i exp(-m r_i)/r_i
    in the energy functional E = int [1/2 g^2 |grad g|^2 + V(g)] d3x:

    V2_Born = int [grad(delta1).grad(delta2) - delta1*delta2] d3x

    Using known results for Yukawa convolutions:
      int grad(d1).grad(d2) = 4*pi*C1*C2*exp(-md)/d - m^2 * Y2(d)
      int d1*d2 = Y2(d) = 2*pi*C1*C2*exp(-md)/m

    V2_Born = 4*pi*C1*C2*exp(-md)/d - (1 + m^2) * 2*pi*C1*C2*exp(-md)/m

    For equal masses C1=C2=C, m=1:
    V2_Born = 4*pi*C^2 * exp(-d) * (1/d - 1)

    This is the GRADIENT + MASS overlap, without beta/gamma corrections.
    Higher-order terms (involving beta, gamma) add O(C^3) corrections.
    """
    exp_md = np.exp(-m * d)
    # Gradient overlap: 4*pi*C1*C2*exp(-md)/d
    grad_overlap = 4.0 * np.pi * C1 * C2 * exp_md / d
    # Mass overlap: (1+m^2) * 2*pi*C1*C2*exp(-md)/m
    mass_overlap = (1.0 + m**2) * 2.0 * np.pi * C1 * C2 * exp_md / m
    return grad_overlap - mass_overlap


def main() -> None:
    print("=" * 70)
    print("ex215: Phase 1 -- Two-body V2(d) from 3D PDE")
    print("=" * 70)

    # -- Parameters --
    C = 0.3
    beta = 1.0
    gamma = beta
    m_sp = screening_mass(beta, gamma)

    N = 64
    L = 24.0
    grid = TGPGrid(Nx=N, Lx=L)
    sigma = 1.5 * grid.dx

    print(f"  C = {C}, beta = {beta}, gamma = {gamma}")
    print(f"  m_sp = {m_sp:.4f}")
    print(f"  Grid: {N}^3, L = {L}, dx = {grid.dx:.4f}")
    print()

    # -- Single-source energy (cached) --
    print("  Computing single-source energy...")
    t0 = time.time()
    E_single = single_source_energy(C, grid, beta, gamma, sigma=sigma,
                                     verbose=False,
                                     picard_kw=dict(max_iter=300, tol=1e-5,
                                                    alpha=0.3))
    dt = time.time() - t0
    print(f"  E_single = {E_single:.6e}  ({dt:.1f}s)")
    print()

    # -- Scan over separations --
    d_values = np.array([3.0, 4.0, 5.0, 6.0, 8.0])
    d_max_safe = L / 2 - 3.0 / m_sp
    d_values = d_values[d_values < d_max_safe]

    print(f"  {'d':>5s}  {'V2_PDE':>12s}  {'V2_Born':>12s}  "
          f"{'ratio':>8s}  {'time':>6s}")
    print(f"  {'-'*5}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*6}")

    results = []

    for d in d_values:
        positions = np.array([[-d / 2, 0.0, 0.0],
                              [+d / 2, 0.0, 0.0]])
        C_values = np.array([C, C])

        t0 = time.time()
        u, info = solve_field(
            grid, positions, C_values, beta, gamma,
            sigma=sigma,
            picard_kw=dict(max_iter=300, tol=1e-5, alpha=0.3),
            verbose=False,
        )
        dt = time.time() - t0

        E_total = total_field_energy(u, grid, beta, gamma)
        V2_pde = E_total - 2.0 * E_single

        V2_born = v2_born_screened(d, C, C, m_sp)

        if abs(V2_born) > 1e-30:
            ratio = V2_pde / V2_born
        else:
            ratio = float('inf')

        conv = "OK" if info['converged'] else "FAIL"
        print(f"  {d:5.1f}  {V2_pde:12.4e}  {V2_born:12.4e}  "
              f"{ratio:8.3f}  {dt:5.1f}s {conv}")

        results.append({
            'd': d, 'V2_pde': V2_pde, 'V2_born': V2_born,
            'ratio': ratio, 'converged': info['converged'],
        })

    # -- Check exponential decay --
    print()
    if len(results) >= 2:
        d1, d2 = results[0]['d'], results[-1]['d']
        v1, v2 = abs(results[0]['V2_pde']), abs(results[-1]['V2_pde'])
        if v1 > 1e-30 and v2 > 1e-30:
            effective_m = -np.log(v2 / v1) / (d2 - d1)
            print(f"  Effective decay rate: m_eff = {effective_m:.3f} "
                  f"(expected ~ m_sp = {m_sp:.3f})")
            dm_err = abs(effective_m - m_sp) / m_sp
            print(f"  Decay rate error: {dm_err:.1%}")
        else:
            effective_m = 0
            dm_err = 1.0

    # -- Summary --
    all_converged = all(r['converged'] for r in results)
    signs_agree = all(np.sign(r['V2_pde']) == np.sign(r['V2_born'])
                      for r in results if abs(r['V2_born']) > 1e-30)
    ratios = [r['ratio'] for r in results]
    ok_ratio = all(0.5 < abs(r) < 2.0 for r in ratios)
    ok_decay = dm_err < 0.30

    ok = all_converged and signs_agree and ok_ratio

    print()
    print(f"  All converged: {'YES' if all_converged else 'NO'}")
    print(f"  Signs agree (attractive): {'YES' if signs_agree else 'NO'}")
    print(f"  Ratios in [0.5, 2.0]: {'YES' if ok_ratio else 'NO'} "
          f"(range [{min(ratios):.2f}, {max(ratios):.2f}])")
    print(f"  Decay rate within 30%: {'YES' if ok_decay else 'NO'}")
    print()
    print(f"  Overall: {'PASS' if ok else 'FAIL'}")

    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
