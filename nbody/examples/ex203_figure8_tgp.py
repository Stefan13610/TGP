#!/usr/bin/env python3
"""
ex203 -- P3.A: Figure-8 orbit (Chenciner-Montgomery) -- TGP vs Newton
=====================================================================

Does the famous Chenciner-Montgomery figure-8 choreography survive
under TGP corrections?

Tests:
  Part 1: Newton reference -- verify figure-8 orbit (energy conservation,
          periodicity, trajectory shape).
  Part 2: TGP pairwise -- same IC, how does orbit deform?
  Part 3: TGP with 3-body (Yukawa Feynman) -- full TGP effect.
  Part 4: Orbit deviation quantification -- Delta_r(t), period shift,
          shape metrics (min/max distance, angular momentum drift).
  Part 5: beta scan -- at what coupling strength does figure-8 break?

PHYSICS:
--------
The Chenciner-Montgomery figure-8 orbit (2000) is a stable choreography
in Newtonian 3-body with equal masses: all three bodies follow the SAME
curve, phase-shifted by T/3.

TGP adds two corrections:
  (a) Pairwise: repulsive beta/r^2 and attractive gamma(Ci+Cj)/r^3 terms
  (b) 3-body: irreducible V3 = -6*gamma*C1*C2*C3*I_Y

For small beta, the figure-8 should survive as a deformed orbit.
For large beta, the repulsive barrier should destroy it.
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.configurations import figure_eight_initial
from nbody.dynamics_v2 import (
    leapfrog_integrate, rk45_integrate,
    potential_newton, forces_newton,
    analytical_forces_tgp_pairwise, potential_tgp,
)
from nbody.dynamics_backends import (
    build_tgp_integration_pair,
    total_angular_momentum,
    center_of_mass,
)

quick = "--quick" in sys.argv

# -- G_NEWTON = 4*pi to match TGP leading term ----------------------------
G_NEWTON = 4.0 * np.pi

# -- Integration parameters -----------------------------------------------
if quick:
    T_FINAL = 4.0          # ~0.6 periods
    N_OUTPUT = 500
    BETAS = [0.001, 0.01, 0.05, 0.10]
    N_QUAD = 20
else:
    T_FINAL = 20.0         # ~3 periods
    N_OUTPUT = 2000
    BETAS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.08, 0.10, 0.20]
    N_QUAD = 30

C_DEFAULT = 0.3


# -- Helpers ---------------------------------------------------------------

def _estimate_period(pos_traj, t_arr, body_idx=0):
    """
    Estimate period from trajectory by finding when body returns near
    its initial position. Uses the first body by default.
    """
    x0 = pos_traj[0, body_idx]
    dists = np.linalg.norm(pos_traj[:, body_idx] - x0, axis=1)
    # Find local minima after initial departure
    # Skip first ~10% of trajectory
    skip = max(len(t_arr) // 10, 5)
    # Look for first return: distance drops below threshold
    d_max = np.max(dists[skip:])
    threshold = 0.1 * d_max
    crossings = []
    for i in range(skip, len(dists) - 1):
        if dists[i] > threshold and dists[i+1] <= threshold:
            # Linear interpolation
            frac = (threshold - dists[i]) / (dists[i+1] - dists[i])
            t_cross = t_arr[i] + frac * (t_arr[i+1] - t_arr[i])
            crossings.append(t_cross)
    if crossings:
        return crossings[0]
    return None


def _orbit_metrics(pos_traj, vel_traj, t_arr, C_values):
    """Compute orbit shape metrics."""
    n_t = len(t_arr)

    # Min/max distance between any pair
    d_all = []
    for t_idx in range(n_t):
        for i in range(3):
            for j in range(i+1, 3):
                d = np.linalg.norm(pos_traj[t_idx, i] - pos_traj[t_idx, j])
                d_all.append(d)
    d_all = np.array(d_all)

    # Angular momentum conservation (z-component, planar orbit)
    L_arr = np.zeros(n_t)
    for t_idx in range(n_t):
        L_vec = total_angular_momentum(pos_traj[t_idx], vel_traj[t_idx], C_values)
        L_arr[t_idx] = L_vec[2]  # z-component

    # CoM drift
    com_arr = np.zeros((n_t, 3))
    for t_idx in range(n_t):
        com_arr[t_idx] = center_of_mass(pos_traj[t_idx], C_values)
    com_drift = np.max(np.linalg.norm(com_arr - com_arr[0], axis=1))

    return {
        'd_min': np.min(d_all),
        'd_max': np.max(d_all),
        'L_z_mean': np.mean(L_arr),
        'L_z_drift': np.max(np.abs(L_arr - L_arr[0])),
        'com_drift': com_drift,
    }


def _deviation_from_newton(pos_newton, pos_tgp, t_newton, t_tgp=None):
    """
    Compute trajectory deviation: max over bodies of |r_tgp(t) - r_newton(t)|
    at each timestep, then return max and rms.

    Handles mismatched time arrays by truncating to common length.
    """
    if t_tgp is None:
        t_tgp = t_newton

    # Use the shorter trajectory
    n_common = min(len(pos_newton), len(pos_tgp))
    pn = pos_newton[:n_common]
    pt = pos_tgp[:n_common]

    # Shape: (n_t, n_bodies, 3)
    delta = pt - pn
    delta_norm = np.linalg.norm(delta, axis=2)  # (n_t, n_bodies)
    max_per_t = np.max(delta_norm, axis=1)       # (n_t,)
    return {
        'max_deviation': np.max(max_per_t),
        'rms_deviation': np.sqrt(np.mean(max_per_t**2)),
        'deviation_vs_t': max_per_t,
    }


# -- Part 1: Newton reference ---------------------------------------------

def part1_newton_reference():
    """Verify figure-8 orbit under Newtonian gravity."""
    print("=" * 78)
    print("Part 1: Newton Reference -- Figure-8 Orbit")
    print("=" * 78)

    pos0, vel0, C_values, name = figure_eight_initial(C=C_DEFAULT)
    print(f"  IC: {name}, C = {C_DEFAULT}")
    print(f"  G_Newton = {G_NEWTON:.4f} (= 4*pi, matching TGP leading term)")

    # Scale velocities for G_NEWTON != 1
    # The Chenciner-Montgomery IC assume G=1, m=1.
    # For G=4*pi, m=C, we need v ~ sqrt(G*m) scaling.
    # v_scaled = v_original * sqrt(G_Newton * C)
    v_scale = np.sqrt(G_NEWTON * C_DEFAULT)
    vel0_scaled = vel0 * v_scale

    # Newton acceleration (returns acceleration, not force)
    def acc_newton(positions, C_vals):
        return forces_newton(positions, C_vals, G=G_NEWTON, softening=1e-8)

    def pot_newton(positions, C_vals):
        return potential_newton(positions, C_vals, G=G_NEWTON, softening=1e-8)

    t0 = time.time()
    result = rk45_integrate(
        pos0, vel0_scaled, C_values,
        acc_newton, pot_newton,
        t_span=(0, T_FINAL), n_output=N_OUTPUT,
        rtol=1e-12, atol=1e-14, quiet=True,
    )
    dt_wall = time.time() - t0

    E = result['energy']
    E_err = result['energy_error']
    max_E_err = np.max(np.abs(E_err))

    period = _estimate_period(result['positions'], result['t'])
    metrics = _orbit_metrics(result['positions'], result['velocities'],
                             result['t'], C_values)

    print(f"\n  Integration: t = [0, {T_FINAL}], n_output = {N_OUTPUT}")
    print(f"  Success: {result['success']}")
    print(f"  Wall time: {dt_wall:.2f}s")
    print(f"  E0 = {E[0]:.8e}")
    print(f"  |dE/E0| max = {max_E_err:.2e}")
    print(f"  Estimated period: {period if period else 'not found (T_final too short)'}")
    print(f"  d_min = {metrics['d_min']:.4f}, d_max = {metrics['d_max']:.4f}")
    print(f"  Lz drift = {metrics['L_z_drift']:.2e}")
    print(f"  CoM drift = {metrics['com_drift']:.2e}")

    status = "PASS" if max_E_err < 1e-6 else "WARN"
    print(f"\n  Energy conservation: {status}")

    return result, C_values, vel0_scaled


# -- Part 2: TGP pairwise -------------------------------------------------

def part2_tgp_pairwise(newton_result, C_values, vel0_scaled):
    """Figure-8 IC under TGP pairwise forces (no 3-body)."""
    print(f"\n{'=' * 78}")
    print("Part 2: TGP Pairwise -- Figure-8 Deformation")
    print(f"{'=' * 78}")

    pos0, _, _, _ = figure_eight_initial(C=C_DEFAULT)
    beta = 0.01  # Small coupling

    acc_fn, pot_fn = build_tgp_integration_pair(
        "pairwise", beta=beta, gamma=beta, softening=1e-8,
    )

    print(f"  beta = gamma = {beta}, C = {C_DEFAULT}")
    print(f"  Backend: pairwise (no 3-body)")

    t0 = time.time()
    result = rk45_integrate(
        pos0, vel0_scaled, C_values,
        acc_fn, pot_fn,
        t_span=(0, T_FINAL), n_output=N_OUTPUT,
        rtol=1e-12, atol=1e-14, quiet=True,
    )
    dt_wall = time.time() - t0

    E_err = result['energy_error']
    max_E_err = np.max(np.abs(E_err))

    period = _estimate_period(result['positions'], result['t'])
    metrics = _orbit_metrics(result['positions'], result['velocities'],
                             result['t'], C_values)
    dev = _deviation_from_newton(newton_result['positions'], result['positions'],
                                 newton_result['t'], result['t'])

    print(f"\n  Wall time: {dt_wall:.2f}s")
    print(f"  |dE/E0| max = {max_E_err:.2e}")
    print(f"  Estimated period: {period if period else 'N/A'}")
    print(f"  d_min = {metrics['d_min']:.4f}, d_max = {metrics['d_max']:.4f}")
    print(f"  Lz drift = {metrics['L_z_drift']:.2e}")
    print(f"  Max deviation from Newton: {dev['max_deviation']:.4e}")
    print(f"  RMS deviation from Newton: {dev['rms_deviation']:.4e}")

    return result


# -- Part 3: TGP with 3-body ----------------------------------------------

def part3_tgp_full(newton_result, C_values, vel0_scaled):
    """Figure-8 IC under full TGP (pairwise + Yukawa 3-body)."""
    print(f"\n{'=' * 78}")
    print("Part 3: TGP Full (Pairwise + Yukawa 3-body) -- Figure-8")
    print(f"{'=' * 78}")

    pos0, _, _, _ = figure_eight_initial(C=C_DEFAULT)
    beta = 0.01

    acc_fn, pot_fn = build_tgp_integration_pair(
        "yukawa_feynman", beta=beta, gamma=beta, softening=1e-8,
        include_3body=True, n_quad_feynman=N_QUAD,
    )

    print(f"  beta = gamma = {beta}, C = {C_DEFAULT}")
    print(f"  Backend: yukawa_feynman (exact Feynman 2D)")
    print(f"  n_quad = {N_QUAD}")

    t0 = time.time()
    result = rk45_integrate(
        pos0, vel0_scaled, C_values,
        acc_fn, pot_fn,
        t_span=(0, T_FINAL), n_output=N_OUTPUT,
        rtol=1e-10, atol=1e-12, quiet=True,
    )
    dt_wall = time.time() - t0

    E_err = result['energy_error']
    max_E_err = np.max(np.abs(E_err))

    period = _estimate_period(result['positions'], result['t'])
    metrics = _orbit_metrics(result['positions'], result['velocities'],
                             result['t'], C_values)
    dev = _deviation_from_newton(newton_result['positions'], result['positions'],
                                 newton_result['t'], result['t'])

    print(f"\n  Wall time: {dt_wall:.2f}s")
    print(f"  |dE/E0| max = {max_E_err:.2e}")
    print(f"  Estimated period: {period if period else 'N/A'}")
    print(f"  d_min = {metrics['d_min']:.4f}, d_max = {metrics['d_max']:.4f}")
    print(f"  Lz drift = {metrics['L_z_drift']:.2e}")
    print(f"  Max deviation from Newton: {dev['max_deviation']:.4e}")
    print(f"  RMS deviation from Newton: {dev['rms_deviation']:.4e}")

    # Compare V3/V2 ratio at initial configuration
    V2 = pot_fn(pos0, C_values)  # full TGP
    V2_only = float(potential_tgp(pos0, C_values, beta, beta, softening=1e-8))
    V3 = V2 - V2_only
    print(f"\n  V2 = {V2_only:.6e}")
    print(f"  V3 = {V3:.6e}")
    if abs(V2_only) > 1e-30:
        print(f"  V3/V2 = {abs(V3/V2_only):.4e}")

    return result


# -- Part 4: Orbit deviation quantification --------------------------------

def part4_deviation_analysis(newton_result, C_values, vel0_scaled):
    """Quantify how TGP deforms the figure-8 for different betas."""
    print(f"\n{'=' * 78}")
    print("Part 4: Orbit Deviation vs beta (Pairwise)")
    print(f"{'=' * 78}")

    pos0, _, _, _ = figure_eight_initial(C=C_DEFAULT)

    print(f"\n  {'beta':>8s} {'max_dev':>12s} {'rms_dev':>12s} {'period':>10s} "
          f"{'dE/E0':>10s} {'d_min':>8s} {'d_max':>8s} {'status':>8s}")
    print(f"  {'-'*82}")

    results = []
    for beta in BETAS:
        acc_fn, pot_fn = build_tgp_integration_pair(
            "pairwise", beta=beta, gamma=beta, softening=1e-8,
        )

        res = rk45_integrate(
            pos0, vel0_scaled, C_values,
            acc_fn, pot_fn,
            t_span=(0, T_FINAL), n_output=N_OUTPUT,
            rtol=1e-10, atol=1e-12, quiet=True,
        )

        max_E_err = np.max(np.abs(res['energy_error']))
        dev = _deviation_from_newton(newton_result['positions'],
                                     res['positions'],
                                     newton_result['t'], res['t'])
        period = _estimate_period(res['positions'], res['t'])
        metrics = _orbit_metrics(res['positions'], res['velocities'],
                                 res['t'], C_values)

        period_str = f"{period:.4f}" if period else "N/A"

        # Orbit "intact" if max deviation < characteristic size
        d_char = metrics['d_max']
        intact = dev['max_deviation'] < 0.5 * d_char
        status = "INTACT" if intact else "BROKEN"

        print(f"  {beta:8.4f} {dev['max_deviation']:12.4e} "
              f"{dev['rms_deviation']:12.4e} {period_str:>10s} "
              f"{max_E_err:10.2e} {metrics['d_min']:8.4f} "
              f"{metrics['d_max']:8.4f} {status:>8s}")

        results.append({
            'beta': beta,
            'max_dev': dev['max_deviation'],
            'rms_dev': dev['rms_deviation'],
            'period': period,
            'intact': intact,
            'max_E_err': max_E_err,
        })

    return results


# -- Part 5: beta threshold -- where does figure-8 break? -----------------

def part5_beta_threshold(results_part4):
    """Determine critical beta where figure-8 choreography breaks."""
    print(f"\n{'=' * 78}")
    print("Part 5: Critical beta for Figure-8 Breakdown")
    print(f"{'=' * 78}")

    intact_betas = [r['beta'] for r in results_part4 if r['intact']]
    broken_betas = [r['beta'] for r in results_part4 if not r['intact']]

    if intact_betas:
        print(f"\n  Figure-8 intact for beta <= {max(intact_betas):.4f}")
    if broken_betas:
        print(f"  Figure-8 broken for beta >= {min(broken_betas):.4f}")

    if intact_betas and broken_betas:
        beta_crit_lo = max(intact_betas)
        beta_crit_hi = min(broken_betas)
        print(f"  -> Critical beta in [{beta_crit_lo:.4f}, {beta_crit_hi:.4f}]")
    elif not broken_betas:
        print(f"  -> Figure-8 survives at all tested beta (max = {max(intact_betas):.4f})")
        print(f"     Consider testing larger beta values.")
    else:
        print(f"  -> Figure-8 broken at all tested beta (min = {min(broken_betas):.4f})")

    # Scaling analysis: deviation ~ beta^alpha
    if len(results_part4) >= 3:
        betas_arr = np.array([r['beta'] for r in results_part4])
        devs_arr = np.array([r['max_dev'] for r in results_part4])
        # Filter out zeros/nans
        mask = (devs_arr > 1e-15) & (betas_arr > 0)
        if np.sum(mask) >= 2:
            log_b = np.log(betas_arr[mask])
            log_d = np.log(devs_arr[mask])
            # Linear fit: log(dev) = alpha * log(beta) + const
            coeffs = np.polyfit(log_b, log_d, 1)
            alpha = coeffs[0]
            print(f"\n  Scaling: dr_max ~ beta^{alpha:.2f}")
            if abs(alpha - 1.0) < 0.3:
                print(f"  -> Consistent with linear perturbative regime (alpha ~ 1)")
            elif abs(alpha - 2.0) < 0.3:
                print(f"  -> Consistent with quadratic regime (alpha ~ 2)")


def main():
    print("=" * 78)
    print("ex203 -- P3.A: Figure-8 Orbit (Chenciner-Montgomery) -- TGP vs Newton")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")
    print(f"  T_final = {T_FINAL}, G_Newton = {G_NEWTON:.4f}")
    print(f"  C = {C_DEFAULT}")

    t_total = time.time()

    # Part 1: Newton reference
    newton_res, C_values, vel0_scaled = part1_newton_reference()

    # Part 2: TGP pairwise
    tgp_pw_res = part2_tgp_pairwise(newton_res, C_values, vel0_scaled)

    # Part 3: TGP full
    tgp_full_res = part3_tgp_full(newton_res, C_values, vel0_scaled)

    # Part 4: beta scan
    p4_results = part4_deviation_analysis(newton_res, C_values, vel0_scaled)

    # Part 5: threshold
    part5_beta_threshold(p4_results)

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")
    print(f"  Newton figure-8: {'OK' if newton_res['success'] else 'FAIL'}")
    n_intact = sum(1 for r in p4_results if r['intact'])
    n_total = len(p4_results)
    print(f"  TGP pairwise: {n_intact}/{n_total} beta values preserve figure-8")
    print(f"  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
