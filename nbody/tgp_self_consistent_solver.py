#!/usr/bin/env python3
"""
tgp_self_consistent_solver.py -- Full strong-field TGP BH soliton solver
=========================================================================

PHYSICAL SETUP:
  TGP metric (isotropic coordinates):
      ds^2 = -(c_0^2/psi) dt^2 + psi (dr^2 + r^2 dOmega^2)

  psi = Phi/Phi_0 is the conformal factor (substrate field / vacuum value).

  Field equation (vacuum, flat Laplacian):
      nabla^2 psi + 2(nabla psi)^2/psi + beta*psi^2*(1-psi) = 0

  where beta = gamma*Phi_0 (vacuum condition).

BLACK HOLE SOLITON:
  The BH in TGP is a topological soliton with:
      psi(r -> 0)   -> infinity   (frozen interior: g_tt -> 0)
      psi(r -> inf)  -> 1          (vacuum)
      psi(r) ~ 1 - r_s/r + ...    (Newtonian at large r)

  This is a TWO-POINT BOUNDARY VALUE PROBLEM.
  The exponential ansatz psi = exp(-r_s/r) -> 0 at center is WRONG for the BH.

NEAR-CENTER ASYMPTOTICS:
  For psi >> 1, the equation becomes:
      psi'' + (2/r)psi' + 2(psi')^2/psi - beta*psi^3 = 0

  Power-law ansatz psi = A/r gives:
      A(2 - A^2 beta) = 0  =>  A = sqrt(2/beta)

  So psi ~ sqrt(2/beta)/r near center.

SELF-CONSISTENT CORRECTION:
  The 4D covariant d'Alembertian introduces a factor psi from sqrt(-g):
      nabla^2 psi + 2(nabla psi)^2/psi + psi * beta*psi^2*(1-psi) = 0
  Equivalently: ... + beta*psi^3*(1-psi) = 0

  This makes the interior (psi>>1) MORE strongly repulsive and the
  transition to frozen interior sharper.

Author: Mateusz Serafin (with Claude)
Date: April 2026
"""

import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import interp1d
import sys


# ============================================================
# 1. Field equations (static, spherically symmetric)
# ============================================================

def rhs_flat(r, y, beta):
    """
    Flat TGP vacuum equation (published):
        psi'' + (2/r)psi' + 2(psi')^2/psi + beta*psi^2*(1-psi) = 0

    System: y[0] = psi, y[1] = psi'
    """
    psi, dpsi = y
    if psi < 1e-30:
        psi = 1e-30
    ddpsi = -(2.0/r)*dpsi - 2.0*dpsi**2/psi - beta*psi**2*(1.0 - psi)
    return [dpsi, ddpsi]


def rhs_selfconsistent(r, y, beta):
    """
    Self-consistent TGP equation (with time-dilation from sqrt(-g)):
        psi'' + (2/r)psi' + 2(psi')^2/psi + beta*psi^3*(1-psi) = 0

    The extra factor of psi in the potential comes from sqrt(-g) = c_0*psi.
    """
    psi, dpsi = y
    if psi < 1e-30:
        psi = 1e-30
    ddpsi = -(2.0/r)*dpsi - 2.0*dpsi**2/psi - beta*psi**3*(1.0 - psi)
    return [dpsi, ddpsi]


# ============================================================
# 2. Near-center asymptotics (psi >> 1)
# ============================================================

def inner_bc(r0, A, alpha=1.0):
    """
    Near-center power-law: psi = A * r^(-alpha)
    For the flat equation: alpha = 1, A = sqrt(2/beta)
    For self-consistent: alpha = 1, A = sqrt(2/beta) (same leading order)

    Returns (psi, psi') at r = r0.
    """
    psi = A * r0**(-alpha)
    dpsi = -alpha * A * r0**(-alpha - 1)
    return psi, dpsi


def inner_bc_corrected(r0, A, beta, alpha=1.0):
    """
    Next-order correction to the near-center asymptotics.

    psi = A/r + B + C*r + ...

    Substituting into the flat equation and matching O(1/r^3), O(1/r^2), O(1/r):

    O(1/r^3): A(2 + 2 - A^2*beta) = 0 => A^2*beta = 2 (skip: this is leading order)
                Wait let me redo. psi = A/r, psi' = -A/r^2, psi'' = 2A/r^3.
                (2/r)psi' = -2A/r^3
                2(psi')^2/psi = 2A/r^3
                kinetic = 2A/r^3
                potential = beta*(A/r)^2*(1 - A/r) = beta*A^2/r^2 - beta*A^3/r^3
                Full: 2A/r^3 + beta*A^2/r^2 - beta*A^3/r^3 = 0
                O(1/r^3): 2A - beta*A^3 = 0 => A = sqrt(2/beta) ✓
                O(1/r^2): beta*A^2 ... this comes from the constant B correction.

    For now, just use the leading order.
    """
    psi0 = A / r0
    dpsi0 = -A / r0**2

    # Add O(1) correction: psi = A/r + B
    # This shifts things a bit. B is determined by matching to the outer solution.
    # For now, leave B = 0 (will be adjusted by the shooting parameter).
    return psi0, dpsi0


# ============================================================
# 3. Outer boundary condition
# ============================================================

def outer_bc(r, r_s):
    """
    Weak-field boundary condition at large r.

    Newtonian: psi ~ 1 - r_s/r
    Exponential: psi ~ exp(-r_s/r) (keeps psi > 0)

    NOTE: This is the EXTERIOR solution (psi < 1), NOT the BH soliton.
    For the BH soliton, the outer BC is psi -> 1 with psi' ~ r_s/r^2.

    But wait: in the BH soliton, psi > 1 in the interior (frozen) and
    psi approaches 1 from ABOVE far away? No: psi -> 1 from below
    (Newtonian potential well: psi = 1 - r_s/r < 1).

    Actually: the soliton has psi -> infinity at center, must decrease
    through psi = 1 at some radius, then continue below 1, approaching
    1 from below at infinity.

    Wait: if psi -> 1 from below at infinity, and psi -> inf at center,
    then psi MUST pass through a MINIMUM below 1 and then cross 1 going
    up. This minimum is the "potential well" of the BH.

    Profile: psi(inf) = 1, psi decreases to psi_min < 1, then increases
    back through 1 and continues to infinity.
    """
    phi = np.exp(-r_s / r)  # standard exterior: psi < 1
    dphi = (r_s / r**2) * phi
    return phi, dphi


# ============================================================
# 4. BVP solver: shooting from both ends
# ============================================================

def shoot_from_center(A, beta, r_min=0.01, r_max=50.0, n_points=20000,
                      use_selfconsistent=False):
    """
    Shoot outward from r_min with inner BC: psi ~ A/r.
    """
    psi0, dpsi0 = inner_bc(r_min, A)
    y0 = [psi0, dpsi0]

    rhs = rhs_selfconsistent if use_selfconsistent else rhs_flat

    r_eval = np.linspace(r_min, r_max, n_points)

    def psi_tiny(r, y, *args):
        return y[0] - 1e-12
    psi_tiny.terminal = True
    psi_tiny.direction = -1

    def psi_huge(r, y, *args):
        return y[0] - 1e10
    psi_huge.terminal = True
    psi_huge.direction = 1

    sol = solve_ivp(
        lambda r, y: rhs(r, y, beta),
        [r_min, r_max],
        y0,
        method='Radau',
        t_eval=r_eval,
        rtol=1e-11,
        atol=1e-14,
        max_step=0.05,
        events=[psi_tiny, psi_huge]
    )

    return sol.t, sol.y[0], sol.y[1], sol


def shoot_from_infinity(r_s, beta, r_max=500.0, r_min=0.1, n_points=20000,
                        use_selfconsistent=False):
    """
    Shoot inward from r_max with outer BC: psi ~ exp(-r_s/r).
    """
    psi0, dpsi0 = outer_bc(r_max, r_s)
    y0 = [psi0, dpsi0]

    rhs = rhs_selfconsistent if use_selfconsistent else rhs_flat

    r_eval = np.linspace(r_max, r_min, n_points)

    def psi_tiny(r, y, *args):
        return y[0] - 1e-12
    psi_tiny.terminal = True
    psi_tiny.direction = -1

    sol = solve_ivp(
        lambda r, y: rhs(r, y, beta),
        [r_max, r_min],
        y0,
        method='Radau',
        t_eval=r_eval,
        rtol=1e-11,
        atol=1e-14,
        max_step=0.05,
        events=[psi_tiny]
    )

    return sol.t, sol.y[0], sol.y[1], sol


# ============================================================
# 5. Matching: find A such that inner and outer solutions connect
# ============================================================

def match_solutions(r_s, beta, r_match, r_min=0.005, r_max=500.0,
                    use_selfconsistent=False):
    """
    Find the amplitude A of the inner solution such that:
    - Inner solution at r_match = outer solution at r_match
    - Derivative continuity (smooth matching)

    We scan A and look for the value where psi_inner(r_match) = psi_outer(r_match).
    """
    print(f"  Matching at r_match = {r_match}")

    # First: compute outer solution
    r_out, psi_out, dpsi_out, sol_out = shoot_from_infinity(
        r_s, beta, r_max=r_max, r_min=r_match*0.5,
        n_points=20000, use_selfconsistent=use_selfconsistent
    )

    if len(r_out) < 10:
        print(f"  Outer solution too short ({len(r_out)} points)")
        return None

    # Interpolate outer solution at r_match
    # r_out is DECREASING, so reverse for interp
    try:
        psi_out_interp = np.interp(r_match, r_out[::-1], psi_out[::-1])
        dpsi_out_interp = np.interp(r_match, r_out[::-1], dpsi_out[::-1])
    except Exception:
        print(f"  Outer interpolation failed")
        return None

    print(f"  Outer solution at r_match: psi = {psi_out_interp:.6f}, "
          f"psi' = {dpsi_out_interp:.6e}")

    # Scan A values to find matching
    A_nominal = np.sqrt(2.0 / beta)
    A_values = np.linspace(0.1 * A_nominal, 5.0 * A_nominal, 50)

    psi_inner_at_match = []
    for A in A_values:
        try:
            r_in, psi_in, dpsi_in, sol_in = shoot_from_center(
                A, beta, r_min=r_min, r_max=r_match*2,
                n_points=5000, use_selfconsistent=use_selfconsistent
            )
            if len(r_in) < 10 or r_in[-1] < r_match:
                psi_inner_at_match.append(np.nan)
                continue
            psi_val = np.interp(r_match, r_in, psi_in)
            psi_inner_at_match.append(psi_val)
        except Exception:
            psi_inner_at_match.append(np.nan)

    psi_inner_at_match = np.array(psi_inner_at_match)

    # Find where psi_inner crosses psi_outer
    target = psi_out_interp
    diff = psi_inner_at_match - target

    valid = ~np.isnan(diff)
    if np.sum(valid) < 2:
        print(f"  Not enough valid inner solutions")
        return None

    # Look for sign changes
    A_valid = A_values[valid]
    diff_valid = diff[valid]

    crossings = []
    for i in range(len(diff_valid)-1):
        if diff_valid[i] * diff_valid[i+1] < 0:
            # Refine with bisection
            A_lo, A_hi = A_valid[i], A_valid[i+1]

            def mismatch(A_try):
                r_in, psi_in, _, _ = shoot_from_center(
                    A_try, beta, r_min=r_min, r_max=r_match*2,
                    n_points=5000, use_selfconsistent=use_selfconsistent
                )
                if len(r_in) < 10 or r_in[-1] < r_match:
                    return float('nan')
                return np.interp(r_match, r_in, psi_in) - target

            try:
                A_star = brentq(mismatch, A_lo, A_hi, xtol=1e-8)
                crossings.append(A_star)
            except Exception:
                pass

    if crossings:
        print(f"  Found {len(crossings)} matching amplitude(s):")
        for A_star in crossings:
            print(f"    A* = {A_star:.6f} (nominal: {A_nominal:.6f})")
        return crossings
    else:
        print(f"  No matching found")
        print(f"  psi_inner range at match: [{np.nanmin(psi_inner_at_match):.4f}, "
              f"{np.nanmax(psi_inner_at_match):.4f}]")
        print(f"  Target psi_outer at match: {target:.4f}")
        return None


# ============================================================
# 6. Full soliton construction
# ============================================================

def construct_soliton(A_star, r_s, beta, r_min=0.005, r_max=500.0,
                      n_points=50000, use_selfconsistent=False):
    """
    Given the matched amplitude A*, construct the full soliton profile.
    """
    # Inner solution
    r_in, psi_in, dpsi_in, _ = shoot_from_center(
        A_star, beta, r_min=r_min, r_max=r_max,
        n_points=n_points, use_selfconsistent=use_selfconsistent
    )

    # Find the matching point (where psi starts to deviate from physical)
    # For the inner solution, it should reach psi = 1 at some radius
    # and then... depends on the dynamics

    # Actually, the inner solution IS the full soliton if A is right
    # because it starts with psi -> inf at r_min and evolves outward
    # toward psi -> 1 at large r.

    return r_in, psi_in, dpsi_in


# ============================================================
# 7. Photon sphere and shadow
# ============================================================

def find_photon_sphere(r_arr, psi_arr):
    """
    Photon sphere: minimum of h(r) = r * psi(r).

    Condition: d/dr[r*psi] = 0 => psi + r*psi' = 0 => psi'/psi = -1/r

    For the BH soliton with psi -> inf at center and psi -> 1 at inf:
        h(0+) = 0 * inf (depends on divergence rate)
        h(inf) = inf

    If psi ~ A/r: h = A (constant) at center.
    Then h decreases from A to some minimum and increases to inf.
    """
    h = r_arr * psi_arr

    # Use numerical gradient
    dh = np.gradient(h, r_arr)

    photon_spheres = []
    for i in range(1, len(dh)-1):
        if dh[i-1] < 0 and dh[i+1] > 0:
            r_ph = r_arr[i]
            psi_ph = psi_arr[i]
            h_ph = h[i]
            R_ph = r_ph * np.sqrt(psi_ph)  # areal radius
            photon_spheres.append({
                'r_ph': r_ph,
                'R_ph': R_ph,
                'psi_ph': psi_ph,
                'b_crit': h_ph,  # impact parameter ~ r*psi
            })

    return photon_spheres


def shadow_ratio(ps_list, r_s):
    """Compare shadow with GR Schwarzschild."""
    b_gr = (3.0 * np.sqrt(3.0) / 2.0) * r_s

    for ps in ps_list:
        ratio = ps['b_crit'] / b_gr
        ps['shadow_ratio'] = ratio
        ps['percent_diff'] = (1.0 - ratio) * 100

    return ps_list


# ============================================================
# 8. Main analysis
# ============================================================

def analyze_soliton(beta=1.0, r_min=0.005, r_max=100.0, n_points=30000):
    """
    Full analysis of the TGP BH soliton.
    """
    sep = "=" * 72
    A_nominal = np.sqrt(2.0 / beta)

    print(sep)
    print("TGP BLACK HOLE SOLITON ANALYSIS")
    print(f"  beta = {beta}")
    print(f"  A_nominal = sqrt(2/beta) = {A_nominal:.6f}")
    print(f"  Integration: r in [{r_min}, {r_max}]")
    print(sep)

    # --- Step 1: Shoot from center with nominal A ---
    print("\n--- Step 1: Inner solution (psi ~ A/r at center) ---")

    A_values = [0.5*A_nominal, 0.8*A_nominal, A_nominal,
                1.2*A_nominal, 1.5*A_nominal, 2.0*A_nominal]

    best_A = None
    best_r_max_reached = 0

    for A in A_values:
        r_arr, psi_arr, dpsi_arr, sol = shoot_from_center(
            A, beta, r_min=r_min, r_max=r_max, n_points=n_points
        )

        r_reached = r_arr[-1]
        psi_final = psi_arr[-1] if len(psi_arr) > 0 else float('nan')

        # Find where psi crosses 1 (if it does)
        crossings_1 = []
        for i in range(len(psi_arr)-1):
            if (psi_arr[i] - 1.0) * (psi_arr[i+1] - 1.0) < 0:
                r_cross = r_arr[i]
                crossings_1.append(r_cross)

        # Find minimum of psi
        psi_min = psi_arr.min()
        r_psi_min = r_arr[np.argmin(psi_arr)]

        # h = r*psi analysis
        h = r_arr * psi_arr
        h_min = h.min()
        r_h_min = r_arr[np.argmin(h)]

        status = f"reached r={r_reached:.2f}, psi_final={psi_final:.4f}"
        if crossings_1:
            status += f", crosses 1 at r={crossings_1[0]:.4f}"
        status += f", psi_min={psi_min:.4f} at r={r_psi_min:.4f}"
        status += f", h_min={h_min:.4f} at r={r_h_min:.4f}"

        print(f"  A = {A:.4f}: {status}")

        if r_reached > best_r_max_reached:
            best_r_max_reached = r_reached
            best_A = A

    # --- Step 2: Detailed analysis of best solution ---
    print(f"\n--- Step 2: Detailed analysis (A = {best_A:.4f}) ---")

    r_arr, psi_arr, dpsi_arr, sol = shoot_from_center(
        best_A, beta, r_min=r_min, r_max=r_max*5, n_points=n_points*2
    )

    print(f"  Integration: {len(r_arr)} points, r in [{r_arr[0]:.4f}, {r_arr[-1]:.4f}]")
    print(f"  psi range: [{psi_arr.min():.6f}, {psi_arr.max():.6f}]")

    # Profile description
    h = r_arr * psi_arr
    print(f"  h(r) = r*psi(r) range: [{h.min():.6f}, {h.max():.6f}]")

    # Look for photon sphere
    ps_list = find_photon_sphere(r_arr, psi_arr)

    if ps_list:
        print(f"\n  *** PHOTON SPHERE FOUND ***")
        for i, ps in enumerate(ps_list):
            print(f"    [{i}] r_ph = {ps['r_ph']:.6f}, psi(r_ph) = {ps['psi_ph']:.6f}")
            print(f"         R_ph = {ps['R_ph']:.6f} (areal)")
            print(f"         b_crit = {ps['b_crit']:.6f}")
    else:
        print(f"\n  No photon sphere found")

    # --- Step 3: Scan A more finely ---
    print(f"\n--- Step 3: Fine A scan ---")

    A_fine = np.linspace(0.3*A_nominal, 3.0*A_nominal, 100)
    results_fine = []

    for A in A_fine:
        try:
            r_a, psi_a, _, sol_a = shoot_from_center(
                A, beta, r_min=r_min, r_max=r_max, n_points=10000
            )
            if len(r_a) < 100:
                continue

            # Key metrics
            h_a = r_a * psi_a
            h_min_a = h_a.min()
            r_h_min_a = r_a[np.argmin(h_a)]
            psi_final_a = psi_a[-1]
            r_reached_a = r_a[-1]

            # Does psi approach 1 at large r?
            asymp_deviation = abs(psi_final_a - 1.0)

            results_fine.append({
                'A': A,
                'h_min': h_min_a,
                'r_h_min': r_h_min_a,
                'psi_final': psi_final_a,
                'r_reached': r_reached_a,
                'asymp_dev': asymp_deviation,
            })
        except Exception:
            pass

    if results_fine:
        # Find A that gives psi closest to 1 at large r
        best_asymp = min(results_fine, key=lambda x: x['asymp_dev'])
        print(f"  Best asymptotic match: A = {best_asymp['A']:.6f}, "
              f"psi(r_max) = {best_asymp['psi_final']:.6f}, "
              f"h_min = {best_asymp['h_min']:.6f}")

        # Find A that gives h_min (deepest potential well)
        deepest = min(results_fine, key=lambda x: x['h_min'])
        print(f"  Deepest potential well: A = {deepest['A']:.6f}, "
              f"h_min = {deepest['h_min']:.6f} at r = {deepest['r_h_min']:.4f}")

        # Reconstruct best solution
        A_best = best_asymp['A']
        print(f"\n--- Reconstructing best soliton (A = {A_best:.6f}) ---")

        r_best, psi_best, dpsi_best, _ = shoot_from_center(
            A_best, beta, r_min=r_min, r_max=r_max*5, n_points=50000
        )

        h_best = r_best * psi_best
        ps_best = find_photon_sphere(r_best, psi_best)

        print(f"  Solution: {len(r_best)} points, r in [{r_best[0]:.4f}, {r_best[-1]:.4f}]")
        print(f"  psi(r_min) = {psi_best[0]:.4f}, psi(r_max) = {psi_best[-1]:.6f}")
        print(f"  h_min = {h_best.min():.6f} at r = {r_best[np.argmin(h_best)]:.4f}")

        if ps_best:
            # Estimate r_s from the Newtonian tail
            # At large r: psi ~ 1 - r_s/r => r_s ~ r*(1-psi)
            r_large = r_best[-100:]
            psi_large = psi_best[-100:]
            r_s_est = np.median(r_large * (1.0 - psi_large))
            print(f"  Estimated r_s (from tail): {r_s_est:.4f}")

            ps_best = shadow_ratio(ps_best, r_s_est)

            print(f"\n  *** PHOTON SPHERE ***")
            for ps in ps_best:
                print(f"    r_ph = {ps['r_ph']:.6f}")
                print(f"    psi(r_ph) = {ps['psi_ph']:.6f}")
                print(f"    b_crit = {ps['b_crit']:.6f}")
                print(f"    b_GR = {(3*np.sqrt(3)/2)*r_s_est:.6f}")
                sign = "smaller" if ps['shadow_ratio'] < 1 else "larger"
                print(f"    Shadow: {abs(ps['percent_diff']):.2f}% {sign} than GR")
        else:
            print(f"  No photon sphere in best solution")

        return {
            'r': r_best, 'psi': psi_best, 'dpsi': dpsi_best,
            'photon_spheres': ps_best,
            'A': A_best,
        }

    return None


# ============================================================
# 9. Self-consistent comparison
# ============================================================

def compare_flat_vs_sc(beta=1.0, r_min=0.005, r_max=100.0, n_points=30000):
    """
    Compare flat and self-consistent equations.
    """
    sep = "=" * 72
    A_nom = np.sqrt(2.0 / beta)

    print(f"\n{sep}")
    print("FLAT vs SELF-CONSISTENT COMPARISON")
    print(sep)

    for label, sc_flag in [("FLAT", False), ("SELF-CONSISTENT", True)]:
        print(f"\n--- {label} equation ---")

        # Scan A
        A_values = np.linspace(0.5*A_nom, 3.0*A_nom, 80)
        best_result = None
        best_dev = float('inf')

        for A in A_values:
            try:
                r_a, psi_a, _, _ = shoot_from_center(
                    A, beta, r_min=r_min, r_max=r_max,
                    n_points=10000, use_selfconsistent=sc_flag
                )
                if len(r_a) < 100 or r_a[-1] < r_max * 0.5:
                    continue

                dev = abs(psi_a[-1] - 1.0)
                if dev < best_dev:
                    best_dev = dev
                    best_result = (A, r_a, psi_a)
            except Exception:
                pass

        if best_result:
            A_best, r_best, psi_best = best_result
            h_best = r_best * psi_best
            ps = find_photon_sphere(r_best, psi_best)

            print(f"  Best A = {A_best:.6f}")
            print(f"  psi(r_max) = {psi_best[-1]:.6f} (target: 1.0)")
            print(f"  h_min = {h_best.min():.6f}")

            if ps:
                # Estimate r_s
                r_large = r_best[-50:]
                psi_large = psi_best[-50:]
                r_s_est = np.median(r_large * (1.0 - psi_large))
                ps = shadow_ratio(ps, max(r_s_est, 0.01))

                for p in ps:
                    sign = "smaller" if p['shadow_ratio'] < 1 else "larger"
                    print(f"  PHOTON SPHERE: r_ph = {p['r_ph']:.4f}, "
                          f"shadow {abs(p['percent_diff']):.1f}% {sign}")
            else:
                print(f"  No photon sphere")
        else:
            print(f"  No valid solution found")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    if sys.stdout.encoding != 'utf-8':
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    quick = "--quick" in sys.argv
    compare = "--compare" in sys.argv

    if compare:
        compare_flat_vs_sc(beta=1.0, r_min=0.01, r_max=50.0, n_points=15000)
    elif quick:
        analyze_soliton(beta=1.0, r_min=0.01, r_max=50.0, n_points=10000)
    else:
        # Full analysis
        result = analyze_soliton(beta=1.0, r_min=0.005, r_max=100.0, n_points=30000)

        # Also compare with self-consistent
        compare_flat_vs_sc(beta=1.0, r_min=0.005, r_max=100.0, n_points=20000)

    print("\nDone.")
