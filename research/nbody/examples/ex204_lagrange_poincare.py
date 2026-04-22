#!/usr/bin/env python3
"""
ex204 -- P3.B: Lagrange points & Poincare sections -- TGP vs Newton
====================================================================

Restricted 3-body problem analysis:
  Part 1: Lagrange point locations (L1-L5) in TGP vs Newton
  Part 2: Linear stability of L4, L5 with TGP corrections
  Part 3: Poincare section of test-mass orbit near L4
  Part 4: Effective potential landscape comparison

PHYSICS:
--------
In Newton's restricted 3-body problem (two heavy + one light mass in
co-rotating frame), the Lagrange points L1-L3 (collinear) are always
unstable, while L4, L5 (equilateral) are stable for mass ratio
q = m2/(m1+m2) < q_crit = (1 - sqrt(23/27))/2 ~ 0.0385.

TGP adds:
  - Repulsive beta/r^2 barrier: shifts Lagrange point locations
  - Attractive gamma/r^3: shifts them back (partially)
  - Irreducible 3-body V3: additional coupling

Key question: does TGP widen or narrow the stability window of L4/L5?

POINCARE SECTION:
-----------------
For a test mass orbiting near L4 in the co-rotating frame,
we reduce to Jacobi integral = const, then plot (x, vx) at y=y_L4 crossings.
  - Regular orbits -> closed curves
  - Chaotic orbits -> scattered points
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import potential_tgp, potential_newton
from nbody.dynamics_backends import build_tgp_integration_pair

quick = "--quick" in sys.argv

G_NEWTON = 4.0 * np.pi

if quick:
    N_GRID = 80
    N_OUTPUT = 500
    T_FINAL_POINCARE = 50.0
    BETAS = [0.0, 0.01, 0.05]
    N_ORBITS = 3
else:
    N_GRID = 200
    N_OUTPUT = 2000
    T_FINAL_POINCARE = 200.0
    BETAS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10]
    N_ORBITS = 6


# -- Restricted 3-body effective potential ---------------------------------

def _effective_potential_newton(x, y, q, Omega, G=G_NEWTON):
    """
    Effective potential in co-rotating frame (Newton).

    Two primaries: m1 = (1-q)*M at (-q*a, 0), m2 = q*M at ((1-q)*a, 0)
    where a = separation, M = total mass, Omega = angular velocity.

    We set a=1, M=1 for simplicity, so Omega^2 = G*M/a^3 = G.

    Phi_eff = -G*m1/r1 - G*m2/r2 - 0.5*Omega^2*(x^2 + y^2)
    """
    m1, m2 = 1.0 - q, q
    # Primary positions (in rotating frame, origin at CoM)
    x1, y1 = -q, 0.0
    x2, y2 = 1.0 - q, 0.0

    r1 = np.sqrt((x - x1)**2 + (y - y1)**2 + 1e-10)
    r2 = np.sqrt((x - x2)**2 + (y - y2)**2 + 1e-10)

    return -G * m1 / r1 - G * m2 / r2 - 0.5 * Omega**2 * (x**2 + y**2)


def _effective_potential_tgp(x, y, q, Omega, beta, G=G_NEWTON):
    """
    Effective potential in co-rotating frame (TGP pairwise).

    V_TGP(r) = -4*pi*Ci*Cj/r + 8*pi*beta*Ci*Cj/r^2 - 12*pi*gamma*Ci*Cj*(Ci+Cj)/r^3

    In restricted problem, test mass C3 << C1, C2, so:
    V(r_i) = -4*pi*C3*Ci/ri + 8*pi*beta*C3*Ci/ri^2 - 12*pi*beta*C3*Ci*(C3+Ci)/ri^3

    We use C1 = 1-q, C2 = q, C3 = epsilon (test mass limit).
    """
    gamma = beta  # vacuum condition
    m1, m2 = 1.0 - q, q
    C3 = 1e-4  # test mass

    x1, y1 = -q, 0.0
    x2, y2 = 1.0 - q, 0.0

    r1 = np.sqrt((x - x1)**2 + (y - y1)**2 + 1e-10)
    r2 = np.sqrt((x - x2)**2 + (y - y2)**2 + 1e-10)

    # Newtonian part (leading order)
    V_N = -G * m1 / r1 - G * m2 / r2

    # TGP corrections (beta/r^2 and gamma/r^3)
    V_beta = 8.0 * np.pi * beta * C3 * m1 / r1**2 + 8.0 * np.pi * beta * C3 * m2 / r2**2
    V_gamma = (-12.0 * np.pi * gamma * C3 * m1 * (C3 + m1) / r1**3
               - 12.0 * np.pi * gamma * C3 * m2 * (C3 + m2) / r2**3)

    # Centrifugal
    V_rot = -0.5 * Omega**2 * (x**2 + y**2)

    return V_N + V_beta + V_gamma + V_rot


def _lagrange_points_newton(q, G=G_NEWTON):
    """
    Find approximate Lagrange point locations for mass ratio q.

    L4 = (0.5 - q, sqrt(3)/2) and L5 = (0.5 - q, -sqrt(3)/2) are exact.
    L1, L2, L3 are found numerically on the x-axis.

    Returns dict with L1..L5 as (x, y) tuples.
    """
    from scipy.optimize import brentq

    Omega2 = G  # with a=1, M=1
    m1, m2 = 1.0 - q, q

    # L4, L5 (exact for Newton)
    L4 = (0.5 - q, np.sqrt(3.0) / 2.0)
    L5 = (0.5 - q, -np.sqrt(3.0) / 2.0)

    # For collinear points, we solve dPhi_eff/dx = 0 on y=0
    def dPhi_dx(x):
        r1 = abs(x - (-q))
        r2 = abs(x - (1.0 - q))
        if r1 < 1e-8 or r2 < 1e-8:
            return 1e10
        s1 = np.sign(x - (-q))
        s2 = np.sign(x - (1.0 - q))
        return (-G * m1 * s1 / r1**2 - G * m2 * s2 / r2**2 + Omega2 * x)

    # L1: between the two primaries
    try:
        x_L1 = brentq(dPhi_dx, -q + 0.01, 1.0 - q - 0.01)
    except (ValueError, RuntimeError):
        x_L1 = 0.5 - q  # fallback

    # L2: beyond the smaller primary
    try:
        x_L2 = brentq(dPhi_dx, 1.0 - q + 0.01, 1.0 - q + 2.0)
    except (ValueError, RuntimeError):
        x_L2 = 1.0 - q + 0.5

    # L3: beyond the larger primary (opposite side)
    try:
        x_L3 = brentq(dPhi_dx, -q - 2.0, -q - 0.01)
    except (ValueError, RuntimeError):
        x_L3 = -q - 1.0

    return {
        'L1': (x_L1, 0.0),
        'L2': (x_L2, 0.0),
        'L3': (x_L3, 0.0),
        'L4': L4,
        'L5': L5,
    }


def _lagrange_L4_tgp(q, beta, G=G_NEWTON, tol=1e-8):
    """
    Find L4 in TGP by minimizing gradient of effective potential.
    Start from Newton L4 and refine.
    """
    from scipy.optimize import minimize

    Omega2 = G
    Omega = np.sqrt(Omega2)

    # Start from Newton L4
    x0 = 0.5 - q
    y0 = np.sqrt(3.0) / 2.0

    def neg_potential(xy):
        # We want saddle/extremum, so find zero of gradient
        return _effective_potential_tgp(xy[0], xy[1], q, Omega, beta, G)

    def gradient(xy):
        x, y = xy
        eps = 1e-7
        dpdx = (_effective_potential_tgp(x+eps, y, q, Omega, beta, G)
                - _effective_potential_tgp(x-eps, y, q, Omega, beta, G)) / (2*eps)
        dpdy = (_effective_potential_tgp(x, y+eps, q, Omega, beta, G)
                - _effective_potential_tgp(x, y-eps, q, Omega, beta, G)) / (2*eps)
        return np.array([dpdx, dpdy])

    # L4 is a saddle point of Phi_eff, so we find gradient = 0
    from scipy.optimize import fsolve
    try:
        result = fsolve(gradient, [x0, y0], full_output=True)
        xy_sol = result[0]
        info = result[1]
        return (xy_sol[0], xy_sol[1])
    except Exception:
        return (x0, y0)


def _hessian_at_point(pot_fn, x, y, q, Omega, beta, eps=1e-5):
    """Compute 2x2 Hessian of effective potential at (x, y)."""
    H = np.zeros((2, 2))
    V0 = pot_fn(x, y, q, Omega, beta)

    # d2V/dx2
    H[0, 0] = (pot_fn(x+eps, y, q, Omega, beta)
               - 2*V0
               + pot_fn(x-eps, y, q, Omega, beta)) / eps**2
    # d2V/dy2
    H[1, 1] = (pot_fn(x, y+eps, q, Omega, beta)
               - 2*V0
               + pot_fn(x, y-eps, q, Omega, beta)) / eps**2
    # d2V/dxdy
    H[0, 1] = (pot_fn(x+eps, y+eps, q, Omega, beta)
               - pot_fn(x+eps, y-eps, q, Omega, beta)
               - pot_fn(x-eps, y+eps, q, Omega, beta)
               + pot_fn(x-eps, y-eps, q, Omega, beta)) / (4*eps**2)
    H[1, 0] = H[0, 1]
    return H


def _linear_stability_L4(q, beta, G=G_NEWTON):
    """
    Linear stability of L4 in the co-rotating frame.

    In the rotating frame, linearized equations about L4:
      d2x/dt2 - 2*Omega*dy/dt = -d2Phi/dx2 * dx - d2Phi/dxdy * dy
      d2y/dt2 + 2*Omega*dx/dt = -d2Phi/dxdy * dx - d2Phi/dy2 * dy

    The eigenvalue equation for lambda:
      lambda^4 + (4*Omega^2 - Uxx - Uyy)*lambda^2 + Uxx*Uyy - Uxy^2 = 0

    Stability requires both roots lambda^2 to be negative (pure imaginary).
    """
    Omega2 = G
    Omega = np.sqrt(Omega2)

    if beta == 0.0:
        # Newton case: use exact Newton potential
        L4 = _lagrange_points_newton(q, G)['L4']
        x4, y4 = L4

        def pot_newton_wrap(x, y, q_unused, Omega_unused, beta_unused):
            return _effective_potential_newton(x, y, q, Omega, G)
        H = _hessian_at_point(pot_newton_wrap, x4, y4, q, Omega, 0.0)
    else:
        L4 = _lagrange_L4_tgp(q, beta, G)
        x4, y4 = L4
        H = _hessian_at_point(_effective_potential_tgp, x4, y4, q, Omega, beta)

    Uxx, Uyy, Uxy = H[0, 0], H[1, 1], H[0, 1]

    # Characteristic equation: s^2 + (4*Omega^2 - Uxx - Uyy)*s + Uxx*Uyy - Uxy^2 = 0
    # where s = lambda^2
    a_coeff = 1.0
    b_coeff = 4.0 * Omega2 - Uxx - Uyy
    c_coeff = Uxx * Uyy - Uxy**2

    disc = b_coeff**2 - 4*a_coeff*c_coeff

    if disc < 0:
        # Complex s -> always unstable
        stable = False
        s1 = s2 = None
    else:
        s1 = (-b_coeff + np.sqrt(disc)) / 2.0
        s2 = (-b_coeff - np.sqrt(disc)) / 2.0
        # Stable if both s1, s2 < 0 (lambda = imaginary)
        stable = (s1 < 0) and (s2 < 0)

    return {
        'L4': (x4, y4),
        'Hessian': H,
        'Uxx': Uxx, 'Uyy': Uyy, 'Uxy': Uxy,
        's1': s1, 's2': s2,
        'discriminant': disc,
        'stable': stable,
    }


# -- Poincare section (co-rotating frame) ---------------------------------

def _integrate_corotating(x0, y0, vx0, vy0, q, Omega, beta, t_final, n_output, G=G_NEWTON):
    """
    Integrate test-mass motion in co-rotating frame.

    EOM in rotating frame:
      x'' = 2*Omega*y' - dPhi/dx
      y'' = -2*Omega*x' - dPhi/dy
    """
    from scipy.integrate import solve_ivp

    def gradient(x, y):
        eps = 1e-7
        if beta == 0.0:
            Phi = lambda xx, yy: _effective_potential_newton(xx, yy, q, Omega, G)
        else:
            Phi = lambda xx, yy: _effective_potential_tgp(xx, yy, q, Omega, beta, G)
        dpdx = (Phi(x+eps, y) - Phi(x-eps, y)) / (2*eps)
        dpdy = (Phi(x, y+eps) - Phi(x, y-eps)) / (2*eps)
        return dpdx, dpdy

    def ode_rhs(t, state):
        x, y, vx, vy = state
        dpdx, dpdy = gradient(x, y)
        ax = 2.0 * Omega * vy - dpdx
        ay = -2.0 * Omega * vx - dpdy
        return [vx, vy, ax, ay]

    y0_vec = [x0, y0, vx0, vy0]
    t_eval = np.linspace(0, t_final, n_output)

    sol = solve_ivp(ode_rhs, (0, t_final), y0_vec, t_eval=t_eval,
                    method='DOP853', rtol=1e-10, atol=1e-12)
    return sol


def _poincare_section(sol, y_section):
    """
    Extract Poincare section: crossings of y = y_section with dy/dt > 0.

    Returns (x_cross, vx_cross) arrays.
    """
    x_arr = sol.y[0]
    y_arr = sol.y[1]
    vx_arr = sol.y[2]
    vy_arr = sol.y[3]

    x_cross = []
    vx_cross = []

    for i in range(len(y_arr) - 1):
        if (y_arr[i] - y_section) * (y_arr[i+1] - y_section) < 0 and vy_arr[i] > 0:
            # Linear interpolation to find crossing
            frac = (y_section - y_arr[i]) / (y_arr[i+1] - y_arr[i])
            x_c = x_arr[i] + frac * (x_arr[i+1] - x_arr[i])
            vx_c = vx_arr[i] + frac * (vx_arr[i+1] - vx_arr[i])
            x_cross.append(x_c)
            vx_cross.append(vx_c)

    return np.array(x_cross), np.array(vx_cross)


# -- Part 1: Lagrange point locations -------------------------------------

def part1_lagrange_points():
    """Compare Lagrange point locations: Newton vs TGP."""
    print("=" * 78)
    print("Part 1: Lagrange Point Locations -- Newton vs TGP")
    print("=" * 78)

    q = 0.01  # Small mass ratio (Sun-Jupiter-like)

    L_newton = _lagrange_points_newton(q)

    print(f"  Mass ratio q = {q}")
    print(f"  G = {G_NEWTON:.4f}")
    print(f"\n  Newton Lagrange points:")
    for name, (x, y) in sorted(L_newton.items()):
        print(f"    {name}: ({x:+.6f}, {y:+.6f})")

    print(f"\n  TGP shifts (beta -> L4 location):")
    print(f"  {'beta':>8s} {'x_L4':>10s} {'y_L4':>10s} {'dx':>10s} {'dy':>10s}")
    print(f"  {'-'*52}")

    x4_N, y4_N = L_newton['L4']
    results = []
    for beta in BETAS:
        if beta == 0:
            x4, y4 = x4_N, y4_N
        else:
            x4, y4 = _lagrange_L4_tgp(q, beta)
        dx = x4 - x4_N
        dy = y4 - y4_N
        print(f"  {beta:8.4f} {x4:10.6f} {y4:10.6f} {dx:10.2e} {dy:10.2e}")
        results.append({'beta': beta, 'x_L4': x4, 'y_L4': y4, 'dx': dx, 'dy': dy})

    return results, q


# -- Part 2: Linear stability of L4 ---------------------------------------

def part2_stability():
    """Linear stability analysis of L4 for various q and beta."""
    print(f"\n{'=' * 78}")
    print("Part 2: Linear Stability of L4")
    print(f"{'=' * 78}")

    # Newton critical mass ratio
    q_crit_newton = (1.0 - np.sqrt(23.0/27.0)) / 2.0
    print(f"  Newton critical q = {q_crit_newton:.6f}")

    if quick:
        q_values = [0.01, 0.02, 0.03, 0.038]
    else:
        q_values = [0.005, 0.01, 0.02, 0.03, 0.035, 0.038, 0.039, 0.04, 0.05]

    print(f"\n  {'q':>6s} {'beta':>8s} {'stable':>8s} {'s1':>12s} {'s2':>12s} "
          f"{'L4_x':>10s} {'L4_y':>10s}")
    print(f"  {'-'*68}")

    results = []
    for q in q_values:
        for beta in BETAS:
            info = _linear_stability_L4(q, beta)
            s1_str = f"{info['s1']:.4e}" if info['s1'] is not None else "complex"
            s2_str = f"{info['s2']:.4e}" if info['s2'] is not None else "complex"
            stab_str = "STABLE" if info['stable'] else "UNSTABLE"

            print(f"  {q:6.3f} {beta:8.4f} {stab_str:>8s} {s1_str:>12s} {s2_str:>12s} "
                  f"{info['L4'][0]:10.6f} {info['L4'][1]:10.6f}")

            results.append({
                'q': q, 'beta': beta,
                'stable': info['stable'],
                's1': info['s1'], 's2': info['s2'],
            })

    # Summary: does TGP change the stability boundary?
    print(f"\n  Stability boundary analysis:")
    for beta in BETAS:
        beta_results = [r for r in results if r['beta'] == beta]
        stable_qs = [r['q'] for r in beta_results if r['stable']]
        unstable_qs = [r['q'] for r in beta_results if not r['stable']]
        if stable_qs and unstable_qs:
            q_max_stable = max(stable_qs)
            q_min_unstable = min(unstable_qs)
            print(f"    beta = {beta:.4f}: q_crit in [{q_max_stable:.4f}, {q_min_unstable:.4f}]")
        elif stable_qs:
            print(f"    beta = {beta:.4f}: stable for all tested q (max = {max(stable_qs):.4f})")
        else:
            print(f"    beta = {beta:.4f}: unstable for all tested q")

    return results


# -- Part 3: Poincare sections near L4 ------------------------------------

def part3_poincare():
    """Poincare sections of test-mass orbits near L4."""
    print(f"\n{'=' * 78}")
    print("Part 3: Poincare Sections near L4")
    print(f"{'=' * 78}")

    q = 0.01
    Omega = np.sqrt(G_NEWTON)
    L4_newton = _lagrange_points_newton(q)['L4']
    x_L4, y_L4 = L4_newton

    print(f"  q = {q}, L4 = ({x_L4:.4f}, {y_L4:.4f})")
    print(f"  T_final = {T_FINAL_POINCARE}")

    results = {}
    for beta in [0.0, BETAS[-1]]:  # Newton and largest beta
        label = f"Newton" if beta == 0.0 else f"TGP(beta={beta})"
        print(f"\n  === {label} ===")

        all_crossings_x = []
        all_crossings_vx = []

        # Launch orbits with different initial displacements from L4
        displacements = np.linspace(0.02, 0.10, N_ORBITS)

        for i, dx in enumerate(displacements):
            # Start displaced from L4 in x-direction, zero velocity perturbation
            x0 = x_L4 + dx
            y0 = y_L4
            vx0 = 0.0
            vy0 = 0.0

            sol = _integrate_corotating(x0, y0, vx0, vy0, q, Omega, beta,
                                        T_FINAL_POINCARE, N_OUTPUT)

            # Poincare section at y = y_L4
            x_cross, vx_cross = _poincare_section(sol, y_L4)

            n_cross = len(x_cross)
            if n_cross > 0:
                dx_range = np.max(x_cross) - np.min(x_cross) if n_cross > 1 else 0.0
                dvx_range = np.max(vx_cross) - np.min(vx_cross) if n_cross > 1 else 0.0
            else:
                dx_range = dvx_range = 0.0

            print(f"    orbit {i+1}: dx0={dx:.3f}, "
                  f"{n_cross} crossings, "
                  f"x-spread={dx_range:.4f}, vx-spread={dvx_range:.4f}")

            all_crossings_x.extend(x_cross.tolist())
            all_crossings_vx.extend(vx_cross.tolist())

        results[label] = {
            'x': np.array(all_crossings_x),
            'vx': np.array(all_crossings_vx),
            'n_total': len(all_crossings_x),
        }

        # Regularity metric: compute average nearest-neighbor distance
        # in Poincare section — smaller = more regular
        if len(all_crossings_x) > 2:
            pts = np.column_stack([all_crossings_x, all_crossings_vx])
            # Normalize
            x_scale = max(np.std(pts[:, 0]), 1e-10)
            vx_scale = max(np.std(pts[:, 1]), 1e-10)
            pts_norm = pts / np.array([x_scale, vx_scale])

            # Simple regularity check: variance of nearest-neighbor distances
            from scipy.spatial import cKDTree
            tree = cKDTree(pts_norm)
            dists, _ = tree.query(pts_norm, k=2)
            nn_dists = dists[:, 1]  # nearest neighbor (not self)
            nn_std = np.std(nn_dists)
            nn_mean = np.mean(nn_dists)
            regularity = nn_std / max(nn_mean, 1e-10)
            print(f"    Regularity metric (nn_std/nn_mean): {regularity:.4f}")
            print(f"    (< 0.5 = regular, > 1.0 = chaotic)")

    return results


# -- Part 4: Effective potential landscape ---------------------------------

def part4_potential_landscape():
    """Compare effective potential landscape around L4."""
    print(f"\n{'=' * 78}")
    print("Part 4: Effective Potential Landscape around L4")
    print(f"{'=' * 78}")

    q = 0.01
    Omega = np.sqrt(G_NEWTON)
    L4 = _lagrange_points_newton(q)['L4']
    x_L4, y_L4 = L4

    # Grid around L4
    dx_range = 0.3
    x_arr = np.linspace(x_L4 - dx_range, x_L4 + dx_range, N_GRID)
    y_arr = np.linspace(y_L4 - dx_range, y_L4 + dx_range, N_GRID)
    X, Y = np.meshgrid(x_arr, y_arr)

    results = {}
    for beta in [0.0, BETAS[-1]]:
        label = f"Newton" if beta == 0.0 else f"TGP(beta={beta})"

        if beta == 0.0:
            Z = np.array([[_effective_potential_newton(x, y, q, Omega)
                           for x in x_arr] for y in y_arr])
        else:
            Z = np.array([[_effective_potential_tgp(x, y, q, Omega, beta)
                           for x in x_arr] for y in y_arr])

        V_L4 = Z[N_GRID//2, N_GRID//2]
        Z_rel = Z - V_L4

        # Curvature at L4
        i_c, j_c = N_GRID//2, N_GRID//2
        dx = x_arr[1] - x_arr[0]
        dy = y_arr[1] - y_arr[0]
        d2V_dx2 = (Z[i_c, j_c+1] - 2*Z[i_c, j_c] + Z[i_c, j_c-1]) / dx**2
        d2V_dy2 = (Z[i_c+1, j_c] - 2*Z[i_c, j_c] + Z[i_c-1, j_c]) / dy**2

        print(f"\n  {label}:")
        print(f"    V(L4) = {V_L4:.8e}")
        print(f"    d2V/dx2 at L4 = {d2V_dx2:.6e}")
        print(f"    d2V/dy2 at L4 = {d2V_dy2:.6e}")
        print(f"    V range around L4: [{np.min(Z_rel):.4e}, {np.max(Z_rel):.4e}]")

        # Depth of potential well/saddle
        results[label] = {
            'V_L4': V_L4,
            'd2V_dx2': d2V_dx2,
            'd2V_dy2': d2V_dy2,
            'V_min_rel': np.min(Z_rel),
            'V_max_rel': np.max(Z_rel),
        }

    # Compare curvatures
    if len(results) == 2:
        labels = list(results.keys())
        r0, r1 = results[labels[0]], results[labels[1]]
        print(f"\n  Curvature change (TGP vs Newton):")
        print(f"    d(d2V/dx2) = {r1['d2V_dx2'] - r0['d2V_dx2']:.6e}")
        print(f"    d(d2V/dy2) = {r1['d2V_dy2'] - r0['d2V_dy2']:.6e}")

    return results


def main():
    print("=" * 78)
    print("ex204 -- P3.B: Lagrange Points & Poincare Sections -- TGP vs Newton")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t_total = time.time()

    p1, q_test = part1_lagrange_points()
    p2 = part2_stability()
    p3 = part3_poincare()
    p4 = part4_potential_landscape()

    # Summary
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")

    # L4 shift
    if len(p1) > 1:
        max_shift = max(np.sqrt(r['dx']**2 + r['dy']**2) for r in p1 if r['beta'] > 0)
        print(f"  Max L4 shift (q={q_test}): {max_shift:.4e}")

    # Stability
    newton_stable = [r for r in p2 if r['beta'] == 0.0 and r['stable']]
    if newton_stable:
        q_max_N = max(r['q'] for r in newton_stable)
        print(f"  Newton: L4 stable up to q = {q_max_N:.4f}")

    beta_max = max(b for b in BETAS if b > 0)
    tgp_stable = [r for r in p2 if r['beta'] == beta_max and r['stable']]
    if tgp_stable:
        q_max_T = max(r['q'] for r in tgp_stable)
        print(f"  TGP (beta={beta_max}): L4 stable up to q = {q_max_T:.4f}")

    print(f"\n  Total time: {time.time() - t_total:.1f}s")


if __name__ == "__main__":
    main()
