#!/usr/bin/env python3
"""
collapse_exponent_multi_v47b.py -- Verify gamma_lim = 1/(2*alpha+3) conjecture.

Previous finding (collapse_exponent_v47b.py):
  For alpha=2, d=3: local gamma -> 1/7 = 1/(2*2+3) as delta -> 0.

Test: does gamma_lim = 1/(2*alpha+3) hold for other alpha values?
Also test d-dependence.
"""
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2


def make_solver(alpha, d=3, r_max=300):
    def solver(g0):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-10)
            source = g**2 * (1.0 - g)
            cross = (alpha / g) * gp**2
            if d == 1:
                return [gp, source - cross]
            if r < 1e-10:
                return [gp, (source - cross) / float(d)]
            return [gp, source - cross - float(d - 1) * gp / r]
        sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                        rtol=1e-12, atol=1e-14, max_step=0.02)
        return sol.t, sol.y[0]
    return solver


def A_tail(solver, g0, r_min=50, r_max=200):
    r, g = solver(g0)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


def measure_gamma_fine(alpha, d=3, n_pts=35):
    """Measure limiting gamma using finest delta resolution."""
    gc = (2*alpha + 4) / (2*alpha + 4 - d)
    solver = make_solver(alpha, d)

    # Use very fine deltas near collapse
    deltas = np.logspace(-3.5, -1.0, n_pts)
    g0_vals = gc - deltas

    prev_ld, prev_lA = None, None
    local_gammas = []

    for i, g0 in enumerate(g0_vals):
        if g0 <= 1.001:
            continue
        A = A_tail(solver, g0)
        if A > 1e-8:
            ld = np.log(deltas[i])
            lA = np.log(A)
            if prev_ld is not None:
                lg = -(lA - prev_lA) / (ld - prev_ld)
                local_gammas.append((deltas[i], lg))
            prev_ld, prev_lA = ld, lA

    if len(local_gammas) < 3:
        return None, []

    # Extrapolate: gamma at smallest deltas
    finest = [(d_val, g_val) for d_val, g_val in local_gammas if d_val < 0.003]
    if len(finest) >= 2:
        gamma_lim = np.mean([g for _, g in finest])
    else:
        gamma_lim = local_gammas[0][1]  # smallest delta available

    return gamma_lim, local_gammas


# ================================================================
print("=" * 70)
print("COLLAPSE EXPONENT CONJECTURE: gamma_lim = 1/(2*alpha + 3)")
print("=" * 70)

# Test 1: varying alpha, d=3
print("\n1. gamma_lim vs alpha (d=3)")
print("-" * 50)
print(f"  {'alpha':>6s} {'gc':>8s} {'gamma_lim':>10s} {'1/(2a+3)':>10s} {'ratio':>8s} {'err':>8s}")

alpha_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

for alpha in alpha_vals:
    gc = (2*alpha + 4) / (2*alpha + 1)
    predicted = 1.0 / (2*alpha + 3)
    gamma_lim, all_gammas = measure_gamma_fine(alpha, d=3)
    if gamma_lim is not None:
        ratio = gamma_lim / predicted
        err = abs(gamma_lim - predicted) / predicted
        print(f"  {alpha:6.1f} {gc:8.4f} {gamma_lim:10.6f} {predicted:10.6f} "
              f"{ratio:8.4f} {err:8.4f}")

        # Show trend of local gamma
        if len(all_gammas) > 0:
            ds = [d for d, g in all_gammas]
            gs = [g for d, g in all_gammas]
            print(f"         range: delta [{min(ds):.5f}, {max(ds):.5f}], "
                  f"gamma [{min(gs):.4f}, {max(gs):.4f}]")
    else:
        print(f"  {alpha:6.1f} {gc:8.4f} {'FAILED':>10s} {predicted:10.6f}")


# Test 2: varying d, alpha=2
print("\n\n2. gamma_lim vs d (alpha=2)")
print("-" * 50)
print(f"  {'d':>3s} {'gc':>8s} {'gamma_lim':>10s} {'1/(2a+3)':>10s} {'d/(N*(N-1))':>12s}")

for d in [2, 3, 4, 5]:
    N = 8  # 2*2+4
    if d >= N:
        continue
    gc = N / (N - d)
    pred_1 = 1.0 / (2*2 + 3)  # 1/7 (alpha-only)
    pred_2 = d / (N * (N - 1))  # d-dependent guess
    gamma_lim, all_gammas = measure_gamma_fine(2.0, d=d)
    if gamma_lim is not None:
        print(f"  {d:3d} {gc:8.4f} {gamma_lim:10.6f} {pred_1:10.6f} {pred_2:12.6f}")
    else:
        print(f"  {d:3d} {gc:8.4f} {'FAILED':>10s}")


# Test 3: Check if gamma_lim depends on N=2a+4 or a alone
print("\n\n3. Pattern analysis")
print("-" * 50)
print("  Does gamma_lim = 1/(N-1) = 1/(2a+3)?")
print("  Or gamma_lim = 1/(N-d) = 1/(2a+1) for d=3?")
print("  Or something else?")

print(f"\n  {'alpha':>6s} {'N':>4s} {'gamma_lim':>10s} {'1/(N-1)':>10s} {'1/(N-d)':>10s} {'d/N^2':>10s}")
for alpha in alpha_vals:
    N = 2*alpha + 4
    gc = N / (N - 3)
    gamma_lim, _ = measure_gamma_fine(alpha, d=3)
    if gamma_lim is not None:
        pred_Nm1 = 1/(N-1)
        pred_Nd = 1/(N-3)
        pred_dN2 = 3/N**2
        print(f"  {alpha:6.1f} {N:4.0f} {gamma_lim:10.6f} {pred_Nm1:10.6f} "
              f"{pred_Nd:10.6f} {pred_dN2:10.6f}")


# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
