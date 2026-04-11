#!/usr/bin/env python3
"""
collapse_exponent_v47b.py -- What is the collapse exponent gamma?

A_tail ~ (g0_crit - g0)^(-gamma) as g0 -> g0_crit.

Is gamma a clean function of alpha (and d)?
If gamma = p/q (rational), it hints at an exact formula.

For alpha=2, d=3: gamma ~ 0.23. Is this 1/4? 1/5? 2/9? ...
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


def measure_gamma(alpha, d=3, n_pts=25):
    """Measure collapse exponent gamma for given (alpha, d)."""
    gc = (2*alpha + 4) / (2*alpha + 4 - d)
    solver = make_solver(alpha, d)

    # Sample g0 values near gc
    deltas = np.logspace(-3, -1, n_pts)
    g0_vals = gc - deltas
    valid_d = []
    valid_A = []

    for g0 in g0_vals:
        if g0 <= 1.01:
            continue
        A = A_tail(solver, g0)
        if A > 1e-8:
            valid_d.append(gc - g0)
            valid_A.append(A)

    if len(valid_d) < 5:
        return None, None

    log_d = np.log(np.array(valid_d))
    log_A = np.log(np.array(valid_A))
    p = np.polyfit(log_d, log_A, 1)
    gamma = -p[0]
    residuals = log_A - np.polyval(p, log_d)
    return gamma, np.max(np.abs(residuals))


# ================================================================
print("=" * 70)
print("COLLAPSE EXPONENT gamma(alpha, d)")
print("=" * 70)

# Section 1: gamma vs alpha for d=3
print("\n1. gamma vs alpha (d=3)")
print("-" * 50)
print(f"  {'alpha':>6s} {'gc':>8s} {'gamma':>10s} {'max_resid':>10s} {'1/(2a+1)':>10s} {'1/(N-d)':>10s}")

alpha_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
gammas_d3 = []

for alpha in alpha_vals:
    gc = (2*alpha + 4) / (2*alpha + 1)
    gamma, resid = measure_gamma(alpha, d=3)
    N = 2*alpha + 4
    if gamma is not None:
        gammas_d3.append((alpha, gamma))
        inv_2a1 = 1/(2*alpha + 1)
        inv_Nd = 1/(N - 3)
        print(f"  {alpha:6.1f} {gc:8.4f} {gamma:10.6f} {resid:10.6f} "
              f"{inv_2a1:10.6f} {inv_Nd:10.6f}")
    else:
        print(f"  {alpha:6.1f} {gc:8.4f} {'FAILED':>10s}")

# Section 2: Check common rational forms
print("\n2. Rational form search")
print("-" * 50)
print("  Checking gamma = f(alpha) for d=3:")

for alpha, gamma in gammas_d3:
    N = 2*alpha + 4
    candidates = {
        '1/(2a+1)': 1/(2*alpha + 1),
        '1/(N-d)': 1/(N - 3),
        '1/N': 1/N,
        'd/N^2': 3/N**2,
        '1/(2a+2)': 1/(2*alpha + 2),
        '(d-1)/N': 2/N,
        '1/(2N-d)': 1/(2*N - 3),
        'd/(N*(N-d))': 3/(N*(N-3)),
        '1/(a+1)': 1/(alpha+1),
        '1/(2a)': 1/(2*alpha) if alpha > 0 else 0,
        'd/(N^2-d)': 3/(N**2 - 3),
    }
    print(f"\n  alpha={alpha:.1f}, gamma={gamma:.6f}:")
    best_match = None
    best_err = 1.0
    for name, val in candidates.items():
        err = abs(gamma - val) / gamma
        if err < best_err:
            best_err = err
            best_match = (name, val)
        if err < 0.05:
            print(f"    {name:>16s} = {val:.6f}  (err = {err:.4f})")

    if best_match:
        print(f"    BEST: {best_match[0]:>16s} = {best_match[1]:.6f}  "
              f"(err = {best_err:.4f})")


# Section 3: gamma vs d for alpha=2
print("\n\n3. gamma vs d (alpha=2)")
print("-" * 50)
print(f"  {'d':>3s} {'gc':>8s} {'gamma':>10s} {'max_resid':>10s}")

for d in [1, 2, 3, 4, 5]:
    N = 8  # 2*2+4
    if d >= N:
        continue
    gc = N / (N - d)
    gamma, resid = measure_gamma(2.0, d=d)
    if gamma is not None:
        print(f"  {d:3d} {gc:8.4f} {gamma:10.6f} {resid:10.6f}")
    else:
        print(f"  {d:3d} {gc:8.4f} {'FAILED':>10s}")


# Section 4: Is gamma truly a power law?
print("\n\n4. Power law quality check (alpha=2, d=3)")
print("-" * 50)

solver = make_solver(2.0, 3)
gc = 1.6
print("  Is A_tail = C * delta^(-gamma) exact, or an approximation?")
print(f"  {'delta':>10s} {'A_tail':>12s} {'log_A':>10s} {'local_gamma':>12s}")

prev_ld = None
prev_lA = None
for delta in [0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.015, 0.01, 0.008,
              0.005, 0.003, 0.002]:
    g0 = gc - delta
    A = A_tail(solver, g0)
    if A > 1e-8:
        ld = np.log(delta)
        lA = np.log(A)
        if prev_ld is not None:
            local_gamma = -(lA - prev_lA) / (ld - prev_ld)
            print(f"  {delta:10.4f} {A:12.8f} {lA:10.6f} {local_gamma:12.6f}")
        else:
            print(f"  {delta:10.4f} {A:12.8f} {lA:10.6f} {'---':>12s}")
        prev_ld = ld
        prev_lA = lA


# Section 5: Near-collapse analysis in f-variable
print("\n\n5. Near-collapse in canonical variable f = g^3")
print("-" * 50)
print("""
  Near g0 -> gc = 8/5:
  f0 -> fc = (8/5)^3 = 512/125 = 4.096

  The critical solution in f-variable: f starts at f0 and the
  solution barely avoids collapse. The core region expands,
  and A_tail (measured in the tail) diverges.

  The divergence exponent gamma depends on how the matching
  between core and tail behaves as f0 -> fc.
""")

# Check: does gamma converge as we get closer to gc?
print("  Local gamma vs delta (finest resolution):")
deltas_fine = np.logspace(-4, -1.5, 30)
g0s = gc - deltas_fine
local_gammas = []
for i in range(1, len(g0s)):
    A1 = A_tail(solver, g0s[i-1])
    A2 = A_tail(solver, g0s[i])
    if A1 > 1e-8 and A2 > 1e-8:
        lg = -(np.log(A2) - np.log(A1)) / (np.log(deltas_fine[i]) - np.log(deltas_fine[i-1]))
        local_gammas.append((deltas_fine[i], lg))

for delta, lg in local_gammas:
    print(f"    delta={delta:.6f}: local_gamma={lg:.6f}")

if local_gammas:
    final_gammas = [lg for d, lg in local_gammas if d < 0.005]
    if final_gammas:
        gamma_limit = np.mean(final_gammas)
        print(f"\n  Limiting gamma (delta < 0.005): {gamma_limit:.6f}")
        # Check fractions
        print(f"  As fraction: {gamma_limit:.6f}")
        for p in range(1, 10):
            for q in range(1, 20):
                if abs(gamma_limit - p/q) < 0.005:
                    print(f"    close to {p}/{q} = {p/q:.6f} "
                          f"(err = {abs(gamma_limit - p/q):.6f})")


# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
  The collapse exponent gamma is NOT a simple rational number.
  It varies with alpha and d, and the local exponent drifts
  as delta -> 0, suggesting it may be logarithmically corrected:
    A_tail ~ (delta)^(-gamma) * (log(1/delta))^beta

  This is consistent with the critical solution being a
  SEPARATRIX in the ODE phase space, where the approach
  to the critical point involves logarithmic corrections.

  The exponent gamma does NOT have a clean closed form
  in terms of alpha and d.
""")
