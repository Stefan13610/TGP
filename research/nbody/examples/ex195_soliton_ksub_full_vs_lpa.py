#!/usr/bin/env python3
"""
ex195 — Soliton profile: K_sub = g^2 (full) vs f = 1+4*ln(g) (LPA)
====================================================================

Solves the radial soliton ODE in both kinetic coupling forms and
compares the resulting profiles, energies, and tail amplitudes.

FULL form (manuscript, sek10_N0_wyprowadzenie.tex):
  K_sub(g) = g^2
  ODE: g^2*g'' + g*(g')^2 + (2/r)*g^2*g' = g^2*(1 - g)

  Rewritten as first-order system:
    g' = u
    u' = (1 - g) - u^2/g - 2*u/r

LPA form (one-loop ERG truncation, legacy nbody):
  f(g) = 1 + 4*ln(g)
  ODE: f(g)*g'' + (2/r)*f(g)*g' + (2/g)*(g')^2 = V'(g)

  where V'(g) = g^2 - g^3 (for beta = gamma = 1).
  Rewritten:
    g' = u
    u' = [V'(g) - (2/g)*u^2 - (2/r)*f(g)*u] / f(g)

Both solved with boundary conditions: g(0) = g_0, g'(0) = 0, g(inf) -> 1.
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp

# ── Parameters ──────────────────────────────────────────────────────────────
G0_VALUES = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.869, 0.95]
R_MAX = 60.0
R_EVAL_N = 6000
ALPHA = 2.0

quick = "--quick" in sys.argv
if quick:
    G0_VALUES = [0.1, 0.5, 0.869]
    R_MAX = 40.0
    R_EVAL_N = 2000


# ── ODE: K_sub = g^2 (FULL) ────────────────────────────────────────────────

def ode_full(r, y):
    """g^2*g'' + g*(g')^2 + (2/r)*g^2*g' = g^2*(1-g)
    => g' = u
    => u' = (1-g) - u^2/g - 2*u/r
    """
    g, u = y
    g = max(g, 1e-30)  # prevent division by zero

    if r < 1e-12:
        # L'Hopital for 2*u/r at r=0: limit is 2*u'(0)
        # g''(0) = (1-g) - u^2/g - 2*g''(0)
        # 3*g''(0) = (1-g) - u^2/g
        # But u(0) = 0, so g''(0) = (1-g)/3
        up = (1.0 - g) / 3.0
    else:
        up = (1.0 - g) - u**2 / g - 2.0 * u / r
    return [u, up]


# ── ODE: f = 1 + 4*ln(g) (LPA) ─────────────────────────────────────────────

def ode_lpa(r, y):
    """f(g)*g'' + (2/r)*f(g)*g' + (alpha/g)*(g')^2 = V'(g)
    where f(g) = 1 + 2*alpha*ln(g), alpha=2 => f = 1 + 4*ln(g)
    V'(g) = g^2 - g^3  (for beta=gamma=1)

    => g' = u
    => u' = [V'(g) - (alpha/g)*u^2 - (2/r)*f(g)*u] / f(g)
    """
    g, u = y
    g = max(g, 1e-30)

    f_g = 1.0 + 2.0 * ALPHA * np.log(g)
    Vp = g**2 - g**3

    if abs(f_g) < 1e-15:
        # f(g) = 0 is a singularity of LPA form
        return [u, 0.0]

    if r < 1e-12:
        # At r=0, u=0, so only V'(g)/(3*f(g)) survives (same L'Hopital as full)
        up = Vp / (3.0 * f_g)
    else:
        up = (Vp - (ALPHA / g) * u**2 - (2.0 / r) * f_g * u) / f_g
    return [u, up]


# ── Solver ──────────────────────────────────────────────────────────────────

def solve_soliton(g0, ode_func, r_max=R_MAX, n_eval=R_EVAL_N):
    """Solve soliton ODE from r=eps to r=r_max."""
    r_start = 1e-6
    r_eval = np.linspace(r_start, r_max, n_eval)
    y0 = [g0, 0.0]

    sol = solve_ivp(ode_func, [r_start, r_max], y0, t_eval=r_eval,
                    method='DOP853', rtol=1e-11, atol=1e-13,
                    max_step=0.05)

    if not sol.success:
        return r_eval[:len(sol.t)], sol.y[0], sol.y[1], False
    return sol.t, sol.y[0], sol.y[1], True


def tail_amplitude(r, g, r_min=15.0, r_max=None):
    """Extract tail amplitude A from g(r) ~ 1 + A*sin(r+delta)/r.

    For the full form, the oscillatory tail has period 2*pi.
    Returns max|g-1|*r in the tail region.
    """
    if r_max is None:
        r_max = r[-1]
    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 10:
        return np.nan
    deviation = np.abs(g[mask] - 1.0) * r[mask]
    return np.max(deviation)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("ex195 -- Soliton profile: K_sub = g^2 (FULL) vs f = 1+4*ln(g) (LPA)")
    print("=" * 72)

    print(f"\n{'g_0':>8s} | {'FULL':^28s} | {'LPA':^28s} | {'Ratio':^12s}")
    print(f"{'':>8s} | {'g(rmax)':>10s} {'A_tail':>8s} {'OK':>6s} | "
          f"{'g(rmax)':>10s} {'A_tail':>8s} {'OK':>6s} | {'A_F/A_L':>10s}")
    print("-" * 98)

    results = []

    for g0 in G0_VALUES:
        # Solve FULL
        r_f, g_f, gp_f, ok_f = solve_soliton(g0, ode_full)
        A_f = tail_amplitude(r_f, g_f) if ok_f else np.nan
        gend_f = g_f[-1] if ok_f else np.nan

        # Solve LPA
        r_l, g_l, gp_l, ok_l = solve_soliton(g0, ode_lpa)
        A_l = tail_amplitude(r_l, g_l) if ok_l else np.nan
        gend_l = g_l[-1] if ok_l else np.nan

        ratio = A_f / A_l if (A_l > 0 and np.isfinite(A_f)) else np.nan

        ok_str_f = "PASS" if ok_f else "FAIL"
        ok_str_l = "PASS" if ok_l else "FAIL"

        print(f"{g0:8.3f} | {gend_f:10.6f} {A_f:8.4f} {ok_str_f:>6s} | "
              f"{gend_l:10.6f} {A_l:8.4f} {ok_str_l:>6s} | {ratio:10.4f}")

        results.append({
            'g0': g0, 'ok_full': ok_f, 'ok_lpa': ok_l,
            'A_full': A_f, 'A_lpa': A_l, 'ratio': ratio,
            'gend_full': gend_f, 'gend_lpa': gend_l
        })

    # ── Key comparison: g_0 = 0.869 (the TGP particle sector calibration) ──
    print("\n" + "=" * 72)
    print("KEY COMPARISON: g_0 = 0.869 (phi-FP calibration point)")
    print("=" * 72)

    key = [r for r in results if abs(r['g0'] - 0.869) < 0.01]
    if key:
        k = key[0]
        print(f"  FULL: A_tail = {k['A_full']:.6f}, g(r_max) = {k['gend_full']:.8f}")
        print(f"  LPA:  A_tail = {k['A_lpa']:.6f}, g(r_max) = {k['gend_lpa']:.8f}")
        print(f"  Ratio A_full/A_lpa = {k['ratio']:.6f}")
        pct = abs(k['ratio'] - 1) * 100
        print(f"  Difference: {pct:.2f}%")

        if pct < 1.0:
            print("  -> SMALL difference: LPA is good approximation at this g_0")
        elif pct < 10.0:
            print("  -> MODERATE difference: LPA corrections non-negligible")
        else:
            print("  -> LARGE difference: LPA form is NOT reliable at this g_0")

    # ── Ghost check ──
    print("\n--- Ghost check (f(g) < 0 in LPA) ---")
    g_ghost = np.exp(-1.0 / (2.0 * ALPHA))  # f(g) = 0 when ln(g) = -1/(2*alpha)
    print(f"  LPA f(g) = 0 at g = exp(-1/(2*alpha)) = {g_ghost:.6f}")
    print(f"  For g < {g_ghost:.4f}, LPA has kinetic ghosts (f < 0)")
    print(f"  FULL K_sub = g^2 > 0 for ALL g > 0 (ghost-free)")

    ghost_affected = [r for r in results if r['g0'] < g_ghost]
    if ghost_affected:
        print(f"  WARNING: g_0 values {[r['g0'] for r in ghost_affected]} "
              f"are in the ghost regime of LPA!")
    else:
        print(f"  All tested g_0 values are above ghost threshold.")

    # ── Summary ──
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)

    n_full_ok = sum(1 for r in results if r['ok_full'])
    n_lpa_ok = sum(1 for r in results if r['ok_lpa'])

    print(f"  FULL form: {n_full_ok}/{len(results)} converged")
    print(f"  LPA form:  {n_lpa_ok}/{len(results)} converged")
    print(f"  Ghost threshold (LPA): g = {g_ghost:.4f}")
    print(f"  Manuscript uses FULL form (K_sub = g^2, ghost-free)")

    # Check overall agreement in weak-field regime
    weak_field = [r for r in results if r['g0'] > 0.5
                  and r['ok_full'] and r['ok_lpa']
                  and np.isfinite(r['ratio'])]
    if weak_field:
        max_diff = max(abs(r['ratio'] - 1) * 100 for r in weak_field)
        print(f"  Max A_tail difference (g_0 > 0.5): {max_diff:.2f}%")
        if max_diff < 5.0:
            print(f"  -> Weak-field regime: LPA ≈ FULL (OK)")
        else:
            print(f"  -> Even in weak field: significant LPA vs FULL deviation!")

    strong_field = [r for r in results if r['g0'] < 0.3
                    and r['ok_full'] and r['ok_lpa']
                    and np.isfinite(r['ratio'])]
    if strong_field:
        max_diff_sf = max(abs(r['ratio'] - 1) * 100 for r in strong_field)
        print(f"  Max A_tail difference (g_0 < 0.3): {max_diff_sf:.2f}%")


if __name__ == "__main__":
    main()
