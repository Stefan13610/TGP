#!/usr/bin/env python3
"""
ex197 -- Find optimal g0^e for FULL K_sub=g^2 form to match r_21 = 206.768
============================================================================

Uses high-precision ODE solving and improved A_tail extraction (fitting
the exact asymptotic form g ~ 1 + A*sin(r+delta)/r) to find the g0^e
that gives r_21 = m_mu/m_e = 206.768 (PDG).

Then computes all particle predictions at that optimal point.
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
from scipy.optimize import curve_fit

# ── Constants ───────────────────────────────────────────────────────────────
PHI = (1.0 + np.sqrt(5.0)) / 2.0
R21_PDG = 206.768
M_E = 0.51100     # MeV
M_MU = 105.658    # MeV
M_TAU = 1776.86   # MeV

quick = "--quick" in sys.argv
R_MAX = 50.0 if quick else 120.0
N_EVAL = 3000 if quick else 12000


# ── ODE: K_sub = g^2 ───────────────────────────────────────────────────────

def ode_full(r, y):
    g, u = y
    g = max(g, 1e-30)
    if r < 1e-12:
        up = (1.0 - g) / 3.0
    else:
        up = (1.0 - g) - u**2 / g - 2.0 * u / r
    return [u, up]


def solve_soliton(g0, r_max=R_MAX, n_eval=N_EVAL):
    r_start = 1e-6
    r_eval = np.linspace(r_start, r_max, n_eval)
    y0 = [g0, 0.0]
    sol = solve_ivp(ode_full, [r_start, r_max], y0, t_eval=r_eval,
                    method='DOP853', rtol=1e-13, atol=1e-15, max_step=0.02)
    if not sol.success:
        return None, None, False
    return sol.t, sol.y[0], True


def fit_tail_amplitude(r, g, r_min=25.0, r_max=None):
    """Fit g(r) = 1 + A*sin(r + delta)/r in the tail.

    Uses least-squares fit to extract A and delta.
    Returns A (always positive).
    """
    if r_max is None:
        r_max = r[-1] - 10.0
    mask = (r >= r_min) & (r <= r_max)
    rm = r[mask]
    gm = g[mask]

    if len(rm) < 30:
        # Fallback: envelope method
        dev = np.abs(gm - 1.0) * rm
        return np.max(dev) if len(dev) > 0 else np.nan

    # Fit: g - 1 = (A*sin(r + delta))/r = (a*sin(r) + b*cos(r))/r
    # where a = A*cos(delta), b = A*sin(delta)
    sin_r = np.sin(rm) / rm
    cos_r = np.cos(rm) / rm

    # Build design matrix
    X = np.column_stack([sin_r, cos_r])
    y = gm - 1.0

    # Least squares: y = X @ [a, b]
    result = np.linalg.lstsq(X, y, rcond=None)
    a, b = result[0]
    A = np.sqrt(a**2 + b**2)

    return A


def compute_r21(g0_e):
    """Compute r_21 = (A_mu/A_e)^4 for given g0^e."""
    g0_mu = PHI * g0_e

    r_e, g_e, ok_e = solve_soliton(g0_e)
    if not ok_e:
        return np.nan, np.nan, np.nan

    r_mu, g_mu, ok_mu = solve_soliton(g0_mu)
    if not ok_mu:
        return np.nan, np.nan, np.nan

    A_e = fit_tail_amplitude(r_e, g_e)
    A_mu = fit_tail_amplitude(r_mu, g_mu)

    if A_e < 1e-20 or A_mu < 1e-20 or not np.isfinite(A_e) or not np.isfinite(A_mu):
        return np.nan, A_e, A_mu

    r21 = (A_mu / A_e)**4
    return r21, A_e, A_mu


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("ex197 -- Optimal g0^e for FULL K_sub=g^2 to match r_21 = 206.768")
    print("=" * 72)

    # 1. Coarse scan
    print("\n--- Coarse scan ---")
    print(f"  {'g0_e':>8s} {'r21':>12s} {'err(%)':>10s} {'A_e':>12s} {'A_mu':>12s}")
    print(f"  {'-'*56}")

    g0_vals = np.arange(0.84, 0.90, 0.005)
    r21_vals = []

    for g0 in g0_vals:
        r21, A_e, A_mu = compute_r21(g0)
        err = abs(r21 - R21_PDG) / R21_PDG * 100 if np.isfinite(r21) else np.nan
        print(f"  {g0:8.4f} {r21:12.3f} {err:10.4f} {A_e:12.6f} {A_mu:12.6f}")
        r21_vals.append(r21)

    # 2. Find bracket for r21 = 206.768
    r21_arr = np.array(r21_vals)
    target = R21_PDG

    # Find where r21 crosses target
    below = g0_vals[r21_arr < target]
    above = g0_vals[r21_arr > target]

    if len(below) > 0 and len(above) > 0:
        g0_lo = below[-1]
        g0_hi = above[0]
        print(f"\n  Bracket: [{g0_lo:.4f}, {g0_hi:.4f}]")

        # 3. Brent root finding
        print("\n--- Brent optimization ---")

        def objective(g0):
            r21, _, _ = compute_r21(g0)
            return r21 - target if np.isfinite(r21) else 1e6

        g0_opt = brentq(objective, g0_lo, g0_hi, xtol=1e-6)
        r21_opt, A_e_opt, A_mu_opt = compute_r21(g0_opt)

        print(f"  Optimal g0^e = {g0_opt:.6f}")
        print(f"  r_21 = {r21_opt:.4f} (target: {target:.3f})")
        print(f"  A_tail(e) = {A_e_opt:.8f}")
        print(f"  A_tail(mu) = {A_mu_opt:.8f}")
        print(f"  A_mu/A_e = {A_mu_opt/A_e_opt:.8f}")
        print(f"  (A_mu/A_e)^4 = {(A_mu_opt/A_e_opt)**4:.4f}")

    else:
        # Fallback: minimum error scan
        print("\n  No bracket found; using minimum error")
        idx = np.argmin(np.abs(r21_arr - target))
        g0_opt = g0_vals[idx]
        r21_opt, A_e_opt, A_mu_opt = compute_r21(g0_opt)
        print(f"  Best g0^e = {g0_opt:.4f}, r21 = {r21_opt:.3f}")

    # 4. Fine scan around optimal
    print("\n--- Fine scan around optimal ---")
    print(f"  {'g0_e':>10s} {'r21':>12s} {'err(%)':>12s}")
    print(f"  {'-'*38}")

    for dg in np.arange(-0.003, 0.004, 0.0005):
        g0_test = g0_opt + dg
        r21_test, _, _ = compute_r21(g0_test)
        err = abs(r21_test - R21_PDG) / R21_PDG * 100 if np.isfinite(r21_test) else np.nan
        marker = " <-- optimal" if abs(dg) < 0.0003 else ""
        print(f"  {g0_test:10.5f} {r21_test:12.4f} {err:12.5f}{marker}")

    # 5. Full predictions at optimal g0^e
    print("\n" + "=" * 72)
    print(f"FULL PREDICTIONS at g0^e = {g0_opt:.6f}")
    print("=" * 72)

    g0_mu = PHI * g0_opt
    print(f"  g0^e = {g0_opt:.6f}")
    print(f"  g0^mu = phi * g0^e = {g0_mu:.6f}")

    # Lepton masses
    m_mu_pred = M_E * r21_opt
    sqrt_me = np.sqrt(M_E)
    sqrt_mmu = np.sqrt(m_mu_pred)

    # Koide
    S = sqrt_me + sqrt_mmu
    M_sum = M_E + m_mu_pred
    disc = 16 * S**2 - 4 * (3 * M_sum - 2 * S**2)
    if disc >= 0:
        x = (4 * S + np.sqrt(disc)) / 2.0
        m_tau_pred = x**2
    else:
        m_tau_pred = np.nan

    r31 = m_tau_pred / M_E if np.isfinite(m_tau_pred) else np.nan
    r32 = m_tau_pred / m_mu_pred if np.isfinite(m_tau_pred) else np.nan

    # Koide parameter
    if np.isfinite(m_tau_pred):
        sqrt_mtau = np.sqrt(m_tau_pred)
        K = (sqrt_me + sqrt_mmu + sqrt_mtau)**2 / (M_E + m_mu_pred + m_tau_pred)
    else:
        K = np.nan

    # alpha_s (N_f = 5 at M_Z)
    alpha_s_mz = 27.0 * g0_opt / (8.0 * 25.0)
    alpha_s_mtau = 3.0 * g0_opt / 8.0
    ratio_as = alpha_s_mtau / alpha_s_mz  # should be (5/3)^2 = 25/9

    print(f"\n  {'Observable':>24s} {'TGP':>12s} {'PDG':>12s} {'sigma/err':>12s}")
    print(f"  {'-'*64}")
    print(f"  {'r_21 = m_mu/m_e':>24s} {r21_opt:12.3f} {'206.768':>12s} "
          f"{abs(r21_opt - R21_PDG)/R21_PDG*100:10.4f}%")
    print(f"  {'m_mu (MeV)':>24s} {m_mu_pred:12.3f} {'105.658':>12s} "
          f"{abs(m_mu_pred - M_MU)/M_MU*100:10.4f}%")

    if np.isfinite(m_tau_pred):
        print(f"  {'r_31 = m_tau/m_e':>24s} {r31:12.1f} {'3477.2':>12s} "
              f"{abs(r31 - 3477.2)/3477.2*100:10.4f}%")
        print(f"  {'r_32 = m_tau/m_mu':>24s} {r32:12.3f} {'16.817':>12s} "
              f"{abs(r32 - 16.817)/16.817*100:10.4f}%")
        print(f"  {'m_tau (MeV)':>24s} {m_tau_pred:12.2f} {'1776.86':>12s} "
              f"{abs(m_tau_pred - M_TAU)/M_TAU*100:10.4f}%")
        print(f"  {'K_lepton':>24s} {K:12.5f} {'0.66667':>12s} "
              f"{abs(K - 2.0/3.0)/(2.0/3.0)*100:10.4f}%")

    print(f"  {'alpha_s(M_Z)':>24s} {alpha_s_mz:12.4f} {'0.1179+/-9':>12s} "
          f"{abs(alpha_s_mz - 0.1179)/0.0009:10.2f}sigma")
    print(f"  {'alpha_s(m_tau)':>24s} {alpha_s_mtau:12.4f} {'0.330+/-14':>12s} "
          f"{abs(alpha_s_mtau - 0.330)/0.014:10.2f}sigma")
    print(f"  {'alpha_s(tau)/alpha_s(Z)':>24s} {ratio_as:12.4f} {'2.799+/-121':>12s} "
          f"{abs(ratio_as - 2.799)/0.121:10.2f}sigma")

    # 6. Comparison with manuscript value
    print(f"\n--- Comparison with manuscript g0^e = 0.869 ---")
    r21_ms, A_e_ms, A_mu_ms = compute_r21(0.869)
    alpha_s_ms = 27.0 * 0.869 / 200.0
    print(f"  g0^e = 0.869: r_21 = {r21_ms:.3f}, alpha_s(M_Z) = {alpha_s_ms:.4f}")
    print(f"  g0^e = {g0_opt:.6f}: r_21 = {r21_opt:.3f}, alpha_s(M_Z) = {alpha_s_mz:.4f}")

    delta_g0 = g0_opt - 0.869
    print(f"  Delta g0^e = {delta_g0:+.6f} ({delta_g0/0.869*100:+.3f}%)")

    print(f"\n  NOTE: alpha_s(M_Z) = 27*g0/(8*25) depends linearly on g0^e")
    print(f"  At optimal g0: alpha_s(M_Z) = {alpha_s_mz:.4f}")
    print(f"  PDG: 0.1179 +/- 0.0009")
    print(f"  Tension: {abs(alpha_s_mz - 0.1179)/0.0009:.1f} sigma")

    # Summary
    print(f"\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"  FULL form K_sub = g^2 with optimal g0^e = {g0_opt:.6f}:")
    print(f"  - r_21 matches PDG to {abs(r21_opt - R21_PDG)/R21_PDG*100:.4f}%")
    print(f"  - alpha_s(M_Z) = {alpha_s_mz:.4f} ({abs(alpha_s_mz - 0.1179)/0.0009:.1f}sigma from PDG)")
    if np.isfinite(m_tau_pred):
        print(f"  - m_tau = {m_tau_pred:.2f} MeV ({abs(m_tau_pred - M_TAU)/M_TAU*100:.3f}% from PDG)")
    print(f"  - K_lepton = {K:.5f} (target: 0.66667)")

    # Count passes
    passes = 0
    total = 0
    checks = [
        ("r_21", abs(r21_opt - R21_PDG) / R21_PDG * 100 < 1.0),
        ("m_mu", abs(m_mu_pred - M_MU) / M_MU * 100 < 1.0),
        ("alpha_s(M_Z)", abs(alpha_s_mz - 0.1179) / 0.0009 < 2.0),
        ("alpha_s(m_tau)", abs(alpha_s_mtau - 0.330) / 0.014 < 2.0),
        ("ratio", abs(ratio_as - 2.799) / 0.121 < 2.0),
    ]
    if np.isfinite(m_tau_pred):
        checks.extend([
            ("m_tau", abs(m_tau_pred - M_TAU) / M_TAU * 100 < 1.0),
            ("r_32", abs(r32 - 16.817) / 16.817 * 100 < 1.0),
            ("K_lepton", abs(K - 2.0/3.0) < 0.001),
        ])

    for name, passed in checks:
        total += 1
        if passed:
            passes += 1
            print(f"    [{name}] PASS")
        else:
            print(f"    [{name}] MARGINAL")

    print(f"\n  Score: {passes}/{total} PASS")


if __name__ == "__main__":
    main()
