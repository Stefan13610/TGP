#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex233_Kg2_both_fp_tau.py
===========================
OBA φ-FP DLA K=g² I PEŁNY TEST τ SOLITONU

KONTEKST (ex232):
  K=g² u-solver znalazł g₀* = 0.8695 (= g₀ᵉ), r₂₁ ✓
  ALE τ soliton kolapsa do g=0 (u < 0)
  Skan sugeruje DRUGĄ φ-FP przy g₀ ≈ 1.25 (stary wynik z ex106)

PLAN:
  §1. Pełny skan r₂₁(g₀) dla K=g² — znajdź WSZYSTKIE φ-FP
  §2. Użyj g-solver (z g'²/g osobliwością) zamiast u-solver
  §3. Dla każdej φ-FP: oblicz τ, r₃₁, Koide
  §4. Porównaj obie φ-FP

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2
R21_PDG = 206.768
R31_PDG = 3477.48

TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")

def koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2


# ============================================================
# K=g² g-SPACE SOLVER (with physical g'²/g friction)
# ============================================================

def solve_K2_g(g0, R_max=200.0, N_pts=8000):
    """
    K=g² soliton in g-space:
      g'' = (1-g) - g'²/g - 2g'/r

    The g'²/g term diverges at g=0, PREVENTING collapse.
    This is PHYSICAL: g=0 means metric degeneracy.
    """
    def rhs(r, y):
        g, gp = y
        g_safe = max(g, 1e-12)
        r_safe = max(r, 1e-10)
        gpp = (1.0 - g_safe) - gp**2 / g_safe - 2.0 * gp / r_safe
        return [gp, gpp]

    r0 = 1e-4
    acc0 = 1.0 - g0  # g''(0) = (1-g₀)
    g_init = g0 + acc0 * r0**2 / 6.0
    gp_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.02, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    g_arr = sol.sol(r_arr)[0]
    return r_arr, g_arr, sol.status


def fit_tail(r_arr, g_arr, r_L=20.0, r_R=60.0):
    """Fit (g-1)·r = B cos(r) + C sin(r)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, 0.0
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)

    y_pred = X @ coefs
    ss_res = np.sum((y_fit - y_pred)**2)
    ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return A, np.arctan2(B, C), R2


def get_A(g0, R_max=120.0, r_L=20.0, r_R=60.0):
    """Get A_tail from g-solver."""
    try:
        r, g, status = solve_K2_g(g0, R_max=R_max, N_pts=6000)
        if status < 0:
            return np.nan, 0.0
        A, _, R2 = fit_tail(r, g, r_L=r_L, r_R=r_R)
        return A, R2
    except:
        return np.nan, 0.0


# ============================================================
# §1. PEŁNY SKAN r₂₁(g₀) — FINDING ALL φ-FPs
# ============================================================
print("=" * 72)
print("§1. PEŁNY SKAN r₂₁(g₀) DLA K=g²")
print("=" * 72)

g0_scan = np.concatenate([
    np.linspace(0.50, 0.95, 19),  # below vacuum
    np.linspace(1.05, 2.00, 20),  # above vacuum
])

print(f"\n  {'g₀':>8s}  {'A(g₀)':>12s}  {'A(φg₀)':>12s}  {'r₂₁':>12s}  {'R²':>8s}")
print("  " + "-" * 56)

ratios = []
for g0 in g0_scan:
    A1, R2_1 = get_A(g0)
    A2, R2_2 = get_A(phi * g0)

    if np.isfinite(A1) and A1 > 1e-15 and np.isfinite(A2):
        r21 = (A2/A1)**4
        ratios.append((g0, r21, A1, A2, R2_1))
        if abs(r21) < 1e8:
            print(f"  {g0:8.4f}  {A1:12.6f}  {A2:12.6f}  {r21:12.2f}  {R2_1:8.4f}")
    else:
        ratios.append((g0, np.nan, A1, A2, 0))
        if g0 < 1.6:  # Don't print all NaNs for large g₀
            print(f"  {g0:8.4f}  {A1 if np.isfinite(A1) else 'NaN':>12}  {'NaN':>12s}")

# Find ALL zero crossings of (r₂₁ - R21_PDG)
phi_fps = []
for i in range(len(ratios)-1):
    g1, r1, _, _, _ = ratios[i]
    g2, r2, _, _, _ = ratios[i+1]
    if np.isfinite(r1) and np.isfinite(r2):
        v1 = r1 - R21_PDG
        v2 = r2 - R21_PDG
        if v1 * v2 < 0:
            try:
                def rf(g0):
                    A1, _ = get_A(g0)
                    A2, _ = get_A(phi * g0)
                    if np.isnan(A1) or A1 < 1e-15:
                        return 1e10
                    return (A2/A1)**4 - R21_PDG
                g0_fp = brentq(rf, g1, g2, xtol=1e-7)
                phi_fps.append(g0_fp)
                print(f"\n  ★ φ-FP at g₀* = {g0_fp:.8f}")
            except Exception as e:
                print(f"\n  brentq failed between {g1:.4f} and {g2:.4f}: {e}")

print(f"\n  Znaleziono {len(phi_fps)} φ-FP(s).")
for i, fp in enumerate(phi_fps):
    print(f"    FP{i+1}: g₀* = {fp:.8f}")


# ============================================================
# §2. DLA KAŻDEJ φ-FP: PEŁNY TEST SPEKTRUM
# ============================================================
print("\n" + "=" * 72)
print("§2. TEST KAŻDEJ φ-FP")
print("=" * 72)

for fp_idx, g0_star in enumerate(phi_fps):
    g0_e = g0_star
    g0_mu = phi * g0_star
    g0_tau = phi**2 * g0_star

    print(f"\n  ═══ FP{fp_idx+1}: g₀* = {g0_star:.8f} ═══")
    print(f"    g₀_e  = {g0_e:.6f}")
    print(f"    g₀_μ  = {g0_mu:.6f}")
    print(f"    g₀_τ  = {g0_tau:.6f}")

    # e soliton
    A_e, R2_e = get_A(g0_e)

    # μ soliton
    A_mu, R2_mu = get_A(g0_mu)

    r21 = (A_mu/A_e)**4
    print(f"    A_e = {A_e:.8f}  (R² = {R2_e:.4f})")
    print(f"    A_μ = {A_mu:.8f}  (R² = {R2_mu:.4f})")
    print(f"    r₂₁ = {r21:.4f}")

    # τ soliton — try multiple configs
    print(f"\n    τ soliton (g₀ = {g0_tau:.4f}):")

    tau_results = []
    for R_max in [100, 150, 200]:
        for r_L, r_R in [(20, 55), (25, 70), (30, 85)]:
            A_tau, R2_tau = get_A(g0_tau, R_max=R_max, r_L=r_L, r_R=r_R)
            tau_results.append((A_tau, R2_tau, R_max, r_L, r_R))

    # Also try stiff solver
    try:
        sol_stiff = solve_ivp(
            lambda r, y: [(gp := y[1]),
                          (1.0 - max(y[0], 1e-12)) - gp**2 / max(y[0], 1e-12) - 2.0*gp/max(r, 1e-10)],
            [1e-4, 200.0],
            [g0_tau + (1-g0_tau)*1e-8/6, (1-g0_tau)*1e-4/3],
            method='Radau', rtol=1e-8, atol=1e-10,
            max_step=0.05, dense_output=True
        )
        if sol_stiff.status >= 0:
            r_stiff = np.linspace(sol_stiff.t[0], min(sol_stiff.t[-1], 200.0), 8000)
            g_stiff = sol_stiff.sol(r_stiff)[0]
            for r_L, r_R in [(20, 55), (25, 70)]:
                A_s, _, R2_s = fit_tail(r_stiff, g_stiff, r_L=r_L, r_R=r_R)
                tau_results.append((A_s, R2_s, 200, r_L, r_R))
    except:
        pass

    # Show results
    for A_t, R2_t, Rm, rL, rR in tau_results:
        status = "OK" if np.isfinite(A_t) and A_t > 0 and R2_t > 0.9 else "FAIL"
        if np.isfinite(A_t):
            print(f"      R={Rm} [{rL},{rR}]: A_τ = {A_t:.6f}, R² = {R2_t:.4f}  [{status}]")

    # Pick best
    valid = [(A, R2, Rm, rL, rR) for A, R2, Rm, rL, rR in tau_results
             if np.isfinite(A) and A > 0 and R2 > 0.5]
    if valid:
        best = max(valid, key=lambda x: x[1])
        A_tau, R2_tau = best[0], best[1]
    elif any(np.isfinite(x[0]) and x[0] > 0 for x in tau_results):
        # Take highest R², even if < 0.5
        valid2 = [(A, R2, Rm, rL, rR) for A, R2, Rm, rL, rR in tau_results
                  if np.isfinite(A) and A > 0]
        best = max(valid2, key=lambda x: x[1])
        A_tau, R2_tau = best[0], best[1]
    else:
        A_tau, R2_tau = np.nan, 0.0

    # Diagnostics: show τ profile
    r_diag, g_diag, st = solve_K2_g(g0_tau, R_max=100.0, N_pts=2000)
    g_min = np.min(g_diag)
    g_max = np.max(g_diag)
    g_end = g_diag[-1]
    print(f"\n    τ profile: g ∈ [{g_min:.4f}, {g_max:.4f}], g(end) = {g_end:.4f}")
    print(f"    τ solver status: {st}")

    # Show key radii
    print(f"    {'r':>8s}  {'g_τ':>10s}")
    for ri in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 20.0, 50.0]:
        idx = np.searchsorted(r_diag, ri)
        if idx < len(r_diag):
            print(f"    {ri:8.1f}  {g_diag[idx]:10.6f}")

    # Results
    if np.isfinite(A_tau) and A_tau > 0:
        r31 = (A_tau / A_e)**4
        K_val = koide(1.0, r21, r31)

        print(f"\n    ★ FP{fp_idx+1} RESULTS:")
        print(f"      A_τ = {A_tau:.8f}  (R² = {R2_tau:.4f})")
        print(f"      r₃₁ = {r31:.2f}  (PDG: {R31_PDG})")
        print(f"      Error r₃₁: {abs(r31-R31_PDG)/R31_PDG*100:.2f}%")
        print(f"      Koide = {K_val:.6f}  (2/3 = {2/3:.6f})")

        record(f"T{2*fp_idx+1}: FP{fp_idx+1} r₃₁",
               abs(r31 - R31_PDG)/R31_PDG < 0.10,
               f"r₃₁ = {r31:.2f}, err = {abs(r31-R31_PDG)/R31_PDG*100:.1f}%")

        record(f"T{2*fp_idx+2}: FP{fp_idx+1} Koide",
               abs(K_val - 2/3)/(2/3) < 0.05,
               f"K = {K_val:.6f}")
    else:
        print(f"\n    ✗ τ solver FAILED for FP{fp_idx+1}")
        record(f"T{2*fp_idx+1}: FP{fp_idx+1} r₃₁", False, "τ NaN")
        record(f"T{2*fp_idx+2}: FP{fp_idx+1} Koide", False, "τ NaN")


# ============================================================
# §3. A_tail(g₀) — PEŁNE MAPOWANIE (g-solver)
# ============================================================
print("\n" + "=" * 72)
print("§3. A_tail(g₀) — PEŁNE MAPOWANIE (g-solver, K=g²)")
print("=" * 72)

g0_full_map = np.concatenate([
    np.linspace(0.3, 0.95, 14),
    np.linspace(1.05, 3.0, 20)
])

print(f"\n  {'g₀':>8s}  {'A_tail':>12s}  {'g_min':>8s}  {'g(end)':>8s}  {'R²':>8s}")
print("  " + "-" * 52)

for g0 in g0_full_map:
    r, g, st = solve_K2_g(g0, R_max=120.0, N_pts=4000)
    A, _, R2 = fit_tail(r, g, r_L=20.0, r_R=60.0)
    g_min = np.min(g)
    g_end = g[-1]

    if np.isfinite(A):
        print(f"  {g0:8.3f}  {A:12.6f}  {g_min:8.4f}  {g_end:8.4f}  {R2:8.4f}")
    else:
        print(f"  {g0:8.3f}  {'NaN':>12s}  {g_min:8.4f}  {g_end:8.4f}")


# ============================================================
# §4. DIRECT COMPARISON WITH OLD K_sub=g² (ex106 form)
# ============================================================
print("\n" + "=" * 72)
print("§4. PORÓWNANIE Z STARĄ FORMĄ ODE")
print("=" * 72)

print("""
  Stara ODE (ex106, V = ½(1-g)², K_sub = g²):
    g'' = (1-g)/g² - g'²/g - 2g'/r       [silna siła w g=0!]

  Nowa ODE (V_natural, K = g²):
    g'' = (1-g) - g'²/g - 2g'/r           [słaba siła w g=0]

  Różnica: (1-g)/g² vs (1-g)
""")

def solve_K2_old(g0, R_max=200.0, N_pts=8000):
    """OLD form: g'' = (1-g)/g² - g'²/g - 2g'/r"""
    def rhs(r, y):
        g, gp = y
        g_safe = max(g, 1e-12)
        r_safe = max(r, 1e-10)
        gpp = (1.0 - g_safe) / g_safe**2 - gp**2 / g_safe - 2.0 * gp / r_safe
        return [gp, gpp]

    r0 = 1e-4
    acc0 = (1.0 - g0) / g0**2
    g_init = g0 + acc0 * r0**2 / 6.0
    gp_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.02, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    g_arr = sol.sol(r_arr)[0]
    return r_arr, g_arr, sol.status


def get_A_old(g0, R_max=120.0, r_L=20.0, r_R=60.0):
    try:
        r, g, st = solve_K2_old(g0, R_max=R_max, N_pts=6000)
        if st < 0:
            return np.nan, 0.0
        A, _, R2 = fit_tail(r, g, r_L=r_L, r_R=r_R)
        return A, R2
    except:
        return np.nan, 0.0


# Find φ-FP with old ODE
def ratio_old(g0):
    A1, _ = get_A_old(g0)
    A2, _ = get_A_old(phi * g0)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

print("\n  Skan φ-FP (STARA ODE):")
g0_scan_old = np.linspace(0.7, 1.5, 17)
for g0 in g0_scan_old:
    A1, _ = get_A_old(g0)
    A2, _ = get_A_old(phi * g0)
    if np.isfinite(A1) and A1 > 1e-15 and np.isfinite(A2):
        r = (A2/A1)**4
        if r < 1e6:
            print(f"    g₀ = {g0:.4f}: r₂₁ = {r:.2f}")

# Find old φ-FP
old_fps = []
for gl in np.arange(0.7, 1.5, 0.05):
    try:
        vl = ratio_old(gl)
        vr = ratio_old(gl + 0.05)
        if np.isfinite(vl) and np.isfinite(vr) and vl * vr < 0:
            g0_old = brentq(ratio_old, gl, gl + 0.05, xtol=1e-7)
            old_fps.append(g0_old)
    except:
        continue

print(f"\n  Stara ODE φ-FP(s): {[f'{x:.6f}' for x in old_fps]}")

# Test old FPs
for fp_idx, g0_star_old in enumerate(old_fps):
    g0_tau_old = phi**2 * g0_star_old

    A_e_old, R2e = get_A_old(g0_star_old)
    A_mu_old, R2mu = get_A_old(phi * g0_star_old)

    r21_old = (A_mu_old/A_e_old)**4
    print(f"\n  Old FP{fp_idx+1}: g₀* = {g0_star_old:.8f}")
    print(f"    r₂₁ = {r21_old:.4f}")
    print(f"    g₀_τ = {g0_tau_old:.4f}")

    # τ
    print(f"    τ soliton:")
    for R_max in [100, 150, 200]:
        for r_L, r_R in [(20, 55), (25, 70), (30, 85)]:
            A_tau_old, R2_tau_old = get_A_old(g0_tau_old, R_max=R_max, r_L=r_L, r_R=r_R)
            if np.isfinite(A_tau_old) and A_tau_old > 0:
                r31_old = (A_tau_old/A_e_old)**4
                K_old = koide(1.0, r21_old, r31_old)
                print(f"      R={R_max} [{r_L},{r_R}]: A_τ={A_tau_old:.6f}, r₃₁={r31_old:.1f}, K={K_old:.4f}, R²={R2_tau_old:.4f}")
            else:
                pass  # don't print failures

    # Profile
    r_tau_old, g_tau_old, st_old = solve_K2_old(g0_tau_old, R_max=100.0, N_pts=2000)
    print(f"    τ profile: g ∈ [{np.min(g_tau_old):.4f}, {np.max(g_tau_old):.4f}], g(end)={g_tau_old[-1]:.4f}")


# ============================================================
# SCORECARD
# ============================================================
print("\n" + "=" * 72)
print("SCORECARD")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")
print(f"\n  {n_pass}/{n_total} testów przeszło.")


# ============================================================
# PODSUMOWANIE
# ============================================================
print("\n" + "=" * 72)
print("PODSUMOWANIE ex233")
print("=" * 72)
print("""
  KLUCZOWE PYTANIE: Czy K=g² daje r₃₁ = 3477 (PDG)?

  Dwie φ-FP dla K=g²:
  1. g₀* ≈ 0.87 (nowa, identyczna z g₀ᵉ z α_s)
  2. g₀* ≈ 1.25 (stara, z ex106)

  Stara ODE (V=½(1-g)², K_sub=g²):
    g'' = (1-g)/g² - g'²/g - 2g'/r
    Ma SILNĄ siłę w g=0 → soliton się odbija

  Nowa ODE (V_natural, K=g²):
    g'' = (1-g) - g'²/g - 2g'/r
    Ma SŁABĄ siłę w g=0 → soliton kolapsa

  RÓŻNICA W POTENCJALE, NIE W K!
""")
