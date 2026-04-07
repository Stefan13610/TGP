#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex231_full_tgp_potential_soliton.py
======================================
SOLITON Z PEŁNYM POTENCJAŁEM TGP: V = (γ/7)g⁷ - (γ/8)g⁸

PRZEŁOMOWA OBSERWACJA:
  Dotychczasowe analizy (ex226-230) używały UPROSZCZONEGO potencjału
  V ≈ ½(1-g)² → daje ODE z osobliwością (1-g)/g² dla K=g⁴.

  PEŁNY potencjał TGP: V(g) = (γ/7)g⁷ - (γ/8)g⁸ (β=γ at vacuum)
  V'(g) = γg⁶(1-g)
  V'(g)/g⁴ = γg²(1-g)  ← REGULARNE w g=0!

  ODE w ψ = g³/3 kanoniczny:
    ψ'' + (2/r)ψ' = γ·(3ψ)^{4/3}·(1 - (3ψ)^{1/3})   [FULL TGP]
  vs
    ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}                     [SIMPLIFIED]

  Kluczowa różnica:
  - Przy g=0 (ψ=0): FULL → 0 (brak siły), SIMPLIFIED → 1 (maks siła)
  - Przy g=1 (ψ=1/3): oba → 0 (vacuum)
  - Soliton z pełnym V "miękko" przechodzi przez g<1 zamiast "odbijać się"

  HIPOTEZA: Pełny potencjał daje ln(A) WKLĘSŁE → poprawne r₃₁!

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
# FULL TGP ψ-SOLVER
# ============================================================

def solve_psi_full(g0, gamma=1.0, R_max=200.0, N_pts=8000):
    """
    Solve FULL TGP soliton in canonical ψ = g³/3.

    ODE: ψ'' + (2/r)ψ' = γ·(3ψ)^{4/3}·(1 - (3ψ)^{1/3})

    At vacuum (ψ=1/3, g=1): RHS = 0
    At g=0 (ψ=0): RHS = 0 (REGULAR! No singularity!)
    For g>1: RHS < 0 (restoring)
    For 0<g<1: RHS > 0 (restoring toward g=1)
    """
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        psi_safe = max(psi, 0.0)
        r_safe = max(r, 1e-10)

        if psi_safe < 1e-30:
            # At ψ≈0: RHS ≈ γ·(3ψ)^{4/3} → 0
            force = 0.0
        else:
            g = (3.0 * psi_safe) ** (1.0/3.0)
            force = gamma * g**4 * (1.0 - g)  # γ·(3ψ)^{4/3}·(1-g) = γ·g⁴·(1-g)

        psi_pp = force - 2.0 * psip / r_safe
        return [psip, psi_pp]

    r0 = 1e-4
    # Initial acceleration: ψ''(0) = γ·g₀⁴·(1-g₀) (since ψ'/r → ψ''(0)/3 at r=0)
    acc0 = gamma * g0**4 * (1.0 - g0)

    psi_init = psi0 + acc0 * r0**2 / 6.0
    psip_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    return r_arr, y_arr[0], y_arr[1]


def solve_psi_simplified(g0, R_max=200.0, N_pts=8000):
    """
    Simplified potential: ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}
    (from ex228, for comparison)
    """
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        psi_safe = max(psi, 1e-30)
        r_safe = max(r, 1e-10)
        g = (3.0 * psi_safe) ** (1.0/3.0)
        psi_pp = (1.0 - g) - 2.0 * psip / r_safe
        return [psip, psi_pp]

    r0 = 1e-4
    acc0 = 1.0 - g0
    psi_init = psi0 + acc0 * r0**2 / 6.0
    psip_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    return r_arr, y_arr[0], y_arr[1]


def fit_tail(r_arr, psi_arr, r_L=20.0, r_R=60.0):
    """Fit (ψ - 1/3)·r = B cos(r) + C sin(r)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, 0.0
    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, res, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)
    delta = np.arctan2(B, C)

    y_pred = X @ coefs
    ss_res = np.sum((y_fit - y_pred)**2)
    ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return A, delta, R2


def fit_tail_full(r_arr, psi_arr, gamma=1.0, r_L=20.0, r_R=60.0):
    """
    Tail fit for FULL TGP potential.

    Linearized near g=1: ψ'' + 2ψ'/r + γ·δψ = 0
    Solution: δψ = A·sin(√γ·r + δ)/r

    So fit (ψ-1/3)·r = B·cos(ω·r) + C·sin(ω·r) with ω = √γ.
    """
    omega = np.sqrt(gamma)
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, 0.0
    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit
    X = np.column_stack([np.cos(omega * r_fit), np.sin(omega * r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)

    y_pred = X @ coefs
    ss_res = np.sum((y_fit - y_pred)**2)
    ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return A, np.arctan2(B, C), R2


# ============================================================
# §1. PORÓWNANIE POTENCJAŁÓW
# ============================================================
print("=" * 72)
print("§1. PORÓWNANIE: PEŁNY TGP vs UPROSZCZONY POTENCJAŁ")
print("=" * 72)

print("""
  Pełny TGP:       V(g) = (γ/7)g⁷ - (γ/8)g⁸
                    V'(g) = γg⁶(1-g)
                    V'/g⁴ = γg²(1-g)

  Uproszczony:      V(g) ≈ ½(1-g)²
                    V'(g) = -(1-g)
                    V'/g⁴ = -(1-g)/g⁴ → ∞ w g=0

  W ψ-space (ψ = g³/3):
    Full:       RHS = γg⁴(1-g) → 0 w g=0, 0 w g=1
    Simplified: RHS = (1-g)    → 1 w g=0, 0 w g=1
""")

# Show force profiles
g_range = np.linspace(0.01, 2.5, 200)

print(f"  {'g':>6s}  {'F_full (γ=1)':>14s}  {'F_simplified':>14s}  {'ratio':>10s}")
print("  " + "-" * 50)
for g in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 2.5]:
    F_full = g**4 * (1 - g)  # γ=1
    F_simpl = (1 - g)
    ratio = F_full / F_simpl if abs(F_simpl) > 1e-10 else np.nan
    print(f"  {g:6.2f}  {F_full:14.6f}  {F_simpl:14.6f}  {ratio:10.4f}")


# ============================================================
# §2. PROFIL SOLITONU: FULL vs SIMPLIFIED
# ============================================================
print("\n" + "=" * 72)
print("§2. PROFILE SOLITONOWE: FULL vs SIMPLIFIED")
print("=" * 72)

# Test with g₀ = 2.17 (τ soliton)
g0_test = 2.17

r_full, psi_full, psip_full = solve_psi_full(g0_test, gamma=1.0, R_max=100.0)
r_simpl, psi_simpl, psip_simpl = solve_psi_simplified(g0_test, R_max=100.0)

g_full = (3.0 * np.maximum(psi_full, 0.0)) ** (1.0/3.0)
g_simpl = (3.0 * np.maximum(psi_simpl, 0.0)) ** (1.0/3.0)

print(f"\n  Profil g(r) dla g₀ = {g0_test}:")
print(f"  {'r':>8s}  {'g_full':>10s}  {'g_simpl':>10s}")
print("  " + "-" * 32)
for ri in [0.001, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0]:
    idx_f = np.searchsorted(r_full, ri)
    idx_s = np.searchsorted(r_simpl, ri)
    gf = g_full[min(idx_f, len(g_full)-1)]
    gs = g_simpl[min(idx_s, len(g_simpl)-1)]
    print(f"  {ri:8.3f}  {gf:10.6f}  {gs:10.6f}")

# A_tail comparison
A_full, _, R2_full = fit_tail_full(r_full, psi_full, gamma=1.0, r_L=20.0, r_R=60.0)
A_simpl, _, R2_simpl = fit_tail(r_simpl, psi_simpl, r_L=20.0, r_R=60.0)

print(f"\n  A_tail (g₀ = {g0_test}):")
print(f"    Full TGP:    A = {A_full:.6f}, R² = {R2_full:.4f}")
print(f"    Simplified:  A = {A_simpl:.6f}, R² = {R2_simpl:.4f}")


# ============================================================
# §3. ★ φ-FP Z PEŁNYM POTENCJAŁEM TGP
# ============================================================
print("\n" + "=" * 72)
print("§3. ★ φ-FIXED POINT Z PEŁNYM POTENCJAŁEM TGP")
print("=" * 72)

def get_A_full(g0, gamma=1.0, R_max=100.0):
    """Get A_tail from FULL TGP potential."""
    try:
        r, psi, psip = solve_psi_full(g0, gamma=gamma, R_max=R_max, N_pts=5000)
        A, _, R2 = fit_tail_full(r, psi, gamma=gamma, r_L=20.0, r_R=60.0)
        return A, R2
    except:
        return np.nan, 0.0


def ratio_func_full(g0, gamma=1.0):
    """(A(φg₀)/A(g₀))⁴ - R21_PDG with FULL TGP potential."""
    A1, _ = get_A_full(g0, gamma)
    A2, _ = get_A_full(phi * g0, gamma)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG


# Scan for γ = 1
print("\n  Skan φ-FP (FULL TGP, γ=1):")
g0_scan = np.linspace(0.5, 1.5, 21)
vals = []
for g0 in g0_scan:
    A1, _ = get_A_full(g0)
    A2, _ = get_A_full(phi * g0)
    if np.isfinite(A1) and A1 > 1e-15 and np.isfinite(A2):
        r = (A2/A1)**4
        vals.append((g0, r))
        print(f"    g₀ = {g0:.4f}: r₂₁ = {r:.2f}")
    else:
        print(f"    g₀ = {g0:.4f}: A₁={A1}, A₂={A2}")

# Find zero crossing
g0_star_full = None
for i in range(len(vals)-1):
    g1, r1 = vals[i]
    g2, r2 = vals[i+1]
    if (r1 - R21_PDG) * (r2 - R21_PDG) < 0:
        try:
            g0_star_full = brentq(ratio_func_full, g1, g2, xtol=1e-6)
        except:
            pass
        break

if g0_star_full is not None:
    print(f"\n  ★ g₀*(FULL TGP) = {g0_star_full:.8f}")

    A_e_full, R2e = get_A_full(g0_star_full)
    A_mu_full, R2mu = get_A_full(phi * g0_star_full)
    A_tau_full, R2tau = get_A_full(phi**2 * g0_star_full, R_max=200.0)

    r21_full = (A_mu_full / A_e_full)**4
    print(f"    A_e = {A_e_full:.8f} (R²={R2e:.4f})")
    print(f"    A_μ = {A_mu_full:.8f} (R²={R2mu:.4f})")
    print(f"    A_τ = {A_tau_full:.8f} (R²={R2tau:.4f})")
    print(f"    r₂₁ = {r21_full:.4f}  (PDG: {R21_PDG})")

    if np.isfinite(A_tau_full) and A_tau_full > 0:
        r31_full = (A_tau_full / A_e_full)**4
        r32_full = (A_tau_full / A_mu_full)**4
        K_full = koide(1.0, r21_full, r31_full)

        print(f"    r₃₁ = {r31_full:.2f}  (PDG: {R31_PDG})")
        print(f"    r₃₂ = {r32_full:.2f}  (PDG: {R31_PDG/R21_PDG:.2f})")
        print(f"    Koide = {K_full:.6f}  (2/3 = {2/3:.6f})")
        print(f"    Error r₃₁: {abs(r31_full - R31_PDG)/R31_PDG*100:.2f}%")

        record("T1: r₂₁ FULL TGP",
               abs(r21_full - R21_PDG)/R21_PDG < 0.01,
               f"r₂₁ = {r21_full:.2f}")

        record("T2: r₃₁ FULL TGP",
               abs(r31_full - R31_PDG)/R31_PDG < 0.10,
               f"r₃₁ = {r31_full:.2f} (PDG: {R31_PDG})")

        record("T3: Koide FULL TGP",
               abs(K_full - 2/3)/(2/3) < 0.10,
               f"K = {K_full:.6f}")

        # Shape analysis
        ratio_log = np.log(A_mu_full/A_e_full) / np.log(A_tau_full/A_e_full)
        print(f"\n    ln(A_μ/A_e)/ln(A_τ/A_e) = {ratio_log:.4f}")
        print(f"    Required (PDG): 0.6539")
        print(f"    K=g⁴ simplified: 0.4032")

        record("T4: ln(A) shape",
               abs(ratio_log - 0.654) < 0.10,
               f"ratio = {ratio_log:.4f} (req: 0.654)")
    else:
        print("    τ solver failed!")
        record("T2: r₃₁", False, "τ NaN")
        record("T3: Koide", False, "No τ")
        record("T4: shape", False, "No τ")

    # Compare with g₀ᵉ(α_s)
    g0e = 0.86941
    print(f"\n    g₀*(FULL) = {g0_star_full:.6f}")
    print(f"    g₀ᵉ(α_s) = {g0e}")
    print(f"    diff: {abs(g0_star_full - g0e)/g0e*100:.2f}%")

    record("T5: g₀* = g₀ᵉ identification",
           abs(g0_star_full - g0e)/g0e < 0.05,
           f"g₀* = {g0_star_full:.6f}, g₀ᵉ = {g0e}")
else:
    print("  φ-FP nie znaleziony!")
    record("T1: φ-FP", False, "Not found")


# ============================================================
# §4. A_tail(g₀) MAPPING — FULL vs SIMPLIFIED
# ============================================================
print("\n" + "=" * 72)
print("§4. A_tail(g₀) — FULL TGP vs SIMPLIFIED")
print("=" * 72)

g0_map = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.1, 1.3, 1.5, 1.7, 2.0, 2.17, 2.5]
print(f"\n  {'g₀':>8s}  {'A_full':>12s}  {'A_simpl':>12s}  {'ratio F/S':>10s}  {'R²_full':>8s}")
print("  " + "-" * 56)

for g0 in g0_map:
    A_f, R2_f = get_A_full(g0, R_max=150.0)

    r_s, psi_s, _ = solve_psi_simplified(g0, R_max=150.0, N_pts=5000)
    A_s, _, R2_s = fit_tail(r_s, psi_s, r_L=20.0, r_R=60.0)

    ratio = A_f / A_s if np.isfinite(A_f) and np.isfinite(A_s) and A_s > 0 else np.nan
    print(f"  {g0:8.3f}  {A_f:12.6f}  {A_s:12.6f}  {ratio:10.4f}  {R2_f:8.4f}" if np.isfinite(A_f) else
          f"  {g0:8.3f}  {'NaN':>12s}  {A_s:12.6f}  {'—':>10s}")


# ============================================================
# §5. SKAN γ — CZY γ TUNE'UJE r₃₁?
# ============================================================
print("\n" + "=" * 72)
print("§5. SKAN γ — CZY PARAMETR γ TUNE'UJE r₃₁?")
print("=" * 72)

print("""
  W pełnym ODE: ψ'' + 2ψ'/r = γ·g⁴·(1-g)
  γ kontroluje SIŁĘ potencjału. Zmiana γ zmienia:
  - Częstość ogona: ω = √γ
  - Profil rdzenia: silniejszy potencjał → węższy soliton
  - Stosunek A_tail(g₀_τ)/A_tail(g₀_e)
""")

# For different γ, find φ-FP and compute r₃₁
gammas_to_test = [0.5, 1.0, 2.0, 5.0, 10.0]

print(f"\n  {'γ':>6s}  {'g₀*':>10s}  {'r₂₁':>10s}  {'A_τ/A_e':>10s}  {'r₃₁':>12s}  {'err_r₃₁':>10s}")
print("  " + "-" * 64)

for gam in gammas_to_test:
    def rf_gamma(g0):
        A1, _ = get_A_full(g0, gamma=gam)
        A2, _ = get_A_full(phi * g0, gamma=gam)
        if np.isnan(A1) or A1 < 1e-15:
            return 1e10
        return (A2/A1)**4 - R21_PDG

    # Find φ-FP
    g0_fp = None
    for gl in np.arange(0.3, 2.0, 0.1):
        try:
            vl = rf_gamma(gl)
            vr = rf_gamma(gl + 0.1)
            if np.isfinite(vl) and np.isfinite(vr) and vl * vr < 0:
                g0_fp = brentq(rf_gamma, gl, gl + 0.1, xtol=1e-5)
                break
        except:
            continue

    if g0_fp is not None:
        A_e_g, _ = get_A_full(g0_fp, gamma=gam)
        A_mu_g, _ = get_A_full(phi * g0_fp, gamma=gam)
        A_tau_g, _ = get_A_full(phi**2 * g0_fp, gamma=gam, R_max=200.0)

        r21_g = (A_mu_g/A_e_g)**4

        if np.isfinite(A_tau_g) and A_tau_g > 0:
            ratio_31 = A_tau_g / A_e_g
            r31_g = ratio_31**4
            err = abs(r31_g - R31_PDG) / R31_PDG * 100
            print(f"  {gam:6.1f}  {g0_fp:10.6f}  {r21_g:10.2f}  {ratio_31:10.4f}  {r31_g:12.2f}  {err:9.2f}%")
        else:
            print(f"  {gam:6.1f}  {g0_fp:10.6f}  {r21_g:10.2f}  {'NaN':>10s}")
    else:
        print(f"  {gam:6.1f}  {'NOT FOUND':>10s}")


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
print("PODSUMOWANIE ex231")
print("=" * 72)
print("""
  TESTOWANA HIPOTEZA:
  Pełny potencjał TGP V = (γ/7)g⁷ - (γ/8)g⁸ daje V'/g⁴ = γg²(1-g),
  co jest REGULARNE w g=0 (w przeciwieństwie do uproszczonego (1-g)/g²).
  To może zmienić kształt A(g₀) z wypukłego na wklęsły → poprawne r₃₁.
""")
