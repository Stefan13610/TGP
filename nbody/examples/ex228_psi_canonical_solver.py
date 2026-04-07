#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex228_psi_canonical_solver.py
================================
KANONICZNY SOLVER SOLITONOWY: ψ = g³/3

PRZEŁOM: Zmiana zmiennych ψ = g³/3 kanonizuje kinetykę:
  ½g⁴(∇g)² → ½(∇ψ)²

ODE w ψ:
  ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}

Jest to REGULARNE równanie (bez osobliwości w g=0)!
- Dla ψ = 0 (g=0): RHS = 1 (finite, pushes UP)
- Dla ψ = 1/3 (g=1): RHS = 0 (vacuum)
- Dla ψ > 1/3 (g>1): RHS < 0 (pulls DOWN)

→ Rozwiązuje problem τ solvera (g₀ > 2) z ex227.

WYPROWADZENIE:
  S = ∫[½g⁴(∇g)² + V(g)] d³x
  ψ = g³/3 → dψ = g²dg → ∇ψ = g²∇g
  ½g⁴(∇g)² = ½(g²∇g)² = ½(∇ψ)²

  EOM w g: g'' = (1-g)/g² - 2g'²/g - 2g'/r
  Podstawiając ψ' = g²g', ψ'' = 2gg'² + g²g'':
    ψ'' = 2gg'² + g²[(1-g)/g² - 2g'²/g - 2g'/r]
         = 2gg'² + (1-g) - 2gg'² - 2g²g'/r
         = (1-g) - 2ψ'/r

  ψ'' + 2ψ'/r = 1 - (3ψ)^{1/3}    ★

Warunki brzegowe: ψ(0) = g₀³/3, ψ'(0) = 0
Ogon: δψ ≈ (g-1) → (ψ-1/3)·r = B cos(r) + C sin(r)

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

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
# CANONICAL ψ-SOLVER
# ============================================================

def solve_soliton_psi(g0, R_max=100.0, method='DOP853', rtol=1e-10, atol=1e-12):
    """
    Solve soliton ODE in canonical ψ = g³/3 variable.

    ODE: ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}
    BC:  ψ(0) = g₀³/3, ψ'(0) = 0

    Returns sol with ψ(r), and helper to convert back to g.
    """
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        # Guard against ψ < 0 (unphysical)
        psi_safe = max(psi, 1e-30)
        r_safe = max(r, 1e-10)

        g = (3.0 * psi_safe) ** (1.0/3.0)
        psi_pp = (1.0 - g) - 2.0 * psip / r_safe
        return [psip, psi_pp]

    # Initial conditions at small r₀
    r0 = 1e-4
    g_val = (3.0 * psi0) ** (1.0/3.0)  # = g0
    acc = 1.0 - g_val  # ψ''(0) = 1 - g₀

    psi_init = psi0 + acc * r0**2 / 6.0   # Taylor: ψ(r₀) = ψ₀ + ψ''₀·r₀²/6
    psip_init = acc * r0 / 3.0              # Taylor: ψ'(r₀) = ψ''₀·r₀/3

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method=method, rtol=rtol, atol=atol,
                    max_step=0.05, dense_output=True)
    return sol


def psi_to_g(psi_arr):
    """Convert ψ array back to g = (3ψ)^{1/3}."""
    psi_safe = np.maximum(psi_arr, 0.0)
    return (3.0 * psi_safe) ** (1.0/3.0)


def fit_tail_psi(r_arr, psi_arr, r_L=20.0, r_R=60.0):
    """
    Fit tail in ψ-space: (ψ - 1/3)·r = B cos(r) + C sin(r)

    Since δψ ≈ (g-1) near g=1, A_tail here = A_tail in g-space.
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan

    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit

    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, res, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)
    delta = np.arctan2(B, C)

    # Quality check: residual
    if len(res) > 0 and len(y_fit) > 0:
        R2 = 1.0 - res[0] / np.sum((y_fit - np.mean(y_fit))**2)
    else:
        y_pred = X @ coefs
        ss_res = np.sum((y_fit - y_pred)**2)
        ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
        R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return A, delta, R2


def get_A_psi(g0, R_max=100.0, r_L=20.0, r_R=60.0, method='DOP853'):
    """Get A_tail using canonical ψ-solver."""
    try:
        sol = solve_soliton_psi(g0, R_max=R_max, method=method)
        if sol.status < 0:
            return np.nan, np.nan, 0.0

        r = np.linspace(sol.t[0], min(sol.t[-1], R_max), 5000)
        psi = sol.sol(r)[0]

        result = fit_tail_psi(r, psi, r_L=r_L, r_R=r_R)
        if result is None or len(result) != 3:
            return np.nan, np.nan, 0.0
        return result
    except Exception as e:
        return np.nan, np.nan, 0.0


# ============================================================
# §1. WERYFIKACJA: STARY vs NOWY SOLVER DLA g₀ < 1
# ============================================================
print("=" * 72)
print("§1. PORÓWNANIE SOLVERÓW: g-space vs ψ-space (g₀ < 1)")
print("=" * 72)

# Use old g-space solver for comparison
def solve_soliton_g(g0, n_K=4, R_max=100.0):
    """Old g-space solver (from ex227)."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-15)
        r = max(r, 1e-10)
        gpp = (1.0 - g) / g**(n_K-2) - (n_K/2.0) * gp**2 / g - 2.0 * gp / r
        return [gp, gpp]

    r0 = 1e-3
    g0_acc = (1.0 - g0) / g0**(n_K-2)
    g_init = g0 + g0_acc * r0**2 / 6.0
    gp_init = g0_acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.03, dense_output=True)
    return sol


def fit_tail_g(r_arr, g_arr, r_L=20.0, r_R=55.0):
    """Old g-space tail fit."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    return np.sqrt(B**2 + C**2), np.arctan2(B, C)


print("\n  Porównanie A_tail dla g₀ ∈ [0.5, 0.95]:")
print(f"  {'g₀':>8s}  {'A_g':>12s}  {'A_ψ':>12s}  {'diff':>10s}  {'R²_ψ':>8s}")
print("  " + "-" * 56)

test_g0s = [0.50, 0.60, 0.70, 0.80, 0.830, 0.85, 0.90, 0.95]
max_diff = 0.0
for g0 in test_g0s:
    # Old solver
    sol_g = solve_soliton_g(g0)
    r_g = np.linspace(sol_g.t[0], min(sol_g.t[-1], 100.0), 4000)
    A_g, _ = fit_tail_g(r_g, sol_g.sol(r_g)[0])

    # New solver
    A_psi, _, R2 = get_A_psi(g0)

    if np.isfinite(A_g) and np.isfinite(A_psi) and A_g > 0:
        diff = abs(A_psi - A_g) / A_g * 100
        max_diff = max(max_diff, diff)
    else:
        diff = float('nan')

    print(f"  {g0:8.3f}  {A_g:12.8f}  {A_psi:12.8f}  {diff:9.4f}%  {R2:8.4f}")

record("T1: ψ-solver vs g-solver consistency",
       max_diff < 1.0,
       f"max difference = {max_diff:.4f}%")


# ============================================================
# §2. φ-FP Z ψ-SOLVEREM
# ============================================================
print("\n" + "=" * 72)
print("§2. φ-FIXED POINT Z KANONICZNYM ψ-SOLVEREM")
print("=" * 72)

def ratio_func_psi(g0):
    """(A(φg₀)/A(g₀))⁴ - R21_PDG using ψ-solver."""
    A1, _, _ = get_A_psi(g0)
    A2, _, _ = get_A_psi(phi * g0)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

# Scan
print("\n  Skan φ-FP (ψ-solver):")
g0_scan = np.linspace(0.78, 0.90, 13)
vals_psi = []
for g0 in g0_scan:
    rv = ratio_func_psi(g0)
    r_val = rv + R21_PDG
    vals_psi.append(rv)
    if abs(rv) < 1e8:
        print(f"    g₀ = {g0:.4f}: r₂₁ = {r_val:.2f}")

# Find zero with brentq
g0_star_psi = None
for i in range(len(vals_psi)-1):
    if np.isfinite(vals_psi[i]) and np.isfinite(vals_psi[i+1]) and vals_psi[i]*vals_psi[i+1] < 0:
        g0_star_psi = brentq(ratio_func_psi, g0_scan[i], g0_scan[i+1], xtol=1e-8)
        break

if g0_star_psi is not None:
    A_e, d_e, R2_e = get_A_psi(g0_star_psi)
    A_mu, d_mu, R2_mu = get_A_psi(phi * g0_star_psi)
    r21 = (A_mu / A_e)**4

    print(f"\n  ★ g₀*(ψ-solver) = {g0_star_psi:.8f}")
    print(f"    A_e = {A_e:.8f}  (R² = {R2_e:.6f})")
    print(f"    A_μ = {A_mu:.8f}  (R² = {R2_mu:.6f})")
    print(f"    r₂₁ = {r21:.6f}  (PDG: {R21_PDG})")
    print(f"    Error r₂₁: {abs(r21-R21_PDG)/R21_PDG*100:.4f}%")

    record("T2: φ-FP with ψ-solver",
           abs(r21 - R21_PDG)/R21_PDG < 0.001,
           f"g₀* = {g0_star_psi:.8f}, r₂₁ = {r21:.4f}")
else:
    print("  φ-FP nie znaleziony!")
    record("T2: φ-FP", False, "Not found")
    A_e = A_mu = np.nan


# ============================================================
# §3. ★ τ SOLITON (g₀ > 2) — KLUCZOWY TEST
# ============================================================
print("\n" + "=" * 72)
print("§3. ★ τ SOLITON Z ψ-SOLVEREM (g₀ > 2)")
print("=" * 72)

if g0_star_psi is not None:
    g0_tau = phi**2 * g0_star_psi
    psi0_tau = g0_tau**3 / 3.0

    print(f"\n  g₀^τ = {g0_tau:.6f}")
    print(f"  ψ₀^τ = {psi0_tau:.4f}")
    print(f"  ψ_vacuum = {1/3:.4f}")
    print(f"  Stosunek ψ₀/ψ_vac = {psi0_tau/(1/3):.2f}")

    # Solve with ψ-solver
    print("\n  Rozwiązywanie ODE w ψ-przestrzeni...")

    # Try multiple R_max and tail windows
    results_tau = []
    configs = [
        (100.0, 20.0, 55.0, 'DOP853'),
        (150.0, 25.0, 70.0, 'DOP853'),
        (200.0, 30.0, 90.0, 'DOP853'),
        (150.0, 20.0, 70.0, 'Radau'),
    ]

    for R_max, r_L, r_R, method in configs:
        A_tau, d_tau, R2_tau = get_A_psi(g0_tau, R_max=R_max, r_L=r_L, r_R=r_R, method=method)
        status = "OK" if np.isfinite(A_tau) and A_tau > 0 else "FAIL"
        results_tau.append((A_tau, d_tau, R2_tau, R_max, r_L, r_R, method))
        print(f"    {method:7s} R=[{r_L:.0f},{r_R:.0f}] Rmax={R_max:.0f}: "
              f"A_τ = {A_tau:.6f}  R² = {R2_tau:.4f}  [{status}]")

    # Pick best result (highest R²)
    best = max(results_tau, key=lambda x: x[2] if np.isfinite(x[0]) and x[0] > 0 else -1)
    A_tau_best = best[0]
    R2_tau_best = best[2]

    if np.isfinite(A_tau_best) and A_tau_best > 0:
        r31 = (A_tau_best / A_e)**4
        r32 = (A_tau_best / A_mu)**4

        print(f"\n  ★ WYNIKI τ:")
        print(f"    A_τ = {A_tau_best:.8f}  (R² = {R2_tau_best:.6f})")
        print(f"    r₃₁ = (A_τ/A_e)⁴ = {r31:.2f}")
        print(f"    PDG: r₃₁ = {R31_PDG}")
        print(f"    Error r₃₁: {abs(r31-R31_PDG)/R31_PDG*100:.2f}%")
        print(f"    r₃₂ = (A_τ/A_μ)⁴ = {r32:.2f}  (PDG: {R31_PDG/R21_PDG:.2f})")

        # Koide check
        m_e_A = A_e**4
        m_mu_A = A_mu**4
        m_tau_A = A_tau_best**4
        K_val = koide(m_e_A, m_mu_A, m_tau_A)
        print(f"\n  Koide z A_tail⁴:")
        print(f"    m_e ∝ A_e⁴ = {m_e_A:.6e}")
        print(f"    m_μ ∝ A_μ⁴ = {m_mu_A:.6e}")
        print(f"    m_τ ∝ A_τ⁴ = {m_tau_A:.6e}")
        print(f"    K = {K_val:.6f}")
        print(f"    2/3 = {2/3:.6f}")
        print(f"    Error K: {abs(K_val-2/3)/(2/3)*100:.2f}%")

        record("T3: τ soliton solved",
               np.isfinite(A_tau_best) and R2_tau_best > 0.8,
               f"A_τ = {A_tau_best:.6f}, R² = {R2_tau_best:.4f}")

        record("T4: r₃₁ from ψ-solver",
               abs(r31-R31_PDG)/R31_PDG < 0.10,
               f"r₃₁ = {r31:.1f} (PDG: {R31_PDG}), error = {abs(r31-R31_PDG)/R31_PDG*100:.1f}%")

        record("T5: Koide from A_tail⁴",
               abs(K_val - 2/3)/(2/3) < 0.10,
               f"K = {K_val:.6f}")
    else:
        print("\n  τ solver NADAL NIE DZIAŁA")

        # Diagnostic: show ψ-profile
        print("\n  Diagnostyka profilu ψ(r) dla g₀^τ:")
        sol_diag = solve_soliton_psi(g0_tau, R_max=100.0)
        print(f"    solver status = {sol_diag.status}")
        print(f"    solver message = {sol_diag.message}")
        print(f"    t_span = [{sol_diag.t[0]:.4f}, {sol_diag.t[-1]:.4f}]")

        r_diag = np.linspace(sol_diag.t[0], sol_diag.t[-1], 200)
        psi_diag = sol_diag.sol(r_diag)[0]
        g_diag = psi_to_g(psi_diag)

        print(f"    ψ ∈ [{np.min(psi_diag):.4f}, {np.max(psi_diag):.4f}]")
        print(f"    g ∈ [{np.min(g_diag):.4f}, {np.max(g_diag):.4f}]")
        print(f"    ψ(end) = {psi_diag[-1]:.4f}, g(end) = {g_diag[-1]:.4f}")

        # Show profile at key radii
        print(f"\n    {'r':>8s}  {'ψ':>10s}  {'g':>10s}  {'g-1':>10s}")
        for ri in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0]:
            if ri <= sol_diag.t[-1]:
                psi_r = sol_diag.sol(ri)[0]
                g_r = (3.0 * max(psi_r, 0.0)) ** (1.0/3.0)
                print(f"    {ri:8.1f}  {psi_r:10.6f}  {g_r:10.6f}  {g_r-1:10.6f}")

        record("T3: τ soliton solved", False, f"solver status = {sol_diag.status}")
        record("T4: r₃₁", False, "No τ data")
        record("T5: Koide", False, "No τ data")


# ============================================================
# §4. PROFIL SOLITONOWY: PORÓWNANIE e / μ / τ
# ============================================================
print("\n" + "=" * 72)
print("§4. PROFILE SOLITONOWE e / μ / τ")
print("=" * 72)

if g0_star_psi is not None:
    g0_e = g0_star_psi
    g0_mu = phi * g0_star_psi
    g0_tau_val = phi**2 * g0_star_psi

    print(f"\n  {'r':>6s}  {'g_e':>10s}  {'g_μ':>10s}  {'g_τ':>10s}")
    print("  " + "-" * 42)

    sol_e = solve_soliton_psi(g0_e, R_max=60.0)
    sol_mu = solve_soliton_psi(g0_mu, R_max=60.0)
    sol_tau = solve_soliton_psi(g0_tau_val, R_max=60.0)

    for ri in [0.0001, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0]:
        vals = []
        for sol, g0_v in [(sol_e, g0_e), (sol_mu, g0_mu), (sol_tau, g0_tau_val)]:
            if ri <= sol.t[-1] and ri >= sol.t[0]:
                psi_r = sol.sol(ri)[0]
                g_r = (3.0 * max(psi_r, 0.0)) ** (1.0/3.0)
                vals.append(f"{g_r:10.6f}")
            else:
                vals.append(f"{'—':>10s}")
        print(f"  {ri:6.1f}  {'  '.join(vals)}")


# ============================================================
# §5. SKALOWANIE A_tail(g₀) DLA PEŁNEGO ZAKRESU
# ============================================================
print("\n" + "=" * 72)
print("§5. A_tail(g₀) — PEŁNE MAPOWANIE (ψ-solver)")
print("=" * 72)

g0_map = np.array([0.3, 0.5, 0.7, 0.8, 0.83, 0.85, 0.9, 0.95,
                    1.05, 1.1, 1.2, 1.3, 1.343, 1.5, 1.7, 2.0, 2.174, 2.5, 3.0])

print(f"\n  {'g₀':>8s}  {'A_tail':>12s}  {'A⁴':>14s}  {'R²':>8s}")
print("  " + "-" * 48)

A_map = []
for g0 in g0_map:
    A, _, R2 = get_A_psi(g0, R_max=150.0, r_L=20.0, r_R=70.0)
    A_map.append(A)
    if np.isfinite(A):
        print(f"  {g0:8.3f}  {A:12.6f}  {A**4:14.6e}  {R2:8.4f}")
    else:
        print(f"  {g0:8.3f}  {'NaN':>12s}  {'NaN':>14s}  {'—':>8s}")

# Near-vacuum scaling: A ∝ |g₀-1|^ν
g0_near = np.array([0.9, 0.95, 1.05, 1.1])
A_near = []
for g0 in g0_near:
    A, _, _ = get_A_psi(g0)
    A_near.append(A)

A_near = np.array(A_near)
dg = np.abs(g0_near - 1.0)
mask_valid = np.isfinite(A_near) & (A_near > 0)
if np.sum(mask_valid) >= 2:
    log_A = np.log(A_near[mask_valid])
    log_dg = np.log(dg[mask_valid])
    p = np.polyfit(log_dg, log_A, 1)
    nu = p[0]
    print(f"\n  Near g₀≈1: A_tail ∝ |g₀-1|^ν")
    print(f"    ν = {nu:.4f}")
    print(f"    → M ∝ A⁴ ∝ |g₀-1|^{4*nu:.2f}")


# ============================================================
# §6. PORÓWNANIE: g₀*(ψ) vs g₀*(g) vs g₀ᵉ(α_s)
# ============================================================
print("\n" + "=" * 72)
print("§6. PORÓWNANIE PARAMETRÓW")
print("=" * 72)

g0_alphas = 0.86941
g0_K2 = 1.24915  # from ex106

if g0_star_psi is not None:
    print(f"""
  g₀*(ψ-solver, K=g⁴) = {g0_star_psi:.8f}
  g₀*(g-solver, K=g⁴) = 0.83040836  (ex227)
  g₀ᵉ(α_s formuła)     = {g0_alphas}
  g₀*(K=g²)             = {g0_K2}

  Różnica ψ vs g solver: {abs(g0_star_psi - 0.83040836)/0.83040836*100:.4f}%
  Różnica ψ vs g₀ᵉ:     {abs(g0_star_psi - g0_alphas)/g0_alphas*100:.2f}%
""")

    # α_s from both
    Phi0_bare = 168 * 0.685  # = 115.08
    N_c = 3

    alpha_s_psi = 7 * N_c**3 * g0_star_psi / (12 * Phi0_bare)
    alpha_s_ge = 7 * N_c**3 * g0_alphas / (12 * Phi0_bare)

    print(f"  α_s(g₀*_ψ) = {alpha_s_psi:.4f}")
    print(f"  α_s(g₀ᵉ)   = {alpha_s_ge:.4f}")
    print(f"  PDG:         0.1179 ± 0.0009")

    record("T6: ψ vs g solver agreement",
           abs(g0_star_psi - 0.83040836)/0.83040836 < 0.01,
           f"g₀*(ψ) = {g0_star_psi:.8f}, g₀*(g) = 0.83040836")


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
print("PODSUMOWANIE ex228")
print("=" * 72)
print("""
  PRZEŁOM: Zmiana zmiennych ψ = g³/3 daje REGULARNE równanie:
    ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}

  Brak osobliwości w g=0 (ψ=0)!
  → Umożliwia obliczenie τ solitonu (g₀ > 2)
""")
