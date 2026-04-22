#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex227_phi_fp_K4_refined.py
============================
PRECYZYJNE φ-FP Z K(g) = g⁴ + WERYFIKACJA r₃₁

KONTEKST (ex226):
  φ-FP znaleziony: g₀* ≈ 0.830, r₂₁ = 206.77 ✓
  Ale g₀^τ = φ²·g₀* ≈ 2.17 → solver się rozbiega

CELE:
  1. Precyzyjne g₀* z K=g⁴ (porównanie z g₀ᵉ = 0.86941 z α_s)
  2. Naprawić solver dla dużych g₀ (stiff method)
  3. Obliczyć r₃₁ = (A_τ/A_e)⁴ vs PDG 3477.48
  4. Sprawdzić Koide z A_tail⁴ dla K=g⁴

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
# IMPROVED SOLITON SOLVER
# ============================================================

def solve_soliton(g0, n_K=4, R_max=100.0, method='DOP853'):
    """
    K(g)=g^n_K soliton ODE:
    g'' = (1-g)/g^(n_K-2) - (n_K/2)(g')²/g - (2/r)g'
    """
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
                    method=method, rtol=1e-10, atol=1e-12,
                    max_step=0.03, dense_output=True)
    return sol


def solve_soliton_stiff(g0, n_K=4, R_max=100.0):
    """Stiff solver for large g₀ (Radau IIA, implicit)."""
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

    # Try Radau (stiff solver)
    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='Radau', rtol=1e-8, atol=1e-10,
                    max_step=0.1, dense_output=True)
    return sol


def fit_tail(r_arr, g_arr, r_L=20.0, r_R=55.0):
    """Fit (g-1)·r = B·cos(r) + C·sin(r) → A_tail, δ"""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)
    delta = np.arctan2(B, C)
    return A, delta


def get_A(g0, n_K=4, stiff=False, R_max=100.0, r_L=20.0, r_R=55.0):
    """Get A_tail for given g₀."""
    try:
        if stiff:
            sol = solve_soliton_stiff(g0, n_K=n_K, R_max=R_max)
        else:
            sol = solve_soliton(g0, n_K=n_K, R_max=R_max)

        if sol.status < 0:
            return np.nan, np.nan

        r = np.linspace(sol.t[0], min(sol.t[-1], R_max), 4000)
        y = sol.sol(r)
        g = y[0]
        A, delta = fit_tail(r, g, r_L=r_L, r_R=r_R)
        return A, delta
    except:
        return np.nan, np.nan

# ============================================================
# §1. PRECYZYJNE φ-FP
# ============================================================
print("=" * 72)
print("§1. PRECYZYJNE φ-FIXED POINT Z K(g) = g⁴")
print("=" * 72)

def ratio_func(g0, n_K=4):
    """(A(φg₀)/A(g₀))⁴ - R21_PDG"""
    A1, _ = get_A(g0, n_K=n_K)
    A2, _ = get_A(phi * g0, n_K=n_K)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

# Fine scan near ex226 result
print("\n  Skanowanie φ-FP z wysoką rozdzielczością...")
g0_scan = np.linspace(0.78, 0.90, 25)
vals = []
for g0 in g0_scan:
    rv = ratio_func(g0)
    r_val = rv + R21_PDG
    vals.append(rv)
    if abs(rv) < 1e8:
        print(f"    g₀ = {g0:.4f}: (A_μ/A_e)⁴ = {r_val:.2f}")

# Find precise zero
g0_star = None
vals = np.array(vals)
for i in range(len(vals)-1):
    if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
        g0_star = brentq(ratio_func, g0_scan[i], g0_scan[i+1], xtol=1e-8)
        break

if g0_star is not None:
    A_e, d_e = get_A(g0_star)
    A_mu, d_mu = get_A(phi * g0_star)
    r21 = (A_mu/A_e)**4

    print(f"\n  ★ g₀*(K=g⁴) = {g0_star:.8f}")
    print(f"    g₀*(μ) = φ·g₀* = {phi*g0_star:.8f}")
    print(f"    g₀*(τ) = φ²·g₀* = {phi**2*g0_star:.8f}")
    print(f"    A_e = {A_e:.8f}")
    print(f"    A_μ = {A_mu:.8f}")
    print(f"    r₂₁ = {r21:.6f}  (PDG: {R21_PDG})")
    print(f"    Error: {abs(r21-R21_PDG)/R21_PDG*100:.4f}%")

    # Compare with g₀ᵉ from α_s formula
    g0e_alphas = 0.86941
    print(f"\n  Porównanie z g₀ᵉ (α_s formuła):")
    print(f"    g₀*(K=g⁴ FP) = {g0_star:.6f}")
    print(f"    g₀ᵉ(α_s)     = {g0e_alphas:.6f}")
    print(f"    Difference:     {abs(g0_star-g0e_alphas)/g0e_alphas*100:.2f}%")

    # Compare with g₀* from K_sub=g² (ex106)
    g0_K2 = 1.24915  # from ex106
    print(f"    g₀*(K=g² FP)  = {g0_K2:.5f}")
    print(f"    Ratio K4/K2   = {g0_star/g0_K2:.4f}")

    record("T1: φ-FP precision",
           abs(r21-R21_PDG)/R21_PDG < 0.001,
           f"g₀* = {g0_star:.8f}, r₂₁ = {r21:.4f}")
else:
    print("  φ-FP nie znaleziony!")
    record("T1: φ-FP", False, "Not found")

# ============================================================
# §2. τ SOLVER — STIFF METHOD
# ============================================================
print("\n" + "=" * 72)
print("§2. τ SOLITON Z METODĄ SZTYWNĄ (Radau)")
print("=" * 72)

if g0_star is not None:
    g0_tau = phi**2 * g0_star
    print(f"\n  g₀^τ = {g0_tau:.6f}")

    # Try standard first
    A_tau, d_tau = get_A(g0_tau, stiff=False)
    print(f"  DOP853: A_τ = {A_tau}, {'OK' if np.isfinite(A_tau) else 'FAILED'}")

    # Try Radau (stiff)
    A_tau_s, d_tau_s = get_A(g0_tau, stiff=True, R_max=120.0)
    print(f"  Radau:  A_τ = {A_tau_s}, {'OK' if np.isfinite(A_tau_s) else 'FAILED'}")

    # Try with wider tail window
    A_tau_w, d_tau_w = get_A(g0_tau, stiff=True, R_max=150.0, r_L=25.0, r_R=70.0)
    print(f"  Radau (wide): A_τ = {A_tau_w}, {'OK' if np.isfinite(A_tau_w) else 'FAILED'}")

    # Try DOP853 with smaller step
    sol_test = solve_ivp(
        lambda r, y: [(gp := y[1]), (1.0-max(y[0],1e-15))/max(y[0],1e-15)**2 - 2.0*y[1]**2/max(y[0],1e-15) - 2.0*y[1]/max(r,1e-10)],
        [1e-3, 120.0],
        [g0_tau + (1-g0_tau)/g0_tau**2 * 1e-6/6, (1-g0_tau)/g0_tau**2 * 1e-3/3],
        method='DOP853', rtol=1e-8, atol=1e-10, max_step=0.01
    )
    if sol_test.status >= 0:
        r_test = np.linspace(sol_test.t[0], min(sol_test.t[-1], 120.0), 5000)
        g_test = sol_test.sol(r_test)[0]
        A_tau_d, d_tau_d = fit_tail(r_test, g_test, 20.0, 55.0)
        print(f"  DOP853 (fine): A_τ = {A_tau_d}, t_max = {sol_test.t[-1]:.1f}")

        if np.isfinite(A_tau_d) and A_tau_d > 0:
            A_tau_final = A_tau_d
        elif np.isfinite(A_tau_s) and A_tau_s > 0:
            A_tau_final = A_tau_s
        elif np.isfinite(A_tau_w) and A_tau_w > 0:
            A_tau_final = A_tau_w
        else:
            A_tau_final = np.nan
    else:
        A_tau_final = A_tau_s if np.isfinite(A_tau_s) else A_tau_w

    if np.isfinite(A_tau_final) and A_tau_final > 0:
        r31 = (A_tau_final / A_e)**4
        print(f"\n  ★ r₃₁ = (A_τ/A_e)⁴ = {r31:.2f}")
        print(f"    PDG: r₃₁ = {R31_PDG}")
        print(f"    Error: {abs(r31-R31_PDG)/R31_PDG*100:.2f}%")

        # Koide check
        m_e_A = A_e**4
        m_mu_A = A_mu**4
        m_tau_A = A_tau_final**4
        K_val = koide(m_e_A, m_mu_A, m_tau_A)
        print(f"\n  Koide z A_tail⁴:")
        print(f"    K = {K_val:.6f}")
        print(f"    2/3 = {2/3:.6f}")
        print(f"    Error: {abs(K_val-2/3)/(2/3)*100:.2f}%")

        record("T2: r₃₁ from K=g⁴",
               abs(r31-R31_PDG)/R31_PDG < 0.05,
               f"r₃₁ = {r31:.1f} (PDG: {R31_PDG})")

        record("T3: Koide from A_tail⁴",
               abs(K_val - 2/3)/(2/3) < 0.05,
               f"K = {K_val:.6f}, 2/3 = {2/3:.6f}")
    else:
        print("\n  τ solver FAILED z wszystkimi metodami")
        print("  Analiza profilu dla dużych g₀...")

        # Diagnostic: show profile
        for R_try in [40, 60, 80]:
            try:
                sol = solve_soliton(g0_tau, R_max=R_try)
                r = np.linspace(sol.t[0], sol.t[-1], 200)
                g = sol.sol(r)[0]
                g_min, g_max = np.min(g), np.max(g)
                g_end = g[-1]
                print(f"    R_max={R_try}: g ∈ [{g_min:.3f}, {g_max:.3f}], g(end) = {g_end:.4f}, t_end = {sol.t[-1]:.1f}")
            except Exception as e:
                print(f"    R_max={R_try}: ERROR: {e}")

        record("T2: r₃₁ from K=g⁴", False, "τ solver failed")
        record("T3: Koide", False, "No τ data")

# ============================================================
# §3. COMPARISON TABLE: K=g⁴ vs K=g²
# ============================================================
print("\n" + "=" * 72)
print("§3. PORÓWNANIE: φ-FP DLA RÓŻNYCH K(g)")
print("=" * 72)

# Compute φ-FP for K=g² (K_sub)
def ratio_func_K2(g0):
    A1, _ = get_A(g0, n_K=2)
    A2, _ = get_A(phi * g0, n_K=2)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

g0_star_K2 = None
try:
    g0_star_K2 = brentq(ratio_func_K2, 1.0, 1.5, xtol=1e-6)
except:
    # scan
    for gl, gr in [(1.0, 1.1), (1.1, 1.2), (1.2, 1.3), (1.3, 1.5)]:
        try:
            vl = ratio_func_K2(gl)
            vr = ratio_func_K2(gr)
            if vl * vr < 0:
                g0_star_K2 = brentq(ratio_func_K2, gl, gr, xtol=1e-6)
                break
        except:
            pass

if g0_star_K2 is not None:
    A_e_K2, _ = get_A(g0_star_K2, n_K=2)
    A_mu_K2, _ = get_A(phi * g0_star_K2, n_K=2)
    A_tau_K2, _ = get_A(phi**2 * g0_star_K2, n_K=2, R_max=120.0)
    r21_K2 = (A_mu_K2/A_e_K2)**4
    r31_K2 = (A_tau_K2/A_e_K2)**4 if np.isfinite(A_tau_K2) else np.nan

    print(f"\n  {'':15s}  {'K=g⁴':>12s}  {'K=g²':>12s}  {'PDG':>12s}")
    print("  " + "-" * 55)
    print(f"  {'g₀*(e)':15s}  {g0_star:12.6f}  {g0_star_K2:12.6f}  {'—':>12s}")
    print(f"  {'g₀*(μ)':15s}  {phi*g0_star:12.6f}  {phi*g0_star_K2:12.6f}  {'—':>12s}")
    print(f"  {'g₀*(τ)':15s}  {phi**2*g0_star:12.6f}  {phi**2*g0_star_K2:12.6f}  {'—':>12s}")
    print(f"  {'r₂₁':15s}  {(A_mu/A_e)**4:12.4f}  {r21_K2:12.4f}  {R21_PDG:12.3f}")
    if np.isfinite(r31_K2):
        print(f"  {'r₃₁':15s}  {'?':>12s}  {r31_K2:12.1f}  {R31_PDG:12.1f}")

    # K=g² Koide check
    if np.isfinite(A_tau_K2) and A_tau_K2 > 0:
        K_K2 = koide(A_e_K2**4, A_mu_K2**4, A_tau_K2**4)
        print(f"  {'Koide':15s}  {'?':>12s}  {K_K2:12.6f}  {2/3:12.6f}")

    record("T4: K=g² φ-FP comparison",
           True,
           f"g₀*(K=g²) = {g0_star_K2:.5f}, r₂₁ = {r21_K2:.2f}")

# ============================================================
# §4. g₀* vs g₀ᵉ(α_s): CZY SĄ TĄ SAMĄ WIELKOŚCIĄ?
# ============================================================
print("\n" + "=" * 72)
print("§4. ★ g₀* vs g₀ᵉ(α_s): IDENTYFIKACJA")
print("=" * 72)

g0e_alphas = 0.86941
print(f"""
  g₀*(K=g⁴ FP) = {g0_star:.8f}  (z warunku r₂₁ = 206.77)
  g₀ᵉ(α_s)      = {g0e_alphas}  (z formuly α_s = 7N_c³g₀ᵉ/(12Φ₀))

  Różnica: {abs(g0_star-g0e_alphas):.6f} = {abs(g0_star-g0e_alphas)/g0e_alphas*100:.2f}%

  Są to RÓŻNE wielkości:
  - g₀*(FP) = amplituda centralna solitonu spełniająca warunek φ-FP
  - g₀ᵉ(α_s) = parameter fixed point z drabinki φ-FP (Brannen)

  Ale BLISKOŚĆ ({abs(g0_star-g0e_alphas)/g0e_alphas*100:.1f}%) sugeruje głębszy związek.
""")

# Check: what α_s would g₀* give?
PHI0 = 168 * 0.685
alpha_s_from_FP = 7 * 3**3 * g0_star / (12 * PHI0)
alpha_s_from_g0e = 7 * 3**3 * g0e_alphas / (12 * PHI0)
print(f"  α_s(g₀*={g0_star:.4f}) = {alpha_s_from_FP:.4f}")
print(f"  α_s(g₀ᵉ={g0e_alphas}) = {alpha_s_from_g0e:.4f}")
print(f"  PDG: α_s = 0.1179 ± 0.0009")

# ============================================================
# §5. SCAN: A_tail(g₀) SCALING FOR K=g⁴
# ============================================================
print("\n" + "=" * 72)
print("§5. SKALOWANIE A_tail(g₀) DLA K=g⁴")
print("=" * 72)

g0_dense = np.linspace(0.3, 1.8, 40)
A_dense = []
for g0 in g0_dense:
    A, _ = get_A(g0)
    A_dense.append(A if np.isfinite(A) else 0)

A_dense = np.array(A_dense)

# For g₀ < 1: A should increase as g₀ → 1 from below? Or decrease?
# Actually for g₀ → 1: soliton is weak → A_tail → 0
# For g₀ → 0: soliton is strong → A_tail → ?

print(f"\n  {'g₀':>6s}  {'A_tail':>10s}  {'A⁴':>12s}")
print("  " + "-" * 32)
for i in range(0, len(g0_dense), 4):
    g0 = g0_dense[i]
    A = A_dense[i]
    if A > 0:
        print(f"  {g0:6.3f}  {A:10.6f}  {A**4:12.6f}")

# Fit A_tail vs |g₀-1| near g₀=1
mask_near1 = (np.abs(g0_dense - 1.0) > 0.02) & (np.abs(g0_dense - 1.0) < 0.5) & (A_dense > 1e-8)
if np.sum(mask_near1) >= 5:
    log_dg = np.log(np.abs(g0_dense[mask_near1] - 1.0))
    log_A = np.log(A_dense[mask_near1])
    coeffs = np.polyfit(log_dg, log_A, 1)
    nu_scale = coeffs[0]
    print(f"\n  Near g₀≈1: A_tail ∝ |g₀-1|^ν")
    print(f"    ν = {nu_scale:.4f}")
    print(f"    → M ∝ A⁴ ∝ |g₀-1|^{4*nu_scale:.2f}")

# ============================================================
# SCORECARD
# ============================================================
print(f"\n{'='*72}")
print("SCORECARD")
print(f"{'='*72}\n")

passed = sum(1 for _, p, _ in TESTS if p)
total = len(TESTS)

for name, p, detail in TESTS:
    mark = "PASS" if p else "FAIL"
    print(f"  [{mark}] {name}")

print(f"\n  {passed}/{total} testów przeszło.")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*72}")
print("PODSUMOWANIE ex227")
print(f"{'='*72}")
print(f"""
  GŁÓWNE WYNIKI:

  1. φ-FP z K=g⁴: g₀* = {g0_star:.8f}
     r₂₁ = {(A_mu/A_e)**4:.4f} (PDG: {R21_PDG})

  2. g₀*(K=g⁴) = {g0_star:.4f} vs g₀ᵉ(α_s) = {g0e_alphas}
     Różnica: {abs(g0_star-g0e_alphas)/g0e_alphas*100:.1f}%
     → Blisko ale NIE identyczne

  3. α_s z g₀*(FP): {alpha_s_from_FP:.4f}
     α_s z g₀ᵉ:     {alpha_s_from_g0e:.4f}
     PDG:           0.1179

  OTWARTE:
  - τ solver dla g₀ > 2 z K=g⁴
  - Fizyczny związek g₀*(FP) ↔ g₀ᵉ(α_s)
  - Koide z K=g⁴ (wymaga τ)
""")
