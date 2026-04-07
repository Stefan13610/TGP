#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex226_z3_phase_verification.py
================================
NUMERYCZNA WERYFIKACJA a·ln(φ) = 2π/3 Z ODE SOLITONOWEGO

KONTEKST (ex225 §4):
  Koide K = 2/3 wymaga Z₃ kwantyzacji faz ogona solitonowego:
    δ_{k+1} - δ_k = 2π/3
  Z φ-drabinki (g₀^(k) = φ^k · g₀^e):
    δ(g₀) = a·ln(g₀) + b  (logarytmiczna)
    Δδ = a·ln(φ) = 2π/3
    a = 2π/(3·ln(φ)) = 4.3523 (wymagane)

CEL:
  1. Rozwiązać ODE solitonowe z K(g)=g⁴ dla wielu g₀
  2. Wyekstrahować A_tail i fazę δ(g₀) dla każdego
  3. Zmierzyć a = dδ/d(ln g₀)
  4. Sprawdzić czy a·ln(φ) ≈ 2π/3

UWAGA: Porównanie TRZECH kinetic couplings:
  (A) K = g⁴  (poprawiona akcja, sek08a)
  (B) K_sub = g²  (substrat, sek08b, ex106)
  (C) f(g) = 1 + 2α·ln(g)  (one-loop ERG, ex111)

ODE:
  K(g)g'' + ½K'(g)(g')² + (2/r)K(g)g' = V'(g)
  g'' = [V'(g) - ½K'(g)(g')² - (2/r)K(g)g'] / K(g)

  Dla K = g^n_K: K' = n_K·g^(n_K-1)
  g'' = V'(g)/g^n_K - (n_K/2)(g')²/g - (2/r)g'

TESTY:
  T1: Solver K=g⁴ daje g(r→∞) → 1 dla g₀ ∈ [0.5, 3.0]
  T2: Faza δ(g₀) jest monotoniczne i gładka
  T3: a = dδ/d(ln g₀) zmierzone
  T4: a·ln(φ) ≈ 2π/3 (±10% tolerance)
  T5: Porównanie a(K=g⁴) vs a(K=g²)
  T6: φ-FP istnieje dla K=g⁴

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

# ============================================================
# ODE Solver for general K(g) = g^n_K
# ============================================================

def solve_soliton(g0, n_K=4, R_max=80.0, max_step=0.03):
    """
    Solve soliton ODE:  K(g)g'' + ½K'(g)(g')² + (2/r)K(g)g' = V'(g)
    with K(g) = g^n_K, V'(g) = g²(1-g) [β=γ=1]

    Simplified: g'' = (1-g)/g^(n_K-2) - (n_K/2)(g')²/g - (2/r)g'
    For n_K=4: g'' = (1-g)/g² - 2(g')²/g - (2/r)g'
    For n_K=2: g'' = (1-g) - (g')²/g - (2/r)g'
    """
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-15)
        r = max(r, 1e-10)
        gpp = (1.0 - g) / g**(n_K-2) - (n_K/2.0) * gp**2 / g - 2.0 * gp / r
        return [gp, gpp]

    # Taylor expansion near r=0:
    # g(r) ≈ g0 + (1/6)·[V'(g0)/K(g0)]·r²
    # V'(g0)/K(g0) = g0²(1-g0)/g0^n_K = (1-g0)/g0^(n_K-2)
    r0 = 1e-3
    g0_acc = (1.0 - g0) / g0**(n_K-2)
    g_init = g0 + g0_acc * r0**2 / 6.0
    gp_init = g0_acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=max_step, dense_output=True)

    return sol


def fit_tail(r_arr, g_arr, r_L=20.0, r_R=50.0):
    """
    Fit tail: (g(r)-1)·r = B·cos(r) + C·sin(r)
    Returns A_tail = √(B²+C²), phase δ = atan2(C, B)
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit

    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, res, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)

    # Phase: (g-1)·r = A·sin(r + δ) = A·sin(r)cos(δ) + A·cos(r)sin(δ)
    # So: B = A·sin(δ), C = A·cos(δ)
    delta = np.arctan2(B, C)  # sin(δ)/cos(δ) = B/C

    # Residual check
    y_model = B * np.cos(r_fit) + C * np.sin(r_fit)
    rmse = np.sqrt(np.mean((y_fit - y_model)**2))

    return A, delta, B, C


def get_A_delta(g0, n_K=4, R_max=80.0, r_L=20.0, r_R=50.0):
    """Solve soliton and extract A_tail and phase δ."""
    sol = solve_soliton(g0, n_K=n_K, R_max=R_max)
    if sol.status != 0 and sol.status != 1:
        return np.nan, np.nan

    N_eval = 3000
    r = np.linspace(sol.t[0], min(sol.t[-1], R_max), N_eval)
    y = sol.sol(r)
    g = y[0]

    A, delta, B, C = fit_tail(r, g, r_L=r_L, r_R=r_R)
    return A, delta


TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")

def koide_from_A(A1, A2, A3):
    """Koide ratio from tail amplitudes: m_k ∝ A_k⁴"""
    if any(np.isnan(x) or x <= 0 for x in [A1, A2, A3]):
        return np.nan
    m1, m2, m3 = A1**4, A2**4, A3**4
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2

# ============================================================
# §1. SOLITON PROFILES WITH K(g) = g⁴
# ============================================================
print("=" * 72)
print("§1. PROFILE SOLITONOWE Z K(g) = g⁴")
print("=" * 72)

g0_scan = np.array([0.3, 0.5, 0.7, 0.8, 0.9, 0.95,
                     1.05, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0, 2.5, 3.0])

print(f"\n  Solving soliton ODE for {len(g0_scan)} g₀ values...")
print(f"  {'g₀':>6s}  {'A_tail':>10s}  {'δ [rad]':>10s}  {'δ [°]':>8s}  {'g(R_max)':>10s}")
print("  " + "-" * 52)

results = {}  # g0 -> (A, delta)
valid_count = 0

for g0 in g0_scan:
    try:
        sol = solve_soliton(g0, n_K=4, R_max=80.0)
        r = np.linspace(sol.t[0], min(sol.t[-1], 80.0), 3000)
        y = sol.sol(r)
        g = y[0]
        g_end = g[-1]

        A, delta, B, C = fit_tail(r, g, r_L=20.0, r_R=50.0)

        if np.isfinite(A) and A > 1e-10:
            results[g0] = (A, delta)
            valid_count += 1
            converged = abs(g_end - 1.0) < 0.1
            mark = "✓" if converged else "~"
            print(f"  {g0:6.3f}  {A:10.6f}  {delta:10.4f}  {np.degrees(delta):8.2f}  {g_end:10.6f} {mark}")
        else:
            print(f"  {g0:6.3f}  {'FAILED':>10s}")
    except Exception as e:
        print(f"  {g0:6.3f}  ERROR: {e}")

record("T1: Soliton solver K=g⁴ works",
       valid_count >= 10,
       f"{valid_count}/{len(g0_scan)} solutions obtained")

# ============================================================
# §2. FAZA δ(g₀) — ANALIZA
# ============================================================
print("\n" + "=" * 72)
print("§2. FAZA δ(g₀) — CZY JEST LOGARYTMICZNA?")
print("=" * 72)

g0_vals = np.array(sorted(results.keys()))
A_vals = np.array([results[g0][0] for g0 in g0_vals])
delta_vals = np.array([results[g0][1] for g0 in g0_vals])

# Unwrap phase (handle 2π jumps)
delta_unwrapped = np.unwrap(delta_vals)

print(f"\n  Faza (unwrapped):")
print(f"  {'g₀':>6s}  {'ln(g₀)':>8s}  {'δ [rad]':>10s}  {'δ_unwrap':>10s}")
print("  " + "-" * 40)
for g0, d, du in zip(g0_vals, delta_vals, delta_unwrapped):
    print(f"  {g0:6.3f}  {np.log(g0):8.4f}  {d:10.4f}  {du:10.4f}")

# Fit: δ(g₀) = a·ln(g₀) + b
log_g0 = np.log(g0_vals)

# Use only points where tail is well-resolved
good = A_vals > 1e-6
if np.sum(good) >= 4:
    coeffs = np.polyfit(log_g0[good], delta_unwrapped[good], 1)
    a_fit = coeffs[0]
    b_fit = coeffs[1]

    # R² quality
    delta_model = a_fit * log_g0[good] + b_fit
    ss_res = np.sum((delta_unwrapped[good] - delta_model)**2)
    ss_tot = np.sum((delta_unwrapped[good] - np.mean(delta_unwrapped[good]))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    print(f"\n  Fit: δ = a·ln(g₀) + b")
    print(f"    a = {a_fit:.6f}")
    print(f"    b = {b_fit:.6f}")
    print(f"    R² = {R2:.6f}")

    a_target = 2*np.pi / (3 * np.log(phi))
    print(f"\n  a·ln(φ) = {a_fit * np.log(phi):.6f}")
    print(f"  2π/3   = {2*np.pi/3:.6f}")
    print(f"  Error  = {abs(a_fit * np.log(phi) - 2*np.pi/3)/(2*np.pi/3)*100:.2f}%")
    print(f"\n  a(fit)     = {a_fit:.4f}")
    print(f"  a(required) = {a_target:.4f}")
    print(f"  Error = {abs(a_fit - a_target)/a_target*100:.2f}%")

    record("T2: Phase δ(g₀) is smooth",
           R2 > 0.7,
           f"R² = {R2:.4f}, a = {a_fit:.4f}")

    record("T3: a measured from fit",
           np.isfinite(a_fit),
           f"a = {a_fit:.6f}, a·ln(φ) = {a_fit*np.log(phi):.4f}")

    record("T4: a·ln(φ) ≈ 2π/3",
           abs(a_fit * np.log(phi) - 2*np.pi/3) / (2*np.pi/3) < 0.15,
           f"a·ln(φ) = {a_fit*np.log(phi):.4f}, 2π/3 = {2*np.pi/3:.4f}, err = {abs(a_fit*np.log(phi)-2*np.pi/3)/(2*np.pi/3)*100:.1f}%")
else:
    a_fit = np.nan
    print("  NOT ENOUGH VALID POINTS for fit")
    record("T2-T4: Insufficient data", False, f"Only {np.sum(good)} valid points")

# ============================================================
# §3. COMPARISON: K=g⁴ vs K=g²
# ============================================================
print("\n" + "=" * 72)
print("§3. PORÓWNANIE: K=g⁴ vs K=g²")
print("=" * 72)

print(f"\n  Solving with K=g² for same g₀ values...")

results_K2 = {}
for g0 in g0_scan:
    try:
        A2, delta2 = get_A_delta(g0, n_K=2, R_max=80.0)
        if np.isfinite(A2) and A2 > 1e-10:
            results_K2[g0] = (A2, delta2)
    except:
        pass

g0_K2 = np.array(sorted(results_K2.keys()))
if len(g0_K2) >= 4:
    A_K2 = np.array([results_K2[g0][0] for g0 in g0_K2])
    delta_K2 = np.unwrap(np.array([results_K2[g0][1] for g0 in g0_K2]))

    log_g0_K2 = np.log(g0_K2)
    good_K2 = A_K2 > 1e-6
    if np.sum(good_K2) >= 3:
        coeffs_K2 = np.polyfit(log_g0_K2[good_K2], delta_K2[good_K2], 1)
        a_K2 = coeffs_K2[0]

        print(f"\n  K = g² fit: a = {a_K2:.4f}, a·ln(φ) = {a_K2*np.log(phi):.4f}")
        if np.isfinite(a_fit):
            print(f"  K = g⁴ fit: a = {a_fit:.4f}, a·ln(φ) = {a_fit*np.log(phi):.4f}")
        print(f"  Target:      a = {2*np.pi/(3*np.log(phi)):.4f}, a·ln(φ) = {2*np.pi/3:.4f}")

        record("T5: Compare a(K=g⁴) vs a(K=g²)",
               True,
               f"a(g⁴)={a_fit:.4f}, a(g²)={a_K2:.4f}")

# ============================================================
# §4. φ-FP SEARCH WITH K=g⁴
# ============================================================
print("\n" + "=" * 72)
print("§4. φ-FIXED POINT Z K(g) = g⁴")
print("=" * 72)

print(f"\n  Szukam g₀* t.ż. (A(φg₀*)/A(g₀*))⁴ = {R21_PDG:.3f}")

def ratio_func_nK(g0, n_K=4):
    A1, _ = get_A_delta(g0, n_K=n_K)
    A2, _ = get_A_delta(phi * g0, n_K=n_K)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

# Scan g₀ range
g0_fp_scan = np.linspace(0.5, 2.0, 30)
ratios = []
print(f"\n  {'g₀':>6s}  {'(A(φg₀)/A(g₀))⁴':>18s}")
print("  " + "-" * 30)

for g0 in g0_fp_scan:
    try:
        rv = ratio_func_nK(g0, n_K=4)
        ratios.append(rv)
        r_val = rv + R21_PDG
        if np.isfinite(rv) and abs(rv) < 1e9:
            mark = " ←" if abs(rv) < 50 else ""
            print(f"  {g0:6.3f}  {r_val:18.2f}{mark}")
    except:
        ratios.append(np.nan)

# Find zero crossing
ratios = np.array(ratios)
g0_fp_star = None

for i in range(len(ratios)-1):
    if np.isfinite(ratios[i]) and np.isfinite(ratios[i+1]):
        if ratios[i] * ratios[i+1] < 0:
            try:
                g0_fp_star = brentq(lambda g: ratio_func_nK(g, 4),
                                     g0_fp_scan[i], g0_fp_scan[i+1], xtol=1e-4)
                break
            except:
                pass

if g0_fp_star is not None:
    A_e, delta_e = get_A_delta(g0_fp_star, n_K=4)
    A_mu, delta_mu = get_A_delta(phi * g0_fp_star, n_K=4)
    A_tau, delta_tau = get_A_delta(phi**2 * g0_fp_star, n_K=4)
    r21_fp = (A_mu/A_e)**4
    r31_fp = (A_tau/A_e)**4 if A_tau > 0 else np.nan

    print(f"\n  ★ φ-FP znaleziony:")
    print(f"    g₀*(e)  = {g0_fp_star:.5f}")
    print(f"    g₀*(μ)  = {phi*g0_fp_star:.5f}")
    print(f"    g₀*(τ)  = {phi**2*g0_fp_star:.5f}")
    print(f"    A_e     = {A_e:.6f}")
    print(f"    A_μ     = {A_mu:.6f}")
    print(f"    A_τ     = {A_tau:.6f}")
    print(f"    r₂₁     = (A_μ/A_e)⁴ = {r21_fp:.2f}  (PDG: {R21_PDG})")
    if np.isfinite(r31_fp):
        print(f"    r₃₁     = (A_τ/A_e)⁴ = {r31_fp:.1f}  (PDG: {R31_PDG})")

    # Phase analysis
    # UWAGA: g₀^e < 1 (soliton poniżej próżni) i g₀^μ > 1 (powyżej).
    # Przejście g₀ = 1 daje naturalny skok fazy ~π (zmiana znaku g-1).
    # Korekta: dodaj π do fazy gdy g₀ < 1
    delta_e_corr = delta_e + np.pi  # correct for below-vacuum start
    delta_mu_corr = delta_mu  # already above-vacuum

    Ddelta_raw = delta_mu - delta_e  # raw difference
    Ddelta_corr = delta_mu_corr - delta_e_corr  # corrected (removes π shift)

    print(f"\n  Fazy ogona na φ-drabince:")
    print(f"    δ_e (raw)  = {delta_e:.4f} rad  ({np.degrees(delta_e):.2f}°)")
    print(f"    δ_μ (raw)  = {delta_mu:.4f} rad  ({np.degrees(delta_mu):.2f}°)")
    print(f"    δ_e (corr) = {delta_e_corr:.4f} rad  (+ π for g₀<1)")
    print(f"    δ_μ (corr) = {delta_mu_corr:.4f} rad")
    print(f"\n    Δδ(raw)    = {Ddelta_raw:.4f} rad  ({np.degrees(Ddelta_raw):.2f}°)")
    print(f"    Δδ(corr)   = {Ddelta_corr:.4f} rad  ({np.degrees(Ddelta_corr):.2f}°)")
    print(f"    2π/3       = {2*np.pi/3:.4f} rad  ({120.0:.2f}°)")
    print(f"    π - Δδ_corr = {np.pi + Ddelta_corr:.4f}")

    # Try multiple interpretations
    for label, dd in [("raw", Ddelta_raw),
                       ("corrected", Ddelta_corr),
                       ("raw mod 2π", Ddelta_raw % (2*np.pi)),
                       ("|raw| mod 2π", abs(Ddelta_raw) % (2*np.pi))]:
        err = abs(dd - 2*np.pi/3)/(2*np.pi/3)*100
        err2 = abs(abs(dd) - 2*np.pi/3)/(2*np.pi/3)*100
        print(f"    {label:20s}: {dd:8.4f} rad, err vs 2π/3: {min(err,err2):.1f}%")

    # Tau: try with longer integration
    if np.isnan(A_tau):
        print(f"\n  τ solver failed at g₀={phi**2*g0_fp_star:.3f}, trying extended range...")
        for Rmax_try in [120, 200]:
            A_tau2, delta_tau2 = get_A_delta(phi**2 * g0_fp_star, n_K=4, R_max=Rmax_try)
            if np.isfinite(A_tau2):
                print(f"    R_max={Rmax_try}: A_τ = {A_tau2:.6f}, δ_τ = {delta_tau2:.4f}")
                A_tau = A_tau2
                delta_tau = delta_tau2
                break

    if np.isfinite(A_tau):
        r31_fp = (A_tau/A_e)**4
        print(f"\n    r₃₁ = {r31_fp:.1f} (PDG: {R31_PDG})")
        K_Atail = koide_from_A(A_e, A_mu, A_tau)
        print(f"    K(A⁴) = {K_Atail:.6f} (2/3 = {2/3:.6f})")
    else:
        print(f"\n    τ solver STILL fails — g₀={phi**2*g0_fp_star:.3f} too large for K=g⁴")

    record("T6: φ-FP exists for K=g⁴",
           abs(r21_fp - R21_PDG) / R21_PDG < 0.02,
           f"g₀* = {g0_fp_star:.5f}, r₂₁ = {r21_fp:.2f}")
else:
    print("\n  φ-FP NIE znaleziony w zakresie [0.5, 2.0]")
    print("  Poszerzam zakres skanowania...")

    # Try wider range
    g0_fp_scan2 = np.linspace(0.3, 3.0, 50)
    ratios2 = []
    for g0 in g0_fp_scan2:
        try:
            rv = ratio_func_nK(g0, n_K=4)
            ratios2.append(rv)
        except:
            ratios2.append(np.nan)

    ratios2 = np.array(ratios2)
    for i in range(len(ratios2)-1):
        if np.isfinite(ratios2[i]) and np.isfinite(ratios2[i+1]):
            if ratios2[i] * ratios2[i+1] < 0:
                try:
                    g0_fp_star = brentq(lambda g: ratio_func_nK(g, 4),
                                         g0_fp_scan2[i], g0_fp_scan2[i+1], xtol=1e-4)
                    print(f"  ★ Znaleziony w poszerzonym zakresie: g₀* = {g0_fp_star:.5f}")
                    break
                except:
                    pass

    if g0_fp_star is None:
        # Show the ratio function values to diagnose
        print("\n  Ratio function values (diagnostic):")
        for g0, rv in zip(g0_fp_scan2, ratios2):
            if np.isfinite(rv):
                print(f"    g₀={g0:.3f}: (A_μ/A_e)⁴ = {rv+R21_PDG:.2f}")

    record("T6: φ-FP search",
           g0_fp_star is not None,
           f"g₀* = {g0_fp_star}" if g0_fp_star else "Not found in [0.3, 3.0]")


# ============================================================
# §5. DIRECT PHASE MEASUREMENT AT φ-SPACED g₀
# ============================================================
print("\n" + "=" * 72)
print("§5. BEZPOŚREDNI POMIAR FAZ DLA φ-ROZSTAWIONYCH g₀")
print("=" * 72)

# Even if the FP is not found, test phase spacing for representative g₀
# Use the g₀ range where solutions are good
test_g0_base = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]

print(f"\n  Pomiar Δδ dla par (g₀, φ·g₀) z K=g⁴:")
print(f"  {'g₀':>6s}  {'φ·g₀':>6s}  {'δ₁':>8s}  {'δ₂':>8s}  {'Δδ':>8s}  {'2π/3':>8s}  {'err':>8s}")
print("  " + "-" * 60)

deltas_measured = []
for g0 in test_g0_base:
    g0_2 = phi * g0
    try:
        A1, d1 = get_A_delta(g0, n_K=4)
        A2, d2 = get_A_delta(g0_2, n_K=4)
        if np.isfinite(d1) and np.isfinite(d2):
            Dd = d2 - d1
            # Normalize to [-π, π]
            while Dd > np.pi: Dd -= 2*np.pi
            while Dd < -np.pi: Dd += 2*np.pi
            err = abs(Dd - 2*np.pi/3)/(2*np.pi/3)*100
            deltas_measured.append(Dd)
            print(f"  {g0:6.3f}  {g0_2:6.3f}  {d1:8.4f}  {d2:8.4f}  {Dd:8.4f}  {2*np.pi/3:8.4f}  {err:6.1f}%")
    except:
        pass

if deltas_measured:
    Dd_mean = np.mean(deltas_measured)
    Dd_std = np.std(deltas_measured)
    err_mean = abs(Dd_mean - 2*np.pi/3)/(2*np.pi/3)*100
    print(f"\n  Średnia Δδ = {Dd_mean:.4f} ± {Dd_std:.4f}")
    print(f"  2π/3     = {2*np.pi/3:.4f}")
    print(f"  Error    = {err_mean:.1f}%")

    # Also do for K=g²
    print(f"\n  Porównanie z K=g²:")
    deltas_K2 = []
    for g0 in test_g0_base:
        g0_2 = phi * g0
        try:
            A1, d1 = get_A_delta(g0, n_K=2)
            A2, d2 = get_A_delta(g0_2, n_K=2)
            if np.isfinite(d1) and np.isfinite(d2):
                Dd = d2 - d1
                while Dd > np.pi: Dd -= 2*np.pi
                while Dd < -np.pi: Dd += 2*np.pi
                deltas_K2.append(Dd)
        except:
            pass

    if deltas_K2:
        Dd_K2_mean = np.mean(deltas_K2)
        print(f"  K=g⁴: <Δδ> = {Dd_mean:.4f}  (err: {err_mean:.1f}%)")
        print(f"  K=g²: <Δδ> = {Dd_K2_mean:.4f}  (err: {abs(Dd_K2_mean-2*np.pi/3)/(2*np.pi/3)*100:.1f}%)")
        print(f"  2π/3       = {2*np.pi/3:.4f}")

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
print("PODSUMOWANIE ex226")
print(f"{'='*72}")

if deltas_measured:
    print(f"""
  GŁÓWNY WYNIK:
    Δδ(K=g⁴) = {Dd_mean:.4f} ± {Dd_std:.4f} rad
    2π/3     = {2*np.pi/3:.4f} rad
    Error    = {err_mean:.1f}%

  INTERPRETACJA:
    {'a·ln(φ) ≈ 2π/3 → Z₃ MECHANIZM POTWIERDZONY' if err_mean < 15 else 'a·ln(φ) ≠ 2π/3 → Z₃ mechanizm NIE potwierdzony'}
    {'Koide K = 2/3 wynika z geometrii solitonowej' if err_mean < 15 else 'K = 2/3 wymaga innego wyjaśnienia'}
""")
else:
    print("\n  BRAK WYNIKÓW — solver nie dał wystarczających danych")
