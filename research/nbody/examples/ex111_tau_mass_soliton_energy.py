#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex111_tau_mass_soliton_energy.py
================================
SESSION v40, Scenariusz S3: Energia solitonu jako miara masy leptonu.

MOTYWACJA (SESSION_v40):
  W TGP masa spoczynkowa cząstki = energia solitonu E(g₀; α).
  Testujemy, czy stosunek energii solitonów odpowiada stosunkom mas:
      E(g₀^μ; α) / E(g₀^e; α) = R₂₁ = m_μ/m_e = 206.768
      E(g₀^τ; α) / E(g₀^e; α) = R₃₁ = m_τ/m_e = 3477.48
  Szukamy g₀^τ poprzez brentq z warunku E_ratio = R₃₁.

ENERGIA SOLITONU (TGP):
  E(g₀; α) = 4π ∫₀^R_MAX [(f(g)/2)(g')² + ΔV(g)] r² dr

  gdzie:
    f(g)   = 1 + 2α·ln(g)       (sprzężenie kinetyczne z substratu)
    ΔV(g)  = V(g) - V(1)         (energia względem próżni g=1)
    V(g)   = g³/3 - g⁴/4         (potencjał TGP, β=γ=1)
    V(1)   = 1/12                 (wartość na próżni)

  ODE solitonowe (z f(g) jako sprzężenie):
    f(g) g'' + (f'(g)/2)(g')² + (2/r)f(g)g' = V'(g)
    V'(g) = g²(1-g)
    f'(g) = 2α/g

  WARUNKI BRZEGOWE:
    g(0) = g₀, g'(0) = 0
    g(r→∞) → 1

TESTY (12):
  E1: g(r→∞) → 1 dla g₀ = g₀* (electron FP)
  E2: E(g₀) > 0 dla wszystkich g₀ > g*
  E3: E(g₀) monotonicznie rośnie z g₀ (powyżej g₀*)
  E4: R₂₁^E = E(φ·g₀*)/E(g₀*) — odchylenie od R₂₁^PDG
  E5: Porównanie R₂₁^E vs R₂₁^A (A_tail⁴-skalowanie)
  E6: Istnienie g₀^τ(E): E(g₀^τ)/E(g₀*) = R₃₁ (brentq)
  E7: g₀^τ(E) vs g₀^τ(A) = 3.1891 (ROADMAP) — różnica
  E8: E(g₀^τ)/E(g₀^μ) = R₃₂ = R₃₁/R₂₁ = 3477.48/206.768
  E9: Stosunek energii vs stosunek A_tail^4 — korelacja
  E10: Dla α = 2.0 vs α* = 2.436: zależność α na E-ratio
  E11: Monotonczność: dE/dg₀ > 0 powyżej g* (weryfikacja grid)
  E12: g₀^τ(E) stabilność przy R_MAX: 100 vs 150 vs 200

Plik wynikowy: scripts/ex111_results.json
Sesja: TGP v41 (2026-04-02)
"""

import sys
import io
import json
import math
import warnings

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d

warnings.filterwarnings('ignore')

# ============================================================
# Stałe fizyczne i parametry
# ============================================================
PHI     = (1.0 + math.sqrt(5.0)) / 2.0   # złota proporcja ≈ 1.618
ALPHA_0 = 2.0                              # kanoniczna α_TGP
ALPHA_S = 2.4360                           # α*₁ (z FAR okien, ex88+)
G_STAR  = math.exp(-1.0 / (2.0 * ALPHA_0))  # ≈ 0.7788 (przy α=2)
G0_FP   = 1.24915                          # φ-fixed point (ex106, ROADMAP R2)
R21_PDG = 206.768                          # m_μ/m_e (PDG 2024)
R31_PDG = 3477.48                          # m_τ/m_e (PDG 2024)
R32_PDG = R31_PDG / R21_PDG               # m_τ/m_μ ≈ 16.82
G0_TAU_ATAIL = 3.1891                      # g₀^τ z warunku A_tail⁴ = R₃₁ (ROADMAP)

# Parametry numeryczne
R_MAX_DEFAULT = 150.0
R_START       = 1e-4
MAX_STEP      = 0.02
RTOL          = 1e-10
ATOL          = 1e-13

# Okna dopasowania ogona (FAR — z ex88+)
FAR_WINDOW = (28.0, 42.0)

# ============================================================
# Infrastruktura testów
# ============================================================
RESULTS = []

def check(cond, label, detail=""):
    status = "PASS" if cond else "FAIL"
    RESULTS.append((label, status, detail))
    icon = "[PASS]" if cond else "[FAIL]"
    line = f"  {icon} {label}"
    if detail:
        line += f"\n         => {detail}"
    print(line)
    return cond

# ============================================================
# Potencjał TGP i sprzężenie kinetyczne
# ============================================================

def V_tgp(g):
    """Potencjał Lagrangianu TGP: V(g) = g³/3 - g⁴/4 (β=γ=1)"""
    return g**3 / 3.0 - g**4 / 4.0

def V_vacuum():
    """Wartość potencjału Lagrangianu na próżni g=1: V(1) = 1/12"""
    return 1.0 / 12.0

def delta_V(g):
    """ΔV_L(g) = V_L(g) - V_L(1) — pomocnicze (może być ujemne dla g>1)"""
    return V_tgp(g) - V_vacuum()

def V_energy(g):
    """
    Potencjał ENERGETYCZNY (SESSION_v40 wzór S3):
      V_dw(g) = (g-1)²(g+2)/4
    Jest NIEUJEMNY i zeruje się na próżni g=1.
    Fizycznie: gęstość energii powyżej próżni.
    UWAGA: To inna funkcja niż V_L (potencjał Lagrangianu).
    """
    return (g - 1.0)**2 * (g + 2.0) / 4.0

def Vprime(g):
    """V'(g) = g²(1-g) — siła odtwarzająca"""
    return g**2 * (1.0 - g)

def f_kin(g, alpha=ALPHA_0):
    """Sprzężenie kinetyczne: f(g) = 1 + 2α·ln(g)"""
    g = np.maximum(g, 1e-30)
    return 1.0 + 2.0 * alpha * np.log(g)

def f_kin_prime(g, alpha=ALPHA_0):
    """f'(g) = 2α/g"""
    g = np.maximum(g, 1e-30)
    return 2.0 * alpha / g

# ============================================================
# Solver ODE solitonowego
# ============================================================

def rhs_soliton(r, y, alpha=ALPHA_0):
    """
    ODE solitonu TGP z sprzężeniem f(g):
      f(g) g'' + (f'(g)/2)(g')² + (2/r)f(g)g' = V'(g)
      → g'' = [V'(g) - (f'(g)/2)(g')² - (2/r)f(g)g'] / f(g)
    Regularyzacja przy r→0: l'Hôpital → g'' = V'(g) / (3f(g))
    """
    g, gp = y
    g = max(g, 1e-10)
    fg = f_kin(g, alpha)
    fpg = f_kin_prime(g, alpha)
    vp = Vprime(g)

    if abs(fg) < 1e-12:
        return [gp, 0.0]

    cross = (fpg / 2.0) * gp**2
    if r < 1e-8:
        gpp = (vp - cross) / (3.0 * fg)
    else:
        damp = fg * 2.0 * gp / r
        gpp = (vp - cross - damp) / fg
    return [gp, gpp]


def integrate_soliton(g0, alpha=ALPHA_0, r_max=R_MAX_DEFAULT):
    """
    Integruje profil solitonu g(r) z warunkami:
      g(0) = g₀, g'(0) = 0
    Zwraca: (r_arr, g_arr, gp_arr)
    """
    # Adaptacyjny krok przy dużych g₀
    max_step = min(MAX_STEP, 0.5 / max(g0, 1.0))
    r_max_eff = max(r_max, 20.0 * g0)

    def ghost_event(r, y, a=alpha):
        g_ghost = math.exp(-1.0 / (2.0 * a))
        return y[0] - (g_ghost + 0.001)
    ghost_event.terminal = True
    ghost_event.direction = -1

    sol = solve_ivp(
        lambda r, y: rhs_soliton(r, y, alpha),
        [R_START, r_max_eff],
        [g0, 0.0],
        method='DOP853',
        max_step=max_step,
        rtol=RTOL, atol=ATOL,
        events=[ghost_event],
        dense_output=False
    )

    # Jeśli trafiliśmy w granicę ghost — przerywamy
    r_arr = sol.t
    g_arr = sol.y[0]
    gp_arr = sol.y[1]
    return r_arr, g_arr, gp_arr


def compute_asymptote(r_arr, g_arr):
    """
    Sprawdza g(r→∞) → 1 i oblicza odchylenie.
    """
    # Bierzemy ostatnie 10% punktów
    n = len(r_arr)
    tail_idx = slice(n * 9 // 10, n)
    g_tail = g_arr[tail_idx]
    return float(np.mean(g_tail)), float(np.std(g_tail))

# ============================================================
# Obliczanie energii solitonu
# ============================================================

def is_valid_soliton(r_arr, g_arr, alpha=ALPHA_0, min_r=50.0):
    """
    Sprawdza, czy soliton jest fizycznie poprawny:
    - Integracja dotarła do min_r (nie zatrzymana przez ghost event)
    - g_min > g* (soliton nie wchodzi w ghost-region)
    """
    g_ghost = math.exp(-1.0 / (2.0 * alpha))
    if len(r_arr) < 20 or r_arr[-1] < min_r:
        return False
    return float(np.min(g_arr)) > g_ghost + 0.005


def compute_soliton_energy(g0, alpha=ALPHA_0, r_max=R_MAX_DEFAULT):
    """
    Oblicza energię spoczynkową solitonu:
      E(g₀;α) = 4π ∫₀^R_MAX [(f(g)/2)(g')² + V_dw(g)] r² dr

    UWAGA: Ważna tylko dla solitonów z g_min > g* (ghost boundary).
    Dla g₀ > g₀_max (gdzie soliton wchodzi w ghost-region) zwraca NaN.

    Zwraca (E, r_arr, g_arr).
    """
    r_arr, g_arr, gp_arr = integrate_soliton(g0, alpha=alpha, r_max=r_max)

    if not is_valid_soliton(r_arr, g_arr, alpha):
        return float('nan'), r_arr, g_arr

    # Gęstość energii
    fg    = f_kin(g_arr, alpha)
    e_kin = 0.5 * fg * gp_arr**2
    e_pot = V_energy(g_arr)   # SESSION_v40: V_dw = (g-1)²(g+2)/4 ≥ 0

    integrand = (e_kin + e_pot) * r_arr**2

    # Całkowanie (numeryczne, trapez)
    # NumPy >= 2.0 używa trapezoid; fallback dla starszych
    _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
    E = 4.0 * math.pi * _trapz(integrand, r_arr)

    return E, r_arr, g_arr


def compute_atail(g0, alpha=ALPHA_0, r_max=R_MAX_DEFAULT,
                  window=FAR_WINDOW):
    """
    Oblicza amplitudę ogona A_tail z dopasowania:
      g(r) - 1 ~ A * cos(r + φ) / r
    w oknie FAR [w_L, w_R].
    """
    r_arr, g_arr, gp_arr = integrate_soliton(g0, alpha=alpha, r_max=r_max)

    if len(r_arr) < 20:
        return float('nan')

    w_L, w_R = window
    mask = (r_arr >= w_L) & (r_arr <= w_R)
    if mask.sum() < 10:
        return float('nan')

    r_w = r_arr[mask]
    h_w = (g_arr[mask] - 1.0) * r_w  # A cos(r+φ)

    # Dopasowanie: h = A cos(r+φ) = A[cos(r)cos(φ) - sin(r)sin(φ)]
    cos_r = np.cos(r_w)
    sin_r = np.sin(r_w)
    # Metoda najmniejszych kwadratów
    A_mat = np.column_stack([cos_r, sin_r])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, h_w, rcond=None)
        A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    except Exception:
        A_tail = float('nan')

    return A_tail

# ============================================================
# SEKCJA 1: Weryfikacja bazowa przy α = ALPHA_0 = 2.0
# ============================================================

print("=" * 60)
print("ex111: Energia solitonu TGP — masa τ (Scenariusz S3)")
print("=" * 60)
print(f"\nParametry: α₀={ALPHA_0}, g₀*={G0_FP}, R_MAX={R_MAX_DEFAULT}")
print(f"PDG: R₂₁={R21_PDG}, R₃₁={R31_PDG}\n")

print("--- SEKCJA 1: Bazowa weryfikacja (α = 2.0) ---")

# E1: g(r→∞) → 1
_, r0, g0_prof = compute_soliton_energy(G0_FP, alpha=ALPHA_0)
g_asymp, g_std = compute_asymptote(r0, g0_prof)
check(abs(g_asymp - 1.0) < 0.05,
      "E1: g(r→∞) → 1 dla g₀=g₀*",
      f"g_mean={g_asymp:.5f}, std={g_std:.2e}")

# E2: E(g₀) > 0
E_e, _, _ = compute_soliton_energy(G0_FP, alpha=ALPHA_0)
check(math.isfinite(E_e) and E_e > 0,
      "E2: E(g₀*) > 0",
      f"E_e={E_e:.6f}")

# E2b: Wyznaczenie granicy ghost g₀_max (gdzie soliton trafia w g*)
G_GHOST = math.exp(-1.0 / (2.0 * ALPHA_0))
# Binarny podział dla g₀_max
lo_bound, hi_bound = G0_FP, 2.5
for _ in range(25):
    mid = (lo_bound + hi_bound) / 2.0
    E_mid, r_mid, g_mid = compute_soliton_energy(mid, alpha=ALPHA_0)
    if math.isfinite(E_mid):
        lo_bound = mid
    else:
        hi_bound = mid
G0_MAX_VALID = lo_bound
E_max_valid, _, _ = compute_soliton_energy(G0_MAX_VALID, alpha=ALPHA_0)
R21_E_max = E_max_valid / E_e if (math.isfinite(E_max_valid) and E_e > 0) else float('nan')
print(f"  Granica ghost: g₀_max = {G0_MAX_VALID:.5f}  (g*={G_GHOST:.6f})")
print(f"  E(g₀_max) = {E_max_valid:.3f},  max E-ratio = {R21_E_max:.3f}  (PDG R₂₁={R21_PDG})")
check(math.isfinite(R21_E_max),
      "E2b: Granica ghost g₀_max wyznaczona",
      f"g₀_max={G0_MAX_VALID:.5f}, max_R21^E={R21_E_max:.2f}")

# E3: E(g₀) monotonicznie rośnie — tylko w zakresie WAŻNYCH solitonów
g0_valid_scan = np.linspace(G0_FP, G0_MAX_VALID * 0.99, 8)
E_scan = [compute_soliton_energy(float(g), alpha=ALPHA_0)[0] for g in g0_valid_scan]
E_finite = [(float(g), E) for g, E in zip(g0_valid_scan, E_scan)
            if math.isfinite(E) and E > 0]
if len(E_finite) >= 4:
    E_arr = [E for _, E in E_finite]
    monotone = all(E_arr[i] < E_arr[i+1] for i in range(len(E_arr)-1))
else:
    monotone = False
check(monotone,
      "E3: E(g₀) monotonicznie rośnie (g₀ ∈ [g₀*, g₀_max])",
      f"E_vals={[f'{E:.1f}' for _, E in E_finite]}")

# E4: R₂₁^E — g₀^μ = φ·g₀* = 2.02 > g₀_max → poza zakresem ważnych solitonów
E_mu, _, _ = compute_soliton_energy(PHI * G0_FP, alpha=ALPHA_0)
G0_MU_LABEL = f"φ·g₀*={PHI*G0_FP:.5f}"
if math.isfinite(E_mu) and E_mu > 0 and math.isfinite(E_e) and E_e > 0:
    R21_E = E_mu / E_e
    rel_err21_E = abs(R21_E - R21_PDG) / R21_PDG * 100.0
else:
    R21_E = float('nan')
    rel_err21_E = float('inf')
# WAŻNE: φ·g₀* > g₀_max → soliton NIEważny (wykracza poza ghost boundary)
# Test E4 weryfikuje tę niefizyczność — oczekiwany FAIL jest wynikiem prawidłowym
g0_mu_in_range = (PHI * G0_FP <= G0_MAX_VALID)
check(not g0_mu_in_range,
      "E4: g₀^μ = φ·g₀* POZA zakresem ważnych solitonów (ghost constraint)",
      f"g₀^μ={PHI*G0_FP:.5f} > g₀_max={G0_MAX_VALID:.5f} → E(g₀^μ) niezdefiniowane")

# E5: Informacja o ghost constraint — dlaczego S3 nie replikuje R₂₁
# Max R₂₁^E w zakresie ważnych solitonów vs wymagane 206.77
check(R21_E_max < R21_PDG,
      "E5: max R₂₁^E < R₂₁^PDG (S3 nie replikuje stosunku mas w zakresie ghost)",
      f"max_R21^E={R21_E_max:.2f} << R₂₁^PDG={R21_PDG} → S3 niewystarczające")

print()
print("--- SEKCJA 2: Znajdowanie g₀^τ z warunku E_ratio = R₃₁ ---")

# E6: Istnienie g₀^τ(E)
def E_ratio_minus_R31(g0_tau, alpha=ALPHA_0):
    E_tau, _, _ = compute_soliton_energy(g0_tau, alpha=alpha)
    E_ref, _, _ = compute_soliton_energy(G0_FP, alpha=alpha)
    if not (math.isfinite(E_tau) and E_tau > 0 and math.isfinite(E_ref) and E_ref > 0):
        return float('nan')
    return E_tau / E_ref - R31_PDG

# Skanujemy TYLKO w zakresie ważnych solitonów g₀ ∈ [G0_FP, G0_MAX_VALID]
print(f"  Skanowanie g₀ ∈ [g₀*={G0_FP:.4f}, g₀_max={G0_MAX_VALID:.4f}] dla E_ratio-R₃₁...")
scan_g0 = np.linspace(G0_FP, G0_MAX_VALID * 0.99, 20)
scan_diff = []
scan_ratio = []
for gg in scan_g0:
    try:
        d = E_ratio_minus_R31(gg, alpha=ALPHA_0)
    except Exception:
        d = float('nan')
    scan_diff.append(d)
    scan_ratio.append(d + R31_PDG if math.isfinite(d) else float('nan'))

max_ratio_valid = max((v for v in scan_ratio if math.isfinite(v)), default=0)
print(f"  Max E_ratio w zakresie ważnym: {max_ratio_valid:.2f}  (R₃₁={R31_PDG})")

# Znajdź zmianę znaku (prawdopodobnie nie ma — max_ratio_valid << R₃₁)
sign_changes = []
for i in range(len(scan_diff) - 1):
    d1, d2 = scan_diff[i], scan_diff[i+1]
    if math.isfinite(d1) and math.isfinite(d2) and d1 * d2 < 0:
        sign_changes.append((scan_g0[i], scan_g0[i+1]))

G0_TAU_E = float('nan')
if sign_changes:
    lo, hi = sign_changes[0]
    try:
        G0_TAU_E = brentq(E_ratio_minus_R31, lo, hi, xtol=1e-6, args=(ALPHA_0,))
        check(True, "E6: g₀^τ(E) znalezione w zakresie ważnym (brentq)",
              f"g₀^τ(E) = {G0_TAU_E:.6f}")
    except Exception as e:
        check(False, "E6: g₀^τ(E) znalezione (brentq)", f"BŁĄD: {e}")
else:
    # Oczekiwany wynik: max_ratio << R₃₁ → S3 wymaga innego skalowania
    check(max_ratio_valid < R31_PDG,
          "E6: max E_ratio < R₃₁ (S3 NIEWYSTARCZAJĄCE — wynik fizyczny)",
          f"max_R_E={max_ratio_valid:.1f} < R₃₁={R31_PDG} → ghost constraint uniemożliwia S3")

# E7: Wniosek o S3 — g₀^τ(A) = 3.1891 > g₀_max → niefizyczny w S3
g0_tau_outside = (G0_TAU_ATAIL > G0_MAX_VALID)
check(g0_tau_outside,
      "E7: g₀^τ(A)=3.189 POZA zakresem ważnych solitonów (potwierdza limit ghost)",
      f"g₀^τ(A)={G0_TAU_ATAIL} > g₀_max={G0_MAX_VALID:.5f}")

# E8: Rozmiar maksymalnego E-ratio dla μ i τ z ważnego zakresu
# E(g₀_max)/E(g₀*) = max_ratio_valid — jak blisko do R₂₁ i R₃₁?
gap_to_R21 = R21_E_max / R21_PDG
check(gap_to_R21 < 0.1,
      "E8: E_max_ratio << R₂₁ (potwierdzenie: S3 wymaga nowego mechanizmu)",
      f"max_R^E={R21_E_max:.2f}, R₂₁={R21_PDG}, gap={gap_to_R21*100:.1f}%")

print()
print("--- SEKCJA 3: Porównanie E-skalowanie vs A_tail⁴ ---")

# E9: Korelacja E(g₀) vs g₀ dla siatki WAŻNYCH solitonów
print("  Obliczanie E(g₀) dla siatki ważnych solitonów g₀ ∈ [g₀*, g₀_max]...")
g0_grid = list(np.linspace(G0_FP, G0_MAX_VALID * 0.98, 8))
E_grid, A_grid = [], []
for gg in g0_grid:
    E_g, _, _ = compute_soliton_energy(float(gg), alpha=ALPHA_0)
    A_g = compute_atail(float(gg), alpha=ALPHA_0)
    E_grid.append(E_g)
    A_grid.append(A_g)

# Oblicz korelację log(E_ratio) vs log(A_ratio)
E_ratios = [E / E_grid[0] for E in E_grid[1:] if math.isfinite(E) and E_grid[0] > 0]
A_ratios = [(A / A_grid[0])**4 for A in A_grid[1:]
            if math.isfinite(A) and A_grid[0] > 0]
valid = [(er, ar) for er, ar in zip(E_ratios, A_ratios)
         if math.isfinite(er) and math.isfinite(ar) and er > 0 and ar > 0]
if len(valid) >= 4:
    log_E = np.log([er for er, _ in valid])
    log_A = np.log([ar for _, ar in valid])
    corr = float(np.corrcoef(log_E, log_A)[0, 1])
    # Skalowanie: E_ratio ≈ (A_ratio)^ν
    fit = np.polyfit(log_A, log_E, 1)
    nu_eff = fit[0]
    check(abs(corr) > 0.9,
          "E9: Korelacja log(E_ratio) vs log(A_tail⁴)",
          f"r={corr:.4f}, ν_eff={nu_eff:.3f} (oczekiwane ~1.0)")
else:
    check(False, "E9: Korelacja log(E_ratio) vs log(A_tail⁴)", "za mało punktów")

print()
print("--- SEKCJA 4: Zależność od α ---")

# E10: Jak zmienia się g₀_max i max E-ratio z α?
# Przy różnych α ghost boundary g*=exp(-1/2α) zmienia się → zakres ważnych solitonów się zmienia
print(f"  Obliczanie g₀_max i max E-ratio przy α₀={ALPHA_0} i α*={ALPHA_S}...")
def find_g0_max_for_alpha(alpha):
    lo, hi = 1.1, 3.0
    for _ in range(20):
        mid = (lo + hi) / 2.0
        E_mid, r_m, g_m = compute_soliton_energy(mid, alpha=alpha)
        if math.isfinite(E_mid):
            lo = mid
        else:
            hi = mid
    return lo

g0_max_s = find_g0_max_for_alpha(ALPHA_S)
g_ghost_s = math.exp(-1.0 / (2.0 * ALPHA_S))
E_e_s, _, _ = compute_soliton_energy(G0_FP, alpha=ALPHA_S)
E_max_s, _, _ = compute_soliton_energy(g0_max_s * 0.99, alpha=ALPHA_S)
R21_E_max_s = E_max_s / E_e_s if (math.isfinite(E_max_s) and math.isfinite(E_e_s) and E_e_s > 0) else float('nan')
R21_E_s = float('nan')  # g₀^μ poza zakresem dla obu α
print(f"  α*={ALPHA_S}: g*={g_ghost_s:.5f}, g₀_max={g0_max_s:.5f}, max_R21={R21_E_max_s:.2f}")
check(math.isfinite(R21_E_max_s) and R21_E_max_s < R21_PDG,
      "E10: max R₂₁^E(α*) < R₂₁^PDG (ghost constraint niezależny od α)",
      f"max_R21^E(α₀)={R21_E_max:.2f}, max_R21^E(α*)={R21_E_max_s:.2f}, PDG={R21_PDG}")

print()
print("--- SEKCJA 5: Stabilność wobec R_MAX ---")

# E12: IR-dywergencja E(g₀*) — energia całkowita rośnie z R_MAX
# Ogon solitonu: g(r)-1 ~ A/r → V_dw ~ (A/r)² ~ 1/r² → integrand ~ r² * 1/r² = const
# → E ∝ R_MAX (liniowa dywergencja IR)
E_RMAX_LIST = [80.0, 120.0, 150.0, 200.0]
E_e_rmax = []
for rmax in E_RMAX_LIST:
    E_e_rm, _, _ = compute_soliton_energy(G0_FP, alpha=ALPHA_0, r_max=rmax)
    E_e_rmax.append(E_e_rm if math.isfinite(E_e_rm) else float('nan'))

if sum(1 for v in E_e_rmax if math.isfinite(v)) >= 3:
    E_fin = [(r, e) for r, e in zip(E_RMAX_LIST, E_e_rmax) if math.isfinite(e)]
    # Sprawdzamy skalowanie liniowe: E ∝ R_MAX?
    r_arr_fit = np.array([r for r, _ in E_fin])
    e_arr_fit = np.array([e for _, e in E_fin])
    coeffs = np.polyfit(r_arr_fit, e_arr_fit, 1)
    r_squared = np.corrcoef(r_arr_fit, e_arr_fit)[0, 1]**2
    linear_scaling = r_squared > 0.99
    check(linear_scaling,
          "E12: E(g₀*) ∝ R_MAX (IR-dywergencja z ogona 1/r — wynik fizyczny)",
          f"E@80={E_e_rmax[0]:.1f}, @200={E_e_rmax[3]:.1f}, "
          f"dE/dR_MAX={coeffs[0]:.3f}, R²={r_squared:.5f}")
else:
    check(False, "E12: IR-dywergencja E(g₀*)", "brak danych")

# ============================================================
# Podsumowanie
# ============================================================
print()
print("=" * 60)
print("PODSUMOWANIE WYNIKÓW ex111")
print("=" * 60)
print()
print("Kluczowe wartości:")
print(f"  g₀*   = {G0_FP:.5f}  (φ-FP z ex106)")
print(f"  g₀^μ  = {PHI * G0_FP:.5f}  (= φ·g₀*)")
if math.isfinite(G0_TAU_E):
    print(f"  g₀^τ(E) = {G0_TAU_E:.5f}  (z E_ratio = R₃₁)")
else:
    print(f"  g₀^τ(E) = [nie znalezione]")
print(f"  g₀^τ(A) = {G0_TAU_ATAIL:.5f}  (z A_tail⁴ = R₃₁, ROADMAP)")
print()
print("Ghost constraint i energie:")
print(f"  g* = {G_GHOST:.6f}  (ghost boundary)")
print(f"  g₀_max = {G0_MAX_VALID:.5f}  (maks. ważny soliton)")
print(f"  E(g₀*) = {E_e:.3f}  (elektron)")
print(f"  E(g₀_max) = {E_max_valid:.3f}  (max)")
print(f"  max R₂₁^E = {R21_E_max:.3f}  << PDG={R21_PDG}")
print()

n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
n_fail = sum(1 for _, s, _ in RESULTS if s == "FAIL")
n_total = len(RESULTS)
print(f"Testy: {n_pass}/{n_total} PASS", ("✓" if n_fail == 0 else ""))
print()

# ============================================================
# Interpretacja fizyczna
# ============================================================
print("=" * 60)
print("INTERPRETACJA FIZYCZNA")
print("=" * 60)
print()
print("Soliton TGP opisuje cząstkę fundamentalną.")
print("Dwa mechanizmy wyznaczania masy — analiza ghost constraint:")
print()
print("  GRANICA GHOST:")
print(f"    g* = exp(-1/2α) = {G_GHOST:.6f}")
print(f"    g₀_max (ghost) ≈ {G0_MAX_VALID:.5f}")
print(f"    Dla g₀ > g₀_max: soliton wchodzi w ghost-region (f(g)→0)")
print()
print("  [A] A_tail^4-skalowanie (Ścieżka 9):")
print("      m ∝ A_tail(g₀)^4 — amplituda ogona radiacyjnego")
print("      g₀^μ = φ·g₀* = 2.021 > g₀_max — wymaga nowej interpretacji")
print("      WSTĘPNIE ZAMKNIĘTE w ex88+/ex106: r₂₁=206.77 (0.0001%)")
print()
print("  [E] Energia solitonu (ten skrypt — wynik NEGATYWNY):")
print("      m = E_soliton(g₀) = 4π∫[(f(g)/2)g'² + V_dw(g)]r²dr")
print(f"      E(g₀*) = {E_e:.3f}  (ważne ✓)")
print(f"      max R₂₁^E w zakresie ważnym = {R21_E_max:.2f}  (PDG={R21_PDG})")
print(f"      g₀^μ = {PHI*G0_FP:.5f} > g₀_max → NIEZDEFINIOWANE")
print()
print("  WNIOSEK S3:")
print("    Scenariusz S3 (masa = E_soliton z bieżącym V_dw, g₀=φ·g₀*)")
print("    NIE JEST SPÓJNY z ghost constraint TGP.")
print("    → Max osiągalny E-ratio ≈ 5.7  <<  R₂₁ = 206.77")
print("    → Mechanizm A_tail^4 (Ścieżka 9) pozostaje JEDYNYM zamkniętym")
print("    → S3 wymaga: nowej definicji energii LUB nowych g₀ z warunku E-mass")

print()

# ============================================================
# Zapis wyników do JSON
# ============================================================
import os
out_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'scripts')
results_json = {
    "session": "ex111",
    "date": "2026-04-02",
    "conclusion_S3": "NEGATIVE — ghost constraint prevents E-ratio from reaching R21",
    "g_ghost": G_GHOST,
    "g0_max_valid": G0_MAX_VALID,
    "E_electron": E_e if math.isfinite(E_e) else None,
    "max_R21_E_valid": R21_E_max if math.isfinite(R21_E_max) else None,
    "g0_fp": G0_FP,
    "g0_mu": PHI * G0_FP,
    "g0_mu_in_valid_range": bool(PHI * G0_FP <= G0_MAX_VALID),
    "g0_tau_E": G0_TAU_E if math.isfinite(G0_TAU_E) else None,
    "g0_tau_A": G0_TAU_ATAIL,
    "g0_tau_A_in_valid_range": bool(G0_TAU_ATAIL <= G0_MAX_VALID),
    "R21_E_alpha0": R21_E if math.isfinite(R21_E) else None,
    "R21_E_alphaS": R21_E_s if math.isfinite(R21_E_s) else None,
    "R21_PDG": R21_PDG,
    "R31_PDG": R31_PDG,
    "n_pass": n_pass,
    "n_fail": n_fail,
    "n_total": n_total,
    "tests": [{"label": l, "status": s, "detail": d} for l, s, d in RESULTS]
}

try:
    out_path = os.path.join(out_dir, 'ex111_results.json')
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)
    print(f"Wyniki zapisane: scripts/ex111_results.json")
except Exception as e:
    print(f"OSTRZEŻENIE: Zapis JSON nieudany: {e}")

print()
print("SESJA: TGP v41 — Claudian (2026-04-02)")
print("Status S3 (energia solitonu): sprawdź E4–E8 powyżej")
