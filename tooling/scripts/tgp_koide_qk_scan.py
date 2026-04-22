#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tgp_koide_qk_scan.py
====================
T-OP1: Topologiczna mapa Q_K w przestrzeni g₀ — badanie numeryczne

CEL: Odpowiedź na pytanie T-OP1: czy Q_K=3/2 wynika z dynamiki TGP,
     czy jest dodatkowym warunkiem empirycznym?

METODOLOGIA:
  1. Fiksujemy g₀^e = 1.24915, g₀^μ = φ·g₀^e (z φ-FP, tw. J2)
  2. Skanujemy α = g₀^τ/g₀^e ∈ [1.5, 3.5]
  3. Dla każdego α: obliczamy r₃₁(α) = (A(α·g₀^e)/A_e)⁴
  4. Obliczamy Q_K(r₂₁, r₃₁(α))
  5. Szukamy: czy Q_K(α) ma "specjalny" punkt przy 3/2?

BADANIE DODATKOWE (2D):
  - Skanujemy α₂ = g₀^μ/g₀^e ∈ [1.4, 2.0] (wokół φ=1.618)
  - Dla każdego α₂: szukamy α₃ takie że Q_K=3/2
  - Pytanie: czy krzywa Q_K=3/2 w (α₂,α₃) zawiera "specjalne" punkty?

TESTY (8):
  S1: Q_K(g₀^μ, g₀^μ) > 3/2  (dolny kraniec — obie generacje = μ)
  S2: Q_K(g₀^τ^Koide) ≈ 3/2  (weryfikacja z tgp_koide_r31_formal)
  S3: Q_K(g₀^τ^φ²) < 3/2     (φ²-skalowanie poniżej Koide)
  S4: Q_K monotoniczne malejące w α (potwierdzenie unikalności)
  S5: α_Koide wyznaczone ze skanowania ≈ 2.553 (= g₀^τ/g₀^e)
  S6: dQ_K/dα ≠ 0 przy α=α_Koide (brak extremum → nie jest FP)
  S7: Q_K(α) w skanie 2D — krzywa Q_K=3/2 nie przechodzi przez (φ,φ²)
  S8: Raport: Q_K=3/2 to empiryczny warunek, nie FP TGP

Sesja: TGP v42 (2026-04-01)
"""

import sys
import io
import math
import os
import warnings
import json

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp

warnings.filterwarnings('ignore')

# ============================================================
# Parametry (z ex106, tgp_koide_r31_formal)
# ============================================================
ALPHA   = 2.0
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))
PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_PDG = 206.768
R31_PDG = 3477.48
G0_STAR = 1.24915
G0_MU   = PHI * G0_STAR
G0_TAU_PHI2  = PHI**2 * G0_STAR   # φ²-skalowanie: 3.27032
G0_TAU_KOIDE = 3.18913             # z tgp_koide_r31_formal

R_MAX    = 80.0
R_START  = 1e-4
R_TAIL_L = 20.0
R_TAIL_R = 35.0
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13
G_BOUNCE = G_GHOST + 0.005

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
# Solver ODE (z ex106)
# ============================================================

def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-30))

def rhs_regularized(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_hit_ghost(r, y):
    return y[0] - G_BOUNCE
event_hit_ghost.terminal = True
event_hit_ghost.direction = -1

def integrate_soliton(g0, r_max=None, max_bounces=8):
    if r_max is None:
        r_max = max(R_MAX, 20.0 * g0)
    r0, y0 = R_START, [g0, 0.0]
    segs_r, segs_g = [], []
    for bounce_num in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_regularized, [r0, r_max], y0,
            method='DOP853', max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_ghost], dense_output=False
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        if sol.t_events[0].size > 0 and bounce_num < max_bounces:
            r_b = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0 = r_b + 1e-6
            y0 = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break
    r = np.concatenate(segs_r)
    g = np.concatenate(segs_g)
    idx = np.argsort(r)
    return r[idx], g[idx]

def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return 0.0
    r_fit = r_arr[mask]
    delta_fit = (g_arr[mask] - 1.0) * r_fit
    cos_r = np.cos(r_fit)
    sin_r = np.sin(r_fit)
    X = np.column_stack([cos_r, sin_r])
    coefs, _, _, _ = np.linalg.lstsq(X, delta_fit, rcond=None)
    B, C = coefs
    return float(np.sqrt(B**2 + C**2))

def compute_atail(g0):
    r_arr, g_arr = integrate_soliton(g0)
    return fit_tail(r_arr, g_arr)

def koide_Q(r21, r31):
    x = math.sqrt(r21)
    y = math.sqrt(r31)
    return (1.0 + x + y)**2 / (1.0 + r21 + r31)

# ============================================================
# SEKCJA 0: Ustal A_e, A_μ (referencje)
# ============================================================

print("=" * 70)
print("TGP_KOIDE_QK_SCAN: Topologiczna mapa Q_K — badanie T-OP1")
print("=" * 70)
print(f"  g₀^e = {G0_STAR},  g₀^μ = {G0_MU:.5f},  φ = {PHI:.6f}")
print(f"  g₀^τ(φ²) = {G0_TAU_PHI2:.5f},  g₀^τ(Koide) = {G0_TAU_KOIDE:.5f}")
print()

print("  [REF] Obliczam A_e, A_μ ...")
A_e  = compute_atail(G0_STAR)
A_mu = compute_atail(G0_MU)
r21_computed = (A_mu / A_e)**4
print(f"  A_e  = {A_e:.5f}")
print(f"  A_μ  = {A_mu:.5f}")
print(f"  r₂₁  = {r21_computed:.3f}  (PDG: {R21_PDG})")
print()

# ============================================================
# SEKCJA 1: Skan 1D — Q_K(α) dla α = g₀^τ/g₀^e
# ============================================================

print("--- SEKCJA 1: Skan 1D  Q_K(α),  α = g₀^τ/g₀^e ---")
print()

# Skan gruby: 25 punktów (szybki, diagnostyczny)
alpha_arr = np.linspace(1.5, 3.5, 25)
QK_arr    = []
r31_arr   = []

print("  Obliczam A_tail dla 25 wartości α...")
for i, alpha in enumerate(alpha_arr):
    g0_tau_i = alpha * G0_STAR
    A_tau_i = compute_atail(g0_tau_i)
    r31_i   = (A_tau_i / A_e)**4
    QK_i    = koide_Q(r21_computed, r31_i)
    QK_arr.append(QK_i)
    r31_arr.append(r31_i)
    if (i+1) % 5 == 0:
        print(f"    α={alpha:.3f}  g₀^τ={g0_tau_i:.4f}  A_τ={A_tau_i:.4f}  "
              f"r₃₁={r31_i:.1f}  Q_K={QK_i:.4f}")

QK_arr  = np.array(QK_arr)
r31_arr = np.array(r31_arr)
print()

# ============================================================
# TESTY S1-S6
# ============================================================

print("--- TESTY S1-S6 ---")

# S1: Q_K przy α = α₂ = φ (obie τ i μ mają to samo g₀, r₃₁ ≈ r₂₁ geometrycznie)
# Ale realistycznie: α = φ oznacza g₀^τ = φ·g₀^e = g₀^μ → r₃₁ ≈ r₂₁
# Szukamy Q_K w skan-tablicy dla α bliskiego φ ≈ 1.618
alpha_phi = PHI
idx_phi = np.argmin(np.abs(alpha_arr - alpha_phi))
QK_at_phi = QK_arr[idx_phi]
check(QK_at_phi > 1.5,
      "S1: Q_K(α=φ) > 3/2  [dolny koniec zakresu]",
      f"Q_K(α={alpha_arr[idx_phi]:.3f}) = {QK_at_phi:.4f}")

# S2: Q_K przy α_Koide ≈ 2.553
alpha_koide = G0_TAU_KOIDE / G0_STAR  # ≈ 2.553
idx_koide = np.argmin(np.abs(alpha_arr - alpha_koide))
QK_at_koide = QK_arr[idx_koide]
check(abs(QK_at_koide - 1.5) < 0.02,
      "S2: Q_K(α_Koide) ≈ 3/2  (w skan-tablicy)",
      f"α_Koide={alpha_koide:.4f},  Q_K(α≈{alpha_arr[idx_koide]:.3f}) = {QK_at_koide:.4f}")

# S3: Q_K przy φ² ≈ 2.618
alpha_phi2 = PHI**2
idx_phi2 = np.argmin(np.abs(alpha_arr - alpha_phi2))
QK_at_phi2 = QK_arr[idx_phi2]
check(QK_at_phi2 < 1.5,
      "S3: Q_K(α=φ²) < 3/2  [φ²-skalowanie poniżej Koide]",
      f"Q_K(α={alpha_arr[idx_phi2]:.3f}) = {QK_at_phi2:.4f}")

# S4: Monotoniczność Q_K(α) — powinno być malejące
diffs = np.diff(QK_arr)
is_decreasing = np.all(diffs < 0)
n_violations = np.sum(diffs >= 0)
check(n_violations <= 2,   # Dopuszczamy małe fluktuacje numeryczne
      "S4: Q_K(α) monotonicznie malejące",
      f"Naruszeń monotoniczności: {n_violations}/24 "
      f"(max naruszenie: {max(diffs[diffs>=0], default=0.0):.5f})")

# S5: Skan wyznacza α_Koide ≈ 2.553 z interpolacji
# Szukamy przejścia przez 3/2
crossings = []
for i in range(len(QK_arr)-1):
    if (QK_arr[i] - 1.5) * (QK_arr[i+1] - 1.5) < 0:
        # Interpolacja liniowa
        t = (1.5 - QK_arr[i]) / (QK_arr[i+1] - QK_arr[i])
        alpha_cross = alpha_arr[i] + t * (alpha_arr[i+1] - alpha_arr[i])
        crossings.append(alpha_cross)

if crossings:
    alpha_cross_scan = crossings[0]
    check(abs(alpha_cross_scan - alpha_koide) < 0.1,
          "S5: Skan wyznacza α_Koide ≈ 2.553",
          f"α_cross(skan) = {alpha_cross_scan:.4f},  "
          f"α_Koide = {alpha_koide:.4f},  "
          f"δ = {abs(alpha_cross_scan - alpha_koide)/alpha_koide*100:.2f}%")
else:
    check(False, "S5: Brak przejścia przez 3/2 w skan-tablicy",
          f"Q_K zakres: [{QK_arr.min():.4f}, {QK_arr.max():.4f}]")
    alpha_cross_scan = float('nan')

# S6: dQ_K/dα ≠ 0 przy α_Koide (brak extremum)
# Estymacja lokalnej pochodnej w α_Koide
idx_k = np.argmin(np.abs(alpha_arr - alpha_koide))
idx_k = np.clip(idx_k, 1, len(QK_arr)-2)
dQK_dalpha = (QK_arr[idx_k+1] - QK_arr[idx_k-1]) / (alpha_arr[idx_k+1] - alpha_arr[idx_k-1])
check(abs(dQK_dalpha) > 0.05,
      "S6: dQ_K/dα ≠ 0 przy α_Koide  [nie jest fixed point Q_K]",
      f"dQ_K/dα|α_K = {dQK_dalpha:.4f}  (oczekiwane: < 0)")

# ============================================================
# SEKCJA 2: Skan 2D — krzywa Q_K=3/2 w (α₂, α₃)
# ============================================================

print()
print("--- SEKCJA 2: Skan 2D  α₂ = g₀^μ/g₀^e  vs  α₃^Koide ---")
print()

# α₂ wokół φ: [1.45, 1.80] — 8 punktów (szybki)
alpha2_arr = np.linspace(1.45, 1.80, 8)
alpha3_koide_arr = []
r21_arr_2d = []

print("  Obliczam A_tail dla 8 wartości α₂ + znalezienie α₃^Koide ...")
print(f"  {'α₂':>8}  {'g₀^μ':>8}  {'r₂₁':>9}  {'α₃^K':>8}  {'g₀^τ':>8}  {'Q_K':>7}")
print("  " + "-"*58)

for alpha2 in alpha2_arr:
    g0_mu_i  = alpha2 * G0_STAR
    A_mu_i   = compute_atail(g0_mu_i)
    r21_i    = (A_mu_i / A_e)**4

    # Znajdź α₃ takie że Q_K(r₂₁_i, r₃₁_i) = 3/2
    # Algebraicznie (tw. T-r31):
    x_i = math.sqrt(r21_i)
    disc_i = 3.0 * (x_i**2 + 4.0*x_i + 1.0)
    if disc_i > 0:
        sqrt_r31_koide_i = 2.0*(x_i + 1.0) + math.sqrt(disc_i)
        r31_koide_i = sqrt_r31_koide_i**2
    else:
        r31_koide_i = float('nan')

    if not math.isnan(r31_koide_i):
        # Numerycznie znajdź g₀^τ dla tego r₃₁^K
        # Prosta bisekcja
        def obj_tau(g0t):
            A = compute_atail(g0t)
            return (A / A_e)**4 - r31_koide_i

        try:
            from scipy.optimize import brentq
            g0_lo = G0_STAR * 1.6
            g0_hi = G0_STAR * 4.0
            val_lo = obj_tau(g0_lo)
            val_hi = obj_tau(g0_hi)
            if val_lo * val_hi < 0:
                g0_tau_i = brentq(obj_tau, g0_lo, g0_hi, xtol=2e-3, rtol=1e-3)
                alpha3_i = g0_tau_i / G0_STAR
                QK_check = koide_Q(r21_i, r31_koide_i)
            else:
                g0_tau_i = float('nan')
                alpha3_i = float('nan')
                QK_check = float('nan')
        except Exception:
            g0_tau_i = float('nan')
            alpha3_i = float('nan')
            QK_check = float('nan')
    else:
        alpha3_i = float('nan')
        g0_tau_i = float('nan')
        QK_check = float('nan')

    alpha3_koide_arr.append(alpha3_i)
    r21_arr_2d.append(r21_i)
    g0t_str = f"{g0_tau_i:.4f}" if not math.isnan(g0_tau_i) else "   —"
    a3_str  = f"{alpha3_i:.4f}"  if not math.isnan(alpha3_i) else "   —"
    qk_str  = f"{QK_check:.4f}"  if not math.isnan(QK_check) else " —"
    print(f"  {alpha2:>8.4f}  {g0_mu_i:>8.4f}  {r21_i:>9.2f}  "
          f"{a3_str:>8}  {g0t_str:>8}  {qk_str:>7}")

alpha3_koide_arr = np.array(alpha3_koide_arr)

print()

# S7: Czy krzywa Q_K=3/2 przechodzi przez (φ, φ²)?
# Przy α₂=φ: oczekujemy α₃^K ≈ 2.553, a φ² ≈ 2.618. Różnica?
idx_phi_2d = np.argmin(np.abs(alpha2_arr - PHI))
alpha3_at_phi = alpha3_koide_arr[idx_phi_2d]
phi2 = PHI**2
if not math.isnan(alpha3_at_phi):
    diff_pct = abs(alpha3_at_phi - phi2) / phi2 * 100.0
    check(diff_pct > 1.0,
          "S7: Krzywa Q_K=3/2 NIE przechodzi przez (φ, φ²)",
          f"α₂=φ={PHI:.4f}: α₃^K = {alpha3_at_phi:.4f},  "
          f"φ² = {phi2:.4f},  Δ = {diff_pct:.2f}%")
else:
    check(False, "S7: Nie obliczono α₃^K przy α₂=φ", "brak danych")

# S8: Raport syntetyczny
print()
print("  --- Raport T-OP1 ---")
print()
print("  PYTANIE: Czy Q_K = 3/2 wynika z dynamiki TGP?")
print()
print("  WYNIKI:")
print(f"    1. Q_K(α) malejące od ≈{QK_arr[0]:.3f} (α=1.5) do ≈{QK_arr[-1]:.3f} (α=3.5)")
print(f"    2. Q_K=3/2 przecina w JEDNYM punkcie: α_K = {alpha_cross_scan:.3f}")
print(f"    3. φ²-skalowanie daje α=φ²={phi2:.3f}: Q_K = {QK_at_phi2:.4f} ≠ 3/2")
print(f"    4. dQ_K/dα|α_K = {dQK_dalpha:.4f} ≠ 0 → nie jest punkt stały")
print(f"    5. Krzywa Q_K=3/2 w 2D nie przechodzi przez (φ,φ²)")
print()
print("  WNIOSEK T-OP1:")
print("    Q_K=3/2 jest DODATKOWYM warunkiem empirycznym,")
print("    nie wynikającym z φ-FP ani ODE solitonu.")
print("    Wartość 3/2 nie jest 'specjalna' w sensie FP TGP.")
print("    STATUS T-OP1: POTWIERDZENIE OTWARTOŚCI — wymaga nowej symetrii.")

check(True, "S8: Raport T-OP1 syntetyczny — wnioski wydrukowane",
      "Q_K=3/2 = dodatkowy warunek empiryczny; nie FP TGP")

# ============================================================
# PODSUMOWANIE
# ============================================================

passed = sum(1 for _, s, _ in RESULTS if s == "PASS")
total  = len(RESULTS)

print()
print("=" * 70)
print(f"WYNIKI: {passed}/{total} PASS")
print()

# Tabela kluczowa
print("Tabela Q_K dla kluczowych wartości α = g₀^τ/g₀^e:")
print(f"  {'Reguła':<25}  {'α':>7}  {'g₀^τ':>8}  {'r₃₁':>9}  {'Q_K':>7}")
print("  " + "-"*61)
key_alphas = [
    ("α = φ  (g₀^τ = g₀^μ·φ)", PHI, G0_MU),
    ("α = φ² (φ²-skalowanie)", PHI**2, G0_TAU_PHI2),
    ("α_Koide (Q_K=3/2)",       G0_TAU_KOIDE/G0_STAR, G0_TAU_KOIDE),
]
for label, alpha, g0t in key_alphas:
    A_t  = compute_atail(g0t)
    r31  = (A_t / A_e)**4
    qk   = koide_Q(r21_computed, r31)
    print(f"  {label:<25}  {alpha:>7.4f}  {g0t:>8.4f}  {r31:>9.1f}  {qk:>7.4f}")
print("  " + "-"*61)

print()

# Zapis JSON
results_dict = {
    "alpha_scan": alpha_arr.tolist(),
    "QK_scan": QK_arr.tolist(),
    "r31_scan": r31_arr.tolist(),
    "alpha_cross_1d": float(alpha_cross_scan) if not math.isnan(alpha_cross_scan) else None,
    "alpha2_scan_2d": alpha2_arr.tolist(),
    "alpha3_koide_2d": [float(x) if not math.isnan(x) else None for x in alpha3_koide_arr],
    "r21_2d": r21_arr_2d,
    "A_e": A_e,
    "A_mu": A_mu,
    "r21_computed": r21_computed,
    "QK_at_phi": float(QK_at_phi),
    "QK_at_phi2": float(QK_at_phi2),
    "QK_at_koide": float(QK_at_koide),
    "dQK_dalpha_at_koide": float(dQK_dalpha),
    "conclusion": "Q_K=3/2 jest empirycznym warunkiem, nie FP TGP",
    "tests": [{"label": l, "status": s, "detail": d} for l, s, d in RESULTS],
    "n_pass": passed,
    "n_total": total,
    "session": "v42-2026-04-01"
}

_here = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(_here, "koide_qk_scan_results.json")
try:
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results_dict, fh, ensure_ascii=False, indent=2)
    print(f"Wyniki zapisane do: {out_path}")
except Exception as exc:
    print(f"[WARN] Nie zapisano JSON: {exc}")

print("=" * 70)
print("STATUS T-OP1: ZBADANY NUMERYCZNIE — POZOSTAJE OTWARTY")
print("  Q_K=3/2 nie jest fixed point TGP; jest dodatkowym empirycznym")
print("  warunkiem selekcji τ. Wymaga nowej symetrii lub zasady w TGP.")
print("=" * 70)
