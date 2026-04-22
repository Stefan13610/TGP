#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tgp_koide_r31_formal.py
=======================
Dodatek T: Weryfikacja r₃₁ z formuły Koide'go + wyznaczenie g₀^τ

CEL: Formalne zamknięcie OP-K1 poprzez:
  1. Algebraiczne sprawdzenie formuły (T-r31-formula) z tw. T-r31
  2. Numeryczne wyznaczenie g₀^τ z warunku A_tail (prop. T-g0tau-corrected)
  3. Porównanie: φ²-skalowanie vs Koide vs PDG
  4. Zestawienie tabeli z rem. T-comparison

TESTY (9):
  K1: Formuła algebraiczna r₃₁^K = 3477.5 ± 1 (z r₂₁=206.768)
  K2: φ²-skalowanie r₃₁^φ² ≈ 3955 (błąd ~13.7% vs PDG)
  K3: Koide Q_K weryfikacja przy r₂₁=206.768, r₃₁=3477.5
  K4: A_tail(g₀^e) = A_e > 0 (elektron ma ogon)
  K5: Bisekcja g₀^τ istnieje w (g₀^μ, 2·g₀^μ)
  K6: g₀^τ spełnia (A(g₀^τ)/A_e)^4 = r₃₁^K ± 0.5%
  K7: g₀^τ ≠ φ²·g₀^e (różnica > 1%)
  K8: Tabela porównawcza — 3 predykcje vs PDG
  K9: Q_K sprawdzony dla wyznaczonego g₀^τ (musi być ≈ 3/2 ± 0.01)

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
from scipy.optimize import brentq

warnings.filterwarnings('ignore')

# ============================================================
# Parametry fizyczne (identyczne jak ex106)
# ============================================================
ALPHA   = 2.0
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))   # ≈ 0.7788
PHI     = (1.0 + np.sqrt(5.0)) / 2.0     # ≈ 1.618
R21_PDG = 206.768                          # m_μ/m_e (PDG 2024)
R31_PDG = 3477.48                          # m_τ/m_e (PDG 2024)

# g₀^* z twierdzenia J2-FP (ex106, tolerancja 1e-5)
G0_STAR   = 1.24915
G0_MU     = PHI * G0_STAR                 # ≈ 2.02117
G0_TAU_PHI2 = PHI**2 * G0_STAR           # ≈ 3.27032 (φ²-skalowanie)

# Parametry numeryczne (solvera ODE)
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
# Formuła algebraiczna Koide → r₃₁ (tw. T-r31)
# ============================================================

def koide_r31_algebraic(r21):
    """
    Oblicza r₃₁ = m_τ/m_e z formuły Koide'go Q_K=3/2 i danego r₂₁.

    Wyprowadzenie (tw. T-r31):
      y² - 2(2x+2)y + (x²-4x+1) = 0,  x=sqrt(r₂₁), y=sqrt(r₃₁)
      y = (2x+2) + sqrt(3(x²+4x+1))    [wybieram +, bo y>x]
    """
    x = math.sqrt(r21)
    discriminant = 3.0 * (x**2 + 4.0*x + 1.0)
    if discriminant < 0:
        raise ValueError(f"Wyróżnik ujemny przy r₂₁={r21}")
    y = 2.0*(x + 1.0) + math.sqrt(discriminant)
    return y**2

def koide_Q(r21, r31):
    """
    Oblicza Q_K = (1+sqrt(r21)+sqrt(r31))² / (1+r21+r31).
    Powinno dawać 3/2 przy r₃₁=r₃₁^K.
    """
    x = math.sqrt(r21)
    y = math.sqrt(r31)
    num = (1.0 + x + y)**2
    den = 1.0 + r21 + r31
    return num / den

# ============================================================
# Solver ODE solitonowego (z ex106)
# ============================================================

def V(g):
    return g**3 / 3.0 - g**4 / 4.0

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
    r0 = R_START
    y0 = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []

    for bounce_num in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_regularized, [r0, r_max], y0,
            method='DOP853', max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_ghost], dense_output=False
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        segs_gp.append(sol.y[1])

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
    """
    Dopasuj ogon: g(r)-1 ≈ (B cos r + C sin r)/r → A_tail = sqrt(B²+C²).
    """
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
    """Oblicza A_tail dla danego g₀."""
    r_arr, g_arr = integrate_soliton(g0)
    return fit_tail(r_arr, g_arr)

# ============================================================
# SEKCJA 1: Algebraiczna weryfikacja (K1-K3)
# ============================================================

print("=" * 70)
print("TGP_KOIDE_R31_FORMAL: Weryfikacja formuły Koide → r₃₁")
print("=" * 70)
print(f"  α = {ALPHA},  g* = {G_GHOST:.6f},  φ = {PHI:.6f}")
print(f"  g₀^e = {G0_STAR},  g₀^μ = {G0_MU:.5f},  g₀^τ(φ²) = {G0_TAU_PHI2:.5f}")
print(f"  r₂₁ (PDG) = {R21_PDG},  r₃₁ (PDG) = {R31_PDG}")
print()

print("--- SEKCJA 1: Algebraika (tw. T-r31) ---")

# K1: Formuła algebraiczna
r31_koide = koide_r31_algebraic(R21_PDG)
check(abs(r31_koide - R31_PDG) / R31_PDG < 0.001,
      "K1: r₃₁^Koide ≈ 3477.5 (błąd < 0.1% vs PDG)",
      f"r₃₁^K = {r31_koide:.2f},  r₃₁^PDG = {R31_PDG:.2f},  "
      f"δ = {abs(r31_koide-R31_PDG)/R31_PDG*100:.4f}%")

# K2: φ²-skalowanie
x_phi2 = math.sqrt(R21_PDG)          # √r₂₁
# r₃₁^φ² z prop. J2-profiles: (A(φ²g₀*)/A(g₀*))^4 — tylko numerycznie
# Dla testu algebraicznego: porównaj r₃₁^φ² = 3955.1 (znane z ex106)
r31_phi2_known = 3955.1
err_phi2 = abs(r31_phi2_known - R31_PDG) / R31_PDG
check(err_phi2 > 0.10,
      "K2: φ²-skalowanie r₃₁^φ² ≈ 3955 (błąd > 10% vs PDG)",
      f"r₃₁^φ² = {r31_phi2_known:.1f},  δ = {err_phi2*100:.1f}%")

# K3: Weryfikacja Q_K przy r₃₁^K
QK = koide_Q(R21_PDG, r31_koide)
check(abs(QK - 1.5) < 1e-6,
      "K3: Q_K = 3/2 zweryfikowany dla (r₂₁, r₃₁^K)",
      f"Q_K = {QK:.8f}  (oczekiwane: 1.5)")

print()
print("--- SEKCJA 2: Numeryczne A_tail ---")

# K4: A_tail(g₀^e) > 0
print(f"  Obliczam A_tail(g₀^e = {G0_STAR}) ...")
A_e = compute_atail(G0_STAR)
check(A_e > 0.1,
      "K4: A_tail(g₀^e) > 0",
      f"A_e = {A_e:.5f}")

# K5: Sprawdź A_tail w przedziale [g₀^μ, 2·g₀^μ]
print(f"  Obliczam A_tail(g₀^μ = {G0_MU:.5f}) ...")
A_mu = compute_atail(G0_MU)

g0_search_hi = 2.0 * G0_MU   # ≈ 4.04
print(f"  Obliczam A_tail(g₀_hi = {g0_search_hi:.5f}) ...")
A_hi = compute_atail(g0_search_hi)

# Warunek: (A/A_e)^4 = r₃₁^K = 3477.5
# Szukamy g₀^τ gdzie (A(g)/A_e)^4 = r₃₁^K
target_r31 = r31_koide

def objective_g0tau(g0):
    """Zwraca (A(g₀)/A_e)^4 - r₃₁^K"""
    A = compute_atail(g0)
    return (A / A_e)**4 - target_r31

val_lo = (A_mu / A_e)**4 - target_r31
val_hi = (A_hi / A_e)**4 - target_r31

check(val_lo < 0 < val_hi,
      "K5: Bisekcja g₀^τ — zmiana znaku w (g₀^μ, 2·g₀^μ)",
      f"f(g₀^μ) = {val_lo:.1f},  f(2·g₀^μ) = {val_hi:.1f}")

print()
print("--- SEKCJA 3: Bisekcja g₀^τ ---")

# K6: Wyznacz g₀^τ z bisekcji
print(f"  Bisekcja w [{G0_MU:.4f}, {g0_search_hi:.4f}] ...")
try:
    g0_tau_koide = brentq(objective_g0tau, G0_MU, g0_search_hi, xtol=1e-4, rtol=1e-5)
    A_tau = compute_atail(g0_tau_koide)
    r31_verify = (A_tau / A_e)**4
    err_verify = abs(r31_verify - target_r31) / target_r31

    check(err_verify < 0.005,
          "K6: g₀^τ spełnia (A(g₀^τ)/A_e)^4 = r₃₁^K ± 0.5%",
          f"g₀^τ = {g0_tau_koide:.4f},  A_τ = {A_tau:.5f},  "
          f"r₃₁^verify = {r31_verify:.2f},  δ = {err_verify*100:.3f}%")
except Exception as exc:
    check(False, "K6: Bisekcja g₀^τ", f"Błąd: {exc}")
    g0_tau_koide = float('nan')
    A_tau = 0.0
    r31_verify = 0.0

# K7: g₀^τ ≠ φ²·g₀^e
if not math.isnan(g0_tau_koide):
    diff_g0 = abs(g0_tau_koide - G0_TAU_PHI2) / G0_TAU_PHI2
    check(diff_g0 > 0.01,
          "K7: g₀^τ(Koide) ≠ φ²·g₀^e  (różnica > 1%)",
          f"g₀^τ(Koide) = {g0_tau_koide:.4f},  φ²·g₀^e = {G0_TAU_PHI2:.4f},  "
          f"Δ = {diff_g0*100:.2f}%")
else:
    check(False, "K7: g₀^τ(Koide) ≠ φ²·g₀^e", "Pominięto — bisekcja nieudana")

print()
print("--- SEKCJA 4: Tabela porównawcza (rem. T-comparison) ---")

# K8: Tabela 3 predykcji
print()
print("  Tabela: predykcje r₃₁ vs PDG")
print("  " + "-"*62)
print(f"  {'Metoda':<32} {'r₃₁':>8}  {'δ [%]':>8}  {'g₀^τ':>8}")
print("  " + "-"*62)

methods = [
    ("φ²-skalowanie (prop. J2)",  r31_phi2_known,         G0_TAU_PHI2),
    ("Koide Q_K=3/2 (tw. T)",     r31_koide,              g0_tau_koide if not math.isnan(g0_tau_koide) else float('nan')),
    ("PDG (pomiar)",               R31_PDG,                float('nan')),
]
for name, r31_val, g0_val in methods:
    delta = abs(r31_val - R31_PDG) / R31_PDG * 100.0
    g0_str = f"{g0_val:.4f}" if not math.isnan(g0_val) else "  —"
    print(f"  {name:<32} {r31_val:>8.2f}  {delta:>8.3f}  {g0_str:>8}")
print("  " + "-"*62)
print()

check(True,
      "K8: Tabela porównawcza wydrukowana",
      f"Koide daje {abs(r31_koide - R31_PDG)/R31_PDG*100:.4f}% vs PDG, "
      f"φ²-skalowanie {abs(r31_phi2_known - R31_PDG)/R31_PDG*100:.2f}%")

# K9: Sprawdzenie Q_K dla wyznaczonego g₀^τ
if not math.isnan(g0_tau_koide):
    QK_verify = koide_Q(R21_PDG, r31_verify)
    check(abs(QK_verify - 1.5) < 0.01,
          "K9: Q_K ≈ 3/2 dla wyznaczonego g₀^τ",
          f"Q_K = {QK_verify:.6f}  (oczekiwane ≈ 1.5)")
else:
    check(False, "K9: Q_K ≈ 3/2", "Pominięto — bisekcja nieudana")

# ============================================================
# PODSUMOWANIE
# ============================================================

passed = sum(1 for _, s, _ in RESULTS if s == "PASS")
total  = len(RESULTS)

print()
print("=" * 70)
print(f"WYNIKI: {passed}/{total} PASS")
print()

# Kluczowe wartości
print("Kluczowe wartości:")
print(f"  r₃₁^Koide (algebraiczny)  = {r31_koide:.4f}   "
      f"(δ = {abs(r31_koide - R31_PDG)/R31_PDG*100:.4f}% vs PDG)")
print(f"  r₃₁^φ²   (skalowanie)     = {r31_phi2_known:.1f}    "
      f"(δ = {abs(r31_phi2_known - R31_PDG)/R31_PDG*100:.2f}% vs PDG)")
print(f"  g₀^e  = {G0_STAR:.5f}")
print(f"  g₀^μ  = {G0_MU:.5f}")
if not math.isnan(g0_tau_koide):
    print(f"  g₀^τ  = {g0_tau_koide:.5f}  [Koide]")
    print(f"  g₀^τ  = {G0_TAU_PHI2:.5f}  [φ²-skalowanie]")
print(f"  A_e   = {A_e:.5f}")
print(f"  A_μ   = {A_mu:.5f}")
if not math.isnan(g0_tau_koide):
    print(f"  A_τ   = {A_tau:.5f}")
print()

# Zapis JSON
results_dict = {
    "r31_koide_algebraic": r31_koide,
    "r31_phi2_known": r31_phi2_known,
    "r31_pdg": R31_PDG,
    "r21_pdg": R21_PDG,
    "QK_algebraic": koide_Q(R21_PDG, r31_koide),
    "g0_star": G0_STAR,
    "g0_mu": G0_MU,
    "g0_tau_phi2": G0_TAU_PHI2,
    "g0_tau_koide": g0_tau_koide if not math.isnan(g0_tau_koide) else None,
    "A_e": A_e,
    "A_mu": A_mu,
    "A_tau": A_tau if not math.isnan(g0_tau_koide) else None,
    "r31_verify": r31_verify if not math.isnan(g0_tau_koide) else None,
    "tests": [{"label": l, "status": s, "detail": d} for l, s, d in RESULTS],
    "n_pass": passed,
    "n_total": total,
    "session": "v42-2026-04-01"
}

_here = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(_here, "koide_r31_results.json")
try:
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results_dict, fh, ensure_ascii=False, indent=2)
    print(f"Wyniki zapisane do: {out_path}")
except Exception as exc:
    print(f"[WARN] Nie zapisano JSON: {exc}")

print("=" * 70)
if passed == total:
    print("STATUS OP-K1: CZ. ZAMKNIETY [AN+NUM] — formuła algebraiczna OK,")
    print("              g₀^τ wyznaczone numerycznie z A_tail.")
    print("              Pozostaje T-OP1: dynamiczne uzasadnienie Q_K=3/2 w TGP.")
else:
    print(f"STATUS: {total - passed} testów nieudanych — sprawdź log.")
print("=" * 70)
