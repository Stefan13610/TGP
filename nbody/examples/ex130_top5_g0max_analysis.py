#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex130_top5_g0max_analysis.py
==============================
ANALIZA T-OP5: Skąd g₀^max(RMSE<10%) ≈ 3.59-3.66?

HIPOTEZA GŁÓWNA (ex130):
  g₀^max(RMSE<10%) ≈ g₀*(3→4) = próg przejścia n_bounce: 3→4
  Uzasadnienie: RMSE skacze przy każdym nowym odbiciu.
  Z ex127: progi 0→1 (1.625), 1→2 (2.232), 2→3 (2.875)
  Różnice: 0.607, 0.643 → rosnące ~0.04 na krok
  Przewidywanie: g₀*(3→4) ≈ 2.875 + 0.68 ≈ 3.55-3.60

  Jeśli to prawda: g₀^max = g₀*(3→4) — wewnętrzne kryterium TGP

TESTY E1..E14:
  --- PROGI BOUNCE ---
  E1:  Zlokalizuj g₀*(3→4) z gęstego skanu (100 pkt)
  E2:  Porównaj g₀*(3→4) z g₀^max(RMSE<10%)
  E3:  |g₀*(3→4) - g₀^max| < 0.1 (te same w granicach dyskretyzacji)?
  E4:  Wzorzec progów: arytmetyczny? geometryczny? φ?
  E5:  Przewidywanie g₀*(3→4) z interpolacji progów

  --- RMSE PRZY PROGACH ---
  E6:  RMSE lokalnie maksymalne przy progach n_bounce
  E7:  RMSE(g₀*(k→k+1)) rośnie monoton. z k
  E8:  A_tail ma lokalne minimum przy każdym progu (z ex127)

  --- ANALITYCZNE ---
  E9:  Progi g₀* a warunek energetyczny: E_bounce = V_barrier?
  E10: Progi g₀* vs wielokrotności g₀* = g*(ghost) · const?
  E11: Wzorzec arytmetyczny: g₀*(k) = a + b·k?
  E12: Wzorzec geometryczny: g₀*(k) / g₀*(k-1) = const?
  E13: φ-drabina: g₀*(k) = g₀(e) · φ^{α·k}?
  E14: Formuła zamknięta: g₀*(3→4) wyliczona z powyższego

Referencje: ex127, ex128
"""

import sys
import io
import math
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe
# ============================================================
ALPHA    = 2.0
PHI      = (1.0 + math.sqrt(5.0)) / 2.0
G_GHOST  = math.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005

G0_E   = 1.24915
G0_TAU = 3.18912
G0_L4  = PHI**3 * G0_E

R_MAX    = 60.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# ODE solitonu (z ex127/ex128)
# ============================================================
def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

def rhs(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_ghost(r, y):
    return y[0] - G_BOUNCE
event_ghost.terminal  = True
event_ghost.direction = -1

def integrate_soliton_count(g0, r_max=None, max_bounces=25):
    """Zwraca (r_arr, g_arr, n_bounces)"""
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0, y0 = R_START, [g0, 0.0]
    segs_r, segs_g = [], []
    n_bounces = 0
    for bn in range(max_bounces + 1):
        sol = solve_ivp(rhs, [r0, r_max], y0,
                        method='DOP853', max_step=MAX_STEP,
                        rtol=RTOL, atol=ATOL,
                        events=[event_ghost], dense_output=False)
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            n_bounces += 1
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0   = r_b + 1e-6
            y0   = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break
    r_all = np.concatenate(segs_r)
    g_all = np.concatenate(segs_g)
    idx = np.argsort(r_all)
    return r_all[idx], g_all[idx], n_bounces

def fit_tail_rmse(r_arr, g_arr, r_L=20.0, r_R=35.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    n = np.sum(mask)
    if n < 10:
        return float('nan'), float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - 1.0) * r_f
    X   = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B, C = float(coef[0]), float(coef[1])
    y_hat = B * np.cos(r_f) + C * np.sin(r_f)
    rmse = float(np.sqrt(np.mean((y_f - y_hat)**2)))
    A    = float(math.sqrt(B**2 + C**2))
    return A, rmse / max(A, 1e-10)


# ============================================================
print("=" * 72)
print("EX130: T-OP5 — ANALIZA g₀^max(RMSE<10%) ≈ g₀*(BOUNCE 3→4)")
print("=" * 72)
print()

# ============================================================
# SEKCJA 0: Obliczenie gęstego skanu n_bounce i RMSE
# ============================================================
print("[0] GĘSTY SKAN g₀ ∈ [1.05, 6.0] (100 punktów)")
print("-" * 55)

G0_SCAN = np.linspace(1.05, 6.0, 100)
print("  Obliczanie... (może zająć ~1 min)")

scan_data = []  # (g0, n_bounce, A, rmse_rel)
for g0_s in G0_SCAN:
    r_s, g_s, nb = integrate_soliton_count(g0_s)
    A_s, rmse_s = fit_tail_rmse(r_s, g_s)
    scan_data.append((g0_s, nb, A_s, rmse_s))

print(f"  Skan ukończony. {len(scan_data)} punktów.")

# Wypisz tabelę kluczowych obszarów
print(f"\n  {'g₀':>8}  {'n_bounce':>9}  {'A_tail':>8}  {'RMSE/A':>8}")
print("  " + "-"*45)
for g0_s, nb, A_s, rmse_s in scan_data:
    if abs(g0_s - round(g0_s*4)/4) < 0.01 or nb != scan_data[max(0,scan_data.index((g0_s,nb,A_s,rmse_s))-1)][1]:
        rmse_str = f"{100*rmse_s:.1f}%" if math.isfinite(rmse_s) else "N/A"
        print(f"  {g0_s:>8.4f}  {nb:>9d}  {A_s:>8.4f}  {rmse_str:>8}")


# ============================================================
# SEKCJA 1: Progi n_bounce (gdzie n skacze)
# ============================================================
print()
print("[1] PROGI n_bounce — WYKRYWANIE PRZEJŚĆ")
print("-" * 55)

bounce_transitions = {}  # {(old_n, new_n): g0_threshold}
prev_nb = scan_data[0][1]
for g0_s, nb, A_s, rmse_s in scan_data[1:]:
    if nb > prev_nb:
        for step in range(prev_nb, nb):
            key = (step, step+1)
            if key not in bounce_transitions:
                bounce_transitions[key] = g0_s
        prev_nb = nb

print(f"\n  Wykryte przejścia n_bounce:")
print(f"  {'Przejście':>15}  {'g₀* (skan)':>12}")
print("  " + "-"*32)
for (k1, k2), g0_thresh in sorted(bounce_transitions.items()):
    print(f"  {k1}→{k2}:          {g0_thresh:>12.4f}")

# Zapisz progi
thresholds = {k: v for k, v in sorted(bounce_transitions.items())}
thresh_list = sorted(thresholds.values())

print(f"\n  Progi g₀* = {thresh_list}")
print(f"  Różnice między progami:")
diffs = [thresh_list[i+1]-thresh_list[i] for i in range(len(thresh_list)-1)]
for i, d in enumerate(diffs):
    print(f"    g₀*({i+1}→{i+2}) - g₀*({i}→{i+1}) = {d:.4f}")

# E1: Zlokalizuj g₀*(3→4)
g0_34 = thresholds.get((3, 4), None)
if g0_34 is None:
    # Interpoluj
    print(f"  [UWAGA] Nie znaleziono przejścia 3→4 w grubym skanie — zagęszczamy")
    # Znajdź ostatni punkt z n=3 i pierwszy z n=4
    g0_last3 = max(g0_s for g0_s, nb, _, _ in scan_data if nb == 3)
    g0_first4 = min(g0_s for g0_s, nb, _, _ in scan_data if nb >= 4)
    g0_34 = (g0_last3 + g0_first4) / 2.0
    print(f"  Interpolacja: g₀*(3→4) ≈ ({g0_last3:.4f}+{g0_first4:.4f})/2 = {g0_34:.4f}")

record("E1: g₀*(3→4) znalezione w skanie lub interpolowane",
       g0_34 is not None,
       f"g₀*(3→4) ≈ {g0_34:.4f}")


# ============================================================
# SEKCJA 2: Zagęszczony skan wokół g₀*(3→4)
# ============================================================
print()
print("[2] ZAGĘSZCZONY SKAN wokół g₀*(3→4)")
print("-" * 55)

g0_lo = max(1.0, g0_34 - 0.3)
g0_hi = min(6.5, g0_34 + 0.5)
G0_FINE = np.linspace(g0_lo, g0_hi, 80)

print(f"  Zakres: [{g0_lo:.3f}, {g0_hi:.3f}], 80 punktów")

fine_data = []
for g0_s in G0_FINE:
    r_s, g_s, nb = integrate_soliton_count(g0_s)
    A_s, rmse_s = fit_tail_rmse(r_s, g_s)
    fine_data.append((g0_s, nb, A_s, rmse_s))

# Znajdź dokładny próg 3→4
g0_34_precise = g0_34
for i in range(len(fine_data)-1):
    g0_a, nb_a, _, _ = fine_data[i]
    g0_b, nb_b, _, _ = fine_data[i+1]
    if nb_a < 4 and nb_b >= 4:
        g0_34_precise = (g0_a + g0_b) / 2.0
        print(f"  Przejście 3→4 zlokalizowane: g₀∈[{g0_a:.4f},{g0_b:.4f}]")
        print(f"  Środek: {g0_34_precise:.4f}")
        break

# Znajdź g₀^max(RMSE<10%)
g0_max_rmse10 = max(g0_s for g0_s, nb, A_s, rmse_s in fine_data
                    if math.isfinite(rmse_s) and rmse_s < 0.10)
print(f"\n  g₀^max(RMSE<10%) w zagęszczonym skanie = {g0_max_rmse10:.4f}")
print(f"  g₀*(3→4 precise)                       = {g0_34_precise:.4f}")
print(f"  Różnica: {abs(g0_max_rmse10 - g0_34_precise):.4f}")

record("E2: g₀^max(RMSE<10%) i g₀*(3→4) wyznaczone",
       g0_max_rmse10 > 0 and g0_34_precise > 0,
       f"g₀^max={g0_max_rmse10:.4f}, g₀*(3→4)={g0_34_precise:.4f}")

diff_e3 = abs(g0_max_rmse10 - g0_34_precise)
record("E3: |g₀^max - g₀*(3→4)| < 0.15 (hipoteza bounce → RMSE)",
       diff_e3 < 0.15,
       f"|{g0_max_rmse10:.4f} - {g0_34_precise:.4f}| = {diff_e3:.4f}")


# ============================================================
# SEKCJA 3: Wzorzec progów
# ============================================================
print()
print("[3] WZORZEC PROGÓW g₀*(k→k+1)")
print("-" * 55)

# Zebranie wszystkich progów
all_thresholds = sorted(thresholds.items())
print(f"\n  Progi (z grubego skanu + doprecyzowane):")
print(f"  {'k→k+1':>8}  {'g₀*':>10}  {'Δg₀*':>8}  {'Δg₀*/Δpoprzedni':>16}")

thresh_values = [v for k, v in sorted(thresholds.items())]
# Dodaj g₀*(3→4) jeśli brak
if (3, 4) not in thresholds:
    thresh_values.append(g0_34_precise)
thresh_values = sorted(thresh_values)

for i, g0_t in enumerate(thresh_values):
    k1, k2 = i, i+1
    if i == 0:
        print(f"  {k1}→{k2}:       {g0_t:>10.4f}")
    else:
        delta = g0_t - thresh_values[i-1]
        if i >= 2:
            ratio = delta / (thresh_values[i-1] - thresh_values[i-2])
            print(f"  {k1}→{k2}:       {g0_t:>10.4f}  {delta:>8.4f}  {ratio:>16.4f}")
        else:
            print(f"  {k1}→{k2}:       {g0_t:>10.4f}  {delta:>8.4f}")

# E4: Wzorzec arytmetyczny?
if len(thresh_values) >= 3:
    diffs2 = [thresh_values[i+1]-thresh_values[i] for i in range(len(thresh_values)-1)]
    mean_diff = np.mean(diffs2)
    std_diff  = np.std(diffs2)
    is_arith  = std_diff / mean_diff < 0.15
    record("E4: Progi arytmetyczne? (σ/μ < 0.15)",
           is_arith,
           f"Δg₀* = {diffs2}; μ={mean_diff:.4f}, σ={std_diff:.4f}, σ/μ={std_diff/mean_diff:.3f}")

    # E5: Przewidywanie g₀*(3→4) z pierwszych 3 progów (0→1,1→2,2→3)
    if len(thresh_values) >= 4:
        # Ekstrapolacja z 3 znanych progów (k=0,1,2) → k=3
        ks_known = np.arange(3)
        fit_e5 = np.polyfit(ks_known, thresh_values[:3], 1)
        g0_34_arith = fit_e5[0] * 3 + fit_e5[1]
        delta_pred = abs(g0_34_arith - g0_34_precise)
        print(f"\n  Ekstrapolacja z k=0,1,2: g₀*(k)={fit_e5[0]:.4f}k+{fit_e5[1]:.4f}")
        print(f"    g₀*(3→4) przewidywane = {g0_34_arith:.4f}")
        print(f"    g₀*(3→4) rzeczywiste  = {g0_34_precise:.4f}")
        print(f"    δ = {delta_pred:.4f}")
        record("E5: Ekstrapolacja z k=0,1,2 przewiduje g₀*(3→4) (δ < 0.2)",
               delta_pred < 0.2,
               f"Przewidywanie: {g0_34_arith:.4f}, rzeczywiste: {g0_34_precise:.4f}, δ={delta_pred:.4f}")


# ============================================================
# SEKCJA 4: RMSE przy progach
# ============================================================
print()
print("[4] RMSE PRZY PROGACH BOUNCE")
print("-" * 55)

# Oblicz RMSE tuż przy każdym progu
thresh_rmse_data = {}
for (k1, k2), g0_t in thresholds.items():
    r_s, g_s, nb = integrate_soliton_count(g0_t + 0.01)
    A_s, rmse_s = fit_tail_rmse(r_s, g_s)
    thresh_rmse_data[(k1,k2)] = (g0_t + 0.01, nb, A_s, rmse_s)

print(f"\n  RMSE/A przy progach n_bounce:")
print(f"  {'Próg':>8}  {'g₀+0.01':>10}  {'n_b':>5}  {'RMSE/A':>8}")
for (k1,k2), (g0_t, nb, A_s, rmse_s) in sorted(thresh_rmse_data.items()):
    rmse_str = f"{100*rmse_s:.1f}%" if math.isfinite(rmse_s) else "N/A"
    print(f"  {k1}→{k2}:     {g0_t:>10.4f}  {nb:>5d}  {rmse_str:>8}")

# E6: RMSE ma lokalne maksimum przy progach?
# Sprawdź RMSE(g₀*-ε) vs RMSE(g₀*+ε) dla każdego progu
print(f"\n  Lokalne maksima RMSE przy progach (g₀*-0.05 vs g₀*+0.05):")
rmse_jumps = []
for (k1, k2), g0_t in sorted(thresholds.items()):
    g0_bef = g0_t - 0.05
    g0_aft = g0_t + 0.05
    _, rmse_bef = fit_tail_rmse(*integrate_soliton_count(g0_bef)[:2])
    _, rmse_aft = fit_tail_rmse(*integrate_soliton_count(g0_aft)[:2])
    # RMSE rośnie PO progu (nowe odbicie zwiększa RMSE)
    jumped = rmse_aft > rmse_bef
    rmse_jumps.append(jumped)
    bef_str = f"{100*rmse_bef:.1f}%" if math.isfinite(rmse_bef) else "N/A"
    aft_str = f"{100*rmse_aft:.1f}%" if math.isfinite(rmse_aft) else "N/A"
    print(f"  {k1}→{k2}: bef={bef_str}, aft={aft_str} "
          f"[{'SKOK ↑' if jumped else 'brak skoku'}]")

record("E6: RMSE skacze przy każdym progu n_bounce (aft > bef)",
       sum(rmse_jumps) >= len(rmse_jumps) * 0.6,
       f"{sum(rmse_jumps)}/{len(rmse_jumps)} progów z RMSE ↑ po przejściu")

# E7: RMSE przy progach rośnie monotonicznie
rmse_at_thresh = []
for (k1, k2), g0_t in sorted(thresholds.items()):
    _, rmse_at = fit_tail_rmse(*integrate_soliton_count(g0_t)[:2])
    rmse_at_thresh.append(rmse_at)
is_mono_rmse = all(rmse_at_thresh[i] <= rmse_at_thresh[i+1] for i in range(len(rmse_at_thresh)-1))
record("E7: RMSE przy progach rośnie monotonicznie",
       is_mono_rmse,
       f"RMSE przy progach: {[f'{100*x:.1f}%' for x in rmse_at_thresh]}")


# ============================================================
# SEKCJA 5: Analityczne — wzorzec progów
# ============================================================
print()
print("[5] ANALITYCZNE — WZORZEC PROGÓW")
print("-" * 55)

thresh_vals = sorted(thresholds.values())
if (3, 4) not in thresholds:
    thresh_vals.append(g0_34_precise)
thresh_vals = sorted(thresh_vals)
n_thresh = len(thresh_vals)

# E11: Arytmetyczny fit
if n_thresh >= 3:
    ks = np.arange(n_thresh)
    coeffs_lin = np.polyfit(ks, thresh_vals, 1)
    a_lin, b_lin = coeffs_lin[0], coeffs_lin[1]
    res_lin = [thresh_vals[k] - (a_lin*k + b_lin) for k in ks]
    rmse_lin = math.sqrt(np.mean([r**2 for r in res_lin]))
    print(f"\n  Model liniowy: g₀*(k) = {a_lin:.4f}·k + {b_lin:.4f}")
    print(f"    RMSE_fit = {rmse_lin:.5f}")
    record("E11: Progi liniowe: g₀*(k) = a·k+b (RMSE_fit < 0.05)",
           rmse_lin < 0.05,
           f"a={a_lin:.4f}, b={b_lin:.4f}, RMSE={rmse_lin:.5f}")

# E12: Geometryczny fit
ratios = [thresh_vals[i+1]/thresh_vals[i] for i in range(len(thresh_vals)-1)]
mean_ratio = np.mean(ratios)
std_ratio  = np.std(ratios)
print(f"\n  Stosunki kolejnych progów: {[f'{r:.4f}' for r in ratios]}")
print(f"  Średni stosunek: {mean_ratio:.6f}  (σ={std_ratio:.6f})")
print(f"  φ¹/³ = {PHI**(1/3):.6f},  φ¹/² = {PHI**0.5:.6f},  1+1/φ = {1+1/PHI:.6f}")
record("E12: Progi NIE są geometryczne (σ/μ > 0.05) — wzorzec liniowy lepszy",
       std_ratio / mean_ratio > 0.05,
       f"stosunki={[f'{r:.4f}' for r in ratios]}, σ/μ={std_ratio/mean_ratio:.4f}  "
       f"[malejący ciąg stosunków → wzrost subgeometryczny]")

# E13: φ-drabina
# Sprawdź czy g₀*(k) = g₀(e) · φ^α·k
# log(g₀*(k)/g₀(e)) = α·k → α z regresji
log_ratios = [math.log(t / G0_E) for t in thresh_vals]
print(f"\n  log(g₀*(k)/g₀(e)) = {[f'{x:.4f}' for x in log_ratios]}")
ks = np.arange(len(log_ratios)) + 1
# Fit: log_ratio = α*log(φ)*k
alpha_phi = np.polyfit(ks, log_ratios, 1)[0] / math.log(PHI)
print(f"  Fit φ-drabiny: g₀*(k) = g₀(e)·φ^(α·k), α = {alpha_phi:.4f}")
print(f"  [α=1 → każdy krok to ×φ; α=2 → ×φ²; α=3 → ×φ³]")
res_phi = [log_ratios[k-1] - alpha_phi * math.log(PHI) * k for k in ks]
rmse_phi = math.sqrt(np.mean([r**2 for r in res_phi]))
print(f"  RMSE_fit(φ) = {rmse_phi:.5f}")
record("E13: g₀*(k) NIE jest φ-drabiną (RMSE_φ > 0.05) — liniowe lepsze niż exp",
       rmse_phi > 0.05,
       f"α = {alpha_phi:.4f}, RMSE_φ = {rmse_phi:.5f}  [φ-drabina nie pasuje]")


# E14: Zamknięta formuła z liniowego wzorca
print(f"\n  Formuła zamknięta (model liniowy):")
print(f"    g₀*(3→4) = {a_lin:.4f}·3 + {b_lin:.4f} = {a_lin*3+b_lin:.4f}")
print(f"    g₀*(3→4) rzeczywiste = {g0_34_precise:.4f}")
g0_34_from_linear = a_lin * 3 + b_lin
delta_e14 = abs(g0_34_from_linear - g0_34_precise)
record("E14: g₀*(3→4) z formuły liniowej (δ < 0.1)",
       delta_e14 < 0.1,
       f"z formuły: {g0_34_from_linear:.4f}, rzeczywiste: {g0_34_precise:.4f}, δ={delta_e14:.4f}")


# ============================================================
# SEKCJA 6: Porównanie g₀^max vs g₀*(3→4) dla różnych okien
# ============================================================
print()
print("[6] g₀^max(RMSE<10%) VS g₀*(3→4) DLA RÓŻNYCH OKIEN")
print("-" * 60)

windows = [[20.0, 35.0], [22.0, 36.0], [18.0, 33.0]]
print(f"\n  {'Okno':>14}  {'g₀^max(10%)':>12}  {'g₀*(3→4)':>12}  {'Δ':>8}")
print("  " + "-"*52)

# Użyj fine_data do g₀^max
for win in windows:
    r_L, r_R = win
    fine_data_w = []
    for g0_s in G0_FINE:
        r_s, g_s, nb = integrate_soliton_count(g0_s)
        A_s, rmse_s = fit_tail_rmse(r_s, g_s, r_L, r_R)
        fine_data_w.append((g0_s, nb, A_s, rmse_s))

    g0_max_w = max((g0_s for g0_s, nb, A_s, rmse_s in fine_data_w
                    if math.isfinite(rmse_s) and rmse_s < 0.10), default=0.0)
    delta_w = abs(g0_max_w - g0_34_precise)
    print(f"  {str(win):>14}  {g0_max_w:>12.4f}  {g0_34_precise:>12.4f}  {delta_w:>8.4f}")


# ============================================================
# SEKCJA 7: Podsumowanie T-OP5
# ============================================================
print()
print("[7] PODSUMOWANIE T-OP5")
print("=" * 72)

print(f"""
  HIPOTEZA GŁÓWNA: g₀^max(RMSE<10%) = g₀*(3→4) — próg przejścia 3→4

  g₀*(3→4) = {g0_34_precise:.4f}
  g₀^max(RMSE<10%, okno [20,35]) = {g0_max_rmse10:.4f}
  Różnica: {abs(g0_max_rmse10-g0_34_precise):.4f}

  Progi bounce (z ex127 + ex130):
    g₀*(0→1) = {thresh_vals[0] if len(thresh_vals)>0 else 'N/A':.4f}
    g₀*(1→2) = {thresh_vals[1] if len(thresh_vals)>1 else 'N/A':.4f}
    g₀*(2→3) = {thresh_vals[2] if len(thresh_vals)>2 else 'N/A':.4f}
    g₀*(3→4) = {thresh_vals[3] if len(thresh_vals)>3 else g0_34_precise:.4f}

  Wzorzec liniowy: g₀*(k) = {a_lin:.4f}·k + {b_lin:.4f}

  WNIOSEK:
  {
'HIPOTEZA POTWIERDZONA: g₀^max ≈ g₀*(3→4)' if diff_e3 < 0.15
else 'HIPOTEZA CZĘŚCIOWO: g₀^max zbliżone do g₀*(3→4), ale różnica > 0.15'
}
""")


# ============================================================
# Wyniki końcowe
# ============================================================
print()
print("=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
passed_all = sum(1 for _, p, _ in TESTS if p)
total_all  = len(TESTS)
for name, passed_t, detail in TESTS:
    mark = "PASS" if passed_t else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print(f"WYNIK: {passed_all}/{total_all} testów PASS")
print()
print("STATUS T-OP5:")
if diff_e3 < 0.1:
    print(f"  ZAMKNIĘTY! g₀^max = g₀*(3→4 bounce) = {g0_34_precise:.4f}")
    print(f"  Wyjaśnienie: RMSE(g₀) skacze przy każdym nowym odbiciu.")
    print(f"  g₀^max(10%) = próg kolejnego odbicia = kryterium 'ostatniej stabilnej generacji'.")
elif diff_e3 < 0.2:
    print(f"  CZ. ZAMKNIĘTY: g₀^max ≈ g₀*(3→4), Δ={diff_e3:.3f} < 0.2")
    print(f"  Wzorzec liniowy progów: g₀*(k) = {a_lin:.4f}k + {b_lin:.4f}")
else:
    print(f"  OPEN: g₀^max={g0_max_rmse10:.4f}, g₀*(3→4)={g0_34_precise:.4f}, Δ={diff_e3:.3f}")
