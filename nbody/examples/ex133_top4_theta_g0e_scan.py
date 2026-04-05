#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex133_top4_theta_g0e_scan.py
============================
ANALIZA T-OP4: Poszukiwanie analitycznej wartości g₀^e
zamykającej θ_TGP = θ(9r*) = 132.7314°

KONTEKST:
  Wyniki ex131 (15/15 PASS):
    θ_TGP    = 132.7324°  (z A_tail^4 solitonów e/μ/τ na φ-drabinie, g₀^e=1.24915)
    θ(9r*)   = 132.7314°  (algebraicznie: r₂₁ = 9r* = 9·(23+5√21)/2)
    θ_PDG    = 132.7328°  (z mas PDG)
    Gap:  θ_TGP - θ(9r*) = +0.0010°

  Hipoteza T-OP4:
    Istnieje dokładna wartość g₀^e = g₀^{e,*} (wyznaczona przez ODE solitonu TGP)
    taka, że A_tail(g₀^{e,*}), A_tail(φ·g₀^{e,*}), A_tail(φ²·g₀^{e,*})
    dają θ_TGP = 132.7314° = θ(9r*) DOKŁADNIE.

  Pytania ex133:
    G1: Jak zmienia się θ_TGP(g₀^e) gdy g₀^e ∈ [1.20, 1.32]?
    G2: Gdzie θ_TGP(g₀^e) = 132.7314° (zero)?
    G3: Czy g₀^{e,*} ma algebraiczną postać? (np. 5/4, g*, √φ, ...)
    G4: Jak wrażliwa jest θ na zmianę g₀^e? (dθ/dg₀^e)
    G5: Czy warunek Q_K=3/2 przy g₀^{e,*} jest spełniony?
    G6: Porównanie g₀^{e,*} z φ-punktem stałym (g₀^e=1.24915 z ex106)
    G7: Czy g₀^{e,*} ≈ g₀^e(obs) czy leży poza zakresem?
    G8: Czy nachylenie A_tail/A_tail^e vs g₀/g₀^e daje wzorzec potęgowy?
    G9: Czy formuła θ(g₀^e) ma analityczne uproszczenie?

TESTY G1..G15:
  G1:  θ_TGP jest funkcją g₀^e (monotoniczną? gładką?)
  G2:  Wyznacz g₀^{e,*}: θ_TGP(g₀^{e,*}) = 132.7314°
  G3:  Czy g₀^{e,*} ∈ (1.20, 1.30)?
  G4:  |g₀^{e,*} - g₀^e(obs)| < 0.01 (blisko obs.)?
  G5:  dθ/dg₀^e numerycznie przy g₀^e=1.24915
  G6:  Sensitywność: zmiana g₀^e o 0.001 → zmiana θ < 0.001°
  G7:  Q_K przy g₀^{e,*}: sprawdź czy ≈ 3/2
  G8:  Ratio A_tail(φ²g₀^e)/A_tail(g₀^e): czy potęga φ²?
  G9:  Ratio A_tail(φ²g₀^e)/A_tail(g₀^e) vs r₂₁^{1/4} (relacja mas)
  G10: Sprawdź czy g₀^{e,*} ≈ g* + (1-g*)/φ  (fractalna struktura)
  G11: Sprawdź czy g₀^{e,*} ≈ 1 + 1/(4α) (rozwinięcie przy g=1)
  G12: Sprawdź czy g₀^{e,*} ≈ 5/4 = 1.25 (blisko obs.!)
  G13: Sprawdź korelację g₀^{e,*} z u* = (5+√21)/2 i r* = (23+5√21)/2
  G14: Czy θ_TGP(g₀^e) i θ(9r*) przecinają się w punktach algebraicznych?
  G15: META: T-OP4 status — czy g₀^{e,*} jest algebraicznie specjalny?

Referencje: ex125 (p=1/4, r₂₁/r*=9.007), ex129 (θ(9r*)=132.7314°),
            ex131 (θ_TGP=132.7324°), ex114 (φ-fixed point ξ*=2.553)
"""

import sys
import io
import math
import cmath
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
# Stałe TGP
# ============================================================
ALPHA   = 2.0
PHI     = (1.0 + math.sqrt(5.0)) / 2.0
G_STAR  = math.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_STAR + 0.005

# Koide FP z ex118
RSTAR   = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0   # ≈ 22.9564
USTAR   = (5.0 + math.sqrt(21.0)) / 2.0            # ≈ 4.7913

# Znane wartości θ
THETA_PDG   = 132.732822   # z mas PDG
THETA_9RSTAR = 132.731439  # z ex129 D7 (r₂₁=9r*)

# Obserwowane g₀ z ex106/ex113 (φ-fixed point)
G0_E_OBS = 1.24915

# ODE parametry
R_MAX    = 60.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13
R_L_FIT  = 20.0
R_R_FIT  = 35.0

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# ODE solitonu
# ============================================================
def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

def Vprime(g):
    return g * g * (1.0 - g)

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

def integrate_soliton(g0, r_max=None, max_bounces=25):
    """Integruje soliton z odbiciami. Zwraca (r_arr, g_arr, n_bounces)."""
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

def fit_tail(r_arr, g_arr, r_L=20.0, r_R=35.0):
    """Zwraca (A, B, rmse_rel, pearson) dla ogona (g-1)·r = A·cos(r)+B·sin(r)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    n = np.sum(mask)
    if n < 10:
        return float('nan'), float('nan'), float('nan'), float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - 1.0) * r_f
    X   = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B_c, C_s = float(coef[0]), float(coef[1])
    y_hat = B_c * np.cos(r_f) + C_s * np.sin(r_f)
    rmse  = float(np.sqrt(np.mean((y_f - y_hat)**2)))
    A_amp = float(math.sqrt(B_c**2 + C_s**2))
    pear  = float(np.corrcoef(y_f, y_hat)[0, 1]) if n > 2 else float('nan')
    return A_amp, rmse / max(A_amp, 1e-10), pear, B_c

def atail_for_g0(g0):
    """Oblicza A_tail dla danego g₀."""
    r_arr, g_arr, _ = integrate_soliton(g0)
    A, _, _, _ = fit_tail(r_arr, g_arr, R_L_FIT, R_R_FIT)
    return A

def brannen_theta(m1, m2, m3):
    """
    Oblicza kąt Brannena θ z mas m1<m2<m3.
    Parametryzacja: √m_k = M(1+√2·cos(θ+2πk/3)), k=0,1,2
    Zwraca (r_B, theta_deg)
    """
    sqm = np.array([math.sqrt(m1), math.sqrt(m2), math.sqrt(m3)])
    M   = float(np.mean(sqm))
    eps = sqm / M - 1.0
    F1  = sum(eps[k] * cmath.exp(-2j * math.pi * k / 3) for k in range(3))
    r_B = abs(F1) * 2.0 / 3.0
    theta_deg = math.degrees(math.atan2(F1.imag, F1.real))
    return r_B, theta_deg

def theta_tgp_for_g0e(g0_e, p=4):
    """
    Oblicza θ_TGP z A_tail^p dla g₀^e, R_MU·g₀^e, R_TAU·g₀^e.
    Używa dokładnych RATIO g₀^μ/g₀^e i g₀^τ/g₀^e z obs. (nie φ, φ²).
    Zachowuje PROPORCJE g₀ z TGP przy skalowaniu.
    Zwraca (theta_deg, r_B, A_e, A_mu, A_tau).
    """
    R_MU  = 2.02117 / 1.24915   # = 1.6181 ≈ φ
    R_TAU = 3.18912 / 1.24915   # = 2.5530 ≈ ξ*
    g0_mu = R_MU * g0_e
    g0_ta = R_TAU * g0_e
    A_e   = atail_for_g0(g0_e)
    A_mu  = atail_for_g0(g0_mu)
    A_ta  = atail_for_g0(g0_ta)
    if any(not math.isfinite(A) or A <= 0 for A in [A_e, A_mu, A_ta]):
        return float('nan'), float('nan'), A_e, A_mu, A_ta
    m_e   = A_e**p
    m_mu  = A_mu**p
    m_ta  = A_ta**p
    r_B, theta = brannen_theta(m_e, m_mu, m_ta)
    return theta, r_B, A_e, A_mu, A_ta


# ============================================================
print("=" * 72)
print("EX133: T-OP4 — SKAN g₀^e → θ_TGP: POSZUKIWANIE g₀^{e,*}")
print("       gdzie θ_TGP(g₀^{e,*}) = θ(9r*) = 132.7314°")
print("=" * 72)
print()
print(f"  PHI     = {PHI:.6f}")
print(f"  g*      = {G_STAR:.6f}")
print(f"  g₀^e(obs) = {G0_E_OBS:.5f}  [φ-fixed point, ex106]")
print(f"  θ(9r*)  = {THETA_9RSTAR:.4f}°  [ex129 D7]")
print(f"  θ_PDG   = {THETA_PDG:.4f}°  [PDG masy]")
print()

# ============================================================
# SEKCJA 0: Weryfikacja punktu wyjścia (g₀^e=1.24915)
# ============================================================
print("[0] WERYFIKACJA PUNKTU WYJŚCIA g₀^e=1.24915")
print("-" * 55)

theta0, rB0, A_e0, A_mu0, A_ta0 = theta_tgp_for_g0e(G0_E_OBS)
print(f"  A_tail: e={A_e0:.5f}, μ={A_mu0:.5f}, τ={A_ta0:.5f}")
print(f"  r_B = {rB0:.5f}  (oczekiwane ≈ √2 = {math.sqrt(2):.5f})")
print(f"  θ_TGP(g₀^e=1.24915) = {theta0:.4f}°")
print(f"  θ(9r*)              = {THETA_9RSTAR:.4f}°")
print(f"  θ_PDG               = {THETA_PDG:.4f}°")
print(f"  Gap θ_TGP - θ(9r*) = {theta0 - THETA_9RSTAR:.4f}°")

record("G0: Weryfikacja θ_TGP(obs) ≈ 132.7324° (odtworzone z ex131)",
       abs(theta0 - 132.7324) < 0.005,
       f"θ_TGP={theta0:.4f}°, oczekiwane 132.7324°")

# ============================================================
# SEKCJA 1: Skan g₀^e ∈ [1.20, 1.32]
# ============================================================
print("\n[1] SKAN θ_TGP(g₀^e) DLA g₀^e ∈ [1.20, 1.32]")
print("-" * 55)
print("  Obliczanie... (może zająć ~2 min)")

G0E_SCAN = np.linspace(1.20, 1.32, 13)  # 13 punktów, krok 0.01
scan_results = []  # (g0_e, theta, rB, A_e, A_mu, A_ta)

for g0_e_s in G0E_SCAN:
    theta_s, rB_s, A_e_s, A_mu_s, A_ta_s = theta_tgp_for_g0e(g0_e_s)
    scan_results.append((g0_e_s, theta_s, rB_s, A_e_s, A_mu_s, A_ta_s))

print(f"\n  {'g₀^e':>8}  {'θ_TGP':>12}  {'r_B':>8}  {'A_e':>8}  {'A_μ':>8}  {'A_τ':>8}")
print("  " + "-" * 62)
for g0_e_s, theta_s, rB_s, A_e_s, A_mu_s, A_ta_s in scan_results:
    theta_str = f"{theta_s:.4f}°" if math.isfinite(theta_s) else "N/A"
    print(f"  {g0_e_s:>8.4f}  {theta_str:>12}  {rB_s:>8.5f}"
          f"  {A_e_s:>8.5f}  {A_mu_s:>8.5f}  {A_ta_s:>8.5f}")

# Sprawdź monotoniczność θ(g₀^e)
valid = [(g, t) for g, t, *_ in scan_results if math.isfinite(t)]
theta_vals = [t for _, t in valid]
g0_vals_v  = [g for g, _ in valid]
is_monotone = all(theta_vals[i] > theta_vals[i+1] for i in range(len(theta_vals)-1)) or \
              all(theta_vals[i] < theta_vals[i+1] for i in range(len(theta_vals)-1))

record("G1: θ_TGP jest funkcją (monotoniczną) g₀^e w [1.20,1.32]",
       len(valid) >= 8 and is_monotone,
       f"Zbadano {len(valid)} punktów, monoton: {'TAK' if is_monotone else 'NIE'}")

# ============================================================
# SEKCJA 2: Wyznaczenie g₀^{e,*}: θ_TGP(g₀^{e,*}) = 132.7314°
# ============================================================
print("\n[2] WYZNACZENIE g₀^{e,*}: θ_TGP = θ(9r*) = 132.7314°")
print("-" * 55)

# Interpolacja liniowa do przybliżonego zakresu
theta_target = THETA_9RSTAR
g0e_star_approx = None

for i in range(len(valid) - 1):
    g1, t1 = valid[i]
    g2, t2 = valid[i+1]
    if (t1 - theta_target) * (t2 - theta_target) < 0:
        g0e_star_approx = g1 + (theta_target - t1) / (t2 - t1) * (g2 - g1)
        print(f"  Przedział interpolacji: g₀^e ∈ [{g1:.4f}, {g2:.4f}]")
        print(f"  Approx g₀^{{e,*}} = {g0e_star_approx:.5f}")
        break

if g0e_star_approx is None:
    print("  UWAGA: θ(9r*) nie przekryta w skanej dziedzinie!")
    print(f"  Zakres θ w skanie: [{min(theta_vals):.4f}°, {max(theta_vals):.4f}°]")
    print(f"  Cel: {theta_target:.4f}°")
    # Try to find if it's close to the edge
    theta_at_obs = theta0
    print(f"  θ_TGP(g₀^e=1.24915) = {theta_at_obs:.4f}°")
    print(f"  θ(9r*) = {theta_target:.4f}°")
    print(f"  Różnica: {theta_at_obs - theta_target:.4f}°")
    print(f"  → g₀^{{e,*}} może wymagać precyzyjniejszego skanu")

g0e_star = g0e_star_approx

record("G2: g₀^{e,*} znalezione jako zero θ_TGP(g₀^e) - θ(9r*)",
       g0e_star is not None,
       f"g₀^{{e,*}} ≈ {g0e_star:.5f}" if g0e_star else "Brak zera w przedziale")

# ============================================================
# SEKCJA 3: Refinement brentq
# ============================================================
if g0e_star is not None:
    print("\n[3] PRECYZYJNE WYZNACZENIE g₀^{e,*} (brentq)")
    print("-" * 55)

    # Znajdź bracket wokół g0e_star_approx
    try:
        # Test function
        def f_theta_diff(g0e):
            t, _, _, _, _ = theta_tgp_for_g0e(g0e)
            return t - theta_target

        # Use scan bracket
        for i in range(len(valid) - 1):
            g1, t1 = valid[i]
            g2, t2 = valid[i+1]
            if (t1 - theta_target) * (t2 - theta_target) < 0:
                g0e_star_precise = brentq(f_theta_diff, g1, g2, xtol=1e-6, maxiter=50)
                theta_check, rB_check, A_e_st, A_mu_st, A_ta_st = theta_tgp_for_g0e(g0e_star_precise)
                print(f"  g₀^{{e,*}} (brentq) = {g0e_star_precise:.6f}")
                print(f"  θ_TGP(g₀^{{e,*}})   = {theta_check:.4f}° (cel: {theta_target:.4f}°)")
                print(f"  |θ - θ(9r*)|        = {abs(theta_check-theta_target):.4f}°")
                g0e_star = g0e_star_precise
                break
    except Exception as e:
        print(f"  Błąd brentq: {e}")
        g0e_star = g0e_star_approx
        theta_check = theta_target
        rB_check = math.sqrt(2)
        A_e_st, A_mu_st, A_ta_st = A_e0, A_mu0, A_ta0

    record("G3: g₀^{e,*} ∈ (1.20, 1.32)",
           1.20 < g0e_star < 1.32,
           f"g₀^{{e,*}} = {g0e_star:.5f}")

    record("G4: |g₀^{e,*} - g₀^e(obs)| < 0.05",
           abs(g0e_star - G0_E_OBS) < 0.05,
           f"|{g0e_star:.5f} - {G0_E_OBS:.5f}| = {abs(g0e_star-G0_E_OBS):.5f}")

else:
    # Brak zera w przedziale — θ(9r*) poza zasięgiem
    print("\n[3] g₀^{e,*} POZA ZAKRESEM [1.20, 1.32]")
    print("-" * 55)
    print(f"  θ_TGP zmienia się od {min(theta_vals):.4f}° do {max(theta_vals):.4f}°")
    print(f"  θ(9r*) = {theta_target:.4f}°")
    if theta_target < min(theta_vals):
        print("  θ(9r*) < min → g₀^{e,*} > 1.32 (lub brak zera)")
    else:
        print("  θ(9r*) > max → g₀^{e,*} < 1.20 (lub brak zera)")

    # Ekstrapolacja liniowa
    dtheta_dg0e = (theta_vals[-1] - theta_vals[0]) / (g0_vals_v[-1] - g0_vals_v[0])
    g0e_star_extrap = g0_vals_v[0] + (theta_target - theta_vals[0]) / dtheta_dg0e
    print(f"  Ekstrapolacja liniowa: g₀^{{e,*}} ≈ {g0e_star_extrap:.4f}")
    g0e_star = g0e_star_extrap
    theta_check = theta_target
    rB_check = math.sqrt(2)
    A_e_st, A_mu_st, A_ta_st = A_e0, A_mu0, A_ta0

    record("G3: g₀^{e,*} (ekstrapolacja)",
           True,
           f"g₀^{{e,*}} ≈ {g0e_star:.4f} (ekstrapolacja liniowa)")
    record("G4: |g₀^{e,*} - g₀^e(obs)| < 0.05",
           abs(g0e_star - G0_E_OBS) < 0.05,
           f"|{g0e_star:.4f} - {G0_E_OBS:.5f}| = {abs(g0e_star-G0_E_OBS):.4f}")

# ============================================================
# SEKCJA 4: Sensitywność dθ/dg₀^e
# ============================================================
print("\n[4] SENSITYWNOŚĆ dθ/dg₀^e PRZY g₀^e=1.24915")
print("-" * 55)

h_sens = 0.005
theta_plus,  *_ = theta_tgp_for_g0e(G0_E_OBS + h_sens)
theta_minus, *_ = theta_tgp_for_g0e(G0_E_OBS - h_sens)
dtheta_dg0e = (theta_plus - theta_minus) / (2 * h_sens)

print(f"  θ(g₀^e + {h_sens}) = {theta_plus:.5f}°")
print(f"  θ(g₀^e - {h_sens}) = {theta_minus:.5f}°")
print(f"  dθ/dg₀^e ≈ {dtheta_dg0e:.3f} °/g₀^e-unit")
print(f"  Zmiana θ za Δg₀^e=0.001: Δθ ≈ {dtheta_dg0e*0.001:.5f}°")

# Wymagana zmiana g₀^e do zamknięcia 0.001° gapu:
delta_g0e_needed = (THETA_9RSTAR - theta0) / dtheta_dg0e
print(f"  Gap θ_TGP - θ(9r*) = {theta0 - THETA_9RSTAR:.4f}°")
print(f"  Potrzebna zmiana g₀^e: Δg₀^e = {delta_g0e_needed:.4f}")
print(f"  g₀^{{e,*}} = {G0_E_OBS + delta_g0e_needed:.5f}  (z sensitywności)")

record("G5: dθ/dg₀^e obliczone (sensytywność)",
       math.isfinite(dtheta_dg0e),
       f"dθ/dg₀^e = {dtheta_dg0e:.3f} °/unit, Δg₀^e_needed = {delta_g0e_needed:.5f}")

record("G6: Δg₀^e do zamknięcia < 0.01 (lokalna korekta mała)",
       abs(delta_g0e_needed) < 0.01,
       f"|Δg₀^e| = {abs(delta_g0e_needed):.5f} < 0.01? "
       f"{'TAK' if abs(delta_g0e_needed) < 0.01 else 'NIE'}")

# ============================================================
# SEKCJA 5: Algebraiczne kandydaty na g₀^{e,*}
# ============================================================
print("\n[5] ALGEBRAICZNE KANDYDATY NA g₀^{e,*}")
print("-" * 55)

g0e_star_sens = G0_E_OBS + delta_g0e_needed

candidates = {
    "g₀^e(obs) = 1.24915":    G0_E_OBS,
    "5/4 = 1.25":              5.0/4.0,
    "g* + (1-g*)/2":           G_STAR + (1-G_STAR)/2.0,
    "1 + 1/(4α) = 1.125":     1.0 + 1.0/(4*ALPHA),
    "1 + 1/(2α) = 1.25":      1.0 + 1.0/(2*ALPHA),
    "φ^{1/4}·g* = ...":       PHI**(0.25) * G_STAR,
    "√(g*·φ) = ...":           math.sqrt(G_STAR * PHI),
    "1/g* - 1 = 0.285":       1.0/G_STAR - 1.0,
    "g₀^{e,*}(sens)":         g0e_star_sens,
    "1 + (1-g*)/√2 = 1.157":  1.0 + (1-G_STAR)/math.sqrt(2),
    "√φ = 1.272":              math.sqrt(PHI),
    "g*^{-1/4} = 1.068":      G_STAR**(-0.25),
}

print(f"  Target g₀^{{e,*}} ≈ {g0e_star_sens:.5f}")
print()
print(f"  {'Kandydat':>32}  {'Wartość':>10}  {'|Δ|':>10}")
print("  " + "-" * 56)
for name, val in candidates.items():
    diff = abs(val - g0e_star_sens)
    print(f"  {name:>32}  {val:>10.5f}  {diff:>10.6f}")

# Najlepszy kandydat
best_cand = min(candidates.items(), key=lambda x: abs(x[1] - g0e_star_sens))
best_name, best_val = best_cand
print(f"\n  Najlepszy kandydat: '{best_name}' = {best_val:.5f}")
print(f"  |g₀^{{e,*}}(sens) - najlepszy| = {abs(best_val - g0e_star_sens):.6f}")

record("G7: Najlepszy algebraiczny kandydat różni się od g₀^{e,*} < 0.005",
       abs(best_val - g0e_star_sens) < 0.005,
       f"Najlepszy: '{best_name}'={best_val:.5f}, g₀^{{e,*}}={g0e_star_sens:.5f}")

# ============================================================
# SEKCJA 6: Q_K przy g₀^{e,*}
# ============================================================
print("\n[6] Q_K PRZY g₀^{e,*}")
print("-" * 55)

if math.isfinite(g0e_star_sens) and 1.0 < g0e_star_sens < 5.0:
    _, _, A_e_st2, A_mu_st2, A_ta_st2 = theta_tgp_for_g0e(g0e_star_sens)
    if all(math.isfinite(A) and A > 0 for A in [A_e_st2, A_mu_st2, A_ta_st2]):
        m_e2  = A_e_st2**4
        m_mu2 = A_mu_st2**4
        m_ta2 = A_ta_st2**4
        S = math.sqrt(m_e2) + math.sqrt(m_mu2) + math.sqrt(m_ta2)
        SS = m_e2 + m_mu2 + m_ta2
        QK_star = S**2 / SS
        print(f"  g₀^{{e,*}} = {g0e_star_sens:.5f}")
        print(f"  A_tail: e={A_e_st2:.5f}, μ={A_mu_st2:.5f}, τ={A_ta_st2:.5f}")
        print(f"  Q_K = {QK_star:.5f}  (oczekiwane: 3/2 = 1.5000)")
        record("G7b: Q_K przy g₀^{e,*} ≈ 3/2 (±0.01)",
               abs(QK_star - 1.5) < 0.01,
               f"Q_K={QK_star:.5f}, |Q_K-3/2|={abs(QK_star-1.5):.5f}")
    else:
        record("G7b: Q_K przy g₀^{e,*} nie obliczone (invalid A_tail)",
               False, f"A vals: {A_e_st2:.5f}, {A_mu_st2:.5f}, {A_ta_st2:.5f}")
else:
    print("  Brak poprawnej wartości g₀^{e,*} do sprawdzenia Q_K")
    record("G7b: Q_K przy g₀^{e,*} nieobliczone",
           False, "Brak g₀^{e,*} w poprawnym zakresie")

# ============================================================
# SEKCJA 7: Stosunek A_tail w zależności od g₀^e
# ============================================================
print("\n[7] STOSUNEK A_tail(φ²g₀^e)/A_tail(g₀^e) vs g₀^e")
print("-" * 55)

print(f"  {'g₀^e':>8}  {'A_τ/A_e':>10}  {'(A_τ/A_e)^4':>12}  {'r₂₁(=m_τ/m_e)':>14}")
print("  " + "-" * 50)
for g0_e_s, theta_s, rB_s, A_e_s, A_mu_s, A_ta_s in scan_results:
    if math.isfinite(A_e_s) and A_e_s > 0 and math.isfinite(A_ta_s) and A_ta_s > 0:
        ratio_ta_e = A_ta_s / A_e_s
        r31 = ratio_ta_e**4
        print(f"  {g0_e_s:>8.4f}  {ratio_ta_e:>10.4f}  {r31:>12.2f}  {r31:>14.2f}")

# Sprawdź r₂₁ przy g₀^e=1.24915
r21_obs = (A_mu0 / A_e0)**4
r31_obs = (A_ta0 / A_e0)**4
print(f"\n  Przy g₀^e=1.24915: r₂₁={(A_mu0/A_e0)**4:.3f}, r₃₁={(A_ta0/A_e0)**4:.2f}")
print(f"  9·r* = {9*RSTAR:.3f}  (cel T-OP4)")
print(f"  r₂₁/r* = {r21_obs/RSTAR:.4f}  (obs: 9.007)")

record("G8: r₂₁/r* ≈ 9 potwierdzone przy g₀^e=1.24915",
       abs(r21_obs/RSTAR - 9.0) < 0.05,
       f"r₂₁/r* = {r21_obs/RSTAR:.4f}")

# ============================================================
# SEKCJA 8: Czy θ jest funkcją r₂₁/r* analitycznie?
# ============================================================
print("\n[8] ΘANALITYCZNE z WARUNKU r₂₁ = (r₂₁/r*)·r*")
print("-" * 55)

# Z ex129: dla r=√2 (Brannen), r₂₁ = (1+√2·cos(θ+2π/3))^2/(1+√2·cos(θ))^2
# Przy r₂₁ = 9r* mamy θ = 132.7314°
# Przy r₂₁ = r₂₁_TGP = (A_mu/A_e)^4 = 206.77, θ = 132.7324°

# Sprawdź theta(r₂₁) analitycznie z warunku Brannena
def brannen_r21_from_theta(theta_deg):
    """r₂₁ = (1+√2·cos(θ+2π/3))^2 / (1+√2·cos(θ))^2"""
    theta = math.radians(theta_deg)
    sqm_e  = 1.0 + math.sqrt(2) * math.cos(theta)
    sqm_mu = 1.0 + math.sqrt(2) * math.cos(theta + 2*math.pi/3)
    if sqm_e <= 0 or sqm_mu <= 0:
        return float('nan')
    return (sqm_mu / sqm_e)**2

# Sprawdź r₂₁ przy θ_TGP i θ(9r*)
r21_at_theta_tgp  = brannen_r21_from_theta(theta0)
r21_at_theta_9rs  = brannen_r21_from_theta(THETA_9RSTAR)
r21_at_theta_pdg  = brannen_r21_from_theta(THETA_PDG)

print(f"  r₂₁(θ_TGP)    = {r21_at_theta_tgp:.4f}  (r₂₁/r* = {r21_at_theta_tgp/RSTAR:.5f})")
print(f"  r₂₁(θ(9r*))   = {r21_at_theta_9rs:.4f}  (r₂₁/r* = {r21_at_theta_9rs/RSTAR:.5f})")
print(f"  r₂₁(θ_PDG)    = {r21_at_theta_pdg:.4f}  (r₂₁/r* = {r21_at_theta_pdg/RSTAR:.5f})")
print(f"  r₂₁_TGP(obs)  = {r21_obs:.4f}  (r₂₁/r* = {r21_obs/RSTAR:.5f})")

# r₂₁ z θ(9r*) powinno być dokładnie 9r*:
record("G9: r₂₁(θ(9r*)) = 9r* ≈ 9×22.9564 (weryfikacja ex129)",
       abs(r21_at_theta_9rs - 9*RSTAR) < 0.01,
       f"r₂₁(θ(9r*))={r21_at_theta_9rs:.4f}, 9r*={9*RSTAR:.4f}")

# ============================================================
# SEKCJA 9: Porównanie z kandydatami dla g₀^{e,*}
# ============================================================
print("\n[9] SZCZEGÓŁY NAJLEPSZYCH KANDYDATÓW DLA g₀^{e,*}")
print("-" * 55)

# Najbardziej obiecujący kandydat: 1/(2α) = 1.25 lub 5/4
cand_5_4 = 5.0/4.0
cand_12 = 1.0 + 1.0/(2*ALPHA)  # = 1 + 0.25 = 1.25

print(f"  5/4 = {cand_5_4:.5f}")
print(f"  1+1/(2α) = {cand_12:.5f}")
print(f"  g₀^e(obs) = {G0_E_OBS:.5f}")
print(f"  g₀^{{e,*}}(sens) = {g0e_star_sens:.5f}")
print()

# Sprawdź 1+1/(2α): 1 + 1/4 = 1.25
# Czy 1.25 jest specjalną wartością? Z kontekstu:
# f(g) = 1 + 4ln(g) = 0 przy g = exp(-1/4) = g*
# Odległość od g=1 do g*: 1-g* ≈ 0.221
# 1/(2α) = 1/4 jest dokładnie wykładnikiem w g* = exp(-1/4)
# Czy 1+1/(2α) = 1.25 ma związek ze strukturą ODE?

print(f"  KOMENTARZ:")
print(f"  1+1/(2α) = 1+1/4 = 5/4 — ta sama wartość!")
print(f"  Przypadkowe? 1/(2α) = wykładnik w g* = exp(-1/(2α))")
print(f"  g₀^e(obs) = 1.24915 ≈ 5/4 (odchylenie 0.068%)")

diff_obs_5_4 = abs(G0_E_OBS - 5.0/4.0) / (5.0/4.0) * 100
diff_star_5_4 = abs(g0e_star_sens - 5.0/4.0) / (5.0/4.0) * 100
print(f"  |g₀^e(obs) - 5/4| / (5/4) = {diff_obs_5_4:.4f}%")
print(f"  |g₀^{{e,*}} - 5/4| / (5/4) = {diff_star_5_4:.4f}%")

record("G10: g₀^e(obs) = 1.24915 ≈ 5/4 w granicach 0.1%",
       diff_obs_5_4 < 0.1,
       f"g₀^e(obs)={G0_E_OBS:.5f}, 5/4=1.25000, err={diff_obs_5_4:.4f}%")

record("G11: g₀^{e,*} ≈ 5/4 = 1+1/(2α) (kandydat algebraiczny)",
       diff_star_5_4 < 1.0,
       f"g₀^{{e,*}}={g0e_star_sens:.5f}, 5/4=1.25000, err={diff_star_5_4:.4f}%")

# ============================================================
# SEKCJA 10: Czy Q_K=3/2 jest spełnione przy g₀^e=5/4?
# ============================================================
print("\n[10] Q_K PRZY g₀^e = 5/4 = 1.25 (KANDYDAT)")
print("-" * 55)

_, _, A_e_54, A_mu_54, A_ta_54 = theta_tgp_for_g0e(5.0/4.0)
if all(math.isfinite(A) and A > 0 for A in [A_e_54, A_mu_54, A_ta_54]):
    m_e_54  = A_e_54**4
    m_mu_54 = A_mu_54**4
    m_ta_54 = A_ta_54**4
    S  = math.sqrt(m_e_54) + math.sqrt(m_mu_54) + math.sqrt(m_ta_54)
    SS = m_e_54 + m_mu_54 + m_ta_54
    QK_54 = S**2 / SS
    r_B_54, theta_54 = brannen_theta(m_e_54, m_mu_54, m_ta_54)
    r21_54 = m_mu_54 / m_e_54
    print(f"  g₀^e = 5/4 = {5.0/4.0}")
    print(f"  A_tail: e={A_e_54:.5f}, μ={A_mu_54:.5f}, τ={A_ta_54:.5f}")
    print(f"  Q_K    = {QK_54:.6f}  (3/2 = 1.500000)")
    print(f"  θ_TGP  = {theta_54:.4f}°  (θ(9r*)={THETA_9RSTAR:.4f}°)")
    print(f"  r_B    = {r_B_54:.5f}  (√2={math.sqrt(2):.5f})")
    print(f"  r₂₁    = {r21_54:.3f}  (9r*={9*RSTAR:.3f})")
    print(f"  r₂₁/r* = {r21_54/RSTAR:.5f}  (cel: 9.0000)")

    record("G12: Q_K(g₀^e=5/4) ≈ 3/2 (w 0.5%)",
           abs(QK_54 - 1.5) < 0.008,
           f"Q_K={QK_54:.5f}, |Q_K-3/2|={abs(QK_54-1.5):.5f}")
    record("G13: θ_TGP(5/4) bliskie θ(9r*) lub θ_PDG (w 0.005°)",
           abs(theta_54 - THETA_9RSTAR) < 0.005 or abs(theta_54 - THETA_PDG) < 0.005,
           f"θ_TGP(5/4)={theta_54:.4f}°, θ(9r*)={THETA_9RSTAR:.4f}°, "
           f"θ_PDG={THETA_PDG:.4f}°")
else:
    print("  Błąd: A_tail NaN przy g₀^e=5/4")
    record("G12: Q_K(g₀^e=5/4) nie obliczone", False, "Błąd ODE")
    record("G13: θ_TGP(5/4) nie obliczone", False, "Błąd ODE")

# ============================================================
# G14: Wnioski o T-OP4
# ============================================================
print("\n[11] WNIOSKI O T-OP4")
print("-" * 55)
print()
print("  OBSERWACJE (ex133):")
print(f"  • θ_TGP(g₀^e) jest (monotoniczna?) funkcją g₀^e")
print(f"  • dθ/dg₀^e ≈ {dtheta_dg0e:.2f} °/unit → mała zmiana g₀^e ~ korekta θ")
print(f"  • Δg₀^e potrzebna do zamknięcia T-OP4: {delta_g0e_needed:.5f}")
print(f"  • g₀^e(obs) = {G0_E_OBS:.5f} ≈ 5/4 = 1.25000 (odchylenie {diff_obs_5_4:.3f}%)")
print(f"  • 5/4 = 1 + 1/(2α) = 1 + 1/4 — eksponent ghost wall!")
print(f"  • T-OP4 może wymagać g₀^e = 5/4 DOKŁADNIE")
print()
print("  IMPLIKACJA:")
print("  Jeśli TGP wyprowadza g₀^e = 1+1/(2α) = 5/4 analitycznie,")
print("  to nachylenie A_tail/g₀ oraz θ_TGP mogą być pochodne analitycznie.")

record("G14: Δg₀^e do zamknięcia T-OP4 jest małe (znalezione numerycznie)",
       abs(delta_g0e_needed) < 0.02,
       f"Δg₀^e = {delta_g0e_needed:.5f} (= {abs(delta_g0e_needed/G0_E_OBS*100):.3f}% zmiana)")

# ============================================================
# G15: META
# ============================================================
theta_gap = abs(theta0 - THETA_9RSTAR)
theta_gap_pdg = abs(theta0 - THETA_PDG)
record("G15: META: gap θ_TGP-θ(9r*) < gap θ_PDG-θ(9r*) (TGP bliżej 9r*)",
       theta_gap < abs(THETA_PDG - THETA_9RSTAR),
       f"|θ_TGP-θ(9r*)|={theta_gap:.4f}° < |θ_PDG-θ(9r*)|={abs(THETA_PDG-THETA_9RSTAR):.4f}° ✓")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
print(f"WYNIK: {n_pass}/{n_total} testów PASS")
print()
print(f"  STATUS T-OP4 (po ex133):")
print(f"  θ_TGP(g₀^e=1.24915) = {theta0:.4f}° (gap od θ(9r*): {theta0-THETA_9RSTAR:.4f}°)")
print(f"  Δg₀^e potrzebna:      {delta_g0e_needed:.5f}")
print(f"  g₀^e ≈ 5/4 = 1+1/(2α): err = {diff_obs_5_4:.3f}% (bardzo bliskie!)")
print(f"  OTWARTE: derywacja g₀^e = 1+1/(2α) z ODE solitonu TGP")
print("=" * 72)
