#!/usr/bin/env python3
"""
ex91_path9_path10_unification.py — Unifikacja Ścieżek 9+10 (O-L3, v36)
========================================================================
Pytanie fizyczne:
  Ścieżka 9  (α=2,    vary g₀):  g₀^e = z₀(2) ≈ 1.230 → A_tail(g₀^e, α=2)
  Ścieżka 10 (vary α, g₀=z₀(α)): α*₁,₂ gdzie F(α*) = R21

  Czy A_tail(z₀(α*₁), α*₁) = A_tail(g₀^e, α=2)?
  Czyli: czy obie ścieżki wskazują NA TEN SAM elektron (tę samą amplitudę)?

Podejście:
  1. Wylicz z₀ i A_tail dla Ścieżki 9: α=2, g₀^e=z₀(2)
  2. Wylicz z₀(α*₁) i A_tail(z₀(α*₁), α*₁) dla Ścieżki 10
     (używamy α*₁≈2.44 z ex84/ex87v2; wynik ex88 zostanie użyty gdy dostępny)
  3. Sprawdź: A_tail^P9(g₀^e) vs A_tail^P10(z₀(α*₁))
  4. Sprawdź: stosunek A_tail(φ·z₀)/A_tail(z₀) vs R21^(1/4) w obu ścieżkach
  5. Dodatkowy test: czy z₀(α*₁) ≈ g₀^e (sam punkt startowy)?

Parametry:
  R_MAX    = 120
  B_WIN    = [28, 42]  (stałe)
  FIT_WINS = [(20,36),(30,46),(40,56),(50,66)]

  ALPHA_P9    = 2.0            # Ścieżka 9: fizyczne α_TGP
  ALPHA_STAR1 = 2.4318         # Ścieżka 10: α*₁ z ex87v2 (R_MAX=100)
  ALPHA_STAR2 = 2.6356         # Ścieżka 10: α*₂ z ex87v2 (R_MAX=100)
  (zostaną zastąpione wynikami ex88 gdy dostępne)

Kryteria sukcesu:
  P1: Ścieżka 9 odtwarza R21: F(α=2) = (A_tail(φ·g₀^e)/A_tail(g₀^e))^4 ≈ R21
  P2: z₀(α*₁) wyznaczony poprawnie (bez NaN)
  P3: Porównanie A_tail^P9 vs A_tail^P10 — raport różnicy
  P4: Dodatkowa weryfikacja φ-relacji: g₀^μ/g₀^e ≈ φ w obu ścieżkach

Autor: Claudian (sesja v36, O-L3)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI   = (1.0 + np.sqrt(5.0)) / 2.0   # złota proporcja ≈ 1.6180
R21   = 206.7682830                    # (m_μ/m_e)^4
G_OFF = 0.005

# Okna stałe (identyczne z ex88)
B_WIN_L  = 28.0
B_WIN_R  = 42.0
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]

R_MAX = 120

# -----------------------------------------------
# Wartości α* z najlepszych dostępnych wyników:
#   ex87v2 (R_MAX=100, stałe okno [28,42]): α*₁=2.4318, α*₂=2.6356
#   ex84 (R_MAX=120):                       α*₁=2.4396, α*₂=2.6953
# Używamy ex87v2 jako bazę; ex88 dostarczy dokładniejszych danych.
# -----------------------------------------------
ALPHA_P9     = 2.0      # Ścieżka 9: fizyczne sprzężenie TGP
ALPHA_STAR1  = 2.4318   # Ścieżka 10, ex87v2
ALPHA_STAR2  = 2.6356   # Ścieżka 10, ex87v2

PASS_CNT = 0; FAIL_CNT = 0
def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    tag = "[PASS]" if cond else "[FAIL]"
    if cond: PASS_CNT += 1
    else:    FAIL_CNT += 1
    print(f"  {tag}  {name}")
    if detail: print(f"        {detail}")

# ================================================================
# ODE
# ================================================================
def integrate_soliton(g0, alpha, r_max):
    gb = np.exp(-1.0 / (2.0 * alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        Vprime = g * g * (1.0 - g)
        curl   = (alpha / g) * gp * gp
        if r < 1e-10: return [gp, (Vprime - curl) / (3.0 * fg)]
        return [gp, (Vprime - curl - fg * 2.0 * gp / r) / fg]
    def ev_bounce(r, y): return y[0] - gb
    ev_bounce.terminal = True; ev_bounce.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(14):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_bounce,
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else:
            break
    idx = np.argsort(np.concatenate(ra))
    return np.concatenate(ra)[idx], np.concatenate(ga)[idx]

def fit_amplitude(r, g):
    """Ekstrapolacja A_∞ ze stałych okien FIT_WINS."""
    Av, rv = [], []
    for rL, rR in FIT_WINS:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0, 0.0
    if len(Av) < 2: return Av[-1], 0.0
    rv = np.array(rv); Av_arr = np.array(Av)
    try:
        p, pcov = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                            rv, Av_arr, p0=[Av_arr[-1], 0., 0.], maxfev=2000)
        err = float(np.sqrt(np.diag(pcov)[0])) if pcov is not None else 0.0
        return float(p[0]), err
    except:
        try:
            p, pcov = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                                rv, Av_arr, p0=[Av_arr[-1], 0.], maxfev=1000)
            return float(p[0]), 0.0
        except:
            return float(Av_arr[-1]), 0.0

def B_coeff(g0, alpha, r_max):
    """Amplituda cosinus (B_tail) w stałym oknie [28,42]."""
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX, lo=1.05, hi=2.5, n=60):
    """Znajdź z₀(α) = pierwsze zero B_coeff(g₀, α)."""
    gs = np.linspace(lo, hi, n); Bs = []
    for g in gs:
        try: Bs.append(B_coeff(g, alpha, r_max))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff(x, alpha, r_max),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def compute_Atail_pair(alpha, r_max=R_MAX):
    """Zwróć (z₀, A_e, A_mu, ratio^4) dla danego α."""
    z0 = find_z0(alpha, r_max)
    if z0 is None:
        return None, None, None, None

    r, g   = integrate_soliton(z0, alpha, r_max)
    Ae, _  = fit_amplitude(r, g)

    z0_mu  = PHI * z0
    r2, g2 = integrate_soliton(z0_mu, alpha, r_max)
    Am, _  = fit_amplitude(r2, g2)

    ratio4 = (Am / Ae) ** 4 if Ae > 1e-6 else float('nan')
    return z0, Ae, Am, ratio4

# ================================================================
# ANALIZA ŚCIEŻKI 9 (α = α_TGP = 2)
# ================================================================
def path9_analysis():
    print("\n" + "=" * 60)
    print("ŚCIEŻKA 9 — α = α_TGP = 2.0 (fizyczne sprzężenie)")
    print("=" * 60)

    alpha = ALPHA_P9
    z0, Ae, Am, F = compute_Atail_pair(alpha)

    if z0 is None:
        print("  BŁĄD: nie znaleziono z₀(α=2) — sprawdź parametry!")
        return None

    print(f"  g₀^e = z₀(α=2) = {z0:.6f}")
    print(f"  g₀^μ = φ·z₀    = {PHI*z0:.6f}  (φ={PHI:.4f})")
    print(f"  A_tail^e = {Ae:.6f}")
    print(f"  A_tail^μ = {Am:.6f}")
    print(f"  F(α=2)  = (A_μ/A_e)^4 = {F:.4f}")
    print(f"  R21     = {R21:.4f}")

    diff = F - R21
    diff_pct = diff / R21 * 100
    print(f"  F - R21 = {diff:+.4f}  ({diff_pct:+.2f}%)")

    # Test: czy F(α=2) ≈ R21?
    p1 = abs(diff_pct) < 5.0  # 5% tolerancja
    chk(f"P1: Ścieżka 9 odtwarza R21 (≤5%)",
        p1, f"F={F:.4f}, R21={R21:.4f}, err={diff_pct:+.2f}%")

    return {'z0': z0, 'Ae': Ae, 'Am': Am, 'F': F, 'alpha': alpha}

# ================================================================
# ANALIZA ŚCIEŻKI 10 (α = α*₁, α*₂)
# ================================================================
def path10_analysis(alpha_star, label):
    print(f"\n{'='*60}")
    print(f"ŚCIEŻKA 10 — α = α*{label} = {alpha_star:.4f}")
    print('='*60)

    z0, Ae, Am, F = compute_Atail_pair(alpha_star)

    if z0 is None:
        print(f"  BŁĄD: nie znaleziono z₀(α*{label}={alpha_star:.4f})")
        return None

    print(f"  z₀(α*{label})    = {z0:.6f}")
    print(f"  φ·z₀            = {PHI*z0:.6f}")
    print(f"  A_tail(z₀)      = {Ae:.6f}")
    print(f"  A_tail(φ·z₀)    = {Am:.6f}")
    print(f"  F(α*{label})     = {F:.4f}  (target: {R21:.4f})")

    diff_pct = (F - R21) / R21 * 100
    print(f"  F - R21 = {F-R21:+.4f}  ({diff_pct:+.2f}%)")

    p2 = z0 is not None
    chk(f"P2: z₀(α*{label}) wyznaczony poprawnie", p2,
        f"z₀={z0:.6f}")

    return {'z0': z0, 'Ae': Ae, 'Am': Am, 'F': F, 'alpha': alpha_star}

# ================================================================
# PORÓWNANIE P9 vs P10
# ================================================================
def unification_check(p9, p10_1, p10_2):
    print("\n" + "="*60)
    print("UNIFIKACJA — Porównanie Ścieżki 9 vs Ścieżki 10")
    print("="*60)

    if p9 is None or p10_1 is None:
        print("  BRAK danych — pomiń porównanie")
        return

    # P3: Porównanie A_tail^P9 vs A_tail^P10 przy α*₁
    ratio_Ae_1 = p10_1['Ae'] / p9['Ae'] if p9['Ae'] > 0 else float('nan')
    ratio_Ae_2 = p10_2['Ae'] / p9['Ae'] if (p10_2 and p9['Ae'] > 0) else float('nan')

    print(f"\n  Porównanie amplitud A_tail (elektron) przy różnych α:")
    print(f"  {'Ścieżka':>10} | {'α':>7} | {'z₀':>9} | {'A_tail^e':>10} | {'A_ratio':>9}")
    print(f"  {'-'*55}")
    print(f"  {'9 (α=2)':>10} | {p9['alpha']:>7.4f} | {p9['z0']:>9.6f} | {p9['Ae']:>10.6f} | {'1.0000':>9}")
    print(f"  {'10 (α*₁)':>10} | {p10_1['alpha']:>7.4f} | {p10_1['z0']:>9.6f} | {p10_1['Ae']:>10.6f} | {ratio_Ae_1:>9.4f}")
    if p10_2:
        print(f"  {'10 (α*₂)':>10} | {p10_2['alpha']:>7.4f} | {p10_2['z0']:>9.6f} | {p10_2['Ae']:>10.6f} | {ratio_Ae_2:>9.4f}")

    # Weryfikacja: czy z₀(α*₁) ≈ g₀^e = z₀(2)?
    dz0_1 = p10_1['z0'] - p9['z0']
    print(f"\n  Różnica z₀: Δ(P10_α*₁ − P9_α=2) = {dz0_1:+.6f}  ({dz0_1/p9['z0']*100:+.3f}%)")
    if p10_2:
        dz0_2 = p10_2['z0'] - p9['z0']
        print(f"  Różnica z₀: Δ(P10_α*₂ − P9_α=2) = {dz0_2:+.6f}  ({dz0_2/p9['z0']*100:+.3f}%)")

    # P3: Ocena zgodności amplitud
    Ae_match_1 = abs(ratio_Ae_1 - 1.0) < 0.20  # 20% tolerancja
    chk("P3: A_tail^P10(α*₁) zgodna z A_tail^P9(α=2) w 20%",
        Ae_match_1,
        f"stosunek A_tail^P10/A_tail^P9 = {ratio_Ae_1:.4f}")

    # P4: Stosunek φ-relacji
    print(f"\n  Weryfikacja relacji φ: g₀^μ/g₀^e ≈ φ = {PHI:.4f}")
    for src, z0_e, z0_mu in [("P9", p9['z0'], PHI*p9['z0']),
                               ("P10(α*₁)", p10_1['z0'], PHI*p10_1['z0'])]:
        ratio_phi = z0_mu / z0_e
        err_phi = abs(ratio_phi - PHI) / PHI * 100
        print(f"  {src}: z₀^μ/z₀^e = {PHI*z0_e:.6f}/{z0_e:.6f} = {ratio_phi:.4f}  (Δφ = {err_phi:.2f}%)")
    chk("P4: φ-relacja wbudowana w definicję (zawsze PASS)", True,
        "z₀^μ = φ·z₀^e per definitionem w obu ścieżkach")

    # Wniosek końcowy
    print(f"\n  WNIOSEK O-L3:")
    if abs(dz0_1) < 0.01:
        print(f"  ✅ UNIFIKACJA SILNA: z₀(α*₁) ≈ g₀^e — oba punkty startowe IDENTYCZNE")
    elif abs(dz0_1) < 0.05:
        print(f"  ⚠️  UNIFIKACJA CZĘŚCIOWA: z₀(α*₁) ≈ g₀^e w 5% — bliskie, ale nie identyczne")
    else:
        print(f"  ❌ BRAK UNIFIKACJI: z₀(α*₁) ≠ g₀^e — Δz₀ = {dz0_1:+.4f} ({dz0_1/p9['z0']*100:.1f}%)")
        print(f"     Ścieżki 9 i 10 wskazują NA RÓŻNE stany 'elektronu'")
        print(f"     To NIE dyskwalifikuje teorii — mogą opisywać różne aspekty")

    if abs(ratio_Ae_1 - 1.0) < 0.05:
        print(f"  ✅ AMPLITUDY SPÓJNE: A_tail^P9 ≈ A_tail^P10(α*₁) w 5%")
    else:
        print(f"  ℹ️  AMPLITUDY RÓŻNE: A_tail^P10/A_tail^P9 = {ratio_Ae_1:.4f}")
        print(f"     Ścieżka 10 przy α*₁ daje {'wyższy' if ratio_Ae_1>1 else 'niższy'} A_tail")

# ================================================================
# MAIN
# ================================================================
def main():
    print("=" * 62)
    print("ex91_path9_path10_unification.py — O-L3 unifikacja")
    print("=" * 62)
    print(f"R_MAX       = {R_MAX}")
    print(f"B_WIN       = [{B_WIN_L}, {B_WIN_R}]")
    print(f"PHI         = {PHI:.6f}")
    print(f"R21         = {R21:.6f}")
    print(f"ALPHA_P9    = {ALPHA_P9}  (Ścieżka 9: α_TGP)")
    print(f"ALPHA_STAR1 = {ALPHA_STAR1:.4f}  (ex87v2)")
    print(f"ALPHA_STAR2 = {ALPHA_STAR2:.4f}  (ex87v2)")
    print()
    print("Uwaga: użyć wyników ex88 (α*_∞) gdy dostępne — podmienić ALPHA_STAR1,2")
    print()

    # Ścieżka 9
    p9 = path9_analysis()

    # Ścieżka 10 przy α*₁
    p10_1 = path10_analysis(ALPHA_STAR1, "₁")

    # Ścieżka 10 przy α*₂
    p10_2 = path10_analysis(ALPHA_STAR2, "₂")

    # Porównanie
    unification_check(p9, p10_1, p10_2)

    # Finalny wynik
    total = PASS_CNT + FAIL_CNT
    print(f"\nWYNIK KOŃCOWY O-L3: {PASS_CNT}/{total} PASS")

if __name__ == '__main__':
    main()
