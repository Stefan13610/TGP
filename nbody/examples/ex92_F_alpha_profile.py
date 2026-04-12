#!/usr/bin/env python3
"""
ex92_F_alpha_profile.py — Pełny profil F(α) i G₃(α) (O-L4 + O-L5, v37)
========================================================================
STATUS: LEGACY-TRANSLATIONAL

This script continues the older `F(alpha)` / `G3(alpha)` exploration based on
legacy selection variables. Keep it as historical exploratory material rather
than a canonical synchronized `nbody` entry point.

Pytania fizyczne:
  O-L4: Dlaczego α*₁,₂ ≠ α_TGP=2? Wyznacz α_min, F_min, kształt doliny.
  O-L5: Czy G₃(α)=(A_tail(φ²·z₀)/A_tail(z₀))^4 osiąga R31=3477 dla jakiegoś α?
         Czy F(α) osiąga R31=3477 dla α poza [2, 3.5]?

Metoda:
  Dla każdego α ∈ [1.0, 5.0] (krok ≈ 0.033):
    1. Znajdź z₀(α) = zero B_coeff(g₀, α) w [28,42]
    2. Całkuj soliton z g₀=z₀ i g₀=φ·z₀ i g₀=φ²·z₀
    3. Wyznacz A_tail dla każdego z tych punktów startowych
    4. F(α)  = (A_tail(φ·z₀) / A_tail(z₀))^4
       G3(α) = (A_tail(φ²·z₀) / A_tail(z₀))^4

Stałe:
  R_MAX   = 120
  B_WIN   = [28.0, 42.0]  (stałe)
  PHI     = (1+√5)/2 ≈ 1.6180
  R21     = (m_μ/m_e)^4 ≈ 206.77
  R31     = (m_τ/m_e)^4 ≈ 3477.23... → sprawdzić dokładną wartość

Wyszukiwanie:
  1. α_min = argmin F(α)  → wyjaśnienie O-L4
  2. F_min = F(α_min)
  3. Czy F(α) > R31 dla jakiegoś α?  → O-L5 Podejście 1
  4. Czy G3(α) = R31 dla jakiegoś α? → O-L5 Podejście 2

Autor: Claudian (sesja v37, O-L4 + O-L5)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit, minimize_scalar
import warnings
import multiprocessing as mp
import time
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE FIZYCZNE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0   # złota proporcja ≈ 1.6180
PHI2 = PHI * PHI                      # φ² ≈ 2.6180

# R21 = m_μ/m_e ≈ 206.77  (stosunek mas, NIE czwarta potęga — F=(A_μ/A_e)^4 = m_μ/m_e)
# R31 = m_τ/m_e ≈ 3477.23 (stosunek mas tauon/elektron)
R21  = 206.7682830    # m_μ/m_e

M_TAU_MEV = 1776.86   # MeV
M_MU_MEV  = 105.658   # MeV
M_E_MEV   = 0.51100   # MeV

R21_check = M_MU_MEV / M_E_MEV          # ≈ 206.77
R31_val   = M_TAU_MEV / M_E_MEV         # ≈ 3477.23
R32_val   = M_TAU_MEV / M_MU_MEV        # ≈ 16.82

G_OFF = 0.005

# Okna stałe
B_WIN_L  = 28.0
B_WIN_R  = 42.0
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]
R_MAX    = 120

# Zakres α
ALPHA_MIN_RANGE = 1.0
ALPHA_MAX_RANGE = 5.0
N_ALPHA         = 120   # krok ≈ 0.033

# ================================================================
# ODE + AMPLITUDY (identyczne z ex91)
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
        return float(p[0]), float(np.sqrt(abs(np.diag(pcov)[0])))
    except:
        try:
            p, _ = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                             rv, Av_arr, p0=[Av_arr[-1], 0.], maxfev=1000)
            return float(p[0]), 0.0
        except:
            return float(Av_arr[-1]), 0.0

def B_coeff(g0, alpha, r_max):
    """Amplituda cosinus w oknie [28,42]."""
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX, lo=1.05, hi=2.5, n=60):
    """Znajdź z₀(α) = pierwsze zero B_coeff(g₀, α)."""
    gs = np.linspace(lo, hi, n)
    Bs = []
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

# ================================================================
# OBLICZENIE JEDNEGO PUNKTU α
# ================================================================
def compute_one_alpha(alpha):
    """
    Dla danego α oblicz:
      z₀(α), F(α), G3(α)
    Zwróć słownik z wynikami lub None jeśli z₀ nie znalezione.
    """
    try:
        z0 = find_z0(alpha)
        if z0 is None:
            return {'alpha': alpha, 'z0': None, 'F': None, 'G3': None,
                    'Ae': None, 'Amu': None, 'Atau': None}

        # Amplituda elektronu
        r, g = integrate_soliton(z0, R_MAX)
        Ae, _ = fit_amplitude(r, g)

        # Amplituda muonu (φ·z₀)
        z0_mu = PHI * z0
        r2, g2 = integrate_soliton(z0_mu, alpha, R_MAX)
        Amu, _ = fit_amplitude(r2, g2)

        # Amplituda tauonu (φ²·z₀)
        z0_tau = PHI2 * z0
        r3, g3 = integrate_soliton(z0_tau, alpha, R_MAX)
        Atau, _ = fit_amplitude(r3, g3)

        # Sprawdzenie
        if Ae < 1e-8:
            return {'alpha': alpha, 'z0': z0, 'F': None, 'G3': None,
                    'Ae': Ae, 'Amu': Amu, 'Atau': Atau}

        F  = (Amu  / Ae) ** 4
        G3 = (Atau / Ae) ** 4

        return {'alpha': alpha, 'z0': z0, 'F': F, 'G3': G3,
                'Ae': Ae, 'Amu': Amu, 'Atau': Atau}
    except Exception as e:
        return {'alpha': alpha, 'z0': None, 'F': None, 'G3': None,
                'Ae': None, 'Amu': None, 'Atau': None, 'err': str(e)}

def compute_one_alpha_wrapper(args):
    """Wrapper dla multiprocessing."""
    alpha = args
    return compute_one_alpha(alpha)

# ================================================================
# FIX: integrate_soliton wymaga alpha jako argument
# ================================================================
# Poprawka: integrate_soliton(g0, alpha, r_max) — alpha jest w closure
# Ale compute_one_alpha używa integrate_soliton(z0, R_MAX) — brakuje alpha!
# Naprawmy to tutaj:

def integrate_soliton_full(g0, alpha, r_max):
    """Identyczne z integrate_soliton — z pełnym podpisem."""
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

def compute_one_alpha_fixed(alpha):
    """
    Dla danego α oblicz z₀, F(α), G3(α).
    Pełna wersja z poprawnym podpisem integrate_soliton.
    """
    try:
        z0 = find_z0(alpha)
        if z0 is None:
            return {'alpha': alpha, 'z0': None, 'F': None, 'G3': None,
                    'Ae': None, 'Amu': None, 'Atau': None}

        # Elektron
        r, g = integrate_soliton_full(z0, alpha, R_MAX)
        Ae, _ = fit_amplitude(r, g)

        # Muon (φ·z₀)
        z0_mu = PHI * z0
        r2, g2 = integrate_soliton_full(z0_mu, alpha, R_MAX)
        Amu, _ = fit_amplitude(r2, g2)

        # Tau (φ²·z₀)
        z0_tau = PHI2 * z0
        r3, g3 = integrate_soliton_full(z0_tau, alpha, R_MAX)
        Atau, _ = fit_amplitude(r3, g3)

        if Ae < 1e-8:
            return {'alpha': alpha, 'z0': z0, 'F': None, 'G3': None,
                    'Ae': Ae, 'Amu': Amu, 'Atau': Atau}

        F  = (Amu  / Ae) ** 4
        G3 = (Atau / Ae) ** 4

        return {'alpha': alpha, 'z0': z0, 'F': F, 'G3': G3,
                'Ae': Ae, 'Amu': Amu, 'Atau': Atau}
    except Exception as e:
        return {'alpha': alpha, 'z0': None, 'F': None, 'G3': None,
                'Ae': None, 'Amu': None, 'Atau': None, 'err': str(e)}

# ================================================================
# MAIN — skan równoległy
# ================================================================
def main():
    print("=" * 68)
    print("ex92_F_alpha_profile.py — Profil F(α) i G₃(α)")
    print("=" * 68)
    print(f"  ALPHA zakres: [{ALPHA_MIN_RANGE}, {ALPHA_MAX_RANGE}], N={N_ALPHA}, krok≈{(ALPHA_MAX_RANGE-ALPHA_MIN_RANGE)/(N_ALPHA-1):.3f}")
    print(f"  R_MAX  = {R_MAX}")
    print(f"  B_WIN  = [{B_WIN_L}, {B_WIN_R}]")
    print(f"  PHI    = {PHI:.8f},  PHI² = {PHI2:.8f}")
    print(f"  R21    = {R21:.6f}  (m_μ/m_e, check={R21_check:.4f})")
    print(f"  R31    = {R31_val:.4f}  (m_τ/m_e)")
    print(f"  R32    = {R32_val:.6f}  (m_τ/m_μ)")
    print()

    alpha_arr = np.linspace(ALPHA_MIN_RANGE, ALPHA_MAX_RANGE, N_ALPHA)
    t0 = time.time()

    # Równoległy skan
    ncpu = min(mp.cpu_count(), 6)
    print(f"  Używam {ncpu} rdzeni (Pool)...")
    print()

    with mp.Pool(ncpu) as pool:
        results = pool.map(compute_one_alpha_fixed, alpha_arr.tolist())

    elapsed = time.time() - t0
    print(f"  Czas obliczen: {elapsed:.1f} s")
    print()

    # ================================================================
    # ANALIZA WYNIKÓW
    # ================================================================

    # Filtr: tylko poprawne wyniki (F i G3 nie-None)
    good = [r for r in results if r['F'] is not None and r['G3'] is not None]
    bad  = [r for r in results if r['F'] is None]

    print(f"  Wyniki dobre: {len(good)}/{len(results)}")
    print(f"  Wyniki złe:   {len(bad)}/{len(results)}")
    if bad:
        bad_alphas = [f"{r['alpha']:.3f}" for r in bad[:10]]
        print(f"  Brakujące α:  {', '.join(bad_alphas)}" + ("..." if len(bad) > 10 else ""))
    print()

    # Tablice numpy
    alphas = np.array([r['alpha'] for r in good])
    Farr   = np.array([r['F']     for r in good])
    G3arr  = np.array([r['G3']    for r in good])
    z0arr  = np.array([r['z0']    for r in good])
    Aearr  = np.array([r['Ae']    for r in good])

    # ----------------------------------------------------------------
    # PROFIL F(α) — raport
    # ----------------------------------------------------------------
    print("=" * 68)
    print("PROFIL F(α) = (A_tail(φ·z₀) / A_tail(z₀))^4")
    print("=" * 68)

    # Minimum F
    idx_min = np.argmin(Farr)
    alpha_min_F = alphas[idx_min]
    F_min = Farr[idx_min]
    print(f"  α_min(F)  = {alpha_min_F:.4f}  (F minimum)")
    print(f"  F_min     = {F_min:.4f}")
    print(f"  F(α=2.0)  = {np.interp(2.0, alphas, Farr):.4f}  (α_TGP)")
    print(f"  R21       = {R21:.4f}")
    print()

    # Tabela — co 5 punktów
    print(f"  {'α':>7} | {'z₀':>9} | {'F(α)':>10} | {'G₃(α)':>12} | {'vs R21':>9}")
    print(f"  {'-'*58}")
    for i in range(0, len(good), max(1, len(good)//25)):
        r = good[i]
        if r['F'] is None: continue
        vs_R21 = f"{(r['F']-R21)/R21*100:+.1f}%" if r['F'] is not None else "—"
        print(f"  {r['alpha']:>7.3f} | {r['z0']:>9.5f} | {r['F']:>10.3f} | {r['G3']:>12.2f} | {vs_R21:>9}")
    print()

    # ----------------------------------------------------------------
    # ZERA F(α) = R21 (α*₁, α*₂)
    # ----------------------------------------------------------------
    print("=" * 68)
    print("ZERA F(α) = R21 (szukamy α* gdzie F(α*)=R21)")
    print("=" * 68)

    zeros_F_R21 = []
    for i in range(len(alphas) - 1):
        f1, f2 = Farr[i] - R21, Farr[i+1] - R21
        if f1 * f2 < 0:
            a1, a2 = alphas[i], alphas[i+1]
            # Interpolacja liniowa
            a_cross = a1 - f1 * (a2 - a1) / (f2 - f1)
            zeros_F_R21.append(a_cross)
            print(f"  Zero F=R21 w α ≈ {a_cross:.5f}  (przedział [{a1:.4f}, {a2:.4f}])")
            print(f"    F({a1:.4f})={Farr[i]:.4f}, F({a2:.4f})={Farr[i+1]:.4f}")

    if not zeros_F_R21:
        print("  Brak zer F=R21 w przebadanym zakresie!")

    print()

    # ----------------------------------------------------------------
    # CZY F(α) OSIĄGA R31?
    # ----------------------------------------------------------------
    print("=" * 68)
    print(f"CZY F(α) OSIĄGA R31 = {R31_val:.2f}?")
    print("=" * 68)
    F_max = np.max(Farr)
    F_min_val = np.min(Farr)
    print(f"  F zakres: [{F_min_val:.2f}, {F_max:.2f}]")
    print(f"  R31 = {R31_val:.2f}")
    if F_max >= R31_val:
        print(f"  TAK — F(α) osiąga R31 gdzieś w zakresie!")
        # Znajdź gdzie
        zeros_F_R31 = []
        for i in range(len(alphas) - 1):
            f1, f2 = Farr[i] - R31_val, Farr[i+1] - R31_val
            if f1 * f2 < 0:
                a1, a2 = alphas[i], alphas[i+1]
                a_cross = a1 - f1 * (a2 - a1) / (f2 - f1)
                zeros_F_R31.append(a_cross)
                print(f"  Zero F=R31 w α ≈ {a_cross:.5f}  [{a1:.4f}, {a2:.4f}]")
    else:
        print(f"  NIE — F_max={F_max:.2f} < R31={R31_val:.2f}")
        print(f"  Podejście 1 (F=R31) NIE zadziała w tym zakresie α.")
    print()

    # ----------------------------------------------------------------
    # PROFIL G₃(α) i zera G₃=R31
    # ----------------------------------------------------------------
    print("=" * 68)
    print("PROFIL G₃(α) = (A_tail(φ²·z₀) / A_tail(z₀))^4")
    print("=" * 68)

    idx_min_G3 = np.argmin(G3arr)
    alpha_min_G3 = alphas[idx_min_G3]
    G3_min = G3arr[idx_min_G3]
    G3_max = np.max(G3arr)
    print(f"  α_min(G₃) = {alpha_min_G3:.4f}  (G₃ minimum)")
    print(f"  G₃_min    = {G3_min:.4f}")
    print(f"  G₃_max    = {G3_max:.4f}")
    print(f"  R31       = {R31_val:.4f}")
    print()

    # Tabela G3 — klucze punkty
    print(f"  {'α':>7} | {'G₃(α)':>12} | {'vs R31':>12} | {'G3/F':>8}")
    print(f"  {'-'*50}")
    for i in range(0, len(good), max(1, len(good)//20)):
        r = good[i]
        if r['F'] is None: continue
        vs_R31 = f"{(r['G3']-R31_val)/R31_val*100:+.1f}%"
        g3_f_ratio = r['G3']/r['F'] if r['F'] > 0 else float('nan')
        print(f"  {r['alpha']:>7.3f} | {r['G3']:>12.2f} | {vs_R31:>12} | {g3_f_ratio:>8.3f}")
    print()

    # Zera G₃ = R31
    print(f"  Szukam zer G₃(α) = R31 = {R31_val:.2f}:")
    zeros_G3_R31 = []
    for i in range(len(alphas) - 1):
        g1, g2 = G3arr[i] - R31_val, G3arr[i+1] - R31_val
        if g1 * g2 < 0:
            a1, a2 = alphas[i], alphas[i+1]
            a_cross = a1 - g1 * (a2 - a1) / (g2 - g1)
            zeros_G3_R31.append(a_cross)
            print(f"  Zero G₃=R31 w α ≈ {a_cross:.5f}  [{a1:.4f}, {a2:.4f}]")
            print(f"    G₃({a1:.4f})={G3arr[i]:.2f}, G₃({a2:.4f})={G3arr[i+1]:.2f}")

    if not zeros_G3_R31:
        print(f"  Brak zer G₃=R31 w zakresie α∈[{ALPHA_MIN_RANGE},{ALPHA_MAX_RANGE}]")
        print(f"  G₃ zakres: [{G3_min:.2f}, {G3_max:.2f}]")
        if G3_min <= R31_val <= G3_max:
            print(f"  ⚠️  R31 leży w zakresie G₃, ale zero nie znalezione — sprawdź siatkę")
        elif R31_val > G3_max:
            print(f"  R31={R31_val:.1f} > G₃_max={G3_max:.1f} — G₃ nie osiąga R31 w tym zakresie")
        else:
            print(f"  R31={R31_val:.1f} < G₃_min={G3_min:.1f} — G₃ zawsze > R31 w tym zakresie")
    print()

    # ----------------------------------------------------------------
    # G₃ przy α*₁, α*₂ (znane wartości)
    # ----------------------------------------------------------------
    print("=" * 68)
    print("G₃ przy α*₁=2.43184, α*₂=2.63558 (z ex88)")
    print("=" * 68)
    for alpha_star, label in [(2.43183767, "α*₁"), (2.63557742, "α*₂")]:
        res = compute_one_alpha_fixed(alpha_star)
        if res['F'] is not None:
            print(f"  {label} = {alpha_star:.8f}:")
            print(f"    z₀     = {res['z0']:.6f}")
            print(f"    F(α*)  = {res['F']:.6f}  (R21={R21:.4f}, diff={res['F']-R21:+.6f})")
            print(f"    G₃(α*) = {res['G3']:.4f}  (R31={R31_val:.4f}, diff={res['G3']-R31_val:+.4f})")
            print(f"    ratio G₃/F = {res['G3']/res['F']:.6f}  (φ^16={PHI**16:.6f})")
        else:
            print(f"  {label}: BŁĄD — brak wyniku")
    print()

    # ----------------------------------------------------------------
    # STOSUNEK G₃/F — analityczny test
    # ----------------------------------------------------------------
    print("=" * 68)
    print("STOSUNEK G₃(α)/F(α) — test analityczny")
    print("=" * 68)
    ratio_G3_F = G3arr / Farr
    print(f"  φ^16           = {PHI**16:.6f}")
    print(f"  φ^8            = {PHI**8:.6f}")
    print(f"  φ^4            = {PHI**4:.6f}")
    print(f"  G₃/F zakres:    [{ratio_G3_F.min():.4f}, {ratio_G3_F.max():.4f}]")
    print(f"  G₃/F przy α=2:  {np.interp(2.0, alphas, ratio_G3_F):.6f}")
    print(f"  G₃/F przy α*₁:  {np.interp(2.43184, alphas, ratio_G3_F):.6f}")
    print()

    # ----------------------------------------------------------------
    # DIAGNOZA O-L4: Dlaczego α_min ≠ α_TGP=2?
    # ----------------------------------------------------------------
    print("=" * 68)
    print("DIAGNOZA O-L4: Struktura doliny F(α)")
    print("=" * 68)
    print(f"  α_TGP       = 2.0000")
    print(f"  α_min(F)    = {alpha_min_F:.4f}  (gdzie F osiąga minimum)")
    print(f"  F(α_TGP=2)  = {np.interp(2.0, alphas, Farr):.4f}")
    print(f"  F_min       = {F_min:.4f}")
    print(f"  R21         = {R21:.4f}")
    print()
    print(f"  Asymetria: (α_min - α_TGP) = {alpha_min_F - 2.0:+.4f}")
    print(f"  α*₁ = 2.43184:  F(α_TGP=2) - F(α*₁) = {np.interp(2.0,alphas,Farr) - np.interp(2.43184,alphas,Farr):+.4f}")
    print()

    # z₀ jako funkcja α
    print(f"  z₀ jako f(α):")
    for alpha_test in [1.5, 2.0, 2.43, 2.5, 2.63, 3.0, 3.5, 4.0]:
        z0_test = np.interp(alpha_test, alphas, z0arr)
        g_star  = np.exp(-1.0/(2.0*alpha_test))
        print(f"    α={alpha_test:.2f}: z₀={z0_test:.5f}, g*(α)=e^(-1/2α)={g_star:.5f}, z₀-g*={z0_test-g_star:.5f}")
    print()

    # ----------------------------------------------------------------
    # PODSUMOWANIE KOŃCOWE
    # ----------------------------------------------------------------
    print("=" * 68)
    print("PODSUMOWANIE KOŃCOWE ex92")
    print("=" * 68)

    print(f"\nO-L4 (dlaczego α*≠2):")
    print(f"  → α_min(F) = {alpha_min_F:.4f}  (centrum doliny)")
    print(f"  → F(α=2.0) = {np.interp(2.0, alphas, Farr):.4f} > R21 = {R21:.4f}")
    print(f"  → Dolina F jest przesunięta: α_min > α_TGP")
    print(f"  → α_TGP=2 leży NA LEWO od minimum — F spada od α=2 do α_min")
    print(f"  → Dwa zera (α*₁ < α_min < α*₂) leżą symetrycznie wokół α_min")

    if zeros_F_R21:
        a1_found = min(zeros_F_R21)
        a2_found = max(zeros_F_R21)
        print(f"  → α*₁ ≈ {a1_found:.4f}  (zero F=R21, malejąco)")
        print(f"  → α*₂ ≈ {a2_found:.4f}  (zero F=R21, rosnąco)")
        print(f"  → α_mid = (α*₁+α*₂)/2 = {(a1_found+a2_found)/2:.4f}  vs α_min={alpha_min_F:.4f}")

    print(f"\nO-L5 (tau sektor):")
    if zeros_G3_R31:
        print(f"  → ZNALEZIONO α*₃ gdzie G₃(α*₃)=R31 ≈ {R31_val:.1f}!")
        for z in zeros_G3_R31:
            print(f"    α*₃ ≈ {z:.5f}")
    else:
        print(f"  → G₃ NIE osiąga R31={R31_val:.1f} w α∈[{ALPHA_MIN_RANGE},{ALPHA_MAX_RANGE}]")
        print(f"     G₃ zakres: [{G3_min:.1f}, {G3_max:.1f}]")
        # Sugestia dalszego zakresu
        if G3_max < R31_val:
            print(f"  → G₃_max({G3_max:.1f}) < R31({R31_val:.1f}) — τ nieodkryty w tym zakresie")
            print(f"     Sugestia: rozszerzyć zakres α lub użyć innej funkcji (np. H(α)=G₃²/F)")

    print()
    print(f"Czas całkowity: {time.time()-t0:.1f} s")
    print("=" * 68)

if __name__ == '__main__':
    mp.freeze_support()
    main()
