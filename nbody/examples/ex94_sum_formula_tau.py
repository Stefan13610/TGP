#!/usr/bin/env python3
"""
ex94_sum_formula_tau.py — Weryfikacja S₃=8 i G₃/F=R32 (O-L5, v37)
====================================================================
Cele:
  1. Precyzyjne α*τμ gdzie G₃(α)/F(α) = R32 = m_τ/m_μ ≈ 16.817
     (Z ex93: crossing między α=2.840 a α=2.850)
  2. Precyzyjne α*₃ przy R_MAX=300 dla najwyższej dokładności
  3. Analiza algebraiczna sumy S₃ = α*₁+α*₂+α*₃:
     - S₃ vs 8 (= 4·α_TGP)
     - S₃ vs inne kombinacje TGP
  4. Testy kandydatów algebraicznych α*₃:
     - α*₃ = 8 − α*₁ − α*₂?
     - α*₃ = (α*₁+α*₂)·k dla prostego k?
     - α*₃ relacja z φ, π, e?

Z ex93 (wartości finalne):
  α*₁ = 2.43183767 (ex88, σ=0)
  α*₂ = 2.63557742 (ex88, σ=0)
  α*₃ ≈ 2.92952433 (ex93, σ=0 dla R_MAX=120,150,200)
  G₃/F przy α=2.840: 17.342 (+3.12%)
  G₃/F przy α=2.850: 16.800 (−0.10%) ← crossing tu!

Autor: Claudian (sesja v37, O-L5)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import multiprocessing as mp
import time
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0
PHI2 = PHI * PHI

M_TAU_MEV = 1776.86
M_MU_MEV  = 105.658
M_E_MEV   = 0.51100

R21  = M_MU_MEV / M_E_MEV    # ≈ 206.767
R31  = M_TAU_MEV / M_E_MEV   # ≈ 3477.221
R32  = M_TAU_MEV / M_MU_MEV  # ≈ 16.817

ALPHA_TGP = 2.0

# Znane wartości z ex88/ex93
ALPHA_STAR1 = 2.43183767
ALPHA_STAR2 = 2.63557742
ALPHA_STAR3_EX93 = 2.92952433  # ex93, R_MAX∈{120,150,200}, σ=0

G_OFF    = 0.005
B_WIN_L  = 28.0
B_WIN_R  = 42.0
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]

# ================================================================
# ODE + AMPLITUDY (identyczne z ex91/ex92/ex93)
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
    for _ in range(20):
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
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, lo=1.05, hi=2.0, n=60):
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

def compute_all(alpha, r_max):
    z0 = find_z0(alpha, r_max)
    if z0 is None:
        return None

    r,  g  = integrate_soliton(z0,         alpha, r_max)
    Ae, _  = fit_amplitude(r, g)

    r2, g2 = integrate_soliton(PHI*z0,     alpha, r_max)
    Amu, _ = fit_amplitude(r2, g2)

    r3, g3 = integrate_soliton(PHI2*z0,    alpha, r_max)
    Atau, _ = fit_amplitude(r3, g3)

    if Ae < 1e-8:
        return None

    F   = (Amu  / Ae) ** 4
    G3  = (Atau / Ae) ** 4
    G3F = G3 / F if F > 0 else float('nan')

    return {'alpha': alpha, 'r_max': r_max, 'z0': z0,
            'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
            'F': F, 'G3': G3, 'G3F': G3F}

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    print("=" * 68)
    print("ex94_sum_formula_tau.py — S₃=8? i G₃/F=R32 (O-L5 v37)")
    print("=" * 68)
    print(f"  R21={R21:.6f}, R31={R31:.6f}, R32={R32:.6f}")
    print(f"  α*₁={ALPHA_STAR1}, α*₂={ALPHA_STAR2}, α*₃_ex93={ALPHA_STAR3_EX93}")
    print(f"  S₂=α*₁+α*₂ = {ALPHA_STAR1+ALPHA_STAR2:.8f}")
    print(f"  S₃_ex93    = {ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93:.8f}")
    print(f"  8-S₃_ex93  = {8-(ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93):+.8f}  "
          f"({(8-(ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93))/8*1e6:+.1f} ppm)")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Precyzyjne α*τμ gdzie G₃/F = R32
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 1: brentq dla G₃(α)/F(α) = R32 = m_τ/m_μ")
    print("=" * 68)
    print(f"  R32 = {R32:.8f}")
    print(f"  Z ex93: crossing między α=2.840 (G₃/F=17.342) a α=2.850 (G₃/F=16.800)")
    print()

    def G3F_minus_R32(alpha, r_max=120):
        res = compute_all(alpha, r_max)
        if res is None: return float('nan')
        return res['G3F'] - R32

    # Weryfikacja bracketa
    f840 = G3F_minus_R32(2.840)
    f850 = G3F_minus_R32(2.850)
    print(f"  G₃/F(2.840)-R32 = {f840:+.6f}")
    print(f"  G₃/F(2.850)-R32 = {f850:+.6f}")

    alpha_taumu = None
    if f840 * f850 < 0:
        alpha_taumu = brentq(G3F_minus_R32, 2.840, 2.850, xtol=1e-9, maxiter=200)
        res_taumu = compute_all(alpha_taumu, 120)
        print(f"\n  α*τμ = {alpha_taumu:.9f}")
        print(f"  G₃/F(α*τμ) = {res_taumu['G3F']:.8f}  (R32={R32:.6f})")
        print(f"  F(α*τμ)    = {res_taumu['F']:.6f}  (R21={R21:.4f})")
        print(f"  G₃(α*τμ)   = {res_taumu['G3']:.4f}  (R31={R31:.4f})")
        print(f"  z₀(α*τμ)   = {res_taumu['z0']:.6f}")

        # Weryfikacja R_MAX
        print(f"\n  Stabilność α*τμ vs R_MAX:")
        rmax_taumu = []
        for rmax in [120, 150, 200]:
            def f_rmax(a): return G3F_minus_R32(a, rmax)
            try:
                f_lo = f_rmax(alpha_taumu - 0.02)
                f_hi = f_rmax(alpha_taumu + 0.02)
                if f_lo * f_hi < 0:
                    a_r = brentq(f_rmax, alpha_taumu-0.02, alpha_taumu+0.02, xtol=1e-9)
                    res_r = compute_all(a_r, rmax)
                    print(f"    R_MAX={rmax}: α*τμ={a_r:.9f}, G₃/F={res_r['G3F']:.6f}")
                    rmax_taumu.append(a_r)
                else:
                    print(f"    R_MAX={rmax}: brak bracketa — G₃/F({alpha_taumu-0.02:.3f})={f_lo:+.4f}, "
                          f"G₃/F({alpha_taumu+0.02:.3f})={f_hi:+.4f}")
            except Exception as e:
                print(f"    R_MAX={rmax}: błąd — {e}")

        if len(rmax_taumu) >= 2:
            sigma_tm = np.std(rmax_taumu)
            print(f"  Spread α*τμ: {max(rmax_taumu)-min(rmax_taumu):.2e}  (σ={sigma_tm:.2e})")
    else:
        print(f"  Brak bracketa! f840*f850 > 0 — sprawdź zakres")
        # Gęsty skan pomocniczy [2.83, 2.87]
        print(f"  Skan pomocniczy [2.83, 2.87], krok=0.002:")
        for a in np.arange(2.83, 2.88, 0.002):
            v = G3F_minus_R32(a)
            print(f"    α={a:.3f}: G₃/F-R32 = {v:+.4f}")

    print()

    # ----------------------------------------------------------------
    # FAZA 2: Precyzyjne α*₃ przy R_MAX=300
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 2: Precyzyjne α*₃ przy R_MAX=300")
    print("=" * 68)

    def G3_minus_R31_r(alpha, r_max):
        res = compute_all(alpha, r_max)
        if res is None: return float('nan')
        return res['G3'] - R31

    alpha_star3_300 = None
    f920_300 = G3_minus_R31_r(2.920, 300)
    f935_300 = G3_minus_R31_r(2.935, 300)
    print(f"  G₃(2.920,300)-R31 = {f920_300:+.4f}")
    print(f"  G₃(2.935,300)-R31 = {f935_300:+.4f}")

    if not np.isnan(f920_300) and not np.isnan(f935_300) and f920_300 * f935_300 < 0:
        alpha_star3_300 = brentq(lambda a: G3_minus_R31_r(a, 300),
                                 2.920, 2.935, xtol=1e-9, maxiter=200)
        res3_300 = compute_all(alpha_star3_300, 300)
        print(f"\n  α*₃ (R_MAX=300) = {alpha_star3_300:.9f}")
        print(f"  G₃(α*₃)         = {res3_300['G3']:.6f}  (R31={R31:.4f}, diff={res3_300['G3']-R31:+.4f})")
        print(f"  F(α*₃)          = {res3_300['F']:.6f}")
        print(f"  z₀              = {res3_300['z0']:.6f}")
        print(f"\n  Porównanie α*₃:")
        print(f"    ex93 (R_MAX=120,150,200): {ALPHA_STAR3_EX93:.9f}")
        print(f"    ex94 (R_MAX=300):         {alpha_star3_300:.9f}")
        print(f"    Różnica: {alpha_star3_300-ALPHA_STAR3_EX93:+.2e}")
    else:
        print(f"  Brak bracketa przy R_MAX=300! Skan pomocniczy:")
        for a in np.arange(2.90, 2.96, 0.005):
            v = G3_minus_R31_r(a, 300)
            print(f"    α={a:.3f}: G₃-R31={v:+.2f}")
        alpha_star3_300 = ALPHA_STAR3_EX93  # fallback

    print()

    # ----------------------------------------------------------------
    # FAZA 3: Analiza algebraiczna S₃
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 3: Analiza algebraiczna sumy S₃ = α*₁+α*₂+α*₃")
    print("=" * 68)

    # Użyj najlepszej wartości α*₃
    best_alpha3 = alpha_star3_300 if alpha_star3_300 else ALPHA_STAR3_EX93
    S3 = ALPHA_STAR1 + ALPHA_STAR2 + best_alpha3
    S2 = ALPHA_STAR1 + ALPHA_STAR2
    D2 = ALPHA_STAR2 - ALPHA_STAR1
    D3 = best_alpha3 - ALPHA_STAR2

    print(f"  Wartości wejściowe:")
    print(f"    α*₁ = {ALPHA_STAR1:.8f}")
    print(f"    α*₂ = {ALPHA_STAR2:.8f}")
    print(f"    α*₃ = {best_alpha3:.8f}")
    print(f"    S₂  = α*₁+α*₂ = {S2:.8f}")
    print(f"    S₃  = α*₁+α*₂+α*₃ = {S3:.8f}")
    print(f"    D₂  = α*₂−α*₁ = {D2:.8f}")
    print(f"    D₃  = α*₃−α*₂ = {D3:.8f}")
    print(f"    D₃/D₂ = {D3/D2:.8f}")
    print()

    # Kandydaci dla S₃
    candidates_S3 = [
        ("4·α_TGP = 8",          8.0),
        ("3π − 1",               3*np.pi - 1),
        ("2π + 3/2",             2*np.pi + 1.5),
        ("2π + φ",               2*np.pi + PHI),
        ("4π − 5",               4*np.pi - 5),
        ("e²",                   np.e**2),
        ("π + e + 1/2",          np.pi + np.e + 0.5),
        ("5φ − 1/10",            5*PHI - 0.1),
        ("3·α*₁ + 3/10",         3*ALPHA_STAR1 + 0.3),
        ("α*₁·φ + α*₂ + α*₂",   ALPHA_STAR1*PHI + 2*ALPHA_STAR2),
    ]

    print(f"  Kandydaci dla S₃ = {S3:.8f}:")
    print(f"  {'Formuła':>30} | {'Wartość':>12} | {'Odchylenie':>14} | {'ppm':>10}")
    print(f"  {'-'*75}")
    for name, val in candidates_S3:
        diff = S3 - val
        ppm  = diff / val * 1e6
        marker = " ←" if abs(ppm) < 1000 else ""
        print(f"  {name:>30} | {val:>12.8f} | {diff:>+14.8f} | {ppm:>+10.1f}{marker}")

    print()

    # Kandydaci dla α*₃
    print(f"  Kandydaci dla α*₃ = {best_alpha3:.8f}:")
    print(f"  {'Formuła':>35} | {'Wartość':>12} | {'ppm':>10}")
    print(f"  {'-'*65}")
    candidates_a3 = [
        ("8 − S₂ = 8 − α*₁ − α*₂",    8.0 - S2),
        ("3 − 1/10",                    2.9),
        ("3π/α_TGP − 1",               3*np.pi/2 - 1),
        ("(S₂ + 1)/2",                  (S2+1)/2),
        ("φ² − 1/3",                    PHI2 - 1.0/3.0),
        ("α*₁ · φ²/φ",                  ALPHA_STAR1 * PHI),
        ("α*₂ + 1/3",                   ALPHA_STAR2 + 1.0/3.0),
        ("α*₁ + α*₁/(φ²)",             ALPHA_STAR1 + ALPHA_STAR1/PHI2),
        ("(α*₁+α*₂) · φ/√5",           S2 * PHI / np.sqrt(5)),
        ("π − 1/4",                     np.pi - 0.25),
    ]
    for name, val in candidates_a3:
        diff = best_alpha3 - val
        ppm  = diff / val * 1e6
        marker = " ←" if abs(ppm) < 1000 else ""
        print(f"  {name:>35} | {val:>12.8f} | {ppm:>+10.1f}{marker}")

    print()

    # Stosunek D₃/D₂
    print(f"  Analiza różnic:")
    print(f"    D₂ = α*₂−α*₁ = {D2:.8f}")
    print(f"    D₃ = α*₃−α*₂ = {D3:.8f}")
    print(f"    D₃/D₂ = {D3/D2:.8f}")
    cands_ratio = [
        ("φ",   PHI),
        ("3/2", 1.5),
        ("√2",  np.sqrt(2)),
        ("4/3", 4.0/3.0),
        ("φ²",  PHI2),
        ("e/2", np.e/2),
    ]
    print(f"  Kandydaci D₃/D₂:")
    for name, val in cands_ratio:
        ppm = (D3/D2 - val)/val*1e6
        marker = " ←" if abs(ppm) < 5000 else ""
        print(f"    D₃/D₂ vs {name:>5} = {val:.6f}: {ppm:>+10.1f} ppm{marker}")

    print()

    # ----------------------------------------------------------------
    # FAZA 4: Podsumowanie
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 4: Podsumowanie ex94")
    print("=" * 68)

    if alpha_taumu:
        print(f"\n  G₃/F = R32 = m_τ/m_μ:")
        print(f"    α*τμ = {alpha_taumu:.9f}")
        print(f"    Interpretacja: stosunek amplitud τ/μ odtwarza stosunek mas")
        print(f"    F(α*τμ) = {compute_all(alpha_taumu,120)['F']:.4f}  vs R21={R21:.4f}  "
              f"(nie jest zero F=R21)")

    print(f"\n  Suma S₃ = α*₁+α*₂+α*₃:")
    print(f"    S₃ = {S3:.8f}")
    best_S3 = min(candidates_S3, key=lambda x: abs(S3-x[1]))
    print(f"    Najlepszy kandydat: {best_S3[0]} = {best_S3[1]:.8f}")
    print(f"    Odchylenie: {(S3-best_S3[1])/best_S3[1]*1e6:+.1f} ppm")

    best_a3 = min(candidates_a3, key=lambda x: abs(best_alpha3-x[1]))
    print(f"\n  Najlepszy kandydat algebraiczny dla α*₃:")
    print(f"    α*₃ ≈ {best_a3[0]} = {best_a3[1]:.8f}")
    print(f"    Odchylenie: {(best_alpha3-best_a3[1])/best_a3[1]*1e6:+.1f} ppm")

    print(f"\nCzas całkowity: {time.time()-t0:.1f} s")
    print("=" * 68)


if __name__ == '__main__':
    mp.freeze_support()
    main()
