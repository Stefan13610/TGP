#!/usr/bin/env python3
"""
ex95_G3_precision.py — Lepsza precyzja G₃ przez rozszerzone okna (O-L6, v38)
==============================================================================
STATUS: LEGACY-TRANSLATIONAL

This file remains in the older tau-sector `G3` / `alpha*` exploration built on
legacy selection language. Keep it as historical exploratory material rather
than a canonical current path.

Problem z ex93/ex94:
  G₃(α*₃) = 3495.2 vs R31=3477.2 → różnica +18 j. (+0.52%)
  Szum G₃ ≈ ±20 j. → niepewność α*₃ ≈ ±0.007 → szum S₃ ≈ ±875 ppm
  Formuła S₃=8 ma odchylenie 383 ppm — niemożliwa do potwierdzenia przy obecnym szumie.

Hipoteza:
  Obecne FIT_WINS = [(20,36),(30,46),(40,56),(50,66)] kończą się na r=66.
  Dla φ²·z₀≈3.3 (duże g₀) ekstrapolacja A∞ może być niestabilna przy małych oknach.
  Rozszerzone okna do r≈100 powinny dać lepszą ekstrapolację.

Plan:
  FAZA 1: Porównanie A_tail dla gałęzi elektron/muon/tau przy różnych zestawach okien
  FAZA 2: G₃ na siatce α∈[2.90, 2.95] (krok 0.001) z oknami standardowymi vs rozszerzonymi
  FAZA 3: Precyzyjny brentq dla α*₃ z rozszerzonymi oknami + R_MAX=150
  FAZA 4: Obliczenie S₃ z nową precyzją i test S₃=8

Autor: Claudian (sesja v38, O-L6)
Data: 2026-03-29
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

ALPHA_STAR1 = 2.43183767
ALPHA_STAR2 = 2.63557742
ALPHA_STAR3_EX93 = 2.929524330   # ex93, R_MAX=120,150,200
ALPHA_STAR3_EX94 = 2.929524307   # ex94, R_MAX=300

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0
R_MAX   = 150   # wyższy niż 120 aby okna do r=106 zmieściły się

# Zestawy okien do porównania
FIT_WINS_STD = [(20, 36), (30, 46), (40, 56), (50, 66)]
FIT_WINS_EXT = [(20, 36), (30, 46), (40, 56), (50, 66),
                (60, 76), (70, 86), (80, 96), (90, 106)]
FIT_WINS_FAR = [(60, 76), (70, 86), (80, 96), (90, 106),
                (100, 116), (110, 126), (120, 136)]  # tylko dalekie okna (wymaga R_MAX=150)

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

def fit_amplitude_windows(r, g, fit_wins):
    """Ekstrapolacja A_∞ z zadanych okien."""
    Av, rv = [], []
    for rL, rR in fit_wins:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0, 0.0, []
    if len(Av) < 2: return Av[-1], 0.0, list(zip(rv, Av))
    rv = np.array(rv); Av_arr = np.array(Av)
    try:
        p, pcov = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                            rv, Av_arr, p0=[Av_arr[-1], 0., 0.], maxfev=2000)
        err = float(np.sqrt(abs(np.diag(pcov)[0])))
        return float(p[0]), err, list(zip(rv.tolist(), Av_arr.tolist()))
    except:
        try:
            p, _ = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                             rv, Av_arr, p0=[Av_arr[-1], 0.], maxfev=1000)
            return float(p[0]), 0.0, list(zip(rv.tolist(), Av_arr.tolist()))
        except:
            return float(Av_arr[-1]), 0.0, list(zip(rv.tolist(), Av_arr.tolist()))

def B_coeff_func(g0, alpha, r_max):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX, lo=1.05, hi=2.0, n=60):
    gs = np.linspace(lo, hi, n)
    Bs = []
    for g in gs:
        try: Bs.append(B_coeff_func(g, alpha, r_max))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff_func(x, alpha, r_max),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def compute_G3_compare(alpha, r_max=R_MAX):
    """Oblicz G₃ z obydwoma zestawami okien."""
    z0 = find_z0(alpha, r_max)
    if z0 is None:
        return None

    # Solitony
    r,  g  = integrate_soliton(z0,       alpha, r_max)
    r3, g3 = integrate_soliton(PHI2*z0,  alpha, r_max)
    r2, g2 = integrate_soliton(PHI*z0,   alpha, r_max)

    results = {}
    for label, wins in [('std', FIT_WINS_STD), ('ext', FIT_WINS_EXT), ('far', FIT_WINS_FAR)]:
        Ae,  dAe,  wAe  = fit_amplitude_windows(r,  g,  wins)
        Atau, dAtau, wAtau = fit_amplitude_windows(r3, g3, wins)
        Amu, dAmu, wAmu = fit_amplitude_windows(r2, g2, wins)
        if Ae > 1e-8:
            F  = (Amu  / Ae) ** 4
            G3 = (Atau / Ae) ** 4
        else:
            F, G3 = float('nan'), float('nan')
        results[label] = {'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
                          'F': F, 'G3': G3, 'dAe': dAe, 'dAtau': dAtau}

    return {'alpha': alpha, 'z0': z0, 'r_max': r_max, **results}

def worker_compare(args):
    alpha, r_max = args
    return compute_G3_compare(alpha, r_max)

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    print("=" * 72)
    print("ex95_G3_precision.py — Lepsza precyzja G₃ (O-L6, v38)")
    print("=" * 72)
    print(f"  R_MAX     = {R_MAX}")
    print(f"  R21={R21:.4f}, R31={R31:.4f}, R32={R32:.4f}")
    print(f"  α*₁={ALPHA_STAR1}, α*₂={ALPHA_STAR2}")
    print(f"  α*₃(ex93)={ALPHA_STAR3_EX93}, α*₃(ex94)={ALPHA_STAR3_EX94}")
    print(f"  S₃(ex93)  = {ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93:.8f}")
    print(f"  8−S₃(ex93)= {8-(ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93):+.8f} ({(8-(ALPHA_STAR1+ALPHA_STAR2+ALPHA_STAR3_EX93))/8*1e6:+.1f} ppm)")
    print()
    print(f"  FIT_WINS_STD = {FIT_WINS_STD}")
    print(f"  FIT_WINS_EXT = {FIT_WINS_EXT}")
    print(f"  FIT_WINS_FAR = {FIT_WINS_FAR}")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Diagnostyka okien przy α*₃
    # ----------------------------------------------------------------
    print("=" * 72)
    print(f"FAZA 1: Diagnostyka okien dla amplitudy przy α*₃={ALPHA_STAR3_EX93:.6f}")
    print("=" * 72)

    res_diag = compute_G3_compare(ALPHA_STAR3_EX93, R_MAX)
    if res_diag:
        print(f"  z₀(α*₃) = {res_diag['z0']:.6f}")
        print(f"  φ·z₀    = {PHI*res_diag['z0']:.6f}")
        print(f"  φ²·z₀   = {PHI2*res_diag['z0']:.6f}")
        print()
        print(f"  {'Okna':>4} | {'Ae':>10} | {'Amu':>10} | {'Atau':>10} | {'F':>9} | {'G3':>10} | {'G3-R31':>10}")
        print(f"  {'-'*72}")
        for lab in ['std','ext','far']:
            d = res_diag[lab]
            g3r31 = d['G3'] - R31 if not np.isnan(d['G3']) else float('nan')
            print(f"  {lab:>4} | {d['Ae']:>10.6f} | {d['Amu']:>10.6f} | {d['Atau']:>10.6f} | "
                  f"{d['F']:>9.4f} | {d['G3']:>10.4f} | {g3r31:>+10.4f}")
        print()

        # Profil A_tail(r) dla Atau (tau sector)
        z0_v = res_diag['z0']
        r3, g3 = integrate_soliton(PHI2*z0_v, ALPHA_STAR3_EX93, R_MAX)
        print(f"  Profil A_tau per okno (φ²·z₀={PHI2*z0_v:.5f}):")
        for rL, rR in FIT_WINS_EXT + [(100,116),(110,126),(120,136)]:
            if rR > r3[-1]: break
            mask = (r3 >= rL) & (r3 <= rR)
            if mask.sum() < 15: continue
            rf = r3[mask]; df = (g3[mask]-1.0)*rf
            M = np.column_stack([np.cos(rf), np.sin(rf)])
            bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
            A_win = float(np.sqrt(bc[0]**2+bc[1]**2))
            print(f"    [{rL:3d},{rR:3d}]: A_tau={A_win:.6f}  (points={mask.sum()})")
    print()

    # ----------------------------------------------------------------
    # FAZA 2: Gęsty skan G₃ na α∈[2.90, 2.95] krok=0.001
    # ----------------------------------------------------------------
    print("=" * 72)
    print("FAZA 2: Gęsty skan G₃(α) na [2.90, 2.95], krok=0.001")
    print("=" * 72)

    alpha_arr = np.arange(2.900, 2.951, 0.001)
    args_list = [(a, R_MAX) for a in alpha_arr]
    ncpu = min(mp.cpu_count(), 6)
    print(f"  {len(alpha_arr)} punktów, {ncpu} rdzeni...")

    with mp.Pool(ncpu) as pool:
        scan_res = pool.map(worker_compare, args_list)

    good = [r for r in scan_res if r is not None]
    print(f"  Wyniki: {len(good)}/{len(scan_res)}")
    print()

    print(f"  {'α':>7} | {'z₀':>8} | {'F_std':>9} | {'G3_std':>10} | "
          f"{'G3_ext':>10} | {'G3_far':>10} | {'G3std-R31':>10}")
    print(f"  {'-'*75}")
    for r in good:
        g3s   = r['std']['G3']
        g3e   = r['ext']['G3']
        g3f   = r['far']['G3']
        fs    = r['std']['F']
        diff  = g3s - R31 if not np.isnan(g3s) else float('nan')
        marker = " ←" if abs(diff) < 25 else ""
        print(f"  {r['alpha']:>7.3f} | {r['z0']:>8.5f} | {fs:>9.4f} | {g3s:>10.4f} | "
              f"{g3e:>10.4f} | {g3f:>10.4f} | {diff:>+10.4f}{marker}")
    print()

    # ----------------------------------------------------------------
    # FAZA 3: brentq z rozszerzonymi oknami
    # ----------------------------------------------------------------
    print("=" * 72)
    print("FAZA 3: brentq α*₃ z oknami EXT i FAR")
    print("=" * 72)

    def G3_minus_R31(alpha, wins_label, r_max=R_MAX):
        res = compute_G3_compare(alpha, r_max)
        if res is None: return float('nan')
        g3 = res[wins_label]['G3']
        return g3 - R31 if not np.isnan(g3) else float('nan')

    results_astar3 = {}
    for wins_label in ['std', 'ext', 'far']:
        print(f"\n  --- Okna: {wins_label} ---")
        # Znajdź bracket z poprzednich skanów
        f920 = G3_minus_R31(2.920, wins_label)
        f935 = G3_minus_R31(2.935, wins_label)
        print(f"  G3_{wins_label}(2.920)-R31 = {f920:+.4f}")
        print(f"  G3_{wins_label}(2.935)-R31 = {f935:+.4f}")

        if not np.isnan(f920) and not np.isnan(f935) and f920*f935 < 0:
            try:
                astar = brentq(lambda a: G3_minus_R31(a, wins_label),
                               2.920, 2.935, xtol=1e-8, maxiter=100)
                res_v = compute_G3_compare(astar, R_MAX)
                g3_v  = res_v[wins_label]['G3'] if res_v else float('nan')
                f_v   = res_v['std']['F'] if res_v else float('nan')
                print(f"  α*₃ ({wins_label}) = {astar:.9f}")
                print(f"  G3  = {g3_v:.4f}  (R31={R31:.4f}, diff={g3_v-R31:+.4f})")
                print(f"  F   = {f_v:.4f}  (R21={R21:.4f})")
                results_astar3[wins_label] = {'alpha': astar, 'G3': g3_v}
            except Exception as e:
                print(f"  brentq błąd: {e}")
        else:
            # Spróbuj szerszego zakresu
            print(f"  Brak bracketa [2.920,2.935] — próba [2.91,2.96]...")
            f910 = G3_minus_R31(2.910, wins_label)
            f960 = G3_minus_R31(2.960, wins_label)
            if not np.isnan(f910) and not np.isnan(f960) and f910*f960 < 0:
                astar = brentq(lambda a: G3_minus_R31(a, wins_label),
                               2.910, 2.960, xtol=1e-8, maxiter=100)
                res_v = compute_G3_compare(astar, R_MAX)
                g3_v  = res_v[wins_label]['G3'] if res_v else float('nan')
                print(f"  α*₃ ({wins_label}) = {astar:.9f}")
                print(f"  G3  = {g3_v:.4f}  (R31={R31:.4f}, diff={g3_v-R31:+.4f})")
                results_astar3[wins_label] = {'alpha': astar, 'G3': g3_v}
            else:
                print(f"  Brak bracketa! G3_{wins_label}∈[{f910:.1f}…{f960:.1f}], R31={R31:.1f}")
                results_astar3[wins_label] = None

    print()

    # ----------------------------------------------------------------
    # FAZA 4: Analiza S₃ z nową precyzją
    # ----------------------------------------------------------------
    print("=" * 72)
    print("FAZA 4: Analiza S₃ = α*₁ + α*₂ + α*₃ z różnymi oknami")
    print("=" * 72)
    print(f"  α*₁ = {ALPHA_STAR1:.8f}")
    print(f"  α*₂ = {ALPHA_STAR2:.8f}")
    print(f"  S₂  = {ALPHA_STAR1+ALPHA_STAR2:.8f}")
    print()
    print(f"  {'Okna':>4} | {'α*₃':>12} | {'S₃':>12} | {'8-S₃':>12} | {'ppm vs 8':>10}")
    print(f"  {'-'*58}")

    # ex93/ex94 baseline
    for label_b, a3_b in [('ex93', ALPHA_STAR3_EX93), ('ex94', ALPHA_STAR3_EX94)]:
        S3 = ALPHA_STAR1 + ALPHA_STAR2 + a3_b
        diff8 = 8 - S3
        ppm8  = diff8/8*1e6
        print(f"  {label_b:>4} | {a3_b:>12.9f} | {S3:>12.9f} | {diff8:>+12.9f} | {ppm8:>+10.1f}")

    for wins_l, res_a3 in results_astar3.items():
        if res_a3 is None:
            print(f"  {wins_l:>4} | {'BRAK':>12} | {'—':>12} | {'—':>12} | {'—':>10}")
            continue
        a3 = res_a3['alpha']
        S3 = ALPHA_STAR1 + ALPHA_STAR2 + a3
        diff8 = 8 - S3
        ppm8  = diff8/8*1e6
        marker = " ←" if abs(ppm8) < 500 else ""
        print(f"  {wins_l:>4} | {a3:>12.9f} | {S3:>12.9f} | {diff8:>+12.9f} | {ppm8:>+10.1f}{marker}")

    print()

    # Ocena konsystencji między oknami
    valid_a3 = [v['alpha'] for v in results_astar3.values() if v is not None]
    if len(valid_a3) >= 2:
        spread = max(valid_a3) - min(valid_a3)
        sigma  = np.std(valid_a3)
        mean   = np.mean(valid_a3)
        print(f"  Konsystencja α*₃ między zestawami okien:")
        print(f"    Spread = {spread:.2e},  σ = {sigma:.2e},  mean = {mean:.9f}")
        S3_mean = ALPHA_STAR1 + ALPHA_STAR2 + mean
        diff_mean = 8 - S3_mean
        print(f"    S₃(mean) = {S3_mean:.8f},  8-S₃ = {diff_mean:+.8f}  ({diff_mean/8*1e6:+.1f} ppm)")
        if spread < 0.001:
            print(f"  ✅ Okna konsystentne (spread < 0.001)")
        elif spread < 0.01:
            print(f"  ⚠️  Umiarkowana konsystencja (spread = {spread:.4f})")
        else:
            print(f"  ❌ Niekonsystentne (spread = {spread:.4f}) — szum G₃ dominuje")

    print()

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print("=" * 72)
    print("PODSUMOWANIE ex95 — O-L6")
    print("=" * 72)
    print(f"\nCzas całkowity: {time.time()-t0:.1f} s")

    # Wniosek o S₃
    print(f"\n  WNIOSEK S₃=8:")
    if valid_a3:
        all_S3 = [ALPHA_STAR1+ALPHA_STAR2+a for a in valid_a3+[ALPHA_STAR3_EX93]]
        all_ppm = [(8-s)/8*1e6 for s in all_S3]
        max_ppm = max(abs(p) for p in all_ppm)
        min_ppm = min(abs(p) for p in all_ppm)
        print(f"    Zakres ppm od 8: [{min(all_ppm):+.1f}, {max(all_ppm):+.1f}]")
        if max_ppm < 500:
            print(f"  ✅ S₃=8 POTWIERDZONE (wszystkie okna <500 ppm od 8)")
        elif max_ppm < 2000:
            print(f"  ⚠️  S₃≈8 PRAWDOPODOBNE (spread < 2000 ppm)")
        else:
            print(f"  ❌ S₃=8 WĄTPLIWE (max odchylenie = {max_ppm:.0f} ppm)")
    print("=" * 72)


if __name__ == '__main__':
    mp.freeze_support()
    main()
