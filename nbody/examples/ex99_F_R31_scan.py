#!/usr/bin/env python3
"""
ex99_F_R31_scan.py — Czy F_far(α) = R₃₁ ma rozwiązanie? (O-L8, v38)
======================================================================
STATUS: LEGACY-TRANSLATIONAL

This file continues the older `F(alpha)` / tau-target exploration in legacy
variables. It should be read as historical exploratory context, not as a
canonical current `nbody` source.

Hipoteza:
  Z ex96/ex98: F_far rośnie monotonicznie dla α>3 (233→672 na [3,5]).
  Interpolacja sugeruje F_far(α)=R₃₁=3477 gdzieś w α≈8–12.

  Warunek F_far=R₃₁ oznaczałby:
    (A_tail(φ·z₀, α) / A_tail(z₀, α))^4 = R₃₁
  — tę samą formułę co dla muona, ale z innym targetem.

  Jeśli α*₃_F ≈ 8 = 4·α_TGP = 4·2: potencjalnie czysta algebraicznie.
  Jeśli α*₃_F = coś innego: wnioski inne.

Plan:
  FAZA 1: Coarse skan F_far(α) na α∈[5, 15], krok=0.5
  FAZA 2: brentq α*₃_F gdzie F_far=R₃₁
  FAZA 3: Dense skan wokół α*₃_F, weryfikacja G₃_far(φ²·z₀) przy tej α
  FAZA 4: Analiza algebraiczna α*₁_far, α*₂_far, α*₃_F

Autor: Claudian (sesja v38, O-L8)
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

R21 = M_MU_MEV / M_E_MEV
R31 = M_TAU_MEV / M_E_MEV
R32 = M_TAU_MEV / M_MU_MEV

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0
R_MAX   = 150

FIT_WINS_FAR = [(60, 76), (70, 86), (80, 96), (90, 106),
                (100, 116), (110, 126), (120, 136)]

ALPHA_STAR1_FAR = 2.436011168
ALPHA_STAR2_FAR = 2.753824880

# ================================================================
# ODE
# ================================================================
def integrate_soliton(g0, alpha, r_max=R_MAX):
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
    for _ in range(30):
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

def amplitude_far(r, g):
    Av, rv = [], []
    for rL, rR in FIT_WINS_FAR:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2)))
        rv.append(float(rL))
    if not Av: return 0.0
    Av = np.array(Av)
    cv = np.std(Av) / (np.mean(Av) + 1e-12)
    if cv < 0.005 or len(Av) < 3:
        return float(np.median(Av))
    try:
        p, _ = curve_fit(lambda x, ai, a, b: ai*(1+a/x+b/x**2),
                         np.array(rv), Av, p0=[Av[-1],0.,0.], maxfev=2000)
        return float(p[0])
    except:
        return float(np.median(Av))

def B_coeff_func(g0, alpha, r_max=R_MAX):
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
    for gg in gs:
        try: Bs.append(B_coeff_func(gg, alpha, r_max))
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

def compute_F_G3_full(alpha, r_max=R_MAX):
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0s = integrate_soliton(z0,      alpha, r_max)
    r1, g1s = integrate_soliton(PHI*z0,  alpha, r_max)
    r2, g2s = integrate_soliton(PHI2*z0, alpha, r_max)
    Ae   = amplitude_far(r0, g0s)
    Amu  = amplitude_far(r1, g1s)
    Atau = amplitude_far(r2, g2s)
    if Ae < 1e-8: return None
    return {'alpha': alpha, 'z0': z0, 'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
            'F':  (Amu/Ae)**4,
            'G3': (Atau/Ae)**4}

def worker(args):
    alpha, r_max = args
    try: return compute_F_G3_full(alpha, r_max)
    except: return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6

    print("=" * 68)
    print("ex99_F_R31_scan.py — F_far(α) = R₃₁? (O-L8, v38)")
    print("=" * 68)
    print(f"  R21={R21:.6f}, R31={R31:.6f}, R32={R32:.6f}")
    print(f"  α*₁_far={ALPHA_STAR1_FAR:.7f}, α*₂_far={ALPHA_STAR2_FAR:.7f}")
    print(f"  S₂_far = {ALPHA_STAR1_FAR+ALPHA_STAR2_FAR:.7f}")
    print(f"  R_MAX={R_MAX}")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Coarse skan F_far(α) na α∈[5, 15], krok=0.5
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 1: Coarse skan F_far(α) na α∈[5,15], krok=0.5")
    print("=" * 68)

    alphas_c = np.arange(5.0, 15.5, 0.5)
    args_c = [(a, R_MAX) for a in alphas_c]
    print(f"  {len(alphas_c)} punktów, {ncpu} rdzeni...")
    with mp.Pool(ncpu) as pool:
        coarse = pool.map(worker, args_c)

    good_c = [r for r in coarse if r]
    print(f"  Wyniki: {len(good_c)}/{len(alphas_c)}")
    print()
    print(f"  {'α':>7} | {'z₀':>8} | {'F_far':>10} | {'F-R21':>9} | {'F-R31':>10} | {'G3_far':>10}")
    print(f"  {'-'*68}")

    sign_changes_F31 = []
    prev_sign = None; prev_alpha = None
    for r in sorted(good_c, key=lambda x: x['alpha']):
        a = r['alpha']
        f_diff21 = r['F'] - R21
        f_diff31 = r['F'] - R31
        marker = ""
        cur_sign = np.sign(f_diff31)
        if prev_sign is not None and cur_sign != prev_sign:
            sign_changes_F31.append((prev_alpha, a))
            marker = " ←"
        print(f"  {a:7.2f} | {r['z0']:8.5f} | {r['F']:10.4f} | {f_diff21:+9.4f} | "
              f"{f_diff31:+10.4f} | {r['G3']:10.2f}{marker}")
        prev_sign = cur_sign; prev_alpha = a

    # ----------------------------------------------------------------
    # FAZA 2: brentq α*₃_F = F_far = R₃₁
    # ----------------------------------------------------------------
    print()
    print("=" * 68)
    print("FAZA 2: brentq α*₃_F — F_far=R₃₁")
    print("=" * 68)

    alpha_star3_F = None

    if sign_changes_F31:
        for a_lo, a_hi in sign_changes_F31:
            def F_minus_R31(alpha):
                res = compute_F_G3_full(alpha)
                return (res['F'] - R31) if res else float('nan')
            try:
                a_star = brentq(F_minus_R31, a_lo, a_hi, xtol=1e-7, maxiter=200)
                res = compute_F_G3_full(a_star)
                alpha_star3_F = a_star
                print(f"  ✓ α*₃_F = {a_star:.9f}")
                if res:
                    print(f"    z₀={res['z0']:.7f}, φ·z₀={PHI*res['z0']:.7f}")
                    print(f"    F_far = {res['F']:.6f} (R₃₁={R31:.6f})")
                    print(f"    G₃_far(φ²·z₀) = {res['G3']:.4f}")
                    print()
                    # Sumy z α*₁,₂_far
                    print(f"  Relacje z α*₁,₂_far:")
                    S3_F = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + a_star
                    P3_F = ALPHA_STAR1_FAR * ALPHA_STAR2_FAR * a_star
                    print(f"    α*₁+α*₂+α*₃ = {S3_F:.9f}")
                    print(f"    α*₁·α*₂·α*₃ = {P3_F:.9f}")
                    print(f"    α*₃/α*₁ = {a_star/ALPHA_STAR1_FAR:.6f}")
                    print(f"    α*₃/α*₂ = {a_star/ALPHA_STAR2_FAR:.6f}")
            except Exception as e:
                print(f"  brentq failed: {e}")
    else:
        print("  Brak zmiany znaku F_far−R₃₁ na [5,15].")
        # Sprawdź czy F jest wszędzie poniżej lub powyżej
        f_vals = [(r['alpha'], r['F']) for r in good_c]
        max_f = max(f_vals, key=lambda x: x[1])
        min_f = min(f_vals, key=lambda x: x[1])
        print(f"  F_far ∈ [{min_f[1]:.2f} (α={min_f[0]:.1f}), {max_f[1]:.2f} (α={max_f[0]:.1f})]")
        print(f"  R₃₁ = {R31:.4f}")
        if max_f[1] < R31:
            print(f"  → F_far < R₃₁ wszędzie na [5,15]. Rozszerzam skan...")
            # Rozszerzony skan do α=25
            alphas_ext = np.arange(15, 25.5, 0.5)
            args_ext = [(a, R_MAX) for a in alphas_ext]
            print(f"  Rozszerzony skan α∈[15,25], krok=0.5...")
            with mp.Pool(ncpu) as pool:
                ext_results = pool.map(worker, args_ext)
            good_ext = [r for r in ext_results if r]

            print(f"  {'α':>7} | {'F_far':>10} | {'F-R31':>10} | {'z₀':>8}")
            print(f"  {'-'*45}")
            prev_sign2 = None; prev_alpha2 = None; sign2 = []
            for r in sorted(good_ext, key=lambda x: x['alpha']):
                a = r['alpha']
                fd = r['F'] - R31
                cur_s = np.sign(fd)
                marker = " ←" if (prev_sign2 and cur_s != prev_sign2) else ""
                print(f"  {a:7.2f} | {r['F']:10.4f} | {fd:+10.4f} | {r['z0']:8.5f}{marker}")
                if prev_sign2 and cur_s != prev_sign2:
                    sign2.append((prev_alpha2, a))
                prev_sign2 = cur_s; prev_alpha2 = a

            if sign2:
                for a_lo2, a_hi2 in sign2:
                    def F_minus_R31_ext(alpha):
                        res = compute_F_G3_full(alpha)
                        return (res['F'] - R31) if res else float('nan')
                    try:
                        a_star2 = brentq(F_minus_R31_ext, a_lo2, a_hi2, xtol=1e-6, maxiter=100)
                        res2 = compute_F_G3_full(a_star2)
                        alpha_star3_F = a_star2
                        print(f"\n  ✓ α*₃_F (ext) = {a_star2:.7f}")
                        if res2:
                            print(f"    F_far={res2['F']:.4f} (R₃₁={R31:.4f})")
                            S3 = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + a_star2
                            print(f"    S₃_F = α*₁+α*₂+α*₃ = {S3:.7f}")
                    except Exception as e:
                        print(f"  brentq ext failed: {e}")

    # ----------------------------------------------------------------
    # FAZA 3: Dense skan α wokół znalezionego α*₃_F ± 0.3, krok=0.05
    # ----------------------------------------------------------------
    if alpha_star3_F is not None:
        print()
        print("=" * 68)
        print(f"FAZA 3: Dense skan wokół α*₃_F={alpha_star3_F:.4f} ± 0.5, krok=0.05")
        print("=" * 68)

        a_lo3 = max(2.0, alpha_star3_F - 0.5)
        a_hi3 = alpha_star3_F + 0.5
        alphas_d = np.arange(a_lo3, a_hi3+0.01, 0.05)
        args_d = [(a, R_MAX) for a in alphas_d]
        with mp.Pool(ncpu) as pool:
            dense = pool.map(worker, args_d)
        good_d = [r for r in dense if r]

        print(f"  {'α':>8} | {'z₀':>8} | {'F_far':>10} | {'F-R31':>11} | {'G3_far':>12}")
        print(f"  {'-'*62}")
        for r in sorted(good_d, key=lambda x: x['alpha']):
            marker = " ←" if abs(r['F']-R31) < 20 else ""
            print(f"  {r['alpha']:8.4f} | {r['z0']:8.5f} | {r['F']:10.4f} | "
                  f"{r['F']-R31:+11.4f} | {r['G3']:12.4f}{marker}")

    # ----------------------------------------------------------------
    # FAZA 4: Analiza algebraiczna
    # ----------------------------------------------------------------
    print()
    print("=" * 68)
    print("FAZA 4: Analiza algebraiczna α*₁_far, α*₂_far, α*₃_F")
    print("=" * 68)

    a1 = ALPHA_STAR1_FAR
    a2 = ALPHA_STAR2_FAR
    a3 = alpha_star3_F

    if a3 is not None:
        S3 = a1+a2+a3
        P3 = a1*a2*a3
        print(f"\n  α*₁ = {a1:.9f}")
        print(f"  α*₂ = {a2:.9f}")
        print(f"  α*₃ = {a3:.9f}")
        print(f"\n  Sumy i iloczyny:")
        print(f"  S₂ = α*₁+α*₂        = {a1+a2:.9f}")
        print(f"  S₃ = α*₁+α*₂+α*₃   = {S3:.9f}")
        print(f"  P₂ = α*₁·α*₂        = {a1*a2:.9f}")
        print(f"  P₃ = α*₁·α*₂·α*₃   = {P3:.9f}")
        print(f"  α*₃/α*₁             = {a3/a1:.9f}")
        print(f"  α*₃/α*₂             = {a3/a2:.9f}")
        print(f"  α*₃/(α*₁+α*₂)      = {a3/(a1+a2):.9f}")

        # Kandydaci dla S₃
        candidates_S3 = [
            ("4·α_TGP = 4·2 = 8",     8.0),
            ("3φ²",                    3*PHI2),
            ("5φ",                     5*PHI),
            ("2e+R21^(1/4)",           2*np.e + R21**0.25),
            ("3√6",                    3*np.sqrt(6)),
            ("e+φ³+1",                 np.e + PHI**3 + 1),
            ("π²",                     np.pi**2),
            ("4π−R32^(1/3)",           4*np.pi - R32**(1/3)),
        ]
        print(f"\n  Kandydaci S₃ = {S3:.6f}:")
        print(f"  {'Kandydat':30s} | {'Wartosc':>10} | {'delta':>10} | {'ppm':>8}")
        print(f"  {'-'*62}")
        for name, val in sorted(candidates_S3, key=lambda x: abs(x[1]-S3)):
            d = S3 - val
            ppm = d/S3*1e6
            print(f"  {name:30s} | {val:10.6f} | {d:+10.6f} | {ppm:+8.1f}")

        # Kandydaci dla α*₃ samego
        candidates_a3 = [
            ("4·α_TGP−S₂ = 8−5.19",   8.0-(a1+a2)),
            ("φ⁴",                     PHI**4),
            ("e²",                     np.e**2),
            ("3φ",                     3*PHI),
            ("2π",                     2*np.pi),
            ("5",                      5.0),
            ("2φ²",                    2*PHI2),
            ("4φ−1",                   4*PHI-1),
            ("φ²+e",                   PHI2+np.e),
            ("R32",                    R32),
        ]
        print(f"\n  Kandydaci α*₃_F = {a3:.6f}:")
        print(f"  {'Kandydat':30s} | {'Wartosc':>10} | {'delta':>10} | {'ppm':>8}")
        print(f"  {'-'*62}")
        for name, val in sorted(candidates_a3, key=lambda x: abs(x[1]-a3)):
            d = a3 - val
            ppm = d/a3*1e6
            print(f"  {name:30s} | {val:10.6f} | {d:+10.6f} | {ppm:+8.1f}")
    else:
        print("  α*₃_F nie znaleziona.")
        print(f"  α*₁_far = {a1:.9f}")
        print(f"  α*₂_far = {a2:.9f}")

    # Podsumowanie
    print()
    print("=" * 68)
    print("PODSUMOWANIE ex99 — F_far = R₃₁?")
    print("=" * 68)
    print(f"Czas: {time.time()-t0:.1f} s")
    if alpha_star3_F is not None:
        print(f"  α*₃_F (F_far=R₃₁) = {alpha_star3_F:.7f}")
        S3 = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + alpha_star3_F
        print(f"  S₃_F = {S3:.7f}  (vs 8.000000)")
        print(f"  S₃_F − 8 = {S3-8:+.7f} ({(S3-8)/8*1e6:+.1f} ppm)")
    else:
        print("  F_far = R₃₁ nie znaleziona w α∈[5,25].")
    print("=" * 68)

if __name__ == '__main__':
    mp.freeze_support()
    main()
