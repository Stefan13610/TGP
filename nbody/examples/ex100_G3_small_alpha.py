#!/usr/bin/env python3
"""
ex100_G3_small_alpha.py — G₃_far(α) na α∈[1.5, 2.45] (S1, v39)
=================================================================
Pytanie: czy G₃_far maleje poniżej R₃₁ przy małych α < 2.4?

Z ex96: G₃_far(2.40)=4534, G₃_far(2.50)=4667 — minimum gdzieś koło α≈2.4.
Ale skanu α<2.4 nie było. Sprawdzamy teraz.

Oczekiwania:
  a) G₃_far rośnie dla α<2 (brak zera)
  b) G₃_far maleje i przekracza R₃₁ od góry dla α<2.4 (ZERO FIZYCZNE!)
  c) z₀ przestaje istnieć dla małych α → G₃ nieokreślone

Dodatkowe: α_min_far (minimum F_far) z FAR oknami — FAZA 2.

Autor: Claudian (sesja v39, S1)
Data: 2026-03-29
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit, minimize_scalar
import multiprocessing as mp
import time
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0
PHI2 = PHI * PHI
PHI3 = PHI2 * PHI

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
    for _ in range(25):
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
    if len(Av) < 3 or np.std(Av)/(np.mean(Av)+1e-12) < 0.005:
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

def find_z0(alpha, r_max=R_MAX, lo=1.00, hi=2.5, n=80):
    """Szerszy przedział dla małych α."""
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

def compute_all(alpha, r_max=R_MAX):
    """Oblicz F, G₃, G₃_phi3 z FAR oknami."""
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0s = integrate_soliton(z0,      alpha, r_max)
    r1, g1s = integrate_soliton(PHI*z0,  alpha, r_max)
    r2, g2s = integrate_soliton(PHI2*z0, alpha, r_max)

    Ae   = amplitude_far(r0, g0s)
    Amu  = amplitude_far(r1, g1s)
    Atau = amplitude_far(r2, g2s)
    if Ae < 1e-8: return None

    # φ³·z₀ soliton (S4 preview)
    try:
        r3, g3s = integrate_soliton(PHI3*z0, alpha, r_max)
        Atau3 = amplitude_far(r3, g3s)
        G3_phi3 = (Atau3/Ae)**4 if Ae > 1e-8 else float('nan')
    except:
        Atau3 = 0.0; G3_phi3 = float('nan')

    return {
        'alpha': alpha, 'z0': z0,
        'Ae': Ae, 'Amu': Amu, 'Atau': Atau, 'Atau3': Atau3,
        'F':  (Amu/Ae)**4,
        'G3': (Atau/Ae)**4,
        'G3_phi3': G3_phi3,
        'phi2_z0': PHI2*z0,
        'phi3_z0': PHI3*z0,
    }

def worker(args):
    alpha, r_max = args
    try: return compute_all(alpha, r_max)
    except: return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6

    print("=" * 70)
    print("ex100_G3_small_alpha.py — G₃_far(α<2.4)? (S1, v39)")
    print("=" * 70)
    print(f"  R21={R21:.6f}, R31={R31:.6f}, R32={R32:.6f}")
    print(f"  PHI={PHI:.8f}, PHI²={PHI2:.8f}, PHI³={PHI3:.8f}")
    print(f"  R_MAX={R_MAX}")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Skan α∈[1.50, 2.55], krok=0.05
    # ----------------------------------------------------------------
    print("=" * 70)
    print("FAZA 1: Skan G₃_far(α) na α∈[1.50, 2.55], krok=0.05")
    print("=" * 70)

    alphas_1 = np.arange(1.50, 2.56, 0.05)
    args_1 = [(a, R_MAX) for a in alphas_1]
    print(f"  {len(alphas_1)} punktów, {ncpu} rdzeni...")
    with mp.Pool(ncpu) as pool:
        results_1 = pool.map(worker, args_1)

    good_1 = [r for r in results_1 if r]
    print(f"  Wyniki: {len(good_1)}/{len(alphas_1)}")
    print()

    print(f"  {'α':>6} | {'z₀':>8} | {'φ²·z₀':>7} | {'F_far':>9} | "
          f"{'G3_far':>10} | {'G3-R31':>11} | {'G3_φ³':>10}")
    print(f"  {'-'*80}")

    sign_G3 = []; prev_sign = None; prev_a = None
    sign_F  = []; prev_sign_F = None
    for r in sorted(good_1, key=lambda x: x['alpha']):
        a = r['alpha']
        g3d = r['G3'] - R31
        fd  = r['F'] - R21
        cur_s  = np.sign(g3d)
        cur_sf = np.sign(fd)
        mk = " ←" if (prev_sign and cur_s != prev_sign) else ""
        if prev_sign  and cur_s  != prev_sign:  sign_G3.append((prev_a, a))
        if prev_sign_F and cur_sf != prev_sign_F: sign_F.append((prev_a, a))
        g3phi3_str = f"{r['G3_phi3']:10.2f}" if not np.isnan(r['G3_phi3']) else "      N/A"
        print(f"  {a:6.3f} | {r['z0']:8.5f} | {r['phi2_z0']:7.4f} | {r['F']:9.4f} | "
              f"{r['G3']:10.4f} | {g3d:+11.4f} | {g3phi3_str}{mk}")
        prev_sign = cur_s; prev_sign_F = cur_sf; prev_a = a

    print()
    print(f"  Zmiany znaku G₃−R₃₁: {sign_G3}")
    print(f"  Zmiany znaku F−R₂₁:  {sign_F}")

    # ----------------------------------------------------------------
    # FAZA 2: brentq jeśli jest zmiana znaku G₃
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 2: brentq α*₃_new (jeśli zmiana znaku G₃−R₃₁)")
    print("=" * 70)

    alpha_star3_new = None
    if sign_G3:
        for a_lo, a_hi in sign_G3:
            def g3_func(alpha):
                res = compute_all(alpha)
                return (res['G3'] - R31) if res else float('nan')
            try:
                a_star = brentq(g3_func, a_lo, a_hi, xtol=1e-9, maxiter=200)
                res = compute_all(a_star)
                if res:
                    alpha_star3_new = a_star
                    print(f"  ✓ α*₃_new = {a_star:.9f}")
                    print(f"    z₀={res['z0']:.7f}, φ²·z₀={PHI2*res['z0']:.7f}")
                    print(f"    G₃_far={res['G3']:.6f} (R₃₁={R31:.6f})")
                    print(f"    F_far={res['F']:.6f} (R₂₁={R21:.6f})")
                    print(f"    G₃_φ³_far={res['G3_phi3']:.4f}")
                    print()
                    S3 = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + a_star
                    print(f"    S₃ = α*₁_far+α*₂_far+α*₃_new = {S3:.9f}")
                    print(f"    S₃ − 8 = {S3-8:+.9f}  ({(S3-8)/8*1e6:+.1f} ppm)")
            except Exception as e:
                print(f"  brentq failed [{a_lo:.3f},{a_hi:.3f}]: {e}")
    else:
        print("  Brak zmiany znaku G₃−R₃₁ na [1.50, 2.55].")
        sorted_r = sorted(good_1, key=lambda x: x['G3'])
        if sorted_r:
            mn = sorted_r[0]
            print(f"  Min G₃_far = {mn['G3']:.4f} przy α={mn['alpha']:.3f} "
                  f"(G₃−R₃₁={mn['G3']-R31:+.4f})")
            # Trend przy małych α
            small = [r for r in good_1 if r['alpha'] <= 1.8]
            if small:
                print(f"\n  Profil dla małych α (trend G₃_far):")
                for r in sorted(small, key=lambda x: x['alpha']):
                    print(f"    α={r['alpha']:.3f}: G₃_far={r['G3']:.2f}, F_far={r['F']:.4f}, z₀={r['z0']:.5f}")

    # ----------------------------------------------------------------
    # FAZA 3: α_min_far — minimum F_far (S5)
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 3: α_min_far — minimum F_far z FAR oknami")
    print("=" * 70)

    # Szukamy minimum F w zakresie [2.0, 3.0] z dobrze rozwiązanymi α
    good_F = [(r['alpha'], r['F']) for r in good_1 if r and r['F'] > 0]
    # Uzupełnij z ex96/ex98: α∈[2.50, 3.00]
    alphas_F2 = np.arange(2.55, 3.05, 0.05)
    args_F2 = [(a, R_MAX) for a in alphas_F2]
    with mp.Pool(ncpu) as pool:
        res_F2 = pool.map(worker, args_F2)
    good_F += [(r['alpha'], r['F']) for r in res_F2 if r and r['F'] > 0]
    good_F.sort(key=lambda x: x[0])

    # Znajdź minimum
    amin_idx = min(range(len(good_F)), key=lambda i: good_F[i][1])
    a_min_approx, F_min_approx = good_F[amin_idx]
    print(f"  Przybliżone minimum F_far: α≈{a_min_approx:.3f}, F_far≈{F_min_approx:.4f}")
    print(f"  R₂₁ = {R21:.4f}")
    print(f"  F_min_far vs R₂₁: {F_min_approx - R21:+.4f}")

    # Dense skan wokół minimum
    a_lo_min = max(1.50, a_min_approx - 0.3)
    a_hi_min = min(3.10, a_min_approx + 0.3)
    alphas_min = np.arange(a_lo_min, a_hi_min+0.001, 0.01)
    args_min = [(a, R_MAX) for a in alphas_min]
    print(f"\n  Dense skan F_far na [{a_lo_min:.2f},{a_hi_min:.2f}], krok=0.01...")
    with mp.Pool(ncpu) as pool:
        res_min = pool.map(worker, args_min)
    good_min = [(r['alpha'], r['F']) for r in res_min if r and r['F'] > 0]

    amin2_idx = min(range(len(good_min)), key=lambda i: good_min[i][1])
    a_min2, F_min2 = good_min[amin2_idx]
    print(f"  Minimum F_far (krok=0.01): α={a_min2:.4f}, F_far={F_min2:.6f}")

    # Brentq dF/dα = 0 (zero pochodnej)
    def F_func(alpha):
        res = compute_all(alpha)
        return res['F'] if res else float('nan')

    # Numeryczna pochodna
    eps = 0.005
    def dF_dalpha(alpha):
        fp = F_func(alpha + eps)
        fm = F_func(alpha - eps)
        if np.isnan(fp) or np.isnan(fm): return float('nan')
        return (fp - fm) / (2*eps)

    # Bracket dla minimum
    a_bracket_lo = max(a_min2 - 0.15, 1.5)
    a_bracket_hi = min(a_min2 + 0.15, 3.0)
    d_lo = dF_dalpha(a_bracket_lo)
    d_hi = dF_dalpha(a_bracket_hi)
    print(f"\n  dF/dα na [{a_bracket_lo:.3f},{a_bracket_hi:.3f}]:")
    print(f"    dF/dα({a_bracket_lo:.3f}) = {d_lo:+.6f}")
    print(f"    dF/dα({a_bracket_hi:.3f}) = {d_hi:+.6f}")

    if d_lo * d_hi < 0:
        try:
            a_min_precise = brentq(dF_dalpha, a_bracket_lo, a_bracket_hi,
                                   xtol=1e-6, maxiter=200)
            res_min_p = compute_all(a_min_precise)
            print(f"\n  ✓ α_min_far (brentq) = {a_min_precise:.7f}")
            if res_min_p:
                print(f"    F_min_far = {res_min_p['F']:.7f}")
                print(f"    F_min_far − R₂₁ = {res_min_p['F']-R21:+.7f}")
                print(f"    z₀={res_min_p['z0']:.7f}")
                print(f"    (α*₁_far+α*₂_far)/2 = {(ALPHA_STAR1_FAR+ALPHA_STAR2_FAR)/2:.7f}")
                print(f"    α_min_far − (α*₁+α*₂)/2 = "
                      f"{a_min_precise-(ALPHA_STAR1_FAR+ALPHA_STAR2_FAR)/2:+.7f}")
        except Exception as e:
            print(f"  brentq dF/dα failed: {e}")
    else:
        print(f"  Brak zmiany znaku dF/dα w [{a_bracket_lo:.3f},{a_bracket_hi:.3f}].")

    # ----------------------------------------------------------------
    # FAZA 4: G₃_φ³ profil — czy φ³-mnożnik daje R₃₁?
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 4: G₃_φ³=(A(φ³·z₀)/A(z₀))^4 vs R₃₁ (preview S4)")
    print("=" * 70)

    # Użyj wyników z FAZY 1
    phi3_results = [(r['alpha'], r['G3_phi3']) for r in good_1
                    if r and not np.isnan(r['G3_phi3'])]
    if phi3_results:
        print(f"\n  {'α':>6} | {'G3_φ³':>12} | {'G3_φ³−R₃₁':>12} | {'φ³·z₀':>8}")
        print(f"  {'-'*50}")
        sign_phi3 = []; prev_s3=None; prev_a3=None
        for r in sorted(good_1, key=lambda x: x['alpha']):
            if np.isnan(r['G3_phi3']): continue
            diff = r['G3_phi3'] - R31
            cur_s3 = np.sign(diff)
            mk = " ←" if (prev_s3 and cur_s3 != prev_s3) else ""
            if prev_s3 and cur_s3 != prev_s3: sign_phi3.append((prev_a3, r['alpha']))
            print(f"  {r['alpha']:6.3f} | {r['G3_phi3']:12.4f} | {diff:+12.4f} | {r['phi3_z0']:8.5f}{mk}")
            prev_s3 = cur_s3; prev_a3 = r['alpha']

        print()
        print(f"  Zmiany znaku G₃_φ³−R₃₁: {sign_phi3}")
        if sign_phi3:
            for a_lo3, a_hi3 in sign_phi3:
                def g3phi3_func(alpha):
                    res = compute_all(alpha)
                    return (res['G3_phi3'] - R31) if res and not np.isnan(res['G3_phi3']) else float('nan')
                try:
                    a_s3 = brentq(g3phi3_func, a_lo3, a_hi3, xtol=1e-8, maxiter=200)
                    res3 = compute_all(a_s3)
                    print(f"  ✓ α*₃_φ³ = {a_s3:.9f}")
                    if res3:
                        print(f"    G₃_φ³={res3['G3_phi3']:.6f}, F_far={res3['F']:.6f}")
                        S3 = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + a_s3
                        print(f"    S₃_φ³ = {S3:.9f}  (S₃−8={S3-8:+.6f}, {(S3-8)/8*1e6:+.1f} ppm)")
                except Exception as e:
                    print(f"  brentq φ³ failed: {e}")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("PODSUMOWANIE ex100 — S1: G₃_far na małych α")
    print("=" * 70)
    print(f"Czas: {time.time()-t0:.1f} s")
    print()

    # Trend G₃ dla małych α
    small_a = [(r['alpha'], r['G3']) for r in good_1 if r and r['alpha'] <= 2.0]
    if small_a:
        print(f"  Trend G₃_far dla α≤2.0:")
        for a, g in sorted(small_a):
            print(f"    α={a:.2f}: G₃_far={g:.2f} (G₃−R₃₁={g-R31:+.2f})")
    print()
    if alpha_star3_new:
        print(f"  ✅ ZNALEZIONO: α*₃_new = {alpha_star3_new:.7f}")
        S3 = ALPHA_STAR1_FAR + ALPHA_STAR2_FAR + alpha_star3_new
        print(f"     S₃ = {S3:.7f}  (S₃−8={S3-8:+.6f})")
    else:
        mn = min(good_1, key=lambda r: r['G3'])
        print(f"  ❌ Brak zera G₃_far=R₃₁ na α∈[1.5,2.55].")
        print(f"     Min G₃_far={mn['G3']:.2f} przy α={mn['alpha']:.3f} "
              f"({mn['G3']-R31:+.2f} nad R₃₁)")
    print("=" * 70)

if __name__ == '__main__':
    mp.freeze_support()
    main()
