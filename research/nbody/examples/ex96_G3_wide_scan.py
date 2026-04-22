#!/usr/bin/env python3
"""
ex96_G3_wide_scan.py — Skan G₃(α) z oknami FAR [60,136] na α∈[2.5, 5.0]
===========================================================================
STATUS: LEGACY-TRANSLATIONAL

This script continues the older tau-sector `G3(alpha)` search chain in legacy
variables. It is useful as exploratory history, but it is not part of the
canonical synchronized `nbody` layer.

Cel:
  Po wykryciu artefaktu w ex95 (okno [20,36] w near-field dla φ²·z₀≈3.31):
  szukamy PRAWDZIWEGO zera G₃=R₃₁ na szerokim zakresie α∈[2.5, 5.0]
  używając WYŁĄCZNIE okien dalekich [60,136], które są w pełni asymptotyczne.

Metoda:
  FAZA 1: Coarse skan α∈[2.5, 5.0], krok=0.05, FAR windows only
            → identyfikacja przedziałów ze zmianą znaku G₃−R₃₁
  FAZA 2: Refinement brentq w każdym znalezionym przedziale
  FAZA 3: Dense skan wokół zera/zer (krok=0.005)
  FAZA 4: Analiza F(α) z FAR oknami — weryfikacja α*₁, α*₂

Okna FAR (wyłączne):
  [(60,76),(70,86),(80,96),(90,106),(100,116),(110,126),(120,136)]
  R_MAX = 150 (okna sięgają do r=136)

Autor: Claudian (sesja v38, O-L5)
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

R21 = M_MU_MEV / M_E_MEV    # ≈ 206.767
R31 = M_TAU_MEV / M_E_MEV   # ≈ 3477.221
R32 = M_TAU_MEV / M_MU_MEV  # ≈ 16.817

ALPHA_STAR1 = 2.43183767
ALPHA_STAR2 = 2.63557742

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0
R_MAX   = 150

# WYŁĄCZNIE okna dalekie — wyeliminowanie near-field bias
FIT_WINS_FAR = [(60, 76), (70, 86), (80, 96), (90, 106),
                (100, 116), (110, 126), (120, 136)]

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

def fit_amplitude_far(r, g):
    """Ekstrapolacja A_∞ z okien FAR [60,136]."""
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
    if not Av: return 0.0, 0.0, []
    if len(Av) == 1: return Av[0], 0.0, list(zip(rv, Av))
    rv = np.array(rv); Av_arr = np.array(Av)
    # Sprawdź rozrzut — jeśli plateau (std<0.5%), użyj mediany
    cv = np.std(Av_arr) / (np.mean(Av_arr) + 1e-12)
    if cv < 0.005:
        # Plateau: wartości są praktycznie stałe, wróć mediany
        return float(np.median(Av_arr)), float(np.std(Av_arr)), list(zip(rv.tolist(), Av_arr.tolist()))
    # Ekstrapolacja z modelem A_∞(1 + a/r + b/r²)
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
            return float(np.median(Av_arr)), float(np.std(Av_arr)), list(zip(rv.tolist(), Av_arr.tolist()))

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

def compute_one(alpha, r_max=R_MAX):
    """Oblicz F i G₃ z oknami FAR dla danego α."""
    z0 = find_z0(alpha, r_max)
    if z0 is None:
        return None
    # Integracje trzech solitonów
    r0, g0_sol = integrate_soliton(z0,       alpha, r_max)
    r1, g1_sol = integrate_soliton(PHI*z0,   alpha, r_max)
    r2, g2_sol = integrate_soliton(PHI2*z0,  alpha, r_max)

    Ae,   dAe,   _ = fit_amplitude_far(r0, g0_sol)
    Amu,  dAmu,  _ = fit_amplitude_far(r1, g1_sol)
    Atau, dAtau, _ = fit_amplitude_far(r2, g2_sol)

    if Ae < 1e-8:
        return None

    F  = (Amu  / Ae) ** 4
    G3 = (Atau / Ae) ** 4

    return {
        'alpha': alpha, 'z0': z0,
        'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
        'dAe': dAe, 'dAmu': dAmu, 'dAtau': dAtau,
        'F': F, 'G3': G3,
        'G3_minus_R31': G3 - R31,
        'F_minus_R21': F - R21,
    }

def worker(args):
    alpha, r_max = args
    try:
        return compute_one(alpha, r_max)
    except Exception as e:
        return {'alpha': alpha, 'error': str(e)}

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6

    print("=" * 72)
    print("ex96_G3_wide_scan.py — Prawdziwe α*₃ z oknami FAR (O-L5, v38)")
    print("=" * 72)
    print(f"  R_MAX     = {R_MAX}")
    print(f"  R21={R21:.4f}, R31={R31:.4f}, R32={R32:.4f}")
    print(f"  PHI={PHI:.6f}, PHI²={PHI2:.6f}")
    print(f"  FIT_WINS_FAR = {FIT_WINS_FAR}")
    print(f"  Okna FAR: tylko [60,136] — eliminacja near-field bias")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Coarse skan α∈[2.5, 5.0], krok=0.05
    # ----------------------------------------------------------------
    print("=" * 72)
    print("FAZA 1: Coarse skan α∈[2.5, 5.0], krok=0.05 (FAR okna)")
    print("=" * 72)

    alphas_coarse = np.arange(2.50, 5.01, 0.05)
    args_list = [(a, R_MAX) for a in alphas_coarse]

    print(f"  {len(alphas_coarse)} punktów, {ncpu} rdzeni...")
    with mp.Pool(ncpu) as pool:
        coarse_results = pool.map(worker, args_list)

    good_coarse = [r for r in coarse_results if r and 'error' not in r]
    print(f"  Wyniki: {len(good_coarse)}/{len(alphas_coarse)}")
    print()

    print(f"  {'α':>7} | {'z₀':>8} | {'F_far':>9} | {'G3_far':>10} | {'G3-R31':>11} | {'F-R21':>9}")
    print("  " + "-" * 72)

    # Szukaj zmiany znaku G₃−R₃₁ i F−R₂₁
    sign_changes_G3 = []
    sign_changes_F  = []
    prev_G3_sign = None
    prev_F_sign  = None
    prev_alpha   = None

    for r in sorted(good_coarse, key=lambda x: x['alpha']):
        a = r['alpha']
        g3diff = r['G3_minus_R31']
        fdiff  = r['F_minus_R21']
        print(f"  {a:7.3f} | {r['z0']:8.5f} | {r['F']:9.4f} | {r['G3']:10.4f} | "
              f"{g3diff:+11.4f} | {fdiff:+9.4f}")

        cur_G3_sign = np.sign(g3diff)
        cur_F_sign  = np.sign(fdiff)
        if prev_G3_sign is not None and cur_G3_sign != prev_G3_sign and cur_G3_sign != 0:
            sign_changes_G3.append((prev_alpha, a))
        if prev_F_sign is not None and cur_F_sign != prev_F_sign and cur_F_sign != 0:
            sign_changes_F.append((prev_alpha, a))

        prev_G3_sign = cur_G3_sign
        prev_F_sign  = cur_F_sign
        prev_alpha   = a

    print()
    print(f"  Zmiany znaku G₃−R₃₁: {sign_changes_G3}")
    print(f"  Zmiany znaku F−R₂₁:  {sign_changes_F}")

    # ----------------------------------------------------------------
    # FAZA 2: Refinement brentq dla G₃ = R₃₁
    # ----------------------------------------------------------------
    print()
    print("=" * 72)
    print("FAZA 2: brentq α*₃ dla G₃=R₃₁ (FAR okna)")
    print("=" * 72)

    alpha_star3_candidates = []

    if not sign_changes_G3:
        print("  ⚠ Brak zmiany znaku G₃−R₃₁ na [2.5, 5.0] przy oknach FAR!")
        print("  → Prawdziwe α*₃ może leżeć poza skanowanym zakresem.")
        # Sprawdź skrajne wartości
        sorted_g = sorted(good_coarse, key=lambda x: x['alpha'])
        if sorted_g:
            g3_vals = [(r['alpha'], r['G3_minus_R31']) for r in sorted_g]
            print(f"  G₃−R₃₁ na krańcach: α={g3_vals[0][0]:.2f}: {g3_vals[0][1]:+.1f}, "
                  f"α={g3_vals[-1][0]:.2f}: {g3_vals[-1][1]:+.1f}")
            # Minimum G₃
            min_r = min(sorted_g, key=lambda x: x['G3'])
            print(f"  Min G₃_far = {min_r['G3']:.4f} przy α={min_r['alpha']:.3f} "
                  f"(G₃−R₃₁={min_r['G3_minus_R31']:+.4f})")
    else:
        for a_lo, a_hi in sign_changes_G3:
            print(f"\n  Bracket G₃=R₃₁: [{a_lo:.3f}, {a_hi:.3f}]")
            # Zagęszczony skan w brackecie
            alphas_fine = np.linspace(a_lo, a_hi, 11)
            args_fine = [(a, R_MAX) for a in alphas_fine]
            with mp.Pool(ncpu) as pool:
                fine_results = pool.map(worker, args_fine)
            good_fine = {r['alpha']: r for r in fine_results if r and 'error' not in r}

            # Próba brentq
            def g3_func(a):
                if a in good_fine:
                    return good_fine[a]['G3_minus_R31']
                res = compute_one(a)
                return res['G3_minus_R31'] if res else float('nan')

            # Sprawdź czy bracket jest spójny
            va_lo = compute_one(a_lo)
            va_hi = compute_one(a_hi)
            if va_lo and va_hi:
                diff_lo = va_lo['G3_minus_R31']
                diff_hi = va_hi['G3_minus_R31']
                print(f"  G₃−R₃₁ na krańcach: [{diff_lo:+.4f}, {diff_hi:+.4f}]")
                if diff_lo * diff_hi < 0:
                    try:
                        a_star = brentq(g3_func, a_lo, a_hi, xtol=1e-8, maxiter=200)
                        res = compute_one(a_star)
                        if res:
                            print(f"  ✓ α*₃ (brentq) = {a_star:.9f}")
                            print(f"    G₃={res['G3']:.4f}, R₃₁={R31:.4f}, G₃−R₃₁={res['G3_minus_R31']:+.6f}")
                            print(f"    z₀={res['z0']:.6f}, φ²·z₀={PHI2*res['z0']:.6f}")
                            S3 = ALPHA_STAR1 + ALPHA_STAR2 + a_star
                            print(f"    S₃=α*₁+α*₂+α*₃ = {S3:.8f}")
                            print(f"    S₃−8 = {S3-8:+.8f} ({(S3-8)/8*1e6:+.1f} ppm)")
                            alpha_star3_candidates.append(a_star)
                    except Exception as e:
                        print(f"  brentq failed: {e}")
                else:
                    print(f"  ⚠ Brak bracketu (ten sam znak) — prawdopodobnie szum G₃")

    # ----------------------------------------------------------------
    # FAZA 3: Dense skan w okolicach znalezionych zer (krok=0.005)
    # ----------------------------------------------------------------
    print()
    print("=" * 72)
    print("FAZA 3: Dense skan w rejonach interesujących (krok=0.005)")
    print("=" * 72)

    # Identyfikuj rejony o małym |G₃−R₃₁| lub ze zmianą znaku
    sorted_g = sorted(good_coarse, key=lambda x: x['alpha'])
    close_regions = []
    for r in sorted_g:
        if abs(r['G3_minus_R31']) < 200:  # blisko R₃₁ — interesujące
            close_regions.append(r['alpha'])

    if not close_regions:
        # Pokaż minimum G₃ jako region
        min_g3 = min(sorted_g, key=lambda x: abs(x['G3_minus_R31']))
        close_regions = [min_g3['alpha']]
        print(f"  Brak bliskich regionów — dense skan wokół minimum |G₃−R₃₁|: α={min_g3['alpha']:.3f}")

    dense_regions = []
    for a_center in close_regions:
        a_lo = max(2.40, a_center - 0.10)
        a_hi = min(5.05, a_center + 0.10)
        dense_regions.append((a_lo, a_hi))

    # Usuń duplikaty/nakładania
    dense_regions = list(set(dense_regions))
    dense_regions.sort()

    # Spłaszcz nakładające się regiony
    merged = []
    for lo, hi in dense_regions:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append([lo, hi])

    for a_lo, a_hi in merged:
        alphas_dense = np.arange(a_lo, a_hi + 0.001, 0.005)
        args_dense = [(a, R_MAX) for a in alphas_dense]
        print(f"\n  Dense skan α∈[{a_lo:.3f},{a_hi:.3f}], krok=0.005, {len(alphas_dense)} punktów...")
        with mp.Pool(ncpu) as pool:
            dense_results = pool.map(worker, args_dense)
        good_dense = [r for r in dense_results if r and 'error' not in r]
        print(f"  Wyniki: {len(good_dense)}/{len(alphas_dense)}")

        print(f"\n  {'α':>7} | {'z₀':>8} | {'F_far':>9} | {'G3_far':>10} | {'G3-R31':>11}")
        print("  " + "-" * 56)

        prev_sign = None
        for r in sorted(good_dense, key=lambda x: x['alpha']):
            cur_sign = np.sign(r['G3_minus_R31'])
            marker = " ←" if (prev_sign is not None and cur_sign != prev_sign) else ""
            print(f"  {r['alpha']:7.4f} | {r['z0']:8.5f} | {r['F']:9.4f} | "
                  f"{r['G3']:10.4f} | {r['G3_minus_R31']:+11.4f}{marker}")
            prev_sign = cur_sign

    # ----------------------------------------------------------------
    # FAZA 4: F(α) z FAR oknami — weryfikacja α*₁, α*₂
    # ----------------------------------------------------------------
    print()
    print("=" * 72)
    print("FAZA 4: F(α) z FAR oknami — weryfikacja α*₁, α*₂")
    print("=" * 72)

    if sign_changes_F:
        print(f"  Zmiany znaku F−R₂₁: {sign_changes_F}")
        for a_lo, a_hi in sign_changes_F:
            def f_func(a):
                res = compute_one(a)
                return (res['F_minus_R21'] if res else float('nan'))
            try:
                a_star_F = brentq(f_func, a_lo, a_hi, xtol=1e-8, maxiter=200)
                res = compute_one(a_star_F)
                if res:
                    print(f"  ✓ α* (F=R₂₁, FAR) = {a_star_F:.9f}  "
                          f"F={res['F']:.6f} (R₂₁={R21:.4f})")
            except Exception as e:
                print(f"  brentq F failed: {e}")
    else:
        print("  Brak zmian znaku F−R₂₁ w coarse skanie.")
        # Sprawdź okolicę znanych α*₁,₂
        print(f"  Sprawdzam F w okolicach α*₁={ALPHA_STAR1:.5f}, α*₂={ALPHA_STAR2:.5f} (z ex93):")
        for a_test in [2.43, ALPHA_STAR1, 2.44, 2.63, ALPHA_STAR2, 2.64]:
            res = compute_one(a_test)
            if res:
                print(f"    α={a_test:.5f}: F_far={res['F']:.5f}, F−R₂₁={res['F_minus_R21']:+.5f}")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    elapsed = time.time() - t0
    print()
    print("=" * 72)
    print("PODSUMOWANIE ex96 — O-L5 (prawdziwe α*₃ z FAR oknami)")
    print("=" * 72)
    print(f"Czas całkowity: {elapsed:.1f} s")
    print()
    if alpha_star3_candidates:
        for a in alpha_star3_candidates:
            S3 = ALPHA_STAR1 + ALPHA_STAR2 + a
            print(f"  ✓ α*₃ (FAR) = {a:.9f}")
            print(f"    S₃ = {S3:.9f}")
            print(f"    S₃−8 = {S3-8:+.9f} ({(S3-8)/8*1e6:+.1f} ppm)")
    else:
        print("  ⚠ Brak zera G₃=R₃₁ na α∈[2.5, 5.0] przy oknach FAR.")
        # Pokaż zakres G₃
        sorted_g = sorted(good_coarse, key=lambda x: x['alpha'])
        if sorted_g:
            g3_min = min(sorted_g, key=lambda x: x['G3'])
            g3_max = max(sorted_g, key=lambda x: x['G3'])
            print(f"  Min G₃_far = {g3_min['G3']:.2f} przy α={g3_min['alpha']:.3f}")
            print(f"  Max G₃_far = {g3_max['G3']:.2f} przy α={g3_max['alpha']:.3f}")
            print(f"  R₃₁ = {R31:.4f}")
            all_above = all(r['G3'] > R31 for r in sorted_g)
            all_below = all(r['G3'] < R31 for r in sorted_g)
            if all_above:
                print(f"  → G₃_far > R₃₁ wszędzie na [2.5, 5.0]")
                print(f"  → α*₃ nie istnieje w tym zakresie (lub metoda FAR nieodpowiednia)")
            elif all_below:
                print(f"  → G₃_far < R₃₁ wszędzie na [2.5, 5.0]")
            else:
                print(f"  → G₃_far przekracza R₃₁ gdzieś w zakresie")
    print()
    print("  Weryfikacja R₃₁ = m_τ/m_e = {R31:.4f} (NIE czwarta potęga)")
    print("=" * 72)

if __name__ == '__main__':
    mp.freeze_support()
    main()
