#!/usr/bin/env python3
"""
ex102_hierarchy.py — Inne hierarchie mnożnikowe + α_min_far precyzja (S4+S5, v39)
==================================================================================
S4: Czy inny mnożnik m (zamiast φ²) może dać G₃=R₃₁ w fizycznym zakresie α?
    Test: m = φ, 1.5, 2.0, e, √5, √3, 3/2 przy α∈{α*₁_far, α*₂_far, 2.5, 3.0}
    Nowy twist: używamy TYLKO okien [80,140] (wyeliminowanie [60,76] biasu)

S5: Precyzyjne α_min_far z brentq dF_far/dα=0
    Z ex100: α_min_far≈2.5595, F_min_far=196.85 (< R₂₁=206.77)

Autor: Claudian (sesja v39, S4+S5)
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

# Okna VFAR — bez [60,76] który jest w tranzycie dla φ²·z₀≈3.1–3.3
FIT_WINS_VFAR = [(80,96),(90,106),(100,116),(110,126),(120,136),(130,146),(140,156)]
FIT_WINS_FAR  = [(60,76),(70,86),(80,96),(90,106),(100,116),(110,126),(120,136)]
FIT_WINS_STD  = [(20,36),(30,46),(40,56),(50,66)]

ALPHA_STAR1_FAR = 2.436011168
ALPHA_STAR2_FAR = 2.753824880

# Mnożniki do testu
MULTIPLIERS = [
    ("φ",    PHI),
    ("φ²",   PHI2),
    ("φ³",   PHI3),
    ("√2",   np.sqrt(2)),
    ("√3",   np.sqrt(3)),
    ("√5",   np.sqrt(5)),
    ("e",    np.e),
    ("2",    2.0),
    ("2.5",  2.5),
    ("3",    3.0),
    ("1+φ",  1+PHI),
]

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

def amplitude_from_wins(r, g, wins):
    Av, rv = [], []
    for rL, rR in wins:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2)))
        rv.append(float(rL))
    if not Av: return 0.0
    Av_arr = np.array(Av)
    if len(Av_arr) < 3 or np.std(Av_arr)/(np.mean(Av_arr)+1e-12) < 0.005:
        return float(np.median(Av_arr))
    try:
        p, _ = curve_fit(lambda x,ai,a,b: ai*(1+a/x+b/x**2),
                         np.array(rv), Av_arr, p0=[Av_arr[-1],0.,0.], maxfev=2000)
        return float(p[0])
    except:
        return float(np.median(Av_arr))

def B_coeff_func(g0, alpha, r_max=R_MAX):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX, lo=1.00, hi=2.5, n=80):
    gs = np.linspace(lo, hi, n)
    Bs = []
    for gg in gs:
        try: Bs.append(B_coeff_func(gg, alpha, r_max))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs)-1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i]*Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff_func(x, alpha, r_max),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def amplitude_vfar(r, g):
    return amplitude_from_wins(r, g, FIT_WINS_VFAR)

def amplitude_far(r, g):
    return amplitude_from_wins(r, g, FIT_WINS_FAR)

def compute_F_vfar(alpha, r_max=R_MAX):
    """F z oknami VFAR [80,156] — wyeliminowanie biasu [60,76]."""
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0s = integrate_soliton(z0,     alpha, r_max)
    r1, g1s = integrate_soliton(PHI*z0, alpha, r_max)
    Ae  = amplitude_vfar(r0, g0s)
    Amu = amplitude_vfar(r1, g1s)
    if Ae < 1e-8: return None
    F = (Amu/Ae)**4
    return {'alpha': alpha, 'z0': z0, 'Ae': Ae, 'Amu': Amu, 'F': F}

def compute_G3_multiplier(alpha, multiplier, r_max=R_MAX):
    """G₃ z mnożnikiem m·z₀, okna VFAR."""
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0s = integrate_soliton(z0,             alpha, r_max)
    rm, gms = integrate_soliton(multiplier*z0,  alpha, r_max)
    Ae = amplitude_vfar(r0, g0s)
    Am = amplitude_vfar(rm, gms)
    if Ae < 1e-8: return None
    G3 = (Am/Ae)**4
    return {'alpha': alpha, 'z0': z0, 'mult': multiplier, 'g0_tau': multiplier*z0,
            'Ae': Ae, 'Am': Am, 'G3': G3, 'G3_minus_R31': G3-R31}

def worker_mult(args):
    alpha, multiplier = args
    try: return compute_G3_multiplier(alpha, multiplier)
    except: return None

def worker_F(args):
    alpha, = args
    try: return compute_F_vfar(alpha)
    except: return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6

    print("=" * 70)
    print("ex102_hierarchy.py — Inne hierarchie + α_min_far (S4+S5, v39)")
    print("=" * 70)
    print(f"  R21={R21:.6f}, R31={R31:.6f}, R32={R32:.6f}")
    print(f"  Okna VFAR = {FIT_WINS_VFAR} (bez biasowanego [60,76])")
    print()

    test_alphas = [ALPHA_STAR1_FAR, ALPHA_STAR2_FAR, 2.50, 2.75, 3.00]

    # ----------------------------------------------------------------
    # FAZA 1: G₃ z różnymi mnożnikami przy wybranych α
    # ----------------------------------------------------------------
    print("=" * 70)
    print("FAZA 1: G₃(m·z₀) z oknami VFAR — różne mnożniki m")
    print("=" * 70)
    print(f"  {'α':>9} | {'m·z₀':>7} | {'Mnoznik':>8} | {'Am_vfar':>10} | {'G3':>12} | {'G3-R31':>12}")
    print(f"  {'-'*72}")

    for alpha in test_alphas:
        z0 = find_z0(alpha)
        if z0 is None: continue
        r0, g0s = integrate_soliton(z0, alpha)
        Ae = amplitude_vfar(r0, g0s)
        print(f"\n  α={alpha:.7f}  z₀={z0:.7f}  Ae_vfar={Ae:.6f}")
        for name, mult in MULTIPLIERS:
            res = compute_G3_multiplier(alpha, mult)
            if res:
                mk = " ← BINGO!" if abs(res['G3_minus_R31']) < 50 else ""
                mk2 = " ← bliskie" if 50 <= abs(res['G3_minus_R31']) < 500 else ""
                print(f"  {'':>9}   {res['g0_tau']:7.4f} | {name:>8} | {res['Am']:10.6f} | "
                      f"{res['G3']:12.4f} | {res['G3_minus_R31']:+12.4f}{mk}{mk2}")

    # ----------------------------------------------------------------
    # FAZA 2: Skan G₃(φ²·z₀) z VFAR oknami na α∈[2.3, 5.0]
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 2: G₃_vfar(α) na α∈[2.3, 5.0], krok=0.1 (VFAR okna)")
    print("=" * 70)
    print(f"  Okna VFAR [80,156] — kontrola biasu [60,76]")
    print()

    alphas_2 = np.arange(2.30, 5.05, 0.10)
    args_2 = [(a, PHI2) for a in alphas_2]
    print(f"  {len(alphas_2)} punktów, {ncpu} rdzeni...")
    with mp.Pool(ncpu) as pool:
        res_2 = pool.map(worker_mult, args_2)

    good_2 = [r for r in res_2 if r]
    print(f"  Wyniki: {len(good_2)}/{len(alphas_2)}")
    print()
    print(f"  {'α':>6} | {'z₀':>8} | {'G3_vfar':>12} | {'G3_vfar-R31':>13} | {'Am_vfar':>10}")
    print(f"  {'-'*62}")
    min_diff = float('inf'); min_r = None
    sign_chg = []; prev_s=None; prev_a=None
    for r in sorted(good_2, key=lambda x: x['alpha']):
        a = r['alpha']
        d = r['G3_minus_R31']
        cur_s = np.sign(d)
        mk = " ←" if (prev_s and cur_s != prev_s) else ""
        if prev_s and cur_s != prev_s: sign_chg.append((prev_a, a))
        if abs(d) < abs(min_diff): min_diff = d; min_r = r
        print(f"  {a:6.2f} | {r['z0']:8.5f} | {r['G3']:12.4f} | {d:+13.4f} | {r['Am']:10.6f}{mk}")
        prev_s = cur_s; prev_a = a

    print(f"\n  Zmiany znaku G₃_vfar−R₃₁: {sign_chg}")
    if min_r:
        print(f"  Min |G₃_vfar−R₃₁| = {abs(min_diff):.4f} przy α={min_r['alpha']:.2f}")

    # ----------------------------------------------------------------
    # FAZA 3: Precyzyjny α_min_far (S5) z VFAR oknami
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 3: Precyzyjny α_min_far z oknami VFAR (S5)")
    print("=" * 70)

    def F_vfar_func(alpha):
        res = compute_F_vfar(alpha)
        return res['F'] if res else float('nan')

    eps = 0.003
    def dF_dalpha(alpha):
        fp = F_vfar_func(alpha + eps)
        fm = F_vfar_func(alpha - eps)
        if np.isnan(fp) or np.isnan(fm): return float('nan')
        return (fp - fm) / (2*eps)

    # Dense skan F_vfar na [2.4, 2.8]
    alphas_F = np.arange(2.40, 2.81, 0.01)
    print(f"  Dense skan F_vfar na [2.40, 2.80], krok=0.01...")
    F_vals = []
    for a in alphas_F:
        res = compute_F_vfar(a)
        if res: F_vals.append((a, res['F']))

    if F_vals:
        amin_idx = min(range(len(F_vals)), key=lambda i: F_vals[i][1])
        a_min_v, F_min_v = F_vals[amin_idx]
        print(f"  Min F_vfar (krok=0.01): α={a_min_v:.4f}, F_vfar={F_min_v:.6f}")
        print(f"  R₂₁ = {R21:.6f}")
        print(f"  F_min_vfar − R₂₁ = {F_min_v-R21:+.6f}")

        # Porównanie FAR vs VFAR alpha_min
        print(f"\n  Porównanie α_min:")
        print(f"    α_min_far  (ex100)  = 2.5595170")
        print(f"    α_min_vfar (ex102)  = {a_min_v:.7f}")
        print(f"    α*₁_far    (ex98)   = {ALPHA_STAR1_FAR:.7f}")
        print(f"    α*₂_far    (ex98)   = {ALPHA_STAR2_FAR:.7f}")
        print(f"    midpoint   α*₁₊₂/2 = {(ALPHA_STAR1_FAR+ALPHA_STAR2_FAR)/2:.7f}")

        # Brentq dF_vfar/dα = 0
        a_lo_b = max(2.40, a_min_v - 0.10)
        a_hi_b = min(2.80, a_min_v + 0.10)
        d_lo = dF_dalpha(a_lo_b)
        d_hi = dF_dalpha(a_hi_b)
        print(f"\n  dF_vfar/dα na [{a_lo_b:.3f},{a_hi_b:.3f}]:")
        print(f"    dF/dα({a_lo_b:.3f}) = {d_lo:+.5f}")
        print(f"    dF/dα({a_hi_b:.3f}) = {d_hi:+.5f}")
        if d_lo * d_hi < 0:
            try:
                a_min_p = brentq(dF_dalpha, a_lo_b, a_hi_b, xtol=1e-6, maxiter=100)
                res_p = compute_F_vfar(a_min_p)
                print(f"\n  ✓ α_min_vfar (brentq) = {a_min_p:.7f}")
                if res_p:
                    print(f"    F_min_vfar = {res_p['F']:.7f}")
                    print(f"    F_min_vfar − R₂₁ = {res_p['F']-R21:+.7f}")
                    print(f"    z₀ = {res_p['z0']:.7f}")
            except Exception as e:
                print(f"  brentq failed: {e}")

    # ----------------------------------------------------------------
    # FAZA 4: F_vfar vs F_far — porównanie przy α∈[2.3,2.8]
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 4: F_vfar vs F_far — czy okna zmieniają α*₁,₂?")
    print("=" * 70)
    print(f"\n  {'α':>6} | {'F_far':>10} | {'F_vfar':>10} | {'ΔF':>8} | {'F-R21(vfar)':>12}")
    print(f"  {'-'*58}")

    prev_sf=None; sign_vfar=[]
    for a in np.arange(2.35, 2.82, 0.05):
        res_far  = compute_G3_multiplier(a, PHI)  # F = G3 z mnoznikiem phi
        res_vfar = compute_F_vfar(a)
        # Oblicz F_far też
        z0 = find_z0(a)
        if z0 is None: continue
        r0, g0s = integrate_soliton(z0, a)
        r1, g1s = integrate_soliton(PHI*z0, a)
        Ae_f  = amplitude_far(r0, g0s)
        Amu_f = amplitude_far(r1, g1s)
        F_far_v = (Amu_f/Ae_f)**4 if Ae_f > 1e-8 else float('nan')
        F_vfar_v = res_vfar['F'] if res_vfar else float('nan')
        delta = F_vfar_v - F_far_v if not (np.isnan(F_far_v) or np.isnan(F_vfar_v)) else float('nan')
        cur_sf = np.sign(F_vfar_v-R21) if not np.isnan(F_vfar_v) else None
        mk = " ←" if (prev_sf and cur_sf and cur_sf != prev_sf) else ""
        if prev_sf and cur_sf and cur_sf != prev_sf: sign_vfar.append(a)
        print(f"  {a:6.3f} | {F_far_v:10.4f} | {F_vfar_v:10.4f} | {delta:+8.4f} | {F_vfar_v-R21:+12.4f}{mk}")
        prev_sf = cur_sf

    if sign_vfar:
        print(f"\n  Zmiany znaku F_vfar−R₂₁ przy α ≈ {sign_vfar}")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("PODSUMOWANIE ex102 — S4+S5")
    print("=" * 70)
    print(f"Czas: {time.time()-t0:.1f} s")
    print()
    if sign_chg:
        print(f"  G₃_vfar zmiana znaku przy α∈{sign_chg} → możliwe α*₃_vfar!")
    else:
        print(f"  Brak zmiany znaku G₃_vfar−R₃₁ na [2.3,5.0] z VFAR oknami.")
        if min_r:
            print(f"  Min G₃_vfar = {min_r['G3']:.2f} przy α={min_r['alpha']:.2f} "
                  f"(+{min_diff:.2f} nad R₃₁)")
    print("=" * 70)

if __name__ == '__main__':
    mp.freeze_support()
    main()
