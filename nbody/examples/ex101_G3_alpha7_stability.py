#!/usr/bin/env python3
"""
ex101_G3_alpha7_stability.py — Rezonans czy fizyka? G₃_far≈R₃₁ przy α≈7.1 (S2, v39)
======================================================================================
Z ex99: G₃_far(7.0)=781 << R₃₁, G₃_far(7.5)=13640 >> R₃₁ — skok o czynnik 17.5.
Pytanie: czy to rezonans fazowy pomiaru (okna FAR [60,136] trafiają w węzeł A(r)),
         czy fizyczne przejście przez R₃₁?

Test rezonansu: per-okno A_tau przy α=6.8, 7.0, 7.1, 7.2, 7.5
  - Jeśli różne okna dają drastycznie różne A → rezonans fazowy
  - Jeśli plateau → fizyczne

Test stabilności R_MAX:
  - G₃_far(α, R_MAX=150) vs G₃_far(α, R_MAX=300) w okolicach α≈7.1
  - Jeśli α* zmienia się o >0.1 przy zmianie R_MAX → niestabilne (rezonans)

Autor: Claudian (sesja v39, S2)
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

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0

# Trzy zestawy okien do porównania
FIT_WINS_FAR  = [(60,76),(70,86),(80,96),(90,106),(100,116),(110,126),(120,136)]
FIT_WINS_MID  = [(30,46),(40,56),(50,66),(60,76),(70,86),(80,96)]
FIT_WINS_VFAR = [(100,116),(110,126),(120,136),(130,146),(140,156),(150,166),(160,176)]  # wymaga R_MAX≥180

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
    if not Av: return 0.0, []
    Av_arr = np.array(Av)
    if len(Av_arr) < 3 or np.std(Av_arr)/(np.mean(Av_arr)+1e-12) < 0.005:
        return float(np.median(Av_arr)), list(zip(rv, Av))
    try:
        p, _ = curve_fit(lambda x,ai,a,b: ai*(1+a/x+b/x**2),
                         np.array(rv), Av_arr, p0=[Av_arr[-1],0.,0.], maxfev=2000)
        return float(p[0]), list(zip(rv, Av))
    except:
        return float(np.median(Av_arr)), list(zip(rv, Av))

def amplitude_per_window(r, g, wins=None):
    """Zwróć listę (rL, rR, A) dla każdego okna z osobna."""
    if wins is None:
        wins = [(r0, r0+16) for r0 in range(20, 271, 10)]
    res = []
    for rL, rR in wins:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        A = float(np.sqrt(bc[0]**2 + bc[1]**2))
        res.append((rL, rR, A))
    return res

def B_coeff_func(g0, alpha, r_max):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, lo=1.00, hi=2.5, n=80):
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

def compute_G3_multiwins(alpha, r_max):
    """G₃ z FAR, MID, VFAR oknami + per-okno profil."""
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0s = integrate_soliton(z0,      alpha, r_max)
    r2, g2s = integrate_soliton(PHI2*z0, alpha, r_max)

    # Per-window profil (co 10 r, szerokość 16)
    probe_wins = [(rr, rr+16) for rr in range(20, min(int(r_max)-20, 281), 10)]
    pw_e   = amplitude_per_window(r0, g0s, probe_wins)
    pw_tau = amplitude_per_window(r2, g2s, probe_wins)

    results = {}
    for label, wins in [('far',FIT_WINS_FAR),('mid',FIT_WINS_MID),('vfar',FIT_WINS_VFAR)]:
        Ae,  _ = amplitude_from_wins(r0, g0s, wins)
        Atau,_ = amplitude_from_wins(r2, g2s, wins)
        G3 = (Atau/Ae)**4 if Ae > 1e-8 else float('nan')
        results[label] = {'Ae': Ae, 'Atau': Atau, 'G3': G3}

    return {'alpha': alpha, 'z0': z0, 'r_max': r_max,
            'pw_e': pw_e, 'pw_tau': pw_tau,
            **results}

def worker_scan(args):
    alpha, r_max = args
    try:
        z0 = find_z0(alpha, r_max)
        if z0 is None: return None
        r0, g0s = integrate_soliton(z0,      alpha, r_max)
        r2, g2s = integrate_soliton(PHI2*z0, alpha, r_max)
        Ae,  _ = amplitude_from_wins(r0, g0s, FIT_WINS_FAR)
        Atau,_ = amplitude_from_wins(r2, g2s, FIT_WINS_FAR)
        G3 = (Atau/Ae)**4 if Ae > 1e-8 else float('nan')
        return {'alpha': alpha, 'z0': z0, 'G3': G3, 'Ae': Ae, 'Atau': Atau}
    except: return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6
    R_MAX_1 = 150
    R_MAX_2 = 300

    print("=" * 70)
    print("ex101_G3_alpha7_stability.py — Rezonans czy fizyka? (S2, v39)")
    print("=" * 70)
    print(f"  R31={R31:.6f}, R21={R21:.6f}")
    print(f"  PHI²={PHI2:.8f}")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Dense skan α∈[6.5, 8.0], krok=0.05, R_MAX=150
    # ----------------------------------------------------------------
    print("=" * 70)
    print("FAZA 1: Dense skan G₃_far(α) na [6.5, 8.0], krok=0.05, R_MAX=150")
    print("=" * 70)

    alphas_1 = np.arange(6.50, 8.05, 0.05)
    args_1 = [(a, R_MAX_1) for a in alphas_1]
    print(f"  {len(alphas_1)} punktów, {ncpu} rdzeni...")
    with mp.Pool(ncpu) as pool:
        res_1 = pool.map(worker_scan, args_1)

    good_1 = [r for r in res_1 if r]
    print(f"  Wyniki: {len(good_1)}/{len(alphas_1)}")
    print()
    print(f"  {'α':>6} | {'z₀':>8} | {'Ae':>8} | {'Atau':>8} | {'G3_far':>12} | {'G3-R31':>11}")
    print(f"  {'-'*68}")

    sign_changes = []; prev_s=None; prev_a=None
    for r in sorted(good_1, key=lambda x: x['alpha']):
        a = r['alpha']
        g3d = r['G3'] - R31 if not np.isnan(r['G3']) else float('nan')
        cur_s = np.sign(g3d) if not np.isnan(g3d) else None
        mk = " ←" if (prev_s and cur_s and cur_s != prev_s) else ""
        if prev_s and cur_s and cur_s != prev_s: sign_changes.append((prev_a, a))
        g3_str = f"{r['G3']:12.4f}" if not np.isnan(r['G3']) else "         N/A"
        g3d_str = f"{g3d:+11.4f}" if not np.isnan(g3d) else "          N/A"
        print(f"  {a:6.2f} | {r['z0']:8.5f} | {r['Ae']:8.6f} | {r['Atau']:8.6f} | "
              f"{g3_str} | {g3d_str}{mk}")
        prev_s = cur_s; prev_a = a

    print(f"\n  Zmiany znaku G₃−R₃₁: {sign_changes}")

    # ----------------------------------------------------------------
    # FAZA 2: Per-okno profil A_tau przy α∈{6.8, 7.0, 7.1, 7.2, 7.5, 8.0}
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 2: Per-okno A_tau — test rezonansu fazowego")
    print("=" * 70)

    probe_alphas = [6.8, 7.0, 7.05, 7.1, 7.15, 7.2, 7.5, 8.0]
    for alpha in probe_alphas:
        print(f"\n  α={alpha:.3f}, R_MAX={R_MAX_1}")
        z0 = find_z0(alpha, R_MAX_1)
        if z0 is None:
            print(f"    z₀ nie znaleziona!"); continue
        print(f"    z₀={z0:.6f}, φ²·z₀={PHI2*z0:.6f}")
        r2, g2s = integrate_soliton(PHI2*z0, alpha, R_MAX_1)
        r0, g0s = integrate_soliton(z0,      alpha, R_MAX_1)
        Ae_med = np.median([A for _,_,A in amplitude_per_window(r0, g0s,
                            [(rr,rr+16) for rr in range(60,130,10)])])
        pw_tau = amplitude_per_window(r2, g2s, [(rr,rr+16) for rr in range(20, 141, 10)])
        print(f"    Ae(FAR median)={Ae_med:.6f}")
        print(f"    {'Okno':12s} | {'A_tau':>10} | {'(A/Ae)^4':>12} | Uwaga")
        print(f"    {'-'*55}")
        prev_A = None
        for rL, rR, A_t in pw_tau:
            G3_w = (A_t/Ae_med)**4 if Ae_med > 1e-8 else float('nan')
            jump = ""
            if prev_A is not None and abs(A_t-prev_A)/max(prev_A,1e-6) > 0.5:
                jump = " ← SKOK!"
            elif prev_A is not None and abs(A_t-prev_A)/max(prev_A,1e-6) > 0.1:
                jump = " ← zmiana >10%"
            print(f"    [{rL:3d},{rR:3d}]    | {A_t:10.6f} | {G3_w:12.4f} |{jump}")
            prev_A = A_t

    # ----------------------------------------------------------------
    # FAZA 3: Weryfikacja R_MAX=300 w okolicach zmiany znaku
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 3: Weryfikacja R_MAX=300 — czy G₃ przy α≈7.1 stabilne?")
    print("=" * 70)

    if sign_changes:
        a_check = [(a_lo+a_hi)/2 for a_lo,a_hi in sign_changes]
        # Dodaj kilka punktów wokół przejścia
        alpha_check = set()
        for a_lo, a_hi in sign_changes:
            for a in np.arange(a_lo-0.1, a_hi+0.15, 0.05):
                alpha_check.add(round(a, 3))
    else:
        # Sprawdź okolice α=7.0 gdzie był skok z ex99
        alpha_check = {6.8, 6.9, 7.0, 7.05, 7.1, 7.15, 7.2, 7.3, 7.4, 7.5}

    alphas_3 = sorted(alpha_check)
    print(f"  Sprawdzam α∈{alphas_3} z R_MAX=300, {ncpu} rdzeni...")
    args_3 = [(a, R_MAX_2) for a in alphas_3]
    with mp.Pool(ncpu) as pool:
        res_3 = pool.map(worker_scan, args_3)
    good_3 = {r['alpha']: r for r in res_3 if r}

    print(f"\n  Porównanie R_MAX=150 vs R_MAX=300:")
    print(f"  {'α':>6} | {'G3(150)':>12} | {'G3(300)':>12} | {'ΔG3':>10} | Stabilne?")
    print(f"  {'-'*60}")

    g3_150 = {r['alpha']: r['G3'] for r in good_1}
    sign_300 = []; prev_s3=None; prev_a3=None
    for a in sorted(good_3.keys()):
        r3 = good_3[a]
        g3_150v = g3_150.get(a, float('nan'))
        g3_300v = r3['G3']
        delta = g3_300v - g3_150v if not np.isnan(g3_150v) else float('nan')
        stable = "✓" if (not np.isnan(delta) and abs(delta)/max(abs(g3_150v),1) < 0.05) else "⚠"
        g3_150_str = f"{g3_150v:12.4f}" if not np.isnan(g3_150v) else "         N/A"
        delta_str  = f"{delta:+10.4f}" if not np.isnan(delta) else "        N/A"
        cur_s3 = np.sign(g3_300v - R31) if not np.isnan(g3_300v) else None
        if prev_s3 and cur_s3 and cur_s3 != prev_s3: sign_300.append((prev_a3, a))
        print(f"  {a:6.3f} | {g3_150_str} | {g3_300v:12.4f} | {delta_str} | {stable}")
        prev_s3 = cur_s3; prev_a3 = a

    print(f"\n  Zmiany znaku G₃−R₃₁ (R_MAX=300): {sign_300}")

    # Brentq z R_MAX=300 jeśli jest bracket
    if sign_300:
        for a_lo3, a_hi3 in sign_300:
            def g3_300_func(alpha):
                res = worker_scan((alpha, R_MAX_2))
                return (res['G3'] - R31) if res and not np.isnan(res['G3']) else float('nan')
            try:
                a_star_300 = brentq(g3_300_func, a_lo3, a_hi3, xtol=1e-7, maxiter=100)
                print(f"  ✓ α*₃(R_MAX=300) = {a_star_300:.7f}")
                # Porównaj z R_MAX=150
                if sign_changes:
                    for a_lo15, a_hi15 in sign_changes:
                        def g3_150_func(alpha):
                            res = worker_scan((alpha, R_MAX_1))
                            return (res['G3'] - R31) if res and not np.isnan(res['G3']) else float('nan')
                        try:
                            a_star_150 = brentq(g3_150_func, a_lo15, a_hi15, xtol=1e-7, maxiter=100)
                            print(f"  ✓ α*₃(R_MAX=150) = {a_star_150:.7f}")
                            print(f"  Δα* = {a_star_300-a_star_150:+.6f} "
                                  f"({'NIESTABILNE! rezonans' if abs(a_star_300-a_star_150)>0.05 else 'stabilne <0.05'})")
                        except: pass
            except Exception as e:
                print(f"  brentq R_MAX=300 failed: {e}")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("PODSUMOWANIE ex101 — S2: rezonans czy fizyka?")
    print("=" * 70)
    print(f"Czas: {time.time()-t0:.1f} s")
    print()
    print(f"  R_MAX=150 sign changes G₃−R₃₁: {sign_changes}")
    print(f"  R_MAX=300 sign changes G₃−R₃₁: {sign_300}")
    if sign_changes or sign_300:
        consistent = bool(sign_changes) and bool(sign_300)
        print(f"  Przejście spójne w obu R_MAX: {'TAK' if consistent else 'NIE'}")
        if consistent:
            print(f"  → Potencjalnie FIZYCZNE (wymaga dalszej weryfikacji)")
        else:
            print(f"  → NIESTABILNE — prawdopodobnie rezonans fazowy")
    else:
        print(f"  Brak przejścia G₃=R₃₁ w α∈[6.5,8.0] przy krok=0.05.")
        print(f"  → Przejście z ex99 (krok=0.5) mogło być artefaktem interpolacji")
    print("=" * 70)

if __name__ == '__main__':
    mp.freeze_support()
    main()
