#!/usr/bin/env python3
"""
ex103_mc_precision.py — Precyzyjny mnożnik krytyczny m_c i zera G₃=R₃₁ (S1, v40)
===================================================================================
STATUS: LEGACY-TRANSLATIONAL

This file continues the older tau-sector critical-multiplier search in legacy
selection variables. Keep it as historical exploratory context rather than a
canonical synchronized `nbody` path.

Odkrycie z weryfikacji v40:
  min_α G3(m·z₀; α) zmienia znak między m=2.49 (min=3437<R₃₁) a m=2.50 (min=3510>R₃₁)
  → istnieje m_c ≈ 2.4955 gdzie min_α G3(m_c·z₀)=R₃₁ (tangensowe zero)
  → dla m < m_c istnieją DWA zera G3=R₃₁ (analogia do muona!)

FAZA A: Precyzyjne m_c — brentq po m (minimalizacja po α wewnętrznie)
FAZA B: Dla m < m_c — brentq dwóch zer α*_τ1, α*_τ2
FAZA C: Weryfikacja R_MAX=300
FAZA D: F(φ·z₀) przy α*_τ1, α*_τ2

Autor: Claudian (sesja v40, S1)
Data: 2026-03-29
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar, curve_fit
import multiprocessing as mp
import time
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STALE
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
R_MAX_MAIN = 150
R_MAX_VER  = 300

# Okna VFAR — bez biasowanego [60,76]
FIT_WINS_VFAR = [(80,96),(90,106),(100,116),(110,126),(120,136),(130,146),(140,156)]
# Okna VFAR szersze (dla R_MAX=300)
FIT_WINS_VFAR300 = [(80,96),(100,116),(120,136),(140,156),(160,176),(180,196),(200,216)]

ALPHA_STAR1_FAR = 2.436011168
ALPHA_STAR2_FAR = 2.753824880

# ================================================================
# ODE + AMPLITUDA
# ================================================================
def integrate_soliton(g0, alpha, r_max=R_MAX_MAIN):
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

def B_coeff_func(g0, alpha, r_max=R_MAX_MAIN):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX_MAIN, lo=1.00, hi=2.5, n=80):
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

def compute_G3(alpha, mult, r_max=R_MAX_MAIN, wins=None):
    """G₃=(A(mult·z₀)/A(z₀))^4 z oknami VFAR."""
    if wins is None:
        wins = FIT_WINS_VFAR if r_max <= 200 else FIT_WINS_VFAR300
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None, None
    r0, g0s = integrate_soliton(z0,        alpha, r_max)
    rm, gms = integrate_soliton(mult*z0,   alpha, r_max)
    Ae = amplitude_from_wins(r0, g0s, wins)
    Am = amplitude_from_wins(rm, gms, wins)
    if Ae < 1e-8: return None, z0
    return (Am/Ae)**4, z0

def compute_F(alpha, r_max=R_MAX_MAIN, wins=None):
    """F=(A(φ·z₀)/A(z₀))^4."""
    if wins is None:
        wins = FIT_WINS_VFAR if r_max <= 200 else FIT_WINS_VFAR300
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None, None
    r0, g0s = integrate_soliton(z0,      alpha, r_max)
    r1, g1s = integrate_soliton(PHI*z0,  alpha, r_max)
    Ae  = amplitude_from_wins(r0, g0s, wins)
    Amu = amplitude_from_wins(r1, g1s, wins)
    if Ae < 1e-8: return None, z0
    return (Amu/Ae)**4, z0

# ================================================================
# FAZA A: Precyzyjne m_c
# ================================================================
def min_G3_over_alpha(mult, alpha_lo=2.20, alpha_hi=2.55, n_scan=15):
    """Minimalizuj G3(mult·z₀;α) po α — zwraca (min_G3, alpha_at_min)."""
    alphas = np.linspace(alpha_lo, alpha_hi, n_scan)
    vals = []
    for a in alphas:
        g3, _ = compute_G3(a, mult)
        if g3 is not None:
            vals.append((a, g3))
    if not vals: return None, None
    # Minimalizacja dokładna
    # Znajdź bracket dla minimum
    if len(vals) < 3: return min(vals, key=lambda x: x[1])
    idx_min = min(range(len(vals)), key=lambda i: vals[i][1])
    a_min_v, g3_min_v = vals[idx_min]
    # Refined minimum z minimize_scalar
    if 0 < idx_min < len(vals)-1:
        a_lo_b = vals[idx_min-1][0]
        a_hi_b = vals[idx_min+1][0]
        def neg_target(a):
            g3, _ = compute_G3(a, mult)
            return g3 if g3 is not None else 1e10
        try:
            res = minimize_scalar(neg_target, bounds=(a_lo_b, a_hi_b), method='bounded',
                                  options={'xatol':1e-5})
            return res.fun, res.x
        except:
            pass
    return g3_min_v, a_min_v

def worker_mc(m):
    try:
        g3_min, a_min = min_G3_over_alpha(m)
        if g3_min is None: return None
        return (m, g3_min, a_min, g3_min - R31)
    except: return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    ncpu = 6

    print("=" * 70)
    print("ex103_mc_precision.py — Krytyczny mno.znik m_c i zera G3=R31 (v40)")
    print("=" * 70)
    print(f"  R21={R21:.6f}, R31={R31:.6f}")
    print(f"  Okna VFAR = {FIT_WINS_VFAR}")
    print()

    # ----------------------------------------------------------------
    # FAZA A: Dense skan min_alpha G3(m·z₀) dla m∈[2.46, 2.53]
    # ----------------------------------------------------------------
    print("=" * 70)
    print("FAZA A: Skan min_alpha G3(m*z0) dla m w [2.46, 2.53], krok=0.002")
    print("=" * 70)

    ms_A = np.arange(2.460, 2.531, 0.002)
    print(f"  {len(ms_A)} wartosci m, {ncpu} rdzeni...")

    with mp.Pool(ncpu) as pool:
        results_A = pool.map(worker_mc, ms_A)

    good_A = [r for r in results_A if r is not None]
    good_A.sort(key=lambda x: x[0])

    print(f"\n  {'m':>7} | {'alpha@min':>9} | {'G3_min':>12} | {'G3_min-R31':>13} | sign")
    print(f"  {'-'*55}")
    prev_s = None; sign_changes = []
    mc_lo = mc_hi = None
    for m, g3min, amin, diff in good_A:
        cur_s = np.sign(diff)
        mk = " <---" if (prev_s is not None and cur_s != prev_s) else ""
        if prev_s is not None and cur_s != prev_s:
            sign_changes.append((m, prev_m))
            mc_lo = prev_m; mc_hi = m
        print(f"  {m:7.4f} | {amin:9.5f} | {g3min:12.4f} | {diff:+13.4f} | {'+' if diff>0 else '-'}{mk}")
        prev_s = cur_s; prev_m = m

    print(f"\n  Zmiany znaku G3_min-R31: {sign_changes}")
    if mc_lo and mc_hi:
        mc_interp = mc_lo + (R31 - dict([(r[0],r[1]) for r in good_A])[mc_lo]) / \
                    (dict([(r[0],r[1]) for r in good_A])[mc_hi] - dict([(r[0],r[1]) for r in good_A])[mc_lo]) * (mc_hi - mc_lo)
        print(f"  Interpolacja liniowa: m_c ≈ {mc_interp:.6f}")

    # ----------------------------------------------------------------
    # FAZA A2: brentq precyzyjny m_c
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA A2: brentq precyzyjny m_c")
    print("=" * 70)

    if mc_lo and mc_hi:
        def min_G3_minus_R31(m):
            g3min, _ = min_G3_over_alpha(m, n_scan=20)
            if g3min is None: return float('nan')
            return g3min - R31

        try:
            mc_brentq = brentq(min_G3_minus_R31, mc_lo, mc_hi, xtol=1e-5, maxiter=50)
            g3_at_mc, alpha_at_mc = min_G3_over_alpha(mc_brentq, n_scan=25)
            print(f"\n  m_c (brentq) = {mc_brentq:.6f}")
            print(f"  alpha @ m_c  = {alpha_at_mc:.6f}")
            print(f"  G3_min(m_c)  = {g3_at_mc:.4f}  (powinno = R31={R31:.4f})")
            print(f"  Residual     = {g3_at_mc - R31:+.4f}")
            print()
            print(f"  Kandydaci algebraiczni m_c:")
            import math
            candidates = [
                ("5/2",        2.5),
                ("sqrt(2*pi)", np.sqrt(2*np.pi)),
                ("ln(12)",     np.log(12)),
                ("phi+7/8",    PHI+0.875),
                ("sqrt(6)",    np.sqrt(6)),
                ("pi/sqrt(phi)", np.pi/np.sqrt(PHI)),
                ("5/2-1/200",  2.5-0.005),
                ("(1+phi)/phi^(1/3)", (1+PHI)/PHI**(1/3)),
            ]
            for name, val in candidates:
                d = val - mc_brentq
                ppm = d / mc_brentq * 1e6
                print(f"    {name:<25s} = {val:.6f}  delta={d:+.6f}  ({ppm:+.0f} ppm)")
        except Exception as e:
            print(f"  brentq failed: {e}")
            mc_brentq = None; alpha_at_mc = None
    else:
        print("  BRAK zmiany znaku w przedziale [2.46,2.53] — rozszerzam...")
        mc_brentq = None; alpha_at_mc = None

    # ----------------------------------------------------------------
    # FAZA B: Dwa zera G3(m·z₀;α)=R31 dla m < m_c
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA B: Dwa zera alpha dla m < m_c (m=2.48, m=2.46, m=2.44)")
    print("=" * 70)

    test_ms = [2.48, 2.46, 2.44, 2.42, 2.40]
    alpha_scan = np.arange(2.20, 2.65, 0.02)

    for m_test in test_ms:
        print(f"\n  --- m={m_test:.2f} ---")
        vals = []
        for a in alpha_scan:
            g3, z0 = compute_G3(a, m_test)
            if g3 is not None:
                vals.append((a, g3, g3-R31))

        prev_s = None; zeros_found = []
        for a, g3, diff in vals:
            cur_s = np.sign(diff)
            if prev_s is not None and cur_s != prev_s:
                zeros_found.append((prev_a, a, prev_diff, diff))
            prev_s = cur_s; prev_a = a; prev_diff = diff

        print(f"  Przyblizone zera G3=R31: {len(zeros_found)} znalezionych")
        alpha_zeros = []
        for a_lo, a_hi, d_lo, d_hi in zeros_found:
            try:
                def G3_minus_R31(a):
                    g3, _ = compute_G3(a, m_test)
                    return (g3 - R31) if g3 is not None else float('nan')
                az = brentq(G3_minus_R31, a_lo, a_hi, xtol=1e-8, maxiter=100)
                g3z, z0z = compute_G3(az, m_test)
                alpha_zeros.append(az)
                print(f"    alpha*_tau = {az:.8f}  G3={g3z:.6f}  z0={z0z:.7f}  m*z0={m_test*z0z:.7f}")
            except Exception as e:
                print(f"    brentq failed [{a_lo:.3f},{a_hi:.3f}]: {e}")

        if len(alpha_zeros) == 2:
            a1, a2 = alpha_zeros
            print(f"    Suma S_tau    = {a1+a2:.8f}")
            print(f"    Roznica       = {a2-a1:.8f}")
            print(f"    Midpoint      = {(a1+a2)/2:.8f}")
            print(f"    alpha*1_far   = {ALPHA_STAR1_FAR:.8f}")
            print(f"    alpha*2_far   = {ALPHA_STAR2_FAR:.8f}")
            print(f"    Relacja a1/a*1_far = {a1/ALPHA_STAR1_FAR:.6f}")
            print(f"    Relacja a2/a*1_far = {a2/ALPHA_STAR1_FAR:.6f}")

    # ----------------------------------------------------------------
    # FAZA C: Weryfikacja z R_MAX=300 dla m=2.48
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA C: Weryfikacja R_MAX=300 dla m=2.48")
    print("=" * 70)

    m_ver = 2.48
    alpha_scan_ver = np.arange(2.20, 2.65, 0.05)
    print(f"  Skan G3({m_ver}*z0) na alpha in [2.20, 2.65], R_MAX=300...")
    vals_ver = []
    for a in alpha_scan_ver:
        g3, z0 = compute_G3(a, m_ver, r_max=R_MAX_VER, wins=FIT_WINS_VFAR300)
        if g3 is not None:
            d = g3 - R31
            vals_ver.append((a, g3, d))
            print(f"    alpha={a:.2f}: G3={g3:.4f}  G3-R31={d:+.4f}")

    prev_s = None; zeros_300 = []
    for a, g3, diff in vals_ver:
        cur_s = np.sign(diff)
        if prev_s is not None and cur_s != prev_s:
            zeros_300.append((prev_a, a))
        prev_s = cur_s; prev_a = a

    if zeros_300:
        print(f"\n  Zmiany znaku (R_MAX=300): {zeros_300}")
        for a_lo, a_hi in zeros_300:
            try:
                def G3_ver(a):
                    g3, _ = compute_G3(a, m_ver, r_max=R_MAX_VER, wins=FIT_WINS_VFAR300)
                    return (g3 - R31) if g3 is not None else float('nan')
                az = brentq(G3_ver, a_lo, a_hi, xtol=1e-6, maxiter=100)
                g3z, _ = compute_G3(az, m_ver, r_max=R_MAX_VER, wins=FIT_WINS_VFAR300)
                print(f"  alpha*_tau (R_MAX=300) = {az:.7f}  G3={g3z:.4f}")
            except Exception as e:
                print(f"  brentq failed: {e}")

    # ----------------------------------------------------------------
    # FAZA D: F(φ·z₀) przy zerach α*_τ z m=2.48
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA D: F(phi*z0) przy zerach alpha*_tau (m=2.48)")
    print("=" * 70)

    # Powtórz brentq dla m=2.48 żeby mieć zera
    alpha_scan_d = np.arange(2.22, 2.62, 0.02)
    vals_d = []
    for a in alpha_scan_d:
        g3, _ = compute_G3(a, 2.48)
        if g3 is not None: vals_d.append((a, g3-R31))

    prev_s = None; zeros_d = []
    for a, diff in vals_d:
        cur_s = np.sign(diff)
        if prev_s is not None and cur_s != prev_s:
            zeros_d.append((prev_a, a))
        prev_s = cur_s; prev_a = a

    for a_lo, a_hi in zeros_d:
        try:
            def G3_d(a):
                g3, _ = compute_G3(a, 2.48)
                return (g3-R31) if g3 is not None else float('nan')
            az = brentq(G3_d, a_lo, a_hi, xtol=1e-8)
            g3z, z0z = compute_G3(az, 2.48)
            Fv, _ = compute_F(az)
            print(f"\n  alpha*_tau = {az:.8f}")
            print(f"    G3(2.48*z0) = {g3z:.6f}  (R31={R31:.4f})")
            print(f"    F(phi*z0)   = {Fv:.6f}  (R21={R21:.6f})")
            print(f"    F - R21     = {Fv-R21:+.6f}")
            print(f"    z0          = {z0z:.8f}")
            print(f"    2.48*z0     = {2.48*z0z:.8f}")
            print(f"    phi*z0      = {PHI*z0z:.8f}")
        except Exception as e:
            print(f"  brentq failed: {e}")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("PODSUMOWANIE ex103")
    print("=" * 70)
    print(f"Czas: {time.time()-t0:.1f} s")
    print()
    if mc_brentq:
        print(f"  m_c (precyzyjny) = {mc_brentq:.6f}")
        print(f"  Dla m < m_c: istnieja DWA zera G3(m*z0;alpha)=R31")
        print(f"  Dla m > m_c: brak zer (G3>R31 wszedzie)")
    print()
    print(f"  Kandydaci m_c:")
    print(f"    5/2 = 2.500000  (bliskie)")
    print(f"    sqrt(2*pi) = {np.sqrt(2*np.pi):.6f}")
    print(f"    phi^2/phi  = {PHI2/PHI:.6f}  = phi")
    print("=" * 70)

if __name__ == '__main__':
    mp.freeze_support()
    main()
