#!/usr/bin/env python3
"""
ex88_fixed_window.py — Poprawiona ekstrapolacja α*₁,₂ do R_MAX→∞ (O-L1b, v36)
==============================================================================
Naprawa błędu ex87v3: okno B_coeff stałe [28,42], niezależne od R_MAX.

Strategia:
  1. Szeroki skan α∈[2.0, 3.5] dla każdego R_MAX → znajdź WSZYSTKIE skrzyżowania
     F(α)=R21 (może być wiele, bo F jest oscylacyjna)
  2. Klasyfikuj które skrzyżowania są "trwałe" (obecne dla ≥3 wartości R_MAX)
  3. Ekstrapoluj α*(R_MAX) → α*(∞) dla trwałych skrzyżowań

Okno B_coeff: [28, 42] — STAŁE dla wszystkich R_MAX.
Okna fit_amplitude: [(20,36),(30,46),(40,56),(50,66)] — STAŁE.

R_MAX_LIST = [100, 120, 150, 200, 300]

Autor: Claudian (sesja v36, O-L1b)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI        = (1.0 + np.sqrt(5.0)) / 2.0
R21        = 206.7682830
S_FORMULA  = 2.0 * np.pi - 1.1
G_OFF      = 0.005

# Okno B_coeff: STAŁE, niezależne od R_MAX
B_WIN_L = 28.0
B_WIN_R = 42.0

# Okna fit_amplitude: STAŁE
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]

RMAX_LIST = [100, 120, 150, 200, 300]

PASS_CNT = 0; FAIL_CNT = 0
def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    tag = "[PASS]" if cond else "[FAIL]"
    if cond: PASS_CNT += 1
    else:    FAIL_CNT += 1
    print(f"  {tag}  {name}")
    if detail: print(f"        {detail}")

# ================================================================
# ODE — identyczne z ex87v2 (tylko r_max jako parametr)
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
    for _ in range(14):   # więcej iteracji dla dużego R_MAX
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
    """Ekstrapolacja A_inf z stałych okien FIT_WINS."""
    Av, rv = [], []
    for rL, rR in FIT_WINS:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0
    if len(Av) < 2: return Av[-1]
    rv = np.array(rv); Av = np.array(Av)
    try:
        p, _ = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                         rv, Av, p0=[Av[-1], 0., 0.], maxfev=2000)
        return float(p[0])
    except:
        try:
            p, _ = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                             rv, Av, p0=[Av[-1], 0.], maxfev=1000)
            return float(p[0])
        except:
            return float(Av[-1])

def B_coeff(g0, alpha, r_max):
    """Amplituda cosinus w STAŁYM oknie [B_WIN_L, B_WIN_R]."""
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, lo=1.05, hi=2.5, n=40):
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

def fz0(alpha, r_max):
    z0 = find_z0(alpha, r_max)
    if z0 is None: return np.nan
    r, g   = integrate_soliton(z0, alpha, r_max)
    Ae = fit_amplitude(r, g)
    r2, g2 = integrate_soliton(PHI * z0, alpha, r_max)
    Am = fit_amplitude(r2, g2)
    if Ae < 1e-6: return np.nan
    return (Am / Ae) ** 4

# ================================================================
# SKAN: znajdź WSZYSTKIE przejścia przez R21 w [a_lo, a_hi]
# ================================================================
def scan_crossings(a_lo, a_hi, n_pts, r_max):
    """Skan α → znajdź wszystkie przejścia F(α)=R21."""
    alphas = np.linspace(a_lo, a_hi, n_pts)
    fvals  = []
    for a in alphas:
        try: fvals.append(fz0(a, r_max))
        except: fvals.append(np.nan)
    fvals = np.array(fvals)
    crossings = []
    for i in range(len(fvals) - 1):
        if np.isnan(fvals[i]) or np.isnan(fvals[i+1]): continue
        if (fvals[i] - R21) * (fvals[i+1] - R21) < 0:
            try:
                a_star = brentq(lambda a: fz0(a, r_max) - R21,
                                alphas[i], alphas[i+1],
                                xtol=1e-11, rtol=1e-12, maxiter=300)
                crossings.append(a_star)
            except: pass
    return alphas, fvals, crossings

def worker_scan(args):
    a_lo, a_hi, n_pts, r_max = args
    alphas, fvals, crossings = scan_crossings(a_lo, a_hi, n_pts, r_max)
    return r_max, alphas, fvals, crossings

# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("  ex88_fixed_window — Ekstrapolacja α*₁,₂ (O-L1b, stałe okno)")
    print("=" * 70)
    print(f"  R21       = {R21:.8f}")
    print(f"  S_formula = 2π−11/10 = {S_FORMULA:.10f}")
    print(f"  B_WIN     = [{B_WIN_L}, {B_WIN_R}]  (STAŁE)")
    print(f"  FIT_WINS  = {FIT_WINS}  (STAŁE)")
    print(f"  R_MAX_LIST= {RMAX_LIST}")

    N_CPU = min(5, cpu_count())
    print(f"  Rdzenie   = {N_CPU}")

    # ----------------------------------------------------------------
    # Krok 1: szeroki skan [2.0, 3.5] dla każdego R_MAX
    # ----------------------------------------------------------------
    SCAN_LO = 2.0
    SCAN_HI = 3.5
    N_SCAN  = 60   # dość gęsto, by złapać wszystkie przejścia

    print(f"\n=== KROK 1: Szeroki skan [{SCAN_LO}, {SCAN_HI}] ({N_SCAN} pkt) ===")

    tasks = [(SCAN_LO, SCAN_HI, N_SCAN, rmax) for rmax in RMAX_LIST]
    with Pool(N_CPU) as pool:
        scan_results = pool.map(worker_scan, tasks)

    # Wyświetl wyniki skanu i zbierz przejścia
    all_crossings = {}   # r_max → lista α*
    for r_max, alphas, fvals, crossings in sorted(scan_results, key=lambda x: x[0]):
        all_crossings[r_max] = crossings
        print(f"\n  R_MAX={r_max}: znaleziono {len(crossings)} przejść:")
        for a in crossings:
            f_v = fz0(a, r_max)
            print(f"    α* = {a:.8f},  F={f_v:.4f},  F-R21={f_v-R21:+.4f}")

    # ----------------------------------------------------------------
    # Krok 2: Klasyfikacja przejść — które są "trwałe"?
    # ----------------------------------------------------------------
    print("\n=== KROK 2: Klasyfikacja trwałości przejść ===")
    print("  (trwałe = znalezione dla ≥3 wartości R_MAX, zbieżne do wspólnej wartości)")

    # Grupuj przejścia w "rodziny" (różniące się o < 0.1 między R_MAX)
    families = []
    for r_max in sorted(all_crossings.keys()):
        for a in all_crossings[r_max]:
            # Czy pasuje do istniejącej rodziny?
            matched = False
            for fam in families:
                # Porównaj z dowolnym elementem rodziny
                for rmax_fam, a_fam in fam.items():
                    if abs(a - a_fam) < 0.12:
                        fam[r_max] = a
                        matched = True
                        break
                if matched: break
            if not matched:
                families.append({r_max: a})

    print(f"\n  Łącznie {len(families)} rodzin przejść:")
    for i, fam in enumerate(families):
        count = len(fam)
        vals  = sorted(fam.items())
        a_min = min(fam.values()); a_max = max(fam.values())
        print(f"  Rodzina {i+1}: {count} pkt, α∈[{a_min:.4f}, {a_max:.4f}]  ",
              end="")
        for rmax, a in vals:
            print(f"  R={rmax}:{a:.4f}", end="")
        print()

    # Wybierz trwałe rodziny (≥3 R_MAX)
    stable = [fam for fam in families if len(fam) >= 3]
    stable.sort(key=lambda f: min(f.values()))
    print(f"\n  Trwałe rodziny (≥3 R_MAX): {len(stable)}")
    for i, fam in enumerate(stable):
        print(f"  Stable {i+1}: R_MAX={sorted(fam.keys())} → "
              f"α={sorted(fam.values())}")

    # ----------------------------------------------------------------
    # Krok 3: Ekstrapolacja dla trwałych rodzin
    # ----------------------------------------------------------------
    print("\n=== KROK 3: Ekstrapolacja α*(R_MAX) → α*(∞) ===")

    def extrapolate_family(fam, label):
        rmax_arr = np.array(sorted(fam.keys()), dtype=float)
        a_arr    = np.array([fam[r] for r in rmax_arr.astype(int)], dtype=float)
        print(f"\n  {label}:")
        for rm, av in zip(rmax_arr, a_arr):
            print(f"    R_MAX={rm:.0f}: α*={av:.8f}")
        if len(rmax_arr) < 3:
            print(f"  Za mało punktów.")
            return None, None
        try:
            p1, cov1 = curve_fit(lambda x, ai, c: ai + c/x,
                                  rmax_arr, a_arr, p0=[a_arr[-1], 0.])
            a_inf = p1[0]; err  = np.sqrt(cov1[0, 0])
            print(f"  Fit 1/R: α*_∞ = {a_inf:.8f} ± {err:.3e}")
            if len(rmax_arr) >= 4:
                p2, cov2 = curve_fit(lambda x, ai, c, d: ai + c/x + d/x**2,
                                      rmax_arr, a_arr, p0=[a_arr[-1], 0., 0.])
                a_inf2 = p2[0]; err2 = np.sqrt(cov2[0, 0])
                print(f"  Fit 1/R²: α*_∞ = {a_inf2:.8f} ± {err2:.3e}")
                return a_inf2, err2
            return a_inf, err
        except Exception as e:
            print(f"  Ekstrapolacja FAIL: {e}")
            return None, None

    extrap_results = []
    for i, fam in enumerate(stable):
        a_inf, err = extrapolate_family(fam, f"Rodzina {i+1}")
        if a_inf is not None:
            extrap_results.append((a_inf, err, i+1))

    # ----------------------------------------------------------------
    # Krok 4: Analiza sumy
    # ----------------------------------------------------------------
    print("\n" + "=" * 70)
    print("  WYNIKI KOŃCOWE")
    print("=" * 70)

    extrap_results.sort(key=lambda x: x[0])
    for a_inf, err, idx in extrap_results:
        print(f"  Rodzina {idx}: α*_∞ = {a_inf:.8f} ± {err:.3e}")

    chk("P1  Znaleziono ≥2 trwałe rodziny", len(stable) >= 2,
        f"Znaleziono: {len(stable)}")

    if len(extrap_results) >= 2:
        a1, e1, _ = extrap_results[0]
        a2, e2, _ = extrap_results[1]
        S     = a1 + a2
        dS    = S - S_FORMULA
        sigma = np.sqrt(e1**2 + e2**2)
        print(f"\n  α*₁(∞) = {a1:.8f} ± {e1:.3e}")
        print(f"  α*₂(∞) = {a2:.8f} ± {e2:.3e}")
        print(f"\n  S(∞)   = {S:.8f} ± {sigma:.3e}")
        print(f"  S_form = {S_FORMULA:.8f}  (2π−11/10)")
        print(f"  S−Sf   = {dS:+.5e}  ({dS/S_FORMULA*1e6:+.1f} ppm)")

        chk("P2  Ekstrapolacja zbieżna (σ_α < 0.003)",
            e1 < 0.003 and e2 < 0.003,
            f"σ₁={e1:.4f}, σ₂={e2:.4f}")
        chk("P3  S(∞) ≠ 2π−11/10: odchylenie > 500 ppm",
            abs(dS/S_FORMULA) * 1e6 > 500,
            f"{abs(dS/S_FORMULA)*1e6:.0f} ppm")
        chk("P4  |α*₁(∞) − α*₂(∞)| wyznaczone z precyzją 0.01",
            e1 < 0.01 and e2 < 0.01, f"e1={e1:.4f}, e2={e2:.4f}")

        # Kandydaci algebraiczni
        print("\n  Kandydaci algebraiczni dla S(∞):")
        cands = [
            ('2π − 11/10',  2*np.pi - 1.1),
            ('2π − 1',      2*np.pi - 1.0),
            ('π + 2',       np.pi + 2.0),
            ('5.10',        5.10),
            ('5.12',        5.12),
            ('5.14',        5.14),
            ('5.15',        5.15),
        ]
        for name, val in cands:
            dev = (S - val) / val * 1e6
            print(f"    {name:16s} = {val:.8f}  ({dev:+.0f} ppm)")

        # Kandydaci dla D = α*₂ − α*₁
        D = a2 - a1
        print(f"\n  D = α*₂ − α*₁ = {D:.8f}")
        for name, val in [('3/10=0.3', 0.3), ('π/10', np.pi/10),
                           ('1/4', 0.25), ('1/3', 1/3)]:
            print(f"    {name:12s} = {val:.8f}  (Δ={D-val:+.6f})")

    else:
        chk("P2  ekstrapolacja", len(extrap_results) >= 2)
        chk("P3  obalenie S", len(extrap_results) >= 2)
        chk("P4  precyzja", len(extrap_results) >= 2)

    print(f"\n  WYNIK: {PASS_CNT}/{PASS_CNT+FAIL_CNT}  PASS")
    print("=" * 70)
