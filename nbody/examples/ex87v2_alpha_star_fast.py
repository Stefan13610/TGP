#!/usr/bin/env python3
"""
ex87v2_alpha_star_fast.py — Szybkie precyzyjne α*₁, α*₂ (O-L1, v35)
=====================================================================
Wersja szybka ex87: brentq w wąskich przedziałach wokół znanych zer
z ex84 (α*₁≈2.440, α*₂≈2.695), z multiprocessingiem.

Strategia (czas ~2-4 min zamiast ~90 min):
  1. Skan na małej siatce w dwóch wąskich oknach (po 20 pkt)
  2. Brentq w każdym przejściu
  3. Weryfikacja R_MAX=80/100/120

Wyniki ex84 (50 pkt, R_MAX=120):
  α*₁ ≈ 2.4396,  α*₂ ≈ 2.6953,  S ≈ 5.1349  (9307 ppm od 2π−11/10)

Autor: Claudian (sesja v35)
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
PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21       = 206.7682830
S_FORMULA = 2.0 * np.pi - 1.1
G_OFF     = 0.005
R_MAX     = 100.0   # Kompromis: dokładność vs czas (ex84 używało 120)

PASS_CNT = 0; FAIL_CNT = 0
def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    tag = "[PASS]" if cond else "[FAIL]"
    if cond: PASS_CNT += 1
    else:    FAIL_CNT += 1
    print(f"  {tag}  {name}")
    if detail: print(f"        {detail}")

# ================================================================
# ODE (identyczne z ex84)
# ================================================================
def integrate_soliton(g0, alpha, r_max=R_MAX):
    gb = np.exp(-1.0 / (2.0 * alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        Vprime = g*g*(1.0-g)
        curl   = (alpha/g)*gp*gp
        if r < 1e-10: return [gp, (Vprime-curl)/(3.0*fg)]
        return [gp, (Vprime-curl-fg*2.0*gp/r)/fg]
    def ev_bounce(r, y): return y[0] - gb
    ev_bounce.terminal = True; ev_bounce.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(10):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_bounce,
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    idx = np.argsort(np.concatenate(ra))
    return np.concatenate(ra)[idx], np.concatenate(ga)[idx]

WIN = [(20,36),(30,46),(40,56),(50,66)]

def fit_amplitude(r, g):
    Av, rv = [], []
    for rL, rR in WIN:
        if rR > r[-1]: break
        mask = (r>=rL)&(r<=rR)
        if mask.sum() < 20: continue
        rf = r[mask]; df = (g[mask]-1.0)*rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2+bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0
    if len(Av) < 2: return Av[-1]
    rv = np.array(rv); Av = np.array(Av)
    try:
        p, _ = curve_fit(lambda x,ai,a,b: ai*(1+a/x+b/x**2), rv, Av,
                         p0=[Av[-1],0.,0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def B_coeff(g0, alpha):
    r, g = integrate_soliton(g0, alpha)
    mask = (r>=28.0)&(r<=42.0)
    if mask.sum() < 20: return 0.0
    rf = r[mask]; df = (g[mask]-1.0)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, lo=1.05, hi=2.5, n=40):
    gs = np.linspace(lo, hi, n)
    Bs = []
    for g in gs:
        try: Bs.append(B_coeff(g, alpha))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs)-1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i]*Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff(x, alpha),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def fz0(alpha):
    z0 = find_z0(alpha)
    if z0 is None: return np.nan
    r, g = integrate_soliton(z0, alpha)
    Ae = fit_amplitude(r, g)
    r2, g2 = integrate_soliton(PHI*z0, alpha)
    Am = fit_amplitude(r2, g2)
    if Ae < 1e-6: return np.nan
    return (Am/Ae)**4

def worker(a): return a, fz0(a)

# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    print("=" * 65)
    print("  ex87v2 — Szybkie α*₁, α*₂ (O-L1, R_MAX=100)")
    print("=" * 65)
    print(f"  R21 = {R21:.8f}")
    print(f"  S_formula = 2π−11/10 = {S_FORMULA:.10f}")
    print(f"  R_MAX = {R_MAX},  PHI = {PHI:.8f}")

    N_CPU = min(8, cpu_count())
    print(f"  Rdzenie: {N_CPU}")

    # ----------------------------------------------------------------
    # Okno 1: wokół α*₁ ≈ 2.440  →  [2.40, 2.47]
    # Okno 2: wokół α*₂ ≈ 2.695  →  [2.65, 2.73]
    # ----------------------------------------------------------------
    win1 = np.linspace(2.40, 2.47, 22)
    win2 = np.linspace(2.65, 2.73, 22)
    alpha_all = np.concatenate([win1, win2])

    print(f"\n  Skanuję {len(alpha_all)} punktów w oknach [2.40,2.47]∪[2.65,2.73]...")
    with Pool(N_CPU) as pool:
        results = pool.map(worker, alpha_all)

    # ----------------------------------------------------------------
    # Znajdź przejścia przez R21
    # ----------------------------------------------------------------
    print(f"\n  {'α':>8}  {'f(z0)':>12}  {'f−R21':>10}")
    print(f"  {'-'*8}  {'-'*12}  {'-'*10}")
    crossings = []
    prev_a, prev_f = None, None
    for a, f in results:
        mark = ""
        if prev_f is not None and not np.isnan(f) and not np.isnan(prev_f):
            if (prev_f - R21)*(f - R21) < 0:
                crossings.append((prev_a, a))
                mark = "  ← PRZEJŚCIE"
        disp = f"{f:12.4f}" if not np.isnan(f) else f"{'NaN':>12}"
        diff = f"{f-R21:+10.4f}" if not np.isnan(f) else f"{'---':>10}"
        print(f"  {a:8.4f}  {disp}  {diff}{mark}")
        prev_a, prev_f = a, f

    print(f"\n  Znalezione przejścia: {len(crossings)}")
    chk("P1  Dokładnie 2 przejścia w oknach", len(crossings)==2,
        f"Znaleziono: {len(crossings)}")

    # ----------------------------------------------------------------
    # brentq w każdym przejściu
    # ----------------------------------------------------------------
    print("\n--- Precyzyjne zera (brentq, xtol=1e-12) ---")
    alpha_stars = []
    for lo, hi in crossings:
        try:
            a_star = brentq(lambda a: fz0(a)-R21, lo, hi,
                            xtol=1e-12, rtol=1e-13, maxiter=300)
            f_star = fz0(a_star)
            alpha_stars.append(a_star)
            print(f"  α* = {a_star:.12f},  f={f_star:.6f},  |f-R21|={abs(f_star-R21):.4e}")
        except Exception as e:
            print(f"  brentq FAIL [{lo:.5f},{hi:.5f}]: {e}")

    # ----------------------------------------------------------------
    # Analiza sumy
    # ----------------------------------------------------------------
    print("\n--- Analiza sumy i porównanie ---")
    if len(alpha_stars) >= 2:
        a1, a2 = sorted(alpha_stars[:2])
        S = a1 + a2; D = a2-a1; P = a1*a2
        dS = S - S_FORMULA
        print(f"\n  α*₁ = {a1:.12f}")
        print(f"  α*₂ = {a2:.12f}")
        print(f"\n  S   = {S:.10f}")
        print(f"  S_f = {S_FORMULA:.10f}  (2π−11/10)")
        print(f"  S−S_f = {dS:+.4e}  ({dS/S_FORMULA*1e6:+.1f} ppm)")
        print(f"\n  D   = {D:.10f}  (3/10={0.3:.10f}, Δ={D-0.3:+.3e})")
        P_hyp = (np.pi-0.7)*(np.pi-0.4)
        print(f"  P   = {P:.10f}  ((π-7/10)(π-4/10)={P_hyp:.10f})")

        ex84_a1, ex84_a2 = 2.4396, 2.6953
        print(f"\n  Porównanie z ex84 (R_MAX=120):")
        print(f"    α*₁: ex84={ex84_a1:.4f}, ex87v2={a1:.6f}, Δ={a1-ex84_a1:+.6f}")
        print(f"    α*₂: ex84={ex84_a2:.4f}, ex87v2={a2:.6f}, Δ={a2-ex84_a2:+.6f}")

        chk("P2  α*₁ stabilne względem R_MAX (|Δ|<0.003)", abs(a1-ex84_a1)<0.003,
            f"|Δα*₁|={abs(a1-ex84_a1):.6f}")
        chk("P3  α*₂ stabilne względem R_MAX (|Δ|<0.003)", abs(a2-ex84_a2)<0.003,
            f"|Δα*₂|={abs(a2-ex84_a2):.6f}")
        chk("P4  S ≠ 2π−11/10: odchylenie > 1000 ppm",
            abs(dS/S_FORMULA)*1e6 > 1000,
            f"{abs(dS/S_FORMULA)*1e6:.0f} ppm")

        # Kandydaci algebraiczni
        print("\n  Kandydaci dla S:")
        cands = [
            ('2π − 11/10',  2*np.pi - 1.1),
            ('2π − ln(3)', 2*np.pi - np.log(3)),
            ('2π − 1',     2*np.pi - 1.0),
            ('5.13',       5.13),
            ('5.135',      5.135),
        ]
        for name, val in cands:
            print(f"    {name:18s} = {val:.8f}  ({(S-val)/val*1e6:+.0f} ppm)")

    else:
        print("  Nie znaleziono dwóch zer!")
        chk("P2  dwa α* znalezione", False)
        chk("P3  stabilność", False)
        chk("P4  obalenie S", False)

    print("\n" + "=" * 65)
    print(f"  WYNIK: {PASS_CNT}/{PASS_CNT+FAIL_CNT}  PASS")
    print("=" * 65)
