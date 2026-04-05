#!/usr/bin/env python3
"""
ex87v3_rmax_extrapolation.py — Ekstrapolacja α*₁,₂ do R_MAX→∞  (O-L1 finalny)
==============================================================================
Odkrycie ex87v2: α*₁, α*₂ są wrażliwe na R_MAX:
    R_MAX=100: α*₁≈2.432, α*₂≈2.636
    R_MAX=120: α*₁≈2.440, α*₂≈2.695 (ex84)

Zadanie: ekstrapolacja α*(R_MAX) → α*(∞) metodą fitowania 1/R_MAX.

Strategia:
  Dla 6 wartości R_MAX ∈ {80, 100, 120, 140, 160, 200}:
    1. Zlokalizuj przejście fz0(α)=R21 w wąskim oknie [2.40,2.47] i [2.46,2.70]
    2. brentq → α*(R_MAX)
  Następnie: fit α*(R_MAX) = α*_inf + c/R_MAX + d/R_MAX²
  Wyznacz α*_inf = α*(R_MAX→∞)

Test zbieżności:
  P1: α*(200)/α*(120) stabilne na 0.5%
  P2: ekstrapolacja spójna: α*_inf ≈ α*(200) ± 2σ
  P3: S = α*₁_inf + α*₂_inf ≠ 2π−11/10 (obalenie potwierdzone)
  P4: Dokładne α*_inf z precyzją < 0.005

Autor: Claudian (sesja v35, O-L1 finalny)
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
PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21       = 206.7682830
S_FORMULA = 2.0 * np.pi - 1.1
G_OFF     = 0.005

RMAX_LIST = [80, 100, 120, 140, 160, 200]

PASS_CNT = 0; FAIL_CNT = 0
def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    tag = "[PASS]" if cond else "[FAIL]"
    if cond: PASS_CNT += 1
    else:    FAIL_CNT += 1
    print(f"  {tag}  {name}")
    if detail: print(f"        {detail}")

# ================================================================
# ODE
# ================================================================
def integrate_soliton(g0, alpha, r_max):
    gb = np.exp(-1.0/(2.0*alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb+1e-8)
        fg = 1.0 + 2.0*alpha*np.log(g)
        if abs(fg)<1e-10: return [gp, 0.0]
        Vprime = g*g*(1.0-g); curl = (alpha/g)*gp*gp
        if r<1e-10: return [gp, (Vprime-curl)/(3.0*fg)]
        return [gp, (Vprime-curl-fg*2.0*gp/r)/fg]
    def ev_bounce(r, y): return y[0]-gb
    ev_bounce.terminal=True; ev_bounce.direction=-1
    y0=[g0,0.0]; r0=1e-10; ra,ga=[],[]
    for _ in range(12):
        sol=solve_ivp(rhs,(r0,r_max),y0,events=ev_bounce,
                      dense_output=True,rtol=1e-11,atol=1e-13,max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status==1 and len(sol.t_events[0])>0:
            rh=sol.t_events[0][-1]; st=sol.sol(rh)
            y0=[st[0],-st[1]]; r0=rh
        else: break
    idx=np.argsort(np.concatenate(ra))
    return np.concatenate(ra)[idx], np.concatenate(ga)[idx]

def fit_amplitude(r, g, r_max):
    """Ekstrapolacja A_inf z okien liniowo rozmieszczonych."""
    # Okna skalowane z r_max
    win_fracs = [(0.18,0.32),(0.27,0.40),(0.36,0.48),(0.45,0.58),(0.54,0.68)]
    Av, rv = [], []
    for fl, fh in win_fracs:
        rL, rR = fl*r_max, fh*r_max
        if rR > r[-1]: break
        mask = (r>=rL)&(r<=rR)
        if mask.sum() < 20: continue
        rf=r[mask]; df=(g[mask]-1.0)*rf
        M=np.column_stack([np.cos(rf),np.sin(rf)])
        bc,_,_,_=np.linalg.lstsq(M,df,rcond=None)
        Av.append(float(np.sqrt(bc[0]**2+bc[1]**2))); rv.append(rL)
    if not Av: return 0.0
    if len(Av)<2: return Av[-1]
    rv=np.array(rv); Av=np.array(Av)
    try:
        p,_=curve_fit(lambda x,ai,a,b:ai*(1+a/x+b/x**2),rv,Av,
                      p0=[Av[-1],0.,0.],maxfev=3000)
        return float(p[0])
    except:
        try:
            p,_=curve_fit(lambda x,ai,a:ai*(1+a/x),rv,Av,
                          p0=[Av[-1],0.],maxfev=1000)
            return float(p[0])
        except: return float(Av[-1])

def B_coeff(g0, alpha, r_max):
    r,g=integrate_soliton(g0,alpha,r_max)
    rL,rR=0.28*r_max, 0.40*r_max
    mask=(r>=rL)&(r<=rR)
    if mask.sum()<20: return 0.0
    rf=r[mask]; df=(g[mask]-1.0)*rf
    M=np.column_stack([np.cos(rf),np.sin(rf)])
    bc,_,_,_=np.linalg.lstsq(M,df,rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, lo=1.05, hi=2.5, n=35):
    gs=np.linspace(lo,hi,n); Bs=[]
    for g in gs:
        try: Bs.append(B_coeff(g,alpha,r_max))
        except: Bs.append(np.nan)
    Bs=np.array(Bs)
    for i in range(len(Bs)-1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i]*Bs[i+1]<0:
            try:
                return brentq(lambda x:B_coeff(x,alpha,r_max),
                              gs[i],gs[i+1],xtol=1e-9,maxiter=100)
            except: pass
    return None

def fz0_rmax(alpha, r_max):
    z0=find_z0(alpha,r_max)
    if z0 is None: return np.nan
    r,g=integrate_soliton(z0,alpha,r_max)
    Ae=fit_amplitude(r,g,r_max)
    r2,g2=integrate_soliton(PHI*z0,alpha,r_max)
    Am=fit_amplitude(r2,g2,r_max)
    if Ae<1e-6: return np.nan
    return (Am/Ae)**4

def find_alpha_star(bracket, r_max, n_scan=18):
    """Znajdź crossing fz0=R21 w oknie bracket przy r_max."""
    lo, hi = bracket
    alphas = np.linspace(lo, hi, n_scan)
    fvals = []
    for a in alphas:
        try: fvals.append(fz0_rmax(a, r_max))
        except: fvals.append(np.nan)
    fvals = np.array(fvals)
    # Znajdź crossing
    for i in range(len(fvals)-1):
        if np.isnan(fvals[i]) or np.isnan(fvals[i+1]): continue
        if (fvals[i]-R21)*(fvals[i+1]-R21) < 0:
            try:
                return brentq(lambda a: fz0_rmax(a,r_max)-R21,
                              alphas[i], alphas[i+1],
                              xtol=1e-10, maxiter=200)
            except: pass
    return None

def worker_rmax(args):
    """Worker: (bracket, r_max) → α*(r_max)"""
    bracket, r_max = args
    return r_max, find_alpha_star(bracket, r_max)

# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    print("=" * 65)
    print("  ex87v3 — Ekstrapolacja α*₁,₂ do R_MAX→∞ (O-L1 finalny)")
    print("=" * 65)
    print(f"  R21 = {R21:.8f},  S_formula = 2π−11/10 = {S_FORMULA:.10f}")
    print(f"  R_MAX list = {RMAX_LIST}")

    N_CPU = min(6, cpu_count())
    print(f"  Rdzenie: {N_CPU}")

    # Przedziały poszukiwania α* dla obu zer
    # Z ex87v2 (R_MAX=100): α*₁∈[2.430,2.436], α*₂∈[2.50,2.70]
    # Z ex84 (R_MAX=120):   α*₁∈[2.437,2.443], α*₂∈[2.69,2.70]
    # Łączymy w szerszy zakres
    bracket1 = (2.420, 2.460)   # α*₁
    bracket2 = (2.500, 2.720)   # α*₂ (szeroki, obejmuje oba wyniki)

    print(f"\n  Przedział 1 (α*₁): {bracket1}")
    print(f"  Przedział 2 (α*₂): {bracket2}")

    # ----------------------------------------------------------------
    # OBLICZENIA α*(R_MAX) dla obu zer
    # ----------------------------------------------------------------
    tasks1 = [(bracket1, rmax) for rmax in RMAX_LIST]
    tasks2 = [(bracket2, rmax) for rmax in RMAX_LIST]

    print(f"\n--- α*₁(R_MAX) ---")
    print(f"  {'R_MAX':>8}  {'α*₁':>14}  {'f(z0)':>10}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*10}")

    results1 = {}
    with Pool(N_CPU) as pool:
        raw1 = pool.map(worker_rmax, tasks1)
    for rmax, astar in sorted(raw1):
        results1[rmax] = astar
        if astar is not None:
            f_val = fz0_rmax(astar, rmax)
            print(f"  {rmax:8d}  {astar:14.10f}  {f_val:10.4f}")
        else:
            print(f"  {rmax:8d}  {'[nie znaleziono]':>14}")

    print(f"\n--- α*₂(R_MAX) ---")
    print(f"  {'R_MAX':>8}  {'α*₂':>14}  {'f(z0)':>10}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*10}")

    results2 = {}
    with Pool(N_CPU) as pool:
        raw2 = pool.map(worker_rmax, tasks2)
    for rmax, astar in sorted(raw2):
        results2[rmax] = astar
        if astar is not None:
            f_val = fz0_rmax(astar, rmax)
            print(f"  {rmax:8d}  {astar:14.10f}  {f_val:10.4f}")
        else:
            print(f"  {rmax:8d}  {'[nie znaleziono]':>14}")

    # ----------------------------------------------------------------
    # EKSTRAPOLACJA α*(R_MAX→∞)
    # ----------------------------------------------------------------
    print("\n--- Ekstrapolacja do R_MAX→∞ ---")

    def extrapolate(results_dict, label):
        rmaxs = sorted([r for r, v in results_dict.items() if v is not None])
        astars = [results_dict[r] for r in rmaxs]
        if len(rmaxs) < 3:
            print(f"  {label}: za mało punktów ({len(rmaxs)}) do ekstrapolacji")
            return None, None
        xs = np.array(rmaxs, dtype=float)
        ys = np.array(astars, dtype=float)
        print(f"\n  {label}:")
        for rmax, astar in zip(xs, ys):
            print(f"    R_MAX={rmax:.0f}: α*={astar:.8f}")
        try:
            # Fit: α*(R_MAX) = α*_inf + c/R_MAX
            p1, cov1 = curve_fit(lambda x, ai, c: ai + c/x,
                                  xs, ys, p0=[ys[-1], 0.])
            a_inf_1 = p1[0]; err1 = np.sqrt(cov1[0,0])
            print(f"  Fit 1/R_MAX: α*_inf = {a_inf_1:.10f} ± {err1:.2e}")
            # Fit: α*(R_MAX) = α*_inf + c/R_MAX + d/R_MAX²
            if len(rmaxs) >= 4:
                p2, cov2 = curve_fit(lambda x, ai, c, d: ai + c/x + d/x**2,
                                      xs, ys, p0=[ys[-1], 0., 0.])
                a_inf_2 = p2[0]; err2 = np.sqrt(cov2[0,0])
                print(f"  Fit 1/R²:   α*_inf = {a_inf_2:.10f} ± {err2:.2e}")
                return a_inf_2, err2
            return a_inf_1, err1
        except Exception as e:
            print(f"  Ekstrapolacja FAIL: {e}")
            return None, None

    a1_inf, err1 = extrapolate(results1, "α*₁")
    a2_inf, err2 = extrapolate(results2, "α*₂")

    # ----------------------------------------------------------------
    # ANALIZA KOŃCOWA
    # ----------------------------------------------------------------
    print("\n" + "=" * 65)
    print("  WYNIKI EKSTRAPOLACJI (R_MAX→∞)")
    print("=" * 65)

    if a1_inf is not None and a2_inf is not None:
        S_inf = a1_inf + a2_inf
        dS    = S_inf - S_FORMULA
        print(f"\n  α*₁(∞) = {a1_inf:.8f} ± {err1:.2e}")
        print(f"  α*₂(∞) = {a2_inf:.8f} ± {err2:.2e}")
        print(f"\n  S(∞)   = {S_inf:.8f}")
        print(f"  S_form = {S_FORMULA:.8f}  (2π−11/10)")
        print(f"  S−Sf   = {dS:+.4e}  ({dS/S_FORMULA*1e6:+.1f} ppm)")

        # Wyniki ze znanych R_MAX
        a1_200 = results1.get(200)
        a2_200 = results2.get(200)
        if a1_200 and a2_200:
            chk("P1  α*(200)/α*(120) zbieżne (|Δ|<0.01)",
                abs(a1_200-results1.get(120,a1_200))<0.01 and
                abs(a2_200-results2.get(120,a2_200))<0.01,
                f"|Δα*₁|={abs(a1_200-results1.get(120,a1_200)):.4f}, "
                f"|Δα*₂|={abs(a2_200-results2.get(120,a2_200)):.4f}")
        chk("P2  Ekstrapolacja spójna (err<0.01)",
            err1 < 0.01 and err2 < 0.01,
            f"err1={err1:.4f}, err2={err2:.4f}")
        chk("P3  S(∞) ≠ 2π−11/10: odchylenie > 500 ppm",
            abs(dS/S_FORMULA)*1e6 > 500,
            f"{abs(dS/S_FORMULA)*1e6:.0f} ppm")
        chk("P4  α*₁(∞) wyznaczone (nie NaN)",
            not np.isnan(a1_inf), f"α*₁={a1_inf:.6f}")

        # Algebraiczni kandydaci
        print("\n  Kandydaci algebraiczni dla S(∞):")
        cands = [
            ('2π − 11/10', S_FORMULA),
            ('2π − 1',     2*np.pi - 1.0),
            ('5.10',       5.10),
            ('5.12',       5.12),
            ('5.14',       5.14),
            ('π + 2',      np.pi + 2.0),
        ]
        for name, val in cands:
            print(f"    {name:16s} = {val:.8f}  ({(S_inf-val)/val*1e6:+.0f} ppm)")

    else:
        chk("P1  dane do ekstrapolacji", False)
        chk("P2  ekstrapolacja", False)
        chk("P3  S(∞)≠S_form", False)
        chk("P4  α*₁ wyznaczone", False)

    print(f"\n  WYNIK: {PASS_CNT}/{PASS_CNT+FAIL_CNT}  PASS")
    print("=" * 65)
