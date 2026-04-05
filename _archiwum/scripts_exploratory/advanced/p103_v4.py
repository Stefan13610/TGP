#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p103_v4.py  --  TGP v1
======================
Skan wykladnika n w V_mod = gamma/3*phi^3 - gamma/4*phi^4 + lam/n*(phi-1)^n
POPRAWIONA METODOLOGIA (po diagnozie p103_diag3):

BLAD p101/p102: PSI2=[3.50,7.00] wykluczalo fizyczna galaz K*_2!
  Prawdziwy K*_2 (gal. B) ma psi_0~2.77, K~0.100  =>  K*_2/K*_1~9.73
  p102 findowal galaz C (psi~4.1, g_U=3.38 != 0) — niesamoistna.

POPRAWKA:
  - PSI2=[1.80, 3.20] z filtrem mono (phi monoton. malej.)
  - K*_2 skan: K in [0.060, 0.150], step~0.003
  - Pierwszy zero F(psi) w [1.80,3.20] z mono=True => galaz B

Data: 2026-03-26
"""
import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA_K = 8.5616; GAM = 1.0; V1 = GAM/3.0 - GAM/4.0
LAM = 5.501357e-06; AGAM = 0.040; R_MAX = 80.0
R21_TARGET = 206.77; K1_REF = 0.010281; K2_TRUE = 0.1000

PSI1_LO, PSI1_HI = 1.10, 1.80
PSI2_LO, PSI2_HI = 1.80, 3.20

K1_LO = K1_REF * 0.70; K1_HI = K1_REF * 1.50
K2_LO = 0.060;  K2_HI = 0.150

N_PSI1 = 50; N_PSI2 = 100; N_K1 = 20; N_K2 = 30
N_EVAL = 600; RTOL = 1e-8; ATOL = 1e-10

N_VALUES = [4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20]


def dV_mod_n(phi, n):
    delta = phi - 1.0
    try:
        corr = LAM * float(delta)**(n - 1)
        if not np.isfinite(corr): corr = np.sign(delta)**(n-1) * 1e15
    except: corr = np.sign(delta)**(n-1) * 1e15
    return GAM*phi**2 - GAM*phi**3 + corr

def V_mod_n(phi, n):
    delta = phi - 1.0
    try:
        stab = (LAM/n)*float(delta)**n
        if not np.isfinite(stab): stab = 1e15
    except: stab = 1e15
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + stab

def ode_rhs_n(r, y, n):
    phi, dphi = y; phi = max(phi, 1e-10); kfac = 1.0 + ALPHA_K/phi
    return [dphi, dV_mod_n(phi,n)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def solve_full(psi, K, n, r_max=R_MAX):
    dphi0 = -K/AGAM**2
    r_ev = AGAM*(r_max/AGAM)**np.linspace(0,1,N_EVAL)
    for method in ('DOP853','Radau'):
        try:
            sol = solve_ivp(lambda r,y: ode_rhs_n(r,y,n), [AGAM,r_max], [psi,dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] < r_max*0.99:
                if method=='DOP853': continue
                return None
            return sol
        except:
            if method=='DOP853': continue
            return None
    return None

def phi_at_rmax(psi, K, n):
    sol = solve_full(psi, K, n)
    return float(sol.y[0,-1]) if sol else np.nan

def is_monotone(sol, tol=0.01):
    phi = sol.y[0]
    return bool(np.all(np.diff(phi) <= tol))

def energy_soliton(sol, K, n):
    r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))): return np.nan
    kfac = 1.0 + ALPHA_K/phi
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vpot = np.clip(np.array([V_mod_n(float(p),n) for p in phi]), -1e12, 1e12)
    Ep = 4*np.pi*np.trapezoid((Vpot - V1)*r**2, r)
    E = Ek+Ep
    return E if np.isfinite(E) else np.nan

def find_psi_zero_mono(K, n, psi_lo, psi_hi, n_psi):
    """Znajdz pierwszy MONOTONNY zero F(psi) w [psi_lo, psi_hi]."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv = [phi_at_rmax(p, K, n) - 1.0 for p in psis]
    for i in range(len(Fv)-1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1) and fi*fi1 < 0): continue
        try:
            pz = brentq(lambda p: phi_at_rmax(p,K,n)-1.0, psis[i], psis[i+1], xtol=1e-7, maxiter=50)
            sol = solve_full(pz, K, n)
            if sol and is_monotone(sol):
                return pz, sol
        except: pass
    return np.nan, None

def g_U(K, n, psi_lo, psi_hi, n_psi):
    pz, sol = find_psi_zero_mono(K, n, psi_lo, psi_hi, n_psi)
    if np.isnan(pz) or sol is None: return np.nan, np.nan
    E = energy_soliton(sol, K, n)
    g = E/(4*np.pi*K)-1.0 if np.isfinite(E) else np.nan
    return g, pz

def find_Kstar(n, K_lo, K_hi, n_K, psi_lo, psi_hi, n_psi):
    Ks = np.linspace(K_lo, K_hi, n_K)
    gv, psv = [], []
    for K in Ks:
        g, pz = g_U(K, n, psi_lo, psi_hi, n_psi)
        gv.append(g); psv.append(pz)
    for i in range(len(gv)-1):
        gi, gi1 = gv[i], gv[i+1]
        if np.isfinite(gi) and np.isfinite(gi1) and gi*gi1 < 0:
            _last = [np.nan]
            def gfunc(K):
                g2, pz2 = g_U(K, n, psi_lo, psi_hi, n_psi)
                _last[0] = pz2
                return g2 if np.isfinite(g2) else np.nan
            try:
                ga = gfunc(Ks[i]); gb = gfunc(Ks[i+1])
                if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                    Kstar = brentq(gfunc, Ks[i], Ks[i+1], xtol=1e-5, maxiter=40)
                    psi_star = _last[0]
                    if not np.isfinite(psi_star):
                        psi_star, _ = find_psi_zero_mono(Kstar, n, psi_lo, psi_hi, n_psi)
                    sol_star = solve_full(psi_star, Kstar, n) if np.isfinite(psi_star) else None
                    E_star = energy_soliton(sol_star, Kstar, n) if sol_star else np.nan
                    g_star = E_star/(4*np.pi*Kstar)-1.0 if np.isfinite(E_star) else np.nan
                    return Kstar, psi_star, g_star
            except:
                Kstar = Ks[i] - gv[i]*(Ks[i+1]-Ks[i])/(gv[i+1]-gv[i])
                psi_star, sol_star = find_psi_zero_mono(Kstar, n, psi_lo, psi_hi, n_psi)
                E_star = energy_soliton(sol_star, Kstar, n) if sol_star else np.nan
                g_star = E_star/(4*np.pi*Kstar)-1.0 if np.isfinite(E_star) else np.nan
                return Kstar, psi_star, g_star
    return np.nan, np.nan, np.nan


t0 = time.time()
print("="*72)
print("p103_v4: Skan wykladnika n  (POPRAWIONA METODOLOGIA)")
print(f"V_mod = g/3*phi^3 - g/4*phi^4 + lam/n*(phi-1)^n")
print(f"lambda={LAM:.3e}, a_Gam={AGAM}, r_max={R_MAX}")
print(f"K*_2 PRAWDZIWE (gal.B): K~0.100, psi_0~2.77, K*_2/K*_1~9.73")
print(f"K*_1 okno psi: [{PSI1_LO},{PSI1_HI}], K:[{K1_LO:.4f},{K1_HI:.4f}]")
print(f"K*_2 okno psi: [{PSI2_LO},{PSI2_HI}] (MONO filter), K:[{K2_LO},{K2_HI}]")
print("="*72)

results = []

for n in N_VALUES:
    tn = time.time()
    print(f"\n=== n = {n} ===", flush=True)

    print(f"  K*_1...", end=' ', flush=True)
    K1, psi1, g1 = find_Kstar(n, K1_LO, K1_HI, N_K1, PSI1_LO, PSI1_HI, N_PSI1)
    if np.isfinite(K1):
        print(f"K*_1={K1:.6f}  psi_0={psi1:.4f}  g={g1:.2e}")
    else:
        print("BRAK")

    print(f"  K*_2...", end=' ', flush=True)
    K2, psi2, g2 = find_Kstar(n, K2_LO, K2_HI, N_K2, PSI2_LO, PSI2_HI, N_PSI2)
    if np.isfinite(K2):
        print(f"K*_2={K2:.6f}  psi_0={psi2:.4f}  g={g2:.2e}")
    else:
        print("BRAK")

    r21 = K2/K1 if (np.isfinite(K1) and np.isfinite(K2)) else np.nan
    if np.isfinite(r21):
        frac = r21/R21_TARGET
        print(f"  K*_2/K*_1 = {r21:.4f}  ({frac:.1%} celu)")
    dt = time.time()-tn
    print(f"  [czas: {dt:.0f}s]", flush=True)
    results.append((n, K1, psi1, g1, K2, psi2, g2, r21))


print("\n"+"="*72)
print("TABELA WYNIKOW")
print(f"  {'n':>4}  {'K*_1':>10}  {'psi1':>6}  {'K*_2':>10}  {'psi2':>6}  {'r21':>8}  g1      g2")
print("  "+"-"*70)
for n, K1, psi1, g1, K2, psi2, g2, r21 in results:
    K1s = f"{K1:.6f}" if np.isfinite(K1) else "   --   "
    K2s = f"{K2:.6f}" if np.isfinite(K2) else "   --   "
    p1s = f"{psi1:.3f}" if np.isfinite(psi1) else "  -- "
    p2s = f"{psi2:.3f}" if np.isfinite(psi2) else "  -- "
    r21s = f"{r21:.4f}" if np.isfinite(r21) else "   --  "
    g1s = f"{g1:.2e}" if np.isfinite(g1) else "  nan "
    g2s = f"{g2:.2e}" if np.isfinite(g2) else "  nan "
    print(f"  {n:>4}  {K1s}  {p1s}  {K2s}  {p2s}  {r21s}  {g1s}  {g2s}")

print("\n"+"="*72)
valid = [(n,r) for n,K1,p1,g1,K2,p2,g2,r in results if np.isfinite(r)]
if valid:
    r_vals = [r for _,r in valid]
    n_max, r_max = max(valid, key=lambda x:x[1])
    n_min, r_min = min(valid, key=lambda x:x[1])
    print(f"K*_2/K*_1 zakres: [{r_min:.4f}, {r_max:.4f}]")
    print(f"Srednia: {np.mean(r_vals):.4f},  Std: {np.std(r_vals):.4f}")
    mono = all(r_vals[i+1]-r_vals[i]>-0.5 for i in range(len(r_vals)-1))
    print(f"Trend: {'monot.rosn' if mono else 'niemonot.'}")
    const = np.std(r_vals)/np.mean(r_vals) < 0.05
    print(f"Stalost: {'STALY (odch<5%)' if const else 'ZMIENNY'}")
    print(f"\nWniosek: K*_2/K*_1 = {np.mean(r_vals):.3f} (+/-{np.std(r_vals):.3f})")
    if const:
        print(f"  => K*_2/K*_1 NIEZALEZY od n (lambda-term zaniedbywalny)")
        print(f"  => Zmiana wykladnika n NIE pomaga osiagnac celu r21={R21_TARGET}")

print("\nPASS/FAIL:")
r_all = [r for _,_,_,_,_,_,_,r in results if np.isfinite(r)]
K1_n6 = next((K1 for n,K1,*_ in results if n==6 and np.isfinite(K1)), np.nan)
P1 = np.isfinite(K1_n6)
P2 = P1 and abs(K1_n6/K1_REF-1.0)<0.05
P3 = len(r_all)>=5
P4 = len(r_all)>=2 and np.std(r_all)/np.mean(r_all)<0.10
P5 = P1 and len(r_all)>0 and abs(r_all[0]-9.73)<1.0  # r21~9.73
P6 = any(r>=R21_TARGET for r in r_all)
for i,(P,d) in enumerate([(P1,"K*_1 znaleziony dla n=6"),
    (P2,f"K*_1(n=6) stabilny (ref={K1_REF})"),
    (P3,"K*_2 znaleziony >=5 n"),
    (P4,"K*_2/K*_1 stalY (std<10%)"),
    (P5,"K*_2/K*_1~9.73 (prawdziwa wartosc)"),
    (P6,f"Cel r21={R21_TARGET} osiagniety")],1):
    print(f"  P{i}: {'PASS' if P else 'FAIL'} -- {d}")

total = time.time()-t0
print(f"\nWynik: {sum([P1,P2,P3,P4,P5,P6])}/6 PASS  (czas: {total/60:.1f} min)")
