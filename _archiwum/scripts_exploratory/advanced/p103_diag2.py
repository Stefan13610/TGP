#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p103_diag2.py — Diagnostyka spurious zero w p103.

Cel: wyjasnij dlaczego psi_0=3.5001 dla wszystkich n.
1. Pokaz F(psi) = phi_at_rmax(psi,K,n) - 1 dla n=6 przy K=K*_2, K*_2+0.01, K*_2-0.01
2. Sprawdz czy PSI2_LO=3.90 naprawia problem (dla n=6 cross-check z p102)
3. Pokaz F(psi) dla n=14 — co sie dzieje gdy lambda term jest ogromny
"""
import sys, io, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA_K = 8.5616; GAM = 1.0; V1 = GAM/3.0 - GAM/4.0
LAM = 5.501357e-06; AGAM = 0.040; R_MAX = 80.0
K2_REF = 0.171272; K1_REF = 0.010281
N_EVAL = 600; RTOL = 1e-8; ATOL = 1e-10

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

def phi_at_rmax(psi, K, n, r_max=R_MAX):
    dphi0 = -K/AGAM**2
    r_ev = AGAM*(r_max/AGAM)**np.linspace(0,1,N_EVAL)
    for method in ('DOP853','Radau'):
        try:
            sol = solve_ivp(lambda r,y: ode_rhs_n(r,y,n), [AGAM,r_max], [psi,dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] < r_max*0.99:
                if method == 'DOP853': continue
                return np.nan
            return float(sol.y[0,-1])
        except:
            if method == 'DOP853': continue
            return np.nan
    return np.nan

def energy_soliton(psi, K, n):
    dphi0 = -K/AGAM**2
    r_ev = AGAM*(R_MAX/AGAM)**np.linspace(0,1,N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs_n(r,y,n), [AGAM,R_MAX], [psi,dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
        if sol.t[-1] < R_MAX*0.99: return np.nan
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))): return np.nan
        kfac = 1.0 + ALPHA_K/phi
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vpot = np.clip(np.array([V_mod_n(float(p),n) for p in phi]), -1e12, 1e12)
        Ep = 4*np.pi*np.trapezoid((Vpot - V1)*r**2, r)
        E = Ek+Ep
        return E if np.isfinite(E) else np.nan
    except: return np.nan

# ── Test 1: F(psi) dla n=6 przy roznych K ─────────────────────────────────────
print("=== Test 1: F(psi) = phi_at_rmax(psi,K,n=6) - 1 ===")
print("psi range: [3.30, 4.50] step=0.05")
psis = np.arange(3.30, 4.55, 0.05)

for K_label, K_val in [("K*_2=0.171272", K2_REF),
                        ("K*_2-0.01=0.161", K2_REF-0.01),
                        ("K*_2+0.01=0.181", K2_REF+0.01),
                        ("K*_2-0.05=0.121", K2_REF-0.05)]:
    print(f"\n  K={K_label}:")
    Fv = [phi_at_rmax(p, K_val, 6) - 1.0 for p in psis]
    for p, f in zip(psis, Fv):
        marker = " <-- ZERO" if abs(f) < 0.05 else ""
        sign = "+" if f > 0 else "-"
        fs = f"{f:+.4f}" if np.isfinite(f) else "  nan  "
        print(f"    psi={p:.2f}  F={fs}{marker}")

# ── Test 2: Czy PSI2_LO=3.90 naprawia problem? ────────────────────────────────
print("\n\n=== Test 2: find_psi_zero z PSI2_LO=3.90 dla n=6 ===")
def find_psi_zero_fixed(K, n, psi_lo=3.90, psi_hi=7.00, n_psi=60):
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv = [phi_at_rmax(p, K, n) - 1.0 for p in psis]
    for i in range(len(Fv)-1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,n)-1.0, psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                E = energy_soliton(pz, K, n)
                g = E/(4*np.pi*K)-1.0 if np.isfinite(E) else np.nan
                return pz, g
            except: pass
    return np.nan, np.nan

K_scan = np.linspace(0.12, 0.23, 12)
print(f"  {'K':>8}  {'psi_0':>8}  {'g_U':>10}  status")
for K in K_scan:
    pz, g = find_psi_zero_fixed(K, 6)
    pz_s = f"{pz:.4f}" if np.isfinite(pz) else "  nan  "
    g_s  = f"{g:.4f}" if np.isfinite(g) else "  nan  "
    ok   = "OK" if (np.isfinite(g) and abs(g)<0.5) else ("bad_g" if np.isfinite(g) else "BRAK")
    print(f"  {K:>8.5f}  {pz_s:>8}  {g_s:>10}  {ok}")

# ── Test 3: g_U scan z PSI2_LO=3.90 — szukamy sign change ────────────────────
print("\n\n=== Test 3: g_U(K) z PSI2_LO=3.90 — poszukiwanie K*_2 dla n=6 ===")
K_fine = np.linspace(0.043, 0.685, 50)
g_fine = []
psi_fine = []
for K in K_fine:
    pz, g = find_psi_zero_fixed(K, 6)
    g_fine.append(g)
    psi_fine.append(pz)

# Pokaz tylko non-nan
valid = [(K, g, pz) for K, g, pz in zip(K_fine, g_fine, psi_fine) if np.isfinite(g)]
print(f"  Znaleziono {len(valid)} valid g_U z {len(K_fine)} K wartosci")
if valid:
    print("  Pierwsze 8 valid:")
    for K, g, pz in valid[:8]:
        print(f"    K={K:.5f}  psi_0={pz:.4f}  g={g:.5f}")
    print("  Ostatnie 8 valid:")
    for K, g, pz in valid[-8:]:
        print(f"    K={K:.5f}  psi_0={pz:.4f}  g={g:.5f}")

    # Szukaj sign change
    gv = [g for _,g,_ in valid]
    Kv = [K for K,_,_ in valid]
    for i in range(len(gv)-1):
        if gv[i]*gv[i+1] < 0:
            print(f"\n  ZMIANA ZNAKU: K=[{Kv[i]:.5f},{Kv[i+1]:.5f}]  g=[{gv[i]:.5f},{gv[i+1]:.5f}]")
            # Wyznacz K*_2
            try:
                def gfunc(K):
                    pz2, g2 = find_psi_zero_fixed(K, 6)
                    return g2 if np.isfinite(g2) else np.nan
                ga = gfunc(Kv[i]); gb = gfunc(Kv[i+1])
                if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                    Kstar = brentq(gfunc, Kv[i], Kv[i+1], xtol=1e-5, maxiter=30)
                    pstar, gstar = find_psi_zero_fixed(Kstar, 6)
                    print(f"  ==> K*_2 = {Kstar:.6f}  psi_0 = {pstar:.4f}  g = {gstar:.6f}")
                    print(f"  ==> K*_2/K*_1 = {Kstar/K1_REF:.4f}")
            except Exception as e:
                print(f"  brentq error: {e}")

# ── Test 4: F(psi) dla n=14 przy K=K*_2 ─────────────────────────────────────
print("\n\n=== Test 4: F(psi) dla n=14 przy K=K*_2=0.171 (czy spurious zniklo?) ===")
psis14 = np.arange(3.30, 7.10, 0.20)
print(f"  K={K2_REF}  n=14:")
for p in psis14:
    f = phi_at_rmax(p, K2_REF, 14)
    fs = f"{f-1.0:+.4f}" if np.isfinite(f) else "  nan  "
    print(f"    psi={p:.2f}  F={fs}")

# ── Test 5: Gdzie jest K*_2 dla n=14? ────────────────────────────────────────
print("\n\n=== Test 5: g_U(K) dla n=14 z PSI2=[3.90,9.00] — gdzie K*_2? ===")
def find_psi_zero_wide(K, n, psi_lo=3.90, psi_hi=9.00, n_psi=60):
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv = [phi_at_rmax(p, K, n) - 1.0 for p in psis]
    for i in range(len(Fv)-1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,n)-1.0, psis[i], psis[i+1], xtol=1e-6, maxiter=50)
                E = energy_soliton(pz, K, n)
                g = E/(4*np.pi*K)-1.0 if np.isfinite(E) else np.nan
                return pz, g
            except: pass
    return np.nan, np.nan

K_scan14 = np.concatenate([np.linspace(0.04, 0.69, 20), np.linspace(0.7, 5.0, 20)])
print(f"  {'K':>8}  {'psi_0':>8}  {'g_U':>10}")
for K in K_scan14:
    pz, g = find_psi_zero_wide(K, 14)
    if np.isfinite(pz):
        g_s = f"{g:.4f}" if np.isfinite(g) else "nan"
        print(f"  {K:>8.5f}  {pz:.4f}  {g_s}")

print("\n=== Diagnostyka zakonczona ===")
