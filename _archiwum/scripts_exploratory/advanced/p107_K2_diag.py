#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p107_K2_diag.py  --  TGP v1
Diagnostyka: precyzyjne K*₂ przez skan g_U(K) przy RTOL=1e-11, N_EVAL=2000.
Odtwarza stary task bb504odkp.
"""
import sys, io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA_K = 8.5616; GAM = 1.0; LAM = 5.501357e-06; AGAM = 0.040; R_MAX = 80.0
V1 = GAM/3 - GAM/4
N_EVAL = 2000; RTOL = 1e-11; ATOL = 1e-13

def V_mod(p):  return GAM/3*p**3 - GAM/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAM*p**2 - GAM*p**3 + LAM*(p-1)**5

def rhs(r, y):
    ph, dp = y; ph = max(ph, 1e-10); kf = 1 + ALPHA_K/ph
    return [dp, dV_mod(ph)/kf + ALPHA_K*dp**2/(2*ph**2*kf) - 2/r*dp]

def gU(K):
    def f(psi):
        dp0 = -K/AGAM**2
        rev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
        sol  = solve_ivp(rhs, [AGAM, R_MAX], [psi, dp0],
                         method='DOP853', rtol=RTOL, atol=ATOL, t_eval=rev)
        return sol.y[0,-1]-1.0 if sol.t[-1] > R_MAX*0.99 else float('nan')
    try:
        psi_star = brentq(f, 2.70, 2.83, xtol=1e-10, maxiter=100)
    except:
        return float('nan'), float('nan'), float('nan')
    dp0 = -K/AGAM**2
    rev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    sol  = solve_ivp(rhs, [AGAM, R_MAX], [psi_star, dp0],
                     method='DOP853', rtol=RTOL, atol=ATOL, t_eval=rev)
    if sol.t[-1] < R_MAX*0.99: return float('nan'), float('nan'), float('nan')
    r  = sol.t; ph = np.maximum(sol.y[0], 1e-10); dp = sol.y[1]
    kf = 1 + ALPHA_K/ph
    Vv = np.array([V_mod(float(q)) for q in ph])
    Ek = 4*np.pi*np.trapezoid(0.5*dp**2*kf*r**2, r)
    Ep = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
    E  = Ek + Ep
    return E/(4*np.pi*K) - 1.0, psi_star, E

print("Diagnostyka K*2 — RTOL=1e-11, N_EVAL=2000, brentq w [2.70,2.83]")
print("  " + "="*60)
print(f"  {'K':>8}  {'psi*':>12}  {'E':>12}  {'4piK':>12}  {'gU':>12}")
print("  " + "-"*60)
for K in [0.09900, 0.09950, 0.10000, 0.10028, 0.10040, 0.10050, 0.10100, 0.10150, 0.10200]:
    g, psi, E = gU(K)
    if not (isinstance(g, float) and g == g):
        print(f"  {K:.5f}  FAIL")
    else:
        print(f"  {K:.5f}  {psi:.8f}  {E:.8f}  {4*np.pi*K:.8f}  {g:.4e}")
print()
print("WNIOSEK: K*2 = K gdzie gU=0 (zmiana znaku)")
