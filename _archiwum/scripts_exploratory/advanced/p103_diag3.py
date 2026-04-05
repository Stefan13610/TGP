#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p103_diag3.py — Weryfikacja czy K*_2=0.171272 z p102 naprawde ma g_U=0.

Pytania:
1. g_U przy (psi=1.2419, K=0.010281) — znane K*_1, powinno byc ~0
2. g_U przy (psi=4.10, K=0.171272) — K*_2 z p102, czy g=0?
3. Poszukiwanie g_U=0 dla wiekszego K (K*_2 moze byc gdzie indziej)
4. Profil phi(r) dla K*_1 i K*_2-kandydatow: czy fizyczne?
"""
import sys, io, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA_K=8.5616; GAM=1.0; V1=GAM/3.0-GAM/4.0; LAM=5.501357e-06; AGAM=0.040; R_MAX=80.0
N_EVAL=600; RTOL=1e-8; ATOL=1e-10

def dV_n6(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.0)**5
def V_n6(phi):  return GAM/3*phi**3 - GAM/4*phi**4 + LAM/6*(phi-1.0)**6

def ode_rhs(r, y):
    phi, dphi = y; phi = max(phi, 1e-10); kfac = 1.0+ALPHA_K/phi
    return [dphi, dV_n6(phi)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def solve_ode(psi, K):
    dphi0 = -K/AGAM**2
    r_ev = AGAM*(R_MAX/AGAM)**np.linspace(0,1,N_EVAL)
    sol = solve_ivp(ode_rhs, [AGAM,R_MAX], [psi,dphi0],
                    method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
    if sol.t[-1] < R_MAX*0.99: return None
    return sol

def energy(sol):
    r = sol.t; phi = np.maximum(sol.y[0],1e-10); dphi = sol.y[1]
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))): return np.nan
    kfac = 1.0+ALPHA_K/phi
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vpot = np.array([V_n6(float(p)) for p in phi])
    Ep = 4*np.pi*np.trapezoid((Vpot-V1)*r**2, r)
    return Ek+Ep

# ── Test 1: K*_1 cross-check ──────────────────────────────────────────────────
print("=== Test 1: K*_1 = 0.010281, psi = 1.2419 ===")
K1 = 0.010281; psi1 = 1.2419
sol1 = solve_ode(psi1, K1)
if sol1:
    E1 = energy(sol1)
    g1 = E1/(4*np.pi*K1)-1
    phi_end1 = sol1.y[0,-1]
    print(f"  phi(r_max) = {phi_end1:.6f}  (powinno byc ~1.0)")
    print(f"  E = {E1:.6f}")
    print(f"  4*pi*K = {4*np.pi*K1:.6f}")
    print(f"  g_U = {g1:.6f}  (powinno byc ~0)")
    # profil
    phi = sol1.y[0]
    print(f"  phi_min = {phi.min():.4f}, phi_max = {phi.max():.4f}")
    print(f"  Czy monoton. malejacy? {np.all(np.diff(phi) <= 0.01)}")

# ── Test 2: K*_2 z p102 ───────────────────────────────────────────────────────
print("\n=== Test 2: K*_2 = 0.171272, psi = 4.10 (z p102) ===")
K2 = 0.171272; psi2 = 4.10
sol2 = solve_ode(psi2, K2)
if sol2:
    E2 = energy(sol2)
    g2 = E2/(4*np.pi*K2)-1
    phi_end2 = sol2.y[0,-1]
    print(f"  phi(r_max) = {phi_end2:.6f}  (powinno byc ~1.0)")
    print(f"  E = {E2:.6f}")
    print(f"  4*pi*K = {4*np.pi*K2:.6f}")
    print(f"  g_U = {g2:.6f}  (z p102 powinno byc ~0)")
    phi = sol2.y[0]
    print(f"  phi_min = {phi.min():.4f}, phi_max = {phi.max():.4f}")
    mono = np.all(np.diff(phi) <= 0.01)
    print(f"  Czy monoton. malejacy? {mono}")
    if not mono:
        non_mono = np.where(np.diff(phi) > 0.01)[0]
        print(f"  Pierwsze naruszenie przy i={non_mono[0]}, r={sol2.t[non_mono[0]]:.2f}")

# ── Test 3: K*_2 przy roznych psi (dla K=0.171272) ───────────────────────────
print("\n=== Test 3: g_U vs psi dla K=0.171272 ===")
print(f"  {'psi':>6}  {'phi_end':>9}  {'E':>12}  {'g_U':>10}  mono")
for psi_t in [1.24, 2.0, 3.0, 3.5, 3.9, 4.0, 4.1, 4.2, 4.5, 5.0, 6.0]:
    sol_t = solve_ode(psi_t, K2)
    if sol_t:
        E_t = energy(sol_t)
        g_t = E_t/(4*np.pi*K2)-1 if np.isfinite(E_t) else np.nan
        phi_t = sol_t.y[0]
        mono_t = np.all(np.diff(phi_t) <= 0.01)
        g_s = f"{g_t:.4f}" if np.isfinite(g_t) else "  nan "
        print(f"  {psi_t:>6.2f}  {sol_t.y[0,-1]:>9.4f}  {E_t:>12.4f}  {g_s:>10}  {mono_t}")
    else:
        print(f"  {psi_t:>6.2f}  FAIL")

# ── Test 4: Skan g_U(K) dla KAZDEGO K, psi_0 z [1.0, 7.0] ───────────────────
print("\n=== Test 4: Skan g_U(K) dla K in [0.001, 0.30] ===")
print("Dla kazdego K: znajdz WSZYSTKIE zera phi_at_rmax=1 w [1.001,7.0]")
print(f"  {'K':>8}  {'psi_0':>7}  {'phi_end':>9}  {'g_U':>10}  mono")

def phi_end(psi, K):
    sol = solve_ode(psi, K)
    return sol.y[0,-1] if sol else np.nan

K_scan = np.concatenate([
    np.linspace(0.001, 0.020, 20),   # K*_1 region
    np.linspace(0.020, 0.200, 30),   # K*_2 region
])

for K_t in K_scan:
    # Skanuj psi w [1.001, 7.0] (szerokie okno)
    psis = np.linspace(1.001, 7.0, 80)
    Fv = [phi_end(p, K_t) - 1.0 for p in psis]
    zeros = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_end(p,K_t)-1.0, psis[i], psis[i+1], xtol=1e-6, maxiter=30)
                sol_z = solve_ode(pz, K_t)
                if sol_z:
                    E_z = energy(sol_z)
                    g_z = E_z/(4*np.pi*K_t)-1 if np.isfinite(E_z) else np.nan
                    mono_z = np.all(np.diff(sol_z.y[0]) <= 0.01)
                    zeros.append((pz, sol_z.y[0,-1], E_z, g_z, mono_z))
            except: pass
    for pz, pe, E_z, g_z, mono_z in zeros:
        g_s = f"{g_z:.4f}" if np.isfinite(g_z) else "  nan "
        print(f"  {K_t:>8.5f}  {pz:>7.4f}  {pe:>9.4f}  {g_s:>10}  {mono_z}")

print("\n=== Koniec diagnostyki ===")
