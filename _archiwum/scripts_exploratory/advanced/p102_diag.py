#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostyka p102 — dlaczego K*_1 nie jest znajdowane."""
import sys, io, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA_K = 8.5616
GAM     = 1.0
V1      = GAM/3.0 - GAM/4.0
LAM_PHYS = 5.501357e-06
K1_REF   = 0.010281   # z p101
PSI1_REF = 1.2419     # z p101
AGAM     = 0.040
R_MAX    = 80.0
N_EVAL   = 600

def dV_mod_lam(phi, lam):
    return GAM*phi**2 - GAM*phi**3 + lam*(phi - 1.0)**5

def V_mod_lam(phi, lam):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + lam/6.0*(phi - 1.0)**6

def ode_rhs_lam(r, y, lam):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    phi2 = max(phi, 1e-10)**2
    return [dphi,
            dV_mod_lam(phi, lam)/kfac
            + ALPHA_K*dphi**2/(2.0*phi2*kfac)
            - 2.0/r*dphi]

def phi_at_rmax(psi, K, lam, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs_lam(r, y, lam),
                        [AGAM, r_max], [psi, dphi0],
                        method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan

def energy_soliton(psi, K, lam, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs_lam(r, y, lam),
                        [AGAM, r_max], [psi, dphi0],
                        method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(dphi)):
            return np.nan
        kfac = 1.0 + ALPHA_K / phi
        Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * kfac * r**2, r)
        Ep   = 4*np.pi * np.trapezoid((V_mod_lam(phi, lam) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

print("=== Test 1: phi_at_rmax z wartosciami referencyjnymi (p101) ===")
val = phi_at_rmax(PSI1_REF, K1_REF, LAM_PHYS)
print(f"phi_at_rmax(psi={PSI1_REF}, K={K1_REF}, lam={LAM_PHYS:.3e}) = {val:.6f}  (oczek. ~1.0)")

E = energy_soliton(PSI1_REF, K1_REF, LAM_PHYS)
g = E / (4*np.pi*K1_REF) - 1.0
print(f"energy = {E:.6f},  g_U = {g:.6f}  (oczek. ~0.0)")

print("\n=== Test 2: skan psi w [1.001,1.80] dla K=K1_REF ===")
psis = np.linspace(1.001, 1.80, 15)
Fv = [phi_at_rmax(p, K1_REF, LAM_PHYS) - 1.0 for p in psis]
for p, f in zip(psis, Fv):
    print(f"  psi={p:.4f}  F={f:.6f}")

print("\n=== Test 3: g_U vs K przy stalym lambda_phys ===")
Ks = np.linspace(0.006, 0.016, 9)
for K in Ks:
    # znajdz psi_0
    psis2 = np.linspace(1.001, 1.80, 20)
    Fv2 = [phi_at_rmax(p, K, LAM_PHYS) - 1.0 for p in psis2]
    psi_zero = np.nan
    for i in range(len(Fv2)-1):
        if np.isfinite(Fv2[i]) and np.isfinite(Fv2[i+1]) and Fv2[i]*Fv2[i+1] < 0:
            try:
                psi_zero = brentq(lambda p: phi_at_rmax(p, K, LAM_PHYS) - 1.0,
                                   psis2[i], psis2[i+1], xtol=1e-6, maxiter=30)
            except Exception:
                pass
            break
    if np.isfinite(psi_zero):
        E2 = energy_soliton(psi_zero, K, LAM_PHYS)
        g2 = E2/(4*np.pi*K) - 1.0 if np.isfinite(E2) else np.nan
        print(f"  K={K:.5f}: psi_0={psi_zero:.4f}, E={E2:.5f}, g_U={g2:.4f}")
    else:
        Fmin = np.nanmin([f for f in Fv2 if np.isfinite(f)]) if any(np.isfinite(f) for f in Fv2) else np.nan
        Fmax = np.nanmax([f for f in Fv2 if np.isfinite(f)]) if any(np.isfinite(f) for f in Fv2) else np.nan
        print(f"  K={K:.5f}: BRAK zera psi, F in [{Fmin:.4f},{Fmax:.4f}]")

print("\n=== Test 4: phi(r) profil dla (psi=1.2419, K=0.01028, lam=lam_phys) ===")
K_test = K1_REF; psi_test = PSI1_REF; lam_test = LAM_PHYS
dphi0 = -K_test / AGAM**2
r_ev = AGAM * (R_MAX / AGAM)**np.linspace(0, 1, 200)
sol2 = solve_ivp(lambda r, y: ode_rhs_lam(r, y, lam_test),
                 [AGAM, R_MAX], [psi_test, dphi0],
                 method='DOP853', rtol=1e-8, atol=1e-10, t_eval=r_ev)
phi_profile = sol2.y[0]
print(f"phi(r_min={r_ev[0]:.4f}) = {phi_profile[0]:.6f}")
print(f"phi(r_max={r_ev[-1]:.2f}) = {phi_profile[-1]:.6f}")
print(f"phi min = {np.min(phi_profile):.6f}, phi max = {np.max(phi_profile):.6f}")
print(f"phi monotone decreasing: {np.all(np.diff(phi_profile) <= 0.001)}")
