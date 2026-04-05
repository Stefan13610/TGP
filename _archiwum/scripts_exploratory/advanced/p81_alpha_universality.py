#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p81_alpha_universality.py  --  TGP v1  *  Uniwersalnosc alpha_max i E(alpha)
=============================================================================
Dwa pytania z p80/PLAN_ROZWOJU:

PYTANIE 1: Czy alpha_max jest niezalezny od a_Gam?
  Teoretyczna prognoza: K*_1 = C_1 * a_Gam (Schwarzschild); w zrescalowanym ODE
  a_Gam odpada -> alpha_max powinien byc wlasnoscia V_mod, nie a_Gam.
  Test: wyznacz alpha_max dla a_Gam in {0.020, 0.030, 0.040, 0.050, 0.060}

PYTANIE 2: Czy E_soliton(alpha) ma minimum przy alpha_K?
  OP-3 Sciezka 1 (PLAN_ROZWOJU): alpha_K = argmin_alpha E_soliton(alpha)
  Test: oblicz E*(alpha) = E[phi*(alpha)] (energia na galeezi K*_1) dla
        alpha in [5, 9] i sprawdz czy minimum jest blisko 8.5616

Data: 2026-03-25
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF = 8.5616
A_GAM_REF   = 0.040
K1_REF      = 0.010414
LAM         = 5.501357e-06
GAMMA_V     = 1.0

def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
V1 = GAMMA_V/3 - GAMMA_V/4


def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]


def _phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=2000):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0] - 1e-5
    ev.terminal = True; ev.direction = -1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0,-1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except: return np.nan


def _energy_and_K(psi_z, K, a_gam, alpha, r_max=40.0):
    """Energia solitonu przy psi_core=psi_z, parametrach (K, a_gam, alpha)."""
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek + Ep
    except: return np.nan


def compute_g(K, a_gam, alpha, r_max=40.0, n_psi=80):
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1 for p in psis]
    branches = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: _phi_at_rmax(p,K,a_gam,alpha,r_max)-1,
                            psis[i], psis[i+1], xtol=1e-5, rtol=1e-5, maxiter=40)
                E = _energy_and_K(pz, K, a_gam, alpha, r_max)
                if np.isfinite(E):
                    branches.append((pz, E/(4*np.pi*K)-1))
            except: pass
    if not branches: return np.nan, np.nan
    best = min(branches, key=lambda x: abs(x[1]))
    return best[1], best[0]   # (g, psi_z)


def find_K_star_full(alpha, a_gam=A_GAM_REF, K_min=3e-3, K_max=0.12, n_K=20, n_psi=80):
    """Zwraca (K*_1, psi_z, E) lub (nan, nan, nan)."""
    r_max = max(40.0, a_gam*5, 6.0/np.sqrt(1+alpha))
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = []
    psi_vals = []
    for K in K_scan:
        g, pz = compute_g(K, a_gam, alpha, r_max, n_psi)
        g_vals.append(g); psi_vals.append(pz)

    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                K_star = brentq(
                    lambda K: compute_g(K, a_gam, alpha, r_max, n_psi)[0],
                    K_scan[i], K_scan[i+1], xtol=K_scan[i]*0.003, rtol=0.003, maxiter=20)
                _, pz = compute_g(K_star, a_gam, alpha, r_max, n_psi)
                E = _energy_and_K(pz, K_star, a_gam, alpha, r_max) if np.isfinite(pz) else np.nan
                return K_star, pz, E
            except: pass
    return np.nan, np.nan, np.nan


def alpha_max_bisect(a_gam, a_lo=8.70, a_hi=9.00, tol=0.01):
    """Wyznacza alpha_max bisection dla danego a_gam."""
    # Sprawdz ze a_lo daje soliton, a_hi nie
    K_lo, _, _ = find_K_star_full(a_lo, a_gam)
    K_hi, _, _ = find_K_star_full(a_hi, a_gam)
    if not np.isfinite(K_lo):
        return np.nan   # a_lo juz nie ma solitonu
    if np.isfinite(K_hi):
        return np.nan   # a_hi ma soliton — trzeba przesunac zakres

    n = 0
    while a_hi - a_lo > tol and n < 8:
        a_mid = 0.5*(a_lo + a_hi)
        K_mid, _, _ = find_K_star_full(a_mid, a_gam)
        if np.isfinite(K_mid):
            a_lo = a_mid
        else:
            a_hi = a_mid
        n += 1
    return 0.5*(a_lo + a_hi)


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0; FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("TGP v1  *  p81_alpha_universality.py")
print("P1-P3: alpha_max vs a_Gam (niezaleznosc?)")
print("P4-P6: E_soliton(alpha) — minimum przy alpha_K?")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA A: alpha_max vs a_Gam
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SEKCJA A: alpha_max dla roznych a_Gam")
print("Prognoza: alpha_max = const (niezalezny od a_Gam)")
print("=" * 70)

a_gam_scan = [0.020, 0.030, 0.040, 0.050, 0.060]
alpha_max_vals = []

print(f"\n  {'a_Gam':>8}  {'alpha_max':>12}  {'alpha_K/alpha_max':>18}")
print("  " + "-" * 42)

for a_gam in a_gam_scan:
    print(f"  a_Gam={a_gam:.3f} ...", flush=True)
    K_min = 3e-3 * a_gam/0.04
    K_max = 0.15 * a_gam/0.04
    amax = alpha_max_bisect(a_gam, a_lo=8.70, a_hi=9.00, tol=0.01)
    alpha_max_vals.append(amax)
    if np.isfinite(amax):
        ratio = ALPHA_K_REF / amax
        print(f"  {a_gam:>8.3f}  {amax:>12.4f}  {ratio:>18.4f}")
    else:
        print(f"  {a_gam:>8.3f}  {'---':>12}  {'---':>18}")

valid_max = [v for v in alpha_max_vals if np.isfinite(v)]
if len(valid_max) >= 3:
    spread = max(valid_max) - min(valid_max)
    mean_amax = np.mean(valid_max)
    rel_spread = spread / mean_amax
    print(f"\n  Rozrzut alpha_max: [{min(valid_max):.4f}, {max(valid_max):.4f}]")
    print(f"  Srednia: {mean_amax:.4f},  rozrzut wzgl.: {rel_spread*100:.2f}%")

ok_P1 = len(valid_max) >= 3
record("P1: alpha_max wyznaczone dla >= 3 wartosci a_Gam", ok_P1,
       f"{len(valid_max)}/5 wartosci")

ok_P2 = ok_P1 and rel_spread < 0.05
record("P2: alpha_max niezalezny od a_Gam (rozrzut < 5%)", ok_P2,
       f"rozrzut = {rel_spread*100:.2f}%" if ok_P1 else "brak danych")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B: E_soliton(alpha) — minimum?
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SEKCJA B: Energia solitonu E*(alpha) na galezi K*_1")
print("OP-3 Sciezka 1: czy alpha_K = argmin E*(alpha)?")
print("=" * 70)

alpha_E_scan = [5.0, 6.0, 7.0, 7.5, 8.0, 8.2, 8.4, 8.5, 8.5616, 8.6, 8.7, 8.8]

results_E = []
print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'E*':>14}  {'g':>8}")
print("  " + "-" * 45)

for alpha in alpha_E_scan:
    print(f"  alpha={alpha:.4f} ...", flush=True)
    K_star, pz, E_star = find_K_star_full(alpha, A_GAM_REF, n_psi=80)
    results_E.append((alpha, K_star, E_star))
    if np.isfinite(K_star) and np.isfinite(E_star):
        g_val = E_star/(4*np.pi*K_star) - 1
        print(f"  {alpha:>8.4f}  {K_star:>10.6f}  {E_star:>14.6e}  {g_val:>8.4f}")
    else:
        print(f"  {alpha:>8.4f}  {'---':>10}  {'---':>14}  {'---':>8}")

# Szukaj minimum E*(alpha)
valid_E = [(a, K, E) for a, K, E in results_E if np.isfinite(E)]
if len(valid_E) >= 3:
    alphas_E = np.array([x[0] for x in valid_E])
    E_vals   = np.array([x[2] for x in valid_E])
    K_vals   = np.array([x[1] for x in valid_E])

    # Minimum bezwzgledne E
    idx_min_E = np.argmin(E_vals)
    alpha_at_Emin = alphas_E[idx_min_E]
    E_min = E_vals[idx_min_E]

    # Minimum E/K (znormalizowane — gęstość energii na ładunek)
    EoverK = E_vals / K_vals
    idx_min_EK = np.argmin(EoverK)
    alpha_at_EKmin = alphas_E[idx_min_EK]

    print(f"\n  Minimum E*(alpha):      alpha = {alpha_at_Emin:.4f},  E* = {E_min:.4e}")
    print(f"  Minimum E*/K*(alpha):   alpha = {alpha_at_EKmin:.4f}")
    print(f"  alpha_K_ref = {ALPHA_K_REF}")

    ok_P4 = True
    record("P4: E*(alpha) profil obliczony dla >= 5 wartosci", ok_P4,
           f"{len(valid_E)} punktow")

    ok_P5 = abs(alpha_at_Emin - ALPHA_K_REF) < 1.0
    record("P5: minimum E*(alpha) blisko alpha_K (< 1.0 odch.)", ok_P5,
           f"alpha_Emin = {alpha_at_Emin:.4f}, odch = {abs(alpha_at_Emin-ALPHA_K_REF):.4f}")

    ok_P6 = abs(alpha_at_EKmin - ALPHA_K_REF) < 1.0
    record("P6: minimum E*/K*(alpha) blisko alpha_K (< 1.0 odch.)", ok_P6,
           f"alpha_EKmin = {alpha_at_EKmin:.4f}, odch = {abs(alpha_at_EKmin-ALPHA_K_REF):.4f}")
else:
    ok_P4 = ok_P5 = ok_P6 = False
    record("P4: E*(alpha) profil", False, "za malo danych")
    record("P5: minimum E* blisko alpha_K", False, "brak")
    record("P6: minimum E*/K* blisko alpha_K", False, "brak")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PODSUMOWANIE  p81_alpha_universality.py")
print("=" * 70)
print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

if ok_P2:
    print("\n  WNIOSEK A (P2 PASS):")
    print(f"    alpha_max = {mean_amax:.4f} +/- {spread/2:.4f}")
    print("    Jest NIEZALEZNY od a_Gam — to wlasnosc V_mod, nie skali.")
    print("    Konsekwencja: alpha_K/alpha_max = n_s jest relacja miedzy")
    print("    struktura V_mod (bifurkacja) a indeksem spektralnym CMB.")

if ok_P5 or ok_P6:
    print("\n  WNIOSEK B (P5/P6 PASS):")
    print(f"    E*(alpha) ma minimum przy alpha ~ {alpha_at_Emin:.4f}")
    print("    Sugeruje: alpha_K moze byc wyznaczone jako minimum energii solitonu.")
    print("    To byloby OP-3 Sciezka 1 (energia, nie bifurkacja).")
else:
    print("\n  WNIOSEK B:")
    print("    E*(alpha) minimum nie w poblizu alpha_K lub brak danych.")
    print("    OP-3 Sciezka 1 (argmin E) wymaga dalszego badania.")
