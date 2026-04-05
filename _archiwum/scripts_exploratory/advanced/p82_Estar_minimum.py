#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p82_Estar_minimum.py  --  TGP v1  *  Precyzyjny argmin E*(alpha)
=================================================================
Cel: Czy minimum energii solitonu E*(alpha) = E_sol[K*_1(alpha), alpha]
     wypada dokladnie przy alpha_K = 8.5616?

Kontekst (p81):
  - Grube skanowanie (krok 0.1) daje minimum E* przy alpha=8.60
  - Roznica E*(8.5616) vs E*(8.60) = 0.14% — zbyt mala by rozstrzygac
  - Potrzebny krok 0.01 w przedziale [8.50, 8.70]

OP-3 Sciezka 1:
  alpha_K = argmin E*(alpha) ?
  Jesli tak -> alpha_K jest wyznaczony energetycznie bez r_21_PDG

Testy:
  P1: E*(alpha) obliczone dla >= 10 wartosci w [8.50, 8.70]
  P2: Minimum E*(alpha) znalezione (parabola fit lub skan)
  P3: |argmin E* - alpha_K| < 0.10 (wskazuje Sciezke 1)
  P4: |argmin E* - alpha_K| < 0.02 (mocne potwierdzenie)
  P5: Minimum K*_1(alpha) wspolbiezne z minimum E* (spodziewane z g=0)

Data: 2026-03-25
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import warnings

warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF = 8.5616
A_GAM       = 0.040
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


def compute_g_and_psi(K, a_gam, alpha, r_max=40.0, n_psi=100):
    """Zwraca (g_best, psi_best)."""
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1 for p in psis]
    branches = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: _phi_at_rmax(p,K,a_gam,alpha,r_max)-1,
                            psis[i], psis[i+1], xtol=1e-6, rtol=1e-6, maxiter=50)
                # Energia
                dphi0 = -K/a_gam**2
                r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, 4000)
                sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                                [pz, dphi0], method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_ev)
                r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
                Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
                Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
                E = Ek + Ep
                if np.isfinite(E):
                    branches.append((pz, E/(4*np.pi*K)-1, E))
            except: pass
    if not branches: return np.nan, np.nan, np.nan
    best = min(branches, key=lambda x: abs(x[1]))
    return best[1], best[0], best[2]   # (g, psi_z, E)


def find_K_star_E(alpha, a_gam=A_GAM, K_min=4e-3, K_max=0.12, n_K=25, n_psi=100):
    """Zwraca (K*_1, E*) lub (nan, nan)."""
    r_max = 40.0
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_list = []
    for K in K_scan:
        g, _, _ = compute_g_and_psi(K, a_gam, alpha, r_max, n_psi)
        g_list.append(g)

    K_star = np.nan
    for i in range(len(g_list)-1):
        gi, gj = g_list[i], g_list[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                def g_func(K):
                    gv, _, _ = compute_g_and_psi(K, a_gam, alpha, r_max, n_psi)
                    return gv if np.isfinite(gv) else 0.0
                K_star = brentq(g_func, K_scan[i], K_scan[i+1],
                                xtol=K_scan[i]*0.001, rtol=0.001, maxiter=25)
                break
            except: pass

    if not np.isfinite(K_star):
        return np.nan, np.nan

    # Energia przy K_star
    g, pz, E = compute_g_and_psi(K_star, a_gam, alpha, r_max, n_psi)
    return K_star, E


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0; FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 68)
print("TGP v1  *  p82_Estar_minimum.py  --  argmin E*(alpha)")
print("Precyzyjny skan alpha w [8.50, 8.70] z krokiem 0.01")
print("=" * 68)
print(f"\n  alpha_K_ref = {ALPHA_K_REF}")
print(f"  a_Gam = {A_GAM},  K_min=4e-3, K_max=0.12, n_K=25, n_psi=100")

# ─────────────────────────────────────────────────────────────────────────────
# Gesty skan alpha w [8.50, 8.70]
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("SKAN alpha w [8.50, 8.70] z krokiem ~0.01")
print("=" * 68)

alpha_scan = np.arange(8.50, 8.71, 0.01)

results = []
print(f"\n  {'alpha':>7}  {'K*_1':>10}  {'E*':>14}  {'E*/K*':>10}")
print("  " + "-" * 47)

for alpha in alpha_scan:
    K_star, E_star = find_K_star_E(alpha)
    results.append((alpha, K_star, E_star))
    if np.isfinite(K_star) and np.isfinite(E_star):
        EoK = E_star/K_star if K_star > 0 else np.nan
        print(f"  {alpha:>7.4f}  {K_star:>10.6f}  {E_star:>14.8f}  {EoK:>10.6f}")
    else:
        print(f"  {alpha:>7.4f}  {'---':>10}  {'---':>14}  {'---':>10}")

# ─────────────────────────────────────────────────────────────────────────────
# Analiza minimum
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("ANALIZA MINIMUM E*(alpha)")
print("=" * 68)

valid = [(a, K, E) for a, K, E in results if np.isfinite(E) and np.isfinite(K)]
ok_P1 = len(valid) >= 10
record("P1: E*(alpha) obliczone dla >= 10 wartosci", ok_P1, f"{len(valid)} punktow")

if len(valid) >= 5:
    a_arr = np.array([x[0] for x in valid])
    K_arr = np.array([x[1] for x in valid])
    E_arr = np.array([x[2] for x in valid])

    # Minimum bezposrednie
    idx_min = np.argmin(E_arr)
    alpha_Emin_direct = a_arr[idx_min]
    E_min_direct = E_arr[idx_min]
    print(f"\n  Minimum bezposrednie: alpha = {alpha_Emin_direct:.4f},  E* = {E_min_direct:.8f}")

    # Fit paraboliczny (jesli wystarczy punktow)
    alpha_Emin_fit = np.nan
    if len(valid) >= 5:
        # Fit na 5 punktach wokol minimum
        i0 = max(1, idx_min-2)
        i1 = min(len(valid)-2, idx_min+2)
        a_fit = a_arr[i0:i1+1]
        E_fit = E_arr[i0:i1+1]
        if len(a_fit) >= 3:
            coeffs = np.polyfit(a_fit, E_fit, 2)
            if coeffs[0] > 0:
                alpha_Emin_fit = -coeffs[1] / (2*coeffs[0])
                E_min_fit = np.polyval(coeffs, alpha_Emin_fit)
                print(f"  Fit paraboliczny:     alpha = {alpha_Emin_fit:.5f},  E* = {E_min_fit:.8f}")

    # Minimum K*_1
    idx_Kmin = np.argmin(K_arr)
    alpha_Kmin = a_arr[idx_Kmin]
    print(f"  Minimum K*_1(alpha):  alpha = {alpha_Kmin:.4f}")
    print(f"\n  alpha_K_ref = {ALPHA_K_REF}")

    best_alpha_min = alpha_Emin_fit if np.isfinite(alpha_Emin_fit) else alpha_Emin_direct
    odch_P3 = abs(best_alpha_min - ALPHA_K_REF)
    odch_P4 = abs(best_alpha_min - ALPHA_K_REF)

    print(f"\n  Odchylenie argmin E* od alpha_K: {odch_P3:.5f}")

    ok_P2 = np.isfinite(best_alpha_min)
    record("P2: Minimum E*(alpha) znalezione", ok_P2,
           f"alpha_min = {best_alpha_min:.5f}")

    ok_P3 = ok_P2 and odch_P3 < 0.10
    record("P3: |argmin E* - alpha_K| < 0.10 (Sciezka 1 sugerowana)", ok_P3,
           f"odch = {odch_P3:.5f}")

    ok_P4 = ok_P2 and odch_P4 < 0.02
    record("P4: |argmin E* - alpha_K| < 0.02 (mocne potwierdzenie)", ok_P4,
           f"odch = {odch_P4:.5f}")

    ok_P5 = abs(alpha_Kmin - best_alpha_min) < 0.05
    record("P5: minimum K*_1 wspolbiezne z minimum E*", ok_P5,
           f"|alpha_Kmin - alpha_Emin| = {abs(alpha_Kmin - best_alpha_min):.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("PODSUMOWANIE  p82_Estar_minimum.py")
print("=" * 68)
print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

if ok_P4:
    print(f"""
  WNIOSEK P4 PASS: argmin E*(alpha) = {best_alpha_min:.5f}
  (odchylenie {odch_P4:.5f} = {odch_P4/ALPHA_K_REF*100:.3f}% od alpha_K)
  -> OP-3 Sciezka 1 POTWIERDZONA: alpha_K minimalizuje energie solitonu!
  -> alpha_K nie wymaga dopasowania do r_21 — wynika z warunku E -> min.
""")
elif ok_P3:
    print(f"""
  WNIOSEK P3 PASS: argmin E*(alpha) ~ {best_alpha_min:.5f}
  (odchylenie {odch_P3:.5f} od alpha_K={ALPHA_K_REF})
  -> OP-3 Sciezka 1 BLISKA (odch < 0.10), ale nie mocno potwierdzona.
  -> Minimum energii blisko alpha_K, moze byc numeryczna niedokladnosc.
""")
else:
    print("""
  WNIOSEK: argmin E*(alpha) odbiega od alpha_K o > 0.10.
  -> OP-3 Sciezka 1 (argmin E) nie potwierdzona w tym zakresie.
""")
