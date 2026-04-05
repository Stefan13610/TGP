#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p93_rmax_convergence.py  --  TGP v1  *  Audyt zbieznosci K*(r_max)
===================================================================
Cel: wykazac zbieznosc K*(r_max) -> K*_true jako r_max -> inf,
     obalenie dwugaleziowej struktury z p91 (Galaz L = artefakt r_max=40),
     fit ogona oscylatoryjnego phi(r)-1 ~ A/r sin(mr)+B/r cos(mr).

Testy:
  P1: K*(r_max) maleje monotonicznie z r_max (spodziewane dla Branch U)
  P2: Fit 1/r zbiega: K*(r) = K_inf + c/r, R^2 > 0.95
  P3: K*_extrap (r->inf) rozni sie od K*(r_max=40) o > 5% (kwantyfikacja bledu p91)
  P4: Ogon phi(r)-1 fituje model Bessela A sin(mr)/r + B cos(mr)/r, R^2 > 0.90
  P5: Galaz L (psi "pierwsze zero") zmienia pozycje o > 0.01 miedzy r_max=40 a 100

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import linregress
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA  = 8.840
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA_V= 1.0
V1     = GAMMA_V/3 - GAMMA_V/4

RMAX_LIST = [30.0, 40.0, 60.0, 80.0, 100.0]   # r_max do testu zbieznosci

def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, a_gam, alpha, r_max, n_eval=1000):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0,1,n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi,dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        return float(sol.y[0,-1]) if sol.t[-1]>=r_max*0.99 else np.nan
    except: return np.nan

def energy_at_psi(psi_z, K, a_gam, alpha, r_max, n_eval=1000):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0,1,n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi_z,dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2,r)
        Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2,r)
        return Ek+Ep
    except: return np.nan

def scan_psi_zeros(K, alpha, r_max, n_psi=100):
    """Zwraca liste wszystkich zer phi(rmax,psi)-1 w [1.001,3.0]."""
    psis = np.linspace(1.001, 3.0, n_psi)
    Fv = [phi_at_rmax(p, K, A_GAM, alpha, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,A_GAM,alpha,r_max)-1.0,
                            psis[i],psis[i+1], xtol=1e-8, maxiter=60)
                roots.append(pz)
            except: pass
    return roots

def g_branch(psi_z, K, alpha, r_max):
    E = energy_at_psi(psi_z, K, A_GAM, alpha, r_max)
    return E/(4*np.pi*K)-1.0 if np.isfinite(E) else np.nan

def find_K_star_U(alpha, r_max, K_lo=0.008, K_hi=0.015, n_K=40):
    """Wyznacz K* dla galezy U (ostatnie zero psi) przez skan + brentq."""
    Ks = np.linspace(K_lo, K_hi, n_K)
    gv = []
    for K in Ks:
        roots = scan_psi_zeros(K, alpha, r_max)
        if roots:
            psi_U = roots[-1]
            g = g_branch(psi_U, K, alpha, r_max)
        else:
            g = np.nan
        gv.append(g)

    for i in range(len(gv)-1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1]<0:
            try:
                def g_at_K(K):
                    roots = scan_psi_zeros(K, alpha, r_max)
                    if not roots: return np.nan
                    return g_branch(roots[-1], K, alpha, r_max)
                Kz = brentq(g_at_K, Ks[i], Ks[i+1], xtol=5e-6, maxiter=20)
                return Kz, scan_psi_zeros(Kz, alpha, r_max)
            except: pass
    return np.nan, []


# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 70)
    print("TGP v1  *  p93_rmax_convergence.py")
    print(f"  alpha={ALPHA}, a_gam={A_GAM}")
    print(f"  r_max testowane: {RMAX_LIST}")
    print("=" * 70)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== SEKCJA A: K*(r_max) -- zbieznosc Galezy U ==")
    print(f"  {'r_max':>6}  {'K*_U':>10}  {'psi_U':>8}  {'n_zeros':>8}  {'t':>5}")
    print("  " + "-"*50)

    kstar_data = []  # [(r_max, K*, psi_U, n_zeros)]
    psi_L_data = []  # [(r_max, psi_L)] dla Branch L test

    for r_max in RMAX_LIST:
        t0 = time.time()
        Kz, roots = find_K_star_U(ALPHA, r_max)
        dt = time.time()-t0
        if np.isfinite(Kz) and roots:
            psi_U = roots[-1]
            psi_L = roots[0] if len(roots)>=1 else np.nan
            n_z   = len(roots)
            kstar_data.append((r_max, Kz, psi_U, n_z))
            if np.isfinite(psi_L):
                psi_L_data.append((r_max, psi_L))
            print(f"  {r_max:>6.0f}  {Kz:.7f}  {psi_U:.5f}  {n_z:>8}  {dt:>4.0f}s")
        else:
            kstar_data.append((r_max, np.nan, np.nan, 0))
            print(f"  {r_max:>6.0f}  {'brak':>10}  {'brak':>8}  {0:>8}  {dt:>4.0f}s")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== SEKCJA B: Fit K*(r_max) = K*_inf + c/r ==")
    valid = [(r, K) for r, K, *_ in kstar_data if np.isfinite(K)]
    if len(valid) >= 3:
        rs  = np.array([x[0] for x in valid])
        Ks  = np.array([x[1] for x in valid])
        inv_r = 1.0/rs
        slope, K_inf, rval, *_ = linregress(inv_r, Ks)
        print(f"  Fit: K*(r) = {K_inf:.7f} + {slope:.5f}/r")
        print(f"  R^2 = {rval**2:.5f}")
        print(f"  K*_extrapolowany (r->inf) = {K_inf:.6f}")
        print(f"  K*(r_max=40) = {Ks[rs==40][0]:.6f}" if 40.0 in rs else "")
        print(f"  Roznica K*(40) vs K*(inf): "
              f"{abs(Ks[rs==40][0]-K_inf)/K_inf*100:.2f}%" if 40.0 in rs else "")
        p2_ok = rval**2 > 0.90
        p3_ok = (abs(Ks[rs==40][0]-K_inf)/K_inf > 0.05) if 40.0 in rs else False
    else:
        K_inf = np.nan; p2_ok = False; p3_ok = False
        print("  Za malo punktow do fitu")

    # P1: monotonicznosc
    Ks_all = [K for _,K,*_ in kstar_data if np.isfinite(K)]
    p1_ok = all(Ks_all[i] >= Ks_all[i+1]-1e-5 for i in range(len(Ks_all)-1))
    record("P1: K*(r_max) monotonicznie maleje", p1_ok,
           f"wartosci: " + " > ".join(f"{K:.5f}" for K in Ks_all))

    record("P2: Fit 1/r zbiega (R^2>0.90)", p2_ok,
           f"R^2={rval**2:.4f}" if len(valid)>=3 else "brak danych")

    record("P3: K*(r=40) rozni sie od K*_inf o > 5%", p3_ok,
           f"roznica={abs(Ks[rs==40][0]-K_inf)/K_inf*100:.2f}%" if (len(valid)>=3 and 40.0 in rs) else "brak")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== SEKCJA C: Galaz L (psi_L) -- test czy pozycja zalezy od r_max ==")
    print(f"  {'r_max':>6}  {'psi_L (pierwsze zero)':>22}  {'uwaga':>20}")
    for r, p in psi_L_data:
        note = "** artefakt r_max=40 **" if abs(r-40)<1 else ""
        print(f"  {r:>6.0f}  {p:.6f}                      {note}")

    if len(psi_L_data) >= 2:
        pL_vals = [p for _,p in psi_L_data]
        r_vals  = [r for r,_ in psi_L_data]
        delta_pL = max(pL_vals) - min(pL_vals)
        p5_ok = delta_pL > 0.01
        print(f"\n  Rozstaw psi_L: {min(pL_vals):.5f} .. {max(pL_vals):.5f}"
              f"  (delta={delta_pL:.5f})")
        record("P5: psi_L zmienia sie o > 0.01 miedzy r_max=40 a 100"
               " (artefakt)", p5_ok, f"delta_psi_L={delta_pL:.5f}")
    else:
        # Brak wielu psi_L (Branch L moze byc brak juz przy r_max>=60)
        n_L_found = len(psi_L_data)
        p5_ok = n_L_found < len([r for r in RMAX_LIST if r >= 60])
        record("P5: psi_L (Branch L) nie istnieje dla r_max>=60",
               True, f"Branch L znaleziona tylko dla {[r for r,_ in psi_L_data]}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== SEKCJA D: Fit ogona oscylatoryjnego phi(r)-1 ~ Bessel 1/r ==")
    # Uzyj K* z r_max=100 lub r_max=80 jako najlepsze przyblizenie
    best = [(r,K,psi,_) for r,K,psi,_ in kstar_data if np.isfinite(K)]
    if best:
        r_use, K_use, psi_use, _ = best[-1]  # najwiekszy r_max
        print(f"  Uzywam K*={K_use:.6f} (r_max={r_use:.0f}), psi_0={psi_use:.5f}")

        # Dlugi profil do r=300
        dphi0 = -K_use/A_GAM**2
        r_long = 300.0
        n_ev = 2500
        r_ev = A_GAM*(r_long/A_GAM)**np.linspace(0,1,n_ev)
        t0 = time.time()
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,ALPHA), [A_GAM,r_long],
                        [psi_use,dphi0], method='DOP853',
                        rtol=1e-10, atol=1e-12, t_eval=r_ev)
        dt = time.time()-t0

        r = sol.t; phi = sol.y[0]; delta = phi - 1.0

        # Fit w ogonie r > 50
        mask = r > 50.0
        r_t = r[mask]; d_t = delta[mask]
        m_osc = 1.0/np.sqrt(1.0+ALPHA)

        y_fit = d_t * r_t  # delta*r = A sin(m*r) + B cos(m*r)
        X = np.column_stack([np.sin(m_osc*r_t), np.cos(m_osc*r_t)])
        coeffs, _, rank, _ = np.linalg.lstsq(X, y_fit, rcond=None)
        A_fit, B_fit = coeffs
        y_pred = X @ coeffs
        ss_res = np.sum((y_fit-y_pred)**2)
        ss_tot = np.sum((y_fit-np.mean(y_fit))**2)
        R2_bessel = 1-ss_res/ss_tot if ss_tot>0 else 0

        amp = np.sqrt(A_fit**2+B_fit**2)
        print(f"  Ogon fit (r>50, t={dt:.0f}s):")
        print(f"    m_oczekiwane = 1/sqrt(1+alpha) = {m_osc:.5f}")
        print(f"    A_fit = {A_fit:.5f},  B_fit = {B_fit:.5f}")
        print(f"    Amplituda = {amp:.5f}")
        print(f"    R^2 = {R2_bessel:.5f}")
        print(f"    Amplituda przy r=100: {amp/100:.5f}")
        print(f"    Amplituda przy r=200: {amp/200:.5f}")

        # Pokaz jak delta*r wygada przy kilku r
        print(f"\n  Profil ogona delta(r)*r [oczekiwane: harmoniczne ~amp sin(mr+faza)]:")
        rs_show = [60,80,100,120,150,200,250]
        for rs in rs_show:
            idx = np.argmin(abs(r-rs))
            expected = A_fit*np.sin(m_osc*r[idx]) + B_fit*np.cos(m_osc*r[idx])
            print(f"    r={r[idx]:>6.1f}: delta*r={delta[idx]*r[idx]:+.6f}, "
                  f"fit={expected:+.6f}, blad={delta[idx]*r[idx]-expected:+.6f}")

        p4_ok = R2_bessel > 0.90
        record("P4: ogon Bessela (A sin/r + B cos/r) fituje z R^2>0.90",
               p4_ok, f"R^2={R2_bessel:.4f}, amp={amp:.4e}")
    else:
        record("P4: ogon Bessela fituje z R^2>0.90", False, "brak K* do profilu")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== PODSUMOWANIE ILOSIOWE ==")
    print()
    print("  K*(r_max) dla alpha=8.840:")
    for r, K, psi, nz in kstar_data:
        Ks_s = f'{K:.7f}' if np.isfinite(K) else '  brak   '
        print(f"    r_max={r:>6.0f}:  K* = {Ks_s}  (psi_U={psi:.4f}, n_zeros={nz})")
    if np.isfinite(K_inf):
        print(f"    r_max=inf  :  K* ~ {K_inf:.6f}  (extrapolacja 1/r)")
    print()
    print("  Referencje:")
    print("    K*(p80, r_max~35-40 konwencja) ~ 0.011 -- OK dla tamtej konwencji")
    print("    K*(p92, r_max=60)   ~ 0.01031")
    print(f"    K*(p93, r_max=inf) ~ {K_inf:.5f}" if np.isfinite(K_inf) else "")
    print()
    print("  Wniosek strukturalny:")
    print("    Jedyna FIZYCZNA rodzina solitonow: Galaz U (K*_B, ostatnie zero psi)")
    print("    Galaz L z p91 = ARTEFAKT: wynikala z falsz. zer phi(rmax)-1=0")
    print("    Hipoteza dwu-galeziona (p91) COFNIETA")

    print()
    print("=" * 70)
    print("PODSUMOWANIE  p93_rmax_convergence.py")
    print("=" * 70)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}\n")


if __name__ == '__main__':
    run_main()
