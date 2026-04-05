#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p88_cascade_bifurcation.py  --  TGP v1  *  Kaskada bifurkacji solitonow
========================================================================
Motywacja (wyniki p86-p87):
  g(K,alpha)=0 ma 4 zera dla alpha<8.77, potem kaskada saddle-node:
    - alpha_SN34 ~ 8.77:  K*_3 i K*_4 lacza sie
    - alpha_SN23 ~ 8.86:  K*_2 i K*_3 lacza sie
    - alpha_max = 8.8734: K*_1 znika (bifurkacja graniczna)

  Pytanie: czy wykladnik gamma jest UNIVERSALNY dla calej kaskady?
    gamma_34  (fit delta_K34 ~ (alpha_SN34-alpha)^gamma)
    gamma_23  (fit delta_K23 ~ (alpha_SN23-alpha)^gamma)
    gamma_12  (z p87: gamma ~ 0.55)
  Jesli gamma_34 ~ gamma_23 ~ gamma_12 ~ 0.5 -> jeden wykladnik uniwersalny.

PLAN:
  A. Skan alpha in [8.72, 8.80] (10 pkt) + K grid [0.010, 0.045]
     -> znajdz K*_3 i K*_4 (tam gdzie istnieja)
     -> fit delta_K34(alpha)
  B. Skan alpha in [8.84, 8.866] (12 pkt) + K grid [0.010, 0.020]
     -> znajdz K*_2 i K*_3 (tam gdzie istnieja)
     -> fit delta_K23(alpha)
  C. Porowniez gamma_34, gamma_23, gamma_12=0.55

Testy:
  P1: K*_3,K*_4 znalezione dla >= 5 alpha w [8.72,8.80]
  P2: K*_2,K*_3 znalezione dla >= 5 alpha w [8.84,8.866]
  P3: Oba fity zbiezne (RMSE < 5%)
  P4: |gamma_34 - gamma_23| < 0.15  (ta sama klasa universalnosci?)
  P5: Wszystkie gamma in [0.35, 0.75]  (blisko 0.5)

Workers: 8, maxtasksperchild=3.
Data: 2026-03-25
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA_V = 1.0
V1     = GAMMA_V/3 - GAMMA_V/4

# Wyniki z poprzednich skryptow
ALPHA_MAX_REF = 8.8734
GAMMA_12      = 0.551   # z p87, staly alpha_max

def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=1400):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0]-1e-5
    ev.terminal=True; ev.direction=-1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0,1,n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi,dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0,-1]) if sol.t[-1]>=r_max*0.99 else np.nan
    except: return np.nan

def energy_at_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0,1,1400)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi_z,dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek+Ep
    except: return np.nan

def g_at_K(K, a_gam, alpha, n_psi=80):
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [phi_at_rmax(p,K,a_gam,alpha)-1 for p in psis]
    best = None
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,a_gam,alpha)-1,
                            psis[i],psis[i+1], xtol=1e-7, maxiter=50)
                E = energy_at_psi(pz,K,a_gam,alpha)
                if np.isfinite(E):
                    gv = E/(4*np.pi*K)-1
                    if best is None or abs(gv)<abs(best[0]):
                        best = (gv, pz, E)
            except: pass
    return best[0] if best else np.nan

def scan_zeros(alpha, a_gam, K_arr, n_psi=80):
    """Znajdz wszystkie zera g(K) w K_arr. Zwraca posortowana liste."""
    g_vals = [g_at_K(K, a_gam, alpha, n_psi) for K in K_arr]
    zeros = []
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                Kz = brentq(lambda K: g_at_K(K, a_gam, alpha, n_psi),
                            K_arr[i], K_arr[i+1],
                            xtol=K_arr[i]*3e-4, rtol=3e-4, maxiter=25)
                zeros.append(Kz)
            except: pass
    return sorted(zeros)

def worker_A(args):
    """Sekcja A: alpha in [8.72,8.80], K grid dla K*_3, K*_4."""
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()
    zeros = scan_zeros(alpha, a_gam, K_arr)
    dt = time.time()-t0
    z_str = ', '.join(f'{z:.5f}' for z in zeros)
    print(f'  A[{idx+1:2d}/{total}] a={alpha:.4f}  zeros=[{z_str}]  t={dt:.1f}s',
          flush=True)
    return (alpha, zeros)

def worker_B(args):
    """Sekcja B: alpha in [8.84,8.866], K grid dla K*_2, K*_3 (bliskie)."""
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()
    zeros = scan_zeros(alpha, a_gam, K_arr, n_psi=100)
    dt = time.time()-t0
    z_str = ', '.join(f'{z:.6f}' for z in zeros)
    print(f'  B[{idx+1:2d}/{total}] a={alpha:.4f}  zeros=[{z_str}]  t={dt:.1f}s',
          flush=True)
    return (alpha, zeros)


def fit_saddle_node(alphas, deltas, label="", amax_guess=None):
    """
    Fit delta_K = C*(alpha_SN - alpha)^gamma.
    Zwraca (gamma, C, alpha_SN, rmse) lub (nan,nan,nan,inf).
    """
    av = np.array(alphas); dv = np.array(deltas)
    if len(av) < 4: return np.nan, np.nan, np.nan, np.inf

    best = None
    if amax_guess is None:
        amax_guess = av[-1] + 0.01

    def model(a, C, g, amax):
        return C * np.maximum(amax-a, 1e-12)**g

    for g0 in [0.3, 0.5, 0.7, 1.0]:
        for da in [-0.01, 0, 0.01, 0.02]:
            try:
                popt, pcov = curve_fit(model, av, dv,
                                       p0=[0.01, g0, amax_guess+da],
                                       bounds=([0, 0.1, av[-1]], [1, 2, av[-1]+0.05]),
                                       maxfev=3000)
                pred = model(av, *popt)
                rmse = np.sqrt(np.mean((pred-dv)**2))/dv.mean()*100
                if best is None or rmse < best[3]:
                    best = (popt[1], popt[0], popt[2], rmse,
                            np.sqrt(np.diag(pcov)))
            except: pass

    if best:
        g, C, amax, rmse, perr = best
        print(f"  Fit {label}: gamma={g:.4f}±{perr[1]:.4f}  "
              f"alpha_SN={amax:.5f}±{perr[2]:.5f}  RMSE={rmse:.3f}%")
        return g, C, amax, rmse
    return np.nan, np.nan, np.nan, np.inf


def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')
    PC=0; FC=0
    N_WORKERS = min(8, cpu_count())

    def record(label, ok, info=""):
        nonlocal PC, FC
        if ok: PC+=1
        else:  FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    print("=" * 72)
    print("TGP v1  *  p88_cascade_bifurcation.py  --  Kaskada saddle-node")
    print("  gamma_34, gamma_23, gamma_12=0.551 -- czy universalne?")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF}")
    print(f"  gamma_12 (z p87) = {GAMMA_12}")
    print(f"  Workers = {N_WORKERS}\n")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA A: K*_3 i K*_4 (alpha in [8.72, 8.80])
    # ─────────────────────────────────────────────────────────────────────────
    print("=" * 72)
    print("SEKCJA A: K*_3 i K*_4  (alpha in [8.72, 8.80])")
    print("  K grid: [0.010, 0.060] zageszczony tam gdzie K*_3,K*_4 oczekiwane")
    print("=" * 72)

    # K*_3 ~ 0.016-0.022, K*_4 ~ 0.030-0.045 (z p86)
    K_A_low  = np.linspace(0.010, 0.015, 6)   # K*_1 i K*_2 region (tlo)
    K_A_mid  = np.linspace(0.015, 0.050, 28)   # K*_3 i K*_4 region
    K_A_high = np.linspace(0.050, 0.065, 5)
    K_arr_A  = np.unique(np.concatenate([K_A_low, K_A_mid, K_A_high]))
    K_arr_A.sort()

    alpha_A = np.linspace(8.72, 8.80, 10)
    args_A  = [(a, A_GAM, K_arr_A, i, len(alpha_A))
               for i, a in enumerate(alpha_A)]

    t0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        res_A = pool.map(worker_A, args_A, chunksize=1)
    res_A.sort(key=lambda x: x[0])
    print(f"\n  Sekcja A: {time.time()-t0:.1f}s")

    # Wyodrebnij K*_3, K*_4
    print(f"\n  {'alpha':>7}  {'n_zer':>5}  {'K*_3':>10}  {'K*_4':>10}  {'dK_34':>10}")
    print("  " + "-" * 50)
    pts_34 = []
    for alpha, zeros in res_A:
        K3 = zeros[2] if len(zeros)>=3 else np.nan
        K4 = zeros[3] if len(zeros)>=4 else np.nan
        dK = K4-K3 if np.isfinite(K3) and np.isfinite(K4) else np.nan
        K3s = f'{K3:.6f}' if np.isfinite(K3) else '  brak  '
        K4s = f'{K4:.6f}' if np.isfinite(K4) else '  brak  '
        dKs = f'{dK:.7f}' if np.isfinite(dK) else '  brak  '
        print(f"  {alpha:>7.4f}  {len(zeros):>5}  {K3s}  {K4s}  {dKs}")
        if np.isfinite(dK): pts_34.append((alpha, dK))

    n_34 = len(pts_34)
    print(f"\n  Punkty do fitu: {n_34}")

    # Fit gamma_34
    print("\n  --- Fit gamma_34 ---")
    if n_34 >= 4:
        gamma_34, C_34, aSN_34, rmse_34 = fit_saddle_node(
            [x[0] for x in pts_34], [x[1] for x in pts_34],
            label="gamma_34", amax_guess=8.79)
    else:
        gamma_34, C_34, aSN_34, rmse_34 = np.nan, np.nan, np.nan, np.inf
        print("  Za malo punktow.")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA B: K*_2 i K*_3 (alpha in [8.84, 8.866])
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA B: K*_2 i K*_3  (alpha in [8.840, 8.866])")
    print("  K grid: [0.010, 0.018] zageszczony (K*_2~0.011, K*_3~0.014-0.013)")
    print("=" * 72)

    # K*_2 ~ 0.012-0.014, K*_3 ~ 0.013-0.016 (z p86/p87)
    K_arr_B = np.linspace(0.0100, 0.0185, 32)

    alpha_B_main = np.linspace(8.840, 8.860, 9)
    alpha_B_fine = np.array([8.861, 8.862, 8.863, 8.864, 8.865, 8.866])
    alpha_B = np.unique(np.concatenate([alpha_B_main, alpha_B_fine]))
    alpha_B.sort()

    args_B = [(a, A_GAM, K_arr_B, i, len(alpha_B))
              for i, a in enumerate(alpha_B)]

    t0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        res_B = pool.map(worker_B, args_B, chunksize=1)
    res_B.sort(key=lambda x: x[0])
    print(f"\n  Sekcja B: {time.time()-t0:.1f}s")

    # Wyodrebnij K*_2, K*_3 z drugiego i trzeciego zera
    print(f"\n  {'alpha':>8}  {'n_zer':>5}  {'K*_2':>10}  {'K*_3':>10}  "
          f"{'dK_23':>10}  {'dK_12':>10}")
    print("  " + "-" * 62)
    pts_23 = []
    for alpha, zeros in res_B:
        K1 = zeros[0] if len(zeros)>=1 else np.nan
        K2 = zeros[1] if len(zeros)>=2 else np.nan
        K3 = zeros[2] if len(zeros)>=3 else np.nan
        dK23 = K3-K2 if np.isfinite(K2) and np.isfinite(K3) else np.nan
        dK12 = K2-K1 if np.isfinite(K1) and np.isfinite(K2) else np.nan
        K2s = f'{K2:.6f}' if np.isfinite(K2) else '  brak  '
        K3s = f'{K3:.6f}' if np.isfinite(K3) else '  brak  '
        dK23s = f'{dK23:.7f}' if np.isfinite(dK23) else '  brak  '
        dK12s = f'{dK12:.7f}' if np.isfinite(dK12) else '  brak  '
        print(f"  {alpha:>8.4f}  {len(zeros):>5}  {K2s}  {K3s}  {dK23s}  {dK12s}")
        if np.isfinite(dK23): pts_23.append((alpha, dK23))

    n_23 = len(pts_23)
    print(f"\n  Punkty do fitu: {n_23}")

    print("\n  --- Fit gamma_23 ---")
    if n_23 >= 4:
        gamma_23, C_23, aSN_23, rmse_23 = fit_saddle_node(
            [x[0] for x in pts_23], [x[1] for x in pts_23],
            label="gamma_23", amax_guess=8.864)
    else:
        gamma_23, C_23, aSN_23, rmse_23 = np.nan, np.nan, np.nan, np.inf
        print("  Za malo punktow.")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA C: Porowniez gammy
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA C: Porowniez wykladnikow bifurkacji")
    print("=" * 72)
    print(f"""
  Bifurkacja      gamma      alpha_SN      RMSE
  ─────────────────────────────────────────────
  K*_3 ↔ K*_4    {gamma_34:.4f}    {aSN_34:.5f}    {rmse_34:.3f}%
  K*_2 ↔ K*_3    {gamma_23:.4f}    {aSN_23:.5f}    {rmse_23:.3f}%
  K*_1 → boundary {GAMMA_12:.4f}    {ALPHA_MAX_REF:.5f}   (ref, z p87)
  ─────────────────────────────────────────────
  Srednia (34+23): {np.nanmean([gamma_34,gamma_23]):.4f}
  Odch. std:       {np.nanstd([gamma_34,gamma_23]):.4f}
  |gamma_34 - gamma_12| = {abs(gamma_34-GAMMA_12):.4f}
  |gamma_23 - gamma_12| = {abs(gamma_23-GAMMA_12):.4f}
""")

    # ─────────────────────────────────────────────────────────────────────────
    # TESTY
    # ─────────────────────────────────────────────────────────────────────────
    print("=" * 72)
    print("TESTY")
    print("=" * 72)

    record("P1: K*_3,K*_4 znalezione dla >= 5 alpha w [8.72,8.80]",
           n_34 >= 5, f"{n_34} pkt")
    record("P2: K*_2,K*_3 znalezione dla >= 5 alpha w [8.84,8.866]",
           n_23 >= 5, f"{n_23} pkt")
    record("P3: Oba fity zbiezne (RMSE < 5%)",
           rmse_34 < 5.0 and rmse_23 < 5.0,
           f"RMSE_34={rmse_34:.3f}%, RMSE_23={rmse_23:.3f}%")
    record("P4: |gamma_34 - gamma_23| < 0.15  (ta sama klasa)",
           np.isfinite(gamma_34) and np.isfinite(gamma_23) and
           abs(gamma_34-gamma_23) < 0.15,
           f"gamma_34={gamma_34:.4f}, gamma_23={gamma_23:.4f}")
    record("P5: Wszystkie gamma in [0.35, 0.75]",
           all(np.isfinite(g) and 0.35<=g<=0.75
               for g in [gamma_34, gamma_23, GAMMA_12]),
           f"[{gamma_34:.3f}, {gamma_23:.3f}, {GAMMA_12:.3f}]")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p88_cascade_bifurcation.py")
    print("=" * 72)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    all_finite = all(np.isfinite(g) for g in [gamma_34, gamma_23])
    if all_finite:
        gmean = np.nanmean([gamma_34, gamma_23, GAMMA_12])
        gstd  = np.nanstd([gamma_34, gamma_23, GAMMA_12])
        universal = abs(gamma_34-gamma_23)<0.15 and abs(gamma_34-GAMMA_12)<0.2
        print(f"""
  gamma (srednia 3 bifurkacji) = {gmean:.4f}  ± {gstd:.4f}
  Odleglosc od 0.5: {abs(gmean-0.5):.4f}
  {'UNIVERSALNOSC POTWIERDZONA' if universal else 'brak universalnosci w tej tolerancji'}
""")
        if abs(gmean - 0.5) < 0.1:
            print("  WNIOSEK: gamma ~ 0.5 dla calej kaskady.")
            print("  Wszystkie bifurkacje TGP sa quasi-saddle-node.")
            print("  Sugeruje wspolny mechanizm (1D fold) w calej rodzinie solitonow.")
        elif abs(gmean - 0.5) < 0.2:
            print(f"  WNIOSEK: gamma ~ {gmean:.2f}, bliskie 0.5 ale z systematycznym odchyleniem.")
            print("  Mozliwy wplyw nielokalnej geometrii ODE (niestandart. saddle-node).")
        else:
            print(f"  WNIOSEK: gamma = {gmean:.2f} - nie klasyczne saddle-node.")
    else:
        print("  Zbyt malo danych -- potrzeba wiecej punktow lub szerszego K_arr.")


if __name__ == '__main__':
    run_main()
