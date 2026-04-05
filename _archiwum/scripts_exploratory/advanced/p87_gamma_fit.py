#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p87_gamma_fit.py  --  TGP v1  *  Wykladnik bifurkacji gamma
==============================================================
Kontekst (wynik p86):
  g(K,alpha) ma 4 zera dla alpha=8.5-8.7 (K*_1 < K*_2 < K*_3 < K*_4).
  K*_2/K*_1 spada monotoniczne: 2.81 (a=8.5) -> 1.18 (a=8.85).
  K*_2 zanika miedzy alpha=8.85 i alpha=8.87.
  -> Saddle-node K*_1/K*_2 przy alpha_max ~ 8.86-8.87.

Cel: Precyzyjny skan alpha in [8.80, 8.873], step 0.005.
  Dla kazdego alpha: znalezc K*_1 i K*_2.
  Fit: delta_K = K*_2 - K*_1 = C*(alpha_max_SN - alpha)^gamma.
  -> gamma ~ 0.5 (saddle-node)? ~ 1.0 (transcritical)? ~ 0.33 (cusp)?

Testy:
  P1: K*_2 znalezione dla >= 8 wartosci alpha
  P2: delta_K monotonicznie maleje
  P3: Fit zbiezny (RMSE < 5%)
  P4: gamma in [0.3, 0.8]  (gdzies w tej okolicy)
  P5: alpha_max_SN in [8.855, 8.877]  (blisko alpha_max_ref=8.8734)

Optymalizacja: dla kazdego alpha jeden skan g(K) z n_K=45, n_psi=100.
8 workers rownoleglo. Czas: ~60s/alpha -> ~8 min total.
Data: 2026-03-25
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_MAX_REF = 8.8734
A_GAM         = 0.040
LAM           = 5.501357e-06
GAMMA_V       = 1.0
V1            = GAMMA_V/3 - GAMMA_V/4

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

def g_at_K(K, a_gam, alpha, n_psi=100):
    """g(K) = E/(4piK)-1 dla galezie min|g|. Nan jesli brak rozwiazania."""
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
                    g_val = E/(4*np.pi*K)-1
                    if best is None or abs(g_val)<abs(best[0]):
                        best = (g_val, pz, E)
            except: pass
    return best[0] if best else np.nan

def find_K1_K2(alpha, a_gam, K_arr, n_psi=100):
    """
    Skanuj g(K) na K_arr. Znajdz K*_1 i K*_2 (pierwsze dwa zera).
    Zwraca (K1, K2, g_arr).
    """
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
    zeros.sort()
    K1 = zeros[0] if len(zeros)>=1 else np.nan
    K2 = zeros[1] if len(zeros)>=2 else np.nan
    return K1, K2, g_vals

def worker(args):
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()
    K1, K2, g_arr = find_K1_K2(alpha, a_gam, K_arr)
    dt = time.time()-t0
    dK = K2-K1 if np.isfinite(K2) and np.isfinite(K1) else np.nan
    dK_str = f'{dK:.6f}' if np.isfinite(dK) else '  brak'
    K2_str = f'{K2:.6f}' if np.isfinite(K2) else '  brak'
    print(f'  [{idx+1:3d}/{total}] a={alpha:.4f}  K*1={K1:.6f}  '
          f'K*2={K2_str}  dK={dK_str}  t={dt:.1f}s', flush=True)
    return (alpha, K1, K2)


def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')
    PC=0; FC=0
    N_WORKERS = min(8, cpu_count())

    def record(label, ok, info=""):
        nonlocal PC, FC
        if ok: PC+=1;
        else:  FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    print("=" * 72)
    print("TGP v1  *  p87_gamma_fit.py  --  Wykladnik bifurkacji K*_2 - K*_1")
    print("  delta_K = K*_2 - K*_1  ~  C*(alpha_max_SN - alpha)^gamma")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF}")
    print(f"  Workers = {N_WORKERS}\n")

    # K grid dla kazdego alpha:
    # Gestszy w [0.009, 0.025] gdzie K*_1 i K*_2 sa spodziewane
    K_low  = np.logspace(np.log10(0.005), np.log10(0.009), 8)[:-1]
    K_mid  = np.linspace(0.009, 0.025, 25)   # obejmuje K*_1 i K*_2
    K_high = np.logspace(np.log10(0.025), np.log10(0.10), 8)[1:]
    K_arr  = np.unique(np.concatenate([K_low, K_mid, K_high]))
    K_arr.sort()
    print(f"  K grid: {len(K_arr)} pkt w [{K_arr[0]:.4f}, {K_arr[-1]:.4f}]")
    print(f"  (zageszczenie w [0.009, 0.025] gdzie K*_1 i K*_2)\n")

    # Alpha: rowno od 8.80 do 8.873, step 0.005 + kilka punktow gesto
    alpha_coarse = np.arange(8.80, 8.87, 0.005)
    alpha_fine   = np.array([8.852, 8.856, 8.860, 8.864, 8.868, 8.871])
    alpha_scan   = np.unique(np.concatenate([alpha_coarse, alpha_fine]))
    alpha_scan.sort()
    print(f"  alpha scan: {len(alpha_scan)} wartosci od {alpha_scan[0]:.4f} do {alpha_scan[-1]:.4f}")

    args_list = [(a, A_GAM, K_arr, i, len(alpha_scan))
                 for i, a in enumerate(alpha_scan)]

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA A: K*_1 i K*_2 (rownolegly)")
    print("=" * 72)

    t0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        res_raw = pool.map(worker, args_list, chunksize=1)
    t_tot = time.time()-t0
    res_raw.sort(key=lambda x: x[0])
    print(f"\n  Lacznie: {t_tot:.1f}s")

    valid_both = [(a,K1,K2) for a,K1,K2 in res_raw
                  if np.isfinite(K1) and np.isfinite(K2)]
    n_K2 = len(valid_both)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA B: Tabela delta_K(alpha)")
    print("=" * 72)
    print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'K*_2':>10}  "
          f"{'delta_K':>10}  {'K*2/K*1':>8}")
    print("  " + "-" * 54)
    for a,K1,K2 in res_raw:
        if np.isfinite(K1):
            K2s = f'{K2:.6f}' if np.isfinite(K2) else '  brak  '
            dK = f'{K2-K1:.7f}' if np.isfinite(K2) else '  brak  '
            rat = f'{K2/K1:.4f}' if np.isfinite(K2) else '  n/a '
            print(f"  {a:>8.4f}  {K1:>10.6f}  {K2s}  {dK}  {rat}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA C: Fit potegowy delta_K = C*(alpha_max_SN - alpha)^gamma")
    print("=" * 72)

    gamma_fit = np.nan; C_fit = np.nan; rmse_fit = np.inf
    amax_SN   = np.nan

    if n_K2 >= 4:
        av  = np.array([x[0] for x in valid_both])
        dKv = np.array([x[2]-x[1] for x in valid_both])

        # Sprawdz monotonocznosc
        is_mon = bool(np.all(np.diff(dKv) < 0))
        print(f"\n  n_K2 = {n_K2},  delta_K monotonicznie maleje: {'TAK' if is_mon else 'NIE'}")

        # Fit ze swobodnym alpha_max_SN
        def model(a, C, g, amax):
            return C * np.maximum(amax - a, 1e-12)**g

        # Probuj rozne punkty startowe
        best_res = None
        for g0 in [0.3, 0.5, 0.7, 1.0]:
            for amax0 in [8.860, 8.865, 8.870, 8.875]:
                try:
                    popt, pcov = curve_fit(model, av, dKv,
                                           p0=[0.01, g0, amax0],
                                           bounds=([0, 0.1, 8.850], [1, 2, 8.880]),
                                           maxfev=3000)
                    pred = model(av, *popt)
                    rmse = np.sqrt(np.mean((pred-dKv)**2))/dKv.mean()*100
                    if best_res is None or rmse < best_res[3]:
                        best_res = (popt[0], popt[1], popt[2], rmse,
                                    np.sqrt(np.diag(pcov)))
                except: pass

        if best_res:
            C_fit, gamma_fit, amax_SN, rmse_fit, perr = best_res
            print(f"\n  Najlepszy fit (swobodny alpha_max_SN):")
            print(f"    C          = {C_fit:.6f}  ± {perr[0]:.6f}")
            print(f"    gamma      = {gamma_fit:.5f}  ± {perr[1]:.5f}")
            print(f"    alpha_max_SN = {amax_SN:.5f}  ± {perr[2]:.5f}")
            print(f"    RMSE       = {rmse_fit:.4f}%")
            print(f"    |amax_SN - alpha_max_ref| = {abs(amax_SN-ALPHA_MAX_REF):.5f}")

        # Tez fit z ustalonym alpha_max
        def model2(a, C, g): return C*(ALPHA_MAX_REF-a)**g
        x_fix = ALPHA_MAX_REF - av
        if np.all(x_fix > 0):
            try:
                popt2, _ = curve_fit(model2, av, dKv, p0=[0.01, 0.5],
                                     bounds=([0,0.1],[1,2]), maxfev=2000)
                pred2 = model2(av, *popt2)
                rmse2 = np.sqrt(np.mean((pred2-dKv)**2))/dKv.mean()*100
                print(f"\n  Fit (staly alpha_max={ALPHA_MAX_REF}):")
                print(f"    C     = {popt2[0]:.6f}")
                print(f"    gamma = {popt2[1]:.5f}")
                print(f"    RMSE  = {rmse2:.4f}%")
                if rmse2 < rmse_fit:
                    gamma_fit = popt2[1]; C_fit = popt2[0]; rmse_fit = rmse2
            except: pass
    else:
        is_mon = False
        print(f"\n  Za malo punktow z K*_2 ({n_K2}) do fitu.")

    # ─────────────────────────────────────────────────────────────────────────
    # TESTY
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TESTY")
    print("=" * 72)

    record("P1: K*_2 znalezione dla >= 8 wartosci alpha", n_K2 >= 8,
           f"{n_K2} pkt")
    record("P2: delta_K monotonicznie maleje",
           n_K2>=3 and is_mon, "TAK" if (n_K2>=3 and is_mon) else "NIE")
    record("P3: Fit zbiezny (RMSE < 5%)",
           rmse_fit < 5.0, f"RMSE={rmse_fit:.4f}%" if np.isfinite(rmse_fit) else "brak")
    record("P4: gamma in [0.3, 0.8]",
           np.isfinite(gamma_fit) and 0.3 <= gamma_fit <= 0.8,
           f"gamma={gamma_fit:.5f}" if np.isfinite(gamma_fit) else "brak")
    record("P5: alpha_max_SN in [8.855, 8.877]",
           np.isfinite(amax_SN) and 8.855 <= amax_SN <= 8.877,
           f"amax_SN={amax_SN:.5f}" if np.isfinite(amax_SN) else "brak")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p87_gamma_fit.py")
    print("=" * 72)
    print(f"\n  PASS: {PC} / {PC+FC}")
    print(f"  FAIL: {FC} / {PC+FC}")
    print(f"""
  WYNIKI:
    n_K2         = {n_K2}
    gamma        = {gamma_fit:.5f}  (interpretacja: 0.5=saddle-node, 1.0=linear)
    C            = {C_fit:.6f}
    RMSE         = {rmse_fit:.4f}%
    alpha_max_SN = {amax_SN:.5f}
""")

    if np.isfinite(gamma_fit):
        if abs(gamma_fit - 0.5) < 0.08:
            print("  WNIOSEK: gamma ~ 0.5  -->  KLASYCZNE SADDLE-NODE")
            print("  delta_K = K*_2 - K*_1 ~ sqrt(alpha_max - alpha)")
        elif abs(gamma_fit - 1.0) < 0.15:
            print("  WNIOSEK: gamma ~ 1.0  -->  LINIOWE ZANIKANIE")
            print("  delta_K ~ (alpha_max - alpha) -- nie klasyczne saddle-node")
        elif gamma_fit < 0.5:
            print(f"  WNIOSEK: gamma={gamma_fit:.4f} < 0.5 -- ostrzejsze niz saddle-node")
            print("  Mozliwy cusp lub wyzszy rzad.")
        else:
            print(f"  WNIOSEK: gamma={gamma_fit:.4f} -- miedzy saddle-node a liniowe")
    else:
        print("  WNIOSEK: brak gamma -- potrzeba wiecej punktow z K*_2.")


if __name__ == '__main__':
    run_main()
