#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p89_codim2_verify.py  --  TGP v1  *  Weryfikacja bifurkacji kodimension-2
=========================================================================
Hipoteza (z p88):
  Przy alpha_max = 8.8734 zachodza JEDNOCZESNIE:
    (C1)  K*_2(alpha) = K*_3(alpha)    [saddle-node K*_2/K*_3]
    (C2)  K*_1(alpha) = K_nan_left(alpha)  [K*_1 dotyka granicy zapad.]

  Jesli alpha_SN23 = alpha_max -> bifurkacja kodimension-2.

Plan:
  Gesty skan g(K) dla alpha in [8.855, 8.873], krok 0.002.
  K grid: [0.0108, 0.0185], 65 pkt.
  -> znajdz K*_1, K*_2, K*_3 (gdy istnieja)
  -> znajdz K_nan_left = lewy brzeg nan-region w g(K)
  -> delta23 = K*_3 - K*_2,  delta1c = K_nan_left - K*_1
  -> ekstrapolacja: gdzie oba zanikaja?

Testy:
  P1: K*_2 i K*_3 znalezione dla >= 5 alpha
  P2: delta23 < 0.0003 przy alpha bliskim alpha_SN23 (K*_2 i K*_3 bliskie)
  P3: delta1c < 0.0003 blisko alpha_max (K*_1 i K_nan_left bliskie)
  P4: |alpha_SN23 - alpha_max_ref| < 0.003  (oba zbiegaja sie blisko?)
  P5: |alpha_SN23 - alpha_zero_1c| < 0.003  (oba zbiegaja sie ze soba?)

Workers: 8.
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
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA_V = 1.0
V1     = GAMMA_V/3 - GAMMA_V/4

def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=1500):
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
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0,1,1500)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [a_gam,r_max],
                        [psi_z,dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek+Ep
    except: return np.nan

def g_at_K_all_branches(K, a_gam, alpha, n_psi=100):
    """
    Zwraca (g_min, is_nan) gdzie:
      g_min = wartosc g dla galezie min|g| (lub nan jesli brak rozwiazania)
      is_nan = True jesli ODE zapada sie dla wszystkich psi (nan-region)
    """
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [phi_at_rmax(p,K,a_gam,alpha)-1 for p in psis]
    best = None
    any_finite = any(np.isfinite(f) for f in Fv)
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,a_gam,alpha)-1,
                            psis[i],psis[i+1], xtol=1e-7, maxiter=60)
                E = energy_at_psi(pz,K,a_gam,alpha)
                if np.isfinite(E):
                    gv = E/(4*np.pi*K)-1
                    if best is None or abs(gv)<abs(best):
                        best = gv
            except: pass
    is_nan = not any_finite
    return (best if best is not None else np.nan), is_nan


def worker_scan(args):
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()

    g_vals = []; is_nan_flags = []
    for K in K_arr:
        gv, inan = g_at_K_all_branches(K, a_gam, alpha)
        g_vals.append(gv)
        is_nan_flags.append(inan)

    # Znajdz zera g(K)
    zeros = []
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                Kz = brentq(lambda K: g_at_K_all_branches(K, a_gam, alpha)[0],
                            K_arr[i], K_arr[i+1],
                            xtol=K_arr[i]*2e-4, rtol=2e-4, maxiter=25)
                zeros.append(Kz)
            except: pass
    zeros.sort()

    # Znajdz lewy brzeg nan-region (K_nan_left)
    K_nan_left = np.nan
    for i, inan in enumerate(is_nan_flags):
        if inan:
            K_nan_left = K_arr[i]
            break
    # Doprecyzuj bisekcja
    if np.isfinite(K_nan_left) and K_nan_left > K_arr[0]:
        K_lo = K_arr[max(0, is_nan_flags.index(True)-1)]
        K_hi = K_nan_left
        for _ in range(15):
            K_mid = 0.5*(K_lo+K_hi)
            _, inan_mid = g_at_K_all_branches(K_mid, a_gam, alpha)
            if inan_mid: K_hi = K_mid
            else:        K_lo = K_mid
            if K_hi-K_lo < K_lo*3e-4: break
        K_nan_left = 0.5*(K_lo+K_hi)

    K1 = zeros[0] if len(zeros)>=1 else np.nan
    K2 = zeros[1] if len(zeros)>=2 else np.nan
    K3 = zeros[2] if len(zeros)>=3 else np.nan

    dK23 = K3-K2 if np.isfinite(K2) and np.isfinite(K3) else np.nan
    dK1c = K_nan_left-K1 if np.isfinite(K1) and np.isfinite(K_nan_left) else np.nan

    dt = time.time()-t0
    K2s   = f'{K2:.6f}' if np.isfinite(K2)       else '  brak  '
    K3s   = f'{K3:.6f}' if np.isfinite(K3)       else '  brak  '
    Knls  = f'{K_nan_left:.6f}' if np.isfinite(K_nan_left) else '  brak  '
    dK23s = f'{dK23:.6f}' if np.isfinite(dK23) else '  brak '
    dK1cs = f'{dK1c:.6f}' if np.isfinite(dK1c) else '  brak '
    print(f'  [{idx+1:2d}/{total}] a={alpha:.4f}  '
          f'K1={K1:.6f}  K2={K2s}  K3={K3s}  '
          f'Knan={Knls}  d23={dK23s}  d1c={dK1cs}  t={dt:.1f}s', flush=True)
    return (alpha, K1, K2, K3, K_nan_left)


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
    print("TGP v1  *  p89_codim2_verify.py  --  Weryfikacja kodimension-2")
    print("  Warunek (C1): K*_2=K*_3 przy alpha_SN23")
    print("  Warunek (C2): K*_1=K_nan_left przy alpha_max")
    print("  Pytanie: czy alpha_SN23 = alpha_max = 8.8734?")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF},  Workers = {N_WORKERS}\n")

    # K grid: [0.0108, 0.0185] gestszy (K*_1~0.011, K*_2~0.012, K*_3~0.013)
    K_arr = np.linspace(0.0106, 0.0192, 65)

    # Alpha: od 8.855 (gdzie K*_3 jeszcze istnieje) do 8.873 (blisko alpha_max)
    alpha_scan = np.arange(8.855, 8.8745, 0.002)
    alpha_scan = np.unique(np.append(alpha_scan, [ALPHA_MAX_REF]))
    alpha_scan.sort()

    args_list = [(a, A_GAM, K_arr, i, len(alpha_scan))
                 for i, a in enumerate(alpha_scan)]

    print(f"  K grid: {len(K_arr)} pkt w [{K_arr[0]:.4f}, {K_arr[-1]:.4f}]")
    print(f"  alpha scan: {len(alpha_scan)} wartosci [{alpha_scan[0]:.4f}, {alpha_scan[-1]:.4f}]\n")

    print("=" * 72)
    print("SKAN GLOWNY (rownolegly)")
    print("=" * 72)

    t0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        results = pool.map(worker_scan, args_list, chunksize=1)
    results.sort(key=lambda x: x[0])
    print(f"\n  Czas: {time.time()-t0:.1f}s")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TABELA WYNIKOW")
    print("=" * 72)
    print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'K*_2':>10}  {'K*_3':>10}  "
          f"{'K_nan':>10}  {'d23':>9}  {'d1c':>9}")
    print("  " + "-" * 76)

    pts_23 = []; pts_1c = []
    for alpha, K1, K2, K3, Knan in results:
        dK23 = K3-K2   if np.isfinite(K2) and np.isfinite(K3)   else np.nan
        dK1c = Knan-K1 if np.isfinite(K1) and np.isfinite(Knan) else np.nan
        K2s   = f'{K2:.6f}'   if np.isfinite(K2)   else '  brak  '
        K3s   = f'{K3:.6f}'   if np.isfinite(K3)   else '  brak  '
        Knans = f'{Knan:.6f}' if np.isfinite(Knan) else '  brak  '
        d23s  = f'{dK23:.6f}' if np.isfinite(dK23) else '  brak '
        d1cs  = f'{dK1c:.6f}' if np.isfinite(dK1c) else '  brak '
        mrk = ' <--' if abs(alpha-ALPHA_MAX_REF) < 0.0001 else ''
        print(f"  {alpha:>8.4f}  {K1:.6f}  {K2s}  {K3s}  "
              f"{Knans}  {d23s}  {d1cs}{mrk}")
        if np.isfinite(dK23): pts_23.append((alpha, dK23))
        if np.isfinite(dK1c): pts_1c.append((alpha, dK1c))

    n_23 = len(pts_23); n_1c = len(pts_1c)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("ANALIZA ZBIEZNOSCI")
    print("=" * 72)

    # Ekstrapolacja liniowa dla obu warunkow
    alpha_SN23 = np.nan; alpha_zero_1c = np.nan

    if n_23 >= 3:
        av = np.array([x[0] for x in pts_23])
        dv = np.array([x[1] for x in pts_23])
        print(f"\n  Trend delta_23 (K*_3 - K*_2):")
        for a,d in pts_23: print(f"    alpha={a:.4f}  d23={d:.7f}")

        # Fit liniowy -> zero
        try:
            p = np.polyfit(av, dv, 1)
            alpha_SN23 = -p[1]/p[0]
            print(f"  Ekstrapolacja liniowa d23=0 przy alpha_SN23 = {alpha_SN23:.5f}")
        except: pass

        # Fit potegowy
        if n_23 >= 4:
            def model(a, C, g, amax):
                return C * np.maximum(amax-a, 1e-12)**g
            try:
                popt, pcov = curve_fit(model, av, dv,
                                       p0=[0.001, 0.7, alpha_SN23 if np.isfinite(alpha_SN23) else 8.870],
                                       bounds=([0,0.1,av[-1]], [1,2,av[-1]+0.02]),
                                       maxfev=3000)
                pred = model(av, *popt)
                rmse = np.sqrt(np.mean((pred-dv)**2))/dv.mean()*100
                perr = np.sqrt(np.diag(pcov))
                alpha_SN23 = popt[2]
                print(f"  Fit potegowy: gamma={popt[1]:.4f}±{perr[1]:.4f}  "
                      f"alpha_SN23={popt[2]:.5f}±{perr[2]:.5f}  RMSE={rmse:.3f}%")
            except Exception as e:
                print(f"  Fit potegowy fail: {e}")

    if n_1c >= 3:
        av = np.array([x[0] for x in pts_1c])
        dv = np.array([x[1] for x in pts_1c])
        print(f"\n  Trend delta_1c (K_nan_left - K*_1):")
        for a,d in pts_1c: print(f"    alpha={a:.4f}  d1c={d:.7f}")

        try:
            p = np.polyfit(av, dv, 1)
            alpha_zero_1c = -p[1]/p[0]
            print(f"  Ekstrapolacja liniowa d1c=0 przy alpha = {alpha_zero_1c:.5f}")
        except: pass

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TESTY")
    print("=" * 72)

    record("P1: K*_2 i K*_3 znalezione dla >= 5 alpha",
           n_23 >= 5, f"{n_23} pkt")

    dK23_near = min((d for a,d in pts_23 if a > ALPHA_MAX_REF-0.01), default=np.nan)
    record("P2: delta_23 < 0.0003 blisko alpha_max",
           np.isfinite(dK23_near) and dK23_near < 0.0003,
           f"min_d23={dK23_near:.7f}" if np.isfinite(dK23_near) else "brak")

    dK1c_near = min((d for a,d in pts_1c if a > ALPHA_MAX_REF-0.003), default=np.nan)
    record("P3: delta_1c < 0.0003 blisko alpha_max",
           np.isfinite(dK1c_near) and dK1c_near < 0.0003,
           f"min_d1c={dK1c_near:.7f}" if np.isfinite(dK1c_near) else "brak")

    record("P4: |alpha_SN23 - alpha_max_ref| < 0.003",
           np.isfinite(alpha_SN23) and abs(alpha_SN23-ALPHA_MAX_REF) < 0.003,
           f"alpha_SN23={alpha_SN23:.5f}" if np.isfinite(alpha_SN23) else "brak")

    record("P5: |alpha_SN23 - alpha_zero_1c| < 0.003",
           np.isfinite(alpha_SN23) and np.isfinite(alpha_zero_1c) and
           abs(alpha_SN23-alpha_zero_1c) < 0.003,
           f"SN23={alpha_SN23:.5f}, 1c={alpha_zero_1c:.5f}"
           if np.isfinite(alpha_SN23) and np.isfinite(alpha_zero_1c) else "brak jednego")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p89_codim2_verify.py")
    print("=" * 72)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    print(f"""
  alpha_SN23   = {alpha_SN23:.5f if np.isfinite(alpha_SN23) else 'brak'}
  alpha_zero_1c = {alpha_zero_1c:.5f if np.isfinite(alpha_zero_1c) else 'brak'}
  alpha_max_ref = {ALPHA_MAX_REF}
  |SN23 - max|  = {abs(alpha_SN23-ALPHA_MAX_REF):.5f if np.isfinite(alpha_SN23) else 'brak'}
  |SN23 - 1c|   = {abs(alpha_SN23-alpha_zero_1c):.5f if np.isfinite(alpha_SN23) and np.isfinite(alpha_zero_1c) else 'brak'}
""")

    if np.isfinite(alpha_SN23) and np.isfinite(alpha_zero_1c):
        if abs(alpha_SN23-ALPHA_MAX_REF)<0.003 and abs(alpha_SN23-alpha_zero_1c)<0.003:
            print("  WNIOSEK: KODIMENSION-2 POTWIERDZONE.")
            print("  K*_2=K*_3 i K*_1=K_nan_left zachodza przy tym samym alpha_max.")
        elif abs(alpha_SN23-ALPHA_MAX_REF)<0.010:
            print(f"  WNIOSEK: QUASI-KODIMENSION-2: alpha_SN23={alpha_SN23:.5f}")
            print(f"  blisko alpha_max={ALPHA_MAX_REF} (odch={abs(alpha_SN23-ALPHA_MAX_REF):.5f}).")
            print("  Wymaga wiekszej precyzji lub inne alpha_max.")
        else:
            print(f"  WNIOSEK: alpha_SN23={alpha_SN23:.5f} != alpha_max={ALPHA_MAX_REF}.")
            print("  Brak kodimension-2 przy tej precyzji.")
    elif n_23 < 3:
        print("  WNIOSEK: Za malo punktow K*_2/K*_3. Sprawdz K_arr lub range alpha.")
    else:
        print("  WNIOSEK: Niewystarczajace dane do oceny.")


if __name__ == '__main__':
    run_main()
