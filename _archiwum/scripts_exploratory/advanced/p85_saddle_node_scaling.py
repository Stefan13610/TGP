#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p85_saddle_node_scaling.py  --  TGP v1  *  Skalowanie bifurkacji saddle-node
=============================================================================
Cel: Weryfikacja mechanizmu siodlo-wezel przy alpha_max=8.8734.
     Blisko bifurkacji powinno zachodzic:

       K*_1(alpha) ~ K_c - A * (alpha_max - alpha)^beta

     beta=1/2 -> klasyczne siodlo-wezel (norma w bifurkacji 1D).
     Jesli potwierdzony -> alpha_max jest analitycznie wyznaczalny
     z warunku dg/dK = 0  (OP-3 Sciezka 3 podpora numeryczna).

Metoda:
  A. Skan alpha in [8.70, 8.87] (25 pkt) -- K*_1 dokladnie (n_K=40, n_psi=150)
  B. Fit potegowy: K*_1(alpha) = K_c - A*(alpha_max - alpha)^beta
     (nielinowy fit, 4 parametry: K_c, A, beta, alpha_max_free)
  C. Skan a_Gam in {0.020, 0.030, 0.040, 0.050, 0.060}
     -> czy beta jest UNIWERSALNE (niezalezne od a_Gam)?

Testy:
  P1: K*_1(alpha) obliczone dla >= 15 wartosci w [8.70, 8.87]
  P2: Fit potegowy zbiezny (RMSE < 1%)
  P3: beta in [0.40, 0.60] (saddle-node region)
  P4: |beta - 0.5| < 0.05 (klasyczne saddle-node)
  P5: |alpha_max_fit - alpha_max_ref| < 0.01

Optymalizacja: multiprocessing.Pool -- wymaga if __name__=='__main__' na Windows
Data: 2026-03-25
"""

import sys
import io
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_MAX_REF = 8.8734
ALPHA_K_REF   = 8.5616
A_GAM_DEFAULT = 0.040
LAM           = 5.501357e-06
GAMMA_V       = 1.0

V1 = GAMMA_V/3 - GAMMA_V/4


def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6

def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5


def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]


def _phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=2000):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0]-1e-5
    ev.terminal = True; ev.direction = -1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0, -1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except:
        return np.nan


def energy_at_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, 1500)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek + Ep
    except:
        return np.nan


def compute_g(K, a_gam, alpha, n_psi=150, r_max=40.0):
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1 for p in psis]
    branches = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: _phi_at_rmax(p, K, a_gam, alpha, r_max)-1,
                            psis[i], psis[i+1], xtol=1e-7, rtol=1e-7, maxiter=60)
                E = energy_at_psi(pz, K, a_gam, alpha, r_max)
                if np.isfinite(E):
                    branches.append(E/(4*np.pi*K)-1)
            except:
                pass
    if not branches:
        return np.nan
    return min(branches, key=abs)


def find_Kstar(alpha, a_gam, K_min=3e-3, K_max=0.15, n_K=30, n_psi=80, r_max=40.0):
    """Zwraca K*_1 lub nan. n_K=30, n_psi=80 -- balans szybkosc/dokladnosc."""
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = [compute_g(K, a_gam, alpha, n_psi, r_max) for K in K_scan]
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                return brentq(lambda K: compute_g(K, a_gam, alpha, n_psi, r_max),
                              K_scan[i], K_scan[i+1],
                              xtol=K_scan[i]*5e-4, rtol=5e-4, maxiter=30)
            except:
                pass
    return np.nan


def worker_Kstar(args):
    """Worker dla Pool.map — (alpha, a_gam, idx, total) -> (alpha, a_gam, K)."""
    alpha, a_gam, idx, total = args
    t0 = time.time()
    K = find_Kstar(alpha, a_gam)
    dt = time.time()-t0
    tag = f'[{idx+1:3d}/{total}] alpha={alpha:.4f}  a_Gam={a_gam:.3f}'
    if np.isfinite(K):
        print(f'  {tag}  ->  K*={K:.6f}  t={dt:.1f}s', flush=True)
    else:
        print(f'  {tag}  ->  brak   t={dt:.1f}s', flush=True)
    return (alpha, a_gam, K)


# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    PASS_COUNT = 0; FAIL_COUNT = 0
    results = {}

    def record(label, ok, info=""):
        nonlocal PASS_COUNT, FAIL_COUNT
        if ok: PASS_COUNT += 1
        else:  FAIL_COUNT += 1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # UWAGA: kazdy worker to osobny Python + scipy/numpy (~150-200 MB RAM).
    # 32 workerow = crash Windowsa. Bezpieczny limit: 8.
    N_WORKERS = min(8, cpu_count())

    print("=" * 72)
    print("TGP v1  *  p85_saddle_node_scaling.py  --  Bifurkacja saddle-node")
    print("  K*_1(alpha) ~ K_c - A*(alpha_max - alpha)^beta  [beta=0.5?]")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF}")
    print(f"  alpha_K_ref   = {ALPHA_K_REF}")
    print(f"  a_Gam default = {A_GAM_DEFAULT}")
    print(f"  CPU logical   = {cpu_count()},  Workers = {N_WORKERS}\n")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA A: Precyzyjny skan alpha in [8.70, 8.87] dla a_Gam=0.040
    # ─────────────────────────────────────────────────────────────────────────
    print("=" * 72)
    print("SEKCJA A: K*_1(alpha) dla alpha in [8.70, 8.87], a_Gam=0.040")
    print(f"  n_K=30, n_psi=80 | Workers={N_WORKERS} (bezpieczny limit RAM)")
    print("=" * 72)

    alpha_A = np.linspace(8.70, 8.873, 25)
    args_A = [(a, A_GAM_DEFAULT, i, len(alpha_A)) for i, a in enumerate(alpha_A)]

    t_A0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        res_A_raw = pool.map(worker_Kstar, args_A, chunksize=1)
    t_A = time.time()-t_A0

    res_A = [(a, K) for a, _, K in res_A_raw if np.isfinite(K)]
    n_valid_A = len(res_A)
    print(f"\n  Sekcja A: {t_A:.1f}s  |  K*_1 znaleziony dla {n_valid_A}/{len(alpha_A)} pkt")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA B: Fit potegowy
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA B: Fit potegowy K*_1(alpha) ~ K_c - A*(alpha_max - alpha)^beta")
    print("=" * 72)

    a_valid = np.array([x[0] for x in res_A])
    K_valid = np.array([x[1] for x in res_A])

    # Tabela
    print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'alpha_max-alpha':>16}")
    print("  " + "-" * 40)
    for a, K in res_A:
        print(f"  {a:>8.4f}  {K:>10.6f}  {ALPHA_MAX_REF-a:>16.5f}")

    free_fit = None
    best_fit = None
    best_rss = np.inf

    # Fit z ustalonym alpha_max
    print(f"\n  --- Fit z ustalonym alpha_max = {ALPHA_MAX_REF} ---")
    x_data = ALPHA_MAX_REF - a_valid
    mask = x_data > 1e-5
    if mask.sum() >= 5:
        xf = x_data[mask]; Kf = K_valid[mask]
        def model3(x, Kc, A, beta): return Kc - A*x**beta
        try:
            p0 = [Kf[-1]+0.001, 0.003, 0.5]
            popt, _ = curve_fit(model3, xf, Kf, p0=p0,
                                bounds=([Kf.max(), 0, 0.1], [Kf.max()+0.005, 1, 2]),
                                maxfev=5000)
            pred = model3(xf, *popt)
            rmse = np.sqrt(np.mean((pred-Kf)**2))/Kf.mean()*100
            print(f"    K_c  = {popt[0]:.7f}")
            print(f"    A    = {popt[1]:.6f}")
            print(f"    beta = {popt[2]:.5f}")
            print(f"    RMSE = {rmse:.4f}%")
            best_fit = (ALPHA_MAX_REF, popt[0], popt[1], popt[2], rmse)
        except Exception as e:
            print(f"    Fail: {e}")

    # Fit ze swobodnym alpha_max
    print(f"\n  --- Fit ze swobodnym alpha_max ---")
    if n_valid_A >= 8:
        def model4(alpha_arr, Kc, A, beta, amax):
            x = np.maximum(amax - alpha_arr, 1e-10)
            return Kc - A * x**beta
        try:
            p0 = [0.012, 0.005, 0.5, 8.875]
            bounds = ([0.010, 1e-5, 0.1, 8.85], [0.016, 2.0, 2.0, 8.92])
            popt4, pcov4 = curve_fit(model4, a_valid, K_valid, p0=p0,
                                     bounds=bounds, maxfev=10000)
            pred4 = model4(a_valid, *popt4)
            rmse4 = np.sqrt(np.mean((pred4-K_valid)**2))/K_valid.mean()*100
            perr4 = np.sqrt(np.diag(pcov4))
            print(f"    K_c       = {popt4[0]:.7f}  ± {perr4[0]:.7f}")
            print(f"    A         = {popt4[1]:.6f}  ± {perr4[1]:.6f}")
            print(f"    beta      = {popt4[2]:.5f}  ± {perr4[2]:.5f}")
            print(f"    alpha_max = {popt4[3]:.5f}  ± {perr4[3]:.5f}")
            print(f"    RMSE      = {rmse4:.4f}%")
            print(f"    |alpha_max_fit - ref| = {abs(popt4[3]-ALPHA_MAX_REF):.5f}")
            free_fit = popt4
            if rmse4 < best_rss:
                best_rss = rmse4
        except Exception as e:
            print(f"    Fail: {e}")

    # ─────────────────────────────────────────────────────────────────────────
    # SEKCJA C: Universalnosc beta wzgledem a_Gam
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA C: Universalnosc beta(a_Gam)")
    print("  alpha_max(a_Gam) ~ 9.023 - 3.75*a_Gam  (z p81)")
    print("=" * 72)

    a_gam_scan = [0.020, 0.030, 0.040, 0.050, 0.060]
    a_gam_amax = {ag: 9.023 - 3.75*ag for ag in a_gam_scan}
    N_PER_AGAM = 8

    args_C = []
    for ag in a_gam_scan:
        amax_ag = a_gam_amax[ag]
        alphas_c = np.linspace(amax_ag - 0.18, amax_ag - 0.01, N_PER_AGAM)
        for a in alphas_c:
            args_C.append((a, ag))
    total_C = len(args_C)
    args_C = [(a, ag, i, total_C) for i, (a, ag) in enumerate(args_C)]

    print(f"\n  Lacznie {total_C} obliczen ({len(a_gam_scan)} x {N_PER_AGAM})")
    t_C0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=3) as pool:
        res_C_raw = pool.map(worker_Kstar, args_C, chunksize=1)
    t_C = time.time()-t_C0
    print(f"\n  Sekcja C: {t_C:.1f}s")

    print(f"\n  {'a_Gam':>7}  {'alpha_max':>10}  {'beta_fit':>9}  "
          f"{'K_c':>10}  {'RMSE%':>7}  {'n':>4}")
    print("  " + "-" * 54)

    betas_C = []
    for ag in a_gam_scan:
        amax_ag = a_gam_amax[ag]
        pts = [(a, K) for a, ag2, K in res_C_raw if ag2 == ag and np.isfinite(K)]
        if len(pts) < 4:
            print(f"  {ag:>7.3f}  {amax_ag:>10.4f}  za malo pkt ({len(pts)})")
            continue
        ac = np.array([p[0] for p in pts])
        Kc = np.array([p[1] for p in pts])
        xc = amax_ag - ac
        mask = xc > 1e-5
        if mask.sum() < 4:
            continue
        try:
            def m3(x, Kc_, A, beta): return Kc_ - A*x**beta
            popt, _ = curve_fit(m3, xc[mask], Kc[mask],
                                p0=[Kc[mask][-1]+0.0005, 0.003, 0.5],
                                bounds=([Kc.max()-0.0001, 0, 0.1],
                                        [Kc.max()+0.005, 2, 2]),
                                maxfev=3000)
            pred = m3(xc[mask], *popt)
            rmse = np.sqrt(np.mean((pred-Kc[mask])**2))/Kc[mask].mean()*100
            print(f"  {ag:>7.3f}  {amax_ag:>10.4f}  {popt[2]:>9.4f}  "
                  f"{popt[0]:>10.7f}  {rmse:>7.4f}%  {mask.sum():>4}")
            betas_C.append(popt[2])
        except Exception as e:
            print(f"  {ag:>7.3f}  fit fail: {e}")

    beta_mean = np.mean(betas_C) if betas_C else np.nan
    beta_std  = np.std(betas_C)  if len(betas_C) > 1 else np.nan
    if betas_C:
        print(f"\n  beta srednia = {beta_mean:.4f}  ±  {beta_std:.4f}")
        print(f"  |beta - 0.5| = {abs(beta_mean-0.5):.4f}")

    # ─────────────────────────────────────────────────────────────────────────
    # TESTY
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TESTY")
    print("=" * 72)

    record("P1: K*_1 obliczone dla >= 15 wartosci w [8.70, 8.87]",
           n_valid_A >= 15, f"{n_valid_A} pkt")

    rmse_best = best_fit[4] if best_fit else np.inf
    record("P2: Fit potegowy zbiezny (RMSE < 1%)",
           rmse_best < 1.0, f"RMSE={rmse_best:.4f}%" if best_fit else "brak")

    beta_best = (free_fit[2] if free_fit is not None
                 else (best_fit[3] if best_fit else np.nan))
    record("P3: beta in [0.40, 0.60]",
           np.isfinite(beta_best) and 0.40 <= beta_best <= 0.60,
           f"beta={beta_best:.4f}" if np.isfinite(beta_best) else "brak")

    record("P4: |beta - 0.5| < 0.05 (klasyczne saddle-node)",
           np.isfinite(beta_best) and abs(beta_best-0.5) < 0.05,
           f"beta={beta_best:.4f}" if np.isfinite(beta_best) else "brak")

    amax_fit = free_fit[3] if free_fit is not None else np.nan
    record("P5: |alpha_max_fit - alpha_max_ref| < 0.01",
           np.isfinite(amax_fit) and abs(amax_fit-ALPHA_MAX_REF) < 0.01,
           f"fit={amax_fit:.5f}  odch={abs(amax_fit-ALPHA_MAX_REF):.5f}"
           if np.isfinite(amax_fit) else "brak")

    # ─────────────────────────────────────────────────────────────────────────
    # PODSUMOWANIE
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p85_saddle_node_scaling.py")
    print("=" * 72)
    print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
    print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")
    print(f"""
  WYNIKI:
    n_valid_A     = {n_valid_A}
    beta_best     = {beta_best:.4f}   (target: 0.5000)
    beta_C_mean   = {beta_mean:.4f}   ± {beta_std:.4f}  (z {len(betas_C)} a_Gam)
    alpha_max_fit = {amax_fit:.5f}  (ref: {ALPHA_MAX_REF})
""")

    if np.isfinite(beta_best) and abs(beta_best-0.5) < 0.05:
        print(f"""  WNIOSEK P4 PASS:
  beta = {beta_best:.4f} ~ 0.5  ->  KLASYCZNA BIFURKACJA SIODLO-WEZEL
  K*_1(alpha) = K_c - A*(alpha_max - alpha)^(1/2)
  Implikacja (OP-3 Sciezka 3):
    Warunek dg/dK=0 wyznacza alpha_max analitycznie z V_mod.
    alpha_K = alpha_max * n_s byloby predykcja substratowa.
""")
        if np.isfinite(beta_mean) and abs(beta_mean-0.5) < 0.05:
            print(f"  BONUS: beta universalne (± {beta_std:.4f}) -- "
                  f"niezalezne od a_Gam!\n")
    elif np.isfinite(beta_best) and 0.40 <= beta_best <= 0.60:
        print(f"""  WNIOSEK P3 PASS (nie P4):
  beta = {beta_best:.4f} -- blisko 0.5 ale poza tolerancja 0.05.
  Moze saddle-node z poprawkami wyzszego rzedu.
""")
    else:
        print(f"""  WNIOSEK: beta = {beta_best:.4f} -- nie klasyczne saddle-node.
""")


if __name__ == '__main__':
    run_main()
