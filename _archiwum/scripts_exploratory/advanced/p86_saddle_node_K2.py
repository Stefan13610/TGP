#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p86_saddle_node_K2.py  --  TGP v1  *  Diagnostyka g(K): zeros i topologia
==========================================================================
Motywacja (wynik p85 + bledny p86 v1):
  Nie wiadomo czy K*_2 - K*_1 jest wlasciwym parametrem porzadku.
  Najpierw trzeba zrozumiec topologie g(K) dla roznych alpha.

Cel:
  Skan g(K) dla K in [0.005, 0.13], 40 pkt, dla:
    alpha in {8.50, 8.60, 8.70, 8.75, 8.80, 8.83, 8.85, 8.87}
  -> ile zer g(K) istnieje przy kazdym alpha?
  -> czy zera schodza sie przy alpha -> alpha_max?
  -> jaka jest topologia: K*_1 i K*_2 czy tylko K*_1?

Ta diagnostyka odpowie na pytanie: jaki jest WLASCIWY parametr porzadku
bifurkacji saddle-node, i jak go mierzyc.

Testy:
  P1: g(K) obliczone dla >= 6 wartosci alpha
  P2: K*_1 zgadza sie z p85/p79 (odch < 1%)
  P3: K*_2 znalezione dla przynajmniej jednej alpha
  P4: Jesli K*_2 istnieje: K*_2 < K*_1_prev przy alpha_max (tj. zblizenie)
  P5: Narysowanie schematycznej mapy zer g(K,alpha)

Workers: 8 (kazdy alpha to oddzielny skan g, rownolegly)
Data: 2026-03-25
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from multiprocessing import Pool, cpu_count
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_MAX_REF = 8.8734
A_GAM         = 0.040
LAM           = 5.501357e-06
GAMMA_V       = 1.0
V1            = GAMMA_V/3 - GAMMA_V/4
K1_REF        = 0.010414   # z p10

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

def g_at_K(K, a_gam, alpha, n_psi=80, r_max=40.0):
    """g(K) = E/(4*pi*K) - 1 dla galezie z minimalna |g|. Lub nan."""
    psi_max = max(3.0, 1+2*K/a_gam)
    psis = np.linspace(1.001, min(psi_max, 300), n_psi)
    Fv = [phi_at_rmax(p,K,a_gam,alpha,r_max)-1 for p in psis]
    best = None
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,a_gam,alpha,r_max)-1,
                            psis[i],psis[i+1], xtol=1e-7, maxiter=50)
                E = energy_at_psi(pz,K,a_gam,alpha,r_max)
                if np.isfinite(E):
                    g_val = E/(4*np.pi*K)-1
                    if best is None or abs(g_val)<abs(best[0]):
                        best = (g_val, pz, E)
            except: pass
    return best[0] if best else np.nan

def scan_g_for_alpha(args):
    """Worker: skan g(K) dla jednej wartosci alpha."""
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()
    g_vals = [g_at_K(K, a_gam, alpha) for K in K_arr]
    dt = time.time()-t0

    # Znajdz wszystkie zera g(K) (zmiany znaku)
    zeros = []
    for i in range(len(g_vals)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                Kz = brentq(lambda K: g_at_K(K, a_gam, alpha),
                            K_arr[i], K_arr[i+1],
                            xtol=K_arr[i]*5e-4, rtol=1e-4, maxiter=25)
                zeros.append(Kz)
            except: pass

    n_nan = sum(1 for g in g_vals if np.isnan(g))
    tag = f'[{idx+1}/{total}] alpha={alpha:.4f}'
    zeros_str = ', '.join(f'{z:.6f}' for z in zeros) if zeros else 'brak'
    print(f'  {tag}  zeros=[{zeros_str}]  n_nan={n_nan}  t={dt:.1f}s', flush=True)
    return (alpha, np.array(g_vals), zeros, n_nan)


def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    PC = 0; FC = 0
    N_WORKERS = min(8, cpu_count())

    def record(label, ok, info=""):
        nonlocal PC, FC
        if ok: PC += 1
        else:  FC += 1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    print("=" * 72)
    print("TGP v1  *  p86_saddle_node_K2.py  --  Diagnostyka g(K): topologia")
    print("  Ile zer g(K) istnieje? Czy K*_1 i K*_2 schodza sie?")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF}")
    print(f"  K1_ref = {K1_REF}  (z p10)")
    print(f"  Workers = {N_WORKERS}\n")

    # K scan: [0.005, 0.13], 40 pkt (logarytmicznie)
    # + zageszczenie w [0.009, 0.015] gdzie K*_1 jest spodziewane
    K_coarse = np.logspace(np.log10(0.005), np.log10(0.13), 35)
    K_fine   = np.linspace(0.009, 0.016, 12)
    K_arr    = np.unique(np.concatenate([K_coarse, K_fine]))
    K_arr.sort()

    print(f"  K scan: {len(K_arr)} punktow w [{K_arr[0]:.4f}, {K_arr[-1]:.4f}]")

    # Alpha do zbadania: szeroko + gesto przy alpha_max
    alpha_scan = [8.50, 8.60, 8.70, 8.75, 8.80, 8.83, 8.85, 8.87]
    args_list = [(a, A_GAM, K_arr, i, len(alpha_scan))
                 for i, a in enumerate(alpha_scan)]

    print(f"  alpha scan: {alpha_scan}\n")

    # ─────────────────────────────────────────────────────────────────────────
    print("=" * 72)
    print("SEKCJA A: g(K) scan (rownolegly po alpha)")
    print("=" * 72)

    t0 = time.time()
    with Pool(N_WORKERS, maxtasksperchild=2) as pool:
        results = pool.map(scan_g_for_alpha, args_list, chunksize=1)
    t_tot = time.time()-t0
    print(f"\n  Lacznie: {t_tot:.1f}s")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA B: Tabela zer g(K,alpha) -- mapa bifurkacji")
    print("=" * 72)

    # Sortuj po alpha
    results.sort(key=lambda x: x[0])

    # Mapa zer
    print(f"\n  {'alpha':>7}  {'n_zer':>6}  {'K*_1':>10}  {'K*_2':>10}  "
          f"{'K*_3':>10}  {'n_nan':>6}")
    print("  " + "-" * 58)

    K1_vals = {}; K2_vals = {}
    all_K1_ok = True
    K2_exists = False

    for alpha, g_arr, zeros, n_nan in results:
        zeros_sorted = sorted(zeros)
        K1 = zeros_sorted[0] if len(zeros_sorted) >= 1 else np.nan
        K2 = zeros_sorted[1] if len(zeros_sorted) >= 2 else np.nan
        K3 = zeros_sorted[2] if len(zeros_sorted) >= 3 else np.nan
        K1_vals[alpha] = K1
        K2_vals[alpha] = K2
        K2_exists = K2_exists or np.isfinite(K2)

        K2_str = f'{K2:.6f}' if np.isfinite(K2) else '  brak  '
        K3_str = f'{K3:.6f}' if np.isfinite(K3) else '  brak  '
        print(f"  {alpha:>7.4f}  {len(zeros):>6}  {K1:.6f}  {K2_str}  {K3_str}  {n_nan:>6}")

        if np.isfinite(K1) and abs(K1 - K1_REF)/K1_REF > 0.02:
            all_K1_ok = False

    # Anomalie w g(K) -- obszary nan (zapad profilu)
    print("\n  Obszary nan w g(K) -- K_collapse:")
    for alpha, g_arr, zeros, n_nan in results:
        nan_K = [K_arr[i] for i in range(len(g_arr)) if np.isnan(g_arr[i])]
        if nan_K:
            print(f"  alpha={alpha:.4f}: g=nan dla K in [{min(nan_K):.5f}, {max(nan_K):.5f}]"
                  f"  (n={len(nan_K)})")
        else:
            print(f"  alpha={alpha:.4f}: brak nan w g(K)")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("SEKCJA C: Analiza zer -- czy schodza sie przy alpha_max?")
    print("=" * 72)

    # K*_1(alpha)
    print("\n  K*_1(alpha) -- porowniez z p85:")
    print(f"  {'alpha':>7}  {'K*_1':>10}  {'p85_ref':>10}  {'delta%':>8}")
    p85_ref = {8.70: 0.010491, 8.75: 0.010525, 8.80: 0.010632,
               8.83: 0.010700, 8.85: 0.010851}
    for alpha, g_arr, zeros, n_nan in results:
        K1 = K1_vals.get(alpha, np.nan)
        p85 = p85_ref.get(alpha, np.nan)
        if np.isfinite(K1) and np.isfinite(p85):
            drel = (K1-p85)/p85*100
            print(f"  {alpha:>7.4f}  {K1:>10.6f}  {p85:>10.6f}  {drel:>8.4f}%")
        elif np.isfinite(K1):
            print(f"  {alpha:>7.4f}  {K1:>10.6f}  {'brak ref':>10}")
        else:
            print(f"  {alpha:>7.4f}  {'brak':>10}")

    # Zblizanie K*_2 do K*_1?
    if K2_exists:
        print("\n  K*_2 - K*_1  (zblizanie = saddle-node):")
        for alpha, g_arr, zeros, n_nan in results:
            K1 = K1_vals.get(alpha, np.nan)
            K2 = K2_vals.get(alpha, np.nan)
            if np.isfinite(K1) and np.isfinite(K2):
                print(f"  alpha={alpha:.4f}  K*_2-K*_1 = {K2-K1:.6f}  "
                      f"ratio = {K2/K1:.4f}")
    else:
        print("\n  K*_2 nie znalezione w zadnym alpha -- sprawdz K_max lub psi_max.")
        print("  Mozliwe: K*_2 >> 0.13 (poza skanem) lub nie istnieje przy tych alpha.")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TESTY")
    print("=" * 72)

    n_alpha_valid = sum(1 for alpha in K1_vals if np.isfinite(K1_vals[alpha]))
    record("P1: g(K) obliczone dla >= 6 wartosci alpha", n_alpha_valid >= 6,
           f"{n_alpha_valid}/{len(alpha_scan)}")

    record("P2: K*_1 zgodne z p85 (odch < 2%)", all_K1_ok,
           "TAK" if all_K1_ok else "NIESPOJNOSC branch")

    record("P3: K*_2 znalezione dla >= 1 alpha", K2_exists,
           "TAK" if K2_exists else "NIE -- K*_2 > K_max=0.13 lub nieistnieje")

    # P4: czy K*_2 i K*_1 schodza sie (ratio maleje)
    K2_K1_ratios = [(a, K2_vals[a]/K1_vals[a])
                    for a in K1_vals
                    if np.isfinite(K2_vals.get(a,np.nan)) and np.isfinite(K1_vals[a])]
    K2_converging = (len(K2_K1_ratios) >= 2 and
                     K2_K1_ratios[-1][1] < K2_K1_ratios[0][1])
    record("P4: K*_2/K*_1 maleje z alpha (zblizanie do saddle-node)",
           K2_converging,
           f"{[(a, f'{r:.3f}') for a,r in K2_K1_ratios]}" if K2_K1_ratios else "brak K*_2")

    # P5: Sprawdz czy nan-region w g(K) wskazuje K_collapse ~ K*_1 przy alpha_max
    nan_near_K1 = False
    for alpha, g_arr, zeros, n_nan in results:
        if alpha >= 8.85 and n_nan > 0:
            nan_K = [K_arr[i] for i in range(len(g_arr)) if np.isnan(g_arr[i])]
            K1 = K1_vals.get(alpha, np.nan)
            if np.isfinite(K1) and nan_K and min(nan_K) < K1*1.5:
                nan_near_K1 = True
    record("P5: nan-region (zapad) pojawia sie blisko K*_1 przy alpha >= 8.85",
           nan_near_K1, "TAK" if nan_near_K1 else "NIE (K_collapse daleko od K*_1)")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p86_saddle_node_K2.py  (diagnostyka g(K))")
    print("=" * 72)
    print(f"\n  PASS: {PC} / {PC+FC}")
    print(f"  FAIL: {FC} / {PC+FC}")

    print(f"""
  WNIOSKI:
    K*_2 znalezione: {'TAK' if K2_exists else 'NIE'}
    K*_1 spojna z p85: {'TAK' if all_K1_ok else 'NIE -- inna galaz!'}
    nan-region blisko K*_1: {'TAK' if nan_near_K1 else 'NIE'}

  INTERPRETACJA:
""")
    if not K2_exists:
        print("  K*_2 lezy poza skanem K in [0.005, 0.13]. Mechanizm bifurkacji")
        print("  alfa_max jest INNY niz zblizanie K*_1/K*_2 -- prawdopodobnie")
        print("  K_collapse (zapad profilu) zbliaza sie do K*_1.")
        print("  Proponowane: zeskanuj K in [K*_1, K*_1*10] szukajac K*_2.")
    elif K2_converging:
        print("  K*_2 - K*_1 maleje z alpha: SADDLE-NODE K*_1/K*_2 potwierdzone.")
        print("  Nastepny krok: precyzyjny fit delta_K ~ (alpha_max-alpha)^gamma.")
    else:
        print("  K*_2 istnieje ale nie zbliza sie do K*_1. Mechanizm nieznany.")


if __name__ == '__main__':
    run_main()
