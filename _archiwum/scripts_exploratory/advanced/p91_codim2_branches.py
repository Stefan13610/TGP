#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p91_codim2_branches.py  --  TGP v1  *  Weryfikacja kodimension-2 (czyste galęzie)
==================================================================================
Diagnoza z p89/p90:

  g_at_K (min|g|) MIESZA galezi psi_0:
    - Galaz L: pierwszy zero phi(rmax)-1 w [1.001,3.0] -> daje zero g_L w K*_A (wiekszy K)
    - Galaz U: drugi  zero phi(rmax)-1 w [1.001,3.0] -> daje zero g_U w K*_B (mniejszy K)
  Pseudo-zero "K*_2" to ARTEFAKT przelaczenia galezi min|g|, nie prawdziwy soliton.

  Prawdziwe zera samozgodnosci (fizyczne solitony):
    K*_B ~ 0.011  (zero g_U, galaz gorna psi_0)
    K*_A ~ 0.0135 (zero g_L, galaz dolna psi_0)

  Prawa granica istnienia solitonow:
    K_nan_R = najmniejsze K > K*_A gdzie obie galeze daja brak (nan-region)

  Hipoteza kodimension-2 (zrewidowana):
    (C1) K*_B -> K*_A  przy alpha -> alpha_SN (saddle-node obu solitonow)
    (C2) K*_A -> K_nan_R  przy alpha -> alpha_nan (granica zanikania)
    Kodimension-2: alpha_SN = alpha_nan = alpha_max = 8.8734?

Testy:
  P1: K*_A i K*_B znalezione dla >= 5 alpha
  P2: delta_AB = K*_A - K*_B maleje ku 0 (solitony sie zbiegaja)
  P3: delta_An = K_nan_R - K*_A maleje ku 0 (K*_A dotyka granicy)
  P4: Oba extrapolowane zera bliskie alpha_max
  P5: |alpha_SN - alpha_nan| < 0.003

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

def phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0, n_eval=1400):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0] - 1e-5
    ev.terminal = True; ev.direction = -1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0, -1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except: return np.nan

def energy_at_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, 1400)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek + Ep
    except: return np.nan


def g_branches_at_K(K, a_gam, alpha, n_psi=80):
    """
    Dla danego K znajdz g na DWOCH galęziach:
      g_L = g dla PIERWSZEGO (najnizszego psi_0) zera phi(rmax)-1
      g_U = g dla DRUGIEGO zera phi(rmax)-1
    Zwraca (g_L, psi_L, g_U, psi_U).
    Jesli danej galezie brak: (nan, nan, ...).
    """
    psis = np.linspace(1.001, 3.0, n_psi)
    Fv = [phi_at_rmax(p, K, a_gam, alpha) - 1.0 for p in psis]

    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, a_gam, alpha) - 1.0,
                            psis[i], psis[i+1], xtol=1e-8, maxiter=60)
                E = energy_at_psi(pz, K, a_gam, alpha)
                if np.isfinite(E):
                    gv = E / (4*np.pi*K) - 1.0
                    roots.append((gv, pz))
            except: pass

    g_L  = roots[0][0] if len(roots) >= 1 else np.nan
    psi_L = roots[0][1] if len(roots) >= 1 else np.nan
    g_U  = roots[1][0] if len(roots) >= 2 else np.nan
    psi_U = roots[1][1] if len(roots) >= 2 else np.nan
    return g_L, psi_L, g_U, psi_U


def find_zero_brentq(branch, K_lo, K_hi, a_gam, alpha, xtol_rel=2e-4):
    """Szuka zera g_branch(K) miedzy K_lo a K_hi (branch = 'L' lub 'U')."""
    def gfunc(K):
        gL, _, gU, _ = g_branches_at_K(K, a_gam, alpha)
        return gL if branch == 'L' else gU

    glo = gfunc(K_lo); ghi = gfunc(K_hi)
    if not (np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0):
        return np.nan
    try:
        return brentq(gfunc, K_lo, K_hi, xtol=K_lo*xtol_rel, rtol=xtol_rel, maxiter=30)
    except: return np.nan


def worker_scan(args):
    alpha, a_gam, K_arr, idx, total = args
    t0 = time.time()

    # Oblicz g_L i g_U dla wszystkich K
    gL_arr = []; gU_arr = []
    psiL_arr = []; psiU_arr = []
    for K in K_arr:
        gL, pL, gU, pU = g_branches_at_K(K, a_gam, alpha)
        gL_arr.append(gL); gU_arr.append(gU)
        psiL_arr.append(pL); psiU_arr.append(pU)

    gL_arr = np.array(gL_arr)
    gU_arr = np.array(gU_arr)

    # --- Znajdz zera g_L (da K*_A) i g_U (da K*_B) ---
    zeros_L = []
    zeros_U = []

    for i in range(len(K_arr)-1):
        # Galaz L
        if np.isfinite(gL_arr[i]) and np.isfinite(gL_arr[i+1]) and gL_arr[i]*gL_arr[i+1] < 0:
            Kz = find_zero_brentq('L', K_arr[i], K_arr[i+1], a_gam, alpha)
            if np.isfinite(Kz): zeros_L.append(Kz)
        # Galaz U
        if np.isfinite(gU_arr[i]) and np.isfinite(gU_arr[i+1]) and gU_arr[i]*gU_arr[i+1] < 0:
            Kz = find_zero_brentq('U', K_arr[i], K_arr[i+1], a_gam, alpha)
            if np.isfinite(Kz): zeros_U.append(Kz)

    # Najwieksze zero g_U (bliskie alpha_max, K*_B ~ 0.011)
    KstB = zeros_U[-1] if zeros_U else np.nan  # ostatnie zero g_U (najblizsze K*_A)

    # Najmniejsze zero g_L (bliskie alpha_max, K*_A ~ 0.013)
    KstA = zeros_L[0] if zeros_L else np.nan  # pierwsze zero g_L (najblizsze K*_B)

    # K_nan_R: PRAWA granica obszaru solitonu (gdzie obie galeze daja brak)
    # Szukamy LEWEGO brzegu regionu gdzie g_L=nan AND g_U=nan (obie brak)
    K_nan_R = np.nan
    for i in range(len(K_arr)-1):
        L_ok = np.isfinite(gL_arr[i])
        U_ok = np.isfinite(gU_arr[i])
        L_ok_next = np.isfinite(gL_arr[i+1])
        U_ok_next = np.isfinite(gU_arr[i+1])
        # Szukamy: przynajmniej jedna galaz istnieje w K[i], obie brak w K[i+1]
        if (L_ok or U_ok) and not (L_ok_next or U_ok_next):
            # Doprecyzuj
            K_lo_b = K_arr[i]; K_hi_b = K_arr[i+1]
            for _ in range(12):
                K_mid = 0.5*(K_lo_b + K_hi_b)
                gLm, _, gUm, _ = g_branches_at_K(K_mid, a_gam, alpha)
                if np.isfinite(gLm) or np.isfinite(gUm):
                    K_lo_b = K_mid
                else:
                    K_hi_b = K_mid
                if K_hi_b - K_lo_b < K_lo_b * 3e-4: break
            K_nan_R = 0.5*(K_lo_b + K_hi_b)
            break

    # Delty diagnostyczne
    dAB = KstA - KstB  if np.isfinite(KstA) and np.isfinite(KstB) else np.nan
    dAn = K_nan_R - KstA if np.isfinite(KstA) and np.isfinite(K_nan_R) else np.nan

    # Liczniki
    nL_found = sum(1 for g in gL_arr if np.isfinite(g))
    nnan_both = sum(1 for gL, gU in zip(gL_arr, gU_arr)
                    if not np.isfinite(gL) and not np.isfinite(gU))

    dt = time.time() - t0

    # Stringi (FIX Python 3.14)
    KstAs  = f'{KstA:.6f}'  if np.isfinite(KstA)   else '  brak  '
    KstBs  = f'{KstB:.6f}'  if np.isfinite(KstB)   else '  brak  '
    KnanRs = f'{K_nan_R:.6f}' if np.isfinite(K_nan_R) else '  brak  '
    dABs   = f'{dAB:.6f}'  if np.isfinite(dAB)   else '  brak '
    dAns   = f'{dAn:.6f}'  if np.isfinite(dAn)   else '  brak '

    print(f'  [{idx+1:2d}/{total}] a={alpha:.4f}  '
          f'K*_B={KstBs}  K*_A={KstAs}  K_nan_R={KnanRs}  '
          f'd_AB={dABs}  d_An={dAns}  '
          f'nan_both={nnan_both}/{len(K_arr)}  t={dt:.1f}s', flush=True)

    return (alpha, KstA, KstB, K_nan_R)


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
    print("TGP v1  *  p91_codim2_branches.py  --  Kodimension-2 (galęzie)")
    print("  g_L = pierwsza galaz psi_0 -> daje K*_A (wiekszy K)")
    print("  g_U = druga   galaz psi_0 -> daje K*_B (mniejszy K)")
    print("  Hipoteza C1: K*_A -> K*_B przy alpha_SN")
    print("  Hipoteza C2: K*_A -> K_nan_R przy alpha_nan")
    print("  Kodimension-2: alpha_SN = alpha_nan = alpha_max?")
    print("=" * 72)
    print(f"\n  alpha_max_ref = {ALPHA_MAX_REF},  Workers = {N_WORKERS}\n")

    # K grid: obejmuje K*_B~0.010, K*_A~0.013, K_nan_R~0.014
    # Drobna siatka [0.009, 0.0160] z 85 pkt (krok ~8.4e-5)
    K_arr = np.linspace(0.009, 0.0160, 85)

    # Alpha: od 8.840 (gdzie K*_A i K*_B oba istnieja) do alpha_max
    alpha_base = np.arange(8.840, 8.874, 0.003)
    alpha_extra = np.array([8.855, 8.860, 8.863, 8.866, 8.870, ALPHA_MAX_REF])
    alpha_scan = np.unique(np.concatenate([alpha_base, alpha_extra]))
    alpha_scan.sort()

    args_list = [(a, A_GAM, K_arr, i, len(alpha_scan))
                 for i, a in enumerate(alpha_scan)]

    print(f"  K grid: {len(K_arr)} pkt w [{K_arr[0]:.4f}, {K_arr[-1]:.4f}]  "
          f"(krok={K_arr[1]-K_arr[0]:.5f})")
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
    print(f"\n  {'alpha':>8}  {'K*_B':>10}  {'K*_A':>10}  {'K_nan_R':>10}  "
          f"{'d_AB':>9}  {'d_An':>9}")
    print("  " + "-" * 75)

    pts_AB = []; pts_An = []

    for alpha, KstA, KstB, K_nan_R in results:
        dAB = KstA - KstB   if np.isfinite(KstA) and np.isfinite(KstB)   else np.nan
        dAn = K_nan_R - KstA if np.isfinite(KstA) and np.isfinite(K_nan_R) else np.nan

        KstAs  = f'{KstA:.6f}'   if np.isfinite(KstA)   else '  brak  '
        KstBs  = f'{KstB:.6f}'   if np.isfinite(KstB)   else '  brak  '
        KnanRs = f'{K_nan_R:.6f}' if np.isfinite(K_nan_R) else '  brak  '
        dABs   = f'{dAB:.6f}'   if np.isfinite(dAB)   else '  brak '
        dAns   = f'{dAn:.6f}'   if np.isfinite(dAn)   else '  brak '
        mrk = ' <--' if abs(alpha - ALPHA_MAX_REF) < 0.0001 else ''

        print(f"  {alpha:>8.4f}  {KstBs}  {KstAs}  {KnanRs}  "
              f"{dABs}  {dAns}{mrk}")

        if np.isfinite(dAB): pts_AB.append((alpha, dAB))
        if np.isfinite(dAn): pts_An.append((alpha, dAn))

    n_AB = len(pts_AB); n_An = len(pts_An)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("ANALIZA ZBIEZNOSCI")
    print("=" * 72)

    alpha_SN = np.nan; alpha_nan = np.nan

    if n_AB >= 3:
        avAB = np.array([x[0] for x in pts_AB])
        dvAB = np.array([x[1] for x in pts_AB])
        print(f"\n  Trend delta_AB (K*_A - K*_B):  [{n_AB} punktow]")
        for a, d in pts_AB:
            print(f"    alpha={a:.4f}  d_AB={d:.7f}")
        try:
            p = np.polyfit(avAB, dvAB, 1)
            alpha_SN = -p[1]/p[0]
            alpha_SN_s = f'{alpha_SN:.5f}'
            print(f"  Ekstrapolacja liniowa d_AB=0: alpha_SN = {alpha_SN_s}")
        except: pass
        if n_AB >= 4:
            def model_AB(a, C, g, amax):
                return C * np.maximum(amax - a, 1e-12)**g
            try:
                a_g = alpha_SN if np.isfinite(alpha_SN) else 8.875
                popt, pcov = curve_fit(model_AB, avAB, dvAB,
                                       p0=[0.003, 0.6, a_g],
                                       bounds=([0, 0.1, avAB[-1]],
                                               [1, 3, avAB[-1]+0.05]),
                                       maxfev=5000)
                pred = model_AB(avAB, *popt)
                rmse = np.sqrt(np.mean((pred-dvAB)**2))/dvAB.mean()*100
                perr = np.sqrt(np.diag(pcov))
                alpha_SN = popt[2]
                print(f"  Fit potegowy: gamma_AB={popt[1]:.4f}+/-{perr[1]:.4f}  "
                      f"alpha_SN={popt[2]:.5f}+/-{perr[2]:.5f}  RMSE={rmse:.3f}%")
            except Exception as e:
                print(f"  Fit potegowy d_AB fail: {e}")

    if n_An >= 3:
        avAn = np.array([x[0] for x in pts_An])
        dvAn = np.array([x[1] for x in pts_An])
        print(f"\n  Trend delta_An (K_nan_R - K*_A):  [{n_An} punktow]")
        for a, d in pts_An:
            print(f"    alpha={a:.4f}  d_An={d:.7f}")
        try:
            p = np.polyfit(avAn, dvAn, 1)
            alpha_nan = -p[1]/p[0]
            alpha_nan_s = f'{alpha_nan:.5f}'
            print(f"  Ekstrapolacja liniowa d_An=0: alpha_nan = {alpha_nan_s}")
        except: pass
        if n_An >= 4:
            def model_An(a, C, g, amax):
                return C * np.maximum(amax - a, 1e-12)**g
            try:
                a_g = alpha_nan if np.isfinite(alpha_nan) else 8.880
                popt, pcov = curve_fit(model_An, avAn, dvAn,
                                       p0=[0.001, 0.5, a_g],
                                       bounds=([0, 0.1, avAn[-1]],
                                               [1, 3, avAn[-1]+0.05]),
                                       maxfev=5000)
                pred = model_An(avAn, *popt)
                rmse = np.sqrt(np.mean((pred-dvAn)**2))/dvAn.mean()*100
                perr = np.sqrt(np.diag(pcov))
                alpha_nan = popt[2]
                print(f"  Fit potegowy: gamma_An={popt[1]:.4f}+/-{perr[1]:.4f}  "
                      f"alpha_nan={popt[2]:.5f}+/-{perr[2]:.5f}  RMSE={rmse:.3f}%")
            except Exception as e:
                print(f"  Fit potegowy d_An fail: {e}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("TESTY")
    print("=" * 72)

    record("P1: K*_A i K*_B znalezione dla >= 5 alpha",
           n_AB >= 5, f"{n_AB} pkt d_AB, {n_An} pkt d_An")

    dAB_near = min((d for a, d in pts_AB if a > ALPHA_MAX_REF-0.02), default=np.nan)
    dAB_ns = f'{dAB_near:.7f}' if np.isfinite(dAB_near) else 'brak'
    record("P2: delta_AB < 0.001 blisko alpha_max",
           np.isfinite(dAB_near) and dAB_near < 0.001, f"min_dAB={dAB_ns}")

    dAn_near = min((d for a, d in pts_An if a > ALPHA_MAX_REF-0.02), default=np.nan)
    dAn_ns = f'{dAn_near:.7f}' if np.isfinite(dAn_near) else 'brak'
    record("P3: delta_An < 0.0003 blisko alpha_max",
           np.isfinite(dAn_near) and dAn_near < 0.0003, f"min_dAn={dAn_ns}")

    SN_s   = f'{alpha_SN:.5f}'  if np.isfinite(alpha_SN)  else 'brak'
    nan_s  = f'{alpha_nan:.5f}' if np.isfinite(alpha_nan) else 'brak'
    diff_SN_max = f'{abs(alpha_SN-ALPHA_MAX_REF):.5f}' if np.isfinite(alpha_SN) else 'brak'
    diff_nan_max = f'{abs(alpha_nan-ALPHA_MAX_REF):.5f}' if np.isfinite(alpha_nan) else 'brak'
    diff_SN_nan = f'{abs(alpha_SN-alpha_nan):.5f}' \
                  if np.isfinite(alpha_SN) and np.isfinite(alpha_nan) else 'brak'

    record("P4: |alpha_SN - alpha_max_ref| < 0.01 (solitony zbiegaja sie blisko max)",
           np.isfinite(alpha_SN) and abs(alpha_SN-ALPHA_MAX_REF) < 0.01,
           f"alpha_SN={SN_s}, odch={diff_SN_max}")

    record("P5: |alpha_SN - alpha_nan| < 0.005 (oba warunkow przy tym samym alpha)",
           np.isfinite(alpha_SN) and np.isfinite(alpha_nan) and
           abs(alpha_SN-alpha_nan) < 0.005,
           f"SN={SN_s}, nan={nan_s}, roznica={diff_SN_nan}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p91_codim2_branches.py")
    print("=" * 72)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    print(f"""
  alpha_SN    = {SN_s}    (K*_A = K*_B)
  alpha_nan   = {nan_s}    (K*_A = K_nan_R)
  alpha_max_ref = {ALPHA_MAX_REF}
  |SN - max|  = {diff_SN_max}
  |nan - max| = {diff_nan_max}
  |SN - nan|  = {diff_SN_nan}
""")

    if np.isfinite(alpha_SN) and np.isfinite(alpha_nan):
        if abs(alpha_SN-ALPHA_MAX_REF) < 0.005 and abs(alpha_nan-ALPHA_MAX_REF) < 0.005:
            print("  WNIOSEK: KODIMENSION-2 POTWIERDZONE.")
            print("  Oba warunki (K*_A=K*_B oraz K*_A=K_nan_R) zachodza przy alpha_max.")
        elif abs(alpha_SN-alpha_nan) < 0.005:
            print(f"  WNIOSEK: QUASI-KODIMENSION-2.")
            print(f"  alpha_SN ~ alpha_nan = {SN_s}, ale odchylenie od alpha_max={diff_SN_max}")
            print("  Moze poprawic sie przy wiekszej precyzji lub innym alpha_max.")
        else:
            print(f"  WNIOSEK: Dwa progi NIE sa rowne (|SN-nan|={diff_SN_nan}).")
            print("  Brak kodimension-2 lub bledne extrapol.")
    elif n_AB < 3:
        print("  WNIOSEK: Za malo danych K*_A/K*_B. Sprawdz K_arr lub alpha range.")
    else:
        print("  WNIOSEK: Niewystarczajace dane do oceny.")


if __name__ == '__main__':
    run_main()
