#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p95_knan_vs_rmax.py  --  TGP v1
=================================
Cel: Zrozumiec zachowanie K_nan_R(r_max) i wyjasnic wynik p94
     (K_nan_R 'brak' przy r_max=60,100 w zakresie K<0.016).

Pytania:
  Q1: Czy K_nan_R(r_max=60,100) istnieje przy wyzszych K (do 0.040)?
  Q2: Jak K_nan_R zmienia sie z r_max? (Czy to artefakt r_max=40?)
  Q3: Czy phi(r) rzeczywiscie dywerguje (phi->inf) przed r_max,
      czy tez "dywergencja" to artefakt oscylacyjnego ogona?

Metoda:
  Czesc A: Skan K w [0.008, 0.040] dla alpha=8.8734, r_max=40,60,100
           -> znajdz K_nan_R dla kazdego r_max
  Czesc B: Profil F(psi) = phi(r_max) - 1  vs K blisko K_nan_R(r_max=40)
           -> czy brak zerow F(psi) jest z powodu dywergencji czy braku
  Czesc C: Profil phi(r) przy K=K_nan_R dla roznych r_max
           -> obserwacja behawioralna: czy phi->inf w srodku czy na krawedzi
  Czesc D: Powrot do pytania alpha_max: gap(alpha, r_max=60,100) w nowym zakresie

Testy:
  P1: K_nan_R(r_max=40) w p94 odtworzone (sanity check)
  P2: K_nan_R(r_max=60) istnieje w [0.016, 0.040]
  P3: K_nan_R(r_max) monotoniczne rosnace z r_max
  P4: Przy K > K_nan_R(r_max=40), phi(r) dywerguje PRZED r_max=40
       (fizyczna dywergencja, nie artefakt)
  P5: alpha_max(r_max=60) wyznaczone i blizkie alpha_max(r_max=40)

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K    = 8.5616
ALPHA_REF  = 8.8734   # alpha_max z p80 przy r_max=40
A_GAM      = 0.040
LAM        = 5.501357e-06
GAMMA_V    = 1.0
V1         = GAMMA_V/3 - GAMMA_V/4

N_PSI  = 80
N_EVAL = 500

# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, alpha, r_max, n_eval=None):
    """Zwraca phi(r_max). np.nan jesli dywergencja przed r_max."""
    ne = n_eval or N_EVAL
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, ne)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] >= r_max*0.99:
            return float(sol.y[0,-1])
        return np.nan
    except:
        return np.nan

def phi_trajectory(psi, K, alpha, r_max, n_eval=800):
    """Zwraca (r_arr, phi_arr) dla pelnego profilu."""
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        return sol.t, sol.y[0]
    except:
        return np.array([A_GAM]), np.array([psi])

def scan_psi_zeros(K, alpha, r_max):
    """Zwraca liste zer F(psi) = phi(r_max)-1."""
    psis = np.linspace(1.001, 3.0, N_PSI)
    Fv = [phi_at_rmax(p, K, alpha, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,alpha,r_max)-1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=40)
                roots.append(pz)
            except: pass
    n_nan = sum(1 for f in Fv if not np.isfinite(f))
    return roots, n_nan, psis, Fv

def find_K_nan_R(alpha, r_max, K_lo=0.008, K_hi=0.040, n_K=40):
    """
    Szuka pierwszej K gdzie scan_psi_zeros zwraca 0 zerow.
    Zwraca (K_nan_R, last_K_with_roots).
    """
    Ks = np.linspace(K_lo, K_hi, n_K)
    last_K_ok = np.nan
    K_nan_R = np.nan
    for K in Ks:
        roots, n_nan, _, _ = scan_psi_zeros(K, alpha, r_max)
        if roots:
            last_K_ok = K
        elif np.isfinite(last_K_ok):
            # Bisection
            Kl, Kh = last_K_ok, K
            for _ in range(15):
                Km = 0.5*(Kl+Kh)
                rt, _, _, _ = scan_psi_zeros(Km, alpha, r_max)
                if rt: Kl = Km
                else:  Kh = Km
                if Kh-Kl < 1e-5: break
            K_nan_R = 0.5*(Kl+Kh)
            break
    return K_nan_R, last_K_ok

def g_func(K, alpha, r_max):
    """Zwraca g = E/(4pi K) - 1 dla Galezy U (ostatnie zero psi)."""
    roots, _, _, _ = scan_psi_zeros(K, alpha, r_max)
    if not roots: return np.nan, np.nan
    psi_z = roots[-1]
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0,1,N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM,r_max],
                        [psi_z,dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2,r)
        Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2,r)
        return (Ek+Ep)/(4*np.pi*K)-1.0, psi_z
    except: return np.nan, psi_z

def find_K_star(alpha, r_max, K_lo=0.008, K_hi=0.020, n_K=25):
    """Szuka K* (zero g_func) przez brentq."""
    Ks = np.linspace(K_lo, K_hi, n_K)
    gv = []
    for K in Ks:
        g, _ = g_func(K, alpha, r_max)
        gv.append(g)
    K_star = np.nan; psi_star = np.nan
    for i in range(len(gv)-1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1] < 0:
            try:
                Kz = brentq(lambda K: g_func(K,alpha,r_max)[0],
                            Ks[i], Ks[i+1], xtol=5e-6, maxiter=20)
                K_star = Kz
                _, psi_star = g_func(Kz, alpha, r_max)
            except: pass
    return K_star, psi_star

# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 68)
    print("TGP v1  *  p95_knan_vs_rmax.py")
    print(f"  alpha_ref={ALPHA_REF}, alpha_K={ALPHA_K}")
    print(f"  Badanie K_nan_R(r_max) w rozszerzonym zakresie K=[0.008, 0.040]")
    print("=" * 68)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    RMAXES_A = [40.0, 60.0, 80.0, 100.0]
    results_A = {}   # r_max -> K_nan_R

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC A: K_nan_R(r_max) przy alpha=8.8734, K w [0.008, 0.040] ==")
    print(f"  {'r_max':>6}  {'K_nan_R':>12}  {'last_K_ok':>12}  t")
    print("  " + "-"*48)

    for r_max in RMAXES_A:
        t0 = time.time()
        K_nan, last_ok = find_K_nan_R(ALPHA_REF, r_max,
                                       K_lo=0.008, K_hi=0.040, n_K=40)
        dt = time.time()-t0
        Kn_s  = f'{K_nan:.6f}'  if np.isfinite(K_nan)  else '  brak  '
        Ko_s  = f'{last_ok:.6f}' if np.isfinite(last_ok) else '  brak  '
        print(f"  {r_max:>6.0f}  {Kn_s}  {Ko_s}  {dt:.0f}s")
        results_A[r_max] = K_nan

    # P1: Odtworz K_nan_R dla r_max=40 (powinno byc ~0.0108 z p94)
    knan_40 = results_A.get(40.0, np.nan)
    p1_ok = np.isfinite(knan_40) and 0.010 < knan_40 < 0.015
    record("P1: K_nan_R(r_max=40) ~ 0.010-0.015 (sanity check z p94)",
           p1_ok, f"K_nan_R={knan_40:.6f}" if np.isfinite(knan_40) else "brak")

    # P2: K_nan_R dla r_max=60 istnieje (szerzej niz 0.016?)
    knan_60 = results_A.get(60.0, np.nan)
    p2_ok = np.isfinite(knan_60)
    record("P2: K_nan_R(r_max=60) znaleziony w [0.008, 0.040]",
           p2_ok, f"K_nan_R={knan_60:.6f}" if np.isfinite(knan_60) else "brak > 0.040 lub nie istnieje")

    # P3: Monotonicznosc K_nan_R(r_max)
    knan_vals = [results_A[r] for r in RMAXES_A if np.isfinite(results_A.get(r, np.nan))]
    if len(knan_vals) >= 2:
        p3_ok = all(knan_vals[i] <= knan_vals[i+1] for i in range(len(knan_vals)-1))
        record("P3: K_nan_R(r_max) monotoniczne rosnace",
               p3_ok, f"vals={[f'{v:.5f}' for v in knan_vals]}")
    else:
        record("P3: monotonicznosc", False, "za malo punktow")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC B: F(psi) przy K bliskim K_nan_R(r_max=40) ==")
    knan_40 = results_A.get(40.0, np.nan)
    if np.isfinite(knan_40):
        K_test_vals = [knan_40 - 0.001, knan_40, knan_40 + 0.001]
        for r_max in [40.0, 60.0]:
            print(f"\n  r_max={r_max:.0f}:")
            for K_t in K_test_vals:
                roots, n_nan, psis, Fv = scan_psi_zeros(K_t, ALPHA_REF, r_max)
                fin_vals = [f for f in Fv if np.isfinite(f)]
                n_roots = len(roots)
                f_range = f"[{min(fin_vals):.4f},{max(fin_vals):.4f}]" if fin_vals else "brak"
                print(f"    K={K_t:.5f}: n_zeros={n_roots}, n_nan={n_nan}/{N_PSI}, "
                      f"F_range={f_range}")

    # P4: Przy K > K_nan_R(r_max=40), phi dywerguje przed r=40
    print("\n== CZESC C: Profil phi(r) przy K=K_nan_R+0.002, r_max=100 ==")
    if np.isfinite(knan_40):
        K_above = knan_40 + 0.002
        psi_test = 1.5
        r_arr, phi_arr = phi_trajectory(psi_test, K_above, ALPHA_REF, 100.0, n_eval=1000)
        r_max_reached = r_arr[-1]
        phi_max = np.max(np.abs(phi_arr)) if len(phi_arr) > 0 else np.nan
        phi_at_40 = np.interp(40.0, r_arr, phi_arr) if r_arr[-1] >= 40 else np.nan
        print(f"  K={K_above:.5f}, psi={psi_test}")
        print(f"  r_max_reached={r_max_reached:.2f}, phi_max={phi_max:.4f}")
        print(f"  phi(r=40)={phi_at_40:.4f}" if np.isfinite(phi_at_40) else "  phi(r=40)=DYWERGUJE przed r=40")

        # Czy phi dywerguje przed r_max=40?
        phi_at_39 = np.interp(39.0, r_arr, phi_arr) if r_arr[-1] >= 39 else np.nan
        p4_ok = (not np.isfinite(phi_at_40)) or (np.isfinite(phi_at_40) and abs(phi_at_40) > 100)
        record("P4: phi dywerguje (|phi|>100 lub nan) przed r=40 przy K>K_nan_R",
               p4_ok, f"phi(r=40)={phi_at_40:.3f}" if np.isfinite(phi_at_40) else "dywerguje przed r=40")
    else:
        record("P4: brak K_nan_R", False, "nie mozna sprawdzic")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC D: gap(alpha, r_max) przy rozszerzonym zakresie K ==")
    ALPHAS_D = [8.840, 8.860, 8.8734, 8.880, 8.890, 8.900, 8.920]
    RMAXES_D = [40.0, 60.0, 100.0]

    print(f"  {'alpha':>7}  {'r_max':>6}  {'K*':>12}  {'K_nan_R':>12}  {'gap':>12}")
    print("  " + "-"*58)

    results_D = {}
    for alpha in ALPHAS_D:
        for r_max in RMAXES_D:
            t0 = time.time()
            K_nan, _ = find_K_nan_R(alpha, r_max, K_lo=0.008, K_hi=0.040, n_K=35)
            K_star, _ = find_K_star(alpha, r_max, K_lo=0.008, K_hi=0.020, n_K=20)
            gap = K_nan - K_star if (np.isfinite(K_nan) and np.isfinite(K_star)) else np.nan
            dt = time.time()-t0
            ks_s   = f'{K_star:.7f}' if np.isfinite(K_star)  else '   brak   '
            kn_s   = f'{K_nan:.7f}'  if np.isfinite(K_nan)   else '   brak   '
            gap_s  = f'{gap:+.7f}'   if np.isfinite(gap)     else '   brak   '
            sign = ' [+]' if (np.isfinite(gap) and gap > 0) else (' [0]' if (np.isfinite(gap) and gap <= 0) else '')
            print(f"  {alpha:.4f}  {r_max:>6.0f}  {ks_s}  {kn_s}  {gap_s}  ({dt:.0f}s){sign}")
            results_D[(alpha, r_max)] = (K_star, K_nan, gap)

    # alpha_max(r_max) z zerow gap
    print("\n== alpha_max(r_max) z CZESCI D ==")
    alpha_max_D = {}
    for r_max in RMAXES_D:
        gaps = [(a, results_D[(a,r_max)][2]) for a in ALPHAS_D
                if np.isfinite(results_D[(a,r_max)][2])]
        # Znajdz przejscie gap: + -> 0 -> -
        last_pos = np.nan
        for a, g in gaps:
            if g > 0: last_pos = a
        # Ekstrapolacja liniowa ostatnich 2-3 punktow
        pos_pts = [(a, g) for a, g in gaps if g > 0]
        neg_pts = [(a, g) for a, g in gaps if g < 0]
        a_extrap = np.nan
        if len(pos_pts) >= 1 and len(neg_pts) >= 1:
            # Interpolacja miedzy ostatnim pos i pierwszym neg
            ap, gp = pos_pts[-1]
            an, gn = neg_pts[0]
            # Liniowa interpolacja do g=0
            a_extrap = ap - gp*(an-ap)/(gn-gp)
        alpha_max_D[r_max] = (last_pos, a_extrap)
        la_s = f'{last_pos:.4f}' if np.isfinite(last_pos) else 'brak'
        ex_s = f'{a_extrap:.4f}' if np.isfinite(a_extrap) else 'brak'
        ratio_l = ALPHA_K/last_pos if np.isfinite(last_pos) else np.nan
        ratio_e = ALPHA_K/a_extrap if np.isfinite(a_extrap) else np.nan
        rl_s = f'{ratio_l:.5f}' if np.isfinite(ratio_l) else 'brak'
        re_s = f'{ratio_e:.5f}' if np.isfinite(ratio_e) else 'brak'
        print(f"  r_max={r_max:.0f}: alpha_max~{la_s}..{ex_s}, "
              f"ratio_last={rl_s}, ratio_extrap={re_s}")

    # P5: alpha_max porownanie r_max=40 vs r_max=60
    am_40 = alpha_max_D.get(40.0, (np.nan,)*2)[1]
    am_60 = alpha_max_D.get(60.0, (np.nan,)*2)[1]
    if np.isfinite(am_40) and np.isfinite(am_60):
        delta = abs(am_40 - am_60)
        p5_ok = delta < 0.05
        record("P5: |alpha_max(r=60) - alpha_max(r=40)| < 0.05",
               p5_ok, f"delta={delta:.4f}, am40={am_40:.4f}, am60={am_60:.4f}")
    else:
        record("P5: alpha_max z r_max=60 wyznaczone", np.isfinite(am_60),
               f"am60={am_60:.4f}" if np.isfinite(am_60) else "brak")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== WYKRES: K_nan_R i K* vs alpha dla r_max=40,60,100 ==")
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Lewy: K_nan_R vs alpha
        colors = {40.0: 'b', 60.0: 'g', 100.0: 'r'}
        for r_max in RMAXES_D:
            kn_vals = [(a, results_D[(a,r_max)][1]) for a in ALPHAS_D
                       if np.isfinite(results_D[(a,r_max)][1])]
            ks_vals = [(a, results_D[(a,r_max)][0]) for a in ALPHAS_D
                       if np.isfinite(results_D[(a,r_max)][0])]
            c = colors[r_max]
            if kn_vals:
                ax1.plot([x[0] for x in kn_vals], [x[1] for x in kn_vals],
                         'o--', color=c, label=f'K_nan_R r={r_max:.0f}')
            if ks_vals:
                ax1.plot([x[0] for x in ks_vals], [x[1] for x in ks_vals],
                         's-', color=c, alpha=0.5, label=f'K* r={r_max:.0f}')
        ax1.axvline(ALPHA_REF, color='k', ls=':', alpha=0.5, label=f'alpha_ref={ALPHA_REF}')
        ax1.set_xlabel('alpha'); ax1.set_ylabel('K')
        ax1.set_title('K_nan_R(alpha) i K*(alpha) dla roznych r_max')
        ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

        # Prawy: gap vs alpha
        for r_max in RMAXES_D:
            gap_vals = [(a, results_D[(a,r_max)][2]) for a in ALPHAS_D
                        if np.isfinite(results_D[(a,r_max)][2])]
            c = colors[r_max]
            if gap_vals:
                ax2.plot([x[0] for x in gap_vals], [x[1] for x in gap_vals],
                         'o-', color=c, label=f'r_max={r_max:.0f}')
        ax2.axhline(0, color='k', ls='-', lw=0.8)
        ax2.axvline(ALPHA_REF, color='k', ls=':', alpha=0.5)
        ax2.set_xlabel('alpha'); ax2.set_ylabel('gap = K_nan_R - K*')
        ax2.set_title('gap(alpha, r_max): alpha_max gdzie gap=0')
        ax2.legend(); ax2.grid(True, alpha=0.3)

        plt.suptitle(f'p95: K_nan_R vs r_max  (alpha_K={ALPHA_K})', fontsize=12)
        plt.tight_layout()
        out_path = 'p95_knan_vs_rmax.png'
        plt.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Wykres zapisany: {out_path}")
    except Exception as e:
        print(f"  BLAD wykresu: {e}")

    # ═══════════════════════════════════════════════════════════════════════
    print()
    print("=" * 68)
    print("PODSUMOWANIE  p95_knan_vs_rmax.py")
    print("=" * 68)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    # Interpretacja
    print("\n== INTERPRETACJA ==")
    knan_40 = results_A.get(40.0, np.nan)
    knan_60 = results_A.get(60.0, np.nan)
    knan_80 = results_A.get(80.0, np.nan)
    knan_100= results_A.get(100.0,np.nan)
    if all(np.isfinite(v) for v in [knan_40, knan_60]):
        if knan_60 > knan_40 * 1.1:
            print("  -> K_nan_R(r_max) ROSNIE z r_max")
            print("     Oznacza to: graniczny K_nan_R jest artefaktem r_max")
            print("     Definicja alpha_max przez K*=K_nan_R jest r_max-zalezna!")
        else:
            print("  -> K_nan_R(r_max) SLABO zalezne od r_max (dobra wiadomosc)")
            print("     Definicja alpha_max przez K*=K_nan_R jest r_max-odporna")
    elif not np.isfinite(knan_60):
        print("  -> K_nan_R(r_max=60) NIE ISTNIEJE w [0.008, 0.040]")
        print("     Brak fizycznej granicy K_nan_R przy wiekszych r_max!")
        print("     Definicja alpha_max=8.8734 jest ARTEFAKTEM r_max=40.")
        print("     Koincydencja alpha_K/alpha_max = n_s jest artefaktem konwencji.")

    print()


if __name__ == '__main__':
    run_main()
