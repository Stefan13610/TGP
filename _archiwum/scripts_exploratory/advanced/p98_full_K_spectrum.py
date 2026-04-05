#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p98_full_K_spectrum.py  --  TGP v1
====================================
Cel: Zbadac, czy wyzsze K*_n (n=2,3,...) odpowiadajace generacjom mu, tau
     sa FIZYCZNE (odporne na r_max) czy sa artefaktami r_max=40
     (analogicznie do K_nan_R i Galezi L).

Kontekst:
  p77: K*_2 = 0.034003 (K*_2/K*_1 = 3.27) przy alpha_K, r_max=40
  p86: "4 solitony" przy alpha=8.5-8.7, r_max=40
  p95: K_nan_R (artefakt fazy oscylacyjnego ogona przy r=40)
  => Hipoteza: K*_2, K*_3 to takze artefakty r_max=40

Metoda:
  Skan K in [0.006, 0.250] (50 punktow) przy alpha_K
  Dla kazdego K: znajdz WSZYSTKIE zera g(K) przez skan psi i g_U(K)
  Porownaj: r_max = 40, 60, 80
  => Czy K*_2, K*_3 przetrwaja przy r_max=60,80?

  Uwaga: uzywa PELNEGO skanu psi in [1.001, 5.0] (szerokie),
         aby nie przeoczyc wyzszych korzeni F(psi)=phi(r_max)-1

Testy:
  P1: K*_1 z p97 odtworzone (K*(inf)~0.010029, sanity check)
  P2: Liczba zer g(K) w [0.006,0.25] MALEJE z r_max
      (artefakty znikaja przy wiekszym r_max)
  P3: Esli K*_2 przezyje r_max=80: K*_2/K*_1 jest stabilne (< 1% zmiana)
  P4: K*_2(inf)/K*_1(inf) -- prawdziwy stosunek mas mu/e
  P5: Esli K*_2 znika: 3-generacyjna hipoteza wymaga modyfikacji V_mod

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.linalg import lstsq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K  = 8.5616
A_GAM    = 0.040
LAM      = 5.501357e-06
GAMMA_V  = 1.0
V1       = GAMMA_V/3 - GAMMA_V/4

# Parametry skanu
K_LO   = 0.006
K_HI   = 0.250
N_K    = 50
N_PSI  = 100
PSI_LO = 1.001
PSI_HI = 5.0   # Szerszy zakres psi!
N_EVAL = 500
RMAXES = [40.0, 60.0, 80.0]

# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K/phi
    return [dphi, dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, r_max):
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y), [A_GAM, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        return float(sol.y[0,-1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except:
        return np.nan

def scan_all_psi_zeros(K, r_max, psi_lo=PSI_LO, psi_hi=PSI_HI, n_psi=N_PSI):
    """Zwraca WSZYSTKIE zera F(psi)=phi(r_max)-1 w [psi_lo, psi_hi]."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv = [phi_at_rmax(p, K, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,r_max)-1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=40)
                roots.append(pz)
            except: pass
    return roots  # Lista wszystkich zer (rosnaco w psi)

def compute_g_at_psi(K, psi0, r_max):
    """g = E/(4piK) - 1."""
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0,1,N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y), [A_GAM,r_max],
                        [psi0,dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        if sol.t[-1] < r_max*0.99: return np.nan
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA_K/phi)*r**2, r)
        Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return (Ek+Ep)/(4*np.pi*K)-1.0
    except:
        return np.nan

def find_g_of_K(K, r_max, branch='U'):
    """
    Zwraca g dla danego K i r_max.
    branch='U': ostatnie zero psi (Galaz U, fizyczna)
    branch='all': lista (psi, g) dla WSZYSTKICH zer
    """
    roots = scan_all_psi_zeros(K, r_max)
    if not roots:
        return np.nan, np.nan
    if branch == 'U':
        psi_use = roots[-1]
        return compute_g_at_psi(K, psi_use, r_max), psi_use
    elif branch == 'all':
        return [(compute_g_at_psi(K, psi, r_max), psi) for psi in roots]
    return np.nan, np.nan

def full_K_scan(r_max):
    """
    Skan K in [K_LO, K_HI], zwraca:
    - Ks: lista wartosci K
    - g_U: g dla Galezy U (ostatnie zero psi)
    - n_zeros: liczba zer psi dla kazdego K
    - all_zeros_g: lista [(psi, g)] dla wszystkich zer przy kazdym K
    """
    Ks = np.linspace(K_LO, K_HI, N_K)
    g_U_list = []
    n_zeros_list = []
    all_glist = []

    for K in Ks:
        roots = scan_all_psi_zeros(K, r_max)
        n_zeros_list.append(len(roots))
        if roots:
            g_U, _ = find_g_of_K(K, r_max, 'U')
            g_U_list.append(g_U)
            # g dla wszystkich zer
            all_g = [(compute_g_at_psi(K, psi, r_max), psi) for psi in roots]
            all_glist.append(all_g)
        else:
            g_U_list.append(np.nan)
            all_glist.append([])

    return Ks, np.array(g_U_list), np.array(n_zeros_list), all_glist

def find_zeros_of_g(Ks, g_arr, atol=1e-6):
    """Znajdz miejsca zerowe g(K) przez zmiane znaku."""
    zeros = []
    for i in range(len(g_arr)-1):
        if np.isfinite(g_arr[i]) and np.isfinite(g_arr[i+1]) and g_arr[i]*g_arr[i+1] < 0:
            try:
                # Interpolacja liniowa
                Kz = Ks[i] - g_arr[i]*(Ks[i+1]-Ks[i])/(g_arr[i+1]-g_arr[i])
                zeros.append(Kz)
            except: pass
    return zeros

# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 68)
    print("TGP v1  *  p98_full_K_spectrum.py")
    print(f"  alpha_K={ALPHA_K}, K in [{K_LO},{K_HI}], n_K={N_K}")
    print(f"  psi in [{PSI_LO},{PSI_HI}], n_psi={N_PSI}")
    print(f"  r_max = {RMAXES}")
    print("=" * 68)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ─────────────────────────────────────────────────────────────────────────
    results = {}  # r_max -> (Ks, g_U, n_zeros, all_glist)

    for r_max in RMAXES:
        print(f"\n== Skan r_max={r_max:.0f} (K=[{K_LO},{K_HI}], {N_K} pkt) ==")
        t0 = time.time()
        Ks, g_U, n_zeros, all_glist = full_K_scan(r_max)
        dt = time.time()-t0
        print(f"  Czas: {dt:.0f}s")

        # Znajdz zera g_U
        zeros_U = find_zeros_of_g(Ks, g_U)
        print(f"  Zera g_U (Galaz U): n={len(zeros_U)}, K* = {[f'{z:.5f}' for z in zeros_U]}")
        print(f"  n_zeros psi: min={n_zeros.min()}, max={n_zeros.max()}, "
              f"srednia={n_zeros.mean():.1f}")

        # Podsumowanie g_U w kilku punktach
        print(f"  {'K':>8}  {'g_U':>12}  {'n_psi_zeros':>12}")
        step = max(1, N_K//10)
        for i in range(0, N_K, step):
            gu_s = f'{g_U[i]:.5f}' if np.isfinite(g_U[i]) else '   nan  '
            print(f"  {Ks[i]:.5f}  {gu_s}  {n_zeros[i]:>12}")

        results[r_max] = (Ks, g_U, n_zeros, all_glist, zeros_U)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== POROWNANIE LICZBY ZER g_U(K) wedlug r_max ==")
    for r_max in RMAXES:
        z = results[r_max][4]
        print(f"  r_max={r_max:.0f}: {len(z)} zer g_U: {[f'{k:.5f}' for k in z]}")

    # P1: K*_1 z p97 odtworzone
    K_ref = 0.010029
    K_1_40 = results.get(40.0, (None,)*5)[4]
    if K_1_40:
        closest = min(K_1_40, key=lambda k: abs(k-K_ref))
        p1_ok = abs(closest - K_ref) < 0.002
        record("P1: K*_1(r=40) ~ 0.010029 (z p97)", p1_ok,
               f"znaleziono K*={closest:.6f}, ref={K_ref:.6f}")
    else:
        record("P1: K*_1 znaleziona", False, "brak zer")

    # P2: liczba zer maleje z r_max?
    n_zeros_by_rmax = {r: len(results[r][4]) for r in RMAXES}
    p2_ok = n_zeros_by_rmax.get(40.0,0) >= n_zeros_by_rmax.get(80.0,0)
    record("P2: liczba zer g_U maleje (lub rowna) r_max: 40->80",
           p2_ok, f"n(r=40)={n_zeros_by_rmax.get(40.0,'?')}, "
           f"n(r=60)={n_zeros_by_rmax.get(60.0,'?')}, "
           f"n(r=80)={n_zeros_by_rmax.get(80.0,'?')}")

    # P3/P4: Czy K*_2 przezywa przy r_max=80?
    zeros_40 = results.get(40.0,(None,)*5)[4] or []
    zeros_60 = results.get(60.0,(None,)*5)[4] or []
    zeros_80 = results.get(80.0,(None,)*5)[4] or []

    K1_40 = zeros_40[0] if zeros_40 else np.nan
    K1_60 = zeros_60[0] if zeros_60 else np.nan
    K1_80 = zeros_80[0] if zeros_80 else np.nan

    if len(zeros_40) >= 2:
        K2_40 = sorted(zeros_40)[1]
        K2_60 = sorted(zeros_60)[1] if len(zeros_60) >= 2 else np.nan
        K2_80 = sorted(zeros_80)[1] if len(zeros_80) >= 2 else np.nan
        ratio_40 = K2_40/K1_40 if np.isfinite(K1_40) and K1_40>0 else np.nan

        if np.isfinite(K2_80):
            ratio_80 = K2_80/K1_80 if np.isfinite(K1_80) and K1_80>0 else np.nan
            if np.isfinite(ratio_40) and np.isfinite(ratio_80):
                p3_ok = abs(ratio_80/ratio_40 - 1) < 0.01
                record("P3: K*_2/K*_1 stabilne przy r_max: 40->80 (< 1%)",
                       p3_ok, f"ratio(40)={ratio_40:.4f}, ratio(80)={ratio_80:.4f}")
            print(f"\n  K*_2 OCALALO przy r_max=80!")
            print(f"  K*_2(r=40)={K2_40:.5f}, K*_2(r=60)={K2_60:.5f}, K*_2(r=80)={K2_80:.5f}")
            print(f"  K*_2/K*_1(r=40)={ratio_40:.4f}")
            # Ekstrapolacja K*_2(inf)
            pts = [(r, K) for r, K in [(40.0,K2_40),(60.0,K2_60),(80.0,K2_80)]
                   if np.isfinite(K)]
            if len(pts) >= 2:
                rs = np.array([x[0] for x in pts])
                Ks_arr = np.array([x[1] for x in pts])
                X = np.column_stack([np.ones_like(rs), 1.0/rs])
                coef, _, _, _ = lstsq(X, Ks_arr, rcond=None)
                K2_inf = coef[0]
                K1_inf = 0.010029  # z p97
                ratio_inf = K2_inf/K1_inf
                print(f"  K*_2(inf) ~ {K2_inf:.6f}")
                print(f"  K*_2/K*_1(inf) ~ {ratio_inf:.4f}")
                print(f"  Porownanie z r_21=206.77: {ratio_inf:.4f} vs 206.77")
                record("P4: K*_2(inf)/K*_1(inf) wyznaczone",
                       True, f"ratio_inf={ratio_inf:.4f}")
        else:
            print(f"\n  K*_2 z r_max=40 ({K2_40:.5f}) NIE OCALALO przy r_max=80!")
            print(f"  To potwierdza: K*_2 = ARTEFAKT r_max=40")
            record("P3: K*_2/K*_1 stabilne", False, "K*_2 znika przy r_max=80")
            record("P5: K*_2 = artefakt (hipoteza 3-generacyjna wymaga modyfikacji V_mod)",
                   True, f"K*_2 nie przezywa r_max=80 — brak 2. fizycznego solitonu")
    else:
        print(f"\n  Brak K*_2 juz przy r_max=40 (tylko 1 zero g_U)")
        record("P3: brak K*_2", False, "tylko 1 zero g_U")
        record("P4: brak K*_2", False, "")
        record("P5: hipoteza 3-generacyjna", True, "wymaga modyfikacji V_mod")

    # ─────────────────────────────────────────────────────────────────────────
    # Wykres g_U(K) dla trzech r_max
    print("\n== WYKRES g_U(K) dla r_max=40,60,80 ==")
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        colors = {40.0:'b', 60.0:'g', 80.0:'r'}

        for r_max in RMAXES:
            Ks, g_U, n_zeros, _, zeros_U = results[r_max]
            c = colors[r_max]
            # Pelny skan
            ax = axes[0]
            valid = np.isfinite(g_U)
            ax.plot(Ks[valid], g_U[valid], 'o-', color=c, lw=1.5, ms=4,
                    label=f'r_max={int(r_max)}, {len(zeros_U)} zer')
            for kz in zeros_U:
                ax.axvline(kz, color=c, alpha=0.5, ls='--', lw=1)

        ax.axhline(0, color='k', ls='-', lw=0.8)
        ax.set_xlabel('K'); ax.set_ylabel('g_U(K)')
        ax.set_title(f'g_U(K) (Galaz U), alpha_K={ALPHA_K}')
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
        ax.set_xlim(K_LO, K_HI); ax.set_ylim(-3, 3)

        # Zoom na K*_1
        ax = axes[1]
        K_zoom = 0.025
        for r_max in RMAXES:
            Ks, g_U, n_zeros, _, zeros_U = results[r_max]
            c = colors[r_max]
            mask = Ks <= K_zoom
            valid = mask & np.isfinite(g_U)
            ax.plot(Ks[valid], g_U[valid], 'o-', color=c, lw=1.5, ms=4,
                    label=f'r_max={int(r_max)}')
            for kz in zeros_U:
                if kz <= K_zoom:
                    ax.axvline(kz, color=c, alpha=0.5, ls='--', lw=1)
        ax.axhline(0, color='k', ls='-', lw=0.8)
        ax.axvline(0.010029, color='m', ls=':', lw=2, label=f'K*(inf)=0.010029')
        ax.set_xlabel('K'); ax.set_ylabel('g_U(K)')
        ax.set_title(f'Zoom: K in [{K_LO},{K_zoom}]')
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

        plt.suptitle(f'p98: Widmo K* (Galaz U), alpha_K={ALPHA_K}', fontsize=12)
        plt.tight_layout()
        out_path = 'p98_full_K_spectrum.png'
        plt.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Wykres zapisany: {out_path}")
    except Exception as e:
        print(f"  BLAD wykresu: {e}")

    # ─────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 68)
    print("PODSUMOWANIE  p98_full_K_spectrum.py")
    print("=" * 68)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    print("\n== INTERPRETACJA ==")
    n80 = n_zeros_by_rmax.get(80.0, 0)
    n40 = n_zeros_by_rmax.get(40.0, 0)
    if n80 < n40:
        print(f"  -> Liczba zer maleje z r_max: {n40}(r=40) -> {n80}(r=80)")
        print("     Wyzsze K*_n sa artefaktami oscylacyjnego ogona!")
        print("     WNIOSEK: Standardowy V_mod ma tylko 1 fizyczny soliton na alpha_K.")
        print("     Hipoteza 3-generacyjna (3 distinct K*) WYMAGA MODYFIKACJI V_mod.")
    elif n80 == n40:
        print(f"  -> Liczba zer stabilna: {n40} przy obu r_max")
        print("     Wszystkie K*_n sa potencjalnie fizyczne!")
    else:
        print(f"  -> Nieoczekiwany wynik: {n40}(r=40) -> {n80}(r=80)")

    print()


if __name__ == '__main__':
    run_main()
