# -*- coding: utf-8 -*-
"""
P59: ZBIEZNOSC SIATKI -- BADANIE OP-13

Cel: Dla punktu B (a=0.040049, alpha=8.5612, lambda=5.4677e-6)
     obliczyc K1, K2 przy N = 500, 1000, 2000, 3000, 5000, 8000, 10000
     i sprawdzic:
       - zbieznosc K1(N), K2(N) do N->inf (Richardson extrapolation)
       - N* takie, ze blad < 0.01%
       - tabela bledow wzglednych vs N
       - estymata N->inf (Richardson)
       - jaki wplyw na r21, r31, Q

Dodatkowe: to samo dla punktu A (historycznego).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq

# -------------------------------------------------------------------
# PARAMETRY
# -------------------------------------------------------------------
R_MAX = 50.0

def E1(x):
    return -expi(-x)

def V_mod(phi, lam):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha, a, lam, N):
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r)*(-r-1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1.0+alpha/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(phi,lam)-V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a, lam, N):
    return energy_num(K, alpha, a, lam, N) / (4*np.pi*K) - 1.0

def find_K1(alpha, a, lam, N, tol=1e-12):
    try:
        return brentq(lambda K: g_func(K, alpha, a, lam, N),
                      1e-5, 0.45, xtol=tol)
    except Exception:
        return np.nan

def find_K2(alpha, a, lam, N, tol=1e-12):
    try:
        return brentq(lambda K: g_func(K, alpha, a, lam, N),
                      0.4, 5.0, xtol=tol)
    except Exception:
        return np.nan

def K3_ei(a, lam):
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=300)
    return np.sqrt(3*I4/(2*lam*I6)) if I6 > 0 else np.nan

def Q_val(K1, K2, K3):
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    return s**2 / (K1 + K2 + K3)

# PDG
M_E   = 0.51099895
M_MU  = 105.6583755
M_TAU = 1776.86
R21_PDG = M_MU / M_E
R31_PDG = M_TAU / M_E

# Punkty
POINTS = {
    'B (P55/P58 ref)': dict(a=0.040049, alpha=8.5612,  lam=5.4677e-6),
    'A (nominalny)':   dict(a=0.040000, alpha=8.5445,  lam=5.4677e-6),
}

# Siatki do badania
N_LIST = [500, 1000, 2000, 3000, 5000, 8000, 10000]

# -------------------------------------------------------------------
print("="*72)
print("P59: ZBIEZNOSC SIATKI -- BADANIE OP-13")
print("="*72)

for pname, pdict in POINTS.items():
    a, alpha, lam = pdict['a'], pdict['alpha'], pdict['lam']
    K3 = K3_ei(a, lam)

    print(f"\n{'='*72}")
    print(f"  PUNKT: {pname}")
    print(f"  a={a}, alpha={alpha}, lambda={lam:.5e}")
    print(f"  K3_Ei = {K3:.10f}")
    print(f"{'='*72}")
    print(f"  {'N':>7}  {'K1':>18}  {'K2':>18}  {'r21':>14}  {'Q_Ei':>14}")
    print(f"  {'-'*7}  {'-'*18}  {'-'*18}  {'-'*14}  {'-'*14}")

    results_N = {}
    for N in N_LIST:
        K1 = find_K1(alpha, a, lam, N)
        K2 = find_K2(alpha, a, lam, N)
        if np.isnan(K1) or np.isnan(K2):
            print(f"  {N:7d}  {'BLAD':>18}  {'BLAD':>18}")
            continue
        r21 = K2/K1
        Q   = Q_val(K1, K2, K3)
        results_N[N] = dict(K1=K1, K2=K2, r21=r21, Q=Q)
        print(f"  {N:7d}  {K1:.14f}  {K2:.14f}  {r21:.10f}  {Q:.12f}")

    # Richardson extrapolation dla K1, K2
    # Zakladamy K(N) = K_inf + c/N^p
    # Uzywamy trojki punktow (N1, N2, N3) gdzie N3=2*N2=4*N1
    # Tutaj: (2000, 5000, 10000) -- nierowne, ale mozliwe
    if 2000 in results_N and 5000 in results_N and 10000 in results_N:
        # Richardson 2-punktowy z (5000, 10000)
        K1_5k  = results_N[5000]['K1']
        K1_10k = results_N[10000]['K1']
        K2_5k  = results_N[5000]['K2']
        K2_10k = results_N[10000]['K2']
        # zakladamy p=2 (trapezoid rule => O(h^2), h~1/N^2 bo log-grid)
        # Richardson: K_inf = K(2N) + (K(2N)-K(N))/(2^p - 1)
        # Ale N nie podwaja sie dokladnie. Uzyjmy prosta ekstrapolacje:
        # K_inf ~ K_10k + (K_10k - K_5k) * (10000/(10000-5000))
        # = K_10k + (K_10k - K_5k)  -- zbyt agresywne
        # Uzyjmy estymacje rzedu:
        K1_2k = results_N[2000]['K1']
        K2_2k = results_N[2000]['K2']

        # estymata rzedu z trojki 2000, 5000, 10000
        # p taki ze K(N) = K_inf + c/N^p
        # (K(2000)-K(5000))/(K(5000)-K(10000)) = (N3^p - N2^p) / (N2^p - N1^p) ... uproszczenie
        # Prosta metoda: log((K2k-K5k)/(K5k-K10k)) / log(5/2)  [proporcja krokow]
        d1_K1 = K1_2k - K1_5k
        d2_K1 = K1_5k - K1_10k
        d1_K2 = K2_2k - K2_5k
        d2_K2 = K2_5k - K2_10k

        if abs(d2_K1) > 1e-16 and d1_K1*d2_K1 > 0:
            p_K1 = np.log(abs(d1_K1/d2_K1)) / np.log(5000/2000)
        else:
            p_K1 = np.nan
        if abs(d2_K2) > 1e-16 and d1_K2*d2_K2 > 0:
            p_K2 = np.log(abs(d1_K2/d2_K2)) / np.log(5000/2000)
        else:
            p_K2 = np.nan

        print(f"\n  Estymata rzedu zbieznosci:")
        print(f"  p_K1 = {p_K1:.3f}  (oczekiwane: ~2.0 dla trapezoid na log-grid)")
        print(f"  p_K2 = {p_K2:.3f}")

        # Richardson z (5000, 10000) zakladajac p~2
        # h~1/N, wiec K(N) ~ K_inf + c/N^p
        # K_inf = (N3^p * K3 - N2^p * K2) / (N3^p - N2^p)
        N2, N3 = 5000, 10000
        pK1_use = p_K1 if not np.isnan(p_K1) else 2.0
        pK2_use = p_K2 if not np.isnan(p_K2) else 2.0
        K1_inf = (N3**pK1_use * K1_10k - N2**pK1_use * K1_5k) / (N3**pK1_use - N2**pK1_use)
        K2_inf = (N3**pK2_use * K2_10k - N2**pK2_use * K2_5k) / (N3**pK2_use - N2**pK2_use)

        print(f"\n  Richardson (N=5000,10000) ekstrapolacja N->inf:")
        print(f"  K1_inf = {K1_inf:.14f}")
        print(f"  K2_inf = {K2_inf:.14f}")
        r21_inf = K2_inf / K1_inf
        Q_inf   = Q_val(K1_inf, K2_inf, K3)
        print(f"  r21_inf = {r21_inf:.10f}  (PDG: {R21_PDG:.8f},  delta={r21_inf-R21_PDG:+.6f})")
        print(f"  Q_inf   = {Q_inf:.12f}   (Q-3/2 = {(Q_inf-1.5)*1e6:+.3f} ppm)")

        # Bledy wzgledne vs N_inf
        print(f"\n  Bledy wzgledne vs Richardson K_inf [ppm]:")
        print(f"  {'N':>7}  {'dK1 [ppm]':>12}  {'dK2 [ppm]':>12}  {'dr21 [ppm]':>12}  {'dQ [ppm]':>12}")
        print(f"  {'-'*7}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")
        for N in N_LIST:
            if N not in results_N:
                continue
            r = results_N[N]
            dK1  = (r['K1']/K1_inf - 1)*1e6
            dK2  = (r['K2']/K2_inf - 1)*1e6
            dr21 = (r['r21']/r21_inf - 1)*1e6
            Q_N  = Q_val(r['K1'], r['K2'], K3)
            dQ   = (Q_N - Q_inf)*1e6
            print(f"  {N:7d}  {dK1:+12.1f}  {dK2:+12.1f}  {dr21:+12.1f}  {dQ:+12.3f}")

        # Znajdz N* dla bledu < 100 ppm (0.01%) w K2
        print(f"\n  N* dla bledu K2 < 100 ppm (0.01%):")
        for N in N_LIST:
            if N not in results_N:
                continue
            dK2 = abs((results_N[N]['K2']/K2_inf - 1)*1e6)
            if dK2 < 100:
                print(f"  N* = {N}  (blad K2 = {dK2:.1f} ppm)")
                break
        else:
            print(f"  Wymagane N > {N_LIST[-1]}")

        # Sprawdz N=5000 (aktualne uzycie)
        dK1_5k = abs((K1_5k/K1_inf - 1)*1e6)
        dK2_5k = abs((K2_5k/K2_inf - 1)*1e6)
        print(f"\n  Status N=5000 (uzywane w P58):")
        print(f"  dK1 = {dK1_5k:.1f} ppm  ({'OK' if dK1_5k < 100 else 'ZBYT DUZY'})")
        print(f"  dK2 = {dK2_5k:.1f} ppm  ({'OK' if dK2_5k < 100 else 'ZBYT DUZY'})")

# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA: WPLYW SIATKI NA r21 vs PDG")
print(f"{'='*72}")
print("""
  Kluczowe pytanie: czy N=5000 wystarczy aby r21 < 1 ppm od PDG?

  Z powyzszej tabeli bledow mozna ocenic:
  - K1 jest bardzo slabo zalezne od N (zbiega szybko)
  - K2 jest bardziej czule (wieksza amplitude ruchu)
  - r21 ~ K2/K1, wiec glownie blad K2

  Jesli dK2(N=5000) < 100 ppm, to dr21(N=5000) < 100 ppm.
  Przy r21_PDG = 206.768, 100 ppm = 0.021 -- co odpowiada
  dokladnosci ~0.02 w r21, co jest ok na obecnym etapie.
""")

print(f"{'='*72}")
print("P59 ZAKONCZONY")
print(f"{'='*72}")
