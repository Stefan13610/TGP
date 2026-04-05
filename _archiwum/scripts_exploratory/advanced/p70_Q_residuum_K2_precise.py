# -*- coding: utf-8 -*-
"""
P70: Q=3/2 DOKLADNIE — RESIDUUM Z NOWYM K2 (OP-4)

Z P69: K2 z R[2/3]+Newton blad = +0.010 ppm @ B; globalnie max 14.6 ppm.
Z P56/P57: Q residuum = +17 ppm @ B, ale K2_Pade mial blad -0.20% (~2000 ppm).
Oczekiwanie: zastapienie K2_Pade przez K2_R23N powinno drastycznie zredukowac Q residuum.

Cel P70:
  1. Wyliczyc Q z K2_num (precyzja maszynowa) — baseline
  2. Wyliczyc Q z K2_R23N = R[2/3]+1 Newton — ile wynosi residuum?
  3. Wyliczyc Q z K2_Pade[2/2] — ile residuum? (reprodukcja P48/P56)
  4. Analiza: skad pochodzi niezerowe residuum Q?
     - K1 blad? K3 blad? Korelacje miedzy K_i?
  5. Czy istnieje wartosc K2* taka ze Q(K1, K2*, K3) = 3/2 dokladnie?
     K2* vs K2_num: roznica ~ ???
  6. Trajektoria Q(K2) przy stalych K1, K3 — dQ/dK2 = ?
  7. Badanie Q na calej siatce (a, alpha) z K2_R23N

Odpowiada OP-4: Q=3/2 dokladnie.
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq
from numpy.linalg import lstsq

R_MAX  = 50.0
N_GRID = 10000

def E1(x):
    return -expi(-x)

def V_mod(phi, lam):
    return phi**3/3.0 - phi**4/4.0 + lam*(phi-1.0)**6/6.0

def energy_num(K, alpha, a, lam, N=N_GRID):
    t   = np.linspace(0, 1, N)
    r   = a*(R_MAX/a)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi= K*np.exp(-r)*(-r-1.0)/r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi*np.trapezoid(0.5*dphi**2*(1.0+alpha/phi)*r**2, r)
    Ep  = 4*np.pi*np.trapezoid((V_mod(phi,lam)-V1)*r**2, r)
    return Ek+Ep

def g_func(K, alpha, a, lam, N=N_GRID):
    return energy_num(K,alpha,a,lam,N)/(4*np.pi*K) - 1.0

def find_K(bracket, alpha, a, lam, N=N_GRID, tol=1e-12):
    try: return brentq(lambda K: g_func(K,alpha,a,lam,N),
                       bracket[0], bracket[1], xtol=tol)
    except: return np.nan

def koide_Q(K1, K2, K3):
    """Q = (sqrt(K1)+sqrt(K2)+sqrt(K3))^2 / (K1+K2+K3). Koide: Q=3/2."""
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    return s**2 / (K1 + K2 + K3)

def rational_interp(K_pts, g_pts, n, m):
    N_pts = len(K_pts)
    A = np.zeros((N_pts, n+m+1))
    b = np.zeros(N_pts)
    for i, (K, g) in enumerate(zip(K_pts, g_pts)):
        for j in range(n+1): A[i,j] = K**j
        for j in range(1,m+1): A[i,n+j] = -g*K**j
        b[i] = g
    x, _, _, _ = lstsq(A, b, rcond=None)
    p = x[:n+1]; q = np.zeros(m+1); q[0]=1.0; q[1:]=x[n+1:]
    return p, q

def find_zero_poly(p, K_low=0.4, K_high=5.0):
    def poly(K): return sum(p[j]*K**j for j in range(len(p)))
    K_s = np.linspace(K_low, K_high, 10000)
    vs  = np.array([poly(K) for K in K_s])
    zeros = []
    for i in range(len(K_s)-1):
        if vs[i]*vs[i+1] < 0:
            try: zeros.append(brentq(poly, K_s[i], K_s[i+1], xtol=1e-14))
            except: pass
    return zeros

def K2_from_R23_Newton(alpha, a, lam, n_newton=1):
    """K2 z R[2/3] + n_newton krokow Newtona."""
    K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
    K_MID  = np.linspace(1.0, 3.0, 9)
    K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
    K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))
    g_fit  = np.array([g_func(K, alpha, a, lam) for K in K_FIT])
    p, q   = rational_interp(K_FIT, g_fit, 2, 3)
    zeros  = find_zero_poly(p)
    if not zeros: return np.nan
    K_c = min(zeros, key=lambda z: abs(z - 2.0))
    for _ in range(n_newton):
        F_c  = g_func(K_c, alpha, a, lam)
        dF_c = (g_func(K_c+1e-6, alpha, a, lam) -
                g_func(K_c-1e-6, alpha, a, lam)) / 2e-6
        if abs(dF_c) < 1e-20: break
        K_c -= F_c / dF_c
    return K_c

def K2_pade22(alpha, a, lam):
    """K2 z R[2/2] na 21 punktach g(K) (odpowiednik P48 Pade[2/2])."""
    K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
    K_MID  = np.linspace(1.0, 3.0, 9)
    K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
    K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))
    g_fit  = np.array([g_func(K, alpha, a, lam) for K in K_FIT])
    p, q   = rational_interp(K_FIT, g_fit, 2, 2)
    zeros  = find_zero_poly(p)
    if not zeros: return np.nan
    return min(zeros, key=lambda z: abs(z - 2.0))

# Parametry Punkt B
A0   = 0.040049
ALP0 = 8.5612
LAM0 = 5.4677e-6
Q_TARGET = 3.0/2.0

print("="*72)
print("P70: Q=3/2 DOKLADNIE — RESIDUUM Z NOWYM K2 (OP-4)")
print("="*72)
print(f"  Punkt B: a={A0}, alpha={ALP0}, lambda={LAM0:.5e}")

# Wyznacz K1, K2, K3 numerycznie z wysoka precyzja (N=10000)
K1_num = find_K([0.001, 0.1],  ALP0, A0, LAM0)
K2_num = find_K([0.4,   5.0],  ALP0, A0, LAM0)
K3_num = find_K([5.0,  100.0], ALP0, A0, LAM0)
Q_num  = koide_Q(K1_num, K2_num, K3_num)

print(f"\n  K1_num = {K1_num:.14f}")
print(f"  K2_num = {K2_num:.14f}")
print(f"  K3_num = {K3_num:.14f}")
print(f"  Q_num  = {Q_num:.14f}  (residuum = {(Q_num-Q_TARGET)*1e6:+.4f} ppm od 3/2)")

# =========================================================
# SEKCJA 1: Trzy metody K2 -> Q
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Q z roznych przyblizen K2")
print(f"{'='*72}")

# (a) K2 z R[2/3] + 1 Newton
K2_R23N1 = K2_from_R23_Newton(ALP0, A0, LAM0, n_newton=1)
K2_R23N2 = K2_from_R23_Newton(ALP0, A0, LAM0, n_newton=2)
K2_P22   = K2_pade22(ALP0, A0, LAM0)

methods = [
    ('K2_num (referencja)',   K2_num),
    ('K2_R23+1N (+0.010 ppm)', K2_R23N1),
    ('K2_R23+2N (~0 ppm)',     K2_R23N2),
    ('K2_Pade22 (~-2000 ppm)', K2_P22),
]

print(f"\n  {'Metoda':30s}  {'K2':18s}  {'dK2 (ppm)':12s}  {'Q':18s}  {'Q-3/2 (ppm)':12s}")
for name, K2_v in methods:
    if np.isnan(K2_v): continue
    Q_v    = koide_Q(K1_num, K2_v, K3_num)
    dK2    = (K2_v/K2_num - 1)*1e6
    dQ     = (Q_v - Q_TARGET)*1e6
    print(f"  {name:30s}  {K2_v:.12f}  {dK2:+12.4f}  {Q_v:.12f}  {dQ:+12.4f}")

# =========================================================
# SEKCJA 2: Analiza zrodla residuum Q (przy K2_num)
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Skad pochodzi niezerowe Q_num - 3/2?")
print(f"{'='*72}")

Q_res = (Q_num - Q_TARGET)*1e6
print(f"\n  Q_num - 3/2 = {Q_res:+.6f} ppm")
print(f"\n  Analiza: K_i = K_i_num (numeryczna siatka N=10000)")
print(f"  Blad K_i wynika z:")
print(f"  - Dyskretyzacja siatki: ~0.2 ppm w K2 (P59: N*=500 wystarczy)")
print(f"  - Aproksymacja profilu phi = 1 + K*e^(-r)/r (model Yukawa)")

# Zbadaj wrazliwosc Q na K_i
dK = 1e-7
dQ_dK1 = (koide_Q(K1_num+dK, K2_num, K3_num) - koide_Q(K1_num-dK, K2_num, K3_num)) / (2*dK)
dQ_dK2 = (koide_Q(K1_num, K2_num+dK, K3_num) - koide_Q(K1_num, K2_num-dK, K3_num)) / (2*dK)
dQ_dK3 = (koide_Q(K1_num, K2_num, K3_num+dK) - koide_Q(K1_num, K2_num, K3_num-dK)) / (2*dK)

print(f"\n  Wrazliwosc Q na K_i:")
print(f"  dQ/dK1 = {dQ_dK1:.8e}  (K1={K1_num:.6f})")
print(f"  dQ/dK2 = {dQ_dK2:.8e}  (K2={K2_num:.6f})")
print(f"  dQ/dK3 = {dQ_dK3:.8e}  (K3={K3_num:.6f})")

# Jaka korekta K2 redukuje Q do 3/2?
dK2_needed = (Q_TARGET - Q_num) / dQ_dK2
print(f"\n  Korekcja K2 potrzebna do Q=3/2 dokladnie:")
print(f"  dK2* = {dK2_needed:.8e}")
print(f"  K2*  = K2_num + dK2* = {K2_num + dK2_needed:.14f}")
print(f"  dK2* w ppm = {dK2_needed/K2_num*1e6:+.6f} ppm")

# =========================================================
# SEKCJA 3: Czy Q=3/2 przy (a_Gamma, alpha_K, lambda_K)?
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Q jako funkcja lambda przy stalych (a, alpha)")
print(f"{'='*72}")

# Szukamy lambda* takie ze Q=3/2 dokladnie
# (P55/P56: lambda* bliska lambda_K)
lam_vals = LAM0 * np.linspace(0.98, 1.02, 9)
print(f"\n  {'lambda':16s}  {'lambda/lambda_K':14s}  {'K1':14s}  {'K2':14s}  {'K3':14s}  {'Q-3/2 (ppm)':12s}")
for lam_v in lam_vals:
    K1_v = find_K([0.001, 0.1],  ALP0, A0, lam_v)
    K2_v = find_K([0.4,   5.0],  ALP0, A0, lam_v)
    K3_v = find_K([5.0,  100.0], ALP0, A0, lam_v)
    if any(np.isnan([K1_v, K2_v, K3_v])): continue
    Q_v = koide_Q(K1_v, K2_v, K3_v)
    print(f"  {lam_v:.8e}  {lam_v/LAM0:14.8f}  {K1_v:.10f}  {K2_v:.10f}  {K3_v:.10f}  {(Q_v-Q_TARGET)*1e6:+12.4f}")

# Szukamy lambda* = zero Q(lambda)-3/2
try:
    def Q_vs_lam(lam_v):
        K1_v = find_K([0.001, 0.1],  ALP0, A0, lam_v)
        K2_v = find_K([0.4,   5.0],  ALP0, A0, lam_v)
        K3_v = find_K([5.0,  100.0], ALP0, A0, lam_v)
        if any(np.isnan([K1_v, K2_v, K3_v])): return np.nan
        return koide_Q(K1_v, K2_v, K3_v) - Q_TARGET

    lam_star = brentq(Q_vs_lam, LAM0*0.99, LAM0*1.01, xtol=1e-20, rtol=1e-14)
    print(f"\n  lambda* (Q=3/2 dokladnie): {lam_star:.12e}")
    print(f"  lambda_K:                  {LAM0:.12e}")
    print(f"  lambda*/lambda_K - 1 =     {lam_star/LAM0 - 1:.8e}  ({(lam_star/LAM0-1)*1e6:+.4f} ppm)")

    K1_star = find_K([0.001, 0.1],  ALP0, A0, lam_star)
    K2_star = find_K([0.4,   5.0],  ALP0, A0, lam_star)
    K3_star = find_K([5.0,  100.0], ALP0, A0, lam_star)
    Q_star  = koide_Q(K1_star, K2_star, K3_star)
    print(f"  Q(lambda*) = {Q_star:.14f}  (residuum: {(Q_star-Q_TARGET)*1e6:+.4f} ppm)")
    print(f"  K1* = {K1_star:.12f}  (vs K1_num {K1_num:.12f}, delta={(K1_star/K1_num-1)*1e6:+.2f} ppm)")
    print(f"  K2* = {K2_star:.12f}  (vs K2_num {K2_num:.12f}, delta={(K2_star/K2_num-1)*1e6:+.2f} ppm)")
    print(f"  K3* = {K3_star:.12f}  (vs K3_num {K3_num:.12f}, delta={(K3_star/K3_num-1)*1e6:+.2f} ppm)")
    r21_star = K2_star/K1_star
    r31_star = K3_star/K1_star
    print(f"  r21* = {r21_star:.6f}  (vs r21_PDG=206.7683, delta={(r21_star/206.7683-1)*1e6:+.2f} ppm)")
    print(f"  r31* = {r31_star:.4f}")
except Exception as e:
    print(f"  Szukanie lambda*: {e}")

# =========================================================
# SEKCJA 4: Q na siatce (a, alpha) z K2_R23N i K2_num
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Q na siatce (a, alpha) — porownanie K2_num vs K2_R23N")
print(f"{'='*72}")

a_grid   = np.array([0.035, 0.040049, 0.046, 0.055])
alp_grid = np.array([7.0, 8.0, 8.5612, 9.5, 11.0])

print(f"\n  {'a':>8} {'alpha':>7}  {'Q_num (ppm)':>14}  {'Q_R23N1 (ppm)':>16}  {'dQ (ppm)':>10}  {'K2 blad (ppm)':>16}")
for a_v in a_grid:
    for alp_v in alp_grid:
        K1_v = find_K([0.001, 0.1],  alp_v, a_v, LAM0)
        K2_v = find_K([0.4,   5.0],  alp_v, a_v, LAM0)
        K3_v = find_K([5.0,  100.0], alp_v, a_v, LAM0)
        if any(np.isnan([K1_v, K2_v, K3_v])): continue

        K2_apx = K2_from_R23_Newton(alp_v, a_v, LAM0, n_newton=1)
        if np.isnan(K2_apx): continue

        Q_v   = koide_Q(K1_v, K2_v, K3_v)
        Q_apx = koide_Q(K1_v, K2_apx, K3_v)
        dQ    = (Q_apx - Q_v)*1e6
        dK2   = (K2_apx/K2_v - 1)*1e6

        print(f"  {a_v:.5f} {alp_v:.4f}  {(Q_v-Q_TARGET)*1e6:+14.4f}  "
              f"{(Q_apx-Q_TARGET)*1e6:+16.4f}  {dQ:+10.4f}  {dK2:+16.4f}")

# =========================================================
# SEKCJA 5: Trajektoria Q(K2) przy stalych K1, K3
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Trajektoria Q(K2) — dQ/dK2 i K2* dla Q=3/2")
print(f"{'='*72}")

K2_range = np.linspace(K2_num*0.998, K2_num*1.002, 21)
print(f"\n  {'K2':20s}  {'dK2 (ppm)':12s}  {'Q':18s}  {'Q-3/2 (ppm)':12s}")
for K2_v in K2_range:
    Q_v = koide_Q(K1_num, K2_v, K3_num)
    dK2 = (K2_v/K2_num - 1)*1e6
    dQ  = (Q_v - Q_TARGET)*1e6
    print(f"  {K2_v:.14f}  {dK2:+12.4f}  {Q_v:.14f}  {dQ:+12.6f}")

# Analityczne dQ/dK2
# Q = (sqrt(K1)+sqrt(K2)+sqrt(K3))^2 / (3*(K1+K2+K3))
# dQ/dK2 = [2S*(1/(2*sqrt(K2)))*D - S^2] / (3*D^2)
#         gdzie S = sqrt(K1)+sqrt(K2)+sqrt(K3), D = K1+K2+K3
S = np.sqrt(K1_num) + np.sqrt(K2_num) + np.sqrt(K3_num)
D = K1_num + K2_num + K3_num
dQ_dK2_an = (2*S*(1/(2*np.sqrt(K2_num)))*D - S**2) / (3*D**2)
print(f"\n  dQ/dK2 (analityczne) = {dQ_dK2_an:.10e}")
print(f"  dQ/dK2 (numeryczne)  = {dQ_dK2:.10e}")
print(f"  Zgodnosc: delta = {abs(dQ_dK2_an - dQ_dK2)/abs(dQ_dK2)*100:.6f}%")

# K2* dla Q=3/2 (analityczne)
Q_num_val = koide_Q(K1_num, K2_num, K3_num)
K2_star_lin = K2_num + (Q_TARGET - Q_num_val) / dQ_dK2_an
print(f"\n  K2* (linearyzacja)   = {K2_star_lin:.14f}")
print(f"  K2_num               = {K2_num:.14f}")
print(f"  K2* - K2_num (ppm)   = {(K2_star_lin/K2_num - 1)*1e6:+.6f} ppm")
Q_check = koide_Q(K1_num, K2_star_lin, K3_num)
print(f"  Q(K2*)               = {Q_check:.14f}  (residuum: {(Q_check-Q_TARGET)*1e6:+.6f} ppm)")

# =========================================================
# SEKCJA 6: Podsumowanie OP-4
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Podsumowanie — diagnoza OP-4")
print(f"{'='*72}")

print(f"""
  WYNIKI (Punkt B, a={A0}, alpha={ALP0}, lambda={LAM0:.5e}):
  =============================================================

  Q z roznych K2:
    K2_Pade22 (-2000 ppm): Q = {koide_Q(K1_num, K2_P22, K3_num)*3-3:.4f}*1/2 + ?
    K2_R23+1N (+0.010 ppm): Q = {koide_Q(K1_num, K2_R23N1, K3_num):.10f}  (Q-3/2 = {(koide_Q(K1_num, K2_R23N1, K3_num)-Q_TARGET)*1e6:+.4f} ppm)
    K2_R23+2N (~0 ppm):     Q = {koide_Q(K1_num, K2_R23N2, K3_num):.10f}  (Q-3/2 = {(koide_Q(K1_num, K2_R23N2, K3_num)-Q_TARGET)*1e6:+.4f} ppm)
    K2_num (baseline):      Q = {Q_num:.10f}  (Q-3/2 = {(Q_num-Q_TARGET)*1e6:+.4f} ppm)

  Wniosek OP-4:
    Q residuum = {Q_res:+.4f} ppm NIE pochodzi od bledu K2_Pade.
    Nawet z K2 precyzji maszynowej Q-3/2 = {(Q_num-Q_TARGET)*1e6:+.4f} ppm ≠ 0.
    Residuum pochodzi z siatki numerycznej (N=10000, ~0.2 ppm w K2)
    lub z modelu (profil Yukawa phi = 1+K*e^(-r)/r).

  dQ/dK2 = {dQ_dK2:.6e}
  dK2 dla Q=3/2: {dK2_needed/K2_num*1e6:+.4f} ppm — mniejsze niz blad siatki!

  STATUS OP-4:
    OP-4 jest POCHODNA OP-1, ale takze POCHODNA dokladnosci siatki.
    Z K2 z R[2/3]+Newton (0.010 ppm): Q residuum ~= Q_num residuum.
    Zrodlo residuum Q: siatka numeryczna (P59: Richardson daje Q_inf~-1.28 ppm).
    Aby Q=3/2 dokladnie: potrzeba lambda* != lambda_K (inna wartosc!).
""")

print("="*72)
print("P70 ZAKONCZONY")
print("="*72)
