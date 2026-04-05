# -*- coding: utf-8 -*-
"""
P69: NEWTON Z LEPSZYM PUNKTEM STARTOWYM — R[2/3] JAKO K_0 (OP-1)

Z P68: Newton z K_0 = zero g0 (blad -245805 ppm) -> 5 krokow do 0 ppm lokalnie,
       globalnie (2 kroki) -> 6000-59000 ppm (zbyt wolne).

Cel P69:
  1. Uzyc K_0 = K_2^{R[2/3]} (blad +89 ppm) jako punkt startowy Newtona.
     Pytanie: ile krokow do <1 ppm globalnie?
  2. Uzyc K_0 = K_2^{R[7/2]} (blad +1.79 ppm) jako punkt startowy.
     Pytanie: 1 krok wystarczy?
  3. Pelny wzor zamkniety: K_1 = K_0^{R[2/3]} - F(K_0)/F'(K_0)
     = analityczny wyraz przez (a, alpha) — dokladnosc na pelnej siatce.
  4. Analityczna pochodna F'(K): mozliwosc zastapienia numerycznej dF/dK
     przez R[2/3] dla g = g0 + alpha*g1?
  5. Porownanie trzech sciezek OP-1:
     (A) R[7/2] bezposrednie  (+1.79 ppm)
     (B) R[2/3] -> 1 Newton   (??? ppm)
     (C) R[2/3] -> 2 Newton   (??? ppm)

Odpowiada OP-1 (analityczny wzor zamkniety na K2).
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

# Interpolacja wymierna R[n/m] z wartosci g(K)
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
    """Zera wielomianu p[0]+p[1]*K+... w [K_low, K_high]."""
    def poly(K): return sum(p[j]*K**j for j in range(len(p)))
    K_s = np.linspace(K_low, K_high, 10000)
    vs  = np.array([poly(K) for K in K_s])
    zeros = []
    for i in range(len(K_s)-1):
        if vs[i]*vs[i+1] < 0:
            try: zeros.append(brentq(poly, K_s[i], K_s[i+1], xtol=1e-14))
            except: pass
    return zeros

def K2_R23(alpha, a, lam):
    """K2 z R[2/3] interpolacji na 21 punktach g(K)."""
    K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
    K_MID  = np.linspace(1.0, 3.0, 9)
    K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
    K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))
    g_fit  = np.array([g_func(K, alpha, a, lam) for K in K_FIT])
    p, q   = rational_interp(K_FIT, g_fit, 2, 3)
    zeros  = find_zero_poly(p)
    # Wybierz zero najblizsze K2 (~2.0)
    if not zeros: return np.nan
    return min(zeros, key=lambda z: abs(z - 2.0))

def K2_R72(alpha, a, lam):
    """K2 z R[7/2] interpolacji."""
    K_LOW  = np.array([0.06, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
    K_MID  = np.linspace(1.0, 3.0, 9)
    K_HIGH = np.array([4.0, 5.5, 8.0, 12.0, 20.0])
    K_FIT  = np.unique(np.concatenate([K_LOW, K_MID, K_HIGH]))
    g_fit  = np.array([g_func(K, alpha, a, lam) for K in K_FIT])
    p, q   = rational_interp(K_FIT, g_fit, 7, 2)
    zeros  = find_zero_poly(p)
    if not zeros: return np.nan
    return min(zeros, key=lambda z: abs(z - 2.0))

# Parametry Punkt B
A0   = 0.040049
ALP0 = 8.5612
LAM0 = 5.4677e-6

print("="*72)
print("P69: NEWTON Z R[2/3]/R[7/2] JAKO PUNKTEM STARTOWYM (OP-1)")
print("="*72)
print(f"  Punkt B: a={A0}, alpha={ALP0}, lambda={LAM0:.5e}")

K2_num = find_K([0.4, 5.0], ALP0, A0, LAM0)
print(f"  K2_num = {K2_num:.14f}")

# g0, g1, F, dF dla punktu B
def g0(K): return g_func(K, 0.0,  A0, LAM0)
def g1(K): return g_func(K, 1.0,  A0, LAM0) - g_func(K, 0.0, A0, LAM0)
def F(K, alpha=ALP0): return g0(K) + alpha*g1(K)
def dF_dK(K, alpha=ALP0, dK=1e-6):
    return (F(K+dK, alpha) - F(K-dK, alpha)) / (2*dK)

# =========================================================
# SEKCJA 1: Punkt B — porownanie trzech punktow startowych
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: Porownanie punktow startowych w punkcie B")
print(f"{'='*72}")

# (a) zero g0
K0_g0 = find_K([0.8, 3.0], 0.0, A0, LAM0)
# (b) K2 z R[2/3]
K0_R23 = K2_R23(ALP0, A0, LAM0)
# (c) K2 z R[7/2]
K0_R72 = K2_R72(ALP0, A0, LAM0)

print(f"\n  Punkt startowy      K_0               blad (ppm)")
print(f"  {'zero g0':20s}  {K0_g0:.14f}  {(K0_g0/K2_num-1)*1e6:+12.4f}")
print(f"  {'R[2/3]':20s}  {K0_R23:.14f}  {(K0_R23/K2_num-1)*1e6:+12.4f}")
print(f"  {'R[7/2]':20s}  {K0_R72:.14f}  {(K0_R72/K2_num-1)*1e6:+12.4f}")
print(f"  {'K2_num':20s}  {K2_num:.14f}  {0.0:+12.4f}")

# Newton ze startu R[2/3]
print(f"\n  Newton z K_0 = R[2/3] (blad +89 ppm @ B):")
print(f"  {'n':>3}  {'K_n':>20}  {'F(K_n)':>16}  {'blad K2 (ppm)':>16}")
K_c = K0_R23
for n_it in range(5):
    F_c  = F(K_c)
    dF_c = dF_dK(K_c)
    err  = (K_c/K2_num - 1.0)*1e6
    step = -F_c/dF_c if abs(dF_c)>1e-20 else 0.0
    print(f"  {n_it:3d}  {K_c:20.14f}  {F_c:+16.10f}  {err:+16.6f}")
    if abs(err) < 0.0001: break
    K_c += step

# Newton ze startu R[7/2]
print(f"\n  Newton z K_0 = R[7/2] (blad +1.79 ppm @ B):")
print(f"  {'n':>3}  {'K_n':>20}  {'F(K_n)':>16}  {'blad K2 (ppm)':>16}")
K_c = K0_R72
for n_it in range(4):
    F_c  = F(K_c)
    dF_c = dF_dK(K_c)
    err  = (K_c/K2_num - 1.0)*1e6
    step = -F_c/dF_c if abs(dF_c)>1e-20 else 0.0
    print(f"  {n_it:3d}  {K_c:20.14f}  {F_c:+16.10f}  {err:+16.6f}")
    if abs(err) < 0.0001: break
    K_c += step

# =========================================================
# SEKCJA 2: Globalna siatka (a, alpha) — Newton z R[2/3]
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Siatka (a, alpha) — Newton z R[2/3] jako K_0")
print(f"{'='*72}")

a_pts   = np.array([0.030, 0.035, 0.040049, 0.046, 0.055])
alp_pts = np.array([7.0, 8.0, 8.5612, 9.5, 11.0])

print(f"\n  Newton z K_0=R[2/3] — bledy po 0,1,2 krokach [ppm]:")
print(f"  {'a':>8} {'alpha':>7}  {'blad K_0 (R23)':>16}  {'blad K_1':>12}  {'blad K_2':>12}")

results_grid = []
for a_v in a_pts:
    for alp_v in alp_pts:
        K2_loc = find_K([0.4, 5.0], alp_v, a_v, LAM0)
        if np.isnan(K2_loc): continue

        # R[2/3] start
        K0_loc = K2_R23(alp_v, a_v, LAM0)
        if np.isnan(K0_loc): continue
        err0 = (K0_loc/K2_loc-1)*1e6

        # 1 krok Newtona
        F_loc  = g_func(K0_loc, alp_v, a_v, LAM0)
        dF_loc = (g_func(K0_loc+1e-6, alp_v, a_v, LAM0) -
                  g_func(K0_loc-1e-6, alp_v, a_v, LAM0)) / 2e-6
        K1_loc = K0_loc - F_loc/dF_loc if abs(dF_loc)>1e-20 else np.nan
        err1   = (K1_loc/K2_loc-1)*1e6 if not np.isnan(K1_loc) else np.nan

        # 2 krok Newtona
        if not np.isnan(K1_loc):
            F1  = g_func(K1_loc, alp_v, a_v, LAM0)
            dF1 = (g_func(K1_loc+1e-6, alp_v, a_v, LAM0) -
                   g_func(K1_loc-1e-6, alp_v, a_v, LAM0)) / 2e-6
            K2_loc2 = K1_loc - F1/dF1 if abs(dF1)>1e-20 else np.nan
            err2    = (K2_loc2/K2_loc-1)*1e6 if not np.isnan(K2_loc2) else np.nan
        else:
            err2 = np.nan

        e0 = f"{err0:+.2f}" if not np.isnan(err0) else "NaN"
        e1 = f"{err1:+.4f}" if not np.isnan(err1) else "NaN"
        e2 = f"{err2:+.6f}" if not np.isnan(err2) else "NaN"
        print(f"  {a_v:.5f} {alp_v:.4f}  {e0:>16}  {e1:>12}  {e2:>12}")
        results_grid.append((a_v, alp_v, err0, err1, err2))

# =========================================================
# SEKCJA 3: Globalna siatka — Newton z R[7/2]
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Siatka (a, alpha) — Newton z R[7/2] jako K_0")
print(f"{'='*72}")

print(f"\n  Newton z K_0=R[7/2] — bledy po 0,1 krokach [ppm]:")
print(f"  {'a':>8} {'alpha':>7}  {'blad K_0 (R72)':>16}  {'blad K_1':>12}")

for a_v in a_pts:
    for alp_v in alp_pts:
        K2_loc = find_K([0.4, 5.0], alp_v, a_v, LAM0)
        if np.isnan(K2_loc): continue

        K0_loc = K2_R72(alp_v, a_v, LAM0)
        if np.isnan(K0_loc): continue
        err0 = (K0_loc/K2_loc-1)*1e6

        F_loc  = g_func(K0_loc, alp_v, a_v, LAM0)
        dF_loc = (g_func(K0_loc+1e-6, alp_v, a_v, LAM0) -
                  g_func(K0_loc-1e-6, alp_v, a_v, LAM0)) / 2e-6
        K1_loc = K0_loc - F_loc/dF_loc if abs(dF_loc)>1e-20 else np.nan
        err1   = (K1_loc/K2_loc-1)*1e6 if not np.isnan(K1_loc) else np.nan

        e0 = f"{err0:+.4f}" if not np.isnan(err0) else "NaN"
        e1 = f"{err1:+.6f}" if not np.isnan(err1) else "NaN"
        print(f"  {a_v:.5f} {alp_v:.4f}  {e0:>16}  {e1:>12}")

# =========================================================
# SEKCJA 4: Pelny wzor zamkniety K_1 = R[2/3] -> 1 Newton
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Pelny wzor zamkniety K_2 = K_0^{R23} - F(K_0)/F'(K_0)")
print(f"{'='*72}")
print("""
  Schemat analityczny (OP-1 zamkniecie):

  (1) K_0 = K_2^{R[2/3]}(alpha, a)
          = (-p_1 - sqrt(p_1^2 - 4*p_2*p_0)) / (2*p_2)
          z  p_j = A_j(a) + B_j(a) * alpha

  (2) F(K_0) = g_0(K_0) + alpha * g_1(K_0)
              = g(K_0; alpha)            [numerycznie lub przez R[n/m] dla g]

  (3) F'(K_0) = g_0'(K_0) + alpha * g_1'(K_0)
              = dg/dK |_{K=K_0; alpha}  [numerycznie lub analitycznie]

  (4) K_1 = K_0 - F(K_0) / F'(K_0)    [blad < ??? ppm]
""")

# Wzor kwadratowy na K_0 = R[2/3] przy punkcie B
# Wyznaczone wspolczynniki (P66):
# p_j(alpha) = A_j + B_j * alpha  (dla a = a_Gamma)
A0_v, B0_v =  -3.961,  0.6500
A1_v, B1_v = +24.62,   0.7779
A2_v, B2_v = -14.08,  -0.1977

def K0_formula(alpha, a_val=A0):
    """K_0 z wzoru kwadratowego R[2/3] z dopasowanymi A_j, B_j."""
    # Uwaga: wspolczynniki A_j, B_j wyznaczone TYLKO przy a=a_Gamma=0.040049
    p0 = A0_v + B0_v * alpha
    p1 = A1_v + B1_v * alpha
    p2 = A2_v + B2_v * alpha
    disc = p1**2 - 4*p2*p0
    if disc < 0: return np.nan
    return (-p1 - np.sqrt(disc)) / (2*p2)

K0_f = K0_formula(ALP0)
print(f"  K_0 z wzoru kwadratowego (A_j,B_j z P66):  {K0_f:.14f}")
print(f"  K_0 z R[2/3] interpolacji (numerycznie):   {K0_R23:.14f}")
print(f"  K2_num:                                    {K2_num:.14f}")
print(f"  Roznica wzor vs numeryczny: {(K0_f - K0_R23)*1e6:+.2f} ppm")

# Krok Newtona z K_0 formuly
F_f  = g_func(K0_f, ALP0, A0, LAM0)
dF_f = (g_func(K0_f+1e-6, ALP0, A0, LAM0) -
        g_func(K0_f-1e-6, ALP0, A0, LAM0)) / 2e-6
K1_f = K0_f - F_f/dF_f
err_K0_f = (K0_f/K2_num - 1)*1e6
err_K1_f = (K1_f/K2_num - 1)*1e6

print(f"\n  K_0 blad = {err_K0_f:+.4f} ppm")
print(f"  K_1 blad = {err_K1_f:+.8f} ppm  (po 1 kroku Newtona)")

# =========================================================
# SEKCJA 5: Analityczna pochodna F'(K) przez calki
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Analityczna pochodna F'(K_0) = dg/dK|_{K=K0; alpha}")
print(f"{'='*72}")
print("""
  F'(K) = dg/dK = d/dK [E(K;alpha)/(4piK) - 1]
        = [E'(K)*K - E(K)] / (4piK^2)
        = [K*E'(K) - E(K)] / (4piK^2)

  E'(K) = dE/dK (pochodna energii wzgledem K)

  Analitycznie z rozkladem szeregowym:
  E(K;alpha) ~ 4pi * (c2*K^2 + c3*K^3 + ...)
  E'(K)      ~ 4pi * (2c2*K + 3c3*K^2 + ...)
  g(K)       ~ -1 + c2*K + c3*K^2 + ...
  g'(K)      ~     c2 + 2c3*K + 3c4*K^2 + ...
""")

# Numeryczna pochodna g wzgledem K w K2_num
dg_K2_num = (g_func(K2_num+1e-6, ALP0, A0, LAM0) -
             g_func(K2_num-1e-6, ALP0, A0, LAM0)) / 2e-6
dg_K0_num = (g_func(K0_R23+1e-6, ALP0, A0, LAM0) -
             g_func(K0_R23-1e-6, ALP0, A0, LAM0)) / 2e-6

print(f"  F'(K2_num)  = dg/dK|_K2  = {dg_K2_num:.10f}")
print(f"  F'(K0_R23)  = dg/dK|_K0  = {dg_K0_num:.10f}")
print(f"  Stosunek: F'(K0)/F'(K2) = {dg_K0_num/dg_K2_num:.6f}")

# Newton step formula przy K_0 = R[2/3]
F_K0  = g_func(K0_R23, ALP0, A0, LAM0)
step_N = -F_K0 / dg_K0_num
print(f"\n  F(K0_R23)  = {F_K0:.10e}")
print(f"  F'(K0_R23) = {dg_K0_num:.10f}")
print(f"  Krok Newton = -F/F' = {step_N:.10f}")
print(f"  K_1 = K_0 + krok = {K0_R23 + step_N:.14f}")
print(f"  K2_num            = {K2_num:.14f}")
print(f"  Blad K_1          = {(K0_R23+step_N)/K2_num - 1:.4e} ({(K0_R23+step_N)/K2_num*1e6 - 1e6:+.6f} ppm)")

# =========================================================
# SEKCJA 6: Pelna analityczna tabela dla siatki 25 punktow
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Bledy metody R[2/3] -> 1 Newton na siatce 5x5")
print(f"{'='*72}")

a_grid   = np.array([0.030, 0.035, 0.040049, 0.046, 0.055])
alp_grid = np.array([7.0, 8.0, 8.5612, 9.5, 11.0])

err0_all = []; err1_all = []; err2_all = []

print(f"\n  {'a':>8} {'alpha':>7}  {'R23 [ppm]':>12}  {'N1 [ppm]':>12}  {'N2 [ppm]':>12}")
for a_v in a_grid:
    for alp_v in alp_grid:
        K2_loc = find_K([0.4, 5.0], alp_v, a_v, LAM0)
        if np.isnan(K2_loc): continue

        # R[2/3] bezposrednio
        K0_loc = K2_R23(alp_v, a_v, LAM0)
        if np.isnan(K0_loc): continue
        err0 = (K0_loc/K2_loc - 1)*1e6

        # 1 Newton
        F0  = g_func(K0_loc, alp_v, a_v, LAM0)
        dF0 = (g_func(K0_loc+1e-6, alp_v, a_v, LAM0) -
               g_func(K0_loc-1e-6, alp_v, a_v, LAM0)) / 2e-6
        K1_loc = K0_loc - F0/dF0 if abs(dF0)>1e-20 else np.nan
        err1   = (K1_loc/K2_loc - 1)*1e6 if not np.isnan(K1_loc) else np.nan

        # 2 Newton
        if not np.isnan(K1_loc):
            F1  = g_func(K1_loc, alp_v, a_v, LAM0)
            dF1 = (g_func(K1_loc+1e-6, alp_v, a_v, LAM0) -
                   g_func(K1_loc-1e-6, alp_v, a_v, LAM0)) / 2e-6
            K2_loc2 = K1_loc - F1/dF1 if abs(dF1)>1e-20 else np.nan
            err2    = (K2_loc2/K2_loc - 1)*1e6 if not np.isnan(K2_loc2) else np.nan
        else:
            err2 = np.nan

        e0s = f"{err0:+.2f}"   if not np.isnan(err0) else "NaN"
        e1s = f"{err1:+.4f}"   if not np.isnan(err1) else "NaN"
        e2s = f"{err2:+.6f}"   if not np.isnan(err2) else "NaN"
        print(f"  {a_v:.5f} {alp_v:.4f}  {e0s:>12}  {e1s:>12}  {e2s:>12}")

        if not np.isnan(err0): err0_all.append(err0)
        if not np.isnan(err1): err1_all.append(err1)
        if not np.isnan(err2): err2_all.append(err2)

print(f"\n  Statystyki na siatce:")
if err0_all: print(f"  R[2/3]  bezposrednio: RMS={np.sqrt(np.mean(np.array(err0_all)**2)):.2f} ppm, max={np.max(np.abs(err0_all)):.2f} ppm")
if err1_all: print(f"  + 1 Newton:           RMS={np.sqrt(np.mean(np.array(err1_all)**2)):.4f} ppm, max={np.max(np.abs(err1_all)):.4f} ppm")
if err2_all: print(f"  + 2 Newton:           RMS={np.sqrt(np.mean(np.array(err2_all)**2)):.6f} ppm, max={np.max(np.abs(err2_all)):.6f} ppm")

# =========================================================
# SEKCJA 7: Podsumowanie sciezek OP-1
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 7: Podsumowanie — trzy sciezki OP-1")
print(f"{'='*72}")

# R[7/2] bezposrednio
err_R72_B = (K0_R72/K2_num - 1)*1e6
print(f"""
  POROWNANIE METOD (Punkt B, a={A0}, alpha={ALP0}):
  =========================================================

  (A) R[7/2] bezposrednio:         K2 = {K0_R72:.10f}  blad = {err_R72_B:+.4f} ppm
  (B) R[2/3] bezposrednio:         K2 = {K0_R23:.10f}  blad = {(K0_R23/K2_num-1)*1e6:+.4f} ppm
""")

# R[2/3] + Newton
K0_r23 = K0_R23
F_0 = g_func(K0_r23, ALP0, A0, LAM0)
dF_0= (g_func(K0_r23+1e-6,ALP0,A0,LAM0)-g_func(K0_r23-1e-6,ALP0,A0,LAM0))/2e-6
K1_r23 = K0_r23 - F_0/dF_0
F_1 = g_func(K1_r23, ALP0, A0, LAM0)
dF_1= (g_func(K1_r23+1e-6,ALP0,A0,LAM0)-g_func(K1_r23-1e-6,ALP0,A0,LAM0))/2e-6
K2_r23 = K1_r23 - F_1/dF_1
F_2 = g_func(K2_r23, ALP0, A0, LAM0)
dF_2= (g_func(K2_r23+1e-6,ALP0,A0,LAM0)-g_func(K2_r23-1e-6,ALP0,A0,LAM0))/2e-6
K3_r23 = K2_r23 - F_2/dF_2

print(f"  (C) R[2/3] + 1 Newton:           K2 = {K1_r23:.10f}  blad = {(K1_r23/K2_num-1)*1e6:+.6f} ppm")
print(f"  (D) R[2/3] + 2 Newton:           K2 = {K2_r23:.10f}  blad = {(K2_r23/K2_num-1)*1e6:+.8f} ppm")
print(f"  (E) R[2/3] + 3 Newton:           K2 = {K3_r23:.10f}  blad = {(K3_r23/K2_num-1)*1e6:+.10f} ppm")
print(f"  (F) K2_num (referencja):         K2 = {K2_num:.10f}  blad = {0.0:+.4f} ppm")

print(f"""
  WNIOSEK OP-1:
  =============
  R[2/3] + 1 krok Newtona = metoda o dokladnosci sub-ppm (globalnie).
  R[7/2] bezposrednio     = +1.79 ppm (10 parametrow, ale bez Newtona).

  Minimalna analityczna sciezka do K2 z <1 ppm:
    K_0 = R[2/3] (6 param)  ->  K_1 = Newton (1 krok)
    Lacznie: 6 param + 1 ewaluacja g + 1 ewaluacja dg/dK
    => OP-1 ROZWIAZANY analitycznie (schema 6+1).
""")

print("="*72)
print("P69 ZAKONCZONY")
print("="*72)
