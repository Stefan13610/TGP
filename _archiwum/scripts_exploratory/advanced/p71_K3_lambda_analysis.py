# -*- coding: utf-8 -*-
"""
P71: ANALIZA K3 — SKAD K3_num - K3* = 2792 ppm? (OP-4 / OP-12)

Z P70:
  - Q@B = 3/2 - 632 ppm; zrodlo: K3_num = 34.291 za duze
  - K3* = 34.195 = K3_Ei (dla Q=3/2 dokladnie przy a_Gamma, alpha_K)
  - lambda*/lambda_K = +5616 ppm
  - K3* = K3_num * (lambda_K/lambda*)^(1/2) <-- weryfikacja

Cel P71:
  1. Weryfikacja: K3 ~ lambda^(-1/2) (dokładna potega)
  2. Potwierdzenie: lambda*/lambda_K = (K3_num/K3*)^2 = (34.291/34.195)^2
  3. Analityczny wzor K3 = sqrt(3*I4 / (2*lambda*I6)):
       I4 = e^(-4a)/a - 4*E1(4a)
       I6 = e^(-6a)/(3a^3) - e^(-6a)/a^2 + 6*e^(-6a)/a - 36*E1(6a)  (scisly wzor)
     Porownanie K3_an vs K3_num: ile ppm?
  4. I6: blad przybliżenia Yukawa (phi^6 ~ K^6*e^(-6r)/r^6)
     vs calka dokladna (numeryczna)
  5. K3 z formuly K3 = sqrt(3*I4_exact / (2*lambda*I6_exact)):
     czy to K3_Ei czy blizej K3_num?
  6. Trzy wartosci lambda:
       lambda_K   = 5.4677e-6  (uzywana w P68-P70)
       lambda*    = 5.4984e-6  (P70: Q=3/2 przy K_i numerycznych)
       lambda_an  = 3*I4 / (2*K3_num^2*I6_exact)  (z formuly analitycznej)
  7. Czy K3_Ei = K3* (dla Q=3/2) to koincydencja czy wynika z fizyki modelu?
  8. Pelna mapa Q(lambda) przy stalych (a, alpha)

Odpowiada OP-4 (Q=3/2) i OP-12 (K3 dokladnosc).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq

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
    s = np.sqrt(K1)+np.sqrt(K2)+np.sqrt(K3)
    return s**2/(K1+K2+K3)

# -------------------- Calki analityczne -----------------------
def I4_an(a):
    """I4 = integral_a^inf e^{-4r}/r dr = E1(4a)."""
    return E1(4*a)

def I6_an(a):
    """I6 = integral_a^inf e^{-6r}/r dr = E1(6a)."""
    return E1(6*a)

def I4_exact(a):
    """Calka z czlonu phi^4 w potencjale: 4pi * integral (phi-1)^4 r^2 dr.
    Dla phi = 1+K*e^{-r}/r: (phi-1)^4 = K^4*e^{-4r}/r^4.
    I4 uzyty w K3: lambda/6 * 4pi * integral K^6*e^{-6r}/r^6 r^2 dr.

    Z derivacji K3 (zerowanie E wzgledem K dla dominantu czlonu lambda):
    d/dK [lambda * 4pi * integral (phi-1)^6/6 r^2 dr] ~ lambda * K^5 * I6
    Balans: c2*K + ... + lambda*K^5*I6 ~ 0 => K3^2 ~ -c2-... / (lambda*I6)

    Scislej: dla dominantu K3>>1:
    E[K] ~ 4pi*(c2*K^2 + lambda*K^6*I6/6)
    g(K) = E/(4piK) - 1 ~ c2*K + lambda*K^5*I6/6 - 1 = 0
    => K3^2 ~ 1/(c2 + ...) * 6/(lambda*I6)   (dla dominujacego lambda)
    Lub z balansu: lambda*K^5*I6/6 ~ 1 => K3 ~ (6/(lambda*I6))^{1/5}

    Prostszy wzor (P61): K3 = sqrt(3*I4_Ei / (2*lambda*I6_Ei))
    gdzie I4_Ei = e^{-4a}/a - 4*E1(4a), I6_Ei = e^{-6a}/(3a^3) - e^{-6a}/a^2 + 6*e^{-6a}/a - 36*E1(6a)

    Uzywamy obu.
    """
    return np.exp(-4*a)/a - 4*E1(4*a)

def I6_Yukawa(a):
    """I6 = integral_a^inf e^{-6r}/r^4 dr (scisly wzor analityczny).
    Przez czesci trzykrotnie:
      integral e^{-6r}/r^4 dr = e^{-6a}/(3a^3) - e^{-6a}/a^2 + 6*e^{-6a}/a - 36*E1(6a)
    """
    ea = np.exp(-6*a)
    return ea/(3*a**3) - ea/a**2 + 6*ea/a - 36*E1(6*a)

def I6_exact_num(a, K3, N=5000):
    """Dokladna I6 = integral phi^6 / K^6 dr (bez Yukawa).
    phi(r) = 1 + K3*e^{-r}/r, wiec phi^6/K^6 =/= e^{-6r}/r^6.
    Uzyjemy calki numerycznej z dokladnym profilem.
    Ale K3 potrzebujemy -- iteracyjnie.
    Lepiej: I6_exact = integral (phi/K - e^{-r}/r)^... skomplikowane.
    Uzywamy brute force: oblicz E_lambda = 4pi*integral lam*(phi-1)^6/6*r^2 dr
    pochodna po K: dE_lambda/dK = 4pi*lam*integral (phi-1)^5*e^{-r}/r*r^2 dr
    """
    t    = np.linspace(0, 1, N)
    r    = a*(R_MAX/a)**t
    phi  = np.maximum(1.0 + K3*np.exp(-r)/r, 1e-10)
    dphi = np.exp(-r)/r   # dphi/dK = e^{-r}/r
    # (phi-1) = K3*e^{-r}/r
    phi_m1 = K3*np.exp(-r)/r
    # I6_exact: z balansu dg/dK = 0 wzgledem czlonu lambda:
    # d/dK [lambda * 4pi * integral (phi-1)^6/6 * r^2 dr]
    # = lambda * 4pi * integral (phi-1)^5 * e^{-r}/r * r^2 dr
    # = lambda * 4pi * K3^5 * I6_eff
    integrand = phi_m1**5 * np.exp(-r)/r * r**2
    I_num, _ = quad(lambda r_v: (K3*np.exp(-r_v)/r_v)**5 * np.exp(-r_v)/r_v * r_v**2,
                    a, R_MAX, limit=500)
    return I_num / K3**5  # I6_eff

def K3_formula(lam, a, I4_val, I6_val):
    """K3 z formuly K3 = sqrt(3*I4 / (2*lam*I6))."""
    val = 3*I4_val / (2*lam*I6_val)
    if val < 0: return np.nan
    return np.sqrt(val)

# Parametry
A0    = 0.040049
ALP0  = 8.5612
LAM0  = 5.4677e-6
LAM_STAR = 5.4984e-6  # z P70: lambda* dla Q=3/2

print("="*72)
print("P71: ANALIZA K3 — SKAD K3_num - K3* = 2792 ppm?")
print("="*72)
print(f"  a_Gamma = {A0}, alpha_K = {ALP0}")
print(f"  lambda_K = {LAM0:.6e}, lambda* = {LAM_STAR:.6e}")
print(f"  lambda*/lambda_K = {LAM_STAR/LAM0:.8f} (+{(LAM_STAR/LAM0-1)*1e6:.2f} ppm)")

# Wyznacz K1, K2, K3 przy lambda_K i lambda*
K1_lK = find_K([0.001,  0.1],  ALP0, A0, LAM0)
K2_lK = find_K([0.4,    5.0],  ALP0, A0, LAM0)
K3_lK = find_K([5.0,  100.0],  ALP0, A0, LAM0)
K1_ls = find_K([0.001,  0.1],  ALP0, A0, LAM_STAR)
K2_ls = find_K([0.4,    5.0],  ALP0, A0, LAM_STAR)
K3_ls = find_K([5.0,  100.0],  ALP0, A0, LAM_STAR)

print(f"\n  K_i @ lambda_K:  K1={K1_lK:.10f}, K2={K2_lK:.10f}, K3={K3_lK:.10f}")
print(f"  K_i @ lambda*:   K1={K1_ls:.10f}, K2={K2_ls:.10f}, K3={K3_ls:.10f}")
print(f"  Q(lambda_K) = {koide_Q(K1_lK,K2_lK,K3_lK):.10f}  ({(koide_Q(K1_lK,K2_lK,K3_lK)-1.5)*1e6:+.2f} ppm)")
print(f"  Q(lambda*)  = {koide_Q(K1_ls,K2_ls,K3_ls):.10f}  ({(koide_Q(K1_ls,K2_ls,K3_ls)-1.5)*1e6:+.2f} ppm)")

# =========================================================
# SEKCJA 1: K3 ~ lambda^(-1/2) — weryfikacja potegi
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 1: K3 ~ lambda^(-1/2) — weryfikacja dokladna")
print(f"{'='*72}")

lam_range = LAM0 * np.array([0.90, 0.95, 1.00, 1.005616, 1.01, 1.05, 1.10, 1.20])
K3_range  = []
print(f"\n  {'lambda':16s}  {'lambda/lambda_K':14s}  {'K3':18s}  {'K3/K3_ref':14s}  {'log(K3)':12s}")
for lam_v in lam_range:
    K3_v = find_K([5.0, 100.0], ALP0, A0, lam_v)
    K3_range.append(K3_v)
    ratio = K3_v/K3_lK if K3_lK > 0 else np.nan
    print(f"  {lam_v:.8e}  {lam_v/LAM0:14.6f}  {K3_v:.14f}  {ratio:14.10f}  {np.log(K3_v):12.8f}")

# Fit potegowy: log(K3) = beta * log(lambda) + const
lam_arr = lam_range
K3_arr  = np.array(K3_range)
valid   = ~np.isnan(K3_arr)
beta_K3 = np.polyfit(np.log(lam_arr[valid]), np.log(K3_arr[valid]), 1)
print(f"\n  Fit potegowy: K3 ~ lambda^beta,  beta = {beta_K3[0]:.8f}")
print(f"  Oczekiwane (P61): beta = -0.5000")
print(f"  Roznica: {abs(beta_K3[0]+0.5)*1e6:.4f} ppm od -1/2")

# Weryfikacja: lambda*/lambda_K = (K3_lK/K3_ls)^2
ratio_lam_sq = (K3_lK/K3_ls)**2
print(f"\n  Weryfikacja: (K3_num/K3*)^2 = ({K3_lK:.8f}/{K3_ls:.8f})^2")
print(f"  = {ratio_lam_sq:.8f}")
print(f"  lambda*/lambda_K = {LAM_STAR/LAM0:.8f}")
print(f"  Roznica = {abs(ratio_lam_sq - LAM_STAR/LAM0)*1e6:.4f} ppm")

# =========================================================
# SEKCJA 2: Calki I4, I6 — Yukawa vs dokladne
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 2: Calki I4, I6 — Yukawa vs numeryczne")
print(f"{'='*72}")

a = A0
I4_Ei  = I4_exact(a)    # Ei formula (= E1(4a) wg wzoru P61 - sprawdzic)
I4_Eip = np.exp(-4*a)/a - 4*E1(4*a)  # z wzoru P48/P61
I6_Ei_v = I6_Yukawa(a)  # scisly wzor analityczny (P71: blad 0.0000 ppm)

# Dokladna I4: integral_a^inf K^4*e^{-4r}/r^4 * r^2 dr = K^4 * integral e^{-4r}/r^2 dr
I4_num_raw, _ = quad(lambda r: np.exp(-4*r)/r**2, a, R_MAX, limit=500)
# (uzycie przy K3: E_V = 4pi*K^4 * integral (phi-1)^4*r^2 dr, phi-1=K*e^{-r}/r)
# = 4pi*K^4 * integral K^4*e^{-4r}/r^4 * r^2 dr = 4pi*K^8 * integral e^{-4r}/r^2 dr

# Wlasciwe I4, I6 uzywane w K3 (P48): z balansu dE/dK=0 przy K>>1
# E[K] ~ 4pi*(c2*K^2 + c3*K^3 + lam*K^6*I6/6)
# dE/dK = 4pi*(2c2*K + 3c3*K^2 + lam*K^5*I6) = 0 przy K3
# Zerowanie: lam*K^5*I6 ~ -(2c2*K + 3c3*K^2) -- ale to znowu Taylor...
# Prostsza forma (P61 wzor): K3 = sqrt(3*I4 / (2*lam*I6))
# gdzie I4 i I6 zdefiniowane jako:
#   I4 = integral_a^inf e^{-4r}/r^2 dr (numerycznie)
#   I6 = integral_a^inf e^{-6r}/r^4 dr (numerycznie)
#   Analyt: I4_Ei = e^{-4a}/a - 4*E1(4a), I6_Ei = e^{-6a}/(3a^3)-e^{-6a}/a^2+6*e^{-6a}/a-36*E1(6a)

I4_num_full, _ = quad(lambda r: np.exp(-4*r)/r**2, a, R_MAX, limit=500)
I6_num_full, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=500)

print(f"\n  a = {a}")
print(f"\n  Calka I4 = integral_a^inf e^(-4r)/r^2 dr:")
print(f"  I4_Yukawa = e^(-4a)/a - 4*E1(4a) = {I4_Eip:.12f}")
print(f"  I4_numeric                          = {I4_num_full:.12f}")
print(f"  Roznica I4: {(I4_Eip/I4_num_full-1)*1e6:+.4f} ppm")

print(f"\n  Calka I6 = integral_a^inf e^(-6r)/r^4 dr:")
print(f"  I6_analyt = e^(-6a)/(3a^3)-e^(-6a)/a^2+6e^(-6a)/a-36*E1(6a) = {I6_Ei_v:.12f}")
print(f"  I6_numeric                                       = {I6_num_full:.12f}")
print(f"  Roznica I6: {(I6_Ei_v/I6_num_full-1)*1e6:+.4f} ppm")

# K3 z formul
K3_Ei_an   = K3_formula(LAM0, a, I4_Eip, I6_Ei_v)
K3_I4num_I6Ei = K3_formula(LAM0, a, I4_num_full, I6_Ei_v)
K3_I4Ei_I6num = K3_formula(LAM0, a, I4_Eip, I6_num_full)
K3_num_an  = K3_formula(LAM0, a, I4_num_full, I6_num_full)

print(f"\n  K3 z formuly sqrt(3*I4/(2*lam*I6)) @ lambda_K:")
print(f"  K3(I4_Ei,   I6_Ei )  = {K3_Ei_an:.10f}  blad={( K3_Ei_an/K3_lK-1)*1e6:+.2f} ppm")
print(f"  K3(I4_num,  I6_Ei )  = {K3_I4num_I6Ei:.10f}  blad={( K3_I4num_I6Ei/K3_lK-1)*1e6:+.2f} ppm")
print(f"  K3(I4_Ei,   I6_num)  = {K3_I4Ei_I6num:.10f}  blad={( K3_I4Ei_I6num/K3_lK-1)*1e6:+.2f} ppm")
print(f"  K3(I4_num,  I6_num)  = {K3_num_an:.10f}  blad={( K3_num_an/K3_lK-1)*1e6:+.2f} ppm")
print(f"  K3_num (referencja)  = {K3_lK:.10f}  blad=0 ppm")
print(f"  K3* (P70, Q=3/2)     = {K3_ls:.10f}  blad={( K3_ls/K3_lK-1)*1e6:+.2f} ppm")

# =========================================================
# SEKCJA 3: Trzy wartosci lambda
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 3: Trzy wartosci lambda")
print(f"{'='*72}")

# lambda_K (z P70/P68)
# lambda* (P70: Q=3/2 przy K_i num)
# lambda_an (z formuly: lambda = 3*I4 / (2*K3_num^2*I6))
lam_an_EiEi   = 3*I4_Eip    / (2*K3_lK**2*I6_Ei_v)
lam_an_numnum = 3*I4_num_full/ (2*K3_lK**2*I6_num_full)
lam_an_Einumt = 3*I4_num_full/ (2*K3_lK**2*I6_Ei_v)
lam_an_numEi  = 3*I4_Eip    / (2*K3_lK**2*I6_num_full)

print(f"\n  lambda_K           = {LAM0:.10e}  (bazowa, P68-P70)")
print(f"  lambda* (Q=3/2)    = {LAM_STAR:.10e}  (+{(LAM_STAR/LAM0-1)*1e6:.2f} ppm)")
print(f"\n  lambda z formuly 3*I4/(2*K3_num^2*I6) @ K3_num={K3_lK:.6f}:")
print(f"  lambda(I4_Ei, I6_Ei)   = {lam_an_EiEi:.10e}  ({(lam_an_EiEi/LAM0-1)*1e6:+.2f} ppm od lambda_K)")
print(f"  lambda(I4_num,I6_num)  = {lam_an_numnum:.10e}  ({(lam_an_numnum/LAM0-1)*1e6:+.2f} ppm od lambda_K)")
print(f"  lambda(I4_num,I6_Ei)   = {lam_an_Einumt:.10e}  ({(lam_an_Einumt/LAM0-1)*1e6:+.2f} ppm od lambda_K)")
print(f"  lambda(I4_Ei, I6_num)  = {lam_an_numEi:.10e}  ({(lam_an_numEi/LAM0-1)*1e6:+.2f} ppm od lambda_K)")

# Kluczowe: lam_an_EiEi ~ lambda* ?
print(f"\n  Kluczowe porownanie:")
print(f"  lambda(I4_Ei, I6_Ei) = {lam_an_EiEi:.10e}")
print(f"  lambda*              = {LAM_STAR:.10e}")
print(f"  Roznica:               {(lam_an_EiEi/LAM_STAR-1)*1e6:+.4f} ppm")

# =========================================================
# SEKCJA 4: Skad pochodzi K3_num - K3_Ei = -2792 ppm?
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 4: Skad pochodzi K3_lK - K3_Ei = 2792 ppm?")
print(f"{'='*72}")

print(f"\n  K3_num (numeryczne) = {K3_lK:.10f}")
print(f"  K3_Ei  (analityczne)= {K3_Ei_an:.10f}")
print(f"  K3*    (Q=3/2)      = {K3_ls:.10f}")
print(f"  K3_num - K3_Ei = {(K3_lK - K3_Ei_an):.8f} = {(K3_lK/K3_Ei_an-1)*1e6:+.2f} ppm")
print(f"  K3_num - K3*   = {(K3_lK - K3_ls):.8f}  = {(K3_lK/K3_ls-1)*1e6:+.2f} ppm")
print(f"  K3_Ei  - K3*   = {(K3_Ei_an - K3_ls):.8f}  = {(K3_Ei_an/K3_ls-1)*1e6:+.4f} ppm")

# Analiza: K3_Ei(lambda_K) vs K3_Ei(lambda*)
K3_Ei_at_lstar = K3_formula(LAM_STAR, a, I4_Eip, I6_Ei_v)
K3_num_at_lstar = K3_ls  # juz wyznaczone
print(f"\n  K3_Ei(lambda_K) = {K3_Ei_an:.10f}")
print(f"  K3_Ei(lambda*)  = {K3_Ei_at_lstar:.10f}  (blad vs K3* = {(K3_Ei_at_lstar/K3_ls-1)*1e6:+.4f} ppm)")
print(f"  K3_num(lambda*) = {K3_num_at_lstar:.10f}")

# Co powoduje roznice K3_num vs K3_Ei przy tej samej lambda?
# K3_num spełnia: g(K3; alpha, a, lam) = 0 (pelna energia)
# K3_Ei pochodzi z: K3 = sqrt(3*I4/(2*lam*I6)) -- przybliżenie dominantu
# Pytanie: o ile formula K3=sqrt(3I4/(2lam*I6)) jest niedokladna?
print(f"\n  Sprawdzenie formuly K3=sqrt(3*I4/(2*lam*I6)) przy K3=K3_num:")
lam_from_K3num_EiEi   = 3*I4_Eip     / (2*K3_lK**2 * I6_Ei_v)
lam_from_K3num_numnum = 3*I4_num_full / (2*K3_lK**2 * I6_num_full)
print(f"  lam(I4_Ei, I6_Ei) z K3_num  = {lam_from_K3num_EiEi:.8e} ({(lam_from_K3num_EiEi/LAM0-1)*1e6:+.2f} ppm)")
print(f"  lam(I4_num,I6_num) z K3_num = {lam_from_K3num_numnum:.8e} ({(lam_from_K3num_numnum/LAM0-1)*1e6:+.2f} ppm)")
print(f"  lambda_K                    = {LAM0:.8e} (0 ppm)")

# =========================================================
# SEKCJA 5: Odchylenie K3_Ei od K3_num — kwantyfikacja bledu
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 5: Parametryzacja bledu formuly K3=sqrt(3*I4/(2*lam*I6))")
print(f"{'='*72}")

# Na siatce a:
a_vals = np.array([0.030, 0.035, 0.040049, 0.046, 0.055])
print(f"\n  {'a':>8}  {'K3_num':>16}  {'K3_Ei':>16}  {'K3_Ei/K3_num-1 (ppm)':>22}  {'I6_Ei/I6_num-1 (ppm)':>22}")
for a_v in a_vals:
    K3_v = find_K([5.0, 100.0], ALP0, a_v, LAM0)
    I4_v = np.exp(-4*a_v)/a_v - 4*E1(4*a_v)
    I6_Ei_loc = I6_Yukawa(a_v)
    I6_nm_loc, _ = quad(lambda r: np.exp(-6*r)/r**4, a_v, R_MAX, limit=500)
    K3_Ei_loc = K3_formula(LAM0, a_v, I4_v, I6_Ei_loc)
    if np.isnan(K3_v) or np.isnan(K3_Ei_loc): continue
    err_K3 = (K3_Ei_loc/K3_v - 1)*1e6
    err_I6 = (I6_Ei_loc/I6_nm_loc - 1)*1e6
    print(f"  {a_v:8.5f}  {K3_v:16.10f}  {K3_Ei_loc:16.10f}  {err_K3:+22.4f}  {err_I6:+22.4f}")

# =========================================================
# SEKCJA 6: Lambda "prawdziwe Koide" — co wyznacza Q=3/2?
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 6: Lambda Koide — pelna analiza")
print(f"{'='*72}")

# Przy danym (a, alpha), znajdz lambda takie ze Q=3/2
def find_lambda_Koide(a_v, alp_v, lam_guess=LAM0, tol=1e-18):
    def Q_minus_Koide(lam_v):
        K1_v = find_K([0.001, 0.1],  alp_v, a_v, lam_v)
        K2_v = find_K([0.4,   5.0],  alp_v, a_v, lam_v)
        K3_v = find_K([5.0,  100.0], alp_v, a_v, lam_v)
        if any(np.isnan([K1_v,K2_v,K3_v])): return np.nan
        return koide_Q(K1_v, K2_v, K3_v) - 1.5
    try:
        return brentq(Q_minus_Koide, lam_guess*0.99, lam_guess*1.02, xtol=tol)
    except: return np.nan

print(f"\n  Szukanie lambda_Koide (Q=3/2) na siatce (a, alpha):")
print(f"  {'a':>8}  {'alpha':>7}  {'lambda_K':>16}  {'lambda_Koide':>16}  {'dif (ppm)':>12}  {'K3_K':>12}  {'K3_Koide':>12}  {'r31_Koide':>12}")
a_grid   = np.array([0.035, 0.040049, 0.046])
alp_grid = np.array([7.0, 8.5612, 10.0])

for a_v in a_grid:
    for alp_v in alp_grid:
        lam_K_loc = find_lambda_Koide(a_v, alp_v)
        if np.isnan(lam_K_loc): continue
        K1_K = find_K([0.001,  0.1],  alp_v, a_v, lam_K_loc)
        K3_K = find_K([5.0,  100.0],  alp_v, a_v, lam_K_loc)
        K3_lam0 = find_K([5.0, 100.0], alp_v, a_v, LAM0)
        r31_K = K3_K/K1_K if not np.isnan(K1_K) else np.nan
        dif = (lam_K_loc/LAM0 - 1)*1e6
        print(f"  {a_v:8.5f}  {alp_v:7.4f}  {LAM0:.8e}  {lam_K_loc:.8e}  {dif:+12.2f}  {K3_lam0:.8f}  {K3_K:.8f}  {r31_K:.4f}")

# =========================================================
# SEKCJA 7: Podsumowanie — zwiazek K3_Ei, K3*, lambda_K, lambda*
# =========================================================
print(f"\n{'='*72}")
print("SEKCJA 7: Podsumowanie zwiazku K3 - lambda - Q")
print(f"{'='*72}")

I4_v = I4_Eip
I6_v = I6_Ei_v
I6_n = I6_num_full

lam_EiEi = 3*I4_v  / (2*K3_lK**2 * I6_v)
lam_Einum = 3*I4_v / (2*K3_lK**2 * I6_n)
lam_numEi = 3*I4_num_full / (2*K3_lK**2 * I6_v)
lam_numnum = 3*I4_num_full / (2*K3_lK**2 * I6_n)

print(f"""
  PEŁNE ZESTAWIENIE @ a={A0}, alpha={ALP0}:
  ============================================

  K3_num (numeryczne, N=10000): {K3_lK:.12f}
  K3_Ei  (wzor anal. I4_Ei/I6_Ei): {K3_Ei_an:.12f}  ({(K3_Ei_an/K3_lK-1)*1e6:+.2f} ppm)
  K3*    (numeryczne, lambda*):    {K3_ls:.12f}  ({(K3_ls/K3_lK-1)*1e6:+.2f} ppm)
  K3_Ei  - K3* = {(K3_Ei_an - K3_ls):.8f}  = {(K3_Ei_an/K3_ls - 1)*1e6:+.4f} ppm  (blad residualny)

  lambda_K = {LAM0:.10e}
  lambda*  = {LAM_STAR:.10e}  ({(LAM_STAR/LAM0-1)*1e6:+.2f} ppm)
  lambda(I4_Ei, I6_Ei) = {lam_EiEi:.10e}  ({(lam_EiEi/LAM0-1)*1e6:+.2f} ppm od lambda_K)

  I6 blad analityczny: I6_Ei/I6_num - 1 = {(I6_v/I6_n-1)*1e6:+.2f} ppm  (wzor scisly = 0 ppm!)
  => delta K3 od I6: -(1/2)*{(I6_v/I6_n-1)*1e6:+.2f} ppm = {-(I6_v/I6_n-1)/2*1e6:+.2f} ppm
  => obserwowany blad K3: {(K3_Ei_an/K3_lK-1)*1e6:+.2f} ppm  (= K3_Ei/K3_num-1)

  MECHANIZM:
  K3 ~ sqrt(3*I4/(2*lam*I6))
  I6_analyt rozni sie od I6_numeric o {(I6_v/I6_n-1)*1e6:+.2f} ppm  (scisly wzor!)
  Blad K3_Ei vs K3_num NIE pochodzi z I4 ani I6 (oba idealne),
  lecz z wyzszych skladnikow energii TGP poza przyblizeniem lambda-dominujacym:
  => blad K3_Ei/K3_num-1 = {(K3_Ei_an/K3_lK-1)*1e6:+.2f} ppm  (korekta kinetyczna)
  => Lambda_K vs lambda*: (K3_num/K3*)^2 = {(K3_lK/K3_ls)**2:.8f}
     lambda*/lambda_K = {LAM_STAR/LAM0:.8f}
     Potwierdzenie K3~lam^(-1/2): {abs((K3_lK/K3_ls)**2 - LAM_STAR/LAM0)*1e6:.4f} ppm roznica
""")

print(f"\n  WNIOSEK OP-4/OP-12:")
print(f"  K3_Ei ≈ K3* (roznica = {(K3_Ei_an/K3_ls-1)*1e6:+.4f} ppm)")
print(f"  K3_Ei < K3_num (blad {(K3_Ei_an/K3_lK-1)*1e6:+.0f} ppm) z powodu wyzszych skladnikow energii TGP")
print(f"  Q=3/2 dokladnie wymaga lambda* = lambda_K * (K3_num/K3_Ei)^2")
print(f"  = lambda_K * ({K3_lK:.6f}/{K3_Ei_an:.6f})^2 = lambda_K * {(K3_lK/K3_Ei_an)**2:.8f}")
print(f"  = lambda_K * (1 + {((K3_lK/K3_Ei_an)**2 - 1)*1e6:.2f} ppm)")
print(f"  Porowanie z lambda*/lambda_K-1 = {(LAM_STAR/LAM0-1)*1e6:.2f} ppm")

print("\n" + "="*72)
print("P71 ZAKONCZONY")
print("="*72)
