# -*- coding: utf-8 -*-
"""
P57: DIAGNOZA ROZBIEZNOSCI LAMBDA (OP-11)

Dwie metody wyznaczania lambda przy tym samym punkcie (a_Gamma, alpha_K):

  Metoda P50: lambda_analyt takie, ze Q_NNLO(K1_Pade, K2_Pade, K3_Ei, lambda) = 3/2
  Metoda P56: lambda_r31   takie, ze r31(K1_num, K3_Ei(lambda)) = r31_PDG = 3477.2

  Wyniki:
    lambda_analyt / lambda_K = 1.000554  (+554 ppm)
    lambda_r31    / lambda_K = 1.000158  (+158 ppm)
    Roznica: Delta_lambda = 396 ppm

CEL: nie rozwiazac, ale ZDIAGNOZOWAC skad pochodzi 396 ppm.

Hipotezy:
  H1: Blad K1_Pade vs K1_num  (P51: -0.025% => ~250 ppm w Q => ???)
  H2: Blad K2_Pade vs K2_num  (P48: -0.20%  => ~2000 ppm w Q => ???)
  H3: Blad K3_Ei vs K3_num    (K3_Ei=34.195 vs K3_num=34.154, +0.12%)
  H4: Inna definicja lambda_analyt (P50 wyznaczal lambda z warunku Q=3/2
      w przyblizonej formule, P56 z warunku r31 w pelnej numeryce)

Plan:
  A: Wylicz Q przy K1_num, K2_num, K3_num (pelna numeryka)
  B: Wylicz Q przy K1_Pade, K2_Pade, K3_Ei (lancuch analityczny)
  C: Zanalizuj wklad kazdego Ki do rozbnicy Q
  D: Wyznacz lambda_r31 i lambda_analyt przy tym samym (a, alpha, K_i)
  E: Rozloz roznice 396 ppm na skladniki
  F: Sformuluj OP-11 precyzyjnie jako problem
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expi

# -------------------------------------------------------------------
# PARAMETRY
# -------------------------------------------------------------------
A_GAM    = 0.040
ALPHA_K  = 8.5445
LAM_K    = 5.4677e-6
R21_PDG  = 206.770
R31_PDG  = 3477.2
R_MAX    = 50.0
N_GRID   = 3000

# Wartosci numeryczne (z P40/P51/P55)
K1_NUM   = 0.009839    # P40, pelna numeryka
K2_NUM   = 2.0344      # P40, pelna numeryka
K3_NUM   = 34.154      # P40, pelna numeryka

# Wartosci analityczne (lancuch P47-P51)
K1_PADE  = 0.009837    # P51, Pade [2/1], blad -0.025%
K2_PADE  = 2.0304      # P48, Pade [2/2], blad -0.20%
# K3_Ei bedzie obliczone ponizej

def E1(x):
    return -expi(-x)

def V_mod(phi, lam=LAM_K):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha, a, lam=LAM_K, N=N_GRID):
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r)*(-r-1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1.0+alpha/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(phi,lam)-V1)*r**2, r)
    return Ek + Ep

def K3_ei(a=A_GAM, lam=LAM_K):
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=200)
    return np.sqrt(3*I4/(2*lam*I6)) if I6>0 else np.nan

def lambda_from_K3(a, K3):
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=200)
    return 3*I4/(2*I6*K3**2) if (I6>0 and K3>0) else np.nan

def Q_koide(K1, K2, K3):
    s = np.sqrt(K1)+np.sqrt(K2)+np.sqrt(K3)
    return s**2/(K1+K2+K3)

def find_K1_num(alpha=ALPHA_K, a=A_GAM, lam=LAM_K):
    def g(K):
        return energy_num(K, alpha, a, lam)/(4*np.pi*K) - 1.0
    try:
        return brentq(g, 1e-4, 0.4, xtol=1e-10)
    except Exception:
        return np.nan

def find_K2_num(alpha=ALPHA_K, a=A_GAM, lam=LAM_K):
    def g(K):
        return energy_num(K, alpha, a, lam)/(4*np.pi*K) - 1.0
    try:
        return brentq(g, 0.5, 5.0, xtol=1e-10)
    except Exception:
        return np.nan

# -------------------------------------------------------------------
print("="*65)
print("P57: DIAGNOZA ROZBIEZNOSCI LAMBDA (OP-11)")
print("="*65)
print(f"  a_Gamma = {A_GAM},  alpha_K = {ALPHA_K},  lambda_K = {LAM_K:.4e}")
print(f"\n  Rozbieznosc do zbadania:")
print(f"    lambda_analyt (P50) / lambda_K = 1.000554  (+554 ppm)")
print(f"    lambda_r31    (P56) / lambda_K = 1.000158  (+158 ppm)")
print(f"    Roznica: Delta = 396 ppm")

# -------------------------------------------------------------------
print("\n"+"="*65)
print("SEKCJA A: PELNA NUMERYKA -- K_i, Q, lambda")
print("="*65)

K3_ei_val = K3_ei(A_GAM, LAM_K)
K1_n = find_K1_num()
K2_n = find_K2_num()
K3_n = K3_ei_val  # uzywamy K3_Ei jako proxy (K3_num=34.154 z P40)

Q_num = Q_koide(K1_NUM, K2_NUM, K3_NUM)
Q_num_Ei = Q_koide(K1_NUM, K2_NUM, K3_ei_val)

print(f"\n  Wartosci K_i (pelna numeryka, P40):")
print(f"    K1_num = {K1_NUM:.6f}   K2_num = {K2_NUM:.6f}   K3_num = {K3_NUM:.4f}")
print(f"    r21    = {K2_NUM/K1_NUM:.4f}  r31 = {K3_NUM/K1_NUM:.2f}")
print(f"    Q(K1_num, K2_num, K3_num) = {Q_num:.8f}  ({(Q_num-1.5)*1e6:+.1f} ppm)")

print(f"\n  K3_Ei obliczone:")
print(f"    K3_Ei = {K3_ei_val:.6f}  (K3_num = {K3_NUM:.6f}, blad = {(K3_ei_val-K3_NUM)/K3_NUM*100:+.4f}%)")
print(f"    Q(K1_num, K2_num, K3_Ei) = {Q_num_Ei:.8f}  ({(Q_num_Ei-1.5)*1e6:+.1f} ppm)")

# lambda z warunku r31 = r31_PDG przy K1_num
K3_for_r31 = R31_PDG * K1_NUM
lam_r31_num = lambda_from_K3(A_GAM, K3_for_r31)
print(f"\n  Metoda P56 (lambda z r31=r31_PDG przy K1_num):")
print(f"    K3_needed = r31_PDG * K1_num = {K3_for_r31:.6f}")
print(f"    lambda_r31 = {lam_r31_num:.8e}  = lambda_K * {lam_r31_num/LAM_K:.8f}")
print(f"    delta_lambda_r31 = {(lam_r31_num/LAM_K-1)*1e6:+.1f} ppm")

# -------------------------------------------------------------------
print("\n"+"="*65)
print("SEKCJA B: LANCUCH ANALITYCZNY -- K_i_Pade, Q_analyt, lambda_analyt")
print("="*65)

# K3_Ei przy lambda_K
K3_EI = K3_ei_val

Q_analyt_P51 = Q_koide(K1_PADE, K2_PADE, K3_EI)
print(f"\n  Wartosci K_i (analityczne, P47-P51):")
print(f"    K1_Pade = {K1_PADE:.6f}  (blad vs num: {(K1_PADE-K1_NUM)/K1_NUM*100:+.4f}%)")
print(f"    K2_Pade = {K2_PADE:.6f}  (blad vs num: {(K2_PADE-K2_NUM)/K2_NUM*100:+.4f}%)")
print(f"    K3_Ei   = {K3_EI:.6f}  (blad vs num: {(K3_EI-K3_NUM)/K3_NUM*100:+.4f}%)")
print(f"    Q_analyt(K1_Pade, K2_Pade, K3_Ei) = {Q_analyt_P51:.8f}  ({(Q_analyt_P51-1.5)*1e6:+.1f} ppm)")

# Metoda P50: wyznacz lambda takie, ze Q = 3/2 przez zmiane K3
# K3 ~ 1/sqrt(lambda), wiec Q(lambda) = Q_koide(K1_Pade, K2_Pade, K3_Ei(lambda))
# Szukamy lambda_analyt: Q(lambda_analyt) = 3/2

def Q_analyt_of_lam(lam):
    K3 = K3_ei(A_GAM, lam)
    return Q_koide(K1_PADE, K2_PADE, K3) - 1.5

lam_lo = LAM_K * 0.99
lam_hi = LAM_K * 1.01
try:
    lam_analyt = brentq(Q_analyt_of_lam, lam_lo, lam_hi, xtol=1e-14)
    K3_at_lam_analyt = K3_ei(A_GAM, lam_analyt)
    Q_check = Q_koide(K1_PADE, K2_PADE, K3_at_lam_analyt)
    print(f"\n  Metoda P50 (lambda z Q=3/2 w lancuchu analitycznym):")
    print(f"    lambda_analyt = {lam_analyt:.8e}  = lambda_K * {lam_analyt/LAM_K:.8f}")
    print(f"    delta_lambda_analyt = {(lam_analyt/LAM_K-1)*1e6:+.1f} ppm")
    print(f"    K3 przy lambda_analyt = {K3_at_lam_analyt:.6f}")
    print(f"    Q(K1_Pade, K2_Pade, K3@lam_analyt) = {Q_check:.10f}")
except Exception as e:
    print(f"  Blad: {e}")
    lam_analyt = LAM_K * 1.000554
    K3_at_lam_analyt = K3_ei(A_GAM, lam_analyt)
    print(f"  Uzywam P50 wartosci: lambda_analyt = lambda_K * 1.000554")

# -------------------------------------------------------------------
print("\n"+"="*65)
print("SEKCJA C: ROZKLAD ROZBIEZNOSCI NA SKLADNIKI")
print("="*65)

# Definicje:
# lambda_r31   = wyznaczone z r31_PDG * K1_num => K3_needed => lambda
# lambda_analyt = wyznaczone z Q=3/2 w lancuchu Pade

Delta_total = (lam_analyt - lam_r31_num)/LAM_K * 1e6  # w ppm
print(f"\n  lambda_analyt - lambda_r31 = {Delta_total:.1f} ppm")
print(f"  (lambda_analyt = {lam_analyt:.6e}, lambda_r31 = {lam_r31_num:.6e})")

# Hipoteza H1: roznica pochodzi z bledu K1_Pade vs K1_num
# Przy K1_Pade zamiast K1_num: K3_needed = r31_PDG * K1_Pade
K3_for_r31_pade = R31_PDG * K1_PADE
lam_r31_pade = lambda_from_K3(A_GAM, K3_for_r31_pade)
print(f"\n  H1: Blad K1 (Pade vs num)")
print(f"    K1_Pade/K1_num - 1 = {(K1_PADE/K1_NUM-1)*1e6:.1f} ppm  ({(K1_PADE/K1_NUM-1)*100:.4f}%)")
print(f"    lambda_r31(K1_Pade) = {lam_r31_pade:.8e}  (+{(lam_r31_pade/LAM_K-1)*1e6:.1f} ppm)")
print(f"    => wklad bledu K1 w Delta_lambda: {(lam_r31_pade-lam_r31_num)/LAM_K*1e6:.1f} ppm")

# Hipoteza H2: roznica pochodzi z bledu K2_Pade vs K2_num
# lambda_analyt jest wyznaczone z Q(K1_Pade, K2_Pade, K3) = 3/2
# lambda_r31    jest wyznaczone z r31 = K3/K1_num = r31_PDG (nie zalezy od K2!)
# => K2 NIE wchodzi do lambda_r31 => blad K2 w CALOŚCI idzie do lambda_analyt

# Oblicz: ile lambda trzeba zmienic, zeby skorygowac blad K2?
# Zalozmy K1=K1_num, K2=K2_num zamiast K2_Pade w metodzie P50
def Q_corr_of_lam(lam):
    K3 = K3_ei(A_GAM, lam)
    return Q_koide(K1_NUM, K2_NUM, K3) - 1.5

try:
    lam_corr = brentq(Q_corr_of_lam, LAM_K*0.99, LAM_K*1.01, xtol=1e-14)
    print(f"\n  H2: Blad K2 (Pade vs num) w metodzie Q=3/2")
    print(f"    K2_Pade/K2_num - 1 = {(K2_PADE/K2_NUM-1)*1e6:.1f} ppm  ({(K2_PADE/K2_NUM-1)*100:.4f}%)")
    print(f"    lambda(K1_num, K2_num) = {lam_corr:.8e}  (+{(lam_corr/LAM_K-1)*1e6:.1f} ppm)")
    print(f"    lambda(K1_Pade, K2_Pade) = {lam_analyt:.8e}  (+{(lam_analyt/LAM_K-1)*1e6:.1f} ppm)")
    print(f"    => wklad bledu K2 w Delta_lambda: {(lam_analyt-lam_corr)/LAM_K*1e6:.1f} ppm")
    print(f"    => wklad bledu K1 w metodzie P50: {(lam_corr-lam_r31_num)/LAM_K*1e6:.1f} ppm")
except Exception as e:
    print(f"  Blad bisekacji H2: {e}")

# Hipoteza H3: blad K3_Ei vs K3_num
print(f"\n  H3: Blad K3_Ei vs K3_num")
print(f"    K3_Ei/K3_num - 1 = {(K3_ei_val/K3_NUM-1)*1e6:.1f} ppm  ({(K3_ei_val/K3_NUM-1)*100:.4f}%)")
# K3 ~ 1/sqrt(lambda) => delta_lambda/lambda = -2 * delta_K3/K3
delta_lam_from_K3 = -2 * (K3_ei_val/K3_NUM - 1) * 1e6
print(f"    => wklad K3_Ei bledu w lambda (przez K3~1/sqrt(lam)): {delta_lam_from_K3:+.1f} ppm")
print(f"    [Interpretacja: K3_Ei > K3_num => lam_Ei < lam_num]")

# Hipoteza H4: definicyjnie rozne warunki
print(f"\n  H4: Definicyjnie rozne warunki wyznaczania lambda:")
print(f"    P50: Q(K1_Pade, K2_Pade, K3_Ei(lambda)) = 3/2")
print(f"         => lambda wyznaczone przez CALY iloraz (K1,K2,K3)")
print(f"    P56: K3_Ei(lambda) = r31_PDG * K1_num")
print(f"         => lambda wyznaczone TYLKO przez K3 (niezalezne od K2!)")
print(f"    => Roznica = efekt bledu K2 w metodzie P50")

# -------------------------------------------------------------------
print("\n"+"="*65)
print("SEKCJA D: TABELA PODSUMOWUJACA ROZBIEZNOSC")
print("="*65)
print(f"""
  Metoda         Warunek                      lambda*           delta [ppm]
  -----------------------------------------------------------------------
  P50 (Pade)     Q(K1_P,K2_P,K3_Ei)=3/2      {lam_analyt:.6e}     {(lam_analyt/LAM_K-1)*1e6:+.1f}
  P56 (r31)      K3_Ei/K1_num = 3477.2        {lam_r31_num:.6e}     {(lam_r31_num/LAM_K-1)*1e6:+.1f}
  P50+K_num      Q(K1_n,K2_n,K3_Ei)=3/2       {lam_corr:.6e}     {(lam_corr/LAM_K-1)*1e6:+.1f}
  lambda_K       [punkt wyjsciowy]             {LAM_K:.6e}      0.0

  Roznica P50 - P56 = {(lam_analyt-lam_r31_num)/LAM_K*1e6:.1f} ppm

  Rozlozenie wkladu:
    Blad K2_Pade (-0.20%) w formule Q => {(lam_analyt-lam_corr)/LAM_K*1e6:+.1f} ppm
    Blad K1_Pade (-0.025%) w formule Q vs r31 => {(lam_corr-lam_r31_num)/LAM_K*1e6:+.1f} ppm
    Suma = {(lam_analyt-lam_corr+lam_corr-lam_r31_num)/LAM_K*1e6:+.1f} ppm  (vs {(lam_analyt-lam_r31_num)/LAM_K*1e6:.1f} ppm)
""")

# -------------------------------------------------------------------
print("\n"+"="*65)
print("SEKCJA E: PRECYZYJNE SFORMULOWANIE OP-11")
print("="*65)
print(f"""
  OP-11: WEWNETRZNA SPOISTNOSC LANCUCHA ANALITYCZNEGO

  Dwie niezalezne metody wyznaczania lambda daja rozne wartosci:

    lambda_analyt (P50): Q(K1_Pade, K2_Pade, K3_Ei(lambda)) = 3/2
                         => lambda_analyt = lambda_K * (1 + 554 ppm)
    lambda_r31 (P56):    K3_Ei(lambda) / K1_num = r31_PDG
                         => lambda_r31 = lambda_K * (1 + 158 ppm)
    Roznica: 396 ppm

  Glowna przyczyna: blad K2_Pade = -0.20% (vs K2_num).
    Efekt w Q: DQ ≈ dQ/dK2 * DK2 ~ {(Q_koide(K1_PADE,K2_NUM,K3_ei_val)-Q_koide(K1_PADE,K2_PADE,K3_ei_val))*1e6:.0f} ppm
  Blad K1_Pade = -0.025% -> efekt ~ {(lam_corr-lam_r31_num)/LAM_K*1e6:.0f} ppm

  Warunki OP-11 (zebranie w jeden problem):
    (a) Podanie dokladnego K2 z pierwszych zasad (OP-1) automatycznie
        rozwiaze OP-11 -- lambda_analyt zrownalaby sie z lambda_r31.
    (b) Do tego czasu: lambda_r31 (P56, +158 ppm) jest lepsza estymata
        bo nie uzywa K2_Pade (dominujace zrodlo bledu).
    (c) "Prawdziwa" lambda_K (wejsciowa) lezy miedzy:
        lambda_r31 (+158 ppm) a lambda_analyt (+554 ppm).
        Roznica 396 ppm jest miarą bledu K2_Pade w formule Q.

  WNIOSEK: OP-11 jest pochodną OP-1 (K2 z pierwszych zasad).
  Rozwiazanie OP-1 automatycznie rozwiaze OP-11.
""")

print("="*65)
print("P57 ZAKONCZONY")
print("="*65)
