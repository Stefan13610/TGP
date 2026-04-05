# -*- coding: utf-8 -*-
"""
P60: DIAGNOZA I KOREKCJA BLEDU K3_Ei (OP-12)

P58 wykryl: K3_Ei = 34.154, K3_num = 34.249 --> blad -0.28% = -2791 ppm.

Cel:
  1. Rozlozyc g(K) na skladniki: kinetic + V_cubic + V_quartic + V_lambda
     i sprawdzic, ktory z nich jest zrodlem bledu przyblizone K3_Ei.
  2. K3_Ei pochodzi z bilansu: -K^4/4 * I4 + lambda*K^6/6 * I6 = 0
     Pomijane sa: kinetyczny K^2, kubiczny K^3, korekcja (1+K/r>>1).
  3. Wyznaczyc korekcje perturbacyjna: K3 = K3_Ei * (1 + delta)
     i sprawdzic dokladnosc poprawionego wzoru.
  4. Porownac z K3_num z pelnej numeryki.

Dodatkowe: zakres waznosci przybliizenia phi ~ K*exp(-r)/r >> 1.
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq

# -------------------------------------------------------------------
R_MAX  = 50.0
N_GRID = 10000   # pelna precyzja

def E1(x):
    return -expi(-x)

def V_mod(phi, lam):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

# Rozklad energii na skladniki
def energy_components(K, alpha, a, lam, N=N_GRID):
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    dr   = np.diff(r)
    r_m  = 0.5*(r[:-1]+r[1:])   # midpoints dla dr

    phi_full = 1.0 + K*np.exp(-r)/r
    phi_safe = np.maximum(phi_full, 1e-10)
    dphi     = K * np.exp(-r)*(-r-1.0)/r**2

    # Kinetic energy - split alpha=0 part and alpha part
    kin_base = 4*np.pi*np.trapezoid(0.5*dphi**2 * r**2, r)
    kin_alph = 4*np.pi*np.trapezoid(0.5*dphi**2 * alpha/phi_safe * r**2, r)

    # Potential: V(phi)-V(1) = [(phi^3-1)/3 - (phi^4-1)/4 + lam*(phi-1)^6/6]
    V1 = V_mod(1.0, lam)
    dV_cubic  = 4*np.pi*np.trapezoid(((phi_safe**3 - 1.0)/3.0)*r**2, r)
    dV_quart  = 4*np.pi*np.trapezoid((-(phi_safe**4 - 1.0)/4.0)*r**2, r)
    dV_lambda = 4*np.pi*np.trapezoid((lam*(phi_safe-1.0)**6/6.0)*r**2, r)

    E_total = kin_base + kin_alph + dV_cubic + dV_quart + dV_lambda
    g_val   = E_total/(4*np.pi*K) - 1.0

    return dict(
        kin_base=kin_base, kin_alph=kin_alph,
        dV_cubic=dV_cubic, dV_quart=dV_quart, dV_lambda=dV_lambda,
        E_total=E_total, g_val=g_val
    )

def I4_analytic(a):
    """I4 = int_a^inf e^{-4r}/r^2 dr = e^{-4a}/a - 4*E1(4a)"""
    return np.exp(-4*a)/a - 4*E1(4*a)

def I6_numerical(a):
    """I6 = int_a^inf e^{-6r}/r^4 dr  (numeryczny)"""
    val, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=500)
    return val

def K3_ei(a, lam):
    I4 = I4_analytic(a)
    I6 = I6_numerical(a)
    return np.sqrt(3*I4/(2*lam*I6)) if I6 > 0 else np.nan

def find_K3_num(alpha, a, lam, tol=1e-12):
    def g(K, N=N_GRID):
        ec = energy_components(K, alpha, a, lam, N)
        return ec['g_val']
    try:
        return brentq(g, 5.0, 80.0, xtol=tol)
    except:
        return np.nan

# -------------------------------------------------------------------
# Parametry punktu A (nominalnego)
# -------------------------------------------------------------------
A_GAM  = 0.040
ALPHA  = 8.5445
LAM_K  = 5.4677e-6

M_E    = 0.51099895
M_MU   = 105.6583755
M_TAU  = 1776.86
R21    = M_MU/M_E
R31    = M_TAU/M_E

print("="*70)
print("P60: DIAGNOZA BLEDU K3_Ei (OP-12)")
print("="*70)
print(f"  Punkt: a={A_GAM}, alpha={ALPHA}, lambda={LAM_K:.5e}")
print(f"  PDG: r21={R21:.6f}, r31={R31:.6f}")

# -------------------------------------------------------------------
# SEKCJA A: K3_Ei vs K3_num
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA A: K3_Ei vs K3_num")
print(f"{'='*70}")

K3e  = K3_ei(A_GAM, LAM_K)
I4v  = I4_analytic(A_GAM)
I6v  = I6_numerical(A_GAM)
print(f"  I4        = {I4v:.12f}")
print(f"  I6        = {I6v:.12e}")
print(f"  K3_Ei     = {K3e:.12f}")
print(f"  Szukam K3_num (moze chwile potrwac)...")

K3n  = find_K3_num(ALPHA, A_GAM, LAM_K)
blad = (K3e/K3n - 1)*1e6
print(f"  K3_num    = {K3n:.12f}")
print(f"  blad K3_Ei = {blad:+.1f} ppm  ({(K3e/K3n-1)*100:.4f}%)")

# -------------------------------------------------------------------
# SEKCJA B: Rozklad g(K) na skladniki przy K=K3_num
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA B: Rozklad energii g(K) przy K = K3_num i K3_Ei")
print(f"{'='*70}")

for Kval, label in [(K3n, "K3_num"), (K3e, "K3_Ei")]:
    ec = energy_components(Kval, ALPHA, A_GAM, LAM_K)
    E  = ec['E_total']
    print(f"\n  K = {Kval:.6f}  [{label}]")
    print(f"  E_kin_base  = {ec['kin_base']:+.6e}  ({ec['kin_base']/E*100:+.3f}%)")
    print(f"  E_kin_alpha = {ec['kin_alph']:+.6e}  ({ec['kin_alph']/E*100:+.3f}%)")
    print(f"  E_V_cubic   = {ec['dV_cubic']:+.6e}  ({ec['dV_cubic']/E*100:+.3f}%)")
    print(f"  E_V_quart   = {ec['dV_quart']:+.6e}  ({ec['dV_quart']/E*100:+.3f}%)")
    print(f"  E_V_lambda  = {ec['dV_lambda']:+.6e}  ({ec['dV_lambda']/E*100:+.3f}%)")
    print(f"  E_total     = {E:+.6e}")
    print(f"  g(K)        = {ec['g_val']:+.3e}")

# -------------------------------------------------------------------
# SEKCJA C: Analityczne aproksymacje skladnikow przy duzy K
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA C: Aproksymacje skladnikow energii przy duzym K")
print(f"{'='*70}")

# Dla duzego K: phi ~ K*exp(-r)/r + 1
# - dV_lambda ~ lam*K^6/6 * I6  (dominujacy)
# - dV_quart  ~ -K^4/4 * I4     (balansujacy)
# - reszta: korekcje

K = K3n
print(f"\n  Przy K = K3_num = {K:.6f}:")
approx_lam  = LAM_K*K**6/6 * I6v * 4*np.pi
approx_qrt  = -K**4/4       * I4v * 4*np.pi
approx_sum  = approx_lam + approx_qrt

ec = energy_components(K, ALPHA, A_GAM, LAM_K)
print(f"  lam*K^6/6*I6 (approx) = {approx_lam:+.6e}")
print(f"  -K^4/4*I4    (approx) = {approx_qrt:+.6e}")
print(f"  Suma approx           = {approx_sum:+.6e}")
print(f"  E_V_lambda  (pelny)   = {ec['dV_lambda']:+.6e}")
print(f"  E_V_quart   (pelny)   = {ec['dV_quart']:+.6e}")
print(f"  blad lam-approx = {(approx_lam/ec['dV_lambda']-1)*1e6:+.1f} ppm")
print(f"  blad qrt-approx = {(approx_qrt/ec['dV_quart']-1)*1e6:+.1f} ppm")

# Wyznaczyc I4_actual, I6_actual z pelnych calek
# phi = 1 + K*exp(-r)/r, wiec:
# (phi-1)^6 = K^6*exp(-6r)/r^6  -->  dV_lambda = 4pi*lam*K^6/6 * int e^{-6r}/r^4 dr
# (phi^4 - 1)/4 = [phi^4-1]/4 -- pelny wyraz, nie tylko K^4*I4

t  = np.linspace(0, 1, N_GRID)
r  = A_GAM * (R_MAX/A_GAM)**t
ph = 1.0 + K*np.exp(-r)/r
I6_exact = np.trapezoid((ph-1.0)**6 / (LAM_K*K**6) * (6.0/LAM_K/K**6) , r)
# Pokaz roznice miedzy int (phi-1)^6 r^2 dr a K^6*int e^{-6r}/r^4 dr
int_phi6  = np.trapezoid((ph-1.0)**6 * r**2, r)
int_Ke6   = K**6 * np.trapezoid(np.exp(-6*r)/r**4, r)
print(f"\n  int (phi-1)^6 r^2 dr       = {int_phi6:.10e}")
print(f"  K^6 * int e^{{-6r}}/r^4 dr  = {int_Ke6:.10e}")
print(f"  blad aproksymacji (phi-1)^6 ~ K^6*e^{{-6r}}/r^6:")
print(f"     = {(int_phi6/int_Ke6 - 1)*1e6:+.3f} ppm")

int_phi4  = np.trapezoid((ph**4 - 1.0)/4.0 * r**2, r)
int_Ke4   = K**4/4 * np.trapezoid(np.exp(-4*r)/r**2, r)
print(f"\n  int (phi^4-1)/4 * r^2 dr   = {int_phi4:.10e}")
print(f"  K^4/4 * I4                  = {int_Ke4:.10e}")
print(f"  blad aproksymacji phi^4 ~ K^4*e^{{-4r}}/r^4:")
print(f"     = {(int_phi4/int_Ke4 - 1)*1e6:+.3f} ppm")

# -------------------------------------------------------------------
# SEKCJA D: Korekcja perturbacyjna K3
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA D: Korekcja perturbacyjna K3_Ei --> K3_corrected")
print(f"{'='*70}")

# g(K) = 0 przy K = K3. Rozwinmy okolo K3_Ei:
# g(K) = g(K3_Ei) + g'(K3_Ei) * (K3-K3_Ei) + O((K-K3_Ei)^2) = 0
# --> delta = K3 - K3_Ei = -g(K3_Ei) / g'(K3_Ei)

def g_scalar(K, N=N_GRID):
    ec = energy_components(K, ALPHA, A_GAM, LAM_K)
    return ec['g_val']

g_at_K3e  = g_scalar(K3e)
dK        = K3e * 1e-4
g_plus    = g_scalar(K3e + dK)
g_minus   = g_scalar(K3e - dK)
g_prime   = (g_plus - g_minus)/(2*dK)

delta     = -g_at_K3e / g_prime
K3_corr   = K3e + delta

print(f"\n  g(K3_Ei)    = {g_at_K3e:+.8f}")
print(f"  g'(K3_Ei)   = {g_prime:+.8f}")
print(f"  delta       = {delta:+.10f}")
print(f"  K3_Ei       = {K3e:.12f}")
print(f"  K3_corrected= {K3_corr:.12f}")
print(f"  K3_num      = {K3n:.12f}")
print(f"  blad K3_Ei:       {(K3e/K3n-1)*1e6:+.1f} ppm")
print(f"  blad K3_corrected:{(K3_corr/K3n-1)*1e6:+.1f} ppm")

delta_rel = delta/K3e
print(f"\n  Wzgledna korekcja: delta/K3_Ei = {delta_rel:+.6f} = {delta_rel*1e6:+.1f} ppm")
print(f"  Glowna przyczyna: g(K3_Ei) != 0 bo formula K3_Ei pomija")
print(f"    1) kinetyczny skladnik (proporcjonalny do K^2, alpha)")
print(f"    2) kubiczny skladnik V_cubic (proporcjonalny do K^3)")
print(f"    3) korekcje do (phi-1)^6 ~ K^6 e^{{-6r}}/r^6 (phi != K/r)")

# -------------------------------------------------------------------
# SEKCJA E: Analityczne korekcje w rownaniu g(K3)=0
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA E: Analityczne oszacowanie korekcji do K3_Ei")
print(f"{'='*70}")

# g(K) = [Ek + Ep]/(4pi*K) - 1 = 0
# Ep = 4pi[-K^4/4*I4 + lam*K^6/6*I6] + korekcje_kvart + korekcje_lam
# Ek = 4pi*K^2*[J2 + alpha*J_alpha]  gdzie J2, J_alpha to calki kinetyczne

# Wyznaczmy J2, J_alpha numerycznie:
t2   = np.linspace(0, 1, N_GRID)
r2   = A_GAM*(R_MAX/A_GAM)**t2
ph2  = np.maximum(1.0 + K3n*np.exp(-r2)/r2, 1e-10)
dphi = K3n*np.exp(-r2)*(-r2-1.0)/r2**2

# Ek = 4pi K^2 * int (exp(-r)(-r-1)/r^2)^2 * (1+alpha/phi) * r^2 dr
#    = 4pi * [K^2*J2 + K^2*alpha*J_al]
integrand_base  = 0.5 * (np.exp(-r2)*(-r2-1.0)/r2**2)**2 * r2**2
integrand_alpha = 0.5 * (np.exp(-r2)*(-r2-1.0)/r2**2)**2 / ph2 * r2**2

J2    = np.trapezoid(integrand_base, r2)
J_al  = np.trapezoid(integrand_alpha, r2)
Ek_approx = 4*np.pi * (K3n**2 * J2 + K3n**2 * ALPHA * J_al)

ec_test = energy_components(K3n, ALPHA, A_GAM, LAM_K)
Ek_full = ec_test['kin_base'] + ec_test['kin_alph']

print(f"\n  J2   = {J2:.8e}  (kinetic base integral)")
print(f"  J_al = {J_al:.8e}  (kinetic alpha integral)")
print(f"  Ek_approx (K^2) = {Ek_approx:.6e}")
print(f"  Ek_full         = {Ek_full:.6e}")
print(f"  blad Ek-approx  = {(Ek_approx/Ek_full-1)*1e6:+.1f} ppm")

# Wzor korekcji na K3:
# g(K) approx = 0:
# [-K^4/4*I4 + lam*K^6/6*I6 + K^2*(J2 + alpha*J_al)] / K = 0
# lam*K^6/6*I6 - K^4/4*I4 + K^2*(J2+alpha*J_al) = K
# Podstaw K3 = K3_Ei + delta, K3_Ei^2 = 3*I4/(2*lam*I6):
# Dominujace: 2*lam*K3_Ei^5*I6*delta/6 - K3_Ei^3*I4*delta + K3_Ei^2*(J2+al*Jal) = K3_Ei
# delta*(lam*K3_Ei^4*I6/3 - K3_Ei^2*I4) + K3_Ei^2*(J2+al*Jal) = K3_Ei
# delta*(K3_Ei^2 * [lam*K3_Ei^2*I6/3 - I4]) + K3_Ei*(J2+al*Jal) - 1 = 0
# Uzywajac lam*K3_Ei^2*I6 = 3*I4/2:
# lam*K3_Ei^2*I6/3 = I4/2
# delta*(K3_Ei^2 * [I4/2 - I4]) + K3_Ei*(J2+al*Jal) - 1 = 0
# delta*(-K3_Ei^2*I4/2) + K3_Ei*(J2+al*Jal) - 1 = 0
# delta = [K3_Ei*(J2+al*Jal) - 1] / (K3_Ei^2*I4/2)

delta_analyt = (K3e*(J2 + ALPHA*J_al) - 1.0) / (K3e**2 * I4v / 2.0)
K3_corr_an   = K3e * (1.0 + delta_analyt)

print(f"\n  Analityczna korekcja (1-stop, kinetic only):")
print(f"  delta_analyt = {delta_analyt:+.8f}")
print(f"  K3_corrected_an = {K3_corr_an:.12f}")
print(f"  blad vs K3_num  = {(K3_corr_an/K3n-1)*1e6:+.1f} ppm")

# Dlaczego numeryczna korekcja (Sekcja D) jest lepsza?
# Bo uwzglednia WSZYSTKIE skladniki energii, nie tylko kinetyczny.
print(f"\n  Porownanie metod korekcji K3:")
print(f"  K3_Ei              = {K3e:.10f}  (blad {(K3e/K3n-1)*1e6:+.0f} ppm)")
print(f"  K3_corr (Newton-1) = {K3_corr:.10f}  (blad {(K3_corr/K3n-1)*1e6:+.1f} ppm)")
print(f"  K3_corr (analyt)   = {K3_corr_an:.10f}  (blad {(K3_corr_an/K3n-1)*1e6:+.1f} ppm)")
print(f"  K3_num  (brentq)   = {K3n:.10f}  (referencja)")

# -------------------------------------------------------------------
# SEKCJA F: Krok Newtona 2 (pelna dokladnosc)
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA F: Iteracja Newtona na g(K)=0 (dwa kroki)")
print(f"{'='*70}")

K_iter = K3e
for step in range(3):
    g_val  = g_scalar(K_iter)
    dK_    = K_iter * 1e-4
    g_p    = (g_scalar(K_iter+dK_) - g_scalar(K_iter-dK_))/(2*dK_)
    K_iter = K_iter - g_val/g_p
    print(f"  Krok {step+1}: K3 = {K_iter:.12f}  (blad vs K3_num: {(K_iter/K3n-1)*1e6:+.2f} ppm)")

# -------------------------------------------------------------------
# SEKCJA G: Wnioski -- czy K3_Ei wymaga poprawy w lancuchu?
# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("SEKCJA G: WNIOSKI - WPLYW BLEDU K3_Ei NA WYNIKI FIZYCZNE")
print(f"{'='*70}")

# K3 wplywa na: r31 = K3/K1, Q (przez sqrt(K3))
K1_B = 0.009833303   # Punkt B, P58
K2_B = 2.033219233

def Q_val(K1, K2, K3):
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    return s**2 / (K1 + K2 + K3)

Q_ei   = Q_val(K1_B, K2_B, K3e)
Q_num  = Q_val(K1_B, K2_B, K3n)
Q_corr = Q_val(K1_B, K2_B, K3_corr)

r31_ei  = K3e  / K1_B
r31_num = K3n  / K1_B
r31_pdg = M_TAU/M_E

print(f"\n  Wplyw bledu K3_Ei na wyniki fizyczne (Punkt B K1,K2):")
print(f"  {'Wariant':20s}  {'K3':>12}  {'r31':>10}  {'Q-3/2 [ppm]':>12}")
print(f"  {'-'*20}  {'-'*12}  {'-'*10}  {'-'*12}")
print(f"  {'K3_Ei (biezacy)':20s}  {K3e:12.6f}  {r31_ei:10.4f}  {(Q_ei-1.5)*1e6:+12.3f}")
print(f"  {'K3_num (pelny)':20s}  {K3n:12.6f}  {r31_num:10.4f}  {(Q_num-1.5)*1e6:+12.3f}")
print(f"  {'K3_corrected':20s}  {K3_corr:12.6f}  {r31_num:10.4f}  {(Q_corr-1.5)*1e6:+12.3f}")
print(f"  {'PDG':20s}  {'—':>12}  {r31_pdg:10.4f}  {'—':>12}")

print(f"\n  Delta r31 = r31_Ei - r31_num = {r31_ei - r31_num:+.4f}")
print(f"  Delta Q   = Q_Ei - Q_num    = {(Q_ei-Q_num)*1e6:+.3f} ppm")
print(f"\n  Blad K3_Ei (-2791 ppm) przekklada sie na:")
print(f"    dr31 = {(r31_ei-r31_num):+.4f} = {(r31_ei/r31_num-1)*1e6:+.1f} ppm w r31")
print(f"    dQ   = {(Q_ei-Q_num)*1e6:+.3f} ppm  (maly, bo Q slabo zalezne od K3)")
print(f"\n  Wniosek: blad K3_Ei dominuje w r31, ale nie w Q!")
print(f"  Q jest mniej czuly na K3 niz r31.")

# -------------------------------------------------------------------
print(f"\n{'='*70}")
print("P60 ZAKONCZONY")
print(f"{'='*70}")
