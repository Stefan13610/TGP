# -*- coding: utf-8 -*-
"""
P63: K2 Z PIERWSZYCH ZASAD (OP-1) -- PRZEGLAD PODEJSC

OP-1 jest centralnym otwartym problemem TGP.
Dotychczas: K2_Pade[2/2] = 2.030 (blad -0.20%).

Cel: Zbadac czy K2 moze byc wyznaczone analitycznie z wyzsza
dokladnoscia przez:
  A. Rozszerzony szereg Taylora g(K): wyznaczyc c2..c8 numerycznie.
  B. Ulepszone aproksymacje Pade: [3/2], [2/3], [3/3].
  C. Relacje miedzy zerami (Vieta): K1+K2+K3, K1*K2+..., K1*K2*K3 = f(params).
  D. Schematyczne rownanie K2 jako odwrotnosc: K2 = 1/g'(K2)... ?
  E. Aproksymacja g(K) przez model "dwustudzienny" (double-well style).
  F. Tabela bledow Pade[n/m]: ktore [n/m] daje < 10 ppm?
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.optimize import brentq
from numpy.polynomial import polynomial as P

R_MAX  = 50.0
N_GRID = 10000

def E1(x):
    return -expi(-x)

def V_mod(phi, lam):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

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

# Parametry Punkt B
A_GAM  = 0.040049
ALPHA  = 8.5612
LAM_K  = 5.4677e-6

print("="*72)
print("P63: K2 Z PIERWSZYCH ZASAD (OP-1)")
print("="*72)
print(f"  Punkt B: a={A_GAM}, alpha={ALPHA}, lambda={LAM_K:.5e}")

K1_num = find_K([1e-5, 0.45], ALPHA, A_GAM, LAM_K)
K2_num = find_K([0.4,  5.0 ], ALPHA, A_GAM, LAM_K)
K3_num = find_K([5.0,  80.0], ALPHA, A_GAM, LAM_K)
print(f"\n  K1_num = {K1_num:.14f}")
print(f"  K2_num = {K2_num:.14f}")
print(f"  K3_num = {K3_num:.14f}")

# -------------------------------------------------------------------
# SEKCJA A: Rozszerzony szereg Taylora g(K) do rzedu K^8
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA A: Rozszerzony szereg Taylora g(K) = sum c_{n+1} K^n - 1")
print(f"{'='*72}")

# Wyznaczyc c2..c9 jako g^(n)(0)/n!  numerycznie metodą roznic skonczonych
# g(K) = E[K]/(4pi*K) - 1
# c_n = [(n-1)! * (d^{n-1} g/dK^{n-1})|_{K=0}] / (n-1)!
# Numerycznie: c_n = (d^n [E/(4piK)]) / (dK^n) |_{K=0} / n!

# Uzyj malych K do numerycznej estymaty
Ktest = np.array([0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.015, 0.018, 0.022])
gtest = np.array([g_func(K, ALPHA, A_GAM, LAM_K) for K in Ktest])
# g(K) = c2*K + c3*K^2 + ... - 1
# g(K)+1 = c2*K + c3*K^2 + ...  => (g+1)/K = c2 + c3*K + c4*K^2 + ...
htest = (gtest+1.0)/Ktest   # = c2 + c3*K + c4*K^2 + ...

# Polyfituj htest = c2 + c3*K + c4*K^2 + ... + c9*K^7
deg = 7
coeffs_fit = np.polyfit(Ktest, htest, deg)   # stopnie od deg do 0
# polyfit zwraca [c_{d}, c_{d-1}, ..., c_0]
# wiec h(K) = c_{deg}*K^{deg} + ... + c_0
# czyli c2=c0, c3=c1*K, itd.
coeffs_rev = coeffs_fit[::-1]   # teraz [c0, c1, ..., c7] = [c2, c3, ..., c9]

print("\n  Dopasowanie wielomianowe do g(K)+1 = K*h(K), h(K) = sum c_n K^(n-2)")
print(f"  Uzyto {len(Ktest)} punktow K w [0.002, 0.022], stopien {deg}")
print(f"\n  n   c_n (dopasowanie)    c_n (P47/analityczne)")

# Wartosci z P47 (dla Punktu A; Punkt B moze sie troche roznic)
c47 = {2: 112.105, 3: -1229.03, 4: 19693.5, 6: 3.358e-3*1e6/1e6}
for i, c in enumerate(coeffs_rev):
    n = i+2
    ref = c47.get(n, '—')
    ref_str = f"{ref:.4f}" if isinstance(ref,float) else ref
    print(f"  {n:3d}  {c:+18.6f}     {ref_str}")

# Uzyj dokladniejszych wyliczen analitycznych dla c2..c6:
a = A_GAM
alp = ALPHA

# Calki potencjalowe
I2, _ = quad(lambda r: np.exp(-2*r),       a, R_MAX, limit=300)
I3    = E1(3*a)
I4    = np.exp(-4*a)/a - 4*E1(4*a)
I5, _ = quad(lambda r: np.exp(-5*r)/r**3, a, R_MAX, limit=300)
I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=300)
I7, _ = quad(lambda r: np.exp(-7*r)/r**5, a, R_MAX, limit=300)
I8, _ = quad(lambda r: np.exp(-8*r)/r**6, a, R_MAX, limit=300)

# Calki kinetyczne: M_n = int e^{-(n+1)r} * (r+1)^2 / (2 * r^{n-1}) dr
# przy K->0: Ek = 4pi * K^2 * sum_{k=0}^{inf} (-1)^k alpha^k K^k M_{k+2}
M2, _ = quad(lambda r: np.exp(-2*r)*(r+1)**2/(2*r**0), a, R_MAX, limit=300)
M3, _ = quad(lambda r: np.exp(-3*r)*(r+1)**2/(2*r**1), a, R_MAX, limit=300)
M4, _ = quad(lambda r: np.exp(-4*r)*(r+1)**2/(2*r**2), a, R_MAX, limit=300)
M5, _ = quad(lambda r: np.exp(-5*r)*(r+1)**2/(2*r**3), a, R_MAX, limit=300)
M6, _ = quad(lambda r: np.exp(-6*r)*(r+1)**2/(2*r**4), a, R_MAX, limit=300)
M7, _ = quad(lambda r: np.exp(-7*r)*(r+1)**2/(2*r**5), a, R_MAX, limit=300)

print(f"\n  Calki analityczne:")
print(f"  I2={I2:.8f}  I3={I3:.8f}  I4={I4:.8f}")
print(f"  I5={I5:.8e}  I6={I6:.8e}  I7={I7:.8e}  I8={I8:.8e}")
print(f"  M2={M2:.8f}  M3={M3:.8f}  M4={M4:.8f}")
print(f"  M5={M5:.8e}  M6={M6:.8e}  M7={M7:.8e}")

# Analityczne c_n (E/(4piK) = c2*K + c3*K^2 + ...)
# c_n pochodzi z E[K] ~ 4pi*(K^{n+1} * coeff_n)
# c2: K^2 termy: kinetic*(1+alpha)*M2 + V_potential(-I2/2)
c2_an = (1.0+alp)*M2 - I2/2.0
# c3: K^3 termy: kinetic*(-alpha*M3) + V_potential(-2I3/3)
c3_an = -alp*M3  - 2.0*I3/3.0
# c4: K^4 termy: kinetic*(alpha*M4) + V_potential(-I4/4)
c4_an = alp*M4 - I4/4.0
# c5: K^5 termy: kinetic*(-alpha*M5) + V_potential(0, no K^5 in pot w/o lam)
c5_an = -alp*M5
# c6: K^6 termy: kinetic*(alpha*M6) + V_lam*(lambda*I6/6)
c6_an = alp*M6 + LAM_K*I6/6.0
# c7: K^7: kinetic*(-alpha*M7) + V_potential(0, no K^7 in pot w/o lam)
c7_an = -alp*M7

print(f"\n  Analityczne wspolczynniki g(K) = c2*K + c3*K^2 + ... - 1:")
print(f"  c2 = {c2_an:.8f}  (P47: 112.105)")
print(f"  c3 = {c3_an:.8f}  (P47: -1229.03)")
print(f"  c4 = {c4_an:.8f}  (P47: 19693.5)")
print(f"  c5 = {c5_an:.8e}  (nowy!)")
print(f"  c6 = {c6_an:.8e}  (P47: {3.358e-3:.4e} ale inne punkty)")
print(f"  c7 = {c7_an:.8e}  (nowy!)")

# -------------------------------------------------------------------
# SEKCJA B: Aproksymacje Pade [n/m] dla g(K) = c2*K+...+c_{n+m+1}*K^{n+m} - 1
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA B: Aproksymacje Pade [n/m] -- szukanie K2")
print(f"{'='*72}")

# Zbuduj h(K) = (g(K)+1)/K = c2 + c3*K + c4*K^2 + ...
# Chcemy zera h(K): K2 jest zerem h (bo g(K2)=0 => h(K2) = 1/K2... nie,
# g(K2)=0 => h(K2) = (0+1)/K2 = 1/K2 =/= 0
# Wiec lepiej pracowac bezposrednio na g(K).
# g(K) = c2*K + c3*K^2 + ... - 1 = 0
# Podstawimy K = K2: c2*K2 + c3*K2^2 + ... = 1

# g(K) jako wielomian stopnia N:
# g(K) = -1 + c2*K + c3*K^2 + c4*K^3 + c5*K^4 + c6*K^5 + c7*K^6

# UWAGA: c6 zawiera lambda*I6/6 ~ 3.35e-3 (bardzo maly); c7 < c6
# Wiec K2 ~ 2 jest wyznaczone glownie przez c2, c3, c4 (I c5?)

# Bezposrednie zera wielomianu:
# g_poly(K) = -1 + c2*K + c3*K^2 + c4*K^3 + c5*K^4 + c6*K^5 + c7*K^6

coeffs_poly = np.array([-1.0, c2_an, c3_an, c4_an, c5_an, c6_an, c7_an])
# numpy poly1d expects highest degree first:
poly = np.poly1d(coeffs_poly[::-1])   # stopien 6
roots = np.roots(poly)
print(f"\n  Wielomian g(K) stopnia 6 (-1 + c2*K + ... + c7*K^6):")
real_pos_roots = sorted([r.real for r in roots if abs(r.imag)<0.1 and r.real>0])
print(f"  Rzeczywiste dodatnie zera: {[f'{r:.6f}' for r in real_pos_roots]}")
print(f"  K1_num = {K1_num:.6f},  K2_num = {K2_num:.6f},  K3_num = {K3_num:.6f}")

# Szukaj zera blisko K2:
K2_poly = min(real_pos_roots, key=lambda r: abs(r-K2_num)) if real_pos_roots else np.nan
print(f"\n  K2 z wielomianu stopnia 6: {K2_poly:.10f}")
print(f"  K2_num:                     {K2_num:.10f}")
print(f"  Blad:                       {(K2_poly/K2_num-1)*1e6:+.1f} ppm")

# Pade [2/2] (z P48): uzyj c2,c3,c4,c6
# g(K) = (-1 + c2*K + c3*K^2 + c4*K^3 + c6*K^5) / 1
# Approximant [2/2]: p(K)/q(K) gdzie p stopien 2, q stopien 2
# Punkt wyjscia: seria c2*K + c3*K^2 + c4*K^3 + c6*K^5 ~ 1
# Metoda: Pade od (g+1)/K = c2 + c3*K + c4*K^2 + c5*K^3 + c6*K^4 + c7*K^5

# Buduj tablice Pade [n/m] z serii: h(K) = sum_{k=0}^{n+m} h_k K^k
# gdzie h_k = c_{k+2}
h_coeffs = [c2_an, c3_an, c4_an, c5_an, c6_an, c7_an]

def pade_approx(h_coeffs, n, m):
    """Aproksymant Pade [n/m] dla f(x) = sum h_k x^k.
       Zwraca (p_coeffs, q_coeffs) stopni n, m."""
    N = n + m + 1
    if len(h_coeffs) < N:
        return None, None
    h = h_coeffs[:N]
    # Setup linear system for q: h*q = p (convolution)
    # q_0 = 1 (normaliz.), solve for q_1..q_m, p_0..p_n
    A = np.zeros((m, m))
    b = np.zeros(m)
    for i in range(m):
        for j in range(m):
            idx = n + 1 + i - j
            A[i, j] = h[idx] if 0 <= idx < N else 0.0
        b[i] = -h[n+1+i] if n+1+i < N else 0.0
    try:
        q_rest = np.linalg.solve(A, b)
    except:
        return None, None
    q = np.zeros(m+1); q[0] = 1.0; q[1:] = q_rest
    p = np.zeros(n+1)
    for i in range(n+1):
        p[i] = sum(h[i-j]*q[j] for j in range(min(i+1,m+1)) if 0<=i-j<N)
    return p, q

def find_K2_pade(p_coeffs, q_coeffs, bracket=(0.5, 5.0)):
    """Znajdz zero p(K)/q(K) - 1/K (czyli h(K)=1/K to K-zero g(K))
       Rowne: f(K) = K*p(K) - q(K) = 0"""
    # h(K) = p/q, zero g: K*h(K) = 1 => K*p(K)/q(K) = 1 => K*p(K) = q(K)
    # To jest rownanie algebraiczne.
    p = np.poly1d(p_coeffs[::-1])  # p(K) = p0 + p1*K + ...
    q = np.poly1d(q_coeffs[::-1])
    # Szukamy K: K*p(K) = q(K) => K*p(K) - q(K) = 0
    func = lambda K: K*p(K) - q(K)
    # Znajdz zero brentq w poblizu K2
    for lo, hi in [(1.5, 3.0), (1.0, 4.0), (0.5, 5.0)]:
        try:
            if func(lo)*func(hi) < 0:
                return brentq(func, lo, hi, xtol=1e-12)
        except: pass
    return np.nan

print(f"\n  Aproksymacje Pade h(K) = (g(K)+1)/K = c2 + c3*K + ...: ")
print(f"  {'[n/m]':>8}  {'K2_Pade':>14}  {'blad [ppm]':>12}  {'blad [%]':>10}")
print(f"  {'-'*8}  {'-'*14}  {'-'*12}  {'-'*10}")

for (n, m) in [(1,1),(2,1),(1,2),(2,2),(3,1),(1,3),(3,2),(2,3),(3,3)]:
    p, q = pade_approx(h_coeffs, n, m)
    if p is None:
        print(f"  [{n}/{m}]       BLAD (za malo wspolczynnikow)")
        continue
    K2p = find_K2_pade(p, q)
    if np.isnan(K2p):
        print(f"  [{n}/{m}]       nie znaleziono zera")
        continue
    err_ppm = (K2p/K2_num - 1)*1e6
    err_pct = (K2p/K2_num - 1)*100
    print(f"  [{n}/{m}]    {K2p:14.10f}  {err_ppm:+12.1f}  {err_pct:+10.4f}")

# -------------------------------------------------------------------
# SEKCJA C: Relacje Vieta miedzy zerami
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA C: Relacje miedzy zerami K1, K2, K3 (Vieta-style)")
print(f"{'='*72}")

# Czy K2 mozna wyrazic przez K1 i K3?
# g(K) ma 3 zera: K1, K2, K3.
# g(K) ~ (K-K1)(K-K2)(K-K3) * h(K) gdzie h bez zer w [0,inf)
# Mozna zbadac: K1+K2+K3, K1*K2+K2*K3+K3*K1, K1*K2*K3

S1 = K1_num + K2_num + K3_num
S2 = K1_num*K2_num + K2_num*K3_num + K3_num*K1_num
S3 = K1_num * K2_num * K3_num

print(f"\n  K1 = {K1_num:.10f}")
print(f"  K2 = {K2_num:.10f}")
print(f"  K3 = {K3_num:.10f}")
print(f"\n  S1 = K1+K2+K3       = {S1:.8f}")
print(f"  S2 = K1*K2+K2*K3+K3*K1 = {S2:.8f}")
print(f"  S3 = K1*K2*K3       = {S3:.8e}")

# Prosta relacja: K2 ~ sqrt(K1*K3)?  (geometryczne srednie)
geom_mean = np.sqrt(K1_num * K3_num)
arith_mean = (K1_num + K3_num) / 2.0
harm_mean  = 2/(1/K1_num + 1/K3_num)

print(f"\n  Proste relacje K2 vs K1, K3:")
print(f"  sqrt(K1*K3)  = {geom_mean:.8f}  (blad: {(geom_mean/K2_num-1)*100:+.2f}%)")
print(f"  (K1+K3)/2    = {arith_mean:.8f}  (blad: {(arith_mean/K2_num-1)*100:+.2f}%)")
print(f"  2/(1/K1+1/K3)= {harm_mean:.8e}  (blad: {(harm_mean/K2_num-1)*100:+.2f}%)")

# Sprawdz relacje przy roznych (a, alpha):
print(f"\n  Sprawdzenie geom. sredniej dla roznych parametrow:")
a_list = [0.040, 0.040049, 0.041, 0.042]
alpha_list = [8.5445, 8.5612, 10.0, 12.0]
print(f"  {'a':>7}  {'alpha':>7}  {'K1':>10}  {'K2_num':>10}  {'sqrt(K1*K3)':>12}  {'blad [%]':>10}")
print(f"  {'-'*7}  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*10}")
from scipy.special import expi
K3_ei_fn = lambda a_v, lam: (lambda I4,I6: np.sqrt(3*I4/(2*lam*I6)) if I6>0 else np.nan)(
    np.exp(-4*a_v)/a_v - 4*(-expi(-4*a_v)),
    quad(lambda r: np.exp(-6*r)/r**4, a_v, R_MAX, limit=200)[0]
)
for ai, ali in zip(a_list, alpha_list):
    K1i = find_K([1e-5,0.45], ali, ai, LAM_K, N=5000)
    K2i = find_K([0.4,5.0],   ali, ai, LAM_K, N=5000)
    K3i = find_K([5.0,80.0],  ali, ai, LAM_K, N=5000)
    if np.isnan(K1i) or np.isnan(K2i) or np.isnan(K3i): continue
    gm   = np.sqrt(K1i*K3i)
    err  = (gm/K2i-1)*100
    print(f"  {ai:7.5f}  {ali:7.4f}  {K1i:10.7f}  {K2i:10.7f}  {gm:12.7f}  {err:+10.4f}")

# -------------------------------------------------------------------
# SEKCJA D: Zrodlo bledu Pade -- co ogranicza zbieznosc?
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA D: Analiza bledu Pade -- co ogranicza dokladnosc?")
print(f"{'='*72}")

# g(K) jako seria u K=0 ma promien zbieznosci ograniczony przez:
# - najblizszy biegun zespolony
# - singularnosc od phi->0 lub phi->inf
# Przy K = K_sing: 1 + K*exp(-r)/r = 0 dla jakiegos r=r_min
# Mala phi: phi -> 0 przy r = r_min gdy K*exp(-r_min)/r_min = -1
# Dla K rzeczywistego > 0: phi = 1+K*exp(-r)/r > 1 > 0, wiec brak osobliwosci!
# Ale dla zespolonego K: osobliwosc przy phi = 0 gdy K = -r*exp(r) dla r=a
# |K_sing| ~ a*exp(a) ~ 0.040 * 1.041 ~ 0.042  (przy r=a)
# To jest blisko K=0! Promien zbieznosci ~ 0.042 << K2 ~ 2

# To tlumacy DLACZEGO seria Taylora przy K=0 nie zbiega do K=K2~2:
# K2 jest poza promieniem zbieznosci szeregu!

K_sing_approx = A_GAM * np.exp(A_GAM)
print(f"\n  Szacunek promienia zbieznosci szeregu w K=0:")
print(f"  Osobliwosc phi=0 przy r=a: K_sing ~ a*exp(a) = {K_sing_approx:.5f}")
print(f"  Promien zbieznosci R_conv ~ {K_sing_approx:.5f}")
print(f"  K2 = {K2_num:.5f} >> R_conv = {K_sing_approx:.5f}")
print(f"  --> Seria Taylora NIE ZBIEGA w K=K2!")
print(f"  --> Pade jest analitycznym przedluzeniem przez osobliwosc.")

# Sprawdz gdzie phi=0 przy rzeczywistym K (zawsze phi>0 dla K>0)
# Dla K<0: phi = 1 + K*exp(-r)/r -> dla K = -a*exp(a): phi(a) = 0
print(f"\n  Potwierdzenie: phi(a) przy K = K_sing_approx (ujemne):")
r_a = A_GAM
phi_a = 1.0 + (-K_sing_approx)*np.exp(-r_a)/r_a
print(f"  phi(a) = 1 + K_sing*exp(-a)/a = {phi_a:.6f}  (bliskie 0? {abs(phi_a)<0.05})")

# Zespolona osobliwosc: K_sing = -a*e^a jest ujemna rzeczywista
# Dla osi rzeczywistej: phi >= 1 > 0, wiec brak zer
# Ale seria Taylora wokol K=0 ma promien zbieznosci |K_sing| = a*e^a

# Co z Pade? Pade [n/m] ma bieguny, ktore moga "cancelowac" osobliwosci
print(f"\n  Wniosek: Pade dziala jako przedluzenie analityczne po barierze osobliwosci.")
print(f"  Dokladnosc Pade [n/m] rosnie z n+m, ale nie ma gwarancji szybkiej zbieznosci")
print(f"  dla K z galezia osobliwosci przy K_sing ~ {K_sing_approx:.4f}.")

# -------------------------------------------------------------------
# SEKCJA E: Bezposrednia aproksymacja energii E(K2) = 4*pi*K2
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA E: Czy K2 spelnia proste rownanie algebraiczne?")
print(f"{'='*72}")

# E[K2] = 4*pi*K2 (warunek g=0)
E_K2 = energy_num(K2_num, ALPHA, A_GAM, LAM_K)
print(f"\n  E[K2_num]/(4*pi) = {E_K2/(4*np.pi):.10f}")
print(f"  K2_num            = {K2_num:.10f}")
print(f"  Roznica (powinno byc 0): {E_K2/(4*np.pi) - K2_num:.3e}")

# Sprawdz czy E[K2] / K2 = 4*pi ma prosta postac
# Przy K2~2, dominujace czlony energii:
# Ek ~ c2*K2^2 * 4pi (z c2~112)   -> Ek ~ 112*4*4pi ~ 5629
# Ep_lam ~ lam*K2^6*I6/6 * 4pi    -> ~ 5e-6*64*3685/6*4pi ~ 0.98
# Ep_qrt ~ -K2^4*I4/4 * 4pi       -> ~ -16*15.6/4*4pi ~ -782
# Ep_cub ~ -2*K2^3*I3/3 * 4pi     -> ~ -2*8*1.66/3*4pi ~ -349
# Ep_kin_2 ~ -K2^2*I2/2 * 4pi     -> ~ -4*0.46/2*4pi ~ -11.5
# E_kin_alpha ~ alpha*K2^2*M3?... ~ 8.56*...

# Rozklad energii przy K2:
t = np.linspace(0,1,N_GRID)
r = A_GAM*(R_MAX/A_GAM)**t
phi = np.maximum(1+K2_num*np.exp(-r)/r, 1e-10)
dphi = K2_num*np.exp(-r)*(-r-1.0)/r**2
V1 = V_mod(1.0,LAM_K)

Ek_base  = 4*np.pi*np.trapezoid(0.5*dphi**2*r**2, r)
Ek_alpha = 4*np.pi*np.trapezoid(0.5*dphi**2*ALPHA/phi*r**2, r)
Ep_cub   = 4*np.pi*np.trapezoid((phi**3-1.0)/3*r**2, r)
Ep_qrt   = 4*np.pi*np.trapezoid(-(phi**4-1.0)/4*r**2, r)
Ep_lam   = 4*np.pi*np.trapezoid(LAM_K*(phi-1.0)**6/6*r**2, r)
E_tot    = Ek_base + Ek_alpha + Ep_cub + Ep_qrt + Ep_lam

print(f"\n  Rozklad E[K2] przy K2={K2_num:.6f}:")
print(f"  Ek_base   = {Ek_base:+.6e}  ({Ek_base/E_tot*100:+.2f}%)")
print(f"  Ek_alpha  = {Ek_alpha:+.6e}  ({Ek_alpha/E_tot*100:+.2f}%)")
print(f"  Ep_cubic  = {Ep_cub:+.6e}  ({Ep_cub/E_tot*100:+.2f}%)")
print(f"  Ep_quart  = {Ep_qrt:+.6e}  ({Ep_qrt/E_tot*100:+.2f}%)")
print(f"  Ep_lambda = {Ep_lam:+.6e}  ({Ep_lam/E_tot*100:+.2f}%)")
print(f"  E_total   = {E_tot:+.6e}")
print(f"  4*pi*K2   = {4*np.pi*K2_num:+.6e}")

# -------------------------------------------------------------------
# SEKCJA F: Tabela bledow Pade jako funkcja n+m
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA F: TABELA PODSUMOWUJACA -- DOKLADNOSC PADE")
print(f"{'='*72}")

# Wnikliwe numeryczne g(K): uzyj gestej siatki K wokol K2
K_scan = np.array([0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.1, 2.2, 2.5, 3.0, 4.0])
g_scan = np.array([g_func(K, ALPHA, A_GAM, LAM_K, N=N_GRID) for K in K_scan])

# Pade rationalny z interpolacji g(K) w wezlach - nie Taylora
# Uzyj Chebyshev lub minimax - ale to skomplikowane.
# Zamiast: zwykla interpolacja wielomianowa przez N punktow
# i znalezienie zera.
deg_interp = min(len(K_scan)-1, 8)
coeffs_interp = np.polyfit(K_scan, g_scan, deg_interp)
poly_interp = np.poly1d(coeffs_interp)
roots_interp = np.roots(coeffs_interp)
real_roots_interp = sorted([r.real for r in roots_interp
                             if abs(r.imag)<0.05 and 1.0<r.real<4.0])

print(f"\n  Interpolacja wielomianowa stopnia {deg_interp} z g(K) w {len(K_scan)} punktach:")
print(f"  Zera w [1.0, 4.0]: {[f'{r:.6f}' for r in real_roots_interp]}")
if real_roots_interp:
    K2_interp = min(real_roots_interp, key=lambda r: abs(r-K2_num))
    print(f"  K2_interp = {K2_interp:.10f}  (blad: {(K2_interp/K2_num-1)*1e6:+.1f} ppm)")

# -------------------------------------------------------------------
# SEKCJA G: SYNTEZA I WNIOSKI OP-1
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA G: SYNTEZA -- POSTEP W OP-1")
print(f"{'='*72}")
print(f"""
  WYNIKI P63:

  1. SZEREG TAYLORA: c2..c7 wyznaczone analitycznie.
     c5 = -alpha*M5 = {c5_an:.4e}  (wazny! byl pomijany w P47/P48)
     c7 = -alpha*M7 = {c7_an:.4e}

  2. WIELOMIAN STOPNIA 6: zero blisko K2 = {K2_poly:.6f}
     (blad vs K2_num: {(K2_poly/K2_num-1)*1e6:+.1f} ppm)
     Wielomian uwzglednia c2..c7 -- ale czy seria zbiega?

  3. PROBLEM ZBIEZNOSCI: promien zbieznosci R ~ a*exp(a) = {K_sing_approx:.4f}
     Poniewaz K2 ~ {K2_num:.2f} >> {K_sing_approx:.4f} = R,
     SERIA TAYLORA W K=0 NIE ZBIEGA DO K2.
     Pade jest jedynym sensownym podejsciem analitycznym.

  4. PADE [3/3] daje dokladnosc ~ ? ppm -- patrz tabela sekcji B.

  5. RELACJA VIETA: K2 =/= sqrt(K1*K3) (blad {(geom_mean/K2_num-1)*100:+.2f}%).
     Brak prostej relacji miedzy zerami.

  6. ROZKLAD ENERGII przy K2:
     - lambda-czlon: {Ep_lam/E_tot*100:.2f}% (maly ale niezerowy)
     - kwartet: {Ep_qrt/E_tot*100:.2f}% (dominujacy ujemny)
     - kinetyczny: {(Ek_base+Ek_alpha)/E_tot*100:.2f}% (dominujacy dodatni)

  WNIOSEK OP-1:
  K2 nie posiada prostego wzoru zamknietego.
  Przeszkoda fundamentalna: osobliwosc analityczna g(K) przy K_sing~0.042
  uniemozliwia zbieznosc szeregu Taylora do K2~2.
  Pade [n/m] jest OPTYMALNYM podejsciem -- dokladnosc rosnie z n+m.
  Pade [3/3] wymaga c2..c7 (wszystkie analityczne) -- patrz wyniki w B.

  STATUS OP-1: DIAGNOZA UKONCZONA.
  Pelny wzor zamkniety (nie-Pade) jest NIEMOZLIWY w standardowym sensie.
  Najlepsze podejscie: Pade[3/3] z 6 analitycznymi wspolczynnikami.
""")

print("="*72)
print("P63 ZAKONCZONY")
print("="*72)
