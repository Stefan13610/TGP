# -*- coding: utf-8 -*-
"""
P61: MECHANIZM NIEZALEZNOSCI K1 OD LAMBDA (OP-10)

P56 odkryl: zmiana lambda x0.5 .. x2.0 nie zmienia K1 (delta < 0.001%).
Cel: wyjasniec przyczyne strukturalna.

Hipoteza: przy K1 ~ 0.01, czlon lambda*K1^5*I6 jest rzedy ~10^{-14}
          -- calkowicie zaniedbywalny w rownaniu g(K1)=0.
          Stad K1 wyznaczony przez SAME skladniki lambda-niezalezne.

Plan:
  A. Rozloz g(K) na skladniki przy K=K1; pokaz dominacje lambda-niezaleznych.
  B. Przeprowadz analize perturbacyjna:
     g(K) = g0(K) + lambda * g_lam(K)
     K1 = K1_0 + lambda * dK1/dlambda + ...
     Wyznacz K1_0 (zero g0) i dK1/dlambda analitycznie.
  C. Powtorz dla K2 i K3 -- pokaz ze przy K2, K3 czlon lambda jest istotny.
  D. Tabela: wzgledna wrazliwosc dlog(Ki)/dlog(lambda) dla i=1,2,3.
  E. Interpretacja fizyczna: "K1 yje w reziomie solitonu plastkiego".
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
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

# Rozklad energii na skladniki NIEZALEZNE i ZALEZNE od lambda
def g_components(K, alpha, a, lam, N=N_GRID):
    t    = np.linspace(0, 1, N)
    r    = a*(R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r-1.0)/r**2

    kin      = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    v_cubic  = 4*np.pi*np.trapezoid((phi**3-1.0)/3*r**2, r)
    v_quart  = 4*np.pi*np.trapezoid(-(phi**4-1.0)/4*r**2, r)
    v_lambda = 4*np.pi*np.trapezoid(lam*(phi-1.0)**6/6*r**2, r)

    E_nolam  = kin + v_cubic + v_quart          # lambda-niezalezne
    E_lam    = v_lambda                          # lambda-zalezne
    E_total  = E_nolam + E_lam
    g_nolam  = E_nolam / (4*np.pi*K)
    g_lam    = E_lam   / (4*np.pi*K)
    g_val    = g_nolam + g_lam - 1.0

    return dict(
        kin=kin, v_cubic=v_cubic, v_quart=v_quart, v_lambda=v_lambda,
        E_nolam=E_nolam, E_lam=E_lam, E_total=E_total,
        g_nolam=g_nolam, g_lam=g_lam, g_val=g_val
    )

def find_K(alpha, a, lam, bracket, N=N_GRID, tol=1e-12):
    def gfun(K): return g_components(K, alpha, a, lam, N)['g_val']
    try: return brentq(gfun, bracket[0], bracket[1], xtol=tol)
    except: return np.nan

# Parametry
A_GAM  = 0.040049   # Punkt B
ALPHA  = 8.5612
LAM_K  = 5.4677e-6

print("="*72)
print("P61: MECHANIZM NIEZALEZNOSCI K1 OD LAMBDA (OP-10)")
print("="*72)
print(f"  Punkt B: a={A_GAM}, alpha={ALPHA}, lambda={LAM_K:.5e}")

# -------------------------------------------------------------------
# SEKCJA A: Rozklad g(K) przy K=K1, K2, K3
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA A: Rozklad g(K) na skladniki przy K = K1, K2, K3")
print(f"{'='*72}")

K1 = find_K(ALPHA, A_GAM, LAM_K, [1e-5, 0.45])
K2 = find_K(ALPHA, A_GAM, LAM_K, [0.4,  5.0 ])
K3 = find_K(ALPHA, A_GAM, LAM_K, [5.0,  80.0])
print(f"\n  K1 = {K1:.12f}")
print(f"  K2 = {K2:.12f}")
print(f"  K3 = {K3:.12f}  (pelna numeryka)")

header = f"  {'K':>10}  {'g_nolam':>14}  {'g_lam':>14}  {'g_lam/g_val':>14}  {'lam*term/E':>14}"
print(f"\n{header}")
print(f"  {'-'*10}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*14}")

for Kval, label in [(K1,'K1'), (K2,'K2'), (K3,'K3')]:
    ec = g_components(Kval, ALPHA, A_GAM, LAM_K)
    lam_frac = ec['E_lam']/max(abs(ec['E_total']),1e-30)*100 if abs(ec['E_total'])>1e-30 else np.nan
    # udzial g_lam w 1 (bo g_nolam + g_lam = 1 w zerze)
    print(f"  {Kval:10.6f}  {ec['g_nolam']:+14.8f}  {ec['g_lam']:+14.10f}  "
          f"{ec['g_lam']/(ec['g_val']+1.0)*100:+14.6f}%  {lam_frac:+14.6f}%")

# Pokaz jawnie jak maly jest czlon lambda przy K1
ec1 = g_components(K1, ALPHA, A_GAM, LAM_K)
print(f"\n  Przy K1 = {K1:.6f}:")
print(f"    lambda * K1^5 * I6 / 6 ~ {LAM_K * K1**5 * 3685 / 6:.3e}  (zaniedbywalne)")
print(f"    g_lam(K1) / 1       = {ec1['g_lam']:+.3e}  ({ec1['g_lam']*1e9:+.3f} ppb)")
print(f"    g_nolam(K1)         = {ec1['g_nolam']:+.10f}")
print(f"    g_val(K1)           = {ec1['g_val']:+.3e}")

# -------------------------------------------------------------------
# SEKCJA B: Perturbacyjna analiza K1(lambda)
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA B: Analiza perturbacyjna K1(lambda)")
print(f"{'='*72}")

# g(K, lambda) = g0(K) + lambda * g1(K) - 1 = 0
# gdzie g0(K) = E_nolam/(4piK), g1(K) = E_lam/(4piK*lambda)
# W zerze K1: g0(K1) + lambda*g1(K1) = 1
# Pochodna po lambda:
# g0'(K1)*dK1/dlambda + g1(K1) + lambda*g1'(K1)*dK1/dlambda = 0
# dK1/dlambda = -g1(K1) / (g0'(K1) + lambda*g1'(K1))
# Przy malejacej roli lambda: ~ -g1(K1) / g0'(K1)

# Oblicz numerycznie g0'(K1) i g1(K1)
dK = K1 * 1e-5
ec_p = g_components(K1+dK, ALPHA, A_GAM, LAM_K)
ec_m = g_components(K1-dK, ALPHA, A_GAM, LAM_K)

dg0_dK = (ec_p['g_nolam'] - ec_m['g_nolam'])/(2*dK)
dg_dK  = (ec_p['g_val']   - ec_m['g_val'])  /(2*dK)   # pelna pochodna

# g1(K1) = g_lam(K1)/lambda
g1_K1 = ec1['g_lam'] / LAM_K

dK1_dlam       = -g1_K1 / dg_dK
dlogK1_dloglam = (LAM_K / K1) * dK1_dlam

print(f"\n  g0'(K1) = dg_nolam/dK |_K1 = {dg0_dK:+.6f}")
print(f"  g'(K1)  = dg/dK       |_K1 = {dg_dK:+.6f}  (pelna)")
print(f"  g1(K1)  = g_lam(K1)/lambda = {g1_K1:+.6e}")
print(f"\n  Czulosc: dK1/dlambda = {dK1_dlam:.6e}")
print(f"  Wzgledna: dlogK1/dlogLambda = {dlogK1_dloglam:+.6e}")
print(f"  -> zmiana lambda o 100% (x2) zmienia K1 o {abs(dlogK1_dloglam)*100:.4f}%")

# Weryfikacja numeryczna: oblicz K1 przy roznych lambda
print(f"\n  Weryfikacja numeryczna:")
print(f"  {'lambda/lambda_K':>18}  {'K1':>18}  {'delta K1 [ppm]':>16}")
print(f"  {'-'*18}  {'-'*18}  {'-'*16}")

K1_ref = K1
lam_ratios = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0]
for ratio in lam_ratios:
    lam_i = LAM_K * ratio
    K1_i  = find_K(ALPHA, A_GAM, lam_i, [1e-5, 0.45])
    delta  = (K1_i/K1_ref - 1)*1e6 if not np.isnan(K1_i) else np.nan
    print(f"  {ratio:18.3f}  {K1_i:18.14f}  {delta:+16.4f}")

# -------------------------------------------------------------------
# SEKCJA C: Porownanie K1, K2, K3 - wrazliwosc na lambda
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA C: Czulosc Ki na lambda -- porownanie K1, K2, K3")
print(f"{'='*72}")

# Oblicz dlogKi/dlogLambda numerycznie dla wszystkich Ki
print(f"\n  {'Ki':>6}  {'Ki_ref':>14}  {'dlogKi/dlogLam':>18}  {'inter. [ppm/x2lam]':>22}")
print(f"  {'-'*6}  {'-'*14}  {'-'*18}  {'-'*22}")

for Kval, bracket, label in [(K1,[1e-5,0.45],'K1'), (K2,[0.4,5.0],'K2'), (K3,[5.0,80.0],'K3')]:
    lam_hi = LAM_K * 1.01
    lam_lo = LAM_K * 0.99
    K_hi   = find_K(ALPHA, A_GAM, lam_hi, bracket)
    K_lo   = find_K(ALPHA, A_GAM, lam_lo, bracket)
    if np.isnan(K_hi) or np.isnan(K_lo):
        print(f"  {label:>6}  {Kval:14.8f}  {'BLAD':>18}")
        continue
    dlogK_dlogL = (np.log(K_hi) - np.log(K_lo)) / (np.log(lam_hi) - np.log(lam_lo))
    sens_ppm = dlogK_dlogL * np.log(2) * 1e6
    print(f"  {label:>6}  {Kval:14.8f}  {dlogK_dlogL:+18.6f}  {sens_ppm:+22.1f}")

# Teoretyczne: K3 ~ lambda^{-1/2} => dlogK3/dlogLam = -1/2
# K1, K2: dlogKi/dlogLam ~ 0

# -------------------------------------------------------------------
# SEKCJA D: Rownanie K1 bez lambda (K1_0) i korekcja
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA D: Rownanie K1_0 (g0=1, bez lambda) i korekcja pierwszorzedowa")
print(f"{'='*72}")

# g0(K) = E_nolam/(4piK) = 1  =>  K1_0
# Wyznaczyc K1_0 numerycznie
def g0_func(K):
    ec = g_components(K, ALPHA, A_GAM, 0.0)  # lambda=0
    return ec['g_val']  # == g_nolam(K) - 1

print("\n  Szukam K1_0 (zero g0 = 1, bez lambda)...")
try:
    K1_0 = brentq(g0_func, 1e-5, 0.45, xtol=1e-12)
    print(f"  K1_0 (lambda=0) = {K1_0:.14f}")
    print(f"  K1   (lambda=lam_K) = {K1:.14f}")
    print(f"  roznica: {(K1/K1_0-1)*1e6:+.4f} ppm  ({(K1/K1_0-1)*100:.8f}%)")
except Exception as e:
    print(f"  BLAD: {e}")
    K1_0 = np.nan

# Korekcja pierwszorzedowa: K1 ~ K1_0 + lam * dK1/dlam
if not np.isnan(K1_0):
    K1_pert = K1_0 + LAM_K * dK1_dlam
    print(f"\n  Korekcja perturbacyjna:")
    print(f"  K1_pert = K1_0 + lam * dK1/dlam = {K1_pert:.14f}")
    print(f"  K1_num                           = {K1:.14f}")
    print(f"  blad korekty: {(K1_pert/K1-1)*1e6:+.4f} ppm")

# -------------------------------------------------------------------
# SEKCJA E: Interpretacja fizyczna
# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("SEKCJA E: INTERPRETACJA FIZYCZNA -- 'SOLITON PLYTKI' przy K1")
print(f"{'='*72}")

# Ocenmy wzgledna wielkosc czlonu lambda w g(K) dla K=K1,K2,K3
print(f"""
  Czlon lambda w energii: E_lam = 4pi * lam * K^6/6 * I6

  Ocena E_lam/(4pi*K) wzgledem 1 (warunek g(K)=0):

  Dla K=K1={K1:.4f}:
    E_lam/(4pi*K1) ~ lam*K1^5*I6/6
                   ~ {LAM_K:.2e} * {K1:.4f}^5 * 3685 / 6
                   ~ {LAM_K * K1**5 * 3685/6:.3e}  <<  1

  Dla K=K2~2:
    E_lam/(4pi*K2) ~ lam*K2^5*I6/6
                   ~ {LAM_K:.2e} * {K2:.4f}^5 * 3685 / 6
                   ~ {LAM_K * K2**5 * 3685/6:.3e}  (porownywalne z innymi czl.)

  Dla K=K3~34:
    E_lam/(4pi*K3) ~ lam*K3^5*I6/6
                   ~ {LAM_K:.2e} * {K3:.4f}^5 * 3685 / 6
                   ~ {LAM_K * K3**5 * 3685/6:.3e}  (DOMINUJACY!)
""")

ec1_v = g_components(K1, ALPHA, A_GAM, LAM_K)
ec2_v = g_components(K2, ALPHA, A_GAM, LAM_K)
ec3_v = g_components(K3, ALPHA, A_GAM, LAM_K)

print(f"  Frakcja E_lam/E_total przy zerach Ki:")
print(f"  {'Ki':>6}  {'K':>10}  {'E_lam/|E_total|':>18}  {'g_lam [ppb]':>14}")
for ec, label, Kv in [(ec1_v,'K1',K1),(ec2_v,'K2',K2),(ec3_v,'K3',K3)]:
    frac = ec['E_lam']/max(abs(ec['E_total']),1e-30)*100
    glam_ppb = ec['g_lam']*1e9
    print(f"  {label:>6}  {Kv:10.5f}  {frac:+18.6f}%  {glam_ppb:+14.3f}")

print(f"""
  Wniosek strukturalny:
  - K1 ~ 0.010: lambda*K1^6 ~ 5e-6 * 1e-12 = 5e-18 -- calkowicie zaniedbywalny
    K1 wyznaczony przez rownowage czlonow kinetycznego i potencjalowego bez lambda.
  - K2 ~ 2:    lambda*K2^6 ~ 5e-6 * 64    = 3e-4   -- mala korekcja (~0.03%)
  - K3 ~ 34:   lambda*K3^6 ~ 5e-6 * 1.5e9 = 7.6e3  -- calkowicie dominujacy!
    K3 ~ 1/sqrt(lambda) -- silna zaleznosc.

  K1 lezy w reziomie "solitonu plaskego" (K1<<1, phi~1+K1*f gdzie f<<1)
  gdzie potencjal lambda*(phi-1)^6 jest pomijanie malym rzedowo.
  Fizycznie: elektrony odpowiadaja solitonowi platemu, gdzie lambda
  (odpowiadajace masie Higgsa) nie wplywa na skale solitonu.
""")

# -------------------------------------------------------------------
# SEKCJA F: Wzor zamkniety na K1_0 (lambda=0)
# -------------------------------------------------------------------
print(f"{'='*72}")
print("SEKCJA F: WZOR ZAMKNIETY NA K1_0 (lambda->0)")
print(f"{'='*72}")

# g0(K) = E_nolam/(4piK) - 1 = 0
# E_nolam = Ek + V_cubic + V_quart
# Perturbacyjnie w K (maly K):
# E_nolam ~ 4pi [K^2*(1+alpha)*C2 - K^2*I2/2 - 2K^3*I3/3 - K^4*I4/4]
# g0(K) ~ K*(1+alpha)*C2 - K*I2/2 - 2K^2*I3/3 - K^3*I4/4 - 1 = 0

# Oblicz C2, I2, I3, I4 analitycznie/numerycznie:
a = A_GAM
I2, _ = quad(lambda r: np.exp(-2*r),       a, R_MAX, limit=300)  # = e^{-2a}/2
I3    = E1(3*a)                                                    # = E1(3a) -- uwaga: int e^{-3r}/r = E1(3a)
I4    = np.exp(-4*a)/a - 4*E1(4*a)                                # = int e^{-4r}/r^2

# C2 = int e^{-2r}(r+1)^2/(2r^2) dr  (numerycznie)
C2, _ = quad(lambda r: np.exp(-2*r)*(r+1)**2/(2*r**2), a, R_MAX, limit=500)
# J_alpha (alpha/phi term) -- przy phi~1 dla malego K:
Ca, _ = quad(lambda r: np.exp(-2*r)*(r+1)**2/(2*r**2), a, R_MAX, limit=500)  # identyczne z C2 przy K->0

print(f"\n  Calki analityczne:")
print(f"  I2   = {I2:.10f}  (= e^{{-2a}}/2 = {np.exp(-2*a)/2:.10f})")
print(f"  I3   = {I3:.10f}  (= E1(3a))")
print(f"  I4   = {I4:.10f}")
print(f"  C2   = {C2:.10f}  (calka kinetyczna)")

# Rownanie wiodace (K^1 i staly):
# K * [(1+alpha)*C2 - I2/2] = 1 + 2K^2*I3/3 + K^3*I4/4 + ...
# Przyblizone: K1_approx ~ 1 / [(1+alpha)*C2 - I2/2]
A1 = (1.0 + ALPHA)*C2 - I2/2.0
K1_lead = 1.0 / A1

print(f"\n  Wspolczynnik wiodacy A1 = (1+alpha)*C2 - I2/2 = {A1:.10f}")
print(f"  K1_leading = 1/A1 = {K1_lead:.10f}")
print(f"  K1_num     = {K1:.10f}")
print(f"  blad K1_leading vs K1_num: {(K1_lead/K1-1)*1e6:+.1f} ppm  ({(K1_lead/K1-1)*100:.4f}%)")

# NLO korekcja: K1 = K1_lead / (1 + 2*K1_lead*I3/3 * 1/A1 + ...)
# Iteracja: K1_nlo = 1/(A1 + 2K1_lead*I3/3)
K1_nlo = 1.0 / (A1 + 2.0*K1_lead*I3/3.0)
print(f"\n  NLO korekcja: K1_NLO = 1/(A1 + 2*K1_lead*I3/3)")
print(f"  K1_NLO  = {K1_nlo:.10f}")
print(f"  blad K1_NLO vs K1_num: {(K1_nlo/K1-1)*1e6:+.1f} ppm  ({(K1_nlo/K1-1)*100:.4f}%)")

# Iteracja do zbieznosci
K1_iter = K1_lead
for step in range(10):
    K1_new = 1.0 / (A1 + 2.0*K1_iter*I3/3.0 + K1_iter**2*I4/4.0)
    if abs(K1_new/K1_iter - 1) < 1e-14:
        K1_iter = K1_new
        break
    K1_iter = K1_new
print(f"\n  Iteracja analityczna (do zbieznosci):")
print(f"  K1_analyt  = {K1_iter:.14f}")
print(f"  K1_num     = {K1:.14f}")
print(f"  K1_0(brentq)= {K1_0:.14f}")
print(f"  blad vs K1_num:    {(K1_iter/K1-1)*1e6:+.4f} ppm")
print(f"  blad vs K1_0:      {(K1_iter/K1_0-1)*1e6:+.4f} ppm")

# -------------------------------------------------------------------
print(f"\n{'='*72}")
print("P61 ZAKONCZONY")
print(f"{'='*72}")
