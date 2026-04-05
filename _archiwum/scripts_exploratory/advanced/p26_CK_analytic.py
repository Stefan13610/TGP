"""
p26_CK_analytic.py
==================
CEL: Analityczne wyprowadzenie stalej C_K = 2.351 w formule K*1 = C_K * a_gam/(1+alpha).

PYTANIE: Skad sie bierze C_K = 2.351?

ODPOWIEDZ:
  C_K pochodzi z warunku samospojnosci TGP i ma dwie skladowe:
    1) C_K^pert  -- z perturbacyjnego rozwiazania (analitycznie)
    2) Korekcja  -- z nieliniowego czlonu kinetycznego (1+alpha/phi)

PLAN:
  Krok 1: Wyprowadz I_k(a) = integral_a^inf e^{-2r}(r+1)^2/r^2 dr  (dokladnie analitycznie)
  Krok 2: Pokaz ze I_k(a) = e^{-2a}(1/a + 1/2)  [termy E1 znoszace sie!]
  Krok 3: Wstaw do warunku samospojnosci => C_K^pert = 2*e^{2a} / (1 + alpha*a/(2*(1+alpha)))
  Krok 4: Oblicz korekte od nieliniowosci (1+alpha/phi vs 1+alpha)
  Krok 5: Pelna formula C_K = C_K^pert * (1 + delta_nl)
"""

import numpy as np
from scipy.special import exp1
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE I PARAMETRY
# ============================================================
R_MAX = 60.0
GAMMA = 1.0
FAMILY_DATA = [
    (5.9148, 0.025, 2.883e-6),
    (6.8675, 0.030, 3.743e-6),
    (7.7449, 0.035, 4.618e-6),
    (8.5616, 0.040, 5.501e-6),
]

print("P26: Analityczne wyprowadzenie stalej C_K")
print("=" * 65)
print()

# ============================================================
# KROK 1: Calka I_k(a) -- wynik scisly
# ============================================================
print("KROK 1: Scisle analityczne obliczenie I_k(a)")
print("-" * 65)
print()
print("  I_k(a) = integral_a^inf e^{-2r} * (r+1)^2/r^2 * dr")
print()
print("  Rozwijamy licznik: (r+1)^2/r^2 = 1 + 2/r + 1/r^2")
print()
print("  => I_k(a) = I_1 + 2*I_2 + I_3")
print()
print("  gdzie:")
print("    I_1 = integral_a^inf e^{-2r} dr          = e^{-2a}/2")
print("    I_2 = integral_a^inf e^{-2r}/r dr         = E_1(2a)")
print("    I_3 = integral_a^inf e^{-2r}/r^2 dr       = e^{-2a}/a - 2*E_1(2a)")
print()
print("    [I_3 przez calkowanie przez czesci: d/dr[e^{-2r}/r] = -2e^{-2r}/r - e^{-2r}/r^2]")
print("     => e^{-2r}/r^2 = -d/dr[e^{-2r}/r] - 2*e^{-2r}/r")
print("     => I_3 = e^{-2a}/a - 2*E_1(2a)")
print()
print("  Suma:")
print("    I_k = I_1 + 2*I_2 + I_3")
print("        = e^{-2a}/2 + 2*E_1(2a) + e^{-2a}/a - 2*E_1(2a)")
print()
print("    TERMY E_1 ZNOSZĄ SIE DOKLADNIE!")
print()
print("  >>> I_k(a) = e^{-2a} * (1/a + 1/2)  [wynik scisly] <<<")
print()

def I_k_exact(a):
    """Wynik analityczny: e^{-2a}*(1/a + 1/2)"""
    return np.exp(-2*a) * (1.0/a + 0.5)

def I_k_numeric(a):
    from scipy.integrate import quad
    result, _ = quad(lambda r: np.exp(-2*r)*(r+1)**2/r**2, a, np.inf)
    return result

def I_k_via_E1(a):
    """Via E1: I1 + 2*I2 + I3."""
    I1 = np.exp(-2*a)/2.0
    I2 = exp1(2*a)
    I3 = np.exp(-2*a)/a - 2.0*exp1(2*a)
    return I1 + 2*I2 + I3

print("  Weryfikacja numeryczna:")
print(f"  {'a':>8} {'I_k(scisly)':>14} {'I_k(numer.)':>14} {'blad':>12}")
print("  " + "-"*55)
for a in [0.010, 0.025, 0.030, 0.035, 0.040, 0.100, 0.500]:
    exact = I_k_exact(a)
    num   = I_k_numeric(a)
    err   = abs(exact - num)/num * 100
    print(f"  {a:>8.4f} {exact:>14.8f} {num:>14.8f} {err:>12.2e}%")
print()

# ============================================================
# KROK 2: Formula perturbacyjna C_K^pert
# ============================================================
print("KROK 2: Perturbacyjna formula C_K^pert")
print("-" * 65)
print()
print("  Warunek samospojnosci TGP: E(K*1)/(4*pi*K*1) = 1")
print()
print("  Dla malego K: phi = 1 + K*exp(-r)/r ~ 1,  dphi = -K*e^{-r}*(r+1)/r^2")
print()
print("  Energia kinetyczna (0-ty rzad w K):")
print("    E_k ~ 2*pi*(1+alpha)*K^2 * I_k(a_gam)   [dla phi ~ 1]")
print()
print("  Energia potencjalna:")
print("    E_p ~ -2*pi*K^2 * I_p(a_gam)")
print("    I_p(a) = integral_a^inf e^{-2r} dr = e^{-2a}/2")
print()
print("  Warunek E = 4*pi*K*1:")
print("    K*1 = 2 / [(1+alpha)*I_k(a) - I_p(a)]")
print()
print("  Wstawiamy I_k = e^{-2a}(1/a+1/2), I_p = e^{-2a}/2:")
print("    denom = e^{-2a} * [(1+alpha)/a + alpha/2]")
print("          = e^{-2a} * (1+alpha)/a * [1 + alpha*a/(2*(1+alpha))]")
print()
print("    K*1 = 2*a*e^{2a} / [(1+alpha) * (1 + alpha*a/(2*(1+alpha)))]")
print("        = 2*a_gam/(1+alpha) * C_K^pert")
print()
print("  >>> C_K^pert = 2*e^{2*a_gam} / [1 + alpha*a_gam/(2*(1+alpha))] <<<")
print()

def C_K_pert(alpha, a):
    return 2.0 * np.exp(2*a) / (1.0 + alpha*a/(2*(1+alpha)))

print("  Wartosci C_K^pert na rodzinie optymalnej:")
print(f"  {'alpha':>8} {'a_gam':>7}  {'e^(2a)':>8}  {'korektor':>10}  {'C_K^pert':>10}")
print("  " + "-"*55)
CK_pert_vals = []
for alpha, a, lam in FAMILY_DATA:
    exp_factor = np.exp(2*a)
    denom_corr = 1.0 + alpha*a/(2*(1+alpha))
    ckp = C_K_pert(alpha, a)
    CK_pert_vals.append(ckp)
    print(f"  {alpha:>8.4f} {a:>7.4f}  {exp_factor:>8.5f}  {denom_corr:>10.6f}  {ckp:>10.6f}")
print()

# ============================================================
# KROK 3: Korekcja od nieliniowosci
# ============================================================
print("KROK 3: Korekcja od nieliniowego czlonu kinetycznego")
print("-" * 65)
print()
print("  Przyblizenie phi~1 jest nieprecyzyjne: phi(a_gam) = 1 + delta,")
print("  gdzie delta = K*exp(-a)/a -- dla K=K*1 delta ~ 0.24-0.33!")
print()

def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=2000, linearized=False):
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    if linearized:
        Ek = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha) * r**2, r)
    else:
        Ek = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam, linearized=False):
    E = energy_log(K, alpha, a_gam, lam, linearized=linearized)
    return E / (4*np.pi*K) - 1.0

def find_K1(alpha, a_gam, lam, linearized=False, K_lo=0.003, K_hi=0.025):
    try:
        g_lo = g_func(K_lo, alpha, a_gam, lam, linearized)
        g_hi = g_func(K_hi, alpha, a_gam, lam, linearized)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi > 0:
            for K_lo2, K_hi2 in [(0.001, 0.030), (0.001, 0.050)]:
                g1 = g_func(K_lo2, alpha, a_gam, lam, linearized)
                g2 = g_func(K_hi2, alpha, a_gam, lam, linearized)
                if np.isfinite(g1) and np.isfinite(g2) and g1*g2 < 0:
                    return brentq(lambda K: g_func(K, alpha, a_gam, lam, linearized),
                                  K_lo2, K_hi2, xtol=1e-8)
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam, linearized),
                      K_lo, K_hi, xtol=1e-8)
    except Exception:
        return np.nan

print("  Porownanie K*1: pelny vs linearyzacja phi~1 vs formuła pert.:")
print(f"  {'alpha':>8} {'a_gam':>7}  {'K1_pelny':>10}  {'K1_lin':>10}  {'K1_pert':>10}  {'nl/lin':>8}  {'lin/pert':>10}")
print("  " + "-"*78)

CK_full_vals = []
CK_lin_vals  = []
for alpha, a, lam in FAMILY_DATA:
    K1_full = find_K1(alpha, a, lam, linearized=False)
    K1_lin  = find_K1(alpha, a, lam, linearized=True)
    K1_pert_val = C_K_pert(alpha, a) * a / (1+alpha)
    CK_f = K1_full*(1+alpha)/a if not np.isnan(K1_full) else np.nan
    CK_l = K1_lin*(1+alpha)/a  if not np.isnan(K1_lin)  else np.nan
    nl_corr  = CK_f / CK_l    if (not np.isnan(CK_f) and CK_l > 0) else np.nan
    lin_pert = CK_l / C_K_pert(alpha, a) if not np.isnan(CK_l) else np.nan
    CK_full_vals.append(CK_f)
    CK_lin_vals.append(CK_l)
    print(f"  {alpha:>8.4f} {a:>7.4f}  {K1_full:>10.7f}  {K1_lin:>10.7f}  {K1_pert_val:>10.7f}  {nl_corr:>8.5f}  {lin_pert:>10.5f}")

CK_full_mean = np.nanmean(CK_full_vals)
CK_lin_mean  = np.nanmean(CK_lin_vals)
CK_pert_mean = np.mean(CK_pert_vals)
print()
print(f"  C_K (pelny):  {CK_full_mean:.5f}")
print(f"  C_K (lin):    {CK_lin_mean:.5f}")
print(f"  C_K (pert):   {CK_pert_mean:.5f}")
print()

print("  Analiza korekcji nieliniowej -- phi na wewnetrznej granicy:")
print()
for alpha, a, lam in FAMILY_DATA:
    K1 = find_K1(alpha, a, lam, linearized=False)
    if not np.isnan(K1):
        phi_inner = 1.0 + K1*np.exp(-a)/a
        delta_inner = K1*np.exp(-a)/a
        ratio_nl = (1+alpha/phi_inner)/(1+alpha)
        print(f"  alpha={alpha:.4f}, a={a:.4f}: K*1={K1:.6f}")
        print(f"    phi(a_gam) = 1 + {delta_inner:.4f} = {phi_inner:.4f}")
        print(f"    (1+alpha/phi) / (1+alpha) = {ratio_nl:.5f}")
        print(f"    lokalna korekcja: {(ratio_nl-1)*100:+.2f}%")
print()

korekcja_nl  = CK_full_mean / CK_lin_mean  if CK_lin_mean  > 0 else np.nan
korekcja_lp  = CK_lin_mean  / CK_pert_mean if CK_pert_mean > 0 else np.nan
korekcja_tot = CK_full_mean / CK_pert_mean if CK_pert_mean > 0 else np.nan

print(f"  PODSUMOWANIE KOREKCJI:")
print(f"    Korekcja od nieliniowosci (full/lin):     {korekcja_nl:.5f}  ({(korekcja_nl-1)*100:+.2f}%)")
print(f"    Korekcja od potencjalu (lin/pert):         {korekcja_lp:.5f}  ({(korekcja_lp-1)*100:+.2f}%)")
print(f"    Laczna korekcja (full/pert):               {korekcja_tot:.5f}  ({(korekcja_tot-1)*100:+.2f}%)")
print()

# ============================================================
# KROK 4: Fizyczne wytlumaczenie korekcji
# ============================================================
print("KROK 4: Fizyczne wytlumaczenie korekcji ~11.7%")
print("-" * 65)
print()
print("  Linearyzacja phi ~ 1 daje: 1+alpha/phi ~ 1+alpha")
print("  ale dla phi(a_gam) = 1 + delta z delta ~ 0.24:")
print()
print("    1 + alpha/phi  =  1 + alpha/(1+delta)")
print("                   =  (1+alpha+alpha*delta/(1+delta))^{-1} * (1+alpha)")
print("                   ~ (1+alpha) * [1 - alpha*delta/((1+alpha)*(1+delta))]")
print()
print("  Lokalny czynnik jest MNIEJSZY niz (1+alpha), wiec I_k^eff < (1+alpha)*I_k.")
print("  Warunek samospojnosci E = 4*pi*K wymaga wiekszego K*1 aby skompensowac.")
print()
print("  Efektywnie: (1+alpha)*I_k -> integral (1+alpha/phi)*dphi^2*r^2 dr")
print("  gdzie phi(r) = 1 + K*e^{-r}/r, wiec dla r ~ a_gam: phi > 1.")
print()

# ============================================================
# KROK 5: Podsumowanie
# ============================================================
print("=" * 65)
print("KROK 5: PELNA FORMULA I PODSUMOWANIE")
print("=" * 65)
print()
print("  +----------------------------------------------------------+")
print("  |  SKAD SIE BIERZE C_K = 2.351?                           |")
print("  +----------------------------------------------------------+")
print("  |                                                          |")
print("  |  KROK ANALITYCZNY -- scisly wynik:                      |")
print("  |    I_k(a) = integral e^{-2r}(r+1)^2/r^2 dr             |")
print("  |           = e^{-2a}(1/a + 1/2)                         |")
print("  |    [termy E_1 znoszą sie dokladnie!]                    |")
print("  |                                                          |")
print("  |  FORMULA PERTURBACYJNA (zerowy rzad, phi~1):            |")
print("  |    C_K^pert = 2*e^{2*a_gam}                            |")
print("  |               / [1 + alpha*a_gam/(2*(1+alpha))]         |")
print(f"  |             = {CK_pert_mean:.4f}  (na rodzinie optymalnej)     |")
print("  |                                                          |")
print("  |  KOREKCJA NIELINIOWA (phi(a_gam) ~ 1.24, nie 1):       |")
print(f"  |    C_K = C_K^pert * {korekcja_tot:.5f}                          |")
print(f"  |        = {CK_pert_mean:.4f} * {korekcja_tot:.5f}                          |")
print(f"  |        = {CK_full_mean:.5f}  [numery.]                         |")
print("  |                                                          |")
print("  |  FORMULA KONCOWA:                                       |")
print("  |    K*1 = C_K * a_gam/(1+alpha)                          |")
print(f"  |        = {CK_full_mean:.4f} * a_gam/(1+alpha)                   |")
print("  |                                                          |")
print("  |  C_K NIE jest 'magiczna stala' -- jest wyznaczona       |")
print("  |  analitycznie przez geometrie profilu Yukawa             |")
print("  |  i nieliniowy czlon kinetyczny.                         |")
print("  +----------------------------------------------------------+")
print()

# Tabela zbiorczа
print("  Tabela C_K^pert vs C_K^num na rodzinie:")
print()
print(f"  {'alpha':>8} {'a_gam':>7}  {'C_K^pert':>10}  {'C_K^num':>10}  {'ratio':>8}  {'delta[%]':>10}")
print("  " + "-"*65)
for (alpha, a, lam), ckp, ckf in zip(FAMILY_DATA, CK_pert_vals, CK_full_vals):
    ratio = ckf/ckp if ckp > 0 else np.nan
    delta_pct = (ratio - 1)*100
    print(f"  {alpha:>8.4f} {a:>7.4f}  {ckp:>10.5f}  {ckf:>10.5f}  {ratio:>8.5f}  {delta_pct:>10.2f}%")
print()

# Asymptotyka
print("  Asymptotyka C_K^pert:")
print()
print("    a_gam -> 0:   e^{2a} -> 1,  denom -> 1  =>  C_K^pert -> 2")
print("    a_gam > 0:    korekcja ~8% na rodzinie ({:.1f}% od e^{{2a}}, {:.1f}% od mianownika)".format(
    (np.exp(2*0.04)-1)*100, (1/(1+8.56*0.04/(2*9.56))-1)*100))
print()
print("    Wiec C_K = 2 + korekcje geometryczne + korekcja nieliniowa")
print(f"             = 2 * {np.exp(2*0.04)/1:.4f} * {korekcja_tot/np.exp(2*0.04)*np.exp(2*0.04):.4f}")
print()

print("GOTOWE: P26 zakonczone.")
