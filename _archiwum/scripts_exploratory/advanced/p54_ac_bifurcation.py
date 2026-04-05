# -*- coding: utf-8 -*-
"""
P54: WYZNACZENIE a_c -- BIFURKACJA DOSTEPNOSCI KOIDEGO

Pytanie OP-6: jaka jest dokladna wartosc a_c, dla ktorej Q_min(a_c) = 3/2?

Dla a < a_c: Q(alpha) > 3/2 dla wszystkich alpha  =>  Q=3/2 niedostepne
Dla a = a_c: Q_min = 3/2 (punkt dotyku = alpha*)
Dla a > a_c: Q_min < 3/2  =>  dwa zera alpha1, alpha2 istnieja

Podejscie:
  A: Skanowanie Q_min(a) i weryfikacja struktury bifurkacji
  B: Bisekacja na Q_min(a) = 3/2 => a_c (precyzja 1e-8)
  C: Wlasciwosci w punkcie bifurkacji: alpha*, K1*, K2*, K3*, r21*
  D: Badanie czy a_c ma postac analityczna (zwiazki z pi, e, lambda, itp.)
  E: Jak szybko zera alpha1, alpha2 oddalaja sie od alpha* po bifurkacji?
     Wyznaczanie eksponentu delta_alpha ~ (a - a_c)^beta
  F: Predykcja: czy a_Gamma = 0.040 jest blisko a_c? Jak blisko?
  G: Wnioski i nowe hipotezy
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expi

# -------------------------------------------------------------------
# PARAMETRY BAZOWE
# -------------------------------------------------------------------
LAM_K   = 5.4677e-6
R_MAX   = 50.0
N_GRID  = 3000

def E1(x):
    return -expi(-x)

def V_mod(phi, lam=LAM_K):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha, a, lam=LAM_K, N=N_GRID):
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2*(1.0 + alpha/phi)*r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a, lam=LAM_K):
    return energy_num(K, alpha, a, lam) / (4*np.pi*K) - 1.0

def find_K1(alpha, a, lam=LAM_K):
    try:
        return brentq(g_func, 1e-4, 0.4, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def find_K2(alpha, a, lam=LAM_K):
    try:
        return brentq(g_func, 0.5, 5.0, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def K3_ei(a, lam=LAM_K):
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=200)
    return np.sqrt(3*I4 / (2*lam*I6)) if I6 > 0 else np.nan

def Q_from_alpha(alpha, a, lam=LAM_K, K3_fixed=None):
    K1 = find_K1(alpha, a, lam)
    K2 = find_K2(alpha, a, lam)
    if np.isnan(K1) or np.isnan(K2):
        return np.nan
    K3 = K3_fixed if K3_fixed is not None else K3_ei(a, lam)
    if np.isnan(K3):
        return np.nan
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    d = K1 + K2 + K3
    return s**2 / d

def Q_min_val(a, lam=LAM_K, alpha_range=(0.2, 20.0), n_scan=80):
    """Wyznacza minimum Q(alpha) dla danego a, przez grube skanowanie + optymalizacje."""
    K3 = K3_ei(a, lam)
    if np.isnan(K3):
        return np.nan, np.nan

    alpha_scan = np.linspace(alpha_range[0], alpha_range[1], n_scan)
    Q_scan = np.array([Q_from_alpha(al, a, lam, K3) for al in alpha_scan])

    valid = ~np.isnan(Q_scan)
    if valid.sum() == 0:
        return np.nan, np.nan

    # Gruby minimum
    idx_min = np.argmin(Q_scan[valid])
    alpha_min_rough = alpha_scan[valid][idx_min]

    # Dokladny minimum przez minimize_scalar
    try:
        res = minimize_scalar(
            lambda al: Q_from_alpha(al, a, lam, K3),
            bounds=(max(0.1, alpha_min_rough - 2.0),
                    min(alpha_range[1], alpha_min_rough + 2.0)),
            method='bounded', options={'xatol': 1e-9}
        )
        return res.fun, res.x
    except Exception:
        return Q_scan[valid][idx_min], alpha_min_rough


# -------------------------------------------------------------------
print("=" * 65)
print("P54: BIFURKACJA DOSTEPNOSCI KOIDEGO -- WYZNACZENIE a_c")
print("=" * 65)
print(f"  lambda = {LAM_K:.4e}")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA A: SKANOWANIE Q_min(a) -- WERYFIKACJA STRUKTURY")
print("=" * 65)

a_scan = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.032, 0.034,
                   0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.045,
                   0.050, 0.060, 0.080, 0.100])

print(f"\n  {'a':>7}  {'Q_min':>10}  {'alpha*':>8}  {'Q_min-3/2':>12}  {'zera Q=3/2':>10}")
Qmin_arr = []
alpha_star_arr = []
for a in a_scan:
    Qm, am = Q_min_val(a)
    Qmin_arr.append(Qm)
    alpha_star_arr.append(am)
    sign = "brak" if Qm > 1.5 else "dwa"
    print(f"  {a:7.4f}  {Qm:10.6f}  {am:8.4f}  {(Qm-1.5):+12.6f}  {sign:>10}")

Qmin_arr = np.array(Qmin_arr)
alpha_star_arr = np.array(alpha_star_arr)

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA B: BISEKACJA a_c (Q_min(a_c) = 3/2)")
print("=" * 65)

# Znajdz przedzialy gdzie Q_min przekracza 3/2
# Z Sekcji A: Q_min(0.030)>3/2, Q_min(0.040)<3/2
def Qmin_minus_32(a):
    Qm, _ = Q_min_val(a, n_scan=100)
    return Qm - 1.5

# Grube skanowanie z duzym n_scan dla dokladnosci
a_lo, a_hi = 0.030, 0.040

print(f"\n  Szukam a_c w przedziale ({a_lo}, {a_hi})")
print(f"  Q_min({a_lo}) - 3/2 = {Qmin_minus_32(a_lo):+.6f}")
print(f"  Q_min({a_hi}) - 3/2 = {Qmin_minus_32(a_hi):+.6f}")

# Bisekacja
a_c = brentq(Qmin_minus_32, a_lo, a_hi, xtol=1e-8, rtol=1e-10)
Qm_c, alpha_c_star = Q_min_val(a_c, n_scan=200)

print(f"\n  *** WYNIK BISEKACJI ***")
print(f"  a_c     = {a_c:.10f}")
print(f"  Q_min   = {Qm_c:.8f}  (cel: 1.50000000)")
print(f"  alpha*  = {alpha_c_star:.6f}  (alpha przy Q_min = 3/2)")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA C: WLASCIWOSCI W PUNKCIE BIFURKACJI")
print("=" * 65)

K3_c = K3_ei(a_c, LAM_K)
K1_c = find_K1(alpha_c_star, a_c, LAM_K)
K2_c = find_K2(alpha_c_star, a_c, LAM_K)

print(f"\n  Punkt bifurkacji (a_c, alpha*):")
print(f"    a_c    = {a_c:.8f}")
print(f"    alpha* = {alpha_c_star:.6f}")
print(f"    K1*    = {K1_c:.8f}")
print(f"    K2*    = {K2_c:.8f}")
print(f"    K3*    = {K3_c:.6f}")
if not np.isnan(K1_c) and not np.isnan(K2_c):
    r21_c = K2_c / K1_c
    r31_c = K3_c / K1_c
    print(f"    r21*   = K2*/K1* = {r21_c:.4f}")
    print(f"    r31*   = K3*/K1* = {r31_c:.2f}")
    Qc_check = (np.sqrt(K1_c)+np.sqrt(K2_c)+np.sqrt(K3_c))**2 / (K1_c+K2_c+K3_c)
    print(f"    Q*     = {Qc_check:.8f}  (weryfikacja)")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA D: POSZUKIWANIE POSTACI ANALITYCZNEJ a_c")
print("=" * 65)

print(f"\n  a_c = {a_c:.10f}")
print(f"\n  Czy a_c jest prostym ulamkiem?")
for denom in range(1, 51):
    numer = round(a_c * denom)
    if numer > 0:
        approx = numer / denom
        err = abs(approx - a_c) / a_c * 100
        if err < 0.1:
            print(f"    {numer}/{denom} = {approx:.8f}  (blad: {err:.4f}%)")

print(f"\n  Zwiazki z lambda:")
print(f"    lambda^(1/3) = {LAM_K**(1/3):.6e}")
print(f"    lambda^(1/4) = {LAM_K**(1/4):.6e}")
print(f"    a_c / lambda^(1/4) = {a_c / LAM_K**(1/4):.6f}")
print(f"    a_c / lambda^(1/6) = {a_c / LAM_K**(1/6):.6f}")

print(f"\n  Zwiazki z pi, e, sqrt(2):")
candidates = {
    'pi/100': np.pi/100,
    '1/(10*pi)': 1/(10*np.pi),
    'pi/90':  np.pi/90,
    'pi^2/300': np.pi**2/300,
    '1/e^(1/a_c)': None,  # placeholder
    'sqrt(2)/40': np.sqrt(2)/40,
    'sqrt(3)/50': np.sqrt(3)/50,
    '1/30': 1/30,
    '1/31': 1/31,
    '1/32': 1/32,
    '1/33': 1/33,
    '1/34': 1/34,
    '1/35': 1/35,
    'exp(-3.5)': np.exp(-3.5),
    'exp(-3.4)': np.exp(-3.4),
}
for name, val in candidates.items():
    if val is not None:
        err = abs(val - a_c) / a_c * 100
        print(f"    {name:25s} = {val:.8f}  (blad od a_c: {err:.3f}%)")

# Logarytmiczne poszukiwanie
print(f"\n  ln(1/a_c) = {np.log(1/a_c):.6f}")
print(f"  a_c^2 * 1000 = {a_c**2 * 1000:.6f}")
print(f"  a_c * 1/K3_c = {a_c / K3_c:.6e}")
print(f"  a_c * r21*   = {a_c * r21_c:.6f}" if not np.isnan(K1_c) else "")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA E: EKSPONENT BIFURKACJI delta_alpha ~ (a - a_c)^beta")
print("=" * 65)

# Dla a > a_c, dwa zera alpha1 < alpha* < alpha2
# Szerokos przedzialu delta_alpha = alpha2 - alpha1 ~ (a - a_c)^beta
# W punkcie bifurkacji siodlowo: spodziewamy sie beta = 1/2 (fold bifurcation)

epsilon_vals = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.010, 0.015, 0.020])
a_above = a_c + epsilon_vals

print(f"\n  {'a - a_c':>10}  {'alpha1':>9}  {'alpha2':>9}  {'alpha*':>9}  {'delta_a':>10}  {'delta_a/eps^0.5':>15}")

delta_alphas = []
valid_eps = []
for eps, a in zip(epsilon_vals, a_above):
    K3a = K3_ei(a, LAM_K)

    def Qa_m32(alpha, aa=a, K3a=K3a):
        return Q_from_alpha(alpha, aa, LAM_K, K3a) - 1.5

    # Skan
    al_scan = np.linspace(0.5, 15.0, 200)
    Q_s = np.array([Q_from_alpha(al, a, LAM_K, K3a) for al in al_scan])
    zs = []
    for i in range(len(Q_s)-1):
        if not (np.isnan(Q_s[i]) or np.isnan(Q_s[i+1])):
            if (Q_s[i]-1.5)*(Q_s[i+1]-1.5) < 0:
                try:
                    z = brentq(Qa_m32, al_scan[i], al_scan[i+1], xtol=1e-8)
                    zs.append(z)
                except Exception:
                    pass

    if len(zs) >= 2:
        z1, z2 = zs[0], zs[-1]
        da = z2 - z1
        ratio = da / np.sqrt(eps)
        print(f"  {eps:10.4f}  {z1:9.5f}  {z2:9.5f}  {alpha_c_star:9.5f}  {da:10.5f}  {ratio:15.4f}")
        delta_alphas.append(da)
        valid_eps.append(eps)
    else:
        print(f"  {eps:10.4f}  brak zer")

if len(valid_eps) >= 3:
    # Fit log(delta_a) = beta * log(eps) + const
    log_eps = np.log(valid_eps)
    log_da  = np.log(delta_alphas)
    beta, log_c = np.polyfit(log_eps, log_da, 1)
    C_fit = np.exp(log_c)
    print(f"\n  Fit: delta_alpha = {C_fit:.4f} * (a - a_c)^{beta:.4f}")
    print(f"  Eksponent beta = {beta:.4f}  (oczekiwane dla fold bifurcation: 0.5)")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA F: POROWNANIE a_c Z a_GAMMA")
print("=" * 65)

A_GAM = 0.040
print(f"\n  a_c    = {a_c:.8f}")
print(f"  a_Gamma= {A_GAM:.8f}")
print(f"  a_Gamma - a_c = {A_GAM - a_c:.6f}  ({(A_GAM-a_c)/a_c*100:.3f}% powyzej progu)")
print(f"  a_Gamma / a_c = {A_GAM / a_c:.6f}")

# Czy a_Gamma = a_c + male corekcie?
print(f"\n  Jak daleko a_Gamma jest od a_c?")
print(f"    eps = a_Gamma - a_c = {A_GAM - a_c:.6f}")
K3_gamma = K3_ei(A_GAM, LAM_K)
print(f"\n  Przy a_Gamma=0.040:")
print(f"    K3_ei = {K3_gamma:.4f}")
print(f"    Q_min = {Qmin_minus_32(A_GAM)+1.5:.6f}")

# Predykcja zer alpha1, alpha2 z eksponentu bifurkacji
if len(valid_eps) >= 3:
    eps_gamma = A_GAM - a_c
    da_pred = C_fit * eps_gamma**beta
    print(f"\n  Predykcja zer przy a_Gamma z eksponentu bifurkacji:")
    print(f"    delta_alpha (pred) = {da_pred:.4f}")
    print(f"    alpha1 (pred) = alpha* - da/2 = {alpha_c_star - da_pred/2:.4f}")
    print(f"    alpha2 (pred) = alpha* + da/2 = {alpha_c_star + da_pred/2:.4f}")
    print(f"    (P53 wyniki:  alpha1=1.5804, alpha2=8.4731)")

# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA G: WNIOSKI")
print("=" * 65)
print(f"""
  1. BIFURKACJA DOSTEPNOSCI:
     a_c = {a_c:.8f}
     Dla a < a_c: Q(alpha) > 3/2 dla wszystkich alpha  (Koide niedostepne)
     Dla a = a_c: Q_min = 3/2, jeden punkt dotyku alpha* = {alpha_c_star:.4f}
     Dla a > a_c: dwa zera Q=3/2 (alpha1, alpha2 powstajs z dotyku)

  2. PUNKT BIFURKACJI:
     alpha* = {alpha_c_star:.6f}  (alpha przy Q_min = 3/2)
     r21*   = {r21_c:.4f}  (K2*/K1* w punkcie bifurkacji)
     [dla lept.: r21=206.77; r21* jest 'poczatek' hierarchii mas]

  3. EKSPONENT BIFURKACJI:
     delta_alpha ~ (a - a_c)^beta,  beta ~ 0.5  (fold/saddle-node)
     Standard: beta=1/2 dla pit bifurkacji kwadrowej

  4. POLOZENIE a_Gamma:
     a_Gamma = {A_GAM:.4f}  jest  {(A_GAM-a_c)/a_c*100:.2f}% powyzej a_c
     Pytanie: czy istnieje zasada wyznaczajaca a_Gamma - a_c?
     Hipoteza: a_Gamma jest minimalna wartosc spelniajaca r21 >= r21_PDG

  5. NOWE PYTANIA (P55+):
     - Czy beta = 1/2 dokladnie? (kwadratowa osobliwosc Q_min)
     - Czy a_c ma postac analityczna przez I4, I6, lambda?
     - Czy warunek r21(a_c) = r21_PDG wyznacza a_Gamma?
""")

print("=" * 65)
print("P54 ZAKONCZONY")
print("=" * 65)
