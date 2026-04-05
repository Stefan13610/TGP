"""
p37_K1_precise.py
=================
CEL: Precyzyjna formula dla K1(alpha, a_Gamma) — glowne zrodlo bledu r21

KONTEKST (P36):
  Glowny blad r21_NLO (max 15.4%) pochodzi z K1_analytic = C_K*a/(1+alpha),
  gdzie C_K = 2.351 jest skalibrowne tylko przy alpha=8.553.
  Przy alpha=3: blad K1 = +9.4%; przy alpha=12: blad K1 = -1.3%.

  Obliczona wyzej struktura:
    K1 = 2/Delta0 * (1 + NLO_correction)
    Delta0 = (1+alpha)*e^{-2a}*(2+a)/(2a) - e^{-2a}/2
            = e^{-2a}*[(1+alpha)*(2+a)/2a - 1/2]
            = e^{-2a}*(2+3a+alpha*(2+a)) / (2a)
  C_K(alpha,a) = K1_num*(1+alpha)/a  -- empiryczna stala

PLAN P37:
  Czesc A -- Mapa C_K(alpha, a_Gamma)
    Siatka alpha x [1.5, 20] x a x [0.010, 0.060]
    Oblicz K1_num numerycznie, C_K(alpha,a) = K1_num*(1+alpha)/a
    Zidentyfikuj zmienna naturalna dla C_K

  Czesc B -- Analityczne NLO dla K1
    Wzor: K1_LO = 2a / [e^{-2a}*(2+alpha+a*(3+alpha))/2 - a] /... (uprosczone)
    Wzor: K1_NLO = K1_LO * (1 + korekcja(alpha, a))
    Test dokladnosci: blad < ?%

  Czesc C -- Precyzyjne r21 z K1_NLO
    r21_precise = K2_NLO / K1_NLO
    Cel: blad max < 5%

  Czesc D -- Kompletny krajobraz fermionow
    Dla wszystkich 3 rodzin + a_Gamma in [0.005, 0.060]
    Q_TGP dla kazdej rodziny przy jej (alpha_f, a_Gamma)
"""

import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.special import exp1
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX  = 60.0
GAMMA  = 1.0
LAM_REF = 1e-5

print("P37: Precyzyjna formula K1(alpha, a_Gamma) i krajobraz fermionow")
print("=" * 70)
print()

# ============================================================
# NARZEDZIA (jak w P36)
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam=LAM_REF, N=2000):
    t   = np.linspace(0, 1, N)
    r   = a_gam * (R_MAX/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam=LAM_REF):
    e = energy_log(K, alpha, a_gam, lam)
    return e / (4*np.pi*K) - 1.0 if np.isfinite(e) else np.nan

def find_K1_num(alpha, a_gam):
    """Numeryczne K1 (pierwsze zero g(K) przy lambda_ref)."""
    K_lo, K_hi = 0.0005, 0.08
    try:
        glo = g_func(K_lo, alpha, a_gam)
        ghi = g_func(K_hi, alpha, a_gam)
        if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
            return brentq(lambda K: g_func(K, alpha, a_gam),
                          K_lo, K_hi, xtol=1e-13)
    except Exception:
        pass
    return np.nan

def find_K2_num(alpha, a_gam):
    """Numeryczne K2 (drugie zero g(K))."""
    for K_lo, K_hi in [(0.05, 0.5), (0.5, 5.0)]:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam),
                              K_lo, K_hi, xtol=1e-10)
        except Exception:
            pass
    return np.nan

# ============================================================
# FORMULY ANALITYCZNE
# ============================================================

def K1_LO(alpha, a_gam):
    """K1 ze wskaznika LO: K1 = 2a / e^{-2a}*[(1+alpha)*(2+a)/2 - a/2]
    Wyprowadzony z zerowania czlonu K w g(K) ~ K*Delta0/2 - 1 = 0.
    Delta0 = (1+alpha)*I_kin_std - I2_u
    I_kin_std = e^{-2a}*(2+a)/(2a)
    I2_u = e^{-2a}/2
    """
    a = a_gam
    e2a = np.exp(-2*a)
    I_kin_std = e2a * (2+a) / (2*a)
    I2_u      = e2a / 2
    Delta0 = (1+alpha)*I_kin_std - I2_u
    if Delta0 <= 0:
        return np.nan
    return 2.0 / Delta0

# Sprawdz K1_LO dla referncyjnych punktow
print("Test K1_LO vs K1_analytic vs K1_num przy a=0.040:")
print(f"  {'alpha':>6}  {'K1_LO':>10}  {'K1_2.351':>10}  {'K1_num':>10}  "
      f"{'err_LO%':>8}  {'C_K_eff':>8}")
print("  " + "-"*65)

C_K_fixed = 2.351
for alpha_test in [3.0, 5.0, 7.0, 8.553, 10.0, 12.0, 15.0]:
    a = 0.040
    K1_lo = K1_LO(alpha_test, a)
    K1_an = C_K_fixed * a / (1 + alpha_test)
    K1_n  = find_K1_num(alpha_test, a)
    if np.isnan(K1_n): continue
    err_lo = (K1_lo - K1_n)/K1_n*100
    C_K_eff = K1_n*(1+alpha_test)/a
    print(f"  {alpha_test:>6.3f}  {K1_lo:>10.6f}  {K1_an:>10.6f}  {K1_n:>10.6f}  "
          f"  {err_lo:>+7.2f}%  {C_K_eff:>8.5f}")
print()

# ============================================================
# CZESC A: Mapa C_K(alpha, a_Gamma)
# ============================================================
print("CZESC A: Mapa C_K(alpha, a_Gamma) = K1_num*(1+alpha)/a")
print("-" * 70)
print()

alpha_arr = np.array([1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.553,
                       9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0, 20.0])
agam_arr  = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060])

data_A = []
print(f"  {'alpha':>6}  {'a_gam':>6}  {'K1_num':>10}  {'K1_LO':>10}  "
      f"{'err_LO%':>8}  {'C_K_eff':>9}")
print("  " + "-"*58)

for alpha in alpha_arr:
    for a in agam_arr:
        K1n = find_K1_num(alpha, a)
        if np.isnan(K1n): continue
        K1lo = K1_LO(alpha, a)
        err_lo = (K1lo - K1n)/K1n*100 if not np.isnan(K1lo) else np.nan
        C_K = K1n*(1+alpha)/a
        data_A.append({'alpha':alpha, 'a':a, 'K1n':K1n, 'K1lo':K1lo,
                       'err_lo':err_lo, 'C_K':C_K})
        # Drukuj co 3. a_gam dla czytelnosci
        if abs(a - 0.040) < 0.001 or abs(a - 0.020) < 0.001:
            print(f"  {alpha:>6.3f}  {a:>6.3f}  {K1n:>10.7f}  {K1lo:>10.7f}  "
                  f"  {err_lo:>+7.2f}%  {C_K:>9.5f}")

print()
# Statystyki bledu K1_LO
errs_lo = [d['err_lo'] for d in data_A if not np.isnan(d['err_lo'])]
print(f"  K1_LO: blad sred = {np.mean(errs_lo):+.2f}%, max_abs = {np.max(np.abs(errs_lo)):.2f}%")

# Zakres C_K
CK_vals = [d['C_K'] for d in data_A]
print(f"  C_K(alpha,a): zakres [{min(CK_vals):.4f}, {max(CK_vals):.4f}], "
      f"srednia = {np.mean(CK_vals):.4f}")
print()

# ============================================================
# CZESC B: Analityczne NLO dla K1
# ============================================================
print("CZESC B: Analityczne wzory na C_K(alpha, a_Gamma)")
print("-" * 70)
print()

# Obserwacja: C_K zalezy glownie od alpha, slabiej od a
# Sprawdz separowalnosc: czy C_K(alpha, a) = f(alpha) * h(a)?
alpha_u = np.unique([d['alpha'] for d in data_A])
a_u     = np.unique([d['a'] for d in data_A])

print("  C_K(alpha, a) -- tabela wartosci:")
header = "  alpha\\a_gam " + "".join([f"  {a:6.3f}" for a in a_u])
print(header)
print("  " + "-"*70)
for alpha in alpha_u:
    row = [d['C_K'] for d in data_A if abs(d['alpha']-alpha)<0.001]
    if len(row) > 0:
        row_str = "".join([f"  {v:6.4f}" for v in row])
        print(f"  {alpha:>5.2f}  " + row_str)
print()

# Analityczne NLO: K1 = 2/Delta0 * correction
# Correction pochodzi z terminu -K^2 * B2 w rozwinieciu g(K):
# B2 = alpha*L1 + 4*I3_u/3
# gdzie L1 = int_a^inf e^{-3r}*(1+r)^2/r^3 dr
# I3_u = e^{-3a}/3

def I3_u_exact(a):
    """int_a^inf e^{-3r} dr / ... NO: int e^{-3r}*r^2/r^2 dr = int e^{-3r}dr = e^{-3a}/3"""
    # Wait: u = e^{-r}/r, u^3 = e^{-3r}/r^3
    # I3_u = int_a^inf u^3 * r^2 dr = int_a^inf e^{-3r}/r dr = E1(3a)
    return exp1(3*a)

def L1_exact(a):
    """L1 = int_a^inf (u')^2 * u * r^2 dr
    = int_a^inf e^{-2r}*(1+r)^2/r^4 * e^{-r}/r * r^2 dr
    = int_a^inf e^{-3r}*(1+r)^2/r^3 dr
    = int_a^inf [e^{-3r}/r^3 + 2e^{-3r}/r^2 + e^{-3r}/r] dr
    """
    # int e^{-3r}/r dr = E1(3a)
    # int e^{-3r}/r^2 dr = e^{-3a}/a - 3*E1(3a)
    # int e^{-3r}/r^3 dr: IBP: = e^{-3a}/(2a^2) - (3/2)*[e^{-3a}/a - 3*E1(3a)]
    #                         = e^{-3a}/(2a^2) - 3*e^{-3a}/(2a) + (9/2)*E1(3a)
    e3a = np.exp(-3*a)
    E1_3a = exp1(3*a)
    part_r3 = e3a/(2*a**2) - 3*e3a/(2*a) + (9/2)*E1_3a
    part_r2 = 2*(e3a/a - 3*E1_3a)
    part_r1 = E1_3a
    return part_r3 + part_r2 + part_r1

def K1_NLO_analytic(alpha, a_gam):
    """K1 NLO z korekcja drugiego rzedu."""
    a = a_gam
    e2a = np.exp(-2*a)
    I_kin = e2a * (2+a) / (2*a)
    I2u   = e2a / 2
    Delta0 = (1+alpha)*I_kin - I2u
    if Delta0 <= 0:
        return np.nan
    K1_lo = 2.0 / Delta0
    # NLO correction: K1 *= (1 + K1_lo * B2 / Delta0)
    L1  = L1_exact(a)
    I3u = I3_u_exact(a)
    B2  = alpha * L1 + (4.0/3.0) * I3u
    correction = 1.0 + K1_lo * B2 / Delta0
    return K1_lo * correction

# Test K1_NLO
print("  Test K1_NLO_analytic vs K1_num:")
print(f"  {'alpha':>6}  {'a':>6}  {'K1_num':>10}  {'K1_LO':>10}  "
      f"{'K1_NLO':>10}  {'err_LO%':>8}  {'err_NLO%':>9}")
print("  " + "-"*72)

err_lo_all  = []
err_nlo_all = []
data_B = []
for d in data_A:
    al, a = d['alpha'], d['a']
    K1n  = d['K1n']
    K1lo = d['K1lo']
    K1nlo = K1_NLO_analytic(al, a)
    err_lo  = (K1lo - K1n)/K1n*100
    err_nlo = (K1nlo - K1n)/K1n*100 if not np.isnan(K1nlo) else np.nan
    err_lo_all.append(err_lo)
    if not np.isnan(err_nlo):
        err_nlo_all.append(err_nlo)
    data_B.append({**d, 'K1nlo':K1nlo, 'err_nlo':err_nlo})
    if abs(a - 0.040) < 0.001:
        print(f"  {al:>6.3f}  {a:>6.3f}  {K1n:>10.7f}  {K1lo:>10.7f}  "
              f"{K1nlo:>10.7f}  {err_lo:>+7.2f}%  {err_nlo:>+8.2f}%")

print()
print(f"  K1_LO:  blad sred = {np.mean(err_lo_all):+.2f}%, "
      f"max_abs = {np.max(np.abs(err_lo_all)):.2f}%")
print(f"  K1_NLO: blad sred = {np.mean(err_nlo_all):+.2f}%, "
      f"max_abs = {np.max(np.abs(err_nlo_all)):.2f}%")
print()

# Jesli NLO nie jest wystarczajacy, dodaj empiryczną korekcję C_K_fit(alpha,a)
print("  Empiryczne dopasowanie C_K(alpha, a):")

alpha_fit = np.array([d['alpha'] for d in data_A])
a_fit     = np.array([d['a']     for d in data_A])
CK_fit    = np.array([d['C_K']   for d in data_A])

# Model 1: C_K = p0 + p1/alpha (monotonicznie malejace z alpha)
def model_CK_inv(X, p0, p1):
    alpha, a = X
    return p0 + p1/alpha

# Model 2: C_K = p0 + p1/alpha + p2*a (z poprawka a)
def model_CK_inv_a(X, p0, p1, p2):
    alpha, a = X
    return p0 + p1/alpha + p2*a

# Model 3: C_K z analitycznej LO (uzywamy K1_LO) = (2/Delta0)*(1+alpha)/a
def CK_LO_formula(X):
    alpha, a = X
    e2a = np.exp(-2*a)
    I_kin = e2a*(2+a)/(2*a)
    I2u   = e2a/2
    D0 = (1+alpha)*I_kin - I2u
    return 2*(1+alpha)/(a*D0)

# Model 4: C_K = p0 + p1*exp(-p2*alpha) (eksponencjalne)
def model_CK_exp(X, p0, p1, p2):
    alpha, a = X
    return p0 + p1*np.exp(-p2*alpha)

# Model 5: C_K = p0 + p1/alpha + p2/(alpha^2) (rozwinięcie Laurent)
def model_CK_laurent(X, p0, p1, p2):
    alpha, a = X
    return p0 + p1/alpha + p2/alpha**2

fits_CK = {}
for name, func, p0_arr in [
    ("inv_alpha",    model_CK_inv,     [2.2, 1.0]),
    ("inv_alpha_a",  model_CK_inv_a,   [2.2, 1.0, 0.5]),
    ("exp_alpha",    model_CK_exp,     [2.2, 1.5, 0.25]),
    ("laurent",      model_CK_laurent, [2.2, 1.0, 0.5]),
]:
    try:
        popt, _ = curve_fit(func, (alpha_fit, a_fit), CK_fit, p0=p0_arr, maxfev=10000)
        CK_pred = func((alpha_fit, a_fit), *popt)
        rmse = np.sqrt(np.mean((CK_pred - CK_fit)**2))
        maxe = np.max(np.abs(CK_pred - CK_fit))
        fits_CK[name] = (popt, rmse, maxe)
        print(f"  {name:>14}: params={[f'{p:.5f}' for p in popt]}  "
              f"RMSE={rmse:.5f}  maxErr={maxe:.5f}")
    except Exception as ex:
        print(f"  {name:>14}: FAIL ({ex})")

# Sprawdz C_K_LO (analityczne)
CK_LO_vals = CK_LO_formula((alpha_fit, a_fit))
rmse_lo = np.sqrt(np.mean((CK_LO_vals - CK_fit)**2))
maxe_lo = np.max(np.abs(CK_LO_vals - CK_fit))
print(f"  {'CK_LO_anal':>14}: (bez dopasowania)  "
      f"RMSE={rmse_lo:.5f}  maxErr={maxe_lo:.5f}")
print()

best_CK = min(fits_CK, key=lambda k: fits_CK[k][1])
bp_CK = fits_CK[best_CK][0]
print(f"  Najlepszy model: {best_CK}  RMSE={fits_CK[best_CK][1]:.5f}")

def K1_fit(alpha, a_gam):
    """K1 z najlepszego dopasowania C_K."""
    if best_CK == "inv_alpha":
        CK = model_CK_inv((alpha, a_gam), *bp_CK)
    elif best_CK == "inv_alpha_a":
        CK = model_CK_inv_a((alpha, a_gam), *bp_CK)
    elif best_CK == "exp_alpha":
        CK = model_CK_exp((alpha, a_gam), *bp_CK)
    elif best_CK == "laurent":
        CK = model_CK_laurent((alpha, a_gam), *bp_CK)
    else:
        CK = 2.351
    return CK * a_gam / (1 + alpha)

# ============================================================
# CZESC C: Precyzyjne r21 z K1_NLO/fit + K2_NLO
# ============================================================
print()
print("CZESC C: r21_precise = K2_NLO / K1_improved")
print("-" * 70)
print()

# Stale z P36 (model bilinearny)
G_p = [0.01316, 2.62013, 0.10759]   # params dla model_bilin z P36

def K2_cubic_eps(alpha, a_gam):
    c   = a_gam * (2*alpha - 4)
    eps = (4*a_gam/3) * np.log(1/a_gam)
    try:
        return brentq(lambda K: K**3 - eps*K**2 - 2*K - c, 0.5, 4.0, xtol=1e-10)
    except Exception:
        return np.nan

def G_bilin(alpha, a_gam):
    return 1 + G_p[0]*alpha + G_p[1]*a_gam + G_p[2]*alpha*a_gam

# Siatka testowa (z P36)
alpha_test = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.553, 9.0, 10.0, 11.0, 12.0])
agam_test  = np.array([0.020, 0.030, 0.040, 0.050, 0.060])

print(f"  {'alpha':>6}  {'a':>6}  {'r21_num':>8}  {'r21_OLD':>8}  {'r21_NEW_NLO':>12}  "
      f"{'r21_NEW_fit':>12}  {'err_OLD%':>9}  {'err_NLO%':>9}  {'err_fit%':>9}")
print("  " + "-"*90)

err_old_all = []; err_nlo_all2 = []; err_fit_all = []
for alpha in alpha_test:
    for a in agam_test:
        K2n = find_K2_num(alpha, a)
        K1n = find_K1_num(alpha, a)
        if np.isnan(K2n) or np.isnan(K1n): continue

        r21_num = K2n / K1n

        # Stare: K2_eps / K1_analytic(2.351)
        K2e    = K2_cubic_eps(alpha, a)
        K1_old = 2.351 * a / (1 + alpha)
        r21_old = K2e / K1_old if not np.isnan(K2e) else np.nan

        # Nowe: K2_NLO / K1_NLO
        K2_NLO_val = K2e * G_bilin(alpha, a) if not np.isnan(K2e) else np.nan
        K1_nlo     = K1_NLO_analytic(alpha, a)
        r21_NLO    = K2_NLO_val / K1_nlo if not (np.isnan(K2_NLO_val) or np.isnan(K1_nlo)) else np.nan

        # Nowe: K2_NLO / K1_fit
        K1_f    = K1_fit(alpha, a)
        r21_fit = K2_NLO_val / K1_f if not (np.isnan(K2_NLO_val) or np.isnan(K1_f)) else np.nan

        err_old = (r21_old - r21_num)/r21_num*100
        err_nlo = (r21_NLO - r21_num)/r21_num*100 if not np.isnan(r21_NLO) else np.nan
        err_fit = (r21_fit - r21_num)/r21_num*100 if not np.isnan(r21_fit) else np.nan

        err_old_all.append(err_old)
        if not np.isnan(err_nlo): err_nlo_all2.append(err_nlo)
        if not np.isnan(err_fit): err_fit_all.append(err_fit)

        if abs(a - 0.040) < 0.001:
            print(f"  {alpha:>6.3f}  {a:>6.3f}  {r21_num:>8.2f}  {r21_old:>8.2f}  "
                  f"{r21_NLO:>12.2f}  {r21_fit:>12.2f}  "
                  f"{err_old:>+9.2f}%  {err_nlo:>+8.2f}%  {err_fit:>+8.2f}%")

print()
print(f"  r21_OLD (cubic-eps / K1_fixed):   "
      f"sred={np.mean(err_old_all):+.2f}%, max_abs={np.max(np.abs(err_old_all)):.2f}%")
print(f"  r21_NLO (K2_NLO / K1_NLO_anal):  "
      f"sred={np.mean(err_nlo_all2):+.2f}%, max_abs={np.max(np.abs(err_nlo_all2)):.2f}%")
print(f"  r21_fit (K2_NLO / K1_fit):        "
      f"sred={np.mean(err_fit_all):+.2f}%, max_abs={np.max(np.abs(err_fit_all)):.2f}%")
print()

# ============================================================
# CZESC D: Kompletny krajobraz fermionow
# ============================================================
print("CZESC D: Kompletny krajobraz fermionow TGP")
print("-" * 70)
print()

fermions = {
    'e/mu/tau': {'m1': 0.511,   'm2': 105.658,  'm3': 1776.86},
    'u/c/t':    {'m1': 2.16,    'm2': 1270.0,    'm3': 172690.0},
    'd/s/b':    {'m1': 4.67,    'm2': 93.4,      'm3': 4180.0},
}

def Q_koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (m1 + m2 + m3)

# Dane PDG
print("  Dane PDG:")
for name, f in fermions.items():
    r21 = f['m2']/f['m1']; r31 = f['m3']/f['m1']
    Q = Q_koide(f['m1'], f['m2'], f['m3'])
    print(f"    {name:>10}: r21={r21:.1f}, r31={r31:.0f}, Q={Q:.4f}")
print()

def find_alpha_for_r21_precise(target_r21, a_gam, alpha_lo=0.5, alpha_hi=40.0):
    """Znajdz alpha takie ze r21_num(alpha, a_Gamma) = target_r21."""
    def f(al):
        K1 = find_K1_num(al, a_gam)
        K2 = find_K2_num(al, a_gam)
        if np.isnan(K1) or np.isnan(K2) or K1 < 1e-12: return np.nan
        return K2/K1 - target_r21
    try:
        flo = f(alpha_lo); fhi = f(alpha_hi)
        if np.isnan(flo) or np.isnan(fhi) or flo*fhi >= 0: return np.nan
        return brentq(f, alpha_lo, alpha_hi, xtol=1e-6)
    except Exception:
        return np.nan

# Skan a_Gamma i wyznacz alpha_f dla kazdej rodziny
agam_scan = np.array([0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040])

print("  Wyznaczanie alpha_f(a_Gamma) dla kazdej rodziny:")
landscape = {name: {} for name in fermions}

for name, ferm in fermions.items():
    r21_obs = ferm['m2'] / ferm['m1']
    print(f"\n  === {name} (r21={r21_obs:.1f}) ===")
    print(f"  {'a_Gamma':>8}  {'alpha_f':>9}  {'r21_TGP':>9}  {'K1':>10}  "
          f"{'K2':>8}  {'alpha/a':>9}")
    print("  " + "-"*62)
    for a in agam_scan:
        af = find_alpha_for_r21_precise(r21_obs, a)
        if not np.isnan(af):
            K1f = find_K1_num(af, a)
            K2f = find_K2_num(af, a)
            r21_check = K2f/K1f if not (np.isnan(K1f) or np.isnan(K2f)) else np.nan
            ratio = af/a
            landscape[name][a] = af
            print(f"  {a:>8.4f}  {af:>9.4f}  {r21_check:>9.1f}  "
                  f"{K1f:>10.7f}  {K2f:>8.4f}  {ratio:>9.2f}")
        else:
            print(f"  {a:>8.4f}  {'FAIL':>9}")

# Streszczenie alpha/a_Gamma
print()
print("  Podsumowanie alpha_f / a_Gamma dla kazdej rodziny:")
print(f"  {'Rodzina':>10}  {'r21':>7}  {'alpha/a_mean':>13}  {'alpha/a_std':>12}  "
      f"{'r21/sqrt2':>10}")
print("  " + "-"*58)
for name, ferm in fermions.items():
    r21_obs = ferm['m2'] / ferm['m1']
    ratios = [af/a for a, af in landscape[name].items()]
    if ratios:
        print(f"  {name:>10}  {r21_obs:>7.1f}  {np.mean(ratios):>13.2f}  "
              f"{np.std(ratios):>12.2f}  {r21_obs/np.sqrt(2):>10.2f}")
    else:
        print(f"  {name:>10}  {r21_obs:>7.1f}  BRAK DANYCH")

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("P37: Precyzyjna formula K1 i krajobraz fermionow TGP", fontsize=13)

# Panel A: C_K(alpha) przy roznych a_Gamma
ax = axes[0, 0]
ax.set_title("C_K(alpha) = K1*(1+alpha)/a vs alpha")
colors = plt.cm.viridis(np.linspace(0, 0.85, len(a_u)))
for i, a_val in enumerate(a_u):
    subdata = [d for d in data_A if abs(d['a']-a_val) < 0.001]
    if not subdata: continue
    alp_s = [d['alpha'] for d in subdata]
    ck_s  = [d['C_K']   for d in subdata]
    ax.plot(alp_s, ck_s, 'o-', color=colors[i], label=f"a={a_val:.3f}", ms=4)
# Nanies C_K_LO (analityczne)
alp_dense = np.linspace(1.5, 20, 100)
ck_lo_line = [CK_LO_formula((al, 0.040)) for al in alp_dense]
ax.plot(alp_dense, ck_lo_line, 'k--', label="C_K_LO(a=0.04)", lw=2)
ax.axhline(2.351, color='r', ls=':', label="C_K=2.351 (P26)", lw=1.5)
ax.set_xlabel(r"$\alpha$"); ax.set_ylabel(r"$C_K$")
ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

# Panel B: Blad K1_LO i K1_NLO vs alpha
ax = axes[0, 1]
ax.set_title("Blad K1_LO i K1_NLO vs alpha (a=0.040)")
subB = [d for d in data_B if abs(d['a']-0.040) < 0.001]
alp_B = [d['alpha'] for d in subB]
err_lo_B  = [d['err_lo']  for d in subB]
err_nlo_B = [d['err_nlo'] for d in subB]
ax.plot(alp_B, err_lo_B,  'ro-', ms=5, label="K1_LO")
ax.plot(alp_B, err_nlo_B, 'bs-', ms=5, label="K1_NLO")
ax.axhline(0, color='k', lw=0.8)
ax.axhspan(-2, 2, alpha=0.15, color='green', label="±2%")
ax.set_xlabel(r"$\alpha$"); ax.set_ylabel("blad K1 [%]")
ax.legend(); ax.grid(True, alpha=0.3)

# Panel C: Blad r21 OLD vs NLO vs fit (a=0.040)
ax = axes[1, 0]
ax.set_title("Blad r21 vs alpha (a=0.040)")
err_old_C = []; err_nlo_C = []; err_fit_C = []; alp_C = []
for alpha in alpha_test:
    a = 0.040
    K2n = find_K2_num(alpha, a); K1n = find_K1_num(alpha, a)
    if np.isnan(K2n) or np.isnan(K1n): continue
    r21_num = K2n/K1n
    K2e = K2_cubic_eps(alpha, a)
    if np.isnan(K2e): continue
    K1_old = 2.351*a/(1+alpha)
    K2_NLO_v = K2e * G_bilin(alpha, a)
    K1_nlo = K1_NLO_analytic(alpha, a)
    K1_f = K1_fit(alpha, a)
    r21_old = K2e/K1_old
    r21_nlo = K2_NLO_v/K1_nlo
    r21_fit = K2_NLO_v/K1_f
    alp_C.append(alpha)
    err_old_C.append((r21_old-r21_num)/r21_num*100)
    err_nlo_C.append((r21_nlo-r21_num)/r21_num*100)
    err_fit_C.append((r21_fit-r21_num)/r21_num*100)
ax.plot(alp_C, err_old_C, 'ro-', ms=5, label="r21_old (eps/2.351)")
ax.plot(alp_C, err_nlo_C, 'g^-', ms=5, label="r21_NLO (K2_NLO/K1_NLO)")
ax.plot(alp_C, err_fit_C, 'bs-', ms=5, label="r21_fit (K2_NLO/K1_fit)")
ax.axhline(0, color='k', lw=0.8)
ax.axhspan(-5, 5, alpha=0.1, color='green', label="±5% cel")
ax.set_xlabel(r"$\alpha$"); ax.set_ylabel("blad r21 [%]")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel D: Krajobraz alpha_f(a_Gamma) dla kazdej rodziny
ax = axes[1, 1]
ax.set_title("Krajobraz fermionow: alpha_f(a_Gamma)")
colors_f = ['blue', 'red', 'green']
for i, (name, ldata) in enumerate(landscape.items()):
    if not ldata: continue
    a_vals = sorted(ldata.keys())
    af_vals = [ldata[a] for a in a_vals]
    ax.plot(a_vals, af_vals, 'o-', color=colors_f[i], ms=5, label=name)
ax.set_xlabel(r"$a_\Gamma$"); ax.set_ylabel(r"$\alpha_f$")
ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
outpath = "TGP/TGP_v1/scripts/advanced/p37_K1_precise.png"
plt.savefig(outpath, dpi=120, bbox_inches='tight')
print(f"\nWykres zapisany: {outpath}")

print()
print("=" * 70)
print("WYNIKI KLUCZOWE P37")
print("=" * 70)
print()
print(f"  A. C_K(alpha, a): zakres [{min(CK_vals):.4f}, {max(CK_vals):.4f}]")
print(f"     Glowna zaleznosc od alpha (monotonicznie malejace)")
print()
print(f"  B. Dokladnosc formul K1:")
print(f"     K1_LO (2/Delta0):   blad max = {np.max(np.abs(err_lo_all)):.1f}%")
print(f"     K1_NLO (analitycz): blad max = {np.max(np.abs(err_nlo_all)):.1f}%")
print()
print(f"  C. Dokladnosc r21:")
print(f"     r21_old (eps/2.351):         blad max = {np.max(np.abs(err_old_all)):.1f}%")
print(f"     r21_NLO (K2_NLO / K1_NLO):  blad max = {np.max(np.abs(err_nlo_all2)):.1f}%")
print(f"     r21_fit (K2_NLO / K1_fit):   blad max = {np.max(np.abs(err_fit_all)):.1f}%")
print()
print("  D. Krajobraz fermionow:")
for name, ldata in landscape.items():
    r21_obs = fermions[name]['m2']/fermions[name]['m1']
    found = sorted(ldata.items())
    if found:
        a0, af0 = found[0]
        print(f"     {name:>10} (r21={r21_obs:.1f}): "
              f"alpha_f(a={a0:.3f}) = {af0:.4f}")
    else:
        print(f"     {name:>10}: BRAK WYNIKOW")
