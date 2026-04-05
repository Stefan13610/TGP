"""
p35_C_of_a.py  (v2)
====================
CEL: Pelna krzywa C(a_gam) i wyprowadzenie poprawki analitycznej

NAPRAWA v2:
  P35v1 miala FAIL dla a < 0.060 -- przepelnienie energii dla bardzo duzego K
  (przy lambda=1e-7 i a=0.025: K3 ~ 200, phi(r=a) ~ 20000, E_sex -> overflow).
  ROZWIAZANIE: uzyj adaptacyjnego zakresu lambda na podstawie formuly analitycznej
  lambda_anal = (C_LO * a / K1 / r31_K)^2 i szukaj w [0.1*lam_anal, 10*lam_anal].

KONTEKST (P31-P34):
  K3 ~ C * a_gam / sqrt(lambda),  C ~ 2.000 (empirycznie dla a=0.04)
  C jest NIEZALEZNE od alpha (P34), C = C(a_gam) TYLKO

  Znane wartosci (P33-P34):
    a=0.001: C ~ 2.121 (granica LO)
    a=0.025: C = 2.015
    a=0.030: C = 2.009
    a=0.040: C = 2.002
    a=0.050: C = 1.999  (minimum?)
    a=0.060: C = 2.000

PLAN P35v2:
  Czesc A -- Skan C(a_gam) [0.003, 0.120] z adaptacyjnym lambda
    Krok 1: K1 = find_K1(alpha, a) (niezalezne od lambda)
    Krok 2: r21, r31_K -- algebraicznie
    Krok 3: lambda_anal = (C_LO * a / K1 / r31_K)^2 -- punkt startowy
    Krok 4: bisekcja w [0.1*lam_anal, 5*lam_anal] --> lambda_K_num
    Krok 5: K3 = find_K3(alpha, a, lambda_K_num) --> C_num = K3*sqrt(lam)/a

  Czesc B -- C_qs vs C_num: rozklad na wklady
  Czesc C -- Dopasowanie analityczne C(a)
  Czesc D -- Poprawa dokladnosci lambda_Koide
"""

import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.special import exp1
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX = 60.0
GAMMA = 1.0
ALPHA_FIX = 8.553
C_LO = np.sqrt(4.5)

print("P35v2: Pelna krzywa C(a_gam) i analityczna poprawka")
print("=" * 70)
print()

# ============================================================
# NARZEDZIA
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=2000):
    t   = np.linspace(0, 1, N)
    r   = a_gam * (R_MAX/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1  = V_mod(1.0, lam)
    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep  = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam):
    e = energy_log(K, alpha, a_gam, lam)
    if not np.isfinite(e):
        return np.nan
    return e / (4*np.pi*K) - 1.0

def find_zero_K(alpha, a_gam, lam, K_lo, K_hi):
    try:
        glo = g_func(K_lo, alpha, a_gam, lam)
        ghi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(glo) and np.isfinite(ghi)):
            return np.nan
        if glo*ghi >= 0:
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam), K_lo, K_hi,
                      xtol=1e-10, maxiter=80)
    except Exception:
        return np.nan

def find_all_zeros(alpha, a_gam, lam):
    """Znajdz K1, K2, K3 z adaptacyjnymi przedzialami -- unika overflow."""
    K3_est = max(2.1, 2.2 * a_gam / np.sqrt(lam))  # <= unikamy eksplodujacych K
    # Przedzialy dobrze pokrywajace wszystkie trzy zera
    intervals = [
        (0.001, 0.05),
        (0.05,  0.5),
        (0.5,   max(3.0, K3_est*0.15)),
        (max(2.0, K3_est*0.4), max(5.0, K3_est*0.7)),
        (max(3.0, K3_est*0.6), max(20.0, K3_est*1.8)),
    ]
    zeros = []
    for K_lo, K_hi in intervals:
        if K_lo >= K_hi:
            continue
        K = find_zero_K(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            # Deduplikacja
            if not zeros or abs(K - zeros[-1]) / max(zeros[-1], 1e-6) > 0.02:
                zeros.append(K)
    return sorted(zeros)

def Q_koide(K1, K2, K3):
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    d = K1 + K2 + K3
    return s**2 / d

def r31_koide_alg(r21):
    """r31 z warunku Q(1,r21,x)=3/2 -- algebraiczne."""
    def q3(x):
        s = 1.0 + np.sqrt(r21) + np.sqrt(x)
        d = 1.0 + r21 + x
        return s**2/d - 1.5
    try:
        lo = r21 * 1.001
        hi = r21 * 30.0
        if q3(lo)*q3(hi) >= 0:
            hi = r21 * 200.0
        return brentq(q3, lo, hi, xtol=1e-8)
    except Exception:
        return np.nan

# Analityczne calki I4, I6 (z P33)
def I4_exact(a):
    return np.exp(-4*a)/a - 4*exp1(4*a)

def I6_exact(a):
    I4b = np.exp(-6*a)/a - 6*exp1(6*a)
    I5  = np.exp(-6*a)/(2*a**2) + 3*I4b
    return np.exp(-6*a)/(3*a**3) + 2*I5

def C_quartic_sextic(a):
    """C z bilansu quartic=sextic (formuła P33): C = sqrt(3*I4/(2*I6)) / a"""
    return np.sqrt(3*I4_exact(a) / (2*I6_exact(a))) / a

# ============================================================
# CZESC A: Skan C(a_gam) -- adaptacyjna bisekcja po lambda
# ============================================================
print("CZESC A: Skan C(a_gam) w zakresie [0.005, 0.120]")
print("-" * 70)
print()

a_arr = np.array([0.005, 0.008, 0.010, 0.015, 0.020, 0.025, 0.030,
                  0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.070,
                  0.080, 0.090, 0.100, 0.110, 0.120])

alpha = ALPHA_FIX
print(f"  alpha = {alpha} (stale)")
print(f"  Liczba punktow a_gam: {len(a_arr)}")
print()
print(f"  {'a_gam':>7}  {'lam_K_num':>12}  {'K3':>8}  {'C_num':>8}  "
      f"{'C_qs':>8}  {'DeltaC':>8}")
print("  " + "-"*65)

C_num_arr = []
C_qs_arr  = []
lam_K_arr = []
K3_arr    = []
K1_arr    = []
K2_arr    = []

# Staly lambda_ref do wyznaczenia K1, K2 (niezalezne od lambda)
LAM_REF = 1e-5

for a in a_arr:
    # Krok 1: K1, K2 przy lambda_ref (niezalezne od lambda)
    zeros_ref = find_all_zeros(alpha, a, LAM_REF)
    if len(zeros_ref) < 2:
        for lst in [C_num_arr, C_qs_arr, lam_K_arr, K3_arr, K1_arr, K2_arr]:
            lst.append(np.nan)
        print(f"  {a:>7.4f}  {'FAIL(K1K2)':>12}")
        continue
    K1_ref, K2_ref = zeros_ref[0], zeros_ref[1]

    # Krok 2: r21, r31_K -- algebraicznie
    r21 = K2_ref / K1_ref
    r31_K = r31_koide_alg(r21)
    if np.isnan(r31_K):
        for lst in [C_num_arr, C_qs_arr, lam_K_arr, K3_arr, K1_arr, K2_arr]:
            lst.append(np.nan)
        print(f"  {a:>7.4f}  {'FAIL(r31K)':>12}")
        continue

    # Krok 3: punkt startowy lambda z formuła analityczną
    # r31_K = K3/K1 => K3 = K1_ref * r31_K
    # K3 ~ C*a/sqrt(lam) => lam_anal = (C*a/(K1_ref*r31_K))^2
    lam_anal = (C_LO * a)**2 / (K1_ref * r31_K)**2

    # Krok 4: bisekcja w adaptacyjnym przedziale
    lam_lo = lam_anal * 0.05
    lam_hi = lam_anal * 15.0

    def f_r31(lam_try):
        zeros_try = find_all_zeros(alpha, a, lam_try)
        if len(zeros_try) < 3:
            return np.nan
        K3t = zeros_try[2]
        K1t = zeros_try[0]
        return K3t/K1t - r31_K

    try:
        flo = f_r31(lam_lo)
        fhi = f_r31(lam_hi)
        if np.isnan(flo) or np.isnan(fhi) or flo*fhi >= 0:
            # Rozszerz zakres
            lam_lo2 = lam_anal * 0.01
            lam_hi2 = lam_anal * 50.0
            flo = f_r31(lam_lo2); fhi = f_r31(lam_hi2)
            if np.isnan(flo) or np.isnan(fhi) or flo*fhi >= 0:
                raise ValueError("no sign change")
            lam_lo, lam_hi = lam_lo2, lam_hi2

        lam_K = brentq(f_r31, lam_lo, lam_hi,
                        xtol=lam_anal*1e-6, maxiter=80)
    except Exception as ex:
        for lst in [C_num_arr, C_qs_arr, lam_K_arr, K3_arr, K1_arr, K2_arr]:
            lst.append(np.nan)
        print(f"  {a:>7.4f}  {'FAIL(lam)':>12}  ({ex})")
        continue

    # Krok 5: K3 dokladne przy lambda_K
    zeros_K = find_all_zeros(alpha, a, lam_K)
    if len(zeros_K) < 3:
        for lst in [C_num_arr, C_qs_arr, lam_K_arr, K3_arr, K1_arr, K2_arr]:
            lst.append(np.nan)
        lam_K_arr[-1] = lam_K
        print(f"  {a:>7.4f}  {lam_K:>12.4e}  {'<3zeros':>8}")
        continue

    K1f, K2f, K3f = zeros_K[0], zeros_K[1], zeros_K[2]
    C_num = K3f * np.sqrt(lam_K) / a
    C_qs  = C_quartic_sextic(a)

    C_num_arr.append(C_num)
    C_qs_arr.append(C_qs)
    lam_K_arr.append(lam_K)
    K3_arr.append(K3f)
    K1_arr.append(K1f)
    K2_arr.append(K2f)

    print(f"  {a:>7.4f}  {lam_K:>12.4e}  {K3f:>8.4f}  {C_num:>8.5f}  "
          f"{C_qs:>8.5f}  {C_num-C_qs:>+8.5f}")

C_num_arr = np.array(C_num_arr)
C_qs_arr  = np.array(C_qs_arr)
lam_K_arr = np.array(lam_K_arr)
K3_arr    = np.array(K3_arr)
K1_arr    = np.array(K1_arr)
K2_arr    = np.array(K2_arr)

ok = np.isfinite(C_num_arr)
a_ok = a_arr[ok]
C_ok = C_num_arr[ok]
C_qs_ok = C_qs_arr[ok]
dC_ok = C_ok - C_qs_ok

print()
print(f"  C_LO = sqrt(4.5) = {C_LO:.5f}")
print(f"  C_num zakres: {np.nanmin(C_num_arr):.5f} -- {np.nanmax(C_num_arr):.5f}")
if ok.sum() > 0:
    imin = np.nanargmin(C_num_arr)
    print(f"  Minimum C_num: {np.nanmin(C_num_arr):.5f} przy a_gam={a_arr[imin]:.4f}")
print(f"  DeltaC zakres: {dC_ok.min():.5f} -- {dC_ok.max():.5f}")

# ============================================================
# CZESC B: Rozklad energii przy K3
# ============================================================
print()
print("CZESC B: Rozklad E(K3)/(4*pi*K3) na cztery czlony")
print("-" * 70)
print()

def energy_terms(K, alpha, a_gam, lam, N=2000):
    """Zwraca E_kin, E_cub, E_qua, E_sex / (4*pi*K) oddzielnie."""
    t   = np.linspace(0, 1, N)
    r   = a_gam * (R_MAX/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    fac = 4*np.pi*K

    Ek  = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r) / fac
    Ec  = 4*np.pi * np.trapezoid(GAMMA/3*(phi**3 - 1.0) * r**2, r) / fac
    Eq  = 4*np.pi * np.trapezoid(-GAMMA/4*(phi**4 - 1.0) * r**2, r) / fac
    Es  = 4*np.pi * np.trapezoid(lam/6*(phi-1.0)**6 * r**2, r) / fac
    return Ek, Ec, Eq, Es

print(f"  {'a_gam':>7}  {'K3':>7}  {'Ek':>9}  {'Ec':>9}  {'Eq':>12}  {'Es':>12}  {'sum':>7}")
print("  " + "-"*74)

for i, a in enumerate(a_arr):
    if not ok[i]:
        continue
    K3 = K3_arr[i]; lam = lam_K_arr[i]
    if np.isnan(K3): continue
    ek, ec, eq, es = energy_terms(K3, alpha, a, lam)
    s = ek + ec + eq + es
    print(f"  {a:>7.4f}  {K3:>7.3f}  {ek:>9.2f}  {ec:>9.2f}  "
          f"{eq:>12.2f}  {es:>12.2f}  {s:>7.4f}")

# ============================================================
# CZESC C: Dopasowanie analityczne C(a)
# ============================================================
print()
print("CZESC C: Dopasowanie analityczne C(a_gam)")
print("-" * 70)
print()

if len(a_ok) < 4:
    print("  Za malo punktow do dopasowania!")
else:
    def model_lin(a, p):
        return C_LO - p*a

    def model_exp(a, b):
        return C_LO * np.exp(-b*a)

    def model_pade01(a, beta):
        return C_LO / (1 + beta*a)

    def model_poly2(a, p, q):
        return C_LO - p*a + q*a**2

    def model_pade11(a, p, q):
        return C_LO * (1 - p*a) / (1 + q*a)

    fits = {}
    for name, func, p0 in [
        ("lin",    model_lin,    [3.0]),
        ("exp",    model_exp,    [1.5]),
        ("pade01", model_pade01, [1.5]),
        ("poly2",  model_poly2,  [3.0, 5.0]),
        ("pade11", model_pade11, [1.0, 1.0]),
    ]:
        try:
            popt, _ = curve_fit(func, a_ok, C_ok, p0=p0, maxfev=5000)
            C_fit = func(a_ok, *popt)
            rmse  = np.sqrt(np.mean((C_fit - C_ok)**2))
            maxe  = np.max(np.abs(C_fit - C_ok))
            fits[name] = (popt, rmse, maxe)
            print(f"  {name:>8}: params={[f'{p:.4f}' for p in popt]}  "
                  f"RMSE={rmse:.5f}  maxErr={maxe:.5f}")
        except Exception as e:
            print(f"  {name:>8}: FAIL ({e})")
            fits[name] = (None, 999, 999)

    best = min(fits, key=lambda k: fits[k][1])
    print(f"\n  Najlepszy model: {best}  (RMSE={fits[best][1]:.5f})")
    bp = fits[best][0]

    def apply_best(a_val):
        if   best == "lin":    return model_lin(a_val, *bp)
        elif best == "exp":    return model_exp(a_val, *bp)
        elif best == "pade01": return model_pade01(a_val, *bp)
        elif best == "poly2":  return model_poly2(a_val, *bp)
        elif best == "pade11": return model_pade11(a_val, *bp)
        return 2.000

    print()
    print(f"  Predykcja C(a) vs numeryk:")
    print(f"  {'a_gam':>7}  {'C_num':>8}  {'C_best':>8}  {'delta':>8}")
    for i, a in enumerate(a_arr):
        if not ok[i]: continue
        C_pred = apply_best(a)
        print(f"  {a:>7.4f}  {C_num_arr[i]:>8.5f}  {C_pred:>8.5f}  "
              f"{C_pred-C_num_arr[i]:>+8.5f}")

# ============================================================
# CZESC D: Porownanie dokladnosci formul
# ============================================================
print()
print("CZESC D: Blad formuly lambda_Koide: C=2.000 vs C=C_LO vs C=C(a)")
print("-" * 70)
print()

def lambda_K_formula(K1, a_gam, r31_K_val, C_val):
    return (C_val * a_gam)**2 / (K1 * r31_K_val)**2

print(f"  {'a_gam':>7}  {'lam_K_num':>12}  {'err(C=2)%':>10}  "
      f"{'err(C_LO)%':>11}  {'err(C(a))%':>11}")
print("  " + "-"*60)

err_C200 = []; err_CLO = []; err_Ca = []; a_val = []

for i, a in enumerate(a_arr):
    if not ok[i]: continue
    lam_true = lam_K_arr[i]
    K1 = K1_arr[i]
    if np.isnan(K1) or np.isnan(lam_true): continue

    # Oblicz r21, r31_K na podstawie K1, K2
    r21 = K2_arr[i] / K1
    r31_K = r31_koide_alg(r21)
    if np.isnan(r31_K): continue

    lam_200  = lambda_K_formula(K1, a, r31_K, 2.000)
    lam_CLO  = lambda_K_formula(K1, a, r31_K, C_LO)
    lam_Ca   = lambda_K_formula(K1, a, r31_K, apply_best(a))

    e200 = (lam_200 - lam_true)/lam_true*100
    eLO  = (lam_CLO - lam_true)/lam_true*100
    eCa  = (lam_Ca  - lam_true)/lam_true*100

    err_C200.append(e200); err_CLO.append(eLO); err_Ca.append(eCa)
    a_val.append(a)

    print(f"  {a:>7.4f}  {lam_true:>12.4e}  {e200:>+10.3f}  {eLO:>+11.3f}  {eCa:>+11.3f}")

err_C200 = np.array(err_C200)
err_CLO  = np.array(err_CLO)
err_Ca   = np.array(err_Ca)
a_val    = np.array(a_val)

print()
print(f"  Statystyki bledu lambda_K:")
print(f"    C=2.000:  sred={np.mean(err_C200):+.3f}%, max_abs={np.max(np.abs(err_C200)):.3f}%")
print(f"    C=C_LO:   sred={np.mean(err_CLO):+.3f}%, max_abs={np.max(np.abs(err_CLO)):.3f}%")
print(f"    C=C(a):   sred={np.mean(err_Ca):+.3f}%, max_abs={np.max(np.abs(err_Ca)):.3f}%")

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("P35: Pelna krzywa C(a_gam) -- analiza i dopasowanie", fontsize=13)

# Panel 1: C_num, C_qs, modele
ax = axes[0,0]
ax.plot(a_ok, C_ok, 'o-', color='navy', label='C_num (numeryk)', lw=2, ms=7)
ax.plot(a_ok, C_qs_ok, 's--', color='firebrick', label='C_qs (quartic-sextic)', lw=1.5, ms=5)
ax.axhline(C_LO, color='gray', ls=':', lw=1.5, label=f'C_LO=sqrt(4.5)={C_LO:.3f}')
ax.axhline(2.000, color='green', ls=':', lw=1.5, label='C=2.000')
if len(a_ok) >= 4:
    a_d = np.linspace(a_ok.min()*0.8, a_ok.max()*1.1, 200)
    clrs = {'lin':'orange','exp':'purple','pade01':'cyan','poly2':'lime','pade11':'magenta'}
    for nm, (popt, rmse, _) in fits.items():
        if popt is None: continue
        if   nm=="lin":    yf=model_lin(a_d,*popt)
        elif nm=="exp":    yf=model_exp(a_d,*popt)
        elif nm=="pade01": yf=model_pade01(a_d,*popt)
        elif nm=="poly2":  yf=model_poly2(a_d,*popt)
        elif nm=="pade11": yf=model_pade11(a_d,*popt)
        ls = '-' if nm==best else '--'
        ax.plot(a_d, yf, ls, color=clrs.get(nm,'gray'), lw=1.5 if nm==best else 0.8,
                label=f'{nm} (RMSE={rmse:.4f})')
ax.set_xlabel('a_gam')
ax.set_ylabel('C = K3*sqrt(lam)/a')
ax.set_title('C(a_gam): numeryk vs dopasowania')
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

# Panel 2: DeltaC = C_num - C_qs
ax = axes[0,1]
ax.plot(a_ok, dC_ok, 'D-', color='darkorange', lw=2, ms=7, label='DeltaC = C_num - C_qs')
ax.axhline(0, color='black', lw=0.8)
ax.set_xlabel('a_gam')
ax.set_ylabel('DeltaC')
ax.set_title('Wklad E_kin + E_cub: DeltaC = C_num - C_qs')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel 3: log scale
ax = axes[1,0]
ax.semilogx(a_ok, C_ok, 'o-', color='navy', lw=2, ms=7, label='C_num')
ax.semilogx(a_ok, C_qs_ok, 's--', color='firebrick', lw=1.5, ms=5, label='C_qs')
ax.axhline(C_LO, color='gray', ls=':', lw=1.5, label=f'sqrt(4.5)={C_LO:.3f}')
ax.axhline(2.000, color='green', ls=':', lw=1.5, label='C=2.000')
if len(a_ok) >= 4:
    a_dl = np.exp(np.linspace(np.log(a_ok.min()*0.8), np.log(a_ok.max()*1.1), 200))
    yf = apply_best(a_dl)
    ax.semilogx(a_dl, yf, 'r-', lw=2, label=f'fit-{best}')
ax.set_xlabel('a_gam (log)')
ax.set_ylabel('C')
ax.set_title('C(a_gam) w skali logarytmicznej')
ax.legend(fontsize=8)
ax.grid(alpha=0.3, which='both')

# Panel 4: Blad formuly
ax = axes[1,1]
if len(a_val) > 0:
    ax.plot(a_val, err_C200, 'o-', color='red', label='C=2.000', lw=2, ms=6)
    ax.plot(a_val, err_CLO,  's-', color='gray', label=f'C=sqrt(4.5)={C_LO:.3f}', lw=1.5, ms=5)
    ax.plot(a_val, err_Ca,   '^-', color='green', label=f'C=C(a) [{best}]', lw=2, ms=6)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xlabel('a_gam')
    ax.set_ylabel('blad lambda_K [%]')
    ax.set_title('Dokladnosc formuly lambda_Koide')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('p35_C_of_a.png', dpi=120, bbox_inches='tight')
print()
print("Wykres zapisany: p35_C_of_a.png")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 70)
print("WYNIK KLUCZOWY P35: C(a_gam)")
print("=" * 70)
print()
print(f"  C_LO = sqrt(4.5) = {C_LO:.5f}  (teoria: granica a->0)")
print(f"  C_num min = {np.nanmin(C_num_arr):.5f} przy a = {a_arr[np.nanargmin(C_num_arr)]:.4f}")
print(f"  C_num dla a=0.040 = {C_num_arr[np.argmin(np.abs(a_arr-0.040))]:.5f}")
print()
if len(a_ok) >= 4:
    print(f"  Najlepsza formula analityczna: {best}")
    if   best=="lin":    print(f"    C(a) = {C_LO:.5f} - {bp[0]:.4f}*a")
    elif best=="exp":    print(f"    C(a) = {C_LO:.5f} * exp(-{bp[0]:.4f}*a)")
    elif best=="pade01": print(f"    C(a) = {C_LO:.5f} / (1 + {bp[0]:.4f}*a)")
    elif best=="poly2":  print(f"    C(a) = {C_LO:.5f} - {bp[0]:.4f}*a + {bp[1]:.4f}*a^2")
    elif best=="pade11": print(f"    C(a) = {C_LO:.5f} * (1-{bp[0]:.4f}*a)/(1+{bp[1]:.4f}*a)")
    print()
    print(f"  Poprawa dokladnosci lambda_Koide:")
    print(f"    C=2.000:  sred={np.mean(err_C200):+.3f}%, max_abs={np.max(np.abs(err_C200)):.3f}%")
    print(f"    C=C_LO:   sred={np.mean(err_CLO):+.3f}%, max_abs={np.max(np.abs(err_CLO)):.3f}%")
    print(f"    C=C(a):   sred={np.mean(err_Ca):+.3f}%, max_abs={np.max(np.abs(err_Ca)):.3f}%")
print()
print(f"  DeltaC = C_num - C_qs (wklad kin+cub): {dC_ok.min():.4f} do {dC_ok.max():.4f}")
print(f"  Interpretacja: DeltaC > 0 zawsze => wklad E_kin+E_cub przesuwa C w gore")
print(f"  C_qs (quartic-sextic only): {np.nanmin(C_qs_arr):.4f} -- {np.nanmax(C_qs_arr):.4f}")
print(f"  Stosunek C_num/C_qs ~ {np.mean(C_ok/C_qs_ok):.4f} (C_num zawsze > C_qs)")
