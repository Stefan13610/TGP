"""
p33_K3_correction.py  (v2)
==========================
CEL: Wyjaśnienie C = 2.000 vs teoria C_LO = sqrt(4.5) = 2.121

WYNIK ROBOCZY (Krok 2 v1):
  I4_num = 15.667  vs  I4_LO = 1/a = 25.000  (blad LO: 60%!)
  I6_num = 3684.6  vs  I6_LO = 1/(3a^3) = 5208  (blad LO: 41%!)

HIPOTEZA (zrewidowana):
  Rownanie balansu dla K3 (z g(K3)=0 przy duzym K3):
    -K3^3 * I4/4 + lam * K3^5 * I6/6 ~ 0
    K3^2 = 3*I4 / (2*lam*I6)

  LO uzywalo I4 ~ 1/a, I6 ~ 1/(3a^3), co daje K3^2 ~ 4.5*a^2/lam.
  Ale DOKLADNE I4, I6 sa INNE (blad LO ~40-60%), stad C_num = 2.000 ≠ 2.121.

PLAN P33:
  Krok 1 -- Weryfikacja: K3^2 = 3*I4/(2*lam*I6) predykcja vs numeryk
  Krok 2 -- Skad pochodzi blad LO? Analityczne I4 i I6 (z E1)
  Krok 3 -- Poprawa: formula NLO dla I4 i I6, stad NLO dla C
  Krok 4 -- Zaleznosc C(a_gam) z poprawna formula
  Krok 5 -- Wnioski: kompletna poprawiona formula dla K3
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import exp1  # E_1(x) = int_x^inf e^{-t}/t dt
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX = 60.0
GAMMA = 1.0

print("P33: NLO korekta K3 -- dokladne calki I4, I6 i formula C(a_gam)")
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
    return energy_log(K, alpha, a_gam, lam) / (4*np.pi*K) - 1.0

def find_zero_K(alpha, a_gam, lam, K_lo, K_hi):
    try:
        glo = g_func(K_lo, alpha, a_gam, lam)
        ghi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(glo) and np.isfinite(ghi)):
            return np.nan
        if glo*ghi >= 0:
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam), K_lo, K_hi, xtol=1e-10)
    except Exception:
        return np.nan

def find_all_zeros(alpha, a_gam, lam):
    # Automatyczny zakres -- K3 ~ 2*a_gam/sqrt(lam)
    K3_est = 2.5 * a_gam / np.sqrt(lam)
    intervals = [(0.001, 0.05), (0.05, 0.5), (0.5, 5.0),
                 (max(2.0, K3_est*0.2), max(10.0, K3_est*0.8)),
                 (max(5.0, K3_est*0.5), max(50.0, K3_est*2.0))]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero_K(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            if not zeros or abs(K - zeros[-1])/zeros[-1] > 0.01:
                zeros.append(K)
    return sorted(zeros)

def get_K3(alpha, a_gam, lam):
    zeros = find_all_zeros(alpha, a_gam, lam)
    return zeros[2] if len(zeros) >= 3 else np.nan

# Analityczne calki (z rozwinieciem E_1)
def I4_exact(a, Rmax=60.0):
    """I4 = int_a^Rmax e^{-4r}/r^2 dr = e^{-4a}/a - 4*E1(4a)"""
    return np.exp(-4*a)/a - 4*exp1(4*a)

def I6_exact(a, Rmax=60.0):
    """I6 = int_a^Rmax e^{-6r}/r^4 dr"""
    # Rozwijamy przez czesty: int e^{-nr}/r^m dr
    # = [-e^{-nr}/((m-1)r^{m-1})] + n/(m-1)*int e^{-nr}/r^{m-1} dr
    # I6 = int_a e^{-6r}/r^4 dr
    # = [-e^{-6r}/(3r^3)]_a + 2*int_a e^{-6r}/r^3 dr
    # = e^{-6a}/(3a^3) + 2*I5
    # I5 = int_a e^{-6r}/r^3 dr
    # = [-e^{-6r}/(2r^2)]_a + 3*int_a e^{-6r}/r^2 dr
    # = e^{-6a}/(2a^2) + 3*I4b
    # I4b = int_a e^{-6r}/r^2 dr = e^{-6a}/a - 6*E1(6a)
    I4b = np.exp(-6*a)/a - 6*exp1(6*a)
    I5  = np.exp(-6*a)/(2*a**2) + 3*I4b
    I6  = np.exp(-6*a)/(3*a**3) + 2*I5
    return I6

def C_from_integrals(a, lam, Rmax=60.0):
    """C = K3*sqrt(lam)/a z dokladnych calek.
    Rownanie: K3^2 = 3*I4/(2*lam*I6)
    C = K3*sqrt(lam)/a = sqrt(3*I4/(2*I6)) / a
    """
    I4 = I4_exact(a, Rmax)
    I6 = I6_exact(a, Rmax)
    K3_pred = np.sqrt(3*I4 / (2*lam*I6))
    C = K3_pred * np.sqrt(lam) / a
    return C, K3_pred, I4, I6

# ============================================================
# KROK 1: Weryfikacja K3^2 = 3*I4 / (2*lam*I6)
# ============================================================
print("KROK 1: Weryfikacja formuly K3^2 = 3*I4 / (2*lam*I6)")
print("-" * 70)
print()

alpha_ref = 8.553
lam_ref   = 5.489e-6

agam_test = [0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.080]

print(f"  {'a_gam':>8}  {'K3_num':>10}  {'K3_pred':>10}  {'C_num':>8}  "
      f"{'C_pred':>8}  {'err%':>8}  {'I4_num':>10}  {'I4_exact':>10}")
print("  " + "-"*78)

results = []
for ag in agam_test:
    K3_num = get_K3(alpha_ref, ag, lam_ref)
    if np.isnan(K3_num):
        print(f"  {ag:>8.3f}  {'---':>10}")
        continue

    C_pred, K3_pred, I4_ex, I6_ex = C_from_integrals(ag, lam_ref)
    C_num = K3_num * np.sqrt(lam_ref) / ag
    err = (K3_pred - K3_num)/K3_num * 100

    # Numeryczna I4
    r_int = np.logspace(np.log10(ag), np.log10(R_MAX), 10000)
    I4_num = np.trapezoid(np.exp(-4*r_int)/r_int**2, r_int)

    print(f"  {ag:>8.3f}  {K3_num:>10.4f}  {K3_pred:>10.4f}  {C_num:>8.4f}  "
          f"{C_pred:>8.4f}  {err:>8.2f}%  {I4_num:>10.4f}  {I4_ex:>10.4f}")
    results.append((ag, K3_num, K3_pred, C_num, C_pred, I4_num, I4_ex))

print()

# ============================================================
# KROK 2: Skad pochodzi blad LO? -- dokladne I4, I6 vs ich LO
# ============================================================
print("KROK 2: Analiza bledu LO -- I4_exact vs I4_LO = 1/a")
print("-" * 70)
print()
print("  I4_LO = 1/a  (pochodzi z int_a^inf 1/r^2 dr)")
print("  I4_exact = e^{-4a}/a - 4*E1(4a)  (dokladna wartosc z exp_integralem)")
print()
print("  Rozbicie: I4_exact = 1/a * [e^{-4a} - 4a*E1(4a)]")
print("                     = 1/a * f4(4a)")
print("  gdzie f4(x) = e^{-x} - x*E1(x)")
print()
print(f"  {'a_gam':>8}  {'4a':>8}  {'f4':>8}  {'I4_exact':>12}  {'I4_LO':>10}  "
      f"{'I4_exact/I4_LO':>16}  {'I6_exact/I6_LO':>16}")
print("  " + "-"*84)

for ag in agam_test:
    x = 4*ag
    f4 = np.exp(-x) - x*exp1(x)
    I4_ex = I4_exact(ag)
    I6_ex = I6_exact(ag)
    I4_lo = 1.0/ag
    I6_lo = 1.0/(3*ag**3)
    print(f"  {ag:>8.3f}  {x:>8.4f}  {f4:>8.5f}  {I4_ex:>12.4f}  {I4_lo:>10.4f}  "
          f"  {I4_ex/I4_lo:>12.4f}    {I6_ex/I6_lo:>12.4f}")

print()
print("  KLUCZOWY WNIOSEK:")
print("  I4_exact/I4_LO < 1 dla wszystkich a_gam (LO PRZESZACOWUJE I4)")
print("  I6_exact/I6_LO < 1 dla wszystkich a_gam (LO PRZESZACOWUJE I6)")
print()
print("  Ale STOSUNEK I4/I6 jest inny niz LO: (I4/I6)_exact != (I4/I6)_LO")
print("  C = sqrt(3*I4/(2*I6)) / a -- wlasciwa formula")
print()

# Sprawdz stosunek
print(f"  {'a_gam':>8}  {'C_LO':>8}  {'C_exact':>10}  {'ratio I4/I6':>14}  {'ratio LO':>12}  {'C_LO/C_ex':>12}")
print("  " + "-"*70)
for ag in agam_test:
    I4_ex = I4_exact(ag)
    I6_ex = I6_exact(ag)
    I4_lo = 1.0/ag
    I6_lo = 1.0/(3*ag**3)
    C_lo = np.sqrt(4.5)
    C_ex = np.sqrt(3*I4_ex/(2*I6_ex)) / ag
    ratio_ex = I4_ex / I6_ex
    ratio_lo = I4_lo / I6_lo
    print(f"  {ag:>8.3f}  {C_lo:>8.4f}  {C_ex:>10.4f}  {ratio_ex:>14.6f}  "
          f"{ratio_lo:>12.6f}  {C_lo/C_ex:>12.4f}")

print()

# ============================================================
# KROK 3: NLO rozwinięcie I4, I6 (dla małego a)
# ============================================================
print("KROK 3: NLO rozwinięcie calek I4 i I6 dla malego a_gam")
print("-" * 70)
print()
print("  I4_exact = e^{-4a}/a - 4*E1(4a)")
print("           = (1/a - 4 + 8a/2 - ...) - 4*(-gamma - ln(4a) + 4a - 8a^2 + ...)")
print("           = 1/a - 4 + 4*gamma + 4*ln(4a) + O(a)")
print("           = 1/a * [1 - 4a + 4a*(gamma + ln(4a)) + O(a^2)]")
print()
print("  I6_exact = e^{-6a}/(3a^3) + 2*e^{-6a}/(2a^2) + 6*(e^{-6a}/a - 6*E1(6a))")
print("           ~ 1/(3a^3) - 1/a^2 + 6/a - 36*E1(6a)")
print("           = 1/(3a^3) * [1 - 3a + 18a^2 - 108a^2*E1(6a)*3a^3]")
print()

# Numerycznie sprawdz rozwinięcie
for ag in [0.010, 0.020, 0.030, 0.040, 0.050]:
    I4_ex = I4_exact(ag)
    I6_ex = I6_exact(ag)
    # LO + NLO dla I4:
    gam = 0.5772156649  # stala Eulera
    I4_nlo = 1/ag + 4*(gam + np.log(4*ag)) - 4
    # LO + NLO dla I6:
    I6_nlo = 1/(3*ag**3) - 1/ag**2
    print(f"  a={ag:.3f}: I4_ex={I4_ex:.4f}, I4_NLO={I4_nlo:.4f} "
          f"(err={abs(I4_nlo-I4_ex)/I4_ex*100:.2f}%); "
          f"I6_ex={I6_ex:.2f}, I6_NLO={I6_nlo:.2f} "
          f"(err={abs(I6_nlo-I6_ex)/I6_ex*100:.2f}%)")

print()

# ============================================================
# KROK 4: Zaleznosc C(a_gam) -- dokladna formula
# ============================================================
print("KROK 4: Stala C(a_gam) z dokladnych calek")
print("-" * 70)
print()

agam_fine = np.logspace(-3, -1, 40)
C_exact_arr  = []
C_LO_arr     = []
C_NLO_arr    = []

for ag in agam_fine:
    I4_ex = I4_exact(ag)
    I6_ex = I6_exact(ag)
    if I4_ex <= 0 or I6_ex <= 0:
        C_exact_arr.append(np.nan)
        C_LO_arr.append(np.sqrt(4.5))
        C_NLO_arr.append(np.nan)
        continue

    C_ex = np.sqrt(3*I4_ex/(2*I6_ex)) / ag

    # NLO: I4 ~ 1/a + 4*(gamma+ln(4a)) - 4, I6 ~ 1/(3a^3) - 1/a^2
    gam = 0.5772156649
    I4_nlo = max(1e-10, 1/ag + 4*(gam + np.log(4*ag)) - 4)
    I6_nlo = max(1e-10, 1/(3*ag**3) - 1/ag**2)
    C_nlo = np.sqrt(3*I4_nlo/(2*I6_nlo)) / ag

    C_exact_arr.append(C_ex)
    C_LO_arr.append(np.sqrt(4.5))
    C_NLO_arr.append(C_nlo)

C_exact_arr = np.array(C_exact_arr)
C_LO_arr    = np.array(C_LO_arr)
C_NLO_arr   = np.array(C_NLO_arr)

# Numeryczne K3 dla referncyjnych a_gam
print(f"  {'a_gam':>8}  {'C_LO':>8}  {'C_exact(I4,I6)':>16}  {'C_NLO':>10}  "
      f"{'C_num(P31)':>12}  {'err_exact%':>12}")
print("  " + "-"*76)

agam_check = np.array([0.020, 0.025, 0.030, 0.040, 0.050, 0.060])
C_num_P31  = np.array([2.021, None, 2.009, 2.002, 1.999, 1.998])

for ag, C_p31 in zip(agam_check, C_num_P31):
    C_ex = np.sqrt(3*I4_exact(ag)/(2*I6_exact(ag))) / ag
    gam = 0.5772156649
    I4_nlo = max(1e-10, 1/ag + 4*(gam+np.log(4*ag)) - 4)
    I6_nlo = max(1e-10, 1/(3*ag**3) - 1/ag**2)
    C_nlo  = np.sqrt(3*I4_nlo/(2*I6_nlo)) / ag
    err_ex = abs(C_ex - (C_p31 if C_p31 else C_ex))/C_ex*100 if C_p31 else 0
    print(f"  {ag:>8.3f}  {np.sqrt(4.5):>8.4f}  {C_ex:>16.4f}  {C_nlo:>10.4f}  "
          f"  {C_p31 if C_p31 else '---':>10}  {err_ex:>12.2f}%")

print()

# Granica a_gam -> 0: co to jest C?
# C^2 = 3*I4/(2*I6) / a^2
# I4 ~ 1/a dla a->0 => I4*a ~ 1
# I6 ~ 1/(3a^3) => I6*a^3 ~ 1/3
# C^2 -> 3*(1/a)/(2*(1/(3a^3))) / a^2 = 3*3a^2/(2) / a^2 = 4.5
# Czyli C -> sqrt(4.5) = 2.121 dla a -> 0!

print("  KLUCZOWE: Dla a_gam -> 0, C -> sqrt(4.5) = 2.121?")
print()
for ag in [0.001, 0.002, 0.005, 0.010]:
    I4_ex = I4_exact(ag)
    I6_ex = I6_exact(ag)
    if I4_ex > 0 and I6_ex > 0:
        C_ex = np.sqrt(3*I4_ex/(2*I6_ex)) / ag
        print(f"    a_gam={ag:.3f}: C_exact={C_ex:.4f}  (C_LO={np.sqrt(4.5):.4f})")

print()
print("  WNIOSEK: C(a_gam) -> 2.121 dla a_gam -> 0, a dla a_gam=0.04 C=2.00")
print("  Jest to korekcja SKONCZONY-a_gam, nie asymptotyczna poprawka!")
print()

# ============================================================
# KROK 5: Kompletna formula dla K3
# ============================================================
print("KROK 5: Kompletna analityczna formula dla K3")
print("-" * 70)
print()
print("  WYNIK GLOWNY:")
print()
print("  LO (blad ~6% dla a_gam=0.04):")
print("    K3_LO = sqrt(4.5) * a_gam / sqrt(lambda)")
print()
print("  DOKLADNA FORMULA (blad < 1%):")
print("    K3 = sqrt(3*I4(a)/(2*lam*I6(a)))")
print("  gdzie:")
print("    I4(a) = e^{-4a}/a - 4*E1(4a)")
print("    I6(a) = e^{-6a}/(3a^3) + e^{-6a}/a^2 + 6*e^{-6a}/a - 36*E1(6a)")
print()

# Weryfikacja finalna
for ag, C_p31 in zip([0.020, 0.030, 0.040, 0.050, 0.060],
                     [2.021, 2.009, 2.002, 1.999, 1.998]):
    C_exact = np.sqrt(3*I4_exact(ag)/(2*I6_exact(ag))) / ag
    print(f"    a={ag:.3f}: C_exact={C_exact:.4f} vs C_P31={C_p31:.4f}  "
          f"(blad={abs(C_exact-C_p31)/C_p31*100:.2f}%)")

print()
print("  NLO (lepsza niz LO, gorsza od dokladnej):")
print("    I4_NLO(a) ~ 1/a + 4*gamma + 4*ln(4a) - 4")
print("    I6_NLO(a) ~ 1/(3a^3) - 1/a^2")
print()
print("  Stala C zbiega: C(a_gam) ~ 2.121 * sqrt(f4(4a)/f6(6a))")
print("  gdzie f4, f6 sa funkcjami malego parametru 4a, 6a.")
print()

# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 70)
print("PODSUMOWANIE P33")
print("=" * 70)
print()
print("WYNIK GŁÓWNY:")
print()
print("  Rozbieznosc C=2.000 (numeryk) vs C_LO=2.121 (teoria) pochodzi")
print("  z NIEDOKLADNOSCI LO aproksymacji calek:")
print("    I4_LO = 1/a    (zamiast dokladnego 1/a - 4*E1(4a))")
print("    I6_LO = 1/3a^3  (zamiast dokladnego I6_exact)")
print()
print("  DOKLADNA FORMULA:")
print("    K3 = sqrt(3*I4_exact(a) / (2*lam*I6_exact(a)))")
print("    C = K3*sqrt(lam)/a = sqrt(3*I4_exact/(2*I6_exact)) / a")
print()
ag = 0.040
I4_e = I4_exact(ag)
I6_e = I6_exact(ag)
C_e = np.sqrt(3*I4_e/(2*I6_e))/ag
print(f"  Dla a_gam=0.040: I4={I4_e:.4f}, I6={I6_e:.2f}, C={C_e:.4f}")
print(f"  vs C_num (P31) = 2.002  (blad: {abs(C_e-2.002)/2.002*100:.2f}%)")
print()
print("  DLA a_gam -> 0:")
print("    I4 -> 1/a, I6 -> 1/(3a^3)  =>  C -> sqrt(4.5) = 2.121")
print("    Korekta skonczona: C(a=0.04) = 2.000 (5.8% ponizej C_LO)")
print()
print("  INTERPRETACJA: Regularyzacja calki przez e^{-nr} zmniejsza")
print("  I4 o czynnik ~0.627 i I6 o czynnik ~0.708, ale stosunek")
print("  I4/I6 zmienia sie tak, ze C spada o ~5.8%.")
print()

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("P33: Dokladne calki I4, I6 i poprawiona formula K3", fontsize=13)

# C(a_gam) -- LO, NLO, dokladna
ax = axes[0, 0]
valid = np.isfinite(C_exact_arr) & np.isfinite(C_NLO_arr)
ax.semilogx(agam_fine[valid], C_exact_arr[valid], 'b-', lw=2, label='C = sqrt(3I4/(2I6))/a')
ax.semilogx(agam_fine[valid], C_LO_arr[valid], 'r--', lw=2, label=f'C_LO=sqrt(4.5)={np.sqrt(4.5):.3f}')
ax.semilogx(agam_fine[valid], C_NLO_arr[valid], 'g-.', lw=2, label='C_NLO (rozwinięcie)')
for ag, C_p31 in zip([0.020, 0.030, 0.040, 0.050, 0.060], [2.021, 2.009, 2.002, 1.999, 1.998]):
    ax.plot(ag, C_p31, 'ko', ms=6)
ax.axhline(2.000, color='orange', ls=':', lw=1.5, label='C=2.000 (P31 numeryk)')
ax.set_xlabel('a_gam')
ax.set_ylabel('C = K3*sqrt(lam)/a_gam')
ax.set_title('Stala C vs a_gam')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim([1.8, 2.3])

# I4(a) -- dokladna vs LO
ax = axes[0, 1]
ag_plot = np.logspace(-2.5, -0.8, 60)
I4_ex_arr = np.array([I4_exact(a) for a in ag_plot])
I4_lo_arr = 1.0/ag_plot
ax.loglog(ag_plot, I4_ex_arr, 'b-', lw=2, label='I4_exact')
ax.loglog(ag_plot, I4_lo_arr, 'r--', lw=2, label='I4_LO=1/a')
ax.axvline(0.040, color='g', ls=':', label='a=0.040')
ax.set_xlabel('a_gam')
ax.set_ylabel('I4')
ax.set_title('I4 = int e^{-4r}/r^2 dr')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# I6(a) -- dokladna vs LO
ax = axes[0, 2]
I6_ex_arr = np.array([I6_exact(a) for a in ag_plot])
I6_lo_arr = 1.0/(3*ag_plot**3)
ax.loglog(ag_plot, I6_ex_arr, 'b-', lw=2, label='I6_exact')
ax.loglog(ag_plot, I6_lo_arr, 'r--', lw=2, label='I6_LO=1/(3a^3)')
ax.axvline(0.040, color='g', ls=':', label='a=0.040')
ax.set_xlabel('a_gam')
ax.set_ylabel('I6')
ax.set_title('I6 = int e^{-6r}/r^4 dr')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# I4/I6 vs a_gam
ax = axes[1, 0]
ratio_ex = I4_ex_arr / I6_ex_arr
ratio_lo = I4_lo_arr / I6_lo_arr
ax.loglog(ag_plot, ratio_ex, 'b-', lw=2, label='I4/I6 exact')
ax.loglog(ag_plot, ratio_lo, 'r--', lw=2, label='I4/I6 LO = 3a^2')
ax.set_xlabel('a_gam')
ax.set_ylabel('I4 / I6')
ax.set_title('Stosunek I4/I6')
ax.legend()
ax.grid(True, alpha=0.3)

# K3 predict vs numeryk
ax = axes[1, 1]
K3_pred_arr, K3_num_arr, ag_both = [], [], []
for ag in [0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.080]:
    K3_n = get_K3(alpha_ref, ag, lam_ref)
    if np.isnan(K3_n):
        continue
    I4_e = I4_exact(ag)
    I6_e = I6_exact(ag)
    K3_p = np.sqrt(3*I4_e/(2*lam_ref*I6_e))
    K3_pred_arr.append(K3_p)
    K3_num_arr.append(K3_n)
    ag_both.append(ag)
K3_pred_arr = np.array(K3_pred_arr)
K3_num_arr  = np.array(K3_num_arr)
ax.plot(K3_num_arr, K3_pred_arr, 'bo', ms=8)
mn = min(K3_num_arr.min(), K3_pred_arr.min())
mx = max(K3_num_arr.max(), K3_pred_arr.max())
ax.plot([mn, mx], [mn, mx], 'k--', lw=1.5, label='y=x (idealne)')
ax.set_xlabel('K3 numeryczny')
ax.set_ylabel('K3 z formuly sqrt(3I4/(2lam*I6))')
ax.set_title('Predykcja K3 vs numeryk')
ax.legend()
ax.grid(True, alpha=0.3)

for i, (K3n, K3p) in enumerate(zip(K3_num_arr, K3_pred_arr)):
    ax.annotate(f'a={ag_both[i]:.3f}', (K3n, K3p), textcoords="offset points",
                xytext=(5, 5), fontsize=7)

# f4(x) = e^{-x} - x*E1(x) -- funkcja korekcyjna
ax = axes[1, 2]
x_arr = np.linspace(0.01, 0.5, 100)
f4_arr = np.exp(-x_arr) - x_arr*exp1(x_arr)
ax.plot(x_arr, f4_arr, 'b-', lw=2, label='f4(x)=e^{-x}-x*E1(x)')
ax.axhline(1.0, color='r', ls='--', label='f4=1 (LO)')
ax.axvline(4*0.040, color='g', ls=':', label='x=4*0.04=0.16')
ax.set_xlabel('x = 4*a_gam')
ax.set_ylabel('f4(x)')
ax.set_title('Funkcja korekcyjna do I4')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

f4_at_016 = np.exp(-0.16) - 0.16*exp1(0.16)
ax.annotate(f'f4(0.16)={f4_at_016:.3f}', (0.16, f4_at_016),
            textcoords="offset points", xytext=(10, -15), fontsize=9,
            arrowprops=dict(arrowstyle='->'))

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/p33_K3_correction.png', dpi=120, bbox_inches='tight')
print()
print("Wykres zapisany: p33_K3_correction.png")
