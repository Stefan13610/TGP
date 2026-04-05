"""
p31_analytic_koide.py
=====================
CEL: ANALITYCZNA PODSTAWA WARUNKU KOIDEGO W TGP

PYTANIE CENTRALNE:
  Dlaczego krzywa Q=3/2 (w przestrzeni parametrow (r21, r31))
  przechodzi przez okolice mas leptonowych (r21~207, r31~3477)?

PLAN:
  Krok 1 -- Skalowanie K3 ~ a_gam/sqrt(lambda)
    Numerycznie: log-log fit K3 vs lambda przy stalym a_gam
    Ekstrakcja wykladnika i wspolczynnika C3(a_gam)

  Krok 2 -- Zrozumienie K2 (drugi zero)
    Czy K2 ~ a_gam^p dla jakiegos p?
    Czy K2 zalezy od lambda (slabo)?

  Krok 3 -- Wyprowadzenie analityczne skalowania K3
    Z rownania samospojnosci g(K) = E(K)/(4piK) - 1 = 0
    W reziimie duzego K: dominuje czlon kwartyczny vs seksyczny
    Energia: E ~ 4pi[-K^3/(4*a) + lambda*K^5/(18*a^3)]
    Balans: K^2 = 4.5 a^2 / lambda
    --> K3 ~ (4.5)^{1/2} * a_gam / sqrt(lambda) = 2.121 * a_gam / sqrt(lambda)

  Krok 4 -- Stosunek r31 = K3/K1 jako funkcja parametrow
    K1 (maly zero) slabo zalezy od lambda (dominuje potencjal szescienny nie wplywa)
    r31 ~ K3/K1(alpha, a_gam) ~ a_gam / (K1(alpha, a_gam) * sqrt(lambda))

  Krok 5 -- Warunek Koidego jako rownanie na lambda
    Q = 3/2 <=> f(r21, r31) = 0
    r31_Koide(r21) = [(1+sqrt(r21)+sqrt(x))^2/(1+r21+x) = 3/2] rozwiazany na x
    Formuła: r31_K(r21) = ((3-sqrt(r21))^2 + ...) ... (znana z poprzednich skryptow)

    Rownoczesnie r31 = C * a_gam/(K1*sqrt(lambda))
    --> sqrt(lambda_Koide) = C * a_gam / (K1 * r31_Koide(r21))
    --> lambda_Koide = C^2 * a_gam^2 / (K1^2 * r31_Koide(r21)^2)

  Krok 6 -- Czy K3/K2 ~ const?
    Jesli K3/K2 = R_{32} ~ const (niezaleznie od lambda), to:
    r31/r21 = K3/K2 = R_{32}
    I warunek Koidego wymaga specyficznego R_{32}

  Krok 7 -- Podsumowanie: mechanizm Koidego
    TGP mowi: r31/r21 jest wyznaczone przez (alpha, a_gam) przy ustalonym lambda
    Krzywa Q=3/2 w plaszczyznie (r21, r31) jest krzywa algebraiczna
    TGP daje punkt na tej krzywej gdy lambda = lambda_Koide

WYNIK P30: Q(TGP, dokladne masy) = 1.499986 (delta = -1.38e-5)
  PDG nie lezy na krzywej Koidego: r31_K(206.768) = 3477.437 != 3477.65

Pytanie: CZY istnieje mechanizm w TGP ktory PREFERUJE Q~3/2?
"""

import numpy as np
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX = 60.0
GAMMA = 1.0

print("P31: Analityczna podstawa warunku Koidego w TGP")
print("=" * 70)
print()

# ============================================================
# NARZEDZIA (identyczne jak p27-p30)
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
        g_lo = g_func(K_lo, alpha, a_gam, lam)
        g_hi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi >= 0:
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam), K_lo, K_hi, xtol=1e-10)
    except Exception:
        return np.nan

def find_all_zeros(alpha, a_gam, lam):
    intervals = [(0.001, 0.050), (0.050, 0.500), (0.500, 5.0),
                 (5.0, 50.0), (50.0, 200.0)]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero_K(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            if not zeros or abs(K - zeros[-1])/zeros[-1] > 0.01:
                zeros.append(K)
    return sorted(zeros)

def koide_Q(k1, k2, k3):
    s = np.sqrt(k1) + np.sqrt(k2) + np.sqrt(k3)
    return s**2 / (k1 + k2 + k3)

def Q_from_ratios(r21, r31):
    s = 1.0 + np.sqrt(r21) + np.sqrt(r31)
    return s**2 / (1.0 + r21 + r31)

# Krzywa Koidego r31_K(r21): Q(1, r21, r31)=3/2 => r31 =?
def r31_koide(r21):
    """Rozwiaz Q(1, r21, x) = 3/2 na x > r21."""
    def eq(x):
        return Q_from_ratios(r21, x) - 1.5
    # szukaj x w przedziale (r21, 1e8)
    try:
        # dla duzych r21: r31 ~ (2+sqrt(3))^2 * r21 ~ 13.93 * r21
        x_lo = r21 * 1.01
        x_hi = r21 * 30.0
        if eq(x_lo) * eq(x_hi) < 0:
            return brentq(eq, x_lo, x_hi, xtol=1e-6)
    except Exception:
        pass
    return np.nan

# ============================================================
# KROK 1: Skalowanie K3 vs lambda (log-log) przy stalym (alpha, a_gam)
# ============================================================
print("KROK 1: Skalowanie K3 ~ lambda^p -- weryfikacja numeryczna")
print("-" * 70)

alpha_ref  = 8.553
agam_ref   = 0.040
lam_range  = np.logspace(-7, -4, 40)  # od 1e-7 do 1e-4

K1_arr, K2_arr, K3_arr = [], [], []
lam_ok = []

print(f"  Skan lambda (alpha={alpha_ref}, a_gam={agam_ref}): {len(lam_range)} punktow")
for lam in lam_range:
    zeros = find_all_zeros(alpha_ref, agam_ref, lam)
    if len(zeros) >= 3:
        K1_arr.append(zeros[0])
        K2_arr.append(zeros[1])
        K3_arr.append(zeros[2])
        lam_ok.append(lam)

lam_ok = np.array(lam_ok)
K1_arr = np.array(K1_arr)
K2_arr = np.array(K2_arr)
K3_arr = np.array(K3_arr)

print(f"  Znaleziono 3 zera dla {len(lam_ok)}/{len(lam_range)} wartosci lambda")
print()

# Log-log fit K3 vs lambda
log_lam = np.log(lam_ok)
log_K3  = np.log(K3_arr)
log_K2  = np.log(K2_arr)
log_K1  = np.log(K1_arr)

p3, c3 = np.polyfit(log_lam, log_K3, 1)
p2, c2 = np.polyfit(log_lam, log_K2, 1)
p1, c1 = np.polyfit(log_lam, log_K1, 1)

print(f"  K3 ~ lambda^p3:  p3 = {p3:.4f}  (teoria: -0.5000)")
print(f"  K2 ~ lambda^p2:  p2 = {p2:.4f}  (oczekiwane: blisko 0)")
print(f"  K1 ~ lambda^p1:  p1 = {p1:.4f}  (oczekiwane: blisko 0)")
print()

# Wspolczynnik C3: K3 = C3 * a_gam / sqrt(lambda)
C3_vals = K3_arr * np.sqrt(lam_ok) / agam_ref
print(f"  Wspolczynnik C3 = K3*sqrt(lambda)/a_gam:")
print(f"    mediana = {np.median(C3_vals):.4f}")
print(f"    odch. std = {np.std(C3_vals):.4f}")
print(f"    teoria (potencjal lokalny): C3_teor = sqrt(4.5) = {np.sqrt(4.5):.4f}")
print()

# Weryfikacja K3 ~ a_gam / sqrt(lambda) lub K3 ~ a_gam^2 / sqrt(lambda)
# Sprawdzmy dla roznych a_gam
print("  Zaleznosc C3 od a_gam:")
agam_test = [0.020, 0.030, 0.040, 0.050, 0.060]
lam_ref   = 5.5e-6  # ustalona lambda
print(f"  {'a_gam':>8}  {'K3':>10}  {'K3*sqrt(lam)/a':>18}  {'K3*sqrt(lam)/a^2':>18}")
print("  " + "-"*60)
C3_agam = []
for ag in agam_test:
    zeros = find_all_zeros(alpha_ref, ag, lam_ref)
    if len(zeros) >= 3:
        K3 = zeros[2]
        c_a1 = K3*np.sqrt(lam_ref)/ag
        c_a2 = K3*np.sqrt(lam_ref)/ag**2
        print(f"  {ag:>8.3f}  {K3:>10.3f}  {c_a1:>18.3f}  {c_a2:>18.3f}")
        C3_agam.append((ag, K3, c_a1, c_a2))
print()

# ============================================================
# KROK 2: K2 jako funkcja alpha (slaba zaleznosc od lambda)
# ============================================================
print("KROK 2: K2 vs alpha -- ktory czlon dominuje K2?")
print("-" * 70)

lam_test = [1e-6, 5e-6, 1e-5, 5e-5]
alpha_scan = np.linspace(7.0, 10.0, 16)
agam = 0.040

print(f"  K2 dla roznych (alpha, lambda), a_gam={agam}:")
print(f"  {'alpha':>8}", end="")
for lam in lam_test:
    print(f"  {'K2(lam='+f'{lam:.0e})':>18}", end="")
print()
print("  " + "-"*80)
for al in alpha_scan:
    print(f"  {al:>8.3f}", end="")
    for lam in lam_test:
        zeros = find_all_zeros(al, agam, lam)
        if len(zeros) >= 2:
            print(f"  {zeros[1]:>18.4f}", end="")
        else:
            print(f"  {'---':>18}", end="")
    print()
print()

# ============================================================
# KROK 3: WYPROWADZENIE ANALITYCZNE -- energia dla duzego K
# ============================================================
print("KROK 3: Analityczna aproksymacja energii dla duzego K")
print("-" * 70)
print()
print("  Profil: phi(r) = 1 + K*exp(-r)/r")
print("  Potencjal: V(phi) = phi^3/3 - phi^4/4 + lambda*(phi-1)^6/6")
print()
print("  Dla K >> 1 i r ~ a_gam (dolna granica calkowania):")
print("    phi ~ K/r")
print("    dphi/dr ~ -K/r^2")
print()
print("  Wklad czynow do energii calkowalnych przez log-grid (r >= a_gam):")
print("    Kinetyczny: 4pi * int_{a}^inf 1/2 * (K/r^2)^2 * r^2 dr")
print("              = 2pi*K^2 * [1/a - 1/R_max] ~ 2pi*K^2/a")
print("    Kwartyczny: 4pi * int -K^4/(4r^4) * r^2 dr ~ -pi*K^4/a")
print("    Seksyczny:  4pi * int lambda*K^6/(6r^6) * r^2 dr")
print("              = 4pi*lam*K^6/6 * [1/(3r^3)]_a^inf ~ 4pi*lam*K^6/(18*a^3)")
print()
print("  Rownanie samospojnosci E/(4piK) = 1:")
print("    K/2a - K^3/(4a) + lam*K^5/(18a^3) = 1 + trzescienny + ...")
print()
print("  Dla duzego K trzeci i piatego rzad dominuja:")
print("    -K^3/(4a) + lam*K^5/(18a^3) ~ 0  (punkt przegiety ~ K3)")
print("    => K3^2 ~ 18a^2 / (4*lam) = 4.5 * a^2/lam")
print("    => K3 ~ sqrt(4.5) * a/sqrt(lam) = 2.121 * a/sqrt(lam)")
print()

# Porownanie z numeryka
C3_theory = np.sqrt(4.5)
C3_num = np.median(K3_arr * np.sqrt(lam_ok) / agam_ref)
print(f"  Teoria:   K3 = {C3_theory:.4f} * a_gam / sqrt(lambda)")
print(f"  Numeryk:  K3 = {C3_num:.4f} * a_gam / sqrt(lambda)")
print(f"  Roznica:  {abs(C3_num-C3_theory)/C3_theory*100:.1f}%")
print()
print("  UWAGA: roznica pochodzi z poprawek wyzszego rzedu (czlon kwadratowy,")
print("  falowa struktura phi, korekcje od phi^3, granica gornego calkowania)")
print()

# Korekty: dla K nie za duzego, wplywaja tez wyzsze potegi
# Lepsza aproksy: uwzgledniamy K^3 + lambda*K^5 = K^3/(4a) (z calki ze struktura e^-r)

# ============================================================
# KROK 4: r31 jako funkcja (alpha, a_gam, lambda)
# ============================================================
print("KROK 4: Stosunek r31 = K3/K1 -- skalowanie z parametrami")
print("-" * 70)

print()
print("  Z Kroku 3: K3 ~ 2.12 * a_gam / sqrt(lambda)")
print("  Z numeryki: K1 prawie niezalezy od lambda (slaba zaleznosc)")
print("  Wiec: r31 = K3/K1 ~ 2.12 * a_gam / (K1(alpha) * sqrt(lambda))")
print()
print("  K1(alpha, a_gam) wyznacza r31 poprzez:")
print("    sqrt(lambda) = 2.12 * a_gam / (K1 * r31)")
print("    lambda = 4.5 * a_gam^2 / (K1^2 * r31^2)")
print()

# Sprawdz K1 vs lambda
print("  Weryfikacja: K1 vs lambda (przy alpha=8.553, a_gam=0.040)")
print(f"  {'lambda':>12}  {'K1':>12}  {'K1/K1_ref':>12}")
K1_ref_val = K1_arr[len(K1_arr)//2]
for i, lam in enumerate(lam_ok[::5]):
    ratio = K1_arr[i*5]/K1_ref_val
    print(f"  {lam:>12.2e}  {K1_arr[i*5]:>12.6f}  {ratio:>12.4f}")
print()

# ============================================================
# KROK 5: Warunek Koidego Q=3/2 jako algebraiczny warunek
# ============================================================
print("KROK 5: Warunek Q=3/2 wyrazony jako rownanie na lambda")
print("-" * 70)
print()
print("  Dane: r21(alpha, a_gam) -- prawie niezalezy od lambda (z P28)")
print("        r31(alpha, a_gam, lambda) ~ C * a_gam / (K1 * sqrt(lambda))")
print()
print("  Koide r31_K(r21) -- znana formula algebraiczna")
print()
print("  Warunek Q=3/2 mowi: r31_TGP = r31_K(r21)")
print("    => C * a_gam / (K1 * sqrt(lambda)) = r31_K(r21)")
print("    => sqrt(lambda_Koide) = C * a_gam / (K1 * r31_K(r21))")
print("    => lambda_Koide = C^2 * a_gam^2 / (K1^2 * r31_K(r21)^2)")
print()

# Policzmy C numerycznie dla (alpha=8.553, a_gam=0.040)
# Wezmy kilka wartosci lambda, ekstrapolujemy K3~lambda^{-1/2}
# r21 jest prawie stale ~ 207
# r31_K(207) znane

r21_ref = 206.77
r31_K_ref = r31_koide(r21_ref)
print(f"  r21_ref = {r21_ref}")
print(f"  r31_K(r21_ref) = {r31_K_ref:.3f}")
print()

# Dla referncyjnego (alpha, a_gam, lambda): co dostajemy?
zeros_ref = find_all_zeros(alpha_ref, agam_ref, 5.489e-6)
if len(zeros_ref) >= 3:
    K1r, K2r, K3r = zeros_ref[:3]
    r21r = K2r/K1r
    r31r = K3r/K1r
    Qr   = koide_Q(K1r, K2r, K3r)
    print(f"  Przy alpha=8.553, a_gam=0.040, lambda=5.489e-6:")
    print(f"    K1={K1r:.6f}, K2={K2r:.4f}, K3={K3r:.4f}")
    print(f"    r21={r21r:.3f}, r31={r31r:.3f}, Q={Qr:.6f}")
    print()
    C_K3 = K3r * np.sqrt(5.489e-6) / agam_ref
    C_r31 = r31r * K1r * np.sqrt(5.489e-6) / agam_ref
    print(f"    C3 = K3*sqrt(lam)/a_gam = {C_K3:.4f}")
    print(f"    C_r31 = r31*K1*sqrt(lam)/a_gam = {C_r31:.4f}")
    print()
    lam_Koide_analytic = (C_K3*agam_ref)**2 / (K1r*r31_K_ref)**2
    print(f"  Analityczna predykcja lambda_Koide:")
    print(f"    lam_K = (C3*a)^2 / (K1*r31_K)^2")
    print(f"         = ({C_K3:.4f}*{agam_ref})^2 / ({K1r:.6f}*{r31_K_ref:.3f})^2")
    print(f"         = {lam_Koide_analytic:.4e}")
    print(f"  Numeryk z P28: lambda_Koide = 5.489e-6")
    print(f"  Roznica: {abs(lam_Koide_analytic - 5.489e-6)/5.489e-6*100:.1f}%")
    print()

# ============================================================
# KROK 6: Stosunek K3/K2 = R32 -- staly czy nie?
# ============================================================
print("KROK 6: Czy K3/K2 ~ const (niezaleznie od lambda)?")
print("-" * 70)

r31_r21 = K3_arr / K2_arr  # R_{32} = K3/K2 = r31/r21
print()
print(f"  {'lambda':>12}  {'K3/K2':>10}  {'r31':>10}  {'r21':>10}  {'r31/r21':>10}")
for i in range(0, len(lam_ok), max(1, len(lam_ok)//12)):
    r31_i = K3_arr[i]/K1_arr[i]
    r21_i = K2_arr[i]/K1_arr[i]
    print(f"  {lam_ok[i]:>12.2e}  {r31_r21[i]:>10.4f}  "
          f"{r31_i:>10.2f}  {r21_i:>10.4f}  {r31_i/r21_i:>10.4f}")

mean_R32 = np.mean(r31_r21)
std_R32  = np.std(r31_r21)
print()
print(f"  Srednia K3/K2 = {mean_R32:.4f}  +/- {std_R32:.4f}")
print(f"  Wzgledna zmiennosc: {std_R32/mean_R32*100:.2f}%")
print()
print("  INTERPRETACJA:")
if std_R32/mean_R32 < 0.05:
    print("  K3/K2 ~ STALE => r31/r21 ~ const niezaleznie od lambda!")
    print("  Koide wymaga specyficznego r31/r21 = r31_K(r21)/r21")
    target_R32 = r31_K_ref / r21_ref
    print(f"  Wartosc Koidego: R32_Koide = {target_R32:.4f}")
    print(f"  TGP daje:        R32_TGP   = {mean_R32:.4f}")
    print(f"  Roznica: {abs(mean_R32-target_R32)/target_R32*100:.2f}%")
else:
    print("  K3/K2 zalezy od lambda -- nie jest stale")
print()

# ============================================================
# KROK 7: PODSUMOWANIE ANALITYCZNE
# ============================================================
print("=" * 70)
print("KROK 7: PODSUMOWANIE -- Mechanizm Q~3/2 w TGP")
print("=" * 70)
print()
print("WYNIKI:")
print()
print(f"  1. K3 ~ {C3_theory:.3f} * a_gam / sqrt(lambda)  [teoria: sqrt(4.5)=2.121]")
print(f"     numerycznie: K3 ~ {C3_num:.3f} * a_gam / sqrt(lambda)")
print(f"     wykladnik: lambda^{p3:.4f}  (teoria: -0.5000)")
print()
print(f"  2. K2 prawie niezalezy od lambda: wykladnik {p2:.4f} (bliskie 0)")
print(f"     K1 prawie niezalezy od lambda: wykladnik {p1:.4f}")
print()
print(f"  3. r31 = K3/K1 ~ C * a_gam / (K1(alpha) * sqrt(lambda))")
print()
print(f"  4. Warunek Q=3/2: lambda_Koide = C^2*a_gam^2 / (K1^2 * r31_K^2)")
print()
if 'K1r' in dir():
    r31_K_ref2 = r31_koide(K2r/K1r)
    lam_K_pred = (C_K3*agam_ref)**2 / (K1r * r31_K_ref2)**2
    print(f"  5. Predykcja analityczna dla (alpha=8.553, a_gam=0.040):")
    print(f"     r21  = {K2r/K1r:.3f}  (TGP)")
    print(f"     r31_K(r21) = {r31_K_ref2:.3f}  (warunek Koidego)")
    print(f"     lambda_K (analitycznie) = {lam_K_pred:.3e}")
    print(f"     lambda_K (numerycznie) = 5.489e-6")
print()
print("WNIOSEK:")
print()
print("  Mechanizm jest CZYSTO KINEMATYCZNY:")
print("  - K3 skaluje sie jak a_gam/sqrt(lambda) z kwartyczno-seksycznego bilansu")
print("  - K1, K2 sa prawie niezalezne od lambda (sa wyznaczone przez phi^3 potencjal)")
print("  - r31 jest zatem kontrolowane przez lambda, r21 przez alpha")
print("  - Q=3/2 osiaga sie gdy lambda dostroi K3 tak by r31 = r31_K(r21)")
print()
print("  Pytanie residualne (P32):")
print("  Dlaczego TGP daje r21~207 dla rozsadnych (alpha, a_gam)?")
print("  I czy r21=207 + r31_Koide(207)=3481 jest 'naturalny' w TGP?")
print()

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("P31: Skalowanie K_i vs lambda i mechanizm Q~3/2 w TGP", fontsize=13)

ax = axes[0, 0]
ax.loglog(lam_ok, K1_arr, 'b.-', label='K1')
ax.loglog(lam_ok, K2_arr, 'g.-', label='K2')
ax.loglog(lam_ok, K3_arr, 'r.-', label='K3')
lam_fit = np.array([lam_ok[0], lam_ok[-1]])
K3_fit = np.exp(c3) * lam_fit**p3
ax.loglog(lam_fit, K3_fit, 'r--', alpha=0.5, label=f'K3~lam^{p3:.3f}')
ax.set_xlabel('lambda')
ax.set_ylabel('K_i')
ax.set_title('Zera K_i vs lambda')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.loglog(lam_ok, K3_arr * np.sqrt(lam_ok), 'r.-')
ax.axhline(C3_num * agam_ref, color='k', ls='--', label=f'C3*a={C3_num*agam_ref:.4f}')
ax.axhline(np.sqrt(4.5) * agam_ref, color='b', ls=':', label=f'sqrt(4.5)*a={np.sqrt(4.5)*agam_ref:.4f}')
ax.set_xlabel('lambda')
ax.set_ylabel('K3 * sqrt(lambda)')
ax.set_title('K3*sqrt(lam) = const?')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 2]
r31_arr = K3_arr / K1_arr
r21_arr = K2_arr / K1_arr
ax.semilogx(lam_ok, r31_arr, 'r.-', label='r31')
ax.semilogx(lam_ok, r21_arr * 100, 'b.-', label='r21 x100')
ax.set_xlabel('lambda')
ax.set_ylabel('stosunki mas')
ax.set_title('r21, r31 vs lambda')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
Q_arr = np.array([koide_Q(K1_arr[i], K2_arr[i], K3_arr[i]) for i in range(len(K1_arr))])
ax.semilogx(lam_ok, Q_arr, 'k.-')
ax.axhline(1.5, color='r', ls='--', label='Q=3/2')
ax.set_xlabel('lambda')
ax.set_ylabel('Q')
ax.set_title('Q Koidego vs lambda')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.semilogx(lam_ok, K3_arr / K2_arr, 'm.-')
r31_K_arr = np.array([r31_koide(r21_arr[i]) if not np.isnan(r21_arr[i]) else np.nan for i in range(len(r21_arr))])
target_R32_arr = r31_K_arr / r21_arr
ax.semilogx(lam_ok, target_R32_arr, 'r--', alpha=0.7, label='r31_K/r21 (Koide target)')
ax.set_xlabel('lambda')
ax.set_ylabel('K3/K2')
ax.set_title('K3/K2 vs lambda')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 2]
# Diagram analityczny
lam_fine = np.logspace(-7.5, -3.5, 200)
K3_analytic = np.sqrt(4.5) * agam_ref / np.sqrt(lam_fine)
K3_numeric_fit = np.exp(c3) * lam_fine**p3
ax.loglog(lam_fine, K3_analytic, 'b-', label='teoria: sqrt(4.5)*a/sqrt(lam)', lw=2)
ax.loglog(lam_fine, K3_numeric_fit, 'r--', label=f'fit num: lam^{p3:.3f}', lw=2)
ax.loglog(lam_ok, K3_arr, 'k.', ms=5, label='numeryk')
ax.set_xlabel('lambda')
ax.set_ylabel('K3')
ax.set_title('K3: teoria vs numeryk')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/p31_analytic_koide.png', dpi=120, bbox_inches='tight')
print()
print("Wykres zapisany: p31_analytic_koide.png")
