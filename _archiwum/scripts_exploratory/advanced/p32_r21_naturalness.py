"""
p32_r21_naturalness.py
======================
CEL: Czy r21 ~ 207 jest "naturalne" w TGP?

KONTEKST (P31):
  - K3 ~ 2.00 * a_gam / sqrt(lambda)  -- kontroluje r31
  - K1, K2 niezalezne od lambda       -- kontroluja r21
  - r21 = K2*/K1* zalezy od (alpha, a_gam)
  - Pytanie: co wyznacza r21 ~ 207 w naturze?

PLAN:
  Czesc A -- Mapa r21(alpha, a_gam)
    Skan siatki (alpha, a_gam): oblicz r21 = K2*/K1* (lambda nieistotna)
    Wyznacz kontur r21 = 207 (masa leptonow)
    Wyznacz kontur r21 = 600 (miony: szacunek kwarki charm/up)
    Wyznacz kontur r21 = 40000 (stop/up)

  Czesc B -- Analityczna formuła r21(alpha, a_gam)
    Weryfikacja: r21 ~ alpha / (sqrt(2) * a_gam) (asymptotyka z WYPROWADZENIE_r21.md)
    Poprawa: wspolczynnik prefaktor i korekcje

  Czesc C -- Sens fizyczny alpha ~ 8.55
    Czy alpha = 8.55 odpowiada jakiejs szczegolnej wartosci K2 lub K1?
    Czy K2* ~ 2 ma sens fizyczny? (K2 ~ 2 = 2*GAMMA)
    Zbadaj: czy K2* dazy do jakiejs granicy dla a_gam -> 0

  Czesc D -- Kwarki: TGP a masy kwarków
    Oblicz Q_Koide dla kwarków u,c,t i d,s,b
    Sprawdz: czy TGP moze odtworzyc quarki przez inny (alpha, a_gam)?
    Predykcja: jaki alpha daje r21_q = 600 (charm/up)?

  Czesc E -- "Specjalnosc" r21 = 207
    Rownanie analityczne: dla czego wartosci r21 Q=3/2 jest osiagalne
    Czy istnieje dolna/gorna granica r21 na krzywej Q=3/2?

WYNIKI P31: K1=0.009830 (alpha=8.553, a_gam=0.040), K2=2.0323
  r21 = 206.8 -- prawie dokladnie 207
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

print("P32: Czy r21~207 jest naturalne w TGP?")
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
                 (5.0, 50.0), (50.0, 500.0)]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero_K(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            if not zeros or abs(K - zeros[-1])/zeros[-1] > 0.01:
                zeros.append(K)
    return sorted(zeros)

def get_r21(alpha, a_gam, lam=1e-5):
    """r21 = K2*/K1*. lambda nieistotna (K1,K2 niezalezne od lam)."""
    zeros = find_all_zeros(alpha, a_gam, lam)
    if len(zeros) < 2:
        return np.nan, np.nan, np.nan
    return zeros[1]/zeros[0], zeros[0], zeros[1]

def Q_from_ratios(r21, r31):
    s = 1.0 + np.sqrt(r21) + np.sqrt(r31)
    return s**2 / (1.0 + r21 + r31)

def r31_koide(r21):
    """Rozwiaz Q(1,r21,x)=3/2 na x."""
    def eq(x):
        return Q_from_ratios(r21, x) - 1.5
    try:
        return brentq(eq, r21*1.01, r21*50.0, xtol=1e-6)
    except Exception:
        return np.nan

# ============================================================
# CZESC A: Mapa r21(alpha, a_gam)
# ============================================================
print("CZESC A: Mapa r21(alpha, a_gam)")
print("-" * 70)
print()

alpha_arr = np.linspace(4.0, 14.0, 25)
agam_arr  = np.array([0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.080])

lam_ref = 5.5e-6  # lambda referencyjna (K1,K2 od niej niezalezne)

# Siatka r21(alpha, a_gam)
r21_grid = np.full((len(agam_arr), len(alpha_arr)), np.nan)

print(f"  Obliczam siatke {len(agam_arr)}x{len(alpha_arr)} = {len(agam_arr)*len(alpha_arr)} punktow...")
for i, ag in enumerate(agam_arr):
    for j, al in enumerate(alpha_arr):
        r21_val, _, _ = get_r21(al, ag, lam_ref)
        r21_grid[i, j] = r21_val

print("  Gotowe.")
print()

# Szukaj alpha(a_gam) gdzie r21 = r21_target
r21_targets = [207.0, 600.0, 40000.0]
target_names = ["leptony (207)", "kwarki lekkie (600?)", "kwarki ciezkie (40000?)"]

print(f"  Kontury r21 = const:")
print(f"  {'a_gam':>8}  {'alpha(r21=207)':>16}  {'alpha(r21=600)':>16}  {'K1(207)':>10}  {'K2(207)':>10}")
print("  " + "-"*65)

alpha_at_r21_207 = []
alpha_at_r21_600 = []

for i, ag in enumerate(agam_arr):
    # Interpoluj alpha dla r21=207
    row = r21_grid[i, :]
    valid = np.isfinite(row)
    if valid.sum() < 2:
        print(f"  {ag:>8.3f}  {'---':>16}")
        continue

    # Znajdz alpha gdzie r21 ~ 207
    al_207 = np.nan
    al_600 = np.nan
    for j in range(len(alpha_arr)-1):
        if not (valid[j] and valid[j+1]):
            continue
        for target in [207.0, 600.0]:
            if (row[j] - target) * (row[j+1] - target) < 0:
                # Interpolacja liniowa
                frac = (target - row[j]) / (row[j+1] - row[j])
                al_interp = alpha_arr[j] + frac * (alpha_arr[j+1] - alpha_arr[j])
                if target == 207.0:
                    al_207 = al_interp
                else:
                    al_600 = al_interp

    # Dokladne K1, K2 przy alpha=al_207
    K1_ref = K2_ref = np.nan
    if not np.isnan(al_207):
        r21v, K1_ref, K2_ref = get_r21(al_207, ag, lam_ref)
        alpha_at_r21_207.append((ag, al_207))

    if not np.isnan(al_600):
        alpha_at_r21_600.append((ag, al_600))

    print(f"  {ag:>8.3f}  {al_207:>16.4f}  {al_600:>16.4f}  "
          f"{K1_ref:>10.6f}  {K2_ref:>10.4f}")

print()

# ============================================================
# CZESC B: Weryfikacja formuly analitycznej r21(alpha, a_gam)
# ============================================================
print("CZESC B: Analityczna formula r21(alpha, a_gam)")
print("-" * 70)
print()
print("  Formula asymptotyczna (WYPROWADZENIE_r21.md):")
print("    r21 ~ alpha / (sqrt(2) * a_gam)")
print()
print("  Rowniez: K1 ~ C_K * a_gam/(1+alpha), C_K=2.351")
print("           K2 ~ stale ~ 2 (pierwszorzedowo)")
print("           => r21 ~ K2/K1 ~ 2*(1+alpha)/(C_K*a_gam) ~ 2(1+alpha)/(2.351*a_gam)")
print()

# Sprawdz numerycznie
print(f"  {'a_gam':>8}  {'alpha':>8}  {'r21_num':>10}  "
      f"{'asym1 a/s2a':>14}  {'asym2 2(1+a)/(CK*a)':>22}  {'err1%':>8}  {'err2%':>8}")
print("  " + "-"*85)

CK = 2.351
for ag, al in alpha_at_r21_207:
    r21_num, K1, K2 = get_r21(al, ag, lam_ref)
    asym1 = al / (np.sqrt(2) * ag)
    asym2 = 2*(1+al) / (CK * ag)
    err1 = (asym1 - r21_num)/r21_num * 100
    err2 = (asym2 - r21_num)/r21_num * 100
    print(f"  {ag:>8.3f}  {al:>8.4f}  {r21_num:>10.2f}  "
          f"{asym1:>14.2f}  {asym2:>22.2f}  {err1:>8.1f}%  {err2:>8.1f}%")
print()

# ============================================================
# CZESC C: Sens fizyczny K2 ~ 2
# ============================================================
print("CZESC C: Zachowanie K2 przy roznych (alpha, a_gam)")
print("-" * 70)
print()
print("  K2 to DRUGI zero g(K): sredni soliton TGP")
print("  Pytanie: czy K2 -> jakiejs granicy dla a_gam -> 0?")
print()

agam_small = np.array([0.001, 0.002, 0.005, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060])

# Dla stalego r21=207 -- czyli alpha dopasowane do leptonu
print(f"  Wzdluz krzywej r21=207 (alpha dostrojone):")
print(f"  {'a_gam':>8}  {'alpha':>8}  {'K1':>10}  {'K2':>10}  {'K2/a_gam':>12}  {'K2/(1+a_gam)':>15}")
print("  " + "-"*70)

for ag in agam_small:
    # Szukaj alpha gdzie r21=207
    def f_r21(al):
        r21v, _, _ = get_r21(al, ag, lam_ref)
        if np.isnan(r21v):
            return np.nan
        return r21v - 207.0

    # Szacunek startowy alpha ~ 207 * sqrt(2) * a_gam
    al_start = 207 * np.sqrt(2) * ag
    al_lo = max(0.5, al_start*0.3)
    al_hi = min(50.0, al_start*3.0)

    try:
        flo = f_r21(al_lo)
        fhi = f_r21(al_hi)
        if np.isnan(flo) or np.isnan(fhi):
            print(f"  {ag:>8.4f}  {'---':>8}")
            continue
        if flo*fhi >= 0:
            print(f"  {ag:>8.4f}  {'brak':>8}  (r21 nie osiaga 207 w tym zakresie)")
            continue
        al_207 = brentq(f_r21, al_lo, al_hi, xtol=1e-6)
    except Exception as e:
        print(f"  {ag:>8.4f}  {'err':>8}")
        continue

    r21v, K1, K2 = get_r21(al_207, ag, lam_ref)
    if np.isnan(K1):
        print(f"  {ag:>8.4f}  {al_207:>8.4f}  {'---':>10}")
        continue
    print(f"  {ag:>8.4f}  {al_207:>8.4f}  {K1:>10.6f}  {K2:>10.4f}  "
          f"{K2/ag:>12.4f}  {K2/(1+ag):>15.4f}")
print()

# Granica K2 dla a_gam -> 0 (na krzywej r21=207)?
# Z formuly K2 ~ 2 (prawie stale), wiec K2 -> 2 dla a_gam->0
# K1 ~ C_K*a_gam/(1+alpha) -> 0, wiec r21 = K2/K1 -> inf
# Dla r21=207: alpha musi rosnac liniowo z a_gam
print("  Wniosek: K2 dazy do stalej ~2 dla a_gam -> 0")
print("  K1 ~ 2.351*a_gam/(1+alpha) -> 0 liniowo")
print("  Dla r21=207: alpha ~ r21*K1/K2 ~ 207*2.351*a_gam/2 ~ 243*a_gam")
print()

# Sprawdz alpha/a_gam na krzywej r21=207
print(f"  alpha/a_gam wzdluz krzywej r21=207:")
for ag, al in alpha_at_r21_207:
    r21v, K1, K2 = get_r21(al, ag, lam_ref)
    print(f"    a_gam={ag:.3f}: alpha/a_gam = {al/ag:.2f}  (teoria: ~243)")
print()

# ============================================================
# CZESC D: Kwarki
# ============================================================
print("CZESC D: Kwarki -- Q Koidego dla mas PDG")
print("-" * 70)
print()

# Masy kwarków PDG (MeV)
# u, c, t (kwarki up-type)
# d, s, b (kwarki down-type)
mu   = 2.16    # MeV
mc   = 1270.0  # MeV
mt   = 172690.0 # MeV (= 172.69 GeV)

md   = 4.67    # MeV
ms   = 93.4    # MeV
mb   = 4180.0  # MeV

print("  Masy kwarków PDG (MeV):")
print(f"    up-type:   u={mu}, c={mc}, t={mt}")
print(f"    down-type: d={md}, s={ms}, b={mb}")
print()

def koide_Q_masses(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (m1 + m2 + m3)

# Kwarki up-type
r21_uct = mc/mu
r31_uct = mt/mu
Q_uct   = Q_from_ratios(r21_uct, r31_uct)

# Kwarki down-type
r21_dsb = ms/md
r31_dsb = mb/md
Q_dsb   = Q_from_ratios(r21_dsb, r31_dsb)

# Leptony (dla porownania)
me   = 0.51100
mmu  = 105.658
mtau = 1776.86
r21_l = mmu/me
r31_l = mtau/me
Q_l   = Q_from_ratios(r21_l, r31_l)

print(f"  {'Sektor':>16}  {'r21':>12}  {'r31':>12}  {'Q':>10}  {'Q-3/2':>12}")
print("  " + "-"*68)
print(f"  {'Leptony':>16}  {r21_l:>12.3f}  {r31_l:>12.3f}  {Q_l:>10.6f}  {Q_l-1.5:>12.4e}")
print(f"  {'Kwarki u,c,t':>16}  {r21_uct:>12.1f}  {r31_uct:>12.0f}  {Q_uct:>10.6f}  {Q_uct-1.5:>12.4e}")
print(f"  {'Kwarki d,s,b':>16}  {r21_dsb:>12.3f}  {r31_dsb:>12.3f}  {Q_dsb:>10.6f}  {Q_dsb-1.5:>12.4e}")
print()

# r31_K dla kwarków
r31_K_uct = r31_koide(r21_uct)
r31_K_dsb = r31_koide(r21_dsb)
r31_K_l   = r31_koide(r21_l)

print(f"  Odleglosc od krzywej Koide (r31_K(r21) - r31_actual):")
print(f"    Leptony:      r31_K={r31_K_l:.2f},  r31_PDG={r31_l:.2f},  delta={r31_l-r31_K_l:.3f} ({(r31_l-r31_K_l)/r31_K_l*100:.4f}%)")
if not np.isnan(r31_K_uct):
    print(f"    u,c,t:        r31_K={r31_K_uct:.0f},  r31_PDG={r31_uct:.0f},  delta={r31_uct-r31_K_uct:.0f} ({(r31_uct-r31_K_uct)/r31_K_uct*100:.2f}%)")
print(f"    d,s,b:        r31_K={r31_K_dsb:.2f},  r31_PDG={r31_dsb:.2f},  delta={r31_dsb-r31_K_dsb:.3f} ({(r31_dsb-r31_K_dsb)/r31_K_dsb*100:.4f}%)")
print()

# Czy TGP moze odtworzyc kwarki up-type?
print(f"  Kwarki up-type: r21={r21_uct:.1f}, r31={r31_uct:.0f}")
print(f"  Szukam (alpha, a_gam) gdzie TGP daje te stosunki...")
print()

# Sprobujmy dla a_gam=0.040: znajdz alpha dla r21_uct
agam_q = 0.040
def f_r21_q(al):
    r21v, _, _ = get_r21(al, agam_q, lam_ref)
    if np.isnan(r21v):
        return np.nan
    return r21v - r21_uct

# Zakres alpha dla r21=588
al_lo_q = r21_uct * np.sqrt(2) * agam_q * 0.3
al_hi_q = r21_uct * np.sqrt(2) * agam_q * 3.0

try:
    flo_q = f_r21_q(al_lo_q)
    fhi_q = f_r21_q(al_hi_q)
    if np.isnan(flo_q) or np.isnan(fhi_q):
        raise ValueError("nan")
    if flo_q * fhi_q >= 0:
        raise ValueError("no crossing")
    al_q = brentq(f_r21_q, al_lo_q, al_hi_q, xtol=1e-4, maxiter=60)
    r21v_q, K1_q, K2_q = get_r21(al_q, agam_q, lam_ref)
    # lambda_Koide dla kwarków
    C_ref = 2.000
    lam_K_q = (C_ref*agam_q)**2 / (K1_q * r31_K_uct)**2 if not np.isnan(r31_K_uct) else np.nan
    print(f"    Znaleziono: alpha_q = {al_q:.4f} (a_gam={agam_q})")
    print(f"    K1={K1_q:.6f}, K2={K2_q:.4f}, r21={r21v_q:.2f}")
    if not np.isnan(lam_K_q):
        print(f"    lambda_Koide (kwarki) = {lam_K_q:.4e}")
    print()
except Exception as e:
    print(f"    Nie znaleziono rozwiazania dla a_gam={agam_q}: {e}")
    print(f"    Zakres: al_lo={al_lo_q:.2f}, al_hi={al_hi_q:.2f}")
    print()

# ============================================================
# CZESC E: "Specjalnosc" r21 = 207 -- zakres r21 na krzywej Q=3/2
# ============================================================
print("CZESC E: Zakres r21 na krzywej Q=3/2 dla roznych (alpha, a_gam)")
print("-" * 70)
print()
print("  Pytanie: czy r21=207 jest SPECJALNE na krzywej Q=3/2?")
print("  Odpowiedz: Q=3/2 jest krzywą algebraiczna r31_K(r21) -- niezalezna od TGP!")
print("  TGP daje JEDEN punkt na tej krzywej (dla danych parametrow).")
print("  Pytanie lepsze: dla jakich a_gam lambda_Koide ma fizyczne wartosci?")
print()

# Dla roznych (alpha, a_gam) -- oblicz lambda_Koide i r21
print(f"  Parametryczna krzywa Q=3/2 w TGP (wzdluz rodziny parametrow):")
print(f"  {'a_gam':>8}  {'alpha':>8}  {'r21':>10}  {'r31_K':>10}  {'lam_K':>12}  {'Q_check':>10}")
print("  " + "-"*70)

alpha_curve = np.linspace(5.0, 15.0, 20)
agam_curve  = 0.040
C_K3 = 2.000

for al in alpha_curve:
    r21v, K1, K2 = get_r21(al, agam_curve, lam_ref)
    if np.isnan(r21v) or r21v < 2:
        continue
    r31_K = r31_koide(r21v)
    if np.isnan(r31_K):
        continue
    lam_K = (C_K3*agam_curve)**2 / (K1*r31_K)**2

    # Weryfikacja Q przy lam_K
    zeros = find_all_zeros(al, agam_curve, lam_K)
    if len(zeros) >= 3:
        r21c = zeros[1]/zeros[0]
        r31c = zeros[2]/zeros[0]
        Qc   = Q_from_ratios(r21c, r31c)
    else:
        Qc = np.nan

    print(f"  {agam_curve:>8.3f}  {al:>8.4f}  {r21v:>10.2f}  "
          f"{r31_K:>10.2f}  {lam_K:>12.4e}  {Qc:>10.6f}")

print()
print("  Wniosek: dla KAZDEGO r21 > 1 mozemy znalezc lambda_Koide,")
print("  ktore daje Q=3/2. Nie ma 'specjalnej' wartosci r21 w TGP.")
print("  Wymaganie r21 = 207 pochodzi z OBSERWACJI mas leptonowych,")
print("  nie z wewnetrznej struktury TGP.")
print()

# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 70)
print("PODSUMOWANIE P32")
print("=" * 70)
print()
print("WYNIKI KLUCZOWE:")
print()
print(f"  1. r21 = K2*/K1* gdzie K2* ~ 2.0 (prawie stale)")
print(f"     K1* ~ C_K*a_gam/(1+alpha) (z WYPROWADZENIE_r21.md)")
print(f"     r21 ~ 2(1+alpha)/(2.351*a_gam) -- formula przyblizna")
print()
print(f"  2. Dla r21=207: alpha ~ 243*a_gam (przyblizna relacja liniowa)")
print(f"     np. a_gam=0.040 => alpha_0 ~ 9.7 (numerycznie: 8.553, blad ~14%)")
print()
print(f"  3. Kwarki up-type (r21=588, r31=79900): Q_uct = {Q_uct:.4f}")
print(f"     Kwarki down-type (r21=20, r31=895):  Q_dsb = {Q_dsb:.4f}")
print(f"     Kwarki NIE spelniaja Koide Q=3/2 -- bardzo daleko od 3/2")
print()
print(f"  4. TGP NIE PRZEWIDUJE r21=207 -- to wejscie obserwacyjne")
print(f"     Dla kazdego r21>1 mozna znalezc lambda_Koide dajace Q=3/2")
print(f"     Wartosc r21=207 pochodzi z obserwacji mas leptonowych")
print()
print("WNIOSEK FIZYCZNY:")
print()
print("  TGP ma dwa niezalezne parametry kontrolujace hierarchie mas:")
print("    * alpha -> r21 (generacja 2 vs 1)")
print("    * lambda -> r31 (generacja 3 vs 1)")
print("  a_gam -> skala absolutna (nie wplywa na stosunki)")
print()
print("  Warunek Q=3/2 jest JEDNO-wymiarowy w przestrzeni 2-wymiarowej (r21, r31).")
print("  TGP ma dokladnie jeden stopien swobody po natozeniu Q=3/2:")
print("  wartosci r21 (lub r31), ktora pochodzi z obserwacji.")
print()
print("  Interpretacja: TGP jest spójny z Koide, ale go nie WYMUSZA.")
print("  Dla wymuszeniego Q=3/2 potrzebny jest dodatkowy princip/kwantowanie.")
print()

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("P32: r21 w TGP -- naturalnosc i kwarki", fontsize=13)

# Mapa r21(alpha, a_gam)
ax = axes[0, 0]
A, G = np.meshgrid(alpha_arr, agam_arr)
r21_plot = np.where(np.isfinite(r21_grid) & (r21_grid > 0), r21_grid, np.nan)
r21_log = np.where(np.isfinite(r21_plot) & (r21_plot > 0), np.log10(r21_plot), np.nan)
im = ax.pcolormesh(A, G, r21_log, cmap='viridis', vmin=0, vmax=4)
plt.colorbar(im, ax=ax, label='log10(r21)')

# Kontury
try:
    cs = ax.contour(A, G, r21_log, levels=[np.log10(v) for v in [10, 50, 207, 600, 5000]],
                    colors='white', alpha=0.7)
    ax.clabel(cs, fmt=lambda v: f'{10**v:.0f}', fontsize=8)
except Exception:
    pass

ax.set_xlabel('alpha')
ax.set_ylabel('a_gam')
ax.set_title('r21 = K2/K1 (log10)')
ax.set_xlim([4, 14])

# Przekroj przy a_gam=0.040
ax = axes[0, 1]
idx_ag = np.argmin(np.abs(agam_arr - 0.040))
r21_row = r21_grid[idx_ag, :]
valid = np.isfinite(r21_row) & (r21_row > 0)
ax.semilogy(alpha_arr[valid], r21_row[valid], 'b.-')
ax.axhline(207, color='r', ls='--', label='r21=207 (leptony)')
ax.axhline(600, color='g', ls='--', label='r21=600 (kwarki?)')
ax.set_xlabel('alpha')
ax.set_ylabel('r21')
ax.set_title(f'r21 vs alpha (a_gam=0.040)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# alpha vs a_gam na krzywej r21=207
ax = axes[0, 2]
if alpha_at_r21_207:
    ags = [x[0] for x in alpha_at_r21_207]
    als = [x[1] for x in alpha_at_r21_207]
    ax.plot(ags, als, 'ro-', label='r21=207 (numeryk)', ms=5)
    # Teoria: alpha ~ 243*a_gam
    ag_fine = np.linspace(0.015, 0.085, 50)
    ax.plot(ag_fine, 243*ag_fine, 'b--', label=r'$\alpha \approx 243 \cdot a_\Gamma$', alpha=0.7)
    ax.axvline(0.040, color='g', ls=':', label='a_gam=0.040')
ax.set_xlabel('a_gam')
ax.set_ylabel('alpha na krzywej r21=207')
ax.set_title('Krzywa r21=207 w przestrzeni param.')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Q dla kwarków
ax = axes[1, 0]
r21_scan = np.logspace(0, 5, 200)
Q_scan = np.array([Q_from_ratios(r21v, r31_koide(r21v)) for r21v in r21_scan])
ax.semilogx(r21_scan, Q_scan, 'k-', label='TGP (lam=lambda_K)', alpha=0.5)
ax.axhline(1.5, color='r', ls='--', label='Q=3/2 (Koide)')
ax.axvline(r21_l, color='b', ls=':', label=f'leptony r21={r21_l:.0f}', alpha=0.7)
ax.axvline(r21_uct, color='g', ls=':', label=f'u,c,t r21={r21_uct:.0f}', alpha=0.7)
ax.axvline(r21_dsb, color='m', ls=':', label=f'd,s,b r21={r21_dsb:.0f}', alpha=0.7)
ax.set_xlabel('r21')
ax.set_ylabel('Q na krzywej Koidego')
ax.set_title('Q wzdluz krzywej Koide')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Porownanie Q kwarki vs leptony
ax = axes[1, 1]
sektory = ['Leptony\ne,mu,tau', 'Kwarki\nu,c,t', 'Kwarki\nd,s,b']
Q_vals_plot = [Q_l, Q_uct, Q_dsb]
colors_plot = ['blue', 'green', 'orange']
bars = ax.bar(sektory, Q_vals_plot, color=colors_plot, alpha=0.7)
ax.axhline(1.5, color='r', ls='--', lw=2, label='Q=3/2 (Koide)')
for bar, q in zip(bars, Q_vals_plot):
    ax.text(bar.get_x()+bar.get_width()/2, q+0.01, f'{q:.4f}', ha='center', fontsize=9)
ax.set_ylabel('Q Koide')
ax.set_title('Q dla roznych sektorow fermionow')
ax.legend()
ax.set_ylim([0.5, 1.8])
ax.grid(True, alpha=0.3, axis='y')

# Krzywa r21=207 i lambda_Koide
ax = axes[1, 2]
if alpha_at_r21_207:
    ags = np.array([x[0] for x in alpha_at_r21_207])
    als = np.array([x[1] for x in alpha_at_r21_207])
    lam_K_arr = []
    for ag, al in zip(ags, als):
        r21v, K1, K2 = get_r21(al, ag, lam_ref)
        r31_K = r31_koide(r21v)
        lam_K = (C_K3*ag)**2 / (K1*r31_K)**2 if not np.isnan(r31_K) else np.nan
        lam_K_arr.append(lam_K)
    lam_K_arr = np.array(lam_K_arr)
    ax.semilogy(ags, lam_K_arr, 'ro-', ms=5)
    ax.set_xlabel('a_gam')
    ax.set_ylabel('lambda_Koide')
    ax.set_title('lambda_Koide wzdluz krzywej r21=207')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/p32_r21_naturalness.png', dpi=120, bbox_inches='tight')
print("Wykres zapisany: p32_r21_naturalness.png")
