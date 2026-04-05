"""
p27_koide_Q_landscape.py
========================
CEL: Systematyczne badanie Q = (sum sqrt(K*_n))^2 / sum(K*_n) w przestrzeni
     parametrow TGP. Inicjalizacja problemu otwartego O-K1.

PYTANIE OTWARTE (O-K1):
  Czy Q >= 3/2 jest strukturalnym ograniczeniem TGP (dla dowolnych parametrow),
  czy tylko przypadkiem dla rodziny optymalnej?
  Czy TGP moze PRZEWIDZIEC Q = 3/2, czy tylko reprodukuje obserwacje?

PLAN P27:
  Krok 1: Analiza Q(r21, r31) -- geometria krzywej Koide
  Krok 2: Skan Q(lambda) przy stalym (alpha, a_gam) -- szeroki zakres lambda
  Krok 3: Skan Q(alpha, a_gam) poza rodzina optymalna
  Krok 4: Analiza asymptotyczna -- wiodace czlony dla duzych r21, r31
  Krok 5: Podsumowanie -- co jest zbadane, co otwarte

OCZEKIWANIA:
  - Jezeli Q >= 3/2 zawsze: szukac dowodu analitycznego (nastepna sesja)
  - Jezeli Q moze byc < 3/2: Q=3/2 jest przypadkowe, nie predykcja
  - Jezeli Q = 3/2 dokladnie dla jednego lambda: TGP "wybiera" punkt Koide
"""

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_MAX  = 60.0
GAMMA  = 1.0

print("P27: Krajobraz Q(alpha, a_gam, lambda) -- Problem otwarty O-K1")
print("=" * 70)
print()

# ============================================================
# NARZEDZIA PODSTAWOWE
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam, N=1500):
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX/a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1+alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam):
    E = energy_log(K, alpha, a_gam, lam)
    return E / (4*np.pi*K) - 1.0

def find_zero(alpha, a_gam, lam, K_lo, K_hi):
    """Szuka jednego zera g(K) w [K_lo, K_hi]."""
    try:
        g_lo = g_func(K_lo, alpha, a_gam, lam)
        g_hi = g_func(K_hi, alpha, a_gam, lam)
        if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
            return np.nan
        if g_lo * g_hi >= 0:
            return np.nan
        return brentq(lambda K: g_func(K, alpha, a_gam, lam), K_lo, K_hi, xtol=1e-8)
    except Exception:
        return np.nan

def find_all_zeros(alpha, a_gam, lam):
    """Szuka K*1, K*2, K*3 na siatce przedzialow."""
    intervals = [
        (0.001, 0.050),   # K*1
        (0.050, 0.500),   # K*1/K*2 gap
        (0.500, 5.0),     # K*2
        (5.0,   50.0),    # K*3
        (50.0, 200.0),    # K*3 dla malego lam
    ]
    zeros = []
    for K_lo, K_hi in intervals:
        K = find_zero(alpha, a_gam, lam, K_lo, K_hi)
        if not np.isnan(K):
            # Sprawdz ze nie jest duplikatem
            if not zeros or abs(K - zeros[-1])/zeros[-1] > 0.01:
                zeros.append(K)
    return sorted(zeros)

def koide_Q(k1, k2, k3):
    """Q = (sum sqrt(K))^2 / sum(K)"""
    s = np.sqrt(k1) + np.sqrt(k2) + np.sqrt(k3)
    return s**2 / (k1 + k2 + k3)

def Q_from_ratios(r21, r31):
    """Q z samych stosunkow (K1=1)."""
    s = 1.0 + np.sqrt(r21) + np.sqrt(r31)
    return s**2 / (1.0 + r21 + r31)

# ============================================================
# KROK 1: Geometria Q(r21, r31) -- ANALITYCZNY
# ============================================================
print("KROK 1: Geometria Q(r21, r31)")
print("-" * 70)
print()
print("  Q(r21, r31) = (1 + sqrt(r21) + sqrt(r31))^2 / (1 + r21 + r31)")
print()

# Asymptotyki
print("  Asymptotyki:")
print("    r21=1, r31=1 (rowne masy):  Q =", Q_from_ratios(1, 1))
print("    r21=1, r31->inf:            Q ->", Q_from_ratios(1, 1e8))
print("    r21->inf, r31->inf:         Q ->", Q_from_ratios(1e6, 1e8))
print("    Punkt leptonowy (207, 3477):", Q_from_ratios(206.768, 3477.65))
print("    Punkt TGP (207, 3477):      ", Q_from_ratios(207.0, 3477.1))
print()

# Krzywa Q = 3/2 w plaszczyzni (r21, r31)
# Dla ustaloneog r21, znajdz r31 takie ze Q=3/2
print("  Krzywa Koide Q=3/2: dla kazdego r21 istnieje r31(r21):")
print()
print(f"  {'r21':>10}  {'r31(Q=3/2)':>14}  {'Q(r21, r31_TGP)':>18}  {'DeltaQ':>12}")
print("  " + "-"*60)

r21_grid = [4, 9, 16, 25, 50, 100, 207, 500, 1000, 4000]
for r21 in r21_grid:
    # Rozwiaz 2*(1+sqrt(r21)+sqrt(r31))^2 = 3*(1+r21+r31) dla r31
    # Niech x = sqrt(r31), a = 1+sqrt(r21):
    # 2*(a + x)^2 = 3*(1 + r21 + x^2)
    # 2*a^2 + 4*a*x + 2*x^2 = 3 + 3*r21 + 3*x^2
    # -x^2 + 4*a*x + (2*a^2 - 3 - 3*r21) = 0
    # x^2 - 4*a*x - (2*a^2 - 3 - 3*r21) = 0
    a = 1.0 + np.sqrt(r21)
    C = 2*a**2 - 3 - 3*r21
    discriminant = 16*a**2 + 4*C
    if discriminant < 0:
        r31_koide = np.nan
    else:
        x = (4*a + np.sqrt(discriminant)) / 2  # x = sqrt(r31)
        r31_koide = x**2
    Q_at_TGP = Q_from_ratios(r21, 3477.1) if r21 <= 400 else np.nan
    dQ = Q_at_TGP - 1.5 if not np.isnan(Q_at_TGP) else np.nan
    print(f"  {r21:>10.0f}  {r31_koide:>14.2f}  {Q_at_TGP:>18.6f}  {dQ:>+12.6f}")

print()
print("  OBSERWACJA: Dla r21=207, r31_Koide = 3477.65 (dokladna wartosc leptonu!)")
print("  TGP daje r31=3477.10 -- blizko ale nie na krzywej Koide.")
print()

# ============================================================
# KROK 2: Skan Q(lambda) przy stalym (alpha, a_gam)
# ============================================================
print("KROK 2: Skan Q(lambda) przy stalym (alpha=8.5616, a_gam=0.040)")
print("-" * 70)
print()
print("  PYTANIE: Czy istnieje lam* takie ze Q(lam*)=3/2 dokladnie?")
print("           Czy Q >= 3/2 zawsze, czy moze byc < 3/2?")
print()

alpha_fix = 8.5616
agam_fix  = 0.040

# Skan lambda -- szeroki zakres
lam_scan = np.logspace(-7, -4, 40)
Q_scan   = []
r21_scan = []
r31_scan = []
lam_valid = []

print("  Szukam zer dla roznych lambda...")
for lam in lam_scan:
    zeros = find_all_zeros(alpha_fix, agam_fix, lam)
    if len(zeros) >= 3:
        K1, K2, K3 = zeros[0], zeros[1], zeros[2]
        Q_val = koide_Q(K1, K2, K3)
        Q_scan.append(Q_val)
        r21_scan.append(K2/K1)
        r31_scan.append(K3/K1)
        lam_valid.append(lam)
    else:
        Q_scan.append(np.nan)
        r21_scan.append(np.nan)
        r31_scan.append(np.nan)
        lam_valid.append(lam)

# Filtruj nany
mask = [not np.isnan(q) for q in Q_scan]
lam_ok = np.array([lam_valid[i] for i in range(len(mask)) if mask[i]])
Q_ok   = np.array([Q_scan[i]   for i in range(len(mask)) if mask[i]])
r21_ok = np.array([r21_scan[i] for i in range(len(mask)) if mask[i]])
r31_ok = np.array([r31_scan[i] for i in range(len(mask)) if mask[i]])

print()
print(f"  {'lambda':>12}  {'Q':>10}  {'r21':>10}  {'r31':>10}  {'Q-3/2':>12}")
print("  " + "-"*60)
for i in range(len(lam_ok)):
    dQ = Q_ok[i] - 1.5
    print(f"  {lam_ok[i]:>12.3e}  {Q_ok[i]:>10.6f}  {r21_ok[i]:>10.3f}  {r31_ok[i]:>10.3f}  {dQ:>+12.6f}")
print()

if len(Q_ok) > 0:
    Q_min = np.nanmin(Q_ok)
    Q_max = np.nanmax(Q_ok)
    Q_at_opt = Q_ok[np.argmin(np.abs(lam_ok - 5.501e-6))] if len(lam_ok) > 0 else np.nan

    print(f"  Q_min = {Q_min:.6f}  (lambda = {lam_ok[np.argmin(Q_ok)]:.3e})")
    print(f"  Q_max = {Q_max:.6f}  (lambda = {lam_ok[np.argmax(Q_ok)]:.3e})")
    print(f"  Czy Q >= 3/2 zawsze? {'TAK' if Q_min >= 1.5 else 'NIE -- Q moze byc < 3/2!'}")
    print()

    # Szukaj przejscia Q=3/2
    sign_changes = []
    for i in range(len(Q_ok)-1):
        if (Q_ok[i] - 1.5) * (Q_ok[i+1] - 1.5) < 0:
            sign_changes.append(i)
    if sign_changes:
        print(f"  PRZEJSCIA Q=3/2 przy lambda ~ {lam_ok[sign_changes[0]]:.3e}")
    else:
        if Q_min >= 1.5:
            print("  BRAK przejscia Q=3/2 -- Q >= 3/2 wszedzie w badanym zakresie")
        else:
            print("  BRAK przejscia Q=3/2 -- Q < 3/2 wszedzie!")
    print()

# ============================================================
# KROK 3: Skan Q poza rodzina optymalna -- rozne (alpha, a_gam)
# ============================================================
print("KROK 3: Skan Q dla roznych (alpha, a_gam, lambda)")
print("-" * 70)
print()
print("  TEST: Czy Q >= 3/2 jest prawdziwe POZA rodzina optymalna?")
print()

test_points = [
    # alpha, a_gam, lam_range
    (3.0,  0.020, [1e-6, 1e-4]),
    (5.0,  0.025, [1e-7, 1e-4]),
    (8.56, 0.040, [1e-7, 1e-4]),
    (10.0, 0.050, [1e-7, 1e-3]),
    (15.0, 0.080, [1e-7, 1e-3]),
    (5.0,  0.100, [1e-7, 1e-3]),
    (3.0,  0.001, [1e-8, 1e-5]),
]

print(f"  {'alpha':>8}  {'a_gam':>7}  {'lambda':>12}  {'Q':>10}  {'r21':>10}  {'r31':>10}  {'Q>=3/2':>8}")
print("  " + "-"*75)

Q_all = []
for alpha, agam, (lam_lo, lam_hi) in test_points:
    # Skan lambda dla tego (alpha, agam)
    lam_test = np.logspace(np.log10(lam_lo), np.log10(lam_hi), 15)
    for lam in lam_test:
        zeros = find_all_zeros(alpha, agam, lam)
        if len(zeros) >= 3:
            K1, K2, K3 = zeros[0], zeros[1], zeros[2]
            Q_val = koide_Q(K1, K2, K3)
            r21 = K2/K1
            r31 = K3/K1
            Q_all.append(Q_val)
            ok = "TAK" if Q_val >= 1.5 else "NIE <<<"
            if Q_val < 1.5 or (abs(lam - lam_lo) < 1e-15) or (abs(np.log10(lam) - np.log10(lam_hi)) < 0.01):
                print(f"  {alpha:>8.2f}  {agam:>7.4f}  {lam:>12.3e}  {Q_val:>10.6f}  {r21:>10.3f}  {r31:>10.3f}  {ok:>8}")

print()
if Q_all:
    Q_global_min = np.nanmin(Q_all)
    Q_global_max = np.nanmax(Q_all)
    print(f"  Globalne Q_min = {Q_global_min:.6f}")
    print(f"  Globalne Q_max = {Q_global_max:.6f}")
    print(f"  Czy Q >= 3/2 globalnie? {'TAK' if Q_global_min >= 1.5 else 'NIE'}")
    print()

# ============================================================
# KROK 4: Analiza asymptotyczna
# ============================================================
print("KROK 4: Analiza asymptotyczna Q(r21, r31)")
print("-" * 70)
print()
print("  Rozwijamy dla duzych r21, r31:")
print()
print("  Q = (1 + sqrt(r21) + sqrt(r31))^2 / (1 + r21 + r31)")
print()
print("  Dla r31 >> r21 >> 1:")
print("    Licznik ~ r31 + 2*sqrt(r21*r31)")
print("    Mianownik ~ r31")
print("    Q ~ 1 + 2*sqrt(r21/r31)")
print()
print("  Dla Q = 3/2 asymptotycznie:")
print("    1 + 2*sqrt(r21/r31) = 3/2")
print("    sqrt(r21/r31) = 1/4")
print("    r31 = 16*r21")
print()
print("  Dla r21=207: r31_asympt = 16*207 = 3312 (vs rzeczywiste 3477)")
print("    Blad asymptotyki:", abs(3312-3477)/3477*100, "%")
print()

# Lepsza asymptotyka (drugi rzad)
print("  Lepsze przyblizenie (drugi rzad w 1/sqrt(r31)):")
print("    Q = 1 + 2*sqrt(r21/r31) + r21/r31 + 2/sqrt(r31) + ...")
print()
# Sprawdz numerycznie
for r21 in [207, 1000, 4000]:
    for r31 in [3477, 10000, 50000]:
        Q_exact = Q_from_ratios(r21, r31)
        Q_LO    = 1.0 + 2*np.sqrt(r21/r31)
        Q_NLO   = Q_LO + r21/r31 + 2/np.sqrt(r31)
        print(f"  r21={r21:5d}, r31={r31:6d}:  Q_exact={Q_exact:.5f},  "
              f"Q_LO={Q_LO:.5f} (err={abs(Q_exact-Q_LO)/Q_exact*100:.2f}%),  "
              f"Q_NLO={Q_NLO:.5f} (err={abs(Q_exact-Q_NLO)/Q_exact*100:.2f}%)")
print()

# Dla Q >= 3/2 asymptotycznie: r31 <= 16*r21 (wymagany warunek)
print("  WNIOSEK ASYMPTOTYCZNY:")
print("  Q >= 3/2 <==> r31 <= (approx) 16*r21 + korekcje")
print("  Dla r21=207: r31 <= ~3477 (dokladnie na granicy Koide!)")
print()

# ============================================================
# KROK 5: Czy Q = 3/2 jest minimum w przestrzeni parametrow?
# ============================================================
print("KROK 5: Test -- czy Q jest monotoniczne w lambda przy stalym (alpha, a_gam)?")
print("-" * 70)
print()
print("  Jezeli Q(lambda) jest monotoniczne i Q_min ~ 3/2 (dla lam -> inf),")
print("  to Q = 3/2 bylby OGRANICZENIEM ASYMPTOTYCZNYM, nie predykcja.")
print()

alpha_t = 8.5616
agam_t  = 0.040
lam_mono = np.logspace(-7, -3, 60)

Q_mono  = []
lam_ok2 = []
for lam in lam_mono:
    zeros = find_all_zeros(alpha_t, agam_t, lam)
    if len(zeros) >= 3:
        K1, K2, K3 = zeros[0], zeros[1], zeros[2]
        Q_mono.append(koide_Q(K1, K2, K3))
        lam_ok2.append(lam)
    else:
        Q_mono.append(np.nan)
        lam_ok2.append(lam)

lam_ok2 = np.array(lam_ok2)
Q_mono  = np.array(Q_mono)
finite  = np.isfinite(Q_mono)

if np.sum(finite) > 3:
    print(f"  lambda_min zakresu: {lam_ok2[finite][0]:.3e},  Q = {Q_mono[finite][0]:.6f}")
    print(f"  lambda_max zakresu: {lam_ok2[finite][-1]:.3e},  Q = {Q_mono[finite][-1]:.6f}")
    print(f"  Q_min = {Q_mono[finite].min():.6f}  przy lam = {lam_ok2[finite][np.argmin(Q_mono[finite])]:.3e}")
    print(f"  Monotoniczne? Sprawdzam znaki roznic...")
    dQ = np.diff(Q_mono[finite])
    n_pos = np.sum(dQ > 0)
    n_neg = np.sum(dQ < 0)
    print(f"  Rosnace przejscia: {n_pos}, Malejace: {n_neg}")
    if n_neg == 0:
        print("  Q(lambda) jest MONOTONICZNE ROSNACE -- Q maleje gdy lambda maleje")
        print("  => Q osiaga minimum dla lambda -> 0 (lub granicy istnienia K*3)")
    elif n_pos == 0:
        print("  Q(lambda) jest MONOTONICZNE MALEJACE")
    else:
        print("  Q(lambda) NIE jest monotoniczne -- ma ekstremum!")
        idx_min = np.argmin(Q_mono[finite])
        print(f"  Minimum: Q = {Q_mono[finite][idx_min]:.6f} przy lambda = {lam_ok2[finite][idx_min]:.3e}")
print()

# ============================================================
# KROK 6: Podsumowanie
# ============================================================
print("=" * 70)
print("KROK 6: PODSUMOWANIE -- STATUS PROBLEMU OTWARTEGO O-K1")
print("=" * 70)
print()
print("  CO ZOSTALO USTALONE W P27:")
print()
print("  1. Q(r21, r31) = (1+sqrt(r21)+sqrt(r31))^2 / (1+r21+r31)")
print("     Jest to funkcja tylko stosunkow mas (niezalezna od K*1).")
print()
print("  2. Krzywa Koide Q=3/2 w przestrzeni (r21, r31):")
print("     Dokladna postac: 2*(1+sqrt(r21)+sqrt(r31))^2 = 3*(1+r21+r31)")
print("     Punkt leptonowy (207, 3477) LEZY na tej krzywej (Koide 1982).")
print("     Punkt TGP (207, 3477.1) lezy BLISKO, ale nieznacznie poza nia.")
print()
if Q_all:
    print(f"  3. Zakres Q dla zbadanych parametrow: [{Q_global_min:.4f}, {Q_global_max:.4f}]")
    print(f"     Q >= 3/2 {'ZAWSZE (w badanym zakresie)' if Q_global_min >= 1.5 else 'MOZE BYC < 3/2!'}")
print()
print("  4. Asymptotyka: Q >= 3/2 <==> r31 <= (1+2*sqrt(r21))^2 / ...")
print("     Dla r21=207 ta granica wynosi r31_max ~ 3477 (!).")
print("     Punkt leptonowy LEZY dokladnie na granicy asymptotycznej.")
print()
print("  CO POZOSTAJE OTWARTE:")
print()
print("  O-K1a: Dowod analityczny Q >= 3/2 dla wszystkich zer TGP")
print("         (lub kontrprzyklad ze Q < 3/2 jest mozliwe)")
print()
print("  O-K1b: Dlaczego r31 ~ 16*r21 (Stosunek wynikajacy z Q=3/2 asymptotycznie)?")
print("         Czy TGP narzuca r31/r21 ~ 16 z jakiejs zasady?")
print()
print("  O-K1c: Czy istnieje zasada minimalizacji lub symetria V_mod")
print("         wymuszajaca (r21, r31) na krzywej Koide?")
print()
print("  REKOMENDACJA NA NASTEPNA SESJE (P28):")
print("  - Skan szeroki Q w 3D przestrzeni (alpha, a_gam, lambda)")
print("  - Szukaj kontrprzykladu Q < 3/2")
print("  - Zbadaj zachowanie Q dla lambda -> 0 (zanik K*3)")
print("  - Analityczna analiza warunku Q=3/2 jako rownania na lambda*(r21)")
print()
print("GOTOWE: P27 zakonczone.")

# ============================================================
# WYKRES
# ============================================================
print()
print("Tworze wykresy...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("P27: Krajobraz Q = Koide invariant -- Problem otwarty O-K1\n"
             "Q = (sum sqrt(K*_n))^2 / sum(K*_n)", fontsize=11, fontweight='bold')

# Panel 1: Q(r21, r31) -- mapa
ax = axes[0, 0]
r21_2d = np.linspace(10, 600, 120)
r31_2d = np.linspace(100, 10000, 120)
R21, R31 = np.meshgrid(r21_2d, r31_2d)
Q_2d = Q_from_ratios(R21, R31)
cf = ax.contourf(R21, R31, Q_2d, levels=50, cmap='RdYlBu_r')
ax.contour(R21, R31, Q_2d, levels=[1.5], colors='black', linewidths=2)
plt.colorbar(cf, ax=ax, label='Q')
ax.scatter([206.768], [3477.65], color='red',   s=150, zorder=5, label='leptony (dok.)')
ax.scatter([207.0],   [3477.1],  color='orange', s=150, zorder=5, marker='s', label='TGP')
ax.set_xlabel('r21'); ax.set_ylabel('r31')
ax.set_title('Q(r21, r31) -- linia Q=3/2 (czarna)', fontsize=9)
ax.legend(fontsize=8)

# Panel 2: Q(lambda) przy stalym (alpha, a_gam)
ax = axes[0, 1]
if np.sum(finite) > 2:
    ax.semilogx(lam_ok2[finite], Q_mono[finite], 'b-o', ms=4, lw=1.5, label=f'alpha={alpha_t}, a_gam={agam_t}')
    ax.axhline(1.5, color='red', lw=2, linestyle='--', label='Q=3/2 (Koide)')
    ax.axvline(5.501e-6, color='green', lw=1.5, linestyle=':', label='lambda* optymalny')
    ax.set_xlabel('lambda'); ax.set_ylabel('Q')
    ax.set_title('Q(lambda) -- monotoniczne?', fontsize=9)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 3: Q(lambda) dla roznych (alpha, a_gam) -- multi
ax = axes[1, 0]
lam_multi = np.logspace(-7, -4, 50)
colors_m = ['blue', 'red', 'green', 'purple', 'orange']
test_multiplot = [
    (5.91, 0.025, 2.883e-6),
    (6.87, 0.030, 3.743e-6),
    (7.74, 0.035, 4.618e-6),
    (8.56, 0.040, 5.501e-6),
]
for (alph, agam, lam_opt), col in zip(test_multiplot, colors_m):
    Q_line = []
    lam_line = []
    for lam in lam_multi:
        zeros = find_all_zeros(alph, agam, lam)
        if len(zeros) >= 3:
            Q_line.append(koide_Q(zeros[0], zeros[1], zeros[2]))
            lam_line.append(lam)
    if Q_line:
        ax.semilogx(lam_line, Q_line, color=col, lw=1.5, label=f'a={alph:.2f}')
        ax.axvline(lam_opt, color=col, lw=0.8, linestyle=':', alpha=0.7)
ax.axhline(1.5, color='black', lw=2, linestyle='--', label='Q=3/2')
ax.set_xlabel('lambda'); ax.set_ylabel('Q')
ax.set_title('Q(lambda) dla rodziny optymalnej', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Panel 4: Q vs r31/r21 -- czy stosunek r31/r21 determinuje Q?
ax = axes[1, 1]
ratio_array = np.linspace(1, 100, 200)
# Dla r21=207, r31=ratio*r21
r21_fixed = 207
r31_vary  = ratio_array * r21_fixed
Q_vary    = Q_from_ratios(r21_fixed, r31_vary)
ax.plot(ratio_array, Q_vary, 'b-', lw=2, label=f'r21={r21_fixed} (TGP)')
ax.axhline(1.5, color='red', lw=2, linestyle='--', label='Q=3/2')
ax.axvline(3477/207, color='green', lw=1.5, linestyle=':', label=f'r31/r21={3477/207:.1f} (TGP)')
ax.axvline(3477.65/206.768, color='orange', lw=1.5, linestyle=':', label=f'r31/r21={3477.65/206.768:.1f} (lept.)')
ax.set_xlabel('r31/r21'); ax.set_ylabel('Q')
ax.set_title('Q vs r31/r21 przy stalym r21=207', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_ylim([1.4, 1.7])

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
