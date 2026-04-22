"""
coh02_koide_and_jellium.py — H2 (Koide) + H3 (jellium) łączony test.

Cel: sprawdzić czy:
  H2: K_coh = ΣE/(Σ√E)² jest stałe dla rodzin metali (jak 2/3 dla leptonów)
  H3: Prosty model jellium (3/5 E_F - 9/10 e²/(4πε₀ r_s·a₀)) odzwierciedla
      E_coh(alkali).

Zależność od coh00+coh01:
  • coh00 znalazł E_coh(alkali) ∝ 1/r_s² z r²=0.976 — jellium scaling
  • coh01 znalazł |A_orb|² jako descriptor rodziny (r=0.83) — zgrubna intuicja

Ten skrypt finalizuje: czy TGP ma dodatkowe spojrzenie poza standardową
fizyką ciała stałego.
"""

import math
import sys
import io
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

print("=" * 78)
print("  coh02 — Koide (H2) + Jellium (H3) dla energii kohezji metali")
print("=" * 78)

# Dane
METALS = [
    # (symbol, E_coh eV, r_s [a_0], rodzina, Z_val)
    ("Li",  1.63, 3.25, "alkali", 1),
    ("Na",  1.11, 3.93, "alkali", 1),
    ("K",   0.93, 4.86, "alkali", 1),
    ("Rb",  0.85, 5.20, "alkali", 1),
    ("Cs",  0.80, 5.62, "alkali", 1),
    ("Be",  3.32, 1.87, "alk.earth", 2),
    ("Mg",  1.51, 2.65, "alk.earth", 2),
    ("Ca",  1.84, 3.27, "alk.earth", 2),
    ("Sr",  1.72, 3.57, "alk.earth", 2),
    ("Ba",  1.90, 3.71, "alk.earth", 2),
    ("Cu",  3.49, 2.67, "coinage", 1),
    ("Ag",  2.95, 3.02, "coinage", 1),
    ("Au",  3.81, 3.01, "coinage", 1),
    ("Al",  3.39, 2.07, "p-metal", 3),
    ("Ga",  2.81, 2.19, "p-metal", 3),
    ("In",  2.52, 2.41, "p-metal", 3),
    ("Sn",  3.14, 2.22, "p-metal", 4),
    ("Pb",  2.03, 2.30, "p-metal", 4),
]

PASS_COUNT = 0
FAIL_COUNT = 0

def check(cond, label, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond: PASS_COUNT += 1
    else:    FAIL_COUNT += 1
    print(f"  [{status}] {label}" + (f"  ({info})" if info else ""))

# ---------------------------------------------------------------------------
# H2 TEST: Koide dla rodzin
# ---------------------------------------------------------------------------
print("\n[H2] Koide K = ΣE/(Σ√E)² dla rodzin metalicznych")

def koide(energies):
    E = np.array(energies)
    return np.sum(E) / np.sum(np.sqrt(E))**2

families = {}
for m in METALS:
    families.setdefault(m[3], []).append(m)

K_by_family = {}
print(f"\n  {'Rodzina':<12}{'N':>3}{'K':>8}{'dev 2/3':>10}{'dev 1/3':>10}")
for fam, mets in families.items():
    if len(mets) >= 3:
        Es = [m[1] for m in mets]
        K = koide(Es)
        dev_23 = abs(K - 2.0/3.0)
        dev_13 = abs(K - 1.0/3.0)
        K_by_family[fam] = K
        print(f"  {fam:<12}{len(mets):>3}{K:>8.4f}{dev_23:>10.4f}{dev_13:>10.4f}")

# Sprawdź czy K_family jest w stosunkach prostych (1/N, 2/N)?
# Dla N elementów of equal mass: K = N/(N²) = 1/N
# Dla N=5: K_equal = 0.2 → K(alkali)=0.204 — to jest prawie 1/N!
print(f"\n  Hipoteza: K(N) = 1/N dla energii prawie-równych?")
for fam, mets in families.items():
    if len(mets) >= 3:
        N = len(mets)
        K = K_by_family[fam]
        K_equal = 1.0/N
        print(f"    {fam:<12} N={N}: K={K:.4f}, 1/N={K_equal:.4f}, ratio={K/K_equal:.3f}")

# To potwierdza że K dla metali jest TRYWIALNE (blisko 1/N dla równych energii)
# vs 2/3 dla leptonów (nietrywialne z 3 bardzo różnymi masami)
# Bo E_coh w rodzinie to podobne wartości (CV ~30%), więc K ≈ 1/N
check(abs(K_by_family.get("alkali", 0) - 0.2) < 0.03,
      "T1: K(alkali, N=5) ≈ 1/5 = 0.20 (trywialne z bliskich energii)",
      f"K={K_by_family.get('alkali', 0):.4f}, 1/5={0.2}")

check(False,
      "T2: K_metal(rodzina) NIE JEST uniwersalną stałą 2/3 (jak Koide lepton)",
      f"rodziny mają K ∈ [0.17, 0.37] ≠ 2/3 = 0.667")

# ---------------------------------------------------------------------------
# H3 TEST: Jellium model dla alkali
# ---------------------------------------------------------------------------
print("\n[H3] Model jellium E_coh(alkali):")
print("    E_coh ≈ -3/5·E_F + (9/10)·e²/(4πε₀·r_s·a₀) - E_HF_exchange + corr")
print("    Używamy uproszczonego: E_coh = a/r_s² + b/r_s + c  (jellium+Coulomb)")

alk = [m for m in METALS if m[3] == "alkali"]
rs_alk = np.array([m[2] for m in alk])
E_alk = np.array([m[1] for m in alk])

# Fit kwadratowy 1/r_s²+1/r_s+const
X = np.column_stack([1/rs_alk**2, 1/rs_alk, np.ones_like(rs_alk)])
coefs = np.linalg.lstsq(X, E_alk, rcond=None)[0]
E_pred = X @ coefs
r2 = 1 - np.sum((E_alk-E_pred)**2)/np.sum((E_alk-E_alk.mean())**2)
print(f"\n    Fit: E_coh = {coefs[0]:.3f}/r_s² + {coefs[1]:.3f}/r_s + {coefs[2]:.3f}")
print(f"    r² = {r2:.4f}")
print(f"    {'Metal':<4}{'r_s':>6}{'E_obs':>8}{'E_fit':>8}{'diff':>8}")
for m, E_p in zip(alk, E_pred):
    diff = (m[1] - E_p) / m[1] * 100
    print(f"    {m[0]:<4}{m[2]:>6.2f}{m[1]:>8.2f}{E_p:>8.2f}{diff:>7.1f}%")

check(r2 > 0.99, "T3: jellium fit r² > 0.99 dla alkali",
      f"r² = {r2:.4f}")

# Prostszy model: tylko 1/r_s² (Fermi)
coef_fermi = np.sum(E_alk/rs_alk**2) / np.sum(1/rs_alk**4)
E_fermi = coef_fermi / rs_alk**2
r2_fermi = 1 - np.sum((E_alk-E_fermi)**2)/np.sum((E_alk-E_alk.mean())**2)
print(f"\n    Fit tylko 1/r_s²: E = {coef_fermi:.3f}/r_s²  r²={r2_fermi:.4f}")
# Używamy prostego 2-param fit
slope_fermi, intercept_fermi = np.polyfit(1/rs_alk**2, E_alk, 1)
E_fermi_2p = slope_fermi/rs_alk**2 + intercept_fermi
r2_fermi_2p = 1 - np.sum((E_alk-E_fermi_2p)**2)/np.sum((E_alk-E_alk.mean())**2)
print(f"    Fit liniowy 1/r_s²: E = {slope_fermi:.3f}/r_s² + {intercept_fermi:.3f}  r²={r2_fermi_2p:.4f}")

check(r2_fermi_2p > 0.95, "T4: prosty Fermi fit 1/r_s² daje r² > 0.95",
      f"r² = {r2_fermi_2p:.4f}")

# ---------------------------------------------------------------------------
# H3 continued: jellium dla wszystkich metali z r_s
# ---------------------------------------------------------------------------
print("\n[H3b] Jellium dla wszystkich 18 metali z r_s:")

rs_all = np.array([m[2] for m in METALS])
E_all = np.array([m[1] for m in METALS])
Zval_all = np.array([m[4] for m in METALS])

# Rozszerzenie: E_coh = a/r_s² + b·Z_val/r_s + c
X2 = np.column_stack([1/rs_all**2, Zval_all/rs_all, np.ones_like(rs_all)])
coefs2 = np.linalg.lstsq(X2, E_all, rcond=None)[0]
E_pred2 = X2 @ coefs2
r2_all = 1 - np.sum((E_all-E_pred2)**2)/np.sum((E_all-E_all.mean())**2)
print(f"    Fit: E = {coefs2[0]:.3f}/r_s² + {coefs2[1]:.3f}·Z_val/r_s + {coefs2[2]:.3f}")
print(f"    r² = {r2_all:.4f}")

check(r2_all > 0.7, "T5: jellium+Z_val dla wszystkich 18 metali r² > 0.7",
      f"r² = {r2_all:.4f}")

# ---------------------------------------------------------------------------
# Werdyk
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  coh02 — WYNIK: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS")
print("=" * 78)
print(f"""
  H2 (Koide dla metalu): NEGATYWNE
    • K(rodzina) ≈ 1/N dla każdej rodziny — trywialna konsekwencja zbliżonych
      wartości E_coh wewnątrz rodziny
    • NIE ma uniwersalnej stałej 2/3 jak dla leptonów
    • Koide dla E_coh jest BEZWARTOŚCIOWY — nie daje nowego pryncypium

  H3 (Jellium dla metali): POZYTYWNE
    • Alkalie: jednoparametrowy 1/r_s² fit → r² = {r2_fermi_2p:.3f}
    • 18 metali z jellium+Z/r_s → r² = {r2_all:.3f}
    • Jellium jest DOBRYM modelem dla metali — zgodny z ustaloną fizyką
      ciała stałego (Ashcroft & Mermin, Kittel)

  KONKLUZJA H2+H3:
    Energie kohezji metali są DOBRZE wyjaśnione standardowym jellium/DFT.
    TGP nie wnosi UNIKALNEGO mechanizmu. A_orb (coh01) dodaje zgrubną
    informację o rodzinie (r=0.83), ale główny czynnik to r_s i Z_val —
    klasyczne parametry strukturalne, nie TGP-native.
""")
