"""
coh00_baseline.py — tabela E_coh i wzorce rodzin chemicznych.

Cel: zebrać dane o energii kohezji metali (Kittel/CRC), kategoryzować,
znaleźć wzorce do testowania hipotez H1-H4 w coh01-coh04.
"""

import math
import sys
import io
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# Dane (Kittel 8th ed + CRC Handbook 97th)
# (symbol, Z, rodzina, E_coh [eV/atom], walencja, orbital_klasa, r_s [a_0])
# orbital_klasa: s=alkalie/II, sp=p-metal, d=d-block, f=f-block
# r_s: Wigner-Seitz radius w jednostkach a_0 (from n_e = 3/(4π r_s³·a_0³))

METALS = [
    # Alkalie (s-block, waln=1)
    ("Li",  3, "alkali",   1.63, 1, "s",  3.25),
    ("Na", 11, "alkali",   1.11, 1, "s",  3.93),
    ("K",  19, "alkali",   0.93, 1, "s",  4.86),
    ("Rb", 37, "alkali",   0.85, 1, "s",  5.20),
    ("Cs", 55, "alkali",   0.80, 1, "s",  5.62),
    # Metale ziem alkalicznych (s², walencja 2)
    ("Be",  4, "alk.earth", 3.32, 2, "s",  1.87),
    ("Mg", 12, "alk.earth", 1.51, 2, "s",  2.65),
    ("Ca", 20, "alk.earth", 1.84, 2, "s",  3.27),
    ("Sr", 38, "alk.earth", 1.72, 2, "s",  3.57),
    ("Ba", 56, "alk.earth", 1.90, 2, "s",  3.71),
    # Coinage (d¹⁰s¹, walencja 1+)
    ("Cu", 29, "coinage",   3.49, 1, "sd", 2.67),
    ("Ag", 47, "coinage",   2.95, 1, "sd", 3.02),
    ("Au", 79, "coinage",   3.81, 1, "sd", 3.01),
    # Metale przejściowe 3d
    ("Fe", 26, "3d",        4.28, None, "d", 2.12),
    ("Ni", 28, "3d",        4.44, None, "d", 2.07),
    ("Cr", 24, "3d",        4.10, None, "d", 2.00),
    ("Mn", 25, "3d",        2.92, None, "d", 2.14),
    ("Co", 27, "3d",        4.39, None, "d", 2.08),
    # Metale przejściowe 4d/5d ciężkie
    ("Mo", 42, "4d",        6.82, None, "d", 2.17),
    ("Nb", 41, "4d",        7.57, None, "d", 2.07),
    ("W",  74, "5d",        8.90, None, "d", 2.02),
    ("Pt", 78, "5d",        5.84, None, "d", 2.00),
    # sp-metale (p blok)
    ("Al", 13, "p-metal",   3.39, 3, "sp", 2.07),
    ("Ga", 31, "p-metal",   2.81, 3, "sp", 2.19),
    ("In", 49, "p-metal",   2.52, 3, "sp", 2.41),
    ("Sn", 50, "p-metal",   3.14, 4, "sp", 2.22),
    ("Pb", 82, "p-metal",   2.03, 4, "sp", 2.30),
    # Outliery
    ("Hg", 80, "p-metal",   0.67, 2, "sp", 2.65),
    ("Zn", 30, "3d",        1.35, 2, "sd", 2.31),
    ("Cd", 48, "4d",        1.16, 2, "sd", 2.59),
]

print("=" * 78)
print("  coh00 — baseline diagnostic dla energii kohezji metali")
print("=" * 78)

# Grupuj po rodzinie
families = {}
for m in METALS:
    fam = m[2]
    families.setdefault(fam, []).append(m)

print(f"\n[1] Dane: {len(METALS)} metali, {len(families)} rodzin")
print(f"\n{'Symb':>4}{'Z':>4}{'Rodzina':>12}{'E_coh':>8}{'Wal':>5}{'Orb':>5}{'r_s/a_0':>9}")
print("-" * 60)
for s, Z, fam, E, val, orb, rs in METALS:
    val_str = str(val) if val is not None else "–"
    print(f"{s:>4}{Z:>4}{fam:>12}{E:>8.2f}{val_str:>5}{orb:>5}{rs:>9.2f}")

# Statystyka per rodzina
print(f"\n[2] Statystyka per rodzina:")
print(f"{'Rodzina':<15}{'N':>4}{'mean':>8}{'std':>8}{'CV':>8}{'min–max':>14}")
print("-" * 56)
family_stats = {}
for fam, mets in families.items():
    Es = np.array([m[3] for m in mets])
    mean = Es.mean()
    std = Es.std()
    cv = std/mean
    family_stats[fam] = (mean, std, cv, Es.min(), Es.max())
    print(f"{fam:<15}{len(mets):>4}{mean:>8.2f}{std:>8.2f}{cv:>7.2%}{f'{Es.min():.2f}–{Es.max():.2f}':>14}")

# Alkalie: trend monotoniczny Li→Cs?
print(f"\n[3] Alkalie — trend monotoniczny?")
alkalies = [m for m in METALS if m[2] == "alkali"]
alkalies.sort(key=lambda m: m[1])  # by Z
for m in alkalies:
    print(f"    {m[0]:<3} Z={m[1]:>3}  E_coh = {m[3]:.2f} eV   r_s = {m[6]:.2f} a_0")

E_alk = np.array([m[3] for m in alkalies])
Z_alk = np.array([m[1] for m in alkalies])
rs_alk = np.array([m[6] for m in alkalies])

# Test monotoniczności
is_monotonic_decr = all(E_alk[i] > E_alk[i+1] for i in range(len(E_alk)-1))
print(f"    Li → Cs monotoniczny spadek E_coh: {is_monotonic_decr}")

# Korelacja E_coh vs 1/r_s (Fermi wealth scaling)
# E_F ∝ 1/r_s², Coulomb ∝ 1/r_s → E_coh ∝ 1/r_s (naive)
inv_rs = 1.0/rs_alk
slope, intercept = np.polyfit(inv_rs, E_alk, 1)
E_pred = slope*inv_rs + intercept
r2 = 1 - np.sum((E_alk - E_pred)**2) / np.sum((E_alk - E_alk.mean())**2)
print(f"    Fit E_coh = {slope:.3f}·(1/r_s) + {intercept:.3f} — r² = {r2:.4f}")

# Korelacja E_coh vs 1/r_s² (Fermi energy)
inv_rs2 = 1.0/rs_alk**2
slope2, intercept2 = np.polyfit(inv_rs2, E_alk, 1)
E_pred2 = slope2*inv_rs2 + intercept2
r2_2 = 1 - np.sum((E_alk - E_pred2)**2) / np.sum((E_alk - E_alk.mean())**2)
print(f"    Fit E_coh = {slope2:.3f}·(1/r_s²) + {intercept2:.3f} — r² = {r2_2:.4f}")

# ---------------------------------------------------------------------------
# [4] Koide-like dla alkali (wstępny H2)
# ---------------------------------------------------------------------------
print(f"\n[4] Koide-like dla alkali: K = ΣE/(Σ√E)²")

def koide_K(energies):
    E = np.array(energies)
    return np.sum(E) / np.sum(np.sqrt(E))**2

K_alkali_all = koide_K(E_alk)
K_alkali_3 = koide_K(E_alk[:3])  # Li,Na,K
K_leptons = 2.0/3.0

print(f"    K(Li,Na,K,Rb,Cs)   = {K_alkali_all:.5f}")
print(f"    K(Li,Na,K)         = {K_alkali_3:.5f}")
print(f"    K_leptons = 2/3    = {K_leptons:.5f}")
print(f"    K_3element_max     = 1.0 (1 dominant)")
print(f"    Uwaga: 5-element Koide przesunięcie ku niższym wartościom.")

# Per rodzina: K_coh
print(f"\n    K_coh per rodzina (3+ members):")
for fam, mets in families.items():
    if len(mets) >= 3:
        Es = [m[3] for m in mets]
        K = koide_K(Es)
        print(f"      {fam:<12} ({len(mets)} els): K = {K:.5f}  vs  2/3 = {K_leptons:.5f}")

# ---------------------------------------------------------------------------
# [5] Proste oszacowanie Fermi-sea dla alkali (H3 preview)
# ---------------------------------------------------------------------------
print(f"\n[5] Fermi-sea preview dla alkali (H3 test w coh03):")
# E_F = (ℏ²/2m)(3π²n_e)^(2/3), n_e = Z_val/V, V = (4/3)π r_s³
# W atomic units: E_F = 0.5·(9π/4)^(2/3) / r_s² × Z_val^(2/3)  [w Hartree]
# Wartości dla r_s ≈ 3-5 a_0: E_F = 0.5·(2.221)/r_s² = 1.11/r_s² Ha
# 1 Ha = 27.2114 eV
print(f"    {'Metal':<4}{'r_s':>6}{'E_F[eV]':>10}{'-e²/r_s[eV]':>14}{'sum[eV]':>10}{'E_coh obs':>12}")
for m in alkalies:
    rs = m[6]
    # E_F for free electron with n_e set by 1 electron per atom
    # E_F [Ha] = (1/2)·(3π²·n_e)^(2/3), n_e [a_0^-3] = 3/(4π·r_s³)
    # = (1/2)·(3π² · 3/(4π r_s³))^(2/3) = (1/2)·(9π/(4 r_s³))^(2/3)
    E_F_Ha = 0.5 * (9*math.pi/(4*rs**3))**(2/3)
    E_F_eV = E_F_Ha * 27.2114
    # Coulomb (Ewald+Madelung avg): ≈ -e²/(r_s·a_0·4πε₀) = -Hartree/r_s
    E_C_eV = -27.2114/rs
    E_sum = 0.6*E_F_eV + 0.9*E_C_eV   # kinetic 3/5·E_F + correlation fudge 0.9
    print(f"    {m[0]:<4}{rs:>6.2f}{E_F_eV:>10.2f}{E_C_eV:>14.2f}{E_sum:>10.2f}{m[3]:>12.2f}")

# ---------------------------------------------------------------------------
# [6] Hipoteza A_orb → E_coh (H1 preview) - potrzebne A_s, A_d z SC
# ---------------------------------------------------------------------------
print(f"\n[6] H1 preview: A_orb (SC) → E_coh korelacja?")
A_s  = -0.111  # z ps: A_s dla s-band alkalie
A_sp = +0.207  # dla sp-metali
A_d  = +0.310  # dla d-block
A_f  = +2.034  # dla f-block
print(f"    A_s (alkali)    = {A_s:+.3f}  |A|² = {A_s**2:.4f}")
print(f"    A_sp (p-metals) = {A_sp:+.3f}  |A|² = {A_sp**2:.4f}")
print(f"    A_d (d-block)   = {A_d:+.3f}  |A|² = {A_d**2:.4f}")
print(f"    A_f             = {A_f:+.3f}  |A|² = {A_f**2:.4f}")
print(f"    Średnie E_coh:")
for fam, (mean, _, _, _, _) in family_stats.items():
    print(f"      {fam:<12}: <E_coh> = {mean:.2f} eV")

# Czy A_s² vs <E_coh>_alkali daje coś sensownego?
# A_s² = 0.0123; E_coh_alkali avg = 1.06 eV → ratio ≈ 86 eV (bezmyślne skalowanie)
# A_d² = 0.096; E_coh_d avg ≈ 4 eV → ratio ≈ 42 eV
# To nie jest tak proste jak |A|²·Ry; może potrzeba innego skalowania.

# ---------------------------------------------------------------------------
# Werdyk coh00
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("  coh00 — BASELINE VERDICT")
print("=" * 78)
print(f"""
  STAN POCZĄTKOWY:
    • Zebrane: {len(METALS)} metali w {len(families)} rodzinach
    • Alkalie: monotoniczny spadek Li→Cs z rosnącym r_s ✓
    • E_coh(alkali) vs 1/r_s² (Fermi): r² = {r2_2:.3f}
    • E_coh(alkali) vs 1/r_s (Coulomb): r² = {r2:.3f}

  Wstępne obserwacje:
    H1: A_orb² (SC) vs <E_coh> rodziny nie ma prostej liniowej zależności.
        Potrzeba szerszego testu w coh01.
    H2: K_coh(alkali) = {K_alkali_all:.4f} — nie jest 2/3 ani blisko.
        Koide dla E_coh NIE JEST uniwersalną stałą dla rodzin.
    H3: Fermi-sea baseline (tabela wyżej) — rząd wielkości OK dla alkali,
        ale fudge factor wymagany.
    H4: Soliton TGP — wymaga osobnej analizy w coh04.

  NASTĘPNY KROK: coh01 (szczegółowy test H1) i coh02 (Koide dla rodzin).
  Ostrzeżenie: dane pokazują że tę fizykę można DOBRZE wyjaśnić standardowymi
  modelami (jellium + poprawki DFT). TGP musi dodać coś NAD standardem, albo
  po prostu reprodukować — ale wtedy nie jest „testem TGP".
""")
