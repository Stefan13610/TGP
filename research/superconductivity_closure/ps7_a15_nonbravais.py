"""
ps7_a15_nonbravais.py  -  Program P5 #7

Cel:
  Struktura A15 (np. Nb3Sn, V3Si, Nb3Ge) NIE jest siatka Bravais!
  Baza zawiera 8 atomow (2 formuly jednostki) w komorce szescienne a.
  Atomy B (Sn, Si, Ge) tworza BCC (a, z=8).
  Atomy A (Nb, V) tworza 3 prostopadle lancuchy wzdluz krawedzi,
     z odstepem wewnatrz-lancucha = a/2.

  Dla SC wazne sa lancuchy A -> efektywna odleglosc soliton-soliton
     w lancuchu = a/2, a liczba sasiadow w lancuchu = z_chain = 2.
     Ale calkowita k-koordynacja dla A jest hybrydowa.

  W ps5 traktowalismy A15 jak zwykly BCC z a=pełna komorka,
     z=12 (nierealistycznie wysokie - ale dalo dobry fit).
     Teraz sprawdzamy LAN'CUCHY A15 jako rzeczywisty mechanizm.

Plan:
  1. Zmodeluj A15 jako superpozycje LANCUCHOW (z_chain=2, a_chain = a/2)
     plus slabsze sprzezenie 3D (z_bcc = 8, a_bcc = a).
  2. Fit T_c z ANISOTROPOWA formula:
        T_c = (w_chain * k_d(2) * J(a/2) + w_3d * k_d(8) * J(a)) * Lambda_E
  3. Porownaj z ps5 (5c all-free) i z BCC-approx.

Fixed wejscia:
  Soliton amplitudy z P4 (czy fit z ps5 5c).
  Parametry C_0, a_star_tgp.

Uwaga: to test KONCEPTUALNY, nie global fit - sprawdzamy czy A15
  ma sens mechanistyczny, nie tylko numerycznie pasuje.
"""

import numpy as np

# =============================================================
# Stale
# =============================================================

K_B = 8.617333e-5     # eV/K (zgodnie z ps5)
A_BOHR_ANGSTROM = 0.52917721067

A_e_P4  = 0.124587
A_mu_P4 = 0.472198
A_tau_P4 = 0.956027

# Wybor: ps5 5c all-free parametry
Lambda_E_5c = 0.1309  # meV
A_map_5c = {
    "s":  -0.1110,
    "sp":  0.2067,
    "d":   0.3096,
    "f":   2.0336,
}

C_0 = 48.8222
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # ~4.088 A
sigma_a_5c = 2.5856

# =============================================================
# k_d(z)
# =============================================================

def k_d_from_z(z):
    table = {2: 0.30, 4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}
    # z_chain=2: XY 1D nie ma ostrego SC, ale efektywnie niska stala porzadku
    return table.get(z, 2.936)

# =============================================================
# Modulacja gauss
# =============================================================

def nearest_tgp_harmonic_delta(a_A):
    n = round(a_A / a_star_tgp_A)
    n = max(1, n)
    return a_A - n * a_star_tgp_A

def M_gauss(a_A, sigma_a):
    d = nearest_tgp_harmonic_delta(a_A)
    return np.exp(-d**2 / sigma_a**2)

# =============================================================
# A15 model: hybryda lancuch + 3D
# =============================================================

def Tc_bravais(a_A, orb, z, Lambda_E, A_map, sigma_a):
    """Zwykly model Bravais (ps5)."""
    A = A_map[orb]
    J = C_0 * A**2
    M = M_gauss(a_A, sigma_a)
    Tc_substr = k_d_from_z(z) * J * M
    return Tc_substr * (Lambda_E * 1e-3) / K_B

def Tc_a15(a_A, orb, Lambda_E, A_map, sigma_a, w_chain=0.6):
    """Hybryda lancuch (z=2, a/2) + 3D ramka (z=8, a)."""
    A = A_map[orb]
    J = C_0 * A**2
    # Lancuch: odstep a/2 jest bardzo blisko a_star_tgp_A
    a_chain = a_A / 2.0
    M_chain = M_gauss(a_chain, sigma_a)
    # 3D ramka
    M_3d = M_gauss(a_A, sigma_a)
    # Superpozycja: wartosc kluczowa to modulacja lancucha
    w_3d = 1.0 - w_chain
    Tc_substr = (w_chain * k_d_from_z(2) * J * M_chain +
                 w_3d * k_d_from_z(8) * J * M_3d)
    return Tc_substr * (Lambda_E * 1e-3) / K_B

# =============================================================
# A15 materialy (i non-A15 dla kontroli)
# =============================================================

# (nazwa, a[A], orb, struktura, Tc_obs)
a15_materials = [
    ("V3Si",  4.725, "d", "A15", 17.10),
    ("Nb3Sn", 5.290, "d", "A15", 18.30),
    ("Nb3Ge", 5.140, "d", "A15", 23.20),
    ("Nb3Al", 5.187, "d", "A15", 18.80),  # kontrola literaturowa
    ("V3Ga",  4.817, "d", "A15", 14.80),
    ("Mo3Ge", 4.934, "d", "A15",  1.50),  # niski T_c
    ("Cr3Si", 4.540, "d", "A15",  0.45),  # bardzo niski
]

# Kontrola: Nb metaliczne BCC
# Nb a=3.301, z=8, T_c=9.26

print("=" * 78)
print("  ps7_a15_nonbravais.py")
print("=" * 78)
print()
print("  Test hipotezy: A15 jako hybryda lancuch (a/2, z=2) + 3D (a, z=8).")
print("  w_chain = 0.6  (wiekszy wklad lancucha - fizycznie: anizotropia Nb/V)")
print()
print(f"  a*_TGP harmonika = {a_star_tgp_A:.4f} A")
print(f"  Parametry z ps5 5c: Lambda_E = {Lambda_E_5c} meV, sigma_a = {sigma_a_5c}")
print()

# =============================================================
# Porownanie trzech modeli
# =============================================================

print("=" * 78)
print("  Part A. Porownanie: Bravais-approx  vs  A15-chain-hybrid")
print("=" * 78)
print()
print(f"  {'Material':>9} {'a[A]':>7} {'T_obs':>7}  {'T_brav(z=12)':>12}  {'T_A15_hybrid':>12}")
print(f"  {'---------':>9} {'-------':>7} {'-------':>7}  {'------------':>12}  {'------------':>12}")

for name, a_A, orb, struct, Tc_obs in a15_materials:
    # Bravais-approx: jak ps5 (z=12, a=full cell)
    Tc_brav = Tc_bravais(a_A, orb, 12, Lambda_E_5c, A_map_5c, sigma_a_5c)
    # A15-hybrid
    Tc_a15_val = Tc_a15(a_A, orb, Lambda_E_5c, A_map_5c, sigma_a_5c, w_chain=0.6)
    print(f"  {name:>9} {a_A:>7.3f} {Tc_obs:>7.2f}  {Tc_brav:>12.2f}  {Tc_a15_val:>12.2f}")

# =============================================================
# Skanowanie w_chain
# =============================================================

print()
print("=" * 78)
print("  Part B. Skan w_chain (0 do 1) - minimalizacja RMS_log dla 7 A15")
print("=" * 78)
print()

def rms_a15(w_chain):
    resid = []
    for name, a_A, orb, struct, Tc_obs in a15_materials:
        Tc_p = Tc_a15(a_A, orb, Lambda_E_5c, A_map_5c, sigma_a_5c, w_chain)
        resid.append(np.log10(max(Tc_p, 1e-6)) - np.log10(Tc_obs))
    return float(np.sqrt(np.mean(np.array(resid)**2)))

print(f"  {'w_chain':>7}  {'RMS_log':>8}")
print(f"  {'-------':>7}  {'--------':>8}")
for w in np.linspace(0.0, 1.0, 11):
    r_w = rms_a15(w)
    print(f"  {w:>7.2f}  {r_w:>8.4f}")

# Znajdz optymalne
from scipy.optimize import minimize_scalar
res = minimize_scalar(rms_a15, bounds=(0.0, 1.0), method="bounded")
w_opt = res.x
rms_opt = res.fun
print()
print(f"  Optymalne w_chain = {w_opt:.4f},  RMS_log = {rms_opt:.4f}")

# =============================================================
# Korelacja a/2 z a_star_tgp
# =============================================================

print()
print("=" * 78)
print("  Part C. Dlaczego hybryda moze dzialac? Sprawdz a/2 vs a*_tgp")
print("=" * 78)
print()
print(f"  a*_TGP = {a_star_tgp_A:.4f} A")
print()
print(f"  {'Material':>9}  {'a/2':>7}  {'delta_a/2':>10}  {'M_cos(a/2)':>10}  {'M_cos(a)':>10}")
print(f"  {'---------':>9}  {'-------':>7}  {'----------':>10}  {'----------':>10}  {'----------':>10}")

for name, a_A, orb, struct, Tc_obs in a15_materials:
    a_half = a_A / 2.0
    d_half = nearest_tgp_harmonic_delta(a_half)
    M_half = np.exp(-d_half**2 / sigma_a_5c**2)
    d_full = nearest_tgp_harmonic_delta(a_A)
    M_full = np.exp(-d_full**2 / sigma_a_5c**2)
    print(f"  {name:>9}  {a_half:>7.3f}  {d_half:>+10.3f}  {M_half:>10.4f}  {M_full:>10.4f}")

# =============================================================
# Werdykt
# =============================================================

print()
print("=" * 78)
print("  Part D. Werdykt ps7")
print("=" * 78)
print()
# Porownaj z Bravais-approx A15
def rms_brav_a15():
    resid = []
    for name, a_A, orb, struct, Tc_obs in a15_materials:
        Tc_p = Tc_bravais(a_A, orb, 12, Lambda_E_5c, A_map_5c, sigma_a_5c)
        resid.append(np.log10(max(Tc_p, 1e-6)) - np.log10(Tc_obs))
    return float(np.sqrt(np.mean(np.array(resid)**2)))

rms_brav = rms_brav_a15()
print(f"  RMS_log (A15 jako BCC z=12):          {rms_brav:.4f}")
print(f"  RMS_log (A15 hybryda, w_opt={w_opt:.2f}): {rms_opt:.4f}")

if rms_opt < rms_brav * 0.9:
    print(f"\n  Hybryda LEPSZA ({rms_brav/rms_opt:.2f}x) -> A15 chain mechanism MA SENS.")
    print("  -> do P6: formalizacja lancuchowych struktur.")
elif rms_opt < rms_brav * 1.05:
    print(f"\n  Hybryda =  BCC-approx ({rms_brav/rms_opt:.2f}x).")
    print("  A15 numerycznie pasuje przez przypadek (2-a*_TGP dzielenie a).")
else:
    print(f"\n  Hybryda GORSZA niz Bravais-approx.")
    print("  Mechanizm chain nie jest glowny dla A15, tylko zbieznosc a/2 ~ a*_TGP.")

print()
print("  Obserwacja: we wszystkich A15, a/2 ~ 2.3-2.6 A, a a*_TGP = 4.09 A.")
print("  Stad lancuchowy odstep wcale nie jest blisko a*_TGP - to 'polharmonika'.")
print("  Gauss(d/sigma) z sigma~2.6 A daje istotna wartosc dla a/2 -> efekt realny,")
print("  ale nie fundamentalny.")

print()
print("=" * 78)
print("  ps7 complete.")
print("=" * 78)
