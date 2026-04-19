"""
ps14_ambient_roomtemp_hunt.py  -  Program P6 final #

Cel:
  Zidentyfikowac realne drogi do AMBIENT room-temperature SC (T_c >= 293 K)
  uzywajac polaczonego modelu P6.A (d-wave cuprates) + P6.B (phonon coupling).

Strategia:
  1. Zebrac wszystkie istniejace rekordy ambient: Hg1223 138K, FeSe/STO 65K, MgB2 39K.
  2. Dla kazdej klasy materialow sformulowac "TGP boost path":
     - cuprates: wiekszy n_layers, optimal a
     - Fe-SC interface: silniejsza omega substrate
     - boron-rich: strain dla wyzszej omega
  3. Policzyc maksimum osiagalne w kazdej klasie.
  4. Ocenic realizm syntetyczny.

Wynik: tabelka konkretnych propozycji dla roomtemp SC ambient.
"""

import numpy as np

K_B = 8.617333e-5
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A
C_0 = 48.8222
sigma_a = 2.5856

A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

# P6.A (z ps10)
K_dw = 3.498
Lambda_E_cup_base = 0.0513  # meV
A_d = 0.3096
A_p = 0.2067
A_ZR_sq = A_d**2 + 2.0 * A_p**2

# P6.B (z ps12)
alpha_P6B = 1.04
Lambda_0_P6B = 0.0962
omega_0 = 15.0

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma_a**2)

# =============================================================
# Model A: P6.A cuprate d-wave (bez phonon booster)
# =============================================================

def Tc_cuprate_P6A(a_A, n_layers, omega_phonon_cup=55.0):
    """P6.A z opcjonalnym P6.B booster dla cupratow."""
    layer_factor = n_layers ** 0.5
    M = M_gauss(a_A)
    # Dla cupratow: Lambda ma inny skalowanie (electron. pairing, slaby phonon effect)
    # Uzyj alfa_cup = 0.3 (sub-linearne)
    alpha_cup = 0.30
    boost = (omega_phonon_cup / omega_0) ** alpha_cup
    Lambda_eff = Lambda_E_cup_base * boost
    Tc_substr = K_dw * k_d(8) * C_0 * A_ZR_sq * M * layer_factor
    return Tc_substr * (Lambda_eff * 1e-3) / K_B

# =============================================================
# Model B: P6.B phonon coupled (dla hydrydow, interface, MgB2)
# =============================================================

def Tc_P6B(a_A, orb, z, omega_phonon):
    A = A_map[orb]
    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega_phonon / omega_0) ** alpha_P6B
    Lambda_eff = Lambda_0_P6B * boost
    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_eff * 1e-3) / K_B

# =============================================================
# Zestawienie istniejacych rekordow ambient
# =============================================================

print("=" * 78)
print("  ps14_ambient_roomtemp_hunt.py")
print("=" * 78)
print()
print("=" * 78)
print("  Part A. Istniejace rekordy ambient + TGP interpretacja")
print("=" * 78)
print()

existing = [
    # (nazwa, T_obs, mechanizm TGP, parametry, komentarz)
    ("Hg1223",   138.0, "P6.A",  "a=3.855, n=3",              "Rekord ambient od 1993"),
    ("HgBaCa3Cu4O", 125.0, "P6.A", "a=3.855, n=4",             "Hg1234 rzadka synteza"),
    ("BiSCCO",   110.0, "P6.A",  "a=3.820, n=2",              "powszechny"),
    ("YBCO",     92.0,  "P6.A",  "a=3.820, n=2",              "komercyjny"),
    ("FeSe/STO", 65.0,  "P6.B",  "a=3.770, omega_STO=80",     "monolayer UHV"),
    ("LSCO",     38.0,  "P6.A",  "a=3.780, n=1",              "domieszka opt."),
    ("MgB2",     39.0,  "P6.B",  "a=3.086, omega=75",         "B-B ambient"),
    ("NdFeAsO-F",55.0,  "mixed", "a=3.970, 1111",             "optymalny 1111"),
    ("K3C60",    19.0,  "P6.B",  "a=14.24, fulleride",        "intercalated"),
    ("CaC6",     11.5,  "P6.B",  "a=2.610, omega=140",        "graphite intercal."),
]

print(f"  {'Material':>12} {'T_obs':>7} {'Mech':>6}  {'Parametry':>30}  {'Komentarz':>25}")
print(f"  {'------------':>12} {'-------':>7} {'------':>6}  {'------------------------------':>30}  {'-------------------------':>25}")
for name, Tobs, mech, params, comment in existing:
    print(f"  {name:>12} {Tobs:>7.1f} {mech:>6}  {params:>30}  {comment[:25]:>25}")

# =============================================================
# Part B. TGP kandydaci "realistyczni" na ambient high-Tc
# =============================================================

print()
print("=" * 78)
print("  Part B. TGP predykcje dla ambient high-Tc (realizowalne)")
print("=" * 78)
print()
print("  Kategoria: materialy, ktorych synteza jest w zasieglu obecnej technologii.")
print()

realistic = [
    # (nazwa, T_pred_func, T_pred_val, droga syntezy)
    ("Hg1245 (n=5)",                Tc_cuprate_P6A(3.855, 5), "high-P synthesis 3-5 GPa, quench"),
    ("Hg1234 / SrTiO3 epitaxial",   Tc_cuprate_P6A(3.905, 4), "MBE/PLD na STO substrate"),
    ("Hg1223 / NdScO3 strain",      Tc_cuprate_P6A(4.005, 3), "silny strain a=4.005"),
    ("Hg cuprate n=7 hypothet.",    Tc_cuprate_P6A(3.855, 7), "super-stack cuprate"),
    ("FeSe / BaTiO3",               Tc_P6B(3.770, "d", 8, 110), "BTO omega_Fuchs~110 meV"),
    ("FeSe / LaAlO3",               Tc_P6B(3.790, "d", 8, 120), "LAO wyzsza omega"),
    ("FeSe monolayer na diamond",   Tc_P6B(3.770, "d", 8, 165), "diament omega~165meV"),
    ("MgB2 / SrTiO3 interface",     Tc_P6B(3.086, "sp",8, 150), "Mg-B-STO heterostruktura"),
    ("MgB2 / diament substrate",    Tc_P6B(3.086, "sp",8, 165), "C-B interface"),
    ("LiC6 ultra-doped",            Tc_P6B(2.460, "sp",6, 200), "Li-graphene wzmocniony"),
    ("Bilayer graphene / BN",       Tc_P6B(2.460, "sp",6, 170), "vdW heterostruktura"),
    ("Ca-doped graphite multilayer",Tc_P6B(2.610, "sp",6, 160), "doped grafit"),
]

print(f"  {'Kandydat':>32}  {'T_pred[K]':>9}  {'Droga syntezy':>40}")
print(f"  {'--------------------------------':>32}  {'---------':>9}  {'----------------------------------------':>40}")
for name, Tp, route in realistic:
    marker = " !!" if Tp > 150 else ""
    print(f"  {name[:32]:>32}  {Tp:>9.1f}  {route[:40]:>40}{marker}")

# =============================================================
# Part C. Spekulatywne - roomtemp range
# =============================================================

print()
print("=" * 78)
print("  Part C. Spekulatywne - czy TGP dopuszcza T_c = 293 K przy ambient?")
print("=" * 78)
print()

speculative = [
    # (nazwa, a, params, T_pred, rationale)
    ("Hg1245 na SrTiO3",               Tc_cuprate_P6A(3.905, 5), "P6.A z strain"),
    ("Hg cuprate n=10 (infinite-like)",Tc_cuprate_P6A(3.855, 10), "hipotetyczny stack"),
    ("Hg n=5 na DyScO3 a=3.946",       Tc_cuprate_P6A(3.946, 5), "strain extreme"),
    ("Hg n=7 na substracie wysokim",   Tc_cuprate_P6A(4.050, 7), "combined strain+stack"),
    ("Perfect cuprate a=4.088, n=5",   Tc_cuprate_P6A(4.088, 5), "idealny w harmonice"),
    ("Perfect cuprate a=4.088, n=10",  Tc_cuprate_P6A(4.088, 10), "idealny + infinite-layer"),

    # P6.B drogi
    ("Metallic Li-H ambient",          Tc_P6B(3.500, "sp",8, 250), "spekulowana metalizacja"),
    ("BC4 / boron-carbon metal",       Tc_P6B(2.700, "sp",6, 250), "hypothetical"),
    ("Diament doped + boron",          Tc_P6B(3.567, "sp",8, 165), "metallic diamond"),
    ("Li-BC3 intercal.",               Tc_P6B(2.560, "sp",6, 210), "eksperymentalnie"),
    ("Hypothetical C-H metalliczna",   Tc_P6B(3.50,  "sp",8, 300), "metastabilne"),

    # Kombinowane P6.A + P6.B
    ("Hg1245 na interface wysoki omega", Tc_cuprate_P6A(3.855, 5, 120), "cuprate + phonon boost"),
    ("Hg cuprate + FeSe/STO-like interf.",Tc_cuprate_P6A(3.855, 3, 150), "dual mechanism"),
]

print(f"  {'Scenariusz':>38}  {'T_pred[K]':>9}  {'Rationale':>25}")
print(f"  {'--------------------------------------':>38}  {'---------':>9}  {'-------------------------':>25}")
for name, Tp, rationale in speculative:
    marker = " <<< ROOMTEMP" if Tp >= 293 else (" <<< >200K" if Tp >= 200 else "")
    print(f"  {name[:38]:>38}  {Tp:>9.1f}  {rationale[:25]:>25}{marker}")

# =============================================================
# Part D. Skan parametrow: mapa "T_c achievable"
# =============================================================

print()
print("=" * 78)
print("  Part D. Mapa parametrow dla cupratow (P6.A, boost omega_phonon)")
print("=" * 78)
print()

print("  Tabela T_c(n_layers, a_inplane, omega_cup):")
print()
print("  omega_cup = 55 meV (std cuprate)")
print(f"  {'a/n':>7}", end="")
for n in [1, 2, 3, 5, 7, 10, 15]:
    print(f" {n:>6}", end="")
print()
for a in [3.80, 3.85, 3.90, 4.00, 4.088, 4.20]:
    print(f"  {a:>7.3f}", end="")
    for n in [1, 2, 3, 5, 7, 10, 15]:
        Tp = Tc_cuprate_P6A(a, n, 55)
        print(f" {Tp:>6.1f}", end="")
    print()

print()
print("  omega_cup = 150 meV (interface engineered)")
print(f"  {'a/n':>7}", end="")
for n in [1, 2, 3, 5, 7, 10, 15]:
    print(f" {n:>6}", end="")
print()
for a in [3.80, 3.85, 3.90, 4.00, 4.088, 4.20]:
    print(f"  {a:>7.3f}", end="")
    for n in [1, 2, 3, 5, 7, 10, 15]:
        Tp = Tc_cuprate_P6A(a, n, 150)
        print(f" {Tp:>6.1f}", end="")
    print()

# =============================================================
# Part E. Tabela "konkretne eksperymenty"
# =============================================================

print()
print("=" * 78)
print("  Part E. KONKRETNE EKSPERYMENTALNE PROPOZYCJE TGP")
print("=" * 78)
print()

experiments = [
    # (rank, nazwa, setup, expected_Tc, trudnosc, impact)
    (1, "Hg1245 na SrTiO3 (a=3.905)",
        "MBE wzrost film 20nm + oxygen annealing",
        Tc_cuprate_P6A(3.905, 5),
        "SREDNIA: Hg toksyczny, synteza n=5 nieustabilizowana",
        "HIGH: pobija Hg1223 rekord, test P6.A + strain"),
    (2, "FeSe monolayer / BaTiO3",
        "MBE FeSe 1ML na BTO z precyzyjnym tenzylnym strain",
        Tc_P6B(3.770, "d", 8, 110),
        "SREDNIA: znana technika dla FeSe/STO",
        "HIGH: testuje P6.B omega scaling przez nowy substrate"),
    (3, "Hg1223 / diament substrate",
        "CVD diament + epitaxial Hg1223 film",
        Tc_cuprate_P6A(3.855, 3, 165),
        "WYSOKA: niezgodnosc strukturalna Hg-cuprate vs diament",
        "EXTREME: test dual P6.A+P6.B mechanism"),
    (4, "MgB2 / SrTiO3 interface",
        "MgB2 thin film na STO z oksydowana warstwa posrednia",
        Tc_P6B(3.086, "sp", 8, 150),
        "SREDNIA: MgB2 thin films juz zrobione",
        "MEDIUM: test P6.B universalnosci poza Fe-SC"),
    (5, "Hg1223 na NdScO3 (a=4.005)",
        "specialne substrat NdScO3 dla optimal a",
        Tc_cuprate_P6A(4.005, 3),
        "SREDNIA: NdScO3 dostepny crystal-grown",
        "HIGH: maksymalne zblizenie do a*_TGP harmonika"),
    (6, "Li-BC3 syntetyczny",
        "solid-state synthesis + post Li-intercalation",
        Tc_P6B(2.560, "sp", 6, 210),
        "NOWA: niesyntetyzowany material",
        "MEDIUM: test TGP dla layered non-cuprate metal"),
]

for rank, name, setup, Tp, trud, impact in experiments:
    print(f"  [{rank}] {name}")
    print(f"      Setup: {setup}")
    print(f"      T_c pred (TGP) = {Tp:.1f} K")
    print(f"      Trudnosc:       {trud}")
    print(f"      Impact:         {impact}")
    print()

# =============================================================
# Werdykt koncowy
# =============================================================
print("=" * 78)
print("  Werdykt koncowy ps14")
print("=" * 78)
print()
print("  CO TGP MOZE DO UZNANIA:")
print("    - ambient T_c do ~150 K (cuprate multilayer, a optimal)")
print("    - interface-boost T_c 100-200 K (FeSe/STO + engineered substrate)")
print("    - MgB2-klasa ambient do ~80 K (strain + strong-coupling)")
print()
print("  ROOM-TEMPERATURE AMBIENT (T_c >= 293 K):")
print("    - Matematycznie TAK (vide ps13: 699K max przy omega=400meV, ideal a/orb/z)")
print("    - Realistycznie: WYMAGA")
print("        * metaliczny crystal z omega_phonon >= 200 meV")
print("        * a_cell w granicach +/-0.1 A od a*_TGP = 4.088")
print("        * orbital d lub hybryd d+p")
print("        * struktura layered n>=5 (dla cuprate-like)")
print("    - Kandydaci HIPOTETYCZNI:")
print("        * Li-H metalliczny przy ambient (ciagle prawnie metastable)")
print("        * boron-carbon metallic layered Li-BC3")
print("        * cuprate-inf-layer + ultra-strain substrate")
print()
print("  PROZNIA (UHV):")
print("    - NIEZBEDNA dla monolayer (FeSe/STO, nickelates)")
print("    - stabilizuje metastabilne fazy wysokotempowe")
print("    - Nie podnosi T_c bezposrednio, ale umozliwia istnienie fazy 'high-T_c'")
print()
print("  NAJBARDZIEJ OBIECUJACY KIERUNEK EKSPERYMENTALNY:")
print("    Hg1245 na SrTiO3 substrate (combine n=5 + optimal strain)")
print(f"    Przewidywany T_c (TGP P6.A): {Tc_cuprate_P6A(3.905, 5):.1f} K")
print(f"    Jesli P6.B interface coupling: {Tc_cuprate_P6A(3.905, 5, 150):.1f} K")
print()
print("=" * 78)
print("  ps14 complete.")
print("=" * 78)
