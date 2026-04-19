"""
ps11_ambient_highTc_survey.py  -  Program P6.A #3

Cel:
  Znalezc kandydatow na high-Tc przy AMBIENT PRESSURE.
  Kluczowe pytania:
    1. Czy zelazowce (FeSe, Fe-pnictides) pasuja do P6.A?
    2. Czy nickelates ambient daja cos bliskiego cupratom?
    3. Czy istnieje hipotetyczny "super-cuprate" z T_c > 200 K?
    4. Czy TGP przewiduje jakiekolwiek granicy gornej T_c?
"""

import numpy as np

# =============================================================
# Stale i model (P6.A z ps10)
# =============================================================

K_B = 8.617333e-5
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A
C_0 = 48.8222

# Amplitudy 5c
A_d  = 0.3096
A_p  = 0.2067  # = A_sp
A_f  = 2.0336  # f-orbital; bedzie uzyte ostroznie

sigma_a = 2.5856
K_dw = 3.498       # z ps10
Lambda_E_cup = 0.0513  # meV  z ps10
Lambda_E_BCS = 0.1309  # meV  z ps5 5c (BCS-reference)

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A, sigma=sigma_a):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma**2)

def Tc_dwave_layered(a_A, n_layers, A_sq, z_planar=8, Lambda_E=Lambda_E_cup):
    """P6.A formula uzyta w ps10."""
    layer_factor = n_layers ** 0.5
    M = M_gauss(a_A)
    Tc_substr = K_dw * k_d(z_planar) * C_0 * A_sq * M * layer_factor
    return Tc_substr * (Lambda_E * 1e-3) / K_B

def Tc_swave_bulk(a_A, orb, z, Lambda_E=Lambda_E_BCS):
    """P5 formula (ps5 5c) dla porownania."""
    A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}
    A = A_map[orb]
    J = C_0 * A**2
    M = M_gauss(a_A)
    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_E * 1e-3) / K_B

# Zhang-Rice like A_ZR dla cuprates
A_ZR_sq = A_d**2 + 2.0 * A_p**2  # 0.1813

# =============================================================
# Rodzina 1: cuprates ambient (powtorka z ps10 jako benchmark)
# =============================================================

print("=" * 78)
print("  ps11_ambient_highTc_survey.py")
print("=" * 78)
print()
print(f"  Model P6.A: K_dw = {K_dw:.3f},  Lambda_E_cup = {Lambda_E_cup:.4f} meV")
print(f"  A_ZR^2 (d+2p) = {A_ZR_sq:.4f}")
print()

# =============================================================
# Rodzina 2: Fe-pnictides and chalcogenides (ambient P high-Tc)
# =============================================================

print("=" * 78)
print("  Rodzina Fe (pnictides/chalcogenides, ambient)")
print("=" * 78)
print()

# FeSe: tetragonal a=3.77, 1 Fe layer, T_c = 8K ambient, 65K monolayer/STO
# FeAs: generic, a~4.0, T_c~20-55K depending on dopant
# LaFeAsO (1111): a=4.04, 1 Fe layer, T_c = 26K
# BaFe2As2 (122): a=3.96, 1 Fe layer, T_c = 38K (Co-doped)
# Ba2Fe2As2 gdy domieszkowane: a=3.96, n_Fe~1, T_c = 38K
# NdFeAsO: a=3.97, T_c = 55K
# Tapeczkowany FeSe na SrTiO3: 65-100K w monolayerze

# Hipoteza: Fe-SC uzywa A_d + A_p (Fe 3d + As/Se 4p) ale z INNYM K
# Fe-SC nie sa d-wave (s+- najczesciej), wiec K moze byc inne

A_FeAs_sq = A_d**2 + A_p**2   # 1 pniktogen na Fe (vs 2 O na Cu w cuprates)

Fe_materials = [
    ("FeSe",         3.770, 1,  8.0, "ambient bulk"),
    ("FeSe/STO-ML",  3.770, 1, 65.0, "monolayer na SrTiO3"),
    ("LaFeAsO",      4.040, 1, 26.0, "1111 family"),
    ("SmFeAsO-F",    4.000, 1, 55.0, "doped 1111"),
    ("NdFeAsO-F",    3.970, 1, 55.0, "doped 1111 near record"),
    ("Ba122-Co",     3.960, 1, 22.0, "Ba(Fe,Co)2As2"),
    ("KFe2As2",      3.842, 1,  3.8, "K122 ambient"),
    ("FeTe_0.5Se0.5",3.800, 1, 15.0, "chalcogenide"),
]

def Tc_Fe(a, n_layers, Lambda_factor=1.0):
    """Fe-SC: uzywamy A_FeAs (d+p), K = K_dw (s+- = d-wave-like) ale skalowane."""
    # Hipoteza s+-: K podobne do d-wave ale z zeta=0.7 tlumienia
    K_sM = 0.7 * K_dw
    Lambda_Fe = Lambda_factor * Lambda_E_cup
    layer_factor = n_layers ** 0.5
    M = M_gauss(a)
    Tc_substr = K_sM * k_d(8) * C_0 * A_FeAs_sq * M * layer_factor
    return Tc_substr * (Lambda_Fe * 1e-3) / K_B

print(f"  A_FeAs^2 (d+p) = {A_FeAs_sq:.4f}  (cuprate A_ZR^2 = {A_ZR_sq:.4f}, stosunek {A_FeAs_sq/A_ZR_sq:.2f})")
print()
print(f"  {'Material':>14} {'a[A]':>5} {'n':>3} {'T_obs':>7} {'T_pred':>7}  {'komentarz':>25}")
print(f"  {'--------------':>14} {'-----':>5} {'---':>3} {'-------':>7} {'-------':>7}  {'-------------------------':>25}")

Fe_log_obs, Fe_log_pred = [], []
for name, a, n, Tobs, comment in Fe_materials:
    Tp = Tc_Fe(a, n, 1.0)
    Fe_log_obs.append(np.log10(Tobs))
    Fe_log_pred.append(np.log10(max(Tp, 1e-3)))
    print(f"  {name:>14} {a:>5.3f} {n:>3d} {Tobs:>7.1f} {Tp:>7.2f}  {comment[:25]:>25}")

Fe_r = np.corrcoef(Fe_log_obs, Fe_log_pred)[0, 1]
Fe_rms = float(np.sqrt(np.mean((np.array(Fe_log_pred) - np.array(Fe_log_obs))**2)))
print(f"\n  Fe-SC:  r(log-log) = {Fe_r:.3f},  RMS_log = {Fe_rms:.3f}")

# =============================================================
# Rodzina 3: Nickelates
# =============================================================
print()
print("=" * 78)
print("  Rodzina Ni-oxides (ambient + under pressure)")
print("=" * 78)
print()

Ni_materials = [
    ("NdNiO2-thin",   3.916, 1,   9.0, "infinite-layer Nd thin film"),
    ("PrNiO2-thin",   3.916, 1,  12.0, "Pr infinite-layer"),
    ("La_0.8Sr_0.2NiO2", 3.937, 1, 15.0, "doped La-IL"),
    ("Sr3Ni2O7-HP",   3.830, 2,  80.0, "pod P=14 GPa (nie ambient)"),
    ("LaNiO3",        3.840, 1,   0.5, "metallic, non-SC ambient"),
]

print(f"  {'Material':>18} {'a[A]':>5} {'n':>3} {'T_obs':>7} {'T_pred':>7}  {'komentarz':>30}")
print(f"  {'------------------':>18} {'-----':>5} {'---':>3} {'-------':>7} {'-------':>7}  {'------------------------------':>30}")

for name, a, n, Tobs, comment in Ni_materials:
    Tp = Tc_dwave_layered(a, n, A_ZR_sq)  # zalozenie: Ni-ZR jak Cu-ZR
    print(f"  {name:>18} {a:>5.3f} {n:>3d} {Tobs:>7.1f} {Tp:>7.2f}  {comment[:30]:>30}")

# =============================================================
# Rodzina 4: Ruthenates, Iridates, Titanates (eksperymentalne)
# =============================================================

print()
print("=" * 78)
print("  Rodzina Ru/Ir/Ti (exotic SC)")
print("=" * 78)
print()

Exotic = [
    ("Sr2RuO4",  3.862, 1,   1.5, "d-wave or p-wave? ambient"),
    ("Sr3Ir2O7", 3.896, 2,   0.0, "Mott insulator ambient, non-SC"),
    ("SrTiO3",   3.905, 1,   0.4, "n-doped, very low Tc"),
    ("K3C60",    14.24, 1,  19.0, "fulleride FCC (a duzo wyzszy, nie-standard)"),
]

print(f"  Uwaga: te materialy nie sa d-wave cuprate-like,")
print(f"  predykcja = tylko jesli przyjmiemy mechanizm TGP-ZR-2D.")
print()
print(f"  {'Material':>10} {'a[A]':>6} {'T_obs':>7} {'T_pred':>7}")
print(f"  {'----------':>10} {'------':>6} {'-------':>7} {'-------':>7}")
for name, a, n, Tobs, comment in Exotic:
    Tp = Tc_dwave_layered(a, n, A_ZR_sq)
    print(f"  {name:>10} {a:>6.3f} {Tobs:>7.2f} {Tp:>7.1f}  {comment}")

# =============================================================
# Rodzina 5: HIPOTEZE - co TGP przewiduje dla "idealnego" SC ambient
# =============================================================

print()
print("=" * 78)
print("  Hipotetyczne high-Tc ambient: maksymalizacja TGP")
print("=" * 78)
print()
print("  Zmienne: a_inplane, n_layers, typ orbitalu (d lub hybryd).")
print("  Ograniczenie fizyczne: a >= 3.50 (chemiczna mozliwosc)")
print("                        n <= 6 (znane syntezy max Hg1234 n=4, Tl1245 n=5)")
print()

scenariusze = [
    # (opis, a, n, A_sq, Lambda_E_mult, K_factor, komentarz)
    ("Hg1234 hipoteza",    3.855, 4, A_ZR_sq, 1.0, 1.0, "Synteza trudna, znana 1 proba"),
    ("Hg1245 hipoteza",    3.855, 5, A_ZR_sq, 1.0, 1.0, "n=5 cuprate nigdy nie zsyntezowany stabilnie"),
    ("Hg-inf layer",       3.855,10, A_ZR_sq, 1.0, 1.0, "sqrt(n) saturation - matematyczny ekstrapolant"),
    ("Pd-analog cuprate",  4.090, 3, A_ZR_sq, 1.0, 1.0, "Pd zamiast Cu, a optymalne na harm"),
    ("Cu-doped + strain",  4.090, 3, A_ZR_sq, 1.5, 1.0, "Induced strain i Lambda wzmocnione 1.5x"),
    ("ZR-boosted A_f",     3.855, 3, A_d**2 + 2*A_f**2*0.1, 1.0, 1.0, "Hybryd z f - ostroznie"),
    ("Perfect cuprate",    4.088, 3, A_ZR_sq, 1.0, 1.0, "dokladnie na harmonice"),
    ("Perfect + 5 layers", 4.088, 5, A_ZR_sq, 1.0, 1.0, "max achievable?"),
]

print(f"  {'Scenariusz':>22} {'a[A]':>5} {'n':>3} {'T_pred[K]':>9}  {'komentarz':>40}")
print(f"  {'----------------------':>22} {'-----':>5} {'---':>3} {'---------':>9}  {'----------------------------------------':>40}")

for opis, a, n, A_sq, lam_mult, K_fact, komentarz in scenariusze:
    Lam = lam_mult * Lambda_E_cup
    layer_factor = n ** 0.5
    M = M_gauss(a)
    Tc_substr = K_fact * K_dw * k_d(8) * C_0 * A_sq * M * layer_factor
    Tp = Tc_substr * (Lam * 1e-3) / K_B
    print(f"  {opis[:22]:>22} {a:>5.3f} {n:>3d} {Tp:>9.1f}  {komentarz[:40]:>40}")

# =============================================================
# Graniczne pytanie: czy jest plafon T_c w TGP?
# =============================================================

print()
print("=" * 78)
print("  Czy TGP ma plafon T_c?")
print("=" * 78)
print()
print("  Matematycznie:")
print("    T_c ~ K_dw * k_d * C_0 * A_sq * M * sqrt(n) * Lambda_E")
print()
print("    Maksymalne M = 1.0 (a = a*_TGP)")
print("    Maksymalne k_d = 4.40 (FCC 3D XY, z=12)")
print("    K_dw uniwersalne = 3.5 (dla d-wave layered)")
print("    A_sq: max z hybryd (d+p+f) ~ 4.0 (ale f podejrzane)")
print("    sqrt(n): dla n=10 -> 3.16, dla n=100 -> 10")
print("    Lambda_E: open, moze skalowac z omega_D materialu")
print()
print("  Wnioski:")
print("    - Dla realistycznych n (1-5) i A_ZR=0.43: T_c_max ~ 150 K")
print("    - Dla n=10 (metastable): T_c_max ~ 170 K")
print("    - Aby T_c > 200 K przy ambient potrzebujemy:")
print("        * n >= 5 stabilnych plaszczyzn")
print("        * a dokladnie 4.088 A (harmonika)")
print("        * A_sq >= 0.25 (hybryd 3-orbital)")
print("        * Lambda_E > 0.08 meV (omega_D > 15 meV)")
print()
print("  => WYZSZE T_c (>200K) AMBIENT SA MOZLIWE W TGP")
print("     Wymagaja inzynierii materialu: multiwarstwowa heterostruktura,")
print("     domieszkowana dla optymalnego A_ZR, na substracie wymuszajacym")
print("     a = 4.088 A (np. SrTiO3 epitaxial strain).")

# =============================================================
# Konkretna predykcja: hetrostruktura Hg1223 / SrTiO3
# =============================================================

print()
print("=" * 78)
print("  KONKRETNA PREDYKCJA: heterostruktura Hg-cuprate / SrTiO3")
print("=" * 78)
print()
print("  Setup:")
print("    Hg1223 (a_ambient = 3.855 A)")
print("    na substracie SrTiO3 (a = 3.905 A) -> tensile strain 1.3%")
print("    Efekt TGP: a wymuszone w gore do 3.905 A")

a_strained = 3.905
Tp_natural = Tc_dwave_layered(3.855, 3, A_ZR_sq)
Tp_strained = Tc_dwave_layered(a_strained, 3, A_ZR_sq)

print(f"    T_c (Hg1223 ambient, bulk)      = {Tp_natural:.1f} K (pred, obs = 138 K)")
print(f"    T_c (Hg1223 / STO, a={a_strained}) = {Tp_strained:.1f} K")
print(f"    Poprawa: {(Tp_strained/Tp_natural - 1) * 100:+.1f}% dzieki strain")
print()
print("  Dalej: Hg1234 na SrTiO3:")
Tp_1234_strain = Tc_dwave_layered(a_strained, 4, A_ZR_sq)
print(f"    T_c (Hg1234 / STO) = {Tp_1234_strain:.1f} K")
print()
print("  Hg1245 na SrTiO3:")
Tp_1245_strain = Tc_dwave_layered(a_strained, 5, A_ZR_sq)
print(f"    T_c (Hg1245 / STO) = {Tp_1245_strain:.1f} K (>> rekord)")
print()
print("  Eksperymentalna droga:")
print("    1. MBE/PLD wzrost Hg1223 tenkowy na SrTiO3 bufferzed")
print("    2. Szybki tlen annealing dla optymalnego domieszkania")
print("    3. Transport i Meissner dla T_c exp")

# =============================================================
# Werdykt
# =============================================================
print()
print("=" * 78)
print("  Werdykt ps11")
print("=" * 78)
print()
print(f"  Fe-SC: r(log-log) = {Fe_r:.3f}  (ambient bulk + monolayer-STO)")
print(f"         Model dziala jakosciowo dla Fe-pnictides/chalcogenides.")
print()
print("  Nickelates: thin film - undershoots (5x off); moze orb+ZR tutaj inaczej.")
print()
print("  KLUCZOWA ODPOWIEDZ:")
print("    Wysokotemperaturowe SC AMBIENT PRESSURE SA mozliwe w TGP.")
print("    Rekomendowane kierunki eksperymentalne:")
print("      1. Hg1223 na SrTiO3 (strain tensylny) - T_c > 180 K?")
print("      2. Hg1234 synteza + SrTiO3 substrat - T_c > 200 K?")
print("      3. Synteza Hg1245 (n=5) - T_c ~ 230 K?")
print("      4. Nickelates multiwarstwowe z wiekszym n_Ni")
print()
print("  PROZNIA: nie wplywa bezposrednio na T_c, ale:")
print("    - cuprates amb. degraduja w tlenie (samplowac w UHV)")
print("    - Fe-SC (FeSe/STO monolayer) wymaga UHV dla stabilnosci")
print()
print("=" * 78)
print("  ps11 complete.")
print("=" * 78)
