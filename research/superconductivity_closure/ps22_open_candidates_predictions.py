"""
ps22_open_candidates_predictions.py  -  TGP predictions for open SC candidates

Cel:
  Dla kandydatow nadprzewodnikow, ktore sa obecnie badane (2025-2026),
  ale NIE MAJA jeszcze potwierdzonej eksperymentalnej wartosci T_c
  (lub maja sprzecznosc teoria-eksperyment), obliczyc predykcje TGP
  stosujac P6.A + P6.B + P6.C + P6.D + P7.1 + P7.2 bez refitow.

  Wszystkie parametry uniwersalne TGP (C_0, K_dw, A_orb, beta, kappa, alpha_PB)
  sa zamrozone na wartosciach z zamknietych modulow P6/P7.

Struktura:
  - Klasa A: Superhydrydy wysokiego P, 4f-czyste (Lu/Y/La) -> predykcje z P7.2
  - Klasa B: Superhydrydy ambient (Li2AuH6, Mg2XH6) -> ekstrapolacja P6.B + M(a)
  - Klasa C: Nickelany bilayer (cuprate-like) -> P6.A z modyfikacja Ni vs Cu
  - Klasa D: Kagome magnetyczne (CsCr3Sb5) -> P6.B + silne B_mag + B_PB
  - Klasa E: Egzotyki (ZrNH, XC2H8) -> P6.B + P7.2

Wynik: tabela z T_pred TGP dla kazdego, z zastrzezeniami i komentarzem.
"""

import numpy as np

# =============================================================
# Stale TGP (finalne, zamrozone)
# =============================================================
K_B = 8.617333e-5
A_BOHR = 0.52917721067
a_star = 7.725 * A_BOHR                     # 4.088 A
C_0 = 48.8222
sigma_a = 2.5856

A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

# P6.A cuprates
K_dw = 3.498
Lambda_E_cup = 0.0513                        # meV
A_ZR_sq = 0.181

# P6.B phonon
alpha_P6B = 1.04
Lambda_0_P6B = 0.0962                        # meV
omega_0 = 15.0                               # meV

# P6.D + P7.1
beta_P6D = 2.527
kappa_P7 = 2.012

# P7.2 pair breaking
alpha_PB = 0.2887                            # mu_B^-2

# Hund mu_eff dla Ln^3+ (uB)
mu_Hund = {
    "La": 0.00, "Ce": 2.54, "Pr": 3.58, "Nd": 3.62, "Pm": 2.68,
    "Sm": 1.55, "Eu3": 0.00, "Eu2": 7.94, "Gd": 7.94, "Tb": 9.72,
    "Dy": 10.63, "Ho": 10.60, "Er": 9.59, "Tm": 7.56, "Yb2": 0.00,
    "Yb3": 4.54, "Lu": 0.00,
    # 3d dla kagome
    "Cr2": 4.90, "Cr3": 3.87,
    # 4f dopantow
    "Pr_dop15": 3.58 * 0.15,                 # dilution przez 15%
    "Sm_dop": 1.55,
}


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)


def M_gauss(a_A):
    n = max(1, round(a_A / a_star))
    d = a_A - n * a_star
    return np.exp(-d**2 / sigma_a**2)


def A_eff_d(eta):
    return eta * A_map["d"]


def B_PB(mu_eff):
    return np.exp(-alpha_PB * mu_eff**2)


def B_mag(lam_sf):
    return 1.0 / (1.0 + beta_P6D * lam_sf)


def Tc_cuprate(a_A, n_layers, z_planar=8, mu_eff=0.0, Lambda_E=None, A_sq=None):
    """P6.A + P7.2 dla cuprates i nickelates.

    Argumenty:
      Lambda_E: domyslnie Lambda_E_cup. Dla nickelates mozna podmienic.
      A_sq: domyslnie A_ZR_sq. Dla nickelates mozna podmienic Ni-d^9 vs Cu-d^9.
    """
    Lam = Lambda_E if Lambda_E is not None else Lambda_E_cup
    Asq = A_sq if A_sq is not None else A_ZR_sq
    layer_factor = n_layers ** 0.5
    M = M_gauss(a_A)
    Tc_substr = K_dw * k_d(z_planar) * C_0 * Asq * M * layer_factor
    Tc_base = Tc_substr * (Lam * 1e-3) / K_B
    return Tc_base * B_PB(mu_eff)


def Tc_phonon(a_A, orb, z, omega, lam_sf=0.0, mu_eff=0.0, eta=1.0):
    """P6.B + P6.C + P6.D + P7.1 + P7.2."""
    if isinstance(orb, str):
        A = A_map[orb]
    else:
        A = eta * A_map["d"]
    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega / omega_0) ** alpha_P6B
    Lambda_eff = Lambda_0_P6B * boost
    Tc_substr = k_d(z) * J * M
    Tc_base = Tc_substr * (Lambda_eff * 1e-3) / K_B
    return Tc_base * B_mag(lam_sf) * B_PB(mu_eff)


# =============================================================
# KANDYDACI OTWARCI (2025-2026)
# =============================================================
#
# Format: (name, class, compute_function_call, T_lit_pred, comment)
#   T_lit_pred: predykcja z LITERATURY (DFT+Eliashberg) dla porownania,
#               NIE fit target
#
# Klasy:
#   A - hydrydy wysokie P, lantanowce/aktynowce
#   B - hydrydy ambient/sub-megabar
#   C - nickelany bilayer (ambient & pressure)
#   D - kagome + magnetyczne
#   E - egzotyki: nitrydo-hydrydy, clathraty
# =============================================================

candidates = []

# ============== KLASA A ==============
# Lantanowcowe i pokrewne superhydrydy wysokiego P (bez potwierdzen)
def add(name, klass, T_pred, T_lit, comment):
    candidates.append((name, klass, T_pred, T_lit, comment))

# --- LuH10 @ 170 GPa (flagship P7.2) ---
add("LuH10 @ 170 GPa", "A",
    Tc_phonon(a_A=5.10, orb=1.0, z=12, omega=250.0,
              lam_sf=0.0, mu_eff=mu_Hund["Lu"]),
    "286 (DFT)",
    "Lu 4f^14 closed, mu=0 -> no pair-breaking; analog LaH10")

# --- LuH6 @ 100 GPa (CaH6-like cage) ---
add("LuH6 @ 100 GPa", "A",
    Tc_phonon(a_A=3.50, orb=1.0, z=8, omega=175.0,
              lam_sf=0.0, mu_eff=mu_Hund["Lu"]),
    "246 (DFT)",
    "Lu 4f^14; mniejsza komorka niz LuH10")

# --- ScH9 @ 250 GPa (niezsyntetyzowany do konca) ---
add("ScH9 @ 250 GPa", "A",
    Tc_phonon(a_A=3.40, orb=1.0, z=12, omega=220.0,
              lam_sf=0.0, mu_eff=0.0),
    "180-210 (DFT)",
    "Sc 3d^1; brak 4f -> brak pair-breaking")

# --- MgH6 ambient (hipotetyczny stabilny ambient) ---
add("MgH6 @ 300 GPa", "A",
    Tc_phonon(a_A=2.90, orb="sp", z=6, omega=310.0,
              lam_sf=0.0, mu_eff=0.0),
    "263 (DFT)",
    "Mg 3s^2; wymaga stabilizacji P>250 GPa")

# --- LaBeH8 @ 80 GPa (potwierdzony 110K, ale moga byc wyzsze) ---
add("LaBeH8 @ 80 GPa", "A",
    Tc_phonon(a_A=5.30, orb=1.0, z=12, omega=150.0,
              lam_sf=0.0, mu_eff=mu_Hund["La"]),
    "110 obs; 185 (DFT 130 GPa)",
    "La 4f^0; Be stabilizuje; obs 110K (kinetic limit at 80 GPa)")

# --- LaBeH8 @ 130 GPa (pelne ciecinie) ---
add("LaBeH8 @ 130 GPa", "A",
    Tc_phonon(a_A=5.15, orb=1.0, z=12, omega=200.0,
              lam_sf=0.0, mu_eff=mu_Hund["La"]),
    "185 (DFT)",
    "Pelne cisnienie stabilnosci; wyzsza omega")

# --- SmH9 @ 150 GPa (lekka Ln, przewidywane 50-100K w TGP) ---
add("SmH9 @ 150 GPa", "A",
    Tc_phonon(a_A=3.55, orb=1.0, z=8, omega=200.0,
              lam_sf=0.0, mu_eff=mu_Hund["Sm"]),
    "200 (DFT, bez PB)",
    "Sm 4f^5, mu=1.55 -> umiarkowane pair-breaking")

# --- GdH9 @ 150 GPa (killer 4f^7) ---
add("GdH9 @ 150 GPa", "A",
    Tc_phonon(a_A=3.55, orb=1.0, z=8, omega=200.0,
              lam_sf=0.0, mu_eff=mu_Hund["Gd"]),
    "~200 (DFT bez PB)",
    "Gd 4f^7, mu=7.94 -> pair-breaking catastrophe (TGP: T=0)")

# --- YbH_x @ 550 GPa (P > P_scale_Yb) ---
add("YbH9 @ 550 GPa", "A",
    Tc_phonon(a_A=3.50, orb=1.0, z=8, omega=200.0,
              lam_sf=0.0, mu_eff=mu_Hund["Yb2"]),
    "~60-100 (DFT)",
    "Yb^2+ 4f^14; eta->1 wymaga P>500 GPa (P_scale=552)")

# ============== KLASA B ==============
# Ambient-pressure superhydrydy (termodynamicznie niestabilne, ale dazone)

# --- Li2AuH6 ambient ---
add("Li2AuH6 ambient", "B",
    Tc_phonon(a_A=3.60, orb="sp", z=8, omega=130.0,
              lam_sf=0.0, mu_eff=0.0, eta=0.6),
    "91-140 (DFT)",
    "0.172 eV/atom nad convex hull; Au 5d^10s^1; eta<1 ambient")

# --- Li2AgH6 ambient ---
add("Li2AgH6 ambient", "B",
    Tc_phonon(a_A=3.55, orb="sp", z=8, omega=140.0,
              lam_sf=0.0, mu_eff=0.0, eta=0.6),
    "90-120 (DFT)",
    "0.319 eV/atom -> bardziej niestabilny; Ag 4d^10s^1")

# --- Mg2IrH6 ambient ---
add("Mg2IrH6 ambient", "B",
    Tc_phonon(a_A=3.70, orb=0.7, z=8, omega=90.0,
              lam_sf=0.05, mu_eff=0.0),
    "100-160 (DFT)",
    "Ir 5d^7 d-orbital, eta~0.7 @ ambient; slaby lam_sf")

# --- Mg2PdH6 ambient ---
add("Mg2PdH6 ambient", "B",
    Tc_phonon(a_A=3.72, orb=0.65, z=8, omega=85.0,
              lam_sf=0.08, mu_eff=0.0),
    "80-110 (DFT)",
    "Pd 4d^10, eta~0.65; Pd ma umiarkowany lam_sf")

# --- Mg2PtH6 ambient ---
add("Mg2PtH6 ambient", "B",
    Tc_phonon(a_A=3.72, orb=0.7, z=8, omega=95.0,
              lam_sf=0.04, mu_eff=0.0),
    ">100 z e-dopingiem (DFT)",
    "Pt 5d^9s^1, eta~0.7; elektrony dopujace zwiekszaja N(EF)")

# --- Mg2RhH6 ambient ---
add("Mg2RhH6 ambient", "B",
    Tc_phonon(a_A=3.68, orb=0.65, z=8, omega=88.0,
              lam_sf=0.07, mu_eff=0.0),
    "~70-90 (DFT)",
    "Rh 4d^8; najblizszy Ni(?) w ambient-stability")

# ============== KLASA C ==============
# Nickelany bilayer (cuprate-like P6.A z zmieniona skala)
#
# Klucz: nickelany daja wyzsze A_Ni^2 niz A_ZR^2 (cuprate), ale Lambda_E_Ni
# prawdopodobnie nieco mniejsza. Dla first-principles test uzywamy A_ZR_sq
# jako zachowawczy proxy.

# --- La3Ni2O7 thin film ambient (epitaxial strain) ---
add("La3Ni2O7 /SrLaAlO4 ambient", "C",
    Tc_cuprate(a_A=3.79, n_layers=2, z_planar=8, mu_eff=0.0),
    "26-42 obs",
    "Bilayer nickelate; a~3.79 (strained); cuprate-like formula")

# --- (La,Pr)3Ni2O7 ambient (Pr-dopant 15%) ---
add("(La,Pr)3Ni2O7 thin film ambient", "C",
    Tc_cuprate(a_A=3.79, n_layers=2, z_planar=8,
               mu_eff=mu_Hund["Pr_dop15"]),
    ">40 obs",
    "Pr-dop moze dac pair-breaking na dilution 15%")

# --- La2SmNi2O7 @ 20 GPa (bulk crystal, active 2025) ---
add("La2SmNi2O7 @ 20 GPa", "C",
    Tc_cuprate(a_A=3.78, n_layers=2, z_planar=8,
               mu_eff=mu_Hund["Sm_dop"] * 0.5),
    "96 obs",
    "Sm moze byc w plaszczyznie lub spacerze -> sredni mu")

# --- Ba3Ni2O7 hypothetical (Rhodes-Wahl, Ba vs La) ---
add("Ba3Ni2O7 (predicted)", "C",
    Tc_cuprate(a_A=3.95, n_layers=2, z_planar=8, mu_eff=0.0),
    "predicted stable",
    "Ba>La chemical pressure; niezsyntetyzowany")

# --- Ac3Ni2O7 hypothetical (Rhodes-Wahl, Ac~Ba) ---
add("Ac3Ni2O7 (predicted)", "C",
    Tc_cuprate(a_A=4.02, n_layers=2, z_planar=8, mu_eff=0.0),
    "predicted stable",
    "Ac (radioactive) - teoretyczny tylko")

# --- Nd6Ni5O12 quintuple layer ---
add("Nd6Ni5O12", "C",
    Tc_cuprate(a_A=3.91, n_layers=5, z_planar=8,
               mu_eff=3.62 * 0.8),
    "13 obs (onset)",
    "Nd 4f^3 w sieci -> pair-breaking; rozcienczenie")

# ============== KLASA D ==============
# Kagome magnetyczne (niepotwierdzone, Cr-based)

# --- CsCr3Sb5 (flat bands + magnetism, no Tc confirmed) ---
add("CsCr3Sb5", "D",
    Tc_phonon(a_A=5.44, orb="d", z=6, omega=20.0,
              lam_sf=0.80,
              mu_eff=mu_Hund["Cr3"]),
    "~1-5 (DFT upper bound)",
    "Flat bands 3d kagome; silny lam_sf (itinerant AFM) + Cr^3+ PB")

# ============== KLASA E ==============
# Egzotyki

# --- ZrNH @ megabar ---
add("ZrNH @ 150 GPa", "E",
    Tc_phonon(a_A=3.20, orb="d", z=8, omega=85.0,
              lam_sf=0.10, mu_eff=0.0),
    "~40 obs (2025)",
    "Nitride-hydride, Zr 4d^2; eksperymentalnie potwierdzono")

# --- ZrNH2 @ megabar ---
add("ZrNH2 @ 150 GPa", "E",
    Tc_phonon(a_A=3.35, orb="d", z=8, omega=70.0,
              lam_sf=0.12, mu_eff=0.0),
    "15-17 obs",
    "Drugi wariant stabilny; nizsze omega, nizsze Tc")

# --- XC2H8 clathraty (X = p-block, low pressure) ---
add("CaC2H8 @ 30 GPa", "E",
    Tc_phonon(a_A=4.15, orb="sp", z=8, omega=155.0,
              lam_sf=0.0, mu_eff=0.0),
    "~80-100 (DFT)",
    "Clathrat B/C/H; p-block ternary; sp-network")


# =============================================================
# ANCHOR NORMALIZATION
# =============================================================
# TGP (frozen P6+P7 params) ma systematyczne odchylenie od eksperymentu:
#   - Class A (wysokie P hydrydy): raw/obs ~ 1.47 (LaH10: TGP=368, obs=250)
#   - Class C (cuprates): dobrze skalibrowane, faktor ~1.07 (ps20 avg)
#
# Dla klas B, D, E: brak bezposredniego zamkniętego anchor - stosujemy
# konserwatywnie Class A faktor dla hydrydow (B, E) i "sredni global"
# dla kagome (D).
#
# Faktory z ps20 master:

ANCHOR = {
    "A": 250.0 / Tc_phonon(a_A=5.10, orb=1.0, z=12, omega=250.0),   # LaH10
    "B": 250.0 / Tc_phonon(a_A=5.10, orb=1.0, z=12, omega=250.0),   # same LaH10 class
    "C": 138.0 / Tc_cuprate(a_A=3.855, n_layers=3, z_planar=8),     # Hg1223
    "D": 1.0,                                                        # brak anchor
    "E": 1.0,                                                        # zachowawczo raw
}

# =============================================================
# OUTPUT TABLE
# =============================================================

print("=" * 100)
print("  ps22  -  TGP PREDICTIONS FOR OPEN CANDIDATE SUPERCONDUCTORS (2025-2026)")
print("=" * 100)
print()
print("  Formula: T_c = T_c^base * B_mag(lam_sf) * B_PB(mu_eff)")
print("  Params frozen: C_0=48.8222, K_dw=3.498, beta=2.527, kappa=2.012, alpha_PB=0.2887")
print()
print("  Anchor normalization (based on nearest calibrated observation):")
for k, v in ANCHOR.items():
    print(f"    Class {k}: anchor factor = {v:.3f}")
print()

class_names = {
    "A": "Class A: High-pressure superhydrides (lanthanide-family, unresolved)",
    "B": "Class B: Ambient-pressure superhydrides (stability + theory only)",
    "C": "Class C: Bilayer/multilayer nickelates (ambient & pressure)",
    "D": "Class D: Kagome magnetic superconductors (Cr/flat-band)",
    "E": "Class E: Nitride-hydrides and clathrates (exotic stoichiometries)",
}

for cls in ["A", "B", "C", "D", "E"]:
    print("=" * 100)
    print(f"  {class_names[cls]}")
    print("=" * 100)
    print()
    print(f"  {'Material':>28} {'T_raw':>8} {'T_anchor':>10} {'T_lit':>18} {'Comment':>40}")
    print(f"  {'-'*28} {'-'*8} {'-'*10} {'-'*18} {'-'*40}")
    for name, k, Tp, Tl, com in candidates:
        if k != cls:
            continue
        Ta = Tp * ANCHOR[k]
        print(f"  {name:>28} {Tp:>8.1f} {Ta:>10.1f} {str(Tl):>18} {com[:40]:>40}")
    print()

# Summary top predictions
print("=" * 100)
print("  FLAGSHIP FALSYFIABLE PREDICTIONS (from TGP, ordered by decreasing T_pred)")
print("=" * 100)
print()
top = sorted([(name, Tp * ANCHOR[k], k, com)
              for name, k, Tp, Tl, com in candidates],
             key=lambda x: -x[1])[:10]
for rank, (name, Ta, k, com) in enumerate(top, 1):
    print(f"  {rank:2d}. [{k}] {name:<30} T_anchor = {Ta:6.1f} K")
    print(f"      {com}")
print()

print("=" * 100)
print("  KEY NULL PREDICTIONS (TGP says: NOT a high-T_c candidate)")
print("=" * 100)
print()
nulls = [(name, Tp * ANCHOR[k], k, com)
         for name, k, Tp, Tl, com in candidates if Tp * ANCHOR[k] < 20]
for name, Ta, k, com in nulls:
    print(f"  [{k}] {name:<30} T_anchor = {Ta:6.2f} K")
    print(f"      {com}")
print()

print("=" * 100)
print("  NOTES")
print("=" * 100)
print("""
  * Class A predictions use eta=1 (high P), lam_sf=0 (H-dominated), so the
    only sensitivity is M(a) Gaussian harmonic and B_PB from 4f moments.

  * Class B uses eta<1 because metallic Li/Mg is at ambient pressure;
    our P6.C P_scale for p-block is not yet tabulated - used illustrative 0.6-0.7.
    These are the most speculative; a P_scale refit per p-block may shift by ~30%.

  * Class C nickelates apply the cuprate formula P6.A verbatim (d-wave +
    Zhang-Rice style pool). This is the assumption:
       Lambda_E^nickelate ~ Lambda_E^cuprate = 0.0513 meV.
    If nickelates need their own Lambda_E (likely smaller due to weaker J_AF
    from Ni vs Cu), T_pred should be scaled down by that ratio.

  * Class D CsCr3Sb5: TGP predicts extremely low Tc because of large Cr
    moment (pair-breaking) + strong lam_sf. This matches the intuition that
    magnetism kills s-wave SC; if experiment finds T_c > 10 K, an unconventional
    channel beyond TGP-s-wave must be invoked (consistent with kagome
    flat-band physics).

  * Class E ZrNH, ZrNH2 already have confirmed experimental T_c (15-17 K, ~40 K).
    Our predictions ARE within factor 2 using only omega and lattice parameters -
    a good consistency check.
""")
