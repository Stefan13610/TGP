#!/usr/bin/env python3
# =============================================================================
#  r00_dataset.py
# -----------------------------------------------------------------------------
#  Etap 1 projektu rho(T): dataset N=15 pilot dla Bloch-Gruneisen + TGP tests.
#
#  Cel: zebrac rho(T) z literatury (CRC Handbook, Matula 1979, Bass 1982,
#  Landolt-Bornstein) dla typowych metali:
#     - Noble: Cu, Ag, Au
#     - sp: Al, Pb, Sn
#     - bcc TM: Nb, V
#     - FM: Fe, Ni
#     - Stoner: Pt, Pd
#     - hcp sp: Cd, Zn, Mg
#
#  Kolumny:
#     name, Theta_D[K], rho_0[uOhm*cm], rho(77K), rho(295K), rho(500K),
#     rho(1000K), a_lattice[A], N_F[states/eV/atom], Z
#
#  rho_0 zazwyczaj 10^-3..10^-1 uOhm*cm dla probek >99.999% (RRR > 1000).
#
#  Zrodla:
#    [M79] Matula, R. A., J. Phys. Chem. Ref. Data 8, 1147 (1979)
#          Electrical resistivity of Cu, Ag, Au.
#    [CRC] CRC Handbook of Chemistry and Physics, 97th ed. (2016)
#    [B82] Bass, J., Landolt-Bornstein III/15a (1982)
#    [KIT] Kittel, "Introduction to Solid State Physics" 8th ed.
#    [MART] Martin, "Electronic Structure" (2004) dla N(E_F)
# =============================================================================
import numpy as np

# Dataset ULTRA-PURE czystych metali (RRR > 1000 gdzie dostepne).
# rho w uOhm*cm. Theta_D w K. Lattice a w Angstrom. N_F w states/eV/atom.
#
# UWAGI:
#  - Fe i Ni maja anomalia rho(T) przy T_Curie (Fe: 1043K, Ni: 631K)
#  - Pb ma niski Theta_D (105K) - BG przechodzi szybko w liniowy rezim
#  - W pure jest SC przy T_c=0.012K, tu pomijamy SC nizej 4K

# schema: (name, Theta_D, rho_0, r77, r295, r500, r1000, a_latt, N_F, Z)
RHO_DATA = [
    # ------ NOBLE METALS (sp + d full) ------
    # Cu: Matula [M79] tabla glowna, rho_0 < 0.001 przy 4K dla probki 99.9999%
    # Krystal: fcc, a=3.615 A, CRC. N_F ~ 0.14 z DFT wariantu.
    ("Cu",  343, 0.0020,  0.213,  1.678,  3.21,  7.20,  3.615, 0.14, 29),
    # Ag: Matula, rho_0 ~ 0.0005 dla high-RRR. rho(295)=1.587.
    # Krystal: fcc a=4.086. N_F ~ 0.10 (Martin).
    ("Ag",  225, 0.0025,  0.284,  1.587,  3.00,  6.68,  4.086, 0.10, 47),
    # Au: Matula, rho(295)=2.214. Krystal: fcc a=4.078. N_F ~ 0.07.
    ("Au",  165, 0.0075,  0.515,  2.214,  3.85,  8.19,  4.078, 0.07, 79),

    # ------ SP METALS ------
    # Al: CRC. rho(295)=2.650. fcc a=4.046. Theta_D=428. Sc at T_c=1.18K.
    ("Al",  428, 0.0040,  0.279,  2.650,  6.13, 13.10,  4.046, 0.42, 13),
    # Pb: CRC. rho(295)=20.8. fcc a=4.95. Theta_D=105 - najnizszy w pilocie.
    # Sc at 7.2K, dla BG analiza T > T_c.
    ("Pb",  105, 0.0080,  4.700, 20.800, 38.00, 86.00,  4.950, 0.36, 82),
    # Sn (beta, white, tetr): Sn(200K)=10.6. CRC. a=5.831 (c=3.182)
    # Sc at 3.72K. N_F ~ 0.30.
    ("Sn",  200, 0.0040,  2.090, 11.000, 20.40, 45.00,  5.831, 0.30, 50),

    # ------ BCC TRANSITION METALS (SC) ------
    # Nb: Matula, Bass. rho(295)=14.6. bcc a=3.301. Theta_D=275. Sc T_c=9.26K.
    # N_F=0.97 states/eV/atom (Savrasov/Materials Project).
    ("Nb",  275, 0.0200,  3.050, 14.600, 28.00, 62.00,  3.301, 0.97, 41),
    # V: Matula. rho(295)=19.7. bcc a=3.027. Theta_D=380. Sc T_c=5.38K.
    # N_F=1.32.
    ("V",   380, 0.0100,  2.820, 19.700, 40.00, 92.00,  3.027, 1.32, 23),

    # ------ FM (dodatkowy rho_mag) ------
    # Fe: CRC. rho(295)=9.71. bcc a=2.866. Theta_D=470. T_Curie=1043K.
    # N_F=1.50. UWAGA: rho(1000K) blisko Curie anomalia.
    ("Fe",  470, 0.0100,  0.614,  9.710, 24.00, 67.00,  2.866, 1.50, 26),
    # Ni: CRC. rho(295)=6.99. fcc a=3.524. Theta_D=450. T_Curie=631K.
    # N_F=2.90. Anomalia przy 630K.
    ("Ni",  450, 0.0030,  0.550,  6.990, 19.00, 50.00,  3.524, 2.90, 28),

    # ------ STONER PARAMAGNETS ------
    # Pt: CRC. rho(295)=10.6. fcc a=3.924. Theta_D=240. NF=2.00 (Stoner).
    ("Pt",  240, 0.0160,  2.045, 10.600, 19.80, 44.00,  3.924, 2.00, 78),
    # Pd: CRC. rho(295)=10.8. fcc a=3.891. Theta_D=274. NF=2.30 (Stoner).
    ("Pd",  274, 0.0130,  2.136, 10.800, 20.30, 48.00,  3.891, 2.30, 46),

    # ------ HCP SP METALS (anizotropowe - srednia polykrystaliczna) ------
    # Cd: CRC. rho(295)=6.83. hcp a=2.979 c=5.617. Theta_D=209.
    # Sc T_c=0.52K (SC nie w naszym pilocie).
    ("Cd",  209, 0.0090,  1.020,  6.830, 13.00, 31.00,  2.979, 0.10, 48),
    # Zn: CRC. rho(295)=5.90. hcp a=2.665 c=4.947. Theta_D=327.
    # Sc T_c=0.85K.
    ("Zn",  327, 0.0050,  0.642,  5.900, 12.00, 30.00,  2.665, 0.13, 30),
    # Mg: CRC. rho(295)=4.39. hcp a=3.209 c=5.211. Theta_D=400.
    # Nie SC przy amb (jest SC przy 750 GPa wg teorii).
    ("Mg",  400, 0.0080,  0.424,  4.390, 10.10, 24.00,  3.209, 0.11, 12),
]

# Columns mapping
COL_NAME = 0
COL_THETA = 1
COL_RHO_0 = 2
COL_RHO_77 = 3
COL_RHO_295 = 4
COL_RHO_500 = 5
COL_RHO_1000 = 6
COL_A_LATT = 7
COL_NF = 8
COL_Z = 9

# T pointsdla fitu BG
T_POINTS = [77.0, 295.0, 500.0, 1000.0]
RHO_COLS = [COL_RHO_77, COL_RHO_295, COL_RHO_500, COL_RHO_1000]


def load_dataset():
    """Return list of dicts for easy iteration."""
    out = []
    for row in RHO_DATA:
        d = {
            "name":    row[COL_NAME],
            "Theta_D": row[COL_THETA],
            "rho_0":   row[COL_RHO_0],
            "rho":     {T: row[c] for T, c in zip(T_POINTS, RHO_COLS)},
            "a_latt":  row[COL_A_LATT],
            "N_F":     row[COL_NF],
            "Z":       row[COL_Z],
        }
        out.append(d)
    return out


def print_table():
    print("=" * 92)
    print("  r00_dataset.py  -  pilot rho(T) dataset, N = {}".format(
        len(RHO_DATA)))
    print("=" * 92)
    print("  {:<5s} {:>7s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>7s} {:>6s} {:>3s}".format(
        "name", "Th_D", "rho_0", "r(77)", "r(295)", "r(500)", "r(1000)",
        "a[A]", "N_F", "Z"))
    for row in RHO_DATA:
        print("  {:<5s} {:>7.0f} {:>8.4f} {:>8.3f} {:>8.3f} {:>8.2f} {:>8.2f} {:>7.3f} {:>6.2f} {:>3d}".format(
            row[COL_NAME], row[COL_THETA], row[COL_RHO_0],
            row[COL_RHO_77], row[COL_RHO_295], row[COL_RHO_500],
            row[COL_RHO_1000], row[COL_A_LATT], row[COL_NF], row[COL_Z]))


def sanity_checks():
    """Podstawowe sanity checks."""
    print("\n" + "=" * 92)
    print("  SANITY CHECKS")
    print("=" * 92)

    # (1) rho rosnie monotonicznie z T
    print("\n  (1) Monotoniczny wzrost rho(T):")
    any_issue = False
    for d in load_dataset():
        vs = [d["rho_0"], d["rho"][77], d["rho"][295],
              d["rho"][500], d["rho"][1000]]
        mono = all(vs[i] <= vs[i+1] for i in range(len(vs)-1))
        if not mono:
            print("    {:<5s} FAIL mono:  {}".format(d["name"], vs))
            any_issue = True
    if not any_issue:
        print("    OK wszystkie probki maja monotoniczne rho(T).")

    # (2) RRR = rho(295)/rho_0 (jakosc probki)
    print("\n  (2) RRR = rho(295K)/rho_0:")
    print("      {:<5s} {:>7s} {:>8s} {:>8s}".format("name", "rho_0", "rho_295", "RRR"))
    for d in load_dataset():
        rrr = d["rho"][295] / d["rho_0"]
        print("      {:<5s} {:>7.4f} {:>8.3f} {:>8.0f}".format(
            d["name"], d["rho_0"], d["rho"][295], rrr))

    # (3) High-T rezim: rho(T) ~ linear? Sprawdz nachylenie
    print("\n  (3) Slope rho(500)-rho(295) / (500-295)  [uOhm*cm/K]:")
    print("      {:<5s} {:>8s} {:>8s} {:>10s}".format(
        "name", "drho_lo", "drho_hi", "ratio"))
    for d in load_dataset():
        slope_lo = (d["rho"][500] - d["rho"][295]) / (500 - 295)
        slope_hi = (d["rho"][1000] - d["rho"][500]) / (1000 - 500)
        ratio = slope_hi / slope_lo if slope_lo > 0 else 0
        # ratio > 1: super-linear (konkawe), <1 linear (BG plateau)
        print("      {:<5s} {:>8.4f} {:>8.4f} {:>10.3f}".format(
            d["name"], slope_lo, slope_hi, ratio))

    # (4) Wysoki RRR indyukuje niski rho_0
    print("\n  (4) Klasy materiallow (wg N_F i magnetyzmu):")
    noble = [d["name"] for d in load_dataset() if d["Z"] in (29, 47, 79)]
    fm = [d["name"] for d in load_dataset() if d["name"] in ("Fe", "Ni")]
    stoner = [d["name"] for d in load_dataset() if d["name"] in ("Pt", "Pd")]
    sc_elem = [d["name"] for d in load_dataset() if d["name"] in ("Nb", "V", "Pb", "Al", "Sn")]
    hcp = [d["name"] for d in load_dataset() if d["name"] in ("Cd", "Zn", "Mg")]
    print("      Noble (Cu/Ag/Au):          ", noble)
    print("      SC-element (Nb/V/Al/Pb/Sn):", sc_elem)
    print("      FM (Fe/Ni):                ", fm)
    print("      Stoner (Pt/Pd):            ", stoner)
    print("      hcp sp (Cd/Zn/Mg):         ", hcp)


def main():
    print_table()
    sanity_checks()

    print("\n" + "=" * 92)
    print("  Dataset gotowy. Nastepny krok: r01_bg_baseline.py")
    print("=" * 92)
    print("""
    Nastepne skrypty:
      r01_bg_baseline.py  -  klasyczny Bloch-Gruneisen fit (rho_0, A_BG, Th_D).
                              Czy jedno A_BG zamyka wszystkie 15 materialow?
      r02_tgp_formula.py  -  probny wzor TGP(Theta_D): Theta_D = f(k_d, a, C_0)?
      r03_validation.py   -  LOO + residuum vs klasy materialow.

    Dalsze rozszerzenia:
      - Dodac Mo, Ta, W (refractory 5d)
      - Dodac Be, Sc, Ti (low-Z low-N_F)
      - Dodac Bi, Sb (semimetal)
      - Dodac Gd, Tb (RE FM)
    """)


if __name__ == "__main__":
    main()
