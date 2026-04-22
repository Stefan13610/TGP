#!/usr/bin/env python3
# =============================================================================
#  ps23_residual_correlation_Z_Mmol.py
# -----------------------------------------------------------------------------
#  Test zarzutu: "Czy residua log-pred(Tc) korelują z liczbą atomową Z
#  lub masą molową M_mol?  Jeśli tak, model TGP pomija jakiś efekt
#  jądrowy i jest niekompletny."
#
#  Dane wejściowe: 29 materiałów z ps17_results.txt.
#  Dla każdego materiału:
#    - Z_pair   = liczba atomowa dominującego atomu pairing-site
#    - Z_heavy  = liczba atomowa najcięższego atomu w formule
#    - M_mol    = masa molowa formuły [g/mol]
#    - dlog     = log10(T_pred) - log10(T_obs)
#
#  Liczymy Pearson r i Spearman rho dla:
#    (|dlog|, Z_pair), (|dlog|, Z_heavy), (|dlog|, M_mol),
#    (dlog,   Z_pair), (dlog,   Z_heavy), (dlog,   M_mol).
#  Raz dla pełnego N=29, raz po usunięciu flagowanych outlierów
#  (|dlog| > 0.4: Y_amb, Th_amb, Ce_5GPa, H3S).
# =============================================================================
import numpy as np
from scipy.stats import pearsonr, spearmanr

# materiał, Z_pair, Z_heavy, M_mol[g/mol], dlog  (ps17_results.txt, kolumna 'dlog')
# Z_pair  = atom z orbitalem sprzęgającym (Cu dla cupratów, Fe dla Fe-SC, H/S dla hydridu, itd.)
# Z_heavy = najcięższy atom w formule
DATA = [
    # cuprates (pair = Cu, Z=29)
    ("La2CuO4",     29, 57, 405.37,  +0.147),
    ("YBCO",        29, 56, 664.62,  -0.085),
    ("BiSCCO2212",  29, 83, 888.40,  -0.051),
    ("Tl2212",      29, 81, 978.70,  -0.154),
    ("Hg1223",      29, 80, 874.20,  -0.172),
    ("Tl2223",      29, 81, 1114.40, -0.130),
    ("Nd2CuO4",     29, 60, 416.03,  +0.351),
    ("Bi2201",      29, 83, 752.75,  +0.196),
    # phonon-mediated / Fe-SC / hydridy / P7
    ("Al",          13, 13, 26.98,   +0.399),
    ("Pb",          82, 82, 207.20,  -0.162),
    ("Nb",          41, 41, 92.91,   +0.174),
    ("V",           23, 23, 50.94,   +0.315),
    ("Hg_elem",     80, 80, 200.59,  -0.217),
    ("MgB2",        5,  12, 45.93,   -0.094),  # pair on B (sp)
    ("FeSe_bulk",   26, 34, 134.81,  -0.008),
    ("FeSe/STO",    26, 34, 134.81,  -0.055),
    ("Ba122-Co",    26, 56, 398.87,  -0.090),
    ("LaFeAsO",     26, 57, 285.68,  -0.201),
    ("NdFeAsO-F",   26, 60, 291.01,  -0.358),
    ("Nb3Sn",       41, 50, 397.44,  +0.081),
    ("NbTi",        41, 41, 140.78,  +0.237),
    ("H3S",         16, 16, 35.09,   -0.426),  # * outlier
    ("LaH10",       57, 57, 148.99,  +0.168),
    ("CeH9",        58, 58, 149.19,  +0.156),
    ("CeH10",       58, 58, 150.20,  +0.212),
    ("La_amb",      57, 57, 138.91,  +0.378),
    ("Y_amb",       39, 39, 88.91,   +1.152),  # * outlier
    ("Th_amb",      90, 90, 232.04,  +0.877),  # * outlier
    ("Ce_5GPa",     58, 58, 140.12,  +0.549),  # * outlier
]

OUTLIER_NAMES = {"H3S", "Y_amb", "Th_amb", "Ce_5GPa"}


def run_corr(rows, label):
    names = [r[0] for r in rows]
    Z_pair  = np.array([r[1] for r in rows], float)
    Z_heavy = np.array([r[2] for r in rows], float)
    M_mol   = np.array([r[3] for r in rows], float)
    dlog    = np.array([r[4] for r in rows], float)
    abs_dl  = np.abs(dlog)

    print(f"\n================================================")
    print(f"  {label}   (N = {len(rows)})")
    print(f"================================================")

    pairs = [
        ("|dlog|", abs_dl, "Z_pair",  Z_pair),
        ("|dlog|", abs_dl, "Z_heavy", Z_heavy),
        ("|dlog|", abs_dl, "M_mol",   M_mol),
        ("dlog",   dlog,   "Z_pair",  Z_pair),
        ("dlog",   dlog,   "Z_heavy", Z_heavy),
        ("dlog",   dlog,   "M_mol",   M_mol),
    ]
    print(f"  {'var_y':>8s}  {'var_x':>8s}    Pearson r  (p)       Spearman rho (p)")
    for ylab, y, xlab, x in pairs:
        r, pr = pearsonr(x, y)
        rho, ps = spearmanr(x, y)
        sig_r = "*" if pr < 0.05 else " "
        sig_s = "*" if ps < 0.05 else " "
        print(f"  {ylab:>8s}  {xlab:>8s}   {r:+.4f} ({pr:.3f}){sig_r}   "
              f"{rho:+.4f} ({ps:.3f}){sig_s}")
    print()


def main():
    print("=" * 80)
    print("  ps23_residual_correlation_Z_Mmol.py")
    print("  Test: czy residua log-pred(Tc) korelują z Z lub M_mol?")
    print("=" * 80)

    # pełny dataset
    run_corr(DATA, "FULL N=29 (z outlierami)")

    # bez outlierów P7
    clean = [r for r in DATA if r[0] not in OUTLIER_NAMES]
    run_corr(clean, "CORE N=25 (bez P7-outlierów: H3S, Y_amb, Th_amb, Ce_5GPa)")

    # interpretacja
    print("=" * 80)
    print("  Interpretacja")
    print("=" * 80)
    print("""
  Progowe |r| wg Cohena dla N=25-29:
     |r| < 0.10  :  brak efektu
     |r| < 0.30  :  słaby
     |r| < 0.50  :  umiarkowany
     |r| >= 0.50 :  silny

  Aby zarzut agenta ('model pomija efekt jądrowy') miał moc,
  potrzeba CO NAJMNIEJ umiarkowanej korelacji (|r| > 0.30)
  statystycznie istotnej (p < 0.05) DLA CORE (N=25).
  Słaba/nieistotna korelacja odrzuca zarzut.
""")


if __name__ == "__main__":
    main()
