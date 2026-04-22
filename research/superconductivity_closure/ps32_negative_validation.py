#!/usr/bin/env python3
# =============================================================================
#  ps32_negative_validation.py
# -----------------------------------------------------------------------------
#  Negatywna walidacja: co Eq. 5 mowi o materialach ktore NIE sa nadprzewodnikami?
#  Jesli model jest uczciwy, powinien przewidywac T_c << 0.01 K dla Cu, Ag, Au,
#  diamond, NaCl itp. Jesli spuriously daje T_c > 1 K, to jest problem.
#
#  Dataset: ~20 materialow ktore NIE sa SC (lub sa poza klasycznym SC domain):
#    - Szlachetne metale: Cu, Ag, Au, Pt, Pd (< 0.01 K lub nie-SC)
#    - Alkaliczne: Li (0.4 mK, marginalne), Na, K (nie-SC ambient)
#    - Ferromagnet: Fe, Co, Ni (rozmagnesowanie + B_mag -> 0)
#    - Ziemie alk.: Ca, Sr, Ba (nie-SC ambient)
#    - Polmetale: C-graphite, diamond, Si, Ge (ambient, nie-SC)
#    - Jonowe: NaCl, LiF (insulators, T_c = 0)
#    - Grupa 12: Zn, Cd (znane T_c ale bardzo male: Zn=0.85K, Cd=0.52K)
#
#  Kluczowy test: dla kazdego materialu obliczamy T_c^TGP i sprawdzamy czy:
#    (a) T_c^TGP << 1 K dla nie-SC (PASS)
#    (b) T_c^TGP spuriously wysokie dla nie-SC (FAIL)
#
#  Plus: ograniczenia mozemy zaszyc jawnie:
#    - ferromagnet (Fe, Co, Ni): lam_sf -> bardzo wysokie -> B_mag -> 0
#    - insulator (NaCl, LiF, diamond): omega fononowe nisko, ale A*A niskie bo
#      brak Fermi surface -> model zaklada Fermi surf, wiec nie jest pelny
#
#  Uwaga: model NIE ma jawnej zaleznosci od N(E_F). To slabosc. Moze byc
#  subproblem P7.10 (albo osobne badanie do paperu).
# =============================================================================
import numpy as np

K_B = 8.617333e-5
A_BOHR = 0.52917721067
A_STAR = 7.725 * A_BOHR
C_0 = 48.8222
SIGMA_A = 2.5856
ALPHA_P6B = 1.04
LAMBDA_0_P6B = 0.0962
OMEGA_0 = 15.0
BETA_P6D = 2.527
A_MAP = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a):
    n = max(1, round(a / A_STAR))
    d = a - n * A_STAR
    return np.exp(-d ** 2 / SIGMA_A ** 2)

def B_mag(lam):
    return 1.0 / (1.0 + BETA_P6D * lam)

def Tc_phonon(a, orb, z, omega, lam):
    A = A_MAP[orb] if isinstance(orb, str) else float(orb) * A_MAP["d"]
    Lambda_eff = LAMBDA_0_P6B * (omega / OMEGA_0) ** ALPHA_P6B
    return (k_d(z) * C_0 * (A ** 2) * M_gauss(a)
            * (Lambda_eff * 1e-3) / K_B * B_mag(lam))


# ==============================================================================
#  Dataset: nie-nadprzewodniki (lub graniczne ambient)
#  (name, T_obs_upperbound[K], a, orb, z, omega, lam_sf, comment)
#
#  T_obs_upperbound:
#    0.0 = nigdy nie zaobserwowano SC (gorna granica z eksperymentu)
#    kilka milikelvinow = poza praktycznym zasiegiem
# ==============================================================================
NON_SC = [
    # nazwa       T_ub [K]  a[A]   orb  z  omega  lam   komentarz
    ("Cu",         0.00,   3.615, "d",  12, 28.0, 0.0,   "nigdy SC"),
    ("Ag",         0.00035,4.086, "d",  12, 21.5, 0.0,   "< 0.35 mK"),
    ("Au",         0.00005,4.078, "d",  12, 17.0, 0.0,   "< 50 uK"),
    ("Pt",         0.00,   3.924, "d",  12, 23.9, 0.0,   "paramagn., nie SC"),
    ("Pd",         0.00,   3.890, "d",  12, 25.2, 0.0,   "blisko-ferromagn."),
    ("Li",         0.0004, 3.510, "s",   8, 30.0, 0.0,   "marginalne, 0.4 mK"),
    ("Na",         0.00,   4.290, "s",   8, 13.0, 0.0,   "nigdy SC"),
    ("K",          0.00,   5.330, "s",   8,  8.6, 0.0,   "nigdy SC"),
    ("Fe",         0.00,   2.870, "d",   8, 35.0, 3.0,   "FERROMAGNET (lam=3)"),
    ("Co",         0.00,   2.510, "d",  12, 40.0, 3.0,   "FERROMAGNET"),
    ("Ni",         0.00,   3.520, "d",  12, 38.0, 2.0,   "FERROMAGNET"),
    ("Ca_amb",     0.00,   5.588, "s",  12, 19.0, 0.0,   "ambient; SC dopiero >100GPa"),
    ("Sr_amb",     0.00,   6.080, "s",  12, 12.0, 0.0,   "ambient; SC 8 GPa"),
    ("Ba_amb",     0.00,   5.030, "s",  12,  7.0, 0.0,   "ambient; SC 5 GPa"),
    ("Mg",         0.00,   3.210, "s",  12, 40.0, 0.0,   "nigdy SC ambient"),
    ("Zn",         0.85,   2.660, "sp", 12, 26.0, 0.0,   "MARGINAL SC 0.85 K"),
    ("Cd",         0.52,   2.980, "sp", 12, 16.0, 0.0,   "MARGINAL SC 0.52 K"),
    ("Diamond",    0.00,   3.570, "sp",  4, 165.0, 0.0,  "insulator"),
    ("Si_amb",     0.00,   5.431, "sp",  4,  55.0, 0.0,  "insulator ambient"),
    ("Ge_amb",     0.00,   5.658, "sp",  4,  35.0, 0.0,  "insulator ambient"),
    ("Graphite",   0.00,   2.460, "sp",  3, 170.0, 0.0,  "semimetal"),
    ("NaCl",       0.00,   5.640, "sp",  6, 26.0, 0.0,   "insulator (jonowy)"),
    ("LiF",        0.00,   4.030, "sp",  6, 80.0, 0.0,   "insulator (jonowy)"),
]


def main():
    print("=" * 78)
    print("  ps32_negative_validation.py")
    print("  Czy Eq. 5 prawidlowo przewiduje T_c ~ 0 dla nie-nadprzewodnikow?")
    print("=" * 78)

    print(f"\n  {'name':<12s} {'T_obs':>9s} {'T_TGP':>9s} {'ratio':>8s}"
          f" {'verdict':>10s}  komentarz")
    print("  " + "-" * 76)

    pass_count = 0
    fail_count = 0
    borderline = 0
    results = []

    for mat in NON_SC:
        name, Tub, a, orb, z, om, lam, comment = mat
        Tp = Tc_phonon(a, orb, z, om, lam)

        # Kryterium PASS:
        #   (1) Tub = 0 (nigdy SC): T_TGP < 0.1 K = PASS, 0.1-1 = BORDER, >1 = FAIL
        #   (2) Tub > 0 (marginalne SC): T_TGP < 5x T_ub = PASS
        if Tub < 1e-4:
            # nigdy-SC
            if Tp < 0.1:
                verdict = "PASS"
                pass_count += 1
            elif Tp < 1.0:
                verdict = "BORDER"
                borderline += 1
            else:
                verdict = "FAIL"
                fail_count += 1
            ratio_str = "  -"
        else:
            # marginalny SC
            ratio = Tp / Tub
            ratio_str = f"{ratio:.1f}x"
            if ratio < 5.0:
                verdict = "PASS"
                pass_count += 1
            elif ratio < 20.0:
                verdict = "BORDER"
                borderline += 1
            else:
                verdict = "FAIL"
                fail_count += 1

        Tub_str = f"{Tub:.4f}" if Tub > 0 else "~0"
        Tp_str = f"{Tp:.3f}" if Tp < 100 else f"{Tp:.1e}"
        print(f"  {name:<12s} {Tub_str:>9s} {Tp_str:>9s} {ratio_str:>8s} "
              f"{verdict:>10s}  {comment}")
        results.append((name, Tub, Tp, verdict))

    # -------- podsumowanie
    print("\n" + "=" * 78)
    print(f"  Wynik: {pass_count} PASS, {borderline} BORDER, {fail_count} FAIL"
          f"  (z {len(NON_SC)} materialow)")
    print("=" * 78)

    fails = [r for r in results if r[3] == "FAIL"]
    if fails:
        print("\n  FAIL (model spuriously przewiduje SC gdzie go nie ma):")
        for name, Tub, Tp, _ in fails:
            print(f"    {name:<12s}  T_TGP = {Tp:.3f} K  vs  T_obs <= "
                  f"{Tub*1e3:.3f} mK")
        print("\n  Analiza FAIL:")
        print("    Model bez jawnej zaleznosci od N(E_F) moze dawac spuriously")
        print("    wysokie T_c dla materialow bez Fermi surface (insulatory),")
        print("    oraz bez B_ferromagn dla ferromagnetykow bez podanego lam_sf.")
    else:
        print("\n  Zadne FAIL — negatywna walidacja PRZESZLA w pelni.")

    borders = [r for r in results if r[3] == "BORDER"]
    if borders:
        print("\n  BORDER (nie-krytyczne, T_TGP w zakresie 0.1-1 K):")
        for name, Tub, Tp, _ in borders:
            print(f"    {name:<12s}  T_TGP = {Tp:.3f} K")

    # -------- Kluczowe lekcje
    print("\n" + "=" * 78)
    print("  Analiza sub-kategoryczna:")
    print("=" * 78)

    categories = {
        "Szlachetne (Cu,Ag,Au,Pt,Pd)": ["Cu","Ag","Au","Pt","Pd"],
        "Alkaliczne (Li,Na,K)":         ["Li","Na","K"],
        "Ferromagnetyki (Fe,Co,Ni)":    ["Fe","Co","Ni"],
        "Ziemie alk. (Ca,Sr,Ba,Mg)":    ["Ca_amb","Sr_amb","Ba_amb","Mg"],
        "Marginal SC (Zn,Cd)":           ["Zn","Cd"],
        "Insulatory (diamond,Si,Ge,NaCl,LiF)": ["Diamond","Si_amb","Ge_amb","NaCl","LiF"],
        "Graphite":                      ["Graphite"],
    }

    for cat, names in categories.items():
        subres = [r for r in results if r[0] in names]
        if not subres: continue
        passes = sum(1 for r in subres if r[3] == "PASS")
        fails = sum(1 for r in subres if r[3] == "FAIL")
        borders = sum(1 for r in subres if r[3] == "BORDER")
        status = "OK" if fails == 0 else "ISSUE"
        print(f"  {cat:<38s} {passes}/{len(subres)} PASS, "
              f"{borders} BORDER, {fails} FAIL  [{status}]")

    # -------- sugestie dla paperu
    print("\n" + "=" * 78)
    print("  Sugestie dla paperu SC v2 (sekcja 'Domain of validity')")
    print("=" * 78)
    print("""
  Model (Eq. 5) z B_mag(lam_sf) automatycznie tlumi SC dla ferromagnetykow
  gdy lam_sf >> 1. Dla insulatorow model NIE MA bezposredniego mechanizmu
  (brak zaleznosci od N(E_F)) — ale ma M_gauss(a) penalty ktory tlumi
  materialy z nietypowym a.

  Potencjalny P7.10 (na przyszlosc): dodac explicit N(E_F) factor w Eq. 5
    A^2 -> A^2 * g(N_F / N_F_ref)  z g(0) = 0 (insulatory -> T_c = 0)

  Dla paperu v2: rozdzial 'Domain of validity' powinien explicite stwierdzic:
    - Dziala: phonon-mediated SC, cuprates, hydrydy (z d_F > 0)
    - Nie dziala: heavy-fermion (non-phonon), insulators (N_F=0), ferromagnetyki
      (wymagana lam_sf kalibracja), czyste 5d elementy (P7.6 partial fix)
""")


if __name__ == "__main__":
    main()
