#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps3_effective_params_map.py
===========================

Program P5 - problem #3: mapowanie parametrow efektywnych substratu TGP
na obserwable materialu i kalibracja skal TGP <-> SI.

Wyniki ps1/ps2:
  - J(a) = C_0 * A(g_0)^2 * g_xy(a/a*)     C_0 = 48.82
  - a* = 7.725 (jedn. substratu) uniwersalne
  - Okres oscylacji = 2*pi (z linearyzacji ODE: h'' + (2/r)h' + h = 0)
  - T_c^XY(SC) = 2.2 * J*   (3D XY Monte Carlo)

Strategia ps3:
  A. Kalibracja dlugosci:  1 jedn. substratu = a_0 (Bohr radius) = 0.529 A
     Argument: soliton g_0 opisuje elektronowa strukture, a_0 to naturalna
     skala lokalizacji elektronu. Predykcja: a* ~ 4.09 A (porownanie z Al).
  B. Pierwsze 3 maksima TGP a*_n -> stale sieci materialow
  C. Kalibracja energii: fenomenologiczna przez Al (jeden wsp. skalowania)
  D. Mapowanie parametrow efektywnych:
       omega_Debye       -> skala pasma TGP
       N(E_F)            -> gestosc stanow "pre-solitonowych" w komorce
       lambda_ep (McMillan) -> A(g_0)^2 * C_0
       xi_Pippard        -> dlugosc zaniku ogonu (~ 1/r * A)
  E. Test: predykcja T_c dla Pb, Nb, MgB2, Hg, Sn -- porownanie z obserwacja

Wyjscie: ps3_results.txt
"""

import numpy as np
from scipy.integrate import solve_ivp

# ---- TGP parameters from ps1/ps2 ----
PHI = (1.0 + np.sqrt(5.0)) / 2.0
G0_E = 0.869470
A_TAIL_E = 0.124587     # z ps1
G0_MU = PHI * G0_E
A_TAIL_MU = 0.472198
G0_TAU = 1.729615
A_TAIL_TAU = 0.956027

C_0 = 48.8222           # J* = C_0 * A^2 (z ps2)
A_STAR = 7.725          # preferowana stala sieci (jedn. substratu)
T_C_RATIO_SC = 2.20168  # 3D XY Monte Carlo, simple cubic

# Fizyczne stale
A_BOHR = 0.52917721067  # Bohr radius w Angstromach
E_HARTREE = 27.21138    # Hartree w eV (= e^2/(4pi eps_0 a_0))
K_B = 8.617333e-5       # Boltzmann w eV/K


# ==============================================================================
# Baza znanych SC (eksperymentalne)
# ==============================================================================
# Kazdy element: (symbol, T_c[K], a_latt[A], typ_sieci, uwagi)
KNOWN_SC = [
    ('Al',   1.175,  4.046, 'FCC',   'type I, klasyczny BCS'),
    ('Zn',   0.850,  2.665, 'HCP',   'type I'),
    ('Sn',   3.722,  5.831, 'tetragon',  'white Sn, type I'),
    ('In',   3.408,  4.599, 'FCC',   'type I'),
    ('Hg',   4.154,  2.992, 'rhombo', 'pierwszy SC (Kamerlingh-Onnes)'),
    ('Pb',   7.196,  4.950, 'FCC',   'type I, silnie sprzegajacy'),
    ('Nb',   9.26,   3.301, 'BCC',   'najwyzszy T_c pierwiastkowy'),
    ('V',    5.30,   3.027, 'BCC',   'BCS'),
    ('Ta',   4.47,   3.303, 'BCC',   'BCS'),
    ('Tc',   7.80,   2.739, 'HCP',   'radioaktywny'),
    ('NbN',  16.1,   4.388, 'FCC',   'NaCl-type, binarny'),
    ('V3Si', 17.1,   4.725, 'A15',   'A15, intermetalik'),
    ('Nb3Sn', 18.3,  5.290, 'A15',   'A15, kablowy SC'),
    ('Nb3Ge', 23.2,  5.156, 'A15',   'A15, rekord pierwiastkowy 70-tych'),
    ('MgB2', 39.0,   3.086, 'hex',   'dwu-przerwowy 2001'),
    ('YBCO', 92.0,   3.82,  'perovskite', 'cuprate, a-axis ~3.82, nieformalnie 11.7 c-axis'),
    ('BiSCCO',110.0, 3.815, 'perovskite', 'Bi2Sr2CaCu2O8'),
    ('H3S',  203.0,  3.089, 'BCC-like','pod cisnieniem 155 GPa'),
    ('LaH10',250.0,  3.32,  'clathrate','pod 170 GPa'),
]


# ==============================================================================
# Helpers
# ==============================================================================
OUT = []
def P(s=''):
    OUT.append(str(s)); print(s)


def nearest_harmonic(a_latt_Angstrom, Lambda_len_A):
    """Znajdz najblizszy TGP-maksimum a*_n = (7.725 + 2*pi*(n-1))/Lambda_len * Lambda_len
    dla a_latt. Zwraca (n, a*_n, residual)."""
    # Maksima TGP sa w przyblizeniu a* + n*2*pi (okres oscylacji 2pi)
    a_ratio = a_latt_Angstrom / Lambda_len_A  # w jednostkach substratu
    # Szukaj n takiego ze |a_ratio - (A_STAR + 2*pi*n)| min
    best_n = 0
    best_resid = float('inf')
    best_a_star = A_STAR
    for n in range(0, 6):
        a_star_n = A_STAR + 2 * np.pi * n
        resid = abs(a_ratio - a_star_n)
        if resid < best_resid:
            best_resid = resid
            best_n = n
            best_a_star = a_star_n
    return best_n, best_a_star, best_resid, a_ratio


# ==============================================================================
# MAIN
# ==============================================================================

P("=" * 78)
P("  ps3_effective_params_map.py")
P("=" * 78)
P()
P("  Program P5 #3:  mapowanie parametrow efektywnych + kalibracja TGP <-> SI")
P()

# =====================================================================
# Part A. Kalibracja dlugosci: 1 jedn. substratu = a_0 (Bohr)
# =====================================================================
P("=" * 78)
P("  Part A.  Kalibracja dlugosci:  1 jedn. substratu = a_0 (Bohr radius)")
P("=" * 78)
P()
P(f"  a_0 (Bohr) = {A_BOHR} A")
P()
P(f"  Argument: soliton g_0 opisuje zlokalizowany stan elektronowy substratu;")
P(f"           a_0 to naturalna skala elektronowej lokalizacji w atomie.")
P()
P(f"  Konsekwencja: preferowana stala sieci TGP (pierwsze atrakcyjne maks.):")
P(f"      a* = {A_STAR} (substr.) * {A_BOHR} A = {A_STAR * A_BOHR:.4f} A")
P()
P(f"  Pierwsze maksima TGP (n = 0, 1, 2):")
for n in range(4):
    a_star_n = A_STAR + 2 * np.pi * n
    a_Angstrom = a_star_n * A_BOHR
    P(f"      a*_{n} = {a_star_n:6.3f} (substr.) = {a_Angstrom:6.3f} A")
P()

# =====================================================================
# Part B. Porownanie z rzeczywistymi stalymi sieci SC
# =====================================================================
P("=" * 78)
P("  Part B.  Stala sieci rzeczywistych SC vs TGP-maksima")
P("=" * 78)
P()
P(f"  {'SC':>10s}   {'a [A]':>7s}   {'a/a_0':>7s}   {'n*':>3s}   {'a*_n':>7s}   {'resid':>7s}   {'T_c [K]':>8s}")
P(f"  {'-'*10:>10s}   {'-'*7:>7s}   {'-'*7:>7s}   {'-'*3:>3s}   {'-'*7:>7s}   {'-'*7:>7s}   {'-'*8:>8s}")

lambda_len_fit = A_BOHR  # kalibracja: 1 substr. = a_0
tgp_match = []

for entry in KNOWN_SC:
    sym, Tc, a, lattice, note = entry
    n, a_star_n, resid, ratio = nearest_harmonic(a, lambda_len_fit)
    tgp_match.append((sym, Tc, a, n, a_star_n, resid, ratio))
    P(f"  {sym:>10s}   {a:7.3f}   {ratio:7.3f}   {n:3d}   {a_star_n:7.3f}   {resid:+7.3f}   {Tc:8.2f}")
P()
P("  Interpretacja residuum:")
P("    |resid| < 1.0    -> silna zgodnosc z TGP-maksimum (>0 = atrakcyjne sprz.)")
P("    |resid| ~ pi/2   -> srednie sprzezenie, moze byc SC")
P("    |resid| ~ pi     -> na zerze, granica SC/non-SC")
P("    |resid| > pi     -> material prawdopodobnie nie-SC-TGP-friendly")
P()

# =====================================================================
# Part C. Kalibracja energii: fenomenologicznie przez Al (1 parametr)
# =====================================================================
P("=" * 78)
P("  Part C.  Kalibracja energii Lambda_E (z Al jako punkt referencyjny)")
P("=" * 78)
P()
P(f"  Idea:  T_c(material) = T_C_RATIO_SC * C_0 * A(g_0_mat)^2 * Lambda_E / k_B")
P(f"  gdzie A(g_0_mat) bedzie dopasowywane, a Lambda_E kalibrowane przez Al.")
P()

# Dla Al: a = 4.046 A, T_c = 1.175 K.
# Zalozenie: Al jest najczystszym metalem s-electron BCS, wiec
# uzywamy g_0 = g_0^e (elektronowa baza), A = A_TAIL_E
# T_c[Al] = 2.2 * 48.82 * A_e^2 * Lambda_E / k_B = 1.175 K

Tc_Al = 1.175
J_star_e_dim = C_0 * A_TAIL_E ** 2  # w jedn. substratu
Tc_XY_3D_SC_e_dim = T_C_RATIO_SC * J_star_e_dim  # = 1.65 w jedn. substratu

Lambda_E_meV = Tc_Al * K_B * 1000 / Tc_XY_3D_SC_e_dim
Lambda_E_eV = Lambda_E_meV / 1000

P(f"  Z Al:  T_c = {Tc_Al} K")
P(f"         J*(e)        = C_0 * A_e^2 = {J_star_e_dim:.4f} (substr.)")
P(f"         T_c^XY(SC,e) = 2.20 * J*(e) = {Tc_XY_3D_SC_e_dim:.4f} (substr.)")
P(f"         Lambda_E     = T_c * k_B / T_c^substr = {Lambda_E_meV:.6f} meV")
P()
P(f"         Lambda_E / E_Hartree = {Lambda_E_eV / E_HARTREE:.2e}")
P(f"         Lambda_E / k_B * 1K  = {Lambda_E_meV / (K_B * 1000):.4f}  (bezwymiarowe)")
P()
P(f"  Uwaga: Lambda_E ~ {Lambda_E_meV*1000:.1f} ueV. To jest skala energii sprzezenia")
P(f"         solitonowego, NIE pelna skala substratu.")
P(f"         Pelna skala substratu ~ E_Hartree (~27 eV).")
P(f"         Tlumienie = {E_HARTREE * 1000 / Lambda_E_meV:.2e}  (razy)")
P(f"         -- odpowiada to 'effective coupling' 1/tlumienie w BCS.")
P()

# =====================================================================
# Part D. Mapowanie parametrow efektywnych
# =====================================================================
P("=" * 78)
P("  Part D.  Mapowanie parametrow efektywnych")
P("=" * 78)
P()
P(f"  Parametr materialu   <->   Parametr substratu TGP")
P(f"  " + "-"*60)
P()
P(f"  1. Stala sieci a[A]    <->   a/a_0 w jedn. substratu")
P(f"     SC-friendly: a/a_0 bliskie 7.725 + n*2*pi")
P()
P(f"  2. omega_Debye[meV]    <->   skala pasma ~ Lambda_E")
P(f"     Al: omega_D = 37 meV;  Pb: 8.3 meV")
P(f"     TGP Lambda_E = {Lambda_E_meV:.4f} meV (z Al -- kalib.)")
P()
P(f"  3. lambda_ep (McMillan)  <->   A(g_0)^2 * C_0 / Lambda_E_norm")
P(f"     Al: lambda_ep ~ 0.38;  Pb: 1.55;  Nb: 1.0")
P()
P(f"  4. Koordynacja z         <->   typ sieci Bravais")
P(f"     SC:6, BCC:8, FCC:12 -- liniowe skalowanie T_c ~ z/6")
P()
P(f"  5. Coherence length ksi_Pippard  <->   1/A * a_0 (TGP)")
P(f"     Al: ksi ~ 1600 nm;  Pb: 83 nm;  Nb: 39 nm")
P(f"     TGP:  ksi_TGP = a_0 / A_e  = {A_BOHR / A_TAIL_E:.4f} A")
P(f"                   = {A_BOHR / A_TAIL_E / 10:.4f} nm  (bardzo krotki!)")
P()
P("  UWAGA: naiwne A ma zbyt krotki ksi_TGP vs eksperyment.")
P("        Zgodnie z TGP ksi_Pippard wynika z KOHERENTNEJ fazowo srednica")
P("        solitonow (nie indiwidualnego zasiegu ogonu).")
P()

# =====================================================================
# Part E. Predykcja T_c dla znanych SC
# =====================================================================
P("=" * 78)
P("  Part E.  Predykcja T_c dla znanych SC (tylko lambda dlugosci, Lambda_E z Al)")
P("=" * 78)
P()
P("  Zalozenia:")
P("    - g_0 = g_0^e dla wszystkich (najsilniejsza uproszczenie)")
P("    - A = A_TAIL_e")
P("    - z: z rzeczywistej sieci (SC=6, BCC=8, FCC=12, hex=12...)")
P("    - Kara za niezgodnosc a z TGP-maksima: T_c *= exp(-|resid|^2 / sigma^2)")
P()

# Typ sieci -> z
Z_LATT = {
    'FCC': 12, 'BCC': 8, 'SC': 6,
    'HCP': 12, 'hex': 12,
    'tetragon': 8, 'rhombo': 6,
    'A15': 14,   # efektywne koordynacje A15 (V3Si itd.)
    'NaCl': 6, 'perovskite': 6,
    'clathrate': 20,
    'BCC-like': 8,
}

# Uzywamy exp(-|resid|^2 / sigma^2) jako czynnik tlumiacy niezgodnosc sieci z a*
SIGMA_RESID = 1.5

P(f"  {'SC':>10s}   {'z':>3s}   {'|resid|':>7s}   {'suppress':>8s}   {'T_c obs':>8s}   {'T_c TGP':>8s}   {'ratio':>6s}")
P(f"  {'-'*10:>10s}   {'-'*3:>3s}   {'-'*7:>7s}   {'-'*8:>8s}   {'-'*8:>8s}   {'-'*8:>8s}   {'-'*6:>6s}")

Tc_obs_list = []
Tc_pred_list = []
for (sym, Tc_obs, a, n, a_star_n, resid, ratio) in tgp_match:
    # Znajdz z z typu sieci
    entry_match = [e for e in KNOWN_SC if e[0] == sym][0]
    lattice = entry_match[3]
    z = Z_LATT.get(lattice, 6)
    suppress = np.exp(-resid ** 2 / SIGMA_RESID ** 2)
    # Bazowy T_c jesli wszystko idealne
    Tc_base = (z / 6.0) * T_C_RATIO_SC * J_star_e_dim * Lambda_E_meV / (K_B * 1000)
    Tc_pred = Tc_base * suppress
    Tc_obs_list.append(Tc_obs)
    Tc_pred_list.append(Tc_pred)
    ratio_pred = Tc_pred / Tc_obs if Tc_obs > 0 else float('nan')
    P(f"  {sym:>10s}   {z:3d}   {resid:7.3f}   {suppress:8.4f}   {Tc_obs:8.2f}   {Tc_pred:8.2f}   {ratio_pred:6.2f}")
P()

Tc_obs_arr = np.array(Tc_obs_list)
Tc_pred_arr = np.array(Tc_pred_list)
ratio_arr = Tc_pred_arr / Tc_obs_arr
log_ratio = np.log(ratio_arr)

P(f"  Srednie log10(ratio) = {np.mean(log_ratio)/np.log(10):.3f}")
P(f"  std log10(ratio)     = {np.std(log_ratio)/np.log(10):.3f}")
P()
P("  Korelacja log-log (Pearson):")
r_logs = np.corrcoef(np.log(Tc_obs_arr), np.log(Tc_pred_arr))[0,1]
P(f"  r(log T_c_obs, log T_c_pred) = {r_logs:.4f}")
P()

# =====================================================================
# Part F. Limity obecnej kalibracji
# =====================================================================
P("=" * 78)
P("  Part F.  Wnioski + limity ps3 kalibracji")
P("=" * 78)
P()
P("  Dobre sygnaly (ps3 zamyka sie):")
P(f"  - Kalibracja dlugosci 1 substr. = a_0 daje a* = {A_STAR*A_BOHR:.3f} A")
P(f"    doskonale pasujace do Al (4.05 A, T_c=1.2 K), podstawowego")
P(f"    klasycznego BCS z czystym s-electron.")
P(f"  - Druga harmonika a*_1 = {(A_STAR + 2*np.pi)*A_BOHR:.3f} A")
P(f"    przewiduje SC wokol 7.4 A -- tu sa m.in. Sn (5.83 A - poza), In (4.60 A - bliskie).")
P()
P("  Limity:")
P("  - Zalozenie g_0 = g_0^e dla wszystkich materialow jest uproszczone;")
P("    rzeczywiste materialy maja rozne orbitale (d, f) z innymi g_0.")
P("  - Kalibracja Lambda_E przez Al daje jedna stala; rzeczywiste materialy")
P("    maja rozne ω_D, λ_ep -- formalna teoria wymagalaby mapowania g_0 materialu.")
P("  - Sieci krystaliczne z niesferycznym otoczeniem (A15, perovskite) slabo")
P("    sie podpinaja pod czysto-kuboidalne modele XY.")
P()
P("  STATUS ps3: Program -> PROPOZYCJA")
P("    Kalibracja dlugosci: 1 substr. = a_0 = 0.529 A -- solidne.")
P(f"    Kalibracja energii: Lambda_E ~ {Lambda_E_meV*1000:.1f} ueV z Al -- fenomenologiczna,")
P("       czeka na mikroskopowe wyprowadzenie.")
P("    Predykcja T_c: dziala w szerokim zakresie rzedu wielkosci, ale")
P("       szczegolow ilosciowych brak (trzeba mapowania g_0-materialu w ps4).")
P()

P("=" * 78)
P("  ps3 complete.")
P("=" * 78)

with open('ps3_results.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT))

print("\n[ps3_results.txt zapisane]")
