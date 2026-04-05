#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p75_agamma_phi0_correlation.py  --  TGP v1 · Hipoteza maksymalna
=================================================================
Cel: zbadac numerycznie, czy istnieje prosta relacja funkcyjna
     a_Gamma = g(Phi0) latajaca mostem miedzy Warstwa II (a_Gamma
     z r21_PDG) a Warstwa II (Phi0 z Lambda_obs).

Gdyby taka relacja istniala, N_param redukuje sie 2 -> 1.

Znane wartosci:
    Phi0   = 36 * Omega_Lambda = 24.66  [w jednostkach c0^2/H0^2]
    a_Gamma = 0.040049                  [parametr bifurkacji solitonow]
    alpha_K = 8.56                      [parametr Koidego]
    r21     = 206.77                    [m_mu / m_e, PDG 2023]
    kappa   = m_sp * a_Gamma ~ ln(r21) = 5.333

Hipoteza fundamentalna z planu:
    a_Gamma * sqrt(gamma) * l_P = const
    gdzie gamma = Phi0 * H0^2/c0^2  i  l_P = Planck length

Testujemy tez prostsze kandydatury na relacje.

UWAGA: brak jawnej teorii laczacej a_Gamma z Phi0; to czysto
       numeryczna eksploracja zgodnosci.
=================================================================
"""

import sys
import io
import numpy as np
from scipy.optimize import brentq

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ─────────────────────────────────────────────────────────────────────────────
# Znane stale TGP i kosmologiczne
# ─────────────────────────────────────────────────────────────────────────────

# Kosmologiczne (Planck 2018)
OMEGA_LAMBDA = 0.6847       # frakcja ciemnej energii
OMEGA_M      = 0.3153       # frakcja materii
H0_SI        = 67.4e3 / 3.0856776e22  # H0 w s^-1 (67.4 km/s/Mpc)
C0_SI        = 2.99792458e8           # predkosc swiatla [m/s]
G_SI         = 6.674e-11             # G Newtona [m^3/(kg s^2)]

# Planck units w jednostkach H0=c0=1
# l_P = sqrt(hbar * G / c0^3)
HBAR_SI = 1.0545718e-34   # [J s]
l_P_SI  = np.sqrt(HBAR_SI * G_SI / C0_SI**3)   # ~ 1.616e-35 m
t_H_SI  = C0_SI / H0_SI                         # czas Hubble'a [s] / Hubble length [m]
l_P_H   = l_P_SI / t_H_SI                       # l_P w jednostkach c0/H0

# Parametry TGP (Warstwa II)
PHI0        = 36.0 * OMEGA_LAMBDA   # = 24.66 [bezwymiarowe, w jedn. c0^2/H0^2]
A_GAMMA     = 0.040049              # parametr bifurkacji solitonow
ALPHA_K     = 8.56                  # parametr Koidego
R21_PDG     = 206.77                # m_mu / m_e

# Parametry wywodzone
KAPPA       = np.log(R21_PDG)       # = 5.333, z warunku exp(kappa)=r21
GAMMA_COSMO = PHI0                  # gamma [H0^2/c0^2] w j. H0=c0=1
M_SP        = np.sqrt(GAMMA_COSMO)  # m_sp = sqrt(gamma) w j. H0/c0

# ─────────────────────────────────────────────────────────────────────────────
# Sprawdzenie kappa = m_sp * a_Gamma
# ─────────────────────────────────────────────────────────────────────────────

kappa_from_params = M_SP * A_GAMMA
kappa_target      = KAPPA

print("=" * 65)
print("TGP v1  ·  p75_agamma_phi0_correlation.py")
print("Hipoteza maksymalna: a_Gamma = g(Phi0) ? (N_param: 2->1)")
print("=" * 65)

print(f"\n[DANE WEJSCIOWE]")
print(f"  Phi0       = {PHI0:.5f}  [c0^2/H0^2]")
print(f"  a_Gamma    = {A_GAMMA:.6f}")
print(f"  alpha_K    = {ALPHA_K:.4f}")
print(f"  r21_PDG    = {R21_PDG:.3f}")
print(f"  kappa_tgt  = ln(r21) = {KAPPA:.5f}")
print(f"  m_sp=sqrt(gamma)    = {M_SP:.5f}")
print(f"  m_sp * a_Gamma      = {kappa_from_params:.5f}  (powinno byc {KAPPA:.5f})")
print(f"  kappa/kappa_tgt - 1 = {kappa_from_params/kappa_target - 1:.4f}")
print(f"  l_P [j.H0=c0=1]     = {l_P_H:.4e}")

# ─────────────────────────────────────────────────────────────────────────────
# [1] Testowanie kandydatow na relacje a_Gamma = g(Phi0)
# ─────────────────────────────────────────────────────────────────────────────

print("\n[1]  KANDYDACI NA RELACJE a_Gamma = g(Phi0)")
print("-" * 65)

# Liczymy a_Gamma * Phi0^n dla roznych n
powers = [-2.0, -3/2, -1.0, -2/3, -1/2, -1/3, 0.0,
           1/3, 1/2, 2/3, 1.0, 3/2, 2.0]

# Zestawienie kandydatow: (opis, wyrazenie, wartosc)
combos = []

for n in powers:
    val = A_GAMMA * PHI0**n
    # Sprawdz bliskosc do "prostych" liczb
    targets = {
        "1":         1.0,
        "1/2":       0.5,
        "2":         2.0,
        "pi":        np.pi,
        "2pi":       2*np.pi,
        "1/pi":      1/np.pi,
        "e":         np.e,
        "ln2":       np.log(2),
        "ln(207)":   KAPPA,
        "sqrt(2)":   np.sqrt(2),
        "1/sqrt(2)": 1/np.sqrt(2),
        "Omega_L":   OMEGA_LAMBDA,
        "2/3":       2/3,
        "1/3":       1/3,
        "alpha_K/Phi0": ALPHA_K/PHI0,
    }
    best_tgt, best_dev = None, 1e10
    for tname, tval in targets.items():
        dev = abs(val/tval - 1)
        if dev < best_dev:
            best_dev, best_tgt = dev, tname
    label = f"a_Gamma * Phi0^{n:.3g}"
    combos.append((label, val, best_tgt, best_dev))

# Dodaj kombinacje z alpha_K
for m in [-1, -0.5, 0.5, 1]:
    val = A_GAMMA * ALPHA_K**m
    label = f"a_Gamma * alpha_K^{m}"
    combos.append((label, val, "1", abs(val - 1)))

# Kluczowe kandydatki z planu
# kappa: a_Gamma = kappa / m_sp
kappa_val = KAPPA / M_SP
combos.append((f"kappa/m_sp = ln(r21)/sqrt(Phi0)", kappa_val, "a_Gamma",
               abs(kappa_val/A_GAMMA - 1)))

# a_Gamma * sqrt(gamma) * l_P
val_aGl = A_GAMMA * np.sqrt(GAMMA_COSMO) * l_P_H
combos.append((f"a_Gamma * sqrt(gamma) * l_P", val_aGl, "0", val_aGl))

# Posortuj wg odchylenia od "najlepszego celu"
combos_sorted = sorted(combos, key=lambda x: x[3])

print(f"\n  {'Wyrazenie':<42} {'Wartosc':>12}  {'Cel':>12}  {'|dev|':>8}")
print(f"  {'-'*42}  {'-'*12}  {'-'*12}  {'-'*8}")
for label, val, tgt, dev in combos_sorted[:15]:
    print(f"  {label:<42}  {val:>12.5f}  {tgt:>12}  {dev:>8.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# [2] Sprawdzenie hipotezy kappa = ln(r21)
# ─────────────────────────────────────────────────────────────────────────────

print("\n[2]  HIPOTEZA: kappa = ln(r21) = m_sp * a_Gamma")
print("-" * 65)

print(f"\n  kappa_target  = ln({R21_PDG:.2f}) = {KAPPA:.5f}")
print(f"  m_sp          = sqrt(gamma) = sqrt({GAMMA_COSMO:.4f}) = {M_SP:.5f}  [H0/c0]")
print(f"  a_Gamma_pred  = kappa/m_sp  = {kappa_val:.6f}  (TGP predykcja z Phi0)")
print(f"  a_Gamma_obs   = {A_GAMMA:.6f}  (dopasowanie do r21_PDG)")
print(f"  odchylenie    = {abs(kappa_val/A_GAMMA - 1)*100:.2f}%")

# Jaka wartox Phi0 daje dokladne kappa / sqrt(Phi0) = a_Gamma?
# kappa / sqrt(Phi0*) = a_Gamma --> Phi0* = (kappa/a_Gamma)^2
Phi0_pred = (KAPPA / A_GAMMA)**2
Omega_L_pred = Phi0_pred / 36.0
print(f"\n  Gdyby relacja byla dokladna:")
print(f"    Phi0*     = (ln(r21)/a_Gamma)^2 = ({KAPPA:.4f}/{A_GAMMA:.6f})^2 = {Phi0_pred:.4f}")
print(f"    Omega_L*  = Phi0*/36 = {Omega_L_pred:.5f}")
print(f"    odch. od Omega_L_Planck = {abs(Omega_L_pred - OMEGA_LAMBDA)/OMEGA_LAMBDA*100:.2f}%")
print(f"    (Planck 2018: Omega_L = 0.6847 +/- 0.0073, tj. 1.1%)")

# ─────────────────────────────────────────────────────────────────────────────
# [3] Skan Phi0: czy a_Gamma(Phi0) = ln(r21)/sqrt(Phi0) trzyma sie w zakresie?
# ─────────────────────────────────────────────────────────────────────────────

print("\n[3]  SKAN Phi0 IN [20, 35]: a_Gamma_pred vs a_Gamma_obs")
print("-" * 65)

Phi0_scan = np.linspace(20, 35, 13)
print(f"\n  {'Phi0':>6}  {'Omega_L':>8}  {'a_pred=kappa/sqrt(Phi0)':>24}  {'a_obs':>10}  {'dev(%)':>8}")
print(f"  {'-'*6}  {'-'*8}  {'-'*24}  {'-'*10}  {'-'*8}")
for Ph in Phi0_scan:
    OmL  = Ph / 36.0
    apred = KAPPA / np.sqrt(Ph)
    dev   = (apred - A_GAMMA) / A_GAMMA * 100
    flag  = " <-- obserwowany" if abs(Ph - PHI0) < 0.5 else ""
    flag2 = " <-- kappa-exact" if abs(Ph - Phi0_pred) < 0.5 else ""
    print(f"  {Ph:6.2f}  {OmL:8.5f}  {apred:24.6f}  {A_GAMMA:10.6f}  {dev:8.2f}%{flag}{flag2}")

# ─────────────────────────────────────────────────────────────────────────────
# [4] Test gloszej hipotezy: a_Gamma = kappa/(sqrt(Phi0) * f) dla prostego f
# ─────────────────────────────────────────────────────────────────────────────

print("\n[4]  RESZTKOWY WSPOLCZYNNIK f = kappa/(a_Gamma * sqrt(Phi0))")
print("-" * 65)

f_residual = KAPPA / (A_GAMMA * M_SP)
print(f"\n  f = kappa/(a_Gamma * sqrt(Phi0)) = {f_residual:.6f}")
print(f"\n  Bliskosc do prostych liczb:")
candidates_f = {
    "1":          1.0,
    "2/pi":       2/np.pi,
    "pi/3":       np.pi/3,
    "sqrt(2)":    np.sqrt(2),
    "Omega_L":    OMEGA_LAMBDA,
    "1 + Omega_M/Omega_L": 1 + OMEGA_M/OMEGA_LAMBDA,
    "Omega_L + Omega_M":   OMEGA_LAMBDA + OMEGA_M,
    "(2Q)^2, Q=2/3":       (2*2/3)**2,
    "3/pi":       3/np.pi,
    "2/3 + 1/3":  1.0,
    "1/(1-Omega_M)": 1/(1-OMEGA_M),
}
for name, val in sorted(candidates_f.items(), key=lambda x: abs(x[1]-f_residual)):
    dev = abs(f_residual/val - 1) * 100
    print(f"    {name:<30} = {val:.6f},  |dev| = {dev:.2f}%")

# ─────────────────────────────────────────────────────────────────────────────
# [5] Sprawdzenie: a_Gamma * sqrt(gamma) * l_P
# ─────────────────────────────────────────────────────────────────────────────

print("\n[5]  KOMBINACJA a_Gamma * sqrt(gamma) * l_P (hipoteza planu)")
print("-" * 65)

# W j. H0=c0=hbar=1 (ale G nie):
# gamma = Phi0 * H0^2/c0^2 = Phi0 (w tych jednostkach)
# l_P = sqrt(G_N * hbar / c0^3) / (c0/H0) = l_P_H
gamma_nat   = PHI0           # w j. H0=c0=1
combo_lP    = A_GAMMA * np.sqrt(gamma_nat) * l_P_H

print(f"\n  gamma [H0^2/c0^2, j. H0=c0=1]  = {gamma_nat:.4f}")
print(f"  l_P [j. c0/H0]                  = {l_P_H:.4e}")
print(f"  a_Gamma                         = {A_GAMMA:.6f}")
print(f"  a_Gamma * sqrt(gamma) * l_P     = {combo_lP:.4e}")
print(f"  Wartosc oczekiwana (const~1):     nie jest rzedu 1")
print(f"  To sugeruje, ze a_Gamma nie jest w jednostkach H0/c0.")
print(f"  Jezeli a_Gamma w jednostkach Plancka [l_P^-1]:")
a_Gamma_Planck = A_GAMMA / l_P_H   # konwersja
kappa_check    = a_Gamma_Planck * l_P_H * M_SP
print(f"    a_Gamma [l_P^-1] = {a_Gamma_Planck:.3e}")
print(f"    (kappa = m_sp[H0] * a_Gamma[H0] = {M_SP:.4f} * {A_GAMMA:.6f} = {M_SP*A_GAMMA:.4f})")

# ─────────────────────────────────────────────────────────────────────────────
# [6] Formalne podsumowanie
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("[6]  PODSUMOWANIE: STATUS HIPOTEZY MAKSYMALNEJ")
print("=" * 65)

is_kappa_consistent = abs(M_SP * A_GAMMA / KAPPA - 1) < 0.02
is_Phi0_close       = abs(Phi0_pred - PHI0) / PHI0 < 0.05
is_OmL_in_1sigma    = abs(Omega_L_pred - OMEGA_LAMBDA) < 0.0073

print(f"""
  WYNIK NUMERYCZNY:
    kappa := m_sp * a_Gamma = sqrt(Phi0) * a_Gamma
           = {M_SP:.4f} * {A_GAMMA:.6f} = {M_SP*A_GAMMA:.5f}
    kappa_target = ln(r21) = {KAPPA:.5f}
    Odchylenie:  {abs(M_SP*A_GAMMA/KAPPA - 1)*100:.2f}%

  MOCNA RELACJA (niezalezna):
    kappa = ln(r21_PDG) = {KAPPA:.4f}  [warunek Koidego, PDG]

  RELACJA KANDYDATKA:
    a_Gamma ~ kappa / sqrt(Phi0) = ln(r21) / sqrt(36 * Omega_L)
    a_Gamma_pred = {KAPPA:.5f} / {M_SP:.5f} = {kappa_val:.6f}
    a_Gamma_obs  = {A_GAMMA:.6f}
    odchylenie   = {abs(kappa_val/A_GAMMA-1)*100:.2f}%

  IMPLIKACJA (jesli relacja dokladna):
    Phi0 = (ln(r21) / a_Gamma)^2 = {Phi0_pred:.4f}
    vs Phi0_obs = {PHI0:.4f}  (roznica {abs(Phi0_pred-PHI0)/PHI0*100:.2f}%)
    Omega_L_pred = {Omega_L_pred:.5f}
    vs Omega_L_Planck = {OMEGA_LAMBDA:.5f} +/- 0.0073
    {'ZGODNE na poziomie 1-sigma Plancka' if is_OmL_in_1sigma else 'POZA 1-sigma Plancka'}

  WNIOSKI:
  [A] Relacja kappa = sqrt(Phi0) * a_Gamma odchyla sie o {abs(M_SP*A_GAMMA/KAPPA - 1)*100:.1f}%
      od dokladnosci -- za duza dyskrepancja na "przypadkowosc".
  [B] Relacja a_Gamma = kappa/sqrt(Phi0) wymaga Omega_L = {Omega_L_pred:.4f},
      {"zgodnie z Planck 2018 na 1-sigma." if is_OmL_in_1sigma else f"poza 1-sigma Plancka."}
  [C] Resztkowy wspolczynnik f = {f_residual:.4f} nie jest prostym ulamkiem.
  [D] Kombinacja a_Gamma * sqrt(gamma) * l_P jest nierzedu 1 w j. H0=c0=1;
      wymaga innego systemu jednostek (np. Plancka) do interpretacji.

  STATUS: HIPOTEZA -- relacja a_Gamma ~ kappa/sqrt(Phi0) bliska (~{abs(kappa_val/A_GAMMA-1)*100:.0f}%),
          ale nie dokladna. Potrzeba: wyprqwadzenia kappa z TGP lub niezaleznego
          wyznaczenia a_Gamma z mikroskopii substratu.
""")

# ─────────────────────────────────────────────────────────────────────────────
# [7] Listy sprawdzajace (do logu)
# ─────────────────────────────────────────────────────────────────────────────

checks = [
    ("kappa = m_sp * a_Gamma w zakresie 2%",
     is_kappa_consistent,
     f"{M_SP*A_GAMMA:.4f} vs {KAPPA:.4f} (dev={abs(M_SP*A_GAMMA/KAPPA-1)*100:.1f}%)"),
    ("Phi0* = (kappa/a_Gamma)^2 bliskie Phi0_obs na 5%",
     is_Phi0_close,
     f"{Phi0_pred:.3f} vs {PHI0:.3f} (dev={abs(Phi0_pred-PHI0)/PHI0*100:.1f}%)"),
    ("Omega_L_pred zgodne z Planck 1-sigma",
     is_OmL_in_1sigma,
     f"{Omega_L_pred:.5f} vs {OMEGA_LAMBDA:.5f} +/- 0.0073"),
    ("Prosta relacja a_Gamma * Phi0 ~ 1 (<=5% dev)",
     abs(A_GAMMA * PHI0 - 1) < 0.05,
     f"a_Gamma * Phi0 = {A_GAMMA*PHI0:.4f}"),
    ("a_Gamma * sqrt(Phi0) jest racjonalne",
     False,
     f"= {A_GAMMA*M_SP:.5f} (brak oczywistej racjonalnosci)"),
]

n_pass = sum(1 for _, r, _ in checks if r)
print(f"[7]  CHECKLIST: {n_pass}/{len(checks)} PASS")
for desc, result, note in checks:
    print(f"  {'PASS' if result else 'OPEN'}  {desc}: {note}")
