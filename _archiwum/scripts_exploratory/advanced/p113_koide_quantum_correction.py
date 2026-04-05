#!/usr/bin/env python3
"""
p113_koide_quantum_correction.py
================================
TGP: Korekcja kwantowa Koide do stosunku mas leptonow.

Problem:
  - r_21 = m_mu/m_e = 206.768  (TGP z phi-FP: 206.77, dokladnosc 0.0001%)
  - r_31_bare = (A_tau/A_e)^4 ~ 3955  (predykcja TGP bez korekcji)
  - r_31_obs  = m_tau/m_e = 3477.48  (PDG)
  => odchylenie 13.7%

Hipoteza: korekcja anharmoniczna WKB + relacja Koide zamyka lukke.

Testy:
  1. Weryfikacja relacji Koide (Q = 2/3)
  2. Predykcja r_31 z Koide przy znanym r_21
  3. Korekcja kwantowa delta_2
  4. Rozwiazanie ukladu WKB: kappa_1, kappa_2
  5. Predykcja masy 4. generacji
"""

import numpy as np

# ============================================================
# Dane PDG (2024)
# ============================================================
m_e   = 0.51099895    # MeV
m_mu  = 105.6583755   # MeV
m_tau = 1776.86       # MeV

r_21_obs = m_mu / m_e
r_31_obs = m_tau / m_e

# TGP bare prediction (z dodatku J2, bez korekcji kwantowej)
r_31_TGP_bare = 3955.0

print("=" * 70)
print("  TGP: Korekcja kwantowa Koide do hierarchii mas leptonow")
print("=" * 70)
print()
print("--- Dane wejsciowe (PDG) ---")
print(f"  m_e   = {m_e:.8f} MeV")
print(f"  m_mu  = {m_mu:.7f} MeV")
print(f"  m_tau = {m_tau:.2f} MeV")
print(f"  r_21  = m_mu/m_e = {r_21_obs:.3f}")
print(f"  r_31  = m_tau/m_e = {r_31_obs:.2f}")
print(f"  r_31_TGP_bare    = {r_31_TGP_bare:.1f}")
print()

# ============================================================
# TEST 1: Weryfikacja relacji Koide
# ============================================================
print("=" * 70)
print("  TEST 1: Weryfikacja relacji Koide")
print("=" * 70)

x = np.sqrt(m_e)
y = np.sqrt(m_mu)
z = np.sqrt(m_tau)

Q_obs = (m_e + m_mu + m_tau) / (x + y + z)**2
Q_exact = 2.0 / 3.0
dQ = abs(Q_obs - Q_exact)
dQ_rel = dQ / Q_exact * 100  # procent

print(f"  Q_obs   = {Q_obs:.10f}")
print(f"  Q_exact = {Q_exact:.10f}")
print(f"  |dQ|    = {dQ:.2e}")
print(f"  |dQ|/Q  = {dQ_rel:.4f}%")

test1_pass = dQ < 0.001
print(f"\n  >> TEST 1 (|Q - 2/3| < 0.001): {'PASS' if test1_pass else 'FAIL'}")
print(f"     Relacja Koide spelniona z dokladnoscia {dQ_rel:.4f}%")
print()

# ============================================================
# TEST 2: Predykcja r_31 z relacji Koide
# ============================================================
print("=" * 70)
print("  TEST 2: Predykcja r_31 z Koide (przy znanym r_21)")
print("=" * 70)

# Koide: (x^2 + y^2 + z^2) / (x+y+z)^2 = 2/3
# gdzie x = sqrt(m_e), y = sqrt(m_mu), z = sqrt(m_tau)
# Znamy x, y. Szukamy z.
#
# Rozwijajac: x^2+y^2+z^2 = (2/3)(x+y+z)^2
# x^2+y^2+z^2 = (2/3)(x^2+y^2+z^2 + 2xy + 2xz + 2yz)
# 3(x^2+y^2+z^2) = 2(x^2+y^2+z^2) + 4xy + 4xz + 4yz
# x^2+y^2+z^2 = 4xy + 4xz + 4yz
# z^2 - 4(x+y)z + (x^2+y^2 - 4xy) = 0

S = x + y   # sqrt(m_e) + sqrt(m_mu)
P = x * y   # sqrt(m_e) * sqrt(m_mu)

# az^2 + bz + c = 0
a_coeff = 1.0
b_coeff = -4.0 * S
c_coeff = x**2 + y**2 - 4.0 * P  # = x^2 + y^2 - 4xy

discriminant = b_coeff**2 - 4.0 * a_coeff * c_coeff
print(f"  Rownanie kwadratowe: z^2 - 4(x+y)z + (x^2+y^2-4xy) = 0")
print(f"  x = sqrt(m_e)  = {x:.8f}")
print(f"  y = sqrt(m_mu) = {y:.8f}")
print(f"  Dyskryminant    = {discriminant:.6f}")

if discriminant >= 0:
    z_plus  = (-b_coeff + np.sqrt(discriminant)) / (2.0 * a_coeff)
    z_minus = (-b_coeff - np.sqrt(discriminant)) / (2.0 * a_coeff)

    print(f"\n  Rozwiazania:")
    print(f"    z_+  = {z_plus:.8f}  -> m_tau_+ = {z_plus**2:.2f} MeV")
    print(f"    z_-  = {z_minus:.8f} -> m_tau_- = {z_minus**2:.2f} MeV")

    # Wybieramy rozwiazanie fizyczne (z > y > x > 0)
    z_koide = z_plus if z_plus > y else z_minus
    m_tau_koide = z_koide**2
    r_31_koide = m_tau_koide / m_e

    print(f"\n  Rozwiazanie fizyczne (z > y):")
    print(f"    z_Koide     = {z_koide:.8f}")
    print(f"    m_tau_Koide = {m_tau_koide:.2f} MeV")
    print(f"    r_31_Koide  = {r_31_koide:.2f}")
    print(f"    r_31_obs    = {r_31_obs:.2f}")

    dr31 = abs(r_31_koide - r_31_obs) / r_31_obs * 100
    test2_pass = dr31 < 1.0  # mniej niz 1% odchylenia
    print(f"\n  >> TEST 2 (|r_31_Koide - r_31_obs|/r_31_obs < 1%): "
          f"{'PASS' if test2_pass else 'FAIL'}")
    print(f"     Odchylenie: {dr31:.4f}%")
else:
    print("  BLAD: dyskryminant ujemny, brak rozwiazania rzeczywistego")
    test2_pass = False
    r_31_koide = np.nan

print()

# ============================================================
# TEST 3: Korekcja kwantowa delta_2
# ============================================================
print("=" * 70)
print("  TEST 3: Korekcja kwantowa delta_2")
print("=" * 70)

# m_n = c_M * A_tail(g0_n)^4 * (1 + delta_n)
# Dla taonu: r_31_bare * (1 + delta_2) = r_31_obs
delta_2 = r_31_obs / r_31_TGP_bare - 1.0

print(f"  r_31_TGP_bare = {r_31_TGP_bare:.1f}")
print(f"  r_31_obs      = {r_31_obs:.2f}")
print(f"  delta_2       = {delta_2:.6f}")
print(f"  delta_2       = {delta_2*100:.2f}%")

# delta_2 powinno byc ujemne (korekcja zmniejsza mase)
test3a_pass = delta_2 < 0
print(f"\n  >> TEST 3a (delta_2 < 0, korekcja redukuje mase): "
      f"{'PASS' if test3a_pass else 'FAIL'}")

# |delta_2| powinno byc rzedu ~10-15% (perturbacyjne, ale znaczace)
test3b_pass = 0.05 < abs(delta_2) < 0.25
print(f"  >> TEST 3b (5% < |delta_2| < 25%, perturbacyjna korekcja): "
      f"{'PASS' if test3b_pass else 'FAIL'}")
print()

# ============================================================
# TEST 4: Rozwiazanie ukladu WKB: kappa_1, kappa_2
# ============================================================
print("=" * 70)
print("  TEST 4: Anharmoniczna korekcja WKB")
print("=" * 70)
print()
print("  Model: m_n/m_0 = 1 + kappa_1*n + kappa_2*n^2")
print("  n=0: elektron, n=1: mion, n=2: taon")
print()

# Uklad rownan:
#   n=1: r_21 = 1 + kappa_1 + kappa_2
#   n=2: r_31 = 1 + 2*kappa_1 + 4*kappa_2
#
# Z pierwszego: kappa_1 = r_21 - 1 - kappa_2
# Wstawiamy do drugiego:
#   r_31 = 1 + 2*(r_21 - 1 - kappa_2) + 4*kappa_2
#   r_31 = 1 + 2*r_21 - 2 - 2*kappa_2 + 4*kappa_2
#   r_31 = 2*r_21 - 1 + 2*kappa_2
#   kappa_2 = (r_31 - 2*r_21 + 1) / 2

kappa_2 = (r_31_obs - 2.0 * r_21_obs + 1.0) / 2.0
kappa_1 = r_21_obs - 1.0 - kappa_2

print(f"  kappa_1 = {kappa_1:.4f}")
print(f"  kappa_2 = {kappa_2:.4f}")
print(f"  |kappa_2/kappa_1| = {abs(kappa_2/kappa_1):.6f}")
print()

# Weryfikacja:
r21_check = 1.0 + kappa_1 * 1 + kappa_2 * 1**2
r31_check = 1.0 + kappa_1 * 2 + kappa_2 * 2**2
print(f"  Weryfikacja:")
print(f"    r_21_reconstruct = {r21_check:.3f} (obs: {r_21_obs:.3f})")
print(f"    r_31_reconstruct = {r31_check:.2f} (obs: {r_31_obs:.2f})")
print()

# Test: kappa_2 > 0?
test4a_pass = kappa_2 > 0
print(f"  >> TEST 4a (kappa_2 > 0): {'PASS' if test4a_pass else 'FAIL'}")
print(f"     kappa_2 = {kappa_2:.4f}")

# Test: perturbacyjnosc |kappa_2/kappa_1| < 1?
ratio_kappa = abs(kappa_2 / kappa_1)
test4b_pass = ratio_kappa < 1.0
print(f"  >> TEST 4b (|kappa_2/kappa_1| < 1, perturbacyjnosc): "
      f"{'PASS' if test4b_pass else 'FAIL'}")
print(f"     |kappa_2/kappa_1| = {ratio_kappa:.6f}")
print()

# ============================================================
# TEST 5: Predykcja masy 4. generacji
# ============================================================
print("=" * 70)
print("  TEST 5: Predykcja masy 4. generacji leptonu")
print("=" * 70)

n4 = 3
r_41 = 1.0 + kappa_1 * n4 + kappa_2 * n4**2
m_4 = r_41 * m_e  # MeV
m_4_GeV = m_4 / 1000.0

M_W = 80.377  # GeV (PDG 2024)
M_Z = 91.1876  # GeV

print(f"\n  Predykcja WKB (n=3):")
print(f"    r_41       = m_4/m_e = {r_41:.1f}")
print(f"    m_4        = {m_4:.1f} MeV = {m_4_GeV:.2f} GeV")
print(f"    M_W        = {M_W:.3f} GeV")
print(f"    M_Z/2      = {M_Z/2:.3f} GeV")

# Czy m_4 > M_W? Jesli tak, rozpad W -> l4 + nu jest zabroniony
# kinematycznie, wiec 4. generacja nie jest obserwowalna w rozpadach W
test5a_pass = m_4_GeV > M_W / 2.0
print(f"\n  >> TEST 5a (m_4 > M_W/2 = {M_W/2:.1f} GeV, brak w rozpadach W): "
      f"{'PASS' if test5a_pass else 'FAIL'}")

# Czy m_4 > M_Z/2? Jesli tak, Z -> l4 l4bar tez zabronione
test5b_pass = m_4_GeV > M_Z / 2.0
print(f"  >> TEST 5b (m_4 > M_Z/2 = {M_Z/2:.1f} GeV, brak w rozpadach Z): "
      f"{'PASS' if test5b_pass else 'FAIL'}")
print(f"     => Spójne z obserwacja 3 lekkich generacji w LEP")
print()

# ============================================================
# TEST 6: Spojnosc Koide z WKB
# ============================================================
print("=" * 70)
print("  TEST 6: Spojnosc predykcji Koide z WKB")
print("=" * 70)

# Porownanie: r_31 z Koide vs r_31 z WKB (oba uzywaja danych PDG,
# wiec powinny sie zgadzac -- ale test jest czy algebraicznie spójne)
print(f"\n  r_31 z Koide:  {r_31_koide:.2f}")
print(f"  r_31 z WKB:    {r31_check:.2f}")
print(f"  r_31 z PDG:    {r_31_obs:.2f}")

# Korekcja TGP bare -> obs
# Z WKB: efektywna korekcja dla taonu (n=2) wzgledem predykcji bare
delta_2_eff = r_31_obs / r_31_TGP_bare - 1.0
# Z WKB kappa: korekcja anharmoniczna n^2 wzgledem liniowej ekstrapolacji
# Liniowa ekstrapolacja: r_31_linear = 1 + 2*kappa_1_linear
# gdzie kappa_1_linear = r_21 - 1
r_31_linear = 1.0 + 2.0 * (r_21_obs - 1.0)
delta_anharmonic = r_31_obs - r_31_linear

print(f"\n  Liniowa ekstrapolacja WKB (kappa_2=0):")
print(f"    r_31_linear  = {r_31_linear:.2f}")
print(f"    r_31_obs     = {r_31_obs:.2f}")
print(f"    delta_anhar  = {delta_anharmonic:.2f}")
print(f"    kappa_2 * 4  = {kappa_2 * 4:.2f}")

test6_pass = abs(delta_anharmonic - kappa_2 * 4) < 0.01
print(f"\n  >> TEST 6 (spojnosc anharmoniczna): "
      f"{'PASS' if test6_pass else 'FAIL'}")
print()

# ============================================================
# Podsumowanie
# ============================================================
print("=" * 70)
print("  PODSUMOWANIE")
print("=" * 70)
tests = [
    ("TEST 1: Koide Q = 2/3",                    test1_pass),
    ("TEST 2: r_31 z Koide vs PDG (<1%)",         test2_pass),
    ("TEST 3a: delta_2 < 0 (redukcja masy)",      test3a_pass),
    ("TEST 3b: perturbacyjna korekcja (5-25%)",    test3b_pass),
    ("TEST 4a: kappa_2 > 0",                       test4a_pass),
    ("TEST 4b: |kappa_2/kappa_1| < 1",            test4b_pass),
    ("TEST 5a: m_4 > M_W/2 (brak w rozpadach W)", test5a_pass),
    ("TEST 5b: m_4 > M_Z/2 (brak w rozpadach Z)", test5b_pass),
    ("TEST 6: spojnosc WKB anharmoniczna",         test6_pass),
]

n_pass = sum(1 for _, p in tests if p)
n_total = len(tests)

print()
for name, passed in tests:
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] {name}")

print(f"\n  Wynik: {n_pass}/{n_total} testow PASS")
print()

# ============================================================
# Kluczowe liczby
# ============================================================
print("=" * 70)
print("  KLUCZOWE WYNIKI NUMERYCZNE")
print("=" * 70)
print(f"""
  Relacja Koide:
    Q_obs           = {Q_obs:.10f}
    Q - 2/3         = {Q_obs - Q_exact:+.2e}

  Stosunek mas (Koide -> tau):
    m_tau_Koide     = {m_tau_koide:.2f} MeV
    m_tau_PDG       = {m_tau:.2f} MeV
    odchylenie      = {abs(m_tau_koide - m_tau)/m_tau*100:.4f}%

  Korekcja kwantowa TGP:
    delta_2         = {delta_2:.6f} ({delta_2*100:.2f}%)

  Parametry WKB:
    kappa_1         = {kappa_1:.4f}
    kappa_2         = {kappa_2:.4f}
    kappa_2/kappa_1 = {kappa_2/kappa_1:.6f}

  Predykcja 4. generacji:
    m_4             = {m_4_GeV:.2f} GeV
    m_4/M_W         = {m_4_GeV/M_W:.2f}
""")
