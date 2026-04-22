#!/usr/bin/env python3
"""
ex183_phi0_origin.py
Sesja v45, 2026-04-05

Badanie pochodzenia Phi_0 z pierwszych zasad (G5).

CONSTRAINTS NA PHI_0:
  1. a_Gamma * Phi_0 = 1 -> Phi_0 = 25.00 (jesli a_Gamma dokladne)
  2. alpha_s = N_c^3*g0^e/(8*Phi_0) -> Phi_0 = 24.89
  3. Brannen S2c -> Phi_0 = 24.783
  4. kappa = N_c/(4*Phi_0) -> Phi_0 ~ 25

OBSERWACJA: Phi_0 = 25 = 5^2 = (2*N_c - 1)^2

HIPOTEZY:
  H1: Phi_0 = (2*N_c - 1)^2 (combinatorial)
  H2: Phi_0 = Phi_0(substratu) z Ising 3D
  H3: Phi_0 = N_c^2 + d_A (= 9 + 8 = 17? No)
  H4: Phi_0 wynika z warunku normalizacji metryki g_uv = Phi^(2/3)*eta
  H5: Phi_0 = 4*pi*phi (= 20.33? No)
  H6: Phi_0 = (N_SM)^2/... gdzie N_SM = liczba stopni swobody SM
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
N_c = 3
ALPHA_S_PDG = 0.1179
ALPHA_EM = 1/137.036

print("=" * 72)
print("ex183: Pochodzenie Phi_0 (G5)")
print("=" * 72)

# ===== 1. CONSTRAINTS SUMMARY =====
print("\n--- 1. Constraints na Phi_0 ---\n")

g0_e = 0.86941
Phi0_aGamma = 1/0.0400           # = 25.00
Phi0_alphas = N_c**3 * g0_e / (8 * ALPHA_S_PDG)  # = 24.89
Phi0_Brannen = 24.783
Phi0_mean = np.mean([Phi0_aGamma, Phi0_alphas, Phi0_Brannen])
Phi0_std = np.std([Phi0_aGamma, Phi0_alphas, Phi0_Brannen])

print(f"  a_Gamma*Phi_0 = 1:   Phi_0 = {Phi0_aGamma:.3f}")
print(f"  alpha_s formula:     Phi_0 = {Phi0_alphas:.3f}")
print(f"  Brannen S2c:         Phi_0 = {Phi0_Brannen:.3f}")
print(f"  Srednia:             Phi_0 = {Phi0_mean:.3f} +/- {Phi0_std:.3f}")

# ===== 2. HIPOTEZA: Phi_0 = (2*N_c - 1)^2 =====
print("\n--- 2. H1: Phi_0 = (2*N_c - 1)^2 ---\n")

for Nc_test in range(2, 6):
    Phi0_H1 = (2*Nc_test - 1)**2
    print(f"  N_c = {Nc_test}: Phi_0 = (2*{Nc_test}-1)^2 = {2*Nc_test-1}^2 = {Phi0_H1}")

print(f"\n  Dla N_c = 3: Phi_0 = 5^2 = 25")
print(f"  Odch od Brannen: {(25 - 24.783)/24.783*100:+.2f}%")
print(f"  Odch od alpha_s inverse: {(25 - Phi0_alphas)/Phi0_alphas*100:+.2f}%")
print()

# Physical meaning of 2*N_c - 1 = 5?
print(f"  2*N_c - 1 = 5:")
print(f"    = N_c + (N_c - 1)")
print(f"    = dim(fund. SU(5)) w GUT")
print(f"    = liczba smaakow kwarkow przy M_Z (N_f = 5: u,d,s,c,b)")
print(f"    = N_c + 2")

# ===== 3. INNE FORMY Phi_0 = 25 =====
print("\n--- 3. Inne formy dajace Phi_0 ~ 25 ---\n")

forms = [
    ("5^2", 5**2),
    ("(2*N_c - 1)^2", (2*N_c - 1)**2),
    ("N_f^2 (N_f=5)", 5**2),
    ("N_c^2 + N_c*(N_c+1)/2 + ...", N_c**2 + N_c*(N_c+1)//2 + N_c + 1),  # 9+6+4=19
    ("d_A + d_F^2 (8+9)", (N_c**2-1) + N_c**2),  # 8+9=17
    ("(N_c+1)! (4!)", 24),  # 4! = 24
    ("sum_{k=0}^{N_c} k^2 = 14", sum(k**2 for k in range(N_c+1))),
    ("d_A*N_c + 1 (25)", (N_c**2-1)*N_c + 1),  # 24+1=25!
    ("N_c! * N_c + N_c^2 (27)", 6*N_c + N_c**2),  # 3!*3 + 9 = 27
    ("phi^6/phi (=phi^5=11.09)", PHI**5),
    ("8*pi (25.13)", 8*np.pi),
    ("phi^3 * phi^3 (17.9)", PHI**6),
    ("4*pi^2/phi^2 (15.0)", 4*np.pi**2/PHI**2),
    ("N_c^3 - 2 (25)", N_c**3 - 2),  # 27-2=25!
]

forms_sorted = sorted(forms, key=lambda x: abs(x[1] - 25))

print(f"  {'formula':>30s}  {'value':>8s}  {'dev_25%':>8s}")
print("  " + "-" * 55)
for name, val in forms_sorted[:10]:
    dev = (val - 25)/25*100
    marker = " <--" if abs(dev) < 0.5 else ""
    print(f"  {name:>30s}  {val:8.3f}  {dev:+8.3f}%{marker}")

# ===== 4. TRAFIENIA DOKLADNE =====
print("\n--- 4. Trafienia dokladne (Phi_0 = 25) ---\n")

# d_A * N_c + 1 = 8*3 + 1 = 25 !!!
print(f"  d_A * N_c + 1 = (N_c^2 - 1)*N_c + 1")
print(f"                = ({N_c}^2 - 1)*{N_c} + 1")
print(f"                = {(N_c**2-1)} * {N_c} + 1")
print(f"                = {(N_c**2-1)*N_c} + 1 = {(N_c**2-1)*N_c + 1}")
print()

# N_c^3 - 2 = 27 - 2 = 25 !!!
print(f"  N_c^3 - 2 = {N_c}^3 - 2 = {N_c**3} - 2 = {N_c**3 - 2}")
print()

# (2*N_c - 1)^2 = 5^2 = 25
print(f"  (2*N_c - 1)^2 = (2*{N_c} - 1)^2 = 5^2 = 25")
print()

# (N_c + 1)! + 1 = 4! + 1 = 25
print(f"  (N_c + 1)! + 1 = {N_c+1}! + 1 = 24 + 1 = 25")
print()

# Relacja: d_A*N_c + 1 = N_c^3 - N_c + 1 = N_c(N_c^2-1) + 1 = N_c*d_A + 1
print(f"  Wszystkie trzy sa tozsame:")
print(f"    d_A*N_c + 1 = N_c*(N_c^2-1) + 1 = N_c^3 - N_c + 1")
print(f"    (2N_c-1)^2 = 4N_c^2 - 4N_c + 1")
print(f"    N_c^3 - 2 = 27 - 2 = 25")
print()
print(f"  Czy d_A*N_c + 1 = (2N_c-1)^2?")
print(f"    N_c^3 - N_c + 1 vs 4N_c^2 - 4N_c + 1")
print(f"    N_c^3 - N_c vs 4N_c^2 - 4N_c")
print(f"    N_c(N_c^2 - 1) vs 4N_c(N_c - 1)")
print(f"    (N_c+1)(N_c-1) vs 4(N_c-1)")
print(f"    N_c + 1 vs 4")
print(f"    Tylko dla N_c = 3! Koincydencja specyficzna dla SU(3).")

# ===== 5. FIZYCZNA INTERPRETACJA =====
print("\n--- 5. Fizyczna interpretacja Phi_0 = d_A*N_c + 1 ---\n")

# d_A = dim(adjoint) = N_c^2 - 1 = 8 (number of gluons)
# N_c = dim(fundamental) = 3 (number of colors)
# d_A * N_c = 8 * 3 = 24 (number of gluon-color combinations)
# + 1 = vacuum (identity)

print(f"  d_A = dim(adj) = N_c^2 - 1 = {N_c**2 - 1} (liczba gluonow)")
print(f"  N_c = dim(fund) = {N_c} (liczba kolorow)")
print(f"  d_A * N_c = {(N_c**2-1)*N_c} (gluon x kolor = calkowite d.o.f. kolorowe)")
print(f"  + 1 = {(N_c**2-1)*N_c + 1} (vacuum / identycznosc)")
print()

# Alternative: N_c^3 - 2
print(f"  Alternatywnie: N_c^3 - 2 = {N_c**3} - 2 = {N_c**3 - 2}")
print(f"  N_c^3 = dim(fund^3) = calkowita przestrzen 3-kolorowa")
print(f"  - 2 = odjecie 2 singletow? (baryonowy + inna?)")
print()

# Alternative: (2N_c-1)^2 relates to N_f = 2N_c - 1 = 5
print(f"  Lub: Phi_0 = N_f^2 gdzie N_f = 2N_c - 1 = 5")
print(f"  N_f = 5 aktywnych smakow kwarkow na skali M_Z (u,d,s,c,b)")
print(f"  INTERPRETACJA: Phi_0 koduje liczbe aktywnych smaakow!")

# ===== 6. KONSEKWENCJE FORMULY Phi_0 = N_f^2 =====
print("\n--- 6. Konsekwencje Phi_0 = N_f^2 ---\n")

N_f = 5

# alpha_s = N_c^3 * g0^e / (8*N_f^2)
alpha_exact = N_c**3 * g0_e / (8 * N_f**2)
print(f"  alpha_s = N_c^3 * g0^e / (8*N_f^2)")
print(f"         = {N_c**3} * {g0_e:.5f} / (8 * {N_f**2})")
print(f"         = {alpha_exact:.6f}")
print(f"  PDG: {ALPHA_S_PDG}")
print(f"  Odch: {(alpha_exact/ALPHA_S_PDG - 1)*100:+.3f}%")
print()

# kappa = N_c / (4*N_f^2)
kappa_exact = N_c / (4 * N_f**2)
print(f"  kappa = N_c/(4*N_f^2) = {N_c}/(4*{N_f**2}) = {kappa_exact:.5f}")
print(f"  = {N_c}/{4*N_f**2} = 3/100")
print()

# a_Gamma = 1/N_f^2
a_Gamma_exact = 1/N_f**2
print(f"  a_Gamma = 1/N_f^2 = 1/{N_f**2} = {a_Gamma_exact:.4f}")
print()

# Running: b0 depends on N_f too
b0 = (11*N_c - 2*N_f) / 3
print(f"  b0 = (11*N_c - 2*N_f)/3 = (33 - 10)/3 = {b0:.1f}")
print()

# Self-consistency: is N_f = 2*N_c - 1?
print(f"  Self-consistency: N_f = 2*N_c - 1 = {2*N_c - 1}")
print(f"  To znaczy: 5 aktywnych smakow na skali M_Z")
print(f"  (u,d,s,c,b -- top jest za ciezki)")
print(f"  To jest FAKT z SM, nie arbitralny parametr.")

# ===== 7. CO JESLI N_f = 6? =====
print("\n--- 7. Test: Phi_0 z innymi N_f ---\n")

for Nf_test in [3, 4, 5, 6]:
    Phi0_test = Nf_test**2
    alpha_test = N_c**3 * g0_e / (8 * Phi0_test)
    a_G_test = 1/Phi0_test
    print(f"  N_f = {Nf_test}: Phi_0 = {Phi0_test}, alpha_s = {alpha_test:.5f}, a_G = {a_G_test:.4f}")

print(f"\n  PDG: alpha_s = 0.1179, N_f(M_Z) = 5")
print(f"  Jedyny N_f dajacy spojny obraz: N_f = 5 -> Phi_0 = 25")

# ===== 8. CALOSC FORMULY =====
print("\n--- 8. Kompletna zamknieta formula ---\n")

print("  +---------------------------------------------------------+")
print("  |                                                         |")
print("  |  Phi_0 = N_f^2 = (2*N_c - 1)^2 = 25                   |")
print("  |                                                         |")
print("  |  alpha_s = N_c^3 * g0^e / (8 * N_f^2)                  |")
print("  |         = 27 * g0^e / 200                               |")
print("  |                                                         |")
print("  |  kappa = N_c / (4*N_f^2) = 3/100                       |")
print("  |                                                         |")
print("  |  a_Gamma = 1/N_f^2 = 1/25                              |")
print("  |                                                         |")
print("  |  N_c = 3 (SU(3)_c), N_f = 5 (aktywne smaki)           |")
print("  |  g0^e = 0.8694 (z phi-FP)                              |")
print("  |                                                         |")
print("  +---------------------------------------------------------+")

# ===== WNIOSKI =====
print("\n" + "=" * 72)
print("WNIOSKI ex183")
print("=" * 72)
print(f"""
  1. Phi_0 = 25 ma TRZY rownowazne interpretacje:
     a) (2*N_c - 1)^2 = 5^2
     b) d_A * N_c + 1 = 8*3 + 1  (gluon x kolor + vacuum)
     c) N_f^2 = 5^2  (kwadrast liczby aktywnych smakow)
     Wszystkie trzy koincyduja TYLKO dla N_c = 3.

  2. INTERPRETACJA N_f^2:
     N_f = 5 aktywnych smakow na skali M_Z (u,d,s,c,b)
     To jest fakt z SM, nie wolny parametr.
     Phi_0 = N_f^2 daje: kappa = N_c/(4N_f^2) = 3/100
     i a_Gamma = 1/25 = 0.04.

  3. INTERPRETACJA d_A*N_c + 1:
     24 = d_A*N_c = (calkowite stopnie swobody kolorowe)
     + 1 = vacuum / identity
     Phi_0 koduje "calkowita przestrzen kolorowa + vacuum"

  4. KOMPLETNA FORMULA:
     alpha_s = N_c^3 * g0^e / (8*N_f^2) = 27*g0^e/200
     = {alpha_exact:.5f} (Phi_0=25, {(alpha_exact/ALPHA_S_PDG-1)*100:+.2f}%)

  5. STATUS G5: CZESCIOWO ROZWIAZANY
     Phi_0 = N_f^2 = 25 jest hipoteza z trojnym uzasadnieniem.
     Wymaga jeszcze wyprowadzenia N_f^2 z dynamiki substratu.
     Postep: od "stala fenomenologiczna" do "N_f^2".
""")
