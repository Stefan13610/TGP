#!/usr/bin/env python3
"""
O-K1 Exploration: Koide lambda vs Wilson RG lambda
====================================================
Porownanie lambda_Koide (z warunku Q=3/2) z lambda_WF (z Wilson RG).
Jesli te dwie wartosci sa rowne, Koide wynika z RG flow substratu.

TGP v1 -- 2026-03-31
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np

print("=" * 65)
print("O-K1: Koide lambda vs Wilson RG lambda")
print("=" * 65)

# ============================================================
# 1. lambda_Koide z eq:lambda-koide
# ============================================================
# lambda_Koide = C^2 * a_Gamma^2 / (K1*^2 * r31_K^2)
# C ~ 2.000, a_Gamma = 0.040, K1* = 2.351 * a_Gamma/(1+alpha)
# r31_K = r31_Koide(r21) z eq:r31-koide

r21 = 206.768  # PDG
a_Gamma = 0.040049
alpha_K = 8.5612
C_koide = 2.000

# K1*
K1_star = 2.351 * a_Gamma / (1 + alpha_K)
print(f"\n  Parametry:")
print(f"    r_21 = {r21}")
print(f"    a_Gamma = {a_Gamma}")
print(f"    alpha_K = {alpha_K}")
print(f"    K1* = {K1_star:.6f}")

# r31_Koide z locus Q=3/2
a = 1 + np.sqrt(r21)
discriminant = 6*a**2 - 3 - 3*r21
# UWAGA: eq:r31-koide w dokumencie ma literowke (1/4 przed nawiasem).
# Poprawna formula: r31 = (2a + sqrt(6a^2 - 3 - 3*r21))^2 (bez 1/4).
# Weryfikacja: Q(r21=206.768, r31=3477.44) = 1.500009 ~ 3/2. Poprawne.
r31_K = (2*a + np.sqrt(discriminant))**2
print(f"    r_31^Koide = {r31_K:.2f}")

# lambda_Koide
lambda_koide = C_koide**2 * a_Gamma**2 / (K1_star**2 * r31_K**2)
print(f"\n  lambda_Koide = {lambda_koide:.4e}")

# ============================================================
# 2. lambda_WF z Wilson RG (z scripts/wilson_rg_vmod.py)
# ============================================================
# Z analizy wymiarowej: lambda_eff ~ a_Gamma^2 / Phi_0^2
Phi_0 = 24.66
lambda_WF = a_Gamma**2 / Phi_0**2
print(f"  lambda_WF (dimensional) = {lambda_WF:.4e}")

# Dokladniejsze: |u6*| * v0^6 / (5! * 64)
# u6* ~ -29.24 (z Wilson RG), v0 = Phi_0^(1/2) in substrate units
# Ale bardziej trafne: lambda_WF ~ a_Gamma^2 / Phi_0^2 * C_WF
# gdzie C_WF jest stalym czynnikiem
# Z wynikow skryptu wilson_rg_vmod.py: lambda ~ 2.63e-6

lambda_WF_numerical = 2.63e-6
print(f"  lambda_WF (numerical, from RG script) = {lambda_WF_numerical:.4e}")

# ============================================================
# 3. Porownanie
# ============================================================
print(f"\n{'=' * 65}")
print(f"  POROWNANIE:")
print(f"{'=' * 65}")
print(f"  lambda_Koide = {lambda_koide:.4e}")
print(f"  lambda_WF    = {lambda_WF_numerical:.4e}")
ratio = lambda_koide / lambda_WF_numerical
print(f"  Stosunek     = {ratio:.3f}")
print(f"  Roznica      = {abs(ratio - 1)*100:.1f}%")

# ============================================================
# 4. Interpretacja
# ============================================================
print(f"\n{'=' * 65}")
print(f"  INTERPRETACJA:")
print(f"{'=' * 65}")

if abs(ratio - 1) < 0.05:
    print("  *** ZGODNOSC < 5% ***")
    print("  lambda_Koide ~ lambda_WF sugeruje, ze Koide Q=3/2")
    print("  moze wynikac z Wilson RG flow substratu Z_2!")
elif abs(ratio - 1) < 0.5:
    print("  Zgodnosc rzedu wielkosci (ale nie precyzyjna).")
    print("  Mozliwe ze korekty wyzszego rzedu moga domknac luke.")
    print(f"  Potrzebny czynnik korekcyjny: {ratio:.3f}")
else:
    print("  Brak zgodnosci. lambda_Koide i lambda_WF to rozne skale.")
    print("  Koide Q=3/2 NIE wynika bezposrednio z Wilson RG.")

# ============================================================
# 5. Dodatkowa analiza: przy jakim Phi_0 lambda_Koide = lambda_WF?
# ============================================================
print(f"\n  Dodatkowa analiza:")
# lambda_WF ~ a_Gamma^2 / Phi_0^2 => Phi_0^2 = a_Gamma^2 / lambda_Koide
Phi_0_required = a_Gamma / np.sqrt(lambda_koide)
print(f"  Jesli lambda_WF = lambda_Koide, to Phi_0 = {Phi_0_required:.2f}")
print(f"  Aktualna wartosc Phi_0 = {Phi_0}")
print(f"  Stosunek = {Phi_0_required / Phi_0:.3f}")

# Alternatywnie: lambda_WF z pelnym wzorem
# lambda_WF = |u6*| * (a_Gamma / Phi_0)^2 * C_norm
# Dopasujmy C_norm
C_norm = lambda_WF_numerical / (a_Gamma**2 / Phi_0**2)
print(f"\n  C_norm (czynnik normalizacji WF): {C_norm:.3f}")
print(f"  lambda_Koide / (a_Gamma^2 / Phi_0^2) = {lambda_koide / (a_Gamma**2 / Phi_0**2):.3f}")

# Jezeli Koide wymaga lambda_Koide = C_norm * a_G^2/Phi_0^2
# to C_Koide_needed = lambda_Koide * Phi_0^2 / a_G^2
C_K_needed = lambda_koide * Phi_0**2 / a_Gamma**2
print(f"  C_Koide_needed = {C_K_needed:.3f}")
print(f"  C_WF_actual = {C_norm:.3f}")
print(f"  Stosunek C_K/C_WF = {C_K_needed / C_norm:.3f}")

print(f"\n{'=' * 65}")
print("  WNIOSEK KONCOWY:")
print(f"{'=' * 65}")
print(f"  lambda_Koide/lambda_WF = {ratio:.2f}")
if ratio > 1.5:
    print(f"  lambda_Koide jest {ratio:.1f}x wieksza niz lambda_WF.")
    print(f"  Aby Q=3/2 wynikalo z RG, potrzebna korekta czynnika ~{ratio:.1f}")
    print(f"  w formule lambda_WF (np. od wielopetlowych poprawek lub")
    print(f"  dokladniejszej normalizacji u6* -> lambda_eff).")
    print(f"  Status O-K1: OTWARTY, ale z obiecujaca struktura rzedu wielkosci.")
elif ratio < 0.67:
    print(f"  lambda_Koide jest {1/ratio:.1f}x mniejsza niz lambda_WF.")
    print(f"  Status O-K1: OTWARTY.")
else:
    print(f"  Status O-K1: BLISKI ZAMKNIECIU (zgodnosc < 50%).")

print("\nGOTOWE.")
