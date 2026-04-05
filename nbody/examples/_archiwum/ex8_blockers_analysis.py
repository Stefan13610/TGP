"""
ex8_blockers_analysis.py
=========================
Analiza dwoch blokujacych pytan TGP:
  1. Jak wyznaczyc C dla realnej materii?
  2. Jak wyznaczyc m_sp?

Odpowiedzi:
  1. C jest wyznaczone przez rownowaznie z Newtonem => NIE jest wolny.
  2. m_sp jest naprawde wolny, ale daje dwie otwarte drogi.
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# ===========================================================================
# BLOKER 1: C z warunku Newtonowskiego
# ===========================================================================
print("=" * 65)
print("BLOKER 1: C dla realnej materii")
print("=" * 65)
print()
print("Potencjal TGP (czlon dominujacy, r >> 1/m_sp):")
print("  V2(r) = -4*pi*C_i*C_j / r")
print()
print("Prawo Newtona:")
print("  V_N(r) = -G*m_i*m_j / r")
print()
print("Warunek zgodnosci (jeden pomiar G wyznacza C):")
print("  4*pi*C^2 = G*m^2")
print("  C_i = m_i * sqrt(G / 4*pi)")
print()
print("W jednostkach Plancka (G = 1):")
print("  C_i = m_i / (2*sqrt(pi)) ~ 0.282 * m_i")
print()

c_factor = 1.0 / (2.0 * np.sqrt(np.pi))

masses = {
    'elektron':  4.185e-23,
    'proton':    7.685e-20,
    'Ziemia':    2.787e32,
    'Slonce':    9.137e37,
    'gal. M31':  8.0e41,
}

print(f"  {'Czialo':12s}  {'m [m_Planck]':>14s}  {'C':>14s}  {'V3/V2 (Coulomb)':>18s}")
print("-" * 65)
for name, m in masses.items():
    C = m * c_factor
    v3v2 = 4.0 * np.pi * C / 3.0
    N_crit = 3.0 / v3v2 if v3v2 > 0 else float('inf')
    print(f"  {name:12s}  {m:14.3e}  {C:14.3e}  {v3v2:18.3e}  N_crit={N_crit:.1e}")

print()
print("Wniosek nr 1:")
print("  C NIE jest wolnym parametrem.")
print("  Jest wyznaczone przez stala Newtona G i mase ciala.")
print("  Jeden pomiar (G) wiaze C ze wszystkimi masami jednoczesnie.")
print()

# ===========================================================================
# BLOKER 2: m_sp naprawde jest wolny
# ===========================================================================
print("=" * 65)
print("BLOKER 2: m_sp jako wolny parametr")
print("=" * 65)
print()
print("m_sp = sqrt(beta). Teoria nie wyznacza beta z wewnetrznych zasad.")
print()
print("Skale fizyczne odpowiadajace roznym m_sp:")
print()

m_Planck_kg = 2.176e-8
l_Planck_m  = 1.616e-35
AU_m  = 1.496e11
pc_m  = 3.086e16
Mpc_m = 3.086e22
fm_m  = 1.0e-15

scenarios = [
    ("Jadro atomowe",            140.0e6 * 1.602e-19 / (1.956e9 * 1.602e-19)),
    ("Uklad Sloneczny (50 AU)",  l_Planck_m / (50 * AU_m)),
    ("Droga Mleczna (50 kpc)",   l_Planck_m / (50 * 1e3 * pc_m)),
    ("Lokalny Wszechswiat (1 Mpc)", l_Planck_m / Mpc_m),
    ("Horyzont kosmologiczny",   l_Planck_m / (1.3e26)),
]

print(f"  {'Scenariusz':35s}  {'m_sp [m_Planck]':>18s}  {'lambda':>15s}")
print("-" * 75)
for name, m_sp_val in scenarios:
    lam_m = l_Planck_m / m_sp_val
    if lam_m > Mpc_m:
        lam_str = f"{lam_m/Mpc_m:.1e} Mpc"
    elif lam_m > 1e3 * pc_m:
        lam_str = f"{lam_m/(1e3*pc_m):.1e} kpc"
    elif lam_m > pc_m:
        lam_str = f"{lam_m/pc_m:.1e} pc"
    elif lam_m > AU_m:
        lam_str = f"{lam_m/AU_m:.1f} AU"
    elif lam_m > 1.0:
        lam_str = f"{lam_m:.1e} m"
    else:
        lam_str = f"{lam_m/fm_m:.1f} fm"
    print(f"  {name:35s}  {m_sp_val:18.3e}  {lam_str:>15s}")

print()
print("Wniosek nr 2:")
print("  m_sp jest NAPRAWDE wolny. Kazdy scenariusz daje inny charakter sily.")
print("  Potrzebny jest DRUGI pomiar (np. odchylenie od 1/r^2 na jakiejs skali).")
print()

# ===========================================================================
# CO ZOSTAJE PRZEWIDYWALNE BEZ m_sp?
# ===========================================================================
print("=" * 65)
print("CO TGP PRZEWIDUJE NIEZALEZNIE OD m_sp?")
print("=" * 65)
print()
print("1. STOSUNEK V3/V2 w granicy Coulomba (m_sp -> 0):")
print("   V3/V2 = 4*pi*gamma*C/3 = (2*sqrt(pi)/3) * m_i  [Planck]")
print()
print("   Dla N cial: N_crit = 3/(V3/V2) = 3/(2*sqrt(pi)*m_i)")
print()

m_sun_P  = 9.137e37
m_earth_P = 2.787e32

for name, m in [("Slonce", m_sun_P), ("Ziemia", m_earth_P)]:
    C = m * c_factor
    v3v2 = 4.0 * np.pi * C / 3.0
    Nc = 3.0 / v3v2
    print(f"   {name}: N_crit(Coulomb) = {Nc:.2e}")

print()
print("   => Dla cial astrofizycznych N_crit jest astronomicznie duze.")
print("   => Efekty 3-cialowe TGP sa zupelnie pomijalnie male dla")
print("      normalnej materii w granicy Coulomba.")
print()
print("2. CZLON 1/r^2 W POTENCJALE (ZAWSZE):")
print("   V2(r) = -G*m^2/r  +  8*pi*beta*C^2/r^2  +  ...")
print("         = -G*m^2/r  +  2*G*beta*m^2/(4*pi*r^2)  +  ...")
print()
print("   Ten czlon nie zanika. Dla m_sp -> 0 i skonczonego beta:")
print("   delta_V/V_Newton = 2*G*beta*m / (4*pi*r) -> 0 dla r -> inf")
print()
print("   ALE na malych skalach (r ~ 1/m_sp) jest obserwowalne odchylenie")
print("   od prawa odwrotnych kwadratow:")

for name, m in [("Slonce", m_sun_P)]:
    C = m * c_factor
    G = 4.0 * np.pi * c_factor**2 * 4.0 * np.pi  # normalization G = 4*pi in TGP units
    r_test_AU = 1.0   # 1 AU
    r_test_P  = r_test_AU * AU_m / l_Planck_m
    ratio = 8.0 * np.pi * 1.0 * C**2 / r_test_P / (4.0 * np.pi * C**2 / r_test_P)
    print(f"\n   {name} at {r_test_AU} AU: czlon 1/r^2 / czlon 1/r = {ratio:.2e}")
    print(f"   (beta=1 w jednostkach Plancka)")

print()
print("=" * 65)
print("PODSUMOWANIE BLOKUJACE")
print("=" * 65)
print()
print("  Parametr   Status       Co wyznacza")
print("  ---------  -----------  ----------------------------------")
print("  C          WYZNACZONY   C = m*sqrt(G/4pi), z pomiaru G")
print("  m_sp       WOLNY        wymaga 2. obserwacji (odchylenie od 1/r^2)")
print("  beta=gamma WOLNY        wymaga 2. obserwacji")
print()
print("  Bez m_sp: TGP to grawitacja Newtonowska + male poprawki.")
print("  Z m_sp:   TGP daje skonczona range sily + N_crit ~ 1/m_sp.")
print()
print("  Jedyne przewidywanie BEZ zadnego wolnego parametru:")
print("  => Efekty trojcialowe w granicy Coulomba sa proporcjonalne do C")
print("  => Sa pomijalnie male dla zwyklej materii (C ~ 10^{-20} - 10^{38})")
print("  => Ale moga byc duze dla hipotetycznych ciezkich cial TGP (C ~ 1)")

# Czesc ze zmiennymi nazwy
c_factor = c_factor  # fix
C_factor = c_factor  # alias used above
