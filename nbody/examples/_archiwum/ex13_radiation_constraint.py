"""
ex13_radiation_constraint.py
=============================
Ograniczenie na m_sp z obserwacji promieniowania.

PROBLEM (wynik ex12):
  Jesli m_sp -> 0, TGP przewiduje promieniowanie skalarne:
    P_TGP ~ C^2 * a^2 / (6*pi)

  Dla Ziemi na orbicie: P_TGP ~ 1.3e7 W
  Dla pulsara PSR 1913+16: P_TGP ~ 4.6e50 W

  GR pasuje do obserwacji pulsarowych z dokladnoscia 0.2%.
  Jesli TGP dodaje promieniowanie skalarne, to musimy:
  P_TGP << P_GR (quadrupole) = 7.35e24 W  [dla pulsara]

ROZWIAZANIE: progowe TGP
  Promieniowanie TGP: P = C^2*a^2/(6pi) * sqrt(1-(m_sp/omega)^2) * theta(omega-m_sp)
  Jesli m_sp > omega_orbital => P_TGP = 0  (pelna supresja)

KLUCZOWE PYTANIE:
  Jakie m_sp daje omega_orbital < m_sp dla wszystkich obserwowanych ukladow?
"""

import numpy as np

c_SI     = 3.0e8
l_Pl     = 1.616e-35
hbar_SI  = 1.055e-34
E_Pl_J   = 1.956e9 * 1.602e-19  # J
t_Pl     = l_Pl / c_SI

m_Pl_kg  = 2.176e-8
m_sun_P  = 1.989e30 / m_Pl_kg
m_earth_P = 5.972e24 / m_Pl_kg
c_factor = 1.0 / (2.0*np.sqrt(np.pi))

print("=" * 70)
print("OGRANICZENIE NA m_sp Z PROMIENIOWANIA")
print("=" * 70)
print()

# ===========================================================================
# Czestotliwosci orbitalne w jednostkach Plancka
# ===========================================================================
print("Czestotliwosci orbitalne w jednostkach Plancka:")
print()
print(f"  {'Uklad':35s}  {'T [s]':>12s}  {'omega [1/t_Pl]':>16s}")
print("-" * 68)

systems = [
    ("Ziemia/Slonce",             3.156e7),
    ("Merkury/Slonce",            7.6e6),
    ("PSR 1913+16 (Hulse-Taylor)",2.79e4),
    ("PSR J0737-3039A (double)",  8773.0),
    ("LIGO GW150914 (merger)",    0.2),
    ("Planck czas",               t_Pl),
]

for name, T_s in systems:
    omega_P = 2.0 * np.pi * t_Pl / T_s
    print(f"  {name:35s}  {T_s:12.3e}  {omega_P:16.4e}")

print()
print("Wniosek: omega_max dla obserwowanych ukladow ~ 10^-43 [1/t_Pl]")
print("         (to jest omega LIGO, merger in-band ~ 100 Hz)")
print()
print(f"         omega_LIGO = {2*np.pi*100*t_Pl:.3e} [1/t_Pl]")
print()

omega_LIGO = 2.0*np.pi * 100 * t_Pl   # 100 Hz in Planck units
omega_pulsar = 2.0*np.pi / (2.79e4 * c_SI/l_Pl)  # PSR 1913+16

print(f"         omega_pulsar = {omega_pulsar:.3e} [1/t_Pl]")
print()

# ===========================================================================
# Ograniczenie na m_sp
# ===========================================================================
print("=" * 70)
print("OGRANICZENIE: m_sp > omega_max daje P_TGP = 0")
print("=" * 70)
print()
print(f"  Aby supresowac promieniowanie TGP dla WSZYSTKICH obserwowanych ukladow:")
print()
print(f"  m_sp > omega_LIGO = {omega_LIGO:.3e}  [Planck]")
print()

lam_LIGO = l_Pl / omega_LIGO
print(f"  Odpowiadajaca skala: lambda = 1/m_sp < {lam_LIGO:.3e} m")
print(f"                               lambda < {lam_LIGO/1e-15:.1e} fm  (1 fm = 1e-15 m)")
print()

lam_nuclear = 1e-15  # 1 fm
m_sp_nuclear = l_Pl / lam_nuclear
print(f"  Dla lambda ~ 1 fm (jadro atomowe): m_sp = {m_sp_nuclear:.3e}  [Planck]")
print()
print(f"  Dla m_sp > {omega_LIGO:.0e} [Planck]:")
print(f"    - Wszystkie orbitalne czestotliwosci < m_sp")
print(f"    - P_TGP = 0 dla WSZYSTKICH znanych procesow orbitalnych")
print(f"    - TGP promieniuje tylko dla procesow o f > {lam_LIGO:.0e} Hz (Planck-scale)")
print()

# ===========================================================================
# Porownanie z m_sp = O(1) [skala Plancka]
# ===========================================================================
print("=" * 70)
print("SCENARIUSZ: m_sp ~ 1 (skala Plancka)")
print("=" * 70)
print()
print("  m_sp = 1 [Planck] => lambda = l_Planck = 1.616e-35 m")
print("  Zasieg sily TGP = 1 dlugosc Plancka")
print()
print("  Konsekwencje:")
print()
print("  1. PROMIENIOWANIE ZEROWE dla kazdego f < f_Planck = 1.86e43 Hz")
print("     => TGP nie promieniuje w zadnym obserwowalnym eksperymencie")
print("     => GR (tensor) jest jedynym promieniowaniem grawitacyjnym")
print("     => ZGODNOSC z LIGO/pulsarami")
print()

# Dla m_sp ~ 1: sila TGP zanika na skali Plancka
# V_TGP(r) ~ exp(-r/l_Pl) -> 0 dla r >> l_Pl
print("  2. ZASIEG SILY: r >> l_Planck => V_TGP ~ 0")
print("     => TGP nie modyfikuje grawitacji na ZADNYCH obserwowalnych skalach!")
print("     => Caly wplyw TGP zamkniety w sfere Plancka!")
print()
print("  3. KONSEKWENCJA: dla m_sp ~ 1, TGP to teoria KWANTOWEJ PIANKI PRZESTRZENI")
print("     Nie jest modyfikacja klasycznej grawitacji.")
print("     Cale nowe zjawiska TGP sa na skali Plancka.")
print()

# ===========================================================================
# Trzy zakresy m_sp
# ===========================================================================
print("=" * 70)
print("TRZY ZAKRESY m_sp - interpretacja fizyczna")
print("=" * 70)
print()

AU_m  = 1.496e11
pc_m  = 3.086e16
kpc_m = 3.086e19
Mpc_m = 3.086e22

m_sp_AU    = l_Pl / AU_m
m_sp_kpc   = l_Pl / kpc_m
m_sp_Mpc   = l_Pl / Mpc_m
m_sp_Hub   = l_Pl / 1.3e26

print(f"  Range 1: m_sp < {m_sp_Hub:.2e}  [lambda > Hubble radius]")
print(f"    => TGP = grawitacja Yukawa z zasiegiem > wszechswiat")
print(f"    => Nieodroznialne od 1/r^2 w obserwowalnym wszechswiecie")
print(f"    => TGP promieniuje skalarne fale (problem z pulsarami!)")
print()
print(f"  Range 2: {m_sp_Mpc:.2e} < m_sp < {m_sp_AU:.2e}  [lambda: kpc - Gpc]")
print(f"    => TGP modyfikuje grawitacje na skalach kosmologicznych")
print(f"    => Mozliwa: rotacja galaktyk, ciemna materia, anomalie CMB")
print(f"    => Promieniowanie TGP supresowane dla orbit (omega << m_sp)")
print()
print(f"  Range 3: m_sp > {m_sp_AU:.2e}  [lambda < 1 AU]")
print(f"    => TGP nie modyfikuje grawitacji solarnej ani dalszej")
print(f"    => Zasieg: subnuklearny - Planckowski")
print(f"    => Brak JAKIEGOKOLWIEK efektu w standardowej fizyce")
print(f"    => TGP = nowa sila TYLKO w sferze Plancka / jadra atomowego")
print()

# ===========================================================================
# Optymalny scenariusz dla TGP
# ===========================================================================
print("=" * 70)
print("OPTYMALNY SCENARIUSZ DLA TGP (Range 2): galaktyczny m_sp")
print("=" * 70)
print()
print("  m_sp ~ 10^-28 [Planck] => lambda ~ 1 Mpc")
print()
print("  Wtedy TGP:")
print("  1. Nie zmienia prawa Newtona w Ukladzie Slonecznym (zasieg < lambda)")
print("     Sprawdzone do 50 AU (New Horizons)")
print()
print("  2. MODYFIKUJE grawitacje na skalach galaktycznych i wiekszych")
print("     V_TGP(r) ~ -G*m^2*exp(-r/Mpc)/r  zamiast -G*m^2/r")
print("     => Zmiana predkosci orbitalnej w galaktykach!")
print()
print("  3. Promieniowanie ZEROWE dla all orbital systems < Mpc scale")
print("     omega_orbital << m_sp dla orbit gwiezdnych w galaktyce")
print()
print("  4. V3 (efekty trojcialowe) negligibly small (exp(-t) supresja)")
print()

# Oblicz efektywna predkosc krazenia dla V_TGP vs Newton
# V_TGP(r) = -4pi*C^2 * exp(-m_sp*r)/r
# F = dV/dr = -4pi*C^2 * (1/r^2 + m_sp/r) * exp(-m_sp*r)
# v^2/r = F/m => v^2 = 4pi*C^2/m * (1/r + m_sp) * exp(-m_sp*r)
# Dla C = m*c_factor: 4pi*C^2/m = 4pi*c_factor^2*m = G*m
# v_TGP^2 = G*m * (1/r + m_sp) * exp(-m_sp*r)
# vs Newton: v_N^2 = G*M/r (M = total enclosed mass, flat rotation)

m_sp_gal = l_Pl / (1e3 * kpc_m)   # lambda = 1 Mpc

r_kpc_arr = np.logspace(-1, 3, 200)  # kpc
r_m_arr   = r_kpc_arr * kpc_m
r_Pl_arr  = r_m_arr / l_Pl

# Masa galaktyki w kuli promienia r: M(<r) = const * r (flat rotation model)
M0_kg = 1e11 * 1.989e30    # 10^11 Msun enclosed within ~10 kpc
M_enclosed_kg = M0_kg * r_kpc_arr / 10.0  # liniowo dla uproszczenia
M_enclosed_P  = M_enclosed_kg / m_Pl_kg

G_SI    = 6.674e-11
v_N_SI  = np.sqrt(G_SI * M_enclosed_kg / r_m_arr)   # Newton [m/s]
v_N_kms = v_N_SI / 1e3

# TGP modification: G_eff(r) = G * exp(-m_sp * r)
v_TGP_kms = v_N_kms * np.sqrt(np.exp(-m_sp_gal * r_Pl_arr) *
                                (1.0 + m_sp_gal * r_Pl_arr))

print("  Krzywa rotacji: Newton vs TGP (lambda=1 kpc)")
print(f"  {'r [kpc]':>10s}  {'v_Newton [km/s]':>16s}  {'v_TGP [km/s]':>14s}  {'ratio':>8s}")
print("-" * 55)
for idx in [20, 40, 60, 100, 140, 180]:
    if idx < len(r_kpc_arr):
        r  = r_kpc_arr[idx]
        vN = v_N_kms[idx]
        vT = v_TGP_kms[idx]
        print(f"  {r:>10.1f}  {vN:>16.2f}  {vT:>14.2f}  {vT/vN:>8.4f}")

print()
print("  Dla lambda=1 kpc: TGP modyfikuje v_rot dla r >> 1 kpc")
print("  Efekt maleje dla duzych r (supresja eksponencjalna)")
print()

# ===========================================================================
# Podsumowanie konczacze
# ===========================================================================
print("=" * 70)
print("PODSUMOWANIE CALEGO BADANIA TGP")
print("=" * 70)
print()
print("MATEMATYCZNY WYNIK (dokl. numeryczny):")
print("  - Calka trojcialowa I_Y obliczona dokladnie przez Feynman 2D")
print("  - dI_Y/dd_ij wyznaczone analitycznie i zweryfikowane numerycznie")
print("  - Blad przyblizone sp (saddle-point): 160-770% (ZLY eksponent)")
print("  - Newton 3. zasada: |F1+F2+F3| = 2.78e-17 (masz. dokladnosc)")
print()
print("FIZYCZNY WYNIK (Droga B):")
print("  - C_i = m_i/(2*sqrt(pi)) [Planck] - jedyny param wyznaczony z G")
print("  - m_sp = wolny param, wyznacza zasieg sily TGP")
print("  - TGP jest FALSYFIKOWALNA: anomalie G(r) na skali lambda")
print()
print("WYNIK ODE DEFEKTU (ex11):")
print("  - Linearyzacja TGP daje ogon SIN(r)/r (oscylacyjny, nie Yukawa!)")
print("  - Pelne f(g) ODE nie zmienia charakteru ogona")
print("  - Yukawa jest zewn. zrodlem (aksjomat Path B), nie wynika z pola")
print()
print("OTWARTE:")
print("  1. Jaki mechanizm fizyczny daje cialo jako zrodlo Yukawa?")
print("  2. Czy m_sp ~ kpc mozna odroznić od ciemnej materii?")
print("  3. Kwantowa teoria defektow topologicznych TGP (Droga C)")
print("  4. Zwiazek TGP z grawitacja Einsteina (GR) - dwa pola")
print()
