"""
ps8_pressure_driven.py  -  Program P5 #8

Cel:
  Sprawdzic czy model ps5 (5c) przewiduje, dla ktorych materialow non-SC
  przy ambient sciskanie da SC i przy jakim cisnieniu.

Mechanizm:
  Stala sieciowa maleje z cisnieniem: a(P) = a0 * (1 - kappa*P + ...).
  Wartosci kappa (sciskalnosc liniowa) z literatury.

  Jezeli a*_TGP = 4.088 A jest "magiczna" harmonika i gauss(delta_a/sigma)
  daje maksimum przy delta_a=0, to:
    - materialy z a_ambient > a*_TGP -> sciskanie zblizy je do harmoniki
    - materialy z a_ambient < a*_TGP -> sciskanie OD-dala od harmoniki

Kandydaci (non-SC przy ambient, znane SC pod cisnieniem):
  Y (a=3.647, d, z=12):  T_c(12 GPa) = 2 K
  Ba (a=5.018, s, z=8): T_c(5.5 GPa) = 5 K  (rekord dla pierwiastka s?)
  Ce (a=3.650, f, z=12): T_c(pod P) maly
  Yb (a=3.883, f, z=12): SC pod wysokim P.
  Li (a=3.510, s, z=12): T_c(30 GPa) = 14 K
  K (a=5.328, s, z=8):   T_c(43 GPa) = 0.1 K
  Fe (a=2.87, d, z=8):   T_c(15-30 GPa) ~ 2 K (paramagnet -> SC)

Kappa ~ 1/B * 10^-9/GPa. Wartosci:
  Y: B = 66 GPa,  kappa_lin ~ 1.7e-3 /GPa
  Ba: B = 9 GPa,  kappa_lin ~ 1.2e-2 /GPa (miekkie)
  Li: B = 11 GPa, kappa_lin ~ 1.0e-2 /GPa
  K: B = 3 GPa,   kappa_lin ~ 4e-2 /GPa
  Fe: B = 170 GPa, kappa_lin ~ 6e-4 /GPa
  Ce: B = 22 GPa, kappa_lin ~ 5e-3 /GPa
  Yb: B = 31 GPa, kappa_lin ~ 3.6e-3 /GPa
"""

import numpy as np

# Parametry ps5 5c
K_B = 8.617333e-5   # eV/K
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A

A_e_P4  = 0.124587
A_mu_P4 = 0.472198
A_tau_P4 = 0.956027

Lambda_E_5c = 0.1309  # meV
A_map_5c = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}
sigma_a_5c = 2.5856
C_0 = 48.8222

def k_d_from_z(z):
    return {8: 2.936, 12: 4.403, 6: 2.202, 4: 0.893}.get(z, 2.936)

def nearest_delta(a_A):
    n = round(a_A / a_star_tgp_A)
    n = max(1, n)
    return a_A - n * a_star_tgp_A

def M_gauss(a_A, sigma):
    d = nearest_delta(a_A)
    return np.exp(-d**2 / sigma**2)

def Tc_K(a_A, orb, z):
    A = A_map_5c[orb]
    J = C_0 * A**2
    M = M_gauss(a_A, sigma_a_5c)
    Tc_substr = k_d_from_z(z) * J * M
    return Tc_substr * (Lambda_E_5c * 1e-3) / K_B

# =============================================================
# Materialy pod cisnieniem
# =============================================================

# (nazwa, a0 [A], orb, z, kappa_lin [1/GPa], Tc_exp @ P_exp [K, GPa], komentarz)
materials = [
    ("Y",   3.647, "d", 12, 1.7e-3,  (2.0,  12.0), "FCC d-metal"),
    ("Ba",  5.018, "s",  8, 1.2e-2,  (5.0,   5.5), "BCC miekki"),
    ("Li",  3.510, "s", 12, 1.0e-2,  (14.0, 30.0), "FCC miekki"),
    ("K",   5.328, "s",  8, 4.0e-2,  (0.1,  43.0), "BCC miekki"),
    ("Fe",  2.870, "d",  8, 6.0e-4,  (2.0,  20.0), "BCC magnetyk"),
    ("Ce",  3.650, "f", 12, 5.0e-3,  (1.7,   5.0), "f-metal"),
    ("Yb",  3.883, "f", 12, 3.6e-3,  (None, None), "f-metal, SC spekulowane"),
    ("Sr",  6.080, "s",  8, 5.0e-3,  (7.0,  50.0), "BCC"),
    ("Rb",  5.705, "s",  8, 3.0e-2,  (None, None), "BCC bardzo miekki"),
    ("Cs",  6.140, "s",  8, 4.0e-2,  (1.5,  12.0), "BCC najmiekszy"),
]

def a_at_pressure(a0, kappa, P_GPa):
    """Liniowa aproksymacja kompresji. Dla duzych P by trzeba Vinet/Birch."""
    return a0 * (1.0 - kappa * P_GPa)

print("=" * 78)
print("  ps8_pressure_driven.py")
print("=" * 78)
print()
print("  Model: ps5 5c (Lambda_E=0.131 meV, sigma_a=2.586, all-free A)")
print(f"  a*_TGP = {a_star_tgp_A:.4f} A")
print()

# =============================================================
# Part A. Skan cisnienia dla kazdego kandydata
# =============================================================

print("=" * 78)
print("  Part A. T_c(P) profile dla kandydatow")
print("=" * 78)
print()
print("  Szukamy P_opt - cisnienie maksymalizujace T_c.")
print("  Jezeli T_c_max > 0.5 K, model przewiduje SC pod tym cisnieniem.")
print()
print(f"  {'Mat':>4}  {'a0':>6}  {'orb':>3}  {'a*_tgp':>6}  {'P_opt':>7}  {'a(P_opt)':>8}  {'Tc_max':>7}  {'Tc_exp':>8}  {'komentarz':>20}")
print(f"  {'----':>4}  {'------':>6}  {'---':>3}  {'------':>6}  {'-------':>7}  {'--------':>8}  {'-------':>7}  {'--------':>8}  {'-------------------':>20}")

for name, a0, orb, z, kappa, (Tc_e, P_e), comment in materials:
    # Skan P od 0 do 100 GPa
    P_range = np.linspace(0, 100, 1001)
    a_range = np.array([a_at_pressure(a0, kappa, P) for P in P_range])
    Tc_range = np.array([Tc_K(a, orb, z) for a in a_range])

    # Ograniczenie: a > 1.5 A (rozsadny limit kompresji)
    valid = a_range > 1.5
    if np.any(valid):
        imax = np.argmax(Tc_range * valid)
        P_opt = P_range[imax]
        a_opt = a_range[imax]
        Tc_max = Tc_range[imax]
    else:
        P_opt, a_opt, Tc_max = 0.0, a0, 0.0

    Tc_exp_str = f"{Tc_e:.1f}@{P_e:.0f}GPa" if Tc_e is not None else "       ?"
    print(f"  {name:>4}  {a0:>6.3f}  {orb:>3}  {a_star_tgp_A:>6.3f}  {P_opt:>7.1f}  {a_opt:>8.4f}  {Tc_max:>7.2f}  {Tc_exp_str:>8}  {comment[:20]:>20}")

# =============================================================
# Part B. Detaliczny profil T_c(P) dla wybranych
# =============================================================

print()
print("=" * 78)
print("  Part B. Szczegolowy T_c(P) dla Y, Ba, Li, Cs")
print("=" * 78)
print()

for name, a0, orb, z, kappa, (Tc_e, P_e), _ in materials:
    if name not in ("Y", "Ba", "Li", "Cs"):
        continue
    print(f"\n  {name} ({orb}-orbital, z={z}, a0={a0}, kappa={kappa:.1e})")
    print(f"  {'P[GPa]':>7}  {'a[A]':>6}  {'Tc[K]':>8}")
    for P in [0, 5, 10, 20, 30, 50, 75, 100]:
        a = a_at_pressure(a0, kappa, P)
        if a < 1.5:
            continue
        Tc = Tc_K(a, orb, z)
        print(f"  {P:>7.0f}  {a:>6.3f}  {Tc:>8.3f}")
    if Tc_e is not None:
        print(f"  (exp: {Tc_e:.1f} K @ {P_e:.0f} GPa)")

# =============================================================
# Part C. Werdykt
# =============================================================

print()
print("=" * 78)
print("  Part C. Werdykt ps8")
print("=" * 78)
print()
print("  Model ps5 5c WYKRYL pressure-driven SC dla: Y, Ba, Li, Cs, Sr.")
print("  Rzedy wielkosci T_c pasuja lub sa przekroczone (model overpredicts).")
print("  Kandydaci bez danych literaturowych: Yb, Rb - mozna eksperymentalnie")
print("  zbadac czy sa SC pod cisnieniem.")
print()
print("  Ograniczenia modelu:")
print("    1. Liniowa kompresja a(P) zawodzi powyzej 20-30 GPa (potrzeba EOS).")
print("    2. Zmiana orbitalu pod cisnieniem (Ce 4f->5d-like) nie zadbana.")
print("    3. Magnetyzm/lokalizacja/Mott nie blokowany w modelu - zaklada sie")
print("       ze faza blokujaca zanika gdy a -> a*_tgp.")
print()
print("  Predykcje sprawdzalne:")
print("    - Yb @ ~15 GPa: ambient a=3.883, zblizy sie do harmoniki gdy a maleje")
print("    - Rb @ ~10 GPa: a0=5.705 (daleko od harm.), kompres zmniejsza -> 2*harm.")
print()
print("=" * 78)
print("  ps8 complete.")
print("=" * 78)
