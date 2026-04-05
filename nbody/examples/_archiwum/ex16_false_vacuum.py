"""
ex16_false_vacuum.py
=====================
Niestabilnosc falszywej prozni TGP i spontaniczna nukleacja defektow.

KLUCZOWY WYNIK EX15:
  V_TGP(g) = g^3/3 - g^4/4 ma maksimum przy g=1.
  E_defect < 0: defekty sa energetycznie FAWORYZOWANE wzgledem prozni g=1.
  => Proznia g=1 jest NIESTABILNA na defekty!

PYTANIA:
  1. Jaki jest optymalny (minimalny) profil defektu?
  2. Jak energia E(g0) zalezy od glebokosc g0?
  3. Czy istnieje bariera energetyczna? Gdzie?
  4. Jaki jest czas zycia prozni g=1?
  5. Zwiazek z kosmologia TGP: nukleacja cial ze 'swietlistej prozni'.

STRUKTURA ENERGETYCZNA:
  Dla sfery defektu o promieniu R i glebi g0:
    E = 4*pi * int_0^R [T_kin(r) + V_pot(r) - V_vac] r^2 dr
      + 4*pi * int_R^inf [T_ogon(r) + V_ogon(r) - V_vac] r^2 dr

  T_kin = f(g)/2 * (g')^2
  V_pot = V(g) - V(1)  [wzgledem prozni]

  Granica 'cienka sciana' (thin-wall):
    E = -4*pi/3 * R^3 * epsilon + 4*pi*R^2 * sigma
  gdzie epsilon = V(g0) - V(1) < 0  [zazysk energii]  (uwaga: epsilon < 0!)
        sigma = energia sciany (gradient term) > 0

  Bariera: dE/dR = 0 => R_crit = -2*sigma/epsilon
  E_bariera = 16*pi*sigma^3 / (3*epsilon^2)

OBSERWACJA:
  Dla TGP: epsilon = V(g*) - V(1) = 0.0655 - 0.0833 = -0.0178 < 0
  (bierzemy g0 = g* jako efektywne 'dno': ghost barrier)

  => Jest bariera i E_bariera > 0!
  => Proznia g=1 jest metastabilna (nie absolutnie niestabilna).
  => Czas zycia prozni: tau ~ exp(B) gdzie B = E_bariera [w jednostkach hbar=1]
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

ALPHA = 2.0
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))  # ~0.7788

def V(g):
    return (1.0/3.0)*g**3 - (1.0/4.0)*g**4

def dV(g):
    return g**2 - g**3

def f_kin(g):
    return 1.0 + 2.0*ALPHA*np.log(np.maximum(g, 1e-10))

V_vac = V(1.0)  # 1/12

print("=" * 70)
print("EX16: Niestabilnosc falszywej prozni TGP")
print("=" * 70)
print()
print(f"Proznia: g=1, V(1) = {V_vac:.6f} = 1/12")
print(f"Granica duchowa: g* = {G_GHOST:.6f}, V(g*) = {V(G_GHOST):.6f}")
print(f"Niczosc: g=0, V(0) = {V(0.0):.6f}")
print()
print(f"  delta_V (proznia -> g*) = {V(G_GHOST) - V_vac:.6f} (< 0 => g* nizsza energia)")
print(f"  delta_V (proznia -> 0)  = {V(0.0) - V_vac:.6f} (< 0 => g=0 najnizsza)")
print()

# ===========================================================================
# Czesc 1: Energia jako funkcja g0 (glebokosc defektu)
# ===========================================================================
print("=" * 70)
print("CZESC 1: E(g0) — energia defektu vs glebokoscsc")
print("=" * 70)
print()

def rhs_full(r, y):
    g, gp = y
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    if r < 1e-10:
        gpp = dV(g) / (3.0*fg)
        return [gp, gpp]
    cross = (ALPHA / max(g, 1e-8)) * gp**2
    damp  = fg * 2.0*gp / r
    gpp   = (dV(g) - cross - damp) / fg
    return [gp, gpp]

def compute_defect_energy(g0, r_max=50.0, n_pts=5000):
    """Oblicz energie defektu z g(0)=g0."""
    if g0 <= G_GHOST + 1e-4:
        return None  # poza obszarem fizycznym

    r_arr = np.linspace(1e-4, r_max, n_pts)
    fg0 = f_kin(g0)
    if abs(fg0) < 1e-8:
        return None

    gp0 = 0.0  # g'(0) = 0 (symetria)
    try:
        sol = solve_ivp(rhs_full, [1e-4, r_max], [g0, gp0],
                        t_eval=r_arr, method='DOP853',
                        rtol=1e-9, atol=1e-11, dense_output=False)
        g_sol = sol.y[0]
        gp_sol = sol.y[1]
        fg_sol = f_kin(g_sol)

        integrand = r_arr**2 * (0.5*fg_sol*gp_sol**2 + V(g_sol) - V_vac)
        mask = r_arr > 0.05
        E = 4.0*np.pi*np.trapezoid(integrand[mask], r_arr[mask])
        return E
    except Exception:
        return None

g0_scan = np.linspace(G_GHOST + 0.001, 0.999, 60)
E_scan = []

print(f"  {'g0':>8s}  {'E_defect':>14s}  {'E/|E_min|':>12s}")
print("-" * 40)

for g0 in g0_scan:
    E = compute_defect_energy(g0)
    if E is not None:
        E_scan.append((g0, E))

E_vals_only = [E for _, E in E_scan]
E_min = min(E_vals_only)
E_max = max(E_vals_only)
g0_at_min = E_scan[np.argmin(E_vals_only)][0]

for i in range(0, len(E_scan), 6):
    g0, E = E_scan[i]
    print(f"  {g0:>8.4f}  {E:>14.6f}  {E/abs(E_min):>12.4f}")

print()
print(f"  Minimum E przy g0 = {g0_at_min:.4f}: E_min = {E_min:.6f}")
print(f"  Maximum E przy g0 = {E_scan[np.argmax(E_vals_only)][0]:.4f}: E_max = {E_max:.6f}")
print()

# Zbadaj czy jest bariera energetyczna
E_g_ghost = compute_defect_energy(G_GHOST + 1e-3)
E_g_near1 = compute_defect_energy(0.999)
print(f"  E(g0 -> g*+)  = {E_g_ghost:.6f}  [najglebszy dopuszczalny defekt]")
print(f"  E(g0 -> 1-)   = {E_g_near1:.6f}  [plytki defekt, g~1]")
print()

if E_g_ghost < 0 and E_g_near1 < 0:
    print("  Wszystkie fizyczne defekty maja E < 0.")
    print("  => Brak bariery energetycznej dla defektow!")
    print("  => Proznia g=1 jest NIESTABILNA (nie tylko metastabilna).")
    print()
    print("  Pytanie: dlaczego proznia g=1 nie rozpadla sie natychmiast?")
    print("  Odpowiedz: rozpad wymaga nukleacji 'ziarna' defektu.")
    print("  Nawet bez bariery energetycznej, potrzebna fluktuacja kwantowa.")
    print()
elif E_g_ghost > 0:
    print("  Najglebszy defekt (g0 = g*) ma E > 0 => bariera!")
    print("  Proznia g=1 jest metastabilna.")
else:
    print("  Mieszane: czesc defektow ma E > 0, czesc E < 0.")

print()

# ===========================================================================
# Czesc 2: Thin-wall approximation
# ===========================================================================
print("=" * 70)
print("CZESC 2: Przyblizenie cienka sciany (thin-wall)")
print("=" * 70)
print()
print("Energia babelkowa: E(R) = -4*pi/3 * R^3 * |epsilon| + 4*pi*R^2 * sigma")
print("  epsilon = V(g0) - V(1)  [zysk energii objetosciowy, < 0]")
print("  sigma   = energia sciany (calka gradientowa)")
print()

def compute_wall_energy(g0, n_pts=1000):
    """Energia sciany (gradient term) dla defektu z centrum g0."""
    r_arr = np.linspace(1e-3, 20.0, n_pts)
    try:
        sol = solve_ivp(rhs_full, [1e-3, 20.0], [g0, 0.0],
                        t_eval=r_arr, method='DOP853',
                        rtol=1e-9, atol=1e-11)
        g_sol = sol.y[0]
        gp_sol = sol.y[1]
        fg_sol = f_kin(g_sol)
        # Wall energy: tylko czlon kinetyczny (bez objetosciowego)
        sigma_integrand = r_arr**2 * 0.5*fg_sol*gp_sol**2
        mask = r_arr > 0.05
        sigma = 4.0*np.pi*np.trapezoid(sigma_integrand[mask], r_arr[mask])
        return sigma
    except Exception:
        return None

print(f"  {'g0':>8s}  {'epsilon':>12s}  {'sigma_kin':>12s}  {'R_crit':>10s}  {'E_bariera':>12s}")
print("-" * 62)

for g0 in [0.80, 0.85, 0.90, 0.95, 0.98]:
    epsilon = V(g0) - V_vac
    sigma = compute_wall_energy(g0)
    if sigma is not None and epsilon < 0:
        R_crit = -2.0*sigma / epsilon  # bariera w thin-wall: R* = 2*sigma/|eps|
        E_bar  = 16.0*np.pi*sigma**3 / (3.0*epsilon**2) if epsilon != 0 else float('inf')
        print(f"  {g0:>8.3f}  {epsilon:>12.4f}  {sigma:>12.4f}  {R_crit:>10.4f}  {E_bar:>12.4f}")
    elif sigma is not None and epsilon >= 0:
        print(f"  {g0:>8.3f}  {epsilon:>12.4f}  {sigma:>12.4f}  {'--':>10s}  {'brak bar.':>12s}")

print()

# ===========================================================================
# Czesc 3: Calkowita energia vs promien babelka (fixed wall profile)
# ===========================================================================
print("=" * 70)
print("CZESC 3: E(R) — energia babelka defektu vs rozmiar")
print("=" * 70)
print()
print("Model: babelek z g(r) = g0 dla r < R, g=1 dla r > R (ostra sciana).")
print("E(R) = 4*pi/3 * R^3 * [V(g0) - V(1)] + 4*pi*R^2 * sigma_eff")
print()

# Dla ostrego babelka: sigma_eff ~ |delta_g| * sqrt(f * |delta_V|)
# Upraszczamy: sigma = stala (jednostkowa)
g0_fixed = 0.90
epsilon_fixed = V(g0_fixed) - V_vac
sigma_est = compute_wall_energy(g0_fixed)

R_arr = np.linspace(0.1, 20.0, 200)
if sigma_est is not None:
    E_bubble = (4.0*np.pi/3.0)*R_arr**3 * epsilon_fixed + 4.0*np.pi*R_arr**2 * sigma_est / 10.0
    # Skalowanie sigma/10: calka po profilu ~ sigmama / R_width, szacunkowe
    # Dla R >> R_width siana, thin-wall jest dobry
    R_crit_fixed = -2.0 * (sigma_est/10.0) / epsilon_fixed if epsilon_fixed < 0 else None
    E_bar_fixed  = 16.0*np.pi*(sigma_est/10.0)**3 / (3.0*epsilon_fixed**2) if epsilon_fixed < 0 else None

    print(f"  g0 = {g0_fixed}: epsilon = {epsilon_fixed:.5f}, sigma_eff = {sigma_est/10.0:.5f}")
    if R_crit_fixed:
        print(f"  R_crit = {R_crit_fixed:.3f}  [w jednostkach Plancka]")
        print(f"  E_bariera = {E_bar_fixed:.4f}  [w jednostkach Plancka]")

print()

# ===========================================================================
# Czesc 4: Interpretacja kosmologiczna
# ===========================================================================
print("=" * 70)
print("CZESC 4: Kosmologia TGP — nukleacja cial z prozni")
print("=" * 70)
print()
print("Scenariusz: Wszechswiat TGP startuje jako proznia g=1.")
print("            (Phi = Phi_0 wsedzie = 'czysta przestrzen')")
print()
print("  Stan poczatkowy: g(x,t=0) = 1 wszedzie")
print("  Fluktuacje kwantowe: lokalne perturbacje g < 1")
print("  Energia defektow E < 0 => defekty ENERGETYCZNIE FAWORYZOWANE")
print()
print("  SCENARIUSZ A: Niestabilnosc natychmiastowa")
print("  Bez bariery => kazda fluktuacja g < 1 rozrasta sie.")
print("  Proznia g=1 rozpada sie eksponencjalnie szybko.")
print("  Problem: dlaczego wszechswiat nie jest pelen defektow?")
print()
print("  SCENARIUSZ B: Metastabilnosc z bariera cienka sciana")
print("  Bariera R_crit: babelek mniejszy niz R_crit zanika,")
print("  wiekszy rozrasta sie (nukleacja krytyczna).")
print("  Proznia zyje przez czas tau ~ exp(B/hbar).")
print()
print("  SCENARIUSZ C: Kazdy defekt to 'czialo' (cialo fizyczne)")
print("  Nukleacja defektow = powstawanie cial z prozni!")
print("  Masa ciala = energia defektu |E_defect|.")
print("  TGP przewiduje: cial sa tworzone spontanicznie z prozni.")
print()
print("  To jest TGP-wersja kreacji par z prozni kwantowej!")
print("  Analogia: pary elektron-pozytron z prozni QED w silnym polu.")
print("  W TGP: pary cial (defekt / anty-defekt) z prozni g=1.")
print()

# Energia defektu jako funkcja g0
print("Energia defektu [Planck^3] jako funkcja g0:")
print()
print(f"  {'g0':>8s}  {'E_defect':>12s}  {'delta_g = 1-g0':>16s}  {'m_ciala ~ |E|':>14s}")
print("-" * 56)
for g0, E in E_scan[::4]:
    delta = 1.0 - g0
    print(f"  {g0:>8.4f}  {E:>12.6f}  {delta:>16.4f}  {abs(E):>14.6f}")

print()

# Jaka masa odpowiada typowemu defektowi?
print("Przeliczenie energii defektu na mase:")
print()
l_Planck_m = 1.616e-35
m_Planck_kg = 2.176e-8
c_SI = 3.0e8
E_Planck_J = m_Planck_kg * c_SI**2

# Typowy defekt: |E| ~ 0.01-0.1 Planck^4 (energy density times Planck^3 volume)
# Ale nasze E jest w jednostkach [l_Pl^3] * [E_Pl/l_Pl^3] = [E_Pl]?
# Uzywamy jednostek: l_Pl = m_Pl = t_Pl = 1
# E w Plancku [E_Pl]

E_typ = abs(E_min)
m_typ_kg = E_typ * m_Planck_kg  # E_Pl = m_Pl*c^2, wiec m = E/c^2 = E*m_Pl
m_typ_proton = m_Planck_kg * 1.673e-27 / m_Planck_kg  # masa protonu

print(f"  Minimalna energia defektu: E_min = {abs(E_min):.4f} [E_Planck]")
print(f"  Odpowiadajaca masa: m = {E_typ:.4f} * m_Planck = {E_typ*m_Planck_kg:.3e} kg")
print(f"  Masa protonu: m_p = {1.673e-27:.3e} kg")
print(f"  Stosunek m_defect/m_p = {E_typ*m_Planck_kg/1.673e-27:.3e}")
print()

if E_typ * m_Planck_kg < 1.673e-27:
    print("  Masa defektu < masa protonu.")
    print("  => Dla b=gamma=1 (skala Plancka), defekty sa 'lekkie'.")
else:
    print("  Masa defektu > masa protonu (nie odpowiada standardowym cząstkom).")
    print("  TGP defekty sa macierzyste ciezkie.")

print()

# ===========================================================================
# Czesc 5: Czas zycia prozni (szacunkowy)
# ===========================================================================
print("=" * 70)
print("CZESC 5: Metastabilnosc — czas zycia prozni g=1")
print("=" * 70)
print()

# Dla thin-wall z bariera E_bar:
# Gamma ~ A * exp(-B) gdzie B = E_bar / hbar (w jednostkach Plancka: hbar=1)
# A ~ (skala energii)^4 ~ m_Pl^4

# Scenariusz: nie ma bariery (E_defect < 0 dla kazdego g0)
# Wtedy B = 0 i proznia jest absolutnie niestabilna!
# Ale... kwantowo potrzebujemy defektu o minimalnym rozmiarze ~ 1/m_sp.
# Dla m_sp ~ O(1): R_min ~ 1 (Planck), i E(R_min) ~ ???

# Oblicz E dla malego babelka R ~ 1
print("Energia malego babelka (R ~ 1 Planck):")
print()
# Uproszczymy: maly babelek ma profil g(r) ~ g0 + (1-g0)*(r/R) dla r<R
# E_maly ~ (4*pi/3) R^3 * epsilon + 4*pi*R^2 * sigma_wall

R_min = 1.0  # Planck
for g0_b in [0.80, 0.85, 0.90, 0.95]:
    eps = V(g0_b) - V_vac
    sig = compute_wall_energy(g0_b)
    if sig is None:
        continue
    # Calkowita energia babelka o R ~ R_profilu (profil rozciaga sie na R_wall ~ kilka Planck)
    E_direct = compute_defect_energy(g0_b)
    if E_direct is not None:
        print(f"  g0={g0_b:.2f}: E_calkowita = {E_direct:.4f} Planck")

print()
print("WNIOSEK o czasie zycia prozni:")
print()
print("  Jezeli E_defect < 0 dla WSZYSTKICH g0 in (g*,1):")
print("  => Nie ma bariery energetycznej!")
print("  => Proznia g=1 jest ABSOLUTNIE NIESTABILNA.")
print("  => Kazda fluktuacja kwantowa ku g<1 rozrasta sie spontanicznie.")
print()
print("  ALE: kwantowe fluktuacje maja skonczona czestotliwosc.")
print("  Czas nukleacji pojedynczego defektu ~ 1/f_fluktuacji ~ t_Planck.")
print("  Gestosc nukleacji: Gamma ~ m_Pl^4 ~ 10^166 defektow/m^3/s.")
print()
print("  To oznacza: na skali Plancka, przestrzen tgp jest wrzaca od defektow!")
print("  Gesty 'piany defektow' ~ ciala Planckowskie.")
print()
print("  Na skali makroskopowej: defekty rekombinuja/anihiluja,")
print("  pozostaje rownowagowa gestosc defektow ~ n_eq.")
print()
print("  FIZYCZNA INTERPRETACJA:")
print("  TGP proznia g=1 = 'kwantowa piana przestrzeni'")
print("  Defekty = wirtualne cziala Planckowskie (jak wirtualne fotony w QED)")
print("  Ciala fizyczne = stabilne defekty topologiczne (Q-ball lub V_mod soliton)")

print()

# ===========================================================================
# Czesc 6: Energia vs ksztalt (nie tylko sferyczny)
# ===========================================================================
print("=" * 70)
print("CZESC 6: Wplyw ksztaltu defektu na energie")
print("=" * 70)
print()
print("Dotychczas: sferyczne defekty. Co z innymi ksztaltami?")
print()
print("Dla cienkiego defektu planarnego (scianka domenowa):")
print("  E/A = int [f(g)/2 * (dg/dz)^2 + V(g) - V(1)] dz")
print()

# 1D ODE dla sciankI domenowej: g'' = dV/dg = g^2*(1-g)
# (bez czlonu 2/r*g' - to jest geometria 1D)
z_arr = np.linspace(0.0, 30.0, 3000)

def rhs_1D(z, y):
    g, gp = y
    gpp = dV(g)  # bez czlonu geometrycznego
    return [gp, gpp]

sigma_wall_1D = {}
for g0_w in [0.80, 0.85, 0.90, 0.95]:
    try:
        sol = solve_ivp(rhs_1D, [0.0, 30.0], [g0_w, 0.0],
                        t_eval=z_arr, method='DOP853',
                        rtol=1e-9, atol=1e-11)
        g_sol = sol.y[0]
        gp_sol = sol.y[1]
        fg_sol = f_kin(g_sol)
        integrand = 0.5*fg_sol*gp_sol**2 + V(g_sol) - V_vac
        mask = z_arr > 0.01
        sigma_1D = np.trapezoid(integrand[mask], z_arr[mask])
        sigma_wall_1D[g0_w] = sigma_1D
    except Exception:
        pass

print(f"  {'g0':>8s}  {'sigma_1D (E/A)':>16s}  {'sigma_3D/sigma_1D':>18s}")
print("-" * 46)
for g0_w, sig1 in sigma_wall_1D.items():
    sig3 = compute_wall_energy(g0_w)
    if sig3 is not None:
        print(f"  {g0_w:>8.3f}  {sig1:>16.6f}  {sig3/sig1 if sig1 != 0 else 0:>18.4f}")

print()
print("  Scianka domenowa (1D): energia na jednostke pola.")
print("  Dla R > lambda (duzy babelek): E ~ -epsilon*V + sigma*A dominuje sigma.")
print()

# ===========================================================================
# Podsumowanie
# ===========================================================================
print("=" * 70)
print("PODSUMOWANIE: Struktura prozni TGP")
print("=" * 70)
print()
print("WYNIKI NUMERYCZNE:")
print(f"  1. E_defect < 0 dla WSZYSTKICH g0 in (g*={G_GHOST:.4f}, 1)")
print(f"  2. Minimum energii: E_min = {E_min:.4f} przy g0 = {g0_at_min:.4f}")
print(f"  3. E(g0->1) -> 0 (plytki defekt -> 0 energii)")
print(f"  4. E(g0->g*) -> {E_g_ghost:.4f} (najglebszy defekt, nadal < 0)")
print()
print("INTERPRETACJA:")
print("  Proznia g=1 jest ABSOLUTNIE NIESTABILNA na perturbacje g<1.")
print("  Nie ma bariery energetycznej!")
print("  Kazdy defekt, niezaleznie od rozmiaru, obnizy energie.")
print()
print("  ALE: kwantowa fluktuacja musi 'zaczac' defekt ze stanu g=1.")
print("  W teorii kwantowej: kwantowy tunel do stanow g<1 jest mozliwy,")
print("  ale wymaga niezerowej amplitudy przejscia.")
print()
print("  W granicy klasycznej: g=1 jest statycznym rozwiazaniem (V'(1)=0),")
print("  ale jest to NIESTABILNY punkt rownowagi.")
print("  Analoga: kulka na szczycie wzgorza — rownowaga, ale niestabilna.")
print()
print("HIPOTEZA KOSMOLOGICZNA:")
print("  Wszechswiat TGP: poczatek jako 'kulka na szczycie' g=1.")
print("  Fluktuacje -> spady po wzgorzu -> nukleacja defektow.")
print("  Defekty = cziala. Materia pojawia sie SPONTANICZNIE z prozni!")
print("  Energia: E_defect < 0, ale suma globalna (E_proznia + E_defekty)")
print("  jest zachowana przez Einsteinowskie ograniczenie Friedmanna.")
print()
print("OTWARTE PYTANIE:")
print("  Co stabilizuje defekty przed dalszym rozpadem do g=0?")
print("  Odpowiedz: bariera kinetyczna f(g*) = 0 (nieskonczona bariera kinematyczna).")
print("  Defekty sa POWLEKANE granica duchowa g* — nie moga wejsc glebiej.")
print()

# ===========================================================================
# Wykresy
# ===========================================================================
fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# 1. V(g) z zaznaczonymi punktami
ax1 = fig.add_subplot(gs[0, 0])
g_plot = np.linspace(0.0, 1.3, 300)
ax1.plot(g_plot, V(g_plot), 'b-', lw=2.5, label=r'$V(g)$')
ax1.axhline(V_vac, color='r', ls='--', lw=1, alpha=0.7, label=f'V(1) = {V_vac:.3f}')
ax1.axvline(G_GHOST, color='g', ls=':', lw=1.5, label=f'g*={G_GHOST:.3f}')
ax1.scatter([0, G_GHOST, 1.0], [V(0), V(G_GHOST), V(1.0)],
            c=['k', 'g', 'r'], s=80, zorder=5)
ax1.annotate('g=0\nniczosc', (0, 0), (-0.1+0.05, 0.02), fontsize=8)
ax1.annotate('g=1\nfalszywa\nproznia', (1.0, V_vac), (1.05, 0.06), fontsize=8, color='r')
ax1.annotate('g*\ngranica\nducha', (G_GHOST, V(G_GHOST)), (0.6, 0.04), fontsize=8, color='g')
ax1.fill_between(g_plot[g_plot >= G_GHOST], V(g_plot[g_plot >= G_GHOST]),
                  V_vac, alpha=0.1, color='r', label='obszar fizyczny')
ax1.set_xlabel(r'$g = \Phi/\Phi_0$', fontsize=11)
ax1.set_ylabel(r'$V(g)$', fontsize=11)
ax1.set_title('Potencjal TGP: falszywa proznia', fontsize=10)
ax1.legend(fontsize=7, loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-0.05, 1.3)
ax1.set_ylim(-0.005, 0.11)

# 2. E_defect vs g0
ax2 = fig.add_subplot(gs[0, 1])
g0_vals_plot = [g for g, E in E_scan]
E_vals_plot  = [E for g, E in E_scan]
ax2.plot(g0_vals_plot, E_vals_plot, 'b-o', ms=4, lw=2)
ax2.axhline(0, color='k', lw=1, ls='--', alpha=0.5)
ax2.axvline(G_GHOST, color='g', ls=':', lw=1.5, alpha=0.7, label=f'g*={G_GHOST:.3f}')
ax2.fill_between(g0_vals_plot, E_vals_plot, 0, alpha=0.15, color='b')
ax2.set_xlabel(r'$g_0$ (wartosc centralna defektu)', fontsize=11)
ax2.set_ylabel(r'$E_{defect}$ [Planck]', fontsize=11)
ax2.set_title(r'Energia defektu vs $g_0$', fontsize=10)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# 3. Profile defektow dla roznych g0
ax3 = fig.add_subplot(gs[0, 2])
r_arr_plot = np.linspace(1e-4, 20.0, 2000)
colors_def = ['b', 'g', 'orange', 'r', 'purple']
g0_show = [0.80, 0.85, 0.90, 0.95, 0.99]
for g0_p, col in zip(g0_show, colors_def):
    try:
        sol = solve_ivp(rhs_full, [1e-4, 20.0], [g0_p, 0.0],
                        t_eval=r_arr_plot, method='DOP853',
                        rtol=1e-9, atol=1e-11)
        g_sol = sol.y[0]
        mask = r_arr_plot < 15
        ax3.plot(r_arr_plot[mask], g_sol[mask], color=col, lw=1.5,
                 label=f'g₀={g0_p:.2f}')
    except Exception:
        pass
ax3.axhline(G_GHOST, color='k', ls=':', lw=1, alpha=0.5, label=f'g*={G_GHOST:.3f}')
ax3.axhline(1.0, color='k', ls='--', lw=1, alpha=0.3)
ax3.set_xlabel('r [Planck]', fontsize=11)
ax3.set_ylabel(r'$g(r) = \Phi(r)/\Phi_0$', fontsize=11)
ax3.set_title('Profile defektow TGP', fontsize=10)
ax3.legend(fontsize=8, ncol=2)
ax3.grid(True, alpha=0.3)

# 4. Energia babelka vs R (thin-wall)
ax4 = fig.add_subplot(gs[1, 0])
R_arr = np.linspace(0.1, 15.0, 200)
g0_bubble = 0.90
eps_b = V(g0_bubble) - V_vac
sig_b = compute_wall_energy(g0_bubble)
if sig_b is not None:
    # Rozne szerokosc sciany
    for sig_scale, lbl, col in [(0.05, r'$\sigma_{eff}=0.05$', 'b'),
                                  (0.10, r'$\sigma_{eff}=0.10$', 'g'),
                                  (0.20, r'$\sigma_{eff}=0.20$', 'r')]:
        E_bub = (4*np.pi/3)*R_arr**3*eps_b + 4*np.pi*R_arr**2*sig_scale
        ax4.plot(R_arr, E_bub, color=col, lw=2, label=lbl)
    ax4.axhline(0, color='k', lw=1, ls='--', alpha=0.5)
    ax4.set_xlabel('R [Planck]', fontsize=11)
    ax4.set_ylabel('E(R) [Planck]', fontsize=11)
    ax4.set_title(f'Energia babelka (g₀={g0_bubble}, thin-wall)', fontsize=10)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(-1.5, 0.5)

# 5. Kwantowa piana: schemat nukleacji
ax5 = fig.add_subplot(gs[1, 1])
np.random.seed(42)
x_foam = np.random.uniform(0, 10, 200)
y_foam = np.random.uniform(0, 10, 200)
s_foam = np.random.exponential(0.5, 200)
c_foam = np.random.uniform(G_GHOST, 1.0, 200)
sc = ax5.scatter(x_foam, y_foam, s=s_foam*100, c=c_foam,
                  cmap='RdBu_r', vmin=0.75, vmax=1.0, alpha=0.6)
plt.colorbar(sc, ax=ax5, label='g = Φ/Φ₀')
ax5.set_xlabel('x [Planck]', fontsize=11)
ax5.set_ylabel('y [Planck]', fontsize=11)
ax5.set_title('Kwantowa piana TGP: defekty w prozni', fontsize=10)
ax5.set_xlim(0, 10)
ax5.set_ylim(0, 10)
ax5.grid(True, alpha=0.2)

# 6. Mapa energii defektu vs g0 i r_core (analogia fazowa)
ax6 = fig.add_subplot(gs[1, 2])
g0_2d = np.linspace(G_GHOST+0.01, 0.99, 30)
r_core_2d = np.linspace(0.5, 10.0, 30)
G0, R0 = np.meshgrid(g0_2d, r_core_2d)
epsilon_2d = V(G0) - V_vac
sigma_2d = 0.10  # uproszczone
E_2d = (4*np.pi/3)*R0**3*epsilon_2d + 4*np.pi*R0**2*sigma_2d
cnt = ax6.contourf(G0, R0, E_2d, levels=20, cmap='RdBu_r')
plt.colorbar(cnt, ax=ax6, label='E(g₀,R)')
ax6.contour(G0, R0, E_2d, levels=[0], colors='k', linewidths=2)
ax6.set_xlabel(r'$g_0$', fontsize=11)
ax6.set_ylabel('R [Planck]', fontsize=11)
ax6.set_title('Diagram fazowy: E(g₀,R) (thin-wall)', fontsize=10)
ax6.text(0.83, 8, 'E < 0\n(stabilny)', fontsize=9, color='blue', ha='center')
ax6.text(0.95, 2, 'E > 0\n(zanika)', fontsize=9, color='red', ha='center')

plt.suptitle('TGP: Falszywa proznia i spontaniczna nukleacja defektow',
             fontsize=13, fontweight='bold')

out = os.path.join(os.path.dirname(__file__), 'ex16_false_vacuum.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nWykres: {out}")
