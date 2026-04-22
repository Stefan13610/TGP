"""
ex15_tgp_soliton.py
====================
Analiza solitonowa TGP: czy pole Phi moze tworzyc kondensate (FDM)?

MOTYWACJA (wynik ex14):
  Proste Yukawa TGP zmniejsza G_eff, nie daje plasty krzywej rotacji.
  Alternatywa: Phi tworzy soliton (boson star) jako efektywna ciemna materia.

PYTANIE:
  Czy V(Phi) TGP dopuszcza rozwiazanie solitonowe?
  Soliton: rozwiazanie stacjonarne z Phi -> 0 dla r -> inf, Phi != 0 w centrum.

MODEL FDM (fuzzy dark matter):
  Pole skalarne Psi z masa m_b:
    i * d_t Psi = -nabla^2/(2*m_b) * Psi + m_b * V_grav * Psi
  Soliton = rozwiazanie stanu podstawowego: Psi(r,t) = phi(r) * exp(-i*omega*t)

  Profile solitonu (Schive et al 2014):
    rho_sol(r) = rho_c / (1 + 0.091*(r/r_c)^2)^8
  gdzie r_c = core radius (zalezy od m_b i masy halo).

PYTANIE DO TGP:
  Czy V(Phi) = beta/3 * Phi^3 - gamma/4 * Phi^4 moze dac soliton?
  Warunek solitonu statycznego: dE/dPhi = 0 z E_total = T + V + V_grav.
  Analiza energetyczna: czy minimum energii istnieje dla Phi_core != 0?

ANALIZA POTENCJALU V_TGP:
  V(g) = (1/3)g^3 - (1/4)g^4  [dla beta=gamma=1, g = Phi/Phi_0]
  dV/dg = g^2 - g^3 = g^2*(1-g)
  V(0) = 0, V(1) = 1/3 - 1/4 = 1/12  [minimum to V=1/12]
  V'(0) = 0, V'(1) = 0  (dwa ekstrema: g=0 niestabilne, g=1 stabilne)
  V''(0) = 0, V''(1) = 2-3 = -1 < 0  => g=1 jest MAKSIMUM!

KLUCZOWY WYNIK:
  V_TGP ma minimum w g=0 i maksimum w g=1.
  Potencjal jest "odwrocony" wzgledem standardowego phi^4.
  W standardowym phi^4: V = lambda(phi^2-v^2)^2, minimum przy phi=v (SSB).
  W TGP: V = g^3/3 - g^4/4, minimum przy g=0, maksimum przy g=1.

  => TGP nie ma spontanicznego lamania symetrii w konwencjonalnym sensie.
  => Proznia g=1 jest maksimum potencjalu.
  => Soliton ze standardowym V_TGP: phi -> 0 w centrum (g=0, NIC).

  PARADOKS: Cala teoria TGP opiera sie na prozni przy g=1 (maksimum!).
  Rownanie ruchu 2. rzedu ma rozwiazania oscylacyjne wokol g=1 (SAD point).

  Fizyczna interpretacja: TGP opisuje NIESTABILNA prozni.
  Oscylacje wokol g=1 to oscylacje wokol siodla potencjalu.
  Pole jest stabilizowane przez kinetyczny czlon f(g) (blokuje g -> 0).
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp, odeint
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

print("=" * 70)
print("EX15: Analiza solitonowa TGP i potencjal V(g)")
print("=" * 70)
print()

# ===========================================================================
# Czesc 1: Analiza potencjalu V_TGP
# ===========================================================================
print("=" * 70)
print("CZESC 1: Analiza potencjalu V(g) TGP")
print("=" * 70)
print()

g_arr = np.linspace(0.0, 1.5, 300)

def V_TGP(g, beta=1.0, gamma=1.0):
    return (beta/3.0)*g**3 - (gamma/4.0)*g**4

def dV_TGP(g, beta=1.0, gamma=1.0):
    return beta*g**2 - gamma*g**3

def d2V_TGP(g, beta=1.0, gamma=1.0):
    return 2.0*beta*g - 3.0*gamma*g**2

V_arr = V_TGP(g_arr)
dV_arr = dV_TGP(g_arr)
d2V_arr = d2V_TGP(g_arr)

print("Punkty krytyczne V_TGP(g) = (1/3)g^3 - (1/4)g^4:")
print()
print("  g=0: V(0) = 0,  V'(0) = 0,  V''(0) = 0  [infleksja]")
print(f"  g=1: V(1) = {V_TGP(1.0):.4f}, V'(1) = {dV_TGP(1.0):.4f}, V''(1) = {d2V_TGP(1.0):.4f}")
print()
print(f"  V''(1) = {d2V_TGP(1.0):.2f} < 0  => g=1 jest MAKSIMUM V!")
print(f"  V''(0) = {d2V_TGP(0.0):.2f}      => g=0 jest infleksja V")
print()
print("  V_TGP NIE ma minimum dla g > 0.")
print("  Potencjal jest niebounded below dla g < 0 i > 1.")
print()

# Porownanie z phi^4 standardowym
def V_phi4(g, lam=1.0, v=1.0):
    return lam*(g**2 - v**2)**2

print("Porownanie: Standardowe phi^4 vs TGP:")
print()
print("  phi^4:  V = lambda*(g^2-1)^2")
print("    V(0) = 1, V(1) = 0  [minimum w g=1 - SSB]")
print("    V''(1) = 8*lambda > 0  [stabilne minimum]")
print("    Soliton: 'kink' 1D lub 'bubble' 3D")
print()
print("  TGP:    V = g^3/3 - g^4/4")
print("    V(0) = 0, V(1) = 1/12  [maksimum w g=1!]")
print("    V''(1) = -1 < 0  [niestabilny punkt siodlowy]")
print("    Brak solitonu w sensie topologicznym.")
print()

# Ogon efektywny: V_eff(delta) dla delta = 1-g
print("Efektywny potencjal dla oscylacji delta = 1-g wokol g=1:")
print("  V_TGP(1-delta) ~ 1/12 - (1/2)*delta^2 + (2/3)*delta^3 - ...")
print("  V_eff(delta) = -1/2*delta^2 + ...  [odwrocona harmonia!]")
print()
print("  => Rownanie oscylacji: delta'' + delta = source (NIE: delta'' - delta)")
print("  => Nie daje Yukawa (exp decay), daje oscylacje sin(r)/r.")
print("  => To potwierdza wynik ex11 analitycznie!")
print()

# ===========================================================================
# Czesc 2: Warunek solitonu
# ===========================================================================
print("=" * 70)
print("CZESC 2: Warunek istnienia solitonu")
print("=" * 70)
print()
print("Soliton statyczny wymaga (Vakhitov-Kolokolov kryterium):")
print("  dE/dg = 0 (rownanie ODE)")
print("  E_total < E_rozproszony (energetyczna korzystnosc)")
print()
print("  Energia solitonu TGP:")
print("  E = 4*pi * int_0^inf r^2 * [1/2*f(g)*(g')^2 + V(g)] dr")
print()

# Sprawdz, czy E_soliton < E_prozni
g0_vals = np.linspace(0.78, 0.999, 20)
E_vals = []

def rhs_soliton(r, y, alpha=2.0):
    g, gp = y
    fg = 1.0 + 2.0*alpha*np.log(max(g, 1e-8))
    if r < 1e-10:
        gpp = g**2*(1.0-g)/(3.0*fg) if abs(fg) > 1e-10 else 0.0
        return [gp, gpp]
    driving = g**2*(1.0-g)
    cross   = (alpha/max(g, 1e-8))*gp**2
    damp    = fg*2.0*gp/r
    gpp = (driving - cross - damp)/fg if abs(fg) > 1e-10 else 0.0
    return [gp, gpp]

r_arr = np.linspace(1e-4, 30.0, 3000)

for g0 in g0_vals[::5]:
    gp0 = g0**2*(1.0-g0)/(3.0*(1.0+2.0*2.0*np.log(max(g0,1e-8)))) * 1e-4
    sol = solve_ivp(rhs_soliton, [1e-4, 30.0], [g0, gp0],
                    t_eval=r_arr, method='DOP853',
                    rtol=1e-8, atol=1e-10, dense_output=False)
    g_sol = sol.y[0]
    gp_sol = sol.y[1]
    fg_sol = 1.0 + 2.0*2.0*np.log(np.maximum(g_sol, 1e-8))
    integrand = r_arr**2 * (0.5*fg_sol*gp_sol**2 + V_TGP(g_sol) - V_TGP(1.0))
    mask = r_arr > 0.01
    if mask.sum() > 10:
        E = 4.0*np.pi*np.trapezoid(integrand[mask], r_arr[mask])
        E_vals.append((g0, E))

print(f"  {'g0':>8s}  {'E_soliton':>14s}  {'E > 0?':>8s}")
print("-" * 36)
for g0, E in E_vals:
    print(f"  {g0:>8.4f}  {E:>14.4f}  {'TAK (niestab.)' if E > 0 else 'NIE (stabilny)'}")

print()
if all(E > 0 for _, E in E_vals):
    print("  WYNIK: Wszystkie profile TGP maja E_soliton > 0 (wzgledem prozni).")
    print("  => Solitony TGP sa ENERGETYCZNIE NIESTABILNE.")
    print("  => Pole 'chce' byc w prozni g=1, nie tworzy stabilnych solitonow.")
    print()
    print("  Wyjscie: stabilizacja przez obrot (Q-ball), lub modyfikacja V.")
else:
    print("  Znaleziono energetycznie stabilne solitony!")

print()

# ===========================================================================
# Czesc 3: Q-ball (rotating soliton)
# ===========================================================================
print("=" * 70)
print("CZESC 3: Q-ball - soliton z U(1) ladunkiem")
print("=" * 70)
print()
print("Jesli TGP ma kompleksowe pole Phi = phi*exp(i*omega*t),")
print("soliton moze byc stabilizowany przez ladunek U(1) (Q-ball).")
print()
print("  Warunek Q-ball: omega^2 < V(phi)/phi^2  dla pewnego phi")
print("  Czyli: istnieje omega taki, ze V_eff(phi, omega) = V - (omega^2/2)*phi^2")
print("         ma globalne minimum w phi > 0.")
print()
print("  Dla V_TGP(g) = g^3/3 - g^4/4:")
print("  V_eff(g,omega) = g^3/3 - g^4/4 - (omega^2/2)*g^2")
print()
print("  Warunek minimum V_eff:")
print("  dV_eff/dg = g^2 - g^3 - omega^2*g = 0")
print("  => g^2 - g^3 = omega^2*g  [dla g>0]")
print("  => g - g^2 = omega^2  (dla g in (0,1))")
print()
print("  Maksimum g-g^2 dla g=0.5: max = 0.5 - 0.25 = 0.25")
print("  => Q-ball istnieje dla omega^2 < 0.25, tj. omega < 0.5.")
print()

# Rownanie Q-ball: dV_eff/dg = 0
omega_vals = np.linspace(0.0, 0.5, 50)
g_qball = []
for omega in omega_vals:
    # g - g^2 = omega^2 => g^2 - g + omega^2 = 0 => g = (1 +- sqrt(1-4*omega^2))/2
    disc = 1.0 - 4.0*omega**2
    if disc >= 0:
        g_plus  = (1.0 + np.sqrt(disc))/2.0
        g_minus = (1.0 - np.sqrt(disc))/2.0
        g_qball.append((omega, g_minus))  # mniejszy korzen: stabilny Q-ball

print(f"  {'omega':>8s}  {'g_qball':>10s}  {'V_eff(g_Q)':>12s}")
print("-" * 36)
for omega, gq in g_qball[::5]:
    Veff = V_TGP(gq) - (omega**2/2)*gq**2
    print(f"  {omega:>8.3f}  {gq:>10.4f}  {Veff:>12.4f}")

print()
print("  WYNIK: Q-ball TGP istnieje dla omega in (0, 0.5).")
print("  g_qball ~ 0 dla malych omega => soliton slaby.")
print("  g_qball = 0.5 dla omega = 0.5 => max ladunek.")
print()

# ===========================================================================
# Czesc 4: Mapowanie na FDM
# ===========================================================================
print("=" * 70)
print("CZESC 4: Mapowanie Q-ball TGP na FDM")
print("=" * 70)
print()
print("Q-ball TGP z omega ~ Omega_galaktyki:")
print()

kpc_m    = 3.086e19
M_sun    = 1.989e30
G_SI     = 6.674e-11
l_Pl     = 1.616e-35
m_Pl_kg  = 2.176e-8
hbar_SI  = 1.055e-34
c_SI     = 3.0e8

# Skala galaktyczna: r ~ 10 kpc
r_gal_m  = 10.0 * kpc_m
T_gal_s  = 2.5e8 * 3.156e7   # 250 Myr orbital period
omega_gal = 2.0*np.pi / T_gal_s  # [rad/s]

# Masa bozonu: hbar * omega = m * c^2 => m = hbar*omega/c^2
m_boson_SI = hbar_SI * omega_gal / c_SI**2
m_boson_eV = m_boson_SI * c_SI**2 / 1.602e-19

print(f"  T_galaktyki ~ 250 Myr => omega = {omega_gal:.3e} rad/s")
print(f"  m_boson = hbar*omega/c^2 = {m_boson_SI:.3e} kg")
print(f"           = {m_boson_eV:.3e} eV/c^2")
print()

m_sp_FDM = m_boson_SI / m_Pl_kg * c_SI / l_Pl  # m_sp w jednostkach Plancka
print(f"  m_sp [Planck] = {m_sp_FDM:.3e}")
print()

# Wymaganie FDM: m ~ 10^-22 eV
print("  Standardowe FDM: m_boson ~ 1e-22 eV/c^2")
print(f"  TGP galaktyczne Q-ball: m_boson ~ {m_boson_eV:.1e} eV/c^2")
print()
ratio = m_boson_eV / 1e-22
print(f"  Roznica: {ratio:.1e} razy wieksze od standardowego FDM.")
print()
print("  Interpretacja:")
print("  TGP Q-ball na skali galaktycznej odpowiada bozonowi znacznie")
print("  ciezszemu niz ultralight FDM ~ 1e-22 eV.")
print()
print("  Dla m_boson ~ 1e-22 eV: T_obrotu ~ 1e-22/hbar c^2 => astronomicznie dlugi.")
print("  Q-ball TGP z tak malym omega prawie sie nie krecl.")
print()

# ===========================================================================
# Czesc 5: Porownanie modeli ciemnej materii
# ===========================================================================
print("=" * 70)
print("CZESC 5: TGP - gdzie pasuje w panoramie CD?")
print("=" * 70)
print()
print("  Model         m_DM          Skala core     TGP analog?")
print("-" * 65)
dm_models = [
    ("CDM (WIMP)",     "1-1000 GeV",  "< 0.01 pc",  "Nie - inne sily"),
    ("Axion CDM",      "1e-6-1e-3 eV","< 1 pc",     "Blisko (ale m=mikroeV)"),
    ("FDM ultralight", "1e-22 eV",    "~ 1 kpc",    "Moze (Q-ball TGP)"),
    ("WDM (steryl nu)","1-10 keV",    "~ 0.3 kpc",  "Nie - ferm., nie skalar"),
    ("Boson star",     "dowolne m",   "~ M/m^2",    "TGP Q-ball (wymaga omega<0.5)"),
]
for model, mass, core, tgp in dm_models:
    print(f"  {model:20s}  {mass:12s}  {core:14s}  {tgp}")

print()
print("  KONKLUZJA:")
print("  TGP Q-ball (obracajace sie pole Phi z ladunkiem U(1))")
print("  jest KONCEPTUALNIE mozliwym FDM kandydatem,")
print("  ale wymaga:")
print("  1. Rozszerzenia TGP o U(1) globalna symetrie (Phi -> e^{i*theta}*Phi)")
print("  2. Odpowiedniego m_boson (z odrobina dostrajania)")
print("  3. Mechanizmu produkcji Q-ball w kosmologii")
print()

# ===========================================================================
# Czesc 6: Modyfikacja V(Phi) - stabilny soliton
# ===========================================================================
print("=" * 70)
print("CZESC 6: Modyfikacja V dla stabilnych solitonow")
print("=" * 70)
print()
print("Zamiast V_TGP = g^3/3 - g^4/4, rozwazmny:")
print("  V_mod(g) = V_TGP(g) + eps * g^2")
print("           = eps*g^2 + g^3/3 - g^4/4")
print()
print("  dV_mod/dg = 2*eps*g + g^2 - g^3 = 0")
print("  => g=0 lub 2*eps + g - g^2 = 0")
print("  => g = (1 +- sqrt(1-8*eps)) / 2")
print()
print("  Dla eps < 1/8: istnieje minimum V_mod w g_min < 1.")
print()

eps_vals = [0.0, 0.01, 0.05, 0.10, 0.124]
print(f"  {'eps':>8s}  {'g_min':>10s}  {'V_mod(g_min)':>14s}  {'V''(g_min)':>12s}")
print("-" * 50)
for eps in eps_vals:
    disc = 1.0 - 8.0*eps
    if disc >= 0:
        g_min_mod = (1.0 - np.sqrt(disc))/2.0
        Vm = eps*g_min_mod**2 + V_TGP(g_min_mod)
        dVm = 2.0*eps*g_min_mod + dV_TGP(g_min_mod)
        d2Vm = 2.0*eps + d2V_TGP(g_min_mod)
        print(f"  {eps:>8.3f}  {g_min_mod:>10.4f}  {Vm:>14.4f}  {d2Vm:>12.4f}")
    else:
        print(f"  {eps:>8.3f}  {'brak':>10s}  {'---':>14s}  {'---':>12s}")

print()
print("  Dla eps ~ 0.05-0.1: stabilne minimum V_mod przy g ~ 0.1-0.3.")
print("  Soliton z g0 ~ g_min jest energetycznie stabilny!")
print()
print("  Interpretacja fizyczna epsilon:")
print("  eps ~ masa skwadratowana bosonu w sektorze S0")
print("  (masa efektywna pola TGP dla malych g)")
print()
print("  ROZSZERZENIE TGP: dodanie eps*Phi^2 do dzialan => stabilny soliton FDM.")
print("  Cena: dodatkowy parametr eps (ale tylko jeden!).")
print()

# ===========================================================================
# Podsumowanie
# ===========================================================================
print("=" * 70)
print("PODSUMOWANIE: TGP a ciemna materia")
print("=" * 70)
print()
print("1. V_TGP z beta=gamma=1 NIE ma stabilnych solitonow.")
print("   Potencjal ma maksimum (nie minimum) przy g=1.")
print("   Pole oscyluje wokol siodla => sin(r)/r (ex11).")
print()
print("2. Q-ball (obracajacy sie soliton) istnieje dla omega < 0.5.")
print("   Odpowiada baryonowi TGP z ladunkiem U(1).")
print("   Masa Q-ball nie pasuje do FDM ~ 1e-22 eV bez dostrajania.")
print()
print("3. Modyfikacja V_mod = eps*g^2 + V_TGP daje stabilny soliton.")
print("   Wymaga jednego dodatkowego parametru eps ~ 0.05-0.1.")
print("   To jest minimalne rozszerzenie TGP zdolne do FDM.")
print()
print("4. TGP (Droga B, m_sp > 3e-41 Planck) jest zgodne z obserwacjami")
print("   promieniowania, ale NIE wyjasnia plasty krzywej rotacji.")
print("   Wymaga standardowego CDM jako dodatkowego skladnika.")
print()
print("REKOMENDACJA:")
print("  Najbardziej realistyczny scenariusz TGP:")
print("  - Grawitacja = GR (klasyczna)")
print("  - Ciemna materia = standardowe CDM lub FDM (niezalezne od TGP)")
print("  - TGP = nowa sila skalarna na skalach sub-Plancka / subnuklearnych")
print("  - Efekty TGP: trojcialowe korekty, progi promieniowania, N_crit")
print("  - Falsyfikacja: anomalie G(r) na skali lambda=1/m_sp")
print()

# ===========================================================================
# Wykres potencjalow
# ===========================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

g_plot = np.linspace(0.0, 1.5, 300)

# Lewy: V_TGP vs phi^4
ax = axes[0]
ax.plot(g_plot, V_TGP(g_plot), 'b-', lw=2, label=r'$V_{TGP} = g^3/3 - g^4/4$')
ax.plot(g_plot, V_phi4(g_plot), 'r--', lw=2, label=r'$V_{\phi^4} = \lambda(g^2-1)^2$')
ax.axhline(0, color='k', lw=0.5)
ax.axvline(1, color='gray', lw=0.5, ls='--', alpha=0.5)
ax.scatter([0, 1], [V_TGP(0), V_TGP(1)], c=['g', 'r'], s=80, zorder=5)
ax.set_xlabel(r'$g = \Phi/\Phi_0$', fontsize=12)
ax.set_ylabel(r'$V(g)$', fontsize=12)
ax.set_title('Potencjaly: TGP vs $\\phi^4$ SSB', fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.2, 0.8)
ax.set_xlim(0, 1.5)
ax.annotate('min TGP\n(g=0)', (0.0, 0.0), (0.15, 0.15),
            arrowprops=dict(arrowstyle='->', color='b'),
            fontsize=8, color='b')
ax.annotate('max TGP\n(g=1)', (1.0, V_TGP(1.0)), (1.15, 0.25),
            arrowprops=dict(arrowstyle='->', color='r'),
            fontsize=8, color='r')

# Srodek: V_eff dla Q-ball
ax = axes[1]
g_qb = np.linspace(0.0, 1.0, 200)
for omega_v, col, lbl in [(0.0, 'k', r'$\omega=0$'), (0.2, 'b', r'$\omega=0.2$'),
                            (0.4, 'g', r'$\omega=0.4$'), (0.49, 'r', r'$\omega=0.49$')]:
    Veff = V_TGP(g_qb) - 0.5*omega_v**2*g_qb**2
    ax.plot(g_qb, Veff, color=col, lw=2, label=lbl)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel(r'$g = \Phi/\Phi_0$', fontsize=12)
ax.set_ylabel(r'$V_{eff}(g,\omega) = V_{TGP} - \frac{\omega^2}{2}g^2$', fontsize=12)
ax.set_title('Q-ball: efektywny potencjal', fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.15, 0.15)
ax.set_xlim(0, 1.0)

# Prawy: V_mod z eps
ax = axes[2]
eps_plot_vals = [0.0, 0.02, 0.05, 0.10]
cols_mod = ['k', 'b', 'g', 'r']
for eps_v, col in zip(eps_plot_vals, cols_mod):
    Vmod = eps_v*g_plot**2 + V_TGP(g_plot)
    lbl = rf'$\varepsilon={eps_v}$'
    ax.plot(g_plot, Vmod, color=col, lw=2, label=lbl)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel(r'$g = \Phi/\Phi_0$', fontsize=12)
ax.set_ylabel(r'$V_{mod}(g) = \varepsilon g^2 + V_{TGP}$', fontsize=12)
ax.set_title('Modyfikacja: stabilny soliton', fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 0.25)
ax.set_xlim(0, 1.5)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), 'ex15_soliton.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres: {out}")
