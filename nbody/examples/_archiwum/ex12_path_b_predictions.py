"""
ex12_path_b_predictions.py
===========================
DROGA B: TGP jako fenomenologiczna teoria z zewnetrznym zrodlem Yukawa.

UZASADNIENIE (wynik ex11):
  Pelne ODE topologicznego defektu daje ogon OSCYLACYJNY (sin(r)/r),
  nie Yukawa. Oznacza to:
  - TGP swobodne pole NIE emituje Yukawa naturalnie
  - Yukawa musi byc zadane jako zewnetrzne zrodlo (efektywny opis cial)

DROGA B - aksjomat:
  Ciala sa zewnetrznymi zrodlami dla pola TGP.
  Rownanie pola: (Laplacian - m_sp^2) * Phi = -4*pi * C_i * delta^3(r - r_i)
  Rozwiazanie: Phi_i(r) = C_i * exp(-m_sp * r) / r

  Dwa parametry:
    C_i = m_i * sqrt(G / 4*pi)  [z warunku Newtonowskiego - WYZNACZONY]
    m_sp = sqrt(beta)            [z 2. obserwacji - WOLNY]

  Przy tych zalozeniach: TGP jest W PELNI PREDYKTYWNA dla danego m_sp.

PREDYKCJE:
  1. Potencjal parowy: V2(r) = -4*pi*C_i*C_j * exp(-m_sp*r)/r
  2. Potencjal trojcialowy: V3 = -6*gamma*C_i*C_j*C_k * I_Y(d_ij, m_sp)
  3. N_crit(m_sp): po ilu cialach efekty trojcialowe dominuja
  4. Promieniowanie: P_TGP = C^2 * a^2 * omega^2 / (6*pi) * sqrt(1 - m_sp^2/omega^2)
  5. Odchylenie od 1/r^2: delta_V/V = 2*beta*C/r na malych skalach
"""

import sys, os
import numpy as np
from scipy.integrate import dblquad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# Stale fizyczne
G_SI        = 6.674e-11   # m^3 kg^-1 s^-2
c_SI        = 3.0e8       # m/s
hbar_SI     = 1.055e-34   # J*s
m_Planck_kg = 2.176e-8    # kg
l_Planck_m  = 1.616e-35   # m

# Przeliczniki
AU_m   = 1.496e11
pc_m   = 3.086e16
kpc_m  = 3.086e19
Mpc_m  = 3.086e22

# Jednostki Plancka dla mas
m_proton_kg  = 1.673e-27
m_earth_kg   = 5.972e24
m_sun_kg     = 1.989e30
m_neutron_kg = 1.675e-27

m_proton_P  = m_proton_kg  / m_Planck_kg
m_earth_P   = m_earth_kg   / m_Planck_kg
m_sun_P     = m_sun_kg     / m_Planck_kg
m_neutron_P = m_neutron_kg / m_Planck_kg

c_factor = 1.0 / (2.0 * np.sqrt(np.pi))  # C = m * c_factor [Planck]

print("=" * 70)
print("DROGA B: TGP fenomenologiczna z zrodlem Yukawa")
print("=" * 70)
print()
print("Rownanie zrodlowe: (nabla^2 - m_sp^2) * Phi = -4*pi*C_i*delta^3(r-r_i)")
print("Rozwiazanie:       Phi_i(r) = C_i * exp(-m_sp*r) / r")
print()
print("Parametry:")
print("  C_i   = m_i * sqrt(G/4*pi) = m_i/(2*sqrt(pi))  [Planck]  WYZNACZONY")
print("  m_sp  = wolny param - wymaga 2. obserwacji (zasieg sily)")
print("  gamma = 1  (tu: przyjmujemy naturalna skale)")
print()

# ===========================================================================
# Czesc 1: C dla roznych cial
# ===========================================================================
print("=" * 70)
print("CZESC 1: Wyznaczone C dla znanych cial")
print("=" * 70)
print()

bodies = [
    ("proton",         m_proton_P),
    ("neutron",        m_neutron_P),
    ("Ziemia",         m_earth_P),
    ("Slonce",         m_sun_P),
    ("10 M_sun",       10.0 * m_sun_P),
    ("10^6 M_sun (BH)", 1e6 * m_sun_P),
]

print(f"  {'Cialo':22s}  {'m [m_Planck]':>14s}  {'C [m_Planck]':>14s}  {'C/m':>8s}")
print("-" * 65)
for name, m in bodies:
    C = m * c_factor
    print(f"  {name:22s}  {m:14.4e}  {C:14.4e}  {c_factor:.4f}")

print()
print(f"  C/m = 1/(2*sqrt(pi)) = {c_factor:.6f}  dla wszystkich cial.")
print()

# ===========================================================================
# Czesc 2: N_crit jako funkcja m_sp i C
# ===========================================================================
print("=" * 70)
print("CZESC 2: N_crit(m_sp) dla roznych m_sp")
print("=" * 70)
print()
print("N_crit = liczba cial, przy ktorej energia trojcialowa ~= parowa")
print("         N_crit ~ 3 / |V3_triplet / V2_pair|")
print()

# Ladujemy modul obliczen trojcialowych
try:
    from nbody.three_body_force_exact import yukawa_overlap_exact

    def V2_yukawa(C, d, m):
        """Potencjal parowy Yukawa."""
        if m * d < 1e-10:
            return -4.0 * np.pi * C**2 / d
        return -4.0 * np.pi * C**2 * np.exp(-m * d) / d

    def V3_equilateral(C, d, m, gamma=1.0):
        """Trojcialowy potencjal dla trojkata rownobocznego."""
        I_Y = yukawa_overlap_exact(d, d, d, m, n_quad=40)
        return -6.0 * gamma * C**3 * I_Y

    def N_crit_func(C, d, m):
        """Szacunkowe N_crit."""
        v2 = V2_yukawa(C, d, m)
        v3 = V3_equilateral(C, d, m)
        if abs(v2) < 1e-30:
            return float('inf')
        ratio = abs(v3) / abs(v2)
        if ratio < 1e-15:
            return float('inf')
        return 3.0 / ratio + 2.0

    # Parametry dla protonu
    C_p = m_proton_P * c_factor
    d_test = 1.0  # 1 l_Planck jako jednostka

    print(f"  Proton: C = {C_p:.3e}  [Planck]")
    print()
    print(f"  {'m_sp':>10s}  {'lambda':>14s}  {'t=m_sp*d':>10s}  {'N_crit':>10s}")
    print("-" * 50)

    m_sp_vals = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    for m_sp in m_sp_vals:
        t = m_sp * d_test
        lam_Planck = 1.0 / m_sp if m_sp > 0 else float('inf')
        lam_m = lam_Planck * l_Planck_m
        if lam_m > Mpc_m:
            lam_str = f"{lam_m/Mpc_m:.2e} Mpc"
        elif lam_m > kpc_m:
            lam_str = f"{lam_m/kpc_m:.2e} kpc"
        elif lam_m > AU_m:
            lam_str = f"{lam_m/AU_m:.2e} AU"
        elif lam_m > 1.0:
            lam_str = f"{lam_m:.2e} m"
        else:
            lam_str = f"{lam_m*1e15:.2e} fm"

        Nc = N_crit_func(C_p, d_test, m_sp)
        Nc_str = f"{Nc:.2e}" if Nc < 1e10 else "> 1e10"
        print(f"  {m_sp:>10.3f}  {lam_str:>14s}  {t:>10.3f}  {Nc_str:>10s}")

    print()

    # N_crit dla Slonca (C ~ 10^37) jest astronomicznie duze
    C_sun = m_sun_P * c_factor
    d_AU = AU_m / l_Planck_m
    print(f"  Slonce: C = {C_sun:.3e}  [Planck]")
    print(f"  d = 1 AU = {d_AU:.3e}  [Planck]")
    print()
    print("  Dla cial makroskopowych (Slonce, Ziemia):")
    print("  t = m_sp * d_AU jest ASTRONOMICZNIE duze dla dowolnego m_sp > 0.")
    print("  => exp(-t) -> 0  => V3 jest supresowane eksponencjalnie.")
    print("  => N_crit -> inf dla cial makroskopowych.")
    print()
    print("  Efekty trojcialowe TGP sa NIEOBSERWOWALNE dla zwyklej materii.")
    print("  Jedynym wyjatkiem: m_sp -> 0 I malutkie odleglosci.")

    have_exact = True
except ImportError:
    print("  [Modul three_body_force_exact niedostepny - pomijam kalkulacje I_Y]")
    have_exact = False
    C_sun = m_sun_P * c_factor
    C_p   = m_proton_P * c_factor

print()

# ===========================================================================
# Czesc 3: Predykcja odchylenia od 1/r^2
# ===========================================================================
print("=" * 70)
print("CZESC 3: Odchylenie od prawa Newtona (1/r^2 -> poprawka)")
print("=" * 70)
print()
print("Potencjal TGP parowy (roznieza od Newtona):")
print("  V_TGP(r) = -4*pi*C^2 * exp(-m_sp*r)/r")
print("  V_N(r)   = -G*m^2/r  = -4*pi*C^2/r")
print()
print("  Roznica: delta_V = V_N - V_TGP = -4*pi*C^2 * (1 - exp(-m_sp*r))/r")
print()
print("  Na malych skalach (m_sp*r << 1):")
print("  delta_V/V_N ~ m_sp*r/2  (linearne w r)")
print()
print("  Na duzych skalach (m_sp*r >> 1):")
print("  V_TGP ~ 0  (sila zanika!)")
print()
print("  TGP PRZEWIDUJE: pomiar G na roznych skalach da rozne wyniki!")
print("  G_eff(r) = G * exp(-m_sp * r)")
print()

# Oblicz dla roznych m_sp
print(f"  {'m_sp':>10s}  {'lambda':>14s}  {'G_eff(1AU)/G':>14s}  {'G_eff(1kpc)/G':>15s}")
print("-" * 58)
m_sp_scenarios = [
    (1e-28, "Horyzont kosmologiczny"),
    (1e-22, "Skala Mpc"),
    (1e-19, "Skala kpc"),
    (1e-14, "Skala AU"),
]
for m_sp, label in m_sp_scenarios:
    lam_m = l_Planck_m / m_sp
    d_AU_P  = AU_m / l_Planck_m
    d_kpc_P = kpc_m / l_Planck_m
    G_AU  = np.exp(-m_sp * d_AU_P)
    G_kpc = np.exp(-m_sp * d_kpc_P)
    if lam_m > Mpc_m:
        lam_str = f"{lam_m/Mpc_m:.1f} Mpc"
    elif lam_m > kpc_m:
        lam_str = f"{lam_m/kpc_m:.1f} kpc"
    elif lam_m > AU_m:
        lam_str = f"{lam_m/AU_m:.1f} AU"
    else:
        lam_str = f"{lam_m:.2e} m"
    print(f"  {m_sp:>10.2e}  {lam_str:>14s}  {G_AU:>14.6f}  {G_kpc:>15.6f}")

print()
print("  Jesli m_sp jest wystarczajaco male (lambda > obserwowalne skale):")
print("  G_eff ~ G wsedzie  =>  nieodroznialne od Newtona.")
print("  Jedyna obserwowalna sygnatura: anomalia w bardzo duzych skalach.")
print()

# ===========================================================================
# Czesc 4: Promieniowanie TGP - konkretne obliczenia
# ===========================================================================
print("=" * 70)
print("CZESC 4: Promieniowanie skalarne TGP")
print("=" * 70)
print()
print("Moc promieniowania (uogolnienie Larmora dla masywnego skalaru):")
print("  P_TGP(omega) = C^2 * a^2 / (6*pi) * sqrt(1 - (m_sp/omega)^2) * theta(omega-m_sp)")
print()
print("  Dla m_sp = 0 (bezmasowy): P_TGP = C^2*a^2/(6*pi)  [analogia z GR]")
print("  Dla m_sp > omega: P_TGP = 0  (bez promieniowania ponizej progowej czest.)")
print()
print("  Calkowita moc: P_TGP = integral P_TGP(omega) d_omega")
print("  Dla ruchu krazacego z omega_0:")
print("  P_TGP = C^2 * a^2 / (6*pi) * sqrt(1 - (m_sp/omega_0)^2) * theta(omega_0 - m_sp)")
print()

# Uklad Sloneczny
m_sp_AU_P = l_Planck_m / AU_m  # m_sp odpowiadajacy lambda = 1 AU
r_AU_P    = AU_m / l_Planck_m
omega_earth_P = np.sqrt(4.0*np.pi*C_sun / r_AU_P**3) if True else 0  # Kepler
# V_Newton = -4pi*C^2/r => F = -dV/dr = -4pi*C^2/r^2 => a = 4pi*C^2/(m*r^2)
# omega^2 = a/r = 4pi*C_sun^2 / (m_earth * r_AU^3) -- ale to jest troche skomplikowane
# Prostrze: omega_earth ~ 2*pi / T_earth
T_earth_s = 3.156e7   # sekundy
T_earth_P = T_earth_s * c_SI / l_Planck_m
omega_earth_P = 2.0 * np.pi / T_earth_P

print(f"  Ziemia: omega_kepler = {omega_earth_P:.3e}  [l_Planck^-1]")
print()
print(f"  Dla m_sp < omega_kepler = {omega_earth_P:.3e}:")
print(f"    a_earth ~ v^2/r (obliczone z keplera)")
print()

# Promien krazenia, predkosc
v_earth_m_s = 2.0*np.pi * AU_m / T_earth_s
v_earth_P   = v_earth_m_s / c_SI
a_earth_P   = v_earth_P * omega_earth_P  # a = v * omega
C_earth     = m_earth_P * c_factor
P_TGP_earth_maximal = C_earth**2 * a_earth_P**2 / (6.0 * np.pi)
P_Newton_radiated = 0.0  # Newton nie promieniuje

print(f"  a_earth = {a_earth_P:.3e}  [c/l_Planck]")
print(f"  C_earth = {C_earth:.3e}  [Planck]")
print(f"  P_TGP_earth (m_sp->0) = {P_TGP_earth_maximal:.3e}  [E_Planck/t_Planck]")
print()

# Przelicz na watty
E_Planck_J = m_Planck_kg * c_SI**2
t_Planck_s = l_Planck_m / c_SI
P_SI       = P_TGP_earth_maximal * E_Planck_J / t_Planck_s
print(f"  P_TGP_earth (SI) = {P_SI:.3e}  W")
print()
print(f"  Promieniowanie grawitacyjne GR dla Ziemi (krazenie) ~ 200 W")
print(f"  Stosunek P_TGP / P_GR_szacunek ~ {P_SI / 200:.3e}")
print()

# Dla ukladu pulsar
m_psr_P = 1.4 * m_sun_P
r_psr_m = 2e8  # 2000 km
r_psr_P = r_psr_m / l_Planck_m
T_psr_s = 0.0059  # PSR 1913+16
omega_psr_P = 2.0*np.pi / (T_psr_s * c_SI / l_Planck_m)
v_psr_P = np.sqrt(4.0*np.pi * (m_psr_P*c_factor)**2 / (m_psr_P * r_psr_P))  # approx
# Prosciej: v ~ omega*r
v_psr_P = omega_psr_P * r_psr_P * (1.0/3.0)   # gruba ocena

a_psr_P = v_psr_P * omega_psr_P
C_psr   = m_psr_P * c_factor
P_TGP_psr = C_psr**2 * a_psr_P**2 / (6.0*np.pi)
P_TGP_psr_SI = P_TGP_psr * E_Planck_J / t_Planck_s

print(f"  Pulsar (PSR 1913+16):")
print(f"    C_psr = {C_psr:.3e}  [Planck]")
print(f"    P_TGP_psr (SI) ~ {P_TGP_psr_SI:.3e}  W")
print(f"    (Obserw. straty GR: ~7.35e24 W)")
print(f"    Stosunek P_TGP / P_obs ~ {P_TGP_psr_SI / 7.35e24:.3e}")
print()

# ===========================================================================
# Czesc 5: Falsyfikacja TGP - sygnatura obserwacyjna
# ===========================================================================
print("=" * 70)
print("CZESC 5: Jak falsyfikowac Droge B TGP?")
print("=" * 70)
print()
print("TGP (Droga B) jest FALSYFIKOWALNA przez:")
print()
print("1. ANOMALIA G(r): Jesli G_eff(r) = G*exp(-m_sp*r) i m_sp > 0:")
print("   => Odchylenie od 1/r^2 na skali lambda = 1/m_sp")
print("   => Poszukuj: anomalne przyspieszenia pionierow NASA (Pioneer anomaly)")
print("   => Poszukuj: anomalie w rotacji galaktyk")
print()
print("2. PROGI PROMIENIOWANIA: Fale skalarne TGP maja czestotliwosc minimalna:")
print("   => Brak promieniowania TGP dla f < m_sp/(2*pi) [Hz]")
print("   => Binarne pulsary: porownanie TGP vs GR strat energii")
print()
print("3. N_crit: Jesli m_sp ~ 1 (skala Plancka), to N_crit ~ 16 cial")
print("   => Efekty wielocialowe TGP na skali jader atomowych:")
print("   => Poszukuj: odchylenia od addytywnosci w jadro-jadro rozpraszaniu")
print()
print("4. KORELACJA CIEZAR-LADUNEK:")
print("   TGP przewiduje: C_i = m_i * const (DOKLADNIE)")
print("   => Jesli badanie wykazaloby C_i != m_i*const dla jakiejs czastki:")
print("      TGP jest falsyfikowana")
print()

# ===========================================================================
# Czesc 6: Podsumowanie predyktywnosci
# ===========================================================================
print("=" * 70)
print("CZESC 6: Podsumowanie - co TGP (Droga B) przewiduje")
print("=" * 70)
print()
print(f"  {'Wielkosc':40s}  {'Status':20s}  {'Wartosc/Wymaganie'}")
print("-" * 90)
table = [
    ("C_i (ladunek zrodlowy)", "WYZNACZONE", "C = m/(2*sqrt(pi)) [Planck]"),
    ("V2(r) (potencjal parowy)", "WYZNACZONE", "V2 = -G*m^2*exp(-m_sp*r)/r"),
    ("Zasieg sily (lambda)", "WOLNE", "wymaga 2. obserwacji"),
    ("m_sp (masa bosonu Phi)", "WOLNE", "m_sp = l_Planck/lambda"),
    ("N_crit (dla m_sp=1)", "OBLICZALNE", "N_crit ~ 16 dla C~0.15"),
    ("Promieniowanie TGP", "OBLICZALNE", "P ~ C^2*a^2/6pi * sqrt(1-m_sp^2/w^2)"),
    ("Odchylenie od 1/r^2", "OBLICZALNE", "delta_G/G ~ 1-exp(-m_sp*r)"),
    ("Efekty trojcialowe (APC)", "OBLICZALNE", "V3 / V2 ~ C * I_Y / (2*pi*d)"),
    ("Falsyfikacja (lambda<1 AU)", "DOSTEPNA", "Pioneer anomaly / MOND"),
    ("Falsyfikacja (progi)","DOSTEPNA", "Binarne pulsary, LIGO"),
]
for name, status, value in table:
    print(f"  {name:40s}  {status:20s}  {value}")

print()
print("WNIOSEK: Droga B daje dokladnie jednoparametrowa teorie (m_sp).")
print("         Dla m_sp wyznaczonego z obserwacji - wszystkie efekty TGP")
print("         sa OBLICZALNE Z DOKLADNOSCIA NUMERYCZNA (calka Feynmana).")
print()
print("         C_i = m_i/(2*sqrt(pi)) oznacza pelne rownowazenie masy")
print("         inercjalnej i ladunku TGP - ZASADA ROWNOWAZNOSCI z TGP!")
print()

# ===========================================================================
# Wykres: V2 i V3/V2 jako funkcja m_sp*r
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Lewy: G_eff(r) dla roznych m_sp
ax = axes[0]
r_AU_arr = np.logspace(-1, 4, 200)  # w AU
for m_sp_label, m_sp_val_Planck in [
    (r"$\lambda=10^4$ AU", l_Planck_m/(1e4*AU_m)),
    (r"$\lambda=10^3$ AU", l_Planck_m/(1e3*AU_m)),
    (r"$\lambda=100$ AU",  l_Planck_m/(100*AU_m)),
    (r"$\lambda=10$ AU",   l_Planck_m/(10*AU_m)),
    (r"$\lambda=1$ AU",    l_Planck_m/AU_m),
]:
    r_P = r_AU_arr * AU_m / l_Planck_m
    G_eff = np.exp(-m_sp_val_Planck * r_P)
    ax.semilogx(r_AU_arr, G_eff, lw=2, label=m_sp_label)

ax.axhline(1.0, color='k', ls='--', lw=1, alpha=0.5)
ax.axhline(0.9, color='r', ls=':', lw=1, alpha=0.3)
ax.set_xlabel("r [AU]", fontsize=12)
ax.set_ylabel(r"$G_{\rm eff}(r)/G$", fontsize=12)
ax.set_title(r"Efektywna stala Newtona $G_{\rm eff}(r) = G\,e^{-m_{sp}r}$", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.05)

# Prawy: moc promieniowania jako funkcja czestotliwosci
ax = axes[1]
omega_arr = np.linspace(0, 5, 300)
for m_sp_v, lbl in [(0.0, r"$m_{sp}=0$"), (1.0, r"$m_{sp}=1$"),
                     (2.0, r"$m_{sp}=2$"), (3.0, r"$m_{sp}=3$")]:
    mask = omega_arr > m_sp_v
    P_arr = np.zeros_like(omega_arr)
    P_arr[mask] = omega_arr[mask]**4 * np.sqrt(1.0 - (m_sp_v/omega_arr[mask])**2)
    # Normalizuj wzgledem P(omega=5, m=0)
    norm = 5.0**4
    ax.plot(omega_arr, P_arr / norm, lw=2, label=lbl)

ax.set_xlabel(r"$\omega / m_{sp,0}$", fontsize=12)
ax.set_ylabel(r"$P(\omega) / P_0(\omega=5)$  [znorm.]", fontsize=12)
ax.set_title(r"Moc promieniowania TGP: $P \propto \omega^4\sqrt{1-m_{sp}^2/\omega^2}$", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, None)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "ex12_path_b.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres: {out}")
