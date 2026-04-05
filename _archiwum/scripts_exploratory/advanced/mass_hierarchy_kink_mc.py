#!/usr/bin/env python3
"""
mass_hierarchy_kink_mc.py
==========================
Problem O16 TGP: Wspolczynnik amplifikacji hierarchii mas F(chi0)

Oblicza profile kinkowe phi_n(x) dla n=1,2,3,4 w teorii phi^4,
energie kinkowe, warunek braku 4. generacji oraz F(chi0) metoda MC.

Parametry: m_sp=1.0, phi0=1.0, beta=gamma=0.5

Uzycie:
    python mass_hierarchy_kink_mc.py
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

# ──────────────────────────────────────────────────────────────────────────────
# Parametry TGP
# ──────────────────────────────────────────────────────────────────────────────
M_SP   = 1.0   # masa pola przestrzennosci
PHI0   = 1.0   # wartoscia prozni phi_0
BETA   = 0.5   # parametr beta
GAMMA  = 0.5   # parametr gamma
N_MC   = 1000  # liczba prob Monte Carlo
X_MAX  = 30.0  # zasieg siatki przestrzennej
N_X    = 3000  # liczba punktow siatki

# ──────────────────────────────────────────────────────────────────────────────
# 1. Profile kinkowe phi_n(x)
# ──────────────────────────────────────────────────────────────────────────────

def kink_profile_analytic(x, x0=0.0, m=M_SP, phi0=PHI0):
    """
    Analityczny profil kinku phi^4 (n=1):
        phi(x) = phi0 * tanh(m * (x - x0) / sqrt(2))
    Warunki brzegowe: phi(-inf) = -phi0, phi(+inf) = +phi0.
    """
    return phi0 * np.tanh(m * (x - x0) / np.sqrt(2.0))


def multi_kink_profile(x, n, m=M_SP, phi0=PHI0, separation=8.0):
    """
    n-kink: zlozenie n alternujacych kinkow i antikinkow.
    Rozmieszczone rownomiernie z odlegloscia 'separation'.
    n parzyste: phi(-inf)=+phi0, phi(+inf)=+phi0 (kink-antikink-kink...)
    n nieparzyste: phi(-inf)=-phi0, phi(+inf)=+phi0

    Uzywamy przyblizonej superpozycji (dokladne dla dostatecznie odleglych kinkow).
    """
    # Pozycje kinkow
    positions = np.linspace(-(n - 1) * separation / 2.0,
                             (n - 1) * separation / 2.0, n)
    # Znaki: naprzemiennie +1, -1 dla kink / antikink
    signs = np.array([(-1)**k for k in range(n)], dtype=float)

    # Inicjalizacja
    phi = np.zeros_like(x, dtype=float)

    # Superpozycja przez produkt (lepsze dla wielokrotnych kinkow)
    # Alternatywnie: phi(x) = phi0 * prod_k tanh(m*(x-x_k)/sqrt(2))^{sgn_k}
    # Uzywamy addytywnej aproksymacji: phi = -phi0 + sum_k 2*phi0 * sigma(m*(x-xk)/sqrt(2))
    # gdzie sigma to funkcja schodkowa dla kinku.

    # Prostsza metoda: zlozenie przez produkt tanh
    result = np.ones_like(x, dtype=float)
    for k, (x0, sgn) in enumerate(zip(positions, signs)):
        arg = m * (x - x0) / np.sqrt(2.0)
        if sgn > 0:
            result = result * np.tanh(arg)
        else:
            result = result * (-np.tanh(arg))

    # Normalizacja do phi0
    phi = phi0 * result

    # Korekta znaku dla parzystej liczby kinkow
    if n % 2 == 0:
        phi = -phi

    return phi


def compute_kink_energy(phi_arr, x_arr, m=M_SP, phi0=PHI0):
    """
    Energia kinku:
        E = integral( (dphi/dx)^2 + V(phi) ) dx
        V(phi) = (m^2/4) * (phi^2 - phi0^2)^2
    """
    dx = x_arr[1] - x_arr[0]

    # Gradient numeryczny
    dphi_dx = np.gradient(phi_arr, dx)

    # Potencjal
    V = (m**2 / 4.0) * (phi_arr**2 - phi0**2)**2

    # Energia
    integrand = dphi_dx**2 + V
    E = np.trapezoid(integrand, x_arr)
    return E


def analytic_kink_energy(m=M_SP, phi0=PHI0):
    """
    Analityczna energia kinku phi^4:
        E_1 = (4/3) * m * phi0^2
    (znany wynik dla teorii phi^4 z potencjalem V = (m^2/4)(phi^2-phi0^2)^2)
    """
    return (4.0 / 3.0) * m * phi0**2


# ──────────────────────────────────────────────────────────────────────────────
# 2. Obliczenia profili i energii
# ──────────────────────────────────────────────────────────────────────────────

def compute_all_kink_energies():
    """Oblicza energie kinkow n=1,2,3,4 numerycznie i analitycznie."""
    x = np.linspace(-X_MAX, X_MAX, N_X)

    energies = {}
    profiles = {}

    print("=" * 60)
    print("PROFILE KINKOWE I ENERGIE (phi^4, m_sp=1, phi0=1)")
    print("=" * 60)

    E1_analytic = analytic_kink_energy()
    print(f"\nWeryfikacja analityczna:")
    print(f"  E_1 (analityczna) = (4/3)*m*phi0^2 = {E1_analytic:.6f}")

    for n in range(1, 5):
        if n == 1:
            phi = kink_profile_analytic(x)
        else:
            phi = multi_kink_profile(x, n)

        E = compute_kink_energy(phi, x)
        energies[n] = E
        profiles[n] = phi

        ratio_to_E1_analytic = E / E1_analytic if E1_analytic > 0 else 0.0

        print(f"\n  n={n} kink:")
        print(f"    E_{n} (numeryczna) = {E:.6f}")
        print(f"    E_{n}/E_1_analytic = {ratio_to_E1_analytic:.4f}")
        print(f"    E_{n}/E_1_numeryczna = {E/energies[1]:.4f}" if n > 1 else "")

    # Weryfikacja E1
    E1_num = energies[1]
    err_pct = abs(E1_num - E1_analytic) / E1_analytic * 100.0
    print(f"\nWeryfikacja E_1:")
    print(f"  E_1 numeryczna  = {E1_num:.6f}")
    print(f"  E_1 analityczna = {E1_analytic:.6f}")
    print(f"  Blad wzgledny   = {err_pct:.3f}%")

    return energies, profiles, x


# ──────────────────────────────────────────────────────────────────────────────
# 3. Warunek braku 4. generacji
# ──────────────────────────────────────────────────────────────────────────────

def check_generation_threshold(energies, m=M_SP, phi0=PHI0):
    """
    Sprawdza warunek braku 4. generacji:
        Prog: E_prog = m_sp * phi0 * sqrt(2) (bariera miedzy prozniami)
        Warunek: E_3 < E_prog < E_4 => 3 generacje stabilne, 4 niestabilna
    """
    E_threshold = m * phi0 * np.sqrt(2.0)

    print("\n" + "=" * 60)
    print("WARUNEK BRAKU 4. GENERACJI")
    print("=" * 60)
    print(f"\nProg energetyczny:")
    print(f"  E_prog = m_sp * phi0 * sqrt(2) = {E_threshold:.6f}")
    print(f"\nEnergie kinkow:")
    for n in range(1, 5):
        E = energies[n]
        status = "< E_prog (stabilna)" if E < E_threshold else "> E_prog (niestabilna)"
        print(f"  E_{n} = {E:.6f}  {status}")

    E3 = energies[3]
    E4 = energies[4]

    print(f"\nStosunek E_3/E_prog = {E3/E_threshold:.4f}")
    print(f"Stosunek E_4/E_prog = {E4/E_threshold:.4f}")

    if E3 < E_threshold < E4:
        print("\nWYNIK: E_3 < E_prog < E_4 => 3 generacje stabilne, 4. niestabilna [OK]")
        three_gen_condition = True
    elif E3 < E_threshold and E4 < E_threshold:
        print("\nWYNIK: E_4 < E_prog => nawet 4. generacja stabilna [NIEZGODNE z obserwacjami]")
        three_gen_condition = False
    else:
        print("\nWYNIK: E_3 > E_prog => nawet 3. generacja niestabilna [PROBLEM]")
        three_gen_condition = False

    return E_threshold, three_gen_condition


# ──────────────────────────────────────────────────────────────────────────────
# 4. Wspolczynnik F(chi0) metoda Monte Carlo
# ──────────────────────────────────────────────────────────────────────────────

def compute_kink_energy_with_background(chi0, x_arr, m=M_SP, phi0=PHI0):
    """
    Energia kinku przy tle chi0: phi_background = phi0 + chi0.
    Kink przesuwa sie miedzy -(phi0+chi0) a +(phi0+chi0).
    Efektywna phi0_eff = phi0 + chi0.

    F(chi0) = E_kink(chi0) / E_kink(0)
    """
    phi0_eff = phi0 + chi0
    # Profil kinku w zmodyfikowanej prozni
    phi = kink_profile_analytic(x_arr, phi0=phi0_eff, m=m)
    E = compute_kink_energy(phi, x_arr, m=m, phi0=phi0_eff)
    return E


def monte_carlo_F(n_samples=N_MC, x_arr=None, m=M_SP, phi0=PHI0):
    """
    Oblicza F(chi0) metoda Monte Carlo:
        chi0 ~ Uniform(0.01, 0.5) * phi0
        F(chi0) = E_kink(chi0) / E_kink(0)
        Zwraca srednia i odchylenie standardowe.
    """
    if x_arr is None:
        x_arr = np.linspace(-X_MAX, X_MAX, N_X)

    # Energia referencyjna (czysta proznia, chi0=0)
    E0 = compute_kink_energy_with_background(0.0, x_arr, m=m, phi0=phi0)

    # Probkowanie MC
    rng = np.random.default_rng(42)
    chi0_samples = rng.uniform(0.01 * phi0, 0.50 * phi0, n_samples)

    F_values = np.zeros(n_samples)
    for i, chi0 in enumerate(chi0_samples):
        E_chi = compute_kink_energy_with_background(chi0, x_arr, m=m, phi0=phi0)
        F_values[i] = E_chi / E0 if E0 > 0 else 1.0

    F_mean = np.mean(F_values)
    F_std  = np.std(F_values)
    F_min  = np.min(F_values)
    F_max  = np.max(F_values)

    print("\n" + "=" * 60)
    print(f"WSPOLCZYNNIK F(chi0) --- MONTE CARLO (N={n_samples})")
    print("=" * 60)
    print(f"  chi0 ~ Uniform({0.01*phi0:.3f}, {0.50*phi0:.3f}) * phi0")
    print(f"  E_kink(chi0=0) = {E0:.6f}")
    print(f"  <F(chi0)>      = {F_mean:.6f}")
    print(f"  sigma_F        = {F_std:.6f}")
    print(f"  F_min          = {F_min:.6f}")
    print(f"  F_max          = {F_max:.6f}")
    print(f"  95% CI         = [{F_mean - 2*F_std:.6f}, {F_mean + 2*F_std:.6f}]")

    return F_mean, F_std, F_values


# ──────────────────────────────────────────────────────────────────────────────
# 5. Stosunek mas generacji
# ──────────────────────────────────────────────────────────────────────────────

def print_mass_ratios(energies, E_threshold):
    """
    Drukuje stosunki mas generacji i porownuje z obserwacjami.

    W prostym modelu n-kink phi^4: E_n ~ n * E_1
    To daje m2/m1 ~ 2, m3/m1 ~ 3 (za male wzgledem obserwacji).

    Obserwacje:
        mu/e  ~ 207     (mion/elektron)
        tau/e ~ 3477    (tau/elektron)

    Topologiczna interpretacja: stosunek energii kinkow ~ stosunek mas.
    """
    E1 = energies[1]
    E2 = energies[2]
    E3 = energies[3]
    E4 = energies[4]

    print("\n" + "=" * 60)
    print("STOSUNKI MAS GENERACJI (model kinkowy phi^4)")
    print("=" * 60)
    print("\nProsty model n-kink (E_n ~ n*E_1):")
    print(f"  m2/m1 = E_1/E_0 = E_2/E_1 = {E2/E1:.4f}")
    print(f"         (teoria TGP vs exp: mu/e ~ 207)")
    print(f"  m3/m1 = E_2/E_0 = E_3/E_1 = {E3/E1:.4f}")
    print(f"         (teoria TGP vs exp: tau/e ~ 3477)")
    print(f"  E_3/E_prog = {E3/E_threshold:.6f}  (< 1 => 3. generacja stabilna)")
    print(f"  E_4/E_prog = {E4/E_threshold:.6f}  (> 1 => 4. generacja niestabilna)")

    print("\nUwaga: W modelu phi^4 E_n ~ n*E_1 (hierarchia liniowa).")
    print("Obserwowana hierarchia mas (mu/e~207, tau/e~3477) wymaga")
    print("dodatkowego mechanizmu wzmocnienia F(chi0) lub nieliniowego")
    print("sprzezenia profili kinkowych z tlem substratowym.")
    print()
    print("Alternatywne podejscie (nakладka profili):")
    # Ratio wyznaczony z profilu nakладки - wzmocnienie eksponencjalne
    # przy zakladaniu ze energia kinku skaluje sie z phi0^2 * m
    # i ze chi0 modyfikuje efektywne phi0
    chi0_mu  = PHI0 * (np.log(207.0) / np.log(E2/E1 + 1e-10) - 1.0) \
               if E2/E1 > 1 else 0.0
    # Uproszczone: ratio_i = exp(kappa * n_i) gdzie kappa z dopasowania
    kappa_fit = np.log(207.0) / 1.0  # z mu/e dla n=1->2
    print(f"  Jesli F(chi0) ~ exp(kappa*n): kappa ~ ln(207) = {np.log(207.0):.4f}")
    print(f"  Daje: m2/m1 = exp(kappa) = {np.exp(np.log(207.0)):.1f} (= 207 z definicji)")
    print(f"  Daje: m3/m1 = exp(2*kappa) = {np.exp(2*np.log(207.0)):.1f} (teoria)")
    print(f"         (exp: tau/e ~ 3477; ln^2: {np.exp(2*np.log(207.0)):.1f})")
    print(f"  Stosunek exp/teoria (m3/m1): {3477.0/np.exp(2*np.log(207.0)):.4f}")


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "=" * 60)
    print("TGP - PROBLEM O16: HIERARCHIA MAS KINKOWA (MC)")
    print(f"Parametry: m_sp={M_SP}, phi0={PHI0}, beta={BETA}, gamma={GAMMA}")
    print("=" * 60)

    # 1. Profile i energie kinkow
    energies, profiles, x_arr = compute_all_kink_energies()

    # 2. Warunek braku 4. generacji
    E_threshold, three_gen_ok = check_generation_threshold(energies)

    # 3. Stosunki mas
    print_mass_ratios(energies, E_threshold)

    # 4. Monte Carlo F(chi0)
    F_mean, F_std, F_values = monte_carlo_F(n_samples=N_MC, x_arr=x_arr)

    # ── Podsumowanie wynikow ──────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("PODSUMOWANIE WYNIKOW O16")
    print("=" * 60)
    E1_analytic = analytic_kink_energy()
    E1_num = energies[1]
    print(f"  E_1 (analityczna) = {E1_analytic:.6f}  [(4/3)*m*phi0^2]")
    print(f"  E_1 (numeryczna)  = {E1_num:.6f}  [blad: {abs(E1_num-E1_analytic)/E1_analytic*100:.3f}%]")
    print(f"  E_2               = {energies[2]:.6f}")
    print(f"  E_3               = {energies[3]:.6f}")
    print(f"  E_4               = {energies[4]:.6f}")
    print(f"  E_prog            = {E_threshold:.6f}")
    print(f"  m2/m1 = E_2/E_1   = {energies[2]/energies[1]:.4f}  (exp: mu/e ~ 207)")
    print(f"  m3/m1 = E_3/E_1   = {energies[3]/energies[1]:.4f}  (exp: tau/e ~ 3477)")
    print(f"  E_3/E_prog        = {energies[3]/E_threshold:.6f}  (< 1: 3 generacje)")
    print(f"  E_4/E_prog        = {energies[4]/E_threshold:.6f}  (> 1: 4. niegeneracja)")
    print(f"  3 generacje OK?   = {three_gen_ok}")
    print(f"  <F(chi0)>         = {F_mean:.6f} +/- {F_std:.6f}")
    print(f"  [wyniki dla LaTeX]")
    print(f"  E_1={E1_num:.4f}, E_2={energies[2]:.4f}, E_3={energies[3]:.4f}, E_4={energies[4]:.4f}")
    print(f"  m2/m1={energies[2]/energies[1]:.4f}, m3/m1={energies[3]/energies[1]:.4f}")
    print(f"  E_3/E_prog={energies[3]/E_threshold:.4f}, F={F_mean:.4f}+/-{F_std:.4f}")

    # Wykresy opcjonalne
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("TGP O16: Profile kinkowe i hierarchia mas", fontsize=14)

        colors = ['b', 'g', 'r', 'm']
        ax = axes[0, 0]
        for n in range(1, 5):
            ax.plot(x_arr, profiles[n], color=colors[n-1], label=f'n={n}')
        ax.axhline(PHI0, color='k', ls='--', alpha=0.5, label=f'phi0={PHI0}')
        ax.axhline(-PHI0, color='k', ls='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('phi(x)')
        ax.set_title('Profile kinkowe phi_n(x)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-20, 20)

        ax = axes[0, 1]
        ns = list(range(1, 5))
        Es = [energies[n] for n in ns]
        ax.bar(ns, Es, color=colors, alpha=0.7)
        ax.axhline(E_threshold, color='k', ls='--', label=f'E_prog={E_threshold:.3f}')
        ax.set_xlabel('n (numer kinku)')
        ax.set_ylabel('E_n')
        ax.set_title('Energie kinkow vs prog')
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        chi0_range = np.linspace(0.01 * PHI0, 0.50 * PHI0, 50)
        E0_ref = compute_kink_energy_with_background(0.0, x_arr)
        F_curve = [compute_kink_energy_with_background(c, x_arr) / E0_ref
                   for c in chi0_range]
        ax.plot(chi0_range / PHI0, F_curve, 'b-')
        ax.axhline(F_mean, color='r', ls='--', label=f'<F>={F_mean:.4f}')
        ax.set_xlabel('chi0 / phi0')
        ax.set_ylabel('F(chi0)')
        ax.set_title('Wspolczynnik amplifikacji F(chi0)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[1, 1]
        ax.hist(F_values, bins=40, color='steelblue', alpha=0.7, edgecolor='k')
        ax.axvline(F_mean, color='r', ls='--', label=f'<F>={F_mean:.4f}')
        ax.axvline(F_mean + F_std, color='r', ls=':', alpha=0.7)
        ax.axvline(F_mean - F_std, color='r', ls=':', alpha=0.7)
        ax.set_xlabel('F(chi0)')
        ax.set_ylabel('Liczba prob')
        ax.set_title(f'Rozklad MC F(chi0), N={N_MC}')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = "C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/plots/O16_kink_mass_hierarchy.png"
        plt.savefig(out_path, dpi=100)
        print(f"\nWykres zapisany: {out_path}")
        plt.close()
    except Exception as e:
        print(f"\n[INFO] Wykresy niedostepne: {e}")

    return {
        'energies': energies,
        'E_threshold': E_threshold,
        'three_gen_ok': three_gen_ok,
        'F_mean': F_mean,
        'F_std': F_std,
    }


if __name__ == '__main__':
    results = main()
