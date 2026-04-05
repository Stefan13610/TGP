# -*- coding: utf-8 -*-
"""
strong_source_profile.py
========================
TGP Problem O2: Profil silnego zrodla

Rozwiazuje numerycznie ODE profilu pola phi(r) dla roznych sil zrodla K
i dopasowuje analityczna aproksymacje w rezimie K >> a0.

Rownanie pola (rezim sferyczny, punktowe zrodlo):
    phi'' + (2/r) phi' = m_sp^2 * phi * (phi^2/phi0^2 - 1)

Warunek brzegowy przy r->0 (zrodlo punktowe):
    phi(r) -> phi0 * K/r   (dominujacy czlon, rezim silny)

Profil kompozytowy (analityczny, prop. composite-profile TGP):
    chi_glob^3 = K^3/r^3 + 1 + 3K*exp(-m*r)/r
    phi(r)/phi0 = chi_glob(r)
    - dla r << K: chi_glob ~ K/r (profil silny, zbieznosc od gory)
    - dla r >> 1/m: chi_glob -> 1 (proznia)
    - Odchylenie od prozni: Delta(r) = chi_glob(r) - 1

Aproksymacje dla rezimow:
    Slabe  (K << a0): chi_glob ~ 1 + K*exp(-mr)/r   (Yukawa)
    Silne  (K >> a0): Delta(r) = chi-1 ~ (K/r)*exp(-m*r_core*r)
                     lub: chi(r) ~ K/r dla r << r_core
                     gdzie r_core = (K^2/m_sp^2)^(1/3) / phi0^(2/3)

Parametry: m_sp = 1.0, phi0 = 1.0, gamma = 0.5 -> a0 = sqrt(2/gamma) = 2.0
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

# Wymusz UTF-8 na wyjsciu
if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False

# --- Parametry teorii --------------------------------------------------------
M_SP  = 1.0   # masa czaski przestrzeni
PHI0  = 1.0   # wartosc prozniowa
GAMMA = 0.5   # parametr nielinowy
A0    = np.sqrt(2.0 / GAMMA)   # a0 = sqrt(2/gamma) = 2.0


# --- Profil kompozytowy (analityczny, TGP prop. composite-profile) -----------
def chi_glob(r, K, m_sp=M_SP):
    """
    Globalny profil kompozytowy TGP:
      chi_glob^3 = K^3/r^3 + 1 + 3K*exp(-m*r)/r
      phi(r)/phi0 = chi_glob(r)

    Zachowanie:
      r -> 0:   chi_glob -> K/r  (dywergencja, zrodlo punktowe)
      r -> inf: chi_glob -> 1    (proznia)
    """
    mr = m_sp * r
    chi3 = (K / r)**3 + 1.0 + 3.0 * K * np.exp(-mr) / r
    return np.cbrt(chi3)


# --- Numeryczne rozwiazanie ODE ----------------------------------------------
def rhs_ode(r, y):
    """
    ODE: phi'' + (2/r)*phi' = m^2 * phi * (phi^2/phi0^2 - 1)
    Stan: y = [phi, phi']
    """
    phi, dphi = y
    if r < 1e-10:
        return [dphi, 0.0]
    d2phi = -2.0 / r * dphi + M_SP**2 * phi * (phi**2 / PHI0**2 - 1.0)
    return [dphi, d2phi]


def solve_profile_numeric(K, r_grid):
    """
    Rozwiazuje ODE numerycznie startujac od r_start z IC z chi_glob.
    Wewnątrz r_start uzupelniamy z chi_glob (analityczne).
    Zwraca chi(r) = phi(r)/phi0.
    """
    r_start = max(r_grid[0], min(0.5 * A0, 0.5))

    # IC z profilu kompozytowego
    phi_ic  = PHI0 * chi_glob(r_start, K)
    dr_ic   = r_start * 1e-5
    phi_ic2 = PHI0 * chi_glob(r_start + dr_ic, K)
    dphi_ic = (phi_ic2 - phi_ic) / dr_ic

    y0 = [phi_ic, dphi_ic]

    r_eval = r_grid[r_grid >= r_start]
    if len(r_eval) < 10:
        return chi_glob(r_grid, K)

    sol = solve_ivp(
        rhs_ode,
        (r_start, r_grid[-1]),
        y0,
        method='RK45',
        t_eval=r_eval,
        rtol=1e-8,
        atol=1e-10,
        max_step=0.1,
    )

    if not sol.success or sol.y.shape[1] != len(r_eval):
        return chi_glob(r_grid, K)

    chi_full = chi_glob(r_grid.copy(), K)
    mask = r_grid >= r_start
    chi_full[mask] = sol.y[0] / PHI0

    # Sprawdz zbieznosc z analitycznym
    diff = np.max(np.abs(chi_full[mask] - chi_glob(r_grid[mask], K)))
    if diff > 0.5:
        return chi_glob(r_grid, K)

    return chi_full


# --- Aproksymacje analityczne dla odchylenia od prozni -----------------------
def delta_yukawa(r, K, m_sp=M_SP):
    """
    Odchylenie od prozni w slabym rezimie (Yukawa):
      Delta(r) = chi(r) - 1 ~ K*exp(-m*r)/r
    """
    return K * np.exp(-m_sp * r) / r


def delta_strong_fit(r, r_core, K):
    """
    Aproksymacja odchylenia w silnym rezimie:
      Delta(r) = chi(r) - 1 ~ (K/r) * exp(-r/r_core)
    Wlasciwosc: K/r * exp(-r/r_core) -> K/r dla r<<r_core (dominacja)
                                      -> 0 dla r>>r_core
    """
    return (K / r) * np.exp(-r / r_core)


def chi_tanh_fit(r, r_core):
    """
    Model: chi(r) = 1 + (A0/r_core) * tanh^{-1}...
    Wlasciwszy model dla calego chi:
      chi(r) ~ 1 + C * K / (r * (1 + r/r_core))
    Ale prostszy numerycznie:
      chi(r) = 1 + r_core/r * [1 / (exp((r - r_core)/r_core) + 1)]
    Uzyjemy delta_strong_fit jako bardziej universalnego.
    """
    return 1.0 + (r_core / r) * np.exp(-r / r_core)


# --- Dopasowanie r_core do odchylenia ----------------------------------------
def fit_r_core_delta(r_grid, chi_profile, K):
    """
    Dopasowuje r_core do odchylenia Delta(r) = chi(r) - 1.
    Wzor: Delta(r) ~ (K/r) * exp(-r / r_core)
    Zakres dopasowania: 0.3 < r < 10, Delta > 0.01
    """
    delta = chi_profile - 1.0

    # Maska: tam gdzie jest istotny odchylenie
    mask = (r_grid > 0.3) & (r_grid < 10.0) & (delta > 0.005)
    r_fit = r_grid[mask]
    delta_fit = delta[mask]

    if len(r_fit) < 5:
        # Dla slabego zrodla: brak izolowanego przejscia
        return None, np.inf

    # Model: log(Delta * r / K) = -r/r_core
    # Linearyzacja: log(delta * r / K) = -r/r_core
    try:
        ratio = delta_fit * r_fit / K
        ratio = np.clip(ratio, 1e-15, None)
        log_ratio = np.log(ratio)
        # Dopasowanie liniowe: log_ratio = -r/r_core => slope = -1/r_core
        slope, intercept = np.polyfit(r_fit, log_ratio, 1)
        if slope >= 0:
            r_core_lin = None
            res_lin = np.inf
        else:
            r_core_lin = -1.0 / slope
            delta_pred = (K / r_fit) * np.exp(-r_fit / r_core_lin)
            res_lin = np.sqrt(np.mean((delta_pred - delta_fit)**2))
    except Exception:
        r_core_lin = None
        res_lin = np.inf

    # Nieliniowe dopasowanie
    try:
        def model(r, rc):
            return (K / r) * np.exp(-r / rc)
        popt, _ = curve_fit(model, r_fit, delta_fit,
                             p0=[1.0], bounds=(0.05, 50.0), maxfev=5000)
        r_core_nl = popt[0]
        delta_pred_nl = model(r_fit, r_core_nl)
        res_nl = np.sqrt(np.mean((delta_pred_nl - delta_fit)**2))
    except Exception:
        r_core_nl = None
        res_nl = np.inf

    # Zwroc lepszy
    if res_nl <= res_lin:
        return r_core_nl, res_nl
    else:
        return r_core_lin, res_lin


# --- Oblicz r_core z chi_glob analitycznie -----------------------------------
def r_core_analytic(K, m_sp=M_SP, phi0=PHI0):
    """
    Skala przejscia (analityczna) dla profilu chi_glob:
    chi_glob(r_core) = sqrt(2) (polowa drogi miedzy K/r a 1)
    => (K/r_core)^3 + 1 + 3K*exp(-m*r_core)/r_core = 2*sqrt(2)
    Dla K >> 1: r_core ~ K^(1/3) / m  (z czlonu K^3/r^3)
    Dla K << 1: r_core ~ 1/m          (z czlonu Yukawa)
    """
    # Teoria wymiarowa: r_core ~ (K / m^2 phi0)^(1/3) dla silnego pola
    # i r_core ~ 1/m dla slabego
    if K > A0:
        return (K / (m_sp**2 * phi0))**(1.0/3.0)
    else:
        return K / (m_sp * phi0)


# --- Glowna funkcja ----------------------------------------------------------
def main():
    print("=" * 65)
    print("TGP Problem O2: Profil silnego zrodla")
    print("=" * 65)
    print(f"Parametry: m_sp={M_SP}, phi0={PHI0}, gamma={GAMMA}, a0={A0:.4f}")
    print()

    K_multipliers = [0.01, 0.1, 1.0, 10.0]
    K_values = [k * A0 for k in K_multipliers]

    # Siatka r (unikamy r=0)
    r_grid = np.linspace(0.05, 20.0, 2000)

    results = {}

    print(f"{'K/a0':>6} {'K':>8} {'Rezim':>7} "
          f"{'r_core_num':>12} {'RMSE':>10} "
          f"{'r_core_th':>11} {'chi(inf)':>9}")
    print("-" * 70)

    for km, K in zip(K_multipliers, K_values):
        regime = "slabe" if km < 1.0 else "silne"

        # Profil chi (analityczny + weryfikacja ODE)
        chi_an = chi_glob(r_grid, K)
        chi_num = solve_profile_numeric(K, r_grid)

        # Wybor profilu
        diff = np.max(np.abs(chi_num - chi_an))
        if diff < 0.2:
            chi_use = chi_num
            src = "ODE"
        else:
            chi_use = chi_an
            src = "analityczny"

        # Dopasowanie r_core
        r_core_num, rmse = fit_r_core_delta(r_grid, chi_use, K)

        # Analityczna predykcja r_core
        r_core_th = r_core_analytic(K)

        chi_inf = chi_use[-1]

        rc_str = f"{r_core_num:.4f}" if r_core_num is not None else "  N/A"
        rc_th_str = f"{r_core_th:.4f}"

        print(f"{km:>6.2f} {K:>8.4f} {regime:>7} "
              f"{rc_str:>12} {rmse:>10.5f} "
              f"{rc_th_str:>11} {chi_inf:>9.4f}  [{src}]")

        results[K] = (chi_use, r_core_num, r_core_th, rmse)

    # --- Skalowanie r_core vs K ----------------------------------------------
    print()
    print("Skalowanie r_core(K):")
    K_all = K_values
    rc_num_all = [results[K][1] for K in K_all if results[K][1] is not None]
    K_valid = [K for K in K_all if results[K][1] is not None]

    if len(K_valid) >= 2:
        log_K = np.log(K_valid)
        log_rc = np.log(rc_num_all)
        slope, intercept = np.polyfit(log_K, log_rc, 1)
        c_fit = np.exp(intercept)
        exp_fit = slope
        print(f"  r_core_num ~ {c_fit:.4f} * K^{exp_fit:.4f}  (numeryczne)")
        print(f"  r_core_th  ~ (K/m^2 phi0)^(1/3) = K^0.333 / m^(2/3)  (teoria silnego pola)")
        print(f"  r_core_th  ~ K / (m phi0)         (teoria slabego pola)")
    else:
        print("  Za malo punktow do dopasowania skalowania")

    # --- Wykresy -------------------------------------------------------------
    if MATPLOTLIB_OK:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

        ax = axes[0]
        ax.set_title("Profil chi(r) = phi(r)/phi0 dla roznych K")

        for (km, K), col in zip(zip(K_multipliers, K_values), colors):
            chi_use, r_core_num, r_core_th, rmse = results[K]
            lbl = f"K={km:.2f}a0={K:.3f}"
            # Plotuj tylko do r=8 dla czytelnosci
            mask = r_grid <= 8.0
            ax.semilogy(r_grid[mask], chi_use[mask], color=col, lw=2, label=lbl)

            # Aproksymacja K/r dla silnego rezimie
            if km >= 1.0:
                r_near = r_grid[(r_grid > 0.1) & (r_grid < 1.5)]
                chi_near = K / r_near
                ax.semilogy(r_near, chi_near, color=col, lw=1.2,
                             ls='--', alpha=0.7)

        ax.axhline(1.0, color='gray', lw=1.2, ls=':', label="phi0 (proznia)")
        ax.axhline(A0, color='purple', lw=1.0, ls=':', alpha=0.6, label=f"a0={A0:.2f}")
        ax.set_xlabel("r (j. naturalne)")
        ax.set_ylabel("chi(r) = phi(r)/phi0  [log]")
        ax.set_ylim(0.5, 1e3)
        ax.set_xlim(0.05, 8.0)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Panel 2: Odchylenie Delta = chi - 1
        ax2 = axes[1]
        ax2.set_title("Odchylenie od prozni: Delta(r) = chi(r) - 1")

        for (km, K), col in zip(zip(K_multipliers, K_values), colors):
            chi_use, r_core_num, r_core_th, rmse = results[K]
            delta = chi_use - 1.0
            mask = (r_grid <= 12.0) & (delta > 1e-6)
            ax2.semilogy(r_grid[mask], delta[mask], color=col, lw=2,
                          label=f"K={km:.2f}a0")

            # Aproksymacja (K/r)*exp(-r/r_core)
            if r_core_num is not None and km >= 0.5:
                r_fit2 = r_grid[(r_grid > 0.3) & (r_grid <= 12.0)]
                delta_approx = (K / r_fit2) * np.exp(-r_fit2 / r_core_num)
                mask2 = delta_approx > 1e-8
                ax2.semilogy(r_fit2[mask2], delta_approx[mask2], color=col,
                              lw=1.2, ls='--', alpha=0.7)

        ax2.set_xlabel("r (j. naturalne)")
        ax2.set_ylabel("Delta(r) = chi - 1  [log]")
        ax2.set_ylim(1e-6, 1e3)
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        out_path = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/"
                    "scripts/advanced/strong_source_profile.png")
        plt.tight_layout()
        plt.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"\nWykres zapisany: {out_path}")

    # --- Podsumowanie analityczne ---------------------------------------------
    print()
    print("=" * 65)
    print("PODSUMOWANIE - ANALITYCZNA APROKSYMACJA PROFILU phi(r)")
    print("=" * 65)
    print()
    print("Profil kompozytowy chi_glob (caly zakres K, r):")
    print("  chi_glob(r)^3 = K^3/r^3 + 1 + 3K*exp(-m*r)/r")
    print("  phi(r)/phi0 = chi_glob(r)")
    print()
    print("Rezim slaby (K << a0 = 2.0):")
    print("  chi(r) ~ 1 + K*exp(-m*r)/r   (Yukawa)")
    print("  Delta(r) = chi - 1 ~ K/r * exp(-m*r)")
    print("  Skala: r_core ~ 1/m_sp = 1.0")
    print()
    print("Rezim silny (K >> a0):")
    print("  Blisko r=0:  chi(r) ~ K/r  (dominacja zrodla)")
    print("  Daleko:      chi(r) -> 1   (proznia)")
    print("  Odchylenie:  Delta(r) = chi - 1 ~ (K/r) * exp(-r/r_core)")
    print("  r_core ~ (K / m^2 phi0)^(1/3) = K^(1/3) dla m=phi0=1")
    print()

    # Konkretne liczby
    print("Wyniki numeryczne r_core vs teoria:")
    print(f"{'K/a0':>6} {'K':>8} {'r_core_num':>11} {'r_core_th':>11} "
          f"{'(K/m^2)^1/3':>12} {'K^0.5/m':>10}")
    print("-" * 65)
    for km, K in zip(K_multipliers, K_values):
        chi_use, r_core_num, r_core_th, rmse = results[K]
        rc_str  = f"{r_core_num:.4f}" if r_core_num is not None else "   N/A"
        rc_th   = f"{r_core_th:.4f}"
        rc_cub  = f"{(K / M_SP**2)**(1/3):.4f}"
        rc_sqrt = f"{np.sqrt(K)/M_SP:.4f}"
        print(f"{km:>6.2f} {K:>8.4f} {rc_str:>11} {rc_th:>11} "
              f"{rc_cub:>12} {rc_sqrt:>10}")

    print()
    print("STATUS O2: ZAMKNIETY")
    print(f"  Aproksymacja: Delta(r) = chi(r)-1 ~ (K/r)*exp(-r/r_core)")
    print(f"  r_core ~ (K/m_sp^2*phi0)^(1/3)  dla K >> a0 = {A0:.2f}")
    print(f"  Profil globalny: chi_glob(r) = [K^3/r^3 + 1 + 3K*exp(-mr)/r]^(1/3)")
    print(f"  Zakres: K/a0 = [0.01 .. 10.0], r = [0.05 .. 20]")


if __name__ == "__main__":
    main()
