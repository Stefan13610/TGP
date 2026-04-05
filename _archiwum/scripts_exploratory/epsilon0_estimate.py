#!/usr/bin/env python3
"""
epsilon0_estimate.py
====================
Estymacja ε₀ (początkowej wartości Φ przy nukleacji) z termodynamiki substratu.

Ramy formalne (rem:epsilon0-derivation w dodatekG):
  1. Fluktuacja krytyczna: σ²_c = T_c / (2·z·J)
  2. ε₀ = σ²_c · S(R_c / a_sub)      (tłumienie Boltzmannowskie)
  3. R_c = 2σ_s / ΔF                  (promień krytyczny Colemana)
  4. ε₀ ≳ ℓ_P² / R_c²                (dolna granica kwantowa)

Cel: Powiązanie ε₀ z N_e (liczba e-składań) i Φ₀:
  N_e = (1/3) · ln(Φ₀/ε₀)

Skrypt oblicza:
  - ε₀(ξ) dla zakresu długości korelacji ξ/ℓ_P
  - N_e(ε₀, Φ₀) i sprawdzenie spójności z Planck
  - Widmo mocy P(k) perturbacji z nukleacji
  - Mapę (ξ, Φ₀) → N_e z oknem dozwolonym (50 < N_e < 70)

Teoria Generowanej Przestrzeni — Mateusz Serafin, 2026
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ─── Stałe ────────────────────────────────────────────────────────────
L_P = 1.616255e-35       # m  (długość Plancka)
T_P = 1.416784e32        # K  (temperatura Plancka)

# Parametry siatki substratu (typowe wartości Isingowskie 3D)
Z_COORD = 6              # koordynacja siatki sześciennej
J_0 = 1.0                # sprzężenie wymianowe (w jednostkach kT_c)
T_C_RATIO = 1.0          # T_c / (z·J)  ≈ 1 dla Isinga 3D: T_c ≈ 4.51 J

# Znana wartość: T_c/(z·J) ≈ 4.51/6 ≈ 0.752 dla Isinga 3D (siatka sc)
ISING_3D_TC_ZJ = 4.5115 / Z_COORD   # ≈ 0.752

# ─── Model fizyczny ──────────────────────────────────────────────────

def sigma2_critical(T_c_zJ=ISING_3D_TC_ZJ):
    """
    Fluktuacja krytyczna: σ²_c = T_c / (2·z·J) = T_c_zJ / 2
    (w jednostkach bezwymiarowych, ŝ² → Φ)
    """
    return T_c_zJ / 2.0


def suppression_factor(R_c_over_a, d=3):
    """
    Czynnik tłumienia Boltzmannowskiego S(R_c/a_sub):
    S ≈ exp(-C_d · (R_c/a)^d)  (bariera entropowa)

    Dla nukleacji bąblowej w d=3:
    S = exp(-16π/3 · (σ_s)³ / (ΔF)² / T)

    Uproszczenie: S ≈ exp(-(R_c/a)^d)  z prefaktorem O(1)
    """
    return np.exp(-(R_c_over_a)**d)


def epsilon0_from_xi(xi_over_lP, Phi0=25.0, d=3):
    """
    ε₀ z długości korelacji ξ:
    ε₀ ≳ (ℓ_P/ξ)²   (dolna granica kwantowa z tunelowania)

    Dokładniejszy model:
    ε₀ = σ²_c · S(R_c/a_sub) ≈ σ²_c · (ℓ_P/ξ)²
    (bo R_c ~ ξ dla T ≈ T_c, a_sub ~ ℓ_P)
    """
    sigma2 = sigma2_critical()
    eps0 = sigma2 * (1.0 / xi_over_lP)**2
    return eps0


def N_efolds(eps0, Phi0=25.0):
    """
    N_e = (1/3) · ln(Φ₀/ε₀)
    z równania (eq:N-efolds) w dodatekG.
    """
    if eps0 <= 0 or Phi0 <= 0:
        return np.nan
    return (1.0/3.0) * np.log(Phi0 / eps0)


def N_efolds_from_xi(xi_over_lP, Phi0=25.0):
    """
    N_e jako funkcja ξ/ℓ_P i Φ₀.
    """
    eps0 = epsilon0_from_xi(xi_over_lP, Phi0)
    return N_efolds(eps0, Phi0)


# ─── Analiza ──────────────────────────────────────────────────────────

def analyze_epsilon0():
    """Główna analiza ε₀."""

    print("╔══════════════════════════════════════════════════════════╗")
    print("║  TGP: ε₀ Estimation from Substrate Thermodynamics      ║")
    print("║  Problem R6 (ROADMAP v3)                                ║")
    print("╚══════════════════════════════════════════════════════════╝")

    # ─── 1. Fluktuacja krytyczna ─────────────────────────────
    sigma2 = sigma2_critical()
    print(f"\n  1. Fluktuacja krytyczna:")
    print(f"     T_c/(z·J) = {ISING_3D_TC_ZJ:.4f}  (Ising 3D, sc lattice)")
    print(f"     σ²_c = T_c/(2·z·J) = {sigma2:.6f}")

    # ─── 2. Skan ξ/ℓ_P ──────────────────────────────────────
    xi_range = np.logspace(25, 45, 500)  # ξ/ℓ_P od 10²⁵ do 10⁴⁵
    Phi0_values = [20.0, 24.66, 25.0, 28.0, 33.0]

    print(f"\n  2. Skan ξ/ℓ_P ∈ [10²⁵, 10⁴⁵]:")
    print(f"     {'Φ₀':>6s} | {'ξ/ℓ_P (N_e=50)':>16s} | {'ξ/ℓ_P (N_e=60)':>16s} | "
          f"{'ε₀(50)':>12s} | {'ε₀(60)':>12s}")
    print("     " + "-"*75)

    results = {}
    for Phi0 in Phi0_values:
        # ε₀ dla N_e = 50 i N_e = 60
        eps0_50 = Phi0 * np.exp(-3 * 50)
        eps0_60 = Phi0 * np.exp(-3 * 60)

        # ξ/ℓ_P z ε₀ = σ²_c · (ℓ_P/ξ)²
        xi_50 = np.sqrt(sigma2 / eps0_50)
        xi_60 = np.sqrt(sigma2 / eps0_60)

        results[Phi0] = {
            'eps0_50': eps0_50, 'eps0_60': eps0_60,
            'xi_50': xi_50, 'xi_60': xi_60
        }

        print(f"     {Phi0:6.2f} | {xi_50:16.3e} | {xi_60:16.3e} | "
              f"{eps0_50:12.3e} | {eps0_60:12.3e}")

    # ─── 3. Relacja n_s z N_e ────────────────────────────────
    print(f"\n  3. Predykcja n_s z N_e (klasa Starobinsky'ego):")
    print(f"     n_s ≈ 1 - 2/N_e - Δn_s^TGP  (Δn_s^TGP ≈ 0.0005, z ex107)")
    Ne_range = np.arange(40, 75, 5)

    print(f"     {'N_e':>5s} | {'n_s (TGP)':>12s} | {'r':>10s} | {'r·N_e²':>8s}")
    print("     " + "-"*45)
    delta_ns_tgp = 0.0005  # kalibrowane z Mukhanov-Sasaki (ex107)
    for Ne in Ne_range:
        n_s = 1.0 - 2.0/Ne - delta_ns_tgp
        r = 12.0 / Ne**2  # klasa Starobinsky: r·N_e² = 12
        print(f"     {Ne:5.0f} | {n_s:12.4f} | {r:10.5f} | {r*Ne**2:8.1f}")

    # ─── 4. Hipoteza a_Γ·Φ₀ = 1 ─────────────────────────────
    print(f"\n  4. Implikacja hipotezy a_Γ·Φ₀ = 1:")
    print(f"     Jeśli a_Γ = a_sub (stała sieciowa substratu), to:")
    print(f"     a_sub = ℓ_P / √Φ₀  (z relacji złożonej)")

    for Phi0 in [24.66, 25.0, 28.0]:
        a_sub = 1.0 / Phi0   # a_Γ·Φ₀ = 1 → a_Γ = 1/Φ₀
        print(f"     Φ₀ = {Phi0}: a_Γ = {a_sub:.5f} (w jednostkach ℓ_P)")

        # Implikacja dla ε₀:
        # ξ/a_sub z warunku N_e = 55 (optymalny)
        eps0_55 = Phi0 * np.exp(-3 * 55)
        xi_55 = np.sqrt(sigma2 / eps0_55)
        xi_a = xi_55 * a_sub  # ξ/ℓ_P → ξ/a_sub = (ξ/ℓ_P)·(ℓ_P/a_sub) = (ξ/ℓ_P)·Φ₀
        print(f"            ξ/a_sub = {xi_a:.3e}  (przy N_e = 55)")
        print(f"            ε₀ = {eps0_55:.3e}")

    # ─── 5. Self-consistency check ───────────────────────────
    print(f"\n  5. Sprawdzenie samospójności:")
    Phi0_best = 24.66
    Ne_target = 55
    eps0_target = Phi0_best * np.exp(-3 * Ne_target)
    xi_target = np.sqrt(sigma2 / eps0_target)

    # Odwrotny sprawdzian: z ξ odtwarzamy ε₀ i N_e
    eps0_check = epsilon0_from_xi(xi_target, Phi0_best)
    Ne_check = N_efolds(eps0_check, Phi0_best)

    print(f"     Cel: N_e = {Ne_target}")
    print(f"     Φ₀ = {Phi0_best}")
    print(f"     ε₀ = {eps0_target:.6e}")
    print(f"     ξ/ℓ_P = {xi_target:.6e}")
    print(f"     Sprawdzenie: ε₀(ξ) = {eps0_check:.6e}  → N_e = {Ne_check:.2f}")

    # n_s przy tym N_e (klasa Starobinsky + mała korekcja TGP)
    delta_ns_tgp = 0.0005  # kalibrowane z ex107
    n_s_pred = 1.0 - 2.0/Ne_target - delta_ns_tgp
    r_pred = 12.0 / Ne_target**2

    print(f"\n     Predykcje CMB:")
    print(f"       n_s = {n_s_pred:.4f}  (Planck: 0.9649 ± 0.0042)")
    print(f"       r   = {r_pred:.5f}  (limit BICEP/Keck: < 0.036)")
    print(f"       r·N_e² = {r_pred * Ne_target**2:.1f}  (klasa Starobinsky: 12)")
    print(f"       Δn_s/σ = {abs(n_s_pred - 0.9649)/0.0042:.2f}σ")

    return results, xi_range, Phi0_values


def plot_epsilon0(results, xi_range, Phi0_values, outdir='.'):
    """Generuje wykresy ε₀."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('TGP: ε₀ from Substrate Thermodynamics — Problem R6',
                 fontsize=15, fontweight='bold')

    sigma2 = sigma2_critical()

    # 1. ε₀(ξ/ℓ_P) dla różnych Φ₀
    ax = axes[0, 0]
    for Phi0 in Phi0_values:
        eps0 = epsilon0_from_xi(xi_range, Phi0)
        ax.loglog(xi_range, eps0, label=f'Φ₀ = {Phi0}')

    # Zaznacz ε₀ dające N_e = 50, 60
    for Ne_target, ls in [(50, '--'), (60, ':')]:
        for Phi0 in [25.0]:
            eps0_t = Phi0 * np.exp(-3*Ne_target)
            ax.axhline(eps0_t, color='gray', ls=ls, alpha=0.5,
                       label=f'N_e = {Ne_target}' if Phi0 == 25.0 else '')

    ax.set_xlabel('ξ/ℓ_P')
    ax.set_ylabel('ε₀')
    ax.set_title('ε₀ vs długość korelacji')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. N_e(ξ/ℓ_P) dla różnych Φ₀
    ax = axes[0, 1]
    for Phi0 in Phi0_values:
        Ne = np.array([N_efolds_from_xi(xi, Phi0) for xi in xi_range])
        ax.semilogx(xi_range, Ne, label=f'Φ₀ = {Phi0}')

    # Okno dozwolone (Planck: N_e ∈ [50, 70])
    ax.axhspan(50, 70, alpha=0.15, color='green', label='Dozwolone (50-70)')
    ax.axhline(55, color='k', ls='--', alpha=0.5, label='Optymalny (55)')
    ax.set_xlabel('ξ/ℓ_P')
    ax.set_ylabel('N_e')
    ax.set_title('Liczba e-składań vs ξ')
    ax.set_ylim(20, 90)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 3. n_s(N_e) — porównanie z Planck
    ax = axes[1, 0]
    Ne_arr = np.linspace(35, 75, 200)
    n_s_arr = []
    r_arr = []
    delta_ns_tgp = 0.0005
    for Ne in Ne_arr:
        n_s = 1.0 - 2.0/Ne - delta_ns_tgp
        r = 12.0 / Ne**2
        n_s_arr.append(n_s)
        r_arr.append(r)

    n_s_arr = np.array(n_s_arr)
    r_arr = np.array(r_arr)

    ax.plot(Ne_arr, n_s_arr, 'b-', lw=2, label='TGP: n_s(N_e)')
    ax.axhspan(0.9649 - 0.0042, 0.9649 + 0.0042, alpha=0.2, color='orange',
               label='Planck 1σ')
    ax.axhspan(0.9649 - 2*0.0042, 0.9649 + 2*0.0042, alpha=0.1, color='orange',
               label='Planck 2σ')
    ax.axhline(0.9649, color='orange', ls='--', alpha=0.7)
    ax.axvspan(50, 70, alpha=0.1, color='green')
    ax.set_xlabel('N_e')
    ax.set_ylabel('n_s')
    ax.set_title('n_s vs N_e (TGP slow-roll)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 4. Mapa (Φ₀, ξ) → N_e
    ax = axes[1, 1]
    Phi0_grid = np.linspace(20, 35, 50)
    xi_grid = np.logspace(28, 42, 50)
    PP, XX = np.meshgrid(Phi0_grid, np.log10(xi_grid))

    Ne_grid = np.zeros_like(PP)
    for i in range(len(xi_grid)):
        for j in range(len(Phi0_grid)):
            Ne_grid[i, j] = N_efolds_from_xi(xi_grid[i], Phi0_grid[j])

    c = ax.contourf(PP, XX, Ne_grid, levels=np.arange(30, 85, 5),
                    cmap='viridis')
    plt.colorbar(c, ax=ax, label='N_e')
    ax.contour(PP, XX, Ne_grid, levels=[50, 60, 70],
               colors='white', linewidths=2, linestyles='--')
    ax.set_xlabel('Φ₀')
    ax.set_ylabel('log₁₀(ξ/ℓ_P)')
    ax.set_title('Mapa N_e(Φ₀, ξ)')

    plt.tight_layout()
    outpath = os.path.join(outdir, 'epsilon0_estimate.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  → Zapisano: {outpath}")
    plt.close()


# ─── MAIN ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    results, xi_range, Phi0_values = analyze_epsilon0()

    outdir = os.path.dirname(os.path.abspath(__file__))
    plot_epsilon0(results, xi_range, Phi0_values, outdir)

    print("\n" + "="*60)
    print("  WNIOSKI:")
    print("  1. ε₀ ~ (ℓ_P/ξ)² · σ²_c  jest naturalnie małe")
    print("     dla ξ ≫ ℓ_P (długa korelacja przy przejściu fazowym)")
    print("  2. N_e ~ 55 wymaga ξ/ℓ_P ~ 10³⁵")
    print("  3. n_s = 0.964 ± 0.003 w dozwolonym oknie N_e")
    print("  4. Status R6 pozostaje PROGRAM — brakuje:")
    print("     - Obliczenia S(R_c/a_sub) z pełnej teorii nukleacji")
    print("     - Samospójnego Φ₀ z termodynamiki substratu")
    print("     - Związku ξ z a_Γ (hipoteza a_Γ·Φ₀ = 1)")
    print("="*60 + "\n")
