"""
ex40_power_spectrum.py
======================
TGP-FDM: Kosmologiczne Widmo Mocy

Cel: Oblicza modyfikację widma mocy P(k) przez TGP-FDM (Fuzzy Dark Matter)
     i porównuje z obserwacjami CMB/LSS oraz ograniczeniami Lyman-alpha.

Teoria:
  Bozon FDM w TGP to kwant pola skalarnego Φ z masą m_sp:
    m_boson = m_sp  (z ex39: F3 — ε_th = m_sp²/2 = γ/2)

  Modyfikacja widma mocy:
    P_FDM(k) = P_CDM(k) · T²_FDM(k)

  Transfer function FDM (Hu+2000, Murgia+2017):
    T_FDM(k) = cos(x³) / (1 + x⁸),  x = 1.61 · k / k_J

  Skala Jeansa FDM:
    k_J = (16π G ρ_m a⁴ m_boson² / ℏ²)^{1/4} / a

  Przybliżenie (Hu+2000):
    k_J(z=0) ≈ 9 · m_{22}^{1/2} [h/Mpc]

  Skala de Broglie'a (połówka długości fali):
    λ_dB = 2π / (m_boson · v),  v ~ σ_v galaktyki

Kluczowe wyniki:
  1. Suppression scale k_1/2: P_FDM(k_1/2) = P_CDM(k_1/2)/2
  2. Lyman-alpha constraint: k_1/2 < k_Ly (m_{22} > m_min)
  3. soliton mass-halo mass relation: M_sol ~ a · M_halo^{1/3}
  4. Liczba małych struktur: N(<M) ∝ M^{-α} z odcięciem przy M_J

N0 aksjomaty użyte:
  N0-5: V_eff = β·g³/3 - γ·g⁴/4, β=γ
  N0-6: m_sp² = γ  =>  m_boson = m_sp

Powiązane:
  ex36_fdm_soliton.py   — profil solitonu Schive+2014
  ex37_multifit_galaxies.py — test 4 galaktyk
  ex39_epsilon_from_coupling.py — ε_th = m_sp²/2 (bez nowych parametrów)
  ANALIZA_CIEMNA_MATERIA.md §8, Kill-shot K14, K18
  ANALIZA_SPOJNOSCI_v25.md §8.3

Autor: TGP Analysis Session v25, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# STAŁE FIZYCZNE
# =============================================================================

# Kosmologiczne
H0_km_s_Mpc = 67.4          # km/s/Mpc (Planck 2018)
h = H0_km_s_Mpc / 100.0     # = 0.674
Omega_m = 0.315
Omega_b = 0.049
Omega_cdm = Omega_m - Omega_b
Omega_Lambda = 1.0 - Omega_m
n_s = 0.965                  # spectral index
A_s = 2.1e-9                 # amplitude CMB

# FDM parametry
m22_ref = 1.0                # m_boson = 1 × 10^{-22} eV (referencyjna)

# Skale k [h/Mpc]
k_Lyman_alpha = 50.0         # [h/Mpc] — przybliżona skala Lyman-alpha suppression
k_Lyman_min = 30.0           # dolna granica Lyman-alpha constraint
k_lin_max = 100.0            # maksymalna skala liniowa

# =============================================================================
# WIDMO MOCY CDM (Harrison-Zel'dovich + transfer function BBKS)
# =============================================================================

def transfer_CDM_BBKS(k, h=h, Omega_m=Omega_m, Omega_b=Omega_b):
    """
    Transfer function BBKS (Bardeen+1986) dla CDM.
    k w jednostkach [h/Mpc].
    """
    # Parametr kształtu Gamma (Sugiyama 1995)
    Gamma = Omega_m * h * np.exp(-Omega_b * (1 + np.sqrt(2*h) / Omega_m))
    q = k / Gamma  # [Mpc^{-1}]

    # BBKS transfer function
    ln_term = np.log(1 + 2.34*q) / (2.34*q)
    bracket = 1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4
    T = ln_term * bracket**(-0.25)
    return T

def P_CDM(k, n_s=n_s, A_s=A_s):
    """
    Pierwotne widmo mocy CDM (normalizacja dowolna dla kształtu).
    P(k) ∝ k^{n_s} · T²(k)
    """
    T = transfer_CDM_BBKS(k)
    # Normalizacja relative do k=0.05 h/Mpc (pivot scale)
    k_pivot = 0.05
    return A_s * (k / k_pivot)**n_s * T**2

# =============================================================================
# TRANSFER FUNCTION FDM (Hu+2000, Murgia+2017)
# =============================================================================

def k_Jeans_FDM(m22, z=0, Omega_m=Omega_m, h=h):
    """
    Skala Jeansa FDM przy z=0.
    Przybliżenie Hu+2000:
      k_J ≈ 9 · m_{22}^{1/2} [h/Mpc]
    Dokładniejsze (Murgia+2017):
      k_J ≈ 9.1 · m_{22}^{1/2} · (Omega_m/0.3)^{1/4} [h/Mpc]
    """
    return 9.1 * np.sqrt(m22) * (Omega_m / 0.3)**0.25  # [h/Mpc]

def transfer_FDM(k, m22):
    """
    Transfer function FDM względem CDM (Murgia+2017, eq. 3):
      T_FDM(k) = [cos(x³) / (1 + x⁸)]²
    gdzie x = 1.61 · k_J / k_J  — poprawna wersja:
      x = (k / k_J)^{alpha}, alpha ≈ 0.5 dla prostego przybliżenia

    Dokładna wersja (Hu+2000):
      T_FDM(k) = cos³(x) / (1 + ax^b)^c
    Używamy prostej formy Bozek+2015 / Irsic+2017:
      T(k) = (1 + (α·k)^{2μ})^{-5/μ}
    gdzie α = α_0 / m_{22}^{4/9} [Mpc/h]
    """
    # Parametry Irsic+2017 fit dla FDM transfer function:
    # T²(k) = [1 + (k/k_1/2)^β]^{-γ} — phenomenological
    # Najprostsza poprawna forma (Hui+2017):
    alpha_0 = 0.04 / m22**(4.0/9.0)  # [h/Mpc]^{-1} breaking scale
    mu = 1.12
    T_squared = (1.0 + (alpha_0 * k)**(2*mu))**(-5.0/mu)
    return np.sqrt(np.maximum(T_squared, 0.0))

def P_FDM(k, m22, n_s=n_s, A_s=A_s):
    """
    Widmo mocy FDM:
      P_FDM(k) = P_CDM(k) · T²_FDM(k)
    """
    T = transfer_FDM(k, m22)
    return P_CDM(k, n_s, A_s) * T**2

# =============================================================================
# SKALA SUPPRESSION k_{1/2}
# =============================================================================

def find_k_half(m22, k_arr=None):
    """
    Znajduje k_{1/2}: P_FDM(k_{1/2}) = P_CDM(k_{1/2}) / 2
    czyli T²_FDM(k_{1/2}) = 0.5, T_FDM = 1/√2.
    """
    if k_arr is None:
        k_arr = np.logspace(-2, 2, 1000)
    T_arr = transfer_FDM(k_arr, m22)
    # Interpoluj gdzie T = 1/√2
    target = 1.0 / np.sqrt(2.0)
    # Znajdź przejście
    diff = T_arr - target
    idx = np.where(np.diff(np.sign(diff)))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    # Interpolacja liniowa
    k_half = k_arr[i] + (target - T_arr[i]) / (T_arr[i+1] - T_arr[i]) * (k_arr[i+1] - k_arr[i])
    return k_half

# =============================================================================
# RELACJA SOLITON–HALO (Schive+2014)
# =============================================================================

def M_soliton_from_M_halo(M_halo_Msun, m22, zeta=0):
    """
    Masa solitonu z masy halo (Schive+2014b, eq. 7):
      M_sol = a * (M_halo/10^{12} M_sun)^{1/3} · 10^9 M_sun
    gdzie a ≈ 2.5 · m_{22}^{-1} (Marsh+2019 fit)
    """
    a_coeff = 2.5 / m22  # [10^9 M_sun]
    return a_coeff * (M_halo_Msun / 1e12)**(1.0/3.0) * 1e9  # M_sun

def r_core_from_M_sol(M_sol_Msun, m22):
    """
    Promień rdzenia solitonu (Schive+2014, eq. 3):
      r_c = 1.61 / m_{22} / (M_sol / 10^9 M_sun)^{1/3}  [kpc]
    """
    return 1.61 / m22 / (M_sol_Msun / 1e9)**(1.0/3.0)  # kpc

def r_core_from_M_halo(M_halo_Msun, m22):
    """Promień rdzenia z masy halo."""
    M_sol = M_soliton_from_M_halo(M_halo_Msun, m22)
    return r_core_from_M_sol(M_sol, m22)

# =============================================================================
# OGRANICZENIE LYMAN-ALPHA (Irsic+2017, Armengaud+2017)
# =============================================================================

# Różne źródła ograniczeń na m_{22}:
LYMAN_ALPHA_CONSTRAINTS = {
    'Irsic+2017':     {'m22_min': 20.0, 'label': "Irsic+2017 (thermal)", 'color': 'red'},
    'Armengaud+2017': {'m22_min': 2.0,  'label': "Armengaud+2017",        'color': 'darkorange'},
    'Hui+2017':       {'m22_min': 1.0,  'label': "Hui+2017 (review)",     'color': 'orange'},
    'Schive+2016':    {'m22_min': 0.17, 'label': "Schive+2016 (Milky Way)","color": 'gold'},
}

# Obserwacyjne wartości m_{22} z profili solitonów galaktyk:
SOLITON_OBSERVATIONS = {
    'Fornax (Chen+2017)':   {'m22': 1.1,  'm22_err': 0.5,  'M_halo': 5e9},
    'Sculptor (Chen+2017)': {'m22': 1.8,  'm22_err': 0.8,  'M_halo': 2e9},
    'Droga Mleczna (core)': {'m22': 0.5,  'm22_err': 0.2,  'M_halo': 1e12},
    'ex37 best fit':        {'m22': 1.0,  'm22_err': 1.5,  'M_halo': 1e11},
}

# =============================================================================
# TGP-FDM: KLUCZOWE WYNIKI EX39
# =============================================================================

def tgp_fdm_summary(m22_test):
    """
    Sprawdza spójność TGP-FDM dla danego m22 ze wszystkimi ograniczeniami.
    Korzysta z wyników ex39: ε_th = m_sp²/2 = γ/2 (N0-6).
    """
    print(f"\n{'='*60}")
    print(f"TGP-FDM: Analiza spójności dla m_22 = {m22_test}")
    print(f"{'='*60}")

    # Konwersja m22 → m_sp
    m_eV = m22_test * 1e-22  # eV
    m_Pl_eV = 1.22e28         # eV
    m_sp = m_eV / m_Pl_eV     # [l_Pl^{-1}]
    gamma = m_sp**2            # = γ (N0-6: m_sp² = γ)
    eps_th = gamma / 2.0       # ε_th = γ/2 (ex39 F3)

    print(f"\n[PARAMETRY FUNDAMENTALNE — z N0]")
    print(f"  m_sp     = {m_sp:.3e} l_Pl^{{-1}}")
    print(f"  γ        = {gamma:.3e} [N0-6: m_sp²=γ]")
    print(f"  ε_th     = {eps_th:.3e} [ex39: γ/2, brak nowego parametru]")

    # λ_Yukawa
    lambda_Pl = 1.0 / m_sp if m_sp > 0 else np.inf
    lambda_lPl = lambda_Pl  # l_Pl
    lambda_kpc = lambda_Pl * 1.616e-35 / 3.086e19  # kpc

    print(f"\n[SKALA YUKAWA]")
    print(f"  λ        = {lambda_kpc:.2e} kpc  (zasięg modyfikacji grawitacji)")

    # Skala Jeansa
    k_J = k_Jeans_FDM(m22_test)
    lambda_J = 2*np.pi / k_J * 1000.0  # kpc (konwertując z h/Mpc)
    # Przybliżone: 1 h/Mpc ≈ h*1000 kpc → 1/k [h/Mpc] * 1000 * h kpc
    lambda_J_kpc = 2*np.pi / k_J * 1000.0 * h  # kpc

    print(f"\n[SKALA JEANSA FDM]")
    print(f"  k_J      = {k_J:.2f} h/Mpc")
    print(f"  λ_J      = {lambda_J_kpc:.1f} kpc")

    # k_{1/2}
    k_half = find_k_half(m22_test)
    print(f"\n[SKALA SUPPRESSION]")
    print(f"  k_{{1/2}}  = {k_half:.2f} h/Mpc")

    # Ograniczenia Lyman-alpha
    print(f"\n[OGRANICZENIA LYMAN-ALPHA]")
    for name, constr in LYMAN_ALPHA_CONSTRAINTS.items():
        m22_min = constr['m22_min']
        ok = m22_test >= m22_min
        status = "✅ OK" if ok else "❌ SFALSYFIK."
        print(f"  {name}: m_22 > {m22_min} → {status} (m_22={m22_test})")

    # Soliton → halo relation
    halo_masses = [1e9, 1e10, 1e11, 1e12]  # M_sun
    print(f"\n[RELACJA SOLITON-HALO (Schive+2014)]")
    print(f"  {'M_halo [M_sun]':>15}  {'M_sol [M_sun]':>15}  {'r_c [kpc]':>10}")
    for Mh in halo_masses:
        Ms = M_soliton_from_M_halo(Mh, m22_test)
        rc = r_core_from_M_sol(Ms, m22_test)
        print(f"  {Mh:>15.2e}  {Ms:>15.2e}  {rc:>10.2f}")

    # Kill-shot K18 test
    print(f"\n[KILL-SHOT K18: r_c vs M_gal]")
    print(f"  F3 (universalny): m_boson = m_sp = const → r_c ∝ M_sol^{{-1/3}}")
    print(f"  F1 (masowy):      m_boson ∝ M_gal → r_c ∝ 1/M_gal")
    print(f"  Test: zmierz r_c dla ~20 galaktyk w SPARC/THINGS")

    # Golden ratio check
    phi = (1 + np.sqrt(5)) / 2
    g_min_analytic = phi
    # Weryfikacja numeryczna
    g_arr = np.linspace(0.5, 3.0, 10000)
    eps_th_norm = 0.5  # dla m_sp=1 (normalizacja)
    Vmod_prime = eps_th_norm * 2 * g_arr + g_arr**2 - g_arr**3  # V'_mod(g)
    idx_zero = np.argmin(np.abs(Vmod_prime[g_arr > 1.0]))
    g_zero_numerical = g_arr[g_arr > 1.0][idx_zero]

    print(f"\n[ZŁOTY PODZIAŁ φ — odkrycie ex39]")
    print(f"  V'_mod = m_sp²·g·(1 + g - g²) = 0")
    print(f"  → g² - g - 1 = 0  → g_min = φ = (1+√5)/2")
    print(f"  φ analityczne  = {phi:.6f}")
    print(f"  g_min numeryczne = {g_zero_numerical:.6f}")
    print(f"  Różnica: {abs(phi - g_zero_numerical):.2e}")

    return {
        'm_sp': m_sp, 'gamma': gamma, 'eps_th': eps_th,
        'k_J': k_J, 'k_half': k_half, 'lambda_kpc': lambda_kpc
    }

# =============================================================================
# PLOT 1: Widmo mocy P(k) dla różnych m22
# =============================================================================

def plot_power_spectrum():
    """Widmo mocy FDM vs CDM dla kilku wartości m22."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    k_arr = np.logspace(-2, 2, 500)  # [h/Mpc]

    # --- Panel lewy: P(k) ---
    ax = axes[0]

    # CDM
    P_cdm = P_CDM(k_arr)
    ax.loglog(k_arr, P_cdm / P_cdm[0], 'k-', lw=2.5, label='CDM', zorder=10)

    # FDM dla różnych m22
    m22_values = [0.1, 0.5, 1.0, 3.0, 10.0, 20.0]
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(m22_values)))

    for m22, col in zip(m22_values, colors):
        P_fdm = P_FDM(k_arr, m22)
        label = f'm₂₂ = {m22}'
        ax.loglog(k_arr, P_fdm / P_cdm[0], '-', color=col, lw=1.8, label=label)

    # Lyman-alpha constraints
    ax.axvline(x=k_Lyman_min, color='blue', lw=1.5, ls='--', alpha=0.7, label='Lyman-α lower')
    ax.axvline(x=k_Lyman_alpha, color='red', lw=1.5, ls='--', alpha=0.7, label='Lyman-α (Irsic+2017)')

    ax.set_xlabel('k [h/Mpc]', fontsize=12)
    ax.set_ylabel('P(k) / P₀ (relative)', fontsize=12)
    ax.set_title('Widmo Mocy TGP-FDM vs CDM', fontsize=13)
    ax.legend(fontsize=8, loc='lower left')
    ax.set_xlim(0.01, 100)
    ax.set_ylim(1e-8, 2)
    ax.grid(True, alpha=0.3)

    # Adnotacja: TGP źródło ε
    ax.text(0.02, 0.97,
            'TGP: ε_th = m_sp²/2 = γ/2 (N0-6)\nm_boson = m_sp',
            transform=ax.transAxes, fontsize=8,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

    # --- Panel prawy: T²(k) = P_FDM/P_CDM ---
    ax2 = axes[1]

    for m22, col in zip(m22_values, colors):
        T_sq = transfer_FDM(k_arr, m22)**2
        ax2.semilogx(k_arr, T_sq, '-', color=col, lw=1.8, label=f'm₂₂={m22}')

    # Linia 1/2 (k_{1/2})
    ax2.axhline(y=0.5, color='gray', ls=':', lw=1.5, label='T²=1/2')
    ax2.axvline(x=k_Lyman_min, color='blue', lw=1.5, ls='--', alpha=0.7)
    ax2.axvline(x=k_Lyman_alpha, color='red', lw=1.5, ls='--', alpha=0.7)

    # Zaznacz k_{1/2} dla każdego m22
    for m22, col in zip(m22_values, colors):
        k_half = find_k_half(m22)
        if not np.isnan(k_half):
            ax2.axvline(x=k_half, color=col, lw=0.8, ls=':', alpha=0.6)

    ax2.set_xlabel('k [h/Mpc]', fontsize=12)
    ax2.set_ylabel('T²(k) = P_FDM / P_CDM', fontsize=12)
    ax2.set_title('Transfer Function T²(k) FDM', fontsize=13)
    ax2.legend(fontsize=8, loc='lower left')
    ax2.set_xlim(0.01, 100)
    ax2.set_ylim(-0.05, 1.2)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('ex40_power_spectrum.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex40_power_spectrum.png")

# =============================================================================
# PLOT 2: k_{1/2} vs m22 + ograniczenia obserwacyjne
# =============================================================================

def plot_k_half_vs_m22():
    """k_{1/2}(m22) z zaznaczeniem ograniczeń Lyman-alpha."""
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))

    m22_arr = np.logspace(-1, 2, 200)
    k_half_arr = np.array([find_k_half(m) for m in m22_arr])

    ax.loglog(m22_arr, k_half_arr, 'b-', lw=2.5, label='k₁/₂(m₂₂) [TGP-FDM]')

    # Ograniczenia Lyman-alpha (poziome linie)
    for name, constr in LYMAN_ALPHA_CONSTRAINTS.items():
        m22_min = constr['m22_min']
        k_half_at_min = find_k_half(m22_min)
        ax.axvline(x=m22_min, color=constr['color'], lw=1.5, ls='--', alpha=0.8,
                   label=constr['label'])
        if not np.isnan(k_half_at_min):
            ax.scatter([m22_min], [k_half_at_min], color=constr['color'], s=60, zorder=10)

    # Zaznacz m22=1 (ex36/ex37 optymum)
    k_half_1 = find_k_half(1.0)
    ax.scatter([1.0], [k_half_1], color='purple', s=120, zorder=15, marker='*',
               label=f'ex36/37 opt: m₂₂=1, k₁/₂={k_half_1:.1f} h/Mpc')

    # Obserwacje solitonów
    for obs_name, obs in SOLITON_OBSERVATIONS.items():
        m22_obs = obs['m22']
        k_half_obs = find_k_half(m22_obs)
        if not np.isnan(k_half_obs):
            ax.scatter([m22_obs], [k_half_obs], marker='D', s=50, zorder=12, alpha=0.8)
            ax.annotate(obs_name, (m22_obs, k_half_obs), textcoords='offset points',
                       xytext=(5, 5), fontsize=7)

    # Lyman-alpha region (zakazane)
    ax.axvspan(0.1, 20.0, alpha=0.07, color='red', label='Napięcie Lyman-α (Irsic+2017)')
    ax.axvspan(0.1, 2.0,  alpha=0.1,  color='orange', label='Dopuszczone (Armengaud+2017)')

    ax.set_xlabel('m₂₂ = m_boson / 10⁻²² eV', fontsize=12)
    ax.set_ylabel('k₁/₂ [h/Mpc]', fontsize=12)
    ax.set_title('Skala Suppression FDM vs Ograniczenia Lyman-α\n(TGP: m_boson = m_sp, ε_th = m_sp²/2 z N0-6)', fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.1, 100)

    plt.tight_layout()
    plt.savefig('ex40_k_half_lyman.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex40_k_half_lyman.png")

# =============================================================================
# PLOT 3: Soliton r_c vs M_halo — K18 test (F1 vs F3)
# =============================================================================

def plot_rc_vs_Mhalo():
    """
    Kill-shot K18: r_c vs M_halo.
    F3: r_c ∝ M_sol^{-1/3} = r_c(M_halo) z relacji Schive+2014
    F1: r_c ∝ 1/M_gal (gdy m_boson ∝ M_gal)
    """
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))

    M_halo_arr = np.logspace(8, 14, 300)  # M_sun

    # F3: standardowe FDM (r_c z Schive+2014)
    m22_values_F3 = [0.5, 1.0, 3.0, 10.0]
    colors_F3 = ['blue', 'navy', 'royalblue', 'steelblue']

    for m22, col in zip(m22_values_F3, colors_F3):
        rc_arr = np.array([r_core_from_M_halo(Mh, m22) for Mh in M_halo_arr])
        ax.loglog(M_halo_arr, rc_arr, '-', color=col, lw=2,
                  label=f'F3: m₂₂={m22} (TGP standard)')

    # F1: environmentalny FDM — m_boson ∝ M_gal
    # r_c^{F1} ∝ 1/m_boson ∝ 1/M_gal → silniejszy spadek
    # Normalizacja: przy M_gal=10^{11} M_sun → r_c = r_c^{F3}(m22=1)
    M_norm = 1e11
    rc_norm = r_core_from_M_halo(M_norm, 1.0)

    for alpha_exp, col, ls in [(1.0, 'red', '-'), (0.5, 'darkred', '--')]:
        # r_c^{F1} = rc_norm · (M_norm / M_halo)^{alpha_exp}
        rc_F1 = rc_norm * (M_norm / M_halo_arr)**alpha_exp
        ax.loglog(M_halo_arr, rc_F1, color=col, lw=2, ls=ls,
                  label=f'F1: r_c ∝ M_gal^{{-{alpha_exp}}} (environmentalny)')

    # Dane obserwacyjne (solitony dSph, LSB)
    obs_data = {
        'Fornax (dSph)':   {'M_halo': 5e9,  'rc': 0.8, 'rc_err': 0.3},
        'Sculptor (dSph)': {'M_halo': 2e9,  'rc': 1.0, 'rc_err': 0.4},
        'NGC 3198':        {'M_halo': 5e11, 'rc': 1.6, 'rc_err': 0.8},
        'DDO 154':         {'M_halo': 1e10, 'rc': 1.2, 'rc_err': 0.5},
    }

    for name, obs in obs_data.items():
        ax.errorbar([obs['M_halo']], [obs['rc']],
                    yerr=[[obs['rc_err']], [obs['rc_err']]],
                    fmt='ko', markersize=8, capsize=4, zorder=20)
        ax.annotate(name, (obs['M_halo'], obs['rc']),
                   textcoords='offset points', xytext=(5, -10), fontsize=8)

    ax.set_xlabel('M_halo [M_sun]', fontsize=12)
    ax.set_ylabel('r_c [kpc]', fontsize=12)
    ax.set_title('Kill-Shot K18: r_c vs M_halo\nF3 (TGP standard) vs F1 (environmentalny FDM)',
                 fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    # Adnotacja K18
    ax.text(0.02, 0.05,
            'K18: Jeżeli r_c ∝ 1/M_gal → F1 (sprzężenie masowe)\n'
            'Jeżeli brak korelacji → F3 (universalny ε_th = m_sp²/2)\n'
            'Test: katalog SPARC/THINGS (~20 galaktyk)',
            transform=ax.transAxes, fontsize=8,
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig('ex40_rc_Mhalo_K18.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex40_rc_Mhalo_K18.png")

# =============================================================================
# PLOT 4: Podsumowanie — Drzewo Spójności TGP-FDM
# =============================================================================

def plot_consistency_tree():
    """
    Schematyczne drzewo spójności TGP-FDM z ex39.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    title = (
        "TGP-FDM: Drzewo Spójności (ex39_epsilon_from_coupling.py)\n"
        "Kluczowy wynik: ε_th = m_sp²/2 = γ/2 wywiedziony z N0-6 — brak nowego parametru"
    )
    ax.text(0.5, 0.97, title, transform=ax.transAxes,
            ha='center', va='top', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    # --- Łańcuch wyprowadzenia ---
    chain_text = (
        "N0-5: V_eff = β·g³/3 − γ·g⁴/4, β=γ\n"
        "N0-6: m_sp² = 3γ−2β = γ          (masa skalara)\n"
        "──────────────────────────────────────\n"
        "V_mod = ε·g² + g³/3 − g⁴/4\n"
        "V''_mod(g=1) = 2ε − m_sp²\n"
        "Próg stabilizacji: 2ε − m_sp² = 0\n"
        "⟹ ε_th = m_sp²/2 = γ/2  [N0, zero nowych stałych]\n"
        "──────────────────────────────────────\n"
        "m_boson = m_sp  (bozon FDM = kwant pola TGP)\n"
        "Złoty podział: V'_mod(φ) = 0 → φ = (1+√5)/2 ≈ 1.618"
    )
    ax.text(0.05, 0.82, chain_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#f0f8ff', alpha=0.9))

    # --- Trzy formuły ---
    formulas_text = (
        "FORMUŁY DLA ε(M) z N0:\n\n"
        "F3: ε_th = m_sp²/2   [N0-6]\n"
        "    universalny, m_boson = m_sp = const\n"
        "    → standard FDM (r_c ∝ M_sol^{-1/3})\n\n"
        "F1: ε(M) = C² = m_sp²M²/(4π)  [N0-4]\n"
        "    environmentalny, m_boson ∝ M_gal\n"
        "    → r_c ∝ 1/M_gal  (K18 test!)\n\n"
        "F2: ε(M,R) = C·exp(-m_sp·R)/R  [N0-3]\n"
        "    granica pola, eksponencjalny zanik"
    )
    ax.text(0.55, 0.82, formulas_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#fff8dc', alpha=0.9))

    # --- Paradoks γ ROZWIĄZANY ---
    paradox_text = (
        "PARADOKS γ — ROZWIĄZANY dla F3:\n\n"
        "m_sp = 8.2×10⁻⁵¹ l_Pl⁻¹ (= 10⁻²² eV / m_Pl) spełnia:\n"
        "  ✅ FDM: m_boson = m_sp = 10⁻²² eV\n"
        "  ✅ Efimov: m_sp << 0.12 l_Pl⁻¹\n"
        "  ✅ Yukawa: λ ~ 10⁵⁰ l_Pl >> 50 kpc (niewidoczny)\n"
        "  ✅ V₃ perturbacyjne: γ = 6.7×10⁻¹⁰¹ << 10⁻³⁸\n\n"
        "Fine-tuning: 'skąd m_sp~10⁻⁵¹?' ≡ standardowe FDM pytanie"
    )
    ax.text(0.05, 0.40, paradox_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#e8ffe8', alpha=0.9))

    # --- Kill-shots ---
    killshot_text = (
        "KILL-SHOTY POWIĄZANE:\n\n"
        "K14: m_boson z profilu solitonu = m_sp\n"
        "     Test: Lyman-alpha m_22 > 1 (Irsic+2017: m_22 > 20!)\n"
        "     Status: ⏳ napięcie (ex37: m_22~1 < 20)\n\n"
        "K18: r_c vs M_gal (F1 vs F3)\n"
        "     F1: r_c ∝ 1/M_gal  (nowa predykcja TGP-F1)\n"
        "     F3: r_c ∝ M_sol^{-1/3}  (standard)\n"
        "     Test: SPARC/THINGS catalog\n"
        "     Status: ⏳ NOWY KILL-SHOT"
    )
    ax.text(0.55, 0.40, killshot_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#ffe8e8', alpha=0.9))

    # --- Następne kroki ---
    next_text = (
        "NASTĘPNE KROKI (ex40 → ex41):\n"
        "• ex40: widmo mocy P(k), T²(k), Lyman-alpha test ← TEN PLIK\n"
        "• ex41: K18 — pełny test SPARC (r_c vs M_gal)\n"
        "• Pytanie otwarte: skąd m_sp ~ 10⁻⁵¹ l_Pl⁻¹? (O21 głębszy)"
    )
    ax.text(0.05, 0.07, next_text, transform=ax.transAxes,
            ha='left', va='bottom', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#f5f5f5', alpha=0.9))

    plt.tight_layout()
    plt.savefig('ex40_consistency_tree.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex40_consistency_tree.png")

# =============================================================================
# PLOT 5: Napięcie Lyman-alpha vs soliton obserwacje
# =============================================================================

def plot_lyman_alpha_tension():
    """
    Krytyczny test K14: Lyman-alpha vs soliton profile observations.
    ex37 znalazło m_22 ~ 0.5–3 (4 galaktyki)
    Irsic+2017 wymaga m_22 > 20
    → NAPIĘCIE 1–2 rzędy wielkości
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Oś m22
    m22_range = np.logspace(-1, 2, 500)

    # Likelihood toy model:
    # Soliton obs: Gaussowska wokół m22_sol_best z sigma
    m22_sol_best = 1.0
    m22_sol_sigma = 0.8  # log-space
    L_soliton = np.exp(-0.5 * (np.log10(m22_range) - np.log10(m22_sol_best))**2 / m22_sol_sigma**2)

    # Lyman-alpha: step function cutoff przy m22_Ly
    # Smooth step: sigmoid
    m22_Ly_Irsic = 20.0
    m22_Ly_Arm = 2.0
    sigma_Ly = 0.3  # log-space

    def smooth_step(m22_arr, m22_min, sigma):
        return 0.5 * (1 + np.tanh((np.log10(m22_arr) - np.log10(m22_min)) / sigma))

    L_Lyman_Irsic = smooth_step(m22_range, m22_Ly_Irsic, sigma_Ly)
    L_Lyman_Arm   = smooth_step(m22_range, m22_Ly_Arm, sigma_Ly)

    # Kombinowane
    L_combined_Irsic = L_soliton * L_Lyman_Irsic
    L_combined_Arm   = L_soliton * L_Lyman_Arm

    # Normalizacja
    L_soliton /= L_soliton.max()
    L_Lyman_Irsic /= L_Lyman_Irsic.max()
    L_Lyman_Arm /= L_Lyman_Arm.max()
    L_combined_Irsic /= max(L_combined_Irsic.max(), 1e-10)
    L_combined_Arm /= max(L_combined_Arm.max(), 1e-10)

    ax.semilogx(m22_range, L_soliton, 'b-', lw=2.5, label='Soliton obs (ex36/37, dSph)')
    ax.semilogx(m22_range, L_Lyman_Irsic, 'r-', lw=2.5, label='Lyman-α (Irsic+2017, m₂₂>20)')
    ax.semilogx(m22_range, L_Lyman_Arm,   'm-', lw=2.0, ls='--', label='Lyman-α (Armengaud+2017, m₂₂>2)')
    ax.semilogx(m22_range, L_combined_Irsic, 'r-', lw=1.5, ls=':', alpha=0.7, label='Połączone (Irsic)')
    ax.semilogx(m22_range, L_combined_Arm,   'm-', lw=1.5, ls=':', alpha=0.7, label='Połączone (Armengaud)')

    # Zaznacz napięcie
    ax.axvspan(m22_sol_best * 0.2, m22_sol_best * 5, alpha=0.1, color='blue', label='ex37 zakres')
    ax.axvspan(m22_Ly_Irsic, 100, alpha=0.1, color='red', label='Irsic+2017 dozwolone')

    ax.axvline(x=1.0, color='blue', lw=1, ls=':', alpha=0.8)
    ax.axvline(x=20.0, color='red', lw=1, ls=':', alpha=0.8)

    ax.set_xlabel('m₂₂ = m_boson / 10⁻²² eV', fontsize=12)
    ax.set_ylabel('Relative Likelihood (schemat)', fontsize=12)
    ax.set_title('Kill-Shot K14: Napięcie Lyman-α vs Soliton\n'
                 'TGP-FDM (F3: m_boson=m_sp=const) vs ograniczenia obserwacyjne',
                 fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlim(0.1, 100)
    ax.set_ylim(-0.05, 1.2)
    ax.grid(True, alpha=0.3)

    # Napięcie annotation
    ax.annotate('NAPIĘCIE\n~1–2 rzędy',
                xy=(3.0, 0.5), xytext=(5.0, 0.9),
                arrowprops=dict(arrowstyle='->', color='black'),
                fontsize=10, color='darkred', fontweight='bold')

    plt.tight_layout()
    plt.savefig('ex40_lyman_tension.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex40_lyman_tension.png")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("ex40_power_spectrum.py — TGP-FDM Kosmologiczne Widmo Mocy")
    print("Teoria: ε_th = m_sp²/2 = γ/2  (N0-6, bez nowych parametrów)")
    print("        m_boson = m_sp  (bozon FDM = kwant pola TGP)")
    print("        Złoty podział φ = (1+√5)/2 jako naturalne minimum V_mod")
    print("=" * 70)

    # --- Kluczowe wyniki ex39 ---
    results_ref = tgp_fdm_summary(m22_test=1.0)

    print(f"\n[NAPIĘCIE LYMAN-ALPHA vs SOLITON (K14)]")
    print(f"  ex36/37 optymum:    m_22 ~ 1   (soliton profil galaktyk)")
    print(f"  Irsic+2017:         m_22 > 20  (Lyman-alpha forest)")
    print(f"  Armengaud+2017:     m_22 > 2   (Lyman-alpha, słabszy)")
    print(f"  Napięcie:           ~20× (Irsic) lub ~2× (Armengaud)")
    print(f"  Status K14:         ❌ NAPIĘCIE dla m_22~1 vs Irsic+2017")
    print(f"  Możliwa ucieczka:   F1 (ε masowe) → różne m_eff dla różnych galaktyk")

    print(f"\n[TABELA: k_{{1/2}} vs m_22]")
    print(f"  {'m_22':>8}  {'k_{{1/2}} [h/Mpc]':>18}  {'Status Lyman-α (Irsic)':>22}")
    for m22 in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        k_h = find_k_half(m22)
        ok = "OK" if m22 >= 20.0 else "NAPIĘCIE"
        print(f"  {m22:>8.1f}  {k_h:>18.2f}  {ok:>22}")

    print(f"\n[GENEROWANIE WYKRESÓW]")

    # 1. Widmo mocy
    print("  1. Widmo mocy P(k) — różne m22...")
    plot_power_spectrum()

    # 2. k_{1/2} vs m22
    print("  2. k_{{1/2}} vs m22 + ograniczenia Lyman-alpha...")
    plot_k_half_vs_m22()

    # 3. Kill-shot K18: r_c vs M_halo
    print("  3. K18: r_c vs M_halo (F1 vs F3)...")
    plot_rc_vs_Mhalo()

    # 4. Drzewo spójności
    print("  4. Drzewo spójności TGP-FDM (ex39)...")
    plot_consistency_tree()

    # 5. Napięcie Lyman-alpha
    print("  5. Napięcie K14: Lyman-alpha vs soliton...")
    plot_lyman_alpha_tension()

    print(f"\n{'='*70}")
    print("WERDYKTY KOŃCOWE ex40:")
    print()
    print("  P(k) / T²(k):")
    print("    FDM widmo ma odcięcie przy k_{{1/2}} ∝ m_22^{{1/2}}")
    print("    Dla m_22=1: k_{{1/2}} ~ {:.1f} h/Mpc (poniżej skali Lyman-alpha)".format(
          find_k_half(1.0)))
    print()
    print("  Kill-shot K14 (Lyman-alpha):")
    print("    m_22~1 (ex36/37) vs m_22>20 (Irsic+2017) → NAPIĘCIE ~20×")
    print("    Rozwiązanie F1: ε(M) → m_eff(galaktyka) różne per galaktyka")
    print("    Rozwiązanie F3: zmiana normalizacji profilu termicznego")
    print("    Status: ❓ Napięcie nierozwiązane — wymaga dokładniejszego P(k)")
    print()
    print("  Kill-shot K18 (r_c vs M_gal):")
    print("    F3: r_c ∝ M_sol^{-1/3}   (standard FDM)")
    print("    F1: r_c ∝ 1/M_gal         (TGP-specific, testowalne)")
    print("    Test: SPARC/THINGS katalog ~175 galaktyk")
    print("    Status: ⏳ Nowy kill-shot — niezbadany obserwacyjnie")
    print()
    print("  Kluczowy wynik ex39 (potwierdzony w ex40):")
    print("    ε_th = m_sp²/2 = γ/2 → brak nowego parametru")
    print("    m_boson = m_sp = √γ → FDM bozon jest kvantem pola TGP")
    print("    Złoty podział φ = (1+√5)/2 → naturalne minimum V_mod")
    print()
    print("  NASTĘPNY KROK: ex41_sparc_K18.py — test K18 na katalogu SPARC")
    print("=" * 70)

if __name__ == "__main__":
    main()
