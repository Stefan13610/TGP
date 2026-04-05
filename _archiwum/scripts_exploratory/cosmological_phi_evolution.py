#!/usr/bin/env python3
"""
cosmological_phi_evolution.py
==============================
Skrypt dla otwartego problemu O22 TGP:
  „Ewolucja kosmologiczna tła φ_bg(z) przez epoki materii i promieniowania
   przy dynamicznych stałych c(Φ), ℏ(Φ), G(Φ)."

Teoria TGP modyfikuje równania Friedmanna przez dynamiczne stałe:
    c(Φ)  = c₀ · (Φ₀/Φ)^{1/2}
    ℏ(Φ)  = ℏ₀ · (Φ₀/Φ)^{1/2}
    G(Φ)  = G₀ · (Φ₀/Φ)^{1}

co prowadzi do modyfikowanego układu równań ODE dla:
    φ(t) ≡ Φ_bg(t)/Φ₀  (znormalizowane tło)
    a(t)               (czynnik skali)

Zmodyfikowane równania Friedmanna (TGP):
    H² = (8πG(φ)/3) · ρ_tot  →  H² = (8πG₀/3) · ρ_tot / φ
    φ̈ + 3H φ̇ + V'_eff(φ)/Φ₀² = S(ρ)     [równanie pola Φ]

gdzie V'_eff(φ) = (∂/∂φ)(β·φ² − γ·φ³) i S(ρ) jest sprzężeniem ze źródłami.

Gęstości energii (kosmologiczne, z zachowaniem energii w TGP):
    ρ_r(a) = ρ_{r,0} · a^{−4}            [promieniowanie]
    ρ_m(a) = ρ_{m,0} · a^{−3}            [materia zimna]
    ρ_Λ(φ)  ≈ (α/2)(∂_t φ)² + V_eff(φ)  [gęstość ciemnej energii TGP]

─────────────────────────────────────────────────────────────────────────────
ROZSZERZENIE O22: Entropia substratu S_Γ (N0-7, wersja v28)
─────────────────────────────────────────────────────────────────────────────
Substrat Γ=(V,E) posiada entropię termodynamiczną S_Γ zdefiniowaną przez
korelacje spinowe na krawędziach:

    S_Γ[φ] ≡ s₀ · (φ − ln φ − 1)

gdzie s₀ jest gęstością entropiczną substratu (parametr TGP, s₀ > 0).
Własności:
  • S_Γ(φ=1) = 0 (minimum przy próżni — stan równowagi)
  • S_Γ(φ) ≥ 0 dla wszystkich φ > 0 (entropia nieujemna)
  • S_Γ''(φ=1) = s₀ > 0 (stabilność termodynamiczna)

Temperatura substratu T_Γ jest powiązana z kwantowymi fluktuacjami próżni.
Analogia z temperaturą Hawkinga: T_Γ ∝ ℏ·H/(2π). W TGP:
    T_Γ(φ, H) = T₀ · H/H₀ · (Φ₀/Φ)^{1/2} = T₀ · √H² / √φ

Swobodna energia substratu: F_Γ = −T_Γ · S_Γ(φ)

Modyfikuje potencjał efektywny:
    V_eff_total(φ) = V_eff(φ) + F_Γ(φ, H)
                   = (β·φ² − γ·φ³) − T_Γ · s₀ · (φ − ln φ − 1)

Pochodna efektywna (dla równania ruchu):
    dV_eff_total/dφ = (2β·φ − 3γ·φ²) − T_Γ · s₀ · (1 − 1/φ)

Człon entropiczny działa jako dodatkowe pole skalarne sprzężone z T_Γ.

OGRANICZENIE SELF-CONSISTENT (N0-7):
  T_Γ wchodzi przez H² samospójnie, więc układ ODE staje się algebraicznie
  sprzężony. Rozwiązujemy iteracyjnie lub perturbacyjnie (patrz opcja --entropy).

Wyniki z entropią mają status „jakościowy O22" — pełna analiza wymaga:
  (1) Skalowania s₀ z mikroskopowego modelu substratu H_Γ
  (2) Renormalizacji T_Γ przez kwantowe poprawki N0-7
  (3) Perturbacji CMB z δS_Γ (niedostępne w tym skrypcie)

Jednostki: naturalne kosmologiczne (c₀=ℏ₀=1, G₀=1, H₀=1).

Użycie:
    python cosmological_phi_evolution.py [--plot] [--gamma FLOAT] [--w0 FLOAT]
    python cosmological_phi_evolution.py --entropy [--s0 FLOAT] [--T0 FLOAT]
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# ──────────────────────────────────────────────────────────────────────────────
#  Parametry fiducjalne (jednostki: H₀=1, c₀=1)
# ──────────────────────────────────────────────────────────────────────────────

# Parametry kosmologiczne (Planck 2018-like)
OMEGA_M  = 0.315   # gęstość materii
OMEGA_R  = 9.0e-5  # gęstość promieniowania (w tym neutrina)
OMEGA_L  = 0.685   # energia ciemna (Lambda)

# Parametry TGP
GAMMA_TGP = 0.5    # w jednostkach Planck (tutaj bezwymiarowe ~ H₀²)
BETA_TGP  = 0.5    # = γ (warunek próżniowy N0-5)
ALPHA_TGP = 1.0    # współczynnik kinetyczny

# Parametry entropii substratu S_Γ (N0-7 — rozszerzenie O22)
S_GAMMA_0  = 0.01   # gęstość entropii substratu s₀ (bezwymiarowe, małe)
T_GAMMA_0  = 0.10   # temperatura substratu T₀ (w jednostkach H₀; T_Γ(dziś) ≈ T₀)
USE_ENTROPY = False  # domyślnie wyłączone (włącz przez --entropy)

# Amplituda oscylacji φ od próżni: φ(t) = 1 + δφ(t), |δφ| ≪ 1
PHI0_REL_DEVIATION = 1e-4   # δφ(z=0) — dziś pole blisko próżni

# Zakres całkowania (log₁₀(a), a=1 dziś)
LOG_A_INI = -7.0   # a = 10⁻⁷ (epoka promieniowania)
LOG_A_FIN = 0.5    # a = 10^{0.5} ≈ 3.16 (przyszłość)


# ──────────────────────────────────────────────────────────────────────────────
#  Funkcje pomocnicze
# ──────────────────────────────────────────────────────────────────────────────

def H_squared(a, phi, dphi_dN, gamma=GAMMA_TGP, beta=BETA_TGP, alpha=ALPHA_TGP):
    """
    Zmodyfikowany parametr Hubble'a H²(a,φ).
    Używamy zmiennej N = ln(a) jako czasu.
    dφ/dN = dφ/dt · dt/dN = φ̇/H

    H² (TGP) = H²_ΛCDM / φ
    gdzie H²_ΛCDM = H₀² [Ω_m/a³ + Ω_r/a⁴ + Ω_Λ(φ)]
    a Ω_Λ(φ) zawiera poprawki TGP.
    """
    # Gęstości kosmologiczne (w jednostkach H₀²)
    rho_m = OMEGA_M / a**3
    rho_r = OMEGA_R / a**4

    # Gęstość TGP ciemnej energii:
    # ρ_TGP = (α/2)(φ̇)² + V_eff(φ)
    # W zmiennej N: φ̇ = H·(dφ/dN), więc (φ̇)² = H²·(dφ/dN)²
    # Musimy to rozwiązać samospójnie. Na razie: approx φ ≈ 1, ρ_TGP ≈ Ω_Λ
    V_phi = beta * phi**2 - gamma * phi**3
    V_phi0 = beta - gamma  # wartość w φ=1 (próżnia)

    # Poprawka TGP do energii próżniowej (względna)
    delta_V = V_phi - V_phi0

    rho_de = OMEGA_L + delta_V  # ciemna energia = ΛCDM + poprawka TGP

    # Człon kinetyczny φ (H²·(dφ/dN)²·α/2), ale wchodzi do H² samospójnie:
    # H² = (ρ_m + ρ_r + ρ_de + α/2·H²·(dφ/dN)²) / φ
    # → H²·(1 - α/(2φ)·(dφ/dN)²) = (ρ_m + ρ_r + ρ_de)/φ
    kin_coeff = alpha / (2.0 * max(phi, 1e-10))
    denom = 1.0 - kin_coeff * dphi_dN**2

    if denom <= 0:
        denom = 1e-10  # ochrona przed osobliwością

    H2 = (rho_m + rho_r + rho_de) / (max(phi, 1e-10) * denom)
    return max(H2, 0.0)


def entropy_substrate(phi, s0=S_GAMMA_0):
    """
    Entropia substratu Γ (N0-7):
        S_Γ(φ) = s₀ · (φ − ln φ − 1)

    Własności:
      S_Γ(1) = 0  (minimum przy próżni Φ₀)
      S_Γ(φ) ≥ 0  dla φ > 0
      S_Γ''(1) = s₀ > 0  (termodynamicznie stabilna)
    """
    phi_safe = max(float(phi), 1e-12)
    return s0 * (phi_safe - math.log(phi_safe) - 1.0)


def entropy_substrate_arr(phi_arr, s0=S_GAMMA_0):
    """Wektoryzowana S_Γ(φ) dla tablic numpy."""
    phi_safe = np.maximum(phi_arr, 1e-12)
    return s0 * (phi_safe - np.log(phi_safe) - 1.0)


def entropy_substrate_prime(phi, s0=S_GAMMA_0):
    """
    Pochodna S_Γ względem φ:
        dS_Γ/dφ = s₀ · (1 − 1/φ)

    Zero przy φ=1 (extremum entropii przy próżni).
    Dodatnia dla φ>1, ujemna dla φ<1.
    """
    phi_safe = max(float(phi), 1e-12)
    return s0 * (1.0 - 1.0 / phi_safe)


def substrate_temperature(H2, phi, T0=T_GAMMA_0):
    """
    Temperatura substratu Γ (analogia Hawkinga):
        T_Γ(H, φ) = T₀ · √H² / φ^{1/2}

    Interpretacja:
      - T₀: normalizacja (T_Γ(dziś, φ=1) = T₀ · H₀/H₀ = T₀)
      - Rośnie jak H (z ekspansją wstecz) — wyższa T_Γ w przeszłości
      - Maleje jak φ^{1/2} — gęstszy substrat = zimniejszy (więcej stopni swobody)

    Parametr:
      T0: bezwymiarowe (w jednostkach H₀)
    """
    H_val = math.sqrt(max(float(H2), 0.0))
    phi_safe = max(float(phi), 1e-12)
    return T0 * H_val / math.sqrt(phi_safe)


def substrate_temperature_arr(H2_arr, phi_arr, T0=T_GAMMA_0):
    """Wektoryzowana T_Γ dla tablic numpy."""
    H_arr = np.sqrt(np.maximum(H2_arr, 0.0))
    phi_safe = np.maximum(phi_arr, 1e-12)
    return T0 * H_arr / np.sqrt(phi_safe)


def V_eff_prime(phi, beta=BETA_TGP, gamma=GAMMA_TGP,
                H2=0.0, s0=S_GAMMA_0, T0=T_GAMMA_0,
                use_entropy=False):
    """
    Pochodna efektywnego potencjału:
        dV_eff_total/dφ = 2β·φ − 3γ·φ²
                        [− T_Γ·dS_Γ/dφ]  (tylko gdy use_entropy=True)

    Człon entropiczny (N0-7):
        −T_Γ · dS_Γ/dφ = −T_Γ · s₀ · (1 − 1/φ)

    Przy φ>1: człon jest ujemny → odpychanie (pole wypychane od próżni)
    Przy φ<1: człon jest dodatni → przyciąganie (pole ciągnie do φ=1)
    → Entropia stabilizuje pole wokół φ=1 (próżni)
    """
    # Standard TGP
    dvdp = 2 * beta * phi - 3 * gamma * phi**2

    # Człon entropiczny S_Γ
    if use_entropy and s0 > 0:
        T_gamma = substrate_temperature(H2, phi, T0)
        dvdp -= T_gamma * entropy_substrate_prime(phi, s0)

    return dvdp


def sprzezenie_ze_zrodlem(a, phi, beta=BETA_TGP, gamma=GAMMA_TGP):
    """
    Człon źródłowy w równaniu pola Φ z materii i promieniowania.
    W TGP: S(ρ) ∝ (ρ_m + ρ_r/3) · f(φ),
    gdzie f(φ) pochodzi z dynamicznego G(φ) = G₀/φ.
    Przybliżenie: S ≈ (Ω_m/a³ + Ω_r/(3a⁴)) · δG/G
    (zmiana G modyfikuje efektywną masę grawitacyjną).
    """
    rho_m = OMEGA_M / a**3
    rho_r = OMEGA_R / (3 * a**4)  # ciśnienie promieniowania

    # Pochodna G(φ) = G₀/φ: dG/dφ = −G₀/φ²
    # Sprzężenie ∝ T^μ_μ · dG/dφ = −(ρ_m − 3p_m − ρ_r + 3p_r)/φ²
    # Dla zimnej materii: p_m=0, T^μ_μ = −ρ_m
    # Dla promieniowania: ρ_r − 3p_r = 0 (śladowe)
    coupling = -(rho_m) / (phi**2 + 1e-30)
    return coupling


def ode_system(N, y, args_dict):
    """
    Układ ODE dla φ(N) i dφ/dN(N), gdzie N = ln(a).

    y = [φ, u = dφ/dN]

    Równanie pola:
        d²φ/dN² + (3 + d ln H/dN)·dφ/dN + V'_total(φ)/H²(N) = S(a,φ)/H²(N)

    gdzie V'_total zawiera opcjonalny człon entropiczny (gdy use_entropy=True):
        V'_total = V'_TGP − T_Γ · dS_Γ/dφ

    Uproszczenie: d ln H/dN ≈ −(1/2)·(1 + w_eff) · 3
    dla płaskiego wszechświata.
    """
    phi, u = y
    a = np.exp(N)

    beta = args_dict['beta']
    gamma = args_dict['gamma']
    alpha = args_dict['alpha']
    use_entropy = args_dict.get('use_entropy', False)
    s0  = args_dict.get('s0', S_GAMMA_0)
    T0  = args_dict.get('T0', T_GAMMA_0)

    # Oblicz H²
    H2 = H_squared(a, phi, u, gamma, beta, alpha)
    if H2 <= 0:
        H2 = 1e-30

    # d ln(H)/dN ≈ ε_H (parametr slow-roll)
    # ε_H = −Ḣ/H² = −dH/dN / H
    # Dla ΛCDM: ε_H = (3/2)·(ρ_m + 4/3·ρ_r)/(ρ_total)
    rho_m = OMEGA_M / a**3
    rho_r = OMEGA_R / a**4
    rho_tot = H2 * max(phi, 1e-10)  # przybliżone ρ_tot = H²·φ
    if rho_tot > 0:
        eps_H = (3.0/2) * (rho_m + (4.0/3) * rho_r) / rho_tot
    else:
        eps_H = 0.0

    # Efektywny człon tłumienia: (3 − ε_H) · u
    damping = (3.0 - eps_H)

    # Człon potencjałowy (z opcjonalną entropią substratu)
    Vprime = V_eff_prime(phi, beta, gamma,
                         H2=H2, s0=s0, T0=T0,
                         use_entropy=use_entropy)

    # Człon źródłowy
    S = sprzezenie_ze_zrodlem(a, phi, beta, gamma)

    # Równanie: du/dN = −damping·u − Vprime/H² + S/H²
    du_dN = -damping * u - (Vprime - S) / H2

    return [u, du_dN]


# ──────────────────────────────────────────────────────────────────────────────
#  Całkowanie
# ──────────────────────────────────────────────────────────────────────────────

def solve_phi_evolution(gamma=GAMMA_TGP, beta=BETA_TGP, alpha=ALPHA_TGP,
                        phi0_dev=PHI0_REL_DEVIATION, verbose=True,
                        use_entropy=False, s0=S_GAMMA_0, T0=T_GAMMA_0):
    """
    Całkuje ewolucję φ(z) od epoki promieniowania do dziś (i przyszłości).

    Warunki brzegowe (z=0, a=1):
        φ(a=1) = 1 + phi0_dev   [dziś pole blisko próżni]
        dφ/dN  = 0              [pole w spoczynku]

    Całkujemy WSTECZ od a=1 do a=a_ini, a potem do przodu do a=a_fin.

    Parametry entropii (tylko gdy use_entropy=True):
        s0: gęstość entropii substratu s₀ (domyślnie S_GAMMA_0)
        T0: współczynnik temperatury T₀ (domyślnie T_GAMMA_0)
    """
    args_dict = dict(beta=beta, gamma=gamma, alpha=alpha,
                     use_entropy=use_entropy, s0=s0, T0=T0)

    # Warunki początkowe przy a=1 (N=0)
    phi_today = 1.0 + phi0_dev
    u_today = 0.0
    y0 = [phi_today, u_today]

    N_today = 0.0
    N_ini = LOG_A_INI * np.log(10)  # ln(a_ini)
    N_fin = LOG_A_FIN * np.log(10)  # ln(a_fin)

    # Punkty wyjściowe
    N_back = np.linspace(N_today, N_ini, 4000)
    N_fwd  = np.linspace(N_today, N_fin, 1000)

    def rhs(N, y):
        return ode_system(N, y, args_dict)

    if verbose:
        print("  Całkowanie wstecz (epoka promieniowania)...")
    sol_back = solve_ivp(
        rhs, [N_today, N_ini], y0,
        t_eval=N_back, method='RK45',
        rtol=1e-8, atol=1e-10,
        dense_output=False
    )

    if verbose:
        print("  Całkowanie do przodu (przyszłość)...")
    sol_fwd = solve_ivp(
        rhs, [N_today, N_fin], y0,
        t_eval=N_fwd, method='RK45',
        rtol=1e-8, atol=1e-10,
        dense_output=False
    )

    # Łączymy: wstecz (odwrócony) + do przodu
    N_all = np.concatenate([sol_back.t[::-1], sol_fwd.t[1:]])
    phi_all = np.concatenate([sol_back.y[0][::-1], sol_fwd.y[0][1:]])
    dphi_all = np.concatenate([sol_back.y[1][::-1], sol_fwd.y[1][1:]])
    a_all = np.exp(N_all)
    z_all = 1.0 / a_all - 1.0

    return dict(
        N=N_all, a=a_all, z=z_all,
        phi=phi_all, dphi_dN=dphi_all,
        success_back=sol_back.success,
        success_fwd=sol_fwd.success,
    )


def compute_observables(result, gamma=GAMMA_TGP, beta=BETA_TGP, alpha=ALPHA_TGP,
                        use_entropy=False, s0=S_GAMMA_0, T0=T_GAMMA_0):
    """
    Oblicza obserwabla pochodne:
      - δG/G₀ = 1/φ − 1             [zmiana G]
      - w_de(z)                     [równanie stanu ciemnej energii]
      - H(z)/H₀                     [parametr Hubble'a]
      - S_Γ(z), T_Γ(z)             [entropia i temperatura substratu — gdy use_entropy]
      - δw_entropy(z)               [poprawka do w_de od S_Γ]
    """
    phi = result['phi']
    dphi = result['dphi_dN']
    a = result['a']
    z = result['z']

    # Zmiana G
    delta_G_over_G = 1.0 / np.clip(phi, 1e-10, None) - 1.0

    # Parametr Hubble'a
    H2 = np.array([
        H_squared(a[i], phi[i], dphi[i], gamma, beta, alpha)
        for i in range(len(a))
    ])
    H_over_H0 = np.sqrt(np.clip(H2, 0, None))

    # Równanie stanu ciemnej energii (bez entropii)
    V_phi = beta * phi**2 - gamma * phi**3
    kin = alpha / 2.0 * H2 * dphi**2
    rho_de = kin + V_phi
    p_de   = kin - V_phi
    w_de = np.where(np.abs(rho_de) > 1e-20, p_de / rho_de, -1.0)

    obs = dict(
        z=z, a=a,
        delta_G=delta_G_over_G,
        H_over_H0=H_over_H0,
        w_de=w_de,
        rho_de=rho_de,
        phi=phi,
    )

    # Opcjonalne obserwabla entropii substratu S_Γ
    if use_entropy:
        S_Gamma = entropy_substrate_arr(phi, s0)
        T_Gamma = substrate_temperature_arr(H2, phi, T0)

        # Gęstość energii substratu: ρ_Γ = −T_Γ·S_Γ (swobodna energia)
        rho_substrate = -T_Gamma * S_Gamma   # ujemna → ciemna energia-like

        # Poprawka do ciśnienia: p_Γ = −ρ_Γ + T_Γ·φ·dS_Γ/dφ (prawa Gibbsa)
        # dS_Γ/dφ = s₀·(1 − 1/φ)
        dS_dphi = s0 * (1.0 - 1.0 / np.maximum(phi, 1e-12))
        p_substrate = -rho_substrate + T_Gamma * phi * dS_dphi

        # Efektywne w_de z poprawką S_Γ
        rho_de_total = rho_de + rho_substrate
        p_de_total   = p_de   + p_substrate
        w_de_eff = np.where(
            np.abs(rho_de_total) > 1e-20,
            p_de_total / rho_de_total,
            -1.0
        )

        # δw od entropii
        delta_w_entropy = w_de_eff - w_de

        obs.update(dict(
            S_Gamma=S_Gamma,
            T_Gamma=T_Gamma,
            rho_substrate=rho_substrate,
            p_substrate=p_substrate,
            w_de_eff=w_de_eff,
            delta_w_entropy=delta_w_entropy,
        ))

    return obs


# ──────────────────────────────────────────────────────────────────────────────
#  Raporty diagnostyczne
# ──────────────────────────────────────────────────────────────────────────────

def print_diagnostics(result, obs):
    """Drukuje kluczowe wartości z ewolucji."""
    z = obs['z']
    phi = obs['phi']
    dG = obs['delta_G']
    w_de = obs['w_de']

    # Wartości przy kluczowych epochach
    epochs = {
        'Dziś (z=0)': 0.0,
        'z=0.5': 0.5,
        'z=1': 1.0,
        'z=2': 2.0,
        'z=10 (rekombinacja)': 10.0,
        'z=1000 (CMB)': 1000.0,
        'z=10000': 10000.0,
    }

    interp_phi = interp1d(z[::-1], phi[::-1], bounds_error=False, fill_value='extrapolate')
    interp_dG  = interp1d(z[::-1], dG[::-1],  bounds_error=False, fill_value='extrapolate')
    interp_wde = interp1d(z[::-1], w_de[::-1], bounds_error=False, fill_value='extrapolate')

    print("\n" + "=" * 65)
    print("  Ewolucja kosmologiczna φ_bg(z) — TGP (O22)")
    print("=" * 65)
    print(f"  {'Epoka':<30} {'φ':>10} {'δG/G₀':>12} {'w_de':>10}")
    print("  " + "-" * 63)
    for label, z_val in epochs.items():
        try:
            phi_v = float(interp_phi(z_val))
            dG_v  = float(interp_dG(z_val))
            wde_v = float(interp_wde(z_val))
            print(f"  {label:<30} {phi_v:>10.6f} {dG_v:>12.2e} {wde_v:>10.4f}")
        except Exception:
            print(f"  {label:<30} {'N/A':>10}")

    print("=" * 65)

    # Diagnostyka entropii substratu (jeśli obliczone)
    if 'S_Gamma' in obs:
        S_G = obs['S_Gamma']
        T_G = obs['T_Gamma']
        dw  = obs['delta_w_entropy']
        rho_S = obs['rho_substrate']

        interp_SG   = interp1d(z[::-1], S_G[::-1],  bounds_error=False, fill_value='extrapolate')
        interp_TG   = interp1d(z[::-1], T_G[::-1],  bounds_error=False, fill_value='extrapolate')
        interp_dw   = interp1d(z[::-1], dw[::-1],   bounds_error=False, fill_value='extrapolate')
        interp_rhoS = interp1d(z[::-1], rho_S[::-1], bounds_error=False, fill_value='extrapolate')

        print("\n" + "─" * 65)
        print("  Entropia substratu S_Γ (N0-7 rozszerzenie O22)")
        print("─" * 65)
        print(f"  {'Epoka':<28} {'S_Γ':>10} {'T_Γ':>10} {'δw_de':>10} {'ρ_Γ':>12}")
        print("  " + "─" * 63)
        for label, z_val in epochs.items():
            try:
                SG_v  = float(interp_SG(z_val))
                TG_v  = float(interp_TG(z_val))
                dw_v  = float(interp_dw(z_val))
                rS_v  = float(interp_rhoS(z_val))
                print(f"  {label:<28} {SG_v:>10.4e} {TG_v:>10.4e} {dw_v:>10.4e} {rS_v:>12.4e}")
            except Exception:
                print(f"  {label:<28} {'N/A':>10}")
        print("─" * 65)

        # Sprawdź warunek stabilności N0-7: S_Γ(dziś) ≈ 0 (minimalne odchylenie od próżni)
        SG_today = float(interp_SG(0.0))
        TG_today = float(interp_TG(0.0))
        dw_today = float(interp_dw(0.0))
        ok_SG = SG_today < 1e-2
        sym = "✓" if ok_SG else "✗"
        print(f"\n  {sym}  S_Γ(dziś) = {SG_today:.4e} (powinno być ≪ 1)")
        print(f"  ✓  T_Γ(dziś) = {TG_today:.4e} (temperatura substratu dziś)")
        print(f"  ✓  δw_de(dziś) = {dw_today:.4e} (korekta równania stanu od S_Γ)")
        print()

    # Constraint BBN: δG/G₀ < 0.1 przy z ~ 10⁹ (BBN)
    # Sprawdzamy przy z=10⁴ jako proxy
    try:
        dG_bbn = float(interp_dG(1e4))
        ok_bbn = abs(dG_bbn) < 0.1
        sym = "✓ PASS" if ok_bbn else "✗ FAIL"
        print(f"\n  {sym}  Ograniczenie BBN: |δG/G₀| < 0.1 przy z≈10⁴")
        print(f"         δG/G₀ = {dG_bbn:.4e}")
    except Exception:
        pass

    # Constraint CMB: δG/G₀ < 0.05 przy z=1000
    try:
        dG_cmb = float(interp_dG(1000))
        ok_cmb = abs(dG_cmb) < 0.05
        sym = "✓ PASS" if ok_cmb else "✗ FAIL"
        print(f"  {sym}  Ograniczenie CMB: |δG/G₀| < 0.05 przy z=1000")
        print(f"         δG/G₀ = {dG_cmb:.4e}")
    except Exception:
        pass

    print()


# ──────────────────────────────────────────────────────────────────────────────
#  Wykresy
# ──────────────────────────────────────────────────────────────────────────────

def make_plots(obs, output_dir="plots"):
    """Generuje wykresy ewolucji kosmologicznej (z opcjonalnymi wykresami S_Γ)."""
    import os
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib niedostępny — pomijam wykresy.")
        return

    os.makedirs(output_dir, exist_ok=True)

    z = obs['z']
    # Ograniczamy do z > 0 dla osi log
    mask = (z > 0.001) & (z < 2e4)
    z_m = z[mask]
    phi_m = obs['phi'][mask]
    dG_m = obs['delta_G'][mask]
    w_m = obs['w_de'][mask]
    H_m = obs['H_over_H0'][mask]

    has_entropy = 'S_Gamma' in obs
    nrows = 3 if has_entropy else 2

    fig, axes = plt.subplots(nrows, 2, figsize=(13, 5 * nrows))
    if nrows == 2:
        axes = axes.reshape(2, 2)
    fig.suptitle('TGP — Ewolucja kosmologiczna φ_bg(z)  [Problem O22]',
                 fontsize=13)

    # (a) φ(z)
    ax = axes[0, 0]
    ax.semilogx(z_m[::-1], phi_m[::-1], color='royalblue', lw=2)
    ax.axhline(1.0, color='gray', ls='--', lw=0.8, label='φ=1 (próżnia)')
    ax.set_xlabel('z (przesunięcie ku czerwieni)')
    ax.set_ylabel('φ = Φ_bg / Φ₀')
    ax.set_title('(a) Ewolucja pola φ(z)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) δG/G₀
    ax = axes[0, 1]
    ax.semilogx(z_m[::-1], dG_m[::-1] * 100, color='firebrick', lw=2)
    ax.axhline(0, color='gray', ls='--', lw=0.8)
    ax.axhline(10.0, color='orange', ls=':', lw=1.2, label='Limit BBN (10%)')
    ax.axhline(5.0, color='gold', ls=':', lw=1.2, label='Limit CMB (5%)')
    ax.set_xlabel('z')
    ax.set_ylabel('δG/G₀  [%]')
    ax.set_title('(b) Zmiana stałej grawitacji G(z)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (c) w_de(z)
    ax = axes[1, 0]
    w_clipped = np.clip(w_m, -3, 1)
    ax.semilogx(z_m[::-1], w_clipped[::-1], color='seagreen', lw=2)
    ax.axhline(-1.0, color='gray', ls='--', lw=0.8, label='w=−1 (Λ)')
    ax.set_xlabel('z')
    ax.set_ylabel('w_de(z)')
    ax.set_title('(c) Równanie stanu ciemnej energii')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-2, 0.5)

    # (d) H(z)/H₀
    ax = axes[1, 1]
    # Porównanie TGP vs ΛCDM
    H_lcdm = np.sqrt(OMEGA_M / (1 + z_m)**3 + OMEGA_R / (1 + z_m)**4 + OMEGA_L)
    ax.loglog(z_m[::-1], H_m[::-1], color='mediumpurple', lw=2, label='TGP')
    ax.loglog(z_m[::-1], H_lcdm[::-1], color='gray', lw=1.5,
              ls='--', label='ΛCDM')
    ax.set_xlabel('z')
    ax.set_ylabel('H(z) / H₀')
    ax.set_title('(d) H(z): TGP vs ΛCDM')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # (e) Entropia substratu S_Γ(z) — gdy dostępna
    if has_entropy:
        S_G_m = obs['S_Gamma'][mask]
        T_G_m = obs['T_Gamma'][mask]
        dw_m_ent = obs['delta_w_entropy'][mask]
        w_eff_m = obs['w_de_eff'][mask]

        ax = axes[2, 0]
        ax.semilogx(z_m[::-1], S_G_m[::-1],
                    color='royalblue', lw=2, label='S_Γ(φ)')
        ax2 = ax.twinx()
        ax2.semilogx(z_m[::-1], T_G_m[::-1],
                     color='firebrick', lw=1.5, ls='--', label='T_Γ(H,φ)')
        ax.set_xlabel('z')
        ax.set_ylabel('S_Γ (entropia substratu)', color='royalblue')
        ax2.set_ylabel('T_Γ (temperatura substratu)', color='firebrick')
        ax.set_title('(e) Entropia S_Γ i temperatura T_Γ substratu')
        ax.grid(True, alpha=0.3)
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper left')

        ax = axes[2, 1]
        ax.semilogx(z_m[::-1], np.clip(w_eff_m[::-1], -3, 0.5),
                    color='seagreen', lw=2, label='w_de z S_Γ')
        ax.semilogx(z_m[::-1], np.clip(w_m[::-1], -3, 0.5),
                    color='gray', lw=1.5, ls='--', label='w_de bez S_Γ')
        ax.axhline(-1.0, color='black', ls=':', lw=0.8, label='w=−1')
        ax.set_xlabel('z')
        ax.set_ylabel('w_de(z)')
        ax.set_title('(f) Równanie stanu: z/bez entropii S_Γ')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-2.5, 0.5)

    plt.tight_layout()
    suffix = '_entropy' if has_entropy else ''
    out_path = os.path.join(output_dir, f'cosmological_phi_evolution{suffix}.png')
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    print(f"  [PLOT] Zapisano: {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
#  Skan parametryczny (wrażliwość na φ₀_dev)
# ──────────────────────────────────────────────────────────────────────────────

def parametric_scan(plot=False, use_entropy=False, s0=S_GAMMA_0, T0=T_GAMMA_0):
    """
    Skan amplitudy odchyłki φ₀_dev ∈ {10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³, 10⁻²}.
    Pokazuje, jak silnie ewolucja zależy od warunku brzegowego dziś.
    Gdy use_entropy=True: porównanie z/bez entropii substratu S_Γ.
    """
    devs = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    ent_label = f"[s₀={s0:.2e}, T₀={T0:.2e}]" if use_entropy else "[bez S_Γ]"
    print(f"\n  [SKAN] Wrażliwość na warunek brzegowy φ₀_dev {ent_label}:")

    if use_entropy:
        print(f"  {'φ₀_dev':>10} {'δφ(z=1000)':>14} {'δG/G₀(z=1000)':>16} {'δw(z=0)':>10}")
    else:
        print(f"  {'φ₀_dev':>10} {'δφ(z=1000)':>14} {'δG/G₀(z=1000)':>16}")
    print("  " + "─" * 55)

    scan_results = []
    for dev in devs:
        try:
            res = solve_phi_evolution(phi0_dev=dev, verbose=False,
                                      use_entropy=use_entropy, s0=s0, T0=T0)
            obs = compute_observables(res, use_entropy=use_entropy, s0=s0, T0=T0)

            z = obs['z']
            phi = obs['phi']
            dG = obs['delta_G']

            interp_phi = interp1d(z[::-1], phi[::-1],
                                  bounds_error=False, fill_value='extrapolate')
            interp_dG  = interp1d(z[::-1], dG[::-1],
                                  bounds_error=False, fill_value='extrapolate')

            dphi_cmb = float(interp_phi(1000)) - 1.0
            dG_cmb = float(interp_dG(1000))

            if use_entropy and 'delta_w_entropy' in obs:
                interp_dw = interp1d(z[::-1], obs['delta_w_entropy'][::-1],
                                     bounds_error=False, fill_value='extrapolate')
                dw_today = float(interp_dw(0.0))
                print(f"  {dev:>10.0e} {dphi_cmb:>14.4e} {dG_cmb:>16.4e} {dw_today:>10.4e}")
            else:
                print(f"  {dev:>10.0e} {dphi_cmb:>14.4e} {dG_cmb:>16.4e}")

            scan_results.append((dev, dphi_cmb, dG_cmb, res, obs))
        except Exception as e:
            print(f"  {dev:>10.0e}  BŁĄD: {e}")

    return scan_results


def entropy_scan(phi0_dev=PHI0_REL_DEVIATION, s0_range=None, T0_range=None,
                 verbose=False):
    """
    Skan parametrów entropii substratu {s₀, T₀} → efekt na δw_de(dziś) i S_Γ.
    Cel: zbadanie wrażliwości O22 na nieznane parametry N0-7.
    """
    if s0_range is None:
        s0_range = [0.001, 0.005, 0.01, 0.05, 0.1]
    if T0_range is None:
        T0_range = [0.05, 0.1, 0.2, 0.5]

    print("\n" + "─" * 68)
    print("  [SKAN ENTROPII] s₀ × T₀ → δw_de(z=0), S_Γ(z=0), T_Γ(z=0)")
    print("─" * 68)
    print(f"  {'s₀':>8}  {'T₀':>8}  {'S_Γ(z=0)':>12}  {'T_Γ(z=0)':>12}  {'δw(z=0)':>12}")
    print("  " + "─" * 60)

    results_grid = {}
    for s0 in s0_range:
        for T0 in T0_range:
            try:
                res = solve_phi_evolution(phi0_dev=phi0_dev, verbose=False,
                                          use_entropy=True, s0=s0, T0=T0)
                obs = compute_observables(res, use_entropy=True, s0=s0, T0=T0)

                z = obs['z']
                SG = obs['S_Gamma']
                TG = obs['T_Gamma']
                dw = obs['delta_w_entropy']

                interp_SG = interp1d(z[::-1], SG[::-1],
                                     bounds_error=False, fill_value='extrapolate')
                interp_TG = interp1d(z[::-1], TG[::-1],
                                     bounds_error=False, fill_value='extrapolate')
                interp_dw = interp1d(z[::-1], dw[::-1],
                                     bounds_error=False, fill_value='extrapolate')

                SG_0 = float(interp_SG(0.0))
                TG_0 = float(interp_TG(0.0))
                dw_0 = float(interp_dw(0.0))

                print(f"  {s0:>8.3e}  {T0:>8.3e}  {SG_0:>12.4e}  {TG_0:>12.4e}  {dw_0:>12.4e}")
                results_grid[(s0, T0)] = (SG_0, TG_0, dw_0)

            except Exception as e:
                print(f"  {s0:>8.3e}  {T0:>8.3e}  BŁĄD: {e}")

    print("─" * 68)
    print("  UWAGA: δw_de(z=0) < 0.01 → entropia nie zmienia obserwacji dziś")
    print("  Dla δw(z=0) > 0.1 → możliwe napięcia z DESI/Euclid w_de pomiarami")

    return results_grid


# ──────────────────────────────────────────────────────────────────────────────
#  Punkt wejścia
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Ewolucja kosmologiczna φ_bg(z) — TGP O22 (z entropią substratu S_Γ)"
    )
    parser.add_argument('--plot', action='store_true',
                        help='Generuj wykresy (wymaga matplotlib)')
    parser.add_argument('--gamma', type=float, default=GAMMA_TGP,
                        help=f'Parametr γ (domyślnie {GAMMA_TGP})')
    parser.add_argument('--phi0-dev', type=float, default=PHI0_REL_DEVIATION,
                        help=f'Odchyłka φ₀ dziś (domyślnie {PHI0_REL_DEVIATION})')
    parser.add_argument('--scan', action='store_true',
                        help='Wykonaj skan wrażliwości na φ₀_dev')
    # Nowe opcje: entropia substratu S_Γ
    parser.add_argument('--entropy', action='store_true',
                        help='Włącz człon entropii substratu S_Γ (N0-7, O22 pełna)')
    parser.add_argument('--s0', type=float, default=S_GAMMA_0,
                        help=f'Gęstość entropii s₀ (domyślnie {S_GAMMA_0})')
    parser.add_argument('--T0', type=float, default=T_GAMMA_0,
                        help=f'Temperatura substratu T₀/H₀ (domyślnie {T_GAMMA_0})')
    parser.add_argument('--entropy-scan', action='store_true',
                        help='Skan parametrów s₀×T₀ (wymaga --entropy)')
    args = parser.parse_args()

    gamma = args.gamma
    beta  = gamma  # N0-5: β = γ
    use_ent = args.entropy
    s0_val  = args.s0
    T0_val  = args.T0

    print("=" * 65)
    print("  TGP — Problem O22: Ewolucja kosmologiczna φ_bg(z)")
    print(f"  γ={gamma:.3g}, β={beta:.3g} (N0-5), φ₀_dev={args.phi0_dev:.2e}")
    if use_ent:
        print(f"  Entropia S_Γ: s₀={s0_val:.3e}, T₀={T0_val:.3e} [N0-7 aktywne]")
    else:
        print(f"  Entropia S_Γ: WYŁĄCZONA (użyj --entropy by włączyć)")
    print("=" * 65)

    # Główna symulacja
    result = solve_phi_evolution(
        gamma=gamma, beta=beta, phi0_dev=args.phi0_dev,
        use_entropy=use_ent, s0=s0_val, T0=T0_val
    )

    if not result['success_back']:
        print("  [WARN] Całkowanie wstecz nie osiągnęło zbieżności!")
    if not result['success_fwd']:
        print("  [WARN] Całkowanie do przodu nie osiągnęło zbieżności!")

    obs = compute_observables(result, gamma=gamma, beta=beta,
                               use_entropy=use_ent, s0=s0_val, T0=T0_val)
    print_diagnostics(result, obs)

    if args.plot:
        make_plots(obs)

    if args.scan:
        parametric_scan(plot=args.plot,
                        use_entropy=use_ent, s0=s0_val, T0=T0_val)

    if args.entropy_scan:
        entropy_scan(phi0_dev=args.phi0_dev)

    if use_ent:
        print("""
  ─────────────────────────────────────────────────────────
  STATUS O22 (N0-7 CZĘŚCIOWO AKTYWNY):
  Entropia substratu S_Γ uwzględniona w tym przebiegu.
  Parametry s₀, T₀ są wolne — wymagają modelu H_Γ.
  Brakujące elementy:
  (1) Mikrofizyczny model s₀ z hamiltonianu H_Γ substratu
  (2) Pełne równanie Boltzmanna z modyfikowanym c(Φ)
  (3) Perturbacje skalarne CMB z δS_Γ
  (4) Renormalizacja T_Γ przez kwantowe poprawki pętelkowe
  ─────────────────────────────────────────────────────────
""")
    else:
        print("""
  ─────────────────────────────────────────────────────────
  UWAGA O22 (OTWARTY): Wyniki mają charakter jakościowy.
  Użyj --entropy by włączyć człon S_Γ (N0-7 rozszerzenie).
  Brakujące elementy do pełnego rozwiązania:
  (1) Entropia substratu S_Γ i jej sprzężenie z φ [dodane!]
  (2) Pełne równanie Boltzmanna z modyfikowanym c(Φ)
  (3) Perturbacje skalarne CMB z dynamicznym G(Φ)
  (4) Spójność z sekcją sek05 (ciemna energia TGP)
  ─────────────────────────────────────────────────────────
""")

    sys.exit(0)
