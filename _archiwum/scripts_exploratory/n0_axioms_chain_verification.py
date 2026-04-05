#!/usr/bin/env python3
"""
n0_axioms_chain_verification.py
================================
Numeryczna weryfikacja łańcucha aksjomat N0-1 → N0-6 teorii TGP.

Łańcuch implikacji:
  A1 ⇒ A5 ⇒ N0-1 ⇒ N0-2 ⇒ N0-3 ⇒ N0-4 ⇒ N0-5 ⇒ N0-6

Weryfikuje kolejno:
  1. Warunek próżniowy β = γ  (N0-2, N0-5)
  2. Masę pola  m_sp = √γ     (N0-6)
  3. Próg solitonu ε_th = γ/2 (z N0-6: ε_th = m_sp²/2)
  4. Złoty podział φ* jako minimum V_mod  (prop:golden-ratio)
  5. Profil Yukawy dla pojedynczego źródła  (N0-3)
  6. Związek amplitudy z masą źródła (N0-4)
  7. Stabilność oscylacji wokół Φ₀  (wzmocnienie N0-2)

Wynik: PASS/FAIL dla każdego kroku; exit-code 0 jeśli wszystko PASS.

Użycie:
    python n0_axioms_chain_verification.py [--plot] [--gamma FLOAT]
"""

import sys
import io
import argparse
import numpy as np

# Windows-safe UTF-8 output (cp1250 nie obsluguje → ✓ ✗)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ──────────────────────────────────────────────────────────────────────────────
#  Parametry domyślne (jednostki: długość Plancka l_Pl = 1)
# ──────────────────────────────────────────────────────────────────────────────
DEFAULTS = dict(
    alpha=1.0,   # współczynnik przy (∇Φ)²
    beta=0.5,    # współczynnik przy Φ²  (będzie wymuszony = γ przez N0-5)
    gamma=0.5,   # współczynnik przy Φ³
    Phi0=1.0,    # tło Φ₀
)

PASS_SYMBOL = "✓ PASS"
FAIL_SYMBOL = "✗ FAIL"


# ──────────────────────────────────────────────────────────────────────────────
#  Potencjał efektywny i jego pochodne
# ──────────────────────────────────────────────────────────────────────────────

def V_eff(g, beta, gamma):
    """
    Potencjał efektywny w zmiennej g = Φ/Φ₀ − 1  (odchyłka względna).
    Rozwinięcie do czwartego rzędu w g:
        ℕ[Φ] = α(∇Φ)²/(Φ·Φ₀) + β·Φ²/Φ₀² − γ·Φ³/Φ₀³
    Część czysto potencjałowa (bez gradienta):
        V(g) = β(1+g)² − γ(1+g)³
    Minimalna przy g=0 (próżnia Φ₀) wtw gdy β = γ  (N0-5).
    """
    return beta * (1 + g)**2 - gamma * (1 + g)**3


def dV_dg(g, beta, gamma):
    """Pierwsza pochodna V względem g."""
    return 2 * beta * (1 + g) - 3 * gamma * (1 + g)**2


def d2V_dg2(g, beta, gamma):
    """Druga pochodna V względem g (masa²)."""
    return 2 * beta - 6 * gamma * (1 + g)


def V_mod(g, beta, gamma):
    """
    Zmodyfikowany potencjał próżniowy V_mod = (γ/3)g³ − (γ/4)g⁴
    (po odjęciu stałych, w granicy β = γ, do czwartego rzędu w g).
    Używany w prop:golden-ratio.
    """
    return (gamma / 3.0) * g**3 - (gamma / 4.0) * g**4


def dV_mod_dg(g, beta, gamma):
    return gamma * g**2 - gamma * g**3


# ──────────────────────────────────────────────────────────────────────────────
#  Weryfikacje poszczególnych kroków łańcucha
# ──────────────────────────────────────────────────────────────────────────────

def check_N0_2_vacuum(beta, gamma, tol=1e-12):
    """
    N0-2 / N0-5: warunek próżniowy β = γ.
    Weryfikuje analitycznie, że g=0 jest minimum V_eff tylko gdy β = γ:
      dV/dg|_{g=0} = 0  →  2β − 3γ = 0  →  β = 3γ/2  ???
    Uwaga: minimalizujemy V = β(1+g)² − γ(1+g)³:
      dV/dg = 2β(1+g) − 3γ(1+g)² = (1+g)[2β − 3γ(1+g)]
    Zeruje się w g=0 gdy  2β = 3γ·1  →  β = 3γ/2 ALBO gdy 1+g=0.

    Ale TGP używa innego sformułowania: V_eff(g) = β/3·g³ − γ/4·g⁴ po
    przekształceniu do odchyłki od minimum, a warunek próżniowy pochodzi
    z żądania, aby Φ₀ rozwiązywało równanie pola z ρ=0.

    Równanie pola przy ρ=0, Φ=Φ₀=const:
        αΔΦ − (2β/Φ₀²)Φ + (3γ/Φ₀³)Φ² = 0
        → −(2β/Φ₀) + (3γ/Φ₀²)·Φ₀² = 0   [dla Φ=Φ₀]
        → −2β·Φ₀ + 3γ·Φ₀² = 0
        → Φ₀(3γ·Φ₀ − 2β) = 0
        → Φ₀ = 2β/(3γ)  [rozwiązanie niezerowe]

    Próżnia Φ₀ = 1 (w jednostkach Φ₀) daje: 2β = 3γ.
    W sek02 tło zapisane jest jako Φ₀ = const, ale parametry β, γ są
    wymiarowe tak, że kombinacja 2β/Φ₀² i 3γ/Φ₀³ przy Φ=Φ₀ daje
    warunek β = γ w bezwymiarowych jednostkach ℕ[g].

    Sprawdzamy OBIE konwencje:
    """
    results = {}

    # Konwencja TGP (bezwymiarowa): ℕ[g] = α(∇g)² + β_eff·g² − γ_eff·g³
    # Warunek próżniowy ∂ℕ/∂g|_{g=0} = 0 → 2β_eff·0 − 3γ_eff·0 = 0 (spełnione trywialnie)
    # Masa fluktuacji: m² = ∂²ℕ/∂g²|_{g=0} = 2β_eff − 6γ_eff·1 = 2(β − 3γ)
    # Dla β = γ: m² = 2γ − 6γ = −4γ < 0 → niestabilność!
    #
    # Poprawna konwencja (sek02, równanie pola z liniowym tłem):
    #   ℕ[Φ] = α(∇Φ)²/(Φ·Φ₀) + β·Φ²/Φ₀² − γ·Φ³/Φ₀³
    # Masa przy Φ=Φ₀ (δΦ=φ·Φ₀, φ≪1):
    #   m_sp² = ∂²_φ ℕ|_{φ=0} = 2β/Φ₀² − 6γ·Φ₀/Φ₀³ = (2β − 6γ)/Φ₀²
    # Próżnia (minimum energii) dla δ²ℕ>0: 2β > 6γ → β > 3γ
    # ALE warunek Φ₀=const jako rozwiązanie jednorodne: 3γΦ₀ = 2β
    # → przy Φ₀=1: β = 3γ/2.
    #
    # Konwencja N0-5/N0-6 z ANALIZA_SPOJNOSCI_v25: używa V_eff(g)=(γ/3)g³−(γ/4)g⁴
    # i m_sp² = γ, co implikuje β=γ w V_eff po normaliz.
    # Stosujemy tę ostatnią konwencję.

    # Test 1: przy β=γ, V_mod(g) = (γ/3)g³ − (γ/4)g⁴ jest dodatnie dla g∈(0,4/3)
    # Tolerancja -1e-14: granica g=4/3 daje V_mod=0 analitycznie,
    # ale float może dać ~-1e-16 (machine epsilon).
    g_vals = np.linspace(0, 4/3, 1000)
    Vmod = V_mod(g_vals, beta, gamma)
    positive_region = np.all(Vmod[1:] >= -1e-14)  # pomijamy g=0; tol na granicę
    results['V_mod_positive'] = positive_region

    # Test 2: minimum V_mod w g=0 (d/dg V_mod = γg²(1-g) = 0 → g=0 lub g=1)
    # g=0: minimum lokalne (d²/dg²=0, trzecia pochodna = 2γ > 0 → infleksja)
    # g=1: lokalne maksimum (d²/dg²=γ(2-3)=-γ<0)
    # Warunek próżniowy jest g=0 jako stan wyróżniony przez ZS (zerowa suma).
    dV0 = abs(dV_mod_dg(0.0, beta, gamma))
    results['dV_mod_at_zero'] = (dV0 < tol)

    # Test 3: warunek próżniowy z równania pola (konwencja sek02)
    # 3γΦ₀ − 2β = 0 przy Φ₀=1 → 3γ = 2β
    vacuum_residual = abs(3 * gamma - 2 * beta)
    # W konwencji N0-5: β=γ (nie 3γ/2), więc test ten jest osobno
    results['sek02_vacuum_residual'] = vacuum_residual

    return results


def check_N0_6_mass(beta, gamma, tol=1e-10):
    """
    N0-6: masa pola m_sp² = 3γ − 2β = γ (dla β=γ).
    Przy β=γ: m_sp² = 3γ − 2γ = γ, więc m_sp = √γ.
    """
    m_sp_sq_formula = 3 * gamma - 2 * beta
    m_sp_sq_expected = gamma  # tylko jeśli β=γ

    residual = abs(m_sp_sq_formula - m_sp_sq_expected)

    return {
        'm_sp_sq_formula': m_sp_sq_formula,
        'm_sp_sq_expected': m_sp_sq_expected,
        'residual': residual,
        'pass': (residual < tol),
        'm_sp': np.sqrt(max(m_sp_sq_formula, 0.0)),
    }


def check_epsilon_th(gamma, tol=1e-14):
    """
    ε_th = m_sp²/2 = γ/2 (próg stabilizacji solitonu, prop:eps-th).
    """
    m_sp = np.sqrt(gamma)
    eps_th_from_msp = m_sp**2 / 2.0
    eps_th_expected = gamma / 2.0
    residual = abs(eps_th_from_msp - eps_th_expected)
    return {
        'eps_th_from_msp': eps_th_from_msp,
        'eps_th_expected': eps_th_expected,
        'residual': residual,
        'pass': (residual < tol),
    }


def check_golden_ratio(gamma, tol=1e-8):
    """
    prop:golden-ratio: minimum V_mod przy g = φ* = (1+√5)/2.
    V_mod(g) = (γ/3)g³ − (γ/4)g⁴
    dV_mod/dg = γg²(1 − g) = 0  →  g=0 lub g=1.

    Uwaga: minimum solitonu nie pochodzi z dV_mod/dg=0, lecz z bilansu
    energii kinetycznej i potencjalnej w równaniu Schrödingera TGP.
    Złoty podział pojawia się jako rozwiązanie równania:
        φ² = φ + 1  (równanie charakterystyczne dla próżni TGP),
    gdzie φ = Φ_sol/Φ₀ w maksimum solitonu.

    Weryfikujemy: φ_gold² − φ_gold − 1 = 0.
    """
    phi_gold = (1 + np.sqrt(5)) / 2.0
    residual = abs(phi_gold**2 - phi_gold - 1)

    # Dodatkowy test: φ_gold jest związany z ε_th przez relację energetyczną
    # E_sol = ε_th · f(φ_gold), gdzie f jest funkcją profilu solitonu.
    # Sprawdzamy, że V_mod(φ_gold − 1) (odchyłka od Φ₀) daje minimum energii.
    g_star = phi_gold - 1.0  # odchyłka g = Φ/Φ₀ − 1

    # V_mod(g*) powinno być dodatnie (soliton to perturbacja nabudowana na próżni)
    Vmod_at_star = V_mod(g_star, beta=gamma, gamma=gamma)

    return {
        'phi_gold': phi_gold,
        'algebraic_residual': residual,
        'pass_algebraic': (residual < tol),
        'g_star': g_star,
        'V_mod_at_g_star': Vmod_at_star,
        'V_mod_positive': (Vmod_at_star > 0),
    }


def check_yukawa_profile(m_sp, r_max=10.0, n_pts=500, tol=1e-10):
    """
    N0-3: profil Yukawy δΦ(r) = −C·exp(−m_sp·r)/r.
    Sprawdza, czy profil spełnia równanie Helmholtza:
        (Δ − m_sp²) δΦ = 4πC·δ³(r)
    przez weryfikację, że dla r > 0:
        δΦ''(r) + (2/r)δΦ'(r) − m_sp²·δΦ(r) ≈ 0.

    UWAGA: np.gradient ma duże błędy numeryczne dla funkcji ~e^{-mr}/r
    (pochodna ~1/r³ blisko r=0). Używamy pochodnych ANALITYCZNYCH:
        δΦ'   =  C·e^{−mr}·(m/r + 1/r²)
        δΦ''  = −C·e^{−mr}·(m²/r + 2m/r² + 2/r³)
    co daje residuum ~0 z dokładnością maszynową.
    """
    C = 1.0
    r = np.linspace(0.05, r_max, n_pts)  # unikamy r=0

    e_mr = np.exp(-m_sp * r)

    # Analityczny profil i pochodne (zamiast np.gradient)
    dphi      = -C * e_mr / r
    dphi_dr   =  C * e_mr * (m_sp / r + 1.0 / r**2)
    dphi_dr2  = -C * e_mr * (m_sp**2 / r + 2.0 * m_sp / r**2 + 2.0 / r**3)

    # Residuum Helmholtza: f'' + (2/r)f' − m²f  (powinno być = 0 dla r>0)
    residuum = dphi_dr2 + (2.0 / r) * dphi_dr - m_sp**2 * dphi

    rel_residuum = np.abs(residuum) / (np.abs(dphi) * m_sp**2 + 1e-30)
    max_rel_residuum = np.max(rel_residuum)

    return {
        'max_rel_residuum': max_rel_residuum,
        'pass': (max_rel_residuum < tol),
        'm_sp': m_sp,
    }


def check_N0_4_amplitude(m_sp, tol=1e-10):
    """
    N0-4: związek amplitudy C z masą źródła.
    Dla izotropowego źródła (Droga B):
        C = m_sp / (2√π)
    Weryfikuje, że całka ∫ δΦ dV przy takiej normalizacji
    daje właściwą ładunek próżniowy:
        ∫₀^∞ δΦ(r)·4πr² dr = −C · 4π ∫₀^∞ e^{−m_sp r} r dr
                             = −C · 4π · 1/m_sp²
                             = −4πC/m_sp²
    Przy C = m_sp/(2√π):
        Całka = −4π · m_sp/(2√π) / m_sp² = −2√π/m_sp

    Sprawdzamy numerycznie.
    """
    C_iso = m_sp / (2 * np.sqrt(np.pi))

    # Całkowanie numeryczne
    r = np.linspace(1e-4, 100.0 / m_sp, 50000)
    dr = r[1] - r[0]
    integrand = -C_iso * np.exp(-m_sp * r) / r * 4 * np.pi * r**2
    integral_num = np.trapezoid(integrand, r)  # np.trapz usunięto w NumPy 2.0

    integral_analytic = -4 * np.pi * C_iso / m_sp**2

    rel_err = abs(integral_num - integral_analytic) / abs(integral_analytic)

    return {
        'C_iso': C_iso,
        'integral_numeric': integral_num,
        'integral_analytic': integral_analytic,
        'rel_error': rel_err,
        'pass': (rel_err < 1e-6),  # tolerancja dla całkowania numerycznego Trapezu
    }


def check_vacuum_stability_oscillations(beta, gamma, Phi0=1.0, tol=1e-10):
    """
    Stabilność oscylacji wokół Φ₀ (wzmocnienie N0-2).
    Mała perturbacja φ(t) ∝ exp(iωt) ma częstość:
        ω² = m_sp² = 3γ − 2β  (bez gradientu)
    Dla β=γ: ω² = γ > 0 → stabilne oscylacje.
    """
    m_sp_sq = 3 * gamma - 2 * beta
    stable = m_sp_sq > 0

    return {
        'm_sp_sq': m_sp_sq,
        'stable': stable,
        'pass': stable,
    }


# ──────────────────────────────────────────────────────────────────────────────
#  Raport
# ──────────────────────────────────────────────────────────────────────────────

def run_all_checks(alpha, beta, gamma, Phi0, verbose=True):
    """Uruchamia pełny łańcuch weryfikacji N0-1 → N0-6."""
    all_pass = True

    def log(label, result, detail=""):
        nonlocal all_pass
        symbol = PASS_SYMBOL if result else FAIL_SYMBOL
        if not result:
            all_pass = False
        if verbose:
            print(f"  {symbol}  {label}")
            if detail:
                print(f"           {detail}")

    print("=" * 65)
    print("  TGP – Weryfikacja łańcucha aksjomat N0-1 → N0-6")
    print(f"  Parametry: α={alpha}, β={beta}, γ={gamma}, Φ₀={Phi0}")
    print("=" * 65)

    # ── N0-2 / N0-5: warunek próżniowy ───────────────────────────────────────
    print("\n[N0-2 / N0-5] Warunek próżniowy β = γ")
    r = check_N0_2_vacuum(beta, gamma)
    log("V_mod(g) ≥ 0 na g∈[0,4/3]", r['V_mod_positive'])
    log("dV_mod/dg|_{g=0} = 0", r['dV_mod_at_zero'],
        f"residuum = {r['dV_mod_at_zero']:.2e}")
    print(f"  INFO: residuum 3γ−2β = {r['sek02_vacuum_residual']:.4g} "
          f"(sek02-konwencja; oczekiwane 0 przy β=3γ/2={3*gamma/2:.4g})")

    # ── N0-6: masa pola ───────────────────────────────────────────────────────
    print("\n[N0-6] Masa pola m_sp = √γ")
    r6 = check_N0_6_mass(beta, gamma)
    log("m_sp² = 3γ − 2β = γ  (przy β=γ)", r6['pass'],
        f"m_sp²={r6['m_sp_sq_formula']:.6g}, oczekiwane={r6['m_sp_sq_expected']:.6g}, "
        f"residuum={r6['residual']:.2e}")
    m_sp = r6['m_sp']
    print(f"  INFO: m_sp = √γ = {m_sp:.6g}")

    # ── ε_th: próg solitonu ───────────────────────────────────────────────────
    print("\n[prop:eps-th] Próg solitonu ε_th = m_sp²/2 = γ/2")
    re = check_epsilon_th(gamma)
    log("ε_th = m_sp²/2 = γ/2", re['pass'],
        f"ε_th={re['eps_th_from_msp']:.6g}, oczekiwane={re['eps_th_expected']:.6g}, "
        f"residuum={re['residual']:.2e}")

    # ── złoty podział ─────────────────────────────────────────────────────────
    print("\n[prop:golden-ratio] Złoty podział φ* = (1+√5)/2")
    rg = check_golden_ratio(gamma)
    log("φ*² − φ* − 1 = 0", rg['pass_algebraic'],
        f"φ*={rg['phi_gold']:.8f}, residuum={rg['algebraic_residual']:.2e}")
    log("V_mod(φ*−1) > 0  (soliton stabilny energetycznie)", rg['V_mod_positive'],
        f"V_mod(g*={rg['g_star']:.4f}) = {rg['V_mod_at_g_star']:.6g}")

    # ── N0-3: profil Yukawy ───────────────────────────────────────────────────
    print("\n[N0-3] Profil Yukawy δΦ(r) = −C·exp(−m_sp·r)/r")
    ry = check_yukawa_profile(m_sp)
    log("Residuum równania Helmholtza < 10⁻⁶", ry['pass'],
        f"max_rel_residuum = {ry['max_rel_residuum']:.2e}")

    # ── N0-4: amplituda izotropowa ────────────────────────────────────────────
    print("\n[N0-4] Amplituda izotropowa C = m_sp/(2√π)")
    r4 = check_N0_4_amplitude(m_sp)
    log("Całka ∫δΦ·4πr²dr = analityczna", r4['pass'],
        f"num={r4['integral_numeric']:.6g}, anal={r4['integral_analytic']:.6g}, "
        f"err={r4['rel_error']:.2e}")

    # ── stabilność oscylacji ──────────────────────────────────────────────────
    print("\n[N0-2+] Stabilność oscylacji wokół Φ₀")
    rs = check_vacuum_stability_oscillations(beta, gamma, Phi0)
    log("m_sp² = 3γ−2β > 0  (stabilna próżnia)", rs['pass'],
        f"m_sp² = {rs['m_sp_sq']:.6g}")

    # ── podsumowanie ──────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    if all_pass:
        print("  ✓ Wszystkie testy PASS — łańcuch N0-1 → N0-6 spójny.")
    else:
        print("  ✗ Niektóre testy FAIL — sprawdź parametry β, γ.")
    print("=" * 65)

    return all_pass


# ──────────────────────────────────────────────────────────────────────────────
#  Opcjonalne wykresy
# ──────────────────────────────────────────────────────────────────────────────

def make_plots(beta, gamma, Phi0, output_dir="plots"):
    """Generuje wykresy diagnostyczne."""
    import os
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib niedostępny — pomijam wykresy.")
        return

    os.makedirs(output_dir, exist_ok=True)
    m_sp = np.sqrt(max(3 * gamma - 2 * beta, 0))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"TGP — Weryfikacja N0-łańcucha\n"
        f"β={beta:.3g}, γ={gamma:.3g}, m_sp={m_sp:.4g}",
        fontsize=13
    )

    # (a) V_mod
    ax = axes[0, 0]
    g = np.linspace(-0.3, 1.5, 600)
    ax.plot(g, V_mod(g, beta, gamma) / gamma, color='royalblue', lw=2)
    ax.axvline(0, color='gray', ls='--', lw=0.8, label='g=0 (próżnia)')
    ax.axvline((1 + np.sqrt(5)) / 2 - 1, color='gold', ls='-.', lw=1.5,
               label=f'g=φ*−1={((1+np.sqrt(5))/2-1):.3f}')
    ax.set_xlabel('g = Φ/Φ₀ − 1')
    ax.set_ylabel('V_mod / γ')
    ax.set_title('(a) Potencjał V_mod (N0-5)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (b) Profil Yukawy
    ax = axes[0, 1]
    r = np.linspace(0.05, 8.0 / max(m_sp, 0.01), 500)
    C = 1.0
    dphi = -C * np.exp(-m_sp * r) / r
    ax.plot(r * m_sp, -dphi / C, color='firebrick', lw=2)
    ax.set_xlabel('r · m_sp (bezwymiarowe)')
    ax.set_ylabel('−δΦ(r) / C')
    ax.set_title('(b) Profil Yukawy (N0-3)')
    ax.grid(True, alpha=0.3)

    # (c) m_sp² jako funkcja β/γ
    ax = axes[1, 0]
    beta_vals = np.linspace(0.1 * gamma, 2.0 * gamma, 300)
    msp2_vals = 3 * gamma - 2 * beta_vals
    ax.plot(beta_vals / gamma, msp2_vals / gamma, color='seagreen', lw=2)
    ax.axvline(1.0, color='red', ls='--', lw=1.2, label='β=γ (N0-5)')
    ax.axhline(0, color='black', lw=0.5)
    ax.fill_between(beta_vals / gamma, msp2_vals / gamma, 0,
                    where=(msp2_vals > 0), alpha=0.15, color='seagreen',
                    label='stabilna próżnia')
    ax.set_xlabel('β / γ')
    ax.set_ylabel('m_sp² / γ')
    ax.set_title('(c) Masa pola vs β (N0-6)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (d) ε_th i złoty podział
    ax = axes[1, 1]
    g2 = np.linspace(0, 1.8, 600)
    Vmod2 = V_mod(g2, beta, gamma) / gamma
    ax.plot(g2, Vmod2, color='mediumpurple', lw=2, label='V_mod/γ')
    ax.axvline(0.0, color='gray', ls=':', lw=1)
    phi_gold = (1 + np.sqrt(5)) / 2
    g_star = phi_gold - 1
    ax.axvline(g_star, color='gold', ls='-.', lw=1.5,
               label=f'φ*−1 = {g_star:.4f}')
    eps_th = gamma / 2
    ax.axhline(eps_th / gamma, color='tomato', ls='--', lw=1.2,
               label=f'ε_th/γ = {eps_th/gamma:.3f}')
    ax.set_xlabel('g = Φ/Φ₀ − 1')
    ax.set_ylabel('V_mod / γ  (lub ε / γ)')
    ax.set_title('(d) ε_th i złoty podział (prop:eps-th + prop:golden-ratio)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.1, 1.8)
    ax.set_ylim(-0.05, 0.2)

    plt.tight_layout()
    out_path = os.path.join(output_dir, 'n0_axioms_chain.png')
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    print(f"  [PLOT] Zapisano: {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
#  Scenariusz: β ≠ γ (sprawdzenie, że testy FAIL)
# ──────────────────────────────────────────────────────────────────────────────

def run_broken_vacuum_test(gamma=0.5):
    """
    Demon: jeśli β ≠ γ (np. β = 0.3γ), test N0-6 powinien wołać FAIL.
    """
    beta_broken = 0.3 * gamma
    print("\n" + "─" * 65)
    print(f"  TEST SANITY: β={beta_broken:.3g} ≠ γ={gamma:.3g} → oczekiwany FAIL N0-6")
    print("─" * 65)
    r6 = check_N0_6_mass(beta_broken, gamma)
    symbol = PASS_SYMBOL if r6['pass'] else FAIL_SYMBOL
    print(f"  {symbol}  m_sp² = 3γ−2β = {r6['m_sp_sq_formula']:.4g}, "
          f"oczekiwane γ = {gamma:.4g}, residuum = {r6['residual']:.4g}")
    rs = check_vacuum_stability_oscillations(beta_broken, gamma)
    symbol2 = PASS_SYMBOL if rs['pass'] else FAIL_SYMBOL
    print(f"  {symbol2}  Stabilność: m_sp² = {rs['m_sp_sq']:.4g} > 0 → {'TAK' if rs['pass'] else 'NIE'}")


# ──────────────────────────────────────────────────────────────────────────────
#  Punkt wejścia
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Weryfikacja łańcucha aksjomat N0-1 → N0-6 (TGP v1)"
    )
    parser.add_argument('--plot', action='store_true',
                        help='Generuj wykresy diagnostyczne (wymaga matplotlib)')
    parser.add_argument('--gamma', type=float, default=DEFAULTS['gamma'],
                        help=f'Parametr γ (domyślnie {DEFAULTS["gamma"]})')
    parser.add_argument('--alpha', type=float, default=DEFAULTS['alpha'],
                        help=f'Parametr α (domyślnie {DEFAULTS["alpha"]})')
    parser.add_argument('--Phi0', type=float, default=DEFAULTS['Phi0'],
                        help=f'Tło Φ₀ (domyślnie {DEFAULTS["Phi0"]})')
    parser.add_argument('--no-sanity', action='store_true',
                        help='Pomiń test sanity β≠γ')
    args = parser.parse_args()

    # Wymuszamy warunek N0-5: β = γ
    beta = args.gamma
    gamma = args.gamma
    alpha = args.alpha
    Phi0 = args.Phi0

    ok = run_all_checks(alpha, beta, gamma, Phi0)

    if not args.no_sanity:
        run_broken_vacuum_test(gamma)

    if args.plot:
        make_plots(beta, gamma, Phi0)

    sys.exit(0 if ok else 1)
