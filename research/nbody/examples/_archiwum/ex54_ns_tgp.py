"""
ex54_ns_tgp.py
==============
TGP — numeryczne obliczenie indeksu spektralnego n_s perturbacji kosmologicznych.

TEORIA:
  Równanie pola w FRW (z sek08_formalizm.tex, stw. prop:TGP-FRW-full):

    φ'' + 3H φ' + V'_mod(φ) = 0                            (1)
    H² = [φ'²/2 + V_mod(φ)] / (3 M*²)                     (2)
    V_mod(φ) = (γ/3)φ³ − (γ/4)φ⁴                           (3)

  gdzie ' = d/dN (pochodna po e-foldingach N = ln a).

  Parametry powolnego staczania:
    ε_H = −Ḣ/H²  =  (φ'²/2) / M*²H²
    ε_ψ = φ'²/(2H²φ²)    [TGP-specyficzny dodatkowy wkład]

  Indeks spektralny TGP (eq:ns-TGP):
    n_s − 1 = −4ε_H − 4ε_ψ

  (Różni się od standardowego n_s−1 = −2ε_H − η przez czynnik 4 i ε_ψ)

FIZYCZNA INTERPRETACJA:
  Dodatkowy czynnik 4 (zamiast 2) i wkład ε_ψ wynikają z miary K(ψ)=ψ⁴
  w działaniu TGP. Zmienna Mukhanova-Sasakiego: v_k = a·ψ_bg²·δψ_k.

PARAMETRY WARSTWY I + II (tab. param-classification):
  γ:        wyznaczane z Λ_obs ≈ γ/12 → γ ≈ 12 H₀²/c₀²
  Φ₀:       ≈ 24.66 [H₀²/c₀²]⁻¹
  ψ_ini:    = 7/6 (hipoteza GL — atraktor Ginzburg-Landau)

WYNIKI OCZEKIWANE (analiza v33):
  n_s ∈ [0.950, 0.975] przy ψ_ini = 7/6

WERYFIKACJA: porównanie z Planck 2018: n_s = 0.9649 ± 0.0042

Użycie:
  python ex54_ns_tgp.py [--psi-ini FLOAT] [--gamma FLOAT] [--Nefold INT] [--plot]

Autor: TGP v1 sesja v33 (2026-03-27)
"""

import numpy as np
from scipy.integrate import solve_ivp
import argparse
import sys

# ─────────────────────────────────────────────────────────────
# PARAMETRY DOMYŚLNE (Warstwa II TGP)
# ─────────────────────────────────────────────────────────────
GAMMA_DEFAULT = 12.0          # γ w jednostkach H₀²/c₀² (N0-5, N0-6)
PSI_INI_DEFAULT = 7.0 / 6.0  # warunek początkowy φ(0) (hipoteza GL)
NEFOLD = 70                   # liczba e-foldingów do integracji
MSTAR_SQ = 1.0                # M*² = m_P² = 1 (jednostki naturalne)

# ─────────────────────────────────────────────────────────────
# POTENCJAŁ I JEGO POCHODNA
# ─────────────────────────────────────────────────────────────

def V_mod(phi, gamma):
    """V_mod(φ) = (γ/3)φ³ − (γ/4)φ⁴  [warunek N0-5: β=γ]"""
    return (gamma / 3.0) * phi**3 - (gamma / 4.0) * phi**4

def dV_mod(phi, gamma):
    """dV_mod/dφ = γφ² − γφ³"""
    return gamma * phi**2 - gamma * phi**3

def d2V_mod(phi, gamma):
    """d²V_mod/dφ² = 2γφ − 3γφ²"""
    return 2.0 * gamma * phi - 3.0 * gamma * phi**2

# ─────────────────────────────────────────────────────────────
# UKŁAD RÓWNAŃ W ZMIENNEJ N = ln a
# ─────────────────────────────────────────────────────────────
# Zmienne stanu: y = [φ, φ' = dφ/dN]
# Równanie Friedmanna: H² = [φ'²/(2M*²) + V/M*²] / 3
# Równanie pola:  φ'' + (3 − ε_H)φ' + V'/H² = 0
# gdzie ε_H = φ'²/(2M*²H²)

def equations(N, y, gamma):
    """
    Układ ODE w zmiennej N = ln(a):
    dy[0]/dN = y[1]          [φ' = dφ/dN]
    dy[1]/dN = RHS           [φ'' = d²φ/dN²]
    """
    phi = y[0]
    dphi = y[1]   # dφ/dN

    V = V_mod(phi, gamma)
    dV = dV_mod(phi, gamma)

    # H² z równania Friedmanna (M*=1)
    # H² = (φ'²/2 + V) / 3   [jednostki: M*²=1]
    H2 = (dphi**2 / 2.0 + V) / 3.0

    if H2 <= 0:
        # Inflacja skończyła się
        return [dphi, 0.0]

    H = np.sqrt(H2)

    # ε_H = −dH/H²dN = dphi²/(2H²)
    eps_H = dphi**2 / (2.0 * H2)

    # Równanie pola: φ'' + (3 − ε_H)φ' + dV/H² = 0
    d2phi = -(3.0 - eps_H) * dphi - dV / H2

    return [dphi, d2phi]

# ─────────────────────────────────────────────────────────────
# PARAMETRY POWOLNEGO STACZANIA
# ─────────────────────────────────────────────────────────────

def slow_roll_params(phi, dphi, gamma):
    """
    Oblicza ε_H, ε_ψ, η_H na podstawie φ i φ'=dφ/dN.
    Zwraca dict z parametrami i n_s.
    """
    V = V_mod(phi, gamma)
    dV = dV_mod(phi, gamma)
    d2V = d2V_mod(phi, gamma)

    H2 = (dphi**2 / 2.0 + V) / 3.0
    if H2 <= 0:
        return None

    # ε_H = −Ḣ/H² = φ'²/(2H²)
    eps_H = dphi**2 / (2.0 * H2)

    # ε_ψ = φ'²/(2H²φ²)  — specyficzny dla TGP
    # (wynika z K(ψ)=ψ⁴ w zmiennej Mukhanova-Sasakiego)
    if phi**2 > 1e-15:
        eps_psi = dphi**2 / (2.0 * H2 * phi**2)
    else:
        eps_psi = 0.0

    # η_H = dε_H/ε_H dN (standardowy)
    # Uproszczone: η ≈ −V''/(3H²) = −d2V/(3H²) w granicy slow-roll
    eta_H = -d2V / (3.0 * H2)

    # n_s TGP = 1 − 4ε_H − 4ε_ψ  (eq:ns-TGP)
    ns_TGP = 1.0 - 4.0 * eps_H - 4.0 * eps_psi

    # n_s standardowy dla porównania: 1 − 6ε_H + 2η_H
    ns_std = 1.0 - 6.0 * eps_H + 2.0 * eta_H

    return {
        'eps_H': eps_H,
        'eps_psi': eps_psi,
        'eta_H': eta_H,
        'H2': H2,
        'H': np.sqrt(H2),
        'n_s_TGP': ns_TGP,
        'n_s_standard': ns_std,
    }

# ─────────────────────────────────────────────────────────────
# WARUNKI POCZĄTKOWE
# ─────────────────────────────────────────────────────────────

def initial_conditions_TGP(psi_ini, gamma, dpsi_frac=0.01):
    """
    Warunki początkowe dla inflacji TGP.

    psi_ini = 7/6: hipoteza GL (atraktor Ginzburg-Landau)
    Pole zaczyna lekko poniżej równowagi (φ_eq = 4/3) z małą prędkością.

    Warunek próżniowy N0-5: minimum V_mod w φ=0 i max przy φ=1,
    ale V_mod(4/3)=0 (punkt siodłowy, nie minimum!).
    Fizyczne minimum TGP: φ=1 (stabilna próżnia).

    Wartość ψ_ini = 7/6 jest bliska 1 (próżni) z strony φ>1.
    """
    phi0 = psi_ini
    # Prędkość początkowa: small φ'  (slow-roll start)
    # W granicy slow-roll: φ' ≈ -dV/(3H)
    V0 = V_mod(phi0, gamma)
    dV0 = dV_mod(phi0, gamma)

    if V0 > 0:
        H0 = np.sqrt(V0 / 3.0)  # aproksymacja dominacji potencjału
        dphi0 = -dV0 / (3.0 * H0)
    else:
        dphi0 = -dpsi_frac * phi0

    return [phi0, dphi0]

# ─────────────────────────────────────────────────────────────
# GŁÓWNA FUNKCJA OBLICZENIOWA
# ─────────────────────────────────────────────────────────────

def compute_ns_TGP(psi_ini=PSI_INI_DEFAULT, gamma=GAMMA_DEFAULT,
                   N_total=NEFOLD, N_horizon_crossing=50.0,
                   verbose=True):
    """
    Całkuje równania TGP i oblicza n_s przy przekroczeniu horyzontu
    N_* e-foldingów przed końcem inflacji.

    Parametry:
      psi_ini:           φ(0) = warunek początkowy
      gamma:             γ (w jednostkach H₀²/c₀²)
      N_total:           łączna liczba e-foldingów
      N_horizon_crossing: moment wyjścia poza horyzont (typowo N*≈50-60)

    Zwraca: dict z n_s i parametrami slow-roll
    """
    y0 = initial_conditions_TGP(psi_ini, gamma)

    if verbose:
        print(f"\n{'='*60}")
        print(f"TGP: obliczanie n_s")
        print(f"  ψ_ini  = {psi_ini:.6f}  (hipoteza GL: 7/6 = {7/6:.6f})")
        print(f"  γ      = {gamma:.4f}  [jednostki H₀²/c₀²]")
        print(f"  φ(0)   = {y0[0]:.4f},  φ'(0) = {y0[1]:.6f}")
        print(f"  N_tot  = {N_total},  N_*  = {N_horizon_crossing}")
        print(f"{'='*60}")

    # Całkowanie
    N_span = (0, N_total)
    N_eval = np.linspace(0, N_total, 5000)

    def event_end_inflation(N, y, gamma):
        """Koniec inflacji: ε_H = 1."""
        phi, dphi = y
        H2 = (dphi**2 / 2.0 + V_mod(phi, gamma)) / 3.0
        if H2 <= 0:
            return -1.0
        eps_H = dphi**2 / (2.0 * H2)
        return eps_H - 1.0
    event_end_inflation.terminal = True
    event_end_inflation.direction = 1

    sol = solve_ivp(
        equations,
        N_span,
        y0,
        args=(gamma,),
        method='DOP853',
        t_eval=N_eval,
        rtol=1e-10,
        atol=1e-12,
        events=event_end_inflation,
        dense_output=True,
    )

    if not sol.success:
        print(f"BŁĄD: całkowanie nie powiodło się: {sol.message}")
        return None

    N_arr = sol.t
    phi_arr = sol.y[0]
    dphi_arr = sol.y[1]

    # Sprawdzenie czy inflacja zakończyła się
    N_end = N_arr[-1]
    if verbose:
        print(f"\n  Koniec inflacji w N_end = {N_end:.2f}")

    # Obliczenie n_s przy N* e-foldingów przed końcem inflacji
    N_star = N_end - (N_total - N_horizon_crossing)
    if N_star < 0 or N_star > N_end:
        N_star = N_end * 0.7  # fallback
        if verbose:
            print(f"  UWAGA: Korekta N* → {N_star:.2f}")

    # Interpolacja w N_star
    phi_star = float(sol.sol(N_star)[0])
    dphi_star = float(sol.sol(N_star)[1])

    sr = slow_roll_params(phi_star, dphi_star, gamma)
    if sr is None:
        print("BŁĄD: H²≤0 przy N*")
        return None

    # Kompilacja wyników
    result = {
        'n_s_TGP': sr['n_s_TGP'],
        'n_s_standard': sr['n_s_standard'],
        'eps_H': sr['eps_H'],
        'eps_psi': sr['eps_psi'],
        'eta_H': sr['eta_H'],
        'phi_star': phi_star,
        'dphi_star': dphi_star,
        'H_star': sr['H'],
        'N_end': N_end,
        'N_star': N_star,
        'gamma': gamma,
        'psi_ini': psi_ini,
    }

    if verbose:
        print(f"\n  Przy przekroczeniu horyzontu (N* = {N_star:.2f}):")
        print(f"  φ*      = {phi_star:.6f}")
        print(f"  φ'*     = {dphi_star:.6f}")
        print(f"  ε_H     = {sr['eps_H']:.6f}")
        print(f"  ε_ψ     = {sr['eps_psi']:.6f}")
        print(f"  η_H     = {sr['eta_H']:.6f}")
        print(f"\n  {'─'*40}")
        print(f"  n_s (TGP)      = {sr['n_s_TGP']:.6f}")
        print(f"  n_s (standard) = {sr['n_s_standard']:.6f}")
        print(f"  Planck 2018:   n_s = 0.9649 ± 0.0042")
        ns = sr['n_s_TGP']
        sigma = abs(ns - 0.9649) / 0.0042
        print(f"  Odchylenie TGP od Plancka: {sigma:.2f}σ")
        print(f"  {'─'*40}")

    return result

# ─────────────────────────────────────────────────────────────
# SKAN ψ_ini
# ─────────────────────────────────────────────────────────────

def scan_psi_ini(gamma=GAMMA_DEFAULT, psi_range=None, n_pts=20):
    """
    Skan n_s(TGP) jako funkcji ψ_ini.
    Pozwala ocenić wrażliwość na warunek początkowy.
    """
    if psi_range is None:
        psi_range = (1.05, 1.50)

    psi_arr = np.linspace(psi_range[0], psi_range[1], n_pts)
    ns_arr = []

    print(f"\nSKAN ψ_ini ∈ [{psi_range[0]:.2f}, {psi_range[1]:.2f}], γ = {gamma:.4f}")
    print(f"{'ψ_ini':>10} {'n_s(TGP)':>12} {'ε_H':>10} {'ε_ψ':>10}  Planck σ")
    print("─" * 65)

    planck_ns = 0.9649
    planck_err = 0.0042

    for psi in psi_arr:
        res = compute_ns_TGP(psi_ini=psi, gamma=gamma, verbose=False)
        if res is not None:
            ns = res['n_s_TGP']
            eps_H = res['eps_H']
            eps_psi = res['eps_psi']
            sigma = abs(ns - planck_ns) / planck_err
            ns_arr.append(ns)
            marker = " ←★ GL" if abs(psi - 7/6) < 0.01 else ""
            print(f"{psi:>10.4f} {ns:>12.6f} {eps_H:>10.6f} {eps_psi:>10.6f}  "
                  f"{sigma:6.2f}σ{marker}")
        else:
            ns_arr.append(np.nan)
            print(f"{psi:>10.4f} {'FAILED':>12}")

    return psi_arr, np.array(ns_arr)

# ─────────────────────────────────────────────────────────────
# OPTYMALNE ψ_ini → n_s = 0.9649
# ─────────────────────────────────────────────────────────────

def find_psi_for_planck(gamma=GAMMA_DEFAULT, target_ns=0.9649):
    """
    Szuka ψ_ini dającego dokładnie n_s = target_ns.
    Metoda bisekcji.
    """
    from scipy.optimize import brentq

    def ns_minus_target(psi):
        res = compute_ns_TGP(psi_ini=psi, gamma=gamma, verbose=False)
        if res is None:
            return np.nan
        return res['n_s_TGP'] - target_ns

    try:
        psi_opt = brentq(ns_minus_target, 1.05, 1.50, xtol=1e-5, maxiter=50)
        res = compute_ns_TGP(psi_ini=psi_opt, gamma=gamma, verbose=False)
        print(f"\nOptymalne ψ_ini dla n_s = {target_ns:.4f}:")
        print(f"  ψ_ini* = {psi_opt:.6f}")
        print(f"  7/6    = {7/6:.6f}   (hipoteza GL)")
        print(f"  Δψ     = {abs(psi_opt - 7/6):.6f}  ({100*abs(psi_opt-7/6)/(7/6):.2f}%)")
        if res:
            print(f"  n_s    = {res['n_s_TGP']:.6f}")
        return psi_opt
    except Exception as e:
        print(f"Nie znaleziono: {e}")
        return None

# ─────────────────────────────────────────────────────────────
# WIZUALIZACJA
# ─────────────────────────────────────────────────────────────

def plot_results(psi_arr, ns_arr, gamma, res_gl):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib niedostępny — pomijam wykresy")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f'TGP: Indeks spektralny $n_s$  |  $\\gamma = {gamma:.2f}\\, H_0^2/c_0^2$',
        fontsize=13)

    # Lewy: n_s(ψ_ini)
    ax = axes[0]
    valid = ~np.isnan(ns_arr)
    ax.plot(psi_arr[valid], ns_arr[valid], 'b-o', ms=4, lw=1.5,
            label='TGP: $n_s = 1 - 4\\varepsilon_H - 4\\varepsilon_\\psi$')
    ax.axhline(0.9649, color='red', lw=1.5, linestyle='--',
               label='Planck 2018: $n_s = 0.9649$')
    ax.axhspan(0.9649 - 0.0042, 0.9649 + 0.0042, alpha=0.15, color='red',
               label='$1\\sigma$ Planck')
    ax.axvline(7/6, color='green', lw=1.5, linestyle=':', label='$\\psi_{\\rm ini} = 7/6$ (GL)')
    if res_gl:
        ax.plot(7/6, res_gl['n_s_TGP'], 'g*', ms=14, zorder=5,
                label=f'GL: $n_s = {res_gl["n_s_TGP"]:.4f}$')
    ax.set_xlabel('$\\psi_{\\rm ini} = \\phi(0)$', fontsize=12)
    ax.set_ylabel('$n_s$', fontsize=12)
    ax.set_title('Wrażliwość na warunek początkowy')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.92, 1.01)

    # Prawy: ewolucja φ(N) dla ψ_ini = 7/6
    ax2 = axes[1]
    if res_gl:
        y0 = initial_conditions_TGP(7/6, gamma)
        N_span = (0, res_gl['N_end'] + 5)
        N_eval = np.linspace(0, res_gl['N_end'] + 5, 1000)

        def event_end_inflation2(N, y, gamma):
            phi, dphi = y
            H2 = (dphi**2 / 2.0 + V_mod(phi, gamma)) / 3.0
            if H2 <= 0:
                return -1.0
            return dphi**2 / (2.0 * H2) - 1.0
        event_end_inflation2.terminal = True
        event_end_inflation2.direction = 1

        sol = solve_ivp(equations, N_span, y0, args=(gamma,),
                        method='DOP853', t_eval=N_eval,
                        rtol=1e-10, atol=1e-12, events=event_end_inflation2)

        if sol.success and len(sol.t) > 1:
            ax2.plot(sol.t, sol.y[0], 'b-', lw=2, label='$\\phi(N)$')
            ax2.axhline(1.0, color='gray', lw=1, linestyle='--', alpha=0.7,
                        label='Próżnia $\\phi=1$')
            ax2.axvline(res_gl['N_star'], color='orange', lw=1.5, linestyle='-.',
                        label=f'$N_* = {res_gl["N_star"]:.1f}$ (horyzont)')
            ax2.scatter([res_gl['N_star']], [res_gl['phi_star']], color='red', s=80, zorder=5)
            ax2.set_xlabel('$N = \\ln a$  (e-foldingi)', fontsize=12)
            ax2.set_ylabel('$\\phi = \\Phi/\\Phi_0$', fontsize=12)
            ax2.set_title(f'Trajektoria tła TGP ($\\psi_{{\\rm ini}}=7/6$)')
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fname = 'ex54_ns_tgp.png'
    plt.savefig(fname, dpi=130, bbox_inches='tight')
    print(f"\nWykres zapisany: {fname}")
    plt.close()

# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='TGP: Obliczenie indeksu spektralnego n_s')
    parser.add_argument('--psi-ini', type=float, default=PSI_INI_DEFAULT,
                        help=f'Warunek początkowy φ(0) (default: {PSI_INI_DEFAULT:.6f} = 7/6)')
    parser.add_argument('--gamma', type=float, default=GAMMA_DEFAULT,
                        help=f'Parametr γ (default: {GAMMA_DEFAULT})')
    parser.add_argument('--Nefold', type=int, default=NEFOLD,
                        help=f'Łączna liczba e-foldingów (default: {NEFOLD})')
    parser.add_argument('--scan', action='store_true',
                        help='Wykonaj skan n_s(ψ_ini)')
    parser.add_argument('--find-optimal', action='store_true',
                        help='Znajdź ψ_ini dające n_s = 0.9649 (Planck)')
    parser.add_argument('--plot', action='store_true',
                        help='Rysuj wykresy (wymaga matplotlib)')
    args = parser.parse_args()

    print("=" * 65)
    print("  TGP — Indeks spektralny n_s (ex54_ns_tgp.py)")
    print("=" * 65)

    # Obliczenie dla ψ_ini domyślnego (GL)
    res_gl = compute_ns_TGP(
        psi_ini=args.psi_ini,
        gamma=args.gamma,
        N_total=args.Nefold,
        verbose=True)

    # PASS/FAIL
    print("\n" + "─" * 50)
    print("WERYFIKACJA:")
    passed = 0
    tests = 0

    if res_gl:
        ns = res_gl['n_s_TGP']
        # T1: n_s w oknie TGP [0.940, 0.985]
        tests += 1
        ok = 0.940 < ns < 0.985
        print(f"  T1 [n_s ∈ (0.940, 0.985)]:    {ns:.5f}  {'PASS ✓' if ok else 'FAIL ✗'}")
        if ok: passed += 1

        # T2: ε_H > 0 (wolne staczanie aktywne)
        tests += 1
        ok = res_gl['eps_H'] > 0 and res_gl['eps_H'] < 1
        print(f"  T2 [0 < ε_H < 1]:             {res_gl['eps_H']:.6f}  {'PASS ✓' if ok else 'FAIL ✗'}")
        if ok: passed += 1

        # T3: ε_ψ < ε_H (wkład TGP nie dominuje)
        tests += 1
        ok = res_gl['eps_psi'] < res_gl['eps_H']
        print(f"  T3 [ε_ψ < ε_H]:               {res_gl['eps_psi']:.6f} < {res_gl['eps_H']:.6f}  {'PASS ✓' if ok else 'FAIL ✗'}")
        if ok: passed += 1

        # T4: n_s < 1 (czerwone przechylenie)
        tests += 1
        ok = ns < 1.0
        print(f"  T4 [n_s < 1 (czerwone)]:      {'PASS ✓' if ok else 'FAIL ✗'}")
        if ok: passed += 1

    print(f"\n  WYNIK: {passed}/{tests} PASS")
    print("─" * 50)

    # Opcje dodatkowe
    psi_arr = ns_arr = None
    if args.scan:
        psi_arr, ns_arr = scan_psi_ini(gamma=args.gamma)

    if args.find_optimal:
        find_psi_for_planck(gamma=args.gamma)

    if args.plot:
        if psi_arr is None:
            psi_arr, ns_arr = scan_psi_ini(gamma=args.gamma, n_pts=30)
        plot_results(psi_arr, ns_arr, args.gamma, res_gl)

    return passed == tests


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
