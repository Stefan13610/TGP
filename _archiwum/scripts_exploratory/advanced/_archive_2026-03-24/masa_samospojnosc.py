"""
TGP — Samospojnosc masy: wyznaczenie dyskretnych wartosci K i mas generacji
============================================================================

IDEA:
Czastka = lokalne zaburzenie Yukawa nad kosmicznym tlem Phi_0:
    Phi(r) = Phi_0 + K * exp(-m_sp * r) / r

Phi_0 jest NIE absolutna proznia, lecz srednia kosmiczna (punkt wyjscia mierzony).
Masa cząstki = energia tego zaburzenia, samospojnie:

    E(K) = M * c^2 = (4*pi*Phi_0*c0^2/q) * K = Lambda * K

Szukamy przeciec krzywej E(K) z prosta E = Lambda*K.
Kazde przeciecie = jedna fizyczna generacja.

Nieliniowos TGP (potencjal -gamma*Phi^4/4) sprawia ze E(K)/K
rosnie, osiaga maksimum i maleje -> mozliwe 2 przeciecia dla malego Lambda.
Poszukujemy warunkow parametrow dla 3 przeciec (3 generacje).

JEDNOSTKI NATURALNE: Phi_0 = 1, m_sp = 1, c_0 = 1, l_P = 1
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ======================================================================
# PARAMETRY
# ======================================================================
PHI0  = 1.0
M_SP  = 1.0        # masa skalaru TGP = sqrt(gamma)
GAMMA = M_SP**2    # = 1.0
Q     = 1.0        # sprzezenie zrodla

# Stala samospojnosci: E = Lambda * K
# Lambda = 4*pi*Phi_0*c0^2/q
# Uwaga: a_Gamma wchodzi gdy uwzgledniamy smearing zrodla na siatce
# Bez smearingu: Lambda = 4*pi / q

A_GAMMA = 0.3      # krok sieci (1/m_sp) — wyznacza regularyzacje UV

# Parametry nieliniowe TGP
ALPHA = 0.5        # wsp. czlonu nielinijnego kinetycznego
BETA  = GAMMA      # warunek prozni: beta = gamma (domyslnie)

# ======================================================================
# PROFIL POLA — rozwiazanie pelnego rownania nieliniowego BVP
# ======================================================================
def solve_field_profile(K, r_min=None, r_max=30.0, N=300):
    """
    Rozwiazuje radialne rownanie pola TGP:
        Phi'' + (2/r)*Phi' + beta*Phi^2 - gamma*Phi^3 = 0  (bez zrodla)
    z warunkami brzegowymi:
        Phi(r_min) = Phi0 + K/r_min * exp(-m_sp*r_min)   [Yukawa przy zrodle]
        Phi(r_max) = Phi0                                  [proznia na nieskoncznosci]

    Zwraca (r, phi) — profil pola.
    """
    if r_min is None:
        r_min = A_GAMMA

    r = np.linspace(r_min, r_max, N)

    # Warunki brzegowe (Yukawa na lewej, proznia na prawej)
    phi_left  = PHI0 + K * np.exp(-M_SP * r_min) / r_min
    phi_right = PHI0

    def ode(r, y):
        phi, dphi = y
        phi = np.maximum(phi, 1e-10)
        # Rownanie: phi'' = -2/r * phi' - beta*phi^2 + gamma*phi^3
        ddphi = -2.0/r * dphi - BETA*phi**2 + GAMMA*phi**3
        return np.array([dphi, ddphi])

    def bc(ya, yb):
        return np.array([
            ya[0] - phi_left,   # Phi(r_min) = Yukawa
            yb[0] - phi_right   # Phi(r_max) = Phi0
        ])

    # Poczatkowe przyblizenie: Yukawa
    phi_init  = PHI0 + K * np.exp(-M_SP * r) / r
    dphi_init = K * np.exp(-M_SP * r) * (-M_SP * r - 1) / r**2

    y_init = np.vstack([phi_init, dphi_init])

    try:
        sol = solve_bvp(ode, bc, r, y_init, max_nodes=3000, tol=1e-5)
        if sol.success:
            return sol.x, sol.y[0]
    except Exception:
        pass

    # Fallback: profil Yukawa
    return r, phi_init


# ======================================================================
# ENERGIA TGP: E(K) = energia zaburzenia nad tlem Phi_0
# ======================================================================
def energy_TGP(K, r_min=None, r_max=30.0, N=400, use_yukawa=False):
    """
    E(K) = 4*pi * integral_r_min^r_max [
        (1/2)(dPhi/dr)^2 * (1 + alpha/(Phi0*Phi))
        + V(Phi) - V(Phi0)
    ] * r^2 dr

    gdzie V(Phi) = gamma/3 * (Phi/Phi0)^3 * Phi0^4 / Phi0
                 - gamma/4 * (Phi/Phi0)^4 * ...
    W znormalizowanych jednostkach (Phi0=1):
        V(phi) = gamma/3 * phi^3 - gamma/4 * phi^4
        V(Phi0=1) = gamma/3 - gamma/4 = gamma/12
    """
    if r_min is None:
        r_min = A_GAMMA

    if use_yukawa:
        r = np.linspace(r_min, r_max, N)
        phi = PHI0 + K * np.exp(-M_SP * r) / r
        dphi = K * np.exp(-M_SP * r) * (-M_SP * r - 1) / r**2
    else:
        r, phi = solve_field_profile(K, r_min=r_min, r_max=r_max, N=N)
        dphi = np.gradient(phi, r)

    phi = np.maximum(phi, 1e-10)

    # Czlon kinetyczny standardowy
    E_kin_std = 0.5 * dphi**2

    # Czlon kinetyczny TGP nieliniowy
    E_kin_nl  = 0.5 * ALPHA * dphi**2 / (PHI0 * phi)

    # Potencjal (ponad proznia)
    V_phi  = GAMMA/3.0 * phi**3  - GAMMA/4.0 * phi**4
    V_phi0 = GAMMA/3.0 * PHI0**3 - GAMMA/4.0 * PHI0**4
    dV = V_phi - V_phi0

    # Gestosc energii radialnie wazena
    integrand = (E_kin_std + E_kin_nl + dV) * r**2

    E = 4.0 * np.pi * np.trapezoid(integrand, r)
    return E


# ======================================================================
# ROWNANIE SAMOSPOJNOSCI
# ======================================================================
def Lambda(a_gamma=A_GAMMA):
    """
    Nachylenie prostej samospojnosci:
    E = Lambda * K
    Lambda = 4*pi*Phi0*c0^2 / q
    (smearing zrodla na siatce daje dodatkowy czynnik a_Gamma)
    """
    return 4.0 * np.pi * PHI0 / Q  # bez a_Gamma (zrodlo punktowe)


def f_selfconsist(K, lam):
    """f(K) = E(K) - Lambda*K. Zera = masy generacji."""
    return energy_TGP(K) - lam * K


# ======================================================================
# SKAN E(K) i szukanie przeciec
# ======================================================================
def scan_and_find_masses(K_min=0.05, K_max=8.0, N_K=120, lam=None):
    """
    Skanuje E(K) i szuka przeciec z prosta Lambda*K.
    Zwraca liste mas (wartosci K przy przecieciach).
    """
    if lam is None:
        lam = Lambda()

    K_arr = np.linspace(K_min, K_max, N_K)
    E_arr = np.array([energy_TGP(K) for K in K_arr])
    f_arr = E_arr - lam * K_arr

    crossings_K = []
    crossings_E = []

    for i in range(len(f_arr) - 1):
        if np.isfinite(f_arr[i]) and np.isfinite(f_arr[i+1]):
            if f_arr[i] * f_arr[i+1] < 0:
                try:
                    K_cross = brentq(
                        lambda K: energy_TGP(K) - lam * K,
                        K_arr[i], K_arr[i+1],
                        xtol=1e-6, maxiter=50
                    )
                    E_cross = energy_TGP(K_cross)
                    crossings_K.append(K_cross)
                    crossings_E.append(E_cross)
                except Exception:
                    pass

    return K_arr, E_arr, crossings_K, crossings_E


# ======================================================================
# SKAN PARAMETROW: szukanie warunkow dla 3 generacji
# ======================================================================
def scan_parameters():
    """
    Skanuje (alpha, a_Gamma) w poszukiwaniu 3 dyskretnych mas
    z hierarchia 1 : 207 : 3477.
    """
    print("\n" + "="*60)
    print("SKAN PARAMETROW (alpha, a_Gamma)")
    print("="*60)

    alpha_vals  = [0.1, 0.5, 1.0, 2.0, 5.0]
    a_gam_vals  = [0.1, 0.3, 0.5, 1.0, 2.0]

    results = []

    for a_g in a_gam_vals:
        for alp in alpha_vals:
            global ALPHA, A_GAMMA
            ALPHA   = alp
            A_GAMMA = a_g

            lam = Lambda()
            K_arr, E_arr, cK, cE = scan_and_find_masses(
                K_min=0.02, K_max=10.0, N_K=80, lam=lam
            )

            n_cross = len(cK)
            ratios = []
            if n_cross >= 2:
                ratios = [cE[i]/cE[0] for i in range(1, n_cross)]

            results.append({
                'alpha': alp, 'a_gamma': a_g,
                'n_cross': n_cross,
                'masses': cE,
                'ratios': ratios
            })

    return results


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':

    print("="*60)
    print("TGP — Samospojnosc masy (masa z energii zaburzenia)")
    print(f"Parametry poczatkowe: alpha={ALPHA}, a_Gamma={A_GAMMA}")
    print(f"gamma={GAMMA}, beta={BETA}, q={Q}")
    print("="*60)

    lam = Lambda()
    print(f"\nStala samospojnosci Lambda = {lam:.4f}")
    print(f"(E = Lambda*K = {lam:.4f} * K)")

    # --- Skan glowny ---
    print("\nObliczanie E(K)...")
    K_arr, E_arr, cK, cE = scan_and_find_masses(
        K_min=0.05, K_max=8.0, N_K=100, lam=lam
    )

    # --- Wyniki ---
    print(f"\n{'='*60}")
    print(f"Znaleziono {len(cK)} rozwiazania samospojne:")
    for i, (K, E) in enumerate(zip(cK, cE)):
        print(f"  Generacja {i+1}: K = {K:.5f},  m = E/c^2 = {E:.5f}")
    if len(cK) >= 2:
        print(f"\nStosunek mas:")
        for i in range(1, len(cK)):
            print(f"  m{i+1}/m1 = {cE[i]/cE[0]:.2f}   (cel: {[207, 3477][i-1] if i<=2 else '?'})")

    if len(cK) == 0:
        # Diagnostyka
        print("\n[Diagnostyka] E(K)/K w kilku punktach:")
        for K_test in [0.1, 0.5, 1.0, 2.0, 4.0, 8.0]:
            E_test = energy_TGP(K_test)
            print(f"  K={K_test:.1f}:  E={E_test:.4f},  E/K={E_test/K_test:.4f},  Lambda={lam:.4f}")

    # --- Wykres E(K) i prosta ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    valid = np.isfinite(E_arr)
    ax.plot(K_arr[valid], E_arr[valid], 'b-', lw=2, label=r'$E(K)$ [TGP]')
    ax.plot(K_arr, lam * K_arr, 'r--', lw=1.5, label=fr'$\Lambda K = {lam:.2f}\,K$')
    for K, E in zip(cK, cE):
        ax.axvline(K, color='g', lw=0.8, ls=':')
        ax.plot(K, E, 'go', ms=9, zorder=5)
    ax.set_xlabel(r'$K$ (sila zrodla)', fontsize=12)
    ax.set_ylabel(r'$E(K)$', fontsize=12)
    ax.set_title('Samospojnosc masy: E(K) vs \u039bK', fontsize=11)
    ax.legend()
    ax.grid(True, alpha=0.3)
    if len(cK) > 0:
        ax.set_ylim(bottom=0)

    # --- E(K)/K vs K ---
    ax2 = axes[1]
    ratio_arr = np.where(np.abs(K_arr) > 1e-10, E_arr / K_arr, np.nan)
    ax2.plot(K_arr[valid], ratio_arr[valid], 'b-', lw=2, label=r'$E(K)/K$')
    ax2.axhline(lam, color='r', lw=1.5, ls='--', label=fr'$\Lambda = {lam:.2f}$')
    for K in cK:
        ax2.axvline(K, color='g', lw=0.8, ls=':', label=f'K*={K:.3f}')
    ax2.set_xlabel(r'$K$', fontsize=12)
    ax2.set_ylabel(r'$E(K)/K$', fontsize=12)
    ax2.set_title(r'Efektywna stroma — przeciecia = masy', fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        'TGP/TGP_v1/scripts/advanced/masa_samospojnosc.png',
        dpi=150, bbox_inches='tight'
    )
    print("\nWykres zapisany: masa_samospojnosc.png")

    # --- Skan parametrow ---
    print("\nRozpoczecie skanu parametrow (alpha x a_Gamma)...")
    results = scan_parameters()

    print("\n" + "="*60)
    print("WYNIKI SKANU:")
    print(f"{'alpha':>6} {'a_Gam':>6} {'n_cross':>8} {'masy':>30} {'ratio m2/m1':>12}")
    print("-"*70)
    for r in results:
        masy_str = ", ".join(f"{m:.3f}" for m in r['masses']) if r['masses'] else "brak"
        rat_str  = f"{r['ratios'][0]:.1f}" if r['ratios'] else "—"
        print(f"{r['alpha']:>6.1f} {r['a_gamma']:>6.1f} {r['n_cross']:>8d}  {masy_str:>30}  {rat_str:>12}")

    # Podsumowanie najlepszych
    best = [r for r in results if r['n_cross'] >= 2]
    if best:
        print(f"\nNajlepsze dopasowanie do 1:207:")
        for r in best:
            if r['ratios']:
                dist = abs(r['ratios'][0] - 207)
                print(f"  alpha={r['alpha']}, a_Gamma={r['a_gamma']}: "
                      f"m2/m1={r['ratios'][0]:.1f}  (delta={dist:.1f})")

    print("\nGotowe.")
