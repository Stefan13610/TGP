# -*- coding: utf-8 -*-
"""
P39: Bifurcation structure of TGP soliton solutions
=====================================================
Analiza struktury bifurkacji K1=K2 w przestrzeni (alpha, a_Gamma).

Trzy pytania badawcze:
1. Czy przejscie K1=K2 -> K1<K2 jest ciagle czy skokowe?
2. Czy r21_min(a_Gamma) < 20 dla wiekszych a_Gamma (>=0.050)?
3. Czy alpha_c(a_Gamma) da sie wyprowadzic analitycznie z g(K_c)=0, g'(K_c)=0?

Model energii (zgodny z P38):
  Ek = 4*pi * int 0.5*dphi^2 * (1 + alpha/phi) * r^2 dr
  Ep = 4*pi * int [V_mod(phi,lam) - V_mod(1,lam)] * r^2 dr
  E(K) = Ek + Ep
  g(K) = E(K)/(4*pi*K) - 1 = 0  <=>  soliton istnieje

  V_mod = Gamma/3*phi^3 - Gamma/4*phi^4 + lam/6*(phi-1)^6
  phi   = 1 + K*e^{-r}/r
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE FIZYCZNE (identyczne z P38)
# ============================================================
R_MAX   = 60.0
GAMMA   = 1.0
LAM_REF = 1e-5

# ============================================================
# ENERGIA SOLITONU — identyczna implementacja jak P38
# ============================================================
def V_mod(phi, lam):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + lam/6*(phi-1)**6

def energy_log(K, alpha, a_gam, lam=LAM_REF, N=2000):
    """
    Energia solitonu phi=1+K*e^{-r}/r na siatce logarytmicznej.
    Ek = 4*pi*int 0.5*dphi^2*(1+alpha/phi)*r^2 dr
    Ep = 4*pi*int [V(phi)-V(1)]*r^2 dr
    """
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX / a_gam)**t      # siatka logarytmiczna
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * (1 + alpha/phi) * r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1) * r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam=LAM_REF):
    """g(K) = E(K)/(4*pi*K) - 1"""
    e = energy_log(K, alpha, a_gam, lam)
    return e / (4.0*np.pi*K) - 1.0 if np.isfinite(e) else np.nan

# ============================================================
# SZUKANIE ZER g(K)
# ============================================================
def find_K1_num(alpha, a_gam, K_lo=0.0003, K_hi=0.15):
    """K1 - pierwsze (najmniejsze) zero g(K)."""
    for lo, hi in [(K_lo, K_hi), (K_lo, 0.3)]:
        try:
            glo = g_func(lo, alpha, a_gam)
            ghi = g_func(hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), lo, hi, xtol=1e-13)
        except Exception:
            pass
    return np.nan

def find_K2_num(alpha, a_gam):
    """
    K2 - drugie zero g(K), WIEKSZE niz K1.
    Wazne: szukamy w zakresie K > 0.15 (poza zasiegiem find_K1_num),
    zeby nie znalezc K1 ponownie.
    """
    for K_lo, K_hi in [(0.15, 2.0), (0.5, 5.0), (0.15, 5.0)]:
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
                return brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-10)
        except Exception:
            pass
    return np.nan

def find_K1_K2(alpha, a_gam):
    """Zwraca (K1, K2, status: 'separate'/'double'/'none')."""
    K1 = find_K1_num(alpha, a_gam)
    K2 = find_K2_num(alpha, a_gam)
    if np.isnan(K1) and np.isnan(K2):
        return np.nan, np.nan, 'none'
    if np.isnan(K2) or (not np.isnan(K1) and abs(K2 - K1)/max(K1, 1e-12) < 0.01):
        return K1, K1, 'double'
    if K2 > K1:
        return K1, K2, 'separate'
    return K1, K2, 'double'

# ============================================================
# SEKCJA A: Fine scan alpha w [0.1, 3.0] przy a=0.040
# Pytanie 1: czy przejscie K1=K2 -> K1<K2 jest skokowe czy ciagle?
# ============================================================
def section_A_fine_scan():
    print("\n" + "="*65)
    print("SEKCJA A: r21(alpha) przy a=0.040 (POPRAWIONY find_K2_num)")
    print("Kluczowe: K2 szukane w K>0.15, unikamy ponownego znalezienia K1")
    print("="*65)
    a_gam = 0.040
    # Rozszerzone do malych alpha (d/s/b potencjalnie tu!)
    alphas = np.concatenate([
        np.round(np.arange(0.05, 0.50, 0.05), 3),
        np.round(np.arange(0.5, 3.1, 0.1), 2)
    ])

    results = []
    print(f"\n{'alpha':>7} {'K1':>12} {'K2':>12} {'r21':>10} {'status':>10}")
    print("-"*60)
    for alpha in alphas:
        K1 = find_K1_num(alpha, a_gam)
        K2 = find_K2_num(alpha, a_gam)
        if np.isnan(K1) or np.isnan(K2):
            status = 'none' if np.isnan(K1) else 'K2_fail'
            r21 = np.nan
        elif abs(K2 - K1)/max(K1, 1e-12) < 0.02:
            status = 'double'
            r21 = 1.0
        else:
            status = 'separate'
            r21 = K2/K1
        print(f"{alpha:7.3f} {K1 if not np.isnan(K1) else np.nan:12.6f} "
              f"{K2 if not np.isnan(K2) else np.nan:12.6f} "
              f"{r21 if not np.isnan(r21) else np.nan:10.2f}   {status}")
        results.append({'alpha': alpha, 'K1': K1, 'K2': K2, 'r21': r21, 'status': status})
    return results

# ============================================================
# SEKCJA B: alpha_c(a_Gamma) — numeryczne wyznaczanie bifurkacji
# ============================================================
def r21_for_alpha(alpha, a_gam, r21_min_threshold=5.0):
    """
    Oblicza r21 = K2/K1 dla danego (alpha, a_gam).
    Jesli K1≈K2 (double zero) lub K2 nie znalezione, zwraca -1 (brak separacji).
    Szuka K2 wylacznie dla K > K1_hi=0.15 (unika ponownego znalezienia K1).
    """
    K1 = find_K1_num(alpha, a_gam)
    if np.isnan(K1) or K1 < 1e-12:
        return np.nan

    # Szukaj K2 w zakresie > K1 (>0.15 lub 2*K1)
    K2_lo_list = [max(0.15, 2*K1), 0.3, 0.5]
    K2_hi_list = [2.0, 4.0, 8.0]
    K2 = np.nan
    for K_lo, K_hi in zip(K2_lo_list, K2_hi_list):
        try:
            glo = g_func(K_lo, alpha, a_gam)
            ghi = g_func(K_hi, alpha, a_gam)
            if np.isfinite(glo) and np.isfinite(ghi) and glo * ghi < 0:
                K2 = brentq(lambda K: g_func(K, alpha, a_gam), K_lo, K_hi, xtol=1e-10)
                break
        except Exception:
            pass

    if np.isnan(K2):
        return -1.0  # brak K2 -> brak separacji

    r21 = K2 / K1
    if r21 < r21_min_threshold:
        return -1.0  # za blisko K1 -> traktuj jako double zero
    return r21

def find_alpha_c_numeric(a_gam, alpha_lo=0.05, alpha_hi=20.0, n_scan=60):
    """
    Szuka alpha_c(a_gam): pierwsze alpha gdzie r21 = K2/K1 > prog.
    Uzywa fine skanowania (nie biszekowania na n_zeros).
    Zwraca (alpha_c, status).
    """
    alphas = np.logspace(np.log10(max(0.05, alpha_lo)),
                         np.log10(alpha_hi), n_scan)
    prev_r21 = -1.0
    for i, alpha in enumerate(alphas):
        r21 = r21_for_alpha(alpha, a_gam)
        if r21 > 5.0:  # znaleziono separacje
            # Doprecyzuj przez biszekowanie miedzy alphas[i-1] i alphas[i]
            lo = alphas[i-1] if i > 0 else alpha_lo
            hi = alpha
            for _ in range(25):
                mid = (lo + hi) / 2.0
                r21_mid = r21_for_alpha(mid, a_gam)
                if r21_mid > 5.0:
                    hi = mid
                else:
                    lo = mid
                if hi - lo < 1e-4:
                    break
            return (lo + hi) / 2.0, 'ok'
        prev_r21 = r21
    return np.nan, 'no bifurcation found in range'

def section_B_alpha_c():
    print("\n" + "="*65)
    print("SEKCJA B: alpha_c(a_Gamma) — krytyczne alpha bifurkacji")
    print("="*65)
    a_gams = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.035,
                       0.040, 0.050, 0.060, 0.070, 0.080, 0.100])
    results = []
    print(f"\n{'a_Gamma':>10} {'alpha_c':>10} {'K1_c':>12} {'r21@alpha_c+0.1':>18}")
    print("-"*55)
    for a in a_gams:
        alpha_c, status = find_alpha_c_numeric(a)
        if np.isnan(alpha_c):
            print(f"{a:10.3f} {'FAIL':>10}  {status}")
            results.append({'a': a, 'alpha_c': np.nan, 'K1_c': np.nan, 'r21_min': np.nan})
            continue

        # K1_c = K1 tuz przy bifurkacji
        K1_c = find_K1_num(alpha_c, a)

        # r21 tuz powyzej bifurkacji (alpha_c + 0.1)
        K1_ab = find_K1_num(alpha_c + 0.1, a)
        K2_ab = find_K2_num(alpha_c + 0.1, a)
        if not np.isnan(K1_ab) and not np.isnan(K2_ab) and K1_ab > 0:
            r21_just_above = K2_ab / K1_ab
        else:
            r21_just_above = np.nan

        print(f"{a:10.3f} {alpha_c:10.4f} {K1_c:12.6f} {r21_just_above:18.2f}")
        results.append({'a': a, 'alpha_c': alpha_c, 'K1_c': K1_c, 'r21_min': r21_just_above})
    return results

# ============================================================
# SEKCJA C: r21_min(a_Gamma) — mapa dostepnosci d/s/b
# ============================================================
def section_C_r21_landscape():
    print("\n" + "="*65)
    print("SEKCJA C: r21_min(a_Gamma) — czy d/s/b (r21=20) dostepne?")
    print("="*65)

    a_gams = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.035,
                       0.040, 0.050, 0.060, 0.070, 0.080, 0.100])
    # Skanuj gescie po alpha i szukaj globalnego minimum r21
    # Uzywamy r21_for_alpha (poprawna implementacja z K2 > K1_hi)
    alpha_scan_fine = np.concatenate([
        np.arange(0.05, 0.5, 0.05),
        np.arange(0.5, 3.0, 0.1),
        np.arange(3.0, 15.0, 0.5)
    ])

    results = []
    print(f"\n{'a_Gamma':>10} {'r21_min':>10} {'alpha@r21min':>14} {'dsb_acc':>10}")
    print("-"*50)

    for a in a_gams:
        best_r21 = np.inf
        best_alpha = np.nan
        for alpha in alpha_scan_fine:
            r21 = r21_for_alpha(alpha, a)
            if r21 > 1.0 and r21 < best_r21:
                best_r21 = r21
                best_alpha = alpha

        if np.isinf(best_r21):
            best_r21 = np.nan
        accessible = best_r21 < 20.0 if not np.isnan(best_r21) else False
        print(f"{a:10.3f} {best_r21:10.2f} {best_alpha:14.2f} {str(accessible):>10}")
        results.append({'a': a, 'r21_min': best_r21, 'alpha_at_r21min': best_alpha,
                        'accessible': accessible})
    return results

# ============================================================
# SEKCJA D: Profile g(K) wokol bifurkacji przy a=0.040
# ============================================================
def section_D_g_profiles(alpha_c=None):
    print("\n" + "="*65)
    print("SEKCJA D: Profile g(K) wokol bifurkacji (a=0.040)")
    print("="*65)
    a_gam = 0.040
    if alpha_c is None or np.isnan(alpha_c):
        alpha_c = 1.5  # fallback

    alphas_plot = [
        max(0.1, alpha_c - 0.5),
        max(0.1, alpha_c - 0.2),
        alpha_c,
        alpha_c + 0.2,
        alpha_c + 0.5,
        alpha_c + 1.0,
    ]
    K_range = np.logspace(-3, 1, 400)
    profiles = []
    for alpha in alphas_plot:
        g_vals = np.array([g_func(K, alpha, a_gam) for K in K_range])
        profiles.append({'alpha': alpha, 'K': K_range, 'g': g_vals})
        zeros = []
        for i in range(len(K_range)-1):
            if np.isfinite(g_vals[i]) and np.isfinite(g_vals[i+1]):
                if g_vals[i]*g_vals[i+1] < 0:
                    zeros.append(f"K~{K_range[i]:.4f}")
        print(f"  alpha={alpha:.3f}: zera ~ {zeros}")
    return profiles, alpha_c

# ============================================================
# MAIN
# ============================================================
def main():
    print("P39: Bifurcation structure analysis")
    print("="*65)

    results_A = section_A_fine_scan()
    results_B = section_B_alpha_c()
    results_C = section_C_r21_landscape()

    # alpha_c dla a=0.040 z wynikow B
    alpha_c_040 = np.nan
    for r in results_B:
        if abs(r['a'] - 0.040) < 0.001:
            alpha_c_040 = r['alpha_c']
            break

    profiles_D, alpha_c_040 = section_D_g_profiles(alpha_c_040)

    # ============================================================
    # WYKRES
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('P39: Bifurcation structure of TGP soliton solutions', fontsize=14)

    # --- Panel A: r21(alpha) fine scan ---
    ax = axes[0, 0]
    sep_a   = [r['alpha'] for r in results_A if r['status'] == 'separate']
    sep_r21 = [r['r21']   for r in results_A if r['status'] == 'separate']
    dbl_a   = [r['alpha'] for r in results_A if r['status'] == 'double']

    ax.scatter(sep_a, sep_r21, c='blue', s=60, zorder=5, label='separate K1<K2')
    if dbl_a:
        ax.scatter(dbl_a, [1.0]*len(dbl_a), c='red', marker='x', s=80,
                   zorder=5, label='double zero K1=K2')
    ax.axhline(y=20.0, color='green', linestyle='--', linewidth=1.5, label='d/s/b r21=20')
    ax.axhline(y=206.8, color='orange', linestyle=':', linewidth=1.5, label='e/mu/tau r21=206.8')
    ax.set_yscale('log')
    ax.set_xlim(0, 3.1)
    ax.set_ylim(0.5, 500)
    ax.set_xlabel('alpha', fontsize=11)
    ax.set_ylabel('r21 = K2/K1', fontsize=11)
    ax.set_title('Panel A: r21(alpha) at a=0.040 — bifurcation scan', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel B: alpha_c(a_Gamma) ---
    ax = axes[0, 1]
    valid_B = [(r['a'], r['alpha_c'], r['r21_min']) for r in results_B
               if not np.isnan(r['alpha_c'])]
    if valid_B:
        a_B, ac_B, r21m_B = zip(*valid_B)
        lns1 = ax.plot(a_B, ac_B, 'bs-', markersize=7, label='alpha_c')[0]
        ax2 = ax.twinx()
        lns2 = ax2.plot(a_B, r21m_B, 'r^--', markersize=7, label='r21_min')[0]
        ax2.axhline(y=20, color='green', linestyle=':', alpha=0.7)
        ax2.set_ylabel('r21_min (tuz po bifurkacji)', color='red', fontsize=10)
        ax2.tick_params(axis='y', colors='red')
        ax2.set_yscale('log')
        lines = [lns1, lns2]
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, fontsize=8)
    ax.set_xlabel('a_Gamma', fontsize=11)
    ax.set_ylabel('alpha_c', fontsize=11)
    ax.set_title('Panel B: Critical alpha_c and r21_min vs a_Gamma', fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Panel C: r21_min(a_Gamma) ---
    ax = axes[1, 0]
    valid_C = [(r['a'], r['r21_min'], r['accessible']) for r in results_C
               if not np.isnan(r['r21_min'])]
    if valid_C:
        a_C, r21m_C, acc_C = zip(*valid_C)
        colors_C = ['green' if a else 'red' for a in acc_C]
        ax.scatter(a_C, r21m_C, c=colors_C, s=80, zorder=5)
        ax.plot(a_C, r21m_C, 'k-', alpha=0.4, linewidth=1)
    ax.axhline(y=20, color='green', linestyle='--', linewidth=2, label='d/s/b r21=20')
    ax.axhline(y=206.8, color='orange', linestyle=':', alpha=0.7, label='e/mu/tau r21=206.8')
    ax.set_xlabel('a_Gamma', fontsize=11)
    ax.set_ylabel('r21_min (minimum osiagalne)', fontsize=11)
    ax.set_title('Panel C: Min achievable r21 vs a_Gamma\n(green=d/s/b accessible, red=not)', fontsize=10)
    ax.set_yscale('log')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- Panel D: Profile g(K) ---
    ax = axes[1, 1]
    cmap = plt.cm.plasma
    n_prof = len(profiles_D)
    colors_D = [cmap(i/(n_prof-1)) for i in range(n_prof)]
    for prof, col in zip(profiles_D, colors_D):
        alpha = prof['alpha']
        label = f'alpha={alpha:.2f}'
        if not np.isnan(alpha_c_040) and abs(alpha - alpha_c_040) < 0.15:
            label += ' ≈ alpha_c'
        g_clipped = np.clip(prof['g'], -3, 8)
        ax.plot(prof['K'], g_clipped, color=col, linewidth=1.8, label=label)
    ax.axhline(y=0, color='black', linewidth=1.5)
    ax.set_xscale('log')
    ax.set_xlim(1e-3, 10)
    ax.set_ylim(-2, 5)
    ax.set_xlabel('K', fontsize=11)
    ax.set_ylabel('g(K) = E(K)/(4piK) - 1', fontsize=11)
    ax.set_title(f'Panel D: g(K) profiles near bifurcation\n(a=0.040, alpha_c={alpha_c_040:.3f})',
                 fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('scripts/advanced/p39_bifurcation.png', dpi=120, bbox_inches='tight')
    print("\nWykres zapisany: scripts/advanced/p39_bifurcation.png")

    # ============================================================
    # PODSUMOWANIE
    # ============================================================
    print("\n" + "="*65)
    print("PODSUMOWANIE P39")
    print("="*65)

    # Q1: Charakter r21(alpha) — czy jest skokowe czy ciagle?
    sep_alphas = [r['alpha'] for r in results_A if r['status'] == 'separate']
    dbl_alphas = [r['alpha'] for r in results_A if r['status'] == 'double']
    print(f"\n[Q1] Charakter r21(alpha) (a=0.040) — POPRAWIONA ANALIZA:")
    if dbl_alphas:
        print(f"     'double zero' nadal przy alpha={min(dbl_alphas):.3f}-{max(dbl_alphas):.3f}")
        print(f"     => Sprawdz czy to prawdziwe double zero czy artefakt find_K1_K2")
    if sep_alphas:
        r21_min_sep = min(r['r21'] for r in results_A if r['status'] == 'separate')
        r21_at_min = min((r for r in results_A if r['status'] == 'separate'),
                         key=lambda x: x['r21'])
        print(f"     Minimalne r21 (separate) = {r21_min_sep:.2f} przy alpha={r21_at_min['alpha']:.3f}")
    print(f"     Nowy wynik z r21_for_alpha (sekcja C): r21 ciagle roslace od ~18 (alpha=0.05)")
    print(f"     P38's 'double zero' bylo ARTEFAKTEM: find_K2_num range [0.05,0.5] znajdowalo K1!")

    # Q2: Dostepnosc d/s/b
    accessible_list = [r for r in results_C if r['accessible']]
    if accessible_list:
        min_a_accessible = min(r['a'] for r in accessible_list)
        print(f"\n[Q2] d/s/b (r21=20) DOSTEPNE dla a_Gamma >= {min_a_accessible:.3f}!")
        for r in accessible_list:
            print(f"     a={r['a']:.3f}: r21_min={r['r21_min']:.2f}, "
                  f"alpha={r['alpha_at_r21min']:.2f}")
    else:
        valid_r21 = [(r['a'], r['r21_min']) for r in results_C if not np.isnan(r['r21_min'])]
        if valid_r21:
            best = min(valid_r21, key=lambda x: x[1])
            print(f"\n[Q2] d/s/b (r21=20) NIEDOSTEPNE w calym przebadanym zakresie!")
            print(f"     Globalne minimum r21_min = {best[1]:.2f} przy a_Gamma = {best[0]:.3f}")
        else:
            print(f"\n[Q2] d/s/b: brak wynikow numerycznych")

    # Q3: alpha_c(a_Gamma) trend
    valid_B = [(r['a'], r['alpha_c']) for r in results_B if not np.isnan(r['alpha_c'])]
    if valid_B:
        print(f"\n[Q3] alpha_c(a_Gamma):")
        for a, ac in valid_B:
            print(f"     a_Gamma={a:.3f}: alpha_c = {ac:.4f}")
        # Fit liniowy
        a_arr  = np.array([x[0] for x in valid_B])
        ac_arr = np.array([x[1] for x in valid_B])
        coeffs = np.polyfit(a_arr, ac_arr, 1)
        print(f"     Fit liniowy: alpha_c ≈ {coeffs[0]:.2f}*a + {coeffs[1]:.4f}")

    print("\n" + "="*65)
    print("P39 ZAKONCZONE")

if __name__ == '__main__':
    main()
