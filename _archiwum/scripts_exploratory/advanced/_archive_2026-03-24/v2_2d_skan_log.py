"""
v2_2d_skan_log.py

Cel: Znalezc (alpha, a_gam, lambda) takie ze:
  r21 = K2/K1 = 207 (+/-2%)
  r31 = K3/K1 = 3477 (+/-2%)

Poprzedni wynik (log grid): alpha=7.0, a_gam=0.030, lam*=3.882e-6
  => r21=211.4 (+2.1%), r31=3477 OK
  Trzeba zmienic (alpha, a_gam) zeby zmniejszyc r21.

Podejscie:
  1. Dla kazdego (alpha, a_gam) w siatce 2D:
     a) Bisektuj lam* takie ze r31=3477
     b) Oblicz r21 przy tym lam*
     c) Sprawdz czy r21 ~ 207

Oczekiwanie: r21 maleje gdy... sprawdzimy numerycznie.

UWAGA: uzywamy WYLACZNIE siatki logarytmicznej (N=1500) dla szybkosci.
Potwierdzenie zbieznosci N bedzie osobno dla najlepszego punktu.
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0

def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_log(K, alpha, a_gam, lam, N=1500):
    r_max = 60.0
    t   = np.linspace(0, 1, N)
    r   = a_gam * (r_max/a_gam)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    psi  = phi
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam, N=1500):
    E = energy_log(K, alpha, a_gam, lam, N)
    return E/K - 4*np.pi

def find_crossings(alpha, a_gam, lam, K_max=200.0, N_scan=800):
    """Znajdz zera g(K) z siatka log."""
    K_arr = np.concatenate([
        np.linspace(1e-4, 0.5, 150),
        np.linspace(0.5, 10.0, 200),
        np.linspace(10.0, K_max, 450),
    ])
    K_arr = np.unique(K_arr)
    f_arr = np.array([g_func(K, alpha, a_gam, lam) for K in K_arr])
    roots = []
    for i in range(len(f_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                root = brentq(
                    lambda K: g_func(K, alpha, a_gam, lam),
                    K_arr[i], K_arr[i+1], xtol=1e-8, maxiter=60)
                roots.append(root)
            except Exception:
                pass
    return sorted(roots)

def find_lam_for_r31(alpha, a_gam, target_r31=3477.0):
    """
    Bisektuj lam* takie ze K3/K1 = target_r31.
    Najpierw skanujemy lam zeby znalezc przedzialy.
    """
    # Wstepny skan lam
    lam_candidates = [
        1e-7, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6,
        7e-6, 1e-5, 2e-5, 5e-5
    ]
    r31_vals = []
    for lam in lam_candidates:
        roots = find_crossings(alpha, a_gam, lam)
        if len(roots) >= 3:
            r31 = roots[2]/roots[0]
            r31_vals.append((lam, r31))

    if len(r31_vals) < 2:
        return None

    # Znajdz przedzial gdzie r31 przechodzi przez target_r31
    lam_lo, lam_hi = None, None
    for i in range(len(r31_vals)-1):
        la, ra = r31_vals[i]
        lb, rb = r31_vals[i+1]
        if (ra - target_r31) * (rb - target_r31) < 0:
            lam_lo, lam_hi = la, lb
            break

    if lam_lo is None:
        return None

    # Bisekcja
    def residual(lam):
        roots = find_crossings(alpha, a_gam, lam)
        if len(roots) < 3:
            return -target_r31
        return roots[2]/roots[0] - target_r31

    try:
        lam_star = brentq(residual, lam_lo, lam_hi,
                          xtol=lam_hi*2e-3, maxiter=40)
    except Exception:
        return None

    roots = find_crossings(alpha, a_gam, lam_star)
    if len(roots) < 3:
        return None
    return lam_star, roots[0], roots[1], roots[2]

# ============================================================
# KROK 1: Skan 2D (alpha, a_gam)
# ============================================================
TARGET_R21 = 207.0
TARGET_R31 = 3477.0

print("=" * 80)
print("2D SKAN (alpha, a_gam): szukamy r21=207 i r31=3477 (siatka log)")
print(f"{'alpha':>7} {'a_gam':>7} | {'lam*':>12} {'K1':>9} {'K2':>9} {'K3':>9}"
      f" | {'r21':>8} {'r31':>8} | {'blad_r21':>9}")
print("-" * 88)

# Siatka (alpha, a_gam)
# Wiemy ze alpha=7.0, a_gam=0.030 daje r21=211.4 (+2.1%)
# r21 maleje gdy... zbadamy
alpha_vals = [5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0]
a_gam_vals = [0.020, 0.025, 0.030, 0.035, 0.040, 0.050]

results = []

for alpha in alpha_vals:
    for a_gam in a_gam_vals:
        res = find_lam_for_r31(alpha, a_gam, TARGET_R31)
        if res is None:
            print(f"{alpha:>7.1f} {a_gam:>7.3f} | --- brak rozwiazania r31=3477 ---")
            continue
        lam_s, K1, K2, K3 = res
        r21 = K2/K1
        r31 = K3/K1
        err21 = 100*(r21 - TARGET_R21)/TARGET_R21
        print(f"{alpha:>7.1f} {a_gam:>7.3f} | {lam_s:>12.4e} {K1:>9.5f} {K2:>9.5f} {K3:>9.4f}"
              f" | {r21:>8.2f} {r31:>8.2f} | {err21:>+8.3f}%")
        results.append({
            'alpha': alpha, 'a_gam': a_gam, 'lam': lam_s,
            'K1': K1, 'K2': K2, 'K3': K3,
            'r21': r21, 'r31': r31, 'err21': err21
        })

print()

# ============================================================
# KROK 2: Najlepszy wynik
# ============================================================
if results:
    best = sorted(results, key=lambda x: abs(x['err21']))[0]
    print("=" * 80)
    print(f"NAJLEPSZY WYNIK (minimalne |blad r21|):")
    print(f"  alpha   = {best['alpha']}")
    print(f"  a_gam   = {best['a_gam']}")
    print(f"  lam*    = {best['lam']:.6e}")
    print(f"  K1={best['K1']:.5f}, K2={best['K2']:.5f}, K3={best['K3']:.4f}")
    print(f"  r21     = {best['r21']:.2f} (blad: {best['err21']:+.2f}%)")
    print(f"  r31     = {best['r31']:.2f}")
    psi_c = 1.0 + best['K3']*np.exp(-best['a_gam'])/best['a_gam']
    psi_n = 1.0/np.sqrt(best['lam'])
    print(f"  psi_core(M3) ~ {psi_c:.1f}")
    print(f"  psi_new = 1/sqrt(lam) = {psi_n:.1f}")
    print(f"  psi_core/psi_new = {psi_c/psi_n:.3f}")
    print()

    # ============================================================
    # KROK 3: Bisekcja po alpha lub a_gam zeby r21=207 dokladnie
    # ============================================================
    print("=" * 80)
    print("KROK 3: Dokladna bisekcja (trzymamy a_gam, skanujemy alpha dla r21=207)")
    print()

    # Znajdz przedzial alpha gdzie r21 przechodzi przez 207
    # Uzyjemy wynikow z kroku 1
    a_gam_best = best['a_gam']
    alpha_r21 = [(r['alpha'], r['r21']) for r in results
                 if abs(r['a_gam'] - a_gam_best) < 1e-4 and r.get('r21')]
    alpha_r21.sort(key=lambda x: x[0])

    print(f"Dla a_gam={a_gam_best:.3f}: r21 vs alpha:")
    for a, r in alpha_r21:
        print(f"  alpha={a:.1f}: r21={r:.2f}")

    # Znajdz przedzial
    al_lo, al_hi = None, None
    for i in range(len(alpha_r21)-1):
        a1, r1 = alpha_r21[i]
        a2, r2 = alpha_r21[i+1]
        if (r1 - TARGET_R21)*(r2 - TARGET_R21) < 0:
            al_lo, al_hi = a1, a2
            break

    if al_lo is not None:
        print(f"\nZmiana znaku (r21-207) w przedziale alpha=[{al_lo:.2f}, {al_hi:.2f}]")
        print("Bisekcja po alpha...")

        def residual_r21(alpha):
            res = find_lam_for_r31(alpha, a_gam_best, TARGET_R31)
            if res is None:
                return TARGET_R21  # neutralna wartosc
            _, K1, K2, K3 = res
            return K2/K1 - TARGET_R21

        try:
            alpha_star = brentq(residual_r21, al_lo, al_hi,
                                xtol=0.01, maxiter=30)
            res = find_lam_for_r31(alpha_star, a_gam_best, TARGET_R31)
            if res:
                lam_s, K1, K2, K3 = res
                r21 = K2/K1
                r31 = K3/K1
                err21 = 100*(r21-TARGET_R21)/TARGET_R21
                print(f"\n==> alpha* = {alpha_star:.4f}, a_gam = {a_gam_best:.3f}")
                print(f"    lam*   = {lam_s:.6e}")
                print(f"    K1={K1:.5f}, K2={K2:.5f}, K3={K3:.4f}")
                print(f"    r21 = {r21:.3f} (blad: {err21:+.3f}%)")
                print(f"    r31 = {r31:.2f}")
                psi_c = 1.0 + K3*np.exp(-a_gam_best)/a_gam_best
                psi_n = 1.0/np.sqrt(lam_s)
                print(f"    psi_core(M3) ~ {psi_c:.1f}")
                print(f"    psi_new = {psi_n:.1f}")
                print(f"    psi_core/psi_new = {psi_c/psi_n:.3f}")
        except Exception as ex:
            print(f"Bisekcja nie zbiegla: {ex}")
    else:
        print("\nBrak zmiany znaku r21-207 dla a_gam_best — probujemy inna a_gam")
        # Pobierz wszystkie wyniki blisko r21=207
        close = sorted(results, key=lambda x: abs(x['r21']-207))[:5]
        print("\nNajblizsze punkty do r21=207:")
        for r in close:
            print(f"  alpha={r['alpha']:.1f}, a_gam={r['a_gam']:.3f}: "
                  f"r21={r['r21']:.2f} (blad {r['err21']:+.2f}%)")

print()
print("=" * 80)
print("Koniec v2_2d_skan_log.py")
