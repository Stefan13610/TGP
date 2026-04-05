"""
v2_bisekcja_alpha.py

Na podstawie v2_2d_skan_log.py znamy przedziaiy, gdzie r21 przechodzi
przez 207 (przy r31=3477 wymuszonym przez lam*).

Przedzials ze zmiana znaku r21-207 (trzymajac a_gam stale):
  a_gam=0.025: alpha in [5.5, 6.0]  (r21: 191.58 -> 210.17)
  a_gam=0.030: alpha in [6.5, 7.0]  (r21: 194.79 -> 211.40)
  a_gam=0.035: alpha in [7.5, 8.0]  (r21: 199.56 -> 214.75)
  a_gam=0.040: alpha in [8.0, 9.0]  (r21: 191.38 -> 219.39)
  a_gam=0.050: alpha in [10.0, 11.0] (r21: 205.88 -> ?)

Dla kazdego a_gam: bisektujemy alpha* takie ze r21=207.
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
    return energy_log(K, alpha, a_gam, lam, N)/K - 4*np.pi

def find_crossings(alpha, a_gam, lam, K_max=150.0):
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
    """Znajdz lam* takie ze r31=target_r31 (bisekcja z wstepnym skanem)."""
    lam_cands = [1e-7, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 7e-6, 1e-5, 2e-5]
    r31_pts = []
    for lam in lam_cands:
        roots = find_crossings(alpha, a_gam, lam)
        if len(roots) >= 3:
            r31_pts.append((lam, roots[2]/roots[0]))

    lam_lo = lam_hi = None
    for i in range(len(r31_pts)-1):
        la, ra = r31_pts[i]
        lb, rb = r31_pts[i+1]
        if (ra - target_r31)*(rb - target_r31) < 0:
            lam_lo, lam_hi = la, lb
            break

    if lam_lo is None:
        return None

    def res(lam):
        roots = find_crossings(alpha, a_gam, lam)
        if len(roots) < 3:
            return -target_r31
        return roots[2]/roots[0] - target_r31

    try:
        lam_s = brentq(res, lam_lo, lam_hi, xtol=lam_hi*1e-3, maxiter=40)
    except Exception:
        return None

    roots = find_crossings(alpha, a_gam, lam_s)
    if len(roots) < 3:
        return None
    return lam_s, roots[0], roots[1], roots[2]

TARGET_R21 = 207.0
TARGET_R31 = 3477.0

# Przedzials ze zmiana znaku r21-207
brackets = [
    (0.025, 5.5, 6.0),
    (0.030, 6.5, 7.0),
    (0.035, 7.5, 8.0),
    (0.040, 8.0, 9.0),
]

print("=" * 80)
print("DOKLADNA BISEKCJA alpha* dla r21=207, r31=3477 (siatka log, N=1500)")
print()

final_results = []

for a_gam, al_lo, al_hi in brackets:
    print(f"--- a_gam = {a_gam:.3f}, alpha in [{al_lo:.1f}, {al_hi:.1f}] ---")

    def residual_r21(alpha):
        res = find_lam_for_r31(alpha, a_gam, TARGET_R31)
        if res is None:
            print(f"  [warning] alpha={alpha:.4f}: brak r31")
            return TARGET_R21
        _, K1, K2, K3 = res
        r21 = K2/K1
        return r21 - TARGET_R21

    # Weryfikuj wartosci na kraincach
    r_lo = residual_r21(al_lo)
    r_hi = residual_r21(al_hi)
    print(f"  r21-207 @ alpha={al_lo:.1f}: {r_lo:+.2f}")
    print(f"  r21-207 @ alpha={al_hi:.1f}: {r_hi:+.2f}")

    if r_lo * r_hi >= 0:
        print("  [skip] brak zmiany znaku!")
        print()
        continue

    try:
        alpha_star = brentq(residual_r21, al_lo, al_hi, xtol=0.005, maxiter=25)
        res = find_lam_for_r31(alpha_star, a_gam, TARGET_R31)
        if res:
            lam_s, K1, K2, K3 = res
            r21 = K2/K1
            r31 = K3/K1
            err21 = 100*(r21 - TARGET_R21)/TARGET_R21
            psi_c = 1.0 + K3*np.exp(-a_gam)/a_gam
            psi_n = 1.0/np.sqrt(lam_s)
            print(f"  => alpha* = {alpha_star:.4f}, a_gam = {a_gam:.3f}")
            print(f"     lam*   = {lam_s:.6e}")
            print(f"     K1={K1:.5f}, K2={K2:.5f}, K3={K3:.4f}")
            print(f"     r21 = {r21:.3f} (blad: {err21:+.3f}%)")
            print(f"     r31 = {r31:.2f}")
            print(f"     psi_core ~ {psi_c:.1f}, psi_new ~ {psi_n:.1f}")
            print(f"     psi_core/psi_new = {psi_c/psi_n:.3f}")
            final_results.append({
                'alpha': alpha_star, 'a_gam': a_gam, 'lam': lam_s,
                'K1': K1, 'K2': K2, 'K3': K3,
                'r21': r21, 'r31': r31, 'err21': err21,
                'psi_c': psi_c, 'psi_n': psi_n
            })
    except Exception as ex:
        print(f"  [blad bisekcji]: {ex}")
    print()

# Dodatkowy skan dla a_gam=0.050 (potrzeba alpha > 10)
print("--- a_gam = 0.050, szukamy alpha > 10 ---")
for alpha_test in [10.5, 11.0, 11.5, 12.0]:
    res = find_lam_for_r31(alpha_test, 0.050, TARGET_R31)
    if res:
        lam_s, K1, K2, K3 = res
        r21 = K2/K1
        err21 = 100*(r21 - TARGET_R21)/TARGET_R21
        print(f"  alpha={alpha_test:.1f}: r21={r21:.2f} ({err21:+.2f}%)")
print()

print("=" * 80)
print("PODSUMOWANIE: Znalezione rozwiazania (r21=207, r31=3477)")
print(f"{'alpha':>8} {'a_gam':>7} {'lam*':>12} | {'r21':>8} {'r31':>8} | "
      f"{'psi_c':>8} {'psi_n':>8} | {'blad_r21':>9}")
print("-" * 80)
for r in final_results:
    print(f"{r['alpha']:>8.4f} {r['a_gam']:>7.3f} {r['lam']:>12.4e} | "
          f"{r['r21']:>8.3f} {r['r31']:>8.1f} | "
          f"{r['psi_c']:>8.1f} {r['psi_n']:>8.1f} | {r['err21']:>+8.3f}%")

print()
print("=" * 80)
print("Koniec v2_bisekcja_alpha.py")
