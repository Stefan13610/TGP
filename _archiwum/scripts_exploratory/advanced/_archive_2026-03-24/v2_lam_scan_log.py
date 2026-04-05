"""
v2_lam_scan_log.py

Problem: Dla lam=5.514e-7 (z v2_fine_scan) r31~9217 >> 3477.
Musimy znalezc WIEKSZE lam (bo wieksze lam => mniejsze K3 => mniejsze r31).

Skanujemy lam od 5.514e-7 wzwyz, szukamy gdzie r31 spada do 3477.
Uzywamy siatki logarytmicznej.
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0

def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_log(K, alpha, a_gam, lam, N=1500):
    r_max = max(60.0, 20.0)
    t     = np.linspace(0, 1, N)
    r     = a_gam * (r_max/a_gam)**t
    phi   = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi  = K*np.exp(-r)*(-r - 1.0)/r**2
    psi   = phi
    V1    = V_mod(1.0, lam)
    Ek    = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep    = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek + Ep

def find_crossings_log(alpha, a_gam, lam, K_max=300.0, N_scan=800):
    K_arr = np.concatenate([
        np.linspace(1e-4, 0.5, 200),
        np.linspace(0.5, 10.0, 200),
        np.linspace(10.0, K_max, 400),
    ])
    K_arr = np.unique(K_arr)
    f_arr = np.array([energy_log(K,alpha,a_gam,lam)/K - 4*np.pi for K in K_arr])
    roots = []
    for i in range(len(f_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                root = brentq(
                    lambda K: energy_log(K,alpha,a_gam,lam)/K-4*np.pi,
                    K_arr[i], K_arr[i+1], xtol=1e-9)
                roots.append(root)
            except Exception:
                pass
    return sorted(roots)

ALPHA = 7.0
A_GAM = 0.030

# ── Skan lam od malych do duzych ───────────────────────────────────────────
print("=" * 65)
print("SKAN lambda (alpha=7.0, a_gam=0.030, siatka log)")
print(f"{'lam':>12} {'#zer':>5} {'K1':>9} {'K2':>9} {'K3':>9}"
      f" | {'r21':>8} {'r31':>8}")
print("-" * 65)

lam_vals = [
    1e-8, 5e-8,
    1e-7, 2e-7, 5e-7, 5.514e-7,
    1e-6, 2e-6, 5e-6,
    1e-5, 2e-5, 5e-5,
]

results_scan = []
for lam in lam_vals:
    roots = find_crossings_log(ALPHA, A_GAM, lam, K_max=300.0)
    n = len(roots)
    if n < 2:
        print(f"{lam:>12.4e} {n:>5} --- za malo zer ---")
        results_scan.append({'lam': lam, 'n': n, 'r21': None, 'r31': None})
        continue
    K1, K2 = roots[0], roots[1]
    r21 = K2/K1
    if n >= 3:
        K3  = roots[2]
        r31 = K3/K1
        print(f"{lam:>12.4e} {n:>5} {K1:>9.5f} {K2:>9.5f} {K3:>9.4f}"
              f" | {r21:>8.2f} {r31:>8.2f}")
        results_scan.append({'lam': lam, 'n': n, 'r21': r21, 'r31': r31,
                             'K1': K1, 'K2': K2, 'K3': K3})
    else:
        print(f"{lam:>12.4e} {n:>5} {K1:>9.5f} {K2:>9.5f} {'---':>9}"
              f" | {r21:>8.2f} {'---':>8}")
        results_scan.append({'lam': lam, 'n': n, 'r21': r21, 'r31': None})

print()
print("Kluczowe obserwacje:")
r31_vals = [(r['lam'], r['r31']) for r in results_scan if r.get('r31')]
for lam, r31 in r31_vals:
    print(f"  lam={lam:.4e}: r31={r31:.1f}")

print()
print("=" * 65)
print("SZUKAMY lam* takie ze r31=3477 (siatka log)")

# Znajdz przedzial [lam_lo, lam_hi] gdzie r31 przechodzi przez 3477
lam_lo_found, lam_hi_found = None, None
for i in range(len(r31_vals)-1):
    lam_a, r31_a = r31_vals[i]
    lam_b, r31_b = r31_vals[i+1]
    if (r31_a - 3477) * (r31_b - 3477) < 0:
        lam_lo_found, lam_hi_found = lam_a, lam_b
        print(f"Przedzial: lam=[{lam_a:.4e}, {lam_b:.4e}], r31=[{r31_a:.1f}, {r31_b:.1f}]")

if lam_lo_found is not None:
    def residual_r31(lam):
        roots = find_crossings_log(ALPHA, A_GAM, lam)
        if len(roots) < 3:
            return -3477.0
        return roots[2]/roots[0] - 3477.0

    print("Bisekcja...")
    lam_star = brentq(residual_r31, lam_lo_found, lam_hi_found,
                      xtol=lam_hi_found*1e-4, maxiter=50)
    roots = find_crossings_log(ALPHA, A_GAM, lam_star)
    if len(roots) >= 3:
        K1, K2, K3 = roots[0], roots[1], roots[2]
        r21 = K2/K1
        r31 = K3/K1
        err21 = 100*(r21-207)/207
        print(f"  lam* = {lam_star:.6e}")
        print(f"  K1={K1:.5f}, K2={K2:.5f}, K3={K3:.4f}")
        print(f"  r21 = {r21:.2f} (err {err21:+.2f}% vs 207)")
        print(f"  r31 = {r31:.2f}")
        psi_c = 1.0 + K3*np.exp(-A_GAM)/A_GAM
        psi_n = 1.0/np.sqrt(lam_star)
        print(f"  psi_core(M3) = {psi_c:.1f}")
        print(f"  psi_new = {psi_n:.1f}")
        print(f"  psi_core/psi_new = {psi_c/psi_n:.3f}")

        # Test zbieznosci
        print()
        print(f"TEST ZBIEZNOSCI E(K3={K3:.4f}) z siatka log:")
        print(f"{'N':>6} {'Etot':>14} {'Etot/K':>14}")
        for N in [500, 1000, 2000, 4000, 8000]:
            E = energy_log(K3, ALPHA, A_GAM, lam_star, N=N)
            print(f"{N:>6} {E:>14.4e} {E/K3:>14.6f}")
        print(f"  Cel: 4*pi = {4*np.pi:.6f}")
else:
    print("Nie znaleziono przedzials z r31=3477 w zbadanym zakresie lam.")
    print("Potrzeba wiekszyego skanu lam.")

# ── Dodatkowy skan dla roznych (alpha, a_gam) ─────────────────────────────
print()
print("=" * 65)
print("SKAN dla roznych (alpha, a_gam): lam=5.514e-7 (siatka log)")
print(f"{'alpha':>6} {'a_gam':>6} | {'r21':>8} {'r31':>8}")
print("-" * 35)
for alpha, a_gam in [(5.9,0.03),(5.9,0.05),(7.0,0.03),(8.0,0.035),(10.0,0.05)]:
    roots = find_crossings_log(alpha, a_gam, 5.514e-7, K_max=300.0)
    if len(roots) >= 3:
        r21 = roots[1]/roots[0]
        r31 = roots[2]/roots[0]
        print(f"{alpha:>6.1f} {a_gam:>6.3f} | {r21:>8.2f} {r31:>8.2f}")
    elif len(roots) == 2:
        print(f"{alpha:>6.1f} {a_gam:>6.3f} | {roots[1]/roots[0]:>8.2f} {'---':>8}")
    else:
        print(f"{alpha:>6.1f} {a_gam:>6.3f} | ---")
