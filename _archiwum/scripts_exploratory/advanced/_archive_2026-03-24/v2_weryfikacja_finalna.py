"""
v2_weryfikacja_finalna.py

Weryfikacja zbieznosci dla znalezionych rozwiazAN r21=207, r31=3477.

Znalezione rozwiazania (z v2_bisekcja_alpha.py):
  1. alpha=5.9148, a_gam=0.025, lam*=2.883e-6
  2. alpha=6.8675, a_gam=0.030, lam*=3.743e-6
  3. alpha=7.7449, a_gam=0.035, lam*=4.618e-6
  4. alpha=8.5616, a_gam=0.040, lam*=5.501e-6  <= r21=206.999 (najdokl.)

Sprawdzamy:
  a) Zbieznosc g(K) z N dla K1, K2, K3
  b) Stabilnosc zer (czy te same K przy N=500 vs N=8000?)
  c) Wymiar energetyczny: E/K -> 4*pi
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0

def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_log(K, alpha, a_gam, lam, N=2000):
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

def g_func(K, alpha, a_gam, lam, N=2000):
    return energy_log(K, alpha, a_gam, lam, N)/K - 4*np.pi

# Rozwiazania z v2_bisekcja_alpha.py
solutions = [
    {'alpha': 5.9148, 'a_gam': 0.025, 'lam': 2.883172e-06,
     'K1': 0.00853, 'K2': 1.76610, 'K3': 29.6702, 'r21': 206.972, 'r31': 3477.11},
    {'alpha': 6.8675, 'a_gam': 0.030, 'lam': 3.742489e-06,
     'K1': 0.00896, 'K2': 1.85472, 'K3': 31.1596, 'r21': 206.969, 'r31': 3477.12},
    {'alpha': 7.7449, 'a_gam': 0.035, 'lam': 4.618392e-06,
     'K1': 0.00939, 'K2': 1.94375, 'K3': 32.6553, 'r21': 206.967, 'r31': 3477.08},
    {'alpha': 8.5616, 'a_gam': 0.040, 'lam': 5.501357e-06,
     'K1': 0.00982, 'K2': 2.03272, 'K3': 34.1440, 'r21': 206.999, 'r31': 3477.00},
]

N_vals = [500, 1000, 2000, 4000, 8000, 16000]

print("=" * 80)
print("WERYFIKACJA ZBIEZNOSCI g(K) = E(K)/K - 4*pi")
print("Cel: 4*pi =", 4*np.pi)
print()

for s in solutions:
    alpha, a_gam, lam = s['alpha'], s['a_gam'], s['lam']
    K1, K2, K3 = s['K1'], s['K2'], s['K3']
    print("=" * 80)
    print(f"alpha={alpha:.4f}, a_gam={a_gam:.3f}, lam*={lam:.6e}")
    print(f"K1={K1:.5f}, K2={K2:.5f}, K3={K3:.4f}")
    print(f"r21={s['r21']:.3f}, r31={s['r31']:.2f}")
    print()

    for K_label, K in [("K1", K1), ("K2", K2), ("K3", K3)]:
        print(f"  g(K={K:.5f}) = E/K - 4*pi  [{K_label}]:")
        print(f"  {'N':>7} {'Ek':>12} {'Ep':>12} {'g':>14} {'delta':>10}")
        g_prev = None
        for N in N_vals:
            r_max = 60.0
            t   = np.linspace(0, 1, N)
            r   = a_gam * (r_max/a_gam)**t
            phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
            dphi = K*np.exp(-r)*(-r - 1.0)/r**2
            psi  = phi
            V1   = V_mod(1.0, lam)
            Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
            Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
            g = (Ek+Ep)/K - 4*np.pi
            delta = f"{100*(g-g_prev)/abs(g_prev+1e-30):+.3f}%" if g_prev is not None else ""
            print(f"  {N:>7} {Ek:>12.4e} {Ep:>12.4e} {g:>14.6e} {delta:>10}")
            g_prev = g
        print()

    print()

# ============================================================
# Zbiorczy test: czy zera K sa stabilne przy wyzszym N?
# ============================================================
print("=" * 80)
print("STABILNOSC ZER: refine K1, K2, K3 dla alpha=8.5616, a_gam=0.040, lam=5.501e-6")
print("przy N=2000 vs N=8000")
print()

s = solutions[3]  # Najdokladniejszy: alpha=8.5616
alpha, a_gam, lam = s['alpha'], s['a_gam'], s['lam']
K1_approx, K2_approx, K3_approx = s['K1'], s['K2'], s['K3']

def find_crossings_N(alpha, a_gam, lam, N, K_max=100.0):
    K_arr = np.concatenate([
        np.linspace(1e-4, 0.5, 150),
        np.linspace(0.5, 10.0, 200),
        np.linspace(10.0, K_max, 450),
    ])
    K_arr = np.unique(K_arr)
    f_arr = np.array([g_func(K, alpha, a_gam, lam, N) for K in K_arr])
    roots = []
    for i in range(len(f_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                root = brentq(
                    lambda K: g_func(K, alpha, a_gam, lam, N),
                    K_arr[i], K_arr[i+1], xtol=1e-9, maxiter=60)
                roots.append(root)
            except Exception:
                pass
    return sorted(roots)

for N in [2000, 4000, 8000]:
    roots = find_crossings_N(alpha, a_gam, lam, N)
    if len(roots) >= 3:
        K1, K2, K3 = roots[0], roots[1], roots[2]
        r21 = K2/K1
        r31 = K3/K1
        print(f"N={N:5d}: K1={K1:.6f}, K2={K2:.6f}, K3={K3:.5f} | "
              f"r21={r21:.4f}, r31={r31:.2f}")
    else:
        print(f"N={N:5d}: znaleziono {len(roots)} zer (brak 3. zera?)")

print()
print("=" * 80)
print("WYNIK FINALNY:")
print()
print("Rodzina rozwiAzan (r21=207, r31=3477) z lam*(psi-1)^6/6:")
print()
print(f"  {'alpha':>8} {'a_gam':>7} {'lam*':>12} | "
      f"{'r21':>8} {'r31':>8} | {'psi_c/psi_n':>12}")
for s in solutions:
    psi_c = 1.0 + s['K3']*np.exp(-s['a_gam'])/s['a_gam']
    psi_n = 1.0/np.sqrt(s['lam'])
    print(f"  {s['alpha']:>8.4f} {s['a_gam']:>7.3f} {s['lam']:>12.4e} | "
          f"{s['r21']:>8.3f} {s['r31']:>8.2f} | {psi_c/psi_n:>12.3f}")

print()
print("Obserwacja: psi_core/psi_new ~ 1.93-1.97 (prawie stalal!)")
print("  To sugeruje geometryczna relacje: rdzen solitonu ~ 2 x minimum potencjalu")
