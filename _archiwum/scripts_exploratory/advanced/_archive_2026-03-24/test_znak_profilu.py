"""
TGP — Test: profil Phi(r) = Phi0 +/- K*exp(-msp*r)/r
======================================================

Sprawdzamy ktory znak jest fizycznie poprawny i jak zmienia sie r21.

Uzasadnienie znaku MINUS:
  Phi0 to tlo WLACZNIE z wkladem tej czestkि.
  Odizolowany profil czastki = odchylenie od prawdziwego tla (bez niej).
  Czastka "zajmuje" pole => lokalne zageszczenie oznacza ze bez niej bylo go wiecej
  => profil czastki = dip (minus).

Konsekwencja:
  E_kin: identyczna dla obu znakow (zalezy od (dPhi/dr)^2)
  E_pot: rozni sie w czlonie kubicznym psi^3/3 (zmiana znaku)
"""
import numpy as np
from scipy.optimize import brentq

GAMMA = 1.0

def energy(K, Phi0, alpha, a_gam, znak=+1):
    """
    znak = +1: Phi(r) = Phi0 + K*exp(-msp*r)/r  [bump, stary wzor]
    znak = -1: Phi(r) = Phi0 - K*exp(-msp*r)/r  [dip,  poprawny fizycznie]
    """
    msp   = np.sqrt(max(GAMMA / Phi0, 1e-9))
    r_max = max(60.0 / msp, 20.0)
    r     = np.linspace(a_gam, r_max, 4000)

    phi  = Phi0 + znak * K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)          # klip: phi > 0 zawsze
    dphi = znak * K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / (Phi0 * phi)) * r**2, r)
    psi = phi / Phi0
    Ep  = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4) * r**2, r)
    return Ek + Ep


def find_crossings(alpha, a_gam, znak=+1, K_min=1e-4, K_max=5.0, N=800):
    K_arr = np.linspace(K_min, K_max, N)
    f_arr = np.array([
        energy(K, 1.0, alpha, a_gam, znak)/K - 4*np.pi for K in K_arr])
    roots = []
    for i in range(len(f_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                root = brentq(
                    lambda K: energy(K, 1.0, alpha, a_gam, znak)/K - 4*np.pi,
                    K_arr[i], K_arr[i+1], xtol=1e-9)
                roots.append(root)
            except Exception:
                pass
    return roots


# ======================================================================
# 1. KSZTALT E(K)/K DLA OBU ZNAKOW
# ======================================================================
print("=" * 70)
print("E(K)/K dla alpha=5.9, a_gam=0.05 — porownanie znakow")
print(f"  4*pi = {4*np.pi:.4f}")
print("=" * 70)
print(f"{'K':>8} {'E/K (+)':>12} {'E/K (-)':>12} {'roznica':>10}")
print("-" * 50)
alpha, a_gam = 5.9, 0.05
for K in [0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0]:
    ep = energy(K, 1.0, alpha, a_gam, +1) / K
    em = energy(K, 1.0, alpha, a_gam, -1) / K
    print(f"{K:>8.3f} {ep:>12.4f} {em:>12.4f} {em-ep:>+10.4f}")


# ======================================================================
# 2. PRZECIECIA: czy istnieja dwa dla znaku minus?
# ======================================================================
print()
print("=" * 70)
print("PRZECIECIA E(K)/K = 4*pi  (dwa znaki, rozne parametry)")
print("=" * 70)
print(f"{'alpha':>6} {'a_gam':>6} | "
      f"{'K1(+)':>9} {'K2(+)':>9} {'r21(+)':>9} | "
      f"{'K1(-)':>9} {'K2(-)':>9} {'r21(-)':>9}")
print("-" * 75)

cases = [
    (4.0, 0.03), (4.0, 0.05),
    (5.9, 0.03), (5.9, 0.05),
    (6.0, 0.03),
    (7.0, 0.03), (7.0, 0.05),
    (8.0, 0.03),
    (10.0, 0.03),
]

results = []
for alpha, a_gam in cases:
    rp = find_crossings(alpha, a_gam, znak=+1)
    rm = find_crossings(alpha, a_gam, znak=-1)

    def fmt_roots(r):
        if len(r) >= 2:
            K1, K2 = sorted(r)[:2]
            return f"{K1:>9.5f} {K2:>9.5f} {K2/K1:>9.1f}"
        elif len(r) == 1:
            return f"{'?':>9} {r[0]:>9.5f} {'?':>9}"
        return f"{'brak':>9} {'brak':>9} {'brak':>9}"

    print(f"{alpha:>6.1f} {a_gam:>6.3f} | {fmt_roots(rp)} | {fmt_roots(rm)}")
    results.append((alpha, a_gam, rp, rm))


# ======================================================================
# 3. ZMIANA r21 — jak rozni sie znak od znaku?
# ======================================================================
print()
print("=" * 70)
print("POROWNANIE r21: znak (+) vs znak (-)")
print("=" * 70)
print(f"{'alpha':>6} {'a_gam':>6} | {'r21(+)':>9} {'r21(-)':>9} {'roznica':>9}")
print("-" * 50)
for alpha, a_gam, rp, rm in results:
    r21p = sorted(rp)[1]/sorted(rp)[0] if len(rp)>=2 else float('nan')
    r21m = sorted(rm)[1]/sorted(rm)[0] if len(rm)>=2 else float('nan')
    diff = 100*(r21m - r21p)/r21p if (np.isfinite(r21p) and np.isfinite(r21m)) else float('nan')
    s_d = f"{diff:>+9.1f}%" if np.isfinite(diff) else "N/A"
    sp = f"{r21p:>9.1f}" if np.isfinite(r21p) else "N/A"
    sm = f"{r21m:>9.1f}" if np.isfinite(r21m) else "N/A"
    print(f"{alpha:>6.1f} {a_gam:>6.3f} | {sp} {sm} {s_d}")


# ======================================================================
# 4. SKAN: jakie (alpha, a_gam) daja r21(-) ~ 207?
# ======================================================================
print()
print("=" * 70)
print("SKAN: kombinacje (alpha, a_gam) => r21(-) ~ 207  [znak minus]")
print("=" * 70)
print(f"{'alpha':>7} {'a_gam':>7} {'r21(-)':>9}")
print("-" * 30)
for alpha in np.arange(3.0, 12.1, 0.5):
    for a_gam in [0.02, 0.03, 0.05, 0.08, 0.10]:
        rm = find_crossings(alpha, a_gam, znak=-1)
        if len(rm) >= 2:
            r21m = sorted(rm)[1] / sorted(rm)[0]
            if 170 < r21m < 250:
                print(f"{alpha:>7.1f} {a_gam:>7.3f} {r21m:>9.1f}")
