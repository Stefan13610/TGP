"""
TGP — Wyznaczenie dwoch dyskretnych mas z warunkow samospojnosci
Parametry: a_Gam=0.1, alpha=1.0 (daje 2 przeciecia)
Sprawdzamy stosunek m2/m1 i szukamy warunkow dla ratio = 207
"""
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHI0 = 1.0; GAMMA = 1.0; Q = 1.0; M_SP = 1.0
LAM = 5.514e-7  # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)

def E_over_K(K, a_gam, alpha):
    r = np.linspace(a_gam, 25.0, 1500)
    phi  = PHI0 + K * np.exp(-M_SP*r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-M_SP*r) * (-M_SP*r - 1.0) / r**2
    E_k  = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1 + alpha/(PHI0*phi)) * r**2, r)
    V    = GAMMA/3*phi**3 - GAMMA/4*phi**4 - (GAMMA/3 - GAMMA/4) + LAM/6*(phi-1.0)**6
    E_p  = 4*np.pi * np.trapezoid(V * r**2, r)
    return (E_k + E_p) / K

def find_crossings(a_gam, alpha, K_min=0.02, K_max=15.0, N=200):
    Lambda = 4*np.pi / Q
    K_vals = np.linspace(K_min, K_max, N)
    f = np.array([E_over_K(K, a_gam, alpha) - Lambda for K in K_vals])
    crossings = []
    for i in range(len(f)-1):
        if np.isfinite(f[i]) and np.isfinite(f[i+1]) and f[i]*f[i+1] < 0:
            try:
                K_c = brentq(
                    lambda K: E_over_K(K, a_gam, alpha) - Lambda,
                    K_vals[i], K_vals[i+1], xtol=1e-7)
                crossings.append(K_c)
            except Exception:
                pass
    return crossings, K_vals, f

Lambda = 4*np.pi
print("="*65)
print("Szczegolowy skan: masy dla roznych (a_Gam, alpha)")
print(f"Lambda = {Lambda:.4f}")
print("="*65)
print(f"{'a_Gam':>6} {'alpha':>6} {'K1':>8} {'K2':>8} {'m1':>8} {'m2':>8} {'m2/m1':>8}")
print("-"*65)

best_cases = []

for a_gam in [0.05, 0.08, 0.10, 0.12, 0.15, 0.20]:
    for alpha in [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]:
        crossings, _, _ = find_crossings(a_gam, alpha)
        if len(crossings) >= 2:
            K1, K2 = crossings[0], crossings[1]
            m1 = Lambda * K1   # masa = Lambda * K (samospojnosc)
            m2 = Lambda * K2
            ratio = m2 / m1
            print(f"{a_gam:>6.3f} {alpha:>6.1f} {K1:>8.4f} {K2:>8.4f} "
                  f"{m1:>8.3f} {m2:>8.3f} {ratio:>8.1f}")
            best_cases.append((a_gam, alpha, K1, K2, m1, m2, ratio))

if not best_cases:
    print("Brak dwoch przeciec w tym zakresie.")

# Najlepsze dopasowanie do 207
if best_cases:
    print()
    print("Najblizsze m2/m1 = 207:")
    best_cases.sort(key=lambda x: abs(x[6] - 207))
    for case in best_cases[:5]:
        a_gam, alpha, K1, K2, m1, m2, ratio = case
        print(f"  a_Gam={a_gam:.3f} alpha={alpha:.1f}: "
              f"m2/m1={ratio:.1f}  (delta={abs(ratio-207):.1f})")

# Wykres dla najlepszego przypadku
if best_cases:
    best = min(best_cases, key=lambda x: abs(x[6] - 207))
    a_gam, alpha = best[0], best[1]
    print(f"\nWykres dla a_Gam={a_gam}, alpha={alpha}:")
    K_vals = np.linspace(0.02, 15.0, 300)
    EK_ratio = np.array([E_over_K(K, a_gam, alpha) for K in K_vals])
    fig, ax = plt.subplots(figsize=(9, 5))
    mask = np.isfinite(EK_ratio) & (EK_ratio > 0) & (EK_ratio < 200)
    ax.plot(K_vals[mask], EK_ratio[mask], 'b-', lw=2, label=r'$E(K)/K$')
    ax.axhline(Lambda, color='r', ls='--', lw=1.5, label=fr'$\Lambda={Lambda:.2f}$')
    crossings, _, _ = find_crossings(a_gam, alpha)
    for i, Kc in enumerate(crossings):
        ax.axvline(Kc, color='g', ls=':', lw=1)
        ax.plot(Kc, Lambda, 'go', ms=10,
                label=f'Gen {i+1}: K={Kc:.3f}, m={Lambda*Kc:.2f}')
    ax.set_xlabel(r'$K$ (sila zrodla)', fontsize=12)
    ax.set_ylabel(r'$E(K)/K$', fontsize=12)
    ax.set_title(f'TGP samospojnosc: a_Gam={a_gam}, alpha={alpha}', fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, min(200, max(EK_ratio[mask])*1.1))
    plt.tight_layout()
    plt.savefig('TGP/TGP_v1/scripts/advanced/dwa_przeciecia.png',
                dpi=150, bbox_inches='tight')
    print("Wykres: dwa_przeciecia.png")
