"""
TGP — Trzecia generacja przez czlon V5 = +delta*psi^5
======================================================

Potencjal rozszerzony:
    V(psi) = gamma/3*psi^3 - gamma/4*psi^4 + delta/5*psi^5

Czlon delta*psi^5 moze stworzyc lokalne minimum V dla psi > 1,
co daje drugie maksimum E(K)/K -> trzecie przeciecie samospojnosci.

Warunek na lokalne minimum V'(psi)=0 dla psi>1:
    V'(psi) = gamma*psi^2 - gamma*psi^3 + delta*psi^4 = 0
    => psi^2(gamma - gamma*psi + delta*psi^2) = 0
    Dla psi>0: gamma - gamma*psi + delta*psi^2 = 0
    => psi_min = [gamma +/- sqrt(gamma^2 - 4*delta*gamma)] / (2*delta)
    Warunek: delta < gamma/4

Szukamy: 3 przeciecia = 3 generacje z ratios m2/m1=207, m3/m1=3477
"""
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHI0 = 1.0; GAMMA = 1.0; Q = 1.0; M_SP = 1.0
Lambda = 4.0 * np.pi / Q   # stala samospojnosci

# ======================================================================
# POTENCJAL Z CZLONEM V5
# ======================================================================
def V_ext(phi, delta):
    """V(phi) - V(phi0), z czlonem kwintycznym."""
    psi  = phi / PHI0
    psi0 = 1.0
    V    = GAMMA/3*psi**3 - GAMMA/4*psi**4 + delta/5*psi**5
    V0   = GAMMA/3 - GAMMA/4 + delta/5
    return V - V0

def dV_dpsi(psi, delta):
    """Pochodna V wzgledem psi — do analizy struktury potencjalu."""
    return GAMMA*psi**2 - GAMMA*psi**3 + delta*psi**4

def find_V_extrema(delta):
    """Znajdz krytyczne punkty V dla psi > 0."""
    psi_vals = np.linspace(0.01, 10.0, 5000)
    dV = dV_dpsi(psi_vals, delta)
    extrema = []
    for i in range(len(dV)-1):
        if dV[i]*dV[i+1] < 0:
            psi_c = brentq(lambda p: dV_dpsi(p, delta), psi_vals[i], psi_vals[i+1])
            d2V   = GAMMA*2*psi_c - 3*GAMMA*psi_c**2 + 4*delta*psi_c**3
            typ   = 'min' if d2V > 0 else 'max'
            extrema.append((psi_c, typ))
    return extrema

# ======================================================================
# E(K)/K z potencjalem V5
# ======================================================================
def E_over_K(K, a_gam, alpha, delta):
    r    = np.linspace(a_gam, 25.0, 1500)
    phi  = PHI0 + K * np.exp(-M_SP*r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-M_SP*r) * (-M_SP*r - 1.0) / r**2
    E_k  = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1 + alpha/(PHI0*phi)) * r**2, r)
    V    = V_ext(phi, delta)
    E_p  = 4*np.pi * np.trapezoid(V * r**2, r)
    return (E_k + E_p) / K

def find_crossings(a_gam, alpha, delta, K_min=0.005, K_max=20.0, N=300):
    K_vals = np.linspace(K_min, K_max, N)
    f = np.array([E_over_K(K, a_gam, alpha, delta) - Lambda for K in K_vals])
    crossings = []
    for i in range(len(f)-1):
        if np.isfinite(f[i]) and np.isfinite(f[i+1]) and f[i]*f[i+1] < 0:
            try:
                Kc = brentq(
                    lambda K: E_over_K(K, a_gam, alpha, delta) - Lambda,
                    K_vals[i], K_vals[i+1], xtol=1e-8)
                crossings.append(Kc)
            except Exception:
                pass
    return crossings, K_vals, f

# ======================================================================
# ANALIZA STRUKTURY POTENCJALU
# ======================================================================
print("="*60)
print("Struktura potencjalu V(psi) = phi^3/3 - phi^4/4 + delta*phi^5/5")
print("="*60)
print(f"{'delta':>8} {'psi_min':>10} {'V(psi_min)-V(1)':>18} {'warunek':>12}")
print("-"*60)
for delta in [0.0, 0.05, 0.10, 0.15, 0.20, 0.22, 0.24]:
    extrema = find_V_extrema(delta)
    min_pts = [(p,t) for p,t in extrema if t=='min' and p > 1.01]
    if min_pts:
        psi_m = min_pts[0][0]
        dV_val = V_ext(PHI0*psi_m, delta) - V_ext(PHI0, delta)
        print(f"{delta:>8.3f} {psi_m:>10.3f} {dV_val:>18.5f}   drugie minimum")
    else:
        print(f"{delta:>8.3f} {'brak':>10} {'':>18}   brak minimum")

# ======================================================================
# SKAN TROJEK PRZECIEC
# ======================================================================
print()
print("="*70)
print("Skan (a_Gam, alpha, delta) -> liczba przeciec i stosunki mas")
print("="*70)
print(f"{'a_Gam':>6} {'alpha':>6} {'delta':>7} {'n':>4} "
      f"{'m1':>8} {'m2':>8} {'m3':>8} {'m2/m1':>8} {'m3/m1':>8}")
print("-"*70)

three_crossing_cases = []

for a_gam in [0.03, 0.05, 0.08, 0.10]:
    for alpha in [2.0, 5.0, 10.0, 20.0]:
        for delta in [0.05, 0.10, 0.15, 0.18, 0.20, 0.22]:
            if delta >= GAMMA/4:   # warunek konieczny: delta < gamma/4 = 0.25
                continue
            crossings, _, _ = find_crossings(a_gam, alpha, delta)
            n = len(crossings)
            masses = [Lambda * K for K in crossings]
            m1 = masses[0] if n >= 1 else None
            m2 = masses[1] if n >= 2 else None
            m3 = masses[2] if n >= 3 else None
            r21 = m2/m1 if n >= 2 else None
            r31 = m3/m1 if n >= 3 else None

            m1s = f"{m1:.2f}" if m1 else "-"
            m2s = f"{m2:.2f}" if m2 else "-"
            m3s = f"{m3:.2f}" if m3 else "-"
            r21s = f"{r21:.1f}" if r21 else "-"
            r31s = f"{r31:.1f}" if r31 else "-"

            if n >= 3:
                marker = "  <<<  3 generacje!"
                three_crossing_cases.append(
                    (a_gam, alpha, delta, masses, r21, r31))
            else:
                marker = ""

            if n >= 2:
                print(f"{a_gam:>6.3f} {alpha:>6.1f} {delta:>7.3f} {n:>4d} "
                      f"{m1s:>8} {m2s:>8} {m3s:>8} {r21s:>8} {r31s:>8}{marker}")

# ======================================================================
# NAJLEPSZE DOPASOWANIE DO 1:207:3477
# ======================================================================
print()
if three_crossing_cases:
    print("="*60)
    print("Przypadki z 3 generacjami — bliskosc do 1:207:3477")
    print("="*60)
    target21, target31 = 207.0, 3477.0
    scored = []
    for case in three_crossing_cases:
        a_gam, alpha, delta, masses, r21, r31 = case
        score = abs(r21 - target21)/target21 + abs(r31 - target31)/target31
        scored.append((score, case))
    scored.sort(key=lambda x: x[0])
    for score, (a_gam, alpha, delta, masses, r21, r31) in scored[:5]:
        print(f"  a_Gam={a_gam:.3f} alpha={alpha:.1f} delta={delta:.3f}: "
              f"m2/m1={r21:.1f} (cel 207), m3/m1={r31:.1f} (cel 3477), "
              f"score={score:.3f}")

    # Wykres dla najlepszego przypadku
    best = scored[0][1]
    a_gam, alpha, delta, masses, r21, r31 = best
    print(f"\nWykres dla najlepszego przypadku:")
    print(f"  a_Gam={a_gam}, alpha={alpha}, delta={delta}")
    K_vals = np.linspace(0.001, 20.0, 500)
    EK_r   = np.array([E_over_K(K, a_gam, alpha, delta) for K in K_vals])
    crossings, _, _ = find_crossings(a_gam, alpha, delta)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Lewy: E(K)/K
    ax = axes[0]
    mask = np.isfinite(EK_r) & (EK_r > 0) & (EK_r < 3*Lambda)
    ax.plot(K_vals[mask], EK_r[mask], 'b-', lw=2, label=r'$E(K)/K$')
    ax.axhline(Lambda, color='r', ls='--', lw=1.5,
               label=fr'$\Lambda={Lambda:.2f}$')
    colors = ['green', 'orange', 'purple']
    gen_names = ['Gen 1 (elektron)', 'Gen 2 (mion)', 'Gen 3 (tau)']
    for i, Kc in enumerate(crossings[:3]):
        ax.axvline(Kc, color=colors[i], ls=':', lw=1.5)
        ax.plot(Kc, Lambda, 'o', color=colors[i], ms=11,
                label=f'{gen_names[i]}: K={Kc:.3f}')
    ax.set_xlabel(r'$K$', fontsize=12)
    ax.set_ylabel(r'$E(K)/K$', fontsize=12)
    ax.set_title(f'Samospojnosc: 3 generacje\n'
                 f'a_Gam={a_gam}, alpha={alpha}, delta={delta}', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Prawy: potencjal V(psi)
    ax2 = axes[1]
    psi_v = np.linspace(0.1, 6.0, 500)
    V0   = V_ext(PHI0*psi_v, 0.0)
    V_d  = V_ext(PHI0*psi_v, delta)
    ax2.plot(psi_v, V0, 'b--', lw=1.5, label=r'$V(\psi)$, $\delta=0$')
    ax2.plot(psi_v, V_d, 'r-',  lw=2,   label=fr'$V(\psi)$, $\delta={delta}$')
    ax2.axhline(0, color='k', lw=0.5)
    ax2.axvline(1, color='gray', lw=0.5, ls=':')
    ax2.set_xlabel(r'$\psi = \Phi/\Phi_0$', fontsize=12)
    ax2.set_ylabel(r'$V(\psi) - V(1)$', fontsize=12)
    ax2.set_title('Potencjal TGP z czlonem V5', fontsize=10)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.5, 0.3)

    plt.tight_layout()
    plt.savefig('TGP/TGP_v1/scripts/advanced/trzecia_generacja.png',
                dpi=150, bbox_inches='tight')
    print("Wykres: trzecia_generacja.png")

else:
    print("Brak przypadkow z 3 przecieciami w przebadanym zakresie.")
    print("Sprawdz wiekszy zakres alpha lub mniejszy a_Gam.")
