"""
TGP — Wersja 2: Zmodyfikowany potencjal z czlonem stabilizujacym
=================================================================

Problem (wersja 1): V(psi) = psi^3/3 - psi^4/4
  Dla duzego K (duze psi w rdzeniu): -psi^4/4 -> -inf
  => g(M) = E(K) - M ma tylko 2 zera (elektron + mion)
  => brak 3. generacji (taon)

Rozwiazanie (wersja 2): dodanie czlonu stabilizujacego
  V_mod(psi) = V(psi) + lambda * (psi-1)^6 / 6

Wlasnosci modyfikacji:
  - (psi-1)^6|_{psi=1} = 0          => V_mod(1) = V(1) = 1/12 (próznia bez zmian)
  - d/dpsi[(psi-1)^6]|_{psi=1} = 0  => V_mod'(1) = 0 (warunek prozni zachowany)
  - (psi-1)^6 ~ psi^6 dla duzego psi => dominuje nad -psi^4/4 dla lambda > 0
  - Zapobiega katastrofie energetycznej przy duzym K

Crossover (psi duze): lambda*psi^6/6 > psi^4/4 gdy psi^2 > 3/(2*lambda)

Dla typowego psi_core ~ K/Phi_bg/a_gam ~ 2/(0.8*0.03) ~ 83:
  Potrzeba lambda > 3/(2*83^2) ~ 2e-4
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

GAMMA = 1.0


# ======================================================================
# ENERGIA Z ZMODYFIKOWANYM POTENCJALEM
# ======================================================================

def energy_v2(K, Phi_bg, msp, alpha, a_gam, lam):
    """
    E[Phi_bg + dPhi] - E[Phi_bg] z V_mod(psi) = V_orig(psi) + lam*(psi-1)^6/6.

    V_mod(psi) - V_mod(1) = (V_orig(psi) - V_orig(1)) + lam*(psi-1)^6/6
                                                         ^^^^^^^^^^^^^^^^
                                                         czlon dodatkowy (zero przy psi=1)
    """
    if msp < 1e-8 or Phi_bg <= 0:
        return 0.0
    r_max = max(40.0 / msp, 15.0)
    r = np.linspace(a_gam, r_max, 3000)

    phi  = Phi_bg + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    # Ekin: bez zmian (nie zalezny od V)
    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / (Phi_bg * phi)) * r**2, r)

    # Epot: V_orig + czlon stabilizujacy
    psi  = phi / Phi_bg
    V1   = GAMMA/3 - GAMMA/4          # V_orig(psi=1) = 1/12

    Vpsi_orig = GAMMA/3 * psi**3 - GAMMA/4 * psi**4
    Vstab     = lam * (psi - 1.0)**6 / 6.0   # zero przy psi=1

    Ep = 4*np.pi * Phi_bg**2 * np.trapezoid(
        (Vpsi_orig - V1 + Vstab) * r**2, r)

    return Ek + Ep


def g_v2(M, xi, alpha, a_gam, lam, Phi0_total=1.0):
    """
    g(M) = E_v2(K;Phi_bg) - M = 0

    Phi_bg = Phi0_total - xi*M   [MINUS, poprawna samospojnosc]
    K      = M / (4*pi*Phi_bg)
    msp    = sqrt(gamma/Phi0_total)
    """
    if M < 1e-12:
        return 0.0
    Phi_bg = Phi0_total - xi * M
    if Phi_bg <= 1e-6:
        return float('inf')
    msp = np.sqrt(max(GAMMA / Phi0_total, 1e-8))
    K   = M / (4.0 * np.pi * Phi_bg)
    E   = energy_v2(K, Phi_bg, msp, alpha, a_gam, lam)
    return E - M


def find_zeros(xi, alpha, a_gam, lam, Phi0_total=1.0, N=600):
    """Szuka wszystkich zer g(M) w przedziale (0, 0.95*Phi0/xi)."""
    M_max = 0.95 * Phi0_total / xi if xi > 0 else 1e4
    # Gestsze próbkowanie przy malym M (korzenie M1 sa tu)
    M_lo  = np.linspace(1e-5,   M_max*0.03, int(N*0.25))
    M_mid = np.linspace(M_max*0.03, M_max*0.3, int(N*0.35))
    M_hi  = np.linspace(M_max*0.3,  M_max,    int(N*0.4))
    M_arr = np.unique(np.concatenate([M_lo, M_mid, M_hi]))

    g_arr = np.array([g_v2(M, xi, alpha, a_gam, lam, Phi0_total)
                      for M in M_arr])

    roots = []
    for i in range(len(g_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)):
            continue
        if gi * gj < 0:
            try:
                M_root = brentq(
                    lambda M: g_v2(M, xi, alpha, a_gam, lam, Phi0_total),
                    M_arr[i], M_arr[i+1], xtol=1e-9, maxiter=100)
                roots.append(M_root)
            except Exception:
                pass
    return M_arr, g_arr, sorted(roots)


# ======================================================================
# 1. DIAGNOSTYKA: jak lambda zmienia ksztalt g(M)?
# ======================================================================

print("=" * 75)
print("TGP v2: Potencjal zmodyfikowany V_mod = V + lambda*(psi-1)^6/6")
print("  Cel: zapobiec katastrofie psi^4 => uzyskac 3. zero g(M)")
print("=" * 75)

xi, alpha, a_gam, Phi0T = 0.001, 7.0, 0.030, 1.0

print(f"\n--- Liczba zer g(M) vs lambda  (xi={xi}, alpha={alpha}, a_gam={a_gam}, Phi0T={Phi0T}) ---")
print(f"{'lambda':>10} {'n_zer':>7} {'M1':>10} {'M2/M1':>9} {'M3/M1':>9}")
print("-" * 50)

lam_scan = [0.0, 1e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3,
            5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
scan_results = []

for lam in lam_scan:
    _, _, roots = find_zeros(xi, alpha, a_gam, lam, Phi0T)
    n = len(roots)
    M1  = roots[0] if n >= 1 else float('nan')
    r21 = roots[1]/roots[0] if n >= 2 else float('nan')
    r31 = roots[2]/roots[0] if n >= 3 else float('nan')
    s21 = f"{r21:>9.1f}" if np.isfinite(r21) else f"{'—':>9}"
    s31 = f"{r31:>9.1f}" if np.isfinite(r31) else f"{'—':>9}"
    flag = " *** 3 GEN!" if n >= 3 else ""
    print(f"{lam:>10.5f} {n:>7d} {M1:>10.5f} {s21} {s31}{flag}")
    scan_results.append((lam, n, roots))


# ======================================================================
# 2. SKAN: (lambda, alpha, a_gam) -> 3 generacje?
# ======================================================================

print(f"\n{'='*75}")
print("SKAN PARAMETROW: kiedy pojawia sie 3. generacja?")
print(f"{'lambda':>8} {'alpha':>6} {'a_gam':>6} {'n':>4} "
      f"{'M1':>9} {'r21':>8} {'r31':>9} {'bl_r21':>8} {'bl_r31':>8}")
print("-" * 75)

best_3gen = []
cases_scan = []
for lam in [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05]:
    for alpha in [5.0, 5.9, 7.0, 8.0, 10.0]:
        for a_gam in [0.02, 0.03, 0.05, 0.08]:
            cases_scan.append((lam, alpha, a_gam))

xi_fixed = 0.001
Phi0T_fixed = 1.0

for lam, alpha, a_gam in cases_scan:
    _, _, roots = find_zeros(xi_fixed, alpha, a_gam, lam, Phi0T_fixed)
    n = len(roots)
    if n < 2:
        continue
    M1  = roots[0]
    r21 = roots[1]/M1
    r31 = roots[2]/M1 if n >= 3 else float('nan')

    bl21 = 100*(r21 - 207)/207
    bl31 = 100*(r31 - 3477)/3477 if np.isfinite(r31) else float('nan')

    s31   = f"{r31:>9.1f}" if np.isfinite(r31) else f"{'—':>9}"
    sbl31 = f"{bl31:>+8.1f}%" if np.isfinite(bl31) else f"{'—':>8}"
    flag  = " *** 3G!" if n >= 3 else ""

    if n >= 3:
        print(f"{lam:>8.5f} {alpha:>6.1f} {a_gam:>6.3f} {n:>4d} "
              f"{M1:>9.5f} {r21:>8.1f} {s31} {bl21:>+8.1f}% {sbl31}{flag}")
        best_3gen.append(dict(lam=lam, alpha=alpha, a_gam=a_gam,
                              xi=xi_fixed, Phi0T=Phi0T_fixed,
                              roots=roots, M1=M1, r21=r21, r31=r31))
    elif abs(r21 - 207) < 30:
        # Ciekawe 2-generacyjne (blisko 207)
        print(f"{lam:>8.5f} {alpha:>6.1f} {a_gam:>6.3f} {n:>4d} "
              f"{M1:>9.5f} {r21:>8.1f} {s31} {bl21:>+8.1f}% {sbl31}")


# ======================================================================
# 3. NAJLEPSZE 3-GENERACYJNE
# ======================================================================

print(f"\n{'='*75}")
if best_3gen:
    print(f"ZNALEZIONO {len(best_3gen)} przypadkow 3-generacyjnych!")
    print("\nNajlepsze trafienia w 1:207:3477:")
    def score(r):
        d21 = abs(r['r21']-207)/207 if np.isfinite(r['r21']) else 1.0
        d31 = abs(r['r31']-3477)/3477 if np.isfinite(r['r31']) else 1.0
        return d21 + d31
    best_3gen.sort(key=score)
    print(f"{'lambda':>8} {'alpha':>6} {'a_gam':>6} {'r21':>8} {'r31':>9} {'score':>8}")
    print("-" * 55)
    for r in best_3gen[:12]:
        print(f"{r['lam']:>8.5f} {r['alpha']:>6.1f} {r['a_gam']:>6.3f} "
              f"{r['r21']:>8.1f} {r['r31']:>9.1f} {score(r):>8.4f}")

    # Szczegoly najlepszego
    best = best_3gen[0]
    print(f"\nNajlepszy: lambda={best['lam']:.5f}, alpha={best['alpha']:.1f}, "
          f"a_gam={best['a_gam']:.3f}")
    print(f"  M1 = {best['M1']:.6f}")
    for i, M in enumerate(best['roots'][:3]):
        print(f"  M{i+1} = {M:.6f}  (M{i+1}/M1 = {M/best['M1']:.1f})")
    print(f"  r21 = {best['r21']:.2f}  (cel: 207, blad: {100*(best['r21']-207)/207:+.1f}%)")
    print(f"  r31 = {best['r31']:.2f}  (cel: 3477, blad: {100*(best['r31']-3477)/3477:+.1f}%)")

else:
    print("Brak 3 generacji w przeskanym zakresie parametrow.")
    print("Sprawdz diagnostyke lambda powyzej.")


# ======================================================================
# 4. WYKRES: g(M) dla roznych lambda
# ======================================================================

fig, axes = plt.subplots(2, 4, figsize=(18, 9))
axes = axes.flatten()

plot_lambdas = [0.0, 1e-4, 2e-4, 5e-4, 1e-3, 5e-3, 0.01, 0.05]
xi_p, alp_p, ag_p, Phi0T_p = 0.001, 7.0, 0.030, 1.0

for idx, lam in enumerate(plot_lambdas):
    M_arr, g_arr, roots = find_zeros(xi_p, alp_p, ag_p, lam, Phi0T_p)
    ax = axes[idx]

    finite = np.isfinite(g_arr)
    clip   = np.abs(g_arr) < 200
    mask   = finite & clip

    ax.plot(M_arr[mask], g_arr[mask], 'b-', lw=2)
    ax.axhline(0, color='k', lw=1, ls='--')

    colors = ['green', 'orange', 'red']
    for i, Mc in enumerate(roots[:3]):
        ax.axvline(Mc, color=colors[i], lw=1.2, ls=':')
        ax.plot(Mc, 0, 'o', color=colors[i], ms=9,
                label=f'M{i+1}={Mc:.3f}')

    n = len(roots)
    extra = ""
    if n >= 2:
        extra = f"\nr21={roots[1]/roots[0]:.0f}"
    if n >= 3:
        extra += f"  r31={roots[2]/roots[0]:.0f}"

    ax.set_title(f'lambda={lam:.0e}  n={n} gen{extra}', fontsize=9)
    ax.set_xlabel('M', fontsize=9)
    ax.set_ylabel('g(M)', fontsize=9)
    if roots:
        ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    yr = 50 if not mask.sum() else min(200, np.nanmax(np.abs(g_arr[mask])) * 0.3 + 5)
    ax.set_ylim(-yr, yr)
    ax.set_xlim(0, M_arr[-1])

plt.suptitle(
    f'TGP v2: g(M) dla V_mod = V + lambda*(psi-1)^6/6\n'
    f'xi={xi_p}, alpha={alp_p}, a_gam={ag_p}, Phi0T={Phi0T_p}',
    fontsize=11)
plt.tight_layout()
out_plot = 'TGP/TGP_v1/scripts/advanced/v2_lambda.png'
plt.savefig(out_plot, dpi=150, bbox_inches='tight')
print(f"\nWykres: {out_plot}")


# ======================================================================
# 5. JESLI SA 3 GENERACJE: dokladny skan wokol najlepszego
# ======================================================================

if best_3gen:
    best = best_3gen[0]
    print(f"\n{'='*75}")
    print(f"DOKLADNY SKAN wokol najlepszego: "
          f"lambda={best['lam']:.5f}, alpha={best['alpha']:.1f}, "
          f"a_gam={best['a_gam']:.3f}")
    print(f"{'lambda':>10} {'alpha':>7} {'a_gam':>7} {'r21':>8} {'r31':>9} "
          f"{'bl21':>7} {'bl31':>7}")
    print("-" * 60)

    for lam2 in np.linspace(best['lam']*0.3, best['lam']*3.0, 12):
        for alpha2 in np.linspace(max(3., best['alpha']-2), best['alpha']+2, 6):
            for ag2 in [best['a_gam']*0.7, best['a_gam'], best['a_gam']*1.4]:
                _, _, r2 = find_zeros(best['xi'], alpha2, ag2, lam2, best['Phi0T'])
                if len(r2) >= 3:
                    r21_ = r2[1]/r2[0]
                    r31_ = r2[2]/r2[0]
                    bl21 = 100*(r21_-207)/207
                    bl31 = 100*(r31_-3477)/3477
                    print(f"{lam2:>10.5f} {alpha2:>7.2f} {ag2:>7.3f} "
                          f"{r21_:>8.1f} {r31_:>9.1f} "
                          f"{bl21:>+7.1f}% {bl31:>+7.1f}%")

print("\nGotowe.")
