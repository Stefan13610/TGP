"""
TGP — Sprzezony uklad samospojnosci (M, Phi_bg)
================================================

POPRAWNA SAMOSPOJNOSC (per rysunek "Wlasciwa samospojnosc"):

  Phi_bg = Phi0_total - qM / (msp^2 * V)   [tlo BEZ tej czastki]
  K      = qM / (4*pi*Phi_bg)               [sprzezenie do pola]
  M      = E[Phi_bg + dPhi] - E[Phi_bg]     [masa = roznica energii]

Gdzie dPhi(r) = K*exp(-msp*r)/r  (PLUS — czastka to lokalne wybrzuszenie nad Phi_bg).

Phi0_total = calkowite tlo kosmologiczne (STALE dla danej epoki).
xi         = q/(msp^2 * V)  — efektywne sprzezenie (skupia msp, V, q).
Phi_bg(M)  = Phi0_total - xi*M  (MALEJE gdy M rosnie — odejmujemy wklad czastki).

Uklad sprzezony:
  { Phi_bg = Phi0_total - xi*M
  { K      = M / (4*pi*Phi_bg)
  { msp    = sqrt(gamma / Phi0_total)    [masa skalaru z globalnego tla]
  { M      = E(K; Phi_bg, msp) / c0^2

Szukamy: g(M) = E(K(M); Phi_bg(M)) - M = 0

Wiele rozwiazan g(M)=0 => wiele generacji!

Blad oryginalnego skryptu: uzywal Phi0 = PHI0_BASE + xi*M (PLUS zamiast MINUS),
czyli wlaczal wklad czastki do tla zamiast go odejmowac.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

GAMMA     = 1.0
LAM = 5.514e-7    # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)
PHI0_BASE = 1.0   # normalizacja prozni (psi=Phi/Phi0=1 to punkty stacjonarne)


# ======================================================================
# ENERGIA TGP: E[Phi_bg + dPhi] - E[Phi_bg]
# ======================================================================
def energy(K, Phi_bg, msp, alpha, a_gam):
    """
    Energia czastki w tle Phi_bg.
    Profil: Phi(r) = Phi_bg + K*exp(-msp*r)/r   [PLUS — wybrzuszenie]
    Energia = calka (Ekin + Epot) r^2 dr — energia prozni (E[Phi_bg]) odjeta
              implicite przez (V(psi) - V(1)).

    Args:
        K      : amplituda Yukawa = M/(4*pi*Phi_bg)
        Phi_bg : tlo BEZ tej czastki
        msp    : masa skalaru (z globalnego Phi0_total)
        alpha  : wspolczynnik TGP nieliniowy kinetyczny
        a_gam  : krok sieci substratu
    """
    if msp < 1e-8 or Phi_bg <= 0:
        return 0.0
    r_max = max(40.0 / msp, 15.0)
    r = np.linspace(a_gam, r_max, 2500)

    phi  = Phi_bg + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    # Ekin: standardowy + TGP nieliniowy (kinetykcznie generowany uklad)
    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / (Phi_bg * phi)) * r**2, r)

    # Epot: V(phi/Phi_bg) - V(1)  [odejmuje energie prozni]
    psi  = phi / Phi_bg
    Vpsi = GAMMA/3 * psi**3 - GAMMA/4 * psi**4
    V1   = GAMMA/3 - GAMMA/4          # V(psi=1)
    Vstab = LAM / 6.0 * (psi - 1.0)**6   # V_mod stabilizacja
    Ep   = 4*np.pi * Phi_bg**2 * np.trapezoid((Vpsi - V1 + Vstab) * r**2, r)

    return Ek + Ep


# ======================================================================
# POPRAWIONY UKLAD SPRZEZONY
# ======================================================================
def g_coupled(M, xi, alpha, a_gam, Phi0_total=1.0):
    """
    g(M) = E(K(M); Phi_bg(M)) - M  dla poprawnej samospojnosci.

    Phi_bg = Phi0_total - xi*M     [MINUS: odejmujemy wklad czastki]
    K      = M / (4*pi*Phi_bg)
    msp    = sqrt(gamma/Phi0_total)  [z globalnego tla]
    """
    if M < 1e-12:
        return 0.0

    Phi_bg = Phi0_total - xi * M
    if Phi_bg <= 1e-6:
        return float('inf')        # Phi_bg ujemne — poza zakresem

    msp = np.sqrt(max(GAMMA / Phi0_total, 1e-8))
    K   = M / (4.0 * np.pi * Phi_bg)

    E = energy(K, Phi_bg, msp, alpha, a_gam)
    return E - M


# ======================================================================
# SZUKANIE PUNKTOW STALYCH (generacji)
# ======================================================================
def find_fixed_points(xi, alpha, a_gam, Phi0_total=1.0,
                      M_min=1e-5, M_max=None, N=500):
    """
    Skanuje g(M) i zwraca wszystkie zera.
    M_max ograniczone przez Phi_bg > 0: M < Phi0_total/xi.
    """
    if M_max is None:
        M_max = 0.95 * Phi0_total / xi if xi > 0 else 1e4

    M_lo  = np.linspace(M_min, M_max*0.05,  int(N*0.3))
    M_mid = np.linspace(M_max*0.05, M_max*0.4, int(N*0.4))
    M_hi  = np.linspace(M_max*0.4, M_max,   int(N*0.3))
    M_arr = np.unique(np.concatenate([M_lo, M_mid, M_hi]))

    g_arr = np.array([g_coupled(M, xi, alpha, a_gam, Phi0_total) for M in M_arr])

    roots = []
    for i in range(len(g_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)):
            continue
        if gi * gj < 0:
            try:
                M_root = brentq(
                    lambda M: g_coupled(M, xi, alpha, a_gam, Phi0_total),
                    M_arr[i], M_arr[i+1], xtol=1e-9, maxiter=80)
                roots.append(M_root)
            except Exception:
                pass

    return M_arr, g_arr, roots


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':

    print("=" * 75)
    print("TGP — Poprawiony uklad sprzezony (znak MINUS — wlasciwa samospojnosc)")
    print("  Phi_bg = Phi0_total - xi*M   [BEZ wkladu tej czastki]")
    print("  K      = M / (4*pi*Phi_bg)")
    print("  msp    = sqrt(gamma/Phi0_total)  [z globalnego tla]")
    print("  M      = E(K; Phi_bg, msp)")
    print("=" * 75)

    # ------------------------------------------------------------------
    # 1. DIAGNOSTYKA: jak wyglada g(M) dla przykladowych parametrow?
    # ------------------------------------------------------------------
    print("\n--- Diagnostyka g(M) ---")
    print(f"{'M':>8} {'Phi_bg':>8} {'K':>10} {'E':>10} {'g=E-M':>10}")
    xi_d, alpha_d, a_gam_d, Phi0T = 0.01, 5.9, 0.05, 1.0
    for M_d in [0.01, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0, 10.0,
                20.0, 40.0, 60.0, 80.0, 95.0]:
        Phi_bg = Phi0T - xi_d * M_d
        if Phi_bg <= 0:
            print(f"{M_d:>8.2f} {'ujemne':>8}")
            break
        msp = np.sqrt(GAMMA / Phi0T)
        K   = M_d / (4*np.pi*Phi_bg)
        E   = energy(K, Phi_bg, msp, alpha_d, a_gam_d)
        print(f"{M_d:>8.2f} {Phi_bg:>8.4f} {K:>10.5f} {E:>10.4f} {E-M_d:>+10.4f}")

    # ------------------------------------------------------------------
    # 2. SKAN: xi, alpha, a_gam, Phi0_total
    # ------------------------------------------------------------------
    print(f"\n{'='*75}")
    print("SKAN: (xi, alpha, a_gam, Phi0_total) => generacje")
    print(f"{'xi':>7} {'alpha':>6} {'a_gam':>6} {'Phi0T':>6} {'n':>4} "
          f"{'M1':>8} {'M2/M1':>8} {'M3/M1':>9}")
    print("-" * 65)

    all_results = []
    best_3gen   = []

    scan_configs = [
        # (xi,   [alphas],             [a_gams],           [Phi0_totals])
        (0.001, [5.0, 5.9, 7.0],      [0.03, 0.05],       [1.0, 2.0, 5.0]),
        (0.005, [5.0, 5.9, 7.0],      [0.03, 0.05],       [1.0, 2.0, 5.0]),
        (0.010, [4.0, 5.9, 7.0, 10.], [0.03, 0.05, 0.08], [1.0, 2.0, 5.0, 10.]),
        (0.050, [5.9, 7.0, 10.],      [0.03, 0.05],       [2.0, 5.0, 10., 20.]),
        (0.100, [5.9, 7.0, 10.],      [0.03, 0.05],       [5.0, 10., 20., 50.]),
        (0.200, [5.9, 7.0],           [0.03, 0.05],       [10., 20., 50.]),
        (0.500, [5.9, 7.0],           [0.03, 0.05],       [20., 50., 100.]),
        (1.000, [5.9, 7.0],           [0.03, 0.05],       [50., 100., 200.]),
    ]

    for xi, alphas, a_gams, Phi0Ts in scan_configs:
        for alpha in alphas:
            for a_gam in a_gams:
                for Phi0T in Phi0Ts:
                    _, _, roots = find_fixed_points(
                        xi, alpha, a_gam, Phi0_total=Phi0T, N=400)
                    n = len(roots)
                    if n < 2:
                        continue

                    M1  = roots[0]
                    r21 = roots[1]/M1 if n >= 2 else None
                    r31 = roots[2]/M1 if n >= 3 else None

                    flag = "  <-- 3 GEN!" if n >= 3 else ""
                    print(f"{xi:>7.3f} {alpha:>6.1f} {a_gam:>6.3f} {Phi0T:>6.1f} "
                          f"{n:>4d} {M1:>8.4f} "
                          f"{r21:>8.1f} "
                          f"{r31:>9.1f}{flag}" if r31 else
                          f"{xi:>7.3f} {alpha:>6.1f} {a_gam:>6.3f} {Phi0T:>6.1f} "
                          f"{n:>4d} {M1:>8.4f} "
                          f"{r21:>8.1f} "
                          f"{'—':>9}{flag}")

                    res = dict(xi=xi, alpha=alpha, a_gam=a_gam,
                               Phi0T=Phi0T, n=n, roots=roots,
                               M1=M1, r21=r21, r31=r31)
                    all_results.append(res)
                    if n >= 3:
                        best_3gen.append(res)

    # ------------------------------------------------------------------
    # 3. NAJLEPSZE TRAFIENIE w 1:207:3477
    # ------------------------------------------------------------------
    print(f"\n{'='*75}")
    if best_3gen:
        print("PRZYPADKI Z 3 GENERACJAMI — odleglosc od 1:207:3477:")
        scored = []
        for r in best_3gen:
            d21 = abs(r['r21'] - 207)/207 if r['r21'] else 1.0
            d31 = abs(r['r31'] - 3477)/3477 if r['r31'] else 1.0
            scored.append((d21 + d31, r))
        scored.sort(key=lambda x: x[0])
        for score, r in scored[:8]:
            print(f"  xi={r['xi']:.3f} alpha={r['alpha']:.1f} "
                  f"a_gam={r['a_gam']:.3f} Phi0T={r['Phi0T']:.1f}: "
                  f"M2/M1={r['r21']:.1f} (cel 207), "
                  f"M3/M1={r['r31']:.1f} (cel 3477), "
                  f"score={score:.4f}")
    else:
        print("Brak 3 generacji w przeskanym zakresie.")
        two_gen = [r for r in all_results if r['n'] >= 2]
        if two_gen:
            best2 = max(two_gen, key=lambda r: r['r21'] if r['r21'] else 0)
            print(f"  Najwyzsze r21: xi={best2['xi']}, alpha={best2['alpha']}, "
                  f"a_gam={best2['a_gam']}, Phi0T={best2['Phi0T']}: "
                  f"M2/M1={best2['r21']:.1f}")

    # ------------------------------------------------------------------
    # 4. WYKRES: krzywa g(M) dla wybranych xi i Phi0_total
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    plot_cases = [
        (0.010, 5.9, 0.05, 1.0),
        (0.010, 5.9, 0.05, 2.0),
        (0.010, 5.9, 0.05, 5.0),
        (0.050, 5.9, 0.05, 5.0),
        (0.050, 7.0, 0.05, 10.0),
        (0.100, 7.0, 0.03, 20.0),
    ]

    for idx, (xi_p, alp, ag, Phi0T) in enumerate(plot_cases):
        M_arr, g_arr, roots = find_fixed_points(
            xi_p, alp, ag, Phi0_total=Phi0T, N=400)

        ax = axes[idx]
        finite = np.isfinite(g_arr)
        clip   = np.abs(g_arr) < 500
        mask   = finite & clip
        ax.plot(M_arr[mask], g_arr[mask], 'b-', lw=2)
        ax.axhline(0, color='k', lw=1, ls='--')

        colors = ['green', 'orange', 'purple']
        for i, Mc in enumerate(roots[:3]):
            ax.axvline(Mc, color=colors[i], lw=1.2, ls=':')
            ax.plot(Mc, 0, 'o', color=colors[i], ms=9,
                    label=f'M{i+1}={Mc:.3f}')

        n = len(roots)
        title_extra = ""
        if n >= 2:
            title_extra = f"  r21={roots[1]/roots[0]:.0f}"
        if n >= 3:
            title_extra += f"  r31={roots[2]/roots[0]:.0f}"

        ax.set_title(f'xi={xi_p} a={alp} ag={ag} Phi0T={Phi0T}\n'
                     f'n={n} gen{title_extra}', fontsize=9)
        ax.set_xlabel('M', fontsize=9)
        ax.set_ylabel('g(M)=E-M', fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        yr = np.nanmax(np.abs(g_arr[mask])) * 0.25 if mask.sum() > 0 else 5
        ax.set_ylim(-yr, yr)
        ax.set_xlim(0, M_arr[-1])

    plt.suptitle('TGP: g(M) = E(K;Phi_bg) - M\n'
                 'POPRAWIONY: Phi_bg = Phi0_total - xi*M  (znak MINUS)', fontsize=11)
    plt.tight_layout()
    out = 'TGP/TGP_v1/scripts/advanced/sprzezony_uklad.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"\nWykres: {out}")
    print("Gotowe.")
