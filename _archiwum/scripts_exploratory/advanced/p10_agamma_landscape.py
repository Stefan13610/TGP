"""
p10_agamma_landscape.py
=======================
Mapa energetyczna g(K, a_Gamma) dla pelnego ODE TGP.

PYTANIE KLUCZOWE:
  Czy K_2=2.033 i K_3=34.14 moga byc samospojnymi solitonami
  przy odpowiednio dobranym a_Gamma^(n)?

PODEJSCIE:
  Dla kazdego (K, a_Gam): znajdz psi_core taki ze phi(R_max)=1,
  oblicz g(K, a_Gam) = E[phi]/(4*pi*K) - 1.

  Zera g = samospojne solitony TGP.

WYNIKI Z POPRZEDNICH SKRYPTOW:
  K*_1(ODE) = 0.010414, a_Gam=0.040 -> g=-1.9e-8 (potwierdzony)
  K_2=2.033, a_Gam=0.040  -> brak profilu (phi -> N_0)
  K_2=2.033, a_Gam=7.81   -> g=-4.53 (Schwarzschild scaling)

ANALIZA:
  1. Skan g(K_2, a_Gam) dla a_Gam in [0.1 .. 50]  -> czy istnieje a_Gam* ?
  2. Skan g(K_3, a_Gam) dla a_Gam in [0.5 .. 500]
  3. Mapa 2D g(K, a_Gam) na siatce 20x20
  4. Kontur g=0: linia samospojnych solitonow w przestrzeni (K, a_Gam)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# PARAMETRY
# ============================================================
ALPHA  = 8.5616
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)

V1 = GAMMA/3 - GAMMA/4   # V_mod(1) = 1/3 - 1/4 + 0

def V_mod(p):
    return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6

def dV_mod(p):
    return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


def phi_at_rmax(psi_core, K, a_gam, r_max, n_eval=2000):
    """Outward integration, return phi(r_max)."""
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1

    r_eval = a_gam * (r_max / a_gam) ** np.linspace(0, 1, n_eval)
    sol = solve_ivp(
        ode_rhs, [a_gam, r_max], [psi_core, dphi0],
        method='DOP853', rtol=1e-9, atol=1e-11,
        t_eval=r_eval, events=[ev_lo]
    )
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]


def best_g(K, a_gam, r_max=None, n_psi=120, psi_max=None):
    """
    Dla (K, a_gam): skanuj psi_core, znajdz galaz z F=0 najblizsza g=0.
    Zwraca (g_best, psi_core_best) lub (nan, None).
    """
    if r_max is None:
        # R_max: tak by Yukawa decay exp(-m_eff*R) < 1e-6
        r_max = max(40.0, a_gam * 5, 6.0 / M_EFF)

    if psi_max is None:
        psi_max = max(3.0, 1.0 + 2.0 * K / a_gam)

    psi_vals = np.linspace(1.001, min(psi_max, 300.0), n_psi)
    F_vals   = np.array([phi_at_rmax(p, K, a_gam, r_max) - 1.0
                         for p in psi_vals])

    # Zbierz wszystkie galezi (zmiany znaku)
    branches = []
    for i in range(len(F_vals) - 1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi * Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K, a_gam, r_max) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-5, rtol=1e-5, maxiter=40
                )
                # Oblicz energie
                dphi0 = -K / a_gam**2
                r_ev  = a_gam * (r_max / a_gam) ** np.linspace(0, 1, 4000)
                sol   = solve_ivp(ode_rhs, [a_gam, r_max], [psi_z, dphi0],
                                  method='DOP853', rtol=1e-9, atol=1e-11,
                                  t_eval=r_ev)
                r   = sol.t
                phi = np.maximum(sol.y[0], 1e-10)
                dphi = sol.y[1]
                Ek  = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
                Ep  = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
                g   = (Ek + Ep) / (4*np.pi*K) - 1.0
                branches.append((psi_z, g))
            except Exception:
                pass

    if not branches:
        return np.nan, None

    # Zwroc galaz z najmniejszym |g|
    best = min(branches, key=lambda x: abs(x[1]))
    return best[1], best[0]


# ============================================================
# 1. SKAN g(K2, a_Gam)  —  czy istnieje samospojne a_Gam*?
# ============================================================
K_list = {
    'K1*=0.01041': 0.010414,
    'K2=2.033':    2.032728,
    'K3=34.14':    34.14450,
}

print("=" * 70)
print("P10: MAPA g(K, a_Gam) — poszukiwanie samospojnych solitonow ODE")
print("=" * 70)
print()

fig_scans, axes_scans = plt.subplots(1, 3, figsize=(16, 5))
fig_scans.suptitle(
    f'p10: g(K, a_Gamma) — skan 1D  (alpha={ALPHA}, lam={LAM:.2e})',
    fontsize=11, fontweight='bold'
)

scan_results = {}

for ax, (label, K) in zip(axes_scans, K_list.items()):
    print(f"{'='*60}")
    print(f"K = {K:.5e}  ({label})")

    # Ustal zakres a_Gam
    if K < 0.1:
        a_range = np.logspace(-2, 0, 50)   # 0.01 .. 1
    elif K < 5:
        a_range = np.logspace(-2, 2, 60)   # 0.01 .. 100
    else:
        a_range = np.logspace(-1, 3, 70)   # 0.1 .. 1000

    g_arr  = []
    psi_arr = []
    valid_a = []

    for a in a_range:
        g, psi = best_g(K, a)
        g_arr.append(g)
        psi_arr.append(psi if psi else np.nan)
        valid_a.append(a)

    g_arr  = np.array(g_arr)
    valid_a = np.array(valid_a)
    scan_results[label] = (valid_a, g_arr)

    # Znajdz zera g(a_Gam)
    zero_a = []
    for i in range(len(g_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                def func_g(a):
                    val, _ = best_g(K, a)
                    return val if np.isfinite(val) else np.nan
                a_star = brentq(func_g, valid_a[i], valid_a[i+1],
                                xtol=1e-4, rtol=1e-4, maxiter=25)
                g_star, psi_star = best_g(K, a_star)
                zero_a.append((a_star, g_star, psi_star))
            except Exception as e:
                pass

    print(f"  a_Gam zakres: [{a_range[0]:.3f}, {a_range[-1]:.1f}]")
    print(f"  g_min={g_arr[np.isfinite(g_arr)].min():.3e}" if np.isfinite(g_arr).any() else "  (brak wynikow)")
    print(f"  g_max={g_arr[np.isfinite(g_arr)].max():.3e}" if np.isfinite(g_arr).any() else "")
    if zero_a:
        print(f"  Znaleziono {len(zero_a)} samospojnych a_Gam*:")
        for a_star, g_star, psi_star in zero_a:
            print(f"    a_Gam* = {a_star:.4f}  (g={g_star:.3e}, psi_core={psi_star:.4f})")
            print(f"    a_Gam*/K = {a_star/K:.4f}   (K1 referecja: {0.040/0.010414:.4f})")
    else:
        print("  BRAK zer g(K, a_Gam) w przebadanym zakresie!")
        if np.isfinite(g_arr).any():
            fin = g_arr[np.isfinite(g_arr)]
            if fin.min() > 0:
                print("  -> g > 0 wszedzie: E > 4*pi*K (soliton 'za ciezki')")
            elif fin.max() < 0:
                print("  -> g < 0 wszedzie: E < 4*pi*K (soliton 'za lekki')")
            else:
                print("  -> Zmiana znaku g bez bisekcji (nieciaglosc?)")
    print()

    # Wykres
    mask = np.isfinite(g_arr)
    g_clip = np.clip(g_arr, -10, 15)
    ax.plot(valid_a[mask], g_clip[mask], 'b-', lw=1.5)
    ax.axhline(0, color='red', lw=1.2, linestyle='--', label='g=0 (samospojny)')
    for a_star, g_star, _ in zero_a:
        ax.axvline(a_star, color='green', lw=1.5, alpha=0.8,
                   label=f'a*={a_star:.3f}')
    ax.set_xscale('log')
    ax.set_xlabel('a_Gamma')
    ax.set_ylabel('g(K, a_Gamma)')
    ax.set_title(f'{label}', fontsize=9)
    ax.set_ylim(-8, 10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
out1 = __file__.replace('.py', '_scan1D.png')
plt.savefig(out1, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out1}")
print()


# ============================================================
# 2. MAPA 2D g(K, a_Gam)
# ============================================================
print("=" * 70)
print("MAPA 2D g(K, a_Gam):")
print()

# Siatka — logarytmiczna w obu osiach
K_2d   = np.logspace(-3, 2, 22)     # 22 wartosci K
A_2d   = np.logspace(-2, 2, 22)     # 22 wartosci a_Gam

G_2d = np.full((len(A_2d), len(K_2d)), np.nan)

total = len(K_2d) * len(A_2d)
done  = 0
for j, K in enumerate(K_2d):
    for i, a in enumerate(A_2d):
        g, _ = best_g(K, a, n_psi=80)
        G_2d[i, j] = g
        done += 1
    print(f"  K={K:.3e}: g range = [{np.nanmin(G_2d[:,j]):.2e}, {np.nanmax(G_2d[:,j]):.2e}]")

# Mapa kolorow
fig2, ax2 = plt.subplots(figsize=(10, 8))

G_clip = np.clip(G_2d, -5, 5)
# Kontur g=0
try:
    cs = ax2.contourf(K_2d, A_2d, G_clip,
                      levels=np.linspace(-5, 5, 41),
                      cmap='RdBu_r', extend='both')
    plt.colorbar(cs, ax=ax2, label='g(K, a_Gamma) = E/(4*pi*K) - 1')
    ax2.contour(K_2d, A_2d, G_2d,
                levels=[0.0], colors='black', linewidths=2.5,
                linestyles='-')
except Exception as e:
    print(f"  Blad konturu: {e}")

# Linie skalowania
K_line = np.logspace(-3, 2, 200)
# Schwarzschild: a_Gam = C1 * K,  C1 = a_Gam1/K1* = 0.040/0.010414 = 3.840
C1 = 0.040 / 0.010414
ax2.plot(K_line, C1 * K_line, 'g--', lw=1.5, label=f'Schwarzschild: a={C1:.2f}*K')
# Compton: a_Gam = C_c / K (jesli a ~ hbar/(mc) ~ 1/K)
C_c = 0.040 * 0.010414
ax2.plot(K_line, C_c / K_line, 'm--', lw=1.5, label=f'Compton: a={C_c:.4f}/K')
# Stalae a_Gam = 0.040
ax2.axhline(0.040, color='orange', lw=1.5, linestyle=':', label='a_Gam=0.040 (stale)')

# Zaznacz znane punkty
ax2.scatter([0.010414], [0.040], s=120, color='black', zorder=5, label='K*_1(ODE)')
ax2.scatter([2.032728, 34.14450], [0.040, 0.040],
            s=100, color='red', marker='^', zorder=5, label='K_2,K_3 Yukawa')

ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlabel('K', fontsize=12)
ax2.set_ylabel('a_Gamma', fontsize=12)
ax2.set_title(f'Mapa g(K, a_Gamma) — kontur g=0 (czarny) to samospojne solitony TGP\n'
              f'(alpha={ALPHA}, lam={LAM:.2e})', fontsize=10)
ax2.legend(fontsize=9, loc='upper left')
ax2.grid(True, alpha=0.3)

out2 = __file__.replace('.py', '_map2D.png')
plt.savefig(out2, dpi=130, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out2}")
print()


# ============================================================
# 3. SKALOWANIE KONTURU g=0
# ============================================================
print("=" * 70)
print("ANALIZA KONTURU g=0:")
print()

# Dla kazdego K na siatce 2D: znajdz a_Gam* gdzie g=0 (interpolacja)
K_contour = []
A_contour = []

for j, K in enumerate(K_2d):
    g_col = G_2d[:, j]
    for i in range(len(g_col) - 1):
        gi, gj = g_col[i], g_col[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            # Liniowa interpolacja w log-a
            t = -gi / (gj - gi)
            a_star = A_2d[i] * (A_2d[i+1]/A_2d[i])**t
            K_contour.append(K)
            A_contour.append(a_star)
            print(f"  K={K:.4e}  a_Gam*={a_star:.4f}  a_Gam*/K={a_star/K:.4f}  "
                  f"a_Gam*^2/K={a_star**2/K:.4f}")

K_c = np.array(K_contour)
A_c = np.array(A_contour)

if len(K_c) >= 3:
    # Dopasuj potege: a ~ C * K^sigma
    lnK = np.log(K_c)
    lnA = np.log(A_c)
    sigma, lnC = np.polyfit(lnK, lnA, 1)
    C_fit = np.exp(lnC)
    print()
    print(f"  Dopasowanie potegowe a_Gam* ~ C * K^sigma:")
    print(f"  sigma = {sigma:.4f}   C = {C_fit:.4f}")
    print(f"  Porownanie: sigma=1 (Schwarzschild), sigma=0 (stale a_Gam)")
    print(f"  Schwarzschild C1 = {C1:.4f}")
    print()

    # Rysuj kontur osobno
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    ax3.scatter(K_c, A_c, s=50, color='blue', zorder=5, label='g=0 kontur (z mapy 2D)')
    K_fit_line = np.logspace(np.log10(K_c.min()), np.log10(K_c.max()), 100)
    ax3.plot(K_fit_line, C_fit * K_fit_line**sigma, 'r--', lw=2,
             label=f'Dopasowanie: a={C_fit:.3f}*K^{sigma:.3f}')
    ax3.plot(K_fit_line, C1 * K_fit_line, 'g:', lw=1.5,
             label=f'Schwarzschild: a={C1:.3f}*K')
    ax3.scatter([0.010414], [0.040], s=120, color='black', zorder=6,
                label='K*1(ODE) potwierdzony')
    ax3.set_xscale('log'); ax3.set_yscale('log')
    ax3.set_xlabel('K')
    ax3.set_ylabel('a_Gamma* (samospojne)')
    ax3.set_title('Kontur g=0: skalowanie a_Gamma*(K)')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    out3 = __file__.replace('.py', '_contour.png')
    plt.savefig(out3, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Zapisano: {out3}")


# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE p10")
print()
print("Pytanie: Czy K2=2.033 i K3=34.14 moga byc samospojnymi solitonami ODE?")
print()
print("  Wyniki skanu 1D g(K, a_Gam):")
for label, (valid_a, g_arr) in scan_results.items():
    fin = g_arr[np.isfinite(g_arr)]
    if len(fin) > 0:
        crosses = sum(1 for i in range(len(g_arr)-1)
                      if np.isfinite(g_arr[i]) and np.isfinite(g_arr[i+1])
                      and g_arr[i]*g_arr[i+1] < 0)
        print(f"  {label}: g in [{fin.min():.2e}, {fin.max():.2e}], "
              f"zmian znaku: {crosses}")
    else:
        print(f"  {label}: brak wynikow")
print()
print("  Kontur g=0 w 2D:")
if len(K_c) > 0:
    print(f"  Punkty konturu: {len(K_c)}")
    print(f"  K zakres: [{K_c.min():.3e}, {K_c.max():.3e}]")
    print(f"  a_Gam* zakres: [{A_c.min():.4f}, {A_c.max():.3f}]")
    if len(K_c) >= 3:
        print(f"  Skalowanie: a_Gam* ~ {C_fit:.3f} * K^{sigma:.3f}")
        if abs(sigma - 1.0) < 0.15:
            print("  -> PRAWIE Schwarzschild (sigma~1)")
        elif abs(sigma) < 0.15:
            print("  -> PRAWIE stalye a_Gam (sigma~0)")
        else:
            print(f"  -> Niestandardowe skalowanie (sigma={sigma:.3f})")
else:
    print("  Brak punktow konturu — g nie zmienia znaku w przebadanym zakresie!")
    print("  WNIOSEK: g(K, a_Gam) < 0 WSZEDZIE (lub nan)")
    print("  -> Jedyny samospojny soliton to K*_1 przy a_Gam=0.040")
    print("  -> Teoria TGP z jednym potencjalem V_mod i stalym lambda*")
    print("     NIE ma multi-generacyjnych solitonow ODE poza K*_1!")
