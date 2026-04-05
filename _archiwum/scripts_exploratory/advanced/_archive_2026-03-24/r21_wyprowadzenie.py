"""
TGP — Analityczne wyprowadzenie r21 = M2/M1
============================================

Strategia: analiza E(K)/K w trzech rezimach K.
Dwa przeciecia E/K = Lambda = 4pi daja K1 < K2,
stad M1 = 4pi*Phi0*K1,  M2 = 4pi*Phi0*K2,  r21 = K2/K1.

Wyprowadzamy K1 i K2 analitycznie, weryfikujemy numerycznie.
"""
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAMMA = 1.0; PHI0 = 1.0; MSP = 1.0
LAM = 5.514e-7                 # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)
LAMBDA = 4.0 * np.pi          # stala samospojnosci

# ======================================================================
# DOKLADNA E(K)/K (numerycznie, do weryfikacji)
# ======================================================================
def E_over_K_exact(K, a_gam, alpha, Phi0=1.0, msp=1.0, N=3000):
    r_max = max(30.0/msp, 10.0)
    r  = np.linspace(a_gam, r_max, N)
    ph = Phi0 + K * np.exp(-msp*r) / r
    ph = np.maximum(ph, 1e-10)
    dp = K * np.exp(-msp*r) * (-msp*r - 1.0) / r**2
    Ek = 4*np.pi * np.trapezoid(0.5*dp**2*(1 + alpha/(Phi0*ph)) * r**2, r)
    psi = ph / Phi0
    Ep  = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4 + LAM/6*(psi-1.0)**6) * r**2, r)
    return (Ek + Ep) / K

def find_crossings_exact(a_gam, alpha, K_max=4.0, N=600):
    Kv = np.concatenate([
        np.linspace(0.001, 0.1,   int(N*0.3)),
        np.linspace(0.1,   1.0,   int(N*0.4)),
        np.linspace(1.0,   K_max, int(N*0.3))])
    Kv = np.unique(Kv)
    fv = np.array([E_over_K_exact(K, a_gam, alpha) - LAMBDA for K in Kv])
    roots = []
    for i in range(len(fv)-1):
        if np.isfinite(fv[i]) and np.isfinite(fv[i+1]) and fv[i]*fv[i+1] < 0:
            try:
                Kc = brentq(lambda K: E_over_K_exact(K, a_gam, alpha) - LAMBDA,
                            Kv[i], Kv[i+1], xtol=1e-9)
                roots.append(Kc)
            except Exception:
                pass
    return roots, Kv, fv

# ======================================================================
# ANALIZA REZIMU A (male K): E ~ 2pi*K^2*(1+alpha)/a_gam
# Energia kinetyczna UV-dominowana przez rdzen
# ======================================================================
# Ek_std ~ integral K^2*exp(-2s)*(1+s)^2/s^2 ds ze skonczonym limitem a_gam
#        ~ 2*pi*K^2 * sqrt(Phi0)/a_gam  (dominacja UV 1/s^2)
#
# Ek_nl  ~ alpha/Phi0^2 * Ek_std  (dla K<<Phi0*a_gam)
#
# Ep_quad ~ -pi*K^2*sqrt(Phi0)  (czlon V''(1)=-1 podcalkowy)
#
# Dla Phi0=1, msp=1:
#   E/K ~ 2*pi*K*(1+alpha)/a_gam - pi*K = pi*K*(2*(1+alpha)/a_gam - 1)
#
# Pierwsze przejscie: E/K = 4*pi  =>
#   K1_analytic = 4*a_gam / (2*(1+alpha) - a_gam)
#
# Dla alpha >> 1, a_gam << 1:
#   K1 ~ 2*a_gam / (1+alpha)                              [formula K1]

def K1_analytic(a_gam, alpha):
    """Wzor analityczny na K1 (wiodace porzadkowanie)."""
    return 4*a_gam / (2*(1+alpha) - a_gam)

# ======================================================================
# ANALIZA REZIMU C (duze K): E ~ 2pi*K^2/a_gam + 2pi*alpha*K - pi*K^4/a_gam
# Czlon kwartyczny potencjalu dominuje przy duzych K
# ======================================================================
# dPhi/dr ~ -K/r^2 przy r~a_gam
# phi_core ~ K/a_gam >> 1  dla K >> Phi0*a_gam
#
# V(phi_core/Phi0) - V(1) ~ -(K/a_gam)^4/4  (dominacja phi^4)
# Ep_quart ~ -4pi*a_gam^3 * (K/a_gam)^4/4 = -pi*K^4/a_gam
#
# Ek_nl_core (dla phi_core >> Phi0):
#   alpha/(Phi0*phi_core) ~ alpha*a_gam/K  <<  zintegrowane: ~ 2pi*alpha*K
#
# Laczna E/K przy duzych K:
#   E/K ~ 2pi*K/a_gam + 2pi*alpha - pi*K^3/a_gam
#
# Drugie przejscie: E/K = 4*pi  =>
#   2*K/a_gam + 2*alpha - K^3/a_gam = 4
#   K^3 - 2*K = a_gam*(2*alpha - 4)   [rownanie kubiczne na K2]

def K2_analytic(a_gam, alpha):
    """
    Analityczny K2 z rownania kubicznego:
       K^3 - 2*K = a_gam*(2*alpha - 4)

    Rozwiazanie przez wzor Cardano lub numerycznie.
    """
    # K^3 - 2K - c = 0, c = a_gam*(2*alpha-4)
    c = a_gam * (2*alpha - 4)
    # Szukamy korzenia powyzej sqrt(2/3) (minimum K^3-2K)
    try:
        # f(K) = K^3 - 2K - c = 0, K > 1
        return brentq(lambda K: K**3 - 2*K - c, 0.5, 10.0, xtol=1e-10)
    except Exception:
        return None

# Wzor Cardano dla K^3 + pK + q = 0:
# K^3 - 2K - c = 0  =>  p=-2, q=-c
# Wyroznik: Delta = -(4p^3 + 27q^2) = -(4*(-2)^3 + 27*c^2) = 32 - 27c^2
# Jesli Delta > 0: trzy rzeczywiste korzenie
# Jesli Delta < 0: jeden rzeczywisty korzen (duze K rozwiazanie)

def K2_cardano(a_gam, alpha):
    """Wzor Cardano dla K^3 - 2K = c."""
    c = a_gam * (2*alpha - 4)
    # K^3 - 2K - c = 0 <=> K^3 + pK + q = 0, p=-2, q=-c
    p, q = -2.0, -c
    Delta = -(4*p**3 + 27*q**2)  # dyskryminant
    if Delta > 0:
        # 3 rzeczywiste korzenie: forma tryg
        m = 2*np.sqrt(-p/3)
        theta = np.arccos(3*q/(p*m))/3
        # Trzy korzenie: m*cos(theta), m*cos(theta+2pi/3), m*cos(theta+4pi/3)
        roots = [m*np.cos(theta + 2*k*np.pi/3) for k in range(3)]
        positive = [r for r in roots if r > 0]
        return max(positive) if positive else None
    else:
        # 1 rzeczywisty korzen: formula Cardano
        A = np.cbrt(-q/2 + np.sqrt(q**2/4 + p**3/27))
        B = np.cbrt(-q/2 - np.sqrt(q**2/4 + p**3/27))
        return A + B

# ======================================================================
# FORMULA NA r21
# ======================================================================
# r21 = K2/K1 = K2 * (2*(1+alpha) - a_gam) / (4*a_gam)
#
# Dla alpha >> 1:
#   K1 ~ 2*a_gam/(1+alpha)
#   K2 ~ root of K^3 - 2K = 2*alpha*a_gam
#   r21 ~ K2*(1+alpha)/(2*a_gam)

def r21_analytic(a_gam, alpha):
    K1 = K1_analytic(a_gam, alpha)
    K2 = K2_analytic(a_gam, alpha)
    if K1 and K2 and K1 > 0:
        return K2/K1
    return None

# ======================================================================
# WERYFIKACJA NUMERYCZNA
# ======================================================================
print("="*70)
print("WERYFIKACJA: K1_analytic vs K1_numeryczny")
print("="*70)
print(f"{'alpha':>7} {'a_gam':>7} {'K1_num':>10} {'K1_an':>10} {'err%':>8} "
      f"{'K2_num':>10} {'K2_an':>10} {'err%':>8}")
print("-"*70)

for alpha in [4.0, 5.9, 6.0, 8.0, 10.0]:
    for a_gam in [0.03, 0.05, 0.08, 0.10]:
        roots, _, _ = find_crossings_exact(a_gam, alpha)
        if len(roots) < 2:
            continue
        K1n, K2n = roots[0], roots[1]
        K1a = K1_analytic(a_gam, alpha)
        K2a = K2_analytic(a_gam, alpha)
        err1 = 100*(K1a - K1n)/K1n
        err2 = 100*(K2a - K2n)/K2n if K2a else None
        s2a = f"{K2a:.5f}" if K2a else "-"
        se2 = f"{err2:+.1f}%" if err2 else "-"
        print(f"{alpha:>7.1f} {a_gam:>7.3f} {K1n:>10.6f} {K1a:>10.6f} {err1:>+7.1f}% "
              f"{K2n:>10.5f} {s2a:>10} {se2:>8}")

# ======================================================================
# r21 ANALITYCZNY vs NUMERYCZNY
# ======================================================================
print()
print("="*70)
print("r21 = K2/K1: ANALITYCZNY vs NUMERYCZNY")
print("="*70)
print(f"{'alpha':>7} {'a_gam':>7} {'r21_num':>10} {'r21_an':>10} {'err%':>8} "
      f"{'r21_asymp':>12}")
print("-"*70)

results = []
for alpha in [4.0, 5.0, 5.9, 6.0, 7.0, 8.0, 10.0]:
    for a_gam in [0.03, 0.05, 0.08, 0.10]:
        roots, _, _ = find_crossings_exact(a_gam, alpha)
        if len(roots) < 2:
            continue
        K1n, K2n = roots[0], roots[1]
        r21n = K2n/K1n
        r21a = r21_analytic(a_gam, alpha)
        # Asymptotyczna formula: r21 ~ sqrt(2)*(1+alpha)/(2*a_gam)
        r21_asymp = np.sqrt(2)*(1+alpha)/(2*a_gam)
        err = 100*(r21a - r21n)/r21n if r21a else None
        sa = f"{r21a:.2f}" if r21a else "-"
        se = f"{err:+.1f}%" if err else "-"
        results.append({'alpha': alpha, 'a_gam': a_gam,
                        'r21n': r21n, 'r21a': r21a, 'r21_asymp': r21_asymp})
        print(f"{alpha:>7.1f} {a_gam:>7.3f} {r21n:>10.2f} {sa:>10} {se:>8} "
              f"{r21_asymp:>12.1f}")

# ======================================================================
# WZOR ZAMKNIETY NA r21
# ======================================================================
print()
print("="*70)
print("WZOR ZAMKNIETY NA r21:")
print("="*70)
print("""
  K1 ~ 4*a_Gam / (2*(1+alpha) - a_Gam)  ~  2*a_Gam / (1+alpha)

  K2 = wiekszy korzen rownania kubicznego:
       K^3 - 2*K = a_Gam*(2*alpha - 4)

  Wzor Cardano: K2 = 2*sqrt(2/3)*cos[ (1/3)*arccos(3*c/(p*m)) ]
     gdzie c = a_Gam*(2*alpha-4),  m = 2*sqrt(2/3)

  r21 = K2/K1 ~ K2*(1+alpha) / (2*a_Gam)

  Asymptotycznie dla alpha >> 1, a_Gam*alpha >> 1:
       K2 ~ (2*alpha*a_Gam)^(1/3)        [czlon kubiczny dominuje]
       r21 ~ (2*alpha*a_Gam)^(1/3) * (1+alpha) / (2*a_Gam)
           ~ (alpha/(2*a_Gam^2))^(1/3) * alpha/2
           ~ alpha^(4/3) / (2*a_Gam)^(2/3)
""")

print("  Weryfikacja wzoru asymptotycznego:")
print(f"  {'alpha':>7} {'a_gam':>7} {'r21_num':>10} {'asym_kubiczny':>15} {'err%':>8}")
for r in results:
    alpha, a_gam = r['alpha'], r['a_gam']
    c = a_gam*(2*alpha-4)
    # Asymptota dla duzego c: K2 ~ c^(1/3)
    K2_asy = c**(1/3) if c > 0 else None
    if K2_asy:
        K1_asy = 2*a_gam/(1+alpha)
        r21_asy = K2_asy/K1_asy
        err = 100*(r21_asy - r['r21n'])/r['r21n']
        print(f"  {alpha:>7.1f} {a_gam:>7.3f} {r['r21n']:>10.2f} {r21_asy:>15.2f} {err:>+8.1f}%")

# ======================================================================
# JAKA KOMBINACJA (alpha, a_Gam) DAJE r21=207?
# ======================================================================
print()
print("="*70)
print("SKAN: dla jakiego (alpha, a_Gam) r21_analytic = 207?")
print("="*70)
print(f"  {'alpha':>7} {'a_gam':>7} {'r21_an':>10} {'r21_num':>10} {'roznica':>10}")
target = 207.0
for alpha in np.linspace(4.0, 12.0, 9):
    for a_gam in np.linspace(0.02, 0.12, 6):
        r21a = r21_analytic(a_gam, alpha)
        if r21a is None: continue
        # Sprawdz bliskosc do 207
        if abs(r21a - target)/target < 0.15:
            roots, _, _ = find_crossings_exact(a_gam, alpha)
            r21n = roots[1]/roots[0] if len(roots) >= 2 else None
            s_n = f"{r21n:.2f}" if r21n else "-"
            diff = f"{100*(r21a-target)/target:+.1f}%" if r21a else "-"
            print(f"  {alpha:>7.2f} {a_gam:>7.3f} {r21a:>10.2f} {s_n:>10} {diff:>10}")

# ======================================================================
# WYKRES: E(K)/K i dwa przeciecia
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

alpha_plot = 5.9; a_gam_plot = 0.05
roots, Kv, fv = find_crossings_exact(a_gam_plot, alpha_plot, K_max=3.0, N=500)
EKv = fv + LAMBDA

ax = axes[0]
mask = np.isfinite(EKv) & (EKv > -50) & (EKv < 200)
ax.plot(Kv[mask], EKv[mask], 'b-', lw=2.5, label=r'$E(K)/K$ (numeryczne)')
ax.axhline(LAMBDA, color='r', ls='--', lw=2, label=rf'$\Lambda = 4\pi = {LAMBDA:.3f}$')

# Analityczne K1, K2
K1a = K1_analytic(a_gam_plot, alpha_plot)
K2a = K2_analytic(a_gam_plot, alpha_plot)
ax.axvline(K1a, color='green', ls=':', lw=2, label=f'K1_an={K1a:.5f}')
ax.axvline(K2a, color='orange', ls=':', lw=2, label=f'K2_an={K2a:.4f}')

if len(roots) >= 2:
    ax.axvline(roots[0], color='darkgreen', ls='-', lw=1.5, label=f'K1_num={roots[0]:.5f}')
    ax.axvline(roots[1], color='darkorange', ls='-', lw=1.5, label=f'K2_num={roots[1]:.4f}')
    r21n = roots[1]/roots[0]
    ax.set_title(f'alpha={alpha_plot}, a_Gam={a_gam_plot}\n'
                 f'r21_num={r21n:.1f}, r21_an={r21_analytic(a_gam_plot, alpha_plot):.1f}',
                 fontsize=10)

ax.set_xlabel(r'$K$', fontsize=12)
ax.set_ylabel(r'$E(K)/K$', fontsize=12)
ax.set_ylim(-20, 120)
ax.set_xlim(0, 2.5)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 2: r21 vs alpha dla roznych a_Gam
ax2 = axes[1]
alpha_arr = np.linspace(3.0, 15.0, 60)
colors2 = ['b', 'g', 'r', 'purple']
for ci, a_g in enumerate([0.03, 0.05, 0.08, 0.10]):
    r21_arr = [r21_analytic(a_g, a) for a in alpha_arr]
    r21_arr = [r if r else np.nan for r in r21_arr]
    ax2.plot(alpha_arr, r21_arr, '-', color=colors2[ci], lw=2,
             label=f'a_Gam={a_g}')

ax2.axhline(207, color='k', ls='--', lw=1.5, label='r21=207 (cel)')
ax2.set_xlabel(r'$\alpha$', fontsize=12)
ax2.set_ylabel(r'$r_{21} = M_2/M_1$', fontsize=12)
ax2.set_title('r21 analityczny vs alpha', fontsize=10)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, 600)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/r21_wyprowadzenie.png',
            dpi=150, bbox_inches='tight')

# ======================================================================
# DOKLADNA FORMULA: Cardano dla K2
# ======================================================================
print()
print("="*70)
print("WZOR CARDANO (dokladny dla K2):")
print("="*70)
print("""
  Rownanie: K^3 - 2K - c = 0,  c = a_Gam*(2*alpha - 4)

  Dysryminant: Delta = 32 - 27*c^2

  Gdy Delta > 0 (trzy rzeczywiste korzenie — 3 generacje!):
       K = (2/sqrt(3)) * cos[ (1/3)*arccos( (3*sqrt(3)/4)*c ) - 2pi*k/3 ]
         dla k = 0, 1, 2

  K1 (maly): k=2  =>  K1_from_cubic = (2/sqrt(3))*cos[(1/3)*arccos(...) + 4pi/3]
  K2 (duzy): k=0  =>  K2_from_cubic = (2/sqrt(3))*cos[(1/3)*arccos(...)]

  WARUNEK 3 GENERACJI: Delta > 0  =>  c^2 < 32/27
       |a_Gam*(2*alpha-4)| < sqrt(32/27) ~ 1.089

  Maksymalny r21 przy c -> sqrt(32/27) (granica 3 generacji):
       K2_max = 2/sqrt(3) = 2*sqrt(3)/3 ~ 1.155
       K1_min = 4*a_Gam/(2*(1+alpha)-a_Gam) ~ 2*a_Gam/(1+alpha)
       r21_max ~ (1+alpha)*sqrt(3)/a_Gam
""")

for alpha in [5.9, 6.0, 8.0, 10.0]:
    for a_gam in [0.03, 0.05, 0.08]:
        c = a_gam*(2*alpha-4)
        Delta = 32 - 27*c**2
        K2c = K2_cardano(a_gam, alpha)
        K1a = K1_analytic(a_gam, alpha)
        r21c = K2c/K1a if K2c else None
        sc = f"{'3gen' if Delta>0 else '1gen'}"
        sr = f"{r21c:.2f}" if r21c else "-"
        print(f"  alpha={alpha:.1f} a_Gam={a_gam:.3f}: c={c:.4f}  Delta={Delta:.3f}  "
              f"K2_Cardano={K2c:.5f}  r21={sr}  [{sc}]")

# Maksymalne r21 w limicie Delta=0:
print()
print("  Graniczne r21 (Delta=0, maksymalne K2=2/sqrt(3)):")
K2_limit = 2/np.sqrt(3)
for alpha in [5.9, 6.0, 8.0, 10.0]:
    for a_gam in [0.03, 0.05, 0.08]:
        K1a = K1_analytic(a_gam, alpha)
        r21_lim = K2_limit / K1a
        a_Gam_crit = np.sqrt(32/27) / (2*alpha-4) if alpha > 2 else None
        print(f"    alpha={alpha:.1f} a_Gam={a_gam:.3f}: r21_max={r21_lim:.1f}")

print()
print("Wykres: r21_wyprowadzenie.png")
print("Gotowe.")
