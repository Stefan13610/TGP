"""
TGP — Poprawiony wzor analityczny na r21
=========================================

POPRAWA 1 (bledna): V''(1) = -1 (czlon kwadratowy Ep) => -pi*K poprawka
    Wynik: K2 maleje (ZLY KIERUNEK)

POPRAWA 2 (prawidlowa): Czlon kubiczny V = psi^3/3 w Ep jest POZYTYWNY.
    Oryginalnie uwzgledniany tylko kwartyczny czlon (Regime C):
        Ep_quart/K ~ -pi*K^3/a_Gam

    Czlon kubiczny psi^3/3 rowniez wystepuje, dajac wklad:
        Ep_cubic/K ~ +4*pi*K^2 * ln(1/a_Gam) / 3

    (calka log-rozbieza przy malym r, przycita na a_Gam)

Poprawione rownanie (4 wyrazy, numerycznie):
    K^3 - eps*K^2 - 2*K = c
    gdzie:
        c   = a_Gam * (2*alpha - 4)
        eps = (4*a_Gam / 3) * ln(1/a_Gam)    [> 0]

    Rozwiazywane brentq (brak prostego wzoru Cardano dla pelnej cubiki).

Fizyczne uzasadnienie:
    Przy K ~ K2 ~ 1.5-2 (Regime C), pole w rdzeniu phi >> 1.
    Potencjal V(psi) = psi^3/3 - psi^4/4. Dla duzego psi czlon psi^3/3
    jest dodatni i zmniejsza |E_pot|, przesuwajac K2 w gore.
    Calka int_aGam^1 (K/r)^3 r^2 dr = K^3 * ln(1/aGam) daje log-korekcje.
"""

import numpy as np
from scipy.optimize import brentq

GAMMA = 1.0
LAM = 5.514e-7  # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)

# ======================================================================
# ENERGIA NUMERYCZNA (referencyjna)
# ======================================================================

def energy_num(K, Phi0, alpha, a_gam):
    msp = np.sqrt(max(GAMMA / Phi0, 1e-9))
    r_max = max(60.0 / msp, 20.0)
    r = np.linspace(a_gam, r_max, 3000)
    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1.0 + alpha/(Phi0*phi)) * r**2, r)
    psi = phi / Phi0
    Ep  = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4 + LAM/6*(psi-1.0)**6) * r**2, r)
    return Ek + Ep

def find_crossings_num(alpha, a_gam, K_min=1e-4, K_max=5.0, N=800):
    """Oba przeciecia E(K)/K = 4pi numerycznie."""
    K_arr = np.linspace(K_min, K_max, N)
    f_arr = np.array([energy_num(K, 1.0, alpha, a_gam)/K - 4*np.pi for K in K_arr])
    roots = []
    for i in range(len(f_arr)-1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi*fj < 0:
            try:
                root = brentq(
                    lambda K: energy_num(K, 1.0, alpha, a_gam)/K - 4*np.pi,
                    K_arr[i], K_arr[i+1], xtol=1e-9)
                roots.append(root)
            except Exception:
                pass
    return roots


# ======================================================================
# WZORY ANALITYCZNE — K1, K2, r21
# ======================================================================

def K1_analytic(a_gam, alpha):
    """K1 ~ 4*a / (2*(1+alpha) - a)  [dokladnosc ~1%]"""
    denom = 2*(1 + alpha) - a_gam
    return 4*a_gam / denom if denom > 0 else None


# --- ORYGINALNE: K^3 - 2*K = c ---

def K2_original_cardano(a_gam, alpha):
    """Cardano dla K^3 - 2*K = c."""
    c = a_gam * (2*alpha - 4)
    # Viete: K^3 + p*K + q = 0, p=-2, q=-c
    Delta = 32.0 - 27.0*c**2        # 4*|-p|^3 - 27*q^2
    if Delta <= 0:
        return None
    m   = 2.0 * np.sqrt(2.0/3.0)   # 2*sqrt(-p/3)
    arg = 3.0*np.sqrt(1.5)*c / 4.0  # 3q / (p*m) -> 3*sqrt(3/2)*c / 4
    if abs(arg) > 1.0:
        return None
    theta = np.arccos(arg) / 3.0
    # trzy korzenie; K2 = Max dodatni (k=0)
    roots = [m * np.cos(theta - 2*np.pi*k/3) for k in range(3)]
    pos = [r for r in roots if r > 0.3]
    return max(pos) if pos else None


# --- POPRAWKA 1 (bledna): K^3 - (2-a)*K = c ---

def K2_quad_correction(a_gam, alpha):
    """Cardano dla K^3 - (2-a)*K = c.
    Uwzglednia V''(1)=-1 => -pi*K korekcja.
    WYNIK: K2 maleje — zly kierunek!"""
    c = a_gam * (2*alpha - 4)
    p_eff = 2.0 - a_gam
    Delta = 4.0*p_eff**3 - 27.0*c**2
    if Delta <= 0:
        return None
    m   = 2.0 * np.sqrt(p_eff / 3.0)
    arg = 3.0 * np.sqrt(3.0) * c / (2.0 * p_eff**1.5)
    if abs(arg) > 1.0:
        return None
    theta = np.arccos(arg) / 3.0
    roots = [m * np.cos(theta - 2*np.pi*k/3) for k in range(3)]
    pos = [r for r in roots if r > 0.3]
    return max(pos) if pos else None


# --- POPRAWKA 2 (prawidlowa): K^3 - eps*K^2 - 2*K = c ---

def K2_cubic_correction(a_gam, alpha):
    """Rozwiazanie numeryczne K^3 - eps*K^2 - 2*K = c.
    eps = (4*a_Gam/3) * ln(1/a_Gam)  [czlon kubiczny V = psi^3/3]
    Korekcja: Ep_cubic/K = +4*pi*K^2*ln(1/a)/3 (pozytywna => K2 rosnie)."""
    c   = a_gam * (2*alpha - 4)
    eps = (4.0 * a_gam / 3.0) * np.log(1.0 / a_gam)
    def f(K):
        return K**3 - eps*K**2 - 2.0*K - c
    # K2 lezy na galezi rosnacej, typowo K2 in [0.5, 6]
    try:
        return brentq(f, 0.5, 6.0, xtol=1e-10)
    except ValueError:
        return None


def r21(K2_fn, a_gam, alpha):
    K1 = K1_analytic(a_gam, alpha)
    K2 = K2_fn(a_gam, alpha)
    if K1 and K2 and K1 > 0:
        return K2 / K1
    return None


# ======================================================================
# WERYFIKACJA — tabela
# ======================================================================

print("=" * 105)
print("TRZY WERSJE: Oryginal | Poprawka-1 (V''quad) | Poprawka-2 (V-cubic log)")
print("  Oryginal:   K^3 - 2*K        = c    [c = a*(2*alpha-4)]")
print("  Poprawka-1: K^3 - (2-a)*K   = c    [V''(1)=-1 => -pi*K, ZLY kierunek K2]")
print("  Poprawka-2: K^3 - eps*K^2 - 2*K = c [eps=(4a/3)*ln(1/a), czlon psi^3/3]")
print("=" * 105)
header = (f"{'alpha':>6} {'a_gam':>6} | "
          f"{'r21_num':>8} | "
          f"{'r21_orig':>8} {'bl_orig':>7} | "
          f"{'r21_p1':>8} {'bl_p1':>7} | "
          f"{'r21_p2':>8} {'bl_p2':>7} | eps")
print(header)
print("-" * 105)

cases = [
    (4.0, 0.03), (4.0, 0.05),
    (5.9, 0.03), (5.9, 0.05),
    (6.0, 0.03), (6.0, 0.05),
    (7.0, 0.03), (7.0, 0.05),
    (8.0, 0.03), (8.0, 0.05),
    (10.0, 0.03), (10.0, 0.05),
]

for alpha, a_gam in cases:
    roots = find_crossings_num(alpha, a_gam)
    if len(roots) < 2:
        print(f"{alpha:>6.1f} {a_gam:>6.3f} | (brak 2 korzeni)")
        continue
    K1_n, K2_n = sorted(roots)[:2]
    r_num = K2_n / K1_n

    r_orig = r21(K2_original_cardano,  a_gam, alpha)
    r_p1   = r21(K2_quad_correction,   a_gam, alpha)
    r_p2   = r21(K2_cubic_correction,  a_gam, alpha)
    eps    = (4.0*a_gam/3.0) * np.log(1.0/a_gam)

    def fmt(r):
        if r is None: return f"{'N/A':>8} {'N/A':>7}"
        err = 100*(r - r_num)/r_num
        return f"{r:>8.1f} {err:>+7.1f}%"

    print(f"{alpha:>6.1f} {a_gam:>6.3f} | {r_num:>8.1f} | {fmt(r_orig)} | {fmt(r_p1)} | {fmt(r_p2)} | {eps:.3f}")


# ======================================================================
# SZCZEGOLY K2: gdzie lezy blad?
# ======================================================================

print()
print("=" * 95)
print("SZCZEGOLY K2: Poprawka-2 (czlon psi^3/3) vs numeryczny")
print("=" * 95)
header2 = (f"{'alpha':>6} {'a_gam':>6} | "
           f"{'K1_num':>9} {'K1_an':>9} {'bl_K1':>7} | "
           f"{'K2_num':>9} {'K2_orig':>9} {'bl_K2o':>7} | "
           f"{'K2_p2':>9} {'bl_K2p2':>8}")
print(header2)
print("-" * 95)

for alpha, a_gam in cases:
    roots = find_crossings_num(alpha, a_gam)
    if len(roots) < 2:
        continue
    K1_n, K2_n = sorted(roots)[:2]
    K1_a  = K1_analytic(a_gam, alpha)
    K2_o  = K2_original_cardano(a_gam, alpha)
    K2_p2 = K2_cubic_correction(a_gam, alpha)

    bl_K1   = 100*(K1_a  - K1_n)/K1_n if K1_a  else float('nan')
    bl_K2o  = 100*(K2_o  - K2_n)/K2_n if K2_o  else float('nan')
    bl_K2p2 = 100*(K2_p2 - K2_n)/K2_n if K2_p2 else float('nan')

    print(f"{alpha:>6.1f} {a_gam:>6.3f} | "
          f"{K1_n:>9.5f} {K1_a:>9.5f} {bl_K1:>+7.2f}% | "
          f"{K2_n:>9.5f} {K2_o:>9.5f} {bl_K2o:>+7.2f}% | "
          f"{K2_p2:>9.5f} {bl_K2p2:>+8.2f}%")


# ======================================================================
# WZOR ZAMKNIETY NA r21 (poprawiony)
# ======================================================================

print()
print("=" * 80)
print("WZOR POPRAWIONY (Poprawka-2):")
print()
print("  Czlon kubiczny potencjalu V(psi) = psi^3/3 - psi^4/4 daje:")
print("    Ep_cubic/K = +4*pi*K^2 * ln(1/a_Gam) / 3  [calka log-rozbieza]")
print()
print("  Poprawione rownienie:")
print("    K^3 - eps*K^2 - 2*K = c")
print("    eps = (4*a_Gam/3) * ln(1/a_Gam),   c = a_Gam*(2*alpha-4)")
print()
print("  K2 rozwiazywane numerycznie (brentq).")
print("  K1 bez zmian: K1 = 4*a_Gam / (2*(1+alpha) - a_Gam)")
print()
print("  Przy a_Gam -> 0: eps -> 0, formula -> oryginalna.")
print()

for alpha, a_gam in [(5.9, 0.03), (7.0, 0.03), (10.0, 0.03)]:
    c   = a_gam * (2*alpha - 4)
    eps = (4.0*a_gam/3.0)*np.log(1.0/a_gam)
    K2p = K2_cubic_correction(a_gam, alpha)
    K1a = K1_analytic(a_gam, alpha)
    r_p = K2p/K1a if (K2p and K1a) else None
    roots = find_crossings_num(alpha, a_gam)
    r_n = sorted(roots)[1]/sorted(roots)[0] if len(roots)>=2 else None
    s_n = sorted(roots)
    print(f"  alpha={alpha:.1f}, a_gam={a_gam:.3f}: eps={eps:.4f}")
    if roots:
        print(f"    K2_p2={K2p:.5f} (num={s_n[1]:.5f}, bl={100*(K2p-s_n[1])/s_n[1]:+.1f}%)")
        print(f"    K1_an={K1a:.5f} (num={s_n[0]:.5f}, bl={100*(K1a-s_n[0])/s_n[0]:+.1f}%)")
    if r_p and r_n:
        print(f"    r21_p2={r_p:.1f}   (num={r_n:.1f}, bl={100*(r_p-r_n)/r_n:+.1f}%)")
    print()
print()
