"""
v3_regularyzacja.py — Regularyzowany czlon kwartyczny

Pomysl: zamiast dodawac osobny czlon stabilizujacy lambda*(psi-1)^6
(ktory tworzyl katastrofalne nowe minimum i numeryczna niestabilnosc),
REGULARYZUJEMY sama kwartyczne czlon potencjalu:

    V_mod(psi) = psi^3/3 - psi^4 / [4*(1 + eps*(psi-1)^3)]

Wlasnosci:
    V_mod(1) = 1/3 - 1/4 = 1/12           (proznia bez zmian)
    V'_mod(1) = 0                          (minimum prozni zachowane)
    V''_mod(1) = -1                        (masa skalaru bez zmian)

Dla duzych psi:
    1 + eps*(psi-1)^3 ~ eps*psi^3
    V_mod ~ psi^3/3 - psi/(4*eps)  => rosnie jak psi^3!

    => Brak katastrofalnego minimum, V_mod ograniczony od dolu przez szybko
       rosnaca czesc kubiczna (psi^3/3 dominuje dla dostatecznie duzego psi).

Nowe minimum (stacjonarnosc):
    V'(psi) = 0 dla psi > 1 przy eps*(psi-1)^3 >> 1:
    psi^2 - 1/(4*eps) = 0  =>  psi_min ~ 1/(2*sqrt(eps))

Glebokos nowego minimum:
    V_mod(psi_min) - V_mod(1) ~ -1/(12*eps^(3/2))
    Dla eps=1e-7: ~ -1/(12*3.16e-10.5) ...
    Znacznie plyciej niz V_orig(psi_min) ~ -psi_min^4/4  ~ -1/(64*eps^2)

BRAK katastrofalnego odejmowania:
    V_mod(psi) - V_mod(1) to JEDNO wyrazenie, nie suma dwoch duzych liczb.

Cel:
    - Znalezc eps* takie ze r31 = K3/K1 = 3477 (stosunek tau/elektron)
    - Sprawdzic czy r21 = K2/K1 ~ 207 (stosunek mion/elektron)
    - Potwierdzic zbieznosc calki energii (vs v2 gdzie byl artefakt)
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0

# ======================================================================
# POTENCJAL
# ======================================================================
def V_mod(psi, eps, n=3):
    """V_mod = psi^3/3 - psi^4 / (4*(1 + eps*(psi-1)^n))"""
    denom = 1.0 + eps * (psi - 1.0)**n
    return GAMMA/3 * psi**3 - GAMMA/4 * psi**4 / denom

def V_mod_vacuum(eps, n=3):
    """V_mod(1) = 1/12 niezaleznie od eps."""
    return GAMMA/3 - GAMMA/4   # = 1/12

def V_new_minimum(eps, n=3):
    """Przyblizenie polozenia i glebokosci nowego minimum dla duzych psi."""
    psi_min = 1.0 / (2.0 * np.sqrt(eps))
    V_depth = -1.0 / (12.0 * eps**1.5)
    return psi_min, V_depth

# ======================================================================
# ENERGIA SOLITONU (ansatz Yukawa)
# ======================================================================
def energy(K, Phi0, alpha, a_gam, eps, n=3, N=4000):
    msp   = np.sqrt(max(GAMMA / Phi0, 1e-9))
    r_max = max(80.0 / msp, 30.0)
    r     = np.linspace(a_gam, r_max, N)
    phi   = np.maximum(Phi0 + K * np.exp(-msp * r) / r, 1e-10)
    dphi  = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek    = 4*np.pi * np.trapezoid(
                0.5 * dphi**2 * (1.0 + alpha / (Phi0 * phi)) * r**2, r)
    psi   = phi / Phi0
    V1    = V_mod_vacuum(eps, n)
    Ep    = 4*np.pi * Phi0**2 * np.trapezoid(
                (V_mod(psi, eps, n) - V1) * r**2, r)
    return Ek + Ep

# ======================================================================
# ZNAJDOWANIE ZEROW g(M) = E(K) / K - 4*pi = 0
# ======================================================================
def find_crossings(alpha, a_gam, eps, n=3,
                   K_min=1e-4, K_max=300.0, N_scan=1000):
    K_arr = np.linspace(K_min, K_max, N_scan)
    f_arr = np.array([
        energy(K, 1.0, alpha, a_gam, eps, n) / K - 4*np.pi
        for K in K_arr])
    roots = []
    for i in range(len(f_arr) - 1):
        fi, fj = f_arr[i], f_arr[i+1]
        if np.isfinite(fi) and np.isfinite(fj) and fi * fj < 0:
            try:
                root = brentq(
                    lambda K: energy(K,1.0,alpha,a_gam,eps,n)/K - 4*np.pi,
                    K_arr[i], K_arr[i+1], xtol=1e-10)
                roots.append(root)
            except Exception:
                pass
    return sorted(roots)

# ======================================================================
# BISEKCJA po eps*: szukamy eps takie ze r31 = 3477
# ======================================================================
def find_eps_for_r31(alpha, a_gam, target_r31=3477.0, n=3,
                      eps_lo=1e-12, eps_hi=1e-5,
                      tol=1e-3):
    """
    Bisekcja: eps* takie ze K3/K1 = target_r31.
    Zwraca (eps*, K1, K2, K3) lub None jesli brak rozwiazania.
    """
    def residual(eps):
        roots = find_crossings(alpha, a_gam, eps, n)
        if len(roots) < 3:
            return len(roots) - 3   # ujemne jesli za malo zer
        K1, K2, K3 = roots[0], roots[1], roots[2]
        return K3/K1 - target_r31

    # Sprawdz czy w ogole mamy 3 zera przy eps_hi
    r_hi = find_crossings(alpha, a_gam, eps_hi, n)
    if len(r_hi) < 3:
        return None

    try:
        eps_star = brentq(residual, eps_lo, eps_hi, xtol=eps_hi*tol,
                          maxiter=60)
    except Exception:
        return None

    roots = find_crossings(alpha, a_gam, eps_star, n)
    if len(roots) < 3:
        return None
    return eps_star, roots[0], roots[1], roots[2]


# ======================================================================
# PSI_CORE i diagnoza
# ======================================================================
def psi_core(K, Phi0, a_gam, eps):
    msp = np.sqrt(max(GAMMA / Phi0, 1e-9))
    phi_core = Phi0 + K * np.exp(-msp * a_gam) / a_gam
    return phi_core / Phi0

# ======================================================================
# GLOWNA ANALIZA
# ======================================================================
print("=" * 75)
print("V_mod(psi) = psi^3/3 - psi^4 / [4*(1 + eps*(psi-1)^3)]")
print("Szukamy eps* takie ze r31 = 3477 (i sprawdzamy r21)")
print("=" * 75)
print()

# Weryfikacja warunkow prozni
eps_test = 1e-8
print(f"Weryfikacja warunkow prozni (eps={eps_test:.1e}):")
psi_arr = np.array([1.0, 1.0001])
V1   = V_mod(1.0, eps_test)
V1p  = V_mod(1.0 + 1e-5, eps_test)
V1m  = V_mod(1.0 - 1e-5, eps_test)
dV   = (V1p - V1m) / (2e-5)
d2V  = (V1p - 2*V1 + V1m) / (1e-5)**2
print(f"  V_mod(1)   = {V1:.8f}  (oczekiwane {GAMMA/3-GAMMA/4:.8f})")
print(f"  V'_mod(1)  = {dV:.6e}  (oczekiwane 0)")
print(f"  V''_mod(1) = {d2V:.6f}  (oczekiwane -1.0)")
print()

# Nowe minimum
psi_min_approx, V_depth_approx = V_new_minimum(eps_test)
print(f"Nowe minimum dla eps={eps_test:.1e}:")
print(f"  psi_min ~ 1/(2*sqrt(eps)) = {psi_min_approx:.1f}")
print(f"  V_depth ~ -1/(12*eps^1.5) = {V_depth_approx:.3e}")
print(f"  (dla porownania: V_orig(psi_min) ~ {-psi_min_approx**4/4:.3e})")
print(f"  => nowe minimum plytsze od V_orig o czynnik {abs(V_depth_approx/(-psi_min_approx**4/4)):.3f}")
print()

# ======================================================================
# SKAN: (alpha, a_gam) x eps
# ======================================================================
print("=" * 75)
print("SKAN BISEKCJI: r31 = 3477")
print(f"{'alpha':>6} {'a_gam':>6} | {'eps*':>12} {'K1':>9} {'K2':>9} {'K3':>9}"
      f" | {'r21':>8} {'r21 err':>8} | {'psi_core(M3)':>13} {'psi_min':>10}")
print("-" * 95)

cases = [
    (5.9, 0.030), (5.9, 0.050),
    (7.0, 0.025), (7.0, 0.030), (7.0, 0.035),
    (8.0, 0.030), (8.0, 0.035), (8.0, 0.040),
    (10.0, 0.040), (10.0, 0.050),
]

TARGET_R21  = 207.0
TARGET_R31  = 3477.0
results_v3  = []

for alpha, a_gam in cases:
    res = find_eps_for_r31(alpha, a_gam, target_r31=TARGET_R31, n=3)
    if res is None:
        print(f"{alpha:>6.1f} {a_gam:>6.3f} | --- brak 3 zer lub brak zbieznosci ---")
        continue
    eps_star, K1, K2, K3 = res
    r21 = K2 / K1
    r31 = K3 / K1
    err21 = 100 * (r21 - TARGET_R21) / TARGET_R21
    pc3   = psi_core(K3, 1.0, a_gam, eps_star)
    pm    = 1.0 / (2.0 * np.sqrt(eps_star))
    print(f"{alpha:>6.1f} {a_gam:>6.3f} | {eps_star:>12.4e} {K1:>9.5f} {K2:>9.5f} {K3:>9.5f}"
          f" | {r21:>8.1f} {err21:>+7.2f}% | {pc3:>13.1f} {pm:>10.1f}")
    results_v3.append({
        'alpha': alpha, 'a_gam': a_gam, 'eps': eps_star,
        'K1': K1, 'K2': K2, 'K3': K3,
        'r21': r21, 'err21': err21, 'psi_core_M3': pc3, 'psi_min': pm
    })

print()

# ======================================================================
# NAJLEPSZY WYNIK i analiza stabilnosci
# ======================================================================
if results_v3:
    best = sorted(results_v3, key=lambda x: abs(x['err21']))[0]
    print("=" * 75)
    print(f"NAJLEPSZY WYNIK: alpha={best['alpha']}, a_gam={best['a_gam']}, eps*={best['eps']:.4e}")
    print(f"  r21 = {best['r21']:.2f} (err {best['err21']:+.2f}%)")
    print(f"  r31 = {best['K3']/best['K1']:.1f}")
    print(f"  psi_core(M3) = {best['psi_core_M3']:.1f}")
    print(f"  psi_min      = {best['psi_min']:.1f}")
    ratio = best['psi_core_M3'] / best['psi_min']
    print(f"  psi_core/psi_min = {ratio:.3f}")
    if ratio < 1.0:
        print("  => M3 PRZED nowym minimum (stabilnie!)")
    elif ratio < 1.5:
        print("  => M3 blisko nowego minimum (sprawdz zbieznosc!)")
    else:
        print("  => M3 daleko za nowym minimum (mozliwy artefakt)")
    print()

    # Zbieznosc calki energii dla M3
    print("=" * 75)
    print(f"TEST ZBIEZNOSCI calki energii dla M3 (eps={best['eps']:.4e}):")
    print(f"{'N':>6} {'E_tot':>15} {'E_tot/K':>15} {'4*pi':>10}")
    K3   = best['K3']
    eps3 = best['eps']
    alpha3, a_gam3 = best['alpha'], best['a_gam']
    E_prev = None
    for N_pts in [500, 1000, 2000, 4000, 8000]:
        E = energy(K3, 1.0, alpha3, a_gam3, eps3, n=3, N=N_pts)
        ratio_4pi = E / K3 / (4*np.pi)
        delta = f"  delta={100*(E-E_prev)/abs(E_prev):+.2f}%" if E_prev is not None else ""
        print(f"{N_pts:>6} {E:>15.4e} {E/K3:>15.4f} {ratio_4pi:>10.6f}{delta}")
        E_prev = E
    print()
    print(f"  Oczekiwane: E_tot/K = 4*pi = {4*np.pi:.6f}")

    # Skladowe energii dla M3
    print()
    print("SKLADOWE ENERGII dla M3 (N=4000):")
    N_diag = 4000
    msp   = 1.0
    r_max = max(80.0 / msp, 30.0)
    r     = np.linspace(a_gam3, r_max, N_diag)
    phi   = np.maximum(1.0 + K3 * np.exp(-msp * r) / r, 1e-10)
    dphi  = K3 * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek3   = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1 + alpha3 / phi) * r**2, r)
    psi3  = phi / 1.0
    V1_   = V_mod_vacuum(eps3, 3)
    integrand = V_mod(psi3, eps3, 3) - V1_
    Ep3   = 4*np.pi * np.trapezoid(integrand * r**2, r)

    max_integrand = np.max(np.abs(integrand))
    print(f"  Ek  = {Ek3:+.6e}")
    print(f"  Ep  = {Ep3:+.6e}")
    print(f"  Etot= {Ek3+Ep3:+.6e}")
    print(f"  max|V_mod - V1| w calce = {max_integrand:.4e}")
    # Sprawdz czy jest katastrofalne odejmowanie
    # (to jest JEDNA calka, wiec nie ma osobnych Ep_orig i Ep_stab)
    print()
    print("  V_mod jest JEDNYM wyrazeniem => brak katastrofalnego odejmowania!")
    print(f"  Dla porownania v2: Ep_orig=-944e6, Ep_stab=+892e6 (ratio ~1726%)")
    print(f"  Tutaj max|integrand|/|Ep| = {max_integrand / max(abs(Ep3),1e-10):.2f}")

# ======================================================================
# ANALIZA PSI_CORE vs PSI_MIN dla wszystkich wynikow
# ======================================================================
print()
print("=" * 75)
print("DIAGNOSTYKA: psi_core(M3) vs psi_min dla wszystkich znalezionych eps*:")
print(f"  Czy M3 jest przed (ratio<1) czy za (ratio>1) nowym minimum?")
print()
for r in results_v3:
    ratio = r['psi_core_M3'] / r['psi_min']
    flag  = "OK (przed min)" if ratio < 1.0 else f"ZA minimum (x{ratio:.2f})"
    print(f"  alpha={r['alpha']:.1f}, a_gam={r['a_gam']:.3f}: "
          f"psi_core={r['psi_core_M3']:.1f}, psi_min={r['psi_min']:.1f}, "
          f"ratio={ratio:.3f} => {flag}")

print()
print("=" * 75)
print("PODSUMOWANIE PODEJSCIA v3:")
print()
print("  V_mod(psi) = psi^3/3 - psi^4 / [4*(1 + eps*(psi-1)^3)]")
print()
print("  ZALETY vs v2 (lambda*(psi-1)^6):")
print("  1. Jedno wyrazenie => brak katastrofalnego odejmowania")
print("  2. Nowe minimum o wiele plytsze (czynnik ~eps^0.5)")
print("  3. V_mod rosnie jak psi^3 dla duzych psi => stabilna proznia")
print()
print("  KLUCZOWE PYTANIE: czy psi_core(M3) < psi_min (M3 fizyczny)?")
print("  Jesli tak => M3 jest stabilny, calka zbiezna, brak artefaktu.")
