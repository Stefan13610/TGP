"""
p2_lambda_mechanizm.py

Problem P2: Wyznaczenie lambda* z pierwszych zasad.

Rodzina rozwiazań daje psi_core/psi_new ~ 1.93-1.97 dla wszystkich punktow.
Pytanie: czy istnieje DOKLADNY warunek wyznaczajacy lambda*?

Testowane hipotezy:
  H1: E_core(K3) + E_shell(K3) = 0 przy r = r_cross
      (rownanie energii w punkcie crossover)
  H2: psi_core = c * psi_cross
      (rdzeń lezy przy stalym wielokrotnosci psi_cross)
  H3: d/dlam [g(K3; lam)] = 0 (lambda minimalizuje |g|)
  H4: Warunek Koide: r31 = 3481 (dokladne 2/3) vs r31=3477 (TGP)
  H5: Warunek energetyczny: Ek(K3) = Ep(K3) (rownomierny podzial)
  H6: psi_core = sqrt(2) * psi_new (konkretny wspolczynnik)

Dla kazdej hipotezy: oblicz predykowane lambda*, porownaj z numerycznym lambda*.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
ALPHA = 8.5616
A_GAM = 0.040

def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def dV_mod(psi, lam):
    return GAMMA*psi**2 - GAMMA*psi**3 + lam*(psi-1)**5

def energy_components(K, lam, N=4000):
    r_max = 60.0
    t   = np.linspace(0, 1, N)
    r   = A_GAM * (r_max/A_GAM)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    psi  = phi
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)

    # Komponent rdzenia (r < r_cross) vs powloki (r >= r_cross)
    # r_cross: phi(r_cross) = psi_cross
    psi_cross = np.sqrt(3.0/(2.0*lam))
    # Znajdz r_cross numerycznie
    phi_vals = 1.0 + K*np.exp(-r)/r
    r_cross_idx = np.searchsorted(-phi_vals, -psi_cross)
    r_cross_idx = min(r_cross_idx, len(r)-2)
    r_cross = r[r_cross_idx] if r_cross_idx > 0 else A_GAM

    mask_core  = r <= r_cross
    mask_shell = r >  r_cross

    Ep_core  = 4*np.pi*np.trapezoid(
        (V_mod(psi[mask_core],lam)-V1)*r[mask_core]**2, r[mask_core]) if mask_core.sum() > 1 else 0.0
    Ep_shell = 4*np.pi*np.trapezoid(
        (V_mod(psi[mask_shell],lam)-V1)*r[mask_shell]**2, r[mask_shell]) if mask_shell.sum() > 1 else 0.0

    return {'Ek': Ek, 'Ep': Ep, 'E': Ek+Ep,
            'Ep_core': Ep_core, 'Ep_shell': Ep_shell,
            'psi_cross': psi_cross, 'r_cross': r_cross}

def find_K3_for_lam(lam, alpha=ALPHA, a_gam=A_GAM, N=2000):
    """Znajdz K3 (3. zero g) dla danego lam."""
    r_max = 60.0
    V1 = V_mod(1.0, lam)
    def g_func(K):
        if K <= 0: return -4*np.pi
        t   = np.linspace(0, 1, N)
        r   = a_gam * (r_max/a_gam)**t
        phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
        dphi = K*np.exp(-r)*(-r - 1.0)/r**2
        psi  = phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
        return (Ek+Ep)/K - 4*np.pi

    # Skan K
    K_arr = np.concatenate([
        np.linspace(0.001, 0.5, 100),
        np.linspace(0.5, 10.0, 150),
        np.linspace(10.0, 100.0, 250)])
    g_arr = [g_func(K) for K in K_arr]

    roots = []
    for i in range(len(g_arr)-1):
        if np.isfinite(g_arr[i]) and np.isfinite(g_arr[i+1]) and g_arr[i]*g_arr[i+1] < 0:
            try:
                root = brentq(g_func, K_arr[i], K_arr[i+1], xtol=1e-7)
                roots.append(root)
            except Exception:
                pass

    if len(roots) >= 3:
        return roots[0], roots[1], roots[2]
    elif len(roots) == 2:
        return roots[0], roots[1], None
    else:
        return None, None, None

# ================================================================
print("=" * 75)
print("P2: ANALIZA MECHANIZMU WYZNACZENIA lambda*")
print(f"Parametry: alpha={ALPHA}, a_gam={A_GAM}")
print(f"Znane rozwiazanie: lam*=5.501e-6, K3=34.144")
print()

LAM_REF = 5.501357e-6
K3_REF  = 34.14450
K1_REF  = 0.009820
K2_REF  = 2.032728
psi_new_ref = 1.0/np.sqrt(LAM_REF)
psi_core_ref = 1.0 + K3_REF*np.exp(-A_GAM)/A_GAM
psi_cross_ref = np.sqrt(3.0/(2.0*LAM_REF))

print(f"psi_core = {psi_core_ref:.2f}")
print(f"psi_new  = {psi_new_ref:.2f}")
print(f"psi_cross = {psi_cross_ref:.2f}")
print(f"psi_core/psi_new   = {psi_core_ref/psi_new_ref:.4f}")
print(f"psi_core/psi_cross = {psi_core_ref/psi_cross_ref:.4f}")
print(f"psi_cross/psi_new  = {psi_cross_ref/psi_new_ref:.4f} (powinna = sqrt(3/2) = {np.sqrt(3/2):.4f})")
print()

# ================================================================
# H1: Rownanie energetyczne: E_core + E_shell = 0 w punkcie crossover
# ================================================================
print("=" * 75)
print("H1: Rownosc wkladow Ep_core i Ep_shell")
comps = energy_components(K3_REF, LAM_REF, N=6000)
print(f"  Ep_core  = {comps['Ep_core']:.4e}")
print(f"  Ep_shell = {comps['Ep_shell']:.4e}")
print(f"  Ep_core + Ep_shell = {comps['Ep_core']+comps['Ep_shell']:.4e}")
print(f"  Stosunek |Ep_core/Ep_shell| = {abs(comps['Ep_core']/comps['Ep_shell']):.4f}")
print(f"  r_cross = {comps['r_cross']:.5f}")
print()

# Czy istnieje lam takie ze Ep_core = |Ep_shell|?
print("  Skan lam: stosunek |Ep_core|/|Ep_shell| (przy K3 z bledem wyzn.):")
for lam_test in [3e-6, 4e-6, 5e-6, 5.5e-6, 6e-6, 7e-6]:
    K1t, K2t, K3t = find_K3_for_lam(lam_test)
    if K3t is None:
        print(f"  lam={lam_test:.2e}: brak K3")
        continue
    c = energy_components(K3t, lam_test, N=3000)
    ratio = abs(c['Ep_core'])/max(abs(c['Ep_shell']), 1e-10)
    print(f"  lam={lam_test:.2e}: K3={K3t:.3f}, |Ep_core/Ep_shell|={ratio:.4f}")
print()

# ================================================================
# H2: psi_core = c * psi_cross (jaki c?)
# ================================================================
print("=" * 75)
print("H2: Relacja psi_core/psi_cross = const dla calyej rodziny")
print()

rodzina = [
    (5.9148, 0.025, 2.883172e-06, 29.6702),
    (6.8675, 0.030, 3.742489e-06, 31.1596),
    (7.7449, 0.035, 4.618392e-06, 32.6553),
    (8.5616, 0.040, 5.501357e-06, 34.1440),
]

print(f"{'alpha':>8} {'a_gam':>7} {'lam*':>12} {'K3':>8} "
      f"{'psi_c':>8} {'psi_n':>8} {'psi_x':>8} "
      f"{'c/n':>8} {'c/x':>8} {'n/x':>8}")
print("-" * 88)
for alpha, a_gam, lam, K3 in rodzina:
    psi_c = 1.0 + K3*np.exp(-a_gam)/a_gam
    psi_n = 1.0/np.sqrt(lam)
    psi_x = np.sqrt(3.0/(2.0*lam))
    print(f"{alpha:>8.4f} {a_gam:>7.3f} {lam:>12.4e} {K3:>8.4f} "
          f"{psi_c:>8.1f} {psi_n:>8.1f} {psi_x:>8.1f} "
          f"{psi_c/psi_n:>8.4f} {psi_c/psi_x:>8.4f} {psi_n/psi_x:>8.4f}")
print()

# ================================================================
# H3: Warunek energetyczny Ek = |Ep| (rownomierny podzial)
# ================================================================
print("=" * 75)
print("H3: Warunek Ek = |Ep| (ekwipartycja Ek i |Ep|)")
print()
for alpha, a_gam, lam, K3 in rodzina:
    def E_for_alpha_agam_lam(K, N=3000):
        r_max = 60.0
        t   = np.linspace(0, 1, N)
        r   = a_gam * (r_max/a_gam)**t
        phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
        dphi = K*np.exp(-r)*(-r - 1.0)/r**2
        V1 = V_mod(1.0, lam)
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi,lam)-V1)*r**2, r)
        return Ek, Ep
    Ek, Ep = E_for_alpha_agam_lam(K3)
    print(f"  alpha={alpha:.4f}, a_gam={a_gam:.3f}: Ek={Ek:.4e}, Ep={Ep:.4e}, Ek/|Ep|={Ek/abs(Ep):.4f}")
print()

# ================================================================
# H4: Koide: r31 = 3481 vs 3477
# ================================================================
print("=" * 75)
print("H4: Warunek Koide (r31=3481) vs TGP (r31=3477)")
print()
# Zbadaj jak zmienia sie lam* gdy r31=3481 zamiast 3477
from scipy.optimize import brentq

def find_lam_for_r31(target_r31, alpha=8.5616, a_gam=0.040):
    def residual(lam):
        K1t, K2t, K3t = find_K3_for_lam(lam, alpha, a_gam)
        if K3t is None: return -target_r31
        return K3t/K1t - target_r31
    try:
        lam_s = brentq(residual, 4e-6, 7e-6, xtol=1e-3*6e-6, maxiter=30)
        K1t, K2t, K3t = find_K3_for_lam(lam_s, alpha, a_gam)
        return lam_s, K1t, K2t, K3t
    except Exception:
        return None, None, None, None

print("  Dla alpha=8.5616, a_gam=0.040:")
for r31_target, label in [(3477.0, "TGP (r31=3477)"), (3481.0, "Koide (r31=3481)")]:
    lam_s, K1t, K2t, K3t = find_lam_for_r31(r31_target)
    if lam_s:
        psi_c = 1.0 + K3t*np.exp(-A_GAM)/A_GAM
        psi_n = 1.0/np.sqrt(lam_s)
        print(f"  {label}: lam*={lam_s:.6e}, K3={K3t:.4f}, "
              f"r21={K2t/K1t:.3f}, psi_c/psi_n={psi_c/psi_n:.4f}")
print()

# ================================================================
# H5: Warunek psi_core = sqrt(2) * psi_new (wspolczynnik sqrt(2))
# ================================================================
print("=" * 75)
print("H5: Hipoteza psi_core = c_exact * psi_new")
print("    Szukamy c_exact takie ze lam*(c_exact) == lam*_numeryczne")
print()

# Dla danego (alpha, a_gam, K3): z psi_core = K3*exp(-a_gam)/a_gam ~ K3/a_gam,
# psi_new = 1/sqrt(lam) => lam = (a_gam/(c*K3))^2 (przyblizenie dla K3>>1)
# Sprawdzamy czy to jest dokladne
for alpha, a_gam, lam, K3 in rodzina:
    psi_c = 1.0 + K3*np.exp(-a_gam)/a_gam
    psi_n = 1.0/np.sqrt(lam)
    c_obs = psi_c / psi_n
    # Predykcja lam z c=2 i c=sqrt(3):
    for c_hyp, c_name in [(2.0, "2"), (np.sqrt(3), "sqrt(3)"),
                           (np.sqrt(2*np.pi/np.e), "sqrt(2pi/e)"),
                           (np.e*np.sqrt(1/2), "e/sqrt(2)")]:
        lam_pred = (psi_c / c_hyp)**(-2)
        err = 100*(lam_pred - lam)/lam
        print(f"  a_gam={a_gam:.3f}: c_obs={c_obs:.4f}, c={c_name}: "
              f"lam_pred={lam_pred:.4e} (blad {err:+.2f}%)")
print()

# ================================================================
# H6: Crossover geometria: V_mod(psi_core) = n * V_mod(psi_cross)
# ================================================================
print("=" * 75)
print("H6: Relacja V_mod(psi_core) / V_mod(psi_cross)")
print()
for alpha, a_gam, lam, K3 in rodzina:
    psi_c = 1.0 + K3*np.exp(-a_gam)/a_gam
    psi_x = np.sqrt(3.0/(2.0*lam))
    psi_n = 1.0/np.sqrt(lam)
    Vc = V_mod(psi_c, lam)
    Vx = V_mod(psi_x, lam)
    V1v = V_mod(1.0, lam)
    print(f"  a_gam={a_gam:.3f}: V_mod(psi_core)={Vc:.4e}, "
          f"V_mod(psi_cross)={Vx:.4e}, "
          f"ratio={(Vc-V1v)/(Vx-V1v):.4f}")
print()

# ================================================================
# PODSUMOWANIE
# ================================================================
print("=" * 75)
print("PODSUMOWANIE hipotez:")
print()
print("  Kluczowe obserwacje:")
print(f"  psi_core/psi_new     = 1.926 - 1.967 (nie stale!)")
print(f"  psi_core/psi_cross   = ? (sprawdz tabele)")
print(f"  Koide: r31=3481 vs TGP r31=3477 (roznica 0.11%)")
print()
print("  Hipotezy:")
print("  H1: Ep_core = -Ep_shell  => sprawdz stosunek")
print("  H2: psi_core/psi_cross = const  => sprawdz powyzej")
print("  H5: psi_core = c * psi_new dla dokladnego c  => sprawdz bledy")
print()
print("  Wniosek: jesli ZADNA z powyzszych nie daje bledow < 0.1%,")
print("  to lambda* jest naprawde wolnym parametrem (brak zasady)")
