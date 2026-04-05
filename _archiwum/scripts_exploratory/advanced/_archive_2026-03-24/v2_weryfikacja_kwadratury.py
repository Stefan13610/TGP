"""
v2_weryfikacja_kwadratury.py

KLUCZOWE PYTANIE:
Czy 3. zero g(K) = E(K)/K - 4*pi znika gdy uzyjemy siatki logarytmicznej
(zamiast rownoomiern)?

Wynik spodziewany (na podstawie v2_kwadryatura.py):
  - E(K=47, log) = -2.39e8  => g < 0  (nie ma 3. zera przy K=47)
  - E(K=70, log) = -4.43e8  => g < 0  (nie ma 3. zera przy K=70)
  - Siatka logarytmiczna zbiega natychmiast (N=500 = N=8000)
  => 3. zero bylo ARTEFAKTEM rownoomiern siatki
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
LAM   = 5.514e-7

def V_mod(psi, lam=LAM):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

V1 = V_mod(1.0)

def energy_log(K, alpha, a_gam, lam=LAM, N=2000):
    """Energia na siatce logarytmicznej."""
    r_max = max(60.0, 20.0)
    t     = np.linspace(0, 1, N)
    r     = a_gam * (r_max/a_gam)**t
    phi   = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi  = K*np.exp(-r)*(-r - 1.0)/r**2
    psi   = phi  # Phi0=1
    Ek    = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep    = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek + Ep

def energy_uniform(K, alpha, a_gam, lam=LAM, N=4000):
    r_max = max(60.0, 20.0)
    r     = np.linspace(a_gam, r_max, N)
    phi   = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi  = K*np.exp(-r)*(-r - 1.0)/r**2
    psi   = phi
    Ek    = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep    = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek + Ep

ALPHA = 7.0
A_GAM = 0.030

print("=" * 70)
print("g(K) = E(K)/K - 4*pi dla alpha=7.0, a_gam=0.030")
print("Siatka rownoomierna (N=4000) vs logarytmiczna (N=2000)")
print("=" * 70)
print()
print(f"{'K':>8} {'g_uniform':>14} {'g_log':>14} {'psi_core':>10}")
print("-" * 50)
K_arr = np.concatenate([
    np.linspace(0.005, 0.1, 12),
    np.linspace(0.1, 5.0, 20),
    np.linspace(5.0, 30.0, 15),
    np.linspace(30.0, 100.0, 15),
    np.linspace(100.0, 200.0, 8),
])
for K in K_arr:
    gu = energy_uniform(K, ALPHA, A_GAM)/K - 4*np.pi
    gl = energy_log(K, ALPHA, A_GAM)/K - 4*np.pi
    pc = 1.0 + K * np.exp(-A_GAM) / A_GAM
    marker = ""
    if gu * gl < 0 and abs(K) > 1.0:
        marker = "  <<< ROZNICA ZNAKU!"
    print(f"{K:>8.3f} {gu:>14.4e} {gl:>14.4e} {pc:>10.1f}{marker}")

print()
print("=" * 70)
print("SZUKANIE ZER g(K) z siatka logarytmiczna:")
print("=" * 70)
print()

# Szukaj zer uzywajac siatki log
K_scan = np.linspace(0.001, 200.0, 2000)
g_log  = np.array([energy_log(K, ALPHA, A_GAM)/K - 4*np.pi for K in K_scan])

roots_log = []
for i in range(len(g_log)-1):
    gi, gj = g_log[i], g_log[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        try:
            root = brentq(
                lambda K: energy_log(K, ALPHA, A_GAM)/K - 4*np.pi,
                K_scan[i], K_scan[i+1], xtol=1e-9)
            roots_log.append(root)
        except Exception:
            pass

print(f"Liczba zer (siatka log): {len(roots_log)}")
for i, K_root in enumerate(roots_log):
    M_root = 4*np.pi * K_root
    print(f"  Zera {i+1}: K = {K_root:.6f}, M = {M_root:.4f}")
if len(roots_log) == 2:
    r21 = roots_log[1] / roots_log[0]
    print(f"  r21 = K2/K1 = {r21:.2f}  (cel: 207)")
    print()
    print("  => Tylko 2 zera, brak M3!")

print()
print("Szukanie zer z siatka rownoomierna (N=4000):")
K_scan2 = np.linspace(0.001, 200.0, 2000)
g_uni   = np.array([energy_uniform(K, ALPHA, A_GAM)/K - 4*np.pi for K in K_scan2])

roots_uni = []
for i in range(len(g_uni)-1):
    gi, gj = g_uni[i], g_uni[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        try:
            root = brentq(
                lambda K: energy_uniform(K, ALPHA, A_GAM)/K - 4*np.pi,
                K_scan2[i], K_scan2[i+1], xtol=1e-9)
            roots_uni.append(root)
        except Exception:
            pass

print(f"Liczba zer (rownoomierna): {len(roots_uni)}")
for i, K_root in enumerate(roots_uni):
    M_root = 4*np.pi * K_root
    print(f"  Zero {i+1}: K = {K_root:.6f}, M = {M_root:.4f}")
if len(roots_uni) >= 2:
    r21 = roots_uni[1]/roots_uni[0]
    print(f"  r21 = {r21:.2f}")
if len(roots_uni) >= 3:
    r31 = roots_uni[2]/roots_uni[0]
    print(f"  r31 = {r31:.2f}  (<<< to jest ARTEFAKT!)")

print()
print("=" * 70)
print("WNIOSEK:")
print()
if len(roots_log) == 2 and len(roots_uni) >= 3:
    print("  Siatka rownoomierna: FALISZYWY 3. korzenia (artefakt kwadratury)")
    print("  Siatka logarytmiczna: poprawnie 2 korzenie")
    print()
    print("  Potencjal V_mod = psi^3/3 - psi^4/4 + lambda*(psi-1)^6/6")
    print("  NIE produkuje 3. solitonu dla zadnego skoncznego K.")
    print("  M3 z v2_fine_scan byl artefaktem rownoomiern siatki!")
elif len(roots_log) == 3:
    print("  Siatka logarytmiczna tez pokazuje 3 korzenie!")
    print("  M3 JEST prawdziwym solitonem (nie artefaktem kwadratury)")
    print("  Ale czy calka zbiega z N?")
else:
    print(f"  Nieoczekiwany wynik: log={len(roots_log)}, uniform={len(roots_uni)}")
