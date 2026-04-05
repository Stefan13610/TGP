"""
v2_konwergencja_finalna.py

Test zbieznosci energii dla kandydatow na 3. zero (siatka log)
dla lam=5.514e-7, alpha=7.0, a_gam=0.030.

Z v2_weryfikacja_kwadratury.py: zmiana znaku g(K) nastepuje miedzy K=80 i K=85.
Sprawdzamy czy ta zmiana znaku jest stabilna z N.
"""
import numpy as np
from scipy.integrate import quad
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
LAM   = 5.514e-7
ALPHA = 7.0
A_GAM = 0.030

def V_mod(psi):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + LAM/6*(psi-1)**6
V1 = V_mod(1.0)

def energy_log(K, N):
    r_max = 60.0
    t   = np.linspace(0, 1, N)
    r   = A_GAM * (r_max/A_GAM)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    psi  = phi
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(psi)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

def energy_quad_adaptive(K):
    """Adaptywna kwadryatura (nie korzysta z dyskretnej siatki)."""
    def Ep_int(r):
        phi = max(1.0 + K*np.exp(-r)/r, 1e-10)
        psi = phi
        return 4*np.pi*(V_mod(psi)-V1)*r**2

    def Ek_int(r):
        phi = max(1.0 + K*np.exp(-r)/r, 1e-10)
        dphi = K*np.exp(-r)*(-r-1.0)/r**2
        return 4*np.pi*0.5*dphi**2*(1+ALPHA/phi)*r**2

    # Podziel przedzial: [a_gam, r_cross, 10, 60]
    # r_cross ~ K/psi_cross, psi_cross = sqrt(3/(2*lam)) ~ 1650
    psi_cross = np.sqrt(3/(2*LAM))
    r_cross   = min(K/psi_cross, 60.0) if K > A_GAM*psi_cross else A_GAM

    Ek = 0.0; Ep = 0.0
    # Odcinki
    breakpoints = sorted(set([A_GAM, max(A_GAM*1.001, r_cross*0.5),
                              r_cross, min(r_cross*2, 60.0), 60.0]))
    for i in range(len(breakpoints)-1):
        a, b = breakpoints[i], breakpoints[i+1]
        if a >= b: continue
        ek_seg, _ = quad(Ek_int, a, b, limit=200, epsabs=1e-2, epsrel=1e-6)
        ep_seg, _ = quad(Ep_int, a, b, limit=200, epsabs=1e-2, epsrel=1e-6)
        Ek += ek_seg; Ep += ep_seg
    return Ek, Ep, Ek+Ep

# ── Test zbieznosci przy kandydatach na 3. zero ────────────────────────────
K_candidates = [75, 78, 80, 82, 85, 90, 100]

print("=" * 70)
print(f"ZBIEZNOSC g(K) = E(K)/K - 4*pi z N (siatka log, lam=5.514e-7)")
print(f"{'K':>6} {'N':>6} {'Ek':>12} {'Ep':>12} {'E/K':>12} {'g=E/K-4pi':>12}")
print("-" * 65)

for K in K_candidates:
    psi_c = 1.0 + K * np.exp(-A_GAM)/A_GAM
    print(f"K={K}, psi_core={psi_c:.0f}:")
    E_prev = None
    for N in [1000, 2000, 4000, 8000, 16000]:
        Ek, Ep, E = energy_log(K, N)
        g = E/K - 4*np.pi
        delta = f" dlt={100*(E-E_prev)/abs(E_prev+1e-10):+.2f}%" if E_prev is not None else ""
        print(f"  {N:>6} {Ek:>12.4e} {Ep:>12.4e} {E/K:>12.4f} {g:>12.4e}{delta}")
        E_prev = E

    # Adaptywna kwadryatura
    try:
        Ek, Ep, E = energy_quad_adaptive(K)
        g = E/K - 4*np.pi
        print(f"  {'quad':>6} {Ek:>12.4e} {Ep:>12.4e} {E/K:>12.4f} {g:>12.4e}")
    except Exception as ex:
        print(f"  quad error: {ex}")
    print()

print()
print(f"Cel: E/K = 4*pi = {4*np.pi:.6f}")
print()
print("INTERPRETACJA:")
print("  Jesli g(K) zbiega do tego samego znaku dla wszystkich N")
print("  => brak 3. zera (lub zero rzeczywiste)")
print("  Jesli g(K) zmienia znak z N => artefakt kwadratury")
