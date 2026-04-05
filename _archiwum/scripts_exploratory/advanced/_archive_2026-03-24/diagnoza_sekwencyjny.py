"""
Diagnostyka sekwencyjnego modelu — wewnetrzna sprzecznosc
=========================================================

Pytanie: czy model sekwencyjny (Phi0_n rosnie o xi*M_{n-1}) moze
         jednoczesnie dac M2/M1=207 i M3/M1=3477?
"""
import numpy as np
from scipy.optimize import brentq

GAMMA = 1.0; PHI0_BASE = 1.0

def energy(K, Phi0, alpha, a_gam):
    msp = np.sqrt(max(GAMMA * PHI0_BASE / Phi0, 1e-6))
    r_max = max(30.0 / msp, 10.0)
    r = np.linspace(a_gam, r_max, 2000)
    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1.0 + alpha / (Phi0 * phi)) * r**2, r)
    psi  = phi / Phi0
    Ep   = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4) * r**2, r)
    return Ek + Ep

def g(M, Phi0, alpha, a_gam):
    if M < 1e-10: return 0.0
    K = M / (4.0 * np.pi * Phi0)
    return energy(K, Phi0, alpha, a_gam) - M

def first_root(Phi0, alpha, a_gam, M_max=30000):
    M_arr = np.concatenate([np.linspace(1e-4, 2, 120), np.linspace(2, 100, 100),
                            np.linspace(100, 2000, 80), np.linspace(2000, M_max, 60)])
    M_arr = np.unique(M_arr)
    g_arr = np.array([g(M, Phi0, alpha, a_gam) for M in M_arr])
    for i in range(len(g_arr)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                return brentq(lambda M: g(M, Phi0, alpha, a_gam),
                              M_arr[i], M_arr[i+1], xtol=1e-6)
            except Exception:
                pass
    return None

alpha = 5.9; a_gam = 0.05
M1 = first_root(1.0, alpha, a_gam)

print('='*65)
print('DIAGNOSTYKA SEKWENCYJNEGO MODELU (alpha=5.9, a_Gam=0.05)')
print('='*65)
print(f'M1 = {M1:.6f} (w Phi0_base=1)')
print()

print('SKALOWANIE M(Phi0):')
print(f"{'Phi0':>8}  {'M':>10}  {'M/M1':>10}  {'wykladnik':>12}")
for Phi0 in [1.0, 2.0, 5.0, 6.3, 10.0, 20.0, 26.6, 50.0]:
    Mn = first_root(Phi0, alpha, a_gam)
    if Mn:
        exp = np.log(Mn/M1)/np.log(Phi0) if Phi0 > 1 else 0.0
        print(f'{Phi0:>8.1f}  {Mn:>10.4f}  {Mn/M1:>10.2f}  Phi0^{exp:.3f}')

print()
print('='*65)
print('WARUNKI NA Phi0 DLA CELOW M2/M1=207, M3/M1=3477:')
print('='*65)
# Znajdz Phi0 ktore daje M/M1 = 207 i 3477 przez interpolacje
targets = [207.0, 3477.0]
Phi0_needed = []
for target in targets:
    # Binary search on Phi0
    lo, hi = 1.0, 200.0
    for _ in range(40):
        mid = (lo+hi)/2
        Mm = first_root(mid, alpha, a_gam)
        if Mm and Mm/M1 < target:
            lo = mid
        else:
            hi = mid
    Phi0_needed.append((lo+hi)/2)

Phi0_2_need, Phi0_3_need = Phi0_needed
print(f'Phi0_2 potrzebne dla M2/M1=207:  Phi0_2 = {Phi0_2_need:.3f}')
print(f'Phi0_3 potrzebne dla M3/M1=3477: Phi0_3 = {Phi0_3_need:.3f}')
print()

xi_needed = (Phi0_2_need - PHI0_BASE) / M1
M2_at = first_root(Phi0_2_need, alpha, a_gam)
print(f'xi z warunku Phi0_2: xi = (Phi0_2 - 1)/M1 = {xi_needed:.3f}')
print(f'M2 (przy Phi0_2={Phi0_2_need:.2f}) = {M2_at:.4f} = {M2_at/M1:.1f}*M1')
print()

Phi0_3_sequential = PHI0_BASE + xi_needed * (M1 + M2_at)
print(f'Phi0_3 SEKWENCYJNE = 1 + xi*(M1+M2) = 1 + {xi_needed:.2f}*{M1+M2_at:.2f} = {Phi0_3_sequential:.1f}')
print(f'Phi0_3 POTRZEBNE   = {Phi0_3_need:.1f}')
print()
print(f'Roznica: sekwencyjny daje Phi0_3 = {Phi0_3_sequential/Phi0_3_need:.0f}x za duze!')
print()
print('PRZYCZYNA:')
print(f'  Phi0_3/Phi0_2 sekwencyjne = 1 + M2/M1 = 1 + {M2_at/M1:.0f} = {1+M2_at/M1:.0f}')
print(f'  Phi0_3/Phi0_2 potrzebne   = {Phi0_3_need/Phi0_2_need:.2f}')
print()
print('=> Model sekwencyjny NIEMOZLIWY dla tych stosunkow mas.')
print('   Phi0 rosnie za szybko miedzy generacjami.')
print()
print('='*65)
print('ALTERNATYWA: Phi0 jako ZEWNETRZNY parametr kosmologiczny')
print('='*65)
print()
print('Jezeli Phi0_n nie pochodzi z masy poprzedniej generacji,')
print('ale z GLOBALNEJ gestosci kosmicznej (rho_cosmic * xi_cosmo),')
print('to trzy wartosci Phi0 sa NIEZALEZNYMI parametrami TGP.')
print()
print(f'  Phi0_1 = {1.0:.3f}  =>  M1 = {M1:.4f}')
M2_check = first_root(Phi0_2_need, alpha, a_gam)
M3_check = first_root(Phi0_3_need, alpha, a_gam)
print(f'  Phi0_2 = {Phi0_2_need:.3f}  =>  M2 = {M2_check:.4f}  (M2/M1 = {M2_check/M1:.1f})')
print(f'  Phi0_3 = {Phi0_3_need:.3f}  =>  M3 = {M3_check:.4f}  (M3/M1 = {M3_check/M1:.1f})')
print()
print(f'  Phi0_2/Phi0_1 = {Phi0_2_need:.3f}')
print(f'  Phi0_3/Phi0_2 = {Phi0_3_need/Phi0_2_need:.3f}')
print(f'  Phi0_3/Phi0_1 = {Phi0_3_need:.3f}')
print()
print('Pytanie fizyczne: co wyznacza te trzy wartosci Phi0?')
print('  -> Gestosci kosmiczne w trzech epokach?')
print('  -> Trzy stabilne fazy pola Phi?')
print('  -> Trzy rozne topologiczne sektory?')
