"""
v2_weryfikacja_lam_duze.py

Dla lam=5.514e-7: K3~81 jest rzeczywistym solitonem (zbiezna calka).
Pytanie: czy dla wiekszych lam (lam=2e-6, 3.88e-6) tez istnieje
         rzeczywiste 3. zero g(K), czy to artefakty kwadratury?

Test: skanujemy g(K) dla lam=2e-6 i lam=3.88e-6 z N=8000 (log grid),
      sprawdzamy czy jest zbiezna zmiana znaku.
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
ALPHA = 7.0
A_GAM = 0.030

def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_log(K, lam, N=8000):
    r_max = 60.0
    t   = np.linspace(0, 1, N)
    r   = A_GAM * (r_max/A_GAM)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r-1.0)/r**2
    psi  = phi
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

print("=" * 75)
print("g(K) scan z N=8000 (log grid) dla roznych lambda")
print("Szukamy zbieznych zmian znaku => rzeczywiste 3. zera")
print()

for lam in [5.514e-7, 2.0e-6, 3.88e-6, 1.0e-5]:
    psi_new   = 1.0/np.sqrt(lam)
    psi_cross = np.sqrt(3.0/(2.0*lam))   # V_mod(psi_cross) ~ V_mod(1)
    print(f"lambda = {lam:.4e}, psi_new={psi_new:.1f}, psi_cross={psi_cross:.1f}")
    print(f"{'K':>8} {'Ek':>12} {'Ep':>12} {'E/K':>12} {'g':>12} {'psi_c':>8}")
    print("-" * 70)

    # K zakres: od malych az do K gdzie psi_core ~ 3*psi_cross
    K_max_scan = psi_cross * A_GAM * 4  # K takie ze psi_core ~ 4*psi_cross
    K_arr = np.concatenate([
        np.linspace(0.005, 5.0, 30),
        np.linspace(5.0, K_max_scan/2, 20),
        np.linspace(K_max_scan/2, K_max_scan, 20),
    ])
    K_arr = np.unique(K_arr[K_arr > 0])

    g_arr = []
    for K in K_arr:
        Ek, Ep, E = energy_log(K, lam, N=8000)
        g   = E/K - 4*np.pi
        psi_c = 1.0 + K*np.exp(-A_GAM)/A_GAM
        g_arr.append(g)
        if abs(g) < 5e6 or K < 5.0:  # drukuj tylko ciekawe wartosci
            print(f"{K:>8.3f} {Ek:>12.4e} {Ep:>12.4e} {E/K:>12.4f} {g:>12.4e} {psi_c:>8.1f}")

    # Znajdz zmiany znaku
    print()
    sign_changes = []
    for i in range(len(K_arr)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            sign_changes.append((K_arr[i], K_arr[i+1], gi, gj))

    print(f"  Zmiany znaku g(K): {len(sign_changes)}")
    for Ka, Kb, ga, gb in sign_changes:
        Kc = Ka + (Kb-Ka)*abs(ga)/(abs(ga)+abs(gb))  # interpolacja
        psi_c = 1.0 + Kc*np.exp(-A_GAM)/A_GAM
        print(f"    K in [{Ka:.3f}, {Kb:.3f}]: ga={ga:.3e}, gb={gb:.3e}")
        print(f"    => K_zero ~ {Kc:.3f}, psi_core ~ {psi_c:.1f}")
    print()

    if len(sign_changes) >= 3:
        K1_approx = sign_changes[0][0]
        K3_approx = sign_changes[2][0]
        print(f"  WSTEPNE: K1~{K1_approx:.4f}, K3~{K3_approx:.3f}")
        print(f"  r31 ~ {K3_approx/K1_approx:.1f}")
    print()
    print("-" * 75)
    print()
