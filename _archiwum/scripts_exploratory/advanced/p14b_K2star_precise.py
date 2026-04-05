"""
p14b_K2star_precise.py
======================
Precyzyjna bisekcja K*2(ODE) - drugiego samospojnego solitonu TGP.

Z p14: istnieje zero g(K) w okolicach K~0.100, psi_core~2.76.
Cel: wyznaczyc K*2 z dokladnoscia 6 cyfr.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

ALPHA  = 8.5616
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX  = 80.0
V1     = GAMMA/3 - GAMMA/4

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=2500):
    dphi0 = -K / A_GAM**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval = A_GAM * (r_max/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def find_psi_core(K, psi_lo=1.5, psi_hi=5.0, tol=1e-6):
    """Znajdz psi_core z phi(R_max)=1 dla danego K (galaz koolo psi~2.7)."""
    F_lo = phi_at_rmax(psi_lo, K) - 1.0
    F_hi = phi_at_rmax(psi_hi, K) - 1.0
    if not (np.isfinite(F_lo) and np.isfinite(F_hi)):
        return np.nan
    if F_lo * F_hi > 0:
        # Przeszukaj szerzej
        for psi_try in np.linspace(psi_lo, psi_hi, 50):
            F_try = phi_at_rmax(psi_try, K) - 1.0
            if not np.isfinite(F_try): continue
            if F_try * F_lo < 0:
                F_hi = F_try; psi_hi = psi_try
                break
            F_lo = F_try; psi_lo = psi_try
    if F_lo * F_hi > 0: return np.nan
    return brentq(lambda p: phi_at_rmax(p, K) - 1.0,
                  psi_lo, psi_hi, xtol=tol, maxiter=60)

def compute_g(K, psi_core, n_eval=8000):
    """Oblicz g = E/(4*pi*K) - 1."""
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-10, atol=1e-12, t_eval=r_eval)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    E  = Ek + Ep
    g  = E / (4*np.pi*K) - 1.0
    return g, E, Ek, Ep

print("P14b: Precyzyjna bisekcja K*2(ODE)")
print("="*65)

# Weryfikacja punktow brzegowych z p14
print("Weryfikacja punktow z p14:")
for K_test in [0.0848, 0.090, 0.095, 0.100, 0.105]:
    psi = find_psi_core(K_test)
    if not np.isnan(psi):
        g, E, Ek, Ep = compute_g(K_test, psi)
        print(f"  K={K_test:.5f}: psi_core={psi:.5f}, g={g:.6f}")
    else:
        print(f"  K={K_test:.5f}: psi_core=NaN")
print()

# Bisekcja K*2
print("Bisekcja K*2 (g=0 na galezi psi~2.76):")
K_lo, K_hi = 0.084, 0.102

def g_track(K):
    """g(K) na ciaglaej galezi psi_core~2.7."""
    psi = find_psi_core(K, psi_lo=1.8, psi_hi=4.5)
    if np.isnan(psi):
        return np.nan
    g, _, _, _ = compute_g(K, psi)
    return g

# Fine scan do weryfikacji
print("Fine skan K in [0.084, 0.102]:")
K_fine = np.linspace(0.084, 0.102, 20)
g_fine = []
for K in K_fine:
    psi = find_psi_core(K, psi_lo=1.8, psi_hi=4.5)
    if not np.isnan(psi):
        g, E, Ek, Ep = compute_g(K, psi)
        g_fine.append(g)
        marker = ' <-- ZERO!' if abs(g) < 0.05 else ''
        print(f"  K={K:.5f}: psi={psi:.4f}, g={g:.6e}{marker}")
    else:
        g_fine.append(np.nan)
        print(f"  K={K:.5f}: psi=NaN")
print()

# Znajdz przedzial ze zmiana znaku
K_zero_found = False
for i in range(len(K_fine)-1):
    gi, gj = g_fine[i], g_fine[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        K_lo, K_hi = K_fine[i], K_fine[i+1]
        print(f"  Przedzial zmiany znaku: [{K_lo:.5f}, {K_hi:.5f}], g: [{gi:.4e}, {gj:.4e}]")
        K_zero_found = True
        break

if K_zero_found:
    print()
    print("Precyzyjna bisekcja brentq:")
    try:
        K_star2 = brentq(g_track, K_lo, K_hi, xtol=1e-6, maxiter=50,
                         full_output=False)
        psi_star2 = find_psi_core(K_star2, psi_lo=1.8, psi_hi=4.5, tol=1e-7)
        g_star2, E_star2, Ek_star2, Ep_star2 = compute_g(K_star2, psi_star2)

        print(f"  K*2(ODE)   = {K_star2:.8f}")
        print(f"  psi_core*2 = {psi_star2:.8f}")
        print(f"  g          = {g_star2:.4e}")
        print(f"  E          = {E_star2:.8f}")
        print(f"  Ek         = {Ek_star2:.8f}")
        print(f"  Ep         = {Ep_star2:.8f}")
        print(f"  Ek/|Ep|    = {Ek_star2/abs(Ep_star2):.4f}  [Derrick: = 3?]")
        print()

        # Porownanie z K*1
        K1_star = 0.010414
        psi1    = 1.2419
        g1, E1, Ek1, Ep1 = compute_g(K1_star, psi1)

        print(f"Porownanie K*1 vs K*2:")
        print(f"  K*1 = {K1_star:.8f},  K*2 = {K_star2:.8f}")
        print(f"  Stosunek K*2/K*1 = {K_star2/K1_star:.6f}")
        print(f"  (oczekiwany z Yukawa r21=207: K2/K1=207)")
        print(f"  psi_core*1 = {psi1:.4f}, psi_core*2 = {psi_star2:.4f}")
        print(f"  E*1 = {E1:.6e}, E*2 = {E_star2:.6e}")
        print(f"  Ek/|Ep|: Gen1={Ek1/abs(Ep1):.3f}, Gen2={Ek_star2/abs(Ep_star2):.3f}")
        print()
        print(f"WYNIK: K*2/K*1 = {K_star2/K1_star:.4f}")
        print(f"  -> ratio {K_star2/K1_star:.2f} (NIE 207!)")
        print()
        print(f"Dlaczego ODE daje inny stosunek niz Yukawa?")
        print(f"  Yukawa jest przyblizeniem liniowym (phi = 1 + delta_phi)")
        print(f"  Dla K*2/K*1 = {K_star2/K1_star:.1f}, delta_phi(a_Gam) ~ K*2/a_Gam = {K_star2/A_GAM:.2f}")
        print(f"  To NIE jest mala perturbacja! Nieliniowe ODE daje inna hierarchie.")
    except Exception as e:
        print(f"  Blad bisekcji: {e}")
else:
    # Moze zero jest przy K=0.100 wlasnie
    print("Sprawdzam K=0.100 bezposrednio:")
    K_test = 0.100
    psi_t = find_psi_core(K_test, psi_lo=1.8, psi_hi=4.5)
    if not np.isnan(psi_t):
        g_t, E_t, Ek_t, Ep_t = compute_g(K_test, psi_t, n_eval=10000)
        print(f"  K={K_test}: psi={psi_t:.6f}, g={g_t:.8e}")
        print(f"  -> {'SAMOSPOJNY' if abs(g_t)<1e-3 else 'niespojny'}")
