"""Diagnostyka g_ODE(A) — szybki skan."""
import numpy as np
from scipy.integrate import solve_ivp
import warnings; warnings.filterwarnings('ignore')

ALPHA = 8.5616; A_GAM = 0.040; LAM = 5.501357e-06; GAMMA = 1.0
M_EFF = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX = 10.0 / M_EFF  # ~30.9

def V_mod(phi):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + LAM/6*(phi-1)**6

def dV_mod(phi):
    return GAMMA*phi**2 - GAMMA*phi**3 + LAM*(phi-1)**5

V1 = V_mod(1.0)

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]

def shoot_and_g(A):
    eR = np.exp(-M_EFF * R_MAX)
    phi_R  = 1.0 + A * eR / R_MAX
    dphi_R = A * eR * (-1.0/R_MAX**2 - M_EFF/R_MAX)
    sol = solve_ivp(
        ode_rhs, [R_MAX, A_GAM], [phi_R, dphi_R],
        method='DOP853', rtol=1e-8, atol=1e-10,
        max_step=0.05
    )
    r    = np.sort(sol.t)
    idx  = np.argsort(sol.t)
    phi  = sol.y[0][idx]
    dphi = sol.y[1][idx]
    if phi.min() < 1e-5:
        return None, None, None, 'phi<0'
    K = r[0]**2 * abs(dphi[0])
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    E  = Ek + Ep
    if K < 1e-25:
        return None, None, None, 'K=0'
    g = E / (4*np.pi*K) - 1.0
    return K, E, g, 'ok'

print(f"R_max={R_MAX:.2f},  m_eff={M_EFF:.4f}")
print(f"Perturbacja A=1 na R_max: delta_phi = {np.exp(-M_EFF*R_MAX)/R_MAX:.4e}")
print()
print(f"{'A':>10}  {'K_eff':>12}  {'E':>12}  {'g_ODE':>12}  {'status':>8}")
print("-"*58)

A_values = np.logspace(-3, 2, 30)
g_vals = []
for A in A_values:
    K, E, g, status = shoot_and_g(A)
    if g is not None:
        print(f"{A:>10.3e}  {K:>12.4e}  {E:>12.4e}  {g:>12.4e}  {status}")
        g_vals.append(g)
    else:
        print(f"{A:>10.3e}  {'---':>12}  {'---':>12}  {'---':>12}  {status}")
        g_vals.append(np.nan)

g_vals = np.array(g_vals)
print()
finite = np.isfinite(g_vals)
if finite.any():
    print(f"g_min = {g_vals[finite].min():.4e}")
    print(f"g_max = {g_vals[finite].max():.4e}")
    sign_changes = 0
    for i in range(len(g_vals)-1):
        if np.isfinite(g_vals[i]) and np.isfinite(g_vals[i+1]):
            if g_vals[i]*g_vals[i+1] < 0:
                sign_changes += 1
                print(f"Zmiana znaku: A in [{A_values[i]:.3e}, {A_values[i+1]:.3e}]")
    print(f"Laczna liczba zmian znaku: {sign_changes}")
print()
print("UWAGA: Jezeli g_ODE > 0 wszsedzie => brak samospojnych solitonow w pelnym ODE")
print("Jezeli g_ODE < 0 wszsedzie => solitony istnieja TYLKO w Yukawa-aproksymacji")
