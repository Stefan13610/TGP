#!/usr/bin/env python3
"""Quick simplified K-test — no shooting, fixed V_0."""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.integrate import solve_ivp

OMEGA_M0  = 0.315
OMEGA_R0  = 9.1e-5
OMEGA_DE0 = 1.0 - OMEGA_M0 - OMEGA_R0
PHI_0_TGP = 24.65
S_NATURAL = 3.0 / (2.0 * PHI_0_TGP**2)

def V(psi):
    return 4.0*psi**3 - 3.0*psi**4
def dVdpsi(psi):
    return 12.0 * psi**2 * (1.0 - psi)
def Omega_m_a(a):
    return OMEGA_M0 / a**3

V_0 = OMEGA_DE0

def eom_N(N, y, s_c):
    psi, u = y
    a = np.exp(N)
    rho_total = OMEGA_R0/a**4 + OMEGA_M0/a**3 + 0.5*u**2 + V_0*V(psi)
    H = np.sqrt(max(rho_total, 1e-30))
    if H < 1e-30:
        return [0, 0]
    dpsi_dN = u / H
    du_dN = (-3.0*H*u - V_0*dVdpsi(psi) - s_c*Omega_m_a(a)) / H
    return [dpsi_dN, du_dN]

a_init = 1e-5
N_init = np.log(a_init)
N_eval = np.linspace(N_init, 0.0, 5000)

print('K-test scenarios (no shooting, fixed V_0 = Omega_DE):')
print(f'  S_natural (TGP) = {S_NATURAL:.4e}')
print(f'  Required Hubble shift |dH/H| = 0.0837')
print()
print(f'  {"label":<15} {"psi_today":>12} {"psi_recomb":>12} {"dpsi_cosmo":>14} {"dL/L":>14} {"|dH/H|":>10}')
print(f'  {"-"*15} {"-"*12} {"-"*12} {"-"*14} {"-"*14} {"-"*10}')

for label, s in [('s=0', 0.0), ('s=natural', S_NATURAL), ('s=0.01', 0.01),
                 ('s=0.1', 0.1), ('s=0.3', 0.3), ('s=1.0', 1.0), ('s=10', 10.0)]:
    sol = solve_ivp(eom_N, [N_init, 0.0], [1.0 - 1e-3, 0.0],
                    args=(s,), t_eval=N_eval, method='RK45',
                    rtol=1e-8, atol=1e-10)
    if not sol.success:
        print(f'  {label:<15} FAILED: {sol.message}')
        continue
    psi = sol.y[0]
    z = 1.0/np.exp(N_eval) - 1.0
    psi_today = psi[-1]
    idx_recomb = np.argmin(np.abs(z - 1100.0))
    psi_recomb = psi[idx_recomb]
    Vt = V(psi_today)
    Vr = V(psi_recomb)
    dLL = (Vt - Vr) / Vt if Vt != 0 else 0
    dH = 0.5 * abs(dLL) * OMEGA_DE0
    print(f'  {label:<15} {psi_today:>12.6f} {psi_recomb:>12.6f} {psi_today-psi_recomb:>+14.3e} {dLL:>+14.3e} {dH:>10.3e}')

print()
print('Done.')
