#!/usr/bin/env python3
"""ERG (Wetterich) LPA -- minimal version with stiff solver."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Phi0 = 24.66
gamma_TGP = 1.0
a_Gamma = 0.040
k_UV = 1.0 / a_Gamma
k_IR = np.sqrt(gamma_TGP)
t_max = np.log(k_UV / k_IR)

print("=" * 65)
print("  ERG (Wetterich) LPA -- STIFF SOLVER")
print("=" * 65)
print(f"  k_UV={k_UV:.1f}, k_IR={k_IR:.1f}, t_max={t_max:.2f}")

N_grid = 40
psi_min, psi_max = 0.05, 2.5
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

def V_TGP(psi):
    return psi**3 / 3 - psi**4 / 4

def Vpp_num(V, dp):
    r = np.zeros_like(V)
    r[1:-1] = (V[2:] - 2*V[1:-1] + V[:-2]) / dp**2
    r[0] = r[1]; r[-1] = r[-2]
    return r

def rhs(t, V):
    k = k_IR * np.exp(t)
    vpp = Vpp_num(V, dpsi)
    denom = np.maximum(vpp + k**2, 0.01 * k**2)  # stronger floor
    return k**5 / (32 * np.pi**2) / denom

V0 = V_TGP(psi_grid)
print(f"  Grid: {N_grid} pts, dpsi={dpsi:.4f}")
print(f"  Integrating UV->IR (Radau)...", flush=True)

sol = solve_ivp(rhs, (t_max, 0.05), V0,
                method='Radau', rtol=1e-4, atol=1e-6,
                max_step=0.5, dense_output=True)

print(f"  Status: {'OK' if sol.success else 'FAIL: '+sol.message}")
print(f"  Steps: {sol.t.shape[0]}")

if sol.success:
    V_IR = sol.y[:, -1]

    # Analysis
    idx_vac = np.argmin(np.abs(psi_grid - 1.0))
    idx_min = np.argmin(V_IR)

    print(f"\n--- Results ---")
    print(f"  Global min V_IR: psi={psi_grid[idx_min]:.3f}, V={V_IR[idx_min]:.6f}")

    vpp_ir = Vpp_num(V_IR, dpsi)
    print(f"  V''(psi_min) = {vpp_ir[idx_min]:.6f}")

    print(f"  At psi=1: V_UV={V0[idx_vac]:.6f}, V_IR={V_IR[idx_vac]:.6f}")
    print(f"    Delta = {V_IR[idx_vac]-V0[idx_vac]:.6f}")

    # Check structure preservation
    # Look for minimum near psi=1
    near = (psi_grid > 0.5) & (psi_grid < 1.5)
    if np.any(near):
        local_idx = np.argmin(V_IR[near])
        psi_loc = psi_grid[near][local_idx]
        print(f"  Local min near psi=1: psi={psi_loc:.3f}")

    # Dimensionless potential flow
    t_pts = np.linspace(t_max, 0.05, 20)
    print(f"\n  Dimensionless u=V/k^4 at psi=1:")
    for t in [t_max, t_max*0.75, t_max*0.5, t_max*0.25, 0.05]:
        k = k_IR * np.exp(t)
        V_t = sol.sol(t)
        u = V_t[idx_vac] / k**4
        print(f"    t={t:.2f} k={k:.1f}: u={u:.2e}")

    # V'' flow at psi=1
    print(f"\n  V''(psi=1) flow:")
    vpp_uv = Vpp_num(V0, dpsi)[idx_vac]
    vpp_ir_val = Vpp_num(V_IR, dpsi)[idx_vac]
    print(f"    UV: {vpp_uv:.4f}")
    print(f"    IR: {vpp_ir_val:.4f}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].plot(psi_grid, V0, 'b-', lw=2, label='V_UV')
    axes[0].plot(psi_grid, V_IR, 'r-', lw=2, label='V_IR')
    axes[0].axvline(1.0, color='gray', ls=':', alpha=0.5)
    axes[0].set_xlabel(r'$\psi$'); axes[0].set_ylabel(r'$V(\psi)$')
    axes[0].set_title('UV vs IR'); axes[0].legend()
    axes[0].set_ylim(-1.5, 1.5); axes[0].grid(True, alpha=0.3)

    V1_flow = [sol.sol(t)[idx_vac] for t in t_pts]
    axes[1].plot(t_pts, V1_flow, 'g-o', ms=3)
    axes[1].set_xlabel('t=ln(k/k_IR)'); axes[1].set_ylabel(r'$V_k(\psi=1)$')
    axes[1].set_title('V at vacuum'); axes[1].grid(True, alpha=0.3)

    vpp_flow = [Vpp_num(sol.sol(t), dpsi)[idx_vac] for t in t_pts]
    axes[2].plot(t_pts, vpp_flow, 'm-o', ms=3)
    axes[2].set_xlabel('t=ln(k/k_IR)'); axes[2].set_ylabel(r"$V''(\psi=1)$")
    axes[2].set_title('Mass^2 flow'); axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    sd = os.path.dirname(os.path.abspath(__file__))
    plt.savefig(os.path.join(sd, 'erg_phi_run.png'), dpi=150)
    print(f"\n  Plot: erg_phi_run.png")

print(f"\n{'='*65}")
print("DONE.")
