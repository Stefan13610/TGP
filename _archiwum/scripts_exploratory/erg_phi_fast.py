#!/usr/bin/env python3
"""ERG (Wetterich) LPA -- szybka wersja (mniejsza siatka)."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parametry TGP
Phi0 = 24.66
gamma_TGP = 1.0
a_Gamma = 0.040
k_UV = 1.0 / a_Gamma   # = 25
k_IR = np.sqrt(gamma_TGP)  # = 1
t_max = np.log(k_UV / k_IR)

print("=" * 65)
print("  ERG (Wetterich) DLA POLA TGP -- FAST")
print("=" * 65)
print(f"  k_UV = {k_UV:.1f}, k_IR = {k_IR:.1f}, t_max = {t_max:.2f}")

# Siatka -- mniejsza dla szybkosci
N_grid = 60
psi_min, psi_max = 0.01, 3.0
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

def V_TGP(psi):
    return psi**3 / 3 - psi**4 / 4

def V_second_derivative(V_arr, dpsi):
    Vpp = np.zeros_like(V_arr)
    Vpp[1:-1] = (V_arr[2:] - 2*V_arr[1:-1] + V_arr[:-2]) / dpsi**2
    Vpp[0] = Vpp[1]
    Vpp[-1] = Vpp[-2]
    return Vpp

def rg_rhs(t, V_flat, psi_grid, dpsi, k_IR):
    V = V_flat.copy()
    k = k_IR * np.exp(t)
    Vpp = V_second_derivative(V, dpsi)
    denominator = np.maximum(Vpp + k**2, 1e-10)
    dVdt = k**5 / (32 * np.pi**2) / denominator
    return dVdt

V_UV = V_TGP(psi_grid)
print(f"  Siatka: {N_grid} punktow, dpsi = {dpsi:.4f}")

t_span = (t_max, 0.01)
t_eval = np.linspace(t_max, 0.01, 30)

print("  Integracja UV -> IR ...", flush=True)
sol = solve_ivp(
    lambda t, V: rg_rhs(t, V, psi_grid, dpsi, k_IR),
    t_span, V_UV, t_eval=t_eval,
    method='RK45', rtol=1e-5, atol=1e-7, max_step=0.2
)

if not sol.success:
    print(f"  BLAD: {sol.message}")
    sys.exit(1)

print(f"  Integracja: SUKCES ({sol.t.shape[0]} krokow)")
V_IR = sol.y[:, -1]

# Analiza
print(f"\n--- Analiza potencjalu IR ---")
idx_min = np.argmin(V_IR)
psi_min_val = psi_grid[idx_min]
V_min = V_IR[idx_min]
print(f"  Minimum V_IR: psi = {psi_min_val:.4f}, V = {V_min:.6f}")

Vpp_IR = V_second_derivative(V_IR, dpsi)
m2_IR = Vpp_IR[idx_min]
print(f"  V''(psi_min) = m^2 = {m2_IR:.6f}")
if m2_IR > 0:
    print(f"  -> Masa: m = {np.sqrt(m2_IR):.4f} (stabilna)")
else:
    print(f"  -> m^2 < 0: niestabilne!")

psi_vac = 1.0
idx_vac = np.argmin(np.abs(psi_grid - psi_vac))
print(f"\n  Na psi = 1 (proznia TGP):")
print(f"    V_UV(1) = {V_UV[idx_vac]:.6f}")
print(f"    V_IR(1) = {V_IR[idx_vac]:.6f}")
print(f"    Delta   = {(V_IR[idx_vac] - V_UV[idx_vac]):.6f}")

# Czy potencjal zachowuje strukture TGP?
# Szukamy minimum V_IR blisko psi=1
mask_near = (psi_grid > 0.5) & (psi_grid < 1.5)
idx_local_min = np.argmin(V_IR[mask_near])
psi_local = psi_grid[mask_near][idx_local_min]
V_local = V_IR[mask_near][idx_local_min]
print(f"\n  Lokalne minimum V_IR blisko psi=1:")
print(f"    psi_min = {psi_local:.4f}")
print(f"    V_min   = {V_local:.6f}")

# Sprawdz bariere miedzy psi=0 a psi=1
if psi_local > 0.3:
    idx_barrier = np.argmax(V_IR[(psi_grid > 0.1) & (psi_grid < psi_local)])
    psi_barrier_range = psi_grid[(psi_grid > 0.1) & (psi_grid < psi_local)]
    V_barrier_range = V_IR[(psi_grid > 0.1) & (psi_grid < psi_local)]
    if len(V_barrier_range) > 0:
        idx_b = np.argmax(V_barrier_range)
        print(f"    Bariera: psi = {psi_barrier_range[idx_b]:.4f}, V = {V_barrier_range[idx_b]:.6f}")
        print(f"    Wysokosc bariery: {V_barrier_range[idx_b] - V_local:.6f}")

# Ograniczonosc
V_UV_max = np.max(np.abs(V_UV))
V_IR_max = np.max(np.abs(V_IR[psi_grid < 2.0]))
print(f"\n  Ograniczonosc:")
print(f"    max|V_UV| = {V_UV_max:.4f}")
print(f"    max|V_IR| (psi<2) = {V_IR_max:.4f}")
print(f"    Stosunek = {V_IR_max / V_UV_max:.4f}")

# Przeplyw V''(psi=1) -- poszukiwanie UV fixed point
Vpp_flow = []
for i in range(sol.t.shape[0]):
    Vpp_i = V_second_derivative(sol.y[:, i], dpsi)
    Vpp_flow.append(Vpp_i[idx_vac])
Vpp_flow = np.array(Vpp_flow)
print(f"\n  Przeplyw V''(psi=1):")
print(f"    UV (t={sol.t[0]:.2f}): {Vpp_flow[0]:.6f}")
print(f"    IR (t={sol.t[-1]:.2f}): {Vpp_flow[-1]:.6f}")
print(f"    Zmiana: {(Vpp_flow[-1] - Vpp_flow[0]):.6f}")

# Bezwymiarowy potencjal (szukanie AS)
# u_k = V_k / k^4  (bezwymiarowy)
print(f"\n  Bezwymiarowy potencjal u = V/k^4 na psi=1:")
for i in [0, len(sol.t)//4, len(sol.t)//2, 3*len(sol.t)//4, -1]:
    k_i = k_IR * np.exp(sol.t[i])
    u_i = sol.y[idx_vac, i] / k_i**4
    print(f"    t={sol.t[i]:.2f} (k={k_i:.1f}): u = {u_i:.2e}")

# Wykresy
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax1 = axes[0]
ax1.plot(psi_grid, V_UV, 'b-', lw=2, label='V_UV (bare TGP)')
ax1.plot(psi_grid, V_IR, 'r-', lw=2, label='V_IR (po RG flow)')
ax1.axvline(1.0, color='gray', ls=':', alpha=0.5)
ax1.set_xlabel(r'$\psi$')
ax1.set_ylabel(r'$V_k(\psi)$')
ax1.set_title('Potencjal: UV vs IR')
ax1.legend()
ax1.set_ylim(-2, 2)
ax1.grid(True, alpha=0.3)

ax2 = axes[1]
V_at_1 = [sol.y[idx_vac, i] for i in range(sol.t.shape[0])]
ax2.plot(sol.t, V_at_1, 'g-o', markersize=3)
ax2.set_xlabel('t = ln(k/k_IR)')
ax2.set_ylabel(r'$V_k(\psi=1)$')
ax2.set_title('Przeplyw V na prozni')
ax2.grid(True, alpha=0.3)

ax3 = axes[2]
ax3.plot(sol.t, Vpp_flow, 'm-o', markersize=3)
ax3.set_xlabel('t = ln(k/k_IR)')
ax3.set_ylabel(r"$V''_k(\psi=1)$")
ax3.set_title('Masa^2 pola: przeplyw RG')
ax3.grid(True, alpha=0.3)

plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'erg_phi_fast.png'), dpi=150)
print(f"\n  Wykres: scripts/erg_phi_fast.png")

print(f"\n{'='*65}")
print(f"  PODSUMOWANIE ERG (FAST)")
print(f"{'='*65}")
print(f"  Potencjal IR zachowuje strukture TGP: {'TAK' if psi_local > 0.7 and psi_local < 1.3 else 'NIE/NIEJASNE'}")
print(f"  Minimum przesuniete do psi = {psi_local:.4f} (oczek. 1.0)")
print(f"  Masa^2 na IR: {m2_IR:.4f}")
print(f"  Stosunek V_IR/V_UV (max): {V_IR_max / V_UV_max:.2f}")
print(f"\nGOTOWE.")
