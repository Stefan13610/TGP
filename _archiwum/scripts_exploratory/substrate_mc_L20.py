#!/usr/bin/env python3
"""MC Ising 3D -- L=20, skalowanie skonczonyego rozmiaru."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L = 20
J = 1.0
N_therm = 3000
N_meas = 10000
N_skip = 3

T_c_exact = 4.5115
# Gestsza siatka blisko T_c
T_values = np.array([3.5, 4.0, 4.2, 4.35, 4.45, 4.50, 4.55, 4.65, 4.8, 5.0, 5.5])

print("=" * 65)
print(f"  MC ISING 3D -- L={L} (skalowanie)")
print("=" * 65)
print(f"  N_therm={N_therm}, N_meas={N_meas}")

class Ising3D:
    def __init__(self, L, J=1.0):
        self.L = L; self.J = J; self.V = L**3
        self.spins = np.random.choice([-1, 1], size=(L, L, L))
    def magnetization(self):
        return np.mean(self.spins)
    def sweep_wolff(self, beta):
        L = self.L
        p_add = 1 - np.exp(-2 * beta * self.J)
        i0, j0, k0 = np.random.randint(0, L, 3)
        s0 = self.spins[i0, j0, k0]
        cluster = set(); stack = [(i0, j0, k0)]; cluster.add((i0, j0, k0))
        while stack:
            i, j, k = stack.pop()
            for di, dj, dk in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                ni, nj, nk = (i+di)%L, (j+dj)%L, (k+dk)%L
                if (ni, nj, nk) not in cluster and self.spins[ni, nj, nk] == s0:
                    if np.random.rand() < p_add:
                        cluster.add((ni, nj, nk)); stack.append((ni, nj, nk))
        for (i, j, k) in cluster:
            self.spins[i, j, k] *= -1
        return len(cluster)
    def sweep_metropolis(self, beta):
        L = self.L
        for _ in range(self.V):
            i, j, k = np.random.randint(0, L, 3)
            s = self.spins[i, j, k]
            nn = (self.spins[(i+1)%L, j, k] + self.spins[(i-1)%L, j, k] +
                  self.spins[i, (j+1)%L, k] + self.spins[i, (j-1)%L, k] +
                  self.spins[i, j, (k+1)%L] + self.spins[i, j, (k-1)%L])
            dE = 2 * self.J * s * nn
            if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
                self.spins[i, j, k] = -s

results = {}
for T in T_values:
    beta = 1.0 / T
    print(f"  T={T:.2f}...", end="", flush=True)
    lat = Ising3D(L, J)
    use_wolff = abs(T - T_c_exact) < 1.5
    for _ in range(N_therm):
        if use_wolff: lat.sweep_wolff(beta)
        else: lat.sweep_metropolis(beta)
    m_list = []
    for n in range(N_meas):
        for _ in range(N_skip):
            if use_wolff: lat.sweep_wolff(beta)
            else: lat.sweep_metropolis(beta)
        m_list.append(lat.magnetization())
    m_arr = np.array(m_list)
    m_abs = np.mean(np.abs(m_arr))
    m2 = np.mean(m_arr**2); m4 = np.mean(m_arr**4)
    chi = lat.V * (m2 - m_abs**2) / T
    U4 = 1 - m4 / (3 * m2**2) if m2 > 1e-15 else 0
    results[T] = {'m_abs': m_abs, 'm2': m2, 'm4': m4, 'chi': chi, 'U4': U4, 'm_hist': m_arr}
    print(f"  <|m|>={m_abs:.4f}, chi={chi:.1f}, U4={U4:.4f}")

# V_eff blisko T_c
print(f"\n--- Fit V_eff(m) ---")
T_fit = min(T_values, key=lambda T: abs(T - T_c_exact))
m_arr = results[T_fit]['m_hist']
hist, be = np.histogram(m_arr, bins=80, density=True)
bc = 0.5*(be[1:]+be[:-1]); mask = hist > 0
V_eff = -np.log(hist[mask]); V_eff -= np.min(V_eff); m_fit = bc[mask]
m_pos = m_fit[m_fit > 0]; V_pos = V_eff[m_fit > 0]
if len(m_pos) > 5:
    X = np.column_stack([m_pos**2, m_pos**4, m_pos**6])
    a2, a4, a6 = np.linalg.lstsq(X, V_pos, rcond=None)[0]
    print(f"  T_fit={T_fit:.2f}: V_eff = {a2:.1f}*m^2 + {a4:.1f}*m^4 + {a6:.1f}*m^6")
    print(f"  u4={a4:.2f}, u6={a6:.2f}, u6/u4^2={a6/a4**2:.4f}" if abs(a4)>1e-10 else "")

# Skalowanie: porownanie z L=10
print(f"\n--- Skalowanie L={L} ---")
T_arr = sorted(results.keys())
chi_arr = [results[T]['chi'] for T in T_arr]
T_c_mc = T_arr[np.argmax(chi_arr)]
print(f"  T_c(MC, L={L}) = {T_c_mc:.2f}")
print(f"  T_c(exact) = {T_c_exact}")
print(f"  Odchylenie = {abs(T_c_mc-T_c_exact)/T_c_exact*100:.2f}%")
# Binder crossing
U4_arr = [results[T]['U4'] for T in T_arr]
# Find where U4 crosses 0.465 (3D Ising universal value)
U4_target = 0.465
for i in range(len(T_arr)-1):
    if (U4_arr[i] - U4_target) * (U4_arr[i+1] - U4_target) < 0:
        # Linear interpolation
        T_cross = T_arr[i] + (U4_target - U4_arr[i]) / (U4_arr[i+1] - U4_arr[i]) * (T_arr[i+1] - T_arr[i])
        print(f"  T_c(Binder crossing U4=0.465) = {T_cross:.3f}")
        break

# Wykresy
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]
ax.plot(T_arr, [results[T]['m_abs'] for T in T_arr], 'bo-')
ax.axvline(T_c_exact, color='red', ls='--', alpha=0.5)
ax.set_xlabel('T'); ax.set_ylabel(r'$\langle |m| \rangle$'); ax.set_title(f'Magnetyzacja L={L}')

ax = axes[1]
ax.plot(T_arr, chi_arr, 'rs-')
ax.axvline(T_c_exact, color='red', ls='--', alpha=0.5)
ax.set_xlabel('T'); ax.set_ylabel(r'$\chi$'); ax.set_title(f'Susceptibility L={L}')

ax = axes[2]
ax.plot(T_arr, U4_arr, 'go-')
ax.axhline(0.465, color='gray', ls=':', label='U4*=0.465')
ax.axvline(T_c_exact, color='red', ls='--', alpha=0.5)
ax.set_xlabel('T'); ax.set_ylabel('U4'); ax.set_title(f'Binder cumulant L={L}')
ax.legend()

plt.tight_layout()
sd = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(sd, 'substrate_mc_L20.png'), dpi=150, bbox_inches='tight')
print(f"\n  Wykres: substrate_mc_L20.png")
print("GOTOWE.")
