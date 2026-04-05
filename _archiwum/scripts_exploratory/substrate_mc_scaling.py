#!/usr/bin/env python3
"""MC Ising 3D -- skalowanie L=8,12,16 (szybkie, porownawcze)."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T_c_exact = 4.5115
# Temperatury blisko T_c (gestsze)
T_values = np.array([4.0, 4.2, 4.35, 4.45, 4.50, 4.55, 4.60, 4.75, 5.0])

L_values = [8, 12, 16]
J = 1.0

print("=" * 70)
print("  MC ISING 3D -- SKALOWANIE SKONCZONYEGO ROZMIARU")
print("=" * 70)

class Ising3D:
    def __init__(self, L):
        self.L = L; self.V = L**3
        self.spins = np.random.choice([-1, 1], size=(L, L, L))
    def magnetization(self):
        return np.mean(self.spins)
    def sweep_wolff(self, beta):
        L = self.L; p_add = 1 - np.exp(-2*beta*J)
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

all_results = {}

for L in L_values:
    N_therm = 2000
    N_meas = 5000
    N_skip = max(2, L // 4)
    print(f"\n  === L = {L} (V={L**3}) ===")
    results = {}
    for T in T_values:
        beta = 1.0 / T
        lat = Ising3D(L)
        for _ in range(N_therm):
            lat.sweep_wolff(beta)
        m_list = []
        for n in range(N_meas):
            for _ in range(N_skip):
                lat.sweep_wolff(beta)
            m_list.append(lat.magnetization())
        m_arr = np.array(m_list)
        m_abs = np.mean(np.abs(m_arr))
        m2 = np.mean(m_arr**2); m4 = np.mean(m_arr**4)
        chi = lat.V * (m2 - m_abs**2) / T
        U4 = 1 - m4/(3*m2**2) if m2 > 1e-15 else 0
        results[T] = {'m_abs': m_abs, 'chi': chi, 'U4': U4, 'm2': m2}
        print(f"    T={T:.2f}: <|m|>={m_abs:.4f}, chi={chi:.1f}, U4={U4:.4f}")
    all_results[L] = results

# ============================================================
# Analiza skalowania
# ============================================================
print(f"\n{'='*70}")
print("  ANALIZA SKALOWANIA")
print("="*70)

# 1. T_c z max chi
print(f"\n  T_c z max chi:")
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    chi_arr = [all_results[L][T]['chi'] for T in T_arr]
    T_c_mc = T_arr[np.argmax(chi_arr)]
    chi_max = max(chi_arr)
    print(f"    L={L:3d}: T_c = {T_c_mc:.2f}, chi_max = {chi_max:.1f}")

# 2. Binder crossing
print(f"\n  Binder cumulant crossing (U4 = 0.465):")
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    U4_arr = [all_results[L][T]['U4'] for T in T_arr]
    for i in range(len(T_arr)-1):
        if (U4_arr[i] - 0.465) * (U4_arr[i+1] - 0.465) < 0:
            T_cross = T_arr[i] + (0.465 - U4_arr[i])/(U4_arr[i+1]-U4_arr[i])*(T_arr[i+1]-T_arr[i])
            print(f"    L={L:3d}: T_cross = {T_cross:.3f}")
            break
    else:
        print(f"    L={L:3d}: nie znaleziono crossingu")

# 3. chi_max skalowanie: chi_max ~ L^(gamma/nu), gamma/nu = 1.9633 (3D Ising)
print(f"\n  Skalowanie chi_max ~ L^(gamma/nu):")
chi_maxs = []
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    chi_arr = [all_results[L][T]['chi'] for T in T_arr]
    chi_maxs.append(max(chi_arr))
    print(f"    L={L:3d}: chi_max = {max(chi_arr):.2f}")
if len(chi_maxs) >= 2:
    # Fit: ln(chi_max) = (gamma/nu)*ln(L) + const
    log_L = np.log(np.array(L_values, dtype=float))
    log_chi = np.log(np.array(chi_maxs))
    coeffs = np.polyfit(log_L, log_chi, 1)
    gamma_nu_mc = coeffs[0]
    gamma_nu_exact = 1.9633
    print(f"    gamma/nu (fit) = {gamma_nu_mc:.3f}")
    print(f"    gamma/nu (exact 3D Ising) = {gamma_nu_exact:.4f}")
    print(f"    Odchylenie = {abs(gamma_nu_mc - gamma_nu_exact)/gamma_nu_exact*100:.1f}%")

# 4. Porownanie z TGP
print(f"\n  Interpretacja TGP:")
print(f"  - Substrat Z_2 daje przejscie 3D Ising (T_c ~ 4.51)")
print(f"  - Skalowanie: gamma/nu powinno byc ~ 2 (3D Ising)")
print(f"  - TGP V(psi): potencjal efektywny emerguje z coarse-grainingu")
print(f"  - Wykladniki krytyczne: nu=0.6302, gamma=1.2372, eta=0.0364")

# ============================================================
# Wykresy
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
colors = {8: 'blue', 12: 'green', 16: 'red'}

# Magnetyzacja
ax = axes[0]
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    m_arr = [all_results[L][T]['m_abs'] for T in T_arr]
    ax.plot(T_arr, m_arr, 'o-', color=colors[L], label=f'L={L}')
ax.axvline(T_c_exact, color='gray', ls='--', alpha=0.5, label=f'T_c={T_c_exact}')
ax.set_xlabel('T'); ax.set_ylabel(r'$\langle |m| \rangle$')
ax.set_title('Magnetyzacja'); ax.legend()

# chi
ax = axes[1]
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    chi_arr = [all_results[L][T]['chi'] for T in T_arr]
    ax.plot(T_arr, chi_arr, 's-', color=colors[L], label=f'L={L}')
ax.axvline(T_c_exact, color='gray', ls='--', alpha=0.5)
ax.set_xlabel('T'); ax.set_ylabel(r'$\chi$')
ax.set_title('Susceptibility (scaling)'); ax.legend()

# Binder
ax = axes[2]
for L in L_values:
    T_arr = sorted(all_results[L].keys())
    U4_arr = [all_results[L][T]['U4'] for T in T_arr]
    ax.plot(T_arr, U4_arr, 'o-', color=colors[L], label=f'L={L}')
ax.axhline(0.465, color='gray', ls=':', alpha=0.5, label='U4*=0.465')
ax.axvline(T_c_exact, color='gray', ls='--', alpha=0.5)
ax.set_xlabel('T'); ax.set_ylabel('U4')
ax.set_title('Binder crossing'); ax.legend()

plt.tight_layout()
sd = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(sd, 'substrate_mc_scaling.png'), dpi=150, bbox_inches='tight')
print(f"\n  Wykres: substrate_mc_scaling.png")
print("GOTOWE.")
