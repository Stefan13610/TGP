#!/usr/bin/env python3
"""MC Ising 3D -- szybka wersja (mniejsza siec, mniej pomiarow)."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L = 10           # mniejsza siec
J = 1.0
N_therm = 2000
N_meas = 8000
N_skip = 3

T_c_exact = 4.5115
T_values = np.array([3.0, 3.5, 4.0, 4.3, 4.5, 4.6, 4.8, 5.0, 5.5, 6.0])

print("=" * 65)
print("  MC SYMULACJA SUBSTRATU TGP (Ising 3D) -- FAST")
print("=" * 65)
print(f"  L = {L}, N_therm = {N_therm}, N_meas = {N_meas}")

class Ising3D:
    def __init__(self, L, J=1.0):
        self.L = L
        self.J = J
        self.V = L**3
        self.spins = np.random.choice([-1, 1], size=(L, L, L))

    def magnetization(self):
        return np.mean(self.spins)

    def energy(self):
        s = self.spins
        E = 0
        for ax in range(3):
            E -= self.J * np.sum(s * np.roll(s, 1, axis=ax))
        return E

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

    def sweep_wolff(self, beta):
        L = self.L
        p_add = 1 - np.exp(-2 * beta * self.J)
        i0, j0, k0 = np.random.randint(0, L, 3)
        s0 = self.spins[i0, j0, k0]
        cluster = set()
        stack = [(i0, j0, k0)]
        cluster.add((i0, j0, k0))
        while stack:
            i, j, k = stack.pop()
            for di, dj, dk in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                ni, nj, nk = (i+di)%L, (j+dj)%L, (k+dk)%L
                if (ni, nj, nk) not in cluster and self.spins[ni, nj, nk] == s0:
                    if np.random.rand() < p_add:
                        cluster.add((ni, nj, nk))
                        stack.append((ni, nj, nk))
        for (i, j, k) in cluster:
            self.spins[i, j, k] *= -1
        return len(cluster)

results = {}
for T in T_values:
    beta = 1.0 / T
    print(f"  T = {T:.2f}...", end="", flush=True)
    lattice = Ising3D(L, J)
    use_wolff = abs(T - T_c_exact) < 1.5
    for _ in range(N_therm):
        if use_wolff:
            lattice.sweep_wolff(beta)
        else:
            lattice.sweep_metropolis(beta)
    m_list = []
    for n in range(N_meas):
        for _ in range(N_skip):
            if use_wolff:
                lattice.sweep_wolff(beta)
            else:
                lattice.sweep_metropolis(beta)
        m_list.append(lattice.magnetization())
    m_arr = np.array(m_list)
    m_abs = np.mean(np.abs(m_arr))
    m2 = np.mean(m_arr**2)
    m4 = np.mean(m_arr**4)
    chi = lattice.V * (m2 - m_abs**2) / T
    U4 = 1 - m4 / (3 * m2**2) if m2 > 1e-15 else 0
    results[T] = {'m_abs': m_abs, 'm2': m2, 'm4': m4, 'chi': chi, 'U4': U4, 'm_hist': m_arr}
    print(f"  <|m|>={m_abs:.4f}, chi={chi:.1f}, U4={U4:.4f}")

# Potencjal efektywny
print(f"\n--- Potencjal efektywny V_eff(m) ---")
T_plot = [T for T in T_values if 4.0 <= T <= 5.5]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax1 = axes[0, 0]
for T in T_plot:
    ax1.hist(results[T]['m_hist'], bins=60, density=True, alpha=0.5, label=f'T={T:.1f}')
ax1.set_xlabel('m'); ax1.set_ylabel('P(m)'); ax1.set_title('Rozklad magnetyzacji')
ax1.legend(fontsize=8)

ax2 = axes[0, 1]
for T in T_plot:
    hist, bin_edges = np.histogram(results[T]['m_hist'], bins=50, density=True)
    bc = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    mask = hist > 0
    V_eff = -np.log(hist[mask]) / L**3
    V_eff -= np.min(V_eff)
    ax2.plot(bc[mask], V_eff, 'o-', markersize=2, label=f'T={T:.1f}')
ax2.set_xlabel('m'); ax2.set_ylabel(r'$V_{\rm eff}(m)/V$')
ax2.set_title('Potencjal efektywny z MC')
ax2.set_ylim(0, 0.01); ax2.legend(fontsize=8)

ax3 = axes[1, 0]
T_arr = sorted(results.keys())
m_arr_T = [results[T]['m_abs'] for T in T_arr]
chi_arr = [results[T]['chi'] for T in T_arr]
ax3.plot(T_arr, m_arr_T, 'bo-', label=r'$\langle |m| \rangle$')
ax3.axvline(T_c_exact, color='red', ls='--', alpha=0.5, label=f'$T_c$={T_c_exact}')
ax3.set_xlabel('T'); ax3.set_ylabel(r'$\langle |m| \rangle$')
ax3.set_title('Magnetyzacja vs T'); ax3.legend()
ax3b = ax3.twinx()
ax3b.plot(T_arr, chi_arr, 'rs-', alpha=0.7)
ax3b.set_ylabel(r'$\chi$', color='red')

ax4 = axes[1, 1]
U4_arr = [results[T]['U4'] for T in T_arr]
ax4.plot(T_arr, U4_arr, 'go-')
ax4.axhline(2/3, color='gray', ls=':', alpha=0.5, label='U4=2/3')
ax4.axhline(0, color='gray', ls=':', alpha=0.5, label='U4=0')
ax4.axvline(T_c_exact, color='red', ls='--', alpha=0.5)
ax4.set_xlabel('T'); ax4.set_ylabel('U4'); ax4.set_title('Binder cumulant')
ax4.legend(fontsize=8)

plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'substrate_mc_fast.png'), dpi=150, bbox_inches='tight')
print(f"  Wykres: scripts/substrate_mc_fast.png")

# Fit V_eff blisko T_c
print(f"\n--- Fit V_eff(m) blisko T_c ---")
T_fit = min(T_values, key=lambda T: abs(T - T_c_exact))
m_arr = results[T_fit]['m_hist']
hist, bin_edges = np.histogram(m_arr, bins=50, density=True)
bc = 0.5 * (bin_edges[1:] + bin_edges[:-1])
mask = hist > 0
V_eff = -np.log(hist[mask])
V_eff -= np.min(V_eff)
m_fit = bc[mask]

m_pos = m_fit[m_fit > 0]
V_pos = V_eff[m_fit > 0]
if len(m_pos) > 5:
    X = np.column_stack([m_pos**2, m_pos**4, m_pos**6])
    coeffs = np.linalg.lstsq(X, V_pos, rcond=None)[0]
    a2, a4, a6 = coeffs
    print(f"  T_fit = {T_fit:.2f}")
    print(f"  V_eff(m) = {a2:.2f}*m^2 + {a4:.2f}*m^4 + {a6:.2f}*m^6")
    print(f"  u4_eff = {a4:.4f}")
    print(f"  u6_eff = {a6:.4f}")
    if abs(a4) > 1e-10:
        print(f"  u6/u4^2 = {a6/a4**2:.4f}")

# Podsumowanie
print(f"\n{'='*65}")
print(f"  PODSUMOWANIE MC Z2")
print(f"{'='*65}")
T_arr_s = sorted(results.keys())
chi_s = [results[T]['chi'] for T in T_arr_s]
T_c_mc = T_arr_s[np.argmax(chi_s)]
print(f"  T_c (MC, max chi) = {T_c_mc:.2f}")
print(f"  T_c (exact) = {T_c_exact}")
print(f"  Odchylenie = {abs(T_c_mc - T_c_exact)/T_c_exact*100:.1f}%")
print(f"\n  Interpretacja TGP:")
print(f"  - Przejscie fazowe 3D Ising = coarse-graining substratu Gamma")
print(f"  - V_eff(m) z MC ~ V_TGP(psi) po identyfikacji m <-> psi")
print(f"  - Krytyczne skalowanie: nu=0.6302, eta=0.0364 (klasa 3D Ising)")
print(f"  - TGP potrzebuje: nu_TGP ~ 0.63 z V(psi) = psi^3/3 - psi^4/4")
print(f"\nGOTOWE.")
