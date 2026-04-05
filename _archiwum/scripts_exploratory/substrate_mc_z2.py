#!/usr/bin/env python3
"""
Monte Carlo symulacja substratu TGP: model Isinga/GL na sieci 3D
=================================================================
Cel: weryfikacja potencjalu V(psi) z mikrodynamiki substratu.

Model: H = -J * sum_{<ij>} s_i*s_j + h * sum_i s_i^4
       (model phi^4 na sieci kubicznej 3D z Z_2)

Mierzymy:
1. Efektywny potencjal V_eff(m) = -ln(P(m))/V  (m = magnetyzacja)
2. Susceptywalnosc chi(T)
3. Binder cumulant U_4 = 1 - <m^4>/(3*<m^2>^2)
4. Wspolczynniki u4, u6 z rozkladu V_eff(m) wokol minimum

Referencja: sek01 (H_Gamma), sek08 (coarse-graining)
TGP v1 -- 2026-03-31
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

# ============================================================
# 1. Parametry symulacji
# ============================================================
L = 16           # rozmiar sieci (L x L x L)
J = 1.0          # sprzezenie
N_therm = 5000   # termalizacja
N_meas = 20000   # pomiary
N_skip = 5       # odleg. miedzy pomiarami (dekorelacja)

# Temperatura krytyczna 3D Ising: T_c ~ 4.5115 * J (kubiczna)
T_c_exact = 4.5115
T_values = np.array([3.0, 3.5, 4.0, 4.2, 4.4, 4.5, 4.6, 4.8, 5.0, 5.5, 6.0, 7.0])

print("=" * 65)
print("  MC SYMULACJA SUBSTRATU TGP (model Isinga 3D)")
print("=" * 65)
print(f"  L = {L}, J = {J}")
print(f"  T_c(3D Ising) = {T_c_exact}")
print(f"  N_therm = {N_therm}, N_meas = {N_meas}")
print(f"  Temperatury: {T_values}")
print()

# ============================================================
# 2. Algorytm Metropolis na sieci 3D
# ============================================================
class Ising3D:
    def __init__(self, L, J=1.0):
        self.L = L
        self.J = J
        self.V = L**3
        # Inicjalizacja: losowa konfiguracja
        self.spins = np.random.choice([-1, 1], size=(L, L, L))

    def energy(self):
        """Calkowita energia."""
        s = self.spins
        E = 0
        for ax in range(3):
            E -= self.J * np.sum(s * np.roll(s, 1, axis=ax))
        return E

    def magnetization(self):
        """Magnetyzacja na spin."""
        return np.mean(self.spins)

    def sweep_metropolis(self, beta):
        """Jeden sweep Metropolis (N prób)."""
        L = self.L
        for _ in range(self.V):
            i, j, k = np.random.randint(0, L, 3)
            s = self.spins[i, j, k]
            # Suma sasiadow
            nn = (self.spins[(i+1)%L, j, k] + self.spins[(i-1)%L, j, k] +
                  self.spins[i, (j+1)%L, k] + self.spins[i, (j-1)%L, k] +
                  self.spins[i, j, (k+1)%L] + self.spins[i, j, (k-1)%L])
            dE = 2 * self.J * s * nn
            if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
                self.spins[i, j, k] = -s

    def sweep_wolff(self, beta):
        """Jeden klaster Wolffa (szybsza termalizacja blisko T_c)."""
        L = self.L
        p_add = 1 - np.exp(-2 * beta * self.J)
        # Losowy punkt startowy
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

        # Flip cluster
        for (i, j, k) in cluster:
            self.spins[i, j, k] *= -1

        return len(cluster)

# ============================================================
# 3. Pomiary
# ============================================================
results = {}

for T in T_values:
    beta = 1.0 / T
    print(f"  T = {T:.2f} (beta = {beta:.4f})...", end="", flush=True)

    lattice = Ising3D(L, J)

    # Termalizacja (Wolff blisko T_c, Metropolis daleko)
    use_wolff = abs(T - T_c_exact) < 1.5
    for _ in range(N_therm):
        if use_wolff:
            lattice.sweep_wolff(beta)
        else:
            lattice.sweep_metropolis(beta)

    # Pomiary
    m_list = []
    e_list = []
    for n in range(N_meas):
        for _ in range(N_skip):
            if use_wolff:
                lattice.sweep_wolff(beta)
            else:
                lattice.sweep_metropolis(beta)
        m = lattice.magnetization()
        m_list.append(m)
        if n % 5000 == 0:
            e_list.append(lattice.energy() / lattice.V)

    m_arr = np.array(m_list)
    e_arr = np.array(e_list) if e_list else np.array([0])

    # Obserwable
    m_abs = np.mean(np.abs(m_arr))
    m2 = np.mean(m_arr**2)
    m4 = np.mean(m_arr**4)
    chi = lattice.V * (m2 - np.mean(np.abs(m_arr))**2) / T
    U4 = 1 - m4 / (3 * m2**2) if m2 > 1e-15 else 0

    results[T] = {
        'm_abs': m_abs, 'm2': m2, 'm4': m4,
        'chi': chi, 'U4': U4,
        'e_mean': np.mean(e_arr),
        'm_hist': m_arr
    }

    print(f"  <|m|>={m_abs:.4f}, chi={chi:.2f}, U4={U4:.4f}")

# ============================================================
# 4. Potencjal efektywny V_eff(m)
# ============================================================
print(f"\n--- Potencjal efektywny V_eff(m) ---")

# Wybieramy T blisko T_c (faza symetryczna i zlamana)
T_plot = [T for T in T_values if 4.0 <= T <= 5.5]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Histogram magnetyzacji dla roznych T
ax1 = axes[0, 0]
for T in T_plot:
    m_arr = results[T]['m_hist']
    ax1.hist(m_arr, bins=80, density=True, alpha=0.5, label=f'T={T:.1f}')
ax1.set_xlabel('m (magnetyzacja)')
ax1.set_ylabel('P(m)')
ax1.set_title('Rozklad magnetyzacji')
ax1.legend(fontsize=8)

# Panel 2: V_eff(m) = -ln(P(m)) / V
ax2 = axes[0, 1]
for T in T_plot:
    m_arr = results[T]['m_hist']
    hist, bin_edges = np.histogram(m_arr, bins=60, density=True)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    mask = hist > 0
    V_eff = -np.log(hist[mask]) / L**3
    V_eff -= np.min(V_eff)  # normalizacja
    ax2.plot(bin_centers[mask], V_eff, 'o-', markersize=2, label=f'T={T:.1f}')
ax2.set_xlabel('m')
ax2.set_ylabel(r'$V_{\rm eff}(m) / V$')
ax2.set_title('Potencjal efektywny z MC')
ax2.set_ylim(0, 0.005)
ax2.legend(fontsize=8)

# Panel 3: <|m|>(T) i chi(T)
ax3 = axes[1, 0]
T_arr = sorted(results.keys())
m_arr_T = [results[T]['m_abs'] for T in T_arr]
chi_arr = [results[T]['chi'] for T in T_arr]
ax3.plot(T_arr, m_arr_T, 'bo-', label=r'$\langle |m| \rangle$')
ax3.axvline(T_c_exact, color='red', ls='--', alpha=0.5, label=f'$T_c$={T_c_exact}')
ax3.set_xlabel('T')
ax3.set_ylabel(r'$\langle |m| \rangle$')
ax3.set_title('Magnetyzacja vs T')
ax3.legend()

ax3b = ax3.twinx()
ax3b.plot(T_arr, chi_arr, 'rs-', alpha=0.7, label=r'$\chi$')
ax3b.set_ylabel(r'$\chi$ (susceptibility)', color='red')
ax3b.legend(loc='upper left')

# Panel 4: Binder cumulant
ax4 = axes[1, 1]
U4_arr = [results[T]['U4'] for T in T_arr]
ax4.plot(T_arr, U4_arr, 'go-')
ax4.axhline(2/3, color='gray', ls=':', alpha=0.5, label='U4 = 2/3 (ordered)')
ax4.axhline(0, color='gray', ls=':', alpha=0.5, label='U4 = 0 (disordered)')
ax4.axvline(T_c_exact, color='red', ls='--', alpha=0.5, label=f'$T_c$={T_c_exact}')
ax4.set_xlabel('T')
ax4.set_ylabel('U4 (Binder cumulant)')
ax4.set_title('Binder cumulant vs T')
ax4.legend(fontsize=8)

plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'substrate_mc_z2.png'), dpi=150, bbox_inches='tight')
print(f"  Wykres: scripts/substrate_mc_z2.png")

# ============================================================
# 5. Fit potencjalu V_eff(m) blisko T_c
# ============================================================
print(f"\n--- Fit V_eff(m) blisko T_c ---")

T_fit = min(T_values, key=lambda T: abs(T - T_c_exact))
m_arr = results[T_fit]['m_hist']
hist, bin_edges = np.histogram(m_arr, bins=60, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
mask = hist > 0
V_eff = -np.log(hist[mask])
V_eff -= np.min(V_eff)
m_fit = bin_centers[mask]

# Fit: V_eff(m) = a2*m^2 + a4*m^4 + a6*m^6
# (Z_2 symetria: tylko parzyste potegi)
m_pos = m_fit[m_fit > 0]
V_pos = V_eff[m_fit > 0]

if len(m_pos) > 5:
    # Design matrix: [m^2, m^4, m^6]
    X = np.column_stack([m_pos**2, m_pos**4, m_pos**6])
    coeffs = np.linalg.lstsq(X, V_pos, rcond=None)[0]
    a2, a4, a6 = coeffs

    print(f"  T_fit = {T_fit:.2f}")
    print(f"  V_eff(m) = {a2:.2f}*m^2 + {a4:.2f}*m^4 + {a6:.2f}*m^6")
    print(f"  u4_eff = a4 = {a4:.4f}")
    print(f"  u6_eff = a6 = {a6:.4f}")
    if abs(a4) > 1e-10:
        print(f"  u6/u4^2 = {a6/a4**2:.4f}")
        print(f"  u6/u4^(3/2) = {a6/abs(a4)**1.5:.4f}")

    # Porownanie z TGP V(psi) = psi^3/3 - psi^4/4
    # Po przesunieciu psi = 1 + delta: V ~ const + (1/2)*delta^2 + ...
    # V''(1) = gamma (masa^2 pola), V'''(1) = -4gamma (kubiczny)
    print(f"\n  Porownanie z TGP:")
    print(f"  TGP oczekuje: u4 < 0 (zmiana znaku potencjalu w fazie zlamanej)")
    print(f"  MC daje: u4 = {a4:.2f} ({'<0 OK' if a4 < 0 else '>0 faza symetryczna'})")
else:
    print(f"  Za malo danych do fitu")

# ============================================================
# 6. Podsumowanie
# ============================================================
print(f"\n{'='*65}")
print(f"  PODSUMOWANIE MC SUBSTRATU Z2")
print(f"{'='*65}")
print(f"  Siec: {L}^3 = {L**3} spinow")
print(f"  T_c (oczekiwane): {T_c_exact}")

# Wyznacz T_c z max chi
T_arr_sorted = sorted(results.keys())
chi_arr_sorted = [results[T]['chi'] for T in T_arr_sorted]
idx_max_chi = np.argmax(chi_arr_sorted)
T_c_mc = T_arr_sorted[idx_max_chi]
print(f"  T_c (MC, max chi): {T_c_mc:.2f}")
print(f"  Roznica: {abs(T_c_mc - T_c_exact)/T_c_exact*100:.1f}%")

print(f"\n  V_eff(m) zmierzony dla {len(T_values)} temperatur.")
print(f"  Dalsze kroki: L=32,64 (skalowanie), T blizej T_c (precyzyjny fit)")
print(f"\nGOTOWE.")
