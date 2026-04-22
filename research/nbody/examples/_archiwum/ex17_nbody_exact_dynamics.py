"""
ex17_nbody_exact_dynamics.py
=============================
Pelna dynamika N-cialowa TGP z dokladnymi silami 3-cialowymi.

NOVUM wzgledem poprzednich skryptow:
  Poprzednie symulacje (dynamics_v2.py) uzywaly PRZYBLIZONYCH sil
  3-cialowych z bledem 160-770% (saddle-point).

  Ten skrypt uzywa DOKLADNYCH sil z three_body_force_exact.py
  (calka Feynmana, blad 1e-10, Newton 3. zasada: 2.78e-17).

UKLAD:
  3, 4, 5 cial z parametrami TGP.
  Integracja: RK45 (scipy solve_ivp) z kontrola kroku.
  Obserwable: energia calkowita, N-cialowy wklad, orbity.

PYTANIE:
  Jak efekty 3-cialowe TGP zmieniaja orbity w porownaniu z 2-cialowymi?
  Dla jakich C i m_sp sa efekty 3-cialowe obserwowalne?

PARAMETRY:
  C = 0.10-0.50 (testowy zakres dla efektow 3-cialowych)
  m_sp = 0.5, 1.0, 2.0
  gamma = 1.0
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

try:
    from nbody.three_body_force_exact import (
        yukawa_overlap_exact, three_body_forces_exact
    )
    EXACT_AVAILABLE = True
    print("Dokladne sily 3-cialowe: DOSTEPNE")
except ImportError as e:
    print(f"UWAGA: {e}")
    EXACT_AVAILABLE = False

print("=" * 70)
print("EX17: Dynamika N-cialowa TGP z dokladnymi silami 3-cialowymi")
print("=" * 70)
print()

# ===========================================================================
# Sily 2-cialowe (Yukawa)
# ===========================================================================
def force_2body(r_i, r_j, C_i, C_j, m_sp):
    """Sila Yukawa na cialo i od ciala j."""
    dr = r_i - r_j
    d  = np.linalg.norm(dr)
    if d < 1e-10:
        return np.zeros(3)
    # V2 = -4*pi*Ci*Cj*exp(-m*d)/d
    # F_i = -dV2/dr_i = -4*pi*Ci*Cj * (dr/d) * (-m_sp - 1/d) * exp(-m*d)/d
    #      = 4*pi*Ci*Cj * exp(-m*d)/d * (m_sp + 1/d) * dr/d
    # Kierunek: dr/d (od j do i), wartość powinna byc ujemna (przyciaganie)
    # F_i = -dV2/dr_i = -V2/d * (m_sp*d + 1)/d * hat_r
    #      gdzie hat_r = dr/d
    coeff = -4.0*np.pi*C_i*C_j * np.exp(-m_sp*d) / d
    dcoeff_dd = coeff * (-m_sp - 1.0/d)  # d/dd [coeff]
    force = -dcoeff_dd * (dr/d)  # F = -dV/dr = -dV/dd * dr/|dr|
    # Poprawnie:
    # V2 = coeff * d^0 = -4pi*Ci*Cj*exp(-m*d)/d
    # dV2/dr_i = d/dr_i [V2(|r_i-r_j|)] = dV2/dd * d(|r_i-r_j|)/dr_i
    # d|r_i-r_j|/dr_i = (r_i-r_j)/|r_i-r_j| = dr/d
    # dV2/dd = d/dd[-4pi*Ci*Cj*exp(-m*d)/d]
    #        = -4pi*Ci*Cj * [(-m_sp)*exp(-m*d)/d + exp(-m*d)*(-1/d^2)]
    #        = -4pi*Ci*Cj * exp(-m*d) * (-m_sp/d - 1/d^2)
    #        = 4pi*Ci*Cj * exp(-m*d) * (m_sp/d + 1/d^2)
    dV2_dd = 4.0*np.pi*C_i*C_j * np.exp(-m_sp*d) * (m_sp/d + 1.0/d**2)
    F_i = -dV2_dd * (dr/d)  # F_i = -dV/dr_i
    return F_i

def potential_2body(d, C_i, C_j, m_sp):
    """V2(d) = -4*pi*Ci*Cj*exp(-m*d)/d."""
    return -4.0*np.pi*C_i*C_j*np.exp(-m_sp*d)/d

# ===========================================================================
# N-body ODE z dokladnymi silami 3-cialowymi
# ===========================================================================
def nbody_ode(t, state, N, masses, C_vals, m_sp, gamma, include_3body=True):
    """
    Rownania ruchu TGP dla N cial.
    state = [x1,y1,z1, vx1,vy1,vz1, x2,..., vz_N]
    """
    pos = state[:3*N].reshape(N, 3)
    vel = state[3*N:].reshape(N, 3)

    acc = np.zeros((N, 3))

    # Sily 2-cialowe
    for i in range(N):
        for j in range(i+1, N):
            F_ij = force_2body(pos[i], pos[j], C_vals[i], C_vals[j], m_sp)
            acc[i] += F_ij / masses[i]
            acc[j] -= F_ij / masses[j]  # Newton 3.

    # Sily 3-cialowe (dokladne lub pomijane)
    if include_3body and EXACT_AVAILABLE and N >= 3:
        forces_3b = three_body_forces_exact(pos, C_vals, m_sp, gamma, n_quad=30)
        for i in range(N):
            acc[i] += forces_3b[i] / masses[i]

    dstate = np.zeros_like(state)
    dstate[:3*N] = vel.flatten()
    dstate[3*N:] = acc.flatten()
    return dstate

def total_energy(pos, vel, N, masses, C_vals, m_sp, gamma, include_3body=True):
    """Calkowita energia (kinetyczna + potencjalna) ukladu."""
    # Kinetyczna
    T = sum(0.5*masses[i]*np.dot(vel[i], vel[i]) for i in range(N))

    # Potencjalna 2-cialowa
    V2 = 0.0
    for i in range(N):
        for j in range(i+1, N):
            d = np.linalg.norm(pos[i]-pos[j])
            V2 += potential_2body(d, C_vals[i], C_vals[j], m_sp)

    # Potencjalna 3-cialowa (Feynman)
    V3 = 0.0
    if include_3body and EXACT_AVAILABLE and N >= 3:
        for i in range(N):
            for j in range(i+1, N):
                for k in range(j+1, N):
                    d12 = np.linalg.norm(pos[i]-pos[j])
                    d13 = np.linalg.norm(pos[i]-pos[k])
                    d23 = np.linalg.norm(pos[j]-pos[k])
                    I_Y = yukawa_overlap_exact(d12, d13, d23, m_sp, n_quad=30)
                    V3 += -6.0*gamma*C_vals[i]*C_vals[j]*C_vals[k]*I_Y

    return T + V2 + V3, T, V2, V3

# ===========================================================================
# Test 1: Trojkat rownoboczny — stabilnosc
# ===========================================================================
print("=" * 70)
print("TEST 1: Trojkat rownoboczny — orbity z efektami 3-cialowymi")
print("=" * 70)
print()

N3 = 3
C0 = 0.20
m_sp_test = 1.0
gamma_test = 1.0
d_eq = 2.0   # odleglosc rownowagowa

# Masy (uproszczone: rowne 1)
masses3 = np.ones(3)
C_vals3 = np.full(3, C0)

# Pozycje rownoboczne
pos0_3 = np.array([
    [d_eq, 0.0, 0.0],
    [-d_eq/2.0,  d_eq*np.sqrt(3)/2.0, 0.0],
    [-d_eq/2.0, -d_eq*np.sqrt(3)/2.0, 0.0],
])
# Centrum masy: 0
pos0_3 -= pos0_3.mean(axis=0)

# Predkosci orbitalne (rotacja wokol centrum)
# v = sqrt(F_centripetalna * r_cm)
# Dla rownobocznego trojkata: r_cm = d/sqrt(3)
r_cm = d_eq / np.sqrt(3)

# Oblicz sile Yukawa na czialo 0 od cial 1 i 2
F_test = (force_2body(pos0_3[0], pos0_3[1], C0, C0, m_sp_test) +
          force_2body(pos0_3[0], pos0_3[2], C0, C0, m_sp_test))
F_radial = np.linalg.norm(F_test)
# 3 ciala kraza razem: F = m*v^2/r_cm => v = sqrt(F*r_cm/m)
v_orbit = np.sqrt(F_radial * r_cm / masses3[0])

# Predkosci prostopadle do ramienia
vel0_3 = np.zeros((3, 3))
for i in range(3):
    r_vec = pos0_3[i]
    r_perp = np.array([-r_vec[1], r_vec[0], 0.0])  # 90 deg CCW
    if np.linalg.norm(r_perp) > 1e-10:
        vel0_3[i] = v_orbit * r_perp / np.linalg.norm(r_perp)

print(f"  C = {C0}, m_sp = {m_sp_test}, d = {d_eq}")
print(f"  Predkosc orbitalna: v = {v_orbit:.4f}")
print()

# Energia poczatkowa
E0_3, T0, V2_0, V3_0 = total_energy(pos0_3, vel0_3, N3, masses3, C_vals3,
                                     m_sp_test, gamma_test)
print(f"  Energia poczatkowa: E = {E0_3:.6f}")
print(f"    T = {T0:.6f}, V2 = {V2_0:.6f}, V3 = {V3_0:.6f}")
print(f"    V3/V2 = {V3_0/V2_0:.4f}  (efekt 3-cialowy: {100*abs(V3_0/V2_0):.2f}%)")
print()

# Integracja z i bez sil 3-cialowych
state0_3 = np.concatenate([pos0_3.flatten(), vel0_3.flatten()])
t_max = 20.0
t_eval = np.linspace(0, t_max, 2000)

print("  Integracja z silami 2-cialowymi (Newton, 2B only)...")
t0 = time.time()
sol_2B = solve_ivp(
    lambda t, s: nbody_ode(t, s, N3, masses3, C_vals3, m_sp_test,
                            gamma_test, include_3body=False),
    [0, t_max], state0_3, t_eval=t_eval,
    method='RK45', rtol=1e-9, atol=1e-11
)
print(f"    Czas: {time.time()-t0:.2f}s, steps: {sol_2B.t.shape[0]}")

if EXACT_AVAILABLE:
    print("  Integracja z silami 3-cialowymi (TGP exact, 2B+3B)...")
    t0 = time.time()
    sol_3B = solve_ivp(
        lambda t, s: nbody_ode(t, s, N3, masses3, C_vals3, m_sp_test,
                                gamma_test, include_3body=True),
        [0, t_max], state0_3, t_eval=t_eval,
        method='RK45', rtol=1e-9, atol=1e-11
    )
    print(f"    Czas: {time.time()-t0:.2f}s, steps: {sol_3B.t.shape[0]}")
else:
    sol_3B = None

print()

# Sprawdz zachowanie energii
def check_energy_conservation(sol, N, masses, C_vals, m_sp, gamma, include_3body):
    energies = []
    for k in range(0, len(sol.t), 50):
        state = sol.y[:, k]
        pos = state[:3*N].reshape(N, 3)
        vel = state[3*N:].reshape(N, 3)
        E, _, _, _ = total_energy(pos, vel, N, masses, C_vals, m_sp, gamma, include_3body)
        energies.append(E)
    return np.array(energies)

E_2B_arr = check_energy_conservation(sol_2B, N3, masses3, C_vals3, m_sp_test, gamma_test, False)
dE_2B = np.max(np.abs(E_2B_arr - E_2B_arr[0])) / abs(E_2B_arr[0])
print(f"  Zmiana energii (2B only): max |dE/E| = {dE_2B:.2e}")

if sol_3B is not None:
    E_3B_arr = check_energy_conservation(sol_3B, N3, masses3, C_vals3, m_sp_test, gamma_test, True)
    dE_3B = np.max(np.abs(E_3B_arr - E_3B_arr[0])) / abs(E_3B_arr[0])
    print(f"  Zmiana energii (2B+3B): max |dE/E| = {dE_3B:.2e}")
print()

# ===========================================================================
# Test 2: Skanowanie C — kiedy 3B efekty sa wazne?
# ===========================================================================
print("=" * 70)
print("TEST 2: Skanowanie C — V3/V2 i odchylenie orbit")
print("=" * 70)
print()
print(f"  {'C':>8s}  {'V3/V2':>10s}  {'|dpos|(t=5)':>14s}  {'|dE/E|':>10s}")
print("-" * 48)

C_scan = [0.05, 0.10, 0.20, 0.30, 0.50] if EXACT_AVAILABLE else [0.05, 0.10, 0.20]
t_check = 5.0

for C_test in C_scan:
    C_v = np.full(3, C_test)
    pos_c = pos0_3.copy()
    F_c = (force_2body(pos_c[0], pos_c[1], C_test, C_test, m_sp_test) +
           force_2body(pos_c[0], pos_c[2], C_test, C_test, m_sp_test))
    F_r = np.linalg.norm(F_c)
    v_c = np.sqrt(max(F_r * r_cm / 1.0, 1e-10))

    vel_c = np.zeros((3, 3))
    for i in range(3):
        r_vec = pos_c[i]
        r_perp = np.array([-r_vec[1], r_vec[0], 0.0])
        if np.linalg.norm(r_perp) > 1e-10:
            vel_c[i] = v_c * r_perp / np.linalg.norm(r_perp)

    E0c, T0c, V2c, V3c = total_energy(pos_c, vel_c, 3, masses3, C_v,
                                       m_sp_test, gamma_test)
    ratio = abs(V3c/V2c) if abs(V2c) > 1e-15 else 0.0

    state_c = np.concatenate([pos_c.flatten(), vel_c.flatten()])
    t_eval_short = np.linspace(0, t_check, 500)

    # 2B only
    sol2 = solve_ivp(
        lambda t, s: nbody_ode(t, s, 3, masses3, C_v, m_sp_test,
                                gamma_test, include_3body=False),
        [0, t_check], state_c, t_eval=t_eval_short,
        method='RK45', rtol=1e-9, atol=1e-11
    )

    if EXACT_AVAILABLE:
        # 3B
        sol3 = solve_ivp(
            lambda t, s: nbody_ode(t, s, 3, masses3, C_v, m_sp_test,
                                    gamma_test, include_3body=True),
            [0, t_check], state_c, t_eval=t_eval_short,
            method='RK45', rtol=1e-9, atol=1e-11
        )
        pos2_final = sol2.y[:9, -1].reshape(3, 3)
        pos3_final = sol3.y[:9, -1].reshape(3, 3)
        dpos = np.mean([np.linalg.norm(pos3_final[i]-pos2_final[i]) for i in range(3)])

        vel3_final = sol3.y[9:, -1].reshape(3, 3)
        E_final, _, _, _ = total_energy(pos3_final, vel3_final, 3, masses3, C_v,
                                         m_sp_test, gamma_test)
        dE = abs(E_final - E0c) / (abs(E0c) + 1e-20)
        print(f"  {C_test:>8.3f}  {ratio:>10.4f}  {dpos:>14.4e}  {dE:>10.4e}")
    else:
        print(f"  {C_test:>8.3f}  {ratio:>10.4f}  {'(bez 3B)':>14s}  {'---':>10s}")

print()

# ===========================================================================
# Test 3: Uklad 4-cialowy
# ===========================================================================
if EXACT_AVAILABLE:
    print("=" * 70)
    print("TEST 3: Uklad 4-cialowy — suma trojek 3-cialowych")
    print("=" * 70)
    print()

    N4 = 4
    C4 = 0.25
    C_vals4 = np.full(N4, C4)
    masses4 = np.ones(N4)

    # Czworobok kwadratowy
    d4 = 2.0
    pos4 = np.array([
        [ d4,  d4, 0.0],
        [-d4,  d4, 0.0],
        [-d4, -d4, 0.0],
        [ d4, -d4, 0.0],
    ])

    # Predkosci orbitalne
    r4 = d4 * np.sqrt(2.0)
    F4_total = np.zeros(3)
    for j in [1, 2, 3]:
        F4_total += force_2body(pos4[0], pos4[j], C4, C4, m_sp_test)
    v4 = np.sqrt(np.linalg.norm(F4_total) * r4 / masses4[0])

    vel4 = np.zeros((N4, 3))
    for i in range(N4):
        r_vec = pos4[i]
        r_perp = np.array([-r_vec[1], r_vec[0], 0.0])
        if np.linalg.norm(r_perp) > 1e-10:
            vel4[i] = v4 * r_perp / np.linalg.norm(r_perp)

    E0_4, T4, V2_4, V3_4 = total_energy(pos4, vel4, N4, masses4, C_vals4,
                                          m_sp_test, gamma_test)
    print(f"  N=4, C={C4}, m_sp={m_sp_test}")
    print(f"  E = {E0_4:.6f}: T={T4:.4f}, V2={V2_4:.4f}, V3={V3_4:.4f}")
    print(f"  V3/V2 = {V3_4/V2_4:.4f}  ({100*abs(V3_4/V2_4):.2f}%)")
    print()

    # Krótka integracja
    state4 = np.concatenate([pos4.flatten(), vel4.flatten()])
    t4 = np.linspace(0, 5.0, 500)
    sol4_2B = solve_ivp(
        lambda t, s: nbody_ode(t, s, N4, masses4, C_vals4, m_sp_test,
                                gamma_test, include_3body=False),
        [0, 5.0], state4, t_eval=t4, method='RK45', rtol=1e-9, atol=1e-11
    )
    sol4_3B = solve_ivp(
        lambda t, s: nbody_ode(t, s, N4, masses4, C_vals4, m_sp_test,
                                gamma_test, include_3body=True),
        [0, 5.0], state4, t_eval=t4, method='RK45', rtol=1e-9, atol=1e-11
    )

    pos4_2B_f = sol4_2B.y[:12, -1].reshape(4, 3)
    pos4_3B_f = sol4_3B.y[:12, -1].reshape(4, 3)
    dpos4 = np.mean([np.linalg.norm(pos4_3B_f[i]-pos4_2B_f[i]) for i in range(4)])
    print(f"  Odchylenie 2B vs 3B po t=5: |dpos| = {dpos4:.4e}")
    print()

# ===========================================================================
# Podsumowanie
# ===========================================================================
print("=" * 70)
print("PODSUMOWANIE: Kiedy efekty 3-cialowe TGP sa wazne?")
print("=" * 70)
print()
print("  KRYTERIUM: V3/V2 > threshold (np. 1%)")
print()
print("  Dla rownobocznego trojkata d=2, m_sp=1:")
print("  V3/V2 ~ (4*pi/3) * C  [dla t = m_sp*d ~ 2, Yukawa supresja]")
print()
print("  Efekty 3-cialowe sa wazne dla:")
print("  1. Duzych C (silny ladunek TGP) -- nierealistyczne dla zwyklej materii")
print("  2. Malych odleglosci d ~ 1/m_sp (t = m_sp*d ~ 1)")
print("  3. Wielu cial jednoczesnie (suma trojek ~ N^2/2 razy)")
print()
print("  Dla C = m/(2*sqrt(pi)) (z warunku Newtona):")
print("  C_proton ~ 2e-20 => V3/V2 ~ 1e-20 -- nieobserwowalne")
print("  C_planck ~ 0.28  => V3/V2 ~ 1.2   -- DOMINUJACE!")
print()
print("  WNIOSEK: Efekty 3-cialowe TGP sa obserwowalne TYLKO dla:")
print("  - Hipotetycznych objetkow Planckowskich (m ~ m_Planck)")
print("  - Lub dla nowych ciezkich cial z duzym ladunkiem TGP")
print()
print("  W granicy m_sp -> inf (sila zanika na skali Plancka):")
print("  Efekty 3-cialowe sa supresowane eksponencjalnie.")
print("  TGP staje sie praktycznie rownowazne z grawitacja Newtonowska.")

# ===========================================================================
# Wykresy
# ===========================================================================
fig = plt.figure(figsize=(15, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.3)

# 1. Orbity (2B vs 3B) dla trojkata rownobocznego
ax1 = fig.add_subplot(gs[0, 0])
colors_bodies = ['b', 'r', 'g']
for i, col in enumerate(colors_bodies):
    x_2B = sol_2B.y[3*i]
    y_2B = sol_2B.y[3*i+1]
    ax1.plot(x_2B, y_2B, color=col, lw=1, ls='--', alpha=0.7,
             label=f'#{i+1} (2B)' if i == 0 else None)
    if sol_3B is not None:
        x_3B = sol_3B.y[3*i]
        y_3B = sol_3B.y[3*i+1]
        ax1.plot(x_3B, y_3B, color=col, lw=1.5, ls='-',
                 label=f'#{i+1} (2B+3B)' if i == 0 else None)
ax1.set_xlabel('x', fontsize=11)
ax1.set_ylabel('y', fontsize=11)
ax1.set_title(f'Orbity: C={C0}, m_sp={m_sp_test}', fontsize=10)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_aspect('equal')

# 2. Roznica pozycji 2B vs 3B
ax2 = fig.add_subplot(gs[0, 1])
if sol_3B is not None:
    dpos_t = []
    for k in range(len(sol_2B.t)):
        dp = 0.0
        for i in range(3):
            p2 = sol_2B.y[3*i:3*i+3, k]
            p3 = sol_3B.y[3*i:3*i+3, k]
            dp += np.linalg.norm(p3 - p2)**2
        dpos_t.append(np.sqrt(dp/3))
    ax2.semilogy(sol_2B.t, np.maximum(dpos_t, 1e-16), 'b-', lw=2)
    ax2.set_xlabel('t', fontsize=11)
    ax2.set_ylabel('|Δpos| (2B vs 3B)', fontsize=11)
    ax2.set_title('Odchylenie orbit 3-cialowych', fontsize=10)
    ax2.grid(True, alpha=0.3)

# 3. Zachowanie energii
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(range(len(E_2B_arr)), E_2B_arr - E_2B_arr[0],
         'b-o', ms=4, lw=1.5, label='2B only')
if sol_3B is not None:
    ax3.plot(range(len(E_3B_arr)), E_3B_arr - E_3B_arr[0],
             'r-s', ms=4, lw=1.5, label='2B+3B')
ax3.axhline(0, color='k', lw=0.5, ls='--')
ax3.set_xlabel('krok', fontsize=11)
ax3.set_ylabel('ΔE = E(t) - E(0)', fontsize=11)
ax3.set_title('Zachowanie energii', fontsize=10)
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# 4. V3/V2 vs C
ax4 = fig.add_subplot(gs[1, 0])
C_range = np.logspace(-3, 0, 100)
m_sp_range = [0.5, 1.0, 2.0]
d_range = 2.0
for m_sp_v, col, lbl in zip(m_sp_range,
                              ['b', 'g', 'r'],
                              [f'm_sp={v}' for v in m_sp_range]):
    ratio_range = []
    for C_v in C_range:
        if EXACT_AVAILABLE:
            I_Y = yukawa_overlap_exact(d_range, d_range, d_range, m_sp_v, n_quad=30)
            V3_v = -6.0*gamma_test*C_v**3*I_Y
        else:
            # Przyblizenie
            t_v = m_sp_v * d_range
            I_Y_approx = np.exp(-np.sqrt(3)*t_v) / (np.sqrt(3)*t_v)
            V3_v = -6.0*C_v**3*I_Y_approx
        V2_v = potential_2body(d_range, C_v, C_v, m_sp_v) * 3.0  # 3 pary
        ratio_range.append(abs(V3_v/V2_v) if abs(V2_v) > 1e-20 else 0)
    ax4.loglog(C_range, ratio_range, color=col, lw=2, label=lbl)

ax4.axhline(0.01, color='k', ls='--', lw=1, alpha=0.5, label='1% threshold')
ax4.axhline(1.00, color='k', ls=':', lw=1, alpha=0.5, label='100% (dominujace)')
ax4.axvline(0.28, color='orange', ls=':', lw=1.5, alpha=0.7, label='C = m_Pl/2√π')
ax4.set_xlabel('C (ladunek TGP)', fontsize=11)
ax4.set_ylabel('|V3/V2|', fontsize=11)
ax4.set_title(r'$|V_3/V_2|$ vs C (trojkat, d=2)', fontsize=10)
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3, which='both')

# 5. Trajektorie 4 cial
if EXACT_AVAILABLE and 'sol4_2B' in dir():
    ax5 = fig.add_subplot(gs[1, 1])
    c4 = ['b', 'r', 'g', 'orange']
    for i, col in enumerate(c4):
        ax5.plot(sol4_2B.y[3*i], sol4_2B.y[3*i+1], color=col, lw=1, ls='--', alpha=0.5)
        ax5.plot(sol4_3B.y[3*i], sol4_3B.y[3*i+1], color=col, lw=1.5)
        ax5.scatter([pos4[i,0]], [pos4[i,1]], c=col, s=60, zorder=5)
    ax5.set_xlabel('x', fontsize=11)
    ax5.set_ylabel('y', fontsize=11)
    ax5.set_title(f'Orbity 4 cial: C={C4}', fontsize=10)
    ax5.grid(True, alpha=0.3)
    ax5.set_aspect('equal')

# 6. Energia kinetyczna vs potencjalna (wirialni)
ax6 = fig.add_subplot(gs[1, 2])
if EXACT_AVAILABLE and sol_3B is not None:
    T_arr = []
    V_arr = []
    for k in range(len(sol_3B.t)):
        state = sol_3B.y[:, k]
        pos_k = state[:9].reshape(3, 3)
        vel_k = state[9:].reshape(3, 3)
        E_k, T_k, V2_k, V3_k = total_energy(pos_k, vel_k, 3, masses3, C_vals3,
                                              m_sp_test, gamma_test)
        T_arr.append(T_k)
        V_arr.append(V2_k + V3_k)
    ax6.plot(sol_3B.t, T_arr, 'b-', lw=1.5, label='T (kinetyczna)')
    ax6.plot(sol_3B.t, V_arr, 'r-', lw=1.5, label='V2+V3 (potencjalna)')
    ax6.plot(sol_3B.t, np.array(T_arr)+np.array(V_arr), 'k-', lw=2, alpha=0.7, label='E_total')
    ax6.set_xlabel('t', fontsize=11)
    ax6.set_ylabel('Energia', fontsize=11)
    ax6.set_title('Wirial: T vs V (2B+3B)', fontsize=10)
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)

plt.suptitle('TGP: Dynamika N-cialowa z dokladnymi silami 3-cialowymi',
             fontsize=13, fontweight='bold')
out = os.path.join(os.path.dirname(__file__), 'ex17_nbody_exact.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nWykres: {out}")
