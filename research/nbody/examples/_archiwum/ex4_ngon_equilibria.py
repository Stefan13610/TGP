"""
ex4_ngon_equilibria.py
======================
TGP: równowagi regularnych N-kątów dla N = 3, 4, 5, 6, 7.

Teoria:
  Dla regularnego N-kąta wpisanego w okrąg o promieniu R,
  odległość od centrum: r = R.
  Odległości między ciałami: d_ij = 2R·sin(π·|i-j|/N).

  Warunek równowagi: ΣᵢF_net(i) = 0 wzdłuż osi radialnej.
  Dla równych mas C: numeryczne znalezienie R_eq.

Uruchomienie:
    python examples/ex4_ngon_equilibria.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from nbody.pairwise import force_2body
from nbody.dynamics_v2 import (leapfrog_integrate,
                                analytical_forces_tgp_pairwise,
                                potential_tgp)

BETA  = 1.0
GAMMA = 1.0
C     = 0.10


def ngon_positions_3d(N, R):
    """Pozycje regularnego N-kąta wpisanego w okrąg promienia R (3D, z=0)."""
    angles = 2.0 * np.pi * np.arange(N) / N
    xs = R * np.cos(angles)
    ys = R * np.sin(angles)
    zs = np.zeros(N)
    return np.column_stack([xs, ys, zs])


def ngon_positions_2d(N, R):
    """Pozycje 2D dla wykresów."""
    angles = 2.0 * np.pi * np.arange(N) / N
    return np.column_stack([R * np.cos(angles), R * np.sin(angles)])


def radial_force_on_body0(R, N, beta_loc=None, gamma_loc=None):
    """Siła radialna (ku zewnątrz) na ciele 0 w N-kącie promienia R."""
    if beta_loc is None:
        beta_loc = BETA
    if gamma_loc is None:
        gamma_loc = GAMMA
    pos2d = ngon_positions_2d(N, R)
    f_total_x = 0.0
    for j in range(1, N):
        d_vec = pos2d[0] - pos2d[j]
        d = np.linalg.norm(d_vec)
        F_scalar = force_2body(d, C, C, beta_loc, gamma_loc)
        # Składowa radialna (jednostkowy wektor od centrum do ciała 0 to [1,0])
        f_total_x += F_scalar * d_vec[0] / d
    return f_total_x


def find_ngon_equilibrium(N, R_min=0.1, R_max=30.0, n_scan=2000):
    """Szuka promienia równowagi dla N-kąta (stabilna, zewnętrzna).

    Skanuje przedział w poszukiwaniu zmian znaku siły radialnej,
    a następnie używa brentq do dokładnego wyznaczenia zera.
    Zwraca ostatnie (największe R) zero — stabilną równowagę.
    """
    try:
        R_scan = np.linspace(R_min, R_max, n_scan)
        F_scan = np.array([radial_force_on_body0(R, N) for R in R_scan])

        # Znajdź wszystkie zmiany znaku
        sign_changes = []
        for k in range(len(R_scan) - 1):
            if F_scan[k] * F_scan[k + 1] < 0:
                sign_changes.append((R_scan[k], R_scan[k + 1]))

        if not sign_changes:
            return None

        # Zwróć ostatnie (stabilne, zewnętrzne) zero
        Ra, Rb = sign_changes[-1]
        R_eq = brentq(lambda R: radial_force_on_body0(R, N),
                      Ra, Rb, xtol=1e-8)
        return R_eq
    except Exception:
        return None


def simulate_ngon(N, R_eq, t_max=24.0, amplitude=0.03):
    """Symuluje N-kąt z małą perturbacją oddechową (3D, z=0)."""
    pos0 = ngon_positions_3d(N, R_eq * (1 + amplitude))
    vel0 = np.zeros_like(pos0)
    masses   = np.full(N, C)

    def force_func(pos, C_vals):
        F = analytical_forces_tgp_pairwise(pos, C_vals, BETA, GAMMA)
        return F / C_vals[:, None]

    def potential_func(pos, C_vals):
        return potential_tgp(pos, C_vals, BETA, GAMMA)

    result = leapfrog_integrate(
        pos0, vel0, masses,
        force_func, potential_func,
        t_span=(0.0, t_max),
        dt=4e-4,
        save_every=300
    )
    return result['positions'], result['t'], result['energy']


if __name__ == '__main__':
    print("="*60)
    print("EX4: Równowagi regularnych N-kątów TGP")
    print("="*60)
    print(f"Parametry: β={BETA}, γ={GAMMA}, C={C}  (β/C = {BETA/C:.1f})")

    results = {}
    fig_eq, ax_eq = plt.subplots(figsize=(8, 5))

    for N in [3, 4, 5, 6, 7]:
        R_eq = find_ngon_equilibrium(N)
        side = 2.0 * R_eq * np.sin(np.pi / N) if R_eq else None
        if R_eq is not None:
            print(f"\nN={N}: R_eq = {R_eq:.5f},  bok d = {side:.5f}")
        else:
            print(f"\nN={N}: brak równowagi")
        results[N] = R_eq

        if R_eq is not None:
            pos_eq = ngon_positions_2d(N, R_eq)
            ax_eq.plot(pos_eq[:, 0], pos_eq[:, 1], 'o-',
                       label=f'N={N}, R={R_eq:.3f}', markersize=8)
            # Zamknij wielokąt
            pos_closed = np.vstack([pos_eq, pos_eq[0]])
            ax_eq.plot(pos_closed[:, 0], pos_closed[:, 1], '-', alpha=0.4)

    ax_eq.set_aspect('equal')
    ax_eq.legend(fontsize=9)
    ax_eq.set_title(f'Równowagi N-kątów TGP (β={BETA}, C={C})', fontsize=12)
    ax_eq.grid(True, alpha=0.3)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out_eq = os.path.join(plots_dir, 'ex4_ngon_configs.png')
    plt.tight_layout()
    plt.savefig(out_eq, dpi=120)
    plt.close()
    print(f"\nZapisano konfiguracje: {os.path.abspath(out_eq)}")

    # Symulacja N=4 (kwadrat)
    if results.get(4):
        print("\nSymulacja kwadratu (N=4)...")
        pos, t_arr, E_arr = simulate_ngon(4, results[4])
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Trajektorie (x-y)
        colors = ['b', 'r', 'g', 'm']
        for i, col in enumerate(colors):
            axes[0].plot(pos[:, i, 0], pos[:, i, 1], col+'-', linewidth=0.8, alpha=0.8)
            axes[0].plot(pos[0, i, 0], pos[0, i, 1], col+'o', markersize=6)
        axes[0].set_aspect('equal', adjustable='datalim')
        axes[0].set_title('Trajektorie N=4 (kwadrat)', fontsize=11)
        axes[0].grid(True, alpha=0.3)

        # Energia
        dE = np.abs((E_arr - E_arr[0]) / (np.abs(E_arr[0]) + 1e-20))
        axes[1].semilogy(t_arr, dE + 1e-15, 'g-')
        axes[1].set_title('Zachowanie energii', fontsize=11)
        axes[1].set_xlabel('Czas', fontsize=10)
        axes[1].set_ylabel(r'$|\Delta E/E_0|$', fontsize=10)
        axes[1].grid(True, alpha=0.3)

        out_sim = os.path.join(plots_dir, 'ex4_square_dynamics.png')
        plt.tight_layout()
        plt.savefig(out_sim, dpi=120)
        plt.close()
        print(f"Zapisano symulację kwadratu: {os.path.abspath(out_sim)}")

    # Wykres R_eq vs N
    valid = [(N, R) for N, R in results.items() if R is not None]
    if len(valid) >= 2:
        Ns, Rs = zip(*valid)
        fig2, ax2 = plt.subplots(figsize=(7, 4))
        ax2.plot(Ns, Rs, 'bo-', markersize=8, linewidth=2)
        ax2.set_xlabel('Liczba ciał N', fontsize=12)
        ax2.set_ylabel('Promień równowagi $R_{eq}$', fontsize=12)
        ax2.set_title(f'TGP N-kąty: promień vs N (β={BETA}, C={C})', fontsize=12)
        ax2.grid(True, alpha=0.3)
        out2 = os.path.join(plots_dir, 'ex4_R_vs_N.png')
        plt.tight_layout()
        plt.savefig(out2, dpi=120)
        plt.close()
        print(f"Zapisano R vs N: {os.path.abspath(out2)}")

    print("\nGotowe.")
