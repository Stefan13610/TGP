"""
ex5_tgp_vs_newton.py
====================
TGP vs grawitacja Newtonowska: porównanie dynamiki.

Kluczowe różnice:
  TGP:    V(d) = -A/d + B/d² - C_conf/d³  ->  istnieje równowaga
  Newton: V(d) = -G·m₁·m₂/d              ->  brak równowagi (twierdzenie Earnshawa)

Pokazujemy:
  1. Krzywe V(d) dla obu teorii
  2. Dynamikę 3-ciał: trajektorie i energia

Uwaga implementacyjna:
  dynamics_v2.leapfrog_integrate przyjmuje force_func i potential_func
  jako argumenty, dzięki czemu TGP i Newton działają przez ten sam integrator.
  Użyte funkcje z dynamics_v2: analytical_forces_tgp_pairwise, potential_tgp,
  forces_newton, potential_newton.

Uruchomienie:
    python examples/ex5_tgp_vs_newton.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from nbody.pairwise import V_eff_total, force_zeros_2body
from nbody.dynamics_v2 import (leapfrog_integrate,
                                analytical_forces_tgp_pairwise,
                                potential_tgp,
                                forces_newton,
                                potential_newton)

BETA  = 1.0
GAMMA = 1.0
C     = 0.15
# Newtonowska stała dobrana tak, by siła dalekiego zasięgu była porównywalna z TGP.
# TGP: F_grad ~ -4π·C²/d² = -4π·0.0225/d²  →  G_N = 4π·C²
G_N   = 4.0 * np.pi * C**2


# ── Pomocnicze opakowania force/potential ─────────────────────────────

def tgp_force_func(pos, C_vals):
    """Przyspieszenie TGP (siła/masa). Wymaga 3D (n,3) ndarray."""
    F = analytical_forces_tgp_pairwise(pos, C_vals, BETA, GAMMA)
    return F / C_vals[:, None]

def tgp_potential_func(pos, C_vals):
    return potential_tgp(pos, C_vals, BETA, GAMMA)

def newton_force_func(pos, C_vals):
    """Przyspieszenie Newtonowskie (siła/masa). forces_newton zwraca siły."""
    F = forces_newton(pos, C_vals, G=G_N)
    return F / C_vals[:, None]

def newton_potential_func(pos, C_vals):
    return potential_newton(pos, C_vals, G=G_N)


# ── Porównanie potencjałów ────────────────────────────────────────────

def compare_potentials():
    """Porównuje potencjały TGP i Newton."""
    d_vals = np.linspace(0.2, 8.0, 1000)
    V_tgp  = V_eff_total(d_vals, C, C, BETA, GAMMA)
    V_newt = -G_N / d_vals

    fig, ax = plt.subplots(figsize=(9, 5))
    mask_tgp  = (V_tgp  > -8) & (V_tgp  < 3)
    mask_newt = (V_newt > -8) & (V_newt < 3)
    ax.plot(d_vals[mask_tgp],  V_tgp[mask_tgp],   'b-',  linewidth=2, label='TGP $V_2(d)$')
    ax.plot(d_vals[mask_newt], V_newt[mask_newt], 'r--', linewidth=2, label='Newton $-G/d$')
    ax.axhline(0, color='k', linewidth=0.5)

    # Oznacz równowagę TGP
    result = force_zeros_2body(C, BETA, GAMMA)
    if result and len(result) == 2:
        d_in, d_out = result
        ax.axvline(d_out, color='g', linestyle=':', alpha=0.8,
                   label=f'TGP równowaga $d_+=$ {d_out:.2f}')
        ax.plot(d_out, V_eff_total(d_out, C, C, BETA, GAMMA), 'gs', markersize=10)

    ax.set_xlabel('Odległość $d$', fontsize=12)
    ax.set_ylabel('Energia potencjalna $V(d)$', fontsize=12)
    ax.set_title('TGP vs Newton: potencjały dwuciałowe', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_ylim(-6, 2)
    ax.grid(True, alpha=0.3)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex5_potentials.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Zapisano: {os.path.abspath(out)}")


# ── Dynamika 3-ciał ───────────────────────────────────────────────────

def simulate_3body(force_func, potential_func, pos0, vel0, masses,
                   t_max=40.0, dt=5e-4, save_every=300):
    """Uruchamia symulację 3-ciał za pomocą leapfrog."""
    result = leapfrog_integrate(
        pos0.copy(), vel0.copy(), masses,
        force_func, potential_func,
        t_span=(0.0, t_max),
        dt=dt,
        save_every=save_every
    )
    return result['positions'], result['t'], result['energy']


def setup_equilateral(d_eq, amplitude=0.08):
    """Tworzy równoboczny trójkąt o boku d_eq z perturbacją."""
    h = d_eq * np.sqrt(3) / 2.0
    pos0 = np.array([
        [-d_eq/2,  -h/3,   0.0],
        [ d_eq/2,  -h/3,   0.0],
        [ 0.0,      2*h/3, 0.0]
    ])
    pos0 = pos0 * (1.0 + amplitude)
    vel0 = np.zeros((3, 3))
    masses = np.array([C, C, C])
    return pos0, vel0, masses


if __name__ == '__main__':
    print("="*60)
    print("EX5: TGP vs Newton — porównanie dynamiki")
    print("="*60)

    print("\n[1] Krzywe potencjału:")
    compare_potentials()

    result_eq = force_zeros_2body(C, BETA, GAMMA)
    if result_eq and len(result_eq) == 2:
        d_eq = result_eq[1]
        print(f"\n[2] Dynamika 3-ciał (d_eq={d_eq:.3f}, perturbacja 8%):")
        pos0, vel0, masses = setup_equilateral(d_eq, amplitude=0.08)

        print("  Symulacja TGP...")
        pos_tgp, t_arr, E_tgp = simulate_3body(
            tgp_force_func, tgp_potential_func, pos0, vel0, masses
        )
        print("  Symulacja Newton...")
        pos_newt, _, E_newt = simulate_3body(
            newton_force_func, newton_potential_func, pos0, vel0, masses
        )

        # Wykresy
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        # Trajektoria TGP (x-y)
        cols = ['b', 'r', 'g']
        for i, col in enumerate(cols):
            axes[0, 0].plot(pos_tgp[:, i, 0], pos_tgp[:, i, 1],
                            col+'-', linewidth=0.8, alpha=0.9, label=f'Ciało {i+1}')
        axes[0, 0].set_aspect('equal', adjustable='datalim')
        axes[0, 0].set_title('TGP: trajektorie 3-ciał', fontsize=11)
        axes[0, 0].legend(fontsize=8)
        axes[0, 0].grid(True, alpha=0.3)

        # Trajektoria Newton (x-y)
        for i, col in enumerate(cols):
            axes[0, 1].plot(pos_newt[:, i, 0], pos_newt[:, i, 1],
                            col+'-', linewidth=0.8, alpha=0.9, label=f'Ciało {i+1}')
        axes[0, 1].set_aspect('equal', adjustable='datalim')
        axes[0, 1].set_title('Newton: trajektorie 3-ciał', fontsize=11)
        axes[0, 1].legend(fontsize=8)
        axes[0, 1].grid(True, alpha=0.3)

        # Odległość 0-1 vs czas
        d_tgp  = np.linalg.norm(pos_tgp[:, 0, :]  - pos_tgp[:, 1, :],  axis=1)
        d_newt = np.linalg.norm(pos_newt[:, 0, :] - pos_newt[:, 1, :], axis=1)
        axes[1, 0].plot(t_arr, d_tgp,  'b-', linewidth=1.2, label='TGP')
        axes[1, 0].plot(t_arr, d_newt, 'r-', linewidth=1.2, label='Newton', alpha=0.8)
        axes[1, 0].axhline(d_eq, color='g', linestyle='--', alpha=0.7,
                            label=f'$d_{{eq}}={d_eq:.2f}$')
        axes[1, 0].set_xlabel('Czas', fontsize=10)
        axes[1, 0].set_ylabel('$d_{01}$', fontsize=10)
        axes[1, 0].set_title('Odległość ciał 1-2 vs czas', fontsize=11)
        axes[1, 0].legend(fontsize=9)
        axes[1, 0].grid(True, alpha=0.3)

        # Zachowanie energii
        dE_tgp  = np.abs((E_tgp  - E_tgp[0])  / (np.abs(E_tgp[0])  + 1e-20))
        dE_newt = np.abs((E_newt - E_newt[0]) / (np.abs(E_newt[0]) + 1e-20))
        axes[1, 1].semilogy(t_arr, dE_tgp  + 1e-15, 'b-', label='TGP', linewidth=1.2)
        axes[1, 1].semilogy(t_arr, dE_newt + 1e-15, 'r-', label='Newton', linewidth=1.2)
        axes[1, 1].set_xlabel('Czas', fontsize=10)
        axes[1, 1].set_ylabel(r'$|\Delta E/E_0|$', fontsize=10)
        axes[1, 1].set_title('Zachowanie energii', fontsize=11)
        axes[1, 1].legend(fontsize=9)
        axes[1, 1].grid(True, alpha=0.3)

        plt.suptitle('TGP vs Newton: dynamika układu 3-ciał', fontsize=14)
        plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        out = os.path.join(plots_dir, 'ex5_tgp_vs_newton.png')
        plt.tight_layout()
        plt.savefig(out, dpi=120)
        plt.close()
        print(f"Zapisano: {os.path.abspath(out)}")
    else:
        print("Brak równowagi TGP dla podanych parametrów.")

    print("\nGotowe.")
