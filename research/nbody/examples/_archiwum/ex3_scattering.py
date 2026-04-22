"""
ex3_scattering.py
=================
TGP: zderzenie trzech ciał — różne parametry uderzenia.

Scenariusz: dwa ciała spoczywają w równowadze, trzecie nadlatuje.
Badamy: przechwycenie, przelot, rozpad w zależności od prędkości i parametru.

Uruchomienie:
    python examples/ex3_scattering.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from nbody.pairwise import force_zeros_2body
from nbody.dynamics_v2 import (leapfrog_integrate,
                                analytical_forces_tgp_pairwise,
                                potential_tgp)

BETA    = 1.0
GAMMA   = 1.0
C       = 0.15
DT      = 5e-4
T_MAX   = 60.0
SAVE_EVERY = 400


def make_force_func():
    def force_func(pos, C_vals):
        F = analytical_forces_tgp_pairwise(pos, C_vals, BETA, GAMMA)
        return F / C_vals[:, None]
    return force_func

def make_potential_func():
    def potential_func(pos, C_vals):
        return potential_tgp(pos, C_vals, BETA, GAMMA)
    return potential_func


def initial_conditions_scattering(v_inf, impact_b, d_eq=None):
    """Dwa ciała w równowadze, trzecie nadlatuje z nieskończoności.
    Pozycje są trójwymiarowe (z=0), zgodnie z konwencją dynamics_v2.
    """
    if d_eq is None:
        result = force_zeros_2body(C, BETA, GAMMA)
        d_eq = result[1] if result and len(result) == 2 else 4.0

    # Ciała 0 i 1 w równowadze (wzdłuż osi x)
    pos = np.array([
        [-d_eq/2,  0.0,      0.0],
        [ d_eq/2,  0.0,      0.0],
        [ -10.0,   impact_b, 0.0]   # ciało 3 nadlatuje z lewej
    ])
    vel = np.array([
        [0.0,   0.0, 0.0],
        [0.0,   0.0, 0.0],
        [v_inf, 0.0, 0.0]
    ])
    return pos, vel, d_eq


def classify_outcome(pos_final, d_eq, threshold=3.0):
    """Klasyfikuje wynik zderzenia: przechwycenie / przelot."""
    d01 = np.linalg.norm(pos_final[-1, 0] - pos_final[-1, 1])
    d02 = np.linalg.norm(pos_final[-1, 0] - pos_final[-1, 2])
    d12 = np.linalg.norm(pos_final[-1, 1] - pos_final[-1, 2])
    max_d = max(d01, d02, d12)
    if max_d < threshold * d_eq:
        return "przechwycenie"
    else:
        return "przelot"


def run_scenario(v_inf, impact_b, label=""):
    """Uruchamia jeden scenariusz zderzenia."""
    pos0, vel0, d_eq = initial_conditions_scattering(v_inf, impact_b)
    masses = np.array([C, C, C])

    result = leapfrog_integrate(
        pos0, vel0, masses,
        make_force_func(), make_potential_func(),
        t_span=(0.0, T_MAX),
        dt=DT,
        save_every=SAVE_EVERY
    )

    pos   = result['positions']
    t_arr = result['t']
    E_arr = result['energy']

    outcome = classify_outcome(pos, d_eq)
    return pos, t_arr, E_arr, outcome, d_eq


def plot_trajectories(scenarios):
    """Rysuje trajektorie dla kilku scenariuszy (wyświetla tylko x-y)."""
    n = len(scenarios)
    fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
    if n == 1:
        axes = [axes]

    colors = ['b', 'r', 'g']
    for ax, (label, pos, outcome, d_eq) in zip(axes, scenarios):
        for i, col in zip(range(3), colors):
            ax.plot(pos[:, i, 0], pos[:, i, 1], col + '-', linewidth=0.8,
                    alpha=0.8, label=f'Ciało {i+1}')
            ax.plot(pos[0,  i, 0], pos[0,  i, 1], col + 'o', markersize=6)
            ax.plot(pos[-1, i, 0], pos[-1, i, 1], col + 's', markersize=6)
        ax.set_title(f'{label}\n-> {outcome}', fontsize=11)
        ax.set_xlabel('x', fontsize=10)
        ax.set_ylabel('y', fontsize=10)
        ax.legend(fontsize=8)
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(True, alpha=0.3)

    plt.suptitle('TGP 3-body: zderzenia przy różnych prędkościach i parametrach', fontsize=13)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex3_scattering.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Zapisano: {os.path.abspath(out)}")


if __name__ == '__main__':
    print("="*60)
    print("EX3: Zderzenia 3-ciał TGP")
    print("="*60)

    cases = [
        (0.3, 0.0,  "v=0.3, b=0.0 (czolowe)"),
        (0.8, 0.5,  "v=0.8, b=0.5"),
        (1.5, 1.0,  "v=1.5, b=1.0 (szybkie)"),
    ]

    scenarios = []
    for v_inf, b, lbl in cases:
        print(f"\nScenariusz: {lbl}")
        pos, t_arr, E_arr, outcome, d_eq = run_scenario(v_inf, b, lbl)
        dE = abs((E_arr[-1] - E_arr[0]) / (abs(E_arr[0]) + 1e-20))
        print(f"  Wynik: {outcome}")
        print(f"  Zachowanie energii: |ΔE/E₀| = {dE:.2e}")
        scenarios.append((lbl, pos, outcome, d_eq))

    print("\nRysowanie trajektorii...")
    plot_trajectories(scenarios)
    print("\nGotowe.")
