"""
ex2_equilateral_3body.py
========================
TGP: równoboczny układ 3-ciał — równowaga, stabilność, oscylacje.

Fizyka:
  - Potencjał 2-ciałowy V₂(d): gradient(-1/d) + beta(+1/d²) + gamma(-1/d³)
  - Potencjał 3-ciałowy V₃: wkład z Φ⁴, ≈-6γC³·f_△(md)
  - Równanie równowagi: d² - 4β·d + 18β·C = 0  (dla γ=β)
  - Próg istnienia: β/C > 9/2 = 4.5

Uruchomienie:
    python examples/ex2_equilateral_3body.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from nbody.pairwise import V_eff_total, force_2body, force_zeros_2body
from nbody.dynamics_v2 import (leapfrog_integrate,
                                analytical_forces_tgp_pairwise,
                                potential_tgp)

# ── Parametry fizyczne ─────────────────────────────────────────────────

BETA  = 1.0    # sprzężenie kwadratowe (= masa ekranowania² w warunku próżni)
GAMMA = BETA   # warunek próżni: β = γ
C     = 0.15   # masa/ładunek TGP każdego ciała (β/C = 6.67 > 4.5 ✓)

# ── Równowaga ─────────────────────────────────────────────────────────

def find_equilibria():
    """Szuka równodległości równobocznego układu 3-ciał."""
    result = force_zeros_2body(C, BETA, GAMMA)
    if result is None or len(result) < 2:
        print(f"  Brak równowagi dla β/C = {BETA/C:.3f} (próg 4.5)")
        return None, None
    d_inner, d_outer = result
    print(f"  β/C = {BETA/C:.3f}  (próg 4.5)")
    print(f"  d_wewnętrzne = {d_inner:.6f}  (niestabilna)")
    print(f"  d_zewnętrzne = {d_outer:.6f}  (stabilna — studnia potencjału)")
    V_in  = V_eff_total(d_inner, C, C, BETA, GAMMA)
    V_out = V_eff_total(d_outer, C, C, BETA, GAMMA)
    print(f"  V(d_wewn)    = {V_in:.6f}")
    print(f"  V(d_zewn)    = {V_out:.6f}")
    return d_inner, d_outer


# ── Krzywa potencjału ─────────────────────────────────────────────────

def plot_potential(d_in, d_out):
    d_vals = np.linspace(0.05, 6.0 * BETA, 2000)
    V_vals = V_eff_total(d_vals, C, C, BETA, GAMMA)

    fig, ax = plt.subplots(figsize=(9, 5))
    mask = (V_vals > -10) & (V_vals < 5)  # przytnij do sensownego zakresu
    ax.plot(d_vals[mask], V_vals[mask], 'b-', linewidth=2, label=r'$V_2(d)$')
    ax.axhline(0, color='k', linewidth=0.5)
    if d_in is not None:
        ax.axvline(d_in,  color='r', linestyle='--', label=f'd_wewn = {d_in:.3f} (niestabilna)')
        ax.axvline(d_out, color='g', linestyle='--', label=f'd_zewn = {d_out:.3f} (stabilna)')
        ax.plot(d_in,  V_eff_total(d_in,  C, C, BETA, GAMMA), 'ro', markersize=8)
        ax.plot(d_out, V_eff_total(d_out, C, C, BETA, GAMMA), 'gs', markersize=8)
    ax.set_xlabel('Odległość $d$', fontsize=12)
    ax.set_ylabel('$V_2(d)$ [energia potencjalna]', fontsize=12)
    ax.set_title(f'TGP potencjał 2-ciałowy  (β={BETA}, γ={GAMMA}, C={C})', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-3, 2)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex2_potential.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"  Zapisano: {os.path.abspath(out)}")


# ── Symulacja oscylacji ────────────────────────────────────────────────

def simulate_breathing(d_eq, amplitude=0.05):
    """Symuluje oscylacje oddechowe wokół równowagi d_eq."""
    # Geometria: równoboczny trójkąt o boku d_eq + perturbacja (pozycje 3D)
    h = d_eq * np.sqrt(3) / 2.0
    pos0 = np.array([
        [-d_eq/2,    -h/3,   0.0],
        [ d_eq/2,    -h/3,   0.0],
        [ 0.0,        2*h/3, 0.0]
    ])
    # Deformacja oddechowa: powiększenie o amplitude
    pos0 = pos0 * (1.0 + amplitude)

    vel0 = np.zeros_like(pos0)
    masses   = np.array([C, C, C])

    # Opakowania force_func i potential_func zgodne z API leapfrog_integrate
    def force_func(pos, C_vals):
        # leapfrog oczekuje przyspieszenia (nie siły)
        F = analytical_forces_tgp_pairwise(pos, C_vals, BETA, GAMMA)
        return F / C_vals[:, None]

    def potential_func(pos, C_vals):
        return potential_tgp(pos, C_vals, BETA, GAMMA)

    # Parametry całkowania
    dt = 5e-4
    t_max = 40.0
    save_every = 200

    print(f"  Symulacja leapfrog: t_max={t_max}, dt={dt}...")
    result = leapfrog_integrate(
        pos0, vel0, masses,
        force_func, potential_func,
        t_span=(0.0, t_max),
        dt=dt,
        save_every=save_every
    )

    pos   = result['positions']   # (n_save, 3, 3)
    t_arr = result['t']
    E_arr = result['energy']

    # Odległość między ciałem 0 a 1
    d_arr = np.linalg.norm(pos[:, 0, :] - pos[:, 1, :], axis=1)

    return t_arr, d_arr, E_arr


def plot_simulation(t_arr, d_arr, E_arr, d_eq):
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    axes[0].plot(t_arr, d_arr, 'b-', linewidth=1.0)
    axes[0].axhline(d_eq, color='r', linestyle='--', alpha=0.7, label=f'd_eq = {d_eq:.4f}')
    axes[0].set_ylabel('Odległość $d_{01}$', fontsize=11)
    axes[0].set_title('Oscylacje oddechowe równobocznego układu 3-ciał TGP', fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    dE = (E_arr - E_arr[0]) / abs(E_arr[0] + 1e-20)
    axes[1].semilogy(t_arr, np.abs(dE) + 1e-15, 'g-', linewidth=1.0)
    axes[1].set_ylabel(r'$|\Delta E / E_0|$', fontsize=11)
    axes[1].set_xlabel('Czas $t$', fontsize=11)
    axes[1].set_title('Zachowanie energii (leapfrog symplektyczny)', fontsize=12)
    axes[1].grid(True, alpha=0.3)

    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex2_oscillations.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"  Zapisano: {os.path.abspath(out)}")


# ── Skan fazowy: d_eq vs β/C ──────────────────────────────────────────

def phase_scan():
    """Skanuje równowagę dla różnych β/C."""
    ratios = np.linspace(4.51, 12.0, 200)
    d_out_list, d_in_list = [], []
    for r in ratios:
        beta_loc = r * C
        res = force_zeros_2body(C, beta_loc, beta_loc)
        if res and len(res) == 2:
            d_in_list.append(res[0])
            d_out_list.append(res[1])
        else:
            d_in_list.append(np.nan)
            d_out_list.append(np.nan)

    r_arr     = ratios
    d_in_arr  = np.array(d_in_list)
    d_out_arr = np.array(d_out_list)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r_arr, d_out_arr, 'g-', linewidth=2, label='$d_+$ (stabilna)')
    ax.plot(r_arr, d_in_arr,  'r--', linewidth=2, label='$d_-$ (niestabilna)')
    ax.axvline(4.5, color='k', linestyle=':', label='Próg β/C = 9/2')
    ax.set_xlabel(r'$\beta / C$', fontsize=12)
    ax.set_ylabel('Równoległość równowagi $d$', fontsize=12)
    ax.set_title('Diagram fazowy układu 3-ciał TGP (równoboczny)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex2_phase.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"  Zapisano: {os.path.abspath(out)}")


# ── Main ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*60)
    print("EX2: Równoboczny układ 3-ciał TGP")
    print("="*60)

    print(f"\nParametry: β={BETA}, γ={GAMMA}, C={C}")

    print("\n[1] Równowaga:")
    d_in, d_out = find_equilibria()

    if d_out is not None:
        print("\n[2] Krzywa potencjału:")
        plot_potential(d_in, d_out)

        print("\n[3] Skan fazowy β/C:")
        phase_scan()

        print("\n[4] Symulacja oscylacji (amplitude=5%):")
        t_arr, d_arr, E_arr = simulate_breathing(d_out, amplitude=0.05)
        plot_simulation(t_arr, d_arr, E_arr, d_out)

        # Estymacja okresu z FFT
        from numpy.fft import fft, fftfreq
        dt_sim = t_arr[1] - t_arr[0] if len(t_arr) > 1 else 1.0
        N = len(d_arr)
        freqs = fftfreq(N, d=dt_sim)
        power = np.abs(fft(d_arr - d_arr.mean()))**2
        peak_freq = np.abs(freqs[1:N//2][np.argmax(power[1:N//2])])
        T_period = 1.0 / peak_freq if peak_freq > 0 else float('inf')
        print(f"  Częstość własna (FFT): ω ≈ {2*np.pi*peak_freq:.4f}")
        print(f"  Okres oscylacji:       T ≈ {T_period:.4f}")

    print("\nGotowe.")
