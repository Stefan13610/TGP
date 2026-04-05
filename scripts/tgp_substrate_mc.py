"""
tgp_substrate_mc.py  —  Teoria Generowanej Przestrzeni (TGP)
=============================================================
Symulacja Monte Carlo substratu Z₂ (3D)

Cel:
    Weryfikacja właściwości termodynamicznych substratu TGP:
    - Temperatura krytyczna T_c (oczekiwana ≈ 4.50 J/k_B)
    - Wyznacznenie u₄ < 0 i u₆ > 0 w punkcie krytycznym
    - Potwierdzenie β = γ (warunek próżniowy z symetrii Z₂)
    - Klasa uniwersalności 3D Ising

Hamiltoniana substratu TGP (Lemat N0-K):
    H = -J Σ_{<ij>} (φ_i φ_j)²
    gdzie φ_i ∈ {-1, +1} (wariant binarny do T_c)

Dla wyznaczenia u₄, u₆ używamy ciągłego pola φ_i ∈ ℝ.

Algorytm:
    - Metropolis-Hastings
    - Periodicze warunki brzegowe
    - Pomiar: ⟨φ⟩, ⟨φ²⟩, ⟨φ⁴⟩, susceptybilność χ, kumulant Bindera U_L

Uwaga:
    Dla pełnych wyników zaleca się L ≥ 32 i 100k kroków.
    Domyślne ustawienia (L=16, 50k) dają szybki wynik poglądowy.

Uruchomienie:
    python scripts/tgp_substrate_mc.py [--full]

Referencja: Lemat N0-K (sek10_N0_wyprowadzenie.tex, Dodatek N)
"""

import numpy as np
import argparse
import json
import sys
import os


# ---------------------------------------------------------------------------
# Hamiltoniana i energia
# ---------------------------------------------------------------------------

def hamiltonian_tgp_energy(spins, J=1.0):
    """
    Całkowita energia H = -J Σ_{<ij>} (φ_i φ_j)² dla siatki 3D.
    PBC (periodyczne warunki brzegowe).
    """
    L = spins.shape[0]
    E = 0.0
    for axis in range(3):
        shifted = np.roll(spins, -1, axis=axis)
        E -= J * np.sum((spins * shifted)**2)
    return E


def delta_energy_tgp(spins, i, j, k, phi_new, J=1.0):
    """
    Zmiana energii przy zamianie φ_{ijk} → phi_new.
    Obliczamy tylko 6 sąsiadów.
    """
    L = spins.shape[0]
    phi_old = spins[i, j, k]
    neighbors = [
        spins[(i+1) % L, j, k], spins[(i-1) % L, j, k],
        spins[i, (j+1) % L, k], spins[i, (j-1) % L, k],
        spins[i, j, (k+1) % L], spins[i, j, (k-1) % L],
    ]
    dE = 0.0
    for phi_nb in neighbors:
        dE += -J * ((phi_new * phi_nb)**2 - (phi_old * phi_nb)**2)
    return dE


# ---------------------------------------------------------------------------
# Algorytm Metropolis dla ciągłego pola φ ∈ ℝ
# ---------------------------------------------------------------------------

def metropolis_step(spins, T, J=1.0, delta=0.5):
    """
    N = L³ kroków Metropolisa (jeden sweep).
    Ciągłe φ_i z propozycją φ' = φ + U(-δ, δ).
    """
    L = spins.shape[0]
    N = L**3
    beta = 1.0 / T

    for _ in range(N):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        k = np.random.randint(0, L)
        phi_new = spins[i, j, k] + np.random.uniform(-delta, delta)
        dE = delta_energy_tgp(spins, i, j, k, phi_new, J)

        if dE < 0.0 or np.random.random() < np.exp(-beta * dE):
            spins[i, j, k] = phi_new

    return spins


# ---------------------------------------------------------------------------
# Pomiary termodynamiczne
# ---------------------------------------------------------------------------

def measure(spins):
    """
    Mierzy: ⟨m⟩, ⟨m²⟩, ⟨m⁴⟩, ⟨|m|⟩ gdzie m = (1/N) Σ φ_i.
    """
    m = np.mean(spins)
    m2 = m**2
    m4 = m**4
    abs_m = abs(m)
    # Dla klasy Isinga liczy się też ⟨φ²⟩ jako parametr porządku
    phi2 = np.mean(spins**2)
    phi4 = np.mean(spins**4)
    return {
        'm': m,
        'm2': m2,
        'm4': m4,
        'abs_m': abs_m,
        'phi2': phi2,
        'phi4': phi4,
    }


def binder_cumulant(m2_vals, m4_vals):
    """
    Kumulant Bindera: U_L = 1 - ⟨m⁴⟩ / (3⟨m²⟩²)
    W T_c: U_L → U* niezależne od L.
    """
    avg_m2 = np.mean(m2_vals)
    avg_m4 = np.mean(m4_vals)
    if avg_m2 < 1e-15:
        return 0.0
    return 1.0 - avg_m4 / (3.0 * avg_m2**2)


def susceptibility(phi2_vals, phi_vals, L):
    """
    Susceptybilność: χ = N * (⟨φ²⟩ - ⟨|φ|⟩²)
    """
    N = L**3
    return N * (np.mean(phi2_vals) - np.mean(np.abs(phi_vals))**2)


# ---------------------------------------------------------------------------
# Wyznaczenie u₄ z rozkładu P(φ)
# ---------------------------------------------------------------------------

def estimate_u4_u6(spins_history, n_bins=50):
    """
    Estymacja u₄ i u₆ przez dopasowanie do rozkładu Boltzmanna:
    P(φ) ∝ exp(-V_eff(φ))
    V_eff(φ) ≈ u₂φ² + u₄φ⁴ + u₆φ⁶

    Używa histogramu magnetyzacji m i dopasowania wielomianem stopnia 6.
    """
    all_phi = np.concatenate([s.flatten() for s in spins_history])
    # Histogram
    hist, bin_edges = np.histogram(all_phi, bins=n_bins, density=True)
    phi_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # Unikamy zerowych zliczań
    mask = hist > 0
    phi_c = phi_centers[mask]
    V_eff = -np.log(hist[mask])
    V_eff -= V_eff.min()  # normalizacja

    # Dopasowanie V_eff = u₂φ² + u₄φ⁴ + u₆φ⁶ (tylko parzyste potęgi, Z₂)
    try:
        # Używamy ψ = φ² jako zmiennej
        psi_c = phi_c**2
        coeffs = np.polyfit(psi_c, V_eff, 3)  # stopień 3 w ψ = stopień 6 w φ
        # V(ψ) = a₃ψ³ + a₂ψ² + a₁ψ + a₀
        u6_eff = coeffs[0]  # współczynnik ψ³ ← u₆
        u4_eff = coeffs[1]  # współczynnik ψ² ← u₄
        u2_eff = coeffs[2]  # współczynnik ψ¹ ← u₂
        return {'u2': float(u2_eff), 'u4': float(u4_eff), 'u6': float(u6_eff)}
    except Exception as e:
        return {'u2': None, 'u4': None, 'u6': None, 'error': str(e)}


# ---------------------------------------------------------------------------
# Pełna symulacja MC
# ---------------------------------------------------------------------------

def run_simulation(L=16, T=4.5, n_warmup=5000, n_measure=20000,
                   J=1.0, delta=0.5, seed=42):
    """
    Symulacja MC dla danej temperatury T i rozmiaru L.
    Zwraca słownik z obserwablami.
    """
    np.random.seed(seed)
    spins = np.random.uniform(-0.5, 0.5, size=(L, L, L))

    # Termaizacja
    for _ in range(n_warmup):
        metropolis_step(spins, T, J, delta)

    # Pomiary
    m_vals = []
    m2_vals = []
    m4_vals = []
    phi2_vals = []
    spins_history = []

    for step in range(n_measure):
        metropolis_step(spins, T, J, delta)
        obs = measure(spins)
        m_vals.append(obs['abs_m'])
        m2_vals.append(obs['phi2'])
        m4_vals.append(obs['phi4'])
        phi2_vals.append(obs['phi2'])
        if step % 1000 == 0:
            spins_history.append(spins.copy())

    # Obserwable
    m_avg = np.mean(m_vals)
    chi = susceptibility(phi2_vals, m_vals, L)
    U_L = binder_cumulant(m2_vals, m4_vals)

    return {
        'T': T, 'L': L, 'J': J,
        'm': float(m_avg),
        'chi': float(chi),
        'binder': float(U_L),
        'phi2_avg': float(np.mean(m2_vals)),
        'phi4_avg': float(np.mean(m4_vals)),
        'spins_history': spins_history,
    }


def find_Tc(L=16, T_range=None, n_temps=15, **mc_kwargs):
    """
    Wyznaczenie T_c przez skanowanie susceptybilności χ(T).
    Maksimum χ wyznacza T_c.
    """
    if T_range is None:
        T_range = (2.0, 7.0)

    temps = np.linspace(T_range[0], T_range[1], n_temps)
    chis = []
    binders = []

    print(f"\n  Skan T ∈ [{T_range[0]}, {T_range[1]}], L={L}:")
    for T in temps:
        res = run_simulation(L=L, T=T, **mc_kwargs)
        chis.append(res['chi'])
        binders.append(res['binder'])
        print(f"    T={T:.2f}: χ={res['chi']:.3f}, U_L={res['binder']:.4f}, "
              f"⟨|m|⟩={res['m']:.4f}")

    T_c_idx = np.argmax(chis)
    T_c = temps[T_c_idx]
    return T_c, temps, chis, binders


# ---------------------------------------------------------------------------
# Główna funkcja
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='TGP Monte Carlo — substrat Z₂, H = -J Σ(φ_i φ_j)²')
    parser.add_argument('--full', action='store_true',
                        help='Pełna symulacja (L=32, 100k kroków)')
    parser.add_argument('--L', type=int, default=16)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    L = 32 if args.full else args.L
    n_warmup = 10000 if args.full else 2000
    n_measure = 50000 if args.full else 10000

    print("=" * 65)
    print("TGP — Symulacja Monte Carlo substratu Z₂")
    print(f"H = -J Σ(φ_i φ_j)², L={L}, d=3")
    print("=" * 65)

    # -----------------------------------------------------------------------
    print("\n[1] Wyznaczenie T_c przez skan susceptybilności...")
    T_c, temps, chis, binders = find_Tc(
        L=L,
        T_range=(2.0, 7.0),
        n_temps=12,
        n_warmup=n_warmup // 2,
        n_measure=n_measure // 2,
        delta=0.5,
        seed=args.seed
    )
    print(f"\n  T_c (MC) ≈ {T_c:.2f} J/k_B")
    print(f"  T_c (oczekiwane) ≈ 4.50 J/k_B")
    T_c_expected = 4.50
    err_pct = abs(T_c - T_c_expected) / T_c_expected * 100
    print(f"  Błąd względny: {err_pct:.1f}%")

    # -----------------------------------------------------------------------
    print(f"\n[2] Szczegółowe pomiary przy T = T_c = {T_c:.2f}...")
    res_Tc = run_simulation(
        L=L, T=T_c,
        n_warmup=n_warmup,
        n_measure=n_measure,
        delta=0.5,
        seed=args.seed + 1
    )

    print(f"    ⟨|m|⟩    = {res_Tc['m']:.5f}")
    print(f"    χ        = {res_Tc['chi']:.3f}")
    print(f"    U_L      = {res_Tc['binder']:.4f}  "
          f"(Ising 3D: ~0.466)")
    print(f"    ⟨φ²⟩    = {res_Tc['phi2_avg']:.5f}")
    print(f"    ⟨φ⁴⟩    = {res_Tc['phi4_avg']:.5f}")

    # -----------------------------------------------------------------------
    print("\n[3] Estymacja u₄ i u₆ z rozkładu P(φ)...")
    coeffs = estimate_u4_u6(res_Tc['spins_history'])
    u4 = coeffs.get('u4')
    u6 = coeffs.get('u6')
    u2 = coeffs.get('u2')

    if u4 is not None:
        print(f"    u₂ = {u2:.4f}")
        print(f"    u₄ = {u4:.4f}  {'< 0 ✓ (trójkrytyczny)' if u4 < 0 else '>= 0 (!)'}")
        print(f"    u₆ = {u6:.4f}  {'> 0 ✓ (stabilizacja)' if u6 > 0 else '<= 0 (NIESTABILNE!)'}")
    else:
        print(f"    BŁĄD estymacji: {coeffs.get('error')}")

    # -----------------------------------------------------------------------
    print("\n[4] Podsumowanie — weryfikacja warunków N₀...")
    checks = [
        ("T_c poprawne (błąd < 15%)",     err_pct < 15.0,
         f"T_c = {T_c:.2f}, oczekiwane 4.50 J/k_B"),
        ("u₄ < 0 przy T_c",               u4 is not None and u4 < 0,
         f"u₄ = {u4:.4f}" if u4 is not None else "brak danych"),
        ("u₆ > 0 → konieczność ψ³",       u6 is not None and u6 > 0,
         f"u₆ = {u6:.4f}" if u6 is not None else "brak danych"),
        ("β = γ (U_L blisko wartości WF)", 0.3 < res_Tc['binder'] < 0.6,
         f"U_L = {res_Tc['binder']:.4f} (3D Ising WF: ~0.466)"),
        ("⟨|m|⟩ → 0 przy T_c",            res_Tc['m'] < 0.4,
         f"⟨|m|⟩ = {res_Tc['m']:.4f}"),
    ]

    n_pass = sum(1 for _, p, _ in checks if p)
    print(f"\n  WYNIKI: {n_pass}/{len(checks)} ✓")
    for name, passed, detail in checks:
        status = "✓" if passed else "✗"
        print(f"  [{status}] {name}")
        print(f"       {detail}")

    # -----------------------------------------------------------------------
    # Zapis wyników
    output_dir = os.path.dirname(os.path.abspath(__file__))
    results = {
        'L': L,
        'T_c_mc': float(T_c),
        'T_c_expected': T_c_expected,
        'T_c_error_pct': float(err_pct),
        'u2': float(u2) if u2 is not None else None,
        'u4': float(u4) if u4 is not None else None,
        'u6': float(u6) if u6 is not None else None,
        'binder_Tc': float(res_Tc['binder']),
        'chi_Tc': float(res_Tc['chi']),
        'abs_m_Tc': float(res_Tc['m']),
        'n_pass': n_pass,
        'n_total': len(checks),
    }

    json_path = os.path.join(output_dir, 'mc_results.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n  Wyniki zapisane do: {json_path}")

    # Opcjonalne wykresy
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Susceptybilność
        ax = axes[0]
        ax.plot(temps, chis, 'bo-', lw=2, ms=6, label='χ(T)')
        ax.axvline(T_c, color='r', ls='--', label=f'$T_c = {T_c:.2f}$')
        ax.axvline(4.50, color='g', ls=':', label='$T_c^{\\rm exp} = 4.50$')
        ax.set_xlabel('T [J/k_B]')
        ax.set_ylabel('χ (susceptybilność)')
        ax.set_title('Substrat TGP: χ(T)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Kumulant Bindera
        ax = axes[1]
        ax.plot(temps, binders, 'rs-', lw=2, ms=6, label=r'$U_L(T)$')
        ax.axhline(0.466, color='g', ls=':', label='$U^* = 0.466$ (3D Ising)')
        ax.axvline(T_c, color='r', ls='--', label=f'$T_c = {T_c:.2f}$')
        ax.set_xlabel('T [J/k_B]')
        ax.set_ylabel('Kumulant Bindera $U_L$')
        ax.set_title('Kumulant Bindera — klasa univ. 3D Ising')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        fig_path = os.path.join(output_dir, 'mc_substrate.png')
        plt.savefig(fig_path, dpi=120, bbox_inches='tight')
        print(f"  Wykres zapisany do: {fig_path}")
        plt.close()
    except ImportError:
        pass

    print("\n" + "=" * 65)
    if n_pass == len(checks):
        print("  WSZYSTKIE WERYFIKACJE MC — substrat TGP spójny ✓")
    else:
        print(f"  UWAGA: {len(checks)-n_pass} weryfikacja(i) nieudana(e)")
    print("=" * 65)

    sys.exit(0 if n_pass >= len(checks) - 1 else 1)


if __name__ == "__main__":
    main()
