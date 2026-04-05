"""
TGP — Sekwencyjny model generacji
===================================

Kluczowy insight uzytkownika:
  Phi0 to srednia kosmiczna WLACZNIE z czastka ktora liczymy.

Konsekwencja: kazda generacja widzi inne tlo Phi0,
wzmocnione przez wszystkie LZEJSZE generacje:

  Gen 1: Phi0^(1) = Phi0_base                     (brak lzejszych)
  Gen 2: Phi0^(2) = Phi0_base + xi * M1            (gen-1 wzmacnia tlo)
  Gen 3: Phi0^(3) = Phi0_base + xi * (M1 + M2)    (gen-1 i gen-2 wzmacniaja)

Samospojnosc dla kazdej generacji niezaleznie:
  Mn = E(Kn; Phi0^(n)) / c0^2

gdzie Kn = Mn / (4*pi*Phi0^(n))

Dla duzego Phi0: E ~ Phi0^2 * E_tilde(kappa)
=> Mn ~ Phi0^(n)^2 * mu  (skala kwadratem tla!)

To naturalnie tworzy SILNA hierarchie mas.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAMMA = 1.0; PHI0_BASE = 1.0; Q = 1.0

# ======================================================================
# ENERGIA TGP w tle Phi0
# ======================================================================
def energy(K, Phi0, alpha, a_gam):
    """Energia Yukawa wzgledem tla Phi0, z efektywna msp(Phi0)."""
    msp = np.sqrt(max(GAMMA * PHI0_BASE / Phi0, 1e-6))
    r_max = max(30.0 / msp, 10.0)
    r = np.linspace(a_gam, r_max, 2000)

    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / (Phi0 * phi)) * r**2, r)

    psi  = phi / Phi0
    Ep   = 4*np.pi * Phi0**2 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4 - GAMMA/3 + GAMMA/4) * r**2, r)

    return Ek + Ep

def g_fixed_Phi0(M, Phi0, alpha, a_gam):
    """g(M) = E(K(M); Phi0_fixed) - M dla stalego tla."""
    if M < 1e-10: return 0.0
    K = M / (4.0 * np.pi * Phi0)
    return energy(K, Phi0, alpha, a_gam) - M

def find_masses_in_background(Phi0, alpha, a_gam, M_max=500.0, N=400):
    """Znajdz samospojne masy w tle Phi0."""
    M_lo  = np.linspace(1e-4, 2.0,    int(N*0.3))
    M_mid = np.linspace(2.0,  50.0,   int(N*0.4))
    M_hi  = np.linspace(50.0, M_max,  int(N*0.3))
    M_arr = np.unique(np.concatenate([M_lo, M_mid, M_hi]))

    g_arr = np.array([g_fixed_Phi0(M, Phi0, alpha, a_gam) for M in M_arr])

    roots = []
    for i in range(len(g_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                r = brentq(
                    lambda M: g_fixed_Phi0(M, Phi0, alpha, a_gam),
                    M_arr[i], M_arr[i+1], xtol=1e-8)
                roots.append(r)
            except Exception:
                pass
    return roots

# ======================================================================
# SEKWENCYJNY MODEL GENERACJI
# ======================================================================
def sequential_generations(xi, alpha, a_gam):
    """
    Oblicza masy 3 generacji sekwencyjnie.

    Gen 1: tlo = Phi0_base
    Gen 2: tlo = Phi0_base + xi * M1
    Gen 3: tlo = Phi0_base + xi * (M1 + M2)
    """
    results = []

    # --- GENERACJA 1 ---
    Phi0_1 = PHI0_BASE
    roots_1 = find_masses_in_background(Phi0_1, alpha, a_gam)
    if not roots_1:
        return None
    M1 = roots_1[0]   # najmniejszy korzen = najlzejsza czastka gen1
    results.append({'gen': 1, 'Phi0': Phi0_1, 'M': M1, 'roots': roots_1})

    # --- GENERACJA 2 ---
    Phi0_2 = PHI0_BASE + xi * M1
    roots_2 = find_masses_in_background(Phi0_2, alpha, a_gam)
    if not roots_2:
        return results
    M2 = roots_2[0]
    results.append({'gen': 2, 'Phi0': Phi0_2, 'M': M2, 'roots': roots_2})

    # --- GENERACJA 3 ---
    Phi0_3 = PHI0_BASE + xi * (M1 + M2)
    roots_3 = find_masses_in_background(Phi0_3, alpha, a_gam)
    if not roots_3:
        return results
    M3 = roots_3[0]
    results.append({'gen': 3, 'Phi0': Phi0_3, 'M': M3, 'roots': roots_3})

    return results

# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':

    print("="*70)
    print("TGP — Sekwencyjny model generacji")
    print("Gen n widzi tlo Phi0_n wzmocnione przez wszystkie lzejsze generacje")
    print("="*70)

    # --- Skan parametrow ---
    print(f"\n{'xi':>7} {'alpha':>6} {'a_gam':>6} "
          f"{'M1':>9} {'M2':>9} {'M3':>10} "
          f"{'M2/M1':>8} {'M3/M1':>9} {'M3/M2':>8}")
    print("-"*80)

    target_21 = 207.0
    target_31 = 3477.0

    best_score = 1e9
    best_params = None

    for xi in [0.005, 0.01, 0.02, 0.05, 0.10, 0.20]:
        for alpha in [4.0, 5.0, 5.9, 6.0, 7.0, 8.0, 10.0]:
            for a_gam in [0.03, 0.05, 0.08, 0.10]:
                res = sequential_generations(xi, alpha, a_gam)
                if not res or len(res) < 2:
                    continue

                M1 = res[0]['M']
                M2 = res[1]['M'] if len(res) >= 2 else None
                M3 = res[2]['M'] if len(res) >= 3 else None

                r21 = M2/M1 if M2 else None
                r31 = M3/M1 if M3 else None
                r32 = M3/M2 if (M3 and M2) else None

                s_M1  = f"{M1:.4f}"
                s_M2  = f"{M2:.3f}"  if M2  else "-"
                s_M3  = f"{M3:.1f}"  if M3  else "-"
                s_r21 = f"{r21:.1f}" if r21 else "-"
                s_r31 = f"{r31:.0f}" if r31 else "-"
                s_r32 = f"{r32:.1f}" if r32 else "-"

                flag = ""
                if M3:
                    score = (abs(r21-target_21)/target_21 +
                             abs(r31-target_31)/target_31)
                    if score < best_score:
                        best_score = score
                        best_params = (xi, alpha, a_gam, res)
                    flag = "  ***"

                if M2:
                    print(f"{xi:>7.3f} {alpha:>6.1f} {a_gam:>6.3f} "
                          f"{s_M1:>9} {s_M2:>9} {s_M3:>10} "
                          f"{s_r21:>8} {s_r31:>9} {s_r32:>8}{flag}")

    # --- Najlepsze dopasowanie ---
    print(f"\n{'='*70}")
    print(f"Cel: M2/M1 = {target_21}, M3/M1 = {target_31}")
    if best_params:
        xi_b, alp_b, ag_b, res_b = best_params
        M1_b = res_b[0]['M']
        M2_b = res_b[1]['M']
        M3_b = res_b[2]['M']
        print(f"Najlepsze: xi={xi_b}, alpha={alp_b}, a_gam={ag_b}")
        print(f"  M1 = {M1_b:.5f}")
        print(f"  M2 = {M2_b:.4f}  M2/M1 = {M2_b/M1_b:.1f}  (cel: 207)")
        print(f"  M3 = {M3_b:.2f}   M3/M1 = {M3_b/M1_b:.0f}  (cel: 3477)")
        print(f"  Score = {best_score:.4f}")
    else:
        print("Brak 3 generacji — potrzeba szerszego zakresu parametrow.")

    # --- Analityczne skalowanie ---
    print(f"\n{'='*70}")
    print("Analityczna predykcja: E ~ Phi0^2 => Mn ~ Phi0_n^2 * mu")
    print(f"  (Phi0_2/Phi0_1)^2 powinno dawac M2/M1")
    print()
    for xi in [0.005, 0.01, 0.02]:
        for alpha in [5.0, 5.9, 6.0]:
            res = sequential_generations(xi, alpha, 0.05)
            if not res or len(res) < 2:
                continue
            M1 = res[0]['M']; Phi0_1 = res[0]['Phi0']
            M2 = res[1]['M']; Phi0_2 = res[1]['Phi0']
            ratio_Phi0 = Phi0_2/Phi0_1
            ratio_M    = M2/M1
            pred_from_Phi0 = ratio_Phi0**2
            print(f"  xi={xi:.3f} alpha={alpha:.1f}: "
                  f"Phi0_2/Phi0_1={ratio_Phi0:.4f}, "
                  f"(Phi0_2/Phi0_1)^2={pred_from_Phi0:.2f}, "
                  f"M2/M1={ratio_M:.1f}")

    # --- Wykres krzywej g(M) dla kazdej generacji ---
    if best_params:
        xi_b, alp_b, ag_b, res_b = best_params
    else:
        xi_b, alp_b, ag_b = 0.01, 5.9, 0.05
        res_b = sequential_generations(xi_b, alp_b, ag_b)

    if res_b and len(res_b) >= 2:
        fig, axes = plt.subplots(1, min(len(res_b), 3), figsize=(15, 5))
        if len(res_b) == 1:
            axes = [axes]

        colors = ['blue', 'orange', 'red']
        gen_labels = ['Gen 1 (elektron)', 'Gen 2 (mion)', 'Gen 3 (tau)']

        for idx, r in enumerate(res_b[:3]):
            ax = axes[idx]
            Phi0_n = r['Phi0']
            M_lo  = np.linspace(1e-4, 5.0,   150)
            M_mid = np.linspace(5.0,  r['M']*3, 150)
            M_arr = np.unique(np.concatenate([M_lo, M_mid]))
            g_arr = np.array([g_fixed_Phi0(M, Phi0_n, alp_b, ag_b) for M in M_arr])

            mask = np.isfinite(g_arr) & (g_arr > -500) & (g_arr < 500)
            ax.plot(M_arr[mask], g_arr[mask], '-', color=colors[idx], lw=2,
                    label=f'g(M), Phi0={Phi0_n:.3f}')
            ax.axhline(0, color='k', lw=1, ls='--')

            for i, Mc in enumerate(r['roots'][:2]):
                ax.axvline(Mc, color='green', lw=1.2, ls=':')
                ax.plot(Mc, 0, 'go', ms=10, label=f'M={Mc:.3f}')

            ax.set_title(f'{gen_labels[idx]}\nPhi0={Phi0_n:.3f}, msp={np.sqrt(GAMMA/Phi0_n):.4f}',
                         fontsize=10)
            ax.set_xlabel('M', fontsize=10)
            ax.set_ylabel('g(M) = E - M', fontsize=10)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            y_rng = np.nanmax(np.abs(g_arr[mask]))*0.4 if mask.sum()>0 else 5
            ax.set_ylim(-y_rng, y_rng)

        plt.suptitle(f'TGP sekwencyjny: xi={xi_b}, alpha={alp_b}, a_Gam={ag_b}',
                     fontsize=12)
        plt.tight_layout()
        plt.savefig('TGP/TGP_v1/scripts/advanced/sekwencyjny_model.png',
                    dpi=150, bbox_inches='tight')
        print("\nWykres: sekwencyjny_model.png")

    print("\nGotowe.")
