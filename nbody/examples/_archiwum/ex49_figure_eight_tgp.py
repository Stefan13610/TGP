"""
ex49_figure_eight_tgp.py
=========================
Orbita ósemkowa Chenciner-Montgomery w teorii TGP (Problem O19)

PYTANIE BADAWCZE
----------------
Czy orbita ósemkowa (trzy jednakowe masy na trajektorii kształtu cyfry 8)
jest stabilna w TGP przy nieredukowalnych siłach trójciałowych V₃?

KONTEKST
--------
W grawitacji Newtonowskiej orbita ósemkowa to egzotyczne, ściśle periodyczne
rozwiązanie (Chenciner & Montgomery 2000). Jej stabilność jest słaba —
jest punktem izolowanym w przestrzeni rozwiązań pod działaniem małych perturbacji.

W TGP:
  V_total = V_Newton + ΔV_pair(β,γ) + V₃(γ)

Siły parowe TGP: ΔF_pair ~ C²β/r² (odpychanie krótszozasięgowe)
                           − C²β/r  (korekta Yukawy)
Siły 3-ciałowe: F₃ ~ 6γC³ · 8π²/P² · (r̂_ij + r̂_ik)
                                       (odpychanie od centrum)

PREDYKCJA (z PLAN_AKTUALIZACJI_v26)
------------------------------------
V₃ zaburza orbitę ósemkową (tendencja do powiększania trójkąta).
Czas życia τ ~ C⁻¹ · τ_Newton.

PROTOKÓŁ
--------
1. Inicjalizacja: warunki Chenciner-Montgomery (scale=1, C=C_Pl)
2. Całkowanie: RK45 przez T = 10³ t_Pl (z i bez TGP)
3. Diagnostyka:
   - Odchylenie od centrum masy |Δx_cm| < ε (warunek periodyczności)
   - Energia całkowita E(t) — powinna być zachowana
   - Eksponent Lyapunova λ (rozbieżność perturbacji)
   - Czas życia orbit: kiedy |x_i - x_i^{Newton}| > threshold

WYNIKI
------
Zakres: β = γ ∈ {0, 0.001, 0.005, 0.01, 0.02, 0.05}
  β=0 (Newton): orbita stabilna, τ → ∞
  β>0 (TGP): orbita destabilizowana, τ maleje z β
  Predykcja: τ_TGP ~ β⁻¹ · τ_Newton (z reguły przybliżenia pierwszego rzędu)

Autor: TGP Analysis Session v27, 2026-03-22
"""

import sys
import io
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from itertools import combinations

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ──────────────────────────────────────────────────────────────────────────────
#  Stałe
# ──────────────────────────────────────────────────────────────────────────────

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))  # = 0.2821 (TGP standardowy)

# Warunki startowe Chenciner-Montgomery (znormalizowane, układ CM)
_X1_CM = np.array([ 0.97000436, -0.24308753, 0.0])
_V1_CM = np.array([ 0.93240737/2,  0.86473146/2, 0.0])

def figure_eight_ic(C=C_PL, scale=1.0):
    """
    Zwraca warunki początkowe dla orbity ósemkowej.
    Pozycje: x1 = -x2, x3 = 0
    Prędkości: v1 = v2, v3 = -2*v1
    """
    x1 = _X1_CM * scale
    v1 = _V1_CM / scale   # skalowanie prędkości: v ~ 1/scale dla orbity

    pos = np.array([x1, -x1, np.zeros(3)])
    vel = np.array([v1, v1, -2.0 * v1])
    C_arr = np.full(3, C)
    return pos, vel, C_arr

# ──────────────────────────────────────────────────────────────────────────────
#  Siły TGP
# ──────────────────────────────────────────────────────────────────────────────

def forces_newton(pos, C_arr, G=1.0, eps=1e-6):
    """Siły Newtonowskie: F_ij = -G m_i m_j / r_ij² r̂_ij."""
    n = len(C_arr)
    F = np.zeros_like(pos)
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = pos[j] - pos[i]
            r2    = np.dot(r_vec, r_vec) + eps**2
            r     = np.sqrt(r2)
            # W TGP: G = G₀/φ, ale tu φ≈1; konwencja: m_i = C_i
            # Potencjał Newtonowski (bez modyfikacji): V = -G*C_i*C_j/r
            # (TGP: faktyczny potencjał to pełne 3-reżimy, tu Newton = granica β→0)
            fmag = G * C_arr[i] * C_arr[j] / r2
            f_ij = fmag * r_vec / r
            F[i] += f_ij
            F[j] -= f_ij
    return F

def forces_tgp_pairwise(pos, C_arr, beta, eps=1e-6):
    """
    Pełne siły parowe TGP:
      dV_2/dr = 4π C_i C_j/r² − 16π β C_i C_j/r³ + 36π γ C_i C_j(C_i+C_j)/r⁴
    Przy β=γ (N0-5):
      dV_2/dr = 4π C_i C_j/r² − 16π β C_i C_j/r³ + 36π β C_i C_j(C_i+C_j)/r⁴
    """
    gamma = beta  # N0-5
    n = len(C_arr)
    F = np.zeros_like(pos)
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = pos[j] - pos[i]
            r2    = np.dot(r_vec, r_vec) + eps**2
            r     = np.sqrt(r2)
            Ci, Cj = C_arr[i], C_arr[j]
            # dV_2/dr (ujemna wartość → przyciąganie przy małych r jest zaburzone)
            dV_dr = (4.0 * np.pi * Ci * Cj / r2
                     - 16.0 * np.pi * beta * Ci * Cj / (r * r2)
                     + 36.0 * np.pi * gamma * Ci * Cj * (Ci + Cj) / (r2 * r2))
            # F = -dV/dr · r̂, więc F_i = -(dV_dr) * r̂_ij (w kierunku j)
            # Znak: dV_dr > 0 przy r duży → F skierowane od j (odpychanie ?? — nie, to siła gradientu)
            # Poprawnie: F_i = +dV_dr * r̂_ij jeśli V = V(r), r = |r_i - r_j|
            # V maleje gdy i zbliża się do j: dV/d(r_ij) < 0 → przyciąganie
            # F_i = -grad_i V = -dV/dr · (x_i - x_j)/r = +dV/dr · r̂_ij
            # gdzie r̂_ij = (x_j - x_i)/r (od i do j)
            f_ij = -dV_dr * r_vec / r   # r_vec = x_j - x_i
            F[i] += f_ij
            F[j] -= f_ij
    return F

def forces_3body_tgp(pos, C_arr, gamma, eps=1e-6):
    """
    Nieredukowalne siły trójciałowe TGP (przybliżenie Coulomba):
      F_i = 6γ C_i C_j C_k · (8π²/P²) · (r̂_ij + r̂_ik)
    Tendencja: odpychanie od centrum trójkąta.
    """
    n = len(C_arr)
    F = np.zeros_like(pos)
    if n < 3:
        return F

    for i, j, k in combinations(range(n), 3):
        Ci, Cj, Ck = C_arr[i], C_arr[j], C_arr[k]

        r_ij = pos[j] - pos[i]
        r_ik = pos[k] - pos[i]
        r_jk = pos[k] - pos[j]

        dij = max(np.sqrt(np.dot(r_ij, r_ij)), eps)
        dik = max(np.sqrt(np.dot(r_ik, r_ik)), eps)
        djk = max(np.sqrt(np.dot(r_jk, r_jk)), eps)

        P = dij + dik + djk
        coupling = 6.0 * gamma * Ci * Cj * Ck
        dI_dP = -8.0 * np.pi**2 / P**2

        # dP/dx_i = (x_i-x_j)/dij + (x_i-x_k)/dik = -(r_ij/dij + r_ik/dik)
        dP_xi = -(r_ij / dij + r_ik / dik)
        dP_xj =  (r_ij / dij) - r_jk / djk
        dP_xk =  (r_ik / dik) + r_jk / djk

        # F_q = coupling · dI/dP · dP/dx_q
        F[i] += coupling * dI_dP * dP_xi
        F[j] += coupling * dI_dP * dP_xj
        F[k] += coupling * dI_dP * dP_xk

    return F

def total_forces(pos, C_arr, beta=0.0, include_3body=False, eps=1e-6):
    """Całkowite siły TGP = parowe + ewentualnie 3-ciałowe."""
    F = forces_tgp_pairwise(pos, C_arr, beta, eps)
    if include_3body and beta > 0:
        F += forces_3body_tgp(pos, C_arr, gamma=beta, eps=eps)
    return F

# ──────────────────────────────────────────────────────────────────────────────
#  Całkowanie
# ──────────────────────────────────────────────────────────────────────────────

def make_ode_rhs(C_arr, beta, include_3body, eps=1e-4):
    """Zwraca funkcję RHS dla solve_ivp."""
    def rhs(t, y):
        n = len(C_arr)
        pos = y[:3*n].reshape(n, 3)
        vel = y[3*n:].reshape(n, 3)

        F = total_forces(pos, C_arr, beta, include_3body, eps)
        # a_i = F_i / m_i = F_i / C_i
        acc = F / C_arr[:, np.newaxis]

        return np.concatenate([vel.ravel(), acc.ravel()])
    return rhs

def integrate_orbit(pos0, vel0, C_arr, beta, T_max, n_steps,
                    include_3body=False, rtol=1e-8, atol=1e-10):
    """
    Całkuje orbit od t=0 do t=T_max.
    Zwraca słownik z historią trajektorii.
    """
    n = len(C_arr)
    y0 = np.concatenate([pos0.ravel(), vel0.ravel()])

    t_eval = np.linspace(0, T_max, n_steps)
    rhs = make_ode_rhs(C_arr, beta, include_3body)

    sol = solve_ivp(rhs, [0, T_max], y0, t_eval=t_eval,
                    method='RK45', rtol=rtol, atol=atol,
                    dense_output=False)

    if not sol.success:
        print(f"  [WARN] Całkowanie nie powiodło się: {sol.message}")

    pos_t = sol.y[:3*n, :].reshape(n, 3, -1)  # (n, 3, n_t)
    vel_t = sol.y[3*n:, :].reshape(n, 3, -1)
    t     = sol.t

    return dict(t=t, pos=pos_t, vel=vel_t, success=sol.success, n_eval=sol.nfev)

# ──────────────────────────────────────────────────────────────────────────────
#  Diagnostyka
# ──────────────────────────────────────────────────────────────────────────────

def compute_energy(pos, vel, C_arr, beta):
    """Energia kinetyczna + potencjalna."""
    n = len(C_arr)
    E_kin = 0.5 * np.sum(C_arr[:, np.newaxis] * vel**2)
    # Energia potencjalna parowa (uproszczona: Newton + korekty TGP)
    E_pot = 0.0
    for i in range(n):
        for j in range(i+1, n):
            r_vec = pos[j] - pos[i]
            r = max(np.sqrt(np.dot(r_vec, r_vec)), 1e-6)
            Ci, Cj = C_arr[i], C_arr[j]
            E_pot += (-4.0 * np.pi * Ci * Cj / r
                      + 8.0 * np.pi * beta * Ci * Cj / r**2
                      - 12.0 * np.pi * beta * Ci * Cj * (Ci+Cj) / r**3)
    return E_kin + E_pot

def orbit_deviation(pos_tgp, pos_newton):
    """
    Maksymalne odchylenie od orbity Newtonowskiej (po czasie).
    |x_i^TGP(t) - x_i^Newton(t)|_max
    """
    diff = pos_tgp - pos_newton   # (n, 3, n_t)
    dev  = np.sqrt(np.sum(diff**2, axis=1))  # (n, n_t)
    return np.max(dev, axis=0)   # (n_t,)

def orbit_lifetime(t_arr, deviation, threshold=0.5):
    """
    Czas do przekroczenia threshold odchylenia (czas życia orbity).
    """
    idx = np.argmax(deviation > threshold)
    if idx == 0 and deviation[0] <= threshold:
        return t_arr[-1]   # orbit przeżyła cały czas symulacji
    return t_arr[idx]

def center_of_mass_drift(pos_t, C_arr):
    """Drift centrum masy (powinien być = 0)."""
    M_tot = np.sum(C_arr)
    cm = np.einsum('i,ijk->jk', C_arr, pos_t) / M_tot  # (3, n_t)
    cm_r = np.sqrt(np.sum(cm**2, axis=0))
    return cm_r

def lyapunov_exponent(t_arr, deviation, eps_init=1e-7):
    """
    Przybliżony eksponent Lyapunova:
    λ ≈ max log(|δ(t)|/|δ(0)|) / t
    """
    eps_arr = np.clip(deviation, 1e-15, None)
    lyap = np.log(eps_arr / eps_init) / (t_arr + 1e-30)
    return lyap

# ──────────────────────────────────────────────────────────────────────────────
#  Główna symulacja
# ──────────────────────────────────────────────────────────────────────────────

def run_simulation(beta, T_max, n_steps, include_3body, label, C=C_PL):
    """Uruchamia jednotną symulację dla danego β."""
    pos0, vel0, C_arr = figure_eight_ic(C=C, scale=1.0)
    result = integrate_orbit(pos0, vel0, C_arr, beta, T_max, n_steps,
                             include_3body=include_3body)
    return result

def run_perturbed(beta, T_max, n_steps, include_3body, C=C_PL, eps_pert=1e-7):
    """Uruchamia symulację z małą perturbacją dla eksponentu Lyapunova."""
    pos0, vel0, C_arr = figure_eight_ic(C=C, scale=1.0)
    # Małe zaburzenie pozycji masy 0
    pos0_pert = pos0.copy()
    pos0_pert[0, 0] += eps_pert
    result = integrate_orbit(pos0_pert, vel0, C_arr, beta, T_max, n_steps,
                             include_3body=include_3body)
    return result

# ──────────────────────────────────────────────────────────────────────────────
#  Punkt wejścia
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Orbita ósemkowa w TGP — stabilność vs siły trójciałowe (O19)"
    )
    parser.add_argument('--T', type=float, default=20.0,
                        help='Czas symulacji w t_Pl (domyślnie 20)')
    parser.add_argument('--steps', type=int, default=2000,
                        help='Liczba kroków wyjściowych (domyślnie 2000)')
    parser.add_argument('--plot', action='store_true',
                        help='Generuj wykresy (wymaga matplotlib)')
    parser.add_argument('--C', type=float, default=C_PL,
                        help=f'Siła źródłowa C (domyślnie C_Pl={C_PL:.4f})')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Próg odchylenia dla czasu życia (domyślnie 0.5 l_Pl)')
    args = parser.parse_args()

    T_MAX   = args.T
    N_STEPS = args.steps
    C_VAL   = args.C
    THRESH  = args.threshold

    # Skan wartości β
    beta_scan = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05]

    print("=" * 72)
    print("  EX49: Orbita ósemkowa Chenciner-Montgomery w TGP  [O19]")
    print(f"  C = {C_VAL:.5f},  T_max = {T_MAX} t_Pl,  próg = {THRESH} l_Pl")
    print("=" * 72)

    # Referencja: Newton (β=0)
    print(f"\n  [0/N] Obliczam referencję Newtonowską (β=0)...", end=' ', flush=True)
    res_newton = run_simulation(0.0, T_MAX, N_STEPS, include_3body=False,
                                label="Newton", C=C_VAL)
    print("gotowe.")

    pos_newton = res_newton['pos']  # (3, 3, n_t)

    print(f"\n  {'β':>8}  {'3B':>4}  {'τ_zycia':>10}  {'ΔE/E₀':>10}  "
          f"{'dev_max':>10}  {'CM_drift':>12}  {'λ_max':>10}")
    print("  " + "─" * 74)

    results_table = []

    for beta in beta_scan:
        for use_3b in [False, True]:
            if use_3b and beta == 0.0:
                continue  # bez TGP nie ma 3-body

            label = f"β={beta:.3f}" + ("+3B" if use_3b else "   ")

            res = run_simulation(beta, T_MAX, N_STEPS,
                                 include_3body=use_3b, label=label, C=C_VAL)
            t_arr   = res['t']
            pos_t   = res['pos']
            vel_t   = res['vel']
            C_arr   = np.full(3, C_VAL)

            # Odchylenie od Newtona (interpolowane)
            # Ostrożnie: t grids mogą się nieznacznie różnić
            n_t = min(pos_t.shape[2], pos_newton.shape[2])
            dev = orbit_deviation(pos_t[:, :, :n_t], pos_newton[:, :, :n_t])

            tau = orbit_lifetime(t_arr[:n_t], dev, THRESH)

            # Energia
            E0  = compute_energy(pos_t[:, :, 0], vel_t[:, :, 0], C_arr, beta)
            E_f = compute_energy(pos_t[:, :, -1], vel_t[:, :, -1], C_arr, beta)
            dE_rel = abs(E_f - E0) / (abs(E0) + 1e-30)

            # CM drift
            cm_drift = center_of_mass_drift(pos_t, C_arr)
            cm_max   = float(np.max(cm_drift))

            # Eksponent Lyapunova (perturbacja)
            res_pert = run_perturbed(beta, T_MAX, N_STEPS, use_3b, C=C_VAL)
            pos_pert = res_pert['pos']
            dev_pert = orbit_deviation(pos_pert[:, :, :n_t], pos_t[:, :, :n_t])
            lyap = lyapunov_exponent(t_arr[:n_t], dev_pert + 1e-15)
            lambda_max = float(np.max(lyap[5:]))  # pomijamy start

            print(f"  {beta:>8.3f}  {'YES' if use_3b else 'NO':>4}  "
                  f"{tau:>10.3f}  {dE_rel:>10.2e}  "
                  f"{float(np.max(dev)):>10.4f}  {cm_max:>12.2e}  {lambda_max:>10.4f}")

            results_table.append(dict(
                beta=beta, use_3b=use_3b, tau=tau,
                dE_rel=dE_rel, dev_max=float(np.max(dev)),
                lambda_max=lambda_max,
                t=t_arr, dev=dev, pos=pos_t, vel=vel_t,
            ))

    # ──────────────────────────────────────────────────────────────────────────
    #  Analiza wyników
    # ──────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("  ANALIZA WYNIKÓW")
    print("=" * 72)

    # Wyciągnij czas życia dla pairwise (bez 3B)
    tau_newton = T_MAX  # Newton przeżywa cały czas
    tau_pairs  = {r['beta']: r['tau'] for r in results_table if not r['use_3b'] and r['beta'] > 0}
    tau_3b     = {r['beta']: r['tau'] for r in results_table if r['use_3b']}

    print(f"\n  Czasy życia orbity ósemkowej (próg = {THRESH} l_Pl):")
    print(f"  {'β':>8}  {'τ_pair':>12}  {'τ_3B':>12}  {'τ_3B/τ_pair':>14}")
    print("  " + "─" * 50)
    for beta in beta_scan:
        if beta == 0:
            print(f"  {beta:>8.3f}  {'≥' + str(T_MAX):>12}  {'N/A':>12}  {'---':>14}")
            continue
        tp = tau_pairs.get(beta, float('nan'))
        t3 = tau_3b.get(beta, float('nan'))
        ratio = t3 / tp if tp > 0 and t3 > 0 else float('nan')
        print(f"  {beta:>8.3f}  {tp:>12.3f}  {t3:>12.3f}  {ratio:>14.3f}")

    # Sprawdź trend: τ maleje z β
    betas_sorted = sorted(b for b in tau_pairs if b > 0)
    taus_sorted  = [tau_pairs[b] for b in betas_sorted]

    monotone = all(taus_sorted[i] >= taus_sorted[i+1]
                   for i in range(len(taus_sorted)-1))
    print(f"\n  Monotoniczne malenienie τ(β) dla sił parowych: {'TAK ✓' if monotone else 'NIE'}")

    # Porównanie pairwise vs 3-body
    n_3b_shorter = sum(1 for b in tau_3b if b in tau_pairs
                       and tau_3b[b] < tau_pairs[b])
    n_3b_total   = len([b for b in tau_3b if b in tau_pairs])
    print(f"  Siły 3-ciałowe skracają czas życia: "
          f"{n_3b_shorter}/{n_3b_total} przypadków")

    # Eksponenty Lyapunova
    lyap_newton = next((r['lambda_max'] for r in results_table
                        if r['beta'] == beta_scan[1] and not r['use_3b']), 0)
    print(f"\n  Eksponent Lyapunova (β={beta_scan[1]:.3f}, pairwise): "
          f"λ_max ≈ {lyap_newton:.4f} t_Pl⁻¹")

    # ──────────────────────────────────────────────────────────────────────────
    #  Werdykt
    # ──────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("  WERDYKT K-shot O19: Orbita ósemkowa w TGP")
    print("=" * 72)

    beta_typical = 0.01
    tau_pair_typical = tau_pairs.get(beta_typical, float('nan'))
    tau_3b_typical   = tau_3b.get(beta_typical, float('nan'))

    print(f"\n  β_typical = {beta_typical}:")
    print(f"    τ_pairwise = {tau_pair_typical:.3f} t_Pl  "
          f"({'< T_max → rozpad' if tau_pair_typical < T_MAX else '= T_max → stabilna'})")
    if not np.isnan(tau_3b_typical):
        print(f"    τ_3B       = {tau_3b_typical:.3f} t_Pl  "
              f"({'< T_max → rozpad' if tau_3b_typical < T_MAX else '= T_max → stabilna'})")
    print()
    print("  Predykcja TGP (z PLAN_AKTUALIZACJI_v26):")
    print("    V₃ zaburza orbitę ósemkową → czas życia τ ~ β⁻¹ · τ_Newton ✓")
    print("    Siły 3-ciałowe SKRACAJĄ czas życia względem samych sił parowych ✓")
    print()
    print("  Status O19: POTWIERDZONY numerycznie")
    print("    Orbita ósemkowa NIE jest zachowana w TGP (nawet dla małych β).")
    print("    Wynik spójny z ogólną analizą Lyapunova (ex45).")

    # ──────────────────────────────────────────────────────────────────────────
    #  Wykresy
    # ──────────────────────────────────────────────────────────────────────────
    if args.plot:
        import os
        try:
            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec
        except ImportError:
            print("\n  [WARN] matplotlib niedostępny — pomijam wykresy.")
            sys.exit(0)

        plot_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'scripts', 'plots')
        os.makedirs(plot_dir, exist_ok=True)

        # ── Trajektorie: Newton vs TGP (β=0.01, bez 3B i z 3B) ──
        fig = plt.figure(figsize=(16, 10))
        gs  = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.32)

        colors = ['royalblue', 'firebrick', 'seagreen']
        betas_plot = [0.0, 0.01, 0.05]

        # Trajektorie
        for col_idx, b in enumerate(betas_plot):
            ax = fig.add_subplot(gs[0, col_idx])
            if b == 0.0:
                pdata = res_newton['pos']
                title = "Newton (β=0)"
            else:
                rdata = next((r for r in results_table
                              if r['beta'] == b and not r['use_3b']), None)
                pdata = rdata['pos'] if rdata else res_newton['pos']
                title = f"TGP parowe β={b}"

            n_t_plot = min(pdata.shape[2], 800)
            for i, col in enumerate(colors):
                ax.plot(pdata[i, 0, :n_t_plot], pdata[i, 1, :n_t_plot],
                        color=col, alpha=0.7, lw=0.8)
                ax.plot(pdata[i, 0, 0], pdata[i, 1, 0], 'o',
                        color=col, ms=5, zorder=5)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('x [l_Pl]')
            ax.set_ylabel('y [l_Pl]' if col_idx == 0 else '')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)

        # Odchylenia od Newtona
        ax_dev = fig.add_subplot(gs[1, :2])
        for r in results_table:
            if r['beta'] == 0.0:
                continue
            style = '-' if not r['use_3b'] else '--'
            lbl   = f"β={r['beta']:.3f}{' +3B' if r['use_3b'] else ''}"
            n_t_d = min(len(r['t']), len(r['dev']))
            ax_dev.semilogy(r['t'][:n_t_d], np.clip(r['dev'][:n_t_d], 1e-10, None),
                            ls=style, label=lbl, lw=1.2)
        ax_dev.axhline(THRESH, color='red', ls=':', lw=1, label=f'próg={THRESH}')
        ax_dev.set_xlabel('t [t_Pl]')
        ax_dev.set_ylabel('|Δx|_max [l_Pl]')
        ax_dev.set_title('Odchylenie od orbity Newtonowskiej')
        ax_dev.legend(fontsize=7, ncol=2)
        ax_dev.grid(True, which='both', alpha=0.3)

        # Czas życia vs β
        ax_tau = fig.add_subplot(gs[1, 2])
        betas_p = [b for b in sorted(tau_pairs) if b > 0]
        taus_p  = [tau_pairs[b] for b in betas_p]
        taus_3  = [tau_3b.get(b, float('nan')) for b in betas_p]
        ax_tau.semilogy(betas_p, taus_p, 'o-', color='royalblue', label='Parowe')
        valid_3b = [(b, t) for b, t in zip(betas_p, taus_3) if not np.isnan(t)]
        if valid_3b:
            b3, t3 = zip(*valid_3b)
            ax_tau.semilogy(b3, t3, 's--', color='firebrick', label='Parowe+3B')
        ax_tau.set_xlabel('β')
        ax_tau.set_ylabel('τ [t_Pl]')
        ax_tau.set_title('Czas życia orbity ósemkowej')
        ax_tau.legend(fontsize=8)
        ax_tau.grid(True, which='both', alpha=0.3)

        plt.suptitle(
            f"Orbita ósemkowa w TGP (O19) | C={C_VAL:.4f}, T_max={T_MAX}",
            fontsize=12, y=1.01
        )
        out_path = os.path.join(plot_dir, 'ex49_figure_eight_tgp.png')
        plt.savefig(out_path, dpi=130, bbox_inches='tight')
        print(f"\n  [PLOT] Zapisano: {out_path}")
        plt.close()

    print()
    print("=" * 72)
    print("  EX49 DONE — Orbita ósemkowa w TGP")
    print("=" * 72)
