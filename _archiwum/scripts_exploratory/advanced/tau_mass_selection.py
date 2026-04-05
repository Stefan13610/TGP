#!/usr/bin/env python3
"""
tau_mass_selection.py
=====================
Systematyczne poszukiwanie zasady selekcji g_0^tau — problem O-J3.

ODE solitonu TGP (Dodatek J, eq. K-ode):
    f(g) g'' + (2/r) g' = V'(g),
    f(g) = 1 + 2*alpha*ln(g),   alpha = 2,
    V'(g) = g^2 (1 - g).

Warunki brzegowe: g(0)=g0, g'(0)=0, g(inf) -> 1.
Punkt osobliwy kinetyczny: g* = exp(-1/(2*alpha)) ~ 0.779.
Przestrzen rozwiazan: g0 in (g*, inf).

Strategia:
  1. Rozwiazanie ODE solitonu TGP
  2. Ekstrakcja A_tail(g0) z oscylacyjnego ogona g(r)-1 ~ (B cos r + C sin r)/r
  3. Wyznaczenie phi-FP: g0* takie ze (A_tail(phi*g0*)/A_tail(g0*))^4 = r_21
  4. Kandydaci na g0^tau: C2-C8 + Koide test

Teoria Generowanej Przestrzeni — Mateusz Serafin, 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, time

# --- Stale fizyczne ---
PHI = (1.0 + np.sqrt(5)) / 2.0         # zlota proporcja = 1.6180339...
ALPHA = 2.0                             # N0-2
GSTAR = np.exp(-1.0 / (2 * ALPHA))     # ~ 0.7788

R21_PDG = 206.768                       # m_mu/m_e (PDG 2024)
R31_PDG = 3477.15                       # m_tau/m_e (PDG 2024)
R32_PDG = R31_PDG / R21_PDG            # m_tau/m_mu ~ 16.817

M_E  = 0.51099895                       # MeV
M_MU = 105.6583755                      # MeV
M_TAU = 1776.86                         # MeV


# --- ODE solitonu TGP ---

def f_kinetic(g, alpha=ALPHA):
    """f(g) = 1 + 2*alpha*ln(g) — czynnik kinetyczny."""
    return 1.0 + 2.0 * alpha * np.log(max(g, 1e-30))

def Vprime(g):
    """V'(g) = g^2(1-g) — pochodna potencjalu TGP."""
    return g**2 * (1.0 - g)

def soliton_ode(r, y):
    """
    ODE solitonu TGP:
      f(g) g'' + (2/r) g' = V'(g)
      => g'' = [V'(g) - (2/r) g'] / f(g)
    """
    g, gp = y
    if g < 1e-15:
        g = 1e-15

    fg = f_kinetic(g)

    if r < 1e-10:
        # L'Hopital w r=0: (f(g0)+2) g''(0) = V'(g0)
        gpp = Vprime(g) / (fg + 2.0)
    else:
        if abs(fg) < 1e-12:
            # Blisko punktu g* — regularyzacja
            gpp = 0.0
        else:
            gpp = (Vprime(g) - 2.0 * gp / r) / fg

    return [gp, gpp]


def solve_soliton(g0, r_max=100.0, n_points=8000):
    """
    Rozwiazuje ODE solitonu z g(0)=g0, g'(0)=0.
    Zwraca (r, g, gp) lub None.
    """
    if g0 <= GSTAR + 0.001:
        return None  # ponizej bariery kinetycznej

    r_start = 1e-5
    r_eval = np.linspace(r_start, r_max, n_points)

    # Warunki poczatkowe z rozwiniecia Taylora
    fg0 = f_kinetic(g0)
    gpp0 = Vprime(g0) / (fg0 + 2.0)
    g_ini = g0 + 0.5 * gpp0 * r_start**2
    gp_ini = gpp0 * r_start

    # Event: blowup
    def event_blowup(r, y):
        return 50.0 - abs(y[0])
    event_blowup.terminal = True

    # Event: g -> g* (ghost)
    def event_ghost(r, y):
        return y[0] - (GSTAR + 0.001)
    event_ghost.terminal = True

    try:
        sol = solve_ivp(
            soliton_ode,
            [r_start, r_max], [g_ini, gp_ini],
            method='Radau', t_eval=r_eval,
            rtol=1e-9, atol=1e-11, max_step=0.5,
            events=[event_blowup, event_ghost]
        )
    except Exception:
        return None

    if len(sol.t) < 100:
        return None

    return sol.t, sol.y[0], sol.y[1]


def extract_tail(r, g, r_min=30.0, r_max=80.0):
    """
    Ekstrahuje A_tail, B_tail, C_tail z oscylacyjnego ogona:
      g(r) - 1 ~ (B cos(r) + C sin(r)) / r

    Metoda: fit liniowy do u(r)*r = B*cos(r) + C*sin(r) w regionie dalekiego pola.
    Zwraca (A_tail, B, C, n_cross).
    """
    mask = (r >= r_min) & (r <= r_max)
    r_fit = r[mask]
    u_fit = (g[mask] - 1.0)  # odchylenie od prozni

    if len(r_fit) < 20:
        return 0.0, 0.0, 0.0, 0

    # Fit: u(r)*r = B*cos(r) + C*sin(r)
    ur = u_fit * r_fit
    design = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    result = np.linalg.lstsq(design, ur, rcond=None)
    B, C = result[0]
    A_tail = np.sqrt(B**2 + C**2)

    # Przejscia przez zero u(r) (= g(r)-1)
    sign_changes = np.sum(np.abs(np.diff(np.sign(u_fit))) > 1)
    n_cross = int(sign_changes)

    return A_tail, B, C, n_cross


def compute_Atail_scan(g0_range, r_max=100.0):
    """Oblicza A_tail(g0) i B_tail(g0) dla tablicy wartosci g0."""
    N = len(g0_range)
    A_tails = np.full(N, np.nan)
    B_tails = np.full(N, np.nan)
    C_tails = np.full(N, np.nan)
    n_crosses = np.zeros(N, dtype=int)

    t0 = time.time()
    for i, g0 in enumerate(g0_range):
        if (i+1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"    {i+1}/{N}  ({elapsed:.1f}s)", flush=True)

        result = solve_soliton(g0, r_max=r_max)
        if result is None:
            continue
        r, g, gp = result
        A, B, C, nc = extract_tail(r, g)
        if A > 0:
            A_tails[i] = A
            B_tails[i] = B
            C_tails[i] = C
            n_crosses[i] = nc

    elapsed = time.time() - t0
    print(f"    Skan zakonczony: {elapsed:.1f}s", flush=True)
    return A_tails, B_tails, C_tails, n_crosses


# --- Zasady selekcji ---

def find_phi_FP(g0_range, A_tails, target_r21=R21_PDG):
    """Znajduje phi-FP: g0* takie ze (A_tail(phi*g0*)/A_tail(g0*))^4 = r_21."""
    valid = ~np.isnan(A_tails) & (A_tails > 0)
    if np.sum(valid) < 20:
        return None, None

    f_Atail = interp1d(g0_range[valid], A_tails[valid],
                       kind='cubic', fill_value=np.nan, bounds_error=False)

    g0_max_fp = np.max(g0_range[valid]) / PHI - 0.01
    g0_min_fp = np.min(g0_range[valid]) + 0.05

    def residual(g0):
        A_e = f_Atail(g0)
        A_mu = f_Atail(PHI * g0)
        if np.isnan(A_e) or np.isnan(A_mu) or A_e <= 0:
            return 1e10
        return (A_mu / A_e)**4 - target_r21

    g0_test = np.linspace(g0_min_fp, g0_max_fp, 1000)
    residuals = np.array([residual(g) for g in g0_test])

    # Znajdz przejscie przez zero
    finite = np.isfinite(residuals)
    sign_changes = np.where(np.diff(np.sign(residuals[finite])))[0]

    if len(sign_changes) == 0:
        return None, f_Atail

    g0_finite = g0_test[finite]
    idx = sign_changes[0]
    g0_star = brentq(residual, g0_finite[idx], g0_finite[idx+1])

    return g0_star, f_Atail


def find_B_zeros(g0_range, B_tails):
    """Znajduje zera B_tail(g0) = 0."""
    valid = ~np.isnan(B_tails)
    g0v = g0_range[valid]
    Bv = B_tails[valid]

    zeros = []
    for i in range(len(Bv)-1):
        if Bv[i] * Bv[i+1] < 0:
            # Interpolacja liniowa
            g0_zero = g0v[i] - Bv[i] * (g0v[i+1] - g0v[i]) / (Bv[i+1] - Bv[i])
            zeros.append(g0_zero)
    return zeros


def find_tau_candidates(g0_star, f_Atail, g0_range, A_tails, B_tails, n_crosses):
    """Generuje kandydatow na g0^tau."""
    if g0_star is None:
        print("  [WARN] phi-FP nie znaleziony")
        return []

    g0_e = g0_star
    g0_mu = PHI * g0_star
    A_e = f_Atail(g0_e)

    candidates = []

    # C2: g0^mu * phi
    g0_c2 = g0_mu * PHI
    A_c2 = f_Atail(g0_c2)
    if not np.isnan(A_c2) and A_c2 > 0:
        r31 = (A_c2 / A_e)**4
        candidates.append({
            'name': 'C2: g0^mu*phi', 'g0': g0_c2,
            'r31': r31, 'delta_pct': 100*(r31/R31_PDG - 1)
        })

    # C3: g0^e * phi^2
    g0_c3 = g0_e * PHI**2
    A_c3 = f_Atail(g0_c3)
    if not np.isnan(A_c3) and A_c3 > 0:
        r31 = (A_c3 / A_e)**4
        candidates.append({
            'name': 'C3: g0^e*phi^2', 'g0': g0_c3,
            'r31': r31, 'delta_pct': 100*(r31/R31_PDG - 1)
        })

    # C4: numeryczne odwrocenie r_31
    def res_c4(g0):
        A_tau = f_Atail(g0)
        if np.isnan(A_tau) or A_tau <= 0:
            return 1e10
        return (A_tau / A_e)**4 - R31_PDG

    g0_search = np.linspace(g0_mu + 0.01, max(g0_range) - 0.2, 500)
    res_vals = np.array([res_c4(g) for g in g0_search])
    finite = np.isfinite(res_vals) & (np.abs(res_vals) < 1e9)
    if np.sum(finite) > 1:
        sc = np.where(np.diff(np.sign(res_vals[finite])))[0]
        if len(sc) > 0:
            g0f = g0_search[finite]
            g0_c4 = brentq(res_c4, g0f[sc[0]], g0f[sc[0]+1])
            candidates.append({
                'name': 'C4: numeryczne', 'g0': g0_c4,
                'r31': R31_PDG, 'delta_pct': 0.0
            })

    # C6: zera B_tail (warunek kwantowania)
    B_zeros = find_B_zeros(g0_range, B_tails)
    for iz, g0z in enumerate(B_zeros):
        A_z = f_Atail(g0z)
        if not np.isnan(A_z) and A_z > 0 and g0z > g0_mu + 0.05:
            r31 = (A_z / A_e)**4
            candidates.append({
                'name': f'C6: B_tail=0 (#{iz+1})', 'g0': g0z,
                'r31': r31, 'delta_pct': 100*(r31/R31_PDG - 1)
            })

    # C8: granice sektorow topologicznych (zmiana n_cross)
    valid = ~np.isnan(A_tails)
    nc_v = n_crosses[valid]
    g0_v = g0_range[valid]
    for k in range(len(nc_v)-1):
        if nc_v[k] != nc_v[k+1]:
            g0_b = 0.5*(g0_v[k] + g0_v[k+1])
            if g0_b > g0_mu + 0.05:
                A_b = f_Atail(g0_b)
                if not np.isnan(A_b) and A_b > 0:
                    r31 = (A_b / A_e)**4
                    candidates.append({
                        'name': f'C8: sector {nc_v[k]}->{nc_v[k+1]}',
                        'g0': g0_b,
                        'r31': r31, 'delta_pct': 100*(r31/R31_PDG - 1)
                    })

    # Koide check dla kazdego kandydata
    for c in candidates:
        m_tau_pred = M_E * c['r31']
        koide_num = M_E + M_MU + m_tau_pred
        koide_den = (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(m_tau_pred))**2
        koide = koide_num / koide_den
        c['koide'] = koide
        c['koide_delta_pct'] = 100 * abs(koide - 2.0/3.0) / (2.0/3.0)
        c['m_tau'] = m_tau_pred

    return candidates


# --- Wizualizacja ---

def plot_results(g0_range, A_tails, B_tails, n_crosses, g0_star, candidates, outdir='.'):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('TGP: tau Mass Selection Problem O-J3\n'
                 f'ODE: f(g)g" + (2/r)g\' = g^2(1-g),  f(g)=1+4ln(g),  g*={GSTAR:.4f}',
                 fontsize=13, fontweight='bold')

    valid = ~np.isnan(A_tails) & (A_tails > 0)

    # 1. A_tail(g0)
    ax = axes[0, 0]
    ax.semilogy(g0_range[valid], A_tails[valid], 'b-', lw=1.5)
    if g0_star is not None:
        ax.axvline(g0_star, color='g', ls='--', lw=1.5,
                   label=f'g0* = {g0_star:.4f} (e)')
        ax.axvline(PHI*g0_star, color='orange', ls='--', lw=1.5,
                   label=f'phi*g0* = {PHI*g0_star:.4f} (mu)')
    for c in candidates:
        if 'C4' in c['name']:
            ax.axvline(c['g0'], color='r', ls=':', lw=2,
                       label=f"g0^tau = {c['g0']:.3f} (C4)")
    ax.axvline(GSTAR, color='gray', ls='-.', alpha=0.5, label=f'g* = {GSTAR:.3f}')
    ax.set_xlabel('g0')
    ax.set_ylabel('A_tail')
    ax.set_title('A_tail(g0)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. B_tail(g0) — zera = warunek kwantowania
    ax = axes[0, 1]
    valid_B = ~np.isnan(B_tails)
    if np.sum(valid_B) > 0:
        ax.plot(g0_range[valid_B], B_tails[valid_B], 'b-', lw=1)
        ax.axhline(0, color='k', ls='-', lw=0.5)
        B_zeros = find_B_zeros(g0_range, B_tails)
        for iz, gz in enumerate(B_zeros):
            ax.axvline(gz, color='r', ls=':', alpha=0.7,
                       label=f'B=0 #{iz+1}: g0={gz:.3f}')
    ax.set_xlabel('g0')
    ax.set_ylabel('B_tail')
    ax.set_title('B_tail(g0) — zera = warunek kwantowania H1')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 3. Stosunek mas
    ax = axes[1, 0]
    if g0_star is not None:
        f_A = interp1d(g0_range[valid], A_tails[valid], kind='cubic',
                       fill_value=np.nan, bounds_error=False)
        A_e = f_A(g0_star)
        mass_ratio = (A_tails[valid] / A_e)**4
        ax.semilogy(g0_range[valid], mass_ratio, 'b-', lw=1.5)
        ax.axhline(R21_PDG, color='orange', ls='--',
                   label=f'r_21 = {R21_PDG:.1f} (mu/e)')
        ax.axhline(R31_PDG, color='r', ls='--',
                   label=f'r_31 = {R31_PDG:.1f} (tau/e)')
        for c in candidates:
            marker = 's' if 'C4' in c['name'] else 'o'
            ax.plot(c['g0'], c['r31'], marker, ms=8,
                    label=f"{c['name']}: d={c['delta_pct']:.1f}%")
    ax.set_xlabel('g0')
    ax.set_ylabel('(A_tail/A_e)^4 = m/m_e')
    ax.set_title('Stosunek mas vs g0')
    ax.legend(fontsize=7, loc='upper left')
    ax.grid(True, alpha=0.3)

    # 4. Tabela
    ax = axes[1, 1]
    ax.axis('off')
    if candidates:
        col_labels = ['Kandydat', 'g0^tau', 'm_tau [MeV]',
                      'dr31 [%]', 'Koide', 'dK [%]']
        table_data = []
        for c in candidates:
            table_data.append([
                c['name'],
                f"{c['g0']:.4f}",
                f"{c['m_tau']:.1f}",
                f"{c['delta_pct']:.1f}",
                f"{c['koide']:.6f}",
                f"{c['koide_delta_pct']:.3f}"
            ])
        table = ax.table(cellText=table_data, colLabels=col_labels,
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.4)
        ax.set_title('Kandydaci tau + Koide', fontsize=12, pad=20)
    else:
        ax.text(0.5, 0.5, 'Brak kandydatow', ha='center', va='center',
                fontsize=14, color='red')

    plt.tight_layout()
    outpath = os.path.join(outdir, 'tau_mass_selection.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  -> Zapisano: {outpath}")
    plt.close()


# --- MAIN ---

if __name__ == '__main__':
    print("=" * 70)
    print("  TGP: tau Mass Selection -- Problem O-J3")
    print("  ODE: f(g)g'' + (2/r)g' = g^2(1-g)")
    print(f"  f(g) = 1 + 2*{ALPHA}*ln(g),   g* = {GSTAR:.5f}")
    print("=" * 70)

    # --- Quick test ---
    print("\n  Quick test: 5 punktow...")
    for g0 in [0.9, 1.24, 1.5, 2.0, 3.0]:
        result = solve_soliton(g0, r_max=100.0, n_points=5000)
        if result is None:
            print(f"    g0={g0:.2f}: FAILED")
        else:
            r, g, gp = result
            A, B, C, nc = extract_tail(r, g)
            print(f"    g0={g0:.2f}: rmax={r[-1]:.1f}, "
                  f"A_tail={A:.4f}, B={B:.4f}, C={C:.4f}, "
                  f"n_cross={nc}, g_end={g[-1]:.6f}")

    # --- Pelny skan ---
    g0_range = np.linspace(GSTAR + 0.05, 4.0, 600)

    print(f"\n  Skan: g0 in [{g0_range[0]:.3f}, {g0_range[-1]:.3f}], "
          f"{len(g0_range)} punktow")
    print("  Rozwiazywanie ODE solitonu TGP...\n")

    A_tails, B_tails, C_tails, n_crosses = compute_Atail_scan(g0_range, r_max=100.0)

    n_valid = np.sum(~np.isnan(A_tails))
    print(f"\n  Uzyskano {n_valid}/{len(g0_range)} waznych rozwiazan")

    # --- phi-FP ---
    print("\n  Szukam phi-fixed point...")
    g0_star, f_Atail = find_phi_FP(g0_range, A_tails)

    if g0_star is not None:
        A_e = f_Atail(g0_star)
        A_mu = f_Atail(PHI * g0_star)
        r21_check = (A_mu / A_e)**4

        print(f"\n  phi-FP znaleziony:")
        print(f"    g0*      = {g0_star:.5f}   (g0^e)")
        print(f"    g0^mu    = phi*g0* = {PHI*g0_star:.5f}")
        print(f"    A_e      = {A_e:.5f}")
        print(f"    A_mu     = {A_mu:.5f}")
        print(f"    (A_mu/A_e)^4 = {r21_check:.3f}  (PDG: {R21_PDG})")
        print(f"    dr_21    = {100*(r21_check/R21_PDG - 1):.4f}%")
        print(f"    g0^e/g*  = {g0_star/GSTAR:.4f}")

        # Porownanie z danymi teoretycznymi
        print(f"\n  Porownanie z teoria (Dod. J):")
        print(f"    g0^e teoria: 1.24,    tu: {g0_star:.4f}  "
              f"(d = {100*(g0_star/1.24-1):.1f}%)")
        print(f"    g0^mu teoria: 2.00,   tu: {PHI*g0_star:.4f}  "
              f"(d = {100*(PHI*g0_star/2.00-1):.1f}%)")
    else:
        print("  [WARN] phi-FP NIE ZNALEZIONY!")

    # --- Zera B_tail ---
    print("\n  Zera B_tail(g0):")
    B_zeros = find_B_zeros(g0_range, B_tails)
    for iz, gz in enumerate(B_zeros):
        A_z = f_Atail(gz) if f_Atail is not None else 0
        print(f"    #{iz+1}: g0 = {gz:.4f}   A_tail = {A_z:.5f}")

    # --- Kandydaci tau ---
    print("\n  Generuje kandydatow na g0^tau...")
    candidates = find_tau_candidates(
        g0_star, f_Atail, g0_range, A_tails, B_tails, n_crosses)

    # --- Wyniki ---
    print("\n" + "=" * 80)
    print("  WYNIKI -- Kandydaci na g0^tau")
    print("=" * 80)
    header = (f"  {'Kandydat':<25s} {'g0^tau':>8s} {'r_31':>10s} "
              f"{'dr31[%]':>9s} {'m_tau[MeV]':>10s} {'Koide':>8s} {'dK[%]':>7s}")
    print(header)
    print("  " + "-" * 78)

    for c in candidates:
        print(f"  {c['name']:<25s} {c['g0']:8.4f} {c['r31']:10.1f} "
              f"{c['delta_pct']:9.1f} {c['m_tau']:10.1f} "
              f"{c['koide']:8.6f} {c['koide_delta_pct']:7.3f}")

    if candidates:
        best = min(candidates, key=lambda c: abs(c['delta_pct']))
        print(f"\n  -> Najlepszy: {best['name']}  (g0^tau = {best['g0']:.4f}, "
              f"dr_31 = {best['delta_pct']:.1f}%)")

        c4 = [c for c in candidates if 'C4' in c['name']]
        if c4:
            c4 = c4[0]
            print(f"\n  Predykcja bazowa (C4): g0^tau = {c4['g0']:.4f}")
            print(f"    m_tau    = {c4['m_tau']:.2f} MeV  (PDG: {M_TAU} MeV)")
            if g0_star:
                print(f"    g0^tau/g0^mu = {c4['g0']/(PHI*g0_star):.4f}")
                print(f"    g0^tau/g0^e  = {c4['g0']/g0_star:.4f}")
            print(f"    Koide    = {c4['koide']:.6f}  (exact: 0.666667)")

    # Sektory topologiczne
    print("\n  Sektory topologiczne:")
    valid = ~np.isnan(A_tails)
    unique_nc = np.unique(n_crosses[valid])
    for nc in unique_nc:
        mask = (n_crosses == nc) & valid
        if np.sum(mask) > 0:
            g_min, g_max = g0_range[mask].min(), g0_range[mask].max()
            print(f"    n_cross = {nc}: g0 in [{g_min:.3f}, {g_max:.3f}]  "
                  f"({np.sum(mask)} pts)")

    # --- Wykresy ---
    outdir = os.path.dirname(os.path.abspath(__file__))
    plot_results(g0_range, A_tails, B_tails, n_crosses, g0_star, candidates, outdir)

    print(f"\n  DONE.\n")
