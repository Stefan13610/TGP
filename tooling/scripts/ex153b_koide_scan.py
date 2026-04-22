"""
ex153b: Dodatkowe testy Koidego.
1. Wszystkie trojki z 4 zer B_tail
2. Wartosc g0^tau dajaca Q_K = 3/2 (z phi-FP)
3. Czy ta wartosc g0^tau jest blisko jakiegos zera B?
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from itertools import combinations

ALPHA = 2.0
PHI = (1 + np.sqrt(5)) / 2
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005
R_MAX = 100.0


def f_tgp(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-12))

def Vp(g):
    return g**2 * (1.0 - g)

def solve_soliton(g0, r_max=R_MAX):
    def rhs(r, y):
        g, gp = y
        g = max(g, G_BOUNCE + 1e-7)
        fg = f_tgp(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        source = Vp(g)
        cross = (ALPHA / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / (3.0 * fg)]
        return [gp, (source - cross - fg * 2.0 * gp / r) / fg]

    def ev_ghost(r, y):
        return y[0] - G_BOUNCE
    ev_ghost.terminal = True
    ev_ghost.direction = -1

    y0 = [g0, 0.0]
    r0 = 0.0
    r_all, g_all = [], []
    for _ in range(15):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_ghost,
                        dense_output=True, rtol=1e-10, atol=1e-12,
                        max_step=0.04)
        r_all.append(sol.t)
        g_all.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]
            st = sol.sol(rh)
            y0 = [st[0], -st[1]]
            r0 = rh
        else:
            break
    r = np.concatenate(r_all)
    g = np.concatenate(g_all)
    idx = np.argsort(r)
    return r[idx], g[idx]


def extract_A_tail(r, g):
    windows = [(20, 35), (25, 40), (30, 50), (25, 55)]
    A_vals = []
    for rL, rR in windows:
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 30:
            continue
        rf = r[mask]
        df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        A = np.sqrt(bc[0]**2 + bc[1]**2)
        if A > 1e-8:
            A_vals.append(A)
    return float(np.median(A_vals)) if A_vals else 0.0


def koide_QK(m1, m2, m3):
    s1 = np.sqrt(abs(m1)) + np.sqrt(abs(m2)) + np.sqrt(abs(m3))
    s0 = abs(m1) + abs(m2) + abs(m3)
    if s1 < 1e-15:
        return np.nan
    return s0 / s1**2 * 3


def run():
    print("=" * 72)
    print("ex153b: Dodatkowe testy Koidego")
    print("=" * 72)

    # Zera z ex153
    zeros = [1.225978, 2.069449, 2.760735, 3.477218]

    # --- 1. A_tail dla kazdego zera ---
    print("\n--- 1. Amplitudy zer B_tail ---")
    A_zeros = []
    for g0z in zeros:
        r, g = solve_soliton(g0z)
        A = extract_A_tail(r, g)
        A_zeros.append(A)
        print(f"  g0 = {g0z:.6f}: A = {A:.6f}")

    # --- 2. Wszystkie trojki ---
    print("\n--- 2. Q_K dla wszystkich trojek zer ---")
    for combo in combinations(range(len(zeros)), 3):
        i, j, k = combo
        m0, m1, m2 = A_zeros[i]**4, A_zeros[j]**4, A_zeros[k]**4
        QK = koide_QK(m0, m1, m2)
        r21 = m1 / m0 if m0 > 0 else 0
        r31 = m2 / m0 if m0 > 0 else 0
        mark = " <<<" if abs(QK - 1.5) < 0.05 else ""
        print(f"  Trojka ({i+1},{j+1},{k+1}): Q_K = {QK:.4f}, "
              f"r_21 = {r21:.1f}, r_31 = {r31:.1f}{mark}")

    # --- 3. phi-FP: g0_e i g0_mu ---
    print("\n--- 3. phi-FP: skan g0^tau dajacy Q_K = 3/2 ---")
    g0_e_fp = 1.249
    g0_mu_fp = PHI * g0_e_fp  # ~ 2.021

    r_e, g_e = solve_soliton(g0_e_fp)
    A_e = extract_A_tail(r_e, g_e)
    r_mu, g_mu = solve_soliton(g0_mu_fp)
    A_mu = extract_A_tail(r_mu, g_mu)

    print(f"  g0_e  = {g0_e_fp:.4f}, A_e  = {A_e:.6f}")
    print(f"  g0_mu = {g0_mu_fp:.4f}, A_mu = {A_mu:.6f}")
    print(f"  r_21  = (A_mu/A_e)^4 = {(A_mu/A_e)**4:.2f}  (PDG: 206.768)")

    # Skan g0_tau
    g0_tau_scan = np.linspace(2.0, 4.0, 100)
    QK_scan = []
    A_tau_scan = []
    print(f"\n  Skanowanie g0_tau in [{g0_tau_scan[0]:.1f}, {g0_tau_scan[-1]:.1f}]...")

    for g0t in g0_tau_scan:
        r_t, g_t = solve_soliton(g0t)
        A_t = extract_A_tail(r_t, g_t)
        A_tau_scan.append(A_t)
        QK = koide_QK(A_e**4, A_mu**4, A_t**4)
        QK_scan.append(QK)

    QK_arr = np.array(QK_scan)

    # Znajdz g0_tau gdzie Q_K = 3/2
    for i in range(len(QK_arr) - 1):
        if np.isnan(QK_arr[i]) or np.isnan(QK_arr[i+1]):
            continue
        if (QK_arr[i] - 1.5) * (QK_arr[i+1] - 1.5) < 0:
            # Interpolacja liniowa
            f1 = QK_arr[i] - 1.5
            f2 = QK_arr[i+1] - 1.5
            g0_interp = g0_tau_scan[i] - f1 * (g0_tau_scan[i+1] - g0_tau_scan[i]) / (f2 - f1)
            r_t, g_t = solve_soliton(g0_interp)
            A_t = extract_A_tail(r_t, g_t)
            QK_check = koide_QK(A_e**4, A_mu**4, A_t**4)

            print(f"\n  >>> Q_K = 3/2 przy g0_tau = {g0_interp:.4f} <<<")
            print(f"      A_tau = {A_t:.6f}")
            print(f"      Q_K   = {QK_check:.6f}")
            print(f"      r_31  = (A_tau/A_e)^4 = {(A_t/A_e)**4:.2f}")

            # Ile to od zer B?
            for iz, gz in enumerate(zeros):
                dist = abs(g0_interp - gz) / gz * 100
                print(f"      Odleglosc od zera #{iz+1} (g0={gz:.4f}): {dist:.1f}%")

            # Czy blisko phi^2 * g0_e?
            g0_phi2 = PHI**2 * g0_e_fp
            dist_phi2 = abs(g0_interp - g0_phi2) / g0_phi2 * 100
            print(f"      phi^2*g0_e = {g0_phi2:.4f}, odleglosc: {dist_phi2:.1f}%")

            # Czy blisko phi * g0_mu?
            g0_phimu = PHI * g0_mu_fp
            dist_phimu = abs(g0_interp - g0_phimu) / g0_phimu * 100
            print(f"      phi*g0_mu  = {g0_phimu:.4f}, odleglosc: {dist_phimu:.1f}%")

    # Tabela Q_K(g0_tau) co 10 krokow
    print(f"\n  Tabela Q_K(g0_tau):")
    print(f"  {'g0_tau':>8}  {'A_tau':>10}  {'Q_K':>8}  {'r_31':>10}")
    print("  " + "-" * 45)
    for i in range(0, len(g0_tau_scan), 5):
        g0t = g0_tau_scan[i]
        At = A_tau_scan[i]
        QK = QK_arr[i]
        r31 = (At / A_e)**4 if A_e > 0 else 0
        mark = " <<<" if abs(QK - 1.5) < 0.03 else ""
        print(f"  {g0t:8.4f}  {At:10.6f}  {QK:8.4f}  {r31:10.1f}{mark}")

    # --- 4. Stosunki g0 zer ---
    print("\n--- 4. Stosunki g0 zer ---")
    for i in range(len(zeros)):
        for j in range(i+1, len(zeros)):
            ratio = zeros[j] / zeros[i]
            is_phi = abs(ratio - PHI) / PHI * 100
            print(f"  g0_{j+1}/g0_{i+1} = {zeros[j]:.4f}/{zeros[i]:.4f} = {ratio:.4f}"
                  f"  (phi={PHI:.4f}, delta={is_phi:.1f}%)")

    return True


if __name__ == "__main__":
    run()
