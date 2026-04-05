"""
ex154_substrate_ode_Z3.py
==========================
R9: Test Z_3 z ODE substratowego (K_sub = g^2).

ODE substratowe:
  K_sub(g) = g^2  =>  S = integral [g^2/2 * (g')^2 + V(g)] r^2 dr
  E-L: g^2*g'' + g*(g')^2 + 2*g^2*g'/r = V'(g)
  Dzielac przez g^2:
  g'' + (1/g)*(g')^2 + (2/r)*g' = V'(g)/g^2

Kluczowe roznice od kanonicznego ODE:
  - Brak bariery duchowej (g^2 > 0 zawsze)
  - Efektywne alpha = 1 (nie 2)
  - Modyfikowane zrodlo: V'(g)/g^2 = (g^2-g^3)/g^2 = 1-g

PLAN:
1. Rozwiaz ODE substratowe dla zakresu g0
2. Wyznacz A_tail(g0) i B_tail(g0)
3. Znajdz zera B_tail
4. Sprawdz phi-FP: g0_mu = phi * g0_e
5. Test Koidego: Q_K z trzech modow
6. Sprawdz r_21 i r_31
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 120.0


def Vp_over_g2(g):
    """V'(g)/g^2 = (g^2 - g^3)/g^2 = 1 - g"""
    return 1.0 - g


def Vp(g):
    return g**2 * (1.0 - g)


def solve_substrate(g0, r_max=R_MAX):
    """
    ODE: g'' + (1/g)*(g')^2 + (2/r)*g' = V'(g)/g^2 = 1 - g
    """
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)  # brak bariery duchowej!
        source = 1.0 - g   # V'(g)/g^2
        cross = (1.0 / g) * gp**2  # alpha_eff = 1
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]

    y0 = [g0, 0.0]
    sol = solve_ivp(rhs, (0.0, r_max), y0,
                    rtol=1e-10, atol=1e-12, max_step=0.04,
                    dense_output=True)
    return sol.t, sol.y[0]


def extract_BC(r, g, rL=30, rR=80):
    """Wyznacz B i C z fitu: r*(g-1) = B*cos(r) + C*sin(r)"""
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 40:
        return 0.0, 0.0
    rf = r[mask]
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return bc[0], bc[1]


def extract_A_tail(r, g):
    B, C = extract_BC(r, g)
    return np.sqrt(B**2 + C**2)


def B_of_g0(g0):
    try:
        r, g = solve_substrate(g0)
        B, C = extract_BC(r, g)
        return B
    except Exception:
        return np.nan


def koide_QK(m1, m2, m3):
    s1 = np.sqrt(abs(m1)) + np.sqrt(abs(m2)) + np.sqrt(abs(m3))
    s0 = abs(m1) + abs(m2) + abs(m3)
    if s1 < 1e-15:
        return np.nan
    return s0 / s1**2 * 3


def cv_of(arr):
    m = np.mean(arr)
    return np.std(arr) / m if abs(m) > 1e-15 else np.nan


def run():
    print("=" * 72)
    print("ex154: ODE substratowe (K_sub=g^2) — test Z_3 i Koidego")
    print("=" * 72)
    print(f"  ODE: g'' + (1/g)(g')^2 + (2/r)g' = 1-g")
    print(f"  Brak bariery duchowej (g^2 > 0)")
    print(f"  phi = {PHI:.6f}")

    # --- 1. Skan B(g0) ---
    print("\n--- 1. Skan B(g0) ---")
    g0_scan = np.concatenate([
        np.linspace(0.50, 0.90, 20),
        np.linspace(0.90, 1.50, 40),
        np.linspace(1.50, 3.00, 30),
    ])

    B_vals = []
    A_vals = []
    C_vals = []
    print(f"  Skanowanie {len(g0_scan)} wartosci g0...", flush=True)

    for g0 in g0_scan:
        try:
            r, g = solve_substrate(g0)
            B, C = extract_BC(r, g)
            A = np.sqrt(B**2 + C**2)
        except Exception:
            B, C, A = np.nan, np.nan, np.nan
        B_vals.append(B)
        C_vals.append(C)
        A_vals.append(A)

    B_arr = np.array(B_vals)

    # Tabela
    print(f"\n  {'g0':>7}  {'B':>12}  {'C':>12}  {'A':>10}")
    print("  " + "-" * 50)
    for i in range(0, len(g0_scan), 5):
        print(f"  {g0_scan[i]:7.4f}  {B_vals[i]:12.6f}  {C_vals[i]:12.6f}  {A_vals[i]:10.6f}")

    # --- 2. Zera B(g0) ---
    print("\n--- 2. Zera B(g0) ---")
    zeros = []
    for i in range(len(g0_scan) - 1):
        if np.isnan(B_arr[i]) or np.isnan(B_arr[i+1]):
            continue
        if B_arr[i] * B_arr[i+1] < 0:
            try:
                g0z = brentq(B_of_g0, g0_scan[i], g0_scan[i+1],
                             xtol=1e-7, rtol=1e-9)
                zeros.append(g0z)
                r, g = solve_substrate(g0z)
                A = extract_A_tail(r, g)
                B, C = extract_BC(r, g)
                print(f"  Zero #{len(zeros)}: g0 = {g0z:.7f}, A = {A:.6f}, "
                      f"B = {B:.2e}")
            except Exception as e:
                print(f"  Bisekcja nieudana: {e}")

    print(f"\n  Znaleziono {len(zeros)} zer B_tail(g0)")

    # --- 3. Amplitudy i stosunki ---
    if len(zeros) >= 2:
        print("\n--- 3. Amplitudy i stosunki g0 ---")
        A_zeros = []
        for g0z in zeros:
            r, g = solve_substrate(g0z)
            A = extract_A_tail(r, g)
            A_zeros.append(A)

        for i, (gz, Az) in enumerate(zip(zeros, A_zeros)):
            print(f"  Zero #{i+1}: g0 = {gz:.6f}, A = {Az:.6f}")

        print()
        for i in range(len(zeros)):
            for j in range(i+1, len(zeros)):
                ratio = zeros[j] / zeros[i]
                phi_err = abs(ratio - PHI) / PHI * 100
                print(f"  g0_{j+1}/g0_{i+1} = {ratio:.4f} "
                      f"(phi = {PHI:.4f}, delta = {phi_err:.1f}%)")

    # --- 4. phi-FP test ---
    print("\n--- 4. phi-FP test ---")
    if len(zeros) >= 1:
        g0_e = zeros[0]
        g0_mu_fp = PHI * g0_e
        print(f"  g0_e (1. zero) = {g0_e:.6f}")
        print(f"  g0_mu = phi*g0_e = {g0_mu_fp:.6f}")

        r_e, g_e = solve_substrate(g0_e)
        A_e = extract_A_tail(r_e, g_e)
        r_mu, g_mu = solve_substrate(g0_mu_fp)
        A_mu = extract_A_tail(r_mu, g_mu)

        r21 = (A_mu / A_e)**4 if A_e > 0 else 0
        print(f"  A_e  = {A_e:.6f}")
        print(f"  A_mu = {A_mu:.6f}")
        print(f"  r_21 = (A_mu/A_e)^4 = {r21:.2f}  (PDG: 206.768)")

        # Skan g0_tau
        print(f"\n  Skan g0_tau...")
        g0_tau_range = np.linspace(g0_mu_fp * 1.01, min(g0_mu_fp * 3.0, 5.0), 80)
        best_QK = None
        best_g0tau = None

        for g0t in g0_tau_range:
            r_t, g_t = solve_substrate(g0t)
            A_t = extract_A_tail(r_t, g_t)
            if A_t < 1e-8:
                continue
            QK = koide_QK(A_e**4, A_mu**4, A_t**4)
            if best_QK is None or abs(QK - 1.5) < abs(best_QK - 1.5):
                best_QK = QK
                best_g0tau = g0t
                best_At = A_t

        if best_g0tau is not None:
            r31 = (best_At / A_e)**4
            print(f"  Najblizsze Q_K = 3/2:")
            print(f"    g0_tau  = {best_g0tau:.6f}")
            print(f"    A_tau   = {best_At:.6f}")
            print(f"    Q_K     = {best_QK:.6f}")
            print(f"    r_31    = {r31:.2f}  (PDG: 3477.15)")
            print(f"    r_21    = {r21:.2f}  (PDG: 206.768)")

            cv = cv_of(np.array([A_e**2, A_mu**2, best_At**2]))
            print(f"    CV(A^2) = {cv:.4f}  (Z_3: 1.0000)")

            # Gdzie jest g0_tau wzgledem zer?
            for iz, gz in enumerate(zeros):
                d = abs(best_g0tau - gz) / gz * 100
                print(f"    Odleglosc od zera #{iz+1}: {d:.1f}%")

            # phi-based
            g0_phi2 = PHI**2 * g0_e
            d_phi2 = abs(best_g0tau - g0_phi2) / g0_phi2 * 100
            print(f"    phi^2*g0_e = {g0_phi2:.4f}, delta = {d_phi2:.1f}%")

    # --- 5. Test Koidego z zerami B ---
    if len(zeros) >= 3:
        print("\n--- 5. Test Koidego z trojkami zer B ---")
        from itertools import combinations
        A_z = []
        for gz in zeros:
            r, g = solve_substrate(gz)
            A_z.append(extract_A_tail(r, g))

        for combo in combinations(range(len(zeros)), 3):
            i, j, k = combo
            m0, m1, m2 = A_z[i]**4, A_z[j]**4, A_z[k]**4
            QK = koide_QK(m0, m1, m2)
            r21 = m1 / m0 if m0 > 0 else 0
            r31 = m2 / m0 if m0 > 0 else 0
            mark = " <<<" if abs(QK - 1.5) < 0.05 else ""
            print(f"  ({i+1},{j+1},{k+1}): Q_K = {QK:.4f}, "
                  f"r_21 = {r21:.1f}, r_31 = {r31:.1f}{mark}")

    # --- Wnioski ---
    print(f"\n{'='*72}")
    print("WNIOSKI")
    print(f"{'='*72}")
    print("""
  ODE substratowe (K_sub=g^2) rozni sie fundamentalnie od kanonicznego:
  - Brak bariery duchowej => pelny zakres g0
  - alpha_eff = 1 (nie 2)
  - Zrodlo: V'/g^2 = 1-g (nie V' = g^2-g^3)

  Pytanie centralne: czy to ODE daje Z_3 naturalnie?
""")

    return True


if __name__ == "__main__":
    run()
