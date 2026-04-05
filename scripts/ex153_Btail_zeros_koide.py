"""
ex153_Btail_zeros_koide.py
===========================
R9: Czy warunek kwantowania B_tail(g0) = 0 daje Z_3?

IDEA:
Soliton ma ogon: u(r) = (B*cos(r) + C*sin(r))/r
A_tail = sqrt(B^2 + C^2), faza = atan2(B, C)
Warunek kwantowania: B(g0) = 0 (brak fali wychodzacej)
=> ogon to czysta fala stojaca sin(r)/r

Jesli kolejne zera B(g0) daja g0^(0), g0^(1), g0^(2),
to odpowiadajace A_tail^(0), A_tail^(1), A_tail^(2)
sa WYZNACZONE.

PYTANIE: Czy CV(A_n^2) = 1 (rownowazne Q_K = 3/2)?

Jesli TAK -> Z_3 wynika z warunku kwantowania
           -> Koide jest KONSEKWENCJA, nie postulatem!
Jesli NIE -> Z_3 pozostaje postulatem

METODA:
1. Gesty skan g0 in [1.01, 4.0]
2. Dla kazdego g0: oblicz B, C z fitu ogona
3. Znajdz zera B(g0) (zmiana znaku)
4. Dokladnie zlokalizuj zera (bisekcja)
5. Dla kazdego trojkatu zer: oblicz Q_K
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

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


def extract_BC(r, g, rL=25, rR=70):
    """Wyznacz B i C z fitu: r*(g-1) = B*cos(r) + C*sin(r)"""
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 40:
        return 0.0, 0.0
    rf = r[mask]
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return bc[0], bc[1]  # B, C


def extract_A_tail(r, g):
    B, C = extract_BC(r, g)
    return np.sqrt(B**2 + C**2)


def B_of_g0(g0):
    """Zwraca B(g0) — skladowa kosinusowa ogona."""
    try:
        r, g = solve_soliton(g0)
        B, C = extract_BC(r, g)
        return B
    except Exception:
        return np.nan


def koide_QK(m1, m2, m3):
    """Q_K = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2 * 3"""
    s1 = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    s0 = m1 + m2 + m3
    if s1 < 1e-15:
        return np.nan
    return s0 / s1**2 * 3


def cv_of_arr(arr):
    """CV = std/mean"""
    m = np.mean(arr)
    if abs(m) < 1e-15:
        return np.nan
    return np.std(arr) / m


def run():
    print("=" * 72)
    print("ex153: Zera B_tail(g0) i test Koidego")
    print("=" * 72)

    # --- 1. Gesty skan B(g0) ---
    print("\n--- 1. Skan B(g0) ---")
    g0_scan = np.concatenate([
        np.linspace(1.01, 1.10, 20),
        np.linspace(1.10, 1.40, 40),
        np.linspace(1.40, 1.65, 30),  # near ghost barrier limit
    ])

    B_vals = []
    A_vals = []
    print(f"  Skanowanie {len(g0_scan)} wartosci g0...", flush=True)
    for g0 in g0_scan:
        try:
            r, g = solve_soliton(g0)
            B, C = extract_BC(r, g)
            A = np.sqrt(B**2 + C**2)
            B_vals.append(B)
            A_vals.append(A)
        except Exception:
            B_vals.append(np.nan)
            A_vals.append(np.nan)

    B_arr = np.array(B_vals)
    A_arr = np.array(A_vals)

    # Wypisz tabele
    print(f"\n  {'g0':>7}  {'B':>12}  {'C':>12}  {'A':>10}  {'B/A':>8}")
    print("  " + "-" * 55)
    for i in range(0, len(g0_scan), 5):
        g0 = g0_scan[i]
        B = B_arr[i]
        A = A_arr[i]
        ratio = B / A if A > 1e-8 else 0
        r, g = solve_soliton(g0)
        _, C = extract_BC(r, g)
        print(f"  {g0:7.4f}  {B:12.6f}  {C:12.6f}  {A:10.6f}  {ratio:8.4f}")

    # --- 2. Znajdz zera B(g0) ---
    print("\n--- 2. Zera B(g0) ---")
    zeros = []
    for i in range(len(g0_scan) - 1):
        if np.isnan(B_arr[i]) or np.isnan(B_arr[i+1]):
            continue
        if B_arr[i] * B_arr[i+1] < 0:
            # Bisekcja
            try:
                g0_zero = brentq(B_of_g0, g0_scan[i], g0_scan[i+1],
                                 xtol=1e-6, rtol=1e-8)
                zeros.append(g0_zero)
                print(f"  Zero #{len(zeros)}: g0 = {g0_zero:.6f}")
            except Exception as e:
                print(f"  Bisekcja nieudana miedzy {g0_scan[i]:.4f} "
                      f"i {g0_scan[i+1]:.4f}: {e}")

    print(f"\n  Znaleziono {len(zeros)} zer B_tail(g0)")

    if len(zeros) == 0:
        print("  BRAK ZER — nie mozna kontynuowac")
        # Sprawdzmy znaki B
        print("\n  Znaki B(g0):")
        for i in range(0, len(g0_scan), 3):
            sign = "+" if B_arr[i] > 0 else ("-" if B_arr[i] < 0 else "0")
            print(f"    g0={g0_scan[i]:.4f}: B={B_arr[i]:.6f} ({sign})")
        return

    # --- 3. Oblicz A_tail dla zer ---
    print("\n--- 3. Amplitudy przy zerach B ---")
    A_at_zeros = []
    for i, g0z in enumerate(zeros):
        r, g = solve_soliton(g0z)
        A = extract_A_tail(r, g)
        B, C = extract_BC(r, g)
        A_at_zeros.append(A)
        print(f"  Zero #{i+1}: g0 = {g0z:.6f}, A = {A:.6f}, "
              f"B = {B:.2e} (should be ~0), C = {C:.6f}")

    # --- 4. Stosunki mas i Koide ---
    if len(zeros) >= 2:
        print("\n--- 4. Stosunki mas z par zer ---")
        for i in range(len(zeros)):
            for j in range(i+1, len(zeros)):
                ratio = (A_at_zeros[j] / A_at_zeros[i])**4
                print(f"  r({j+1},{i+1}) = (A_{j+1}/A_{i+1})^4 = {ratio:.2f}")

    if len(zeros) >= 3:
        print("\n--- 5. TEST KOIDEGO ---")
        for i in range(len(zeros) - 2):
            A0 = A_at_zeros[i]
            A1 = A_at_zeros[i+1]
            A2 = A_at_zeros[i+2]
            # m_n ~ A_n^4
            m0, m1, m2 = A0**4, A1**4, A2**4
            QK = koide_QK(m0, m1, m2)
            r21 = m1 / m0
            r31 = m2 / m0

            # CV(A^2)
            cv = cv_of_arr(np.array([A0**2, A1**2, A2**2]))

            print(f"  Trojka zer #{i+1},{i+2},{i+3}:")
            print(f"    g0 = ({zeros[i]:.6f}, {zeros[i+1]:.6f}, {zeros[i+2]:.6f})")
            print(f"    A  = ({A0:.6f}, {A1:.6f}, {A2:.6f})")
            print(f"    m/m_0 = (1, {r21:.2f}, {r31:.2f})")
            print(f"    Q_K   = {QK:.6f}  (PDG: 1.500000)")
            print(f"    CV(A^2) = {cv:.4f}  (Z_3: 1.0000)")
            print(f"    r_21  = {r21:.2f}  (PDG: 206.768)")
            print(f"    r_31  = {r31:.2f}  (PDG: 3477.15)")

            if abs(QK - 1.5) < 0.05:
                print(f"    >>> Q_K BLISKO 3/2! Z_3 MOZE WYNIKAC Z KWANTOWANIA <<<")
            else:
                print(f"    >>> Q_K daleko od 3/2 — Z_3 NIE wynika z kwantowania B <<<")

    # --- 6. Rozszerzony skan (do g0 = 4) jesli brakuje zer ---
    if len(zeros) < 3:
        print("\n--- 6. Rozszerzony skan g0 in [1.6, 4.0] ---")
        g0_ext = np.linspace(1.60, 4.0, 60)
        print(f"  Skanowanie {len(g0_ext)} dodatkowych wartosci...", flush=True)
        B_ext = []
        for g0 in g0_ext:
            try:
                r, g = solve_soliton(g0)
                B, C = extract_BC(r, g)
                B_ext.append(B)
            except Exception:
                B_ext.append(np.nan)

        B_ext = np.array(B_ext)
        for i in range(0, len(g0_ext), 4):
            sign = "+" if B_ext[i] > 0 else ("-" if B_ext[i] < 0 else "?")
            bval = B_ext[i] if not np.isnan(B_ext[i]) else 0
            print(f"    g0={g0_ext[i]:.4f}: B={bval:.6f} ({sign})")

        # Szukaj zer
        for i in range(len(g0_ext) - 1):
            if np.isnan(B_ext[i]) or np.isnan(B_ext[i+1]):
                continue
            if B_ext[i] * B_ext[i+1] < 0:
                try:
                    g0z = brentq(B_of_g0, g0_ext[i], g0_ext[i+1],
                                 xtol=1e-5, rtol=1e-7)
                    zeros.append(g0z)
                    r, g = solve_soliton(g0z)
                    A = extract_A_tail(r, g)
                    A_at_zeros.append(A)
                    print(f"  NOWE ZERO #{len(zeros)}: g0 = {g0z:.6f}, A = {A:.6f}")
                except Exception:
                    pass

        # Ponowny test Koidego jesli mamy 3 zera
        if len(zeros) >= 3:
            print("\n--- 7. TEST KOIDEGO (z rozszerzonym skanem) ---")
            A0 = A_at_zeros[0]
            A1 = A_at_zeros[1]
            A2 = A_at_zeros[2]
            m0, m1, m2 = A0**4, A1**4, A2**4
            QK = koide_QK(m0, m1, m2)
            r21 = m1 / m0
            r31 = m2 / m0
            cv = cv_of_arr(np.array([A0**2, A1**2, A2**2]))

            print(f"  g0 = ({zeros[0]:.6f}, {zeros[1]:.6f}, {zeros[2]:.6f})")
            print(f"  A  = ({A0:.6f}, {A1:.6f}, {A2:.6f})")
            print(f"  Q_K = {QK:.6f}  (PDG: 1.500000)")
            print(f"  CV(A^2) = {cv:.4f}  (Z_3: 1.0000)")
            print(f"  r_21 = {r21:.2f}  (PDG: 206.768)")
            print(f"  r_31 = {r31:.2f}  (PDG: 3477.15)")

    # --- Wnioski ---
    print(f"\n--- WNIOSKI ---")
    print(f"  Znaleziono {len(zeros)} zer B_tail(g0)")
    if len(zeros) >= 3:
        print(f"  Test Koidego: Q_K = {koide_QK(A_at_zeros[0]**4, A_at_zeros[1]**4, A_at_zeros[2]**4):.4f}")
    print("""
  Jesli Q_K ~ 3/2: Z_3 WYNIKA z warunku kwantowania B=0
  => Koide jest konsekwencja dynamiki solitonowej, nie postulatem
  => R9 ZAMKNIETY

  Jesli Q_K != 3/2: Z_3 NIE wynika z samego kwantowania B=0
  => potrzebny dodatkowy warunek selekcji
  => R9 pozostaje otwarty
""")

    return True


if __name__ == "__main__":
    run()
