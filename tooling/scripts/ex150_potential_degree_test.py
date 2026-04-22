"""
ex150_potential_degree_test.py
==============================
Test argumentu prop:K-exponent: wykladnik skalowania A_tail(g0)
wynika ze stopnia potencjalu V(g).

Hipoteza (dodatek K, prop:K-exponent):
  A_tail ~ (g0-1)^deg(V)  bo:
  - energia solitonu E_sol ~ (g0-1)^2 (z V''(1) ~ delta^2)
  - amplituda ogona A ~ E^2  (promieniowanie kwadratowe)
  - wiec A ~ (g0-1)^4  gdy deg(V) = 4

Test: rozwiazujemy ODE TGP z roznymi potencjalami V(g):
  V_n(g) = g^3/3 - g^n/n   dla n = 3, 4, 5, 6
i fitujemy A_tail(g0) ~ (g0-1)^k.  Czy k ~ n?

Przewidywanie argumentu heurystycznego:
  k(n=3) ~ 3,  k(n=4) ~ 4,  k(n=5) ~ 5,  k(n=6) ~ 6
Albo jesli A ~ E^2 i E ~ (g0-1)^(n/2):
  k(n) ~ n
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp

ALPHA = 2.0
R_MAX = 80.0


def f_tgp(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-12))


G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))


def make_potential(n_deg):
    """Potencjal V_n(g) = g^3/3 - g^n/n, probka g=1."""
    def V(g):
        return g**3 / 3.0 - g**n_deg / n_deg

    def Vp(g):
        return g**2 - g**(n_deg - 1)

    return V, Vp


def solve_soliton(g0, Vp_func, r_max=R_MAX):
    """Rozwiaz ODE TGP z danym V'(g)."""
    g_bounce = G_GHOST + 0.005

    def rhs(r, y):
        g, gp = y
        g = max(g, g_bounce + 1e-7)
        fg = f_tgp(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        source = Vp_func(g)
        cross = (ALPHA / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / (3.0 * fg)]
        return [gp, (source - cross - fg * 2.0 * gp / r) / fg]

    def ev_ghost(r, y):
        return y[0] - g_bounce
    ev_ghost.terminal = True
    ev_ghost.direction = -1

    y0 = [g0, 0.0]
    r0 = 0.0
    r_all, g_all = [], []
    for _ in range(10):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_ghost,
                        dense_output=True, rtol=1e-9, atol=1e-11,
                        max_step=0.05)
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
    """Ekstrakcja amplitudy ogona z fitu sin/cos."""
    windows = [(16, 28), (20, 32), (24, 36), (28, 40)]
    A_vals = []
    for rL, rR in windows:
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 20:
            continue
        rf = r[mask]
        df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        A = np.sqrt(bc[0]**2 + bc[1]**2)
        if A > 1e-8:
            A_vals.append(A)
    return float(np.median(A_vals)) if A_vals else 0.0


def scan_and_fit(n_deg, g0_vals):
    """Skan g0, oblicz A_tail, fituj wykladnik."""
    _, Vp = make_potential(n_deg)

    data = []
    for g0 in g0_vals:
        try:
            r, g = solve_soliton(g0, Vp)
            A = extract_A_tail(r, g)
            if A > 1e-6:
                data.append((g0, A))
        except Exception:
            pass

    if len(data) < 5:
        return None, None, data

    g0_arr = np.array([d[0] for d in data])
    A_arr = np.array([d[1] for d in data])

    # Fit log(A) vs log(g0 - g*) AND vs log(g0 - 1)
    delta_gstar = g0_arr - G_GHOST
    delta_one = g0_arr - 1.0

    results_fit = {}
    for label, delta in [("g0-g*", delta_gstar), ("g0-1", delta_one)]:
        valid = delta > 0.01
        if np.sum(valid) < 5:
            results_fit[label] = (None, None)
            continue
        logd = np.log(delta[valid])
        logA = np.log(A_arr[valid])
        k, logc = np.polyfit(logd, logA, 1)
        res = np.std(logA - (logc + k * logd))
        results_fit[label] = (k, res)

    return results_fit, data


def run():
    print("=" * 72)
    print("ex150: Test argumentu stopnia potencjalu -> wykladnik A_tail")
    print("=" * 72)

    g0_vals = np.concatenate([
        np.linspace(1.02, 1.10, 5),
        np.linspace(1.10, 1.30, 12),
        np.linspace(1.30, 1.60, 6),
        np.linspace(1.60, 2.50, 10),
    ])

    n_degs = [3, 4, 5, 6]

    print(f"\n{'n_deg':>5}  {'k(g0-g*)':>10}  {'res':>7}  {'k(g0-1)':>10}  {'res':>7}  {'N_pts':>5}")
    print("-" * 65)

    results = {}
    for n in n_degs:
        print(f"\n  Obliczanie n_deg = {n}...", flush=True)
        fits, data = scan_and_fit(n, g0_vals)
        results[n] = (fits, data)

        k_gs = fits.get("g0-g*", (None, None))
        k_one = fits.get("g0-1", (None, None))
        k_gs_v = k_gs[0] if k_gs[0] is not None else float('nan')
        k_one_v = k_one[0] if k_one[0] is not None else float('nan')
        r_gs = k_gs[1] if k_gs[1] is not None else float('nan')
        r_one = k_one[1] if k_one[1] is not None else float('nan')
        print(f"  {n:5d}  {k_gs_v:10.3f}  {r_gs:7.4f}  {k_one_v:10.3f}  {r_one:7.4f}  {len(data):5d}")

    # --- Podsumowanie ---
    print("\n" + "=" * 72)
    print("PODSUMOWANIE")
    print("=" * 72)

    print(f"\n{'n (stopien V)':>15}  {'k(g0-g*)':>10}  {'k(g0-1)':>10}  {'k(g0-g*)/n':>12}")
    print("-" * 55)
    k_gs_all = []
    n_all = []
    for n in n_degs:
        fits, _ = results[n]
        k_gs = fits.get("g0-g*", (None, None))[0]
        k_one = fits.get("g0-1", (None, None))[0]
        if k_gs is not None:
            k_gs_all.append(k_gs)
            n_all.append(n)
            ratio = k_gs / n
            print(f"  {n:13d}  {k_gs:10.3f}  {k_one:10.3f}  {ratio:12.3f}")
        else:
            print(f"  {n:13d}  {'---':>10}  {'---':>10}  {'---':>12}")

    if len(k_gs_all) >= 3:
        slope, intercept = np.polyfit(n_all, k_gs_all, 1)
        print(f"\n  Fit liniowy k(n) = {slope:.3f} * n + {intercept:.3f}")
        print(f"  (idealnie: slope=1, intercept=0)")
        corr = np.corrcoef(n_all, k_gs_all)[0,1]
        print(f"  Korelacja: {corr:.4f}")

    # --- Wnioski ---
    fits4 = results[4][0] if 4 in results else {}
    k4_gs = fits4.get("g0-g*", (None,None))[0]
    k4_gs = k4_gs if k4_gs is not None else float('nan')
    print(f"""
WNIOSKI:

  1. Dla V(g) = g^3/3 - g^4/4 (TGP): k(g0-g*) = {k4_gs:.2f}
     (oczekiwane z prop:K-exponent: k ~ 4)

  2. Test: czy k skaluje sie liniowo z deg(V)?
     Jesli TAK -> argument heurystyczny prop:K-exponent jest POTWIERDZONY
     Jesli NIE -> wykladnik 4 jest specyficzny dla V(g) stopnia 4,
                  ale mechanizm jest inny niz prosty stopien

  3. Argument heurystyczny (Dod. K):
     E_sol ~ (g0-1)^2  (z V''(1) = -1)
     A_tail ~ E_sol^2 ~ (g0-1)^4
     Ale to nie tlmaczy dlaczego E_sol ~ (g0-1)^2 a nie (g0-1)^(n-2)

  4. Alternatywny argument wymiarowy:
     V(g) = g^3/3 - g^4/4, V(1) = 1/12
     V(g0) - V(1) ~ (g0-1)^2 * V''(1)/2 = -(g0-1)^2/2
     Energia solitonu E_sol ~ integral V ~ (g0-1)^2
     Niezaleznie od n_deg (bo V''(1) = 2-n, zawsze kwadratowy!)
     Wiec A ~ (g0-1)^k z k wynikajacym z innego mechanizmu.
""")

    return True


if __name__ == "__main__":
    run()
