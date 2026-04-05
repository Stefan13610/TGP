"""
ex152_tail_energy_A4.py
========================
R10: Nieperturbacyjny mechanizm M ~ A^4.

NOWA HIPOTEZA:
Masa fizyczna leptonu NIE jest M_core (energia jadra solitonu),
lecz E_tail — energia OGONA (r > r_c, gdzie r_c = pierwszy zer u=g-1).

Dlaczego E_tail ~ A^4?
======================
W ogonie: u(r) = g(r)-1 ~ A*sin(r+d)/r (tryb zerowy)
Lagrangian w ogonie (do rzedu 4):
  L = f(g)/2 * (g')^2 + V(g) - V(1)
  = 1/2*(u')^2 + nieliniowe

Energia drugiego rzedu: E^(2) = integral [(u')^2 - u^2] r^2 dr
Dla u = A*sin(r)/r: E^(2) = 0 (tryb zerowy!) - analitycznie

Energia czwartego rzedu: E^(4) = integral [V^(4)(1)/24 * u^4 + ...] r^2 dr
V(g) = g^3/3 - g^4/4, V(1) = 1/12
V''(1) = -1, V'''(1) = 2-6 = -4?, V''''(1) = -6
Wiec E^(4) ~ integral u^4 r^2 dr ~ A^4 * integral sin^4(r)/r^2 dr

ALE: perturbacyjnie E^(3) != 0 (ex148)! Wiec dlaczego nie A^3?

KLUCZOWY TEST:
1. Oblicz E_tail(g0) = integral_r_c^R [f(g)/2*(g')^2 + V(g)-V(1)] r^2 dr
2. Fituj E_tail vs A_tail — jaki wykladnik?
3. Oblicz E_tail_2, E_tail_4 osobno (dekompozycja)
4. Sprawdz E_tail/A^4 = const?
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp, simpson

ALPHA = 2.0
PHI = (1 + np.sqrt(5)) / 2
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005
R_MAX = 80.0


def f_tgp(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-12))


def V_pot(g):
    return g**3 / 3.0 - g**4 / 4.0

V1 = V_pot(1.0)  # = 1/12


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
    r_all, g_all, gp_all = [], [], []
    for _ in range(10):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_ghost,
                        dense_output=True, rtol=1e-9, atol=1e-11,
                        max_step=0.05)
        r_all.append(sol.t)
        g_all.append(sol.y[0])
        gp_all.append(sol.y[1])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]
            st = sol.sol(rh)
            y0 = [st[0], -st[1]]
            r0 = rh
        else:
            break
    r = np.concatenate(r_all)
    g = np.concatenate(g_all)
    gp = np.concatenate(gp_all)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def extract_A_tail(r, g):
    windows = [(16, 28), (20, 32), (24, 36), (28, 40)]
    A_vals = []
    for rL, rR in windows:
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 20:
            continue
        rf, df = r[mask], (g[mask] - 1.0) * r[mask]
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        A = np.sqrt(bc[0]**2 + bc[1]**2)
        if A > 1e-6:
            A_vals.append(A)
    return float(np.median(A_vals)) if A_vals else 0.0


def find_first_crossing(r, g):
    """Znajdz pierwszy punkt gdzie g(r) = 1 (u = 0) po r > 1."""
    u = g - 1.0
    idx0 = np.searchsorted(r, 1.0)
    if idx0 >= len(r) - 1:
        return None
    sc = np.where(np.diff(np.sign(u[idx0:])))[0]
    if len(sc) == 0:
        return None
    return r[idx0 + sc[0]]


def compute_energies(r, g, gp, r_cut):
    """
    Oblicz energie w roznych regionach:
    E = 4*pi * integral [f(g)/2*(g')^2 + V(g) - V(1)] r^2 dr

    E_core: [0, r_cut]
    E_tail: [r_cut, R_max]
    E_total: [0, R_max]
    """
    results = {}
    for label, rL, rR in [("core", 0, r_cut), ("tail", r_cut, R_MAX),
                            ("total", 0, R_MAX)]:
        mask = (r >= rL) & (r <= rR)
        if np.sum(mask) < 10:
            results[label] = 0.0
            continue
        r_m, g_m, gp_m = r[mask], g[mask], gp[mask]
        fg = f_tgp(g_m)
        integrand = (fg / 2.0 * gp_m**2 + V_pot(g_m) - V1) * r_m**2
        results[label] = 4.0 * np.pi * simpson(integrand, x=r_m)
    return results


def compute_tail_orders(r, g, gp, r_cut):
    """
    Dekompozycja energii ogona na rzedy w A:
    u = g - 1, u' = g'

    E^(2) = integral [(u')^2 - u^2] r^2 dr  (tryb zerowy)
    E^(3) = integral [-u^2*(u+1) + ...] r^2 dr  (nieliniowe)
    E^(4) = integral [u^4 * ...] r^2 dr

    Zamiast dekompozycji perturbacyjnej, obliczam:
    E_quad = 4pi integral [1/2*(u')^2 - 1/2*u^2] r^2 dr
    E_nonlin = E_tail - E_quad
    """
    mask = (r >= r_cut) & (r <= R_MAX)
    if np.sum(mask) < 10:
        return {}
    r_m = r[mask]
    u = g[mask] - 1.0
    up = gp[mask]

    # Energia kwadratowa (przyblizona, f(g)~1 w ogonie)
    integrand_quad = (0.5 * up**2 - 0.5 * u**2) * r_m**2
    E_quad = 4.0 * np.pi * simpson(integrand_quad, x=r_m)

    # Pelna energia ogona
    fg = f_tgp(g[mask])
    integrand_full = (fg / 2.0 * up**2 + V_pot(g[mask]) - V1) * r_m**2
    E_full = 4.0 * np.pi * simpson(integrand_full, x=r_m)

    # Nieliniowa czesc
    E_nonlin = E_full - E_quad

    # Korekcja kinetyczna: (f(g)-1)/2 * (g')^2
    integrand_fkin = ((fg - 1.0) / 2.0 * up**2) * r_m**2
    E_fkin = 4.0 * np.pi * simpson(integrand_fkin, x=r_m)

    # Korekcja potencjalowa: V(1+u) - V(1) + u^2/2
    # V(1+u) = 1/12 + ... -u^2/2 + cubic + quartic
    V_exact = V_pot(g[mask]) - V1
    V_quad = -0.5 * u**2  # V''(1) = -1
    V_nonlin = V_exact - V_quad
    integrand_Vnonlin = V_nonlin * r_m**2
    E_Vnonlin = 4.0 * np.pi * simpson(integrand_Vnonlin, x=r_m)

    return {
        "E_quad": E_quad,
        "E_full": E_full,
        "E_nonlin": E_nonlin,
        "E_fkin": E_fkin,
        "E_Vnonlin": E_Vnonlin,
    }


def run():
    print("=" * 72)
    print("ex152: R10 — Masa jako energia ogona solitonu")
    print("=" * 72)

    # Skan g0
    g0_vals = np.concatenate([
        np.linspace(1.02, 1.10, 5),
        np.linspace(1.10, 1.30, 15),
        np.linspace(1.30, 1.60, 8),
    ])

    data = []
    print("\n  Obliczanie solitonow...", flush=True)
    for g0 in g0_vals:
        r, g, gp = solve_soliton(g0)
        A = extract_A_tail(r, g)
        r_c = find_first_crossing(r, g)

        if A > 1e-6 and r_c is not None:
            E = compute_energies(r, g, gp, r_c)
            T = compute_tail_orders(r, g, gp, r_c)
            data.append({
                "g0": g0, "A": A, "r_c": r_c,
                **{f"E_{k}": v for k, v in E.items()},
                **T
            })

    print(f"  Obliczono {len(data)} solitonow\n")

    if len(data) < 10:
        print("  ZA MALO DANYCH")
        return

    g0_arr = np.array([d["g0"] for d in data])
    A_arr = np.array([d["A"] for d in data])

    # --- 1. Tabela energii ---
    print("--- 1. Tabela energii ---")
    print(f"  {'g0':>7}  {'A_tail':>9}  {'r_c':>5}  {'E_core':>11}  {'E_tail':>11}  {'E_total':>11}  {'E_nonlin':>11}")
    print("  " + "-" * 75)
    for d in data:
        print(f"  {d['g0']:7.4f}  {d['A']:9.5f}  {d['r_c']:5.2f}  "
              f"{d['E_core']:11.4e}  {d['E_tail']:11.4e}  {d['E_total']:11.4e}  "
              f"{d.get('E_nonlin', 0):11.4e}")

    # --- 2. Fit E vs A^k ---
    print("\n--- 2. Fit energii vs A^k ---")
    quantities = ["E_core", "E_tail", "E_total", "E_nonlin", "E_quad"]

    for q in quantities:
        vals = np.array([d.get(q, 0) for d in data])
        valid = (A_arr > 1e-4) & (np.abs(vals) > 1e-14)
        if np.sum(valid) < 5:
            print(f"  {q:12s}: za malo danych")
            continue
        logA = np.log(A_arr[valid])
        logE = np.log(np.abs(vals[valid]))
        k, logc = np.polyfit(logA, logE, 1)
        res = np.std(logE - (logc + k * logA))
        sign = "+" if np.mean(vals[valid]) > 0 else "-"
        print(f"  {q:12s}: |E| ~ A^{k:.3f}  (sign={sign}, residual={res:.4f})")

    # --- 3. Stosunek mas z roznych energii ---
    print("\n--- 3. Stosunek mas r_21 z roznych energii ---")
    # phi-FP: g0_mu = phi * g0_e
    g0_star = 1.249
    best_e = min(data, key=lambda d: abs(d["g0"] - g0_star))
    best_mu = min(data, key=lambda d: abs(d["g0"] - PHI * g0_star))

    print(f"  g0_e  = {best_e['g0']:.4f}")
    print(f"  g0_mu = {best_mu['g0']:.4f}")
    print(f"  A_e   = {best_e['A']:.6f}")
    print(f"  A_mu  = {best_mu['A']:.6f}")
    print(f"  (A_mu/A_e)^4 = {(best_mu['A']/best_e['A'])**4:.2f}")
    print()

    for q in quantities:
        e_val = best_e.get(q, 0)
        mu_val = best_mu.get(q, 0)
        if abs(e_val) > 1e-14 and abs(mu_val) > 1e-14:
            ratio = abs(mu_val / e_val)
            print(f"  |{q}(mu)/{q}(e)| = {ratio:.2f}")

    # --- 4. Dekompozycja szczegolowa ---
    print("\n--- 4. Dekompozycja ogona: E_quad vs E_nonlin ---")
    E_quad_arr = np.array([d.get("E_quad", 0) for d in data])
    E_nonlin_arr = np.array([d.get("E_nonlin", 0) for d in data])
    E_fkin_arr = np.array([d.get("E_fkin", 0) for d in data])
    E_Vnonlin_arr = np.array([d.get("E_Vnonlin", 0) for d in data])

    print(f"  {'g0':>7}  {'E_quad':>11}  {'E_nonlin':>11}  {'E_fkin':>11}  {'E_Vnonlin':>11}  {'ratio':>8}")
    print("  " + "-" * 65)
    for i, d in enumerate(data):
        eq = E_quad_arr[i]
        en = E_nonlin_arr[i]
        ef = E_fkin_arr[i]
        ev = E_Vnonlin_arr[i]
        ratio = abs(en / eq) if abs(eq) > 1e-14 else float('nan')
        print(f"  {d['g0']:7.4f}  {eq:11.4e}  {en:11.4e}  {ef:11.4e}  {ev:11.4e}  {ratio:8.3f}")

    # --- 5. E_nonlin / A^4 = const? ---
    print("\n--- 5. E_nonlin / A^4 = const? ---")
    ratios = []
    for i, d in enumerate(data):
        if d["A"] > 1e-4 and abs(E_nonlin_arr[i]) > 1e-14:
            c = E_nonlin_arr[i] / d["A"]**4
            ratios.append(c)
            if i < 5 or i > len(data) - 3:
                print(f"  g0={d['g0']:.4f}: E_nonlin/A^4 = {c:.6f}")

    if ratios:
        print(f"  ...")
        print(f"  Mean = {np.mean(ratios):.6f}")
        print(f"  Std  = {np.std(ratios):.6f}")
        print(f"  CV   = {np.std(ratios)/abs(np.mean(ratios))*100:.1f}%")
        print(f"  (CV < 10% => c_M = const)")

    # --- 6. r_21 z E_nonlin ---
    print("\n--- 6. r_21 z E_nonlin ---")
    if abs(best_e.get("E_nonlin", 0)) > 1e-14:
        r21_En = abs(best_mu["E_nonlin"] / best_e["E_nonlin"])
        print(f"  r_21(E_nonlin) = |E_nonlin(mu)/E_nonlin(e)| = {r21_En:.2f}")
        print(f"  r_21(A^4)      = (A_mu/A_e)^4              = {(best_mu['A']/best_e['A'])**4:.2f}")
        print(f"  r_21(PDG)      = 206.768")

    # --- Wnioski ---
    print(f"\n--- 7. Wnioski ---")
    print("""
  PYTANIE CENTRALNE: Jaka jest FIZYCZNA wielkosc dajaca M ~ A^4?

  Opcje:
  (a) M = E_core (energia jadra) — ex149: M_core ~ A^{1.7} => NIE
  (b) M = E_tail (energia ogona) — testowane tutaj
  (c) M = E_nonlin (nieliniowa czesc energii ogona)
  (d) M = E_total (calkowita energia solitonu)
  (e) M ~ A^4 jest POSTULATEM — masa to "amplituda substratowa do 4"
      bez interpretacji energetycznej
""")

    return True


if __name__ == "__main__":
    run()
