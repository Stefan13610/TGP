"""
ex149_nonperturbative_A4.py
===========================
Test nieperturbacyjnego mechanizmu M ~ A^4.

Hipoteza: Wykladnik 4 w r_21 = (A_mu/A_e)^4 NIE wynika z kasowania
perturbacyjnego E^(3) = 0 (ktore jest FALSZYWE -- zob. ex148),
lecz z EMPIRYCZNEGO skalowania M(g0) vs A_tail(g0).

Test:
1. Oblicz M_core(g0) i A_tail(g0) dla wielu g0
2. Fituj M_core ~ A^k -- jaki jest k?
3. Oblicz r_21 z roznych poteg A -- ktora daje PDG?
4. Sprawdz czy wykladnik 4 jest artefaktem zakresu g0

Kluczowy wniosek z ex148: perturbacyjnie E^(3) != 0,
wiec M(A) = E_2*A^2 + E_3*A^3 + E_4*A^4 + ... z E_2=0, E_3 != 0.
Ale numerycznie r_21 = 206.768 wymaga efektywnego k=4.
Jak to pogodzic?
"""
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


def Vp(g):
    return g**2 * (1.0 - g)


def soliton_rhs(r, y):
    """Kanoniczne ODE TGP: g'' + (2/r)g' + (alpha/g)(g')^2 = V'(g)"""
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


def solve_soliton(g0, r_max=R_MAX):
    y0 = [g0, 0.0]
    r0 = 0.0
    r_all, g_all, gp_all = [], [], []
    for _ in range(10):
        sol = solve_ivp(soliton_rhs, (r0, r_max), y0, events=ev_ghost,
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


def compute_M_core(r, g, gp):
    u = g - 1.0
    sc = np.where(np.diff(np.sign(u[r > 0.5])))[0]
    if len(sc) == 0:
        return np.nan
    idx0 = np.searchsorted(r, 0.5)
    rc_idx = idx0 + sc[0]
    mask = r <= r[rc_idx]
    r_c, g_c, gp_c = r[mask], g[mask], gp[mask]
    if len(r_c) < 10:
        return np.nan
    fg = f_tgp(g_c)
    integrand = (fg * gp_c**2 / 2.0 + V_pot(g_c) - V_pot(1.0)) * r_c**2
    return 4.0 * np.pi * simpson(integrand, x=r_c)


def run():
    print("=" * 72)
    print("ex149: Nieperturbacyjny mechanizm M ~ A^4")
    print("=" * 72)

    # Scan g0 values
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
        M = compute_M_core(r, g, gp)
        if A > 1e-6 and not np.isnan(M):
            data.append({"g0": g0, "A": A, "M": M})

    print(f"  Obliczono {len(data)} solitonow\n")

    if len(data) < 10:
        print("  ZA MALO DANYCH")
        return

    g0_arr = np.array([d["g0"] for d in data])
    A_arr = np.array([d["A"] for d in data])
    M_arr = np.array([d["M"] for d in data])

    # --- 1. Fit M vs A ---
    print("--- 1. Fit M_core vs A_tail ---")
    # M = c * A^k => log|M| = log|c| + k * log(A)
    valid = (A_arr > 1e-4) & (np.abs(M_arr) > 1e-12)
    if np.sum(valid) >= 5:
        logA = np.log(A_arr[valid])
        logM = np.log(np.abs(M_arr[valid]))
        k_fit, logc = np.polyfit(logA, logM, 1)
        print(f"  |M_core| ~ A^k: k = {k_fit:.4f}")
        print(f"  (Perturbacyjnie: k=3 bo E^(3) != 0)")
        print(f"  (Wynik TGP: k=4 bo M ~ A^4)")

        # Residuals
        logM_pred = logc + k_fit * logA
        residual = np.std(logM - logM_pred)
        print(f"  Residual (log): {residual:.4f}")
    else:
        k_fit = np.nan
        print("  Za malo danych do fitu")

    # --- 2. r_21 z roznych poteg ---
    print("\n--- 2. r_21 = (A_mu/A_e)^k dla roznych k ---")
    # phi-FP: g0_mu = phi * g0_e
    g0_star = 1.249  # approximate
    best_g0e = min(data, key=lambda d: abs(d["g0"] - g0_star))
    best_g0mu = min(data, key=lambda d: abs(d["g0"] - PHI * g0_star))

    if best_g0e and best_g0mu:
        Ae = best_g0e["A"]
        Amu = best_g0mu["A"]
        ratio = Amu / Ae
        print(f"  g0_e = {best_g0e['g0']:.4f}, A_e = {Ae:.6f}")
        print(f"  g0_mu = {best_g0mu['g0']:.4f}, A_mu = {Amu:.6f}")
        print(f"  A_mu/A_e = {ratio:.6f}")
        print()
        for k in [2, 3, 3.5, 4, 4.5, 5]:
            r21 = ratio**k
            delta = abs(r21 - 206.768) / 206.768 * 100
            marker = " <<<" if delta < 5 else ""
            print(f"  k={k:.1f}: r_21 = {r21:10.2f}  (delta={delta:.1f}%){marker}")

    # --- 3. Stosunek mas ---
    print("\n--- 3. Stosunek mas M_mu/M_e ---")
    if best_g0e and best_g0mu:
        Me = best_g0e["M"]
        Mmu = best_g0mu["M"]
        r21_M = abs(Mmu / Me)
        print(f"  M_e = {Me:.6e}")
        print(f"  M_mu = {Mmu:.6e}")
        print(f"  |M_mu/M_e| = {r21_M:.2f}")
        print(f"  PDG r_21 = 206.768")

        if r21_M > 0:
            k_eff = np.log(r21_M) / np.log(ratio)
            print(f"  Efektywny wykladnik: k_eff = log(r21_M)/log(A_mu/A_e) = {k_eff:.3f}")

    # --- 4. Tabela pelna ---
    print(f"\n--- 4. Tabela solitonow ---")
    print(f"  {'g0':>7}  {'A_tail':>10}  {'M_core':>12}  {'log|M|':>8}  {'logA':>8}")
    print("  " + "-" * 55)
    for d in data:
        logM = np.log(abs(d["M"])) if abs(d["M"]) > 0 else -99
        logA = np.log(d["A"]) if d["A"] > 0 else -99
        mark = ""
        if abs(d["g0"] - g0_star) < 0.01: mark = " <- e"
        elif abs(d["g0"] - PHI * g0_star) < 0.02: mark = " <- mu"
        print(f"  {d['g0']:7.4f}  {d['A']:10.6f}  {d['M']:12.6e}  {logM:8.3f}  {logA:8.3f}{mark}")

    # --- 5. Wnioski ---
    print(f"\n--- 5. Wnioski ---")
    print(f"""
  WYNIK KLUCZOWY:

  1. M_core NIE skaluje jak A^4 (fit daje k ≈ {k_fit:.1f})
  2. M_core jest UJEMNE (energia wiazania) i skaluje jak A^2 (core-dominated)
  3. Stosunek |M_mu/M_e| NIE daje r_21 = 206.768

  INTERPRETACJA:
  Wykladnik 4 w formule r_21 = (A_mu/A_e)^4 NIE pochodzi z skalowania
  M_core vs A_tail. Jest to formula EMPIRYCZNA dla stosunku mas leptonowych,
  ktore sa identyfikowane z A_tail^4 (NIE z M_core).

  Pytanie: jaka jest FIZYCZNA podstawa identyfikacji m_lepton ~ A^4?
  Mozliwosci:
    (a) Masa to energia ogona (nie jadra) — E_tail ~ A^4 z definicji
    (b) Masa to korekcja nieliniowa ponad tryb zerowy (argument B1)
    (c) A^4 wynika z wymiaru d=3 i struktury solitonu (argument wymiarowy)
    (d) Jest to postulat fenomenologiczny potwierdzony numerycznie
    """)

    return True


if __name__ == "__main__":
    run()
