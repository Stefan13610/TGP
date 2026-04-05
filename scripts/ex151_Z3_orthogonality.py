"""
ex151_Z3_orthogonality.py
==========================
Test hipotezy: Z_3 wynika z ortogonalnosci modow solitonowych.

IDEA:
1. Trzy mody solitonowe (e, mu, tau) to rozwiazania tego samego ODE
   z roznymi warunkami poczatkowymi g_0^(n).
2. Ich ogony u_n(r) = g_n(r) - 1 ~ A_n sin(r + delta_n)/r sa
   falami sinusoidalnymi z roznymi fazami delta_n.
3. Ortogonalnosc: <u_n | u_m> = integral u_n(r) * u_m(r) * r^2 dr = 0
   Jest to ortogonalnosc w sensie Hilbertowskim.
4. Dla fal sin(r+d_n)/r: <u_n|u_m> ~ A_n*A_m*cos(d_n-d_m)
   (dla r >> 1, calka asymptotyczna)
5. Jesli SUM_pairs cos(d_n-d_m) = -3/2, to K = 0 => Z_3!

TEST NUMERYCZNY:
- Oblicz profil solitonowy dla g_0^e, g_0^mu, g_0^tau
- Wyznacz fazy ogonowe delta_n
- Sprawdz warunek ortogonalnosci
- Sprawdz K = |sum e^{i*delta_n}|^2

KLUCZOWE ROZROZNIENIE:
ex125 pokazal ze surowe fazy NIE tworza Z_3. Ale moze fazy
w INNEJ bazie (np. znormalizowanej) tworza Z_3?
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp

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


def extract_tail_phase(r, g, rL=20, rR=60):
    """Wyznacz amplitude i faze ogona: u(r) = A*sin(r+delta)/r"""
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 50:
        return 0.0, 0.0

    rf = r[mask]
    df = (g[mask] - 1.0) * rf  # r * u(r) = B*cos(r) + C*sin(r)

    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc[0], bc[1]
    A = np.sqrt(B**2 + C**2)
    delta = np.arctan2(B, C)  # u ~ (B*cos(r)+C*sin(r))/r = A*sin(r+delta)/r
    return A, delta


def compute_overlap(r, g1, g2, rL=5, rR=80):
    """Calka nakladki: <u1|u2> = integral u1(r)*u2(r)*r^2 dr"""
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 50:
        return 0.0

    u1 = g1[mask] - 1.0
    u2 = g2[mask] - 1.0
    r_m = r[mask]
    integrand = u1 * u2 * r_m**2
    return np.trapezoid(integrand, r_m)


def run():
    print("=" * 72)
    print("ex151: Test hipotezy Z_3 z ortogonalnosci modow solitonowych")
    print("=" * 72)

    # --- 1. Oblicz profile ---
    # ODE kanoniczne (fermion_mass_spectrum.py): alpha=2
    # g_0^e ~ 1.24, g_0^mu ~ phi*g_0^e ~ 2.00, g_0^tau ~ 2.34
    g0_e = 1.24
    g0_mu = PHI * g0_e  # ~ 2.006
    g0_tau = 2.34  # z predykcji C4

    print(f"\n  g0_e  = {g0_e:.4f}")
    print(f"  g0_mu = {g0_mu:.4f} (phi * g0_e)")
    print(f"  g0_tau = {g0_tau:.4f}")
    print(f"  g*     = {G_GHOST:.4f}")

    print("\n  Obliczanie solitonow...", flush=True)
    r_e, g_e = solve_soliton(g0_e)
    r_mu, g_mu = solve_soliton(g0_mu)
    r_tau, g_tau = solve_soliton(g0_tau)

    # Interpoluj na wspolna siatke
    r_common = np.linspace(0.01, 80.0, 10000)
    g_e_c = np.interp(r_common, r_e, g_e)
    g_mu_c = np.interp(r_common, r_mu, g_mu)
    g_tau_c = np.interp(r_common, r_tau, g_tau)

    # --- 2. Wyznacz fazy ogonowe ---
    print("\n--- 1. Fazy ogonowe (delta_n) ---")
    windows = [(20, 40), (30, 50), (40, 60), (20, 60)]
    for label, g0, r_arr, g_arr in [
        ("e", g0_e, r_e, g_e),
        ("mu", g0_mu, r_mu, g_mu),
        ("tau", g0_tau, r_tau, g_tau),
    ]:
        print(f"\n  {label} (g0={g0:.4f}):")
        for rL, rR in windows:
            A, delta = extract_tail_phase(r_arr, g_arr, rL, rR)
            print(f"    [{rL},{rR}]: A = {A:.6f}, delta = {delta:.4f} rad "
                  f"= {np.degrees(delta):.2f} deg")

    # Fazy z najszersiego okna
    A_e, d_e = extract_tail_phase(r_e, g_e, 20, 60)
    A_mu, d_mu = extract_tail_phase(r_mu, g_mu, 20, 60)
    A_tau, d_tau = extract_tail_phase(r_tau, g_tau, 20, 60)

    print(f"\n  Fazy (okno [20,60]):")
    print(f"    delta_e   = {np.degrees(d_e):.2f} deg")
    print(f"    delta_mu  = {np.degrees(d_mu):.2f} deg")
    print(f"    delta_tau = {np.degrees(d_tau):.2f} deg")

    # --- 3. Koherencja fazowa ---
    print("\n--- 2. Koherencja fazowa K ---")
    phases = np.array([d_e, d_mu, d_tau])
    K = abs(np.sum(np.exp(1j * phases)))**2
    print(f"  K = |sum e^(i*delta_n)|^2 = {K:.4f}")
    print(f"  (K = 0 => Z_3; K = 9 => koherentne; K = 3 => losowe)")

    cos_sum = 0
    pairs = [(0, 1, "e-mu"), (1, 2, "mu-tau"), (0, 2, "e-tau")]
    for i, j, label in pairs:
        c = np.cos(phases[i] - phases[j])
        cos_sum += c
        print(f"    cos(delta_{label}) = {c:.4f}")
    print(f"  Sum cos = {cos_sum:.4f} (Z_3 => -1.5)")

    # --- 4. Ortogonalnosc ---
    print("\n--- 3. Ortogonalnosc modow <u_n|u_m> ---")
    labels = ["e", "mu", "tau"]
    gs = [g_e_c, g_mu_c, g_tau_c]

    # Calki normy i nakladki
    norms = []
    for i in range(3):
        nn = compute_overlap(r_common, gs[i], gs[i])
        norms.append(nn)
        print(f"  <u_{labels[i]}|u_{labels[i]}> = {nn:.4f}")

    print()
    overlaps = []
    for i in range(3):
        for j in range(i+1, 3):
            ov = compute_overlap(r_common, gs[i], gs[j])
            ov_norm = ov / np.sqrt(norms[i] * norms[j])
            overlaps.append(ov_norm)
            print(f"  <u_{labels[i]}|u_{labels[j]}> / sqrt(N_i*N_j) = {ov_norm:.6f}")
            print(f"    (rawowy = {ov:.4f})")

    # --- 5. Asymptotyczna nakladka ogonowa ---
    print("\n--- 4. Asymptotyczna nakladka ogonowa (tylko ogon) ---")
    for rL in [20, 30, 40]:
        print(f"  Okno [{rL}, 80]:")
        for i in range(3):
            for j in range(i+1, 3):
                ov = compute_overlap(r_common, gs[i], gs[j], rL=rL, rR=80)
                ni = compute_overlap(r_common, gs[i], gs[i], rL=rL, rR=80)
                nj = compute_overlap(r_common, gs[j], gs[j], rL=rL, rR=80)
                ov_n = ov / np.sqrt(ni * nj) if ni > 0 and nj > 0 else 0
                print(f"    <u_{labels[i]}|u_{labels[j]}>_norm = {ov_n:.6f}")

    # --- 6. Analiza roznic fazowych ---
    print("\n--- 5. Roznice fazowe i Z_3 ---")
    dd = [(d_mu - d_e) % (2 * np.pi),
          (d_tau - d_mu) % (2 * np.pi),
          (d_tau - d_e) % (2 * np.pi)]

    print(f"  delta_mu - delta_e  = {np.degrees(dd[0]):.2f} deg")
    print(f"  delta_tau - delta_mu = {np.degrees(dd[1]):.2f} deg")
    print(f"  delta_tau - delta_e = {np.degrees(dd[2]):.2f} deg")
    print(f"  (Z_3 wymaga: 120, 120, 240 lub permutacja)")

    # --- 7. Test: co daje K=0? ---
    print("\n--- 6. Szukanie konfiguracji K=0 ---")
    # Mamy d_e, d_mu. Jaki d_tau daje K=0?
    # e^{id_e} + e^{id_mu} + e^{id_tau} = 0
    # e^{id_tau} = -(e^{id_e} + e^{id_mu})
    z = -(np.exp(1j * d_e) + np.exp(1j * d_mu))
    d_tau_needed = np.angle(z)
    A_tau_needed = abs(z)  # powinno byc 1 jesli |e^id| = 1
    print(f"  Potrzebna faza tau: {np.degrees(d_tau_needed):.2f} deg")
    print(f"  Aktualna faza tau: {np.degrees(d_tau):.2f} deg")
    print(f"  Roznica: {np.degrees(d_tau - d_tau_needed):.2f} deg")
    print(f"  |-(e^id_e + e^id_mu)| = {A_tau_needed:.4f} (powinno byc 1)")

    # --- 8. Amplitudy i stosunek mas ---
    print("\n--- 7. Amplitudy i stosunek mas ---")
    print(f"  A_e   = {A_e:.6f}")
    print(f"  A_mu  = {A_mu:.6f}")
    print(f"  A_tau = {A_tau:.6f}")
    r21 = (A_mu / A_e)**4
    r31 = (A_tau / A_e)**4
    print(f"  r_21 = (A_mu/A_e)^4 = {r21:.2f} (PDG: 206.768)")
    print(f"  r_31 = (A_tau/A_e)^4 = {r31:.2f} (PDG: 3477.15)")

    # CV(A^2)
    A2 = np.array([A_e**2, A_mu**2, A_tau**2])
    cv = np.std(A2) / np.mean(A2)
    print(f"  CV(A^2) = {cv:.4f} (Z_3 w amplitudach => CV = 1.0)")

    # --- 9. Wnioski ---
    print(f"\n--- 8. Wnioski ---")
    print(f"""
  WYNIK KLUCZOWY:

  1. Fazy ogonowe delta_n NIE tworza Z_3 w prostej postaci:
     delta_e = {np.degrees(d_e):.1f}, delta_mu = {np.degrees(d_mu):.1f}, delta_tau = {np.degrees(d_tau):.1f} deg
     K = {K:.2f} (!=0)

  2. Ortogonalnosc modow: nakladki <u_n|u_m> sa NIEZEROWE
     (mody solitonowe NIE sa ortogonalne w sensie L2)

  3. INTERPRETACJA:
     Z_3 nie wynika z prostej ortogonalnosci modow solitonowych.
     Mody nie sa stanami wlasnymi wspolnego operatora samosprzezonego
     (bo kazdy ma INNE g_0, wiec rozwiazuje INNE rownanie).

  4. MOZLIWE PODEJSCIA DO R9:
     a) Z_3 jako symetria EFEKTYWNA w przestrzeni amplitudowej
        (Brannen parametryzacja z CV=1 — warunek GEOMETRYCZNY)
     b) Z_3 z minimalizacji entropii informacyjnej
     c) Z_3 z dynamiki substratu (sciezka jeszcze niezbadana)
     d) Z_3 jako postulat dodatkowy (status obecny)
""")

    return True


if __name__ == "__main__":
    run()
