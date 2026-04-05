#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p_frw_full_evolution.py
========================
R1 ROADMAP: Pelne rozwiazanie kosmologiczne psi(t) od BBN do dzis.

Rownanie operatorowe TGP w FRW (zmienna N = ln a):

  Oryginalne: psi_tt + 3H psi_t + 2 psi_t^2/psi + gamma(psi^2 - psi^3) = -(3/2) rho_m / Phi_0

  W zmiennej N = ln(a):
    psi'' + (3-eps_H) psi' + 2(psi')^2/psi = [-gamma(psi^2 - psi^3) - kappa*Omega_m/a^3] / H^2

  gdzie:
    psi = Phi / Phi_0  (znormalizowane pole)
    V_self(psi) = gamma*(psi^2 - psi^3) = gamma*psi^2*(1-psi)
      -> V_self jest czlonem potencjalowym BEZPOSREDNIO z operatora TGP
      -> V_self(1) = 0 (warunek prozni, beta=gamma)
      -> V_self(7/6) = gamma*(49/36)*(-1/6) = -49*gamma/216 < 0
    kappa = 3/(2*Phi_0)  (stala sprzezenia z materia; z q = 4*pi*G/c^2)
    H^2 = [Omega_r/a^4 + Omega_m/a^3 + Omega_Lambda] / psi  (zmodyfikowane Friedmann)
    eps_H = -d(ln H)/dN  (parametr slow-roll Hubble'a)

  FIZYKA:
    - V_self na LEWEJ stronie = czlon potencjalowy w operatorze D[Phi]
    - Na RHS mamy: -V_self/H^2 (przeniesiony na prawa strone)
    - Dla psi > 1: -V_self = gamma*psi^2*(psi-1) > 0 -> odpycha OD psi=1
    - Materia: -kappa*rho_m/H^2 < 0 -> przyciaga KU psi=1
    - Rownowaga daje psi(z=0) ~ 1 + delta_psi

  KLUCZOWE WYJSCIA:
    1. delta_psi = |psi(z=0) - 1|  vs  granica LLR: 0.02
    2. |G_dot/G|/H_0 = |psi'(0)|/psi(0)  vs  LLR: 0.02
    3. Krzywe psi(z), G_eff(z)/G_0, c_eff(z)/c_0

Sesja: ROADMAP v3, R1 (2026-03-29)
"""

import sys
import io

# Windows-safe UTF-8
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import argparse
import numpy as np

try:
    from scipy.integrate import solve_ivp
except ImportError:
    print("FATAL: scipy wymagany. Zainstaluj: pip install scipy")
    sys.exit(1)

# ============================================================
#  Parametry
# ============================================================

# --- Kosmologiczne (Planck 2018) ---
OMEGA_M = 0.3153
OMEGA_R = 9.15e-5
OMEGA_L = 0.6847

# --- TGP ---
GAMMA_TGP = 0.5        # gamma = beta (bezwymiarowe, w H_0^2)
PSI_INI   = 7.0 / 6.0  # warunek poczatkowy (atraktor GL)
PHI0      = 24.66       # Phi_0 (skala generowanej przestrzeni)
KAPPA_DEFAULT = 1.5 / PHI0  # = 3/(2*Phi_0) ~ 0.061

# --- Ograniczenia ---
GDOT_G_BOUND = 0.02  # |Gdot/G| < 0.02*H_0 (LLR)

# --- Zakres ---
Z_BBN = 1e9
N_INI = -np.log(1 + Z_BBN)
N_FIN = 0.0


# ============================================================
#  Funkcje fizyczne
# ============================================================

def V_self(psi, gamma=GAMMA_TGP):
    """
    Czlon potencjalowy z operatora TGP:
        V_self(psi) = gamma * (psi^2 - psi^3) = gamma * psi^2 * (1 - psi)

    Wlasnosci:
      V_self(0) = 0
      V_self(1) = 0  (warunek prozni beta=gamma)
      V_self(2/3) = gamma * (4/9 - 8/27) = 4*gamma/27  (maximum)
      V_self -> -inf  dla psi -> +inf

    W rownaniu TGP: V_self jest po LEWEJ stronie operatora:
      D[Phi] = box(psi) + ... + V_self(psi) = -source
    Przeniesiony na prawa:
      box(psi) + ... = -V_self(psi) - source
    Wiec -V_self(psi) = gamma*psi^2*(psi-1) jest SILA NAPEDOWA.
    """
    return gamma * psi**2 * (1.0 - psi)


def H_squared(a, psi):
    """
    Zmodyfikowane Friedmann TGP:
        H^2 = [Omega_r/a^4 + Omega_m/a^3 + Omega_Lambda] / psi

    Pomijamy czlon kinetyczny (maly dla wolno zmiennego psi).
    Omega_Lambda absorbuje energe prozni przy psi=1.
    """
    rho = OMEGA_R / a**4 + OMEGA_M / a**3 + OMEGA_L
    psi_safe = max(psi, 1e-15)
    return max(rho / psi_safe, 1e-30)


def eps_H(a, H2, psi):
    """
    eps_H = -dH/dN / H = (3/2)*(rho_m + (4/3)*rho_r) / (H^2 * psi)

    W erze radiacji: eps_H -> 2
    W erze materii: eps_H -> 3/2
    W erze DE: eps_H -> 0
    """
    rho_m = OMEGA_M / a**3
    rho_r = OMEGA_R / a**4
    denom = H2 * max(psi, 1e-15)
    if denom <= 0:
        return 0.0
    return min((3.0/2.0) * (rho_m + (4.0/3.0) * rho_r) / denom, 4.0)


# ============================================================
#  Uklad ODE
# ============================================================

def ode_rhs(N, y, gamma=GAMMA_TGP, kappa=KAPPA_DEFAULT):
    """
    Uklad ODE dla psi(N), dpsi/dN(N), N = ln(a).

    y = [psi, u]  gdzie u = dpsi/dN

    Rownanie TGP-FRW:
        psi'' + (3 - eps_H)*psi' + 2*(psi')^2/psi
          = [-V_self(psi) - kappa*Omega_m/a^3] / H^2

    KLUCZOWE:
      - V_self = gamma*psi^2*(1-psi) jest BEZPOSREDNI czlon operatora
      - Na RHS mamy -V_self (ze znakiem minus po przeniesieniu)
      - kappa = 3/(2*Phi_0) to stala sprzezenia z materia
    """
    psi, u = y

    if psi <= 1e-10:
        psi = 1e-10

    a = np.exp(N)

    # Friedmann
    H2 = H_squared(a, psi)

    # eps_H
    eps = eps_H(a, H2, psi)

    # Tlumienie
    damping = 3.0 - eps

    # Sila napedowa:
    # -V_self = -gamma*psi^2*(1-psi) = gamma*psi^2*(psi-1)
    # Materia: -kappa*Omega_m/a^3
    force_vacuum = -V_self(psi, gamma)   # = gamma*psi^2*(psi-1)
    force_matter = -kappa * OMEGA_M / a**3

    total_force = (force_vacuum + force_matter) / H2

    # Czlon nieliniowy alpha=2
    nonlinear = 2.0 * u**2 / psi

    # du/dN = -damping*u - nonlinear + total_force
    du_dN = -damping * u - nonlinear + total_force

    return [u, du_dN]


# ============================================================
#  Solver
# ============================================================

def solve_psi_evolution(psi_ini=PSI_INI, dpsi_ini=0.0,
                        z_start=Z_BBN, gamma=GAMMA_TGP,
                        kappa=KAPPA_DEFAULT,
                        n_eval=5000, rtol=1e-10, atol=1e-12):
    """Rozwiaz ewolucje psi(N) od z_start do z=0."""

    N_start = -np.log(1 + z_start)
    N_end = 0.0
    N_eval = np.linspace(N_start, N_end, n_eval)

    y0 = [psi_ini, dpsi_ini]

    def rhs(N, y):
        return ode_rhs(N, y, gamma, kappa)

    sol = solve_ivp(rhs, [N_start, N_end], y0,
                    method='Radau',
                    t_eval=N_eval,
                    rtol=rtol, atol=atol,
                    max_step=0.05)

    if not sol.success:
        print(f"  OSTRZEZENIE: solver problem: {sol.message}")
        sol = solve_ivp(rhs, [N_start, N_end], y0,
                        method='Radau',
                        t_eval=N_eval,
                        rtol=1e-8, atol=1e-10,
                        max_step=0.1)

    N = sol.t
    psi = sol.y[0]
    dpsi_dN = sol.y[1]
    a = np.exp(N)
    z = 1.0 / a - 1.0

    H2 = np.array([H_squared(a[i], psi[i]) for i in range(len(N))])
    G_ratio = 1.0 / np.maximum(psi, 1e-15)
    c_ratio = 1.0 / np.sqrt(np.maximum(psi, 1e-15))

    return {
        'N': N, 'z': z, 'a': a,
        'psi': psi, 'dpsi_dN': dpsi_dN,
        'H2': H2, 'G_ratio': G_ratio, 'c_ratio': c_ratio,
        'success': sol.success, 'message': sol.message,
        'gamma': gamma, 'psi_ini': psi_ini, 'kappa': kappa
    }


# ============================================================
#  Diagnostyka
# ============================================================

RESULTS = []
VERBOSE = True

def test(section, name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    RESULTS.append((section, name, passed, detail))
    if VERBOSE:
        print(f"  [{status}] {name}: {detail}")


def run_diagnostics(res):
    """Pelna diagnostyka wynikow."""

    psi = res['psi']
    dpsi = res['dpsi_dN']
    z = res['z']
    a = res['a']
    G_ratio = res['G_ratio']
    H2 = res['H2']

    psi_0 = psi[-1]
    dpsi_0 = dpsi[-1]
    delta_psi = abs(psi_0 - 1.0)
    Gdot_G = abs(dpsi_0) / max(psi_0, 1e-15)

    print(f"\n{'='*70}")
    print(f"  DIAGNOSTYKA EWOLUCJI KOSMOLOGICZNEJ TGP")
    print(f"{'='*70}")

    # ---- S1: Wartosc dzisiejsza ----
    print(f"\n--- S1: Wartosci dzisiejsze psi(z=0) ---")
    print(f"  psi(z=0) = {psi_0:.8f}")
    print(f"  delta_psi = |psi(0)-1| = {delta_psi:.6f}")
    print(f"  dpsi/dN(z=0) = {dpsi_0:.8e}")
    print(f"  |G_dot/G|/H_0 = {Gdot_G:.6f}")

    test("S1", "T1: delta_psi < 0.05 (liberalna)",
         delta_psi < 0.05, f"delta_psi = {delta_psi:.4f}")
    test("S1", "T2: delta_psi < 0.02 (LLR)",
         delta_psi < 0.02, f"delta_psi = {delta_psi:.4f}, LLR = 0.02")
    test("S1", "T3: |Gdot/G|/H_0 < 0.02 (LLR)",
         Gdot_G < 0.02, f"|Gdot/G|/H_0 = {Gdot_G:.6f}")

    # ---- S2: Zamrozenie w erze radiacji ----
    print(f"\n--- S2: Zamrozenie psi w erze radiacji ---")
    mask_rad = a < 1e-4
    if np.any(mask_rad):
        psi_rad = psi[mask_rad]
        psi_rad_std = np.std(psi_rad)
        psi_rad_mean = np.mean(psi_rad)
        print(f"  <psi>_rad = {psi_rad_mean:.6f} (oczekiwane ~{PSI_INI:.6f})")
        print(f"  std(psi)_rad = {psi_rad_std:.2e}")
        test("S2", "T4: psi zamrozone w erze radiacji",
             psi_rad_std < 0.01, f"std = {psi_rad_std:.2e}")
        test("S2", "T5: <psi>_rad ~ psi_ini",
             abs(psi_rad_mean - PSI_INI) < 0.05,
             f"<psi> = {psi_rad_mean:.4f}")

    # ---- S3: Ewolucja materia ----
    print(f"\n--- S3: Ewolucja w erze materii ---")
    mask_mat = (z > 10) & (z < 3000)
    if np.any(mask_mat):
        psi_mat = psi[mask_mat]
        drift = abs(psi_mat[-1] - psi_mat[0])
        print(f"  psi(z=3000) = {psi_mat[0]:.6f}")
        print(f"  psi(z=10) = {psi_mat[-1]:.6f}")
        print(f"  drift = {drift:.6f}")
        test("S3", "T6: drift materii < 0.2",
             drift < 0.2, f"drift = {drift:.4f}")

    # ---- S4: Ewolucja pozna ----
    print(f"\n--- S4: Ewolucja pozna (z < 5) ---")
    for z_key in [5.0, 2.0, 1.0, 0.5, 0.1, 0.0]:
        idx = np.argmin(np.abs(z - z_key))
        print(f"  psi(z={z_key:.1f}) = {psi[idx]:.6f}")

    # ---- S5: Stale dynamiczne ----
    print(f"\n--- S5: Stale dynamiczne ---")
    G_today = G_ratio[-1]
    print(f"  G_eff(z=0)/G_0 = {G_today:.6f}")
    test("S5", "T7: G_eff(0) bliskie G_0 (10%)",
         abs(G_today - 1.0) < 0.10, f"G/G_0 = {G_today:.4f}")

    mask_bbn = (z > 1e8) & (z < 1e10)
    if np.any(mask_bbn):
        G_bbn = np.mean(G_ratio[mask_bbn])
        print(f"  G_eff(BBN)/G_0 = {G_bbn:.6f} (ocz. {1/PSI_INI:.6f})")
        test("S5", "T8: G_eff(BBN) = G_0/psi_ini",
             abs(G_bbn - 1.0/PSI_INI) / (1.0/PSI_INI) < 0.01,
             f"G_bbn = {G_bbn:.4f}")

    # ---- S6: Wlasnosci globalne ----
    print(f"\n--- S6: Wlasnosci globalne ---")
    test("S6", "T9: psi(0) < psi_ini",
         psi_0 < PSI_INI, f"psi(0)={psi_0:.4f}")
    test("S6", "T10: psi > 0 zawsze",
         np.all(psi > 0), f"min(psi)={np.min(psi):.6f}")
    test("S6", "T11: H^2 > 0 zawsze",
         np.all(H2 > 0), f"min(H^2)={np.min(H2):.2e}")

    return {
        'psi_0': psi_0, 'delta_psi': delta_psi,
        'Gdot_G': Gdot_G, 'dpsi_0': dpsi_0
    }


def run_force_analysis(res):
    """Analiza sil napedowych w roznych epokach."""

    print(f"\n--- Analiza sil napedowych ---")
    psi = res['psi']
    dpsi = res['dpsi_dN']
    a = res['a']
    z = res['z']
    H2 = res['H2']
    gamma = res['gamma']
    kappa = res['kappa']

    epochs = [
        ("Radiacja (z~10^6)", 1e5, 1e7),
        ("Rownowaga (z~3400)", 2000, 5000),
        ("Materia (z~100)", 50, 200),
        ("Materia (z~10)", 5, 20),
        ("Przejscie (z~1)", 0.5, 2.0),
        ("Dzis (z~0)", 0.0, 0.1),
    ]

    print(f"  {'Epoka':<25} {'<psi>':>8} {'F_vac/H2':>10} {'F_mat/H2':>10} {'F_net/H2':>10} {'damping':>8}")
    print(f"  {'-'*25} {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    for name, z_lo, z_hi in epochs:
        mask = (z >= z_lo) & (z <= z_hi)
        if not np.any(mask):
            continue

        psi_ep = psi[mask]
        a_ep = a[mask]
        H2_ep = H2[mask]

        # Srednie sily
        f_vac = np.mean([-V_self(psi_ep[i], gamma) / H2_ep[i] for i in range(len(psi_ep))])
        f_mat = np.mean([-kappa * OMEGA_M / (a_ep[i]**3 * H2_ep[i]) for i in range(len(a_ep))])
        f_net = f_vac + f_mat

        eps_arr = np.array([eps_H(a_ep[i], H2_ep[i], psi_ep[i]) for i in range(len(a_ep))])
        damp = np.mean(3.0 - eps_arr)

        print(f"  {name:<25} {np.mean(psi_ep):8.4f} {f_vac:10.4e} {f_mat:10.4e} {f_net:10.4e} {damp:8.3f}")


def scan_kappa(kappas, psi_ini=PSI_INI, gamma=GAMMA_TGP):
    """Skan delta_psi vs kappa (stala sprzezenia z materia)."""

    print(f"\n{'='*70}")
    print(f"  SKAN KAPPA: delta_psi vs kappa (stala sprzezenia)")
    print(f"{'='*70}")
    print(f"  Phi_0 = {PHI0:.2f}, kappa_fiz = 3/(2*Phi_0) = {KAPPA_DEFAULT:.4f}")
    print(f"\n  {'kappa':>10} | {'Phi_0_eff':>10} | {'psi(0)':>10} | {'delta_psi':>10} | {'|Gdot/G|':>10}")
    print(f"  {'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")

    results = []
    for k in kappas:
        try:
            phi0_eff = 1.5 / k if k > 0 else float('inf')
            res = solve_psi_evolution(psi_ini=psi_ini, gamma=gamma,
                                      kappa=k, n_eval=3000,
                                      rtol=1e-9, atol=1e-11)
            psi_0 = res['psi'][-1]
            dpsi_0 = res['dpsi_dN'][-1]
            delta = abs(psi_0 - 1.0)
            gdot = abs(dpsi_0) / max(psi_0, 1e-15)
            print(f"  {k:10.5f} | {phi0_eff:10.2f} | {psi_0:10.6f} | {delta:10.6f} | {gdot:10.6f}")
            results.append((k, phi0_eff, psi_0, delta, gdot))
        except Exception as e:
            print(f"  {k:10.5f} | FAIL: {e}")
    return results


def scan_gamma(gammas, psi_ini=PSI_INI, kappa=KAPPA_DEFAULT):
    """Skan delta_psi vs gamma."""

    print(f"\n{'='*70}")
    print(f"  SKAN GAMMA: delta_psi vs gamma")
    print(f"{'='*70}")
    print(f"\n  {'gamma':>10} | {'psi(0)':>10} | {'delta_psi':>10} | {'|Gdot/G|':>10}")
    print(f"  {'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")

    for g in gammas:
        try:
            res = solve_psi_evolution(psi_ini=psi_ini, gamma=g,
                                      kappa=kappa, n_eval=3000,
                                      rtol=1e-9, atol=1e-11)
            psi_0 = res['psi'][-1]
            dpsi_0 = res['dpsi_dN'][-1]
            delta = abs(psi_0 - 1.0)
            gdot = abs(dpsi_0) / max(psi_0, 1e-15)
            print(f"  {g:10.4f} | {psi_0:10.6f} | {delta:10.6f} | {gdot:10.6f}")
        except Exception as e:
            print(f"  {g:10.4f} | FAIL: {e}")


# ============================================================
#  Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="TGP: Pelna ewolucja FRW psi(t)")
    parser.add_argument('--gamma', type=float, default=GAMMA_TGP)
    parser.add_argument('--psi-ini', type=float, default=PSI_INI)
    parser.add_argument('--kappa', type=float, default=KAPPA_DEFAULT)
    parser.add_argument('--z-start', type=float, default=Z_BBN)
    parser.add_argument('--scan-kappa', action='store_true',
                        help="Skan delta_psi vs kappa")
    parser.add_argument('--scan-gamma', action='store_true',
                        help="Skan delta_psi vs gamma")
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    global VERBOSE
    VERBOSE = True

    print("="*70)
    print("  TGP: PELNA EWOLUCJA FRW psi(z)")
    print("  Roadmap v3, R1")
    print("="*70)
    print(f"  gamma = beta = {args.gamma}")
    print(f"  psi_ini = {args.psi_ini:.6f}")
    print(f"  kappa = {args.kappa:.5f}  (Phi_0_eff = {1.5/args.kappa:.2f})")
    print(f"  z_start = {args.z_start:.2e}")
    print(f"  Omega_m = {OMEGA_M}, Omega_r = {OMEGA_R:.2e}, Omega_L = {OMEGA_L}")

    # ---- Glowne rozwiazanie ----
    print(f"\n--- Calkowanie glowne ---")
    res = solve_psi_evolution(psi_ini=args.psi_ini, gamma=args.gamma,
                              kappa=args.kappa, z_start=args.z_start,
                              n_eval=5000, rtol=1e-10, atol=1e-12)
    print(f"  Solver: {'SUKCES' if res['success'] else 'PROBLEM'} ({res['message']})")
    print(f"  Punktow: {len(res['N'])}")

    # ---- Diagnostyka ----
    diag = run_diagnostics(res)

    # ---- Analiza sil ----
    run_force_analysis(res)

    # ---- Skan kappa ----
    if args.scan_kappa:
        kappas = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15, 0.20]
        scan_kappa(kappas, psi_ini=args.psi_ini, gamma=args.gamma)

    # ---- Skan gamma ----
    if args.scan_gamma:
        gammas = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
        scan_gamma(gammas, psi_ini=args.psi_ini, kappa=args.kappa)

    # ---- Wykresy ----
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import os

            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            z = res['z']
            psi = res['psi']
            mask = z > 0

            ax = axes[0, 0]
            ax.semilogx(z[mask] + 1, psi[mask], 'b-', lw=1.5)
            ax.axhline(1.0, color='k', ls='--', alpha=0.5, label='psi=1')
            ax.axhline(PSI_INI, color='r', ls=':', alpha=0.5, label=f'psi_ini={PSI_INI:.4f}')
            ax.set_xlabel('1+z'); ax.set_ylabel('psi(z)')
            ax.set_title('(a) Ewolucja pola psi(z)')
            ax.legend(); ax.grid(True, alpha=0.3)

            ax = axes[0, 1]
            ax.semilogx(z[mask] + 1, res['G_ratio'][mask], 'r-', lw=1.5)
            ax.axhline(1.0, color='k', ls='--', alpha=0.5)
            ax.set_xlabel('1+z'); ax.set_ylabel('G_eff/G_0')
            ax.set_title('(b) G_eff(z)/G_0'); ax.grid(True, alpha=0.3)

            ax = axes[1, 0]
            ax.semilogx(z[mask] + 1, res['c_ratio'][mask], 'g-', lw=1.5)
            ax.axhline(1.0, color='k', ls='--', alpha=0.5)
            ax.set_xlabel('1+z'); ax.set_ylabel('c_eff/c_0')
            ax.set_title('(c) c_eff(z)/c_0'); ax.grid(True, alpha=0.3)

            ax = axes[1, 1]
            ax.semilogx(z[mask] + 1, res['dpsi_dN'][mask], 'm-', lw=1.5)
            ax.axhline(0.0, color='k', ls='--', alpha=0.5)
            ax.set_xlabel('1+z'); ax.set_ylabel("dpsi/dN")
            ax.set_title("(d) psi'(N)"); ax.grid(True, alpha=0.3)

            plt.suptitle(f'TGP FRW: psi_ini={PSI_INI:.4f}, gamma={args.gamma}, kappa={args.kappa:.4f}',
                        fontsize=13)
            plt.tight_layout()
            os.makedirs("plots", exist_ok=True)
            outpath = "plots/frw_full_evolution.png"
            plt.savefig(outpath, dpi=150, bbox_inches='tight')
            print(f"\n  Wykres: {outpath}")
        except ImportError:
            print("\n  matplotlib niedostepny")

    # ---- Podsumowanie ----
    n_pass = sum(1 for r in RESULTS if r[2])
    n_total = len(RESULTS)

    print(f"\n{'='*70}")
    print(f"  PODSUMOWANIE: {n_pass}/{n_total} PASS")
    print(f"{'='*70}")
    print(f"  delta_psi = {diag['delta_psi']:.6f}")
    print(f"  |Gdot/G|/H_0 = {diag['Gdot_G']:.6f}")
    print(f"  Granica LLR: 0.02")

    if diag['delta_psi'] < 0.02:
        print(f"\n  *** WYNIK: delta_psi < 0.02 -> ZGODNE z LLR ***")
    elif diag['delta_psi'] < 0.05:
        print(f"\n  *** WYNIK: 0.02 < delta_psi < 0.05 -> NAPIECIE (czynnik {diag['delta_psi']/0.02:.1f}x) ***")
    else:
        print(f"\n  *** WYNIK: delta_psi > 0.05 -> NIEZGODNE z LLR ***")

    return 0 if n_pass == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
