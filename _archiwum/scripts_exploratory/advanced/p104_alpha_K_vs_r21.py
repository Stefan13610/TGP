#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p104_alpha_K_vs_r21.py  --  TGP v1
====================================
Skan alpha_K w [3, 16]: obliczenie K*_1(alpha_K), K*_2(alpha_K),
oraz r21 = K*_2/K*_1(alpha_K).

Pytanie fundamentalne: czy istnieje alpha_K takie, ze K*_2/K*_1 = 206.77?
Jesli tak -> nowe wyprowadzenie alpha_K z warunku r21_PDG.
Jesli nie -> ratio jest stale ~9.75 niezaleznie od alpha_K.

Metodologia (z p103_v4 - POPRAWIONA):
  - Gal. A (K*_1): PSI=[1.0, 3.5], K in [1e-4, 0.50] log, mono=True
  - Gal. B (K*_2): PSI=[1.5, 6.0], K in [1e-4, 5.00] log, mono=True
  - 16 procesow rownoleglych (ProcessPoolExecutor)

Ref (alpha_K=8.5616): K*_1=0.010280, K*_2=0.100276, r21=9.7547

Data: 2026-03-26
"""
import sys
import io
import warnings
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')

# ── Stale globalne (dostepne we wszystkich procesach po re-imporcie) ──────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 500
RTOL   = 1e-8
ATOL   = 1e-10

R21_TARGET = 206.77
N_WORKERS  = 16

# Grid alpha_K: gestszy okolo wartosci referencyjnej 8.5616
ALPHA_VALUES = np.unique(np.round(np.concatenate([
    np.linspace(3.0,  6.0,  7),   # krok 0.50
    np.linspace(6.0,  9.5, 15),   # krok 0.25, pokrywa ref 8.5616
    np.array([8.5616]),            # punkt referencyjny
    np.linspace(9.5,  13.0,  8),  # krok 0.50
    np.linspace(13.0, 16.0,  4),  # krok 1.00
]), 4))

# Okna poszukiwan (szersze niz p103 - dla roznych alpha_K)
PSI1_LO, PSI1_HI, N_PSI1 = 1.0, 3.5,  70
PSI2_LO, PSI2_HI, N_PSI2 = 1.5, 6.0, 100
N_K1 = 25
N_K2 = 30


# ── Funkcje robocze (module-level -> picklable) ───────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6


def compute_one_alpha(alpha_K):
    """
    Oblicz K*_1 i K*_2 dla zadanego alpha_K.
    Wywolywana w osobnym procesie - nie uzywa zmiennych globalnych poza stalymi.
    """
    warnings.filterwarnings('ignore')

    # -- Lokalne definicje ODE (domkniecie nad alpha_K - OK w tym procesie) ----
    def dV(phi):
        return GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1.0)**5

    def ode_rhs(r, y):
        phi, dphi = y
        phi = max(phi, 1e-10)
        kfac = 1.0 + alpha_K / phi
        return [dphi,
                dV(phi)/kfac
                + alpha_K*dphi**2 / (2.0*phi**2*kfac)
                - 2.0/r*dphi]

    def solve_full(psi, K):
        dphi0 = -K / AGAM**2
        r_ev  = AGAM * (R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
        for method in ('DOP853', 'Radau'):
            try:
                sol = solve_ivp(ode_rhs, [AGAM, R_MAX], [psi, dphi0],
                                method=method, rtol=RTOL, atol=ATOL,
                                t_eval=r_ev)
                if sol.t[-1] < R_MAX * 0.99:
                    if method == 'DOP853':
                        continue
                    return None
                return sol
            except Exception:
                if method == 'DOP853':
                    continue
                return None
        return None

    def phi_end(psi, K):
        sol = solve_full(psi, K)
        return float(sol.y[0, -1]) if sol else np.nan

    def is_mono(sol, tol=0.01):
        return bool(np.all(np.diff(sol.y[0]) <= tol))

    def energy_sol(sol, K):
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))):
            return np.nan
        kfac = 1.0 + alpha_K / phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vpot = np.clip(np.array([V_mod(float(p)) for p in phi]), -1e12, 1e12)
        Ep   = 4*np.pi*np.trapezoid((Vpot - V1)*r**2, r)
        E    = Ek + Ep
        return E if np.isfinite(E) else np.nan

    def find_psi_mono(K, psi_lo, psi_hi, n_psi):
        """Pierwszy monotonny zero F(psi) = phi_end(psi,K) - 1."""
        psis = np.linspace(psi_lo, psi_hi, n_psi)
        Fv   = [phi_end(p, K) - 1.0 for p in psis]
        for i in range(len(Fv) - 1):
            fi, fi1 = Fv[i], Fv[i+1]
            if not (np.isfinite(fi) and np.isfinite(fi1) and fi*fi1 < 0):
                continue
            try:
                pz  = brentq(lambda p: phi_end(p, K) - 1.0,
                             psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                sol = solve_full(pz, K)
                if sol and is_mono(sol):
                    return pz, sol
            except Exception:
                pass
        return np.nan, None

    def g_U(K, psi_lo, psi_hi, n_psi):
        pz, sol = find_psi_mono(K, psi_lo, psi_hi, n_psi)
        if np.isnan(pz) or sol is None:
            return np.nan, np.nan
        E = energy_sol(sol, K)
        g = E / (4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan
        return g, pz

    def find_Kstar(K_lo, K_hi, n_K, psi_lo, psi_hi, n_psi):
        Ks = np.logspace(np.log10(K_lo), np.log10(K_hi), n_K)
        gv, psv = [], []
        for K in Ks:
            g, pz = g_U(K, psi_lo, psi_hi, n_psi)
            gv.append(g); psv.append(pz)

        for i in range(len(gv) - 1):
            gi, gi1 = gv[i], gv[i+1]
            if not (np.isfinite(gi) and np.isfinite(gi1) and gi*gi1 < 0):
                continue
            _last = [np.nan]

            def gfunc(K, _l=_last):
                g2, pz2 = g_U(K, psi_lo, psi_hi, n_psi)
                _l[0] = pz2
                return g2 if np.isfinite(g2) else np.nan

            try:
                ga = gfunc(Ks[i]); gb = gfunc(Ks[i+1])
                if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                    Kstar    = brentq(gfunc, Ks[i], Ks[i+1],
                                      xtol=1e-5, maxiter=40)
                    psi_star = _last[0]
                    if not np.isfinite(psi_star):
                        psi_star, _ = find_psi_mono(Kstar, psi_lo, psi_hi, n_psi)
                    sol_s = solve_full(psi_star, Kstar) if np.isfinite(psi_star) else None
                    E_s   = energy_sol(sol_s, Kstar) if sol_s else np.nan
                    g_s   = E_s/(4*np.pi*Kstar) - 1.0 if np.isfinite(E_s) else np.nan
                    return Kstar, psi_star, g_s
            except Exception:
                Kstar    = Ks[i] - gv[i]*(Ks[i+1]-Ks[i]) / (gv[i+1]-gv[i])
                psi_star, sol_s = find_psi_mono(Kstar, psi_lo, psi_hi, n_psi)
                E_s  = energy_sol(sol_s, Kstar) if sol_s else np.nan
                g_s  = E_s/(4*np.pi*Kstar) - 1.0 if np.isfinite(E_s) else np.nan
                return Kstar, psi_star, g_s

        return np.nan, np.nan, np.nan

    # -- Glowne obliczenia -----------------------------------------------------
    tn = time.time()

    K1, psi1, g1 = find_Kstar(1e-4, 0.50, N_K1, PSI1_LO, PSI1_HI, N_PSI1)
    K2, psi2, g2 = find_Kstar(1e-4, 5.00, N_K2, PSI2_LO, PSI2_HI, N_PSI2)

    r21 = K2/K1 if (np.isfinite(K1) and np.isfinite(K2) and K1 > 0) else np.nan
    dt  = time.time() - tn

    return alpha_K, K1, psi1, g1, K2, psi2, g2, r21, dt


# ── Glowny blok (uruchamiany tylko w procesie macierzystym) ───────────────────

if __name__ == '__main__':
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer,
                                      encoding='utf-8', errors='replace')

    from concurrent.futures import ProcessPoolExecutor, as_completed

    t0 = time.time()
    print("="*72)
    print("p104: Skan alpha_K -> K*_2/K*_1(alpha_K)")
    print(f"Pytanie: czy r21 = {R21_TARGET} dla jakiegos alpha_K in [{ALPHA_VALUES[0]:.1f},{ALPHA_VALUES[-1]:.1f}]?")
    print(f"Liczba punktow: {len(ALPHA_VALUES)}  |  Watki: {N_WORKERS}")
    print(f"K*_1 okno psi: [{PSI1_LO},{PSI1_HI}], K log:[1e-4, 0.50], N_K={N_K1}")
    print(f"K*_2 okno psi: [{PSI2_LO},{PSI2_HI}] (MONO), K log:[1e-4, 5.00], N_K={N_K2}")
    print(f"Ref (a=8.5616): K*_1=0.010280, K*_2=0.100276, r21=9.7547")
    print("="*72)
    sys.stdout.flush()

    results = {}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(compute_one_alpha, ak): ak
                   for ak in ALPHA_VALUES}
        for future in as_completed(futures):
            ak = futures[future]
            try:
                res = future.result()
                results[ak] = res
                _, K1, psi1, g1, K2, psi2, g2, r21, dt = res
                K1s  = f"{K1:.6f}"  if np.isfinite(K1)  else "   --   "
                K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "   --   "
                r21s = f"{r21:.4f}" if np.isfinite(r21) else "   --  "
                flag = ""
                if np.isfinite(r21):
                    if r21 >= R21_TARGET * 0.80:
                        flag = f"  *** BLISKI {R21_TARGET}!"
                    elif r21 >= R21_TARGET * 0.50:
                        flag = f"  ** >50% celu"
                print(f"  a={ak:7.4f}  K1={K1s}  psi1={psi1:.3f}  "
                      f"K2={K2s}  psi2={psi2:.3f}  r21={r21s}  [{dt:.0f}s]{flag}")
                sys.stdout.flush()
            except Exception as e:
                print(f"  a={ak:.4f}  BLAD: {e}")
                sys.stdout.flush()

    # ── Tabela posortowana ─────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA WYNIKOW (posortowana po alpha_K)")
    print(f"  {'alpha_K':>8}  {'K*_1':>10}  {'psi1':>5}  "
          f"{'K*_2':>10}  {'psi2':>5}  {'r21':>9}  g1        g2")
    print("  " + "-"*72)

    r21_valid = []
    for ak in sorted(results.keys()):
        _, K1, psi1, g1, K2, psi2, g2, r21, _ = results[ak]
        K1s  = f"{K1:.6f}"  if np.isfinite(K1)  else "   --   "
        K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "   --   "
        p1s  = f"{psi1:.3f}" if np.isfinite(psi1) else "  -- "
        p2s  = f"{psi2:.3f}" if np.isfinite(psi2) else "  -- "
        r21s = f"{r21:.5f}" if np.isfinite(r21) else "    --   "
        g1s  = f"{g1:.2e}"  if np.isfinite(g1)  else "  nan "
        g2s  = f"{g2:.2e}"  if np.isfinite(g2)  else "  nan "
        star = " <-- BLISKI CELU" if (np.isfinite(r21) and r21 >= R21_TARGET*0.5) else ""
        print(f"  {ak:>8.4f}  {K1s}  {p1s}  {K2s}  {p2s}  {r21s}  {g1s}  {g2s}{star}")
        if np.isfinite(r21):
            r21_valid.append((ak, r21, K1, K2))

    # ── Analiza trendu ─────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print(f"ANALIZA: r21(alpha_K)")
    print(f"  Znaleziono r21 dla {len(r21_valid)}/{len(ALPHA_VALUES)} wartosci alpha_K")

    if r21_valid:
        r_vals = [r for _, r, _, _ in r21_valid]
        a_vals = [a for a, _, _, _ in r21_valid]

        print(f"  Zakres r21: [{min(r_vals):.5f}, {max(r_vals):.5f}]")
        print(f"  Srednia: {np.mean(r_vals):.5f},  Std: {np.std(r_vals):.5f}")

        # Trend
        if len(r21_valid) >= 2:
            mono_inc = all(r21_valid[i+1][1] >= r21_valid[i][1]-0.01
                          for i in range(len(r21_valid)-1))
            mono_dec = all(r21_valid[i+1][1] <= r21_valid[i][1]+0.01
                          for i in range(len(r21_valid)-1))
            trend = "ROSNACY" if mono_inc else ("MALEJACY" if mono_dec else "NIEMONOTONICZNY")
            print(f"  Trend: {trend}")

            # Fit liniowy r21(alpha_K)
            try:
                c = np.polyfit(a_vals, r_vals, 1)
                print(f"  Fit liniowy: r21 = {c[0]:.4f}*alpha_K + {c[1]:.4f}")
                alpha_extrap = (R21_TARGET - c[1]) / c[0]
                print(f"  Ekstrapolacja liniowa: r21={R21_TARGET} przy alpha_K = {alpha_extrap:.2f}")
            except Exception:
                pass

            # Fit potegowy r21 ~ alpha_K^beta
            try:
                pos = [(a, r) for a, r in zip(a_vals, r_vals) if a > 0 and r > 0]
                if len(pos) >= 3:
                    loga = np.log([a for a, r in pos])
                    logr = np.log([r for a, r in pos])
                    c2   = np.polyfit(loga, logr, 1)
                    beta, lnA = c2
                    A = np.exp(lnA)
                    print(f"  Fit potegowy: r21 = {A:.4f} * alpha_K^{beta:.4f}")
                    alpha_pow = (R21_TARGET / A)**(1.0/beta)
                    print(f"  Ekstrapolacja potegowa: r21={R21_TARGET} przy alpha_K = {alpha_pow:.2f}")
            except Exception:
                pass

        # Czy cel osiagniety?
        if max(r_vals) >= R21_TARGET:
            print(f"\n  *** CEL {R21_TARGET} OSIAGNIETY! ***")
            for i in range(len(r21_valid)-1):
                a0, r0 = r21_valid[i][0], r21_valid[i][1]
                a1, r1 = r21_valid[i+1][0], r21_valid[i+1][1]
                if (r0 - R21_TARGET) * (r1 - R21_TARGET) < 0:
                    print(f"  Przejscie: alpha_K in [{a0:.4f}, {a1:.4f}]")
                    alpha_cross = a0 + (R21_TARGET - r0)*(a1-a0)/(r1-r0)
                    print(f"  Interpolacja liniowa: alpha_K ~ {alpha_cross:.4f}")
        else:
            frac = max(r_vals) / R21_TARGET
            print(f"\n  Cel {R21_TARGET} NIEDOSTEPNY w zakresie [{ALPHA_VALUES[0]},{ALPHA_VALUES[-1]}]")
            print(f"  Max r21 = {max(r_vals):.5f}  ({frac:.2%} celu)")

    # ── PASS/FAIL ──────────────────────────────────────────────────────────────
    print("\nPASS/FAIL:")
    ref_ak = min(results.keys(), key=lambda a: abs(a - 8.5616))
    ref    = results[ref_ak]
    K1_ref, r21_ref = ref[1], ref[7]

    P1 = np.isfinite(K1_ref)
    P2 = P1 and abs(K1_ref / 0.010280 - 1.0) < 0.02
    P3 = len(r21_valid) >= int(len(ALPHA_VALUES) * 0.7)
    P4 = len(r21_valid) >= 2 and np.std(r_vals)/np.mean(r_vals) > 0.01  # ratio varies
    P5 = len(r21_valid) >= 2 and max(r_vals) / min(r_vals) > 2.0         # wide range
    P6 = len(r21_valid) >= 1 and max([r for _, r, _, _ in r21_valid]) >= R21_TARGET

    for i, (P, d) in enumerate([
        (P1, "K*_1 znaleziony przy alpha_K=8.5616"),
        (P2, f"K*_1(8.5616) stabilny (ref=0.010280, tol=2%)"),
        (P3, f"r21 znaleziony dla >=70% punktow alpha_K"),
        (P4, "r21 zmienia sie z alpha_K (std/mean > 1%)"),
        (P5, "r21 zakres >= 2x (szeroka dynamika)"),
        (P6, f"Cel r21={R21_TARGET} osiagniety"),
    ], 1):
        print(f"  P{i}: {'PASS' if P else 'FAIL'} -- {d}")

    total = time.time() - t0
    passed = sum([P1, P2, P3, P4, P5, P6])
    print(f"\nWynik: {passed}/6 PASS  (czas: {total/60:.1f} min)")
