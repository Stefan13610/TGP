#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p105_alpha_SN_bifurcation.py  --  TGP v1
==========================================
Precyzyjne wyznaczenie punktu bifurkacji saddle-node galezi B:
  alpha_SN = inf{ alpha_K : galaz B (K*_2) nie istnieje }

Hipoteza: alpha_K / alpha_SN = n_s = 0.9649  (spektralny indeks CMB!)
  => alpha_SN = alpha_K / n_s = 8.5616 / 0.9649 = 8.8733

Weryfikacja: jesli tak -> alpha_K jest zdeterminowane przez
  PRODUKT bifurkacji solitonowej (alpha_SN z V_mod) i CMB (n_s).

Metodologia:
  - Kontynuacja galezi B od alpha_K=8.5616 w gore (krok 0.05)
  - Dla kazdego alpha_K: ciasne okno wokol poprzedniego rozwiazania
  - Wykryj alpha_SN = punkt znikniecia K*_2
  - Sprawdz: alpha_K / alpha_SN vs n_s = 0.9649
  - Rownolegle watki (ThreadPoolExecutor, 16) na wewnetrzna petle K-gridu

Data: 2026-03-26
"""
import sys
import io
import warnings
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from concurrent.futures import ThreadPoolExecutor, as_completed

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Stale ────────────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 600          # wiecej dla precyzji
RTOL   = 1e-9
ATOL   = 1e-11
R21_TARGET = 206.77
N_S_PLANCK = 0.9649   # Planck 2018 srednia

# Punkt startowy (p103_v4, potwierdzone)
AK_REF  = 8.5616
K2_REF  = 0.100276
PSI2_REF = 2.7676

# Kontynuacja w gore
DALPHA  = 0.05
AK_SCAN_HI = 9.20

# Okno adaptacyjne
K_FRAC   = 2.5
PSI_DELT = 0.6
N_K      = 30
N_PSI    = 100

N_WORKERS = 16

# Przewidywana wartość alpha_SN
AK_SN_PRED = AK_REF / N_S_PLANCK


# ── Funkcje ODE ───────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6


def make_rhs(alpha_K):
    def rhs(r, y):
        phi, dphi = y
        phi = max(phi, 1e-10)
        kfac = 1.0 + alpha_K / phi
        dV   = GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1.0)**5
        return [dphi,
                dV/kfac + alpha_K*dphi**2 / (2.0*phi**2*kfac) - 2.0/r*dphi]
    return rhs


def solve_full(psi, K, rhs_fn):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(rhs_fn, [AGAM, R_MAX], [psi, dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] < R_MAX * 0.99:
                if method == 'DOP853': continue
                return None
            return sol
        except Exception:
            if method == 'DOP853': continue
            return None
    return None


def phi_end(psi, K, rhs_fn):
    sol = solve_full(psi, K, rhs_fn)
    return float(sol.y[0, -1]) if sol else np.nan


def is_mono(sol, tol=0.01):
    return bool(np.all(np.diff(sol.y[0]) <= tol))


def energy_sol(sol, alpha_K):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))):
        return np.nan
    kfac = 1.0 + alpha_K / phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vpot = np.clip([V_mod(float(p)) for p in phi], -1e12, 1e12)
    Ep   = 4*np.pi*np.trapezoid((np.array(Vpot) - V1)*r**2, r)
    E    = Ek + Ep
    return E if np.isfinite(E) else np.nan


# ── Rownolegle g_U(K) ─────────────────────────────────────────────────────────

def compute_gU_for_K(args):
    K, alpha_K, psi_lo, psi_hi, n_psi, rhs_fn = args
    warnings.filterwarnings('ignore')
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_end(p, K, rhs_fn) - 1.0 for p in psis]
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1) and fi*fi1 < 0):
            continue
        try:
            pz  = brentq(lambda p: phi_end(p, K, rhs_fn) - 1.0,
                         psis[i], psis[i+1], xtol=1e-8, maxiter=60)
            sol = solve_full(pz, K, rhs_fn)
            if sol and is_mono(sol):
                E = energy_sol(sol, alpha_K)
                g = E/(4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan
                return K, pz, g, E
        except Exception:
            pass
    return K, np.nan, np.nan, np.nan


def find_K2star(alpha_K, K2_prev, psi2_prev, rhs_fn):
    """Znajdz K*_2 (galaz B) przy alpha_K, okno adaptacyjne wokol poprzedniego."""
    K_lo   = max(K2_prev / K_FRAC, 1e-5)
    K_hi   = min(K2_prev * K_FRAC, 5.0)
    psi_lo = max(psi2_prev - PSI_DELT, 1.001)
    psi_hi = min(psi2_prev + PSI_DELT, 7.0)

    Ks   = np.logspace(np.log10(K_lo), np.log10(K_hi), N_K)
    args = [(K, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn) for K in Ks]

    gv   = [np.nan] * len(Ks)
    psv  = [np.nan] * len(Ks)
    Ev   = [np.nan] * len(Ks)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(compute_gU_for_K, a): i for i, a in enumerate(args)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                K_out, pz, g, E = fut.result()
                gv[i]  = g
                psv[i] = pz
                Ev[i]  = E
            except Exception:
                pass

    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1) and gi*gi1 < 0):
            continue
        _last = [np.nan, np.nan, np.nan]

        def gfunc(K, _l=_last, _pl=psi_lo, _ph=psi_hi, _r=rhs_fn, _ak=alpha_K):
            _, pz2, g2, E2 = compute_gU_for_K((K, _ak, _pl, _ph, N_PSI, _r))
            _l[0] = pz2; _l[1] = g2; _l[2] = E2
            return g2 if np.isfinite(g2) else np.nan

        try:
            ga = gfunc(Ks[i]); gb = gfunc(Ks[i+1])
            if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                Kstar    = brentq(gfunc, Ks[i], Ks[i+1], xtol=1e-6, maxiter=50)
                psi_star = _last[0]
                E_star   = _last[2]
                g_star   = _last[1]
                return Kstar, psi_star if np.isfinite(psi_star) else psi2_prev, g_star, E_star
        except Exception:
            Kstar = Ks[i] - gv[i]*(Ks[i+1]-Ks[i]) / (gv[i+1]-gv[i])
            _, ps, g2, E2 = compute_gU_for_K((Kstar, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn))
            return Kstar, ps if np.isfinite(ps) else psi2_prev, g2, E2

    return np.nan, np.nan, np.nan, np.nan


# ── Glowna petla ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    t0 = time.time()
    print("="*72)
    print("p105: Wyznaczenie alpha_SN (bifurkacja saddle-node galezi B)")
    print(f"Hipoteza: alpha_K / alpha_SN = n_s = {N_S_PLANCK}")
    print(f"  alpha_K  = {AK_REF}")
    print(f"  alpha_SN = alpha_K / n_s = {AK_SN_PRED:.6f}  (przewidywana)")
    print(f"Skan: alpha_K in [{AK_REF},{AK_SCAN_HI}], krok={DALPHA}")
    print(f"Watki: {N_WORKERS}  N_K={N_K}  N_PSI={N_PSI}  r_max={R_MAX}")
    print("="*72)

    # Punkt startowy
    K2, psi2 = K2_REF, PSI2_REF
    scan_results = [(AK_REF, K2, psi2, 0.0, K2/0.010280)]
    last_alive = AK_REF

    print(f"  a={AK_REF:.4f} [ref]  K2={K2:.6f}  psi2={psi2:.4f}  "
          f"g≈0 [p103_v4]  r21={K2/0.010280:.4f}")

    ak = AK_REF
    consecutive_fails = 0

    while ak < AK_SCAN_HI:
        ak = round(ak + DALPHA, 4)
        tn = time.time()
        rhs_fn = make_rhs(ak)
        K2_new, psi2_new, g2_new, E2_new = find_K2star(ak, K2, psi2, rhs_fn)
        dt = time.time() - tn

        if np.isfinite(K2_new):
            consecutive_fails = 0
            K2   = K2_new
            psi2 = psi2_new
            last_alive = ak

            # K*_1 dla r21 (z prostym regresywnym K1 ~ 0.010280 * skalowanie)
            # Uzywamy K*_1 z referencji (zmienia sie powoli)
            K1_est = 0.010280 * np.exp(-0.060 * (ak - AK_REF))  # gruba estymata
            r21_est = K2 / K1_est

            g2s = f"{g2_new:.3e}" if np.isfinite(g2_new) else "  nan "
            print(f"  a={ak:.4f} [+]  K2={K2:.6f}  psi2={psi2:.4f}  "
                  f"g={g2s}  [{dt:.0f}s]")
            scan_results.append((ak, K2, psi2, g2_new, r21_est))
        else:
            consecutive_fails += 1
            print(f"  a={ak:.4f} [!]  K*_2 BRAK  [{dt:.0f}s]  "
                  f"(kolejny brak: {consecutive_fails})")
            scan_results.append((ak, np.nan, np.nan, np.nan, np.nan))
            if consecutive_fails >= 3:
                print(f"  -> Galaz B znikla definitywnie.")
                break

    # ── Wyznaczenie alpha_SN ───────────────────────────────────────────────────
    print("\n" + "="*72)

    # Ostatni alpha_K z istniejacym K*_2
    alive = [(ak, K2, psi2, g) for ak, K2, psi2, g, _ in scan_results
             if np.isfinite(K2)]
    dead  = [(ak, K2, psi2, g) for ak, K2, psi2, g, _ in scan_results
             if not np.isfinite(K2)]

    if alive and dead:
        ak_last  = alive[-1][0]
        ak_first_dead = dead[0][0]

        # alpha_SN ~ midpoint interpolacji
        # Bardziej precyzyjnie: interpolacja liniowa K2(alpha_K) -> K2=0
        aks_alive  = [a for a, _, _, _ in alive]
        K2s_alive  = [k for _, k, _, _ in alive]

        # Fit liniowy K*_2(alpha_K) ~ a + b*(alpha_K - ak_last)
        if len(alive) >= 2:
            try:
                c = np.polyfit(aks_alive[-3:], K2s_alive[-3:], 1)
                ak_SN_interp = -c[1] / c[0]  # K*_2 -> 0
                print(f"Fit liniowy K*_2(a): K*_2 = {c[0]:.6f}*(a) + {c[1]:.6f}")
                print(f"Interpolacja K*_2=0: alpha_SN ~ {ak_SN_interp:.5f}")
            except Exception:
                ak_SN_interp = (ak_last + ak_first_dead) / 2

        ak_SN_bracket = f"alpha_SN in ({ak_last:.4f}, {ak_first_dead:.4f})"
        ak_SN_mid     = (ak_last + ak_first_dead) / 2.0

        print(f"\nBifurkacja saddle-node galezi B:")
        print(f"  Ostatni alpha_K z K*_2: {ak_last:.4f}")
        print(f"  Pierwszy alpha_K bez K*_2: {ak_first_dead:.4f}")
        print(f"  Przedział: {ak_SN_bracket}")
        print(f"  Midpoint: alpha_SN ~ {ak_SN_mid:.5f}")
        if len(alive) >= 2:
            print(f"  Interpolacja liniowa: alpha_SN ~ {ak_SN_interp:.5f}")

        # ── Test hipotezy alpha_K / alpha_SN = n_s ─────────────────────────────
        print("\n" + "="*72)
        print("TEST HIPOTEZY: alpha_K / alpha_SN = n_s ?")
        print(f"  alpha_K      = {AK_REF}")
        print(f"  n_s (Planck) = {N_S_PLANCK}")
        print(f"  Przewidywane alpha_SN = alpha_K / n_s = {AK_SN_PRED:.6f}")
        print()

        for label, ak_SN_test in [
            ("midpoint",   ak_SN_mid),
            ("interpolacja liniowa", ak_SN_interp if len(alive)>=2 else ak_SN_mid),
            ("przewidywana",          AK_SN_PRED),
        ]:
            ratio = AK_REF / ak_SN_test
            dev   = abs(ratio - N_S_PLANCK) / N_S_PLANCK * 100
            in_bracket = ak_last < ak_SN_test < ak_first_dead
            bracket_ok = "✓ w przedziale" if in_bracket else "✗ poza"
            print(f"  [{label}]  alpha_SN={ak_SN_test:.5f}  "
                  f"alpha_K/alpha_SN={ratio:.6f}  "
                  f"vs n_s={N_S_PLANCK}  "
                  f"delta={dev:.3f}%  {bracket_ok}")

        print()
        ratio_mid  = AK_REF / ak_SN_mid
        delta_mid  = abs(ratio_mid - N_S_PLANCK) / N_S_PLANCK * 100

        if delta_mid < 1.0:
            print(f"  HIPOTEZA POTWIERDZONA (delta < 1%): alpha_K/alpha_SN = n_s")
            print(f"  WNIOSEK: alpha_K = n_s * alpha_SN(V_mod, a_Gam)")
            print(f"  IMPLIKACJA: alpha_K jest zdeterminowane przez bifurkacje (alpha_SN)")
            print(f"              i spektralny indeks CMB (n_s) LACZNIE.")
        elif delta_mid < 5.0:
            print(f"  HIPOTEZA CZESCIOWO POTWIERDZONA (delta < 5%): wymaga precizyj.")
        else:
            print(f"  HIPOTEZA OBALONA (delta >= 5%): alpha_K/alpha_SN != n_s")

    else:
        print("Za mało danych do wyznaczenia alpha_SN.")
        print(f"Galaz B: alive={len(alive)}, dead={len(dead)}")

    # ── Tabela K*_2(alpha_K) ──────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA K*_2(alpha_K)")
    print(f"  {'alpha_K':>8}  {'K*_2':>10}  {'psi2':>6}  {'g_U':>10}")
    print("  " + "-"*42)
    for ak, K2, psi2, g, _ in scan_results:
        K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "   --   "
        p2s  = f"{psi2:.4f}" if np.isfinite(psi2) else "  --  "
        gs   = f"{g:.3e}"   if np.isfinite(g)   else "  nan  "
        print(f"  {ak:>8.4f}  {K2s}  {p2s}  {gs}")

    print(f"\n[Czas: {(time.time()-t0)/60:.1f} min]")
