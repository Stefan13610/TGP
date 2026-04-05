#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p104_v2_continuation.py  --  TGP v1
=====================================
Kontynuacja gałęzi A (K*_1) i B (K*_2) wzgledem alpha_K.

Strategia: startujemy z POTWIERDZONYCH wartosci przy alpha_K=8.5616 (p103_v4):
  K*_1 = 0.010280, psi_0^(1) = 1.2417  (galaz A)
  K*_2 = 0.100276, psi_0^(2) = 2.7676  (galaz B, mono filter)
i kroczymy po alpha_K z krokiem DALPHA, adaptujac okna wokol poprzedniego rozwiazania.

Cel: K*_2/K*_1(alpha_K) - czy moze osiagnac 206.77?

Rownoleglenie: ThreadPoolExecutor(16) wewnatrz petli g_U(K) (wewnetrzny K-grid).
Kontynuacja: sekwencyjna w alpha_K, dwie galęzie (rosnaca i malejaca od 8.5616).

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
N_EVAL = 500
RTOL   = 1e-8
ATOL   = 1e-10
R21_TARGET = 206.77

# Punkt startowy (z p103_v4, potwierdzone)
AK_REF  = 8.5616
K1_REF  = 0.010280;  PSI1_REF = 1.2417
K2_REF  = 0.100276;  PSI2_REF = 2.7676

# Kroki kontynuacji
DALPHA  = 0.5                       # krok alpha_K
AK_MIN  = 3.0;  AK_MAX = 16.0      # zasieg skanu
N_WORKERS = 16                      # watki do rownoleglenia K-gridu

# Okno adaptacyjne (fraction otaczajaca poprzednie rozwiazanie)
K_FRAC   = 3.0    # K in [K_prev/K_FRAC, K_prev*K_FRAC]
PSI_DELT = 0.8    # psi in [psi_prev - DELT, psi_prev + DELT]
N_K   = 30        # punkty K
N_PSI = 80        # punkty psi


# ── ODE / energie ─────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6


def make_rhs(alpha_K):
    """Zwraca funkcje RHS ODE dla danego alpha_K (domkniecie)."""
    def rhs(r, y):
        phi, dphi = y
        phi = max(phi, 1e-10)
        kfac = 1.0 + alpha_K / phi
        dV   = GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1.0)**5
        return [dphi,
                dV/kfac + alpha_K*dphi**2/(2.0*phi**2*kfac) - 2.0/r*dphi]
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
    Vpot = np.clip(np.array([V_mod(float(p)) for p in phi]), -1e12, 1e12)
    Ep   = 4*np.pi*np.trapezoid((Vpot - V1)*r**2, r)
    E    = Ek + Ep
    return E if np.isfinite(E) else np.nan


# ── Szukanie K* przy danym alpha_K i oknie adaptacyjnym ──────────────────────

def compute_gU_for_K(args):
    """Pomocnicza - oblicza g_U dla jednego K. Uruchamiana w watku."""
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
                         psis[i], psis[i+1], xtol=1e-7, maxiter=50)
            sol = solve_full(pz, K, rhs_fn)
            if sol and is_mono(sol):
                E = energy_sol(sol, alpha_K)
                g = E/(4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan
                return K, pz, g
        except Exception:
            pass
    return K, np.nan, np.nan


def find_Kstar_adaptive(alpha_K, K_prev, psi_prev, rhs_fn):
    """
    Znajdz K* przy alpha_K uzywajac okna adaptacyjnego wokol (K_prev, psi_prev).
    Wewnetrzna petla K rownolewa na N_WORKERS watkach.
    """
    K_lo  = max(K_prev / K_FRAC, 1e-5)
    K_hi  = min(K_prev * K_FRAC, 10.0)
    psi_lo = max(psi_prev - PSI_DELT, 1.001)
    psi_hi = min(psi_prev + PSI_DELT, 7.0)

    Ks   = np.logspace(np.log10(K_lo), np.log10(K_hi), N_K)
    args = [(K, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn) for K in Ks]

    # Rownolegle obliczenie g_U(K) na N_WORKERS watkach
    gv  = [np.nan] * len(Ks)
    psv = [np.nan] * len(Ks)
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(compute_gU_for_K, a): i for i, a in enumerate(args)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                K_out, pz, g = fut.result()
                gv[i]  = g
                psv[i] = pz
            except Exception:
                pass

    # Szukaj znakowej zmiany g_U(K)
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1) and gi*gi1 < 0):
            continue
        _last = [np.nan]

        def gfunc(K, _l=_last, _plo=psi_lo, _phi=psi_hi, _rhs=rhs_fn, _ak=alpha_K):
            _, pz2, g2 = compute_gU_for_K((K, _ak, _plo, _phi, N_PSI, _rhs))
            _l[0] = pz2
            return g2 if np.isfinite(g2) else np.nan

        try:
            ga = gfunc(Ks[i]); gb = gfunc(Ks[i+1])
            if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                Kstar    = brentq(gfunc, Ks[i], Ks[i+1], xtol=1e-5, maxiter=40)
                psi_star = _last[0]
                return Kstar, psi_star if np.isfinite(psi_star) else psi_prev
        except Exception:
            Kstar = Ks[i] - gv[i]*(Ks[i+1]-Ks[i]) / (gv[i+1]-gv[i])
            _, pstar, _ = compute_gU_for_K((Kstar, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn))
            return Kstar, pstar if np.isfinite(pstar) else psi_prev

    return np.nan, np.nan


def compute_one_step(alpha_K, K1_prev, psi1_prev, K2_prev, psi2_prev):
    """Jeden krok kontynuacji - oblicz K*_1 i K*_2 dla nowego alpha_K."""
    rhs_fn = make_rhs(alpha_K)
    tn = time.time()

    K1, p1 = find_Kstar_adaptive(alpha_K, K1_prev, psi1_prev, rhs_fn)
    K2, p2 = find_Kstar_adaptive(alpha_K, K2_prev, psi2_prev, rhs_fn)

    r21 = K2/K1 if (np.isfinite(K1) and np.isfinite(K2) and K1 > 0) else np.nan
    dt  = time.time() - tn
    return K1, p1, K2, p2, r21, dt


# ── Glowna petla kontynuacji ──────────────────────────────────────────────────

def run_continuation(direction, alpha_start, K1_start, psi1_start,
                     K2_start, psi2_start, results):
    """
    Kontynuacja galezi w zadanym kierunku (+1 lub -1).
    Modyfikuje slownik results in-place.
    """
    K1, p1 = K1_start, psi1_start
    K2, p2 = K2_start, psi2_start
    ak = alpha_start

    while True:
        ak = round(ak + direction * DALPHA, 4)
        if ak < AK_MIN or ak > AK_MAX:
            break

        K1_new, p1_new, K2_new, p2_new, r21, dt = compute_one_step(ak, K1, p1, K2, p2)

        K1s  = f"{K1_new:.6f}"  if np.isfinite(K1_new)  else "   --   "
        K2s  = f"{K2_new:.6f}"  if np.isfinite(K2_new)  else "   --   "
        r21s = f"{r21:.5f}"     if np.isfinite(r21)     else "   --   "
        flag = ""
        if np.isfinite(r21):
            if r21 >= R21_TARGET * 0.8: flag = "  *** BLISKI CELU!"
            elif r21 >= R21_TARGET * 0.3: flag = "  ** >30% celu"
        print(f"  a={ak:7.4f} [{'+' if direction>0 else '-'}]  "
              f"K1={K1s} psi1={p1_new:.3f}  K2={K2s} psi2={p2_new:.3f}  "
              f"r21={r21s}  [{dt:.0f}s]{flag}")
        sys.stdout.flush()

        results[ak] = (K1_new, p1_new, K2_new, p2_new, r21)

        # Aktualizuj punkty startowe (kontynuacja)
        if np.isfinite(K1_new):
            K1, p1 = K1_new, p1_new
        if np.isfinite(K2_new):
            K2, p2 = K2_new, p2_new

        # Zatrzymaj jesli galaz znikla
        if not (np.isfinite(K1_new) and np.isfinite(K2_new)):
            print(f"    -> Galaz znikla przy alpha_K={ak:.4f}, kontynuacja przerwana.")
            break


# ── Punkt startowy ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    t0 = time.time()
    print("="*72)
    print("p104_v2: Kontynuacja gałęzi A+B wzgl. alpha_K")
    print(f"Pytanie: czy r21 = K*_2/K*_1 osiaga {R21_TARGET} dla jakiegos alpha_K?")
    print(f"Start: alpha_K={AK_REF}, K*_1={K1_REF}, psi1={PSI1_REF}")
    print(f"       K*_2={K2_REF}, psi2={PSI2_REF}")
    print(f"Krok: DALPHA={DALPHA}, zakres: [{AK_MIN},{AK_MAX}]")
    print(f"Rownolegle watki (wew.): {N_WORKERS}  N_K={N_K}  N_PSI={N_PSI}")
    print("="*72)
    sys.stdout.flush()

    results = {}

    # Punkt startowy (potwierdzony przez p103_v4)
    results[AK_REF] = (K1_REF, PSI1_REF, K2_REF, PSI2_REF, K2_REF/K1_REF)
    print(f"  a={AK_REF:7.4f} [*]  K1={K1_REF:.6f} psi1={PSI1_REF:.4f}  "
          f"K2={K2_REF:.6f} psi2={PSI2_REF:.4f}  "
          f"r21={K2_REF/K1_REF:.5f}  [ref p103_v4]")
    sys.stdout.flush()

    print("\n--- Kierunek: alpha_K rosnace (od 8.5616 do 16.0) ---")
    run_continuation(+1, AK_REF, K1_REF, PSI1_REF, K2_REF, PSI2_REF, results)

    print("\n--- Kierunek: alpha_K malejace (od 8.5616 do 3.0) ---")
    run_continuation(-1, AK_REF, K1_REF, PSI1_REF, K2_REF, PSI2_REF, results)

    # ── Tabela i analiza ────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA WYNIKOW (posortowana po alpha_K)")
    print(f"  {'alpha_K':>8}  {'K*_1':>10}  {'psi1':>5}  "
          f"{'K*_2':>10}  {'psi2':>5}  {'r21':>9}")
    print("  " + "-"*65)

    r21_list = []
    for ak in sorted(results.keys()):
        K1, p1, K2, p2, r21 = results[ak]
        K1s  = f"{K1:.6f}"  if np.isfinite(K1)  else "   --   "
        K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "   --   "
        p1s  = f"{p1:.3f}"  if np.isfinite(p1)  else "  -- "
        p2s  = f"{p2:.3f}"  if np.isfinite(p2)  else "  -- "
        r21s = f"{r21:.5f}" if np.isfinite(r21) else "    --   "
        star = " <--" if (np.isfinite(r21) and r21 >= R21_TARGET*0.5) else ""
        print(f"  {ak:>8.4f}  {K1s}  {p1s}  {K2s}  {p2s}  {r21s}{star}")
        if np.isfinite(r21):
            r21_list.append((ak, r21))

    # Analiza
    print("\n" + "="*72)
    if r21_list:
        r_vals = [r for _, r in r21_list]
        a_vals = [a for a, _ in r21_list]
        print(f"Znaleziono r21 dla {len(r21_list)} wartosci alpha_K")
        print(f"Zakres r21: [{min(r_vals):.5f}, {max(r_vals):.5f}]")

        mono_inc = all(r21_list[i+1][1] >= r21_list[i][1]-0.1
                       for i in range(len(r21_list)-1))
        mono_dec = all(r21_list[i+1][1] <= r21_list[i][1]+0.1
                       for i in range(len(r21_list)-1))
        print(f"Trend: {'ROSNACY' if mono_inc else ('MALEJACY' if mono_dec else 'NIEMONOTONICZNY')}")

        if len(r_vals) >= 2:
            try:
                c = np.polyfit(a_vals, r_vals, 1)
                ak_target = (R21_TARGET - c[1]) / c[0]
                print(f"Fit liniowy: r21 = {c[0]:.4f}*a + {c[1]:.4f}")
                print(f"Ekstrapolacja liniowa: r21={R21_TARGET} przy alpha_K = {ak_target:.1f}")
            except Exception:
                pass

        if max(r_vals) >= R21_TARGET:
            print(f"\n*** CEL {R21_TARGET} OSIAGNIETY! ***")
        else:
            print(f"\nCel {R21_TARGET} NIEDOSTEPNY. Max r21 = {max(r_vals):.5f} "
                  f"({max(r_vals)/R21_TARGET:.2%} celu)")

    print(f"\n[Czas calkowity: {(time.time()-t0)/60:.1f} min]")
