#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p105_v2_bisection.py  --  TGP v1
=================================
Precyzyjne wyznaczenie alpha_SN (bifurkacja saddle-node galezi B)
metodą bisekcji binarnej.

Z p105:  alive=8.9616 (K*2=0.016488, psi=1.3115), dead=9.0116
=> szukamy alpha_SN w (8.9616, 9.0116) z precyzja 0.001

Potem: testujemy hipoteze  alpha_K / alpha_SN = n_s = 0.9649
  => przewidywane alpha_SN = 8.5616 / 0.9649 = 8.8733
  => oczekujemy: NIE zostanie potwierdzone (wynik p105 wskazuje ~8.99)

Watki: ThreadPoolExecutor(16) na wewnetrzna petle K-gridu

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

# ── Stale ─────────────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 600
RTOL   = 1e-9
ATOL   = 1e-11

N_S_PLANCK = 0.9649   # Planck 2018
AK_REF     = 8.5616
AK_SN_PRED = AK_REF / N_S_PLANCK  # = 8.8733

# Nasiona z p105 (ostatni alive + pierwsze dead)
AK_ALIVE  = 8.9616
K2_ALIVE  = 0.016488
PSI2_ALIVE = 1.3115
AK_DEAD   = 9.0116

# Parametry wyszukiwania
K_FRAC   = 2.5    # okno K: [K_prev/K_FRAC, K_prev*K_FRAC]
PSI_DELT = 0.5    # okno psi: psi_prev +/- PSI_DELT
N_K      = 40     # wieksza siatka K dla precyzji
N_PSI    = 120    # wieksza siatka psi

N_WORKERS = 16
BISECT_TOL = 0.001  # docelowa precyzja alpha_SN


# ── Funkcje ODE ───────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6


def make_rhs(alpha_K):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
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


# ── Rownolegla g_U(K) ─────────────────────────────────────────────────────────

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


def find_K2star(alpha_K, K2_prev, psi2_prev, rhs_fn, verbose=False):
    """Szuka K*_2 przy alpha_K, adaptacyjne okno wokol poprzedniego."""
    K_lo   = max(K2_prev / K_FRAC, 5e-4)
    K_hi   = min(K2_prev * K_FRAC, 5.0)
    psi_lo = max(psi2_prev - PSI_DELT, 1.001)
    psi_hi = min(psi2_prev + PSI_DELT, 6.0)

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

    if verbose:
        valid = [(Ks[i], gv[i], psv[i]) for i in range(len(Ks))
                 if np.isfinite(gv[i])]
        print(f"    [debug] {len(valid)} valid K points, K_lo={K_lo:.4f} K_hi={K_hi:.4f}")

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
                Kstar = brentq(gfunc, Ks[i], Ks[i+1], xtol=1e-6, maxiter=50)
                return Kstar, _last[0] if np.isfinite(_last[0]) else psi2_prev, _last[1], _last[2]
        except Exception:
            Kstar = Ks[i] - gv[i]*(Ks[i+1]-Ks[i]) / (gv[i+1]-gv[i])
            _, ps, g2, E2 = compute_gU_for_K((Kstar, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn))
            return Kstar, ps if np.isfinite(ps) else psi2_prev, g2, E2

    return np.nan, np.nan, np.nan, np.nan


# ── Bisekcja binarna alpha_SN ──────────────────────────────────────────────────

def test_alive(alpha_K, K2_prev, psi2_prev):
    """Sprawdza, czy galaz B istnieje przy alpha_K. Zwraca (alive, K2, psi2, g)."""
    t0     = time.time()
    rhs_fn = make_rhs(alpha_K)
    K2_new, psi2_new, g2, _ = find_K2star(alpha_K, K2_prev, psi2_prev, rhs_fn)
    dt     = time.time() - t0

    if np.isfinite(K2_new) and abs(g2) < 0.1:
        return True, K2_new, psi2_new, g2, dt
    else:
        return False, K2_prev, psi2_prev, np.nan, dt


if __name__ == '__main__':
    t0_total = time.time()
    print("="*72)
    print("p105_v2: Bisekcja binarna alpha_SN (bifurkacja saddle-node galezi B)")
    print(f"Hipoteza: alpha_K / alpha_SN = n_s = {N_S_PLANCK}")
    print(f"  alpha_K  = {AK_REF}")
    print(f"  Przewidywane alpha_SN = {AK_SN_PRED:.6f}")
    print(f"  Poczatkowy przedział: [{AK_ALIVE}, {AK_DEAD}]")
    print(f"  Docelowa precyzja: {BISECT_TOL}")
    print(f"Watki: {N_WORKERS}  N_K={N_K}  N_PSI={N_PSI}")
    print("="*72)

    # Stan poczatkowy: alive=8.9616, dead=9.0116
    lo     = AK_ALIVE   # ostatni alive
    hi     = AK_DEAD    # pierwsze dead

    K2_lo  = K2_ALIVE
    psi2_lo = PSI2_ALIVE

    print(f"\nNasiona bisekcji:")
    print(f"  ALIVE: alpha_K={lo:.4f}  K2={K2_lo:.6f}  psi2={psi2_lo:.4f}")
    print(f"  DEAD:  alpha_K={hi:.4f}")
    print(f"\nBisekcja:")

    history = []
    step = 0

    while (hi - lo) > BISECT_TOL:
        step += 1
        mid    = round((lo + hi) / 2.0, 6)
        print(f"  Krok {step}: testuję alpha_K={mid:.6f}  (przedział {hi-lo:.5f}) ...", end='', flush=True)

        alive, K2_mid, psi2_mid, g_mid, dt = test_alive(mid, K2_lo, psi2_lo)

        if alive:
            lo      = mid
            K2_lo   = K2_mid
            psi2_lo = psi2_mid
            status  = f"ALIVE  K2={K2_mid:.6f}  psi2={psi2_mid:.4f}  g={g_mid:.2e}"
        else:
            hi      = mid
            status  = "DEAD"

        print(f"  {status}  [{dt:.0f}s]")
        history.append((mid, alive, K2_mid if alive else np.nan,
                        psi2_mid if alive else np.nan, g_mid, dt))

    # Wynik bisekcji
    ak_SN = (lo + hi) / 2.0
    ak_SN_lo = lo
    ak_SN_hi = hi

    print(f"\n{'='*72}")
    print(f"WYNIK BISEKCJI:")
    print(f"  Ostatni ALIVE: alpha_K = {lo:.6f}  K2={K2_lo:.6f}")
    print(f"  Pierwsze DEAD: alpha_K = {hi:.6f}")
    print(f"  alpha_SN = ({lo:.6f}, {hi:.6f})")
    print(f"  Midpoint: alpha_SN ~ {ak_SN:.6f}  +/- {(hi-lo)/2:.6f}")

    # ── Test hipotezy ──────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("TEST HIPOTEZY: alpha_K / alpha_SN = n_s ?")
    print(f"  alpha_K      = {AK_REF}")
    print(f"  n_s (Planck) = {N_S_PLANCK}")
    print(f"  Przewidywane alpha_SN = {AK_SN_PRED:.6f}")
    print(f"  Zmierzone    alpha_SN = {ak_SN:.6f} +/- {(hi-lo)/2:.6f}")
    print()

    ratio_meas = AK_REF / ak_SN
    dev_meas   = abs(ratio_meas - N_S_PLANCK) / N_S_PLANCK * 100
    ratio_lo   = AK_REF / hi
    ratio_hi   = AK_REF / lo
    print(f"  alpha_K/alpha_SN (midpoint)  = {ratio_meas:.6f}  "
          f"vs n_s={N_S_PLANCK}  delta={dev_meas:.3f}%")
    print(f"  Przedzial ratios: [{ratio_lo:.6f}, {ratio_hi:.6f}]  "
          f"(n_s={'w' if ratio_lo <= N_S_PLANCK <= ratio_hi else 'POZA'} przedzialem)")
    print()

    # Alternatywne interpretacje
    print("Alternatywne hipotezy dla zmierzonego ratio:")
    ratio_actual = ratio_meas
    cands = [
        ("n_s Planck 2018", 0.9649, "spektralny indeks CMB"),
        ("n_s Planck +sigma", 0.9670, "CMB +1sigma"),
        ("n_s Planck -sigma", 0.9609, "CMB -1sigma"),
        ("1 - 1/2pi^2",      1.0 - 1.0/(2*np.pi**2), "korekta 1/(2pi^2)"),
        ("1 - 3/(4pi^2)",    1.0 - 3.0/(4*np.pi**2), "korekta 3/(4pi^2)"),
        ("alpha_K/9",        AK_REF/9.0, "alpha_K/9"),
        ("closest_integer",  AK_REF / round(ak_SN * 10) * 10, "zaokr. alpha_SN"),
    ]
    print(f"  Zmierzone alpha_K/alpha_SN = {ratio_actual:.6f}")
    for name, val, desc in cands:
        d = abs(ratio_actual - val) / val * 100
        marker = " <-- NAJBLIZEJ" if d < 0.5 else ""
        print(f"  {name:30s} = {val:.6f}  delta={d:.3f}%  ({desc}){marker}")

    if dev_meas < 1.0:
        verdict = "POTWIERDZONA (delta < 1%)"
    elif dev_meas < 3.0:
        verdict = "CZESCIOWO POTWIERDZONA (delta < 3%)"
    else:
        verdict = "OBALONA (delta >= 3%)"
    print(f"\n  WERDYKT: HIPOTEZA {verdict}")
    print(f"  WNIOSEK: alpha_SN = {ak_SN:.4f}, alpha_K/alpha_SN = {ratio_meas:.4f}, n_s = {N_S_PLANCK}")

    # ── Historia bisekcji ──────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("HISTORIA BISEKCJI:")
    print(f"  {'krok':>4}  {'alpha_K':>10}  {'wynik':>6}  {'K*_2':>10}  {'psi2':>6}  {'g_U':>10}  czas")
    print("  " + "-"*60)
    for i, (ak, al, K2v, psi2v, gv, dtv) in enumerate(history):
        K2s  = f"{K2v:.6f}" if np.isfinite(K2v) else "   --   "
        p2s  = f"{psi2v:.4f}" if np.isfinite(psi2v) else "  --  "
        gs   = f"{gv:.3e}" if np.isfinite(gv) else "  nan  "
        res  = "ALIVE" if al else " DEAD"
        print(f"  {i+1:>4}  {ak:>10.6f}  {res}  {K2s}  {p2s}  {gs}  {dtv:.0f}s")

    print(f"\n[Czas calkowity: {(time.time()-t0_total)/60:.1f} min]")
    print("="*72)
