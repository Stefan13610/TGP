#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p84_branch1_scan_parallel.py  --  TGP v1  *  Parallel scan minimum E* podgalezi B1
====================================================================================
Wersja pod multiprocessing:
  - parallel po alpha
  - brak sekwencyjnego prev_psi w coarse scanie
  - cache lokalny na proces
  - BLAS/OMP ustawione na 1 watek na proces
  - opcjonalny refine wokol minimum tez rownolegly

Uzycie:
  python p84_branch1_scan_parallel.py
  python p84_branch1_scan_parallel.py --workers 12
  python p84_branch1_scan_parallel.py --workers 16 --fine
  python p84_branch1_scan_parallel.py --workers 12 --alpha-min 5.0 --alpha-max 8.6 --alpha-step 0.05

Data: 2026-03-25
"""

import os

# Bardzo wazne przy multiprocessing + numpy/scipy:
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import io
import time
import math
import argparse
import warnings
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings("ignore")

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETRY FIZYCZNE
# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF = 8.5616
A_GAM = 0.040
K1_REF = 0.010414
LAM = 5.501357e-06
GAMMA_V = 1.0

# Zakres B1: E/K in [12.40, 12.65]
B1_EOK_MIN = 12.40
B1_EOK_MAX = 12.65

# Numeryka
RTOL = 1e-9
ATOL = 1e-11
PHI_FLOOR = 1e-10
EVENT_STOP = 1e-5
ROUND_DIGITS = 12

# Lokalne szukanie psi
LOCAL_WIDTH_DEFAULT = 0.20
LOCAL_NPSI_DEFAULT = 32

# Parametry domyslne skanu
DEFAULT_ALPHA_MIN = 5.0
DEFAULT_ALPHA_MAX = 8.60
DEFAULT_ALPHA_STEP = 0.05

DEFAULT_FINE_WINDOW = 0.50
DEFAULT_FINE_STEP = 0.02

DEFAULT_K_MIN = 4e-3
DEFAULT_K_MAX = 0.12
DEFAULT_NK = 25
DEFAULT_NPSI = 120
DEFAULT_RMAX = 40.0

# ─────────────────────────────────────────────────────────────────────────────
# POTENCJAL
# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):
    return GAMMA_V / 3 * p**3 - GAMMA_V / 4 * p**4 + LAM / 6 * (p - 1) ** 6

def dV_mod(p):
    return GAMMA_V * p**2 - GAMMA_V * p**3 + LAM * (p - 1) ** 5

V1 = GAMMA_V / 3 - GAMMA_V / 4

# ─────────────────────────────────────────────────────────────────────────────
# NARZĘDZIA
# ─────────────────────────────────────────────────────────────────────────────
def q(x, nd=ROUND_DIGITS):
    return round(float(x), nd)

def is_finite(x):
    return np.isfinite(x)

# ─────────────────────────────────────────────────────────────────────────────
# ODE
# ─────────────────────────────────────────────────────────────────────────────
def ode_rhs(r, y, alpha):
    phi, dphi = y
    phi = max(phi, PHI_FLOOR)
    kfac = 1.0 + alpha / phi
    return [
        dphi,
        dV_mod(phi) / kfac + alpha * dphi**2 / (2 * phi**2 * kfac) - 2 / r * dphi,
    ]

@lru_cache(maxsize=300000)
def _phi_at_rmax_cached(psi, K, a_gam, alpha, r_max):
    dphi0 = -K / a_gam**2

    def ev(r, y):
        return y[0] - EVENT_STOP

    ev.terminal = True
    ev.direction = -1

    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max],
            [psi, dphi0],
            method="DOP853",
            rtol=RTOL,
            atol=ATOL,
            events=[ev],
            dense_output=False,
        )
        if sol.t[-1] >= r_max * 0.99:
            return float(sol.y[0, -1])
        return np.nan
    except Exception:
        return np.nan

def _phi_at_rmax(psi, K, a_gam, alpha, r_max=40.0):
    return _phi_at_rmax_cached(q(psi), q(K), q(a_gam), q(alpha), q(r_max))

@lru_cache(maxsize=200000)
def _energy_at_psi_cached(psi_z, K, a_gam, alpha, r_max):
    dphi0 = -K / a_gam**2
    r_ev = a_gam * (r_max / a_gam) ** np.linspace(0, 1, 4000)

    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max],
            [psi_z, dphi0],
            method="DOP853",
            rtol=RTOL,
            atol=ATOL,
            t_eval=r_ev,
        )
        r = sol.t
        phi = np.maximum(sol.y[0], PHI_FLOOR)
        dphi = sol.y[1]

        Ek = 4 * np.pi * np.trapezoid(0.5 * dphi**2 * (1 + alpha / phi) * r**2, r)
        Ep = 4 * np.pi * np.trapezoid((V_mod(phi) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

def energy_at_psi(psi_z, K, a_gam, alpha, r_max=40.0):
    return _energy_at_psi_cached(q(psi_z), q(K), q(a_gam), q(alpha), q(r_max))

# ─────────────────────────────────────────────────────────────────────────────
# ROOT FINDING
# ─────────────────────────────────────────────────────────────────────────────
def _root_func_factory(K, a_gam, alpha, r_max):
    return lambda p: _phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0

def find_psi_zeros_in_range(K, a_gam, alpha, psi_lo, psi_hi, n_psi=40, r_max=40.0):
    psi_lo = max(1.001, float(psi_lo))
    psi_hi = max(psi_lo + 1e-6, float(psi_hi))

    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv = [_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0 for p in psis]

    roots = []
    froot = _root_func_factory(K, a_gam, alpha, r_max)

    for i in range(len(Fv) - 1):
        fi, fj = Fv[i], Fv[i + 1]
        if is_finite(fi) and is_finite(fj) and fi * fj < 0:
            try:
                pz = brentq(
                    froot,
                    psis[i],
                    psis[i + 1],
                    xtol=1e-7,
                    rtol=1e-7,
                    maxiter=60,
                )
                if not any(abs(pz - rr) < 1e-6 for rr in roots):
                    roots.append(pz)
            except Exception:
                pass

    return roots

def find_all_psi_zeros(K, a_gam, alpha, n_psi=150, r_max=40.0, psi_hint=None):
    if psi_hint is not None and is_finite(psi_hint):
        roots_local = find_psi_zeros_in_range(
            K,
            a_gam,
            alpha,
            psi_hint - LOCAL_WIDTH_DEFAULT,
            psi_hint + LOCAL_WIDTH_DEFAULT,
            n_psi=LOCAL_NPSI_DEFAULT,
            r_max=r_max,
        )
        if roots_local:
            return roots_local

    psi_max = max(3.0, 1.0 + 2.0 * K / a_gam)
    psis_hi = min(psi_max, 400.0)

    return find_psi_zeros_in_range(
        K,
        a_gam,
        alpha,
        1.001,
        psis_hi,
        n_psi=n_psi,
        r_max=r_max,
    )

# ─────────────────────────────────────────────────────────────────────────────
# BRANCHE
# ─────────────────────────────────────────────────────────────────────────────
def compute_branch_data(K, psi_roots, a_gam, alpha, r_max=40.0):
    branches = []
    for pz in psi_roots:
        E = energy_at_psi(pz, K, a_gam, alpha, r_max)
        if is_finite(E):
            g = E / (4 * np.pi * K) - 1.0
            EoK = E / K
            branches.append(
                {"K": K, "psi": pz, "E": E, "g": g, "EoK": EoK}
            )
    return branches

def select_B1(branches, prev_psi=None):
    b1_cands = [b for b in branches if B1_EOK_MIN <= b["EoK"] <= B1_EOK_MAX]
    if not b1_cands:
        return None

    if prev_psi is None or not is_finite(prev_psi):
        return min(b1_cands, key=lambda b: abs(b["g"]))
    return min(b1_cands, key=lambda b: abs(b["psi"] - prev_psi))

# ─────────────────────────────────────────────────────────────────────────────
# K*
# ─────────────────────────────────────────────────────────────────────────────
def find_Kstar_branches(
    alpha,
    a_gam=A_GAM,
    K_min=DEFAULT_K_MIN,
    K_max=DEFAULT_K_MAX,
    n_K=DEFAULT_NK,
    n_psi=DEFAULT_NPSI,
    r_max=DEFAULT_RMAX,
    psi_hint=None,
):
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)

    g_vals = []
    psi_guess_by_K = psi_hint

    for K in K_scan:
        roots = find_all_psi_zeros(
            K,
            a_gam,
            alpha,
            n_psi=n_psi,
            r_max=r_max,
            psi_hint=psi_guess_by_K,
        )

        if roots:
            branches = compute_branch_data(K, roots, a_gam, alpha, r_max=r_max)
            if branches:
                best_branch = min(branches, key=lambda b: abs(b["g"]))
                g_vals.append(best_branch["g"])
                psi_guess_by_K = best_branch["psi"]
            else:
                g_vals.append(np.nan)
        else:
            g_vals.append(np.nan)

    K_star = np.nan

    for i in range(len(g_vals) - 1):
        gi, gj = g_vals[i], g_vals[i + 1]
        if is_finite(gi) and is_finite(gj) and gi * gj < 0:
            try:
                local_psi_hint = psi_hint

                def g_interp(K):
                    nonlocal local_psi_hint
                    roots = find_all_psi_zeros(
                        K,
                        a_gam,
                        alpha,
                        n_psi=n_psi,
                        r_max=r_max,
                        psi_hint=local_psi_hint,
                    )
                    if not roots:
                        return 0.0

                    branches = compute_branch_data(K, roots, a_gam, alpha, r_max=r_max)
                    if not branches:
                        return 0.0

                    best_branch = min(branches, key=lambda b: abs(b["g"]))
                    local_psi_hint = best_branch["psi"]
                    return float(best_branch["g"])

                K_star = brentq(
                    g_interp,
                    K_scan[i],
                    K_scan[i + 1],
                    xtol=K_scan[i] * 0.001,
                    rtol=0.001,
                    maxiter=25,
                )
                break
            except Exception:
                pass

    if not is_finite(K_star):
        return []

    roots_star = find_all_psi_zeros(
        K_star,
        a_gam,
        alpha,
        n_psi=n_psi,
        r_max=r_max,
        psi_hint=psi_hint,
    )

    if not roots_star:
        return []

    return compute_branch_data(K_star, roots_star, a_gam, alpha, r_max=r_max)

# ─────────────────────────────────────────────────────────────────────────────
# WORKERY
# ─────────────────────────────────────────────────────────────────────────────
def solve_one_alpha(task):
    """
    Worker do pojedynczego alpha.
    """
    alpha = float(task["alpha"])
    psi_hint = task.get("psi_hint", None)
    a_gam = task.get("a_gam", A_GAM)
    K_min = task.get("K_min", DEFAULT_K_MIN)
    K_max = task.get("K_max", DEFAULT_K_MAX)
    n_K = task.get("n_K", DEFAULT_NK)
    n_psi = task.get("n_psi", DEFAULT_NPSI)
    r_max = task.get("r_max", DEFAULT_RMAX)

    t0 = time.perf_counter()
    branches = find_Kstar_branches(
        alpha=alpha,
        a_gam=a_gam,
        K_min=K_min,
        K_max=K_max,
        n_K=n_K,
        n_psi=n_psi,
        r_max=r_max,
        psi_hint=psi_hint,
    )
    b1 = select_B1(branches, prev_psi=psi_hint)
    dt = time.perf_counter() - t0

    return {
        "alpha": alpha,
        "b1": b1,
        "elapsed_s": dt,
    }

def parallel_alpha_scan(alpha_values, workers, psi_hint=None, a_gam=A_GAM,
                        K_min=DEFAULT_K_MIN, K_max=DEFAULT_K_MAX,
                        n_K=DEFAULT_NK, n_psi=DEFAULT_NPSI, r_max=DEFAULT_RMAX,
                        verbose=True):
    tasks = []
    for alpha in alpha_values:
        tasks.append(
            {
                "alpha": float(alpha),
                "psi_hint": psi_hint,
                "a_gam": a_gam,
                "K_min": K_min,
                "K_max": K_max,
                "n_K": n_K,
                "n_psi": n_psi,
                "r_max": r_max,
            }
        )

    results = {}
    times = {}

    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut_map = {ex.submit(solve_one_alpha, task): task["alpha"] for task in tasks}

        done = 0
        total = len(tasks)

        for fut in as_completed(fut_map):
            res = fut.result()
            alpha = res["alpha"]
            results[alpha] = res["b1"]
            times[alpha] = res["elapsed_s"]
            done += 1

            if verbose:
                b1 = res["b1"]
                if b1 is None:
                    print(f"[{done:>3}/{total}] alpha={alpha:.4f}  ->  B1: ---   t={res['elapsed_s']:.2f}s")
                else:
                    print(
                        f"[{done:>3}/{total}] alpha={alpha:.4f}  ->  "
                        f"K*={b1['K']:.6f}, psi={b1['psi']:.4f}, E={b1['E']:.8f}, "
                        f"E/K={b1['EoK']:.5f}, g={b1['g']:.5f}, t={res['elapsed_s']:.2f}s"
                    )

    results = dict(sorted(results.items()))
    times = dict(sorted(times.items()))
    return results, times

# ─────────────────────────────────────────────────────────────────────────────
# ANALIZA WYNIKÓW
# ─────────────────────────────────────────────────────────────────────────────
def analyze_minimum(b1_map):
    valid = [(a, b["E"], b["K"], b["psi"]) for a, b in b1_map.items() if b is not None and is_finite(b["E"])]
    out = {
        "valid": valid,
        "n_B1": len(valid),
        "alpha_min_dir": np.nan,
        "alpha_min_fit": np.nan,
        "best_argmin": np.nan,
        "odch": np.nan,
        "has_inner_min": False,
    }

    if len(valid) < 3:
        return out

    a_arr = np.array([x[0] for x in valid], dtype=float)
    E_arr = np.array([x[1] for x in valid], dtype=float)

    idx_min = np.argmin(E_arr)
    alpha_min_dir = a_arr[idx_min]
    out["alpha_min_dir"] = alpha_min_dir

    alpha_min_fit = np.nan
    if len(valid) >= 5:
        i0 = max(1, idx_min - 2)
        i1 = min(len(valid) - 2, idx_min + 2)
        a_fit = a_arr[i0:i1+1]
        E_fit = E_arr[i0:i1+1]

        if len(a_fit) >= 3:
            coeffs = np.polyfit(a_fit, E_fit, 2)
            if coeffs[0] > 0:
                alpha_min_fit = -coeffs[1] / (2 * coeffs[0])
    out["alpha_min_fit"] = alpha_min_fit

    best = alpha_min_fit if is_finite(alpha_min_fit) else alpha_min_dir
    out["best_argmin"] = best
    out["odch"] = abs(best - ALPHA_K_REF)
    out["has_inner_min"] = (idx_min > 0 and idx_min < len(valid) - 1)

    return out

def print_table(title, b1_map):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)
    print(f"\n  {'alpha':>7}  {'K*':>10}  {'psi_B1':>8}  {'E_B1':>14}  {'E1/K':>9}  {'g':>8}")
    print("  " + "-" * 62)

    for alpha, b1 in sorted(b1_map.items()):
        if b1 is None:
            print(f"  {alpha:>7.4f}  {'---':>10}  {'---':>8}  {'---':>14}  {'---':>9}  {'---':>8}")
        else:
            print(
                f"  {alpha:>7.4f}  {b1['K']:>10.6f}  {b1['psi']:>8.4f}  "
                f"{b1['E']:>14.8f}  {b1['EoK']:>9.5f}  {b1['g']:>8.5f}"
            )

def print_summary(analysis, label="COARSE"):
    print("\n" + "=" * 72)
    print(f"PODSUMOWANIE {label}")
    print("=" * 72)
    print(f"\n  Liczba punktow B1: {analysis['n_B1']}")

    if is_finite(analysis["alpha_min_dir"]):
        print(f"  Minimum bezposrednie: alpha = {analysis['alpha_min_dir']:.5f}")
    if is_finite(analysis["alpha_min_fit"]):
        print(f"  Minimum z fitu:      alpha = {analysis['alpha_min_fit']:.5f}")
    if is_finite(analysis["best_argmin"]):
        print(f"  Najlepszy argmin:    alpha = {analysis['best_argmin']:.5f}")
        print(f"  alpha_K_ref:         {ALPHA_K_REF:.5f}")
        print(f"  Odchylenie:          {analysis['odch']:.5f}")
        print(f"  Minimum wewnatrz:    {analysis['has_inner_min']}")

def print_runtime_stats(times, label="RUNTIME"):
    vals = np.array(list(times.values()), dtype=float) if times else np.array([])
    print("\n" + "=" * 72)
    print(label)
    print("=" * 72)
    if len(vals) == 0:
        print("  Brak danych.")
        return

    print(f"  Punkty:   {len(vals)}")
    print(f"  Min:      {vals.min():.2f} s")
    print(f"  Mean:     {vals.mean():.2f} s")
    print(f"  Median:   {np.median(vals):.2f} s")
    print(f"  Max:      {vals.max():.2f} s")
    print(f"  Sum CPU:  {vals.sum():.2f} s")

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def build_parser():
    p = argparse.ArgumentParser(description="Parallel B1 scan for TGP.")
    p.add_argument("--workers", type=int, default=None, help="Liczba workerow.")
    p.add_argument("--alpha-min", type=float, default=DEFAULT_ALPHA_MIN)
    p.add_argument("--alpha-max", type=float, default=DEFAULT_ALPHA_MAX)
    p.add_argument("--alpha-step", type=float, default=DEFAULT_ALPHA_STEP)
    p.add_argument("--fine", action="store_true", help="Uruchom refine wokol minimum.")
    p.add_argument("--fine-window", type=float, default=DEFAULT_FINE_WINDOW)
    p.add_argument("--fine-step", type=float, default=DEFAULT_FINE_STEP)
    p.add_argument("--nK", type=int, default=DEFAULT_NK)
    p.add_argument("--nPsi", type=int, default=DEFAULT_NPSI)
    p.add_argument("--rmax", type=float, default=DEFAULT_RMAX)
    p.add_argument("--quiet", action="store_true")
    return p

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    args = build_parser().parse_args()

    cpu_total = os.cpu_count() or 4
    workers = args.workers if args.workers is not None else min(16, cpu_total)

    alpha_vals = np.arange(args.alpha_min, args.alpha_max + 0.5 * args.alpha_step, args.alpha_step)

    print("=" * 72)
    print("TGP v1  *  p84_branch1_scan_parallel.py")
    print("=" * 72)
    print(f"\n  CPU logical:    {cpu_total}")
    print(f"  Workers:        {workers}")
    print(f"  alpha range:    [{args.alpha_min}, {args.alpha_max}] step {args.alpha_step}")
    print(f"  a_Gam:          {A_GAM}")
    print(f"  K range:        [{DEFAULT_K_MIN}, {DEFAULT_K_MAX}]")
    print(f"  nK:             {args.nK}")
    print(f"  nPsi:           {args.nPsi}")
    print(f"  rmax:           {args.rmax}")
    print(f"  Fine refine:    {args.fine}")
    print()

    t0 = time.perf_counter()

    # COARSE
    b1_coarse, times_coarse = parallel_alpha_scan(
        alpha_values=alpha_vals,
        workers=workers,
        psi_hint=None,
        a_gam=A_GAM,
        K_min=DEFAULT_K_MIN,
        K_max=DEFAULT_K_MAX,
        n_K=args.nK,
        n_psi=args.nPsi,
        r_max=args.rmax,
        verbose=not args.quiet,
    )

    coarse_analysis = analyze_minimum(b1_coarse)

    if args.quiet:
        print_table("SEKCJA A: Coarse scan alpha — galaz B1", b1_coarse)
    print_summary(coarse_analysis, label="COARSE")
    print_runtime_stats(times_coarse, label="RUNTIME COARSE")

    # FINE
    b1_fine = None
    fine_analysis = None
    times_fine = None

    if args.fine and is_finite(coarse_analysis["best_argmin"]):
        alpha_center = coarse_analysis["best_argmin"]
        a_lo = max(args.alpha_min, alpha_center - args.fine_window)
        a_hi = min(args.alpha_max, alpha_center + args.fine_window)
        alpha_fine = np.arange(a_lo, a_hi + 0.5 * args.fine_step, args.fine_step)

        print("\n" + "=" * 72)
        print(f"SEKCJA B: Fine scan wokol alpha={alpha_center:.5f}")
        print("=" * 72)

        b1_fine, times_fine = parallel_alpha_scan(
            alpha_values=alpha_fine,
            workers=workers,
            psi_hint=None,
            a_gam=A_GAM,
            K_min=DEFAULT_K_MIN,
            K_max=DEFAULT_K_MAX,
            n_K=args.nK,
            n_psi=args.nPsi,
            r_max=args.rmax,
            verbose=not args.quiet,
        )

        fine_analysis = analyze_minimum(b1_fine)

        if args.quiet:
            print_table("SEKCJA B: Fine scan alpha — galaz B1", b1_fine)
        print_summary(fine_analysis, label="FINE")
        print_runtime_stats(times_fine, label="RUNTIME FINE")

    t1 = time.perf_counter()
    wall = t1 - t0

    print("\n" + "=" * 72)
    print("WNIOSKI")
    print("=" * 72)

    best_analysis = fine_analysis if fine_analysis is not None else coarse_analysis

    if is_finite(best_analysis["best_argmin"]):
        odch = best_analysis["odch"]
        best = best_analysis["best_argmin"]

        print(f"\n  argmin E*(B1)   = {best:.5f}")
        print(f"  alpha_K_ref     = {ALPHA_K_REF:.5f}")
        print(f"  odchylenie      = {odch:.5f}")

        if odch < 0.02:
            print("  WNIOSEK: mocne trafienie (< 0.02)")
        elif odch < 0.10:
            print("  WNIOSEK: trafienie umiarkowane (< 0.10)")
        else:
            print("  WNIOSEK: brak potwierdzenia bliskosci do alpha_K")
    else:
        print("\n  Nie znaleziono wiarygodnego minimum B1.")

    print("\n" + "=" * 72)
    print("CZAS CALKOWITY")
    print("=" * 72)
    print(f"  Wall time: {wall:.2f} s")

if __name__ == "__main__":
    main()