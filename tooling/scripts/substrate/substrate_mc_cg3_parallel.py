"""
substrate_mc_cg3_parallel.py -- TGP Monte Carlo with parallel temperatures
===========================================================================
Wrapper around substrate_mc_cg3.py that runs multiple temperatures
in parallel using multiprocessing.

Usage:
    python scripts/substrate/substrate_mc_cg3_parallel.py

Uses 3 processes (one per temperature). Light CPU load (~10%).
Estimated runtime: ~50 min for L=32.
"""

import multiprocessing as mp
import numpy as np
import sys
import os
import time
import io

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from substrate_mc_cg3 import run_cg3_simulation


def run_single_temperature(args):
    """Worker: run MC for one temperature ratio."""
    T_ratio, L, b_list, n_warmup, n_measure, n_skip, seed = args

    # Redirect stdout per-process to avoid interleaving
    result = run_cg3_simulation(
        L=L, T_over_Tc=T_ratio, b_list=b_list,
        n_warmup=n_warmup, n_measure=n_measure,
        n_skip=n_skip, seed=seed + int(T_ratio * 1000)
    )
    return T_ratio, result


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

    print("=" * 72)
    print("  TGP -- Monte Carlo CG-3 (PARALLEL)")
    print("  L=32, 3 temperatures, multiprocessing")
    print("=" * 72)

    # Parameters
    L = 32
    n_warmup = 3000
    n_measure = 5000
    n_skip = 5
    seed = 42
    b_list = [2, 4, 8, 16]  # L=32 allows b up to 16
    T_ratios = [0.70, 0.80, 0.90]

    # Prepare args
    args_list = [
        (T_ratio, L, b_list, n_warmup, n_measure, n_skip, seed)
        for T_ratio in T_ratios
    ]

    print(f"\n  L={L}, blocks={b_list}")
    print(f"  Sweeps: {n_warmup} warmup + {n_measure} measure")
    print(f"  Temperatures: T/T_c = {T_ratios}")
    print(f"  Processes: {len(T_ratios)}")
    print(f"\n  Starting parallel MC...\n")

    t0 = time.time()

    # Run in parallel (3 processes)
    with mp.Pool(processes=len(T_ratios)) as pool:
        results_list = pool.map(run_single_temperature, args_list)

    dt_total = time.time() - t0
    print(f"\n  Total wall time: {dt_total:.0f}s ({dt_total/60:.1f} min)")

    # Collect results
    all_results = {}
    for T_ratio, res in results_list:
        all_results[T_ratio] = res

    # Print per-temperature summary
    for T_ratio in T_ratios:
        res = all_results[T_ratio]
        print(f"\n{'-' * 72}")
        print(f"  T/T_c = {T_ratio:.2f}")
        print(f"  <phi^2> = {res['v_sq']:.4f} +/- {res['v_sq_err']:.4f}")
        print(f"  xi_corr = {res['xi_avg']:.2f} +/- {res['xi_err']:.2f}")

        for b in b_list:
            if b not in res['blocks']:
                continue
            br = res['blocks'][b]
            print(f"\n  b={b} (L_B={br['L_B']}, Lb={br['Lb']}):")
            print(f"    <Phi_B> = {br['Phi_mean']:.4f} +/- {br['Phi_std']:.4f}")
            if br['c_star'] is not None:
                print(f"    c* = {br['c_star']:.6f}")
            else:
                print(f"    c* = N/A")
            if br['K1_Phi_product'] is not None:
                print(f"    K1*Phi = {br['K1_Phi_product']:.6f} (var={br['K1_Phi_relvar']:.2f})")
            else:
                print(f"    K1*Phi = N/A")
            print(f"    alpha=2 test: {'PASS' if br['alpha2_pass'] else 'FAIL'}")
            print(f"    separation: {'OK' if br['separation_ok'] else 'INSUFFICIENT'}")

            vr = br['veff']
            if 'error' not in vr:
                print(f"    Phi_0 = {vr['Phi_0']:.4f}")
                print(f"    m_sp^2 = {vr['m_sp_sq']:.4f}")
                if vr.get('beta_eff') and vr.get('gamma_eff'):
                    ratio = vr['beta_eff'] / vr['gamma_eff'] if vr['gamma_eff'] != 0 else float('inf')
                    print(f"    beta_eff = {vr['beta_eff']:.4f}, gamma_eff = {vr['gamma_eff']:.4f}")
                    print(f"    beta/gamma = {ratio:.3f}")

    # ================================================================
    # AGGREGATE TESTS
    # ================================================================
    print(f"\n{'=' * 72}")
    print("  AGGREGATE TESTS (L=32)")
    print(f"{'=' * 72}")

    # T1: xi scaling
    xi_vals = [(r, all_results[r]['xi_avg']) for r in T_ratios
               if all_results[r]['xi_avg'] > 0]
    if len(xi_vals) >= 2:
        log_eps = [np.log(1 - r) for r, _ in xi_vals]
        log_xi = [np.log(xi) for _, xi in xi_vals]
        coeffs = np.polyfit(log_eps, log_xi, 1)
        nu_eff = -coeffs[0]
        T1_pass = 0.3 < nu_eff < 1.0
        print(f"\n  T1: nu_eff = {nu_eff:.3f} (expected: 0.63)")
        print(f"      {'[PASS]' if T1_pass else '[FAIL]'} nu in [0.3, 1.0]")
    else:
        T1_pass = False
        print("\n  T1: insufficient data")

    # T2: Scale separation
    best_sep = False
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                if all_results[T_ratio]['blocks'][b]['separation_ok']:
                    best_sep = True
                    break
    T2_pass = best_sep
    print(f"\n  T2: Scale separation: {'[PASS]' if T2_pass else '[FAIL]'}")

    # T3: c* > 0 (CRITICAL)
    all_cstars = []
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                cs = all_results[T_ratio]['blocks'][b]['c_star']
                if cs is not None:
                    all_cstars.append(cs)
    T3_pass = len(all_cstars) > 0 and all(c > 0 for c in all_cstars)
    c_star_min = min(all_cstars) if all_cstars else None
    print(f"\n  T3: c* > 0 (Lemma A1): {'[PASS]' if T3_pass else '[FAIL]'}")
    if c_star_min is not None:
        print(f"      c*_min = {c_star_min:.6f}")

    # T4: K1*Phi ~ const
    all_alpha2 = []
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                all_alpha2.append(all_results[T_ratio]['blocks'][b]['alpha2_pass'])
    T4_pass = any(all_alpha2)
    print(f"\n  T4: K1*Phi ~ const (alpha=2): {'[PASS]' if T4_pass else '[FAIL]'}")

    # T5: Phi_0 > 0
    all_phi0 = []
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                vr = all_results[T_ratio]['blocks'][b]['veff']
                if 'Phi_0' in vr:
                    all_phi0.append(vr['Phi_0'])
    T5_pass = len(all_phi0) > 0 and all(p > 0 for p in all_phi0)
    print(f"\n  T5: Phi_0 > 0: {'[PASS]' if T5_pass else '[FAIL]'}")

    # T6: m_sp^2 > 0
    all_msp = []
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                vr = all_results[T_ratio]['blocks'][b]['veff']
                if 'm_sp_sq' in vr:
                    all_msp.append(vr['m_sp_sq'])
    T6_pass = len(all_msp) > 0 and all(m > 0 for m in all_msp)
    print(f"\n  T6: m_sp^2 > 0: {'[PASS]' if T6_pass else '[FAIL]'}")

    # T7: beta/gamma
    bg_ratios = []
    for T_ratio in T_ratios:
        for b in b_list:
            if b in all_results[T_ratio]['blocks']:
                vr = all_results[T_ratio]['blocks'][b]['veff']
                be = vr.get('beta_eff')
                ga = vr.get('gamma_eff')
                if be and ga and ga != 0:
                    bg_ratios.append(abs(be / ga - 1))
    T7_pass = len(bg_ratios) > 0 and min(bg_ratios) < 0.3
    print(f"\n  T7: |beta/gamma - 1| < 0.3: {'[PASS]' if T7_pass else '[FAIL]'}")
    if bg_ratios:
        print(f"      min |beta/gamma - 1| = {min(bg_ratios):.3f}")

    # Summary
    tests = [T1_pass, T2_pass, T3_pass, T4_pass, T5_pass, T6_pass, T7_pass]
    n_pass = sum(tests)
    n_total = len(tests)

    print(f"\n{'=' * 72}")
    print(f"  RESULT: {n_pass}/{n_total} PASS")
    print(f"{'=' * 72}")

    labels = ['T1 (xi scaling)', 'T2 (separation)', 'T3 (c*>0, A1)',
              'T4 (alpha=2, A3)', 'T5 (Phi_0>0)', 'T6 (m_sp^2>0)',
              'T7 (beta=gamma)']
    for label, passed in zip(labels, tests):
        status = '[PASS]' if passed else '[FAIL]'
        critical = ' <-- CRITICAL' if 'A1' in label and not passed else ''
        print(f"  {status} {label}{critical}")

    print(f"\n  Wall time: {dt_total:.0f}s ({dt_total/60:.1f} min)")
    print(f"\n{'=' * 72}")
    print("  DONE -- substrate_mc_cg3_parallel.py (L=32)")
    print(f"{'=' * 72}")

    return n_pass, n_total


if __name__ == "__main__":
    mp.freeze_support()  # Windows compatibility
    n_pass, n_total = main()
    sys.exit(0 if n_pass >= 5 else 1)
