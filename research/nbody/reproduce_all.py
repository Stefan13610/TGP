#!/usr/bin/env python3
"""
reproduce_all.py -- One-command reproducibility check for TGP N-body paper
==========================================================================

Runs all regression suites and key result scripts, reports status.
Exit code 0 = all passed, 1 = failures detected.

Usage:
    python -m nbody.reproduce_all          # full run
    python -m nbody.reproduce_all --quick  # fast (~5 min)
"""

import subprocess
import sys
import time
import os

REPO = os.path.dirname(os.path.abspath(__file__))
REPO_PARENT = os.path.dirname(REPO)

quick = "--quick" in sys.argv
quick_flag = ["--quick"] if quick else []

# Scripts to run in order: (description, module path)
STAGES = [
    # Stage 1: Regression suites
    ("Regression: verify_all (59 tests)", "nbody.examples.verify_all"),
    ("Regression: EOM quick", "nbody.examples.verify_nbody_eom_quick"),
    ("Regression: Lyapunov quick", "nbody.examples.verify_nbody_lyapunov_quick"),
    # Stage 2: Key result scripts
    ("Result 1: Earnshaw (ex201)", "nbody.examples.ex201_full_hessian_stability_scan"),
    ("Result 2: Chaos beta scan (ex200)", "nbody.examples.ex200_lyapunov_beta_scan_yukawa"),
    ("Result 2: Multi-IC (ex207)", "nbody.examples.ex207_lyapunov_multi_ic_beta_scan_v3"),
    ("Result 5: Multipole (ex202)", "nbody.examples.ex202_multipole_vs_feynman_triple_overlap"),
    ("Result 6: V3/V2 regime (ex209)", "nbody.examples.ex209_v3_v2_regime_msp_scan"),
    ("Result 7: Equilibria+Hill (ex210)", "nbody.examples.ex210_analytical_equilibria_hill"),
    ("Result 8: Phase diagram (ex211)", "nbody.examples.ex211_unequal_mass_phase_diagram"),
]

def run_stage(desc, module):
    """Run a module and return (success, elapsed)."""
    cmd = [sys.executable, "-m", module] + quick_flag
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, cwd=REPO_PARENT,
            capture_output=True, text=True, timeout=600,
            encoding="utf-8", errors="replace"
        )
        elapsed = time.time() - t0
        ok = result.returncode == 0
        # For scripts that don't use exit codes, check for FAIL/Error in output
        if ok and ("FAIL" in result.stdout or "Traceback" in result.stdout):
            # Some scripts print FAIL but exit 0 — check context
            # Only flag if FAIL is standalone (not in "FAILED: 0")
            lines = result.stdout.split("\n")
            real_fail = any(
                "FAIL" in l and "FAILED:  0" not in l and "FAILED, 0" not in l
                and "0 FAILED" not in l
                for l in lines
            )
            if real_fail:
                ok = False
        detail = ""
        if not ok:
            detail = (result.stderr or "")[-500:] + "\n" + (result.stdout or "")[-500:]
        return ok, elapsed, detail
    except subprocess.TimeoutExpired:
        return False, 600.0, "TIMEOUT"
    except Exception as e:
        return False, time.time() - t0, str(e)

def main():
    print("=" * 70)
    print("TGP N-body: Reproducibility Check")
    print(f"Mode: {'--quick' if quick else 'full'}")
    print("=" * 70)

    results = []
    total_t0 = time.time()

    for desc, module in STAGES:
        print(f"\n--- {desc} ---")
        ok, elapsed, err = run_stage(desc, module)
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {elapsed:.1f}s")
        if not ok and err:
            # Print last lines of error (stderr + stdout)
            for line in err.strip().split("\n")[-15:]:
                print(f"    {line}")
        results.append((desc, ok, elapsed))

    total_time = time.time() - total_t0

    # Summary
    n_pass = sum(1 for _, ok, _ in results if ok)
    n_fail = sum(1 for _, ok, _ in results if not ok)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for desc, ok, elapsed in results:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {desc} ({elapsed:.1f}s)")

    print(f"\nTotal: {n_pass} PASS, {n_fail} FAIL ({total_time:.0f}s)")

    if n_fail == 0:
        print("\nAll reproducibility checks passed.")
    else:
        print(f"\n{n_fail} check(s) FAILED. See output above.")

    return 0 if n_fail == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
