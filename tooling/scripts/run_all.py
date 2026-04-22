#!/usr/bin/env python3
"""
run_all.py — TGP Reproducibility Runner
=========================================

Runs all key verification scripts for TGP (Teoria Generowanej Przestrzeni).
Each script is a self-contained numerical test; PASS/FAIL is determined by
comparing TGP predictions against PDG/Planck/DESI data.

Usage:
    python run_all.py              # run all core scripts
    python run_all.py --quick      # run only quick-check scripts (~30s)
    python run_all.py --full       # run all scripts including slow ones

Requirements:
    Python 3.9+
    numpy, scipy (pip install numpy scipy)
"""

import subprocess
import sys
import time
import os
from pathlib import Path

# ============================================================
# Script registry: (filename, category, estimated_time_s, description)
# ============================================================
SCRIPTS = [
    # --- CORE: particle sector (the strongest arguments) ---
    ("ex147_substrate_ode_phiFP.py",  "core",  120,
     "Substrate ODE solver: soliton g(r), phi-FP mechanism"),
    ("ex157_one_param_prediction.py", "core",  120,
     "One-parameter prediction: g0^e -> r21, m_mu, m_tau"),
    ("ex174_tgp_predictions_summary.py", "core", 30,
     "Full predictions summary: 11 observables from ODE"),
    ("ex184_discrete_running.py",     "core",  3,
     "Discrete running: alpha_s(N_f) at flavor thresholds"),
    ("ex186_mass_coupling_unification.py", "core", 180,
     "Mass-coupling unification: r21 -> g0^e -> alpha_s (0 free params)"),

    # --- CORE: cosmology ---
    ("ex187_cosmo_data_confrontation.py", "core", 5,
     "Cosmological data confrontation: H(z), w(z), G(z) vs Planck/DESI"),
    ("ex165_slowroll_tgp.py",         "core",  3,
     "Slow-roll inflation: n_s, r predictions"),

    # --- SUPPORTING: soliton & mass hierarchy ---
    ("ex147d_actual_ode.py",          "support", 5,
     "Correct ODE form verification"),
    ("ex150_potential_degree_test.py", "support", 3,
     "Potential degree constraint"),
    ("ex153_Btail_zeros_koide.py",    "support", 5,
     "Tail zeros and Koide relation"),
    ("ex162_phi0_from_substrate.py",  "support", 3,
     "Phi_0 from substrate self-consistency"),

    # --- SUPPORTING: gauge sector ---
    ("ex140_z3_entropy_optimum.py",   "support", 3,
     "Z_3 entropy optimum"),
    ("ex145_B2_qk_z3_bridge.py",     "support", 5,
     "Bridge B2: Q_K from Z_3 symmetry"),
    ("ex178_alphas_new_formula.py",   "support", 3,
     "alpha_s formula derivation and test"),
    ("ex183_phi0_origin.py",          "support", 3,
     "Phi_0 = N_f^2 hypothesis"),

    # --- SUPPORTING: ERG & stability ---
    ("ex141_erg_full_K_flow.py",      "support", 10,
     "ERG flow: full K(phi) running"),
    ("ex142_A4_perturbative.py",      "support", 5,
     "A4 perturbative computation"),
    ("ex146_E3_cancellation_proof.py", "support", 3,
     "E3 cancellation proof"),
    ("ex148_E3_analytical_proof.py",  "support", 3,
     "E3 analytical verification"),

    # --- SUPPORTING: Koide sector ---
    ("ex168_cross_sector_koide.py",   "support", 3,
     "Cross-sector Koide relations"),
    ("ex155_tau_selection.py",        "support", 3,
     "Tau mass selection mechanism"),

    # --- SUPPORTING: PPN & metric ---
    ("ex167_ppn_alpha1.py",           "support", 3,
     "PPN parameters for alpha=1"),
    ("ex166_alpha1_vs_alpha2.py",     "support", 3,
     "alpha=1 vs alpha=2 comparison"),

    # --- EXPLORATORY: running analysis ---
    ("ex185_running_analysis.py",     "explore", 10,
     "Multi-loop QCD running comparison (partial)"),
    ("ex175_alphas_alpha1.py",        "explore", 3,
     "alpha_s at alpha=1"),
    ("ex180_phi0_exact.py",           "explore", 3,
     "Phi_0 exact value exploration"),
]

# Quick-check subset (runs in ~30s)
QUICK_SCRIPTS = [s[0] for s in SCRIPTS if s[1] == "core"]


def run_script(script_name, script_dir):
    """Run a single script and return (success, duration, output)."""
    path = script_dir / script_name
    if not path.exists():
        return False, 0, f"FILE NOT FOUND: {path}"

    t0 = time.time()
    try:
        result = subprocess.run(
            [sys.executable, "-X", "utf8", str(path)],
            capture_output=True, text=True,
            timeout=300, cwd=str(script_dir),
            encoding="utf-8", errors="replace"
        )
        duration = time.time() - t0
        output = result.stdout + result.stderr
        success = result.returncode == 0
        return success, duration, output
    except subprocess.TimeoutExpired:
        return False, time.time() - t0, "TIMEOUT (>120s)"
    except Exception as e:
        return False, time.time() - t0, f"ERROR: {e}"


def main():
    script_dir = Path(__file__).parent

    # Parse args
    mode = "core"  # default: core scripts only
    if "--quick" in sys.argv:
        mode = "quick"
    elif "--full" in sys.argv:
        mode = "full"
    elif "--all" in sys.argv:
        mode = "full"

    if mode == "quick":
        scripts_to_run = [(s, cat, t, desc) for s, cat, t, desc in SCRIPTS
                          if s in QUICK_SCRIPTS]
    elif mode == "full":
        scripts_to_run = SCRIPTS
    else:
        scripts_to_run = [(s, cat, t, desc) for s, cat, t, desc in SCRIPTS
                          if cat in ("core",)]

    print("=" * 70)
    print(f"TGP Reproducibility Runner  |  Mode: {mode}")
    print(f"Scripts to run: {len(scripts_to_run)}")
    print("=" * 70)

    results = []
    total_t0 = time.time()

    for script_name, category, est_time, description in scripts_to_run:
        print(f"\n--- [{category.upper():>7s}] {script_name} ---")
        print(f"    {description}")

        success, duration, output = run_script(script_name, script_dir)
        status = "PASS" if success else "FAIL"
        results.append((script_name, category, status, duration))

        print(f"    Status: {status}  ({duration:.1f}s)")
        if not success:
            # Print last 5 lines of output for debugging
            lines = output.strip().split('\n')
            for line in lines[-5:]:
                print(f"    > {line}")

    total_time = time.time() - total_t0

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    n_pass = sum(1 for _, _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, _, s, _ in results if s == "FAIL")

    print(f"\n{'Script':<45s} {'Cat':>7s} {'Status':>6s} {'Time':>6s}")
    print("-" * 68)
    for name, cat, status, dur in results:
        marker = "+" if status == "PASS" else "X"
        print(f"  [{marker}] {name:<42s} {cat:>7s} {status:>6s} {dur:5.1f}s")

    print("-" * 68)
    print(f"Total: {n_pass} PASS, {n_fail} FAIL  |  "
          f"Time: {total_time:.0f}s")

    if n_fail == 0:
        print("\nAll scripts passed successfully.")
    else:
        print(f"\n{n_fail} script(s) failed. Check output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
