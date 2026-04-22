#!/usr/bin/env python3
"""
verify_nbody_bridge_extended.py
===============================

Extended runner for the synchronized `nbody` bridge:

  1. FULL vs LPA synchronization (`ex195`, `ex196`, `ex197`)
  2. Classical defect -> C_eff -> Yukawa bridge (`ex205`)
  3. Effective N-body EOM regression (`verify_nbody_eom_quick.py`)
  4. Representative Lyapunov / chaos checks (`ex148`, `ex170`, `ex172`, `ex198`)

This is intentionally heavier than `verify_nbody_canonical_quick.py`. Use it
when you want one reproducible entry point spanning the current bridge from the
defect layer into EOM and chaos diagnostics.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


EXAMPLES_DIR = Path(__file__).resolve().parent

EXTENDED_CHAIN = [
    ("ex195_soliton_ksub_full_vs_lpa.py", ["--quick"]),
    ("ex196_phi_fp_full_vs_lpa.py", ["--quick"]),
    ("ex197_optimal_g0_full_form.py", ["--quick"]),
    ("ex205_path_c_yukawa_from_defect.py", ["--quick"]),
    ("verify_nbody_eom_quick.py", []),
    ("ex148_lyapunov_tgp_vs_newton_pairwise.py", ["--quick"]),
    ("ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py", ["--quick"]),
    ("ex172_lyapunov_pairwise_spectrum_sum_leapfrog.py", ["--quick"]),
    ("ex198_lyapunov_p1_closure.py", ["--quick"]),
]


def main() -> None:
    failures: list[str] = []
    for script_name, extra_args in EXTENDED_CHAIN:
        script_path = EXAMPLES_DIR / script_name
        cmd = [sys.executable, str(script_path), *extra_args]
        print(f"==> {script_name}")
        completed = subprocess.run(cmd, cwd=str(EXAMPLES_DIR))
        if completed.returncode != 0:
            failures.append(script_name)

    if failures:
        print("FAILED:", ", ".join(failures))
        raise SystemExit(1)

    print("verify_nbody_bridge_extended: all PASS")


if __name__ == "__main__":
    main()
