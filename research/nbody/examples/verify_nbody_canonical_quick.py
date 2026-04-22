#!/usr/bin/env python3
"""
verify_nbody_canonical_quick.py
===============================

Quick runner for the CURRENT canonical `nbody` path:

  1. FULL vs LPA synchronization (`ex195`, `ex196`, `ex197`)
  2. Defect -> C_eff -> Yukawa bridge (`ex205`)
  3. Effective N-body EOM regression (`verify_nbody_eom_quick.py`)

This script is meant as the shortest reproducible route through the modern
`nbody` layer after the theory-sync cleanup.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


EXAMPLES_DIR = Path(__file__).resolve().parent

QUICK_CHAIN = [
    ("ex195_soliton_ksub_full_vs_lpa.py", ["--quick"]),
    ("ex196_phi_fp_full_vs_lpa.py", ["--quick"]),
    ("ex197_optimal_g0_full_form.py", ["--quick"]),
    ("ex205_path_c_yukawa_from_defect.py", ["--quick"]),
    ("verify_nbody_eom_quick.py", []),
]


def main() -> None:
    failures: list[str] = []
    for script_name, extra_args in QUICK_CHAIN:
        script_path = EXAMPLES_DIR / script_name
        cmd = [sys.executable, str(script_path), *extra_args]
        print(f"==> {script_name}")
        completed = subprocess.run(cmd, cwd=str(EXAMPLES_DIR))
        if completed.returncode != 0:
            failures.append(script_name)

    if failures:
        print("FAILED:", ", ".join(failures))
        raise SystemExit(1)

    print("verify_nbody_canonical_quick: all PASS")


if __name__ == "__main__":
    main()
