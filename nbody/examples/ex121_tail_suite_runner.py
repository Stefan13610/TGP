#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex121_tail_suite_runner.py
============================
Kolejno uruchamia łańcuch skryptów **ogon / okna / linearyzacja** (ex114–ex120, ex122–ex124):

  ex114 → … → ex120 → ex122 → ex123 → ex124

Zapisuje log czasów i kodów wyjścia do `_outputs/ex121_tail_suite_log.txt`.
Kończy się kodem **1**, jeśli którykolwiek skrypt się wyłoży.

Opcje:
  --dry-run   tylko wypisz kolejność, bez uruchamiania

Przykład:
  python ex121_tail_suite_runner.py
"""

from __future__ import annotations

import argparse
import io
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "_outputs"

# Kolejność: ex115 przed ex120 (CSV); ex114 przed reszty (solver); reszta luźno
SUITE = [
    "ex114_tail_phase_map.py",
    "ex115_tail_window_scan.py",
    "ex116_tail_rL_from_epsilon.py",
    "ex117_linear_operator_residual.py",
    "ex118_rmse_vs_linear_residual.py",
    "ex119_tail_pipeline_digest.py",
    "ex120_ex115_scan_overlay.py",
    "ex122_cross_term_ratio.py",
    "ex123_koide_epistemics.py",
    "ex124_dense_g0_solver_compare.py",
]


def main() -> int:
    ap = argparse.ArgumentParser(description="Runner łańcucha tail ex114–ex120 + ex122–ex124")
    ap.add_argument("--dry-run", action="store_true", help="Nie uruchamiaj skryptów")
    args = ap.parse_args()

    print("=" * 70)
    print("EX121: TAIL SUITE RUNNER (ex114 … ex120, ex122–ex124)")
    print("=" * 70)
    for i, name in enumerate(SUITE, 1):
        print(f"  {i}. {name}")
    print()

    if args.dry_run:
        print("  (--dry-run) zakończono.")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = [
        f"# ex121_tail_suite_log  {datetime.now(timezone.utc).isoformat()}Z",
        f"# cwd={HERE}",
        "",
    ]
    failures: list[str] = []

    for script in SUITE:
        path = HERE / script
        if not path.is_file():
            msg = f"BRAK PLIKU: {path}"
            print(f"  [SKIP/FAIL] {msg}")
            log_lines.append(f"FAIL missing {script}")
            failures.append(script)
            continue

        print(f"  >>> {script} …", flush=True)
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(
                [sys.executable, str(path)],
                cwd=str(HERE),
                timeout=300,
            )
        except subprocess.TimeoutExpired:
            dt = time.perf_counter() - t0
            print(f"  [FAIL] timeout po {dt:.1f}s")
            log_lines.append(f"TIMEOUT {dt:.3f}s  {script}")
            failures.append(script)
            continue

        dt = time.perf_counter() - t0
        status = "OK" if proc.returncode == 0 else "FAIL"
        print(f"  [{status}] exit={proc.returncode}  czas={dt:.2f}s")
        log_lines.append(f"{status} exit={proc.returncode}  {dt:.3f}s  {script}")
        if proc.returncode != 0:
            failures.append(script)

    log_lines.append("")
    log_lines.append(f"# failures: {len(failures)}")
    for f in failures:
        log_lines.append(f"#   {f}")

    log_path = OUT_DIR / "ex121_tail_suite_log.txt"
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print()
    print(f"  Log: {log_path}")
    print("=" * 70)

    if failures:
        print(f"  EX121: ZAKOŃCZONO Z BŁĘDAMI ({len(failures)} skryptów)")
        return 1
    print("  EX121: WSZYSTKIE SKRYPTY OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
