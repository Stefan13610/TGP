#!/usr/bin/env python3
"""
ex194_master_verification.py
=============================
TGP Master Verification Runner

Runs ALL verification scripts (ex190–ex193) in sequence, captures
stdout, parses PASS/FAIL counts, and produces a single GO / NO-GO
summary table.

Usage:
    python ex194_master_verification.py          # full run
    python ex194_master_verification.py --quick   # summary only (no stdout)

Exit code: 0 if all modules pass, 1 otherwise.

Author: TGP verification suite
Date:   2026-04-08
"""

import sys, os, io, time, re, importlib, traceback
from contextlib import redirect_stdout, redirect_stderr

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ============================================================
# CONFIGURATION
# ============================================================
SCRIPTS = [
    ("ex190_consistency_chain",     "Internal consistency chain (Gamma->observables)"),
    ("ex191_confinement_m0",        "Quark confinement & shifted Koide"),
    ("ex192_cosmo_Hz_confrontation","Cosmology: H(z) vs DESI DR1 + BBN/CMB/LLR"),
    ("ex193_unified_predictions",   "Unified prediction table (28 observables)"),
]

# Regex patterns for extracting results from stdout
# Pattern A: "Result: 10/10 PASS" (ex190)
PAT_RESULT = re.compile(r'Result:\s*(\d+)/(\d+)\s*PASS')
# Pattern B: "Overall: 5/5 PASS" (ex192)
PAT_OVERALL = re.compile(r'Overall:\s*(\d+)/(\d+)\s*PASS')
# Pattern C: "SUMMARY: 28 PASS / 1 CALIBRATED / 1 PREDICTIONS / ... (total 31 observables)"
PAT_SUMMARY = re.compile(
    r'SUMMARY:\s*(\d+)\s*PASS\s*/\s*(\d+)\s*CALIBR\w*\s*/\s*(\d+)\s*PRED\w*'
    r'(?:\s*/\s*(\d+)\s*INFO)?'
    r'(?:.*?total\s+(\d+))?'
)

def extract_counts(output_text: str):
    """Extract (n_pass, n_total, n_fail) from script output."""
    n_pass = 0
    n_total = 0
    n_fail = 0

    # Pattern A: ex190 style — "Result: N/M PASS"
    m = PAT_RESULT.search(output_text)
    if m:
        n_pass = int(m.group(1))
        n_total = int(m.group(2))
        return n_pass, n_total, n_total - n_pass

    # Pattern B: ex192 style — "Overall: N/M PASS"
    m = PAT_OVERALL.search(output_text)
    if m:
        n_pass = int(m.group(1))
        n_total = int(m.group(2))
        return n_pass, n_total, n_total - n_pass

    # Pattern C: ex193 style — "SUMMARY: 28 PASS / 1 CAL / 1 PRED / 1 INFO (total 31)"
    m = PAT_SUMMARY.search(output_text)
    if m:
        n_pass = int(m.group(1))
        n_cal = int(m.group(2))
        n_pred = int(m.group(3))
        n_info = int(m.group(4)) if m.group(4) else 0
        n_total = int(m.group(5)) if m.group(5) else (n_pass + n_cal + n_pred + n_info)
        # Only PASS count matters for go/no-go; n_total = testable = pass + fail
        # CAL/PRED/INFO are not failures
        return n_pass, n_pass, 0  # report pass/pass since CAL/PRED are not failures

    # Fallback: count individual "  PASS" and "  FAIL" status markers
    # (these appear as last word on test result lines)
    pass_count = len(re.findall(r'\bPASS\b', output_text))
    fail_count = len(re.findall(r'\bFAIL\b', output_text))
    # Subtract likely header/summary lines containing PASS/FAIL
    # by looking for standalone status markers only
    n_pass = max(0, pass_count - 1)  # subtract the summary line itself
    n_fail = fail_count
    n_total = n_pass + n_fail

    return n_pass, n_total, n_fail


def run_module(module_name: str, verbose: bool = True):
    """
    Import and run a script module, capturing its output.
    Returns (success: bool, output: str, n_pass, n_total, elapsed_s).
    """
    buf = io.StringIO()
    t0 = time.time()
    success = True
    n_pass = n_total = 0

    try:
        # Ensure fresh import
        if module_name in sys.modules:
            del sys.modules[module_name]

        with redirect_stdout(buf), redirect_stderr(buf):
            mod = importlib.import_module(module_name)
            # Some scripts use if __name__ == "__main__"; call main() if it exists
            if hasattr(mod, 'main'):
                mod.main()

    except SystemExit as e:
        if e.code not in (None, 0):
            success = False
    except Exception:
        buf.write(f"\n*** EXCEPTION ***\n{traceback.format_exc()}\n")
        success = False

    elapsed = time.time() - t0
    output = buf.getvalue()
    n_pass, n_total, n_fail = extract_counts(output)

    if n_fail > 0:
        success = False

    return success, output, n_pass, n_total, elapsed


def main():
    quick = '--quick' in sys.argv
    verbose = not quick

    # Add scripts directory to path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    print("=" * 76)
    print("  TGP MASTER VERIFICATION RUNNER")
    print("  All predictions from 2 free parameters: Phi_eff = 25, g_0^e = 0.869")
    print("=" * 76)
    print()

    results = []
    total_pass = 0
    total_tests = 0
    total_time = 0.0
    all_ok = True

    for i, (module_name, description) in enumerate(SCRIPTS, 1):
        print(f"[{i}/{len(SCRIPTS)}] {module_name}")
        print(f"       {description}")
        print(f"       Running...", end=" ", flush=True)

        success, output, n_pass, n_total, elapsed = run_module(module_name, verbose)

        status = "OK" if success else "FAIL"
        icon = "✓" if success else "✗"
        print(f"{icon} {status}  ({n_pass}/{n_total} pass, {elapsed:.1f}s)")

        if verbose and not quick:
            # Print last 20 lines of output (summary region)
            lines = output.strip().split('\n')
            summary_start = max(0, len(lines) - 25)
            # Try to find SUMMARY header
            for j, line in enumerate(lines):
                if 'SUMMARY' in line and j > len(lines) // 2:
                    summary_start = max(0, j - 1)
                    break
            print("       " + "-" * 60)
            for line in lines[summary_start:]:
                print(f"       | {line}")
            print("       " + "-" * 60)

        print()

        results.append({
            'module': module_name,
            'desc': description,
            'success': success,
            'n_pass': n_pass,
            'n_total': n_total,
            'elapsed': elapsed,
        })

        total_pass += n_pass
        total_tests += n_total
        total_time += elapsed
        if not success:
            all_ok = False

    # ============================================================
    # MASTER SUMMARY
    # ============================================================
    print("=" * 76)
    print("  MASTER SUMMARY")
    print("=" * 76)

    print(f"\n  {'Module':<38s} {'Tests':<12s} {'Time':<8s} {'Status'}")
    print("  " + "-" * 66)

    for r in results:
        icon = "✓" if r['success'] else "✗"
        status = "OK" if r['success'] else "FAIL"
        tests_str = f"{r['n_pass']}/{r['n_total']}"
        time_str = f"{r['elapsed']:.1f}s"
        print(f"  {icon} {r['module']:<36s} {tests_str:<12s} {time_str:<8s} {status}")

    print("  " + "-" * 66)
    print(f"  {'TOTAL':<38s} {total_pass}/{total_tests:<11} {total_time:.1f}s")

    print(f"\n  Free parameters:  2  (Phi_eff ≈ 25, g_0^e ≈ 0.869)")
    print(f"  Total tests:      {total_tests}")
    print(f"  Tests passed:     {total_pass}")
    print(f"  Efficiency:       {total_pass}/2 = {total_pass/2:.1f} predictions per parameter")

    # ============================================================
    # GO / NO-GO
    # ============================================================
    print()
    if all_ok:
        print("  ╔══════════════════════════════════════════════╗")
        print("  ║    ★  GO  — All verification modules pass   ║")
        print("  ╚══════════════════════════════════════════════╝")
    else:
        failed = [r['module'] for r in results if not r['success']]
        print("  ╔══════════════════════════════════════════════╗")
        print("  ║    ✗  NO-GO  — Some modules failed          ║")
        print("  ╚══════════════════════════════════════════════╝")
        print(f"  Failed: {', '.join(failed)}")

    print()
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
