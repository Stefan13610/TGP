# -*- coding: utf-8 -*-
"""
mass_test_runner.py -- TGP batch test runner
Uruchamia wszystkie skrypty .py z podfolderow i raportuje PASS/FAIL.
"""
import subprocess
import sys
import os

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))

SKIP_PATTERNS = [
    'monte_carlo', 'kink_mc', 'vw_substrate_3d',
    'cmb_tgp', 'cmb_tensor', 'big_bang', 'bbn_timeline',
    'cosmological_phi_evolution', 'deep_consistency',
    'master_verification', 'mass_test_runner',
    'check_refs', 'fix_latex_typos',  # wymagają uruchomienia z TGP_v1/, nie scripts/
]

scripts = []
for root, dirs, files in os.walk(SCRIPTS_DIR):
    dirs[:] = [d for d in dirs if d not in ('archive', '__pycache__')]
    for f in sorted(files):
        if f.endswith('.py') and f != '__init__.py':
            scripts.append(os.path.join(root, f))

passed = []
failed = []
timeouts = []
skipped = []

for script in sorted(scripts):
    rel = os.path.relpath(script, SCRIPTS_DIR).replace('\\', '/')
    skip = any(p in rel for p in SKIP_PATTERNS)
    if skip:
        skipped.append(rel)
        continue
    try:
        r = subprocess.run(
            [sys.executable, script],
            capture_output=True, text=True,
            encoding='utf-8', errors='replace',
            timeout=45,
            cwd=os.path.dirname(script)
        )
        if r.returncode == 0:
            passed.append(rel)
            print(f"  [PASS] {rel}")
        else:
            out = (r.stdout or '') + (r.stderr or '')
            lines = [l for l in out.strip().split('\n') if l.strip()]
            snippet = '\n         '.join(lines[-3:]) if lines else '?'
            failed.append((rel, snippet))
            print(f"  [FAIL] {rel}")
            print(f"         {snippet}")
    except subprocess.TimeoutExpired:
        timeouts.append(rel)
        print(f"  [TOUT] {rel}  (>45s)")

print()
print("=" * 65)
print(f"  PASS:    {len(passed)}")
print(f"  FAIL:    {len(failed)}")
print(f"  TIMEOUT: {len(timeouts)}")
print(f"  SKIP:    {len(skipped)}")
print("=" * 65)

if failed:
    print("\n  FAILURES DETAIL:")
    for name, err in failed:
        print(f"\n  [{name}]")
        for line in err.split('\n'):
            print(f"    {line}")

if timeouts:
    print("\n  TIMEOUTS:")
    for name in timeouts:
        print(f"    {name}")
