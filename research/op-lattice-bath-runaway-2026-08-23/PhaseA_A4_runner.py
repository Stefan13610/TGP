#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
PhaseA A4 helper: uruchamia 8 skryptow rdzenia maszynerii 2 (obiekt audytu A4)
i przechwytuje ich FAKTYCZNE outputy do plikow PhaseA_A4_output_<nazwa>.txt.

Cel (LOCK A4 ii): porownanie SUMMARY skryptu z faktycznym outputem.
Wyniki skryptow NIE sa traktowane jako prawda (LOCK) -- sluza wylacznie
do audytu deskryptywnego (czy skrypt ma osiagalny FAIL, czy SUMMARY
zgadza sie z liczbami, rejestr INPUT/DERIVED).

Timeout per skrypt: 2400 s. Rownolegle: 4 procesy.
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.abspath(os.path.join(
    HERE, "..", "..", "tooling", "scripts"))

SCRIPTS = [
    "ngen_collapse_proof_v47b.py",
    "gcrit_energy_proof_v47b.py",
    "gcrit_pohozaev_v47b.py",
    "atail_asymptotic_v47b.py",
    "atail_functional_v47b.py",
    "ode_koide_formA_exact_v47b.py",
    "collapse_exponent_v47b.py",
    "a3d_soliton_brannen_r.py",
]

TIMEOUT = 2400


def run_one(name):
    src = os.path.join(SCRIPTS_DIR, name)
    out_path = os.path.join(
        HERE, "PhaseA_A4_output_%s.txt" % name.replace(".py", ""))
    t0 = time.time()
    try:
        proc = subprocess.run(
            [sys.executable, src], cwd=SCRIPTS_DIR,
            capture_output=True, text=True, timeout=TIMEOUT)
        body = proc.stdout + ("\n[STDERR]\n" + proc.stderr
                              if proc.stderr else "")
        status = "rc=%d" % proc.returncode
    except subprocess.TimeoutExpired as exc:
        body = ((exc.stdout.decode("utf-8", "replace")
                 if isinstance(exc.stdout, bytes) else (exc.stdout or ""))
                + "\n[TIMEOUT po %d s]" % TIMEOUT)
        status = "TIMEOUT"
    dt = time.time() - t0
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("# %s  (%s, %.0f s)\n" % (name, status, dt))
        fh.write(body)
    return name, status, dt


if __name__ == "__main__":
    todo = sys.argv[1:] or SCRIPTS   # opcjonalnie: lista skryptow z CLI
    print("A4 runner: %d skryptow, timeout %d s" % (len(todo), TIMEOUT))
    with ThreadPoolExecutor(max_workers=4) as ex:
        for name, status, dt in ex.map(run_one, todo):
            print("  %-36s %-8s %.0f s" % (name, status, dt))
    print("A4 runner: DONE")
