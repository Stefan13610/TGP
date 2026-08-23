#!/usr/bin/env python3
# -*- coding: ascii -*-
"""A4: uruchomienie a3d_soliton_brannen_r.py z przechwyceniem BAJTOWYM.

Powod: skrypt rdzenia przestawia stdout na UTF-8 i drukuje znaki
niedostepne w cp1250 -- przechwycenie tekstowe (PhaseA_A4_runner.py,
text=True) pada z UnicodeDecodeError w watku czytajacym. To samo
odnotowane w tabeli A4 jako wlasnosc skryptu (zalezy od kodowania
konsoli)."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(
    HERE, "..", "..", "tooling", "scripts", "a3d_soliton_brannen_r.py"))
OUT = os.path.join(HERE, "PhaseA_A4_output_a3d_soliton_brannen_r.txt")

t0 = time.time()
p = subprocess.run([sys.executable, SRC], cwd=os.path.dirname(SRC),
                   capture_output=True, timeout=2400)  # bajty, nie tekst
body = p.stdout.decode("utf-8", "replace")
if p.stderr:
    body += "\n[STDERR]\n" + p.stderr.decode("utf-8", "replace")
with open(OUT, "w", encoding="utf-8", errors="replace") as fh:
    fh.write("# a3d_soliton_brannen_r.py  (rc=%d, %.0f s)\n"
             % (p.returncode, time.time() - t0))
    fh.write(body)
print("done rc=%d, %.0f s" % (p.returncode, time.time() - t0))
