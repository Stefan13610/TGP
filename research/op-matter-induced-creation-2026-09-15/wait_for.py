#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Pomocnik higieny batchy: AKTYWNE czekanie na marker w pliku
(tura NIE konczy sie w trakcie obliczen -- lekcja poprzednikow).
Uzycie: python wait_for.py <plik> <marker> <timeout_s>"""
import sys
import time

path, marker, tmo = sys.argv[1], sys.argv[2], float(sys.argv[3])
t0 = time.time()
while time.time() - t0 < tmo:
    try:
        with open(path, "r", errors="replace") as fh:
            txt = fh.read()
        if marker in txt:
            print("MARKER '%s' FOUND po %.0f s" % (marker, time.time() - t0))
            sys.exit(0)
    except OSError:
        pass
    time.sleep(5)
print("TIMEOUT po %.0f s (marker '%s' nie znaleziony)" % (tmo, marker))
sys.exit(1)
