#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Pomocnik higieny batchy: aktywne czekanie na marker w pliku.
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
            print("MARKER FOUND po %.0f s" % (time.time() - t0))
            sys.exit(0)
    except OSError:
        pass
    time.sleep(10)
print("TIMEOUT po %.0f s" % tmo)
sys.exit(1)
