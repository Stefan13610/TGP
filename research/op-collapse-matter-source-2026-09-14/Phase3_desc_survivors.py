#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Deskryptywna charakterystyka biegow Q-H2 bez zdarzenia brzegowego
(POST-HOC; zero zmian kryteriow). Stan koncowy, osiadlosc, energie."""
import sys

import numpy as np

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")

DIR = "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14"
RES = DIR + "/Phase3_results"
OUT = []


def w(s=""):
    OUT.append(str(s))


RIDS = ["qh2_qR3_a-0.20__lam0.05__h05", "qh2_qR3_a-0.20__lam0.2__h05",
        "qh2_g_a-0.30_s3__lam0.05__h05", "qh2_g_a-0.30_s3__lam0.2__h05"]
w("Deskryptywnie: biegi Q-H2 bez zdarzenia brzegowego (h05, t=1000)")
for rid in RIDS:
    d = np.load(RES + "/" + rid + ".npz")
    ts, psi0 = d["ts"], d["psi0"]
    ec_, ecf = d["ecore"], d["ecore_field"]
    vmax = d["vmax"]
    pts, profs = d["prof_ts"], d["profs"]
    r = d["r"]
    m80 = r <= 80.0
    mw = pts >= 950.0
    psibar = np.mean(profs[mw], axis=0)
    D = float(np.max(np.abs(psibar[m80] - 1.0)))
    V = float(np.max(vmax[ts >= 950.0]))
    w("%s:" % rid)
    w("  E_core(50)=%.4f E_core(1000)=%.4f ; E_field(50)=%.4f"
      " E_field(1000)=%.4f" % (ec_[np.argmin(np.abs(ts - 50.0))],
                               ec_[-1],
                               ecf[np.argmin(np.abs(ts - 50.0))],
                               ecf[-1]))
    w("  psi(0,1000)=%.6f ; min/max psibar[950,1000]=%.6f/%.6f"
      % (psi0[-1], float(np.min(psibar)), float(np.max(psibar))))
    w("  osiadlosc okno [950,1000]: V=%.3e D=%.3e V<=0.01D: %s"
      % (V, D, V <= 0.01*max(D, 1e-12)))
    w("  min psi(0,t) po t: %.6f @ t=%.1f ; max: %.6f @ t=%.1f"
      % (float(np.min(psi0)), ts[int(np.argmin(psi0))],
         float(np.max(psi0)), ts[int(np.argmax(psi0))]))

with open(DIR + "/Phase3_desc_survivors_output.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
print("desc done")
