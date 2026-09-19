#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Deskryptywny probe lokusa zdarzen brzegowych (POST-HOC, zero zmian
kryteriow/progow; wylacznie do konfrontacji z P1-H3 w FINAL).
Powtarza wybrane krotkie biegi i raportuje, GDZIE (r) pole pierwsze
przekracza pas klasyfikacyjny i jaki jest stan w centrum."""
import sys

import numpy as np

sys.path.insert(0, "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14")
import engine_core as ec

DIR = "TGP/TGP_v1/research/op-collapse-matter-source-2026-09-14"
OUT = []


def w(s=""):
    OUT.append(str(s))


def probe(label, lam, start, h=0.05, dt=0.005, tmax=30.0):
    eng = ec.Engine(h, 200.0, sponge=True, lam=lam)
    if start == "vac":
        g = ec.start_vacuum(eng)
    elif start.startswith("qR3"):
        g = ec.start_quasiR3(eng, float(start.split("_a")[1]))
    else:
        _, a, s = start.split("_")
        g = ec.start_gauss(eng, float(a[1:]), float(s[1:]))
    pi = np.zeros(eng.N)
    prev = g.copy()
    for n in range(int(round(tmax/dt))):
        try:
            g, pi = eng.step(g, pi, dt)
        except ec.NonConvergence:
            w("%s: NonConvergence t=%.3f" % (label, (n + 1)*dt))
            g = prev
            break
        st = eng.band_status(g)
        if st is not None:
            t = (n + 1)*dt
            if not np.all(np.isfinite(g)):
                idx = int(np.argmax(~np.isfinite(g)))
                w("%s: %s t=%.3f pierwszy niefinityczny r=%.3f"
                  % (label, st, t, eng.r[idx]))
            else:
                imx = int(np.argmax(g))
                imn = int(np.argmin(g))
                w("%s: %s t=%.3f" % (label, st, t))
                w("   max psi=%.6f @ r=%.3f (rho_hat=%.3f);"
                  " min psi=%.6e @ r=%.3f (rho_hat=%.3f)"
                  % (float(g[imx]), eng.r[imx],
                     float(np.exp(-eng.r[imx]**2/18.0)),
                     float(g[imn]), eng.r[imn],
                     float(np.exp(-eng.r[imn]**2/18.0))))
                w("   psi(0)=%.6f" % float(g[0]))
            return
        prev = g.copy()
    else:
        w("%s: bez zdarzenia do t=%.1f" % (label, tmax))


w("Deskryptywny probe lokusa zdarzen brzegowych (POST-HOC)")
probe("QH1 lam=0.2  vac", 0.2, "vac")
probe("QH1 lam=0.5  vac", 0.5, "vac")
probe("QH2 qR3+0.20 lam=0.05", 0.05, "qR3_a+0.20")
probe("QH2 qR3-0.20 lam=0.5", 0.5, "qR3_a-0.20")
probe("QH2 g+0.15s6 lam=0.05", 0.05, "g_a+0.15_s6")
probe("QH2 g-0.30s3 lam=0.5", 0.5, "g_a-0.30_s3")

with open(DIR + "/Phase3_desc_locus_output.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
print("locus done")
