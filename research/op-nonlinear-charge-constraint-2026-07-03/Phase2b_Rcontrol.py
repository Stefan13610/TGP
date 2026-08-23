#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-nonlinear-charge-constraint (Phase 2, supplement) -- R-control.

LOCK Phase0_balance.md sec. 5: 'zbieznosc: N in {2000,4000,8000}, R=60
(kontrola R in {40,80})' for BORDERLINE results.  Borderline here: the
P2c localized counts against the omega-dependent edge (mu: drop 2 -> 0
between omega = 0.10 and 0.15; tau: 3 -> 1-2), for omega in the
non-ghosted window {0.05..0.25}.  N scaled with R to keep h = R/N
fixed at 60/4000 = 0.015 (as #62 R-control).

Reuses Phase2_qball_family.py definitions verbatim (import).
"""
import numpy as np

import Phase2_qball_family as P2   # executes the main scan too (stdout)

print("\n" + "=" * 78)
print("[R-CONTROL] mu/tau, omega in {0.05..0.25}, R in {40, 80}, h=0.015")
print("  (LOCK sec. 5: kontrola R dla wynikow granicznych; edge(omega)")
print("   niezalezna od R -- z prozni omega)")
print("    %-4s %-6s %-4s %-6s %-28s %-10s %s"
      % ("tag", "omega", "R", "N", "loc(L+) [lam<edge-1e-3]",
         "lam_min(L-)", "uwagi"))
for tag, g0 in (("mu", P2.G0_MU), ("tau", P2.G0_TAU)):
    for om in (0.05, 0.10, 0.15, 0.20, 0.25):
        ed = P2.edges[om]['edge']
        thresh = ed - P2.LOC_MARGIN
        for R in (40.0, 80.0):
            N = int(R / 60.0 * 4000)
            r, g, gp, d2, nb, gmin, filled = P2.profile_hard(
                g0, om, R=R, N=N)
            if not filled:
                print("    %-4s %-6.2f %-4.0f %-6d PROFIL NIEKOMPLETNY"
                      % (tag, om, R, N))
                continue
            dP, eP = P2.build_L(r, g, gp, d2, om, P2.f_log, P2.fp_log,
                                P2.fpp_log, 'plus')
            dM, eM = P2.build_L(r, g, gp, d2, om, P2.f_log, P2.fp_log,
                                P2.fpp_log, 'minus')
            eigsP = P2.low_eigs(dP, eP)
            eigsM = P2.low_eigs(dM, eM, k=4)
            locP = [v for v in eigsP if v < thresh]
            print("    %-4s %-6.2f %-4.0f %-6d %-28s %-10.4f"
                  % (tag, om, R, N,
                     "%d: %s" % (len(locP),
                                 ["%.4f" % v for v in locP[:4]]),
                     float(eigsM[0])))
print("=" * 78)
