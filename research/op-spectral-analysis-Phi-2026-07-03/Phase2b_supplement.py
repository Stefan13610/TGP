#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase 2b supplement (required by Phase0 LOCK):
 (a) Yukawa F-A nonlinear Newton background pushed to amplitude ~1.0
     (S scan; target amp in {0.3 done, 1.0 missing}).
 (b) Solitons F-S / F-S': count LOCALIZED negative modes, i.e.
     lambda < edge - 1e-3 with edge = -1 (tachyonic continuum of the F-S
     vacuum measured in Phase 2, C5). Continuum count N_box ~ R/pi
     explained analytically: lambda_n = -1 + (n pi / R)^2.
 (c) box-count consistency check: N_neg(vacuum F-S) vs floor(R/pi).
"""
import numpy as np
from Phase2_bvp_spectrum import (newton_background, spectrum_on_background,
                                 soliton_profile, vacuum_case, TOL_NEG,
                                 G0_E, G0_MU, G0_TAU)

print("(a) Yukawa F-A nonlinear, push amplitude toward 1.0:")
for S in (8.0, 12.0, 18.0):
    r, psi, dpsi, d2psi, rho, res = newton_background(S, 60.0, 4000)
    amp = psi.max() - 1
    vals = spectrum_on_background('F-A', r, psi, dpsi, d2psi, extra_rho=rho)
    nneg = int(np.sum(vals < TOL_NEG))
    print("  S=%5.1f amp=%.3f res=%.1e  lam_min=%+.6f  N_neg=%d"
          % (S, amp, res, vals[0], nneg))

print()
print("(b) localized negative modes (lambda < -1 - 1e-3), N=4000, R=60:")
EDGE = -1.0
for form, tag, g0 in (('F-S', 'e', G0_E), ('F-S', 'mu', G0_MU), ('F-S', 'tau', G0_TAU),
                      ('F-Sp', 'sub_e', G0_E), ('F-Sp', 'sub_mu', G0_MU),
                      ('F-Sp', 'sub_tau', G0_TAU)):
    r, g, gp, d2, nb, gmin = soliton_profile(form, g0, R=60.0, N=4000)
    for l in (0, 1):
        vals = spectrum_on_background(form, r, g, gp, d2, l=l)
        nloc = int(np.sum(vals < EDGE - 1e-3))
        below = [round(float(v), 4) for v in vals[vals < EDGE - 1e-3]]
        print("  %-8s %-8s l=%d  N_localized=%d  %s" % (form, tag, l, nloc, below))

print()
print("(c) box-count consistency (vacuum F-S): expect floor(R/pi)")
for R in (40.0, 60.0, 80.0):
    r, vals = vacuum_case('F-S', -1.0, R=R, N=int(R / 60 * 4000))
    nneg = int(np.sum(vals < TOL_NEG))
    print("  R=%.0f: N_neg=%d, floor(R/pi)=%d" % (R, nneg, int(R / np.pi)))
