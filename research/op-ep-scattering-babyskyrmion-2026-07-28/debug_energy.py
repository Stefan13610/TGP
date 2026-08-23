#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Debug: validate 2D energy code against 1D reduction on a SMOOTH analytic profile."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1

# smooth analytic profile, NO relaxation
Rmax = 20.0; M = 800
r_prof = np.linspace(0, Rmax, M + 1)
f_prof = np.pi * np.exp(-r_prof / 3.0)
f_prof[0] = np.pi
dr = r_prof[1] - r_prof[0]

E1d = bs.profile_energy_1d(f_prof, dr, KAPPA2, MU2)
fp = np.gradient(f_prof, dr); s = np.sin(f_prof); rr = r_prof.copy(); rr[0] = 1e-9
Eg1 = 2*np.pi*np.sum(rr[1:-1]*dr*0.5*(fp[1:-1]**2 + s[1:-1]**2/rr[1:-1]**2))
Es1 = 2*np.pi*np.sum(rr[1:-1]*dr*0.5*KAPPA2*fp[1:-1]**2*s[1:-1]**2/rr[1:-1]**2)
Ep1 = 2*np.pi*np.sum(rr[1:-1]*dr*MU2*(1-np.cos(f_prof[1:-1])))
print(f"1D:  E={E1d:.4f}  Egrad={Eg1:.4f}  Eskyrme={Es1:.4f}  Epot={Ep1:.4f}")

print("\n2D energy of the SAME analytic field, several grids (should -> 1D values):")
for N in (128, 256, 384):
    X, Y, dx = bs.make_grid(N, 24.0)
    n = bs.build_from_profile(X, Y, r_prof, f_prof, charge=+1)
    Et, Eg, Es, Ep = bs.total_energy(n, dx, KAPPA2, MU2)
    Q = bs.topological_charge(n, dx)
    ed, edg, eds, edp = bs.energy_density(n, dx, KAPPA2, MU2)
    imax = np.unravel_index(np.argmax(ed), ed.shape)
    print(f"  N={N:4d} dx={dx:.3f}  E={Et:8.3f} Eg={Eg:7.3f} Es={Es:8.3f} Ep={Ep:6.3f}"
          f"  Q={Q:+.4f}  max e-dens={ed.max():.2f} @ {imax}"
          f"  r@max={np.hypot(X[imax],Y[imax]):.3f}")
