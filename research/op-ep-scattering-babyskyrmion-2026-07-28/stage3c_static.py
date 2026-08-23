#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE 3c -- Static interaction energy E_int(d, chi) (clean, no dynamics/init artifacts).
E_int(d) = E(two solitons at separation d) - 2 E_single.
  * Attraction: E_int decreases as d decreases (force pulls together).
  * Repulsion:  E_int increases as d decreases.
Coulomb hypothesis: sign of E_int fixed by Q_A*Q_B, chi-independent, slow (log/1/d) decay.
Dipole/orientation: sign flips with chi, fast (exp) decay (massive field, screened by mu).
This corroborates the dynamic scattering test with a static, momentum-free measurement.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1

def superpose(X, Y, dx, r_prof, f_prof, d, chi, QA, QB):
    vac = np.zeros((3,) + X.shape); vac[2] = 1.0
    nA = bs.build_from_profile(X, Y, r_prof, f_prof, x0=-d/2, y0=0.0, charge=QA, phase=0.0)
    nB = bs.build_from_profile(X, Y, r_prof, f_prof, x0=+d/2, y0=0.0, charge=QB, phase=chi)
    n = bs.normalize(vac + (nA - vac) + (nB - vac))
    return n

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE 3c: static interaction energy (Coulomb vs dipole discriminator)")
    print("=" * 78)
    d = np.load('soliton_profile.npz'); r_prof, f_prof = d['r'], d['f']
    N, L = 256, 64.0
    X, Y, dx = bs.make_grid(N, L)
    Esingle = bs.total_energy(bs.build_from_profile(X, Y, r_prof, f_prof, x0=0, y0=0, charge=1),
                              dx, KAPPA2, MU2)[0]
    print(f"  grid N={N} L={L} dx={dx:.3f}   E_single={Esingle:.4f}\n")
    seps = [4.0, 6.0, 8.0, 10.0, 12.0]
    chis = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    for (QA, QB, tag) in [(-1, +1, 'UNLIKE (-,+)'), (-1, -1, 'LIKE   (-,-)')]:
        print(f"--- {tag}:  E_int(d) = E_pair - 2 E_single ---")
        print(f"    {'d':>5} " + " ".join(f"chi={c:4.2f}" for c in chis))
        for dd in seps:
            row = []
            for chi in chis:
                n = superpose(X, Y, dx, r_prof, f_prof, dd, chi, QA, QB)
                Epair = bs.total_energy(n, dx, KAPPA2, MU2)[0]
                row.append(Epair - 2 * Esingle)
            print(f"    {dd:5.1f} " + " ".join(f"{x:+8.4f}" for x in row))
        print()
    print("Read-off: E_int>0 = repulsive core, and the SIGN of d(E_int)/dd gives the force.")
    print("If sign pattern is chi-dependent (flips across a row) -> orientation/dipole, NOT Coulomb.")
    print("If E_int has same sign & monotonic trend for all chi, set by Q_A*Q_B -> Coulomb-like.")
