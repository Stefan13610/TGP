#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE 3b -- Large-impact-parameter scan: the DECISIVE discriminator.
  * Coulomb (if it emerged): LONG-range, force sign = Q_A*Q_B, chi-INDEPENDENT.
    -> at large b, deflection stays sign-definite by charge product, does NOT flip with chi,
       and falls off slowly (~1/b).
  * Dipole/orientation (baby-Skyrme expectation): SHORT-range (screened by mass mu),
    sign flips with chi, falls off fast (exp) with b.
This is where Coulomb would DOMINATE if present. If sign still flips with chi and decays
fast -> operational charge is NOT Coulombic; the sec.4 falsifier returns NEGATIVE.
"""
import numpy as np
import stage3_scatter as s   # sets utf-8 stdout on import

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE 3b: large-b scan (Coulomb long-range vs dipole short-range)")
    print("=" * 78)
    print("dvAy<0 = A deflects toward partner (attract); >0 = away (repel).")
    print("Coulomb => UNLIKE all<0, LIKE all>0, chi-independent, slow 1/b decay.\n")
    v = 0.4
    chis = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    for b in (5.0, 8.0):
        for (QA, QB, tag) in [(-1, +1, 'UNLIKE'), (-1, -1, 'LIKE  ')]:
            row = []
            for chi in chis:
                r = s.run_collision(QA, QB, chiA=0.0, chiB=chi, b=b, v=v,
                                    N=224, L=64.0, d0=14.0, T=60.0)
                row.append(r['dvAy'])
            signs = ''.join('-' if x < -0.01 else ('+' if x > 0.01 else '0') for x in row)
            print(f"  b={b:.0f} {tag}: dvAy(chi=0,pi/2,pi,3pi/2) = "
                  f"[{row[0]:+.4f} {row[1]:+.4f} {row[2]:+.4f} {row[3]:+.4f}]  signs={signs}")
        print()
    print("Read-off: if UNLIKE signs are not all '-' and LIKE not all '+', OR signs flip")
    print("with chi -> NO Coulomb charge; interaction is orientation/dipole dominated.")
