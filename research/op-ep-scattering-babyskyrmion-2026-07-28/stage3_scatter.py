#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE 3 -- Two-soliton scattering: the RAMY sec.4 falsifier.
Question: does 'charge' (operational momentum transfer) show COULOMB SELECTIVITY,
i.e. sign of the deflection fixed by the product of topological charges
(unlike attract, like repel), INDEPENDENT of internal orientation chi?
  * If sign is Q-product-definite for ALL chi  -> Coulomb-like charge EMERGES.
  * If sign flips with chi  -> orientation-dominated (dipole), NO Coulomb charge.
We INSERTED no pairwise law; interaction is only the O(3)+Skyrme field overlap.

Diagnostic observable: transverse deflection dv_y of projectile A and target B,
plus longitudinal slowdown. Attraction bends trajectories together; repulsion apart.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1

def load_profile():
    d = np.load('soliton_profile.npz')
    return d['r'], d['f']

def two_soliton_init(X, Y, dx, r_prof, f_prof, d0, b, v, QA, QB, chiA, chiB):
    """A: left at (-d0, +b/2) moving +x with v. B: right at (+d0, -b/2) at rest."""
    gamma = 1.0 / np.sqrt(1.0 - v * v)
    XA = -d0 + (X + d0) * gamma
    nA = bs.build_from_profile(XA, Y, r_prof, f_prof, x0=-d0, y0=+b/2, charge=QA, phase=chiA)
    piA = -v * bs.d_dx(nA, dx)
    nB = bs.build_from_profile(X, Y, r_prof, f_prof, x0=+d0, y0=-b/2, charge=QB, phase=chiB)
    # combine via stereographic-like overlay: pick field of whichever soliton is closer;
    # simplest robust superposition on S^2: add tangent deviations from vacuum, renormalize.
    vac = np.zeros_like(nA); vac[2] = 1.0
    n = vac + (nA - vac) + (nB - vac)
    n = bs.normalize(n)
    pi = piA.copy()               # only A carries momentum initially
    n, pi = bs.project(n, pi)
    return n, pi

def track_two(n, X, Y, dx, c1, c2, R=6.0):
    q = np.abs(bs.charge_density(n, dx))
    out = []
    for (cx, cy) in (c1, c2):
        mask = ((X - cx) ** 2 + (Y - cy) ** 2) < R * R
        w = (q * mask)
        W = w.sum()
        if W < 1e-9:
            out.append((cx, cy, 0.0)); continue
        out.append((np.sum(w * X) / W, np.sum(w * Y) / W, W * dx * dx))
    return out

def run_collision(QA, QB, chiA, chiB, b, v=0.4, d0=12.0, N=192, L=52.0, T=52.0,
                  cfl=0.15, verbose=False):
    r_prof, f_prof = load_profile()
    X, Y, dx = bs.make_grid(N, L)
    dt = cfl * dx
    n, pi = two_soliton_init(X, Y, dx, r_prof, f_prof, d0, b, v, QA, QB, chiA, chiB)
    E0 = bs.total_energy_dyn(n, pi, dx, KAPPA2, MU2)
    Px0, Py0 = bs.field_momentum(n, pi, dx)
    cA = (-d0, +b/2); cB = (+d0, -b/2)
    # measure incoming velocity of A over first segment
    hist = []
    nsteps = int(T / dt)
    for it in range(1, nsteps + 1):
        n, pi = bs.rk4_step(n, pi, dt, dx, KAPPA2, MU2)
        if it % 20 == 0 or it == nsteps:
            tr = track_two(n, X, Y, dx, cA, cB)
            cA = (tr[0][0], tr[0][1]); cB = (tr[1][0], tr[1][1])
            hist.append((it * dt, cA[0], cA[1], cB[0], cB[1], tr[0][2], tr[1][2]))
    hist = np.array(hist)
    Efin = bs.total_energy_dyn(n, pi, dx, KAPPA2, MU2)
    Pxf, Pyf = bs.field_momentum(n, pi, dx)
    # velocities from finite differences of trajectories (early vs late)
    t = hist[:, 0]
    def vel(col):
        # linear fit slope over first 25% and last 25%
        k = max(3, len(t) // 4)
        vin = np.polyfit(t[:k], hist[:k, col], 1)[0]
        vout = np.polyfit(t[-k:], hist[-k:, col], 1)[0]
        return vin, vout
    vAx_in, vAx_out = vel(1); vAy_in, vAy_out = vel(2)
    vBx_in, vBx_out = vel(3); vBy_in, vBy_out = vel(4)
    res = dict(QA=QA, QB=QB, chiA=chiA, chiB=chiB, b=b,
               dvAy=vAy_out - vAy_in, dvBy=vBy_out - vBy_in,
               dvAx=vAx_out - vAx_in, dvBx=vBx_out - vBx_in,
               vBy_out=vBy_out, vBx_out=vBx_out,
               Edrift=(Efin - E0) / E0, Pxdrift=Pxf - Px0,
               wA_fin=hist[-1, 5], wB_fin=hist[-1, 6])
    if verbose:
        print(f"    E drift={res['Edrift']:.2e} Px drift={res['Pxdrift']:.2e}"
              f"  |Q|A_fin={res['wA_fin']:.2f} |Q|B_fin={res['wB_fin']:.2f}")
    return res

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE 3: two-soliton scattering -- Coulomb-selectivity test")
    print("=" * 78)
    print("Interpretation key: B target starts at rest. If A ATTRACTS B, target B is")
    print("pulled toward A's incoming line; sign of transverse kick on B (vBy_out) and")
    print("projectile deflection dvAy reveal attract(+together)/repel(+apart).")
    print()
    b = 3.0; v = 0.4
    chis = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    for (QA, QB, tag) in [(-1, +1, 'UNLIKE (-,+)'), (-1, -1, 'LIKE (-,-)')]:
        print(f"--- {tag}  (b={b}, v={v}) ---")
        print(f"    {'chiA':>6} {'chiB':>6} {'dvAy':>9} {'dvBy':>9} {'vBx_out':>9} {'Edrift':>9}")
        for chi in chis:
            r = run_collision(QA, QB, chiA=0.0, chiB=chi, b=b, v=v)
            print(f"    {0.0:6.2f} {chi:6.2f} {r['dvAy']:9.4f} {r['dvBy']:9.4f}"
                  f" {r['vBx_out']:9.4f} {r['Edrift']:9.1e}")
        print()
