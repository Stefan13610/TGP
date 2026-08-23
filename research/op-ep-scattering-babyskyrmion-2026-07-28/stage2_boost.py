#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE 2 -- Boosted single soliton: certify the LOSSLESS dynamics.
A moving soliton must:
  (1) conserve total energy (this MEASURES the 'lossless' claim, not assumes it),
  (2) conserve topological charge Q,
  (3) travel at ~constant velocity ~= v (COM of |q|), shape preserved (low radiation).
If energy drifts a lot or the soliton radiates/decays, the arena is not trustworthy
for scattering. No physics claim here either.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1

def run(N=192, L=40.0, v=0.3, T=30.0, cfl=0.15):
    d = np.load('soliton_profile.npz')
    r_prof, f_prof = d['r'], d['f']
    X, Y, dx = bs.make_grid(N, L)
    dt = cfl * dx
    n, pi = bs.boost_soliton(X, Y, r_prof, f_prof, x0=-L*0.25, y0=0.0, v=v, charge=+1)

    E0 = bs.total_energy_dyn(n, pi, dx, KAPPA2, MU2)
    Q0 = bs.topological_charge(n, dx)
    Px0, Py0 = bs.field_momentum(n, pi, dx)
    xc0, yc0, _ = bs.com_of_charge(n, X, Y, dx)

    print(f"N={N} L={L} dx={dx:.3f} dt={dt:.4f} v={v}  gamma={1/np.sqrt(1-v*v):.4f}")
    print(f"  t=0:  E={E0:.5f}  Q={Q0:+.4f}  Px={Px0:.5f}  xc={xc0:.3f}")
    print(f"  {'t':>6} {'E':>10} {'dE/E0':>10} {'Q':>8} {'xc':>8} {'v_meas':>8} {'Px':>9}")

    nsteps = int(T / dt)
    log_every = max(1, nsteps // 12)
    xc_prev, t_prev = xc0, 0.0
    traj = []
    for it in range(1, nsteps + 1):
        n, pi = bs.rk4_step(n, pi, dt, dx, KAPPA2, MU2)
        t = it * dt
        if it % log_every == 0 or it == nsteps:
            E = bs.total_energy_dyn(n, pi, dx, KAPPA2, MU2)
            Q = bs.topological_charge(n, dx)
            xc, yc, W = bs.com_of_charge(n, X, Y, dx)
            Px, Py = bs.field_momentum(n, pi, dx)
            vmeas = (xc - xc_prev) / (t - t_prev)
            xc_prev, t_prev = xc, t
            traj.append((t, E, Q, xc, vmeas, Px))
            print(f"  {t:6.2f} {E:10.5f} {(E-E0)/E0:10.2e} {Q:8.4f} {xc:8.3f} {vmeas:8.4f} {Px:9.5f}")

    Efin = traj[-1][1]
    print("\n  VERDICT stage 2:")
    print(f"    energy drift |dE/E0| max = {max(abs((e[1]-E0)/E0) for e in traj):.2e}")
    print(f"    Q conserved: {abs(traj[-1][2]-Q0):.2e}")
    vmeans = np.mean([e[4] for e in traj[2:]])
    print(f"    mean measured velocity = {vmeans:.4f}  (injected v={v})")
    ok = (max(abs((e[1]-E0)/E0) for e in traj) < 0.02 and abs(traj[-1][2]-Q0) < 0.02
          and abs(vmeans - v) < 0.05)
    print(f"    DYNAMICS TRUSTWORTHY: {ok}")
    return ok

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE 2: boosted single soliton (lossless dynamics check)")
    print("=" * 78)
    run(N=192, L=40.0, v=0.3, T=30.0, cfl=0.15)
