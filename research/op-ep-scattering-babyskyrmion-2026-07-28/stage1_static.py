#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE 1 (v2) -- Static baby-skyrmion via 1D hedgehog profile, then VALIDATE in 2D.
  (0) Gradient check: analytic dE/df vs finite-difference (trust the code, not the formula).
  (1) 1D profile relaxes; virial E_skyrme == E_potential (2D scale virial theorem).
  (2) Embed in 2D: Q -> +/-1 and converges toward integer with grid refinement.
  (3) Charge conjugation: antiskyrmion mirrors skyrmion energy.
No physics claim -- this only certifies the arena/solver.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2 = 1.0
MU2 = 0.1

def gradient_check():
    print("-- (0) analytic-gradient check (1D profile energy) --")
    rng = np.random.default_rng(0)
    M = 40; dr = 0.3
    f = np.pi * np.exp(-np.arange(M + 1) * dr / 3.0) + 0.1 * rng.standard_normal(M + 1)
    f[0] = np.pi; f[-1] = 0.0
    g = bs.profile_grad_1d(f, dr, KAPPA2, MU2)
    eps = 1e-6
    gnum = np.zeros_like(f)
    for j in range(1, M):
        fp = f.copy(); fp[j] += eps
        fm = f.copy(); fm[j] -= eps
        gnum[j] = (bs.profile_energy_1d(fp, dr, KAPPA2, MU2)
                   - bs.profile_energy_1d(fm, dr, KAPPA2, MU2)) / (2 * eps)
    err = np.max(np.abs(g[1:M] - gnum[1:M])) / (np.max(np.abs(gnum[1:M])) + 1e-30)
    print(f"    max rel error analytic vs numeric = {err:.3e}  -> {'OK' if err < 1e-5 else 'FAIL'}")
    return err < 1e-5

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE 1 v2: static baby-skyrmion  (kappa^2=%.2f, mu^2=%.2f)" % (KAPPA2, MU2))
    print("=" * 78)
    gok = gradient_check()

    print("\n-- (1) relax 1D profile & 2D virial --")
    r_prof, f_prof = bs.relax_profile_1d(M=2000, Rmax=24.0, kappa2=KAPPA2, mu2=MU2)
    dr = r_prof[1] - r_prof[0]
    E1d = bs.profile_energy_1d(f_prof, dr, KAPPA2, MU2)
    # split 1D energy into skyrme vs potential for virial
    fp = np.gradient(f_prof, dr); s = np.sin(f_prof); rr = r_prof.copy(); rr[0] = 1e-9
    Esk = 2*np.pi*np.sum(rr[1:-1]*dr*0.5*KAPPA2*fp[1:-1]**2*s[1:-1]**2/rr[1:-1]**2)
    Ep = 2*np.pi*np.sum(rr[1:-1]*dr*MU2*(1-np.cos(f_prof[1:-1])))
    print(f"    E_1d = {E1d:.5f}   E_skyrme={Esk:.5f}  E_pot={Ep:.5f}  virial={Esk/Ep:.4f} (->1)")
    print(f"    f(0)={f_prof[0]:.4f} (=pi)  f(Rmax)={f_prof[-1]:.4f} (=0)")

    print("\n-- (2) embed in 2D, grid convergence of Q & E --")
    rows = []
    for N in (128, 192, 256):
        X, Y, dx = bs.make_grid(N, 24.0)
        n = bs.build_from_profile(X, Y, r_prof, f_prof, charge=+1)
        Q = bs.topological_charge(n, dx)
        Et, Eg, Es, Ep2 = bs.total_energy(n, dx, KAPPA2, MU2)
        rows.append((N, dx, Q, Et, Es/Ep2))
        print(f"    N={N:4d} dx={dx:.3f}  Q={Q:+.5f}  E2D={Et:.5f}  virial(2D)={Es/Ep2:.3f}")

    print("\n-- (3) charge conjugation --")
    X, Y, dx = bs.make_grid(192, 24.0)
    npos = bs.build_from_profile(X, Y, r_prof, f_prof, charge=+1)
    nneg = bs.build_from_profile(X, Y, r_prof, f_prof, charge=-1)
    Ep_ = bs.total_energy(npos, dx, KAPPA2, MU2)[0]
    En_ = bs.total_energy(nneg, dx, KAPPA2, MU2)[0]
    Qp = bs.topological_charge(npos, dx); Qn = bs.topological_charge(nneg, dx)
    print(f"    Q(+)= {Qp:+.5f}  Q(-)= {Qn:+.5f}   E(+)={Ep_:.5f}  E(-)={En_:.5f}  |dE|={abs(Ep_-En_):.2e}")

    print("\n" + "=" * 78)
    Qs = [abs(r[2]) for r in rows]
    virs = [r[4] for r in rows]
    ok = (gok and abs(Esk/Ep - 1) < 0.1 and all(abs(q-1) < 0.03 for q in Qs)
          and abs(abs(Qp)-1) < 0.03 and abs(Ep_-En_) < 1e-6)
    print(f"VERDICT stage 1:  grad_ok={gok}  virial_1d={Esk/Ep:.3f}  |Q|grids={[f'{q:.4f}' for q in Qs]}")
    print(f"  charge-conj exact: {abs(Ep_-En_) < 1e-6}   ARENA VALIDATED: {ok}")
    if ok:
        np.savez('soliton_profile.npz', r=r_prof, f=f_prof, kappa2=KAPPA2, mu2=MU2, E=E1d)
        print("  saved soliton_profile.npz (clean soliton for stages 2-3)")
