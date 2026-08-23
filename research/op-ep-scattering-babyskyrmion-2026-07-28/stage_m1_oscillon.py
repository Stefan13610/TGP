#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAGE M1a -- GATE for option A: does a long-lived neutral (Q=0) oscillon exist?
Build a Q=0 lump, evolve it freely, measure whether energy stays LOCALIZED for a time
>> the scattering experiment (T~50-100).
  * If a localized breather persists (slow energy loss) -> neutral 'balloon' exists -> A viable.
  * If it disperses fast -> no stable neutral grain in this arena -> A blocked at M1 (a finding).
We scan seed amplitude g0 (bump depth) to find the best oscillon, if any.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import bs2d as bs

KAPPA2, MU2 = 1.0, 0.1

def localized_fraction(n, pi, X, Y, dx, x0, y0, Rloc):
    ed = bs.energy_density(n, dx, KAPPA2, MU2)[0] + 0.5 * np.sum(pi * pi, axis=0)
    cell = dx * dx
    Etot = ed.sum() * cell
    mask = ((X - x0) ** 2 + (Y - y0) ** 2) < Rloc ** 2
    Eloc = (ed * mask).sum() * cell
    return Etot, Eloc, Eloc / Etot

def run_osc(g0, w=2.0, N=192, L=48.0, T=100.0, cfl=0.15, Rloc=6.0):
    X, Y, dx = bs.make_grid(N, L)
    dt = cfl * dx
    n = bs.bump_Q0(X, Y, g0=g0, w=w)
    pi = np.zeros_like(n)
    n, pi = bs.project(n, pi)
    Q = bs.topological_charge(n, dx)
    E0, El0, frac0 = localized_fraction(n, pi, X, Y, dx, 0, 0, Rloc)
    nsteps = int(T / dt)
    log_every = max(1, nsteps // 10)
    print(f"  g0={g0:.2f} w={w}:  Q={Q:+.4f}  E0={E0:.4f}  E_loc/E0(t=0)={frac0:.3f}")
    fr_hist = [frac0]
    for it in range(1, nsteps + 1):
        n, pi = bs.rk4_step(n, pi, dt, dx, KAPPA2, MU2)
        if it % log_every == 0 or it == nsteps:
            Et, El, fr = localized_fraction(n, pi, X, Y, dx, 0, 0, Rloc)
            fr_hist.append(fr)
            ed = bs.energy_density(n, dx, KAPPA2, MU2)[0]
            print(f"    t={it*dt:6.1f}  E_loc/E0={fr:6.3f}  peak_e={ed.max():7.3f}  Q={bs.topological_charge(n,dx):+.4f}")
    return fr_hist[-1]

if __name__ == '__main__':
    print("=" * 78)
    print("STAGE M1a: neutral (Q=0) oscillon existence gate")
    print("=" * 78)
    print("E_loc = energy within Rloc=6 of origin. Oscillon => stays high (localized) at t=100.")
    print("Free dispersion => E_loc/E0 decays toward the 'vacuum' baseline (spreads out).\n")
    results = {}
    for g0 in (1.0, 2.0, 2.6):
        results[g0] = run_osc(g0)
        print()
    print("=" * 78)
    print("GATE VERDICT:")
    for g0, fr in results.items():
        verdict = "LOCALIZED (oscillon)" if fr > 0.6 else ("partial" if fr > 0.3 else "DISPERSED")
        print(f"  g0={g0:.2f}: E_loc/E0(t=100) = {fr:.3f}  -> {verdict}")
    best = max(results.values())
    print(f"\n  Neutral balloon viable in this arena: {best > 0.6}")
    print("  (If not, option A is blocked at M1 -- itself an honest finding about TGP ontology.)")
