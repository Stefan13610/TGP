#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex190_leapfrog_energy_drift_matched_h_family_short.py
=====================================================
P1.C (cd.): **kontrola energii** przy ``leapfrog_integrate`` dla **tej samej** rodziny
IC co ``ex178`` / ``ex186`` (Newton ``v=0`` + trzy TGP z ``H`` dopasowanym do ``H_N``
w swoim potencjale). Te same seede predkosci i ``n_quad`` co w ``ex178``.

UWAGA: jak w ``ex167`` — **bardzo krotkie** ``t_final``; przy dluzszym horyzoncie energia
leapfroga moze silnie dryfowac (Burrau, regulator). To nie jest dowod stabilnosci na ``t``
z Benettina (``ex178`` / ``ex188``).

Obliczenia IC/sil: ``matched_h_family_leapfrog_branches`` w ``ex178_...py``.
Ten sam IC, calka adaptacyjna (referencja Hamiltona): ``ex191_rk45_energy_diag_matched_h_family_short.py``. Obie kolumny naraz: ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``; CSV: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.
Synteza: ``tgp_lyapunov_benettin.tex`` (akapit matched ``H``).

Uruchomienie: ``python ex190_leapfrog_energy_drift_matched_h_family_short.py [--quick] [--n-quad N]``
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_v2 import leapfrog_integrate

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ex178():
    path = os.path.join(_HERE, "ex178_lyapunov_matched_energy_family_table.py")
    spec = importlib.util.spec_from_file_location("_ex178_matched_h", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex178 from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _max_rel_energy_drift(
    pos0: np.ndarray,
    vel0: np.ndarray,
    C: np.ndarray,
    acc,
    pot,
    *,
    t_final: float,
    dt: float,
    save_every: int,
) -> float:
    r = leapfrog_integrate(
        pos0,
        vel0,
        C,
        acc,
        pot,
        t_span=(0.0, t_final),
        dt=dt,
        save_every=save_every,
    )
    return float(np.max(np.abs(r["energy_error"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n-quad", type=int, default=0, help="n_quad Yukawa (0 = domyslnie)")
    args = parser.parse_args()

    ex178 = _load_ex178()
    try:
        branches, meta = ex178.matched_h_family_leapfrog_branches(
            quick=bool(args.quick), n_quad_override=int(args.n_quad)
        )
    except ValueError as e:
        print(f"ex190: FAIL: {e}")
        raise SystemExit(1)

    n_quad = int(meta["n_quad"])
    # Krotki horyzont: przy matched-H predkosciach wieksze niz v=0 (ex167) —
    # dryf rosnie szybko; te parametry utrzymuja max|dE/E0| ponizej progu na calej czworce.
    if args.quick:
        t_final = 0.1
        dt = 0.015
        save_every = 2
        tol = 0.017
    else:
        t_final = 0.11
        dt = 0.014
        save_every = 2
        tol = 0.012

    C = meta["C"]
    print(
        "ex190: leapfrog max |(E-E0)/E0| (ex178 matched-H family)\n"
        f"  t_final={t_final} dt={dt} save_every={save_every} n_quad(Yukawa)={n_quad}"
    )

    all_ok = True
    for lab, pos0, vel, acc, pot in branches:
        d = _max_rel_energy_drift(
            pos0, vel, C, acc, pot, t_final=t_final, dt=dt, save_every=save_every
        )
        ok = d < tol
        flag = "PASS" if ok else "FAIL"
        print(f"  {lab}: max|dE/E0|={d:.6g}  (tol {tol})  {flag}")
        if not ok:
            all_ok = False

    if not all_ok:
        raise SystemExit(1)
    print("ex190: all PASS")


if __name__ == "__main__":
    main()
