#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex191_rk45_energy_diag_matched_h_family_short.py
==============================================
P1.C (cd.): **diagnostyka energii** ``rk45_integrate`` (DOP853) dla **tej samej**
rodziny IC co ``ex178`` / ``ex190`` — cztery galezie (Newton + trzy TGP matched ``H``).
Te same ``(acc, pot)`` co leapfrog w ``matched_h_family_leapfrog_branches`` (``ex178``).

Tryb pelny: ``t_final=0.12`` (``ex192`` pelny uzywa ``0.11`` dla wspolnej kolumny z leapfrog).
Na **krotkim** ``t_final`` solver zwykle zwraca ``success=True`` i bardzo maly
``max|(E-E0)/E0|`` dla TGP (orbita Hamiltonowska w ciaglej dynamice); Newton
ma wiekszy blad na tej samej skali czasu niz TGP — to **nie** jest ten sam obiekt
co trajektoria leapfroga z Benettina (``ex178``, ``ex188``). Por. stdout RK45
w ``ex152`` / ``ex164`` przy dluzszym ``t`` (czesto ``success=False`` mimo
malego ``|dE/E0|`` na ``t_eval``).

Obliczenia IC: ``matched_h_family_leapfrog_branches`` w ``ex178_...py``.
Porownanie z dyskretnym leapfrogiem: ``ex190_leapfrog_energy_drift_matched_h_family_short.py``. Tabela obok siebie (wspolny ``t_final``): ``ex192_matched_h_leapfrog_vs_rk45_energy_table.py``; CSV: ``ex193_matched_h_lf_vs_rk45_energy_csv.py``; odczyt: ``ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py``.
Synteza (Benettin vs LF vs RK45, matched ``H``): ``tgp_lyapunov_benettin.tex``.

Uruchomienie: ``python ex191_rk45_energy_diag_matched_h_family_short.py [--quick] [--n-quad N]``
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

from nbody.dynamics_v2 import rk45_integrate

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ex178():
    path = os.path.join(_HERE, "ex178_lyapunov_matched_energy_family_table.py")
    spec = importlib.util.spec_from_file_location("_ex178_matched_h", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load ex178 from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
        print(f"ex191: FAIL: {e}")
        raise SystemExit(1)

    n_quad = int(meta["n_quad"])
    C = meta["C"]

    if args.quick:
        t_final = 0.1
        n_output = 80
        newton_max_tol = 0.001
    else:
        t_final = 0.12
        n_output = 96
        newton_max_tol = 0.0015

    rtol = 1e-9
    atol = 1e-11

    print(
        "ex191: RK45 (DOP853) energy diag — ex178 matched-H family\n"
        f"  t_final={t_final} n_output={n_output} rtol={rtol:g} atol={atol:g} "
        f"n_quad(Yukawa)={n_quad}"
    )

    all_ok = True
    for lab, pos0, vel, acc, pot in branches:
        r = rk45_integrate(
            pos0,
            vel,
            C,
            acc,
            pot,
            t_span=(0.0, float(t_final)),
            n_output=int(n_output),
            rtol=rtol,
            atol=atol,
            quiet=True,
        )
        ee = r["energy_error"]
        if not np.all(np.isfinite(ee)):
            print(f"  {lab}: FAIL non-finite energy_error")
            all_ok = False
            continue
        mx = float(np.max(np.abs(ee)))
        ok_s = bool(r["success"])
        if not ok_s:
            print(
                f"  {lab}: success=False max|(E-E0)/E0|={mx:.3e}  "
                f"message={r.get('message', '')!r}  FAIL"
            )
            all_ok = False
            continue
        extra = ""
        if lab == "Newton":
            if mx > newton_max_tol:
                print(
                    f"  {lab}: success=True max|(E-E0)/E0|={mx:.3e}  "
                    f"(Newton tol {newton_max_tol})  FAIL"
                )
                all_ok = False
                continue
            extra = f"  (Newton tol {newton_max_tol})"
        print(
            f"  {lab}: success=True max|(E-E0)/E0|={mx:.3e}{extra}"
        )

    if not all_ok:
        raise SystemExit(1)
    print("ex191: PASS (finite energy, success=True, Newton |dE/E0| bound)")


if __name__ == "__main__":
    main()
