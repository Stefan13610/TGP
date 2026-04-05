#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex140_yukawa_feynman_N4_scaling_rk.py
=====================================
Kontynuacja ``ex139``:

1. **Skalowanie kosztu** — ten sam backend ``yukawa_feynman`` i ``n_quad``:
   pomiar czasu wielokrotnej ewaluacji ``acc_fn(pos, C)`` dla $N=3$ vs $N=4$
   (dla $N=4$ jest $\\binom{4}{3}=4$ tripletów vs $1$ przy $N=3$).

2. **Leapfrog** — krótka trajektoria dla $N=4$ (umiarkowany próg energii).

3. **RK45 (DOP853)** — bardzo krótki odcinek czasu z tym samym backendem
   (adaptacyjny krok; mniej kroków jawnych niż przy długim leapfrog).

Uruchomienie z katalogu ``examples/`` (``TGP_v1/`` na ``sys.path``), jak ``ex138``/``ex139``.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody import dynamics_v2


def _bench_acc(acc_fn, pos: np.ndarray, C: np.ndarray, repeats: int) -> float:
    t0 = time.perf_counter()
    for _ in range(repeats):
        _ = acc_fn(pos, C)
    return time.perf_counter() - t0


def main() -> None:
    beta = gamma = 1.0
    eps = 1e-5
    n_quad = 20
    repeats = 60

    # N=3: ta sama geometria co ex139
    a = 3.0
    pos3 = np.array(
        [
            [0.0, 0.0, 0.0],
            [a, 0.02, 0.0],
            [0.5 * a, 0.866 * a, 0.01],
        ],
        dtype=float,
    )
    C3 = np.array([0.25, 0.25, 0.25], dtype=float)

    # N=4: lekko zaburzony czworościan
    pos4 = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.1, 0.0, 0.0],
            [1.0, 1.75, 0.0],
            [1.05, 0.55, 1.45],
        ],
        dtype=float,
    )
    C4 = np.array([0.22, 0.22, 0.22, 0.22], dtype=float)

    acc_y, _ = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=eps,
        include_3body=True,
        n_quad_feynman=n_quad,
    )

    t3 = _bench_acc(acc_y, pos3, C3, repeats)
    t4 = _bench_acc(acc_y, pos4, C4, repeats)
    ratio = t4 / t3 if t3 > 0 else float("nan")

    print("ex140: yukawa_feynman - skalowanie N=3 vs N=4")
    print(f"  n_quad={n_quad}  repeats={repeats}")
    print(f"  wall acc-only: N=3 -> {t3:.4f}s   N=4 -> {t4:.4f}s   ratio t4/t3 = {ratio:.2f}")

    # --- Leapfrog N=4 ---
    vel0 = np.zeros_like(pos4)
    acc_lf, pot_lf = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=eps,
        include_3body=True,
        n_quad_feynman=n_quad,
    )
    t0 = time.perf_counter()
    out_lf = dynamics_v2.leapfrog_integrate(
        pos4,
        vel0,
        C4,
        acc_lf,
        pot_lf,
        t_span=(0.0, 0.08),
        dt=5e-5,
        save_every=25,
    )
    wall_lf = time.perf_counter() - t0
    E = out_lf["energy"]
    dE = float(np.max(E) - np.min(E))
    rel_lf = dE / (abs(np.mean(E)) + 1e-30)
    ok_lf = rel_lf < 2e-4
    print(f"  leapfrog N=4: wall={wall_lf:.2f}s  max|dE|/mean|E|={rel_lf:.3g}  {'PASS' if ok_lf else 'FAIL'}")

    # --- RK45 krótki ---
    acc_rk, pot_rk = acc_lf, pot_lf
    t0 = time.perf_counter()
    out_rk = dynamics_v2.rk45_integrate(
        pos4,
        vel0,
        C4,
        acc_rk,
        pot_rk,
        t_span=(0.0, 0.05),
        n_output=24,
        rtol=1e-9,
        atol=1e-11,
    )
    wall_rk = time.perf_counter() - t0
    Er = out_rk["energy"]
    dEr = float(np.max(Er) - np.min(Er))
    rel_rk = dEr / (abs(np.mean(Er)) + 1e-30)
    ok_rk = bool(out_rk["success"]) and rel_rk < 1e-5
    print(f"  rk45 N=4: wall={wall_rk:.2f}s  max|dE|/mean|E|={rel_rk:.3g}  {'PASS' if ok_rk else 'FAIL'}")

    if not (ok_lf and ok_rk):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
