#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex152_lyapunov_yukawa_feynman_short.py
======================================
P1 (cd.): ``lambda_max`` z backendu ``yukawa_feynman`` (dokladne I_Y + V2),
N=3, IC Burrau — **krotki** horyzont czasu i umiarkowane ``n_quad`` (koszt FD
Jacobiana: O(n_quad) na kazde wywolanie ``acc``).

Porownanie referencyjne: ``coulomb_3b`` przy tych samych ``beta,gamma,soft``.

Po estymacie ``lambda_max`` (RK4+Benettin): **diagnostyka energii** przez
``rk45_integrate`` (DOP853) na tym samym ``[0, t_final]`` — tylko stdout;
Benettin nadal na RK4. ``solve_ivp`` czesto zwraca ``success=False`` przy
bliskich przejsciach Burrau, mimo ze ``max|dE/E0|`` na zapisanych punktach
bywa rzedu 1e-8 (por. ``ex167`` dla leapfroga na krotkim ``t``).

Uruchomienie: ``python ex152_lyapunov_yukawa_feynman_short.py [--quick]``

Oczekiwany czas: ``--quick`` ~0.5–3 min (CPU); pelny tryb ~kilka minut.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.dynamics_backends import build_tgp_integration_pair
from nbody.dynamics_v2 import rk45_integrate
from nbody.lyapunov import largest_lyapunov_exponent_benettin, pythagorean_three_body_burrau


def _print_rk45_energy_diag(
    *,
    label: str,
    pos0: np.ndarray,
    vel0: np.ndarray,
    C: np.ndarray,
    t_final: float,
    beta: float,
    gamma: float,
    soft: float,
    n_quad: int,
) -> None:
    acc, pot = build_tgp_integration_pair(
        label,
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )
    r = rk45_integrate(
        pos0,
        vel0,
        C,
        acc,
        pot,
        t_span=(0.0, float(t_final)),
        n_output=400,
        rtol=1e-9,
        atol=1e-11,
        quiet=True,
    )
    mx = float(np.max(np.abs(r["energy_error"])))
    print(
        f"  rk45 energy diag {label}: success={r['success']}  "
        f"max|(E-E0)/E0|={mx:.3e}  (same t_final as Benettin)"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--n-quad",
        type=int,
        default=0,
        help="n_quad Feynman (0 = domyslne wg trybu quick/full)",
    )
    args = parser.parse_args()

    beta = gamma = 0.08
    soft = 1e-6
    pos0, vel0, C = pythagorean_three_body_burrau()

    if args.quick:
        t_final = 1.35
        dt = 0.05
        renorm = 5
        jac_eps = 1e-4
        n_quad = 12
    else:
        t_final = 2.5
        dt = 0.04
        renorm = 8
        jac_eps = 7e-5
        n_quad = 16
    if args.n_quad > 0:
        n_quad = int(args.n_quad)

    acc_y, _ = build_tgp_integration_pair(
        "yukawa_feynman",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
        n_quad_feynman=n_quad,
    )
    acc_c, _ = build_tgp_integration_pair(
        "coulomb_3b",
        beta=beta,
        gamma=gamma,
        softening=soft,
        include_3body=True,
    )

    print(
        "ex152: lambda_max yukawa_feynman vs coulomb_3b (Burrau, N=3), "
        f"n_quad={n_quad}"
    )
    print(f"  t_final={t_final} dt={dt} renorm={renorm} jac_eps={jac_eps}")

    t0 = time.perf_counter()
    lam_y, sy = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_y,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=np.random.default_rng(1521),
    )
    t1 = time.perf_counter()
    lam_c, sc = largest_lyapunov_exponent_benettin(
        pos0,
        vel0,
        C,
        acc_c,
        t_final=t_final,
        dt=dt,
        renorm_every=renorm,
        jac_eps=jac_eps,
        rng=np.random.default_rng(1522),
    )
    t2 = time.perf_counter()

    print(f"  yukawa_feynman lambda={lam_y:.6f} steps={sy} cpu={t1 - t0:.2f}s")
    print(f"  coulomb_3b     lambda={lam_c:.6f} steps={sc} cpu={t2 - t1:.2f}s")
    print(f"  delta          {lam_y - lam_c:+.6f}")

    _print_rk45_energy_diag(
        label="yukawa_feynman",
        pos0=pos0,
        vel0=vel0,
        C=C,
        t_final=t_final,
        beta=beta,
        gamma=gamma,
        soft=soft,
        n_quad=n_quad,
    )
    _print_rk45_energy_diag(
        label="coulomb_3b",
        pos0=pos0,
        vel0=vel0,
        C=C,
        t_final=t_final,
        beta=beta,
        gamma=gamma,
        soft=soft,
        n_quad=n_quad,
    )

    ok = np.isfinite(lam_y) and np.isfinite(lam_c) and lam_y > 0.02 and lam_c > 0.02
    if not ok:
        raise SystemExit(1)
    print("  PASS (finite, positive chaotic-scale lambda)")


if __name__ == "__main__":
    main()
