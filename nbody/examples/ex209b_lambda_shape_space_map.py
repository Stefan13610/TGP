#!/usr/bin/env python3
"""
ex209b -- Shape-space map of lambda(q1,q2) for Yukawa triple overlap
====================================================================
STATUS: ACTIVE-EXPLORATORY

Maps the large-t exponential rate

    lambda(q1,q2) = min_{alpha in Delta_2} sqrt(Q/Delta)

for the normalized triangle shape-space:

    d_min <= d_mid <= d_max,
    q1 = d_min / d_max,
    q2 = d_mid / d_max,
    0 < q1 <= q2 <= 1,
    q1 + q2 >= 1.

Interpretation:
  - For large t = m_sp * d_max, the exact overlap behaves as
        I_Y(t; q1, q2) ~ exp(-lambda(q1,q2) * t) * [prefactor + corrections].
  - Smaller lambda means weaker exponential suppression and therefore
    relatively larger irreducible 3-body contribution at fixed t.

Outputs:
  - CSV map: `_outputs/ex209b_lambda_shape_space.csv`

Quick mode:
  python -m nbody.examples.ex209b_lambda_shape_space_map --quick
"""

from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.three_body_force_exact import yukawa_overlap_shape_rate

quick = "--quick" in sys.argv


def admissible_q1_values(q2: float, n_q1: int) -> np.ndarray:
    """Uniform q1 samples on the admissible segment for fixed q2."""
    q1_lo = max(1.0 - q2, 1e-6)
    q1_hi = q2
    return np.linspace(q1_lo, q1_hi, n_q1)


def compute_map():
    """Compute lambda(q1,q2) over the admissible shape-space."""
    if quick:
        q2_values = np.linspace(0.55, 1.0, 10)
        n_q1 = 10
    else:
        q2_values = np.linspace(0.52, 1.0, 28)
        n_q1 = 28

    rows = []
    for q2 in q2_values:
        for q1 in admissible_q1_values(float(q2), n_q1):
            rate = yukawa_overlap_shape_rate(float(q1), float(q2))
            rows.append(
                {
                    "q1": float(q1),
                    "q2": float(q2),
                    "lambda": rate["lambda"],
                    "alpha1": rate["alpha1"],
                    "alpha2": rate["alpha2"],
                    "alpha3": rate["alpha3"],
                }
            )
    return rows


def summarize(rows):
    """Print high-level summary of the lambda map."""
    lambdas = np.array([r["lambda"] for r in rows], dtype=float)
    i_min = int(np.argmin(lambdas))
    i_max = int(np.argmax(lambdas))
    r_min = rows[i_min]
    r_max = rows[i_max]

    print("=" * 78)
    print("ex209b -- Shape-space map of lambda(q1,q2)")
    print("=" * 78)
    print(f"  Mode: {'QUICK' if quick else 'FULL'}")
    print(f"  Admissible domain: 0 < q1 <= q2 <= 1, q1 + q2 >= 1")
    print(f"  Samples: {len(rows)}")
    print()
    print(f"  lambda_min = {r_min['lambda']:.6f} at q1={r_min['q1']:.4f}, q2={r_min['q2']:.4f}")
    print(f"    alpha* = ({r_min['alpha1']:.4f}, {r_min['alpha2']:.4f}, {r_min['alpha3']:.4f})")
    print(f"  lambda_max = {r_max['lambda']:.6f} at q1={r_max['q1']:.4f}, q2={r_max['q2']:.4f}")
    print(f"    alpha* = ({r_max['alpha1']:.4f}, {r_max['alpha2']:.4f}, {r_max['alpha3']:.4f})")
    print()

    # Benchmark special shapes.
    specials = [
        ("equilateral", 1.0, 1.0),
        ("isosceles q1=0.5", 0.5, 1.0),
        ("near-degenerate", 0.05, 0.95),
    ]
    print(f"  {'shape':>16s} {'q1':>7s} {'q2':>7s} {'lambda':>10s}")
    print(f"  {'-'*46}")
    for label, q1, q2 in specials:
        if q1 + q2 < 1.0 or q1 > q2:
            continue
        rate = yukawa_overlap_shape_rate(q1, q2)
        print(f"  {label:>16s} {q1:7.3f} {q2:7.3f} {rate['lambda']:10.6f}")


def write_csv(rows):
    """Write the map to CSV."""
    outdir = os.path.join(os.path.dirname(__file__), "_outputs")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "ex209b_lambda_shape_space.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["q1", "q2", "lambda", "alpha1", "alpha2", "alpha3"])
        for r in rows:
            w.writerow(
                [
                    f"{r['q1']:.8f}",
                    f"{r['q2']:.8f}",
                    f"{r['lambda']:.10f}",
                    f"{r['alpha1']:.10f}",
                    f"{r['alpha2']:.10f}",
                    f"{r['alpha3']:.10f}",
                ]
            )
    print()
    print(f"  Written: {path}")


def main():
    t0 = time.time()
    rows = compute_map()
    summarize(rows)
    write_csv(rows)
    print(f"\n  Total time: {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()
