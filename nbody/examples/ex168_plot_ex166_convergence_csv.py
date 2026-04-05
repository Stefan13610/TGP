#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex168_plot_ex166_convergence_csv.py
===================================
P1 (cd.): **odczyt i interpretacja** CSV z ``ex166`` (siatka ``dt`` x ``n_quad``,
``lambda_max`` przy tym samym seedzie). Opcjonalnie zapisuje PNG (wymaga
``matplotlib``).

Domyslne wejscie: ``examples/_outputs/ex166_yukawa_convergence_grid.csv``
(wygeneruj najpierw ``ex166``, albo uruchom ``verify_nbody_lyapunov_quick`` —
``ex166`` jest przed ``ex168``).

``--quick``: tylko walidacja CSV + statystyki na stdout (bez matplotlib).

Uruchomienie:
``python ex168_plot_ex166_convergence_csv.py [--quick] [--csv PATH] [--png PATH]``
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "_outputs", "ex166_yukawa_convergence_grid.csv")
_DEFAULT_PNG = os.path.join(_HERE, "_outputs", "ex168_ex166_lambda_grid.png")


def _load_rows(path: str) -> list[dict[str, float | int]]:
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        need = {
            "t_final",
            "dt",
            "n_quad",
            "renorm_every",
            "seed",
            "lambda_max",
            "steps",
        }
        if r.fieldnames is None or not need.issubset(set(r.fieldnames)):
            raise ValueError(f"CSV {path}: expected columns {sorted(need)}")
        rows: list[dict[str, float | int]] = []
        for row in r:
            rows.append(
                {
                    "t_final": float(row["t_final"]),
                    "dt": float(row["dt"]),
                    "n_quad": int(row["n_quad"]),
                    "renorm_every": int(row["renorm_every"]),
                    "seed": int(row["seed"]),
                    "lambda_max": float(row["lambda_max"]),
                    "steps": int(row["steps"]),
                }
            )
    return rows


def _pivot(
    rows: list[dict[str, float | int]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dts = sorted({float(r["dt"]) for r in rows})
    nqs = sorted({int(r["n_quad"]) for r in rows})
    z = np.full((len(dts), len(nqs)), np.nan, dtype=float)
    key = {(float(r["dt"]), int(r["n_quad"])): float(r["lambda_max"]) for r in rows}
    for i, dt in enumerate(dts):
        for j, nq in enumerate(nqs):
            z[i, j] = key.get((dt, nq), np.nan)
    return np.array(dts), np.array(nqs, dtype=int), z


def _print_interpretation(
    rows: list[dict[str, float | int]],
    dts: np.ndarray,
    nqs: np.ndarray,
    z: np.ndarray,
) -> None:
    lams = np.array([float(r["lambda_max"]) for r in rows], dtype=float)
    mean_l = float(np.mean(lams))
    spread = (float(np.max(lams)) - float(np.min(lams))) / mean_l if mean_l > 0 else float("inf")
    print(f"  rows={len(rows)}  lambda_max in [{np.min(lams):.6g}, {np.max(lams):.6g}]")
    print(f"  global relative spread (max-min)/mean = {spread:.6g}")

    # Per dt: relative spread across n_quad (quadrature sensitivity at fixed dt)
    by_dt: dict[float, list[float]] = defaultdict(list)
    for r in rows:
        by_dt[float(r["dt"])].append(float(r["lambda_max"]))
    quad_spreads: list[float] = []
    for dt in sorted(by_dt):
        v = by_dt[dt]
        if len(v) >= 2:
            m = float(np.mean(v))
            quad_spreads.append((max(v) - min(v)) / m if m > 0 else float("inf"))
    if quad_spreads:
        print(
            "  per-dt spread across n_quad: "
            f"max (max-min)/mean = {max(quad_spreads):.6g}"
        )

    # Per n_quad: range attributable to changing dt only
    by_nq: dict[int, list[float]] = defaultdict(list)
    for r in rows:
        by_nq[int(r["n_quad"])].append(float(r["lambda_max"]))
    for nq in sorted(by_nq):
        v = by_nq[nq]
        if len(v) >= 2:
            m = float(np.mean(v))
            sp = (max(v) - min(v)) / m if m > 0 else float("inf")
            print(f"  n_quad={nq}: (max-min)/mean over dt grid = {sp:.6g}")

    print()
    print("  Interpretation (same seed, analytic J, ex166 setup):")
    print("  - If per-dt n_quad spreads are tiny but dt-spread is large, the")
    print("    time discretization dominates; Feynman quadrature is already")
    print("    close to saturated on this (dt, n_quad) window.")
    print("  - lambda_max is not a converged 'physical' exponent here; it is")
    print("    a finite-time Benettin estimate tied to the trajectory mesh.")


def _edges_midpoints(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return np.array([0.0, 1.0])
    if x.size == 1:
        h = 0.01 * max(abs(x[0]), 1.0)
        return np.array([x[0] - h, x[0] + h])
    inner = 0.5 * (x[:-1] + x[1:])
    left = x[0] - (inner[0] - x[0])
    right = x[-1] + (x[-1] - inner[-1])
    return np.concatenate([[left], inner, [right]])


def _try_plot_png(
    dts: np.ndarray,
    nqs: np.ndarray,
    z: np.ndarray,
    png_path: str,
    t_final: float,
    seed: int,
) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"  (matplotlib not available — PNG skipped: {e})")
        return False

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    dt_e = _edges_midpoints(dts)
    nq_e = _edges_midpoints(nqs.astype(float))
    # Corner grids for pcolormesh: Z shape (ny, nx), X,Y shape (ny+1, nx+1)
    xx = np.broadcast_to(dt_e, (len(nq_e), len(dt_e)))
    yy = np.broadcast_to(nq_e[:, np.newaxis], (len(nq_e), len(dt_e)))
    im = ax0.pcolormesh(xx, yy, z.T, shading="flat")
    ax0.set_xlabel("dt")
    ax0.set_ylabel("n_quad")
    ax0.set_title(f"lambda_max heatmap (t_final={t_final:g}, seed={seed})")
    fig.colorbar(im, ax=ax0, label="lambda_max")

    for j, nq in enumerate(nqs):
        col = z[:, j]
        mask = np.isfinite(col)
        if np.any(mask):
            ax1.plot(dts[mask], col[mask], "o-", label=f"n_quad={nq}")
    ax1.set_xlabel("dt")
    ax1.set_ylabel("lambda_max")
    ax1.set_title("lambda_max vs dt (per n_quad)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(os.path.abspath(png_path)) or ".", exist_ok=True)
    fig.savefig(png_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {png_path}")
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--csv", default="", help=f"default: {_DEFAULT_CSV}")
    parser.add_argument("--png", default="", help=f"default: {_DEFAULT_PNG}")
    args = parser.parse_args()

    csv_path = args.csv or _DEFAULT_CSV
    png_path = args.png or _DEFAULT_PNG

    print("ex168: read ex166 convergence CSV + summary (+ optional PNG)")
    if not os.path.isfile(csv_path):
        print(f"  FAIL: missing {csv_path}")
        print("  Run: python nbody/examples/ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py --quick")
        raise SystemExit(1)

    rows = _load_rows(csv_path)
    if len(rows) < 4:
        print(f"  FAIL: expected at least 4 data rows, got {len(rows)}")
        raise SystemExit(1)

    t0 = float(rows[0]["t_final"])
    sd = int(rows[0]["seed"])
    if not all(float(r["t_final"]) == t0 and int(r["seed"]) == sd for r in rows):
        print("  FAIL: mixed t_final or seed inside CSV (ex168 assumes one experiment)")
        raise SystemExit(1)

    dts, nqs, z = _pivot(rows)
    lams = np.array([float(r["lambda_max"]) for r in rows], dtype=float)
    if not np.all(np.isfinite(lams)) or np.any(lams <= 0.02):
        print("  FAIL: non-finite lambda_max or lambda_max <= 0.02")
        raise SystemExit(1)

    print(f"  csv={csv_path}")
    print(f"  t_final={t0:g}  seed={sd}  unique dt={len(dts)}  unique n_quad={len(nqs)}")
    _print_interpretation(rows, dts, nqs, z)

    spread = (float(np.max(lams)) - float(np.min(lams))) / float(np.mean(lams))
    spread_tol_quick = 0.42
    if args.quick and spread >= spread_tol_quick:
        print(f"  FAIL: global spread {spread:.4f} >= quick tol {spread_tol_quick}")
        raise SystemExit(1)

    if not args.quick:
        _try_plot_png(dts, nqs, z, png_path, t0, sd)

    print("ex168: PASS")


if __name__ == "__main__":
    main()
