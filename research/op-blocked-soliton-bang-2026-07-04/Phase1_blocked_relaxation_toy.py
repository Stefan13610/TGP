#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase1_blocked_relaxation_toy.py
================================

Toy test for the blocked-relaxation bridge:

    single soliton  -> relaxes away
    many solitons   -> overlap blocks the vacuum-relaxation channel

This is intentionally not a GR/MOND solver.  It is a deterministic
2D substrate toy locked by Phase0_balance.md before the first run.

Run:
    python Phase1_blocked_relaxation_toy.py
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import sys
from typing import Dict, Iterable, List, Tuple

import numpy as np


try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass


# Phase0 LOCKED parameters.  Do not tune after seeing output.
N_GRID = 64
L_BOX = 24.0
SEED = 12345
STEPS = 1200
DT = 0.025
KAPPA = 0.18
M2 = 0.55
LAMBDA4 = 0.12
SOURCE_SIGMA = 0.75
BLOCK_SIGMA = 2.40
SOURCE_STRENGTH = 0.55
BLOCK_STRENGTH = 18.0
AMP_RATE = 0.035
DECAY_THRESHOLD = 0.72
ALIVE_THRESHOLD = 0.15
PAIR_DISTANCES = (2.0, 3.0, 4.0, 6.0, 8.0)

# Numerical guard, not a success mechanism.
PHI_MIN = -5.0
PHI_MAX = 8.0
RUNAWAY_LIMIT = 0.98 * PHI_MAX


@dataclass
class Config:
    n_grid: int = N_GRID
    l_box: float = L_BOX
    steps: int = STEPS
    dt: float = DT
    kappa: float = KAPPA
    m2: float = M2
    lambda4: float = LAMBDA4
    source_sigma: float = SOURCE_SIGMA
    block_sigma: float = BLOCK_SIGMA
    source_strength: float = SOURCE_STRENGTH
    block_strength: float = BLOCK_STRENGTH
    amp_rate: float = AMP_RATE
    decay_threshold: float = DECAY_THRESHOLD
    alive_threshold: float = ALIVE_THRESHOLD


@dataclass
class RunResult:
    name: str
    n_initial: int
    n_alive: int
    amp_mean: float
    amp_min: float
    amp_max: float
    e_positive: float
    e_functional: float
    blocked_fraction: float
    mean_overlap_at_solitons: float
    max_phi: float
    min_phi: float
    runaway: bool
    finite: bool


def make_grid(cfg: Config) -> Tuple[np.ndarray, np.ndarray, float]:
    x = np.linspace(-cfg.l_box / 2.0, cfg.l_box / 2.0, cfg.n_grid)
    dx = x[1] - x[0]
    xx, yy = np.meshgrid(x, x, indexing="ij")
    return xx, yy, dx


def laplacian_dirichlet(phi: np.ndarray, dx: float, bc_value: float = 1.0) -> np.ndarray:
    padded = np.pad(phi, 1, mode="constant", constant_values=bc_value)
    return (
        padded[2:, 1:-1]
        + padded[:-2, 1:-1]
        + padded[1:-1, 2:]
        + padded[1:-1, :-2]
        - 4.0 * phi
    ) / (dx * dx)


def source_and_overlap(
    xx: np.ndarray,
    yy: np.ndarray,
    positions: np.ndarray,
    amplitudes: np.ndarray,
    cfg: Config,
) -> Tuple[np.ndarray, np.ndarray]:
    source = np.zeros_like(xx)
    support = np.zeros_like(xx)

    source_norm = 1.0 / (2.0 * math.pi * cfg.source_sigma * cfg.source_sigma)
    for (px, py), amp in zip(positions, amplitudes):
        r2 = (xx - px) ** 2 + (yy - py) ** 2
        source += (
            cfg.source_strength
            * amp
            * source_norm
            * np.exp(-r2 / (2.0 * cfg.source_sigma * cfg.source_sigma))
        )
        support += np.exp(-r2 / (2.0 * cfg.block_sigma * cfg.block_sigma))

    return source, support


def overlap_at_solitons(positions: np.ndarray, cfg: Config) -> np.ndarray:
    n = len(positions)
    if n == 0:
        return np.zeros(0)

    out = np.zeros(n)
    for i in range(n):
        delta = positions[i] - positions
        dist2 = np.sum(delta * delta, axis=1)
        weights = np.exp(-dist2 / (2.0 * cfg.block_sigma * cfg.block_sigma))
        weights[i] = 0.0
        out[i] = float(np.sum(weights))
    return out


def update_amplitudes(
    amplitudes: np.ndarray,
    positions: np.ndarray,
    cfg: Config,
    frozen: bool = False,
) -> np.ndarray:
    if frozen:
        return amplitudes.copy()
    local_overlap = overlap_at_solitons(positions, cfg)
    growth = cfg.amp_rate * amplitudes * (local_overlap - cfg.decay_threshold)
    updated = amplitudes + cfg.dt * growth
    return np.clip(updated, 0.0, 2.5)


def field_energy(
    phi: np.ndarray,
    source: np.ndarray,
    dx: float,
    cfg: Config,
) -> Tuple[float, float]:
    gx, gy = np.gradient(phi, dx, edge_order=1)
    u = phi - 1.0
    e_density = (
        0.5 * cfg.kappa * (gx * gx + gy * gy)
        + 0.5 * cfg.m2 * u * u
        + 0.25 * cfg.lambda4 * u**4
    )
    e_positive = float(np.sum(e_density) * dx * dx)
    e_functional = float(np.sum(e_density - source * u) * dx * dx)
    return e_positive, e_functional


def step_field(
    phi: np.ndarray,
    source: np.ndarray,
    support: np.ndarray,
    dx: float,
    cfg: Config,
) -> Tuple[np.ndarray, np.ndarray]:
    overlap_excess = np.maximum(support - 1.0, 0.0)
    mobility = 1.0 / (1.0 + cfg.block_strength * overlap_excess * overlap_excess)
    u = phi - 1.0
    rhs = (
        cfg.kappa * laplacian_dirichlet(phi, dx)
        - mobility * cfg.m2 * u
        - cfg.lambda4 * u**3
        + source
    )
    new_phi = np.clip(phi + cfg.dt * rhs, PHI_MIN, PHI_MAX)
    return new_phi, mobility


def initial_positions(name: str, cfg: Config) -> np.ndarray:
    if name == "single":
        return np.array([[0.0, 0.0]], dtype=float)

    if name == "cluster":
        radius = 2.35
        angles = np.linspace(0.0, 2.0 * math.pi, 7, endpoint=False)
        ring = np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))
        return np.vstack((np.array([[0.0, 0.0]]), ring))

    raise ValueError(f"unknown scenario: {name}")


def run_dynamic(name: str, cfg: Config, steps: int | None = None) -> RunResult:
    if steps is None:
        steps = cfg.steps

    xx, yy, dx = make_grid(cfg)
    positions = initial_positions(name, cfg)
    amplitudes = np.ones(len(positions), dtype=float)
    phi = np.ones_like(xx)
    mobility = np.ones_like(xx)
    source = np.zeros_like(xx)
    support = np.zeros_like(xx)

    for _ in range(steps):
        source, support = source_and_overlap(xx, yy, positions, amplitudes, cfg)
        phi, mobility = step_field(phi, source, support, dx, cfg)
        amplitudes = update_amplitudes(amplitudes, positions, cfg, frozen=False)

    source, support = source_and_overlap(xx, yy, positions, amplitudes, cfg)
    e_positive, e_functional = field_energy(phi, source, dx, cfg)
    local_overlap = overlap_at_solitons(positions, cfg)
    finite = bool(np.all(np.isfinite(phi)) and np.all(np.isfinite(amplitudes)))
    max_phi = float(np.max(phi))
    min_phi = float(np.min(phi))

    return RunResult(
        name=name,
        n_initial=len(positions),
        n_alive=int(np.sum(amplitudes >= cfg.alive_threshold)),
        amp_mean=float(np.mean(amplitudes)) if len(amplitudes) else 0.0,
        amp_min=float(np.min(amplitudes)) if len(amplitudes) else 0.0,
        amp_max=float(np.max(amplitudes)) if len(amplitudes) else 0.0,
        e_positive=e_positive,
        e_functional=e_functional,
        blocked_fraction=float(np.mean(mobility < 0.5)),
        mean_overlap_at_solitons=float(np.mean(local_overlap)) if len(local_overlap) else 0.0,
        max_phi=max_phi,
        min_phi=min_phi,
        runaway=bool(max_phi >= RUNAWAY_LIMIT or min_phi <= PHI_MIN + 1e-9),
        finite=finite,
    )


def relax_fixed_energy(positions: np.ndarray, cfg: Config, steps: int = 900) -> float:
    xx, yy, dx = make_grid(cfg)
    amplitudes = np.ones(len(positions), dtype=float)
    phi = np.ones_like(xx)

    for _ in range(steps):
        source, support = source_and_overlap(xx, yy, positions, amplitudes, cfg)
        phi, _ = step_field(phi, source, support, dx, cfg)

    source, _ = source_and_overlap(xx, yy, positions, amplitudes, cfg)
    _, e_functional = field_energy(phi, source, dx, cfg)
    return e_functional


def run_pair_scan(cfg: Config) -> List[Tuple[float, float]]:
    single_e = relax_fixed_energy(np.array([[0.0, 0.0]], dtype=float), cfg)
    rows: List[Tuple[float, float]] = []
    for d in PAIR_DISTANCES:
        positions = np.array([[-0.5 * d, 0.0], [0.5 * d, 0.0]], dtype=float)
        two_e = relax_fixed_energy(positions, cfg)
        rows.append((d, two_e - 2.0 * single_e))
    return rows


def pass_fail(name: str, passed: bool, detail: str) -> str:
    return f"[{'OK' if passed else 'FAIL'}] {name} :: {detail}"


def print_result(result: RunResult) -> None:
    print(
        f"{result.name:8s} "
        f"N_alive={result.n_alive:2d}/{result.n_initial:<2d} "
        f"amp_mean={result.amp_mean:.6f} "
        f"amp_range=[{result.amp_min:.6f},{result.amp_max:.6f}] "
        f"E_pos={result.e_positive:.8f} "
        f"E_func={result.e_functional:.8f} "
        f"blocked={result.blocked_fraction:.6f} "
        f"overlap_i={result.mean_overlap_at_solitons:.6f} "
        f"phi=[{result.min_phi:.6f},{result.max_phi:.6f}] "
        f"runaway={result.runaway}"
    )


def main() -> int:
    np.random.default_rng(SEED)  # Locks deterministic environment for future stochastic variants.
    cfg = Config()

    print("=" * 78)
    print("  PHASE 1: BLOCKED RELAXATION TOY")
    print("=" * 78)
    print("LOCK: Phase0_balance.md")
    print(
        "params: "
        f"N={cfg.n_grid}, L={cfg.l_box}, steps={cfg.steps}, dt={cfg.dt}, "
        f"kappa={cfg.kappa}, m2={cfg.m2}, lambda4={cfg.lambda4}, "
        f"source_sigma={cfg.source_sigma}, block_sigma={cfg.block_sigma}, "
        f"block_strength={cfg.block_strength}, amp_rate={cfg.amp_rate}, "
        f"decay_threshold={cfg.decay_threshold}"
    )
    print()

    single = run_dynamic("single", cfg)
    cluster = run_dynamic("cluster", cfg)

    print("DYNAMIC SCENARIOS")
    print_result(single)
    print_result(cluster)
    print()

    print("PAIR SCAN")
    pair_rows = run_pair_scan(cfg)
    for d, v in pair_rows:
        print(f"d={d:4.1f}  V_pair={v:+.10f}")
    v_values = np.array([v for _, v in pair_rows])
    v_span = float(np.max(v_values) - np.min(v_values))
    print(f"V_pair_span={v_span:.10f}")
    print()

    c1 = single.finite and cluster.finite and not single.runaway and not cluster.runaway
    c2 = single.n_alive <= 0 or single.amp_mean < cfg.alive_threshold
    c3 = (
        cluster.n_alive >= 4
        and cluster.amp_mean >= 0.35
        and cluster.blocked_fraction > 0.01
    )
    c4 = v_span > 1e-3

    print("LOCKED CRITERIA")
    print(pass_fail("C1 deterministic finite execution", c1, f"finite={single.finite and cluster.finite}, runaway={single.runaway or cluster.runaway}"))
    print(pass_fail("C2 single soliton decay channel", c2, f"N_alive={single.n_alive}, amp_mean={single.amp_mean:.6f}"))
    print(pass_fail("C3 cluster blocked survival", c3, f"N_alive={cluster.n_alive}, amp_mean={cluster.amp_mean:.6f}, blocked={cluster.blocked_fraction:.6f}"))
    print(pass_fail("C4 nontrivial pair interaction", c4, f"span={v_span:.10f}"))

    n_pass = sum([c1, c2, c3, c4])
    print(f"SUMMARY: PASS={n_pass}/4 FAIL={4 - n_pass}/4")
    print("=" * 78)

    return 0 if c1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
