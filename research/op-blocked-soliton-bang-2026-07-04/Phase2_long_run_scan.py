#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase2_long_run_scan.py
=======================

Exploratory follow-up to Phase 1.  This script does not change the
Phase1 LOCK criteria.  It checks whether the C2 failure at 1200 steps
is a time-scale issue and whether cluster survival depends on the
field-level blocking term or only on the amplitude overlap rule.

Default run is short enough for routine reproduction.  Use:

    python Phase2_long_run_scan.py --profile full

for a larger scan suitable for a longer unattended run.
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from dataclasses import replace
from typing import Iterable, List, Tuple

import numpy as np

import Phase1_blocked_relaxation_toy as phase1


try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass


def c2_single_decay(result: phase1.RunResult, cfg: phase1.Config) -> bool:
    return result.n_alive <= 0 or result.amp_mean < cfg.alive_threshold


def c3_cluster_survival(result: phase1.RunResult) -> bool:
    return (
        result.n_alive >= 4
        and result.amp_mean >= 0.35
        and result.blocked_fraction > 0.01
    )


def print_time_scan(cfg: phase1.Config, steps_list: Iterable[int]) -> None:
    print("TIME SCAN: default parameters")
    print("steps scenario N_alive amp_mean blocked E_func verdict")
    for steps in steps_list:
        single = phase1.run_dynamic("single", cfg, steps=steps)
        cluster = phase1.run_dynamic("cluster", cfg, steps=steps)
        print(
            f"{steps:5d} single  {single.n_alive:2d}/{single.n_initial:<2d} "
            f"{single.amp_mean:.6f} {single.blocked_fraction:.6f} "
            f"{single.e_functional:+.8f} C2={'PASS' if c2_single_decay(single, cfg) else 'FAIL'}"
        )
        print(
            f"{steps:5d} cluster {cluster.n_alive:2d}/{cluster.n_initial:<2d} "
            f"{cluster.amp_mean:.6f} {cluster.blocked_fraction:.6f} "
            f"{cluster.e_functional:+.8f} C3={'PASS' if c3_cluster_survival(cluster) else 'FAIL'}"
        )
    print()


def print_block_scan(
    cfg: phase1.Config,
    block_strengths: Iterable[float],
    decay_thresholds: Iterable[float],
    steps: int,
) -> None:
    print("BLOCK / DECAY SCAN")
    print("steps block_strength decay_threshold scenario N_alive amp_mean blocked E_func status")
    for block_strength, decay_threshold in itertools.product(block_strengths, decay_thresholds):
        local_cfg = replace(
            cfg,
            block_strength=block_strength,
            decay_threshold=decay_threshold,
            steps=steps,
        )
        single = phase1.run_dynamic("single", local_cfg, steps=steps)
        cluster = phase1.run_dynamic("cluster", local_cfg, steps=steps)
        print(
            f"{steps:5d} {block_strength:14.3f} {decay_threshold:15.3f} "
            f"single  {single.n_alive:2d}/{single.n_initial:<2d} "
            f"{single.amp_mean:.6f} {single.blocked_fraction:.6f} "
            f"{single.e_functional:+.8f} C2={'PASS' if c2_single_decay(single, local_cfg) else 'FAIL'}"
        )
        print(
            f"{steps:5d} {block_strength:14.3f} {decay_threshold:15.3f} "
            f"cluster {cluster.n_alive:2d}/{cluster.n_initial:<2d} "
            f"{cluster.amp_mean:.6f} {cluster.blocked_fraction:.6f} "
            f"{cluster.e_functional:+.8f} C3={'PASS' if c3_cluster_survival(cluster) else 'FAIL'}"
        )
    print()


def print_pair_block_scan(cfg: phase1.Config, block_strengths: Iterable[float]) -> None:
    print("PAIR ENERGY VS FIELD BLOCKING")
    print("block_strength V_span V_values")
    for block_strength in block_strengths:
        local_cfg = replace(cfg, block_strength=block_strength)
        rows = phase1.run_pair_scan(local_cfg)
        values = np.array([v for _, v in rows])
        span = float(np.max(values) - np.min(values))
        value_text = ", ".join(f"d={d:.1f}:{v:+.6f}" for d, v in rows)
        print(f"{block_strength:14.3f} {span:.10f} {value_text}")
    print()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile",
        choices=("default", "full"),
        default="default",
        help="default is reproducible and quick; full is for longer unattended scans",
    )
    args = parser.parse_args()

    started = time.time()
    cfg = phase1.Config()

    if args.profile == "full":
        steps_list = [1200, 2400, 3600, 4800, 7200]
        block_strengths = [0.0, 3.0, 6.0, 18.0, 36.0, 72.0]
        decay_thresholds = [0.55, 0.72, 0.95, 1.20]
        block_steps = 3600
    else:
        steps_list = [1200, 2400, 3600]
        block_strengths = [0.0, 6.0, 18.0, 36.0]
        decay_thresholds = [0.55, 0.72, 0.95]
        block_steps = 2400

    print("=" * 78)
    print("  PHASE 2: LONG-RUN / PARAMETER SCAN")
    print("=" * 78)
    print(f"profile={args.profile}")
    print("NOTE: Phase2 is exploratory; it does not rewrite Phase1 C1-C4.")
    print()

    print_time_scan(cfg, steps_list)
    print_block_scan(cfg, block_strengths, decay_thresholds, block_steps)
    print_pair_block_scan(cfg, block_strengths)

    elapsed = time.time() - started
    print("INTERPRETATION FLAGS")
    print("- If single passes C2 only at steps > 1200, Phase1 C2 failure is a locked time-scale miss.")
    print("- If block_strength=0 keeps cluster amplitudes but C3 fails, amplitude-overlap and field-blocking are distinct.")
    print("- If V_span changes with block_strength, pair energy is sensitive to the field-level blocked-relaxation channel.")
    print(f"elapsed_seconds={elapsed:.3f}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
