#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex142_yukawa_overlap_quadrature_convergence.py
==============================================
Zbieżność kwadratury ``n_quad`` dla ``yukawa_overlap_exact`` (jeden triplet,
geometria równoboczna).  Referencja: ``n_quad=56``.

Sens: przy wyborze ``n_quad_feynman`` w ``build_tgp_integration_pair`` widać
kompromis koszt / błąd (warstwa numeryczna, nie analityka EOM).

PASS: względna różnica przy ``n_quad=32`` względem referencji < 1e-6.
"""

from __future__ import annotations

import os
import sys

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.three_body_force_exact import yukawa_overlap_exact


def main() -> None:
    d = 1.55
    m = 1.0
    ref_n = 56
    I_ref = yukawa_overlap_exact(d, d, d, m, n_quad=ref_n)
    print("ex142: yukawa_overlap_exact(d,d,d), m=1, d=1.55")
    print(f"  reference n_quad={ref_n}  I={I_ref:.12f}")
    rows = [12, 16, 20, 24, 28, 32, 40, 48]
    for nq in rows:
        I = yukawa_overlap_exact(d, d, d, m, n_quad=nq)
        rel = abs(I - I_ref) / abs(I_ref)
        print(f"  n_quad={nq:2d}  rel_err={rel:.3e}")
    I32 = yukawa_overlap_exact(d, d, d, m, n_quad=32)
    rel32 = abs(I32 - I_ref) / abs(I_ref)
    ok = rel32 < 1e-6
    print(f"  criterion rel(32 vs {ref_n}) < 1e-6 -> {rel32:.3e}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
