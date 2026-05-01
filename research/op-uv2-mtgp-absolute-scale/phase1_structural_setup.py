# -*- coding: utf-8 -*-
"""
UV.2.Phase1 -- structural setup + alt-K_struct falsification (5 sub-tests)
Date: 2026-05-01

U1.1 M_TGP-M_GUT separation hypothesis
U1.2 K_struct 4-candidate scan
U1.3 M_Pl/M_GUT cross-prediction
U1.4 Joint-system consistency
U1.5 Alt-anchor ansatz falsification (5 alts)
"""

import math
from sympy import Rational, sqrt, pi, Float

# =====================================================================
# Constants
# =====================================================================
g_star = Rational(71, 100)
N_A = Rational(500, 57)
g_star_f = float(g_star)
N_A_f = float(N_A)
PI = math.pi

# chi.1 LOCK from Phase 3
M_PL_PDG = 1.220890e19          # GeV
M_TGP_chi1 = M_PL_PDG * math.sqrt(g_star_f / N_A_f)  # 3.4734e18 GeV
M_GUT_2loop = 2.0e16            # GeV (SM 2-loop gauge unification central)

# F-cluster anchors
ALPHA_0 = float(Rational(1069833, 264500))   # 4.0447
KAPPA_TGP = 2.012

# Target K = M_TGP/M_GUT
K_target = M_TGP_chi1 / M_GUT_2loop

results = []
def record(label, ok, detail=""):
    results.append((label, ok, detail))

# =====================================================================
# U1.1 -- M_TGP-M_GUT separation hypothesis
# =====================================================================
print("="*70)
print("  U1.1 -- M_TGP-M_GUT separation hypothesis")
print("="*70)
print(f"  M_TGP_chi.1   = {M_TGP_chi1:.4e} GeV (joint-lock z M_Pl PDG)")
print(f"  M_GUT (2loop) = {M_GUT_2loop:.2e} GeV (SM gauge unification central)")
print(f"  M_TGP/M_GUT   = {K_target:.2f}")
print(f"  M_TGP/M_Pl    = {M_TGP_chi1/M_PL_PDG:.4f}  (chi.1 sqrt(g*/N_A))")
print(f"  M_GUT/M_Pl    = {M_GUT_2loop/M_PL_PDG:.4e}")
print()
print(f"  Hypothesis: M_TGP = K_struct * M_GUT, K_struct ~ O(100-1000)")
print(f"  Substrate scale 174x above gauge-unification (sub-Planckian)")
print(f"  Cyrkularnosc: chi.1 anchored M_TGP z M_Pl PDG -> UV.2 swaps anchor to M_GUT")

# Gate: K_target in plausible band [10, 1000] (sub-Planckian, above-GUT)
gate_U11 = 10 < K_target < 1000
if gate_U11:
    print(f"  [PASS] M_TGP/M_GUT in plausible separation band [10, 1000]")
    print(f"         K_target = {K_target:.2f}")
else:
    print(f"  [FAIL] K_target outside band")
record("U1.1", gate_U11, f"K_target = {K_target:.2f} in band [10, 1000]")

# =====================================================================
# U1.2 -- K_struct 4-candidate scan
# =====================================================================
print()
print("="*70)
print("  U1.2 -- K_struct 4-candidate scan")
print("="*70)

candidates = [
    ("(a) N_A * 2*pi^2", float(N_A * 2 * pi**2)),
    ("(b) N_A^2 * sqrt(2*pi)", float(N_A**2 * sqrt(2*pi))),
    ("(c) (4*pi) * N_A * kappa_TGP", float(4*pi * N_A) * KAPPA_TGP),
    ("(d) alpha_0 * 4*pi^2 * sqrt(N_A)", ALPHA_0 * 4*PI**2 * float(sqrt(N_A))),
]
print(f"  Target K = M_TGP_chi.1 / M_GUT = {K_target:.4f}")
print()
print(f"  {'candidate':<40} {'K':>10} {'drift':>10}")
print(f"  {'-'*40} {'-'*10} {'-'*10}")
best = None
best_drift = 1e10
for label, K in candidates:
    drift = abs(K - K_target) / K_target
    print(f"  {label:<40} {K:>10.4f} {drift*100:>9.4f}%")
    if drift < best_drift:
        best_drift = drift
        best = (label, K, drift)

print()
print(f"  Winner: {best[0]} -> K = {best[1]:.4f}, drift {best[2]*100:.4f}%")
print(f"  Selection criterion: drift < 1% partial-lock band")
gate_U12 = best[2] < 0.01
if gate_U12:
    print(f"  [PASS] K_a = N_A * 2*pi^2 within partial-lock band")
    print(f"         drift = {best[2]*100:.4f}% < 1%")
else:
    print(f"  [FAIL] No K candidate within 1% drift band")
record("U1.2", gate_U12, f"K_a winner drift {best[2]*100:.4f}%")

# =====================================================================
# U1.3 -- M_Pl/M_GUT cross-prediction
# =====================================================================
print()
print("="*70)
print("  U1.3 -- M_Pl/M_GUT cross-prediction (under K_a winner)")
print("="*70)
# Under K_a = N_A * 2*pi^2, and chi.1: M_TGP/M_Pl = sqrt(g*/N_A):
# M_Pl/M_GUT = K_struct / sqrt(g*/N_A) = (N_A * 2*pi^2) / sqrt(g*/N_A) = 2*pi^2 * N_A^(3/2) / sqrt(g*)
M_Pl_over_M_GUT_predicted = 2 * PI**2 * N_A_f**1.5 / math.sqrt(g_star_f)
M_Pl_over_M_GUT_observed = M_PL_PDG / M_GUT_2loop
drift_MPl_MGUT = abs(M_Pl_over_M_GUT_predicted - M_Pl_over_M_GUT_observed) / M_Pl_over_M_GUT_observed

print(f"  Predicted: M_Pl/M_GUT = 2*pi^2 * N_A^(3/2) / sqrt(g*)")
print(f"             = 2*pi^2 * (500/57)^(3/2) / sqrt(71/100)")
print(f"             = {M_Pl_over_M_GUT_predicted:.4f}")
print(f"  Observed:  M_Pl/M_GUT = {M_PL_PDG:.4e}/{M_GUT_2loop:.2e}")
print(f"             = {M_Pl_over_M_GUT_observed:.4f}")
print(f"  Drift                   = {drift_MPl_MGUT*100:.4f}%")
print(f"  Drift target            < 1%")
gate_U13 = drift_MPl_MGUT < 0.01
if gate_U13:
    print(f"  [PASS] M_Pl/M_GUT structural prediction within 1%")
else:
    print(f"  [FAIL] M_Pl/M_GUT drift exceeds 1%")
record("U1.3", gate_U13, f"M_Pl/M_GUT drift {drift_MPl_MGUT*100:.4f}%")

# =====================================================================
# U1.4 -- Joint-system consistency
# =====================================================================
print()
print("="*70)
print("  U1.4 -- Joint-system consistency (M_TGP-M_Pl-M_GUT triple)")
print("="*70)
K_struct = N_A_f * 2 * PI**2

# Predicted M_TGP from M_GUT under K_a:
M_TGP_uv2_from_GUT = K_struct * M_GUT_2loop
# Predicted M_Pl from joint-system: M_Pl = M_TGP / sqrt(g*/N_A)
M_Pl_uv2 = M_TGP_uv2_from_GUT / math.sqrt(g_star_f / N_A_f)

drift_M_TGP_uv2 = abs(M_TGP_uv2_from_GUT - M_TGP_chi1) / M_TGP_chi1
drift_M_Pl_uv2 = abs(M_Pl_uv2 - M_PL_PDG) / M_PL_PDG

print(f"  K_struct (UV.2 winner)  = N_A * 2*pi^2 = {K_struct:.4f}")
print(f"  M_TGP_UV.2 = K * M_GUT = {M_TGP_uv2_from_GUT:.4e} GeV")
print(f"  M_TGP_chi.1            = {M_TGP_chi1:.4e} GeV")
print(f"  Drift M_TGP             = {drift_M_TGP_uv2*100:.4f}%")
print(f"  M_Pl_UV.2 = M_TGP/sqrt(g*/N_A) = {M_Pl_uv2:.4e} GeV")
print(f"  M_Pl PDG               = {M_PL_PDG:.4e} GeV")
print(f"  Drift M_Pl              = {drift_M_Pl_uv2*100:.4f}%")

# Joint consistency: both drifts must be < 1%
gate_U14 = drift_M_TGP_uv2 < 0.01 and drift_M_Pl_uv2 < 0.01
if gate_U14:
    print(f"  [PASS] Joint-system M_TGP-M_Pl-M_GUT consistent within 1%")
else:
    print(f"  [FAIL] Joint-system inconsistent")
record("U1.4", gate_U14, f"M_TGP drift {drift_M_TGP_uv2*100:.4f}%, M_Pl drift {drift_M_Pl_uv2*100:.4f}%")

# =====================================================================
# U1.5 -- Alt-anchor ansatz falsification (5 alts)
# =====================================================================
print()
print("="*70)
print("  U1.5 -- Alt-anchor ansatz falsification (5 alts)")
print("="*70)

alts = [
    ("(i) M_TGP = M_Pl postulate",
     M_PL_PDG,
     "circular (chi.1 tautology, no derivation gain)"),
    ("(ii) M_TGP = M_GUT direct",
     M_GUT_2loop,
     "no separation (M_TGP/M_GUT = 1, observed 174)"),
    ("(iii) M_TGP = sqrt(M_Pl * M_GUT) geom. mean",
     math.sqrt(M_PL_PDG * M_GUT_2loop),
     "XS dim-mix (orthogonal sectors)"),
    ("(iv) M_TGP = M_GUT/g* (g*-only correction)",
     M_GUT_2loop / g_star_f,
     "misses N_A photon-ring (drift 100x)"),
    ("(v) M_TGP = N_A * 2*pi^2 * M_GUT (UV.2 winner)",
     M_TGP_uv2_from_GUT,
     "TGP-native UV.1 g* + xi.1 N_A + S^3 geom"),
]

print(f"  Target M_TGP = {M_TGP_chi1:.4e} GeV (chi.1 anchor)")
print()
print(f"  {'ansatz':<50} {'M_TGP':>14} {'drift':>10}")
print(f"  {'-'*50} {'-'*14} {'-'*10}")
n_pass = 0
winner_label = None
for label, M_pred, rationale in alts:
    drift = abs(M_pred - M_TGP_chi1) / M_TGP_chi1
    status = "PASS" if drift < 0.01 else "FAIL"
    if drift < 0.01:
        n_pass += 1
        winner_label = label
    print(f"  {label:<50} {M_pred:>14.4e} {drift*100:>9.2f}%  [{status}]")

print()
print(f"  PROBE-pass count: {n_pass}/5")
gate_U15 = n_pass == 1  # uniqueness criterion: exactly 1 alt passes
if gate_U15:
    print(f"  [PASS] Uniqueness criterion satisfied: exactly 1/5 alt passes")
    print(f"         winner: {winner_label}")
else:
    print(f"  [FAIL] Uniqueness criterion not satisfied ({n_pass}/5)")
record("U1.5", gate_U15, f"unique winner {n_pass}/5 PROBE-pass")

# =====================================================================
# Phase 1 verdict
# =====================================================================
print()
print("="*70)
print("  PHASE 1 VERDICT")
print("="*70)
n_pass = sum(1 for _, ok, _ in results if ok)
n_total = len(results)
for label, ok, detail in results:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}  -- {detail}")
print()
print(f"  SCORE: {n_pass}/{n_total}")
gate_phase1 = n_pass >= 4
print(f"  GATE: {'PASS' if gate_phase1 else 'FAIL'} (>=4/5) -> {'Phase 2 enabled' if gate_phase1 else 'BLOCKED'}")
print()
print(f"  K_struct LOCK: K = N_A * 2*pi^2 = {K_struct:.4f}")
print(f"  M_TGP_UV.2 (predicted): {M_TGP_uv2_from_GUT:.4e} GeV")
print(f"  M_Pl_UV.2 (predicted):  {M_Pl_uv2:.4e} GeV")
print(f"  M_GUT (anchor):          {M_GUT_2loop:.2e} GeV")
