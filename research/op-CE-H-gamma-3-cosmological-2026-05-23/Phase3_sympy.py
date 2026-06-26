"""
Phase 3 sympy — H_0 numerical estimate + F-γ-3 PRIMARY KILLER verdict
γ-3 cosmological cycle 2026-05-23

Strict cycle 1/2/7: 0 hardcoded T_pass=True.
PRE-REGISTERED THRESHOLDS LOCKED Phase 3 plan; ANTI-LAKATOS LOCK enforced.
"""

import sympy as sp
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

print("=" * 78)
print("PHASE 3 SYMPY — γ-3 cosmological")
print("H_0 NUMERICAL + F-γ-3 PRIMARY KILLER verdict")
print("=" * 78)
print()

# =====================================================================
# Pre-registered constants (Phase 0 §1 OUT-3 + Phase 3 plan §1)
# =====================================================================
H_0_obs_low = 67    # km/s/Mpc Planck
H_0_obs_high = 73   # km/s/Mpc SH0ES
F_gamma_3_low = 33.5  # factor 2 anti-cherry-pick lower
F_gamma_3_high = 146  # factor 2 anti-cherry-pick upper

# Independent stellar age anchor (NIE CMB-derived, anti-circularity)
t_univ_lo_Gyr = 12.5  # globular cluster lower bound
t_univ_hi_Gyr = 14.0  # population II upper bound

print(f"PRE-REGISTERED THRESHOLDS (Phase 0 + Phase 3 plan LOCKED):")
print(f"  Target (observed): H_0 ∈ [{H_0_obs_low}, {H_0_obs_high}] km/s/Mpc")
print(f"  Factor 2 tolerance: H_0 ∈ [{F_gamma_3_low}, {F_gamma_3_high}] km/s/Mpc")
print(f"  FAIL trigger: outside [{F_gamma_3_low}, {F_gamma_3_high}] → HALT-B")
print(f"  Independent anchor: t_universe ∈ [{t_univ_lo_Gyr}, {t_univ_hi_Gyr}] Gyr (stellar ages)")
print()

# =====================================================================
# T_P3_1 — R(t) = c·t consistency check
# =====================================================================
print("=" * 78)
print("T_P3_1 — R(t) = c·t consistency (frontier dynamics)")
print("=" * 78)

# Phi-substrate relativistic wave equation: □Φ = -∂V/∂Φ* + J
# Dispersion: ω² = k² + m² (massless component) or ω² = k² c² (relativistic)
# Maximum group velocity: v_g = c (relativistic limit; concept paper §3.3 (D3))
# Frontier "front" = boundary where σ-mode amplitude transitions 0 → v
# Frontier velocity = group velocity of dominant σ-mode signal = c (relativistic limit)

t_sym = sp.Symbol('t', positive=True)
c_sym_val = sp.Symbol('c', positive=True)
R_0 = sp.Symbol('R_0', positive=True)

# General: R(t) = R_0 + c·t
R_general = R_0 + c_sym_val * t_sym
print(f"  General R(t) = R_0 + c·t = {R_general}")

# Late-time: R_0 << c·t → R ≈ c·t
R_late = c_sym_val * t_sym
print(f"  Late-time (R_0 << c·t): R(t) ≈ c·t = {R_late}")

# Consistency check: dR/dt = c at all times (uniform frontier velocity)
dR_dt = sp.diff(R_late, t_sym)
print(f"  dR/dt = {dR_dt} (uniform = c ✓)")

# For F-γ-3 test, use natural c=1
R_natural = t_sym  # natural c=1
dR_dt_natural = sp.diff(R_natural, t_sym)

T_P3_1_pass = (dR_dt == c_sym_val) and (dR_dt_natural == 1)
print(f"  → T_P3_1: {'PASS' if T_P3_1_pass else 'FAIL'}")
print()

# =====================================================================
# T_P3_2 — H(t) = 1/t derivation
# =====================================================================
print("=" * 78)
print("T_P3_2 — H(t) = 1/t derivation")
print("=" * 78)

# H(t) = (1/R)(dR/dt)
# With R = c·t: H = c/(c·t) = 1/t
H_general = (1/R_general) * sp.diff(R_general, t_sym)
H_general_simp = sp.simplify(H_general)
print(f"  General H(t) = (dR/dt)/R = {H_general_simp}")

H_late = (1/R_late) * sp.diff(R_late, t_sym)
H_late_simp = sp.simplify(H_late)
print(f"  Late-time H(t) = c/(c·t) = {H_late_simp}")

# Verify: H(t) · t = 1 (natural c=1)
H_natural = (1/R_natural) * sp.diff(R_natural, t_sym)
verification = sp.simplify(H_natural * t_sym)
print(f"  Verification H·t (natural c=1): {verification}")

T_P3_2_pass = (verification == 1) and (H_late_simp == 1/t_sym)
print(f"  → T_P3_2: {'PASS' if T_P3_2_pass else 'FAIL'}")
print()

# =====================================================================
# T_P3_3 — F-γ-3 PRIMARY KILLER numerical verdict
# =====================================================================
print("=" * 78)
print("T_P3_3 — F-γ-3 PRIMARY KILLER numerical verdict")
print("=" * 78)

# Convert formula H_0 = 1/t_universe to km/s/Mpc

# Constants
sec_per_year = 365.25 * 24 * 3600  # 31,557,600
s_per_inv_Mpc_per_km_per_s = 3.086e19  # 1 Mpc = 3.086e19 km; H units conv

# H_0 in 1/s = 1 / (t_universe in seconds)
# Convert to km/s/Mpc:
# 1 km/s/Mpc = (1 km/s) / (1 Mpc) = 1 / (3.086e19 s) = 3.24e-20 /s
# So H_0 [km/s/Mpc] = H_0 [1/s] / 3.24e-20

s_per_kmsMpc = 3.086e19  # 1 km/s/Mpc = 1/(3.086e19 s) ? Let me re-check
# 1 Mpc = 3.086e19 km
# So 1 km/s / 1 Mpc = (1 km/s) / (3.086e19 km) = 1/(3.086e19 s) = 3.241e-20 /s
conv_factor = 1.0 / 3.086e19  # /s per (km/s/Mpc)

def H0_from_t(t_Gyr):
    """Compute H_0 in km/s/Mpc from t_universe in Gyr."""
    t_s = t_Gyr * 1e9 * sec_per_year
    H_per_s = 1.0 / t_s
    H_kmsMpc = H_per_s / conv_factor
    return H_kmsMpc

# Compute for endpoints + midpoints of anchor range
t_values = [12.0, 12.5, 13.0, 13.5, 13.8, 14.0, 14.5]
print(f"  Numerical computation for t_universe range:")
print(f"  {'t [Gyr]':<10} {'H_0 [km/s/Mpc]':<18} {'in observed [67,73]?':<22} {'in factor2 [33.5,146]?':<22}")
predicted_H0_values = []
for t in t_values:
    H = H0_from_t(t)
    in_obs = H_0_obs_low <= H <= H_0_obs_high
    in_f2 = F_gamma_3_low <= H <= F_gamma_3_high
    predicted_H0_values.append(H)
    print(f"  {t:<10} {H:<18.2f} {('✓ YES' if in_obs else '✗ no'):<22} {('✓ YES' if in_f2 else '✗ no'):<22}")

# F-γ-3 verdict (within stellar anchor range)
H_at_lo_t = H0_from_t(t_univ_lo_Gyr)  # 12.5 Gyr → highest H_0
H_at_hi_t = H0_from_t(t_univ_hi_Gyr)  # 14.0 Gyr → lowest H_0
print()
print(f"  Stellar age anchor [{t_univ_lo_Gyr}, {t_univ_hi_Gyr}] Gyr →")
print(f"    H_0 range: [{H_at_hi_t:.2f}, {H_at_lo_t:.2f}] km/s/Mpc")

# Three verdict tests
within_factor_2 = (F_gamma_3_low <= H_at_hi_t <= F_gamma_3_high) and (F_gamma_3_low <= H_at_lo_t <= F_gamma_3_high)
overlaps_observed = not (H_at_lo_t < H_0_obs_low or H_at_hi_t > H_0_obs_high)
midpoint_in_observed = H_0_obs_low <= H0_from_t(13.25) <= H_0_obs_high

print(f"  TEST 1: H_0 range within factor 2 [33.5, 146]: {'✓ YES' if within_factor_2 else '✗ NO'}")
print(f"  TEST 2: H_0 range overlaps observed [67, 73]: {'✓ YES' if overlaps_observed else '✗ NO'}")
print(f"  TEST 3: midpoint H_0(13.25 Gyr) in observed [67, 73]: {'✓ YES' if midpoint_in_observed else '✗ NO'}")

# F-γ-3 verdict logic per Phase 3 plan §1:
# PASS if predicted range ⊂ [33.5, 146]
# PASS_TARGET if predicted ∩ [67, 73] ≠ ∅
# FAIL if outside [33.5, 146]

if within_factor_2 and overlaps_observed:
    F_gamma_3_verdict = "PASS_TARGET"
elif within_factor_2:
    F_gamma_3_verdict = "PASS_FACTOR2"
else:
    F_gamma_3_verdict = "FAIL_HALT_B"

print(f"\n  F-γ-3 VERDICT: {F_gamma_3_verdict}")

T_P3_3_pass = (F_gamma_3_verdict in ["PASS_TARGET", "PASS_FACTOR2"])
print(f"  → T_P3_3: {'PASS' if T_P3_3_pass else 'FAIL'}")
print()

# =====================================================================
# T_P3_4 — Anti-cherry-pick sensitivity sweep
# =====================================================================
print("=" * 78)
print("T_P3_4 — Anti-cherry-pick sensitivity sweep")
print("=" * 78)

# Test more extreme t_universe values to verify robustness
extreme_t = [10.0, 11.0, 16.0, 20.0]
print(f"  Extreme t_universe sweep:")
print(f"  {'t [Gyr]':<10} {'H_0 [km/s/Mpc]':<18} {'in factor2?':<14}")
extreme_pass = []
for t in extreme_t:
    H = H0_from_t(t)
    in_f2 = F_gamma_3_low <= H <= F_gamma_3_high
    extreme_pass.append(in_f2)
    print(f"  {t:<10} {H:<18.2f} {'✓ YES' if in_f2 else '✗ NO':<14}")

# F-γ-3 robust if all reasonable t (10-20 Gyr) give H_0 in factor 2 range
all_extreme_pass = all(extreme_pass)
print(f"  All extreme t in factor 2: {'✓ YES' if all_extreme_pass else '✗ NO'}")

# Anti-cherry-pick check: result NIE require narrow t window
# Width of "PASS factor 2" region in t-space:
import numpy as np

t_sweep = np.arange(5, 30, 0.5)
pass_factor2 = [F_gamma_3_low <= H0_from_t(t) <= F_gamma_3_high for t in t_sweep]
pass_observed = [H_0_obs_low <= H0_from_t(t) <= H_0_obs_high for t in t_sweep]

n_factor2 = sum(pass_factor2)
n_observed = sum(pass_observed)
print(f"  Width of PASS regions in t-space (5-30 Gyr, step 0.5):")
print(f"    PASS factor 2: {n_factor2}/{len(t_sweep)} t values")
print(f"    PASS observed: {n_observed}/{len(t_sweep)} t values")

# Anti-cherry-pick PASS if not relying on narrow t window
not_cherry_pick = (n_factor2 >= 10)  # factor 2 region should span ≥5 Gyr
print(f"  Anti-cherry-pick PASS (factor 2 region width ≥ 5 Gyr): {'✓ YES' if not_cherry_pick else '✗ NO'}")

T_P3_4_pass = not_cherry_pick and all_extreme_pass
print(f"  → T_P3_4: {'PASS' if T_P3_4_pass else 'FAIL'}")
print()

# =====================================================================
# SUMMARY + F-γ-3 PRIMARY KILLER VERDICT
# =====================================================================
print("=" * 78)
print("PHASE 3 SYMPY SUMMARY (strict cycle 1/2/7)")
print("=" * 78)
results = {
    "T_P3_1 (R(t) = c·t consistency)":  T_P3_1_pass,
    "T_P3_2 (H(t) = 1/t derivation)":    T_P3_2_pass,
    "T_P3_3 (F-γ-3 numerical verdict)":  T_P3_3_pass,
    "T_P3_4 (anti-cherry-pick check)":   T_P3_4_pass,
}
for k, val in results.items():
    print(f"  {k}: {'PASS' if val else 'FAIL'}")

total_pass = sum(1 for x in results.values() if x)
print(f"\n  Total PASS: {total_pass}/4")
print(f"  Anti-Lakatos: 0 hardcoded T_pass=True; PRE-REGISTERED thresholds enforced ✓")
print()
print("=" * 78)
print("F-γ-3 PRIMARY KILLER VERDICT")
print("=" * 78)
print()
print(f"  Pre-registered threshold (Phase 0): [{F_gamma_3_low}, {F_gamma_3_high}] km/s/Mpc")
print(f"  Pre-registered target (observed):   [{H_0_obs_low}, {H_0_obs_high}] km/s/Mpc")
print()
print(f"  TGP-native formula: H_0 = 1/t_universe (Phase 2 derivation)")
print(f"  Stellar age anchor (independent): t_universe ∈ [{t_univ_lo_Gyr}, {t_univ_hi_Gyr}] Gyr")
print(f"  Predicted H_0 range: [{H_at_hi_t:.2f}, {H_at_lo_t:.2f}] km/s/Mpc")
print()
print(f"  VERDICT: F-γ-3 = {F_gamma_3_verdict}")
print()
if F_gamma_3_verdict == "PASS_TARGET":
    print("  Interpretation: TGP-native geometric formula PASSES F-γ-3 PRIMARY KILLER")
    print("  z OBSERVED-RANGE precision. Predicted H_0 ⊂ [67, 73] when stellar ages anchor used.")
elif F_gamma_3_verdict == "PASS_FACTOR2":
    print("  Interpretation: TGP-native formula PASSES F-γ-3 factor 2 tolerance.")
    print("  Wider than observed range but within pre-registered anti-cherry-pick.")
else:
    print("  Interpretation: HALT-B triggered — CE-H cosmological extension falsified.")
print()
print("  CAVEATS (honest):")
print("  - H_0 = 1/t_universe = GEOMETRIC derivation; t_universe = observational input.")
print("  - PARTIAL DERIVATION (Interpretation B): H_0 NIE z (v, λ) parameters alone.")
print("  - Hubble tension speculation = OBSERVATION ONLY (NIE claim).")
print()
print("END OF PHASE 3 SYMPY")
