"""
η.1.Phase1 — Wolfenstein (A, ρ̄, η̄) numerical landscape audit
=============================================================

5 sub-tests:
  T1.1: A = 0.790 ± 0.012 PDG + top-5 rational candidates
  T1.2: ρ̄ = 0.141 ± 0.020 PDG + top-5 rational candidates
  T1.3: η̄ = 0.357 ± 0.014 PDG + top-5 rational candidates (5/14 hypothesis)
  T1.4: Unitarity triangle apex closure α + β + γ = π
  T1.5: Jarlskog J = A²·λ_C⁶·η̄ vs PDG J

Inputs:
  PDG 2024 Wolfenstein (A, λ, ρ̄, η̄) = (0.790 ± 0.012, 0.22650 ± 0.00048,
                                         0.141 ± 0.020, 0.357 ± 0.014)
  Measured J = (3.07 ± 0.10) · 10⁻⁵
  ζ.1/θ.1 single-anchor λ_C = 0.22550

Sympy 30-digit precision throughout.
"""

import sympy as sp

# PDG 2024 Wolfenstein parameters
A_PDG = sp.Float("0.790", 30)
A_PDG_ERR = sp.Float("0.012", 30)
LAMBDA_PDG = sp.Float("0.22650", 30)
LAMBDA_PDG_ERR = sp.Float("0.00048", 30)
RHO_PDG = sp.Float("0.141", 30)
RHO_PDG_ERR = sp.Float("0.020", 30)
ETA_PDG = sp.Float("0.357", 30)
ETA_PDG_ERR = sp.Float("0.014", 30)

# ζ.1/θ.1 cross-sector single Cabibbo anchor
LAMBDA_C = sp.Float("0.22550", 30)

# PDG 2024 Jarlskog J
J_PDG = sp.Float("3.07e-5", 30)
J_PDG_ERR = sp.Float("0.10e-5", 30)

PI = sp.N(sp.pi, 30)


def drift_pct(value, target):
    """Percent drift of value vs target."""
    expr = sp.Abs(value - target) / sp.Abs(target) * 100
    return sp.Float(sp.N(expr, 30), 30)


def best_rationals(target, max_denom, top_n=5):
    """Return top_n rationals p/q (q ≤ max_denom) closest to target."""
    target_f = float(target)
    candidates = []
    for q in range(2, max_denom + 1):
        p = round(target_f * q)
        if p < 1:
            continue
        r = sp.Rational(p, q)
        d = drift_pct(sp.Float(r, 30), target)
        candidates.append((r, d, p, q))
    # Deduplicate same value
    seen = set()
    unique = []
    for r, d, p, q in sorted(candidates, key=lambda x: float(x[1])):
        key = (r.p, r.q)
        if key in seen:
            continue
        seen.add(key)
        unique.append((r, d, p, q))
        if len(unique) >= top_n:
            break
    return unique


print("=" * 78)
print("η.1.Phase1 — Wolfenstein (A, ρ̄, η̄) numerical landscape audit")
print("=" * 78)
print()
print("Inputs (PDG 2024):")
print(f"  A      = {A_PDG} ± {A_PDG_ERR}")
print(f"  λ      = {LAMBDA_PDG} ± {LAMBDA_PDG_ERR}")
print(f"  λ_C    = {LAMBDA_C}  (ζ.1/θ.1 single-anchor)")
print(f"  ρ̄      = {RHO_PDG} ± {RHO_PDG_ERR}")
print(f"  η̄      = {ETA_PDG} ± {ETA_PDG_ERR}")
print(f"  J_PDG  = {J_PDG} ± {J_PDG_ERR}")
print()

results = {}

# ------------------------------------------------------------------
# T1.1 — A = 0.790 PDG + top-5 rational candidates
# ------------------------------------------------------------------
print("-" * 78)
print("T1.1 — A = 0.790 PDG + top-5 rational candidates (denom ≤ 100)")
print("-" * 78)
top_A = best_rationals(A_PDG, 100, top_n=5)
print(f"  Target: A_PDG = {A_PDG}")
print(f"  Top-5 rational candidates (denom ≤ 100):")
for i, (r, d, p, q) in enumerate(top_A, 1):
    print(f"    [{i}] {p}/{q} = {sp.Float(r, 12)}  drift = {sp.Float(d, 8)}%")
best_drift_A = top_A[0][1]
T1_1_pass = float(best_drift_A) < 1.0
results["T1.1"] = T1_1_pass
print(f"  Best candidate drift = {sp.Float(best_drift_A, 8)}%")
print(f"  Verdict: {'PASS' if T1_1_pass else 'FAIL'}  (drift < 1% required)")
print()

# ------------------------------------------------------------------
# T1.2 — ρ̄ = 0.141 PDG + top-5 rational candidates
# ------------------------------------------------------------------
print("-" * 78)
print("T1.2 — ρ̄ = 0.141 PDG + top-5 rational candidates (denom ≤ 100)")
print("-" * 78)
top_rho = best_rationals(RHO_PDG, 100, top_n=5)
print(f"  Target: ρ̄_PDG = {RHO_PDG}")
print(f"  Top-5 rational candidates (denom ≤ 100):")
for i, (r, d, p, q) in enumerate(top_rho, 1):
    print(f"    [{i}] {p}/{q} = {sp.Float(r, 12)}  drift = {sp.Float(d, 8)}%")
best_drift_rho = top_rho[0][1]
T1_2_pass = float(best_drift_rho) < 2.0
results["T1.2"] = T1_2_pass
print(f"  Best candidate drift = {sp.Float(best_drift_rho, 8)}%")
print(f"  Verdict: {'PASS' if T1_2_pass else 'FAIL'}  (drift < 2% required)")
print()

# ------------------------------------------------------------------
# T1.3 — η̄ = 0.357 PDG + top-5 rational candidates (5/14 hypothesis)
# ------------------------------------------------------------------
print("-" * 78)
print("T1.3 — η̄ = 0.357 PDG + top-5 rational candidates (5/14 hypothesis)")
print("-" * 78)
top_eta = best_rationals(ETA_PDG, 100, top_n=5)
eta_5_14 = sp.Rational(5, 14)
eta_5_14_drift = drift_pct(sp.Float(eta_5_14, 30), ETA_PDG)
print(f"  Target: η̄_PDG = {ETA_PDG}")
print(f"  Top-5 rational candidates (denom ≤ 100):")
for i, (r, d, p, q) in enumerate(top_eta, 1):
    is_5_14 = (r.p == 5 and r.q == 14)
    marker = " ← 5/14 hypothesis" if is_5_14 else ""
    print(f"    [{i}] {p}/{q} = {sp.Float(r, 12)}  drift = {sp.Float(d, 8)}%{marker}")
print(f"  Explicit 5/14 check: drift = {sp.Float(eta_5_14_drift, 8)}%")
top_is_5_14 = (top_eta[0][0] == eta_5_14)
T1_3_pass = top_is_5_14 and float(eta_5_14_drift) < 0.1
results["T1.3"] = T1_3_pass
print(f"  5/14 ranks #1: {top_is_5_14}")
print(f"  Verdict: {'PASS' if T1_3_pass else 'FAIL'}  (5/14 #1, drift < 0.1% required)")
print()

# ------------------------------------------------------------------
# T1.4 — Unitarity triangle apex closure α + β + γ = π
# ------------------------------------------------------------------
print("-" * 78)
print("T1.4 — Unitarity triangle apex closure α + β + γ = π")
print("-" * 78)

# Compute V_us, V_cb, V_ub, V_td z Wolfenstein O(λ⁴)
# Using PDG λ (not λ_C) for self-consistency check
lam = LAMBDA_PDG
A_v = A_PDG
rho = RHO_PDG
eta = ETA_PDG

# V_us = λ + O(λ^7)
V_us = lam
# V_cb = A·λ²
V_cb = A_v * lam ** 2
# V_ub = A·λ³·(ρ - iη)  (note Wolfenstein uses ρ + iη with sign convention)
V_ub_complex = A_v * lam ** 3 * (rho - sp.I * eta)
# V_ud ≈ 1 - λ²/2
V_ud = 1 - lam ** 2 / 2
# V_cd ≈ -λ
V_cd = -lam
# V_tb ≈ 1
V_tb = sp.Integer(1)
# V_td = A·λ³·(1 - ρ - iη)·(...)  → leading: A·λ³·(1 - ρ + iη)·(-1) =
#  use V_td = A·λ³·(1 - (ρ - iη))
V_td_complex = A_v * lam ** 3 * (1 - (rho - sp.I * eta))

# Conjugates
V_ub_conj = sp.conjugate(V_ub_complex)
V_td_conj = sp.conjugate(V_td_complex)
V_cb_conj = V_cb  # real to this order
V_tb_conj = V_tb
V_ud_conj = V_ud
V_cd_conj = V_cd

# Triangle angles (standard CKM unitarity):
# α = arg(- V_td · V_tb* / (V_ud · V_ub*))
# β = arg(- V_cd · V_cb* / (V_td · V_tb*))
# γ = arg(- V_ud · V_ub* / (V_cd · V_cb*))
num_alpha = -V_td_complex * V_tb_conj
den_alpha = V_ud * V_ub_conj
ratio_alpha = sp.simplify(num_alpha / den_alpha)
alpha_v = sp.arg(ratio_alpha)

num_beta = -V_cd * V_cb_conj
den_beta = V_td_complex * V_tb_conj
ratio_beta = sp.simplify(num_beta / den_beta)
beta_v = sp.arg(ratio_beta)

num_gamma = -V_ud * V_ub_conj
den_gamma = V_cd * V_cb_conj
ratio_gamma = sp.simplify(num_gamma / den_gamma)
gamma_v = sp.arg(ratio_gamma)

alpha_n = sp.N(alpha_v, 30)
beta_n = sp.N(beta_v, 30)
gamma_n = sp.N(gamma_v, 30)
sum_n = sp.N(alpha_n + beta_n + gamma_n, 30)
closure_diff = sp.Abs(sum_n - PI)

def to_deg(rad):
    return sp.Float(sp.N(rad * 180 / PI, 30), 12)

print(f"  α = {sp.Float(alpha_n, 12)} rad  ({to_deg(alpha_n)}°)")
print(f"  β = {sp.Float(beta_n, 12)} rad  ({to_deg(beta_n)}°)")
print(f"  γ = {sp.Float(gamma_n, 12)} rad  ({to_deg(gamma_n)}°)")
print(f"  Σ (α+β+γ) = {sp.Float(sum_n, 12)} rad  ({to_deg(sum_n)}°)")
print(f"  π          = {sp.Float(PI, 12)} rad")
print(f"  |Σ - π|    = {sp.Float(closure_diff, 12)} rad")
T1_4_pass = float(closure_diff) < 1e-3
results["T1.4"] = T1_4_pass
print(f"  Verdict: {'PASS' if T1_4_pass else 'FAIL'}  (closure < 10⁻³ rad required)")
print()

# ------------------------------------------------------------------
# T1.5 — Jarlskog J = A²·λ_C⁶·η̄ z PDG values vs PDG J
# ------------------------------------------------------------------
print("-" * 78)
print("T1.5 — Jarlskog J = A²·λ_C⁶·η̄ vs PDG J")
print("-" * 78)
J_TGP = A_PDG ** 2 * LAMBDA_C ** 6 * ETA_PDG
J_drift = drift_pct(J_TGP, J_PDG)
print(f"  J_TGP (z PDG A, λ_C, η̄) = {sp.Float(J_TGP, 12)}")
print(f"  J_PDG                    = {sp.Float(J_PDG, 12)}")
print(f"  Drift                    = {sp.Float(J_drift, 8)}%")
T1_5_pass = float(J_drift) < 5.0
results["T1.5"] = T1_5_pass
print(f"  Verdict: {'PASS' if T1_5_pass else 'FAIL'}  (drift < 5% required)")
print()

# ------------------------------------------------------------------
# Final verdict
# ------------------------------------------------------------------
print("=" * 78)
print("η.1.Phase1 verdict summary")
print("=" * 78)
for k, v in results.items():
    print(f"  {k}: {'PASS' if v else 'FAIL'}")
n_pass = sum(1 for v in results.values() if v)
n_total = len(results)
print()
print(f"  Total: {n_pass}/{n_total} PASS")
if n_pass == n_total:
    print(f"  → η.1.Phase2 proceeds z best-candidate triple.")
elif n_pass >= 4:
    print(f"  → η.1.Phase2 proceeds z minor gap.")
else:
    print(f"  → η.1.Phase1 reframing required.")
print("=" * 78)
