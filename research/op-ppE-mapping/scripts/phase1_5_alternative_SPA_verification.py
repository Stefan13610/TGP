"""
op-ppE-mapping Phase 1.5 — Alternative SPA verification (orthogonal route)

Date: 2026-05-09
Purpose: Independent cross-check of β_ppE^TGP^(b=-1) = -15/4 finding.
Original Phase 1.5 derivation (phase1_5_G_SPA_derivation.py) used:
  1. Symbolic E(v) = -(η v²/2)(1 + e_1 v² + e_2 v⁴ + ...) with e_n placeholders.
  2. Symbolic F(v) similarly.
  3. Integrate symbolically to get Ψ(v) coefficients.
  4. Extract α_4 = 30·e_2 - 20·e_1·p_1 + 10·p_1² - 10·p_2 (TaylorF2 2PN-phase formula).
  5. Substitute TGP and GR specific e_n, p_n values; compute δα_4 = -40.
  6. β_ppE = (3/(128η)) · δα_4 = -15/4 at η=1/4.

This ALTERNATIVE script:
  1. Substitute concrete e_n^TGP, e_n^GR, p_n values into E(v), F(v) FIRST.
  2. Compute dt/dv = -E'(v)/F(v) for TGP and GR separately, with concrete numbers.
  3. Integrate to get t_TGP(v), t_GR(v); similarly Φ_TGP(v), Φ_GR(v).
  4. Compute ΔΨ(v) = Ψ_TGP(v) - Ψ_GR(v) directly.
  5. Read v^(-1) coefficient → β_ppE^(b=-1).

If this gives the same β = -15/4, then the original derivation is verified
methodologically (no α_4 extraction error). If different, identifies a chain
issue.

Cross-check goal: Δα_4 = -40 (sympy) → β_ppE = -15/4 (η=1/4) — verify or revise.
"""

import sympy as sp

# ============================================================================
# §1 — Concrete e_n^TGP, e_n^GR, p_n values (locked from Phase 1.5 / literature)
# ============================================================================
v = sp.Symbol('v', positive=True)
eta = sp.Rational(1, 4)  # equal-mass binary

# Test-particle e_n (extracted via E_b/m = -x/2(1 + e_1 x + e_2 x² + ...))
e_GR_test = {
    1: sp.Rational(-3, 4),
    2: sp.Rational(-27, 8),
    3: sp.Rational(-675, 64),
}
e_TGP_test = {
    1: sp.Rational(-3, 4),  # 1PN matches GR exactly (β=γ=1)
    2: sp.Rational(-113, 24),  # Phase 1.5 LOCK L3
    3: sp.Rational(-10315, 576),
}

# GR test-particle flux coefficients (Blanchet 2014 LR §11)
p_GR = {
    1: sp.Rational(-1247, 336),
    2: sp.Rational(-44711, 9072),
    3: sp.Rational(0, 1),
}
# LOCK L4: F_TGP = F_GR for M9.1'' "metric-only" leading 2PN-orbital
p_TGP = p_GR

print("=" * 80)
print("op-ppE-mapping Phase 1.5 — ALTERNATIVE SPA verification")
print("Date: 2026-05-09")
print("=" * 80)
print()
print("§1 — Concrete e_n, p_n values (test-particle limit, LOCK L4 F_TGP=F_GR):")
print(f"  GR test-particle:  e_1 = {e_GR_test[1]}, e_2 = {e_GR_test[2]}, e_3 = {e_GR_test[3]}")
print(f"  TGP test-particle: e_1 = {e_TGP_test[1]}, e_2 = {e_TGP_test[2]}, e_3 = {e_TGP_test[3]}")
print(f"  Δe_1 = {sp.simplify(e_TGP_test[1] - e_GR_test[1])}")
print(f"  Δe_2 = {sp.simplify(e_TGP_test[2] - e_GR_test[2])}")
print(f"  Δe_3 = {sp.simplify(e_TGP_test[3] - e_GR_test[3])}")
print(f"  p_1 = {p_GR[1]}, p_2 = {p_GR[2]}, p_3 = {p_GR[3]}")
print()

# ============================================================================
# §2 — Build E(v), F(v) with concrete numerical coefficients
# ============================================================================
def E_concrete(e_dict):
    """E(v) = -(η v²/2)·(1 + e_1 v² + e_2 v⁴ + e_3 v⁶) with concrete e_n"""
    return -(eta * v**2 / 2) * (1 + e_dict[1]*v**2 + e_dict[2]*v**4 + e_dict[3]*v**6)

def F_concrete(p_dict):
    """F(v) = (32 η²/5)·v¹⁰·(1 + p_1 v² + p_2 v⁴ + p_3 v⁶) with concrete p_n"""
    return sp.Rational(32, 5) * eta**2 * v**10 * (1 + p_dict[1]*v**2 + p_dict[2]*v**4 + p_dict[3]*v**6)

E_GR_full = E_concrete(e_GR_test)
E_TGP_full = E_concrete(e_TGP_test)
F_GR_full = F_concrete(p_GR)
F_TGP_full = F_concrete(p_TGP)

print("§2 — E(v), F(v) concrete (numerical) at η=1/4:")
print(f"  E_GR(v) = {sp.expand(E_GR_full)}")
print(f"  E_TGP(v) = {sp.expand(E_TGP_full)}")
print(f"  F(v) = {sp.expand(F_GR_full)} (same for TGP and GR by LOCK L4)")
print()

# ============================================================================
# §3 — Compute dt/dv = -E'(v)/F(v) for TGP and GR separately
# ============================================================================
dt_dv_GR = -sp.diff(E_GR_full, v) / F_GR_full
dt_dv_TGP = -sp.diff(E_TGP_full, v) / F_TGP_full

# Series expand each (dt/dv has leading 1/v^9 in the integrand after cancellations)
# In sympy, we can expand using fractions: dt/dv = (a polynomial in v starting v^1) / (v^10 · polynomial in v starting at 1)
# Multiply numerator times (1/(F polynomial)) expansion.

# Direct approach: simplify and expand
# dt/dv * v^9 should be a regular polynomial in v at v=0
# Let's compute v^9 * dt/dv and expand around v=0:

dt_dv_GR_v9 = sp.series(v**9 * dt_dv_GR, v, 0, 7).removeO()
dt_dv_GR_v9 = sp.expand(dt_dv_GR_v9)
dt_dv_TGP_v9 = sp.series(v**9 * dt_dv_TGP, v, 0, 7).removeO()
dt_dv_TGP_v9 = sp.expand(dt_dv_TGP_v9)

print("§3 — dt/dv · v^9 (regular polynomial near v=0):")
print(f"  GR:  {dt_dv_GR_v9}")
print(f"  TGP: {dt_dv_TGP_v9}")
print(f"  Δ:   {sp.expand(dt_dv_TGP_v9 - dt_dv_GR_v9)}")
print()

# ============================================================================
# §4 — Integrate dt/dv → t(v); compute (v³/M)·dt/dv → dΦ/dv → Φ(v)
# ============================================================================
# We have v^9·(dt/dv) = polynomial p(v) = a_0 + a_2 v^2 + a_4 v^4 + a_6 v^6 + ...
# So dt/dv = (a_0 + a_2 v² + a_4 v⁴ + a_6 v⁶)/v^9
#         = a_0 v^(-9) + a_2 v^(-7) + a_4 v^(-5) + a_6 v^(-3)
# Integral: t(v) = a_0 v^(-8)/(-8) + a_2 v^(-6)/(-6) + a_4 v^(-4)/(-4) + a_6 v^(-2)/(-2)
#                 = -a_0/(8 v^8) - a_2/(6 v^6) - a_4/(4 v^4) - a_6/(2 v^2)

def extract_coeffs(poly_in_v, max_power):
    """Extract coefficients of v^k for k = 0, 2, 4, 6 from polynomial in v."""
    coeffs = {}
    for k in range(0, max_power + 1, 2):
        coeffs[k] = sp.simplify(poly_in_v.coeff(v, k))
    return coeffs

a_GR = extract_coeffs(dt_dv_GR_v9, 6)
a_TGP = extract_coeffs(dt_dv_TGP_v9, 6)

print("§4 — Coefficients a_k of dt/dv · v^9 = a_0 + a_2 v² + a_4 v⁴ + a_6 v⁶:")
print(f"{'k':<4}{'a_k^GR':<25}{'a_k^TGP':<25}{'Δa_k':<25}")
for k in [0, 2, 4, 6]:
    da = sp.simplify(a_TGP[k] - a_GR[k])
    print(f"{k:<4}{str(a_GR[k]):<25}{str(a_TGP[k]):<25}{str(da):<25}")
print()

# t(v) = -[a_0/(8 v^8) + a_2/(6 v^6) + a_4/(4 v^4) + a_6/(2 v^2)]
def t_of_v(a_dict):
    return -(a_dict[0]/(8*v**8) + a_dict[2]/(6*v**6) + a_dict[4]/(4*v**4) + a_dict[6]/(2*v**2))

t_GR = t_of_v(a_GR)
t_TGP = t_of_v(a_TGP)

# dΦ/dv = (v³/M)·dt/dv. In M=1 units: dΦ/dv = v³·dt/dv = v³·v^(-9)·p(v) = p(v)/v^6
# Integrate: Φ(v) = -[a_0/(5 v^5) + a_2/(3 v^3) + a_4/v - a_6 v]
def phi_of_v(a_dict):
    return -(a_dict[0]/(5*v**5) + a_dict[2]/(3*v**3) + a_dict[4]/v) + a_dict[6]*v

phi_GR = phi_of_v(a_GR)
phi_TGP = phi_of_v(a_TGP)

# ============================================================================
# §5 — Ψ(v) = 2 v³ t(v) - 2 Φ(v) - π/4 (in M=1)
# ============================================================================
Psi_GR = 2 * v**3 * t_GR - 2 * phi_GR
Psi_TGP = 2 * v**3 * t_TGP - 2 * phi_TGP

Psi_GR_expanded = sp.expand(Psi_GR)
Psi_TGP_expanded = sp.expand(Psi_TGP)

# ΔΨ(v) = Ψ_TGP - Ψ_GR
delta_Psi = sp.expand(Psi_TGP_expanded - Psi_GR_expanded)
delta_Psi_simplified = sp.simplify(delta_Psi)

print("§5 — Ψ(v) = 2 v³ t(v) - 2 Φ(v) (in M=1, η=1/4):")
print(f"  Ψ_GR(v)  = {Psi_GR_expanded}")
print(f"  Ψ_TGP(v) = {Psi_TGP_expanded}")
print()
print(f"ΔΨ(v) = Ψ_TGP - Ψ_GR = {delta_Psi_simplified}")
print()

# ============================================================================
# §6 — Read v^(-1) coefficient → β_ppE^(b=-1)
# ============================================================================
# In ppE convention: δΨ(v) = β_ppE · v^b at b = -1: δΨ = β/v.
# Identification: β_ppE^(b=-1) = coefficient of v^(-1) in ΔΨ(v).
# Extract Laurent coefficient at v^(-1) via: ΔΨ · v at v=0 → constant = β_ppE^(b=-1).

def laurent_coeff(expr, var, k):
    """Extract coefficient of var^k for any integer k (positive or negative).
    For k < 0: multiply by var^(-k) and take limit at var=0.
    For k >= 0: standard series coefficient."""
    if k < 0:
        shifted = sp.simplify(expr * var**(-k))
        return sp.simplify(shifted.subs(var, 0))
    else:
        return sp.simplify(sp.series(expr, var, 0, k+1).removeO().coeff(var, k))

beta_ppE_alternative = laurent_coeff(delta_Psi_simplified, v, -1)

print("=" * 80)
print("§6 — β_ppE^TGP^(b=-1) extracted directly from ΔΨ(v) at v^(-1)")
print("=" * 80)
print()
print(f"β_ppE^TGP^(b=-1) (alternative SPA route) = {beta_ppE_alternative}")
print(f"                                         ≈ {float(beta_ppE_alternative):.6f}")
print()

# Compare to original Phase 1.5 derivation:
beta_ppE_original = sp.Rational(-15, 4)
print(f"β_ppE^TGP^(b=-1) (original Phase 1.5)    = {beta_ppE_original}")
print(f"                                         ≈ {float(beta_ppE_original):.6f}")
print()

# Verification:
diff = sp.simplify(beta_ppE_alternative - beta_ppE_original)
print(f"Difference (alternative - original) = {diff}")
if sp.simplify(diff) == 0:
    print(f"  → ALTERNATIVE SPA VERIFIES ORIGINAL Phase 1.5 result EXACTLY ✓")
    print(f"  → Possibility (B) [methodological error in Phase 1.5] RULED OUT")
else:
    print(f"  → ALTERNATIVE DISAGREES with original by {diff}")
    print(f"  → Phase 1.5 has a methodological error; needs revision.")
print()

# ============================================================================
# §7 — Also extract v^(-3) coefficient (1PN-phase, b=-3) as sanity check
# ============================================================================
# At b=-3 (1PN-phase), TGP β should be exactly 0 since 1PN matches GR (β=γ=1).
beta_ppE_b3 = laurent_coeff(delta_Psi_simplified, v, -3)
print(f"§7 — Sanity check: β_ppE^TGP^(b=-3) (1PN-phase)")
print(f"  Should be 0 (TGP matches GR exactly at 1PN: β_PPN=γ_PPN=1).")
print(f"  Alternative result: β_ppE^TGP^(b=-3) = {beta_ppE_b3}")
if sp.simplify(beta_ppE_b3) == 0:
    print(f"  → MATCHES expected 0 ✓")
else:
    print(f"  → UNEXPECTED non-zero (1PN should match GR)")
print()

# ============================================================================
# §8 — Also extract v^(-5) coefficient (0PN, leading)
# ============================================================================
# At b=-5 (0PN-phase), this is the LEADING phase. Both TGP and GR should
# give the same coefficient (3/(128η)). So Δ should be 0.
beta_ppE_b5 = laurent_coeff(delta_Psi_simplified, v, -5)
print(f"§8 — Sanity check: ΔΨ(v) at v^(-5) (0PN leading)")
print(f"  Should be 0 (TGP and GR have same 0PN-phase prefactor 3/(128η)).")
print(f"  Alternative result: Δ at v^(-5) = {beta_ppE_b5}")
if sp.simplify(beta_ppE_b5) == 0:
    print(f"  → MATCHES expected 0 ✓")
else:
    print(f"  → UNEXPECTED non-zero")
print()

# Also v^(+1) coefficient (3PN-phase, b=+1):
beta_ppE_p1 = laurent_coeff(delta_Psi_simplified, v, 1)
print(f"§9 — Bonus: β_ppE^TGP^(b=+1) (3PN-phase) — alternative-derived:")
print(f"  β at v^(+1) = {beta_ppE_p1}")
print(f"            ≈ {float(beta_ppE_p1):.6f}")
print(f"  (Compare to TGP's M911-P2 multi-coefficient predicted ratio β_3PN/β_2PN = -23/10.")
print(f"   Phase 1 gave β_3PN at central -23/128. Alternative...)")
beta_ppE_b1_phase1 = sp.Rational(-23, 128)  # Phase 1 OOM, η=1/4, G=1
ratio_alt_to_phase1 = sp.simplify(beta_ppE_p1 / beta_ppE_b1_phase1)
print(f"  Alternative / Phase 1 ratio = {ratio_alt_to_phase1}")

# Also β_3PN/β_2PN ratio:
if beta_ppE_alternative != 0:
    ratio_3PN_to_2PN = sp.simplify(beta_ppE_p1 / beta_ppE_alternative)
    print(f"  Alternative β_3PN/β_2PN ratio = {ratio_3PN_to_2PN}")
    print(f"  Phase 1 predicted ratio = -23/10")
print()

print("=" * 80)
print("ALTERNATIVE SPA VERIFICATION COMPLETE 2026-05-09")
print("=" * 80)
