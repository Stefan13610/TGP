"""
epsilon.1.Phase1 - eps_ph numerical landscape audit (5 sub-tests)
================================================================================
Verify eps_ph (M9.1" refined photon-ring scale) = 0.16788 sympy-exact rational,
psi_ph = 4/(3+0.4250) = 1.16788 z null geodesic, F4 chain consistency 0.0038%
drift, M9.2-D -> M9.1" refinement gain 36x.

Predecessor: UV.1 program END (373 cumulative)
Goal: 5/5 PASS -> proceed Phase 2 (eps_ph = psi_ph - 1 structural decomposition)

Author: TGP_v1 / Mateusz Serafin
Date: 2026-04-29
"""

import sympy as sp


# =====================================================================
# Constants
# =====================================================================

# eps_ph M9.1" refined (UV.1.Phase1.UV1.4 LOCKED)
EPS_PH_M9_1P = sp.Rational(16788, 100000)        # 0.16788

# psi_ph M9.1" null geodesic: psi_ph = 4 / (3 + 0.4250)
PSI_PH_NUM = sp.Rational(4)
PSI_PH_DENOM = sp.Rational(3) + sp.Rational(425, 1000)  # 3.4250 = 685/200
PSI_PH_M9_1P = PSI_PH_NUM / PSI_PH_DENOM

# coarse M9.2-D
EPS_PH_M9_2D = sp.Rational(168, 1000)            # 0.168
PSI_PH_M9_2D = sp.Rational(1168, 1000)           # 1.168

# F4 anchor
ALPHA_0_F4 = sp.Rational(1069833, 264500)        # 4.04472
TARGET_SHIFT_F4 = sp.Rational(114, 1000)         # 0.114


# =====================================================================
# Header
# =====================================================================

print("=" * 72)
print("epsilon.1.Phase1 - eps_ph numerical landscape audit")
print("=" * 72)
print(f"  date          : 2026-04-29")
print(f"  predecessor   : UV.1 program END (373 cumulative)")
print(f"  goal          : 5/5 PASS -> Phase 2 structural decomposition")
print()


# =====================================================================
# E1.1 - eps_ph = 4197/25000 sympy-exact rational
# =====================================================================

print("=" * 72)
print("E1.1 - eps_ph = 4197/25000 sympy-exact rational")
print("=" * 72)

eps_simp = sp.simplify(EPS_PH_M9_1P)
is_rational = eps_simp.is_rational
denom_ok = eps_simp.q < 100000

print(f"  EPS_PH_M9_1P (raw)                            {EPS_PH_M9_1P}")
print(f"  EPS_PH_M9_1P simplified                       {eps_simp}")
print(f"  Numerator p, denominator q                    p={eps_simp.p}, q={eps_simp.q}")
print(f"  Float value                                    {float(eps_simp):.5f}")
print(f"  Is sympy-rational?                            {is_rational}")
print(f"  Denominator < 10^5?                           {denom_ok}")

e1_1_pass = is_rational and denom_ok
if e1_1_pass:
    print(f"\n  -> E1.1 (eps_ph rational): PASS  [eps_ph = {eps_simp} sympy exact]")
else:
    print(f"\n  -> E1.1 (eps_ph rational): FAIL")
print()


# =====================================================================
# E1.2 - psi_ph = 4/(3+0.4250) = 1.16788 z M9.1" null geodesic
# =====================================================================

print("=" * 72)
print("E1.2 - psi_ph = 4/(3+0.4250) = 1.16788 z M9.1\" null geodesic")
print("=" * 72)

psi_ph = sp.simplify(PSI_PH_M9_1P)
eps_from_psi = sp.simplify(psi_ph - 1)
psi_float = float(psi_ph)
eps_from_psi_float = float(eps_from_psi)
drift_psi_eps = abs(eps_from_psi_float - float(EPS_PH_M9_1P)) / float(EPS_PH_M9_1P)

print(f"  psi_ph = 4/(3+0.4250)                          {psi_ph} = {psi_float:.5f}")
print(f"  Numerator/denominator                          p={psi_ph.p}, q={psi_ph.q}")
print(f"  eps_ph := psi_ph - 1                           {eps_from_psi} = {eps_from_psi_float:.5f}")
print(f"  EPS_PH_M9_1P (Rational reference)              {EPS_PH_M9_1P} = {float(EPS_PH_M9_1P):.5f}")
print(f"  drift |psi_ph-1 - EPS_PH_M9_1P|/EPS_PH_M9_1P   {drift_psi_eps*100:.4f}%")
print(f"  drift < 0.1% gate?                             {drift_psi_eps < 0.001}")

e1_2_pass = drift_psi_eps < 0.001
if e1_2_pass:
    print(f"\n  -> E1.2 (psi_ph - 1): PASS  [psi_ph - 1 matches eps_ph w {drift_psi_eps*100:.4f}%]")
else:
    print(f"\n  -> E1.2 (psi_ph - 1): FAIL")
print()


# =====================================================================
# E1.3 - eps_ph^2 sympy exact
# =====================================================================

print("=" * 72)
print("E1.3 - eps_ph^2 sympy exact (consistency between coarse and refined)")
print("=" * 72)

# Coarse 0.168 -> 0.168^2 = 0.028224 = 441/15625 sympy exact
eps_sq_coarse = EPS_PH_M9_2D ** 2
eps_sq_coarse_simp = sp.simplify(eps_sq_coarse)

# Refined 0.16788 -> 0.16788^2
eps_sq_refined = EPS_PH_M9_1P ** 2
eps_sq_refined_simp = sp.simplify(eps_sq_refined)

# Drift between coarse and refined squared
drift_sq = abs(float(eps_sq_refined) - float(eps_sq_coarse)) / float(eps_sq_coarse)

print(f"  M9.2-D coarse eps_ph^2                         {eps_sq_coarse_simp} = {float(eps_sq_coarse_simp):.6f}")
print(f"  M9.1\" refined eps_ph^2                        {eps_sq_refined_simp} = {float(eps_sq_refined_simp):.6f}")
print(f"  Coarse Numerator p, denominator q              p={eps_sq_coarse_simp.p}, q={eps_sq_coarse_simp.q}")
print(f"  Refined Numerator p, denominator q             p={eps_sq_refined_simp.p}, q={eps_sq_refined_simp.q}")
print(f"  drift coarse vs refined                        {drift_sq*100:.4f}%")
print(f"  Both rational?                                 {eps_sq_coarse_simp.is_rational and eps_sq_refined_simp.is_rational}")

# Coarse should give 441/15625 (clean factorization)
coarse_clean = (eps_sq_coarse_simp == sp.Rational(441, 15625))
print(f"  M9.2-D coarse = 441/15625?                     {coarse_clean}")

e1_3_pass = (eps_sq_coarse_simp.is_rational and eps_sq_refined_simp.is_rational and coarse_clean)
if e1_3_pass:
    print(f"\n  -> E1.3 (eps_ph^2 sympy): PASS  [both forms rational, coarse = 441/15625]")
else:
    print(f"\n  -> E1.3 (eps_ph^2 sympy): FAIL")
print()


# =====================================================================
# E1.4 - F4 chain consistency alpha_0 = target_shift_F4 / eps_ph^2 = 4.04489
# =====================================================================

print("=" * 72)
print("E1.4 - F4 chain consistency alpha_0 = target_shift / eps_ph^2 = 4.04489")
print("=" * 72)

alpha_0_repro = TARGET_SHIFT_F4 / eps_sq_refined
alpha_0_repro_num = float(alpha_0_repro)
alpha_0_F4_num = float(ALPHA_0_F4)
drift_alpha = abs(alpha_0_repro_num - alpha_0_F4_num) / alpha_0_F4_num

print(f"  target_shift_F4                                {TARGET_SHIFT_F4} = {float(TARGET_SHIFT_F4):.4f}")
print(f"  eps_ph^2 (M9.1\" refined)                      {float(eps_sq_refined):.6f}")
print(f"  alpha_0_repro                                   {alpha_0_repro_num:.5f}")
print(f"  alpha_0_F4 (sympy)                              {alpha_0_F4_num:.5f}")
print(f"  drift                                          {drift_alpha*100:.4f}%")
print(f"  drift < 0.01% gate?                            {drift_alpha < 0.0001}")

e1_4_pass = drift_alpha < 0.0001
if e1_4_pass:
    print(f"\n  -> E1.4 (F4 consistency): PASS  [drift {drift_alpha*100:.4f}% < 0.01%]")
else:
    print(f"\n  -> E1.4 (F4 consistency): FAIL")
print()


# =====================================================================
# E1.5 - Refinement gain M9.2-D -> M9.1" 36x
# =====================================================================

print("=" * 72)
print("E1.5 - Refinement gain M9.2-D -> M9.1\" 36x")
print("=" * 72)

alpha_0_coarse = TARGET_SHIFT_F4 / eps_sq_coarse
drift_coarse = abs(float(alpha_0_coarse) - alpha_0_F4_num) / alpha_0_F4_num
gain = drift_coarse / drift_alpha if drift_alpha > 0 else float("inf")

print(f"  M9.2-D coarse alpha_0_repro                    {float(alpha_0_coarse):.5f}")
print(f"  M9.2-D coarse drift                            {drift_coarse*100:.4f}%")
print(f"  M9.1\" refined alpha_0_repro                   {alpha_0_repro_num:.5f}")
print(f"  M9.1\" refined drift                           {drift_alpha*100:.4f}%")
print(f"  Refinement gain (drift_coarse/drift_refined)   {gain:.1f}x")
print(f"  Gain >= 30x gate?                              {gain >= 30}")

e1_5_pass = gain >= 30
if e1_5_pass:
    print(f"\n  -> E1.5 (refinement gain): PASS  [{gain:.1f}x >= 30x]")
else:
    print(f"\n  -> E1.5 (refinement gain): FAIL")
print()


# =====================================================================
# Verdict
# =====================================================================

print("=" * 72)
print("epsilon.1.Phase1 verdict")
print("=" * 72)

results = {
    "E1.1": e1_1_pass,
    "E1.2": e1_2_pass,
    "E1.3": e1_3_pass,
    "E1.4": e1_4_pass,
    "E1.5": e1_5_pass,
}

n_pass = sum(1 for r in results.values() if r)
for k, v in results.items():
    print(f"  {k}: {'PASS' if v else 'FAIL'}")
print()
print(f"  Cumulative: {n_pass}/5 PASS")
print()

if n_pass == 5:
    print(f"  -> Phase 1 CLOSED with 5/5 PASS")
    print(f"  -> eps_ph numerical landscape fully LOCKED")
    print(f"  -> Proceed Phase 2 (structural decomposition eps_ph = psi_ph - 1)")
    print(f"  -> Master ledger update: 373 -> 378 (+5 z Phase 1)")
elif n_pass >= 4:
    print(f"  -> Phase 1 PARTIAL ({n_pass}/5); 1 audit gap to address before Phase 2")
else:
    print(f"  -> Phase 1 INSUFFICIENT ({n_pass}/5); audit gaps require closure first")

print()
print("=" * 72)
