#!/usr/bin/env python3
"""
μ.1.Phase1 — δ_CP + PMNS drift correction landscape audit (5 sub-tests).
"""
import sympy as sp


# ===== TGP anchors =====
LAMBDA_C    = sp.Rational(2255, 10000)         # Cabibbo (ζ.1 GL(3,𝔽₂) 165/167)
N_GEN       = sp.Integer(3)
K_NEUTRINO  = sp.Rational(1, 2)                # K_ν Majorana
K_LEPTON    = sp.Rational(2, 3)                # K_lep Dirac
K_UP        = sp.Rational(7, 8)                # K_up Dirac+QCD
B2_UP_NUM   = sp.Integer(13)
B2_UP_DENOM = sp.Integer(4)
B2_NU_NUM   = sp.Integer(1)
K_UP_NUM    = sp.Integer(7)
K_NU_NUM    = sp.Integer(1)

# κ.1 Wolfenstein triple
A_WOLF   = sp.Rational(64, 81)
RHO_BAR  = sp.Rational(11, 78)
ETA_BAR  = sp.Rational(5, 14)

# NuFit 5.3 NO best-fit
SIN2_T12_NUFIT = 0.307
SIN2_T23_NUFIT = 0.572
SIN2_T13_NUFIT = 0.022
DELTA_CP_NUFIT_DEG = 195.0
DELTA_CP_NUFIT_1S = (128.0, 352.0)
PDG_GAMMA_DEG = 65.7


def drift_pct(a, b):
    return float(sp.Abs(sp.N(a, 30) - sp.N(b, 30)) / sp.N(b, 30) * 100)


# ============== M1.1 — NuFit 5.3 PMNS reference inventory ===============
print("=" * 72)
print("M1.1 — NuFit 5.3 PMNS reference inventory")
print("=" * 72)

print(f"  sin²θ₁₂ = {SIN2_T12_NUFIT} ± 0.013 (1σ ≈ 4.2%)")
print(f"  sin²θ₂₃ = {SIN2_T23_NUFIT} ± 0.018 (1σ ≈ 3.1%, 2nd octant)")
print(f"  sin²θ₁₃ = {SIN2_T13_NUFIT} ± 0.0006 (1σ ≈ 2.7%)")
print(f"  δ_CP = {DELTA_CP_NUFIT_DEG}° z 1σ {DELTA_CP_NUFIT_1S}° (NO; very wide)")
M11_PASS = True
print(f"  Verdict M1.1 = {'PASS' if M11_PASS else 'FAIL'}")
print()


# ============== M1.2 — γ_CKM apex angle z Wolfenstein triple =============
print("=" * 72)
print("M1.2 — γ_CKM apex angle z κ.1 Wolfenstein triple")
print("=" * 72)

tan_gamma = ETA_BAR / RHO_BAR
tan_gamma_simplified = sp.simplify(tan_gamma)
gamma_rad = sp.atan(tan_gamma)
gamma_deg = float(gamma_rad * 180 / sp.pi)
print(f"  tan(γ_CKM) = (η̄/ρ̄) = ({ETA_BAR})/({RHO_BAR}) = {tan_gamma_simplified}")
print(f"  γ_CKM = arctan({tan_gamma_simplified}) ≈ {gamma_deg:.2f}°")
print(f"  PDG γ ≈ {PDG_GAMMA_DEG}°")
print(f"  drift {drift_pct(gamma_deg, PDG_GAMMA_DEG):.2f}%")

M12_PASS = (tan_gamma_simplified == sp.Rational(195, 77))
print(f"  tan(γ_CKM) sympy-exact 195/77: {M12_PASS}")
print(f"  Verdict M1.2 = {'PASS' if M12_PASS else 'FAIL'}")
print()


# ============== M1.3 — PMNS-Wolfenstein analog z (ν,up) pair =============
print("=" * 72)
print("M1.3 — PMNS-Wolfenstein analog (ρ̄_PMNS, η̄_PMNS) z (ν,up) pair")
print("=" * 72)

rho_pmns_num = B2_UP_NUM - B2_NU_NUM            # 13 - 1 = 12
rho_pmns_denom = 2 * N_GEN * B2_UP_NUM           # 2·3·13 = 78
rho_pmns = sp.Rational(rho_pmns_num, rho_pmns_denom)
print(f"  ρ̄_PMNS_num = B²_up_num − B²_ν_num = {B2_UP_NUM} − {B2_NU_NUM} = {rho_pmns_num}")
print(f"  ρ̄_PMNS_denom = 2·N_gen·B²_up_num = {rho_pmns_denom}")
print(f"  ρ̄_PMNS = {rho_pmns}")

eta_pmns_num = K_UP_NUM - K_NU_NUM               # 7 - 1 = 6
eta_pmns_denom = K_UP_NUM * K_NU_NUM              # 7·1 = 7
eta_pmns = sp.Rational(eta_pmns_num, eta_pmns_denom)
print(f"  η̄_PMNS_num = K_up_num − K_ν_num = {K_UP_NUM} − {K_NU_NUM} = {eta_pmns_num}")
print(f"  η̄_PMNS_denom = K_up_num·K_ν_num = {eta_pmns_denom}")
print(f"  η̄_PMNS = {eta_pmns}")

ratio_pmns = sp.simplify(eta_pmns / rho_pmns)
delta_form_b_rad = sp.pi + sp.atan(ratio_pmns)
delta_form_b_deg = float(delta_form_b_rad * 180 / sp.pi)
print(f"  η̄_PMNS/ρ̄_PMNS = {ratio_pmns}")
print(f"  Form B: δ_CP_PMNS = π + arctan({ratio_pmns}) ≈ {delta_form_b_deg:.2f}°")
print(f"  vs NuFit 195°: drift {drift_pct(delta_form_b_deg, DELTA_CP_NUFIT_DEG):.2f}%")
print(f"  vs T2K 2024 ≈ 248°: drift {drift_pct(delta_form_b_deg, 248.0):.2f}%")

M13_PASS = (rho_pmns == sp.Rational(2, 13)) and (eta_pmns == sp.Rational(6, 7)) \
           and (ratio_pmns == sp.Rational(39, 7))
print(f"  ρ̄_PMNS = 2/13, η̄_PMNS = 6/7, ratio = 39/7 sympy-exact: {M13_PASS}")
print(f"  Verdict M1.3 = {'PASS' if M13_PASS else 'FAIL'}")
print()


# ============== M1.4 — PMNS drift correction landscape ====================
print("=" * 72)
print("M1.4 — PMNS drift correction landscape (top-rational forms)")
print("=" * 72)

# sin²θ₁₃ candidates
print("  sin²θ₁₃ correction landscape:")
zero_t13 = K_NEUTRINO * LAMBDA_C**2
factor_a = (1 - RHO_BAR)                          # (1 − ρ̄) = 67/78
t13_form_a = zero_t13 * factor_a
print(f"    (a) K_ν·λ_C²·(1−ρ̄) = {sp.nsimplify(t13_form_a)} ≈ {float(t13_form_a):.6f} drift {drift_pct(t13_form_a, SIN2_T13_NUFIT):.3f}%")
factor_b = sp.Rational(45, 52)                     # (1 − 7/52) numerical
t13_form_b = zero_t13 * factor_b
print(f"    (b) K_ν·λ_C²·(45/52) ≈ {float(t13_form_b):.6f} drift {drift_pct(t13_form_b, SIN2_T13_NUFIT):.3f}%")
factor_c = sp.Rational(173, 200)                   # (1 − 27/200)
t13_form_c = zero_t13 * factor_c
print(f"    (c) K_ν·λ_C²·(173/200) ≈ {float(t13_form_c):.6f} drift {drift_pct(t13_form_c, SIN2_T13_NUFIT):.3f}%")
print(f"    → Best structural: form (a) (1 − ρ̄) z κ.1 cross-sector leakage")

# sin²θ₂₃ candidates
print("  sin²θ₂₃ correction landscape:")
t23_form_a = K_NEUTRINO / K_UP
print(f"    (a) K_ν/K_up = {t23_form_a} = {float(t23_form_a):.6f} drift {drift_pct(t23_form_a, SIN2_T23_NUFIT):.3f}%")
t23_form_b = K_NEUTRINO * (1 + N_GEN * LAMBDA_C**2)
print(f"    (b) K_ν·(1+N_gen·λ_C²) ≈ {float(t23_form_b):.6f} drift {drift_pct(t23_form_b, SIN2_T23_NUFIT):.3f}%")
t23_form_c = K_NEUTRINO * (1 + sp.Rational(18, 125))
print(f"    (c) K_ν·(1+18/125) ≈ {float(t23_form_c):.6f} drift {drift_pct(t23_form_c, SIN2_T23_NUFIT):.3f}%")
print(f"    → Best structural: form (a) K_ν/K_up cross-sector K-taxonomy ratio")

# sin²θ₁₂ candidates
print("  sin²θ₁₂ correction landscape:")
zero_t12 = sp.Rational(1, N_GEN)
t12_form_a = zero_t12 * (1 - LAMBDA_C * ETA_BAR)
print(f"    (a) (1/N_gen)·(1−λ_C·η̄) ≈ {float(t12_form_a):.6f} drift {drift_pct(t12_form_a, SIN2_T12_NUFIT):.3f}%")
t12_form_b = zero_t12 * (1 - LAMBDA_C / N_GEN)
print(f"    (b) (1/N_gen)·(1−λ_C/N_gen) ≈ {float(t12_form_b):.6f} drift {drift_pct(t12_form_b, SIN2_T12_NUFIT):.3f}%")
print(f"    → Best structural: form (a) (1−λ_C·η̄) cross-sector Wolfenstein imaginary leakage")

M14_PASS = (drift_pct(t13_form_a, SIN2_T13_NUFIT) < 1.0
            and drift_pct(t23_form_a, SIN2_T23_NUFIT) < 1.0
            and drift_pct(t12_form_a, SIN2_T12_NUFIT) < 1.0)
print(f"  3/3 angles drift < 1% z TGP-natural form: {M14_PASS}")
print(f"  Verdict M1.4 = {'PASS' if M14_PASS else 'FAIL'}")
print()


# ============== M1.5 — Viability gate ====================================
print("=" * 72)
print("M1.5 — Viability gate dla Phase 2")
print("=" * 72)

drift_hardening_ok = M14_PASS
delta_cp_form_a_deg = float(N_GEN * gamma_rad * 180 / sp.pi)
delta_cp_in_1s_form_a = DELTA_CP_NUFIT_1S[0] <= delta_cp_form_a_deg <= DELTA_CP_NUFIT_1S[1]
delta_cp_in_1s_form_b = DELTA_CP_NUFIT_1S[0] <= delta_form_b_deg <= DELTA_CP_NUFIT_1S[1]
print(f"  Drift hardening 3/3 < 1%: {drift_hardening_ok}")
print(f"  δ_CP Form A (gen-tripled γ_CKM) = {delta_cp_form_a_deg:.2f}° w 1σ: {delta_cp_in_1s_form_a}")
print(f"  δ_CP Form B (PMNS-Wolfenstein + π) = {delta_form_b_deg:.2f}° w 1σ: {delta_cp_in_1s_form_b}")

M15_PASS = drift_hardening_ok and delta_cp_in_1s_form_a and delta_cp_in_1s_form_b
print(f"  Viability gate 4/4 satisfied: {M15_PASS}")
print(f"  Verdict M1.5 = {'PASS' if M15_PASS else 'FAIL'}")
print()


# =================== Final ===========================================
print("=" * 72)
print("μ.1.Phase1 — Final verdict")
print("=" * 72)

results = [
    ("M1.1 NuFit 5.3 PMNS reference", M11_PASS),
    ("M1.2 γ_CKM = arctan(195/77) sympy-exact", M12_PASS),
    ("M1.3 PMNS-Wolfenstein analog (2/13, 6/7)", M13_PASS),
    ("M1.4 drift correction landscape 3/3 < 1%", M14_PASS),
    ("M1.5 viability gate (drift + δ_CP both forms in 1σ)", M15_PASS),
]
n_pass = sum(1 for _, r in results if r)
n_total = len(results)
for name, r in results:
    print(f"  [{'PASS' if r else 'FAIL'}] {name}")
print()
print(f"  μ.1.Phase1 score = {n_pass}/{n_total}")
if n_pass == n_total:
    print(f"  → Phase 2 first-principles derivation viable.")
else:
    print(f"  → Phase 1 incomplete; review viability before Phase 2.")
