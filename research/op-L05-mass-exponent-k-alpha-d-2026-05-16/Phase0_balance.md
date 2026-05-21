---
title: "Phase 0 — Balance sheet + 6/6 gate + scope"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟢 BALANCE OK — proceed to Phase 1
---

# Phase 0 — Balance sheet

## §1 — Wprowadzane wielkości (this cycle)

| Symbol | Status | Definicja | Źródło |
|---|---|---|---|
| `k_full(α, d)` | **DERIVED-TARGET** | exponent w M_full ∝ A^k_full (volumetric integral) | Phase 1 T4 (sympy) |
| `k_obs(α, d)` | **DERIVED-TARGET** | exponent w m_obs ∝ A_tail^k_obs (tail coupling) | Phase 1 T9 (sympy) |
| `σ_match(α, d)` | **DERIVED-TARGET** | A_tail ~ A^σ_match exponent (core-tail matching) | Phase 1 T8 (sympy) |

## §2 — Inherited LOCKs (z poprzednich cykli)

| Inheritance | Status | Źródło |
|---|---|---|
| `K(φ) = K_geo · φ^α` z α=2 canonical | LIVE | L04 thm:D-uniqueness (op-L04-ODE-canonicalization-2026-05-04) |
| `V_eff(φ) ~ λ·φ⁴/4` z TGP-FOUNDATIONS §3 | LIVE | TGP_FOUNDATIONS.md §3.5 |
| `p(α) = 5 − α` numerical (α∈{1,2} EXACT) | NUMERICAL_INPUT | why_n3/r3_observable_vs_full_mass.py |
| LP-4 `k = 4` dla K=g², d=3 | LITERATURE_PHASED | mass_scaling_k4/lp4_mass_exponent_verification.py 9/9 PASS |
| Lepton mass ratios m_μ/m_e, m_τ/m_e | LIVE (PDG) | PDG 2024 + r3_alpha2_full_closure.py |

## §3 — 6/6 P-requirements gate (target)

| P# | Requirement | Verification source (Phase 1) |
|---|---|---|
| P1 | k_full(α, d) z virial scaling | T1-T5 sympy |
| P2 | k_obs(α, d) z asymptotic matching | T6-T9 sympy |
| P3 | Reconciliation theorem: LP-4 (M_full) AND R3 (m_obs) | T10-T11 sympy |
| P4 | m_obs vs M_full distinction operationally formalized | T4 vs T9 distinct exponents |
| P5 | L05 Możliwość A confirmed constructively | T10-T11 verify α=2→p=3, α=1→p=4 |
| P6 | S05 single-Φ preserved | T13 declarative |

## §4 — Risks (R-flags)

| R# | Risk | Mitigation |
|---|---|---|
| R1 | Virial scaling depends on V form — alternative V may give different k_full | Use canonical λφ⁴/4; document if testing alternative V |
| R2 | Asymptotic tail expansion — higher-order corrections | Verify next-to-leading O(1/r²) terms don't change k_obs |
| R3 | Sobolev critical exponent connection — flag if structural | PASS-WITH-FLAG annotation if discovered |
| R4 | numerical p(α)=5-α EXACT only dla α∈{1,2}; intermediate α has 2-3% deviation | Phase 1 derivation may explain via higher-order finite-size corrections; document scope |

## §5 — Scope (in/out)

**IN:**
- Radial soliton φ(r) in d=3 (primary) with parametric d-extension
- K(φ) = K_geo · φ^α z α ∈ {1, 2} primary, general α parametric
- V_eff(φ) = canonical quartic ~ φ⁴ (TGP-FOUNDATIONS §3)
- Derrick virial scaling (Phase 1 T2-T4)
- Linearized asymptotic tail (Phase 1 T6-T9)
- Core-tail matching condition (Phase 1 T8)

**OUT:**
- Non-radial profiles (rotation, deformation)
- Higher-order V (φ⁶, φ⁸ etc.) — separate cycle
- Quantum corrections (1-loop renormalization of mass formula)
- Connection to L08 emergent Dirac propagator (downstream cycle scope)
- Numerical N-body verification (inherited from why_n3 numerical data)

## §6 — Native vs Standard physics split

Per CYCLE_KICKOFF_TEMPLATE §0.3 — three-mode taxonomy:

- **Native part:** K(φ) = K_geo · φ^α form jest TGP-specific (S05 single-Φ + TGP-FOUNDATIONS §3 emergent
  scalar substrate); standard Derrick analysis assumes K=const
- **Standard part:** virial scaling + asymptotic tail linearization is universal variational physics
- **Mode:** Native-with-mapping — standard analysis adapted to TGP K(φ)=φ^α; L2 reduction (Sobolev/Derrick)
  is OPTIONAL last stage

## §7 — Validator gate compatibility

Validator `tooling/validate_kickoff.py` checks:
- `contract::L1_native::output_observable` non-empty ✓ (mass ratios + exponents)
- `contract::L1_native::measurement_instrument` non-empty ✓ (PDG 2024 + why_n3 numerical)
- `contract::L1_native::falsification_rule` pre-registered ✓ (with date 2026-05-16)
- `pre_registration_date` matches cycle date ✓
- 6/6 P-requirements declared ✓

**Status:** Phase 0 balance OK; proceed to Phase 1 symbolic derivation.

## Cross-references

- [[./README.md]] — kickoff contract (full)
- [[./Phase1_sympy.py]] — symbolic derivation (next deliverable)
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §2.6 (Phase 0 template)
- [[../../audyt/L05_mass_exponent_drift/README.md]] §Rekomendacja (Phase plan source)
