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
| `C_2-defect` | **DEFINED** | 2-particle configuration space = ((R³ × RP²)² \ Δ)/S_2 | Phase 1 T1 |
| `γ_exchange` | **DEFINED** | exchange path: continuous deformation x₁↔x₂ z R(t) rotation | Phase 1 T6 |
| `[γ_exchange]` | **DERIVED-TARGET** | homotopy class element of π₁(C_2-defect) | Phase 1 T5 + T8 |
| `χ_exchange` | **DERIVED-TARGET** | exchange phase = exp(i·Berry phase) ∈ {+1, -1} | Phase 1 T8 + T9 |
| `Ψ(x_1, x_2)` | **STRUCTURAL** | 2-particle wavefunction, multi-component (spinor ⊗ spinor) | Phase 1 T9 |

## §2 — Inherited LOCKs (z poprzednich cykli)

| Inheritance | Status | Źródło |
|---|---|---|
| RP² ↔ S²/Z₂ projective target | LIVE | why_n3 PHASE3 + tgp_emergent_dirac_propagator §3 |
| π₁(RP²) = Z₂ topological fact | LIVE (standard math) | Hatcher "Algebraic Topology" Thm 1.16 (RP^n) |
| Single-defect Berry phase γ_Berry = π under 2π rotation | LIVE | why_n3 PHASE3 §2.4 (7/7 PASS) |
| Spin-1/2 transformation Ψ(2π) = -Ψ | LIVE | why_n3 PHASE3 §2.5 |
| S05 single-Φ axiom | LIVE | TGP_FOUNDATIONS §1 |
| m_obs vs M_full distinction | LIVE | op-L05-mass-exponent-k-alpha-d-2026-05-16 (CLOSED A−) |

## §3 — 6/6 P-requirements gate (target)

| P# | Requirement | Verification source (Phase 1) |
|---|---|---|
| P1 | C_2-defect configuration space topology explicit | T1-T4 sympy |
| P2 | Lifted γ_exchange ∈ π₁ Z₂ non-trivial loop | T5-T6 sympy |
| P3 | Berry connection additivity for two RP² defects | T7 sympy |
| P4 | Berry phase ∮ A = π along γ_exchange | T8 sympy |
| P5 | 2-particle wavefunction antisymmetric Ψ(x_1,x_2) = -Ψ(x_2,x_1) | T9-T10 sympy |
| P6 | Spin-statistics theorem verified; S05 preserved | T11 + T13 declarative |

## §4 — Risks (R-flags)

| R# | Risk | Mitigation |
|---|---|---|
| R1 | Configuration space topology — particle indistinguishability S_2 quotient | Explicit handling z S_2-orbit quotient construction (T1, T4) |
| R2 | Antipodal map on relative S² + RP² target orientation flip — coherent Z₂ effects | T3 + T4 explicit decomposition |
| R3 | Berry connection additivity — Aharonov-Bohm-like independent sectors | T7 verify additivity in tensor product Hilbert space |
| R4 | Exchange path uniqueness — multiple defect-orientation classes | T5 enumerate full π₁ structure; document if richer than Z₂ |

## §5 — Scope (in/out)

**IN:**
- 2-particle exchange in 3D flat configuration space
- RP² hedgehog defects from why_n3 Phase 3 (RP² topology + Berry phase π)
- Configuration space topology + π₁ computation
- Berry connection additivity + exchange path Berry phase
- Spin-statistics theorem consistency check

**OUT:**
- Full Dirac propagator iε prescription (deferred to multi-session cycle)
- N≥3 particle exchange braid statistics (3D trivially Sym + Antisym only; 2D braid groups outside scope)
- Defect core profile + stabilization radius (separate; m_eff(ψ) derivation outside scope)
- Connection to lepton flavor (e/μ/τ) — separate cycle
- SU(2) substrate extension (audit L08 path D) — different cycle if needed

## §6 — Native vs Standard physics split

Per CYCLE_KICKOFF_TEMPLATE §0.3:

- **Native part:** RP² hedgehog (z S05 single-Φ + Z₂ substrate); 2-defect configuration as
  product of TGP-emergent excitations
- **Standard part:** Configuration space topology (Hatcher); Finkelstein-Rubinstein argument
  (1968); π₁ + Berry phase calculations are universal mathematical tools
- **Mode:** Native-with-mapping — standard FR argument adapted to TGP-specific RP² defect;
  L2 reduction (FR original + spin-statistics theorem) is OPTIONAL last stage

## §7 — Validator gate compatibility

Validator `tooling/validate_kickoff.py` checks (mental verification):
- `contract::L1_native::output_observable` non-empty ✓ (χ_exchange + Pauli + spin-statistics)
- `contract::L1_native::measurement_instrument` non-empty ✓ (standard QM 2-particle Fock)
- `contract::L1_native::falsification_rule` pre-registered ✓ (with date 2026-05-16)
- `pre_registration_date` matches cycle date ✓
- 6/6 P-requirements declared ✓

**Status:** Phase 0 balance OK; proceed to Phase 1 symbolic derivation.

## Cross-references

- [[./README.md]] — kickoff contract
- [[./Phase1_sympy.py]] — symbolic derivation (next deliverable)
- [[../why_n3/PHASE3_RP2_defect_quantization.md]] — single-defect predecessor
- [[../../audyt/L08_kink_fermion_closure/README.md]] §1 problem 1
