---
title: "Phase 1 — F-FCR-A/B/C/D: frontier creation rate — STRUCTURAL_CONDITIONAL"
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-frontier-creation-rate-derivation-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'op-frontier-creation-rate-derivation'"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "8/8 PASS; 0 hardcoded T_pass=True"
falsifiers_resolved: "F-FCR-A A1 DERIVED (S/ρ̄ = H EXACT, conditional on B) / F-FCR-B B1 ε_G = 3/2 EXACT (STRUCTURAL_POSTULATE status) / F-FCR-C PARTIAL_concept_mismatch / F-FCR-D STRUCTURAL_CONDITIONAL"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — frontier creation rate: derivation + prediction matrix

## §0 — Verdict at a glance

**F-FCR-D = STRUCTURAL_CONDITIONAL** (mechanical per Phase 0 §1.3). Key structural results:

| Result | Status |
|---|---|
| **ε_G = (3/2)·Ω_m EXACT skeleton** (4πGρ̄t² with M ∝ t) | DERIVED (symbolic identity) |
| B1: R = c·t ∧ R = 2GM/c² ⇒ M = c³t/(2G); **ρ̄ = 3H²/(8πG) EXACT** (critical-density identity); ε_G = 3/2 EXACT | STRUCTURAL_POSTULATE (consistent, minimal, pre-existing roadmap task — not forced by LOCKED machinery) |
| M(t₀) = 8.80×10⁵² kg vs γ-7 rough input 10⁵³ | INFORMATIONAL consistency (ratio 0.88) |
| **A1: S_matter/ρ̄ = Ṁ/M = H EXACT** — concept paper §10.6 **hyp-Q3 RESOLVED POSITIVELY** (conditional on B) | DERIVED |
| A2: EQ-5 stationarity ⇒ S_Φ = 3H⟨Φ⟩ (∝ H consistent); Φ→matter bridge | GAP DECLARED |
| Native τ_init = 1/(1+z_rec) = 9.17×10⁻⁴ (γ-3 mapping 1+z = t₀/t EXACT) | DERIVED |
| γ-7/R17 τ_init = 2.75×10⁻⁵ provenance | **ΛCDM-borrowed (33×); R1 candidate flag #22** (forward-only; NO retroactive changes) |
| F-FCR-C bulk-form selection | **PARTIAL_concept_mismatch** (boundary-localized creation declared in EQ-5, bulk transport NOT specified) |

## §1 — Prediction matrix (bands LOCKED; G_obs = 10³ comparison-only)

| | C2a-form | C2b-form | **C2c-form (bulk-clean)** |
|---|---|---|---|
| **B1: Ω_m = 1 (zero-energy)**, ε = 3/2 | 10⁰·⁷ FAIL_LOW | no growth | **p = (√7−1)/2 EXACT → 10²·⁵⁰ PASS_BAND** ⭐ |
| **B2: Ω_m = 0.31 (E2 claim)**, ε = 0.465 | no growth | no growth | 10¹·⁰⁵ PARTIAL_BAND |

**The star cell (B1 × C2c):** a **parameter-free** number — log₁₀G = [(√7−1)/2]·log₁₀(1091) = 2.500
— from pure structure (zero-energy condition + source-free bulk + γ-3 mapping), landing inside the
PASS band, factor ~3 below observed 10³. No constant in this number is fitted: it contains only
the integers/roots of the indicial equation and z_rec (observational epoch anchor).

**Circularity guard (FP7):** all p₊ depend only on ε; G_obs absent from every derivation.
**ξ = 1 sensitivity (INFORMATIONAL, not adopted):** B1×C2c → 10³·⁹⁶; B1×C2a → 10²·²².

## §2 — What remains for PREDICTION_REALIZED (the missing-piece list, pre-registered)

1. **B1 derivation**: zero-energy/horizon condition from TGP substrate energetics (concept paper
   roadmap task "Derive Schwarzschild R_s z critical density Ω → 1" — pre-existing). Until then:
   STRUCTURAL_POSTULATE.
2. **C-form derivation**: bulk transport/homogenization of frontier-created matter (does
   boundary-localized creation leave sub-horizon bulk equations source-free? = C2c-form).
3. **A2 bridge**: Φ-substrate creation → matter-sector creation mapping.
4. (B2 alternative) E2-equilibrium Ω_m ≈ 0.31 verification (its own pre-registered F5).

If 1+2 (or 4+2) are derived, the growth factor becomes a genuine pre-registrable prediction →
PR-022 LOCK at THAT point. **PR-022 NOT appended now** (forbidden move #5 honored).

## §3 — R1 candidate #22 (NEW)

γ-7 Phase 3 (and R17 inheriting it) used τ_init = 2.75×10⁻⁵ = ΛCDM age-at-recombination
(1.2×10¹³ s), not the γ-3-native mapping τ = 1/(1+z) (→ 9.17×10⁻⁴; ratio 33×). Forward-looking
flag for ANY future cosmological cycle; γ-7 HALT-B and R17 ARTIFACT_PARTIAL verdicts UNCHANGED
(their internal regression logic is self-consistent; the flag affects native-prediction work only).
Severity: MEDIUM. Candidate §3.6.16 sub-rule: epoch mappings in TGP-native cycles must use γ-3
kinematics, not ΛCDM age tables.

## §4 — Predecessor invariance + anti-Lakatos (Phase 1): COMPLIANT ✓

0 predecessor verdicts modified (R17/P25/γ-7/F8/PR-017..020 PRESERVED). 0/8 hardcoded; routes and
bands LOCKED Phase 0; B-values restricted to pre-declared {1, 0.31}; ξ = 1 reported not adopted;
A2 bridge gap declared not improvised; PR-022 withheld; F8/Ω_DE untouched.

## §5 — Recommended next step

Phase FINAL closure: claim_status `CLOSED-RESOLVED STRUCTURAL_CONDITIONAL` + R1 #22 registration +
follow-up derivation targets (§2 list) as candidate cycles (e.g. `op-zero-energy-condition-derivation`,
`op-bulk-transport-frontier-matter`). Awaiting user reaction.
