---
title: "op-L08-Phase6-FR-antisymmetry — Finkelstein-Rubinstein 2-particle exchange antisymmetry dla RP² hedgehog defects (spin-statistics closure)"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  L1_native:
    output_observable: "Two-particle wavefunction exchange phase χ_exchange [dimensionless ∈ {+1, -1}] for two RP² hedgehog defects; Pauli exclusion validity [boolean]; spin-statistics consistency [boolean] — operational construction emergentnego propagatora Diraca z antykomutacyjnymi własnościami (audit L08 problem #1)"
    measurement_instrument: "Standard QM 2-particle Fock space (Hilbert space tensor product z particle-exchange operator); inherits from why_n3 Phase 3 RP² topology + Berry phase derivation"
    native_coefs_constrained:
      - "χ_exchange = exp(iπ) = -1 dla RP² hedgehogs (FR antisymmetry)"
      - "Loop class [γ_exchange] ∈ π₁(C_2-defect) selected"
      - "Statistics phase = spin phase (spin-statistics consistency)"
    falsification_rule: "Jeśli symbolic derivation w π₁(C_2-defect) configuration space gives [γ_exchange] = trivial (Z₂ identity), wówczas RP² defects ARE BOSONIC, NOT fermionic — TGP kink-as-fermion roszczenie FALSIFIED structurally → wymaga (a) different defect topology (np. SU(2) skyrmion z full spinor target) lub (b) acknowledgment że kink-as-fermion jest niedomknięty (audit L08 path C: status hipotezy)."
    pre_registration_date: "2026-05-16"

  L2_framework_reduction:
    target_frameworks:
      - "Finkelstein-Rubinstein (1968) original argument dla SO(3) defects"
      - "Standard QM particle statistics (Bose/Fermi dichotomy in 3+1D)"
      - "Spin-statistics theorem (Pauli 1940 / Lüders-Zumino 1958)"
    reduction_type: "not-attempted"
    validation_transfer: ""
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "Spin-statistics theorem: spin-1/2 ↔ Fermi statistics (universal in 3+1D)", constrains: "χ_exchange = -1 dla spin-1/2 defects", window: "structural; verified across all relativistic QFTs", status: "pending Phase 1 derivation" }
    - { bound: "Pauli exclusion principle (atomic structure, NS interior, EOS)", constrains: "antisymmetric wavefunction for identical fermions", window: "universal; falsifiable if violated for any spin-1/2 system", status: "preserved if χ_exchange=-1 confirmed" }
    - { bound: "why_n3 Phase 3 Berry phase π (RP² topology) 7/7 PASS", constrains: "single-defect spinor transformation Ψ(2π)=-Ψ", window: "inherited; cycle extends to 2-particle exchange", status: "inherited PASS" }

tgp_status:
  level: L1
  kind: derivation
  output_type: structural-formula
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "Full Dirac propagator iε prescription (downstream cycle scope)"
    - "Three-particle exchange braid statistics (3D trivial; 2D braid; outside scope)"
    - "Generalization to multi-defect SU(N) chiral substrate (audit L08 path D)"
  depends_on:
    - "research/why_n3/PHASE3_RP2_defect_quantization.md (CLOSED 2026-05-01; RP² topology + Berry phase π)"
    - "research/why_n3/tgp_emergent_dirac_propagator.md §3-4 (projective defect + π₁ structure)"
    - "TGP_FOUNDATIONS §4 warstwa 3c (kink-as-fermion hipoteza)"
  impacts:
    - "L08 audit problem #1 (spin-statistics) — operational construction"
    - "L08 audit problem #4 (Dirac algebra) — partial: anticommutation structure for fermionic field operators"
    - "TGP_FOUNDATIONS §4 warstwa 3c — upgrade from (H) hipoteza toward (D) derived"
    - "research/why_n3 Phase 6+ — fundamental closure step"
  source_of_status:
    - "audyt/L08_kink_fermion_closure §1 problem 1 (spin-statistics)"
    - "research/why_n3/tgp_emergent_dirac_propagator.md §14 (central missing derivation)"

predecessors:
  - "[[../why_n3/PHASE3_RP2_defect_quantization.md]] (CLOSED; RP² → π₁=Z₂ → Berry π → spin-1/2)"
  - "[[../why_n3/tgp_emergent_dirac_propagator.md]] (§3-4 RP² defect; §14 central missing derivation)"
  - "[[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (CLOSED A−; m_obs vs M_full distinction LIVE)"

related:
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] §1 problem 1 (spin-statistics)"
  - "[[../../audyt/L05_mass_exponent_drift/README.md]] (closed 2026-05-16; bridge to L08 closure)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-2 (BINDING contract)"

classification: DERIVATION — L08 path A (operational Dirac propagator construction) sub-step
priority: high (P2 OPEN; address audit L08 problem #1 spin-statistics, the deepest open gap)
goal: "First-principles symbolic derivation of Finkelstein-Rubinstein 2-particle exchange antisymmetry for RP² hedgehog defects. Show that exchange path in configuration space C_2-defect generates non-trivial Z₂ loop ⇒ Berry phase π under exchange ⇒ wavefunction antisymmetric ⇒ Pauli exclusion. Closes L08 audit problem #1 (spin-statistics theorem) operationally — kink-as-fermion przesuwa się z 'roszczenia strukturalnego' do 'konstrukcji operacyjnej' (audit §1)."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 symbolic + Phase FINAL closure compressed)"
target_window: "Phase 1: 2-particle configuration space C_2-defect = ((R³ × RP²)² \\ Δ)/S_2 explicit topology; exchange path γ_exchange ⊂ C_2-defect homotopy class; Berry phase along γ_exchange = π via π₁ projection."

six_requirements_target:
  - "P1: Configuration space C_2-defect topology explicit z (R³)² \\ Δ relative position S² × R+, exchange ≃ antipodal map"
  - "P2: Lifted exchange path γ_exchange ∈ π₁(C_2-defect) Z₂ non-trivial loop verified"
  - "P3: Berry connection for two RP² defects compatible z Phase 3 single-defect A_φ (additive)"
  - "P4: Berry phase ∮_{γ_exchange} A = π derived symbolically"
  - "P5: 2-particle wavefunction Ψ(x_1, x_2) = exp(iπ)·Ψ(x_2, x_1) = -Ψ(x_2, x_1) (antisymmetric)"
  - "P6: Spin-statistics theorem verified: spin-1/2 (Phase 3) + antisymmetry (this cycle) ⇒ fermionic Fock space; S05 single-Φ preserved (single profile generates multi-defect via product)"

risk_flags:
  - "R1: Configuration space topology — need to handle particle indistinguishability correctly (S_2 quotient)"
  - "R2: Antipodal map on relative position S² combines z RP² target orientation flip — both Z₂ effects must add coherently (not cancel)"
  - "R3: Berry connection additivity for two defects — verify Aharonov-Bohm-like additivity holds for independent topological sectors"
  - "R4: Exchange path uniqueness — π₁(C_2-defect) might have more than Z₂ structure if multiple defect orientations matter; document if so"

phase_plan:
  Phase_0: "Balance sheet + 6/6 gate + scope (focus on FR antisymmetry; full Dirac propagator iε deferred)"
  Phase_1: "First-principles symbolic: 2-particle config space, exchange path, Berry phase derivation, antisymmetry"
  Phase_FINAL: "Closure + L2 reduction (FR original + spin-statistics theorem comparison) + L3 falsification map check + L08 audit closure note dla problem #1"

tags:
  - L08
  - L08-Phase6
  - emergent-dirac
  - spin-statistics
  - finkelstein-rubinstein
  - antisymmetry
  - first-principles
  - audit-L08-problem-1
  - cycle-scaffold-2026-05-16
---

# op-L08-Phase6-FR-antisymmetry-2026-05-16

> **Cel:** First-principles symbolic derivation of Finkelstein-Rubinstein 2-particle
> exchange antisymmetry dla RP² hedgehog defects. Closes L08 audit problem #1
> (spin-statistics theorem) operationally. Kink-as-fermion: roszczenie strukturalne →
> konstrukcja operacyjna.

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; audit L08 §1 (problem 1 spin-statistics)]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- `χ_exchange` [dimensionless ∈ {+1, -1}] — phase factor under 2-particle exchange
- `[γ_exchange]` [element of π₁(C_2-defect) = Z₂] — homotopy class
- `Ψ(x_1, x_2)` vs `Ψ(x_2, x_1)` — 2-particle wavefunction symmetry property

**Instrument:** standard QM 2-particle Fock space + symbolic computation; inherits why_n3
Phase 3 RP² topology + Berry phase calculation framework

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli symbolic derivation w π₁(C_2-defect) configuration space gives `[γ_exchange] = trivial`
> (Z₂ identity element), then RP² defects ARE BOSONIC, NOT fermionic — TGP kink-as-fermion
> roszczenie FALSIFIED structurally → wymaga (a) different defect topology (np. SU(2)
> skyrmion z full spinor target) lub (b) acknowledgment że kink-as-fermion jest niedomknięty
> (audit L08 path C: status hipotezy permanentny).

```
pre_registration_date: 2026-05-16
recovery_scope:
  allowed_directions:
    - "Subleading corrections to Berry connection (do not change phase modulo 2π)"
    - "Alternative parametrizations of exchange path γ (homotopy-equivalent)"
  forbidden_directions:
    - "Free parameter in χ_exchange (must be exact ±1 from topology)"
    - "Post-hoc lifting of Z₂ to U(1) (would change classification)"
    - "Adding fundamental fermionic field on top (S05 violation)"
  if_recovery_exhausted: "Honest verdict: kink-as-fermion conjecture not derivable; audit L08 path C remains preferred outcome"
```

### §0.3 — TGP-native check (mandatory)

- [x] **Q1 (Pattern coverage):** Pattern 2.8 (topological emergence) relevant
- [x] **Q2 (Red flags):** NONE — standard topology + Berry phase, established mathematics
- [x] **Q3 (Inherited LOCKs):** RP² topology + π₁=Z₂ + Berry phase π (Phase 3) LIVE
- [x] **Q4 (Standard-physics tools):** Finkelstein-Rubinstein argument (1968) is standard;
      native-relevance: applied to TGP-specific RP² hedgehog (S05 single-Φ origin)
- [x] **Q5 (m_Φ usage):** N/A — topological argument, no field mass involved
- [x] **Q6 (GR limit framing):** N/A — flat space configuration space
- [x] **Q7 (ASK-RULE self-check):** methodology cited; gaps in R-flags
- [x] **Q8 (BD-drift audit plan):** Phase FINAL flag if applicable

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2
- [x] Przeczytano [[../../audyt/L08_kink_fermion_closure/README.md]]
- [x] Przeczytano [[../why_n3/PHASE3_RP2_defect_quantization.md]] (predecessor)
- [x] Przeczytano [[../why_n3/tgp_emergent_dirac_propagator.md]] §3-§4, §14

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16

### §0.5 — Sympy substance plan

**Plan testów Phase 1 (target ~10-12 tests, ≥8 first-principles):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | C_2-defect = ((R³ × RP²)² \\ Δ)/S_2 configuration space explicit |
| T2 | **FIRST_PRINCIPLES** | Relative position factorization: R³² \\ Δ ≃ R³_CM × (R³ \\ {0}) ≃ R³_CM × S² × R+ |
| T3 | **FIRST_PRINCIPLES** | Exchange operation (x₁ ↔ x₂) ⇔ antipodal map on relative-position S² |
| T4 | **FIRST_PRINCIPLES** | Combined config space relative direction: S²/Z₂ × RP² × RP² (additional Z₂ from indistinguishability) |
| T5 | **FIRST_PRINCIPLES** | π₁ calculation: π₁(C_2-defect) = Z₂ ⊕ Z₂ ⊕ Z₂ structure (exchange + 2 individual orientations) |
| T6 | **FIRST_PRINCIPLES** | Exchange path γ_exchange: continuous deformation x_1(t)=R(t)e₁, x_2(t)=R(t)e₂ z R(t) rotation matrix |
| T7 | **FIRST_PRINCIPLES** | Berry connection 2-particle = A_1 ⊗ I + I ⊗ A_2 (additive for independent defects) |
| T8 | **FIRST_PRINCIPLES** | Berry phase along γ_exchange: ∮ A = π z half-twist of relative S² |
| T9 | **FIRST_PRINCIPLES** | 2-particle wavefunction transformation: Ψ(x_1,x_2) = exp(iπ)·Ψ(x_2,x_1) = -Ψ(x_2,x_1) |
| T10 | **FIRST_PRINCIPLES** | Pauli exclusion: dla x_1 → x_2, Ψ(x,x) = -Ψ(x,x) ⇒ Ψ(x,x) = 0 (impossibility of same state) |
| T11 | **FIRST_PRINCIPLES** | Spin-statistics theorem verification: spin-1/2 (Phase 3) + antisymmetry (T9) ⇒ Fermi-Dirac statistics consistent z Pauli (1940) |
| T12 | **LITERATURE_ANCHORED** | Comparison z Finkelstein-Rubinstein (1968) "Connection Between Spin Statistics and Kinks" — same SO(3) twist mechanism (sigma model precedent) |
| T13 | **DECLARATIVE** | S05 single-Φ preservation — multi-defect z product of single profiles (separate, NIE w PASS total) |

**Target:** 12/12 PASS sympy (T1-T12) + 1 structural declaration (T13 separate).

**Ratio:** 11 FIRST_PRINCIPLES (91.7%) + 1 LITERATURE_ANCHORED (8.3%) + 1 DECLARATIVE separate.

---

## §1 — Phase 0: balance sheet

[Patrz `Phase0_balance.md`]

## §2 — Phase 1: native derivation

[Patrz `Phase1_sympy.py` + `Phase1_results.md`]

## §FINAL — Closure

[Patrz `Phase_FINAL_close.md`]

---

## Status

🟢 **ACTIVE — opened 2026-05-16** per user authorization "ok działaj z L08 op-why_n3-Phase6-dirac".

Cycle scope: focused FR antisymmetry derivation (audit L08 problem #1, the deepest gap);
full Dirac propagator iε prescription deferred to downstream multi-session cycle.

This session deliverables:
- README.md (this file) z BINDING contract — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (12-test first-principles symbolic)
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED if Phase 1 PASS**

---

**Cross-references:**
- [[../../audyt/L08_kink_fermion_closure/README.md]] §1 problem 1 (spin-statistics)
- [[../why_n3/PHASE3_RP2_defect_quantization.md]] (Berry phase π; single-defect)
- [[../why_n3/tgp_emergent_dirac_propagator.md]] §14 (central missing derivation)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (predecessor closure note)
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]
- [[../../STATE.md]]
