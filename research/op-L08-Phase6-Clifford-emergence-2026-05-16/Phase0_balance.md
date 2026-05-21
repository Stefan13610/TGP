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
| `γ^a` (a=0,1,2,3) | **DEFINED** | 4×4 Dirac gamma matrices, chiral rep | Phase 1 T1 |
| `{γ^a, γ^b}` | **DERIVED-TARGET** | flat Cl(1,3) anticommutator = 2η^ab · 𝟙_4 | Phase 1 T2-T3 |
| `dim(min rep Cl(1,3))` | **DERIVED-TARGET** | = 4 (Dirac spinor) via Clifford classification | Phase 1 T4 |
| `γ^μ(x)` | **DERIVED-TARGET** | curved γ matrices = e_a^μ(x) γ^a (z M9.1'' tetrad) | Phase 1 T6 |
| `{γ^μ, γ^ν}` | **DERIVED-TARGET** | curved-space anticommutator = 2g^μν · 𝟙_4 | Phase 1 T7 |
| `D²` | **DERIVED-TARGET** | Dirac operator squared → KG dispersion | Phase 1 T8-T9 |
| `σ^ab = γ^[a γ^b]/2` | **DERIVED** | Spin(3,1) generators | Phase 1 T11 |

## §2 — Inherited LOCKs (z poprzednich cykli)

| Inheritance | Status | Źródło |
|---|---|---|
| M9.1'' metric A(ψ) = (4-3ψ)/ψ, B = 1/A | LIVE | op-emergent-metric-from-interaction-2026-05-09 |
| M9.1'' tetrad e^0_t = c_0·√A, e^a_i = √A·δ | LIVE | tgp_emergent_dirac_propagator.md §2 |
| RP² hedgehog → spin-1/2 transformation | LIVE | why_n3 PHASE3 (Berry phase π) |
| Antisymmetric 2-particle Fock {Ψ(x_1,x_2) = -Ψ(x_2,x_1)} | LIVE | op-L08-Phase6-FR-antisymmetry-2026-05-16 (sister cycle, A−) |
| S05 single-Φ axiom | LIVE | TGP_FOUNDATIONS §1 |
| Lorentz signature η = diag(-1,+1,+1,+1) | LIVE (via M9.1'') | TGP_FOUNDATIONS §2 + emergent-metric Phase 4 |

## §3 — 6/6 P-requirements gate (target)

| P# | Requirement | Verification source (Phase 1) |
|---|---|---|
| P1 | Cl(1,3) algebra explicit z {γ^a, γ^b} = 2η^ab | T1-T3 sympy |
| P2 | Minimal rep dim = 4 | T4 sympy |
| P3 | M9.1'' tetrad inheritance explicit | T5 sympy |
| P4 | Curved {γ^μ, γ^ν} = 2g^μν | T6-T7 sympy |
| P5 | Dirac² = KG dispersion | T8-T9 sympy |
| P6 | Connection to FR + S05 preserved | T10 + T13 declarative |

## §4 — Risks (R-flags)

| R# | Risk | Mitigation |
|---|---|---|
| R1 | γ representation choice (chiral/Dirac/Majorana) — physical results independent | Use chiral; document equivalence under conjugation; standard fact |
| R2 | Tetrad ψ-dependence: Cl algebra must hold point-wise | Phase 1 T7 verify {γ^μ(x), γ^ν(x)} = 2g^μν(x) for all ψ |
| R3 | Clifford anticommutator vs Fock anticommutator — same structural type, different domains | T10 explicit connection: spinor space vs particle space |
| R4 | Standard Cl construction inherits Lorentz; TGP inherits Lorentz z M9.1'' (geometry, not substrate Z₂) | Flag as inheritance: Cl algebra is consequence of Lorentz signature emergence in M9.1'', NOT separate derivation z Z₂ substrate |

## §5 — Scope (in/out)

**IN:**
- 4×4 Dirac γ matrices explicit (chiral representation)
- Cl(1,3) algebra {γ^a, γ^b} = 2η^ab anticommutator
- Minimal representation dim = 4 (Clifford classification)
- M9.1'' tetrad inheritance + curved γ^μ(x)
- Curved-space anticommutator {γ^μ, γ^ν} = 2g^μν
- Dirac operator D = iγ^μ ∂_μ - m_eff; D² Klein-Gordon
- Spin(3,1) generators σ^ab + spin-1/2 representation
- Consistency z FR Fock anticommutation (structural)

**OUT:**
- Full Dirac propagator iε prescription (downstream cycle)
- Connection to L05 m_obs/M_full (separate cycle)
- Spinor bilinears ψ̄γ^μψ → currents (deferred)
- Majorana/Weyl reduction (alternative reps; deferred)
- Connection to L08 problem #2 (e²/4 in mass) — separate cycle
- Connection to L08 problem #3 (quarks, gauge bosons) — multi-session

## §6 — Native vs Standard physics split

- **Native part:** Lorentz signature emergence z M9.1'' geometry (TGP-specific:
  A(ψ) = (4-3ψ)/ψ); tetrad structure; spin-1/2 from RP² topology
- **Standard part:** Cl(1,3) algebra structure {γ^a, γ^b} = 2η^ab is universal mathematical
  fact for Lorentz signature in 4D; γ representations standard QFT
- **Mode:** Native-with-mapping — standard Cl algebra inherited via TGP's Lorentz signature
  emergence; L2 reduction (Lounesto classification, standard QED limit) optional last stage

## §7 — Validator gate compatibility (mental check)

- `contract::L1_native::output_observable` non-empty ✓
- `contract::L1_native::measurement_instrument` non-empty ✓
- `contract::L1_native::falsification_rule` pre-registered ✓ (2026-05-16)
- 6/6 P-requirements declared ✓

**Status:** Phase 0 balance OK; proceed to Phase 1 symbolic derivation.

## Cross-references

- [[./README.md]] — kickoff contract
- [[./Phase1_sympy.py]] — symbolic derivation (next deliverable)
- [[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] — sister cycle
- [[../why_n3/tgp_emergent_dirac_propagator.md]] §7 — Dirac operator predecessor
