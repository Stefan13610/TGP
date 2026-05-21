---
title: "Phase 0 — Balance sheet: emergent Dirac propagator z FR + Cl + L05 inheritance"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-balance
phase: 0
gate_status: 8/8 ☑ PASS
binding: "Phase 6 ABSOLUTE BINDING gate (CALIBRATION_GATE_ENFORCEMENT)"
tags:
  - balance-sheet
  - phase0
  - pre-derivation
  - L08-Dirac-propagator
  - inheritance-ledger
---

# Phase 0 — Balance sheet (8/8 ☑ gate)

> **MANDATORY per [[../../meta/CALIBRATION_PROTOCOL.md]] §2** + Phase 6 ABSOLUTE BINDING gate.
> **All 8 checks ☑ PASS** wymagane przed Phase 1 entry.

## ☑ Check 1 — External inputs (lit data declared)

**External literature inputs declared (LIT-grade):**

| Input | Source | Used in Phase 1 (planned) |
|---|---|---|
| PDG 2024 m_μ/m_e = 206.7682275 | Particle Data Group | T7 pole mass validation (m_obs identification) |
| PDG 2024 m_τ/m_e = 3477.23 | Particle Data Group | T7 (inherited z L05) |
| Bjorken-Drell convention Dirac propagator S_F = i(γ·p+m)/(p²-m²+iε) | Bjorken-Drell 1964 §15.3 | T1 reference form (L2 reduction target) |
| Finkelstein-Rubinstein 1968 exchange phase | J. Math. Phys. 9, 1762 | already consumed via L08-FR LIVE |
| Cl(1,3) algebra {γ^μ, γ^ν} = 2η^μν | textbook (Peskin-Schroeder §3.2) | already consumed via L08-Clifford LIVE |
| CPT theorem (Lüders-Pauli) | standard | T11 declarative (S05 + Lorentz emergent → CPT) |

**TGP-internal LOCKED inputs (z LIVE inheritance):**

| Inheritance | Source cycle | Status |
|---|---|---|
| FR antisym γ_exchange = π Berry phase | [[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] | A− LIVE |
| Cl(1,3) emergence z M9.1'' Lorentz signature | [[../op-L08-Phase6-Clifford-emergence-2026-05-16/]] | A− LIVE |
| m_obs ≠ M_full distinction; k_obs(α=2,d=3)=3 | [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] | A− LIVE |
| why_n3 Phase 5 universal mass formula m=c·A²·g_0^(e²(1-α/4)) | [[../why_n3/]] (2026-05-01) | CLOSED LIVE |
| S05 single-Φ axiom | [[../../TGP_FOUNDATIONS.md]] §1 | LOCKED |
| L01 ρ ≡ -T^μ_μ/c_0² | [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] | LOCKED |
| M9.1'' canonical form (Lorentz signature inheritance) | sek08a/sek08c + S07 audit caveat | LOCKED |

**Status:** ☑ All external inputs declared explicit; all LIVE inheritance enumerated.

## ☑ Check 2 — Structural axioms TGP-internal LOCKED

**S05 single-Φ axiom:** No new fundamental fields. Emergent fermion = bound state of
2 kinki w Φ field. Verified declaratively w Phase 1 T13.

**ax:metric-coupling (S04, B9 6/6 PASS):** Universal coupling g_eff[Φ] to wszystkich
matter species. Preserved bezwarunkowo (B9 MICROSCOPE 8.3·10¹⁰× safe).

**Z₂ chirality axiom (FOUNDATIONS §1):** Substrate Z₂ symmetry — inherited via FR
+ Cl cycles. Phase 1 T11 verifies CPT operational.

**M9.1'' metric canonical:** Inherited z S07 caveat — Lorentz signature dla Cl algebra
emergence proper.

**Status:** ☑ All TGP-internal axioms enumerated + LOCKED status confirmed.

## ☑ Check 3 — Derived outputs (target observables)

**Phase 1 target outputs:**

1. **S_F^TGP(p)** explicit form: `i(γ·p + m_obs)/(p² - m_obs² + iε)`
2. **m_obs pole identification:** k_obs(α=2, d=3) = 3 LIVE z L05
3. **Lorentz boost covariance:** S_F(p) → S(Λ) S_F(Λp) S(Λ)⁻¹ z Cl emergence
4. **Free-field limit:** lim_{Φ→Φ_0} S_F^TGP = S_F^Dirac (canonical)
5. **CPT consistency:** S_F^†(p) = γ_0 S_F(p) γ_0 (Hermiticity)

**Phase 2 (optional) target outputs:**
- Wilson coefs dla curvature × propagator corrections (Phase 2 deferred)
- PDG a_e precision comparison (Phase 2 deferred)
- Lamb shift δ_TGP[m_obs] (Phase 2 deferred)

**Status:** ☑ All derived outputs operationalized z explicit observables.

## ☑ Check 4 — Tautology test (czy NIE assume conclusion)

**Pytanie:** Czy cykl NIE udowadnia założenia jako wnioskowania?

**Test:**
- ❌ NIE zakładamy explicit Dirac equation a priori; derive z FR + Cl + L05
- ❌ NIE zakładamy {γ^μ, γ^ν} = 2η^μν directly — inherited z L08-Clifford derivation
- ❌ NIE zakładamy antisymmetry directly — inherited z L08-FR Berry phase π derivation
- ❌ NIE zakładamy m_obs as pole mass postulate — derivation z L05 universal mass formula
- ✅ Pole identification (p² = m_obs²) jest wnioskiem z propagator construction +
  L05 inheritance (NIE postulatem)
- ✅ Lorentz covariance jest wnioskiem z Cl algebra inheritance (NIE postulatem)

**Status:** ☑ NO tautology — cycle constructs S_F^TGP z LIVE inheritance, NIE assumes.

## ☑ Check 5 — Falsifiability test (explicit gates)

**Falsification gates per Phase 1 (BINDING):**

| Sub-test | Falsifies if... | Recovery? |
|---|---|---|
| T1 (2-particle Fock) | Antisymmetry FR inheritance NIE applicable | HALT-B: FR scope insufficient |
| T2 (γ^μ inheritance) | Cl algebra NIE consistent z propagator structure | HALT-B: Cl emergence cycle gap |
| T3-T4 (propagator construction) | Pole NIE at p² = m_obs² | HALT-B: structural obstruction |
| T5-T7 (m_obs identification) | k_obs(α=2,d=3)=3 NIE matches pole | HALT-B: L05 distinction insufficient |
| T8-T9 (Lorentz covariance) | Boost NIE consistent z Cl emergence | HALT-B: M9.1'' Lorentz gap |
| T10 (free-field limit) | NIE reduces do canonical Dirac | HALT-B: limit obstruction |
| T11 (CPT) | Substrate Z₂ + Lorentz NIE → CPT | HALT-B: discrete symmetry gap |
| T12 (Hermiticity) | S_F^†(p) ≠ γ_0 S_F(p) γ_0 | HALT-B: unitarity gap |
| T13 DEC (S05) | NIE preserved | CRITICAL: cycle abandon, S05 violation |

**Pre-registered falsification rule (per README §0.2):** LOCKED. No post-hoc weakening.

**Status:** ☑ Explicit falsification gates per sub-test + global rule LOCKED.

## ☑ Check 6 — Independent path test

**Pytanie:** Czy cycle używa ≥3 niezależnych path/anchor sources?

**Independent inheritance sources (NIE z jednego cyklu):**

1. **L08-FR** (antisymmetry; problem #1 closure) — derivation z FR 1968 + Berry transport
2. **L08-Clifford** (algebra; problem #4 closure) — derivation z M9.1'' Lorentz signature
3. **L05** (m_obs ≠ M_full; mass exponent reconciliation) — derivation z Derrick virial + Sobolev
4. **why_n3 Phase 1-5** (kink-fermion strukturalny szkic) — earlier 2026-05-01 closure
5. **PDG 2024** — independent experimental anchor (m_μ/m_e, m_τ/m_e)
6. **Bjorken-Drell convention** — independent QFT textbook reference

**Status:** ☑ 4+ niezależnych inheritance sources + 2 external anchors → robust.

## ☑ Check 7 — Risk assessment (R-flags)

**Identified risks (pre-Phase-1):**

| R# | Risk | Mitigation |
|---|---|---|
| **R1** | FR antisymmetry + Cl algebra MAY not combine smoothly w propagator construction | Phase 1 T1-T4 explicit construction; if obstruction → HALT-B documented honestly |
| **R2** | m_obs identification z k_obs(α=2,d=3)=3 LIVE może wymagać dimensional rescaling | Phase 1 T5-T7 explicit unit analysis; coupling c_M z why_n3 inheritance |
| **R3** | Lorentz boost generator construction z M9.1'' Cl emergence — non-trivial | Inherited z L08-Clifford; if gap → consult predecessor cycle scope |
| **R4** | Free-field limit Φ → Φ_0 might require subtle V_eff regularization | Phase 1 T10 explicit limit check; canonical Dirac as L2 target |
| **R5** | CPT theorem operational claim — discrete Z₂ + emergent Lorentz interaction subtle | Phase 1 T11 DECLARATIVE check; cite Lüders-Pauli theorem reference |
| **R6** | Phase 2 Wilson coefs (a_e correction) may diverge | Phase 2 OUT-OF-SCOPE this cycle; documented honestly |

**Status:** ☑ 6 R-flags identified + per-flag mitigation; NO showstopper pre-Phase-1.

## ☑ Check 8 — Honesty test (no hardcoded T_pass=True)

**Phase 1 sympy code commitment:**

```python
# Each test computes substantive result, then VERIFIES against expected/derived.
# NEVER:
#   T_pass = True
# ALWAYS:
#   T_pass = (sp.simplify(expr_lhs - expr_rhs) == 0) and (some other condition)
```

**Hardcoded `T_pass = True` policy:** **0 allowed** (Phase 6 ABSOLUTE BINDING gate).

**Substance metric target (Phase 1):**
- FIRST_PRINCIPLES: ≥70% (target 11/13 ≈ 84.6%)
- LITERATURE_ANCHORED: ≤30% (LIT inputs: PDG ratios, Bjorken-Drell convention reference)
- DECLARATIVE (separate): 1 (T13 S05)
- Hardcoded `T_pass = True`: **0**

**Pre-registration:**
- All 6 P-requirements stated pre-Phase-1 (README §0.3)
- Falsification rule LOCKED pre-Phase-1 (README §0.2)
- Pre-registration date 2026-05-16 (immutable)

**Status:** ☑ Phase 6 ABSOLUTE BINDING gate enforced; 0 hardcoded `T_pass = True` policy.

---

## §9 — Gate verdict

| Check | Status |
|---|---|
| 1. External inputs declared | ☑ PASS |
| 2. TGP-internal axioms LOCKED | ☑ PASS |
| 3. Derived outputs operationalized | ☑ PASS |
| 4. NO tautology | ☑ PASS |
| 5. Falsifiability gates explicit | ☑ PASS |
| 6. Independent path test (≥3 sources) | ☑ PASS (4+ sources) |
| 7. Risk assessment (R1-R6) | ☑ PASS |
| 8. Honesty test (0 hardcoded) | ☑ PASS |

**Phase 0 gate verdict:** 🟢 **8/8 ☑ PASS** — Phase 1 ENABLED.

**Phase 1 activation pre-condition:** User authorization explicit (e.g.,
"aktywuj L08-Dirac Phase 1" or "/aktywuj fazę 1").

## §10 — Phase 1 sympy plan (preview)

Planowane 13 testów:

| # | Test | Type | Inheritance |
|---|---|---|---|
| T1 | 2-particle Fock state z exchange antisymmetry | FP | L08-FR |
| T2 | γ^μ operator construction z Cl(1,3) inheritance | FP | L08-Clifford |
| T3 | Dirac equation operator (i γ·∂ - m) Ψ = 0 | FP | T1+T2 |
| T4 | Propagator construction S_F = i(γ·p + m_obs)/(p² - m_obs² + iε) | FP | T3 |
| T5 | Pole identification at p² = m_obs² | FP | T4 + L05 |
| T6 | m_obs z why_n3 Phase 5 universal mass formula | FP | L05 + why_n3 |
| T7 | PDG m_μ/m_e = 206.7682 ratio verification | LIT | T6 + PDG |
| T8 | Lorentz boost covariance | FP | L08-Cl + M9.1'' |
| T9 | Hermiticity S_F^†(p) = γ_0 S_F(p) γ_0 | FP | T4 + T8 |
| T10 | Free-field limit lim_{Φ→Φ_0} S_F^TGP = S_F^Dirac | FP | T4 + Bjorken-Drell |
| T11 | CPT theorem operational (declarative) | DEC | Substrate Z₂ + emergent Lorentz |
| T12 | Two-point function ⟨0|Tψ(x)ψ̄(y)|0⟩ = S_F(x-y) | FP | T4 + Wick rotation |
| **T13** | **DECLARATIVE S05 single-Φ axiom preservation** | **DEC** | FOUNDATIONS §1 |

**Substance target:** 11 FP + 1 LIT + 1 DEC = 11/13 FP (84.6%); 0 hardcoded.

## §11 — Cross-references

- [[./README.md]] §0 BINDING contract
- [[../../meta/CALIBRATION_PROTOCOL.md]] §2 Phase 0 mandatory
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] Phase 6 ABSOLUTE BINDING
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] BINDING contract template
- [[../op-L08-Phase6-FR-antisymmetry-2026-05-16/Phase0_balance.md]] (predecessor template)
- [[../op-L08-Phase6-Clifford-emergence-2026-05-16/Phase0_balance.md]] (predecessor template)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/Phase0_balance.md]] (predecessor template)
