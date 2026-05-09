---
title: "Phase 0 — Balance sheet + scope"
date: 2026-05-09
parent: "[[README.md]]"
type: phase0-balance
phase: 0
---

# Phase 0 — Balance sheet + scope

## §0 — TL;DR

**Balance sheet status: BLOCKER ACTIVE — gravity sector predictions w limbo
post-2026-05-09 GWTC-3 falsification of M9.1'' (4−3ψ)/ψ form.**

This Phase 0 establishes:
1. Definitive falsification status (M9.1'' specific form RULED OUT)
2. Hard constraints inventory (C1-C10 sympy formalization plan)
3. Candidate parametrization families enumeration plan
4. Pre-cycle gate criterion: 8/8 PASS

## §1 — Falsification scope (definitive 2026-05-09)

### What IS falsified

| Claim | Source | Status |
|---|---|---|
| f(ψ) = (4−3ψ)/ψ specific form | M9.1'' canonical | **FALSIFIED 5.02σ** |
| g_tt = −c²(4−3ψ)/ψ specific metric | M9.1'' canonical | **FALSIFIED 5.02σ** |
| β_ppE^TGP^(b=-1) = -15/4 prediction | Phase 1.5 G_SPA=48 | **OBSERVATIONAL CONFIRMED**, **theory ruled out** |
| BF_TGP/GR ≈ 1 (Phase 2 original) | pre-Phase 1.5 | **REVISED: BF=3.55·10⁻⁶** |

### What is NOT falsified

| Claim | Status | Reason |
|---|---|---|
| TGP framework (Φ field + Z₂ + substrate) | NOT falsified | Single-field axiom independent of f(ψ) form |
| Dual-V framework (V_grav vs V_matter) | NOT falsified | Structural separation per dual-V clarification 2026-05-09 |
| Mass spectrum (m_μ/m_e, m_τ/m_e) | NOT falsified | A_tail mechanism V-independent (G.0 P22) |
| 1PN exact GR (γ=β=1) | NOT falsified | Holds for ANY f(ψ) satisfying C2 |
| Spin-1/2 (RP² Berry) | NOT falsified | V-independent topology |
| Φ_0 EFT scale-dependent | NOT falsified | EFT free parameter (analog Higgs VEV) |

### Critical distinction

**The FALSIFIED entity is the SPECIFIC f(ψ) = (4−3ψ)/ψ ANSATZ**, NOT the
M9.1'' philosophical framework (gravity emerges via metric ansatz on Φ
substrate). Alternative f(ψ) forms remain viable for exploration provided
they satisfy hard constraints C1-C10.

## §2 — Hard constraints inventory (C1-C10)

### C1 — α=2 vacuum Φ-EOM preserved

**Source:** T-D-uniqueness theorem (TGP foundation, sek08a):
K(ψ) = ψ⁴ from conformal substrate metric

**Mathematical statement:** The kinetic term in S_TGP must scale as ψ⁴ in
some equivalent formulation, dictating α=2 in field equations.

**Sympy testable?** YES — kinetic action variation check.

**Open question:** Czy K(ψ)=ψ⁴ pozostaje invariant pod f(ψ) update? Albo
musi być re-derived dla each candidate.

### C2 — 1PN exact GR

**Source:** G.0 P23 sympy LOCK 5/5 PASS (2026-05-02):
γ_PPN = β_PPN = 1 EXACT (NIE 1+ε perturbative).

**Mathematical statement:** Expanding f(ψ(U)) = 1 - 2U + 2U² + O(U³) at
isotropic radial coordinate U = GM/(rc²), and similarly for h(ψ).

**Constraint:**
- Coefficient of U at f(ψ) = -2 (gives γ_PPN = 1)
- Coefficient of U² at f(ψ) = +2 (gives β_PPN = 1)

**Sympy testable?** YES — Taylor expansion + master formula PPN reading.

### C3 — GWTC-3 2PN constraint

**Source:** [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]

**Constraint:** |β_ppE^(b=-1)| ≤ 0.78 (1σ) at η=1/4

**Implication via SPA chain:**
β_ppE = -(3/(128η)) · Δα_3 · G_SPA

At η=1/4:
β_ppE = -(3/32) · Δα_3 · G_SPA

|β_ppE| ≤ 0.78 ⟹ |Δα_3 · G_SPA| ≤ 0.78·(32/3) ≈ 8.32

**Sympy testable?** YES (SPA chain, modulo η-correction uncertainty).

### C4 — |Δα_3 · G_SPA| ≤ 8.32 (η=1/4)

**Derived from C3.** Provides the actionable bound on f(ψ) third-order
expansion coefficient + chain prefactor.

**Strategy:** if Δα_3 = 0 (1PN exact + 2PN matches GR exactly), then
constraint trivially satisfied regardless of G_SPA.

**Alternative:** if Δα_3 ≠ 0, then G_SPA(f) needs computation per candidate.

### C5 — Newton limit κ = 3/(4·Φ_0) preserved

**Source:** G.0 P32 (Phase 3 sympy LOCK 5/5):
q·c²/Φ_0 = (4/5)πG_0 ⟹ κ = 3/(4Φ_0_new) INVARIANT po re-fit.

**Mathematical statement:** Linear FRW perturbations + matter coupling
give FRW κ_eff = 4πG_0/(3H_0²) ≈ 3/(4Φ_0).

**Sympy testable?** YES — leading-order perturbation derivation.

### C6 — Mass spectrum (V-independent A_tail)

**Source:** G.0 P22 (sympy LOCK 5/5):
m_μ/m_e = 206.766 (PDG -0.0013%)
m_τ/m_e = 3477.40 (PDG +0.0049%)

**Key property:** A_tail mechanism uses ψ-soliton boundary behavior at
infinity; V-independent provided V → 0 at ψ → ψ_vac.

**Sympy testable?** YES (analytical, modulo numerical confirmation).

### C7 — Vacuum stability m_sp² > 0

**Source:** G.0 P21 sympy LOCK 4/5 PASS:
For V_M911 = -γψ²(4-3ψ)²/12, around ψ=1: m_sp² = +γ (stable).

**For alternative f(ψ):** new V_grav(ψ) ansatz must give m_sp² > 0 at
chosen vacuum point ψ_vac.

**Sympy testable?** YES — V''(ψ_vac) > 0 check.

### C8 — Hyperbolicity natural BH cutoff

**Source:** Structural M9.1'' philosophy:
ψ_h = 4/3 horizon coincides z g₀_crit (R3 critical, bariera ≡ horyzont).

**Form constraint:** f(ψ) → 0 (or → ∞ depending on sign convention) at
some finite ψ_h, signaling BH horizon emergence as natural pole/zero.

**Note:** this constraint is SOFT — alternative gravity emergence may
allow different BH horizon mechanism (entropy on Hubble sphere, Padmanabhan-
style). Phase 1 will treat C8 as DESIRABLE but not strictly REQUIRED.

### C9 — √(-g_eff) consistency

**Source:** sek08a (M9.1'' canonical: √(-g) = c·ψ/(4-3ψ)).

**Mathematical statement:** For chosen g_tt = -c²f(ψ), g_rr = h(ψ), the
4-volume element √(-g_eff) = c·√(f·h³) (in spherical coords with isotropic
radial). Constraint: this volume element must lead to UNIQUE V(ψ) from
R3-projection (G.0 derivation pattern).

**Sympy testable?** YES — derive V_grav unique form per f(ψ) candidate.

### C10 — Dual-V matter sector independence

**Source:** [[../op-dual-V-structure-clarification-2026-05-09/]]

**Mathematical statement:** V_grav (gravity sector ansatz) and V_matter
(scalar matter sector) are sector-tagged separately; modifications to
V_grav DO NOT affect V_matter or its derived predictions (Phase 5 Mach
inertia z γ → m_C² mapping unchanged).

**Sympy testable?** YES — sector-tag check via formal independence proof.

## §3 — Candidate parametrization families (Phase 1 enumeration plan)

### Family P (polynomial)

f(ψ) = a₀ + a₁ψ + a₂ψ² + a₃ψ³ + ... + a_N ψ^N

**Constraint anchors:**
- f(1) = 1 (vacuum at infinity matches GR)
- f'(1) = -2 (γ_PPN=1 from U-coefficient)
- f''(1) = 4 (β_PPN=1 from U²-coefficient — to be verified)

**Subcases:**
- P_min: f(ψ) = a₀ + a₁ψ + a₂ψ² (minimal degree-2)
- P_3: degree-3
- P_4: degree-4

### Family R (rational)

f(ψ) = P(ψ) / Q(ψ)

**Subcases:**
- R_(1,1): f = (a + bψ)/(c + dψ) [(4-3ψ)/ψ był specjalnym przypadkiem]
- R_(2,1): f = (a + bψ + cψ²)/(d + eψ)
- R_(2,2): f = quadratic/quadratic

### Family E (exponential)

f(ψ) = exp(g(ψ)) for some g(ψ) (analog Form II z S07 Diagnoza)

**Subcases:**
- E_lin: f = exp(-2(ψ-1)) (linear in ψ)
- E_quad: f = exp(-2(ψ-1) + b(ψ-1)²) (with 2PN tunable correction)

### Family H (hybrid)

f(ψ) = polynomial · exponential, or f = (polynomial)^β

### Family C (constraint-natural)

Start from C1-C10 algebraic system, derive forced form f(ψ) z constraint
solver. NOT empirical fitting — formal solution of constraint algebra.

## §4 — Strategic considerations

### S4.1 — Δα_3 = 0 strategy (preferred)

If candidate gives Δα_3 = 0 EXACTLY (i.e., GR-matching at U³ level), then
C3 (GWTC-3 2PN) trivially satisfied regardless of G_SPA.

**Implication:** Look for f(ψ) with f(ψ(U)) expansion = -2U + 2U² + 0·U³ + ...
(GR's coefficient at U³ is -7/3 hmm wait — actually GR is -2U + 2U² + (4/3)U³ + ...
TODO: clarify. The Δα_3 = -5/6 in M9.1'' was relative to GR, so absolute
α_3 in M9.1'' was -7/3 = -2.333, while GR has α_3 = -3/2 = -1.5? Or...

(To be verified in Phase 1 — GR isotropic Schwarzschild expansion baseline
must be locked via sympy.)

### S4.2 — Higher-order moments freedom

Alternative f(ψ) can deviate from M9.1'' at U⁴ and higher (3PN+) without
violating C3. Only Δα_3 (U³ = 2PN) directly constrained by GWTC-3.

### S4.3 — Constraint redundancy check

Some C_i may be **automatically** satisfied if others hold:
- C1+C9 → constraints on K(ψ), V(ψ) per f(ψ)
- C2+C5 → coupled (PPN + Newton)
- C7 → emergent from V derivation z C9

Phase 1 will identify minimal independent constraint set.

## §5 — Phase 0 gate criteria (8/8 needed)

| # | Criterion | Status |
|---|---|---|
| G1 | README.md scaffolding complete | ✅ |
| G2 | Phase0_balance.md (this document) created | ✅ (in progress) |
| G3 | C1-C10 inventory documented z source/sympy testability | ✅ (this §2) |
| G4 | Falsification scope clearly delineated | ✅ (this §1) |
| G5 | Candidate families enumerated (Phase 1 plan) | ✅ (this §3) |
| G6 | C1-C10 sympy formalization (Phase0_constraints_sympy.py) | PENDING |
| G7 | Strategic considerations documented | ✅ (this §4) |
| G8 | NEEDS.md inventory + cross-references | PENDING |

**Status:** 6/8 → continue with G6 (sympy script) and G8 (NEEDS).

## §6 — Cross-references

- [[README.md]] — cycle overview
- [[../../audyt/S07_M911_derivation/README.md]] — S07 audit source
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — G_SPA=48 LOCK (falsification source)
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] — 5.02σ ruled out
- [[../op-g0-r3-from-canonical-projection/]] — G.0 P21-P33 reference (M9.1'' V_M911 LOCK)
- [[../../TGP_FOUNDATIONS.md]] §3.5 — dual-V structure
