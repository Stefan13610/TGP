---
title: "Phase 2 setup — S07 alternative f(ψ): symbolic Bayesian α-mapping + family distinguishability"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 2
predecessor: "[[./Phase1_results.md]] — 12/12 PASS, β_ppE^poly(α) = (15/16)·α LINEAR SCALING; recovery region α ∈ [-0.832, 0.832]"
authorization: "user '/autoryzuje opcja A' 2026-05-13 conversation"
status: 🟡 ACTIVE — sympy + results pending
---

# Phase 2 setup — S07-reset alternative f(ψ)

## §0 — Pre-flight methodology re-confirmation

Per BINDING workflow (CYCLE_KICKOFF_TEMPLATE.md §1-§2 + §0.4 mandatory pre-flight):

- [x] Phase 1 closed: 12/12 PASS, 10 FP / 2 LIT / 2 DEC, 0 hardcoded; substantive linear-scaling discovery
- [x] PR-010 LOCKED-PENDING-PHASE-1 → przechodzi do LOCKED-PHASE-2-IN-PROGRESS
- [x] Pre-registered falsification rule (PR-010) IMMUTABLE — Phase 2 NIE modyfikuje recovery_scope ani decision_rule
- [x] M9.1'' inheritance reframed per `meta/M9_RESTRUCTURE_NOTE.md` §3 — Path 2 anchor (specific point z {a_n^M911, c_0·κ_σ=4/3}), NIE canonical metric
- [x] Three-layer L1/L2/L3 presentation MANDATORY per `meta/PPN_AS_PROJECTION.md` §3.1 (β_ppE = L2 projection, NIE primary native)
- [x] Anti-BD-drift Triggers A-D executed per `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1.1 — see §0.1 below

## §0.1 — ASK-RULE Triggers A-D execution (anti-BD-drift)

| Trigger | Status | Notes |
|---|---|---|
| **A — TGP analogy visible?** | ✅ OK | Bayesian α-mapping referuje L1 native Taylor coefs {a_n, ξ_n, c_0, κ_σ} z LIGO-3G-native A− cycle; ppE basis cited explicitly jako L2 projection chart. |
| **B — Predecessor LOCK inheritance audit** | ✅ OK z conditional citation | Inheriting: c_0·κ_σ = 4/3 z `op-c0-derivation`+`op-kappa-sigma-2body-PN`+`op-emergent-metric` Phase 4 LOCK; β_ppE = (45/16)·Δe_2_native z `op-LIGO-3G-native-phase-residual` A−; β_ppE^poly = (15/16)·α z Phase 1. **Klasyfikacja per M9_RESTRUCTURE §1.4: kategoria (b)** — Path 2 anchor citation conditional. NIE kategoria (c) ani (d). |
| **C — Reproducing literature without TGP mechanism?** | ✅ OK | GWTC-3 ppE Bayesian posterior (Abbott+2021) używany jako observational input dla α-posterior derivation, NIE jako theoretical mechanism dla TGP. Phase model dla GW phase Δφ(f) = -(15/4)·Δe_2_native/(M·(πMf)^(1/3)) [radians] is L1 native (z LIGO-3G-native A− cycle). |
| **D — Hardcoded T_pass = True?** | ✅ NIE — protocol mandate 0 hardcoded; każdy test musi mieć symbolic verification step. |

**Sign-off:** Claudian @ 2026-05-13 Phase 2 setup. Pre-flight PASS.

## §1 — Phase 2 decision question (Q2)

> **Q2:** Który f(ψ) family member (jeśli którykolwiek) w pre-bounded recovery region
> α ∈ [-0.832, 0.832] jest *maximally compatible* z pełnym GWTC-3 90 BBH ppE posterior, i czy
> quadratic / transcendental families distinguish themselves od polynomial przy higher-PN
> orders (BH5 QNM ringdown, ε.1 photon ring) na poziomie observational?

**Phase 2 zamyka decyzję:**

- **H1a:** α_ML w recovery region z 5σ confidence + family distinguishability identified → Phase
  FINAL closes A− pending observational LIGO-O5 A+ ~2027 verification
- **H1b:** α_ML *outside* recovery region (z 5σ) OR wszystkie family pathologically degenerate
  → S07 INSUFFICIENT, accept M9.1'' framework-level falsification per PR-010 immutable

**Brak H1c/H1d** — anti-Lakatos LOCKED w PR-010.

## §2 — Sub-deliverable scope (B.1-B.4)

### §2.1 — B.1: Bayesian α posterior z GWTC-3 90 BBH (symbolic, NIE MCMC)

GWTC-3 ppE Bayesian combined posterior `p(β_ppE | GWTC-3)` (Abbott+2021) jest
LITERATURE_ANCHORED input. Mapowanie do native α via Phase 1 linear scaling:

```
β_ppE^poly(α) = (15/16)·α
↔  α = (16/15)·β_ppE
```

Posterior transformation Jacobian (linear, constant):

```
p_α(α) = p_β(β = (15/16)·α) · |dβ/dα|
       = p_β((15/16)·α) · (15/16)
```

**α posterior ML estimate:** dla typical GWTC-3 ToGR null result `β_ML ≈ 0`, dostajemy
`α_ML ≈ 0` z 1σ bound `|α_1σ| ≤ 0.832` (per Phase 1 T4 reproduced).

**Future LIGO-O5 A+ projection:** PR-002 LOCKED `SNR_M911 = 15.05σ` z M9.1'' Path 2 anchor
β_M911 = -15/4. Implied measurement precision:
`σ_β^O5 = β_M911 / SNR_M911 = (15/4) / 15.05 ≈ 0.249`. Translates to
`σ_α^O5 = (16/15)·σ_β^O5 ≈ 0.266`. **5σ exclusion bound (LIGO-O5):**
`|α|_5σ^O5 = 5·σ_α^O5 ≈ 1.33`. Wszystkie polynomial-family α in current GWTC-3 1σ recovery
region [-0.832, 0.832] are **WITHIN** O5 5σ window (NIE 5σ-excluded by O5 alone), ALE
posterior precision improves ×3.13.

### §2.2 — B.2: Higher-PN family distinguishability (BH5 QNM + ε.1 photon ring)

Phase 1 T5b/T6b argued że quadratic + transcendental degenerate z polynomial at 2.5PN
leading. Phase 2 derives **first non-degenerate order symbolic** per family via second
derivative `d²f/dψ²(ψ_0)`:

| Family | f(ψ) | df/dψ\|_{ψ_0} | d²f/dψ²\|_{ψ_0} | Phase 2 marker |
|---|---|---|---|---|
| polynomial | 1 + α(ψ-ψ_0) | α | **0** | only α |
| quadratic | 1 + α(ψ-ψ_0) + β_q(ψ-ψ_0)² | α | **2β_q** | α + new β_q |
| transcendental | exp(α(ψ-ψ_0)) | α | **α²** | α + nonlinear α² |

Drugi derivative kontroluje BH5 (QNM ringdown) leading-order modyfikację (przez K_eff(ψ_0)
local curvature) i ε.1 (photon ring) leading asymptotic via b_crit shift.

**Key observational implications:**

- BH5 QNM: `δω_QNM/ω_QNM^GR ∝ d²f/dψ²(ψ_0)` z LIGO-O5+/Cosmic Explorer ~2030+ sensitivity
- ε.1 photon ring: M9.1'' specific value `+14.6%` deviation (op-eht cycle T1) jest Path 2
  anchor data point; family-dependent extrapolation via b_crit(α, d²f/dψ²) symbolic structure

### §2.3 — B.3: Cross-cycle Δe_2_native consistency check

Z LIGO-3G-native A− cycle:

```
Δe_2_native = -4·ξ_3 + 4 - a_3/8 + c_0·κ_σ
β_ppE^TGP   = (45/16) · Δe_2_native
```

Phase 1 polynomial result β_ppE^poly = (15/16)·α implies (substytucja):

```
(15/16)·α = (45/16) · Δe_2_native(α)
→ Δe_2_native(α) = α/3
```

Z c_0·κ_σ = 4/3 LOCKED Path 2 anchor:

```
α/3 = -4·ξ_3(α) + 4 - a_3(α)/8 + 4/3
```

Dla każdego α ∈ recovery region, equation defines **1-parameter family w {ξ_3, a_3} space**
(podstawiona Δe_2 jest jednowymiarowa, dwie unknowny → 1-parametric). Phase 2 sympy verifies:

- α=0 (trivial): constraint `-4ξ_3 + 4 - a_3/8 + 4/3 = 0` → 1-param {ξ_3, a_3}
- α=-4 (M9.1''): reproduces Δe_2 = -4/3 (anchor consistency with Path 2 c_0·κ_σ=4/3)
- α=α_ML (recovery): native parameter freedom audit per PPN_AS_PROJECTION §3.3

### §2.4 — B.4: H1a/H1b verdict draft (Phase FINAL setup)

Decision matrix:

| Scenario | α_ML estimate | Family distinguishability | Verdict |
|---|---|---|---|
| α_ML ≈ 0 z GWTC-3 1σ; quadratic/trans non-degenerate at higher PN | within recovery | YES | **H1a tentative** — pending LIGO-O5 A+ verification |
| α_ML w recovery region, ale family universally degenerate | within recovery | NO | **H1a partial** — recovery exists, but Phase 2 cannot distinguish |
| α_ML excluded z GWTC-3 5σ | outside recovery | n/a | **H1b** — S07 INSUFFICIENT, accept M9.1'' framework-level falsification |

Phase 2 outputs decision draft — Phase FINAL `Phase_FINAL_close.md` (separate session) zamyka
cycle z explicit claim_status.

## §3 — Sympy substance plan (target ≥75% FP, BINDING; aim ≥80%)

**Plan: 12 FP + 3 LIT + 2 DEC** (= 15 PASS-counted; FP fraction = 12/15 = 80.0%; DEC separate).

Per AUDIT_2026-05-11 substance taxonomy:
- **FIRST_PRINCIPLES (FP):** symbolic derivation z TGP axioms (S05, ax:metric-coupling, c_0·κ_σ
  LOCK), no literature constants in PASS condition
- **LITERATURE_ANCHORED (LIT):** observational bounds (GWTC-3, LIGO-O5 A+, ngEHT) explicit-cited
- **DECLARATIVE (DEC):** structural commitments (anti-Lakatos, S05 preservation), separate count

| # | Test | Klasa | Question (substantive) |
|---|---|---|---|
| 1 | T1 | FP | Jacobian transformation `α = (16/15)·β_ppE` symbolic verification (linear, dα/dβ = 16/15) |
| 2 | T2 | FP | α posterior MAP estimate: jeśli β_ML = 0 (typical GWTC-3 ToGR null), α_ML = 0 |
| 3 | T3 | FP | α 1σ bound `|α| ≤ 0.832` derived z |β_ppE| ≤ 0.78 (Phase 1 T4 reproduced via Jacobian) |
| 4 | T4 | FP | α 1σ bound future LIGO-O5 A+ projection: σ_α^O5 = (16/15)·(β_M911/SNR_O5) symbolic |
| 5 | T5 | FP | Polynomial f(ψ) higher-PN structure: d²f/dψ²(ψ_0) = 0 → no β_quad-like contribution |
| 6 | T6 | FP | Quadratic f(ψ) higher-PN structure: d²f/dψ²(ψ_0) = 2β_quad symbolic; new family parameter |
| 7 | T7 | FP | Transcendental f(ψ) higher-PN structure: d²f/dψ²(ψ_0) = α² symbolic; nonlinear α |
| 8 | T8 | FP | BH5 QNM modification structural marker: leading-order δω_QNM ∝ d²f/dψ²(ψ_0) → distinguishes families |
| 9 | T9 | FP | ε.1 photon ring asymptotic: b_crit shift sensitive do f(ψ_photon) tail; M9.1'' Path 2 anchor +14.6% point datum |
| 10 | T10 | FP | Cross-cycle Δe_2_native(α) = α/3 derived z β_ppE^poly = (15/16)·α + β_ppE^TGP = (45/16)·Δe_2 |
| 11 | T11 | FP | Constraint -4ξ_3 + 4 - a_3/8 + 4/3 = α/3 z c_0·κ_σ=4/3 LOCK; 1-param {ξ_3, a_3} family per α |
| 12 | T12 | FP | GR-limit consistency α=0: Δe_2=0, β_ppE=0, f=1 trivial — wszystkie families collapse |
| 13 | T13 | LIT | GWTC-3 |β_ppE| ≤ 0.78 (Abbott+2021); M9.1'' β = -15/4 rejected at 5.02σ |
| 14 | T14 | LIT | LIGO-O5 A+ ~2027 SNR=15.05σ na M9.1'' β = -15/4 (PR-002 LOCKED inheritance) |
| 15 | T15 | LIT | ngEHT photon ring +14.6% M9.1'' (op-eht cycle); cross-channel BH5 QNM Cosmic Explorer ~2030 |
| 16 | T16 | DEC | Anti-Lakatos LOCKED PR-010: brak H1c/H1d backstop; recovery_scope α ∈ [-0.832, 0.832] |
| 17 | T17 | DEC | Three-layer L1/L2/L3 presentation explicit per PPN_AS_PROJECTION §3.1 (mandatory in results.md) |

**0 hardcoded `T_pass = True`. 100% non-trivial.** Każdy test ma explicit pytanie fizyczne
weryfikowane symbolic.

## §4 — Risk register Phase 2

| ID | Risk | Mitigation |
|---|---|---|
| R-P2.1 | Bayesian symbolic mapping nadmiernie uproszczone vs full MCMC | Explicit annotation: pełny PE pipeline (z noise model) jest separate cycle; Phase 2 scope = symbolic Jacobian transformation rigorous, sufficient dla H1a/H1b binary decision |
| R-P2.2 | BH5/ε.1 channel coefficients wymagają deeper SPA chain re-derivation | Phase 2 derives FIRST non-degenerate order symbolic (d²f/dψ²(ψ_0)); pełny SPA chain per family deferred do Phase 3 jeśli verdict wymaga |
| R-P2.3 | **Inheritance trap** — c_0·κ_σ = 4/3 LOCK adoptowany bez re-derivation | Explicit cite jako Path 2 anchor inheritance per M9_RESTRUCTURE §1.4 kategoria (b); dokumentowane w §0.1 ASK-RULE Trigger B; conditional citation `pending Phase 4 emergent-metric LOCK preserved` |
| R-P2.4 | **BD-drift risk** — Bayesian framing slips do "TGP α posterior" jako primary output (L2 mimicry) | L1/L2/L3 layering MANDATORY w Phase2_results.md per PPN_AS_PROJECTION §3.1; primary output (L1) = constraint on native Taylor coefs {a_n, ξ_n, c_0, κ_σ}, β_ppE posterior = L2 projection consistency check |
| R-P2.5 | **Lakatos drift** — kuszenie do dodania H1c "different anchor not M9.1''" lub H1d "BH5 channel rules out family" w trakcie cyklu | PR-010 immutable; jeśli verdict H1b → close-NULL honest per cluster-sterile-ν precedent |
| R-P2.6 | Higher-PN tests T8/T9 nadmiernie ambitne dla 1-session sympy | Tests T8/T9 verify STRUCTURAL marker (d²f/dψ²(ψ_0) per family), NIE numerical photon ring shift; numerical extrapolation deferred do dedicated `op-eht-S07-followup` cycle |

## §5 — Six P-requirements update (Phase 2 contribution)

| P | Phase 1 status | Phase 2 contribution |
|---|---|---|
| P1: f(ψ) family enumeration | ✅ RESOLVED (3 families) | extended z d²f/dψ² distinguishability marker |
| P2: β_ppE prediction per alternative | ✅ RESOLVED (β_ppE = (15/16)·α) | extended z α posterior MAP estimate |
| P3: GR limit recovery | ✅ RESOLVED (T1+T2+T5a+T6a) | reverified T12 sanity check |
| P4: GWTC-3 Bayesian compatibility | 🟡 PARTIAL (analytical range) | **CLOSED** symbolic Bayesian mapping (B.1) |
| P5: Cross-cycle consistency | ✅ RESOLVED (T7) | extended (B.3) Δe_2_native consistency under c_0·κ_σ=4/3 |
| P6: S05 preservation | ✅ RESOLVED (T11) | reverified (declarative) |

**Phase 2 zamyka P4 (GWTC-3 compatibility) symbolic Bayesian.** Pełny MCMC pipeline = separate cycle (out of substance protocol scope).

## §6 — Cross-references

- [[./README.md]] — cycle BINDING contract
- [[./Phase0_balance.md]] — initial scaffold
- [[./Phase1_results.md]] — predecessor (12/12 PASS, linear scaling discovery)
- [[../op-LIGO-3G-native-phase-residual-2026-05-11/]] — Δe_2_native methodology + PR-002
- [[../op-emergent-metric-from-interaction-2026-05-09/]] — Phase 4 zero-β region {A,B,C} + c_0·κ_σ=4/3 LOCK
- [[../op-eht/]] — ngEHT photon ring +14.6% M9.1'' anchor data point
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] — 5.02σ M9.1'' rejection
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 BINDING
- [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 — M9.1'' Path 2 anchor framing
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift Triggers A-D
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] — BINDING contract structure
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-010 — anti-Lakatos LOCKED

---

**Phase 2 setup complete.** Sympy execution + results draft pending. Estymata: 1 sesja sympy
+ results.md; Phase FINAL closure separate session.
