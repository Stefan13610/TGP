---
title: "op-L08-Phase6-Dirac-propagator — Emergent Dirac propagator z FR antisym + Clifford + L05 m_obs"
date: 2026-05-16
pre_registration_date: 2026-05-16
parent: "[[../../audyt/L08_kink_fermion_closure/README.md]]"
cycle: L08 Phase 6 (warstwa 3c kink-fermion closure)
status: CLOSED-RESOLVED — A− claim_status, 13/13 sympy PASS, single-session sprint
folder_status: closed-resolved
claim_status: A-
closure_date: 2026-05-16
may_edit_core: false
pre_registration: "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-### TBD (Phase 0 entry)"
audit_source: "[[../../audyt/L08_kink_fermion_closure/README.md]] §Problem #1 closure aftermath + audit Ścieżka A"
priority: P2 (klaster D ontology)
related:
  - "[[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] (Problem #1 CLOSED A−; FR antisymmetry inheritance)"
  - "[[../op-L08-Phase6-Clifford-emergence-2026-05-16/]] (Problem #4 CLOSED A−; Cl(1,3) algebra inheritance)"
  - "[[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (m_obs ≠ M_full distinction LIVE)"
  - "[[../why_n3/]] Phase 1-5 closed 2026-05-01"
  - "[[../../TGP_FOUNDATIONS.md]] §4 warstwa 3c"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] BINDING contract spec"
  - "[[../../meta/CYCLE_LIFECYCLE.md]] claim_status taxonomy"
  - "[[../../meta/PPN_AS_PROJECTION.md]] §3.1 three-layer L1/L2/L3 mandatory"
tgp_status:
  level: L1
  kind: derivation
  output_type: observable  # emergent Dirac propagator z fermion kinematic predictions
  core_compatibility: extension
  may_edit_core: false
  exports_findings: true
  has_needs_file: false  # to be created Phase 0 if multi-need scope
  has_findings_file: false  # to be created Phase FINAL
  depends_on:
    - "L08-FR antisymmetry (FR 2-particle exchange = Z₂ generator)"
    - "L08-Clifford Cl(1,3) emergence (M9.1'' Lorentz signature)"
    - "L05 m_obs ≠ M_full distinction (k_obs vs k_full)"
    - "why_n3 Phase 1-5 (kink-as-fermion strukturalny szkic)"
    - "S05 single-Φ axiom preserved"
  impacts:
    - "TGP_FOUNDATIONS §4 warstwa 3c (H → partial-(D) upgrade path)"
    - "L08 audit problem disposition (NEW: emergent Dirac propagator operational closure)"
    - "Future cycles: quark/neutrino/boson sector (problem #3) inheritance"
tags:
  - TGP
  - L08
  - Phase6
  - emergent-Dirac
  - propagator
  - kink-fermion
  - warstwa-3c
  - FR-antisymmetry-inheritance
  - Clifford-inheritance
  - m_obs-renormalized-pole
  - pre-registration-locked
---

# op-L08-Phase6-Dirac-propagator — BINDING contract

> **Cykl scaffold:** 2026-05-16
> **Status:** PARKING (Phase 0 scaffold ready; Phase 1 wymaga user authorization)
> **claim_status:** TBD pending Phase 1+ execution
> **Pre-registration date:** 2026-05-16 (immutable, BINDING)

## §0 — BINDING contract (pre-Phase-1 LOCKED)

### §0.1 — Cycle goal (binding)

**Wyprowadzić explicit emergent Dirac propagator** w TGP-substrate framework,
łącząc trzy fundament zamknięte:

1. **L08-FR antisymmetry** (problem #1 CLOSED A− 2026-05-16): 2-particle Fock space
   antisymmetric pod exchange (Z₂ × Z₂ × Z₂ z C_2-defect; γ_exchange = π Berry phase)
2. **L08-Clifford algebra** (problem #4 CLOSED A− 2026-05-16): Cl(1,3) emergence
   z M9.1'' Lorentz signature (NIE z Z₂ samodzielnie — pełna geometria substrate)
3. **L05 m_obs ≠ M_full** (CLOSED-RESOLVED A− 2026-05-16): k_obs(α=2, d=3)=3 LIVE
   jako renormalized pole mass (analog ADM-vs-Komara w GR)

**Target form (predykcja LOCKED pre-Phase 1):**

```
S_F^TGP(p) = i / (γ^μ p_μ − m_obs + iε)
         = i (γ^μ p_μ + m_obs) / (p² − m_obs² + iε)

z m_obs = c_M · A_tail² · g_0^(e²(1−α/4))|_{α=2}  (L05 + why_n3 Phase 5 inheritance)
       ≈ k_obs(α=2, d=3)-fit  (LP-4 + R3 reconciliation)
```

**3-warstwowa specyfikacja (per PPN_AS_PROJECTION §3.1 BINDING):**

- **L1 (native predictions):** kink-pair propagator z explicit antisymmetric Fock state
  → emergent S_F(p) z pole at m_obs; consistency z Hermitian conjugation +
  CPT z substrate-level Z₂ + Lorentz emergent
- **L2 (standard QFT projection):** S_F^TGP must reduce do canonical Dirac propagator
  w limicie Φ → Φ_0 (today's vacuum) + free-field limit (no Φ-mediated interaction)
- **L3 (falsifikator):** PDG fermion mass spectrum (m_μ/m_e, m_τ/m_e); QED elektrowy
  g-2 anomaly; Lamb shift z modyfikacji propagatora w obecności external V_eff[Φ]

### §0.2 — Pre-registered falsification rule (LOCKED, immutable)

**Falsification statement (BINDING, no post-hoc weakening):**

> *Jeśli Phase 1 NIE potwierdzi struktury* `S_F^TGP(p) = i(γ·p + m_obs)/(p² - m_obs² + iε)`
> *w limicie Φ → Φ_0 + free-field, z m_obs identifying jako pole mass od L05 k_obs(α=2, d=3)=3
> via universal mass formula why_n3 Phase 5, TGP emergent Dirac claim jest:*
> - *Strukturalnie nieusunięty* (Phase 1 ujawnia brak natural Dirac structure)
> - *Wymaga substrate-level pivot* (Path D — extension Z₂ → SU(2)/Spin(3); naruszenie S05/§1 FOUNDATIONS)
> - *Lub demoted do (H) hipoteza permanently* (Path C audit Ścieżka C)

**Recovery scope (pre-bounded, anti-Lakatos LOCKED):**

```
allowed (Phase 1+):
  - Refinement Lorentz boost generator construction z M9.1'' (per Cl emergence cycle)
  - Refinement of m_obs pole identification (L05 k_obs LIVE; do precision <1% from PDG)
  - Wilson-coefficient style anomalous magnetic moment correction (Phase 2 deferred)

forbidden:
  - Substrate Z₂ → SU(2)/Spin(3) extension (S05 violation; FOUNDATIONS §1)
  - Post-hoc Φ_0 variation tuning to match observation
  - L01-N1-style trace anomaly retrofit (separate cycle; orthogonal coupling)

if_recovery_exhausted:
  H1c: TGP emergent Dirac propagator insufficient z FR + Cl + L05 alone
  → multi-session extension (γ^5 chirality + axion-like coupling z L06 m_X FREE)
  → OR Path C audit (warstwa 3c permanent hypothesis)
```

### §0.3 — Six P-requirements (declarative pre-Phase-1)

Per [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §2.3:

| # | P-requirement | Resolution mechanism |
|---|---|---|
| **P1** | Native S_F^TGP(p) explicit derivation | Phase 1 sympy T1-T4: 2-particle Fock state + antisym + propagator construction |
| **P2** | Pole mass identification z L05 m_obs | Phase 1 T5-T7: pole at p² = m_obs² verification; m_obs from L05 k_obs(α=2,d=3)=3 LIVE |
| **P3** | Lorentz covariance z Cl(1,3) | Phase 1 T8-T9: γ^μ algebra inheritance from L08-Clifford; boost generators consistent |
| **P4** | Free-field limit recovery (Φ → Φ_0) | Phase 1 T10: limit Φ → Φ_0 reduces S_F^TGP do canonical Dirac S_F^Dirac |
| **P5** | CPT consistency | Phase 1 T11: substrate-level Z₂ + Lorentz emergent → CPT theorem operational |
| **P6** | S05 single-Φ preservation | Phase 1 T13 DECLARATIVE: no new fundamental fields; emergent fermion via kink-Φ |

### §0.4 — Inheritance ledger (Phase 0)

| Source | Element | Status (pre-Phase-1) |
|---|---|---|
| L08-FR | 2-particle antisymmetric Fock space | ✅ LIVE A− (12/12 sympy PASS) |
| L08-FR | γ_exchange = π Berry phase (Z₂ generator) | ✅ LIVE |
| L08-FR | χ_exchange = exp(iπ) = -1 | ✅ LIVE |
| L08-Clifford | Cl(1,3) algebra z M9.1'' Lorentz signature | ✅ LIVE A− |
| L08-Clifford | {γ^μ, γ^ν} = 2g^μν emergence (substrate-derived) | ✅ LIVE |
| L05 | m_obs ≠ M_full distinction (operational) | ✅ LIVE A− |
| L05 | k_obs(α=2, d=3) = 3 EXACT (Sobolev p_crit−α) | ✅ LIVE |
| L05 | Reconciliation theorem LP-4 ⊕ R3 | ✅ LIVE |
| why_n3 | Phase 5 universal mass formula `m_obs = c_M·A²·g_0^(e²(1-α/4))` | ✅ LIVE (numerical anchor) |
| why_n3 | Phase 3 RP² + Berry phase π (spin-1/2) | ✅ LIVE (2026-05-01 closure) |
| L01 | ρ ≡ -T^μ_μ/c_0² formal definition | ✅ LIVE (L01 cykl 2026-05-04) |
| S05 | Single-Φ axiom (no new fundamental fields) | ✅ LOCKED |

### §0.5 — Output type (binding pre-Phase-1)

**`output_type: observable`** — cycle target IS native observable (Dirac propagator
structure z m_obs pole identification → PDG fermion mass spectrum projection).
Per claim_status taxonomy:
- A+ requires L2 reduction analytical FP-grade + validation transfer
- A requires L2 attempted + numerical agreement
- A− requires native L1 derivation + PR entry (no L2 yet)

**Initial target:** **A−** (single-session Phase 0 + Phase 1); upgrade do **A** possible
if Phase 2 includes explicit PDG ratio comparison; **A+** wymaga full L2 QFT mapping
(deferred multi-session).

### §0.6 — Three-layer L1/L2/L3 specification (BINDING)

#### L1 — Native predictions (PRIMARY)

```
1. Two-particle kink Fock state |k_1, k_2⟩_anti z explicit antisymmetry pod exchange:
   |k_1, k_2⟩_anti = (1/√2)(|k_1⟩ ⊗ |k_2⟩ - |k_2⟩ ⊗ |k_1⟩)
   χ_exchange = -1 inheritance from L08-FR T7

2. Emergent γ^μ operators z transport Berry'ego po C_2-defect:
   {γ^μ, γ^ν} = 2η^μν   (η = M9.1'' Lorentz signature)
   L08-Clifford inheritance T6

3. Dirac equation operator (i γ^μ ∂_μ - m_obs) Ψ = 0:
   m_obs = c_M · A_tail² · g_0^(e²(1-α/4))|_{α=2}
         z L05 + why_n3 Phase 5 universal mass formula

4. Propagator:
   S_F^TGP(p) = i(γ·p + m_obs)/(p² - m_obs² + iε)
   Pole at p² = m_obs² (physical mass)
   Hermicity: S_F^†(p) = γ_0 S_F(p) γ_0 (CPT operational)

5. Free-field limit (Φ → Φ_0, no Φ-mediated interaction):
   S_F^TGP|_{Φ→Φ_0} = S_F^Dirac canonical (Bjorken-Drell convention)
```

#### L2 — Standard QFT projection (consistency map)

```
1. S_F^TGP must reduce do canonical Dirac propagator w Φ → Φ_0 limit:
   lim_{Φ→Φ_0} S_F^TGP(p) = S_F^Dirac(p) = (i)(γ·p + m)/(p² - m² + iε)
   
2. Wick rotation consistency: Euclidean propagator G_E(p) = 1/(γ·p_E + m)
   z proper iε prescription dla Feynman rules

3. Anomalous magnetic moment as Wilson coefficient (Phase 2 deferred):
   a_e^TGP = (α/(2π)) · (1 + Wilson_TGP[Φ_0])
   compared with PDG a_e = 0.00115965218073(28)

4. Standard QFT Slavnov-Taylor identities preserved (BRST cohomology bypass —
   TGP NIE ma gauge symmetry breaking issue z emergent Lorentz)
```

#### L3 — Falsification map

```
1. PDG fermion mass spectrum:
   - m_μ/m_e = 206.7682 (PDG): k_obs(α=2,d=3)=3 LIVE z L05 + why_n3 Phase 5
   - m_τ/m_e = 3477.23 (PDG): inherited from L05 reconciliation
   - m_e (absolute) NOT predicted z tego cyklu (substrate scale parametrowanie)
   
2. QED anomalous magnetic moment a_e (precision 0.28 ppb):
   PDG a_e = 0.00115965218073(28)
   TGP correction projected ~(α/(2π))·(Φ_TGP/Φ_0)^k correction — Phase 2 deferred

3. Lamb shift (precision MHz):
   Δν_Lamb,TGP = Δν_Lamb,QED + δ_TGP[m_obs]
   Phase 2 estimate; not directly Phase 1 scope

4. CPT theorem (EOTV: |g_e+ − g_e-|/g < 2·10⁻¹⁰; m_p+−m_p- < 10⁻⁸ z PDG):
   Substrate Z₂ + emergent Lorentz → CPT operational; Phase 1 T11 DECLARATIVE check
```

### §0.7 — What this cycle does NOT do (scope-out)

- **NIE wprowadza** Wilson coefficient corrections z curvature × propagator
  (deferred do Phase 2 extension OR L08-Phase6-precision cycle)
- **NIE rozstrzyga** problem #2 (e² in exponent) — that's L08-e² cycle scope (PARTIAL B+)
- **NIE rozstrzyga** problem #3 (quarks/neutrinos/bosons) — multi-session
- **NIE modyfikuje** core/*.tex files (may_edit_core: false)
- **NIE dotyka** L01-N1 quantum trace anomaly z EM coupling (orthogonal cycle)
- **NIE pisze pełnego** Boltzmann/loop-level QED phenomenology (substance protocol scope)

### §0.8 — Phase plan (program)

| Phase | Goal | Estimated effort | Sympy tests |
|---|---|---|---|
| **Phase 0** | Balance sheet (8/8 ☑) + Phase 1 setup | ~30 min | N/A (scaffold) |
| **Phase 1** | T1-T13 first-principles + declarative S05 | ~2-3 sesje | 13/13 PASS target |
| **Phase 2** (optional) | Wilson coefs + PDG precision comparison | ~1-2 sesje | 5-8 tests |
| **Phase FINAL** | Closure ceremony + cross-cycle propagation | ~30 min | N/A |

**Single-session ambition:** Phase 0 + Phase 1 + Phase FINAL = A− closure dziś
(jak L05 + L08-FR + L08-Clifford precedensy 2026-05-16).

### §0.9 — Pre-registration commitment

**Cycle pre_registration_date:** 2026-05-16 (LOCKED, immutable per PRE_REGISTERED_FALSIFIERS protocol)

**Anti-Lakatos LOCK:** falsification rule §0.2 jest BINDING. Phase 1 outcome A−/B+/HALT-B
musi być honestly reported per cycle-substance metrics; brak post-hoc weakening.

**Audit invariant:** This README jest BINDING contract. Any modifications post-Phase-1 są
restricted do §0.4 (inheritance ledger update z LIVE → CONSUMED) i §0.5 (claim_status
update). Sections §0.1-§0.3, §0.6-§0.8 są IMMUTABLE post-Phase-1.

## §1 — Activation pre-conditions

Cycle aktywuje się przez:
1. ✅ Phase 0 balance sheet 8/8 ☑ PASS (per [[./Phase0_balance.md]])
2. ⏸ Authorization "L08-Dirac Phase 1 start" lub similar explicit user trigger
3. ⏸ WIP slot available (check STATE.md current occupation)

## §2 — Cross-references

### §2.1 — Audit predecessors

- [[../../audyt/L08_kink_fermion_closure/README.md]] §Ścieżka A (cycle proposed)
- [[../../audyt/L05_mass_exponent_drift/README.md]] (m_obs distinction)
- [[../../audyt/PRIORITY_MATRIX.md]] EXT-4 / L08 P2 entry

### §2.2 — Cycle predecessors (closed, LIVE inheritance)

- [[../op-L08-Phase6-FR-antisymmetry-2026-05-16/]] (Problem #1 A−)
- [[../op-L08-Phase6-Clifford-emergence-2026-05-16/]] (Problem #4 A−)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (m_obs distinction A−)
- [[../why_n3/]] Phase 1-5 (2026-05-01 closure)

### §2.3 — Methodological references

- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] (BINDING contract template)
- [[../../meta/CYCLE_LIFECYCLE.md]] (claim_status taxonomy)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 (three-layer L1/L2/L3)
- [[../../meta/CALIBRATION_PROTOCOL.md]] (Phase 6 ABSOLUTE BINDING)
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] (PR entry placeholder)

### §2.4 — Downstream impact (post-closure)

- TGP_FOUNDATIONS §4 warstwa 3c (H → partial-(D) upgrade if A− achieved)
- audyt/L08 problem #1+#4 inherited → emergent Dirac OPERATIONAL
- Future cycles: quark sector (problem #3) inheritance available

## Status

🟢 **CLOSED-RESOLVED** — single-session sprint executed 2026-05-16 (user "aktywuj
L08-Dirac Phase 1").

**Phase 1 verdict:** **13/13 sympy PASS** — 10 FP (76.9%) + 1 LIT + 2 DEC; 0 hardcoded.
**6/6 P-requirements RESOLVED.** **5/6 R-flags closed** (R6 Wilson coefs deferred Phase 2).

**Phase FINAL closure:** claim_status **A−** (STRUCTURAL_DERIVED_NATIVE_PARTIAL).

**Strukturalny milestone:** Emergent Dirac propagator
`S_F^TGP(p) = i(γ·p + m_obs)/(p² − m_obs² + iε)` **OPERACYJNIE SKONSTRUOWANY** z
triple inheritance FR antisym + Cl(1,3) algebra + L05 m_obs distinction. **L08 audit
warstwa 3c (H) → partial-(D) STRENGTHENED**.

**Full closure documentation:** [[./Phase_FINAL_close.md]].
