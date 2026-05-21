---
title: "op-L07-zero-sum-Z2-derivation — Derywacja zasady zerowej sumy ZS1/ZS2 z Z₂-symetrii substratu (L07 audit closure)"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  L1_native:
    output_observable: "Integralny warunek znikania ∫_Σ Δ(x)√h d³x [skalar, energia·objętość]; klasyfikacja statusu (aksjomat vs derywacja Z₂-tożsamość) [boolean+verdict]; konsystencja z prop:Lambda-positive (sek05) [boolean] — operationalna derywacja warunku globalnego z lokalnej Z₂-symetrii substratu (audit L07 EXT-2)"
    measurement_instrument: "Symbolic operator-level analysis z Z₂-invariant Hamiltonian H_Γ; expectation value calculus na ground state |Ω⟩; sympy symbolic verification toksamości"
    native_coefs_constrained:
      - "ZS1 ∫_Σ Δ(x)√h d³x = 0 → derived AS Z₂ identity (Path A from L07 audit)"
      - "ZS2 ∫_Σ (Φ - Φ₀)√h d³x = 0 → partial: linear part follows from ZS1 via Φ = (φ/φ_ref)²·Φ₀; quadratic remainder requires additional boundary condition"
      - "Status: ZS1 promoted from aksjomat → tożsamość; ZS2 partially-derived (gauge/boundary character)"
    falsification_rule: "Jeśli (a) symbolic derivation pokazuje że ZS1 nie wynika z Z₂-invariance H_Γ + Z₂-invariant universe state, LUB (b) relacja Φ = (φ/φ_ref)²·Φ₀ daje sprzeczność z ZS2 dla typowych konfiguracji materii, wówczas Path A z audytu L07 jest OBSTRUCTED — wymagana ścieżka B (Lagrange multiplier), C (φ_eff redefinition) lub D (nonlokalność). Werdykt B+ acceptable jeśli ZS1 derives ale ZS2 wymaga additional condition; HALT-B jeśli ani jedno ani drugie nie działa."
    pre_registration_date: "2026-05-16"

  L2_framework_reduction:
    target_frameworks:
      - "Standard QFT spontaneous symmetry breaking (Goldstone 1961, Nambu 1960)"
      - "Wess-Zumino consistency of global currents (Adler-Bardeen-Bell-Jackiw)"
      - "Gauge fixing on global zero-modes (constraint Hamiltonian, Dirac 1950)"
    reduction_type: "consistency-mapping"
    validation_transfer: "ZS1 maps to chiral order-parameter ⟨q̄γ⁵q⟩ vanishing globally w QCD analogii"
    failure_disposition: "L1-stands"

  L3_falsification_map:
    - { bound: "prop:Lambda-positive (sek05_ciemna_energia §240-293)", constrains: "ZS2 must hold for non-trivial Λ_eff > 0", window: "structural; cosmological", status: "pending Phase 1" }
    - { bound: "Substrate Z₂ symmetry of H_Γ (sek01_ontologia, eq.zero-sum-precise eq.zero-sum-Phi)", constrains: "ZS1 derivable as Z₂ identity", window: "structural; substrate-level", status: "pending Phase 1" }
    - { bound: "Φ(x) > 0 everywhere (sek02_pole §116-117)", constrains: "ZS2 cannot require pointwise Φ < Φ₀ in any region with positive substrate", window: "structural; field-positivity", status: "pending Phase 1 compatibility check" }
    - { bound: "Cosmological consistency H_0² ~ γ/12 (sek05 prop:Lambda-positive)", constrains: "average ⟨(δΦ)²⟩ sets Λ_eff scale; γ-derivation outside L07 scope", window: "cosmological; T-Λ closure preserved", status: "inherited from closure_2026-04-26" }

tgp_status:
  level: L1
  kind: derivation
  output_type: structural-identity
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "ZS2 quadratic remainder full structural origin (likely cosmological boundary condition, deferred)"
    - "Wilson-coefficient correction terms to ZS1 under higher-derivative substrate corrections (separate cycle)"
    - "Generalization to inhomogeneous Z₂-broken phases beyond mean-field (audit L07 path D path)"
  depends_on:
    - "core/sek01_ontologia/sek01_ontologia.tex ax:zero (eq.zero-sum + remark:zero-precyzacja ZS1/ZS2 split)"
    - "core/sek05_ciemna_energia/sek05_ciemna_energia.tex prop:Lambda-positive (eq.zero-sum-phi)"
    - "TGP_FOUNDATIONS §1 (single field Φ z Z₂ substrate)"
    - "research/closure_2026-04-26/Lambda_from_Phi0/ (T-Λ closure preserved)"
  impacts:
    - "audyt/L07_zero_sum_axiom — Path A closure (Z₂ identity for ZS1; partial for ZS2)"
    - "core/sek01_ontologia ax:zero status — aksjomat → derived identity (ZS1) + boundary condition (ZS2)"
    - "core/sek05_ciemna_energia prop:Lambda-positive — strengthened foundation (no longer hangs on raw aksjomat)"
    - "Cosmological constant problem disposition — Λ_eff > 0 mechanism backed by Z₂ structural argument"

predecessors:
  - "[[../closure_2026-04-26/Lambda_from_Phi0/]] (T-Λ closure; γ/12 = M_Pl²·H₀² scale fixed)"
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04/]] (ρ-T^μ_μ bridge; full SM matter coverage)"
  - "[[../op-Q2-vacuum-budget-2026-05-10/]] (substrate-vacuum decoupling from SM vacua; preserves Q2 F1 mechanism)"

related:
  - "[[../../audyt/L07_zero_sum_axiom/README.md]] (this cycle addresses)"
  - "[[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2 (review trigger)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]"
  - "[[../../STATE.md]]"

classification: DERIVATION — L07 audit Path A (Z₂-tożsamość derivation for global zero-sum)
priority: high (P2 OPEN; foundational closure for sek05 dark energy mechanism; cosmological constant problem disposition)
goal: "First-principles symbolic derivation of ZS1 (chiralna, ∫Δ√h = 0) jako Z₂-tożsamość substratu, oraz ZS2 (przestrzenna, ∫(Φ-Φ₀)√h = 0) jako częściowa konsekwencja relacji Φ = (φ/φ_ref)²·Φ₀ + Z₂. Identyfikacja statusu ZS2 quadratic remainder jako boundary condition (gauge fixing) vs additional axiom. Closure for audit L07 (P2 EXT-2)."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 symbolic + Phase FINAL closure compressed; honest partial expectation)"
target_window: "Phase 1: Z₂-tożsamość derivation ZS1 (5-6 tests); ZS2 quadratic-linear decomposition (3-4 tests); compatibility check ZS2 z Φ > 0 + Λ_eff scale (2-3 tests)"

six_requirements_target:
  - "P1: H_Γ Z₂-invariance pod φ → -φ symbolicznie zweryfikowana"
  - "P2: Z₂-invariant ground state |Ω⟩ ⇒ ⟨Δ(x)⟩_Ω = 0 pointwise (symmetric phase) lub global (spontaneous-domain phase)"
  - "P3: ZS1 ∫_Σ Δ(x)√h d³x = 0 derived AS Z₂ identity (Path A audit closure)"
  - "P4: Relacja Φ = (φ/φ_ref)²·Φ₀ expansion δΦ = 2(v/φ_ref²)·Φ₀·δφ + Φ₀/φ_ref² · (δφ)² — explicit linear + quadratic split"
  - "P5: ZS2 linear part ∫(2v/φ_ref²)·Φ₀·δφ√h = 0 — follows from ZS1 (since δφ jest Z₂-odd lokalnie around symmetric vacuum)"
  - "P6: ZS2 quadratic part ∫(Φ₀/φ_ref²)·(δφ)²√h = (Φ₀/φ_ref²)·V·⟨(δφ)²⟩ — NON-ZERO; requires boundary condition lub gauge fixing — explicit status documented"

risk_flags:
  - "R1: Spontaneous Z₂ breaking with ⟨φ⟩ = +v vs -v domains — ZS1 zachodzi tylko jeśli domains balansują się globalnie; bez tego ZS1 wymaga additional condition"
  - "R2: Φ = (φ/φ_ref)²·Φ₀ jest Z₂-EVEN — nie ma simple Z₂ identity dla ZS2; ZS2 to NIE czysta Z₂ tożsamość"
  - "R3: Compatibility Φ > 0 everywhere (sek02 § 116-117) z ZS2 wymagającym δΦ < 0 w pewnych regionach — handle przez (δΦ)_avg = 0 nie point-wise"
  - "R4: Higher-order corrections w expansion δφ/v: poniżej rzędu O((δφ)³/v³) deferred do extension cycle"
  - "R5: Cosmological boundary conditions (FRW topology) — ZS2 quadratic remainder należy do boundary condition class; explicit cosmological time-dependence outside cycle scope"

phase_plan:
  Phase_0: "Balance sheet + 6/6 gate + scope (ZS1 vs ZS2 split; honest partial expectation pre-registered)"
  Phase_1: "First-principles symbolic: Z₂-invariance H_Γ, ZS1 derivation, ZS2 linear/quadratic decomposition, status documentation"
  Phase_FINAL: "Closure + L2 reduction (SSB analog + gauge fixing) + L3 falsification map check + L07 audit closure note z honest partial verdict (B+ expected jeśli ZS2 quadratic requires additional condition)"

tags:
  - L07
  - L07-Phase1
  - zero-sum
  - dark-energy
  - axiom-vs-derived
  - Z2-symmetry
  - substrate-foundations
  - first-principles
  - audit-L07-path-A
  - cycle-scaffold-2026-05-16
---

# op-L07-zero-sum-Z2-derivation-2026-05-16

> **Cel:** First-principles symbolic derivation of zasady zerowej sumy (ax:zero z sek01_ontologia)
> jako Z₂-tożsamości substratu — closure audytu L07 (EXT-2 external review).
> ZS1 (chiralna) jest czystą Z₂-identity; ZS2 (przestrzenna) wymaga additional argument
> dla quadratic part. Cykl honestly partial — pre-registered B+ outcome możliwe.

## §0 — Cel + native-first contract

[CITE: `audyt/L07_zero_sum_axiom/README.md` §A-D; `audyt/EXTERNAL_REVIEW_2026-05-06.md` §EXT-2;
`core/sek01_ontologia/sek01_ontologia.tex` ax:zero + remark:zero-precyzacja;
`core/sek05_ciemna_energia/sek05_ciemna_energia.tex` prop:Lambda-positive]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- `ZS1_LHS` ≡ ∫_Σ Δ(x)√h d³x — chiral zero-sum LHS (Δ = local Z₂ order parameter signed)
- `ZS2_LHS` ≡ ∫_Σ (Φ(x) - Φ₀)√h d³x — spatial zero-sum LHS (Φ-fluctuation around vacuum)
- `ZS1_status` ∈ {derived-identity, axiom-stands, falsified}
- `ZS2_status` ∈ {derived-identity, partially-derived-with-boundary, axiom-stands, falsified}

**Instrument:** Z₂-invariant Hamiltonian H_Γ (sek01_ontologia eq:H-Gamma analog);
symbolic field expansion δφ = φ - v wokół symmetric vacuum lub Z₂-domain-symmetric superposition;
sympy symbolic verification operator identities.

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli (a) symbolic derivation pokazuje że ZS1 NIE wynika z Z₂-invariance H_Γ +
> Z₂-invariant universe state, LUB (b) relacja `Φ = (φ/φ_ref)²·Φ₀` daje SPRZECZNOŚĆ
> z ZS2 dla typowych konfiguracji materii, wówczas Path A z audytu L07 jest OBSTRUCTED
> — wymagana Path B (Lagrange multiplier z osobnym parametrem skali λ), Path C (φ_eff
> redefinition jako odchylenie od średniej), lub Path D (nielokalność na skali
> kosmologicznej).
>
> **Werdykt B+ pre-registered acceptable** jeśli (X1) ZS1 derives cleanly jako Z₂-tożsamość
> ALE (X2) ZS2 wymaga additional condition dla quadratic remainder (boundary/gauge nature).
> Werdykt **A−** wymaga (Y1) ZS1 + (Y2) ZS2 oboje derivable bez additional axioms.
> Werdykt **HALT-B** jeśli ani jedno ani drugie nie zachodzi.

```
pre_registration_date: 2026-05-16
recovery_scope:
  allowed_directions:
    - "Subleading corrections to Δ(x) operator definition (do not change Z₂ transformation)"
    - "Alternative spontaneous-domain parametrizations (homotopy-equivalent w Z₂ orbit)"
    - "Cosmological boundary condition declaration explicit jako separate axiom (ZS2 remainder)"
  forbidden_directions:
    - "Free parameter w ZS1 (jeśli derived from Z₂, prefactor jest exact)"
    - "Post-hoc lifting of Z₂ to U(1) (would change classification of substrate)"
    - "Adding fundamental Z₂-breaking term on top of H_Γ (S05 violation)"
  if_recovery_exhausted: "Honest verdict B+ partial; document ZS2 remainder as boundary condition status; audit L07 Path A partial closure"
```

### §0.3 — TGP-native check (mandatory)

- [x] **Q1 (Pattern coverage):** Pattern 2.5 (chiral symmetry-protected identity) relevant — ZS1 jako analog
      do ⟨q̄γ⁵q⟩=0 w QCD vacuum
- [x] **Q2 (Red flags):** WIELKIE potential red flag — globalny więz bez Lagrange'a multipliera z dynamiki;
      cykl explicit addresses this. NO post-hoc fitting allowed.
- [x] **Q3 (Inherited LOCKs):** T-Λ closure γ/12 = M_Pl²·H₀² LIVE (closure_2026-04-26);
      single-Φ axiom z Z₂ substrate LIVE (TGP_FOUNDATIONS §1);
      Φ > 0 everywhere LIVE (sek02_pole §116-117)
- [x] **Q4 (Standard-physics tools):** SSB Goldstone-Nambu framework standard; gauge-fixing
      on global zero-modes standard (Dirac constraint Hamiltonian); native-relevance:
      applied to TGP-specific substrate Z₂
- [x] **Q5 (m_Φ usage):** N/A — operator-level argument; mass scale enters only via
      consistency check z prop:Lambda-positive (γ = m²·M_Pl² approximate)
- [x] **Q6 (GR limit framing):** ZS2 zawiera √h (metric determinant of spatial slice);
      Φ > 0 ensures √h finite + positive; FRW cosmological slice handled implicitly
- [x] **Q7 (ASK-RULE self-check):** methodology cited; gaps documented w R-flags;
      honest partial expectation pre-registered
- [x] **Q8 (BD-drift audit plan):** Phase FINAL flag jeśli ZS2 remainder reintroduces
      effective Brans-Dicke-style scalar boundary condition

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2
- [x] Przeczytano [[../../audyt/L07_zero_sum_axiom/README.md]]
- [x] Przeczytano [[../../core/sek01_ontologia/sek01_ontologia.tex]] ax:zero + remark zero-precyzacja
- [x] Przeczytano [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]] §240-293 (prop:Lambda-positive)

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16

### §0.5 — Sympy substance plan

**Plan testów Phase 1 (target 10 tests, ≥8 first-principles):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | H_Γ Z₂-invariance: H_Γ[-φ] = H_Γ[φ] explicit dla typical (J·φ²·φ² + ...) Hamiltonian |
| T2 | **FIRST_PRINCIPLES** | Z₂ transformation Δ(x) → -Δ(x) explicit (Δ jest order parameter z definicji) |
| T3 | **FIRST_PRINCIPLES** | Z₂-symmetric ground state \|Ω⟩: P_Z₂\|Ω⟩ = \|Ω⟩ ⇒ ⟨Ω\|Δ(x)\|Ω⟩ = -⟨Ω\|Δ(x)\|Ω⟩ ⇒ = 0 |
| T4 | **FIRST_PRINCIPLES** | ZS1 derivation: ∫_Σ ⟨Δ(x)⟩√h d³x = 0 jako Z₂-tożsamość (Path A closure) |
| T5 | **FIRST_PRINCIPLES** | Φ(φ) = (φ/φ_ref)²·Φ₀ Z₂-EVEN: Φ(-φ) = Φ(φ) explicit (NIE bezpośrednia Z₂ tożsamość dla ZS2) |
| T6 | **FIRST_PRINCIPLES** | δΦ expansion around φ=v: δΦ = (Φ₀/v²)·(2v·δφ + (δφ)²) z δφ = φ - v |
| T7 | **FIRST_PRINCIPLES** | ZS2 linear part: ∫(2v·Φ₀/v²)·δφ√h = (2Φ₀/v)·∫δφ√h — vanishes IF δφ averages to zero (ZS1-like condition na δφ, NIE pełne ZS1 na Δ) |
| T8 | **FIRST_PRINCIPLES** | ZS2 quadratic part: ∫(Φ₀/v²)·(δφ)²√h = (Φ₀/v²)·V·⟨(δφ)²⟩ ≥ 0 — NIE może być zero dla nontrivial fluctuations |
| T9 | **FIRST_PRINCIPLES** | Compatibility check: δΦ = -δΦ_avg implementacja jako boundary condition dla globalnego ⟨δΦ⟩_Σ = 0 (gauge fixing) — NOT pointwise Φ<Φ₀ |
| T10 | **FIRST_PRINCIPLES** | Consistency z prop:Lambda-positive: Λ_eff ~ G·⟨U(δφ)⟩_Σ; ZS2 ensures ⟨δΦ⟩_Σ = 0 ALE ⟨(δΦ)²⟩_Σ > 0 — direct source for Λ_eff > 0 |
| T11 | **LITERATURE_ANCHORED** | Comparison z QCD ⟨q̄γ⁵q⟩ vanishing (chiral analog); SSB framework Goldstone-Nambu (1960-61) |
| T12 | **DECLARATIVE** | S05 single-Φ + Z₂ substrate preserved structurally; no new free parameters introduced (separate, NIE w PASS total) |

**Target:** 11/11 PASS sympy (T1-T11) + 1 structural declaration (T12 separate).

**Ratio:** 10 FIRST_PRINCIPLES (90.9%) + 1 LITERATURE_ANCHORED (9.1%) + 1 DECLARATIVE separate.

---

## §1 — Phase 0: balance sheet

[Patrz `Phase0_balance.md`]

## §2 — Phase 1: native derivation

[Patrz `Phase1_sympy.py` + `Phase1_results.md`]

## §FINAL — Closure

[Patrz `Phase_FINAL_close.md`]

---

## Status

🟢 **ACTIVE — opened 2026-05-16** per user authorization "wybierz kolejny task z research i rozpocznij pracę"
(L07 selected from STATE.md candidate list).

Cycle scope: focused ZS1+ZS2 derivation z Z₂-symetrii substratu (audit L07 Path A, P2);
full nonlocal foundations (L07 Path D) deferred to downstream multi-session cycle.

This session deliverables:
- README.md (this file) z BINDING contract — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (11-test first-principles symbolic)
- Phase1_sympy.txt — **PLANNED** (sympy output)
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED if Phase 1 PASS or HALT-acceptable**

Pre-registered partial outcome: **B+ expected** (ZS1 derived clean; ZS2 quadratic-remainder
status as boundary condition explicit) — analog to L08-e²-derivation cycle (2026-05-16) honest partial.

---

**Cross-references:**
- [[../../audyt/L07_zero_sum_axiom/README.md]] (cycle addresses this audit issue)
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-2 (review trigger)
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] (ax:zero source; remark:zero-precyzacja ZS1/ZS2)
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]] (prop:Lambda-positive depends on ZS2)
- [[../closure_2026-04-26/Lambda_from_Phi0]] (T-Λ closure; γ/12 = H₀² scale)
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] (ρ-T^μ_μ bridge; full SM matter coverage)
- [[../op-Q2-vacuum-budget-2026-05-10/]] (substrate-vacuum decoupling)
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]]
- [[../../STATE.md]]
