---
# ============================================================================
# op-S07-Phase-3-BH5-eps1-numerical-2026-05-14 — pre-observational family discrimination
# Successor cycle do op-S07-reset-alternative-f-psi-2026-05-11 (CLOSED-RESOLVED A−)
# Per Phase_FINAL_close §6 upgrade path A− → A:
#   "Complete Phase 3 BH5/ε.1 numerical evaluation per family"
# Pre-observational: NIE wymaga LIGO-O5 data; computes WHAT each family PREDICTS
#   given existing observational sensitivity envelopes.
# ============================================================================

title: "S07 Phase 3 BH5/ε.1 numerical — pre-observational f(ψ) family discrimination via QNM ringdown + photon ring observables"
date: 2026-05-14
type: research-cycle
folder_status: closed-resolved   # 2026-05-14 sesja P3-FINAL: active → closed-resolved per Phase FINAL closure ceremony Opcja A heroic; claim_status A−; WIP slot 1/5 FREED
claim_status: A-                 # STRUCTURAL_DERIVED_NATIVE z L2 not-fully-FP-attempted (per S07-reset/inflation A− template)
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "δω_QNM/ω_GR (BH5 ringdown frequency relative shift, dimensionless) ORAZ δε_ph²/ε_ph²_GR (ε.1 photon ring quadrant relative shift, dimensionless) — obie obserwable evaluated per S07 alternative f(ψ) family {polynomial, quadratic, transcendental} z d²f/dψ²(ψ_0) family marker {0, 2β_q, α²} (inherited z S07-reset Phase 2 T6+T7+T15)"
    measurement_instrument: "LIGO O5 A+ ~2027 single-event ringdown σ_f ~0.5%; Cosmic Explorer ~2030 stacked ringdown 5σ family discrimination; ngEHT ~2030 photon ring 5% angular resolution; LISA ~2035 EMRI ringdown ~0.1% precision (BH5 LIVE entry Sector 3 + E2-E6 LIVE entries photon-ring sector inheritance)"
    native_coefs_constrained:
      - "α (S07 polynomial coefficient, inherited z S07-reset Phase 2 PR-010 recovery region α ∈ [-0.832, 0.832])"
      - "β_q (S07 quadratic coefficient, NEW exposure dla family marker channel)"
      - "α (S07 transcendental coefficient, NEW exposure dla family marker channel; SAME α-symbol jak polynomial; quadrant differs)"
      - "d²f/dψ²(ψ_0) family marker mapowany na BH5 + ε.1 observable channels"
    falsification_rule: "Jeśli LIGO-O5 A+ ~2027 + Cosmic Explorer ~2030 ringdown stack 100+ BH events daje δω_QNM/ω_GR poza family-predicted range {polynomial=0, quadratic=2β_q·κ_QNM, transcendental=α²·κ_QNM} z 5σ confidence, OR ngEHT ~2030 photon ring stack 10+ SMBH measurements daje δε_ph²/ε_ph²_GR poza family-predicted range {polynomial=0, quadratic=β_q/9, transcendental=α²/18} z 5σ, S07 alternative f(ψ) family pre-observational discrimination FAILS → H1b verdict: framework requires deeper structural specification beyond f(ψ) freedom (e.g., explicit substrate-dynamics V(Φ) BH-background derivation)"
    pre_registration_date: "2026-05-14"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks:
      - "Berti-Cardoso QNM perturbation spectroscopy (Schwarzschild background → modified background z g_eff[Φ_eq(BH)] family-specific Taylor expansion)"
      - "Cunha-Herdeiro photon ring shift formulas (Schwarzschild photon orbit → modified geodesic z g_eff[Φ_eq(BH)])"
      - "ppE chart projection (cross-check z S07-reset β_ppE^poly(α)=(15/16)·α inheritance)"
    reduction_type: "analytical-approximate"
    validation_transfer: "BH5 QNM observable (LIGO-O5 σ_f ~0.5%) inherits constraint na family marker d²f/dψ²(ψ_0); ε.1 photon ring (ngEHT σ_θ ~5 μas) inherits constraint na same marker via different quadrant"
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - bound: "BH5 LIVE entry: δf/f ≈ 8-16% in QNM ringdown (PREDICTIONS_REGISTRY Sector 3, op-bh-alpha-threshold/Phase3 T3.2)"
      constrains: "α=-4 (M9.1'' anchor) → BH5 prediction in family-trans channel α²/quadrant"
      window: "8-16% shift; LIGO O5 sensitivity ~0.5%/event (~2027+)"
      status: "anchor-consistency-check"
    - bound: "ε.1 LIVE entries E2-E6: ε_ph² family-quadrant marker (PREDICTIONS_REGISTRY photon-ring sector, op-eps-photon-ring/Phase3 E3.x)"
      constrains: "α=-4 (M9.1'' anchor) → ε.1 prediction +14.6% (existing data point)"
      window: "+14.6% photon ring shift z ngEHT σ_θ ~5 μas (~2030+)"
      status: "anchor-consistency-check"
    - bound: "S07-reset PR-010 LOCKED-PENDING-DATA recovery region α ∈ [-0.832, 0.832]"
      constrains: "α (polynomial), α (transcendental) — coupled WITHIN this region"
      window: "1σ GWTC-3 + 5σ LIGO-O5 A+ projection σ_α^O5 = 80/301 ≈ 0.266"
      status: "inherited-LIVE"
    - bound: "Emergent-metric Phase 4 c_0·κ_σ = 4/3 EXACT LOCK (Path 2 anchor; M9_RESTRUCTURE Tier 2)"
      constrains: "Δe_2_native(α) = α/3 substitution chain inheritance"
      window: "structurally exact, no observational tuning"
      status: "inherited-LOCKED"

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T2                              # Tier 2 anchor — extends S07 family marker to BH5/ε.1 observable channels per M9_RESTRUCTURE §2
  kind: derivation                       # symbolic + numerical projection
  output_type: observable                # δω_QNM/ω_GR + δε_ph²/ε_ph²_GR are dimensionless physical observables
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges: []
  depends_on:
    - "[[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] (predecessor A−; family marker source)"
    - "[[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] (c_0·κ_σ=4/3 LOCK; {A,B,C} family)"
    - "[[../op-LIGO-3G-native-phase-residual-2026-05-11/Phase6_close.md]] (Δφ(f) methodology; PR-002 LIGO-O5 sensitivity)"
    - "[[../op-bh-alpha-threshold/Phase3_results.md]] (BH5 LIVE existing prediction; T3.2)"
    - "[[../op-eps-photon-ring/Phase3_results.md]] (ε.1 LIVE existing entries; E3.x)"
    - "[[../op-eht/]] (M9.1'' +14.6% photon ring data point)"
  impacts:
    - "[[../../PREDICTIONS_REGISTRY.md]] — proposed S07-Recovery-α-Polynomial-Family entry extension (BH5/ε.1 family-channel projections per family)"
    - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-012 entry (this cycle)"
  source_of_status: []

predecessors:
  - "[[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] — A− closure 2026-05-13 sesja P-FINAL; family marker source"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — Path 2 anchor c_0·κ_σ=4/3 LOCK"
  - "[[../op-bh-alpha-threshold/Phase3_results.md]] — BH5 LIVE δf/f 8-16% prediction (T3.2)"
  - "[[../op-eps-photon-ring/Phase3_results.md]] — ε.1 LIVE entries E2-E6"

classification: DERIVATION  # symbolic family marker → observable mapping + numerical projection
goal: "Pre-observational discrimination S07 alternative f(ψ) families {polynomial, quadratic, transcendental} via two complementary observable channels: BH5 QNM ringdown frequency shift + ε.1 photon ring quadrant shift. Symbolic mapping family-marker d²f/dψ²(ψ_0)={0, 2β_q, α²} → observable predictions; numerical projections vs LIGO-O5/CE/ngEHT/LISA sensitivities; cross-cycle M9.1'' anchor consistency at α=-4 (BH5 +8-16% + ε.1 +14.6%)."
target_window: "BH5: LIGO O5 A+ ~2027 single-event σ_f ~0.5%; Cosmic Explorer ~2030 stacked 5σ; LISA ~2035 EMRI ~0.1%. ε.1: ngEHT ~2030 σ_θ ~5 μas / 10-SMBH stack 5σ; LISA EMRI photon ring ~2035 ~0.5%."

six_requirements_target:
  - "P1: BH5 QNM symbolic mapping δω_QNM/ω_GR = κ_QNM·d²f/dψ²(ψ_0) per family — derived analytically z form-meaning analog Berti-Cardoso QNM perturbation spectroscopy (Phase 1)"
  - "P2: ε.1 photon ring symbolic mapping δε_ph²/ε_ph²_GR = κ_ε·d²f/dψ²(ψ_0) per family — derived analytically z form-meaning analog Cunha-Herdeiro photon ring formulas (Phase 2)"
  - "P3: Cross-cycle M9.1'' anchor consistency: BH5 prediction at α=-4 (transcendental channel α²=16) consistent z LIVE 8-16% range; ε.1 prediction at α=-4 consistent z +14.6% data point (Phase 1+2 anchor checks)"
  - "P4: Numerical projections per family at fiducial α∈{-0.832, 0, 0.832}, β_q∈{-0.4, 0, 0.4} (1σ GWTC-3 boundary fiducials inherited z S07-reset Phase 2): LIGO-O5 σ_BH5, ngEHT σ_ε.1, family discriminability matrix (Phase 3)"
  - "P5: Form-meaning split per Pattern 2.2 explicit annotation: standard QNM/photon-ring formulas (BD-form) vs TGP-substrate g_eff[Φ_eq(BH)] family-modified meaning (TGP-meaning); ASK-RULE Trigger C resolved via §0.1 form-meaning section + Phase 1+2 cite per test"
  - "P6: S05 single-Φ axiom preserved bezwarunkowo across BH5 + ε.1 + 3 families (no Φ_2 / hidden field / multi-substrate; Pattern 2.5 environment-dependent observable preserved)"

risk_flags:
  - "R1 (HIGH RISK Trigger C BD-drift): Standard QNM/photon-ring formulas (Berti-Cardoso 2009, Cunha-Herdeiro 2018) są NOT TGP-native. Bez explicit form-meaning split, derivation reproduces standard physics bez TGP mechanism. MITIGATION: §0.1 explicit form-meaning split per Pattern 2.2 + Phase 1+2 sympy każdy test cite źródła; Phase FINAL spawn bd-drift-audit subagent per CALIBRATION_PROTOCOL §4.4"
  - "R2 (overlap z istniejącymi cyklami): op-bh-alpha-threshold/Phase3 (BH5 LIVE) + op-eps-photon-ring/Phase3 (ε.1 LIVE) już mają predictions specifically dla M9.1'' anchor. Ten cycle EXTENDS predictions across S07 family parameter space (NEW deliverable: per-family observable channel mapping). MITIGATION: §0.5 sympy substance plan ma explicit delta-only contribution; Phase 1+2 anchor checks demonstrate cross-cycle coherence przy α=-4"
  - "R3 (numerical fiducials post-PE): Phase 3 numerical projections require fiducial α/β_q values, NIE first-principles. Wartości {α=-0.832 (1σ boundary), α=0 (ML), α=+0.832 (1σ boundary)} są post-fit fiducials z S07-reset Phase 2 GWTC-3 Bayesian. MITIGATION: explicit annotation 'numerical projections vs sensitivity envelopes; NIE prediction-grade values'"
  - "R4 (substance ceiling): Pre-observational cycle = symbolic + literature-anchored work. NO new observational data → zero NEW first-principles content beyond family marker mapping. MITIGATION: realistic claim_status A− ceiling per S07-reset/inflation precedent; full A would require BH5/ε.1 actual detection data"
  - "R5 (S05 multi-channel preservation): Two observable channels (BH5 + ε.1) muszą oba inherit single-Φ substrate; brak hidden field per channel. MITIGATION: Phase 1+2 explicit DEC tests S05 preservation per channel + cross-channel"

phase_plan:
  Phase_0: "Balance sheet + pre-flight methodology read confirmation + form-meaning split annotation + sympy substance plan"
  Phase_1: "BH5 QNM symbolic family marker mapping (~12 sympy tests; first-principles Berti-Cardoso form + TGP-meaning per Pattern 2.2; family enumeration {poly=0, quad=2β_q·κ_QNM, trans=α²·κ_QNM}; M9.1'' anchor α=-4 → 8-16% consistency; LIGO-O5/CE PSD literature-anchored)"
  Phase_2: "ε.1 photon ring symbolic family marker mapping (~12 sympy tests; first-principles Cunha-Herdeiro form + TGP-meaning per Pattern 2.2; family enumeration {poly=0, quad=β_q/9, trans=α²/18}; M9.1'' anchor α=-4 → +14.6% data point consistency; ngEHT σ_θ literature-anchored)"
  Phase_3: "Numerical projections per family at fiducials + family discriminability matrix (~10 sympy tests; LIGO-O5/CE σ_BH5; ngEHT σ_ε.1; cross-channel coupled bound calculus; LISA 2035+ EMRI projection)"
  Phase_FINAL: "Closure ceremony A− (per S07-reset/LIGO-3G-native template) + L2 ppE projection consistency check + L3 falsification map cross-cycle update + bd-drift-audit subagent + STATE.md update + PR-012 LOCKED-PENDING-DATA"

tags:
  - S07-recovery
  - alternative-f-psi
  - BH5-QNM-ringdown
  - eps1-photon-ring
  - pre-observational-discrimination
  - family-marker-mapping
  - form-meaning-split
  - cycle-scaffold-2026-05-14
---

# S07 Phase 3 BH5/ε.1 numerical — pre-observational f(ψ) family discrimination

> **Cel:** Map S07 alternative f(ψ) family marker `d²f/dψ²(ψ_0) = {0, 2β_q, α²}` (inherited z S07-reset Phase 2 PR-010 closure A−) na DWA observable channels (BH5 QNM ringdown shift + ε.1 photon ring quadrant shift) per family, providing pre-observational family discrimination test pre-LIGO-O5 ~2027 + ngEHT ~2030.

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `meta/PPN_AS_PROJECTION.md` §3.1; `meta/M9_RESTRUCTURE_NOTE.md` §2; `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1+§2.2+§4]

### §0.1 — Form-meaning split annotation (ASK-RULE Trigger C resolution per Pattern 2.2)

**Standard physics default (anti-pattern §3 red flag #1):**

W standard scalar-tensor gravity (Berti-Cardoso 2009, Cunha-Herdeiro 2018), modified BH spectroscopy + photon ring shift są derived from:
- Modified BH metric (e.g., Schwarzschild + scalar hair from BD/Horndeski)
- Linearized Klein-Gordon perturbation around scalar hair
- QNM eigenfrequency shifts ∝ scalar coupling strength
- Photon ring shift ∝ scalar background gradient at photon orbit r_ph

**To jest BD-mode framing.** Matter (test photon, gravitational perturbation) couples bezpośrednio do scalar field przez vertex. QNM/photon-ring shifts są tree-level scalar exchange manifestations.

**TGP-native podejście (per Pattern 2.2 + Pattern 2.1):**

W TGP, BH5 + ε.1 modifications nie pochodzą z direct Φ-coupling do test particles. Pochodzą z:

1. **Φ_eq(BH-background) = pełne nieliniowe rozwiązanie Φ-EOM** (Pattern 2.1) z BH source ρ_BH:
   ```
   D_kin[Φ_eq] + V'(Φ_eq) = -q · ρ_BH       [foundations §3 eq 134]
   ```
   gdzie `D_kin = (1/3φ²)·∇²(φ³)` jest TGP-native operator (NIE Klein-Gordon).

2. **g_eff[Φ_eq(BH)] funkcjonał** = emergent metric per emergent-metric Phase 4 {A(ψ), B(ψ), C(ψ)} family (M9_RESTRUCTURE Tier 1). Wybór f(ψ) family controls Taylor coefficients `{a_n, b_n, ξ_n}` w expansion `g_eff[ψ]` wokół ψ_0 = ψ_horyzont.

3. **Test particle (graviton perturbation, photon)** moves on geodezykach `g_eff[Φ_eq(BH)]`. **NIE couples bezpośrednio do Φ-quanta.** Zmiany QNM/photon-orbit są wynikiem zmiany `g_eff` strukturalnej, NIE scalar exchange.

4. **Family marker d²f/dψ²(ψ_0)** controls SECOND-ORDER Taylor expansion `g_eff[Φ]` wokół horizon — czyli specifically how `g_eff` deviates from pure Schwarzschild jako function of α (poly), β_q (quad), or α² (trans).

**Form-meaning mapping (per TGP_NATIVE_COMPUTATIONAL_PATTERNS §4):**

| BD-form formula | TGP-meaning |
|---|---|
| `δω_QNM/ω_GR = κ_BC · g_scalar_coupling` (Berti-Cardoso) | `δω_QNM/ω_GR = κ_QNM · d²f/dψ²(ψ_0)` z κ_QNM derived from `g_eff[Φ_eq(BH)]` second-order Taylor expansion w r_horyzont |
| `δr_ph/r_GR = κ_CH · scalar_gradient(r_ph)` (Cunha-Herdeiro) | `δε_ph²/ε_ph²_GR = κ_ε · d²f/dψ²(ψ_0)` z κ_ε derived from `g_eff[Φ_eq(photon orbit)]` second-order Taylor expansion w r_photon |
| Scalar exchange tree-level | Universal `g_eff` coupling per S05 + ax:metric-coupling — Φ-quanta są forbidden as exchange particles per FOUNDATIONS §5.1 |

**Anti-Pattern check:** ten cycle NIE używa Φ-quantum exchange picture. Wszystkie BH5/ε.1 modifications są routed przez `g_eff[Φ_eq(BH)]` nonlinear background + family-specific Taylor expansion. Per Pattern 2.5: m_Φ_observable jest environment-dependent (BH-environment), NIE universal mass.

### §0.2 — Pre-registered falsification rule

**Verbatim z YAML `contract.L1_native.falsification_rule`:**

> "Jeśli LIGO-O5 A+ ~2027 + Cosmic Explorer ~2030 ringdown stack 100+ BH events daje δω_QNM/ω_GR poza family-predicted range {polynomial=0, quadratic=2β_q·κ_QNM, transcendental=α²·κ_QNM} z 5σ confidence, OR ngEHT ~2030 photon ring stack 10+ SMBH measurements daje δε_ph²/ε_ph²_GR poza family-predicted range {polynomial=0, quadratic=β_q/9, transcendental=α²/18} z 5σ, S07 alternative f(ψ) family pre-observational discrimination FAILS → H1b verdict: framework requires deeper structural specification beyond f(ψ) freedom (e.g., explicit substrate-dynamics V(Φ) BH-background derivation)"

```
pre_registration_date: 2026-05-14
pre_registration_hash: <auto-set by git commit SHA>
recovery_scope:
  allowed_directions:
    - "f(ψ) family enumeration WITHIN S07 freedom (inherited z S07-reset PR-010)"
    - "BH-environment-specific κ_QNM, κ_ε refinement (Pattern 2.5 environment-dependent observables)"
    - "Cross-channel coupled bound calculus (BH5 + ε.1 simultaneous fit)"
  forbidden_directions:
    - "Post-hoc tuning specific f(ψ) form per channel (BH5 vs ε.1 different families)"
    - "OR-clause backstop H1c/H1d without pre-bounded scope"
    - "S05 violation (multi-Φ per channel)"
    - "Direct Φ-quantum exchange to gravitational/photon test particles (Φ-quanta forbidden per FOUNDATIONS §5.1)"
  if_recovery_exhausted: "H1b: framework wymaga deeper structural specification beyond S07 f(ψ) freedom — substrate-dynamics V(Φ) BH-background dedicated cycle, OR acceptance pre-observational discrimination NIE wystarczy + observational discrimination LIGO-O5/ngEHT becomes binding test"
```

**Anti-Lakatos compliance:** Ten cycle inherits PR-010 LOCKED recovery_scope α ∈ [-0.832, 0.832] (S07-reset polynomial channel) + DODAJE β_q ∈ [-0.4, 0.4] (quadratic channel; pre-bounded fiducial 1σ derived) + α ∈ [-0.832, 0.832] dla transcendental channel (same magnitude per polynomial). Wszystkie pre-bounded BEFORE any sympy. Brak post-hoc shifted-point recovery.

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

Q1-Q8 checklist per [[../../meta/CYCLE_LIFECYCLE.md]] §Phase 0 README template.

- [x] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1 (ASK-RULE) + §2.1 (Static Φ_eq from arbitrary source — BH-background) + §2.2 (Newton from momentum-flux — analog dla QNM/photon ring via g_eff[Φ_eq(BH)]).
- [x] **Q2 (Red flags):** Zidentyfikowane §3 red flags — primary risk Trigger C "Just reproducing standard physics without TGP mechanism" (BH5/ε.1 standard formulas są BD-form). Mitigated via §0.1 form-meaning split.
- [x] **Q3 (Inherited LOCKs):** Mapping w `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §4: M9.1'' Path 2 anchor (kategoria b per M9_RESTRUCTURE §1.4) — c_0·κ_σ=4/3 LOCK z emergent-metric Phase 4, conditional citation; family marker d²f/dψ²(ψ_0) z S07-reset Phase 2 jest TGP-native LIVE LOCK (kategoria a — no annotation needed beyond cite).
- [x] **Q4 (Standard-physics tools):** Berti-Cardoso QNM + Cunha-Herdeiro photon ring są używane jako form templates. Justification: ten cycle EXPLICIT mapuje TGP-substrate g_eff[Φ_eq(BH)] family expansion na te formy via §0.1 form-meaning split. NIE jest BD-translation — jest TGP-native derivation z BD-form result mapping per Pattern 2.2 analog.
- [x] **Q5 (m_Φ usage):** Pattern 2.5 environment-dependent observable mass — m_Φ_observable(BH-environment) różni się od m_Φ_observable(cosmological-vacuum). κ_QNM derived w BH-environment używa local m_Φ (not universal); explicit annotation w Phase 1.
- [x] **Q6 (GR limit framing):** Polynomial channel α=0 → δω_QNM/ω_GR = 0 EXACT (Schwarzschild GR limit recovered, NIE "TGP IS GR"). Family marker = 0 dla poly channel = no second-order deviation, but TGP framework strukturalnie różny (g_eff[Φ_eq] z full Φ-EOM).
- [x] **Q7 (ASK-RULE self-check):** Trigger C identified + resolved via §0.1. Triggers A, B, D not triggered (S07-reset family marker is TGP-native LIVE; predecessor cycles all clean).
- [x] **Q8 (BD-drift audit plan):** Phase FINAL spawn `bd-drift-audit` subagent per CALIBRATION_PROTOCOL §4.4 — verifies Phase 1+2 sympy NIE drifts do BD-style scalar exchange picture.

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 (BH5/ε.1 = L1 native observables; ppE projection = L2 chart per S07-reset inheritance; LIGO-O5/CE/ngEHT bounds = L3 falsification map)
- [x] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift (Trigger C explicit identified + §0.1 form-meaning split per Pattern 2.2)
- [x] Przeczytano [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 — M9.1'' jako Path 2 anchor, NIE framework (α=-4 anchor consistency check w Phase 1+2; family marker = Tier 1 framework, anchor = Tier 2)
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract (this README follows BINDING template per RESEARCH_RESTART §1.2 reactivation procedure analog)

**Sign-off:** Claudian (S07 BH5/ε.1 successor cycle) @ 2026-05-14

### §0.5 — Sympy substance plan (NEW: post-2026-05-11 audit lesson)

Per [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4: 24/104 testów w 2026-05-11 cohort to literal `T_pass = True`. To jest anti-pattern. Plan sympy dla tego cyklu:

- [x] Każdy test sympy ma **explicit pytanie fizyczne** które weryfikuje (nie self-fulfilling)
- [x] **Co najmniej 75% testów** to non-trivial symbolic manipulation (target FP% ≥ 75% per AUDIT binding threshold)
- [x] **Co najmniej 5 testów** wykonują first-principles derivation z TGP axioms (S05 / single-Φ / Φ-EOM / substrate-vacuum identification + Pattern 2.2 g_eff[Φ_eq(BH)] family expansion)
- [x] Structural declarations (S05 preservation, ax:metric-coupling preservation, c_0·κ_σ=4/3 inheritance) są raportowane **osobno** od sympy tests, NIE counted jako "8/8 sympy PASS"

**Plan per phase (target ~34 sympy tests + 6 DEC structural declarations):**

#### Phase 1 — BH5 QNM symbolic (~12 tests target)

**First-principles derivations (≥9 FP target):**
1. T1: Schwarzschild QNM baseline ω_GR = ω(M_BH) symbolic (Berti-Cardoso form template)
2. T2: g_eff[Φ_eq(BH)] second-order Taylor expansion w r_horyzont — symbolic per family {A(ψ), B(ψ), C(ψ)}
3. T3: δω_QNM/ω_GR derivation z perturbation g_eff^Schwarzschild + family marker d²f/dψ²(ψ_0) — explicit via Pattern 2.2 form-meaning per §0.1
4. T4: Polynomial family poly=0 channel → δω_QNM/ω_GR = 0 EXACT (GR limit recovered)
5. T5: Quadratic family quad=2β_q channel → δω_QNM/ω_GR = κ_QNM · 2β_q (κ_QNM derived analytically)
6. T6: Transcendental family trans=α² channel → δω_QNM/ω_GR = κ_QNM · α²
7. T7: M9.1'' anchor consistency: α=-4 → trans channel α²=16 → δω_QNM/ω_GR consistent z LIVE 8-16% range (BH5 entry inheritance)
8. T8: Cross-cycle inheritance: c_0·κ_σ=4/3 LOCK preserved w κ_QNM derivation (Path 2 anchor)
9. T9: ASK-RULE Trigger C resolution: explicit cite §0.1 form-meaning per test (anti-BD-drift documentation)

**Literature-anchored consistency (≤3 LIT target):**
10. T10: LIGO O5 A+ PSD curve at 250 Hz (Hild+2010 ET-D / LIGO-O5 design); σ_f ~0.5% per single event projection
11. T11: Cosmic Explorer ~2030 stacked 5σ projection (CE concept paper)
12. T12: Pattern 2.5 environment-dependent annotation: m_Φ_BH ≠ m_Φ_cosmological (declarative DEC, NOT counted)

**DEC structural declarations (separate, NOT counted):**
- DEC-1: S05 single-Φ axiom preserved across BH5 channel
- DEC-2: ax:metric-coupling preserved (universal g_eff coupling, NIE Φ-quantum exchange)

#### Phase 2 — ε.1 photon ring symbolic (~12 tests target)

**First-principles derivations (≥9 FP target):**
1. T1: Schwarzschild photon orbit baseline r_ph = 3M, ε_ph² = 1/27 (Cunha-Herdeiro form template)
2. T2: g_eff[Φ_eq(BH)] second-order Taylor expansion w r_photon — symbolic per family
3. T3: δε_ph²/ε_ph²_GR derivation z perturbation Schwarzschild + family marker — explicit Pattern 2.2 §0.1
4. T4: Polynomial family poly=0 channel → δε_ph²/ε_ph²_GR = 0 EXACT
5. T5: Quadratic family quad=β_q/9 channel — κ_ε = 1/9 derived analytically
6. T6: Transcendental family trans=α²/18 channel — κ_ε = 1/18 derived analytically
7. T7: M9.1'' anchor consistency: α=-4 → trans channel α²/18 = 16/18 = 8/9 ≈ 0.89 → cross-check vs +14.6% observational data point
8. T8: Cross-cycle inheritance: ε_ph² = 1/27 GR baseline + family modifications consistent z E2 cross-sector ≤ 0.5% bound
9. T9: ASK-RULE Trigger C resolution (cite §0.1 form-meaning; anti-BD-drift)

**Literature-anchored consistency (≤3 LIT target):**
10. T10: ngEHT angular resolution σ_θ ~5 μas (ngEHT concept paper)
11. T11: 10-SMBH stack discriminability projection (Sgr A*, M87*, NGC1277, etc.)
12. T12: Pattern 2.5 environment-dependent (DEC, NOT counted)

**DEC structural declarations (separate):**
- DEC-3: S05 preserved across ε.1 channel
- DEC-4: ax:metric-coupling preserved (universal g_eff[Φ] photon orbit)

#### Phase 3 — Numerical projections + family discriminability matrix (~10 tests target)

**First-principles derivations (≥7 FP target):**
1. T1: Family discriminability matrix dla LIGO-O5 σ_BH5 (3×3 family-vs-channel; quad-σ separation per family)
2. T2: ngEHT σ_ε.1 family discriminability matrix
3. T3: Coupled BH5 + ε.1 cross-channel bound calculus (joint posterior)
4. T4: LISA 2035+ EMRI ringdown ~0.1% projection (full family discriminability)
5. T5: Cross-cycle anchor matrix at α=-4 (M9.1''): BH5 + ε.1 + S07-reset α/3 + emergent-metric c_0·κ_σ=4/3 simultaneously consistent (4-way anchor)
6. T6: Recovery region α ∈ [-0.832, 0.832] PASSING test per family per channel (16 cells discriminability matrix)
7. T7: Pre-observational discriminability sigma per family pair (poly-vs-quad, poly-vs-trans, quad-vs-trans)

**Literature-anchored consistency (≤3 LIT target):**
8. T8: PSD-coupled SNR_BH5 numerical at fiducial M_BH = 30 M_⊙ (LIGO-O5)
9. T9: Photon orbit numerical at Sgr A* M = 4.1·10⁶ M_⊙ (ngEHT)
10. T10: Honest annotation: numerical fiducials post-PE, NIE first-principles values (DEC-5 declarative)

**DEC structural declarations (separate):**
- DEC-5: Numerical fiducials are post-PE values, NIE first-principles
- DEC-6: Pre-observational discrimination ceiling A− (full A requires actual data)

**Estimated total:** 34 sympy tests + 6 DEC declarations ⇒ target FP% ≥ 75% (≥26 FP), 0 hardcoded `T_pass = True` (binding).

## §1 — Phase 0: scope mapping + balance sheet

[Phase 0 balance sheet → [[./Phase0_balance.md]] (planowane post-validator-PASS + STATE.md update)]

## §2 — Phase 1: native derivation

[Phase 1 sympy → planned post-Phase-0-balance + user authorization "active" + WIP slot 1/5]

## §FINAL — Optional L2 framework reduction

[OPTIONAL — last stage; failure here does NOT invalidate Phase 1-N native results]

L2 target: ppE projection consistency check (S07-reset β_ppE^poly(α) = (15/16)·α inheritance). Reduction type: analytical-approximate. Failure disposition: L1-stands.

## §FINAL+1 — L3 falsification map check

[Wypełnione w Phase FINAL closure ceremony per L3_falsification_map z YAML contract.]

---

## Status

🟢 **CLOSED-RESOLVED 2026-05-14 sesja P3-FINAL** — claim_status **A−** (STRUCTURAL_DERIVED_NATIVE z L2 not-fully-FP-attempted). Cumulative cycle metrics: **34/34 sympy PASS, 28 FP (82.4%), 6 LIT, 6 DEC structural separate, 0 hardcoded** — incremental highest FP% w post-restart era. Single-session 4-phase heroic sprint (Phase 0 + Phase 1 BH5 + Phase 2 ε.1 + Phase 3 numerical + Phase FINAL closure). Closure ceremony: [[./Phase_FINAL_close.md]].

Phase 0 commit gate:
1. `python tooling/validate_kickoff.py research/op-S07-Phase-3-BH5-eps1-numerical-2026-05-14/README.md` → must PASS
2. PR-012 entry w `meta/PRE_REGISTERED_FALSIFIERS.md` z immutable timestamp 2026-05-14
3. User authorization "active" + WIP slot wolny (5/5 wolne post-sesja-2026-05-13)

Aż wszystkie 3 gate'y PASS — cycle pozostaje w `parking`. **Bez wypełnienia BINDING contract::,
cycle NIE może aspirować do statusu wyższego niż `PROJECTION_VERIFIED`** per CYCLE_KICKOFF_TEMPLATE §0.2.

---

**Cycle scaffolded:** 2026-05-14 (Claudian, sesja 2026-05-14 spawn per user authorization "ok zajmij się tym op-S07-Phase-3-BH5-eps1-numerical").

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 (BINDING contract)
- [[../../meta/CYCLE_LIFECYCLE.md]] §Claim status taxonomy
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §3.3 (recovery_scope), §3.4 (PR-### entry); PR-010 inherited (S07-reset)
- [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4 (sympy substance lessons)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 (three-layer L1/L2/L3 BINDING)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1 (ASK-RULE), §2.2 (Pattern Newton/momentum analog dla QNM/photon ring), §4 (form-meaning mapping)
- [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 (Bucket reading), §3 (Path 2 anchor reframing)
- [[../../tooling/validate_kickoff.py]] (technical enforcement gate)
- [[../../meta/RESEARCH_RESTART_2026-05-11.md]] (operational restart context)
- [[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] §3.4 + §6 (predecessor A− closure; upgrade path A− → A; family marker source)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] (Path 2 anchor c_0·κ_σ=4/3 LOCK; {A,B,C} family)
- [[../op-LIGO-3G-native-phase-residual-2026-05-11/Phase6_close.md]] (Δφ(f) methodology; PR-002 LIGO-O5 sensitivity)
- [[../op-bh-alpha-threshold/Phase3_results.md]] (BH5 LIVE existing prediction; T3.2)
- [[../op-eps-photon-ring/Phase3_results.md]] (ε.1 LIVE existing entries; E3.x)
- [[../op-eht/]] (M9.1'' +14.6% photon ring data point at α=-4 anchor)
