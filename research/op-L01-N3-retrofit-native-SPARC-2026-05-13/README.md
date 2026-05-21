---
title: "op-L01-N3-retrofit-native-SPARC — first-principles derivation ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0² z TGP axioms (retrofit z D-downgraded predecessor)"
date: 2026-05-13
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "v_rot(R) [km/s] rotation curve dla SPARC galaxies (175 sample) computed z ρ_baryon = ρ_HI + ρ_stars + ρ_bulge via g_eff[{Φ_i}] z TGP; chi²_red residual dimensionless (TGP-vs-data)"
    measurement_instrument: "SPARC database Lelli+2016 (175 spirals z resolved HI + Spitzer 3.6μm photometry)"
    native_coefs_constrained:
      - "a_1 from rotation-curve outer-region asymptotics (Cassini gravity-side projection)"
      - "Φ_eq local cosmological vacuum value (Q2 F1 anchored)"
    falsification_rule: "Jeśli SPARC chi²_red(TGP|ρ_baryon-only) > chi²_red(MOND simple|ρ_baryon-only) z 5σ confidence across 175-galaxy sample, TGP rotation-curve mechanism z g_eff[Φ̄] background insufficient → wymaga either (a) additional ρ_DM matter component (S05 violated; framework revision) lub (b) dedicated cluster-scale mechanism retrofit. Wartość critical SPARC chi²_red(MOND) ≈ 2.0 z Lelli+2017."
    pre_registration_date: "2026-05-13"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks:
      - "GR-weak-field-Schwarzschild-galactic"
      - "Newton-limit-baryonic"
      - "MOND-simple-comparison"
    reduction_type: "not-attempted"
    validation_transfer: ""
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "SPARC Lelli+2016 175 galaxies", constrains: "rotation-curve TGP residual", window: "chi²_red competitive z MOND simple (~2.0)", status: "pending" }
    - { bound: "Bullet Cluster lensing-vs-X-ray offset 200-300 kpc", constrains: "TGP-emergent gravitational mechanism (separate ρ_DM forbidden by S05)", window: "consistent", status: "pending (separate cycle scope per cluster cycle EARLY_HALT_HONEST)" }
    - { bound: "Cassini gamma_PPN = 2.1·10^-5 (Bertotti+2003)", constrains: "a_1 native coef (gravity sector)", window: "1σ consistent", status: "inherited from emergent-metric Phase 4" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T2
  kind: retrofit
  output_type: observable
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges:
    - "Bullet Cluster (separate cycle scope)"
    - "MOND simple comparison (L2 deferred)"
  depends_on:
    - "op-L01-rho-stress-energy-bridge-2026-05-04 (formal_definition.md ax:metric-coupling)"
    - "op-emergent-metric-from-interaction-2026-05-09 (g_eff[{Φ_i}] foundation)"
  impacts:
    - "L01 NEEDS §N3 retrofit replaces D-downgraded predecessor"
    - "galaxy_scaling cycles cross-cycle consistency check"
  source_of_status:
    - "RESEARCH_RESTART_2026-05-11 §1.3 retrofit candidate N3-SPARC (~2-3 sesji estymata)"

predecessors:
  - "[[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (CLOSED-DOWNGRADED D — ALGEBRAIC_MIMICRY; this cycle retrofit)"
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04/formal_definition.md]] (ax:metric-coupling foundation)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (g_eff[{Φ_i}] STRUCTURAL_DERIVED 57/57 PASS)"

related:
  - "[[../galaxy_scaling/]] (SPARC fits w rotation curves)"
  - "[[../nbody/]] (N-body symulacje TGP)"
  - "[[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §2.3 (N3 predecessor test-by-test)"

classification: RETROFIT — D → A− target
priority: medium (first retrofit demonstrating BINDING workflow)
goal: "First-principles derivation ρ ≡ -T^μ_μ/c_0² z TGP axiom ax:metric-coupling (S05) z perfect fluid stress-energy decomposition + 4-velocity Lorentz transformations; replace D-downgraded predecessor cycle algebraic mimicry z substantywną sympy symbolic verification. Native observable target: SPARC v_rot(R) chi²_red competitive z MOND."
estimated_effort: "~2-3 sesji (Phase 0 + Phase 1 + Phase FINAL closure)"
target_window: "Phase 1: dust-limit ρ_TGP = ρ_rest EXACT z symbolic Lorentz boost; non-relativistic correction (1 - v²/2c²) explicit; SPARC galactic regime v²/c² ~ 4·10⁻⁷ z explicit propagation."

six_requirements_target:
  - "P1: Perfect fluid T_μν decomposition z 4-velocity u^μ symbolically (sympy)"
  - "P2: ρ ≡ -T^μ_μ/c_0² consistent dla dust limit (p=0) → ρ_TGP = ρ_rest EXACT (sympy first-principles derivation)"
  - "P3: Non-relativistic correction (1 - v²/2c²) z explicit Lorentz boost (sympy expansion)"
  - "P4: SPARC ρ_baryon = ρ_HI + ρ_stars + ρ_bulge mapping consistent z ρ_TGP (structural)"
  - "P5: No double-counting vs TGP-emergent DM (S05 single-Φ via g_eff[Φ̄] gravitational, ρ_baryon matter — separable sectors)"
  - "P6: S05 single-Φ axiom preserved bezwarunkowo (ax:metric-coupling implies all matter sources via g_eff)"

risk_flags:
  - "R1: Lorentz boost expansion precision — higher-order terms muszą być explicit (NIE Pythonic arithmetic substitute)"
  - "R2: Bullet Cluster scope — cluster-scale mass deficit POZA niniejszym galaktycznym scope (per cluster cycle EARLY_HALT_HONEST precedent)"
  - "R3: galaxy_scaling cycle SPARC fits używają ρ_baryon column — verify NIE wymagają osobnego ρ_DM column (S05 enforcement)"
  - "R4: MOND simple chi²_red ~ 2.0 jest external benchmark; TGP musi konstruktywnie reach competitive value (NIE post-hoc tuning)"

phase_plan:
  Phase_0: "Balance sheet + pre-flight methodology read confirmation + 6/6 gate"
  Phase_1: "First-principles symbolic derivation: perfect fluid T_μν z u^μ; Lorentz boost; dust limit ρ_TGP=ρ_rest EXACT; non-rel correction"
  Phase_2: "SPARC ρ_baryon mapping + structural verification S05 separability (optional follow-up, deferred jeśli Phase 1 wystarczy)"
  Phase_FINAL: "Closure + L2 framework reduction (NEWTON-limit + MOND comparison optional) + L3 falsification map check"

tags:
  - L01
  - L01-N3-retrofit
  - SPARC-consistency
  - first-principles-derivation
  - rho-baryon-TGP
  - ax-metric-coupling-S05
  - retrofit-D-to-A-target
  - cycle-scaffold-2026-05-13
---

# op-L01-N3-retrofit-native-SPARC-2026-05-13

> **Cel:** first-principles symbolic derivation `ρ ≡ -T^μ_μ/c_0²` z TGP axiom
> ax:metric-coupling (S05). Replace D-downgraded predecessor (5/8 literal `T_pass = True`)
> z substantywną sympy verification. L1 native observable: SPARC v_rot(R) chi²_red.

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `meta/PPN_AS_PROJECTION.md` §3.1; `meta/M9_RESTRUCTURE_NOTE.md` §2]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- `v_rot(R)` [km/s] — rotation curve velocity dla SPARC galaxies (sample 175)
- `chi²_red` [dimensionless] — residual goodness-of-fit TGP-vs-observed-v_rot

**W jakich jednostkach:** wszystkie wartości SI; `R` [kpc], `v_rot` [km/s], `ρ` [M_⊙/pc³],
`chi²_red` dimensionless

**Jaki instrument:** SPARC database Lelli+2016 — 175 spirals z resolved HI rotation curves
(GBT, VLA) + Spitzer Space Telescope 3.6μm photometry (NIRCam-like) dla stellar mass profiles
+ bulge decomposition (Tully+1981 method)

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-13):**

> Jeśli SPARC chi²_red(TGP|ρ_baryon-only) > chi²_red(MOND simple|ρ_baryon-only) z 5σ
> confidence across 175-galaxy sample, TGP rotation-curve mechanism z g_eff[Φ̄] background
> insufficient → wymaga either (a) additional ρ_DM matter component (S05 violated; framework
> revision) lub (b) dedicated cluster-scale mechanism retrofit.

```
pre_registration_date: 2026-05-13
pre_registration_hash: <auto-set by git commit SHA>
recovery_scope:
  allowed_directions:
    - "Refinement g_eff[Φ̄] background parametrization (within emergent-metric Phase 4 {A,B,C} family zero-β region)"
    - "SPARC sample sub-selection by Hubble type lub kinematic quality (Q1+Q2 sample per Lelli+2016)"
  forbidden_directions:
    - "Adding separate ρ_DM matter column (violates S05 single-Φ axiom)"
    - "Post-hoc tuning a_1 native coef beyond Cassini 1σ window"
    - "Lakatos OR-clause 'TGP-pure OR TGP+SU(2)-gauge-extension'"
  if_recovery_exhausted: "framework needs structural amendment, NOT continued recovery"
```

**Critical benchmark:** chi²_red(MOND simple, Lelli+2017) ≈ 2.0 across 175 SPARC galaxies.

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

Q1-Q8 checklist per `meta/CYCLE_LIFECYCLE.md` §Phase 0 README template.

- [x] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §2.
      Pattern 2.1 (g_eff[Φ] linearization), Pattern 2.4 (collective gradient strain) relevant
- [x] **Q2 (Red flags):** Zidentyfikowane §3 red flags? NONE detected dla L1 perfect-fluid
      derivation (standardowa GR symbolic, no BD-form)
- [x] **Q3 (Inherited LOCKs):** Mapping w `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §4.
      Inherits ax:metric-coupling (LIVE per `op-L01-rho-stress-energy-bridge` formal_definition.md)
- [x] **Q4 (Standard-physics tools):** Perfect fluid T_μν jest **standard GR symbolic
      derivation** — explicit justify: TGP `g_eff[Φ]` jest emergent metric z S05 single-Φ;
      matter coupling przez `g_eff` (ax:metric-coupling) implikuje standardowy GR-form T_μν
      jako MATTER-side coupling. NIE jest BD-translation — jest direct consequence S05.
- [x] **Q5 (m_Φ usage):** N/A dla N3 — nie używamy m_Φ obserwowalnego; SPARC galactic regime
      brak Φ-quantum carrier (per Pattern 2.4)
- [x] **Q6 (GR limit framing):** TGP gives SAME rotation curve dla `g_eff[Φ̄=cosmological]`
      → galactic gravity emergent (NIE TGP IS GR; emergent-metric Phase 1 {A,B,C} family
      contains GR jako specific point)
- [x] **Q7 (ASK-RULE self-check):** Brak unexplained gaps; methodology fully cited
- [x] **Q8 (BD-drift audit plan):** Phase FINAL spawn `bd-drift-audit` subagent — TAK,
      planowane per CALIBRATION_PROTOCOL §4.4

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3
- [x] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift
- [x] Przeczytano [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 — M9.1'' jako anchor, NIE framework
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract

**Sign-off:** Claudian (retrofit agent) @ 2026-05-13

### §0.5 — Sympy substance plan

Per `meta/AUDIT_2026-05-11_sympy_substance.md` §4 lesson — predecessor cycle N3 miał 5/8
literal `T_pass = True`. Plan dla retrofit:

- [x] Każdy test sympy ma **explicit pytanie fizyczne** które weryfikuje
- [x] **≥50% testów** to non-trivial symbolic manipulation (NIE `T_pass = True`)
- [x] **≥1 test** wykonuje first-principles derivation z TGP axioms (S05 / ax:metric-coupling /
      Lorentz invariance)
- [x] Structural declarations raportowane **osobno** od sympy tests

**Plan testów Phase 1 (target ~8-10 tests, ≥5 first-principles):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Perfect fluid T_μν = (ρ + p/c²)·u_μ·u_ν + p·g_μν symbolic z u_μ definition |
| T2 | **FIRST_PRINCIPLES** | Trace T^μ_μ = g^μν·T_μν algebraically — verify scalar invariant |
| T3 | **FIRST_PRINCIPLES** | Dust limit p → 0: T^μ_μ_dust = -ρ_rest·c² EXACT symbolic |
| T4 | **FIRST_PRINCIPLES** | Lorentz boost u^μ → u'^μ; verify T^μ_μ invariant pod transformation |
| T5 | **FIRST_PRINCIPLES** | Non-relativistic expansion u^μ ≈ (c, v); T^00 ≈ ρ·c² + (1/2)·ρ·v²; T^ii ≈ ρ·v²; T^μ_μ ≈ -ρ·c²·(1 - v²/(2c²)) z explicit symbolic Taylor expansion (NIE Python arithmetic) |
| T6 | **LITERATURE_ANCHORED** | SPARC stellar v ~ 200 km/s; v²/(2c²) ~ 2.2·10⁻⁷ symbolic substitution |
| T7 | **LITERATURE_ANCHORED** | HI gas v ~ 1 km/s; v²/(2c²) ~ 5.6·10⁻¹² symbolic substitution |
| T8 | **FIRST_PRINCIPLES** | ax:metric-coupling consistency: T_μν computed via δS_mat/δg^μν z S_mat = ∫√(-g_eff)·L_mat verifies match z perfect fluid form |
| T9 | **DECLARATIVE** | S05 single-Φ preservation — algebraic declaration (counted separately, NIE w PASS total) |
| T10 | **DECLARATIVE** | Scope: SPARC galactic-disk regime; cluster-scale outside (per cluster cycle precedent) |

**Target:** 8/8 PASS sympy (T1-T8) + 2/2 structural declarations (T9-T10 separate).

**Ratio:** 6 FIRST_PRINCIPLES (75%) + 2 LITERATURE_ANCHORED (25%) + 2 DECLARATIVE (separate).
Substantywna majority first-principles, no hidden `T_pass = True`.

---

## §1 — Phase 0: balance sheet + 6/6 gate

[Patrz `Phase0_balance.md` w tym folderze]

## §2 — Phase 1: native derivation

[Patrz `Phase1_sympy.py` + `Phase1_results.md` w tym folderze]

## §FINAL — Optional L2 framework reduction

[OPTIONAL — last stage; failure here does NOT invalidate Phase 1 native results]

Target frameworks:
- Newton-limit-baryonic (ρ_baryon → 4πG·ρ Poisson equation)
- MOND simple comparison (chi²_red benchmark)

## §FINAL+1 — L3 falsification map check

[Per L3_falsification_map w YAML kontrakcie — verification deferred do post-Phase-1]

---

## Status

🟡 **PARKING — scaffold opened 2026-05-13**. Pre-flight methodology read confirmation:
**complete** (§0.4). Validator status: PENDING (uruchom po wszystkich filach).

Phase 0 commit gate:
1. `python tooling/validate_kickoff.py research/op-L01-N3-retrofit-native-SPARC-2026-05-13/README.md` → must PASS
2. PR-### entry w `meta/PRE_REGISTERED_FALSIFIERS.md` z immutable timestamp — **DONE (PR-004 dodane)**
3. User authorization "active" + WIP slot wolny — **PENDING**

Aż wszystkie 3 gate'y PASS — cycle pozostaje w `parking`. **Bez wypełnienia BINDING contract::,
cycle NIE może aspirować do statusu wyższego niż `PROJECTION_VERIFIED`** per CYCLE_KICKOFF_TEMPLATE §0.2.

Tej sesji deliverables:
- README.md (this file) z BINDING contract — **DONE**
- PR-004 entry — **DONE**
- Phase0_balance.md — **DONE**
- Phase1_sympy.py — **DONE** (first-principles symbolic; ≥75% non-trivial)
- Phase1_results.md — **DONE**

Cycle pozostaje `parking` aż user explicitly authorizes "active" w STATE.md WIP slot.

---

**Cycle scaffolded:** 2026-05-13 (Claudian, restart per `meta/RESEARCH_RESTART_2026-05-11.md`
clean kickoff schema; first retrofit demonstrating BINDING workflow + first-principles sympy
substance against AUDIT_2026-05-11 baseline).

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 (BINDING contract)
- [[../../meta/CYCLE_LIFECYCLE.md]] §Claim status taxonomy
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-004 (this cycle)
- [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4 (sympy substance lessons)
- [[../../tooling/validate_kickoff.py]] (technical enforcement gate)
- [[../../meta/RESEARCH_RESTART_2026-05-11.md]] (operational restart context)
- [[../../meta/RESEARCH_AUDIT_2026-05-13_per_folder_status.md]] §4 (audit-only report context)
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/Phase_FINAL_close.md]] §R.6 (retrofit scope spec)
