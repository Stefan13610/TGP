---
title: "op-L05-mass-exponent-k-alpha-d — formal derivation k(α, d) z finite-energy condition + m_obs vs M_full distinction"
date: 2026-05-16
type: research-cycle
folder_status: parking
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "Lepton mass ratios m_μ/m_e, m_τ/m_e [dimensionless]; soliton total energy M_full [GeV] vs tail-coupling observable mass m_obs [GeV] — distinction structurally formalized z explicit exponents k_full(α,d) and k_obs(α,d) jako functions of kinetic prefactor power α and spatial dimension d"
    measurement_instrument: "PDG 2024 lepton mass ratios (m_μ/m_e = 206.7682, m_τ/m_e = 3477.23); inherits from why_n3 numerical scan r3_observable_vs_full_mass.py (α∈[0.75, 2.5] sweep)"
    native_coefs_constrained:
      - "k_full(α, d) — exponent volumetric (M_full ∝ A^k_full) z virial scaling argument"
      - "k_obs(α, d) — exponent tail-coupling (m_obs ∝ A_tail^k_obs) z asymptotic matching"
      - "Reconciliation map: relationship k_full ↔ k_obs przez wewnętrzną strukturę (core-tail) solitonu"
    falsification_rule: "Jeśli analytical k_obs(α=2, d=3) ≠ 3 (matching empirical p=5-α dla α=2) z mathematical rigor sympy expansion, lub jeśli k_full(α=1, d=3) ≠ 4 (matching LP-4 9/9 PASS), L05 reconciliation FAILS → wymaga (a) extending energy functional beyond canonical K(φ)·(∂φ)² + V(φ) form lub (b) revising m_obs definition. Successful reconciliation = both LP-4 (Możliwość A) AND R3 (Możliwość A) jednocześnie consistent przez m_obs vs M_full distinction."
    pre_registration_date: "2026-05-16"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks:
      - "Standard scalar-field soliton scaling theory (Derrick's theorem)"
      - "Sobolev embedding critical exponent p_crit = (d+2)/(d-2)"
    reduction_type: "not-attempted"
    validation_transfer: ""
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "PDG m_μ/m_e = 206.7682 (8e-7 precision)", constrains: "k_obs(α=2, d=3) for canonical TGP", window: "α=2 → p=3 exact within 0.099%", status: "inherited from r3_alpha2_full_closure.py 6/6 PASS" }
    - { bound: "PDG m_τ/m_e = 3477.23 (8e-7 precision)", constrains: "k_obs(α=2, d=3) + Koide K=2/3", window: "−0.085% diff inherited", status: "inherited from r3_alpha2_full_closure.py" }
    - { bound: "LP-4 9/9 PASS volumetric mass argument (lp4_mass_exponent_verification.py)", constrains: "k_full(α=1, d=3) integer constraint via finite-energy condition", window: "k=4 unique integer in d=3 for K=g² kinetic", status: "inherited but reinterpreted (k_full vs k_obs distinction)" }

# ============== END KICKOFF CONTRACT ==============

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
    - "L08 (kink fermion closure) — uses k_obs for emergent Dirac propagator residue"
    - "L06 (m_X axion mass locked) — separate cycle; uses different scaling"
  depends_on:
    - "op-L04-ODE-canonicalization-2026-05-04 (α=2 canonical decision)"
    - "research/why_n3 Phase 1-5 closed (numerical p(α)=5-α discovery)"
    - "audyt/L05_mass_exponent_drift README.md (problem statement + 3 możliwości A/B/C)"
  impacts:
    - "audyt/L05 → CLOSED-RESOLVED if Możliwość A confirmed"
    - "research/mass_scaling_k4 → renamed/split (k_obs vs k_full)"
    - "core/sek08b_ghost_resolution thm:B1'' → reinterpreted (specific to M_full, not m_obs)"
    - "research/why_n3 CORRECTIONS_2026-05-01.md → analytical backbone added"
  source_of_status:
    - "audyt/PRIORITY_MATRIX.md klaster D L05 P2 OPEN"
    - "audyt/L05_mass_exponent_drift/README.md §Rekomendacja Phase 1"

predecessors:
  - "[[../op-L04-ODE-canonicalization-2026-05-04/]] (α=2 canonical via thm:D-uniqueness)"
  - "[[../why_n3/r3_observable_vs_full_mass.py]] (numerical p(α)=5-α discovery)"
  - "[[../mass_scaling_k4/]] (paused; LP-4 9/9 PASS for K=g² specifically)"

related:
  - "[[../../audyt/L05_mass_exponent_drift/README.md]] (problem statement)"
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] (downstream — fermion sector closure)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-2 (BINDING contract)"

classification: DERIVATION — L05 audit reconciliation
priority: high (P2 OPEN klaster D ontology; bridge between why_n3 numerical and LP-4 formal)
goal: "First-principles formal derivation of mass exponent k(α, d) for radial soliton z K(φ)=K_geo·φ^α + φ⁴ potential via (i) virial scaling argument for full energy M_full ∝ A^k_full, (ii) asymptotic tail matching for observable m_obs ∝ A_tail^k_obs, (iii) explicit reconciliation map showing LP-4 (k=4 for K=g², via M_full) AND R3 (p=5-α via m_obs) jednocześnie consistent przez distinction. Closes Możliwość A from L05 audit."
estimated_effort: "~1-2 sesji (Phase 0 + Phase 1 symbolic derivation + Phase FINAL)"
target_window: "Phase 1: explicit sympy derivation k_full(α,d) z Derrick virial scaling; k_obs(α,d) z asymptotic linearization wokół φ=1; reconciliation theorem (k_full, k_obs related przez core-tail matching). Verification: α=2,d=3 → k_obs=3 EXACT; α=1,d=3 → k_full=4 EXACT."

six_requirements_target:
  - "P1: Symbolic derivation k_full(α, d) z virial scaling argument E_kin = E_pot dla K(φ)=K_geo·φ^α + V_eff(φ)~φ⁴ (sympy explicit Lagrangian + scaling transformation)"
  - "P2: Symbolic derivation k_obs(α, d) z asymptotic tail matching φ→1+δ, A_tail amplitude z linearized EOM"
  - "P3: Reconciliation theorem: k_full(α=1, d=3) = 4 EXACT (matches LP-4); k_obs(α=2, d=3) = 3 EXACT (matches R3 empirical p=5-α dla α=2); structural relationship k_obs = k_full − (d-2)·σ_match dla σ_match z core-tail matching"
  - "P4: m_obs vs M_full distinction operationally formalized — explicit symbolic difference (analogous do ADM vs Komara mass w GR, bare vs renormalized w QFT)"
  - "P5: L05 audit Możliwość A constructively confirmed (LP-4 argument applies do M_full; R3 formula stosuje do m_obs); Możliwości B i C eliminated z analytical evidence"
  - "P6: S05 single-Φ axiom preserved (single soliton profile φ(r); two observables k_full and k_obs są DWA projekcje tej samej struktury)"

risk_flags:
  - "R1: Virial scaling argument requires specifying V_eff(φ) form — use canonical φ⁴ z TGP-FOUNDATIONS §3; document if other V forms give different k_full"
  - "R2: Asymptotic tail expansion order — verify next-to-leading corrections nie zmieniają k_obs"
  - "R3: Sobolev critical exponent connection — IF k_obs = p_crit(d) for some α, may indicate deeper structural reason (note as PASS-WITH-FLAG)"
  - "R4: numerical p(α)=5-α from R3 holds EXACTLY only dla α∈{1,2}; values w between (e.g., α=1.5, p=3.428 vs 5-α=3.5, diff +2.1%) wymagają higher-order correction — document scope"

phase_plan:
  Phase_0: "Balance sheet + 6/6 gate + scope (canonical V=φ⁴, d=3 primary; d-general derivation parametrized)"
  Phase_1: "First-principles symbolic: Derrick virial scaling → k_full(α,d); asymptotic tail linearization → k_obs(α,d); reconciliation theorem"
  Phase_FINAL: "Closure + L2 reduction (Sobolev embedding optional) + L3 falsification check + L05 audit closure note"

tags:
  - L05
  - mass-exponent
  - virial-scaling
  - asymptotic-tail
  - m-obs-vs-M-full
  - first-principles
  - audit-closure-candidate
  - cycle-scaffold-2026-05-16
---

# op-L05-mass-exponent-k-alpha-d-2026-05-16

> **Cel:** First-principles symbolic derivation of mass exponent `k(α, d)` for
> radial soliton z `K(φ) = K_geo·φ^α + V_eff(φ) ~ φ⁴`. Distinction `M_full` (volumetric)
> vs `m_obs` (asymptotic tail-coupling) operationally formalized. Reconciles
> **LP-4** (k=4 for K=g²) AND **R3** (p=5-α empirically) jednocześnie przez
> Możliwość A z audyt/L05_mass_exponent_drift.

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `meta/PPN_AS_PROJECTION.md` §3.1]

### §0.1 — Native observable target

**Co fizycznie liczymy:**

- `k_full(α, d)` [dimensionless] — exponent w `M_full ∝ A^k_full` z volumetric integral
- `k_obs(α, d)` [dimensionless] — exponent w `m_obs ∝ A_tail^k_obs` z tail coupling
- `m_μ/m_e`, `m_τ/m_e` [dimensionless] — PDG lepton mass ratios as benchmark

**Instrument:** PDG 2024 ratios + numerical r3_observable_vs_full_mass.py α-scan

### §0.2 — Pre-registered falsification rule

**Decision rule WRITTEN BEFORE any calculation (2026-05-16):**

> Jeśli analytical k_obs(α=2, d=3) ≠ 3 (matching empirical p=5-α dla α=2) z
> mathematical rigor sympy expansion, lub jeśli k_full(α=1, d=3) ≠ 4 (matching
> LP-4 9/9 PASS), L05 reconciliation FAILS → wymaga (a) extending energy
> functional beyond canonical K·(∂φ)² + V(φ) form lub (b) revising m_obs definition.

```
pre_registration_date: 2026-05-16
pre_registration_hash: <auto-set by git commit SHA>
recovery_scope:
  allowed_directions:
    - "Higher-order corrections in asymptotic tail expansion (sub-leading O(1/r) terms)"
    - "Refinement of A_tail vs A_core matching condition (logarithmic corrections allowed)"
  forbidden_directions:
    - "Adding free fitting parameter do mass formula post-hoc"
    - "Lakatos OR-clause 'TGP-canonical OR alternative K form'"
    - "Redefining 'observable mass' to match desired exponent"
  if_recovery_exhausted: "Możliwość B (R3 fitting artifact) lub Możliwość C (LP-4 specific) — honest verdict"
```

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

Q1-Q8 checklist per `meta/CYCLE_LIFECYCLE.md` §Phase 0 README template.

- [x] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md`.
      Pattern 2.6 (Derrick scaling) + Pattern 2.7 (asymptotic matching) — both relevant
- [x] **Q2 (Red flags):** NONE — standardowy GR/soliton variational analysis, no BD-drift
- [x] **Q3 (Inherited LOCKs):** α=2 canonical via L04 thm:D-uniqueness (LIVE)
- [x] **Q4 (Standard-physics tools):** Derrick virial argument + linearized asymptotic expansion
      jest standard variational physics; native-relevance: K(φ)=K_geo·φ^α form is TGP-specific
      generalization beyond canonical K=const scalar field
- [x] **Q5 (m_Φ usage):** N/A — pracujemy z classical soliton structure
- [x] **Q6 (GR limit framing):** N/A — flat space d=3 spatial (Minkowski background)
- [x] **Q7 (ASK-RULE self-check):** Methodology fully cited; gaps documented w R-flags
- [x] **Q8 (BD-drift audit plan):** Phase FINAL spawn `bd-drift-audit` subagent if applicable

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [x] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3
- [x] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract
- [x] Przeczytano [[../../audyt/L05_mass_exponent_drift/README.md]] — problem statement
- [x] Przeczytano [[../why_n3/CORRECTIONS_2026-05-01.md]] §RESOLUTION — m_obs vs M_full insight

**Sign-off:** Claudian (theoretical physics agent) @ 2026-05-16

### §0.5 — Sympy substance plan

Per `meta/AUDIT_2026-05-11_sympy_substance.md` §4 — substance-first lessons:

- [x] Każdy test sympy ma **explicit pytanie fizyczne** które weryfikuje
- [x] **≥75% testów** to non-trivial symbolic manipulation
- [x] **≥1 test** wykonuje first-principles derivation from scaling/variational principle
- [x] Structural declarations raportowane **osobno** od sympy tests

**Plan testów Phase 1 (target ~10-12 tests, ≥8 first-principles):**

| Test | Klasa | Pytanie fizyczne |
|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Action functional S[φ] = ∫d^d x [½K_geo·φ^α·(∂φ)² + V(φ)] explicit symbolic z V=λφ⁴/4 |
| T2 | **FIRST_PRINCIPLES** | Derrick scaling: φ(x) → φ(x/L); compute E_kin(L), E_pot(L) explicit z α |
| T3 | **FIRST_PRINCIPLES** | Stationarity ∂E/∂L = 0 → L_star^2 z A z α (virial condition) |
| T4 | **FIRST_PRINCIPLES** | M_full = E(L_star) ∝ A^k_full; derive k_full(α, d) symbolically |
| T5 | **FIRST_PRINCIPLES** | Verify k_full(α=1, d=3) = 4 EXACT (matches LP-4) |
| T6 | **FIRST_PRINCIPLES** | Asymptotic tail: φ = 1 + δ z linearized EOM K_geo·∇²δ = V''(1)·δ |
| T7 | **FIRST_PRINCIPLES** | Tail solution δ(r) = A_tail · exp(-m·r)/r^((d-1)/2); m² = V''(1)/K_geo |
| T8 | **FIRST_PRINCIPLES** | Core-tail matching: A_tail ~ A^σ_match(α, d) z saddle point inner-outer |
| T9 | **FIRST_PRINCIPLES** | m_obs ∝ A_tail^k_obs; k_obs derivation z external Yukawa coupling |
| T10 | **FIRST_PRINCIPLES** | Verify k_obs(α=2, d=3) = 3 EXACT (matches R3 empirical p=5-α dla α=2) |
| T11 | **FIRST_PRINCIPLES** | Verify k_obs(α=1, d=3) = 4 EXACT (matches LP-4 dla α=1; consistent with R3 p=5-α dla α=1) |
| T12 | **LITERATURE_ANCHORED** | Compare with Derrick (1964) virial theorem — TGP K(φ)=φ^α extension |
| T13 | **DECLARATIVE** | S05 single-Φ preservation — single soliton profile (separate, NIE w PASS total) |

**Target:** 12/12 PASS sympy (T1-T12) + 1 structural declaration (T13 separate).

**Ratio:** 11 FIRST_PRINCIPLES (91.7%) + 1 LITERATURE_ANCHORED (8.3%) + 1 DECLARATIVE separate.
Substantywna majority first-principles, no hidden `T_pass = True`.

---

## §1 — Phase 0: balance sheet + 6/6 gate

[Patrz `Phase0_balance.md` w tym folderze]

## §2 — Phase 1: native derivation

[Patrz `Phase1_sympy.py` + `Phase1_results.md` w tym folderze]

## §FINAL — Optional L2 framework reduction

[OPTIONAL — last stage; failure here does NOT invalidate Phase 1 native results]

Target frameworks:
- Derrick's theorem original form (K=const) — limit α→0 case (if applicable)
- Sobolev embedding critical exponent comparison

## §FINAL+1 — L3 falsification map check

[Per L3_falsification_map w YAML kontrakcie — verification z PDG lepton mass ratios via inherited why_n3 data]

---

## Status

🟡 **PARKING — scaffold opened 2026-05-16**. Pre-flight methodology read confirmation:
**complete** (§0.4). Validator status: PENDING.

Phase 0 commit gate:
1. `python tooling/validate_kickoff.py research/op-L05-mass-exponent-k-alpha-d-2026-05-16/README.md` → must PASS
2. PR-### entry w `meta/PRE_REGISTERED_FALSIFIERS.md` z immutable timestamp — **PENDING (PR-013 proposed)**
3. User authorization "active" + WIP slot wolny — **GRANTED via session prompt "wybrać kolejny projekt z research i przystapić do jego realizacji w ramach TGP_v1" 2026-05-16**

This session deliverables target:
- README.md (this file) z BINDING contract — **DONE**
- Phase0_balance.md — **PLANNED**
- Phase1_sympy.py — **PLANNED** (first-principles symbolic; ≥75% non-trivial)
- Phase1_results.md — **PLANNED**
- Phase_FINAL_close.md — **PLANNED if Phase 1 PASS**

---

**Cycle scaffolded:** 2026-05-16 (Claudian, theoretical physics expert role; opens after
L01 N1-N5 retrofit cycles all CLOSED-RESOLVED A− 2026-05-13; addresses next P2 OPEN audit
item from klaster D ontology).

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 (BINDING contract)
- [[../../meta/CYCLE_LIFECYCLE.md]] §Claim status taxonomy
- [[../../audyt/L05_mass_exponent_drift/README.md]] (problem statement)
- [[../../audyt/PRIORITY_MATRIX.md]] klaster D L05 P2 OPEN
- [[../op-L04-ODE-canonicalization-2026-05-04/]] (α=2 canonical predecessor)
- [[../why_n3/CORRECTIONS_2026-05-01.md]] (numerical p(α)=5-α discovery)
