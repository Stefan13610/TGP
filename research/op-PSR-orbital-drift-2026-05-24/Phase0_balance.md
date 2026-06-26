---
title: "Phase 0 — balance sheet + pre-registration (op-PSR-orbital-drift)"
type: phase_balance
status: PHASE0_DRAFT
phase: 0
cycle: op-PSR-orbital-drift-2026-05-24
parent_motivation: "sek08a §3838-3840 native O(U³) Schwarzschild deviation"
created_date: 2026-05-24
authorization_pending: "User to LOCK Phase 0 before Phase 1 execution"
methodology_binding: "CALIBRATION_PROTOCOL.md §3.6 BINDING (incl. §3.6.13 THIRD practical application)"
anti_lakatos_lock: PRESERVED
independent_of_F8: YES (different observable, different mechanism category)
---

# Phase 0 — op-PSR-orbital-drift pre-registration

## §1 — Scope declaration

### §1.1 — Primary objective

Derive TGP-native prediction of orbital-drift residual in binary pulsars from sek08a Schwarzschild O(U³) correction, then compare to observed timing residuals (Hulse-Taylor B1913+16 + J0737-3039) to either:
- **Detect** TGP signal above observational precision (positive falsification target)
- **Bound** TGP signal below observational precision (consistency constraint)

### §1.2 — Out-of-scope (anti-Lakatos)

This cycle is **independent of F8** (dark energy acceleration). Per [[../../meta/F8_FORENSIC_2026-05-24.md]] §9 discipline:

- ❌ NOT testing F8 mechanism
- ❌ NOT continuation of γ-3/γ-3'/γ-5/γ-7
- ❌ NOT inheriting F8-related falsifier thresholds (factor 10 from γ-7)
- ✓ Independent native prediction from sek08a (concept paper)
- ✓ Independent observable (pulsar timing, NOT cosmological d_L)
- ✓ Own pre-registered threshold (observational precision-based)

### §1.3 — Mechanism category

**This cycle tests:** modifications of Schwarzschild metric at O(U³) order in TGP's emergent gravity.

**OTHER categories OUT OF SCOPE for this cycle:**
- Vacuum stress-energy (separate cycle: op-LAM-vacuum-substrate, QUEUED)
- Independent γ derivation (separate cycle: op-G-substrate-derivation, QUEUED)
- Emergent time formalism (separate cycle: op-EMT-emergent-time, DEFERRED)

---

## §2 — Mandatory reading (Phase 1 start prerequisites)

Before Phase 1 sympy work:

1. **core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex**
   - §3838-3840: Schwarzschild O(U³) deviation formula
   - thm:einstein-emergence §3493 (FRW exact, O(U²) general)
   - thm:einstein-emergence-linearized §10856 (general U(t,x))
   - thm:general-emergence §10533

2. **core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex**
   - PPN γ = β = 1 derivation (line 587-589)
   - V_M911 potential (vacuum stability)

3. **core/sek08_formalizm/sek08_formalizm.tex**
   - Linearized Einstein tensor for general U(t,x,y,z) (thm:einstein-emergence-linearized)

4. **meta/CALIBRATION_PROTOCOL.md**
   - §3.6.1-§3.6.14 BINDING
   - Especially §3.6.13 THIRD practical application (constants classification)

5. **meta/PRE_REGISTERED_FALSIFIERS.md**
   - Confirm O(U³) orbital drift NOT already pre-registered
   - Identify PR-### number for this cycle

6. **External observational references** (to be detailed in Phase 1):
   - Damour & Taylor 1992 — Hulse-Taylor binary pulsar
   - Kramer et al. 2006/2021 — double pulsar J0737-3039 tests
   - PPN constraints summary (Will 2014 review)

---

## §3 — Pre-registered falsifiers

### F-PSR-A — TGP O(U³) orbital-residual magnitude derivation

**Hypothesis:** TGP sek08a O(U³) Schwarzschild correction → specific predicted residual in periastron advance, orbital decay, or Shapiro delay (Phase 1 to derive which observable shows leading O(U³) signature).

**Pre-registered acceptance criteria:**
- **PASS**: explicit numerical prediction Δ_TGP for chosen observable, derived strictly from sek08a equations + binary pulsar U-magnitude (U = GM/(c²r) ~ 10⁻⁶ for typical NS-NS), 0 hardcoded factors
- **FAIL**: TGP prediction cannot be derived from concept paper (would indicate sek08a is incomplete for this observable)

**Computation route (Phase 1):**
1. Identify leading O(U³) observable for binary pulsars (periastron advance Δω̇, orbital decay ΔṖ_b, Shapiro delay residual Δr)
2. Compute Δ_TGP/Δ_GR ratio symbolically using sek08a metric + g_eff^μν
3. Numerical: substitute B1913+16 + J0737-3039 orbital parameters

### F-PSR-B — Comparison to observational bound

**Hypothesis:** TGP prediction either lies within observed Δ_obs error bars (consistency) or exceeds them (TGP RULED OUT).

**Pre-registered acceptance criteria:**
- **PASS** (consistency): |Δ_TGP - Δ_obs| ≤ σ_obs (1σ observational precision)
- **FAIL_DETECTED** (positive falsification): |Δ_TGP - Δ_obs| > 3σ_obs → TGP O(U³) prediction ruled out
- **FAIL_TINY** (signal below precision): |Δ_TGP| < 0.1·σ_obs → consistent BUT not a strong test (registered as future-test target)

**No factor-10 inheritance.** Threshold = observational precision (σ_obs from referenced papers).

### F-PSR-C — Cross-system independence check

**Hypothesis:** TGP prediction must be consistent between B1913+16 and J0737-3039 (different masses, eccentricities, periods).

**Pre-registered acceptance criteria:**
- **PASS**: TGP predictions for both systems independently pass F-PSR-B
- **FAIL**: prediction fails for one system but not the other (would indicate sek08a approximation breakdown)

---

## §4 — Methodology (BINDING)

### §4.1 — Strict cycle 1/2/7
- 0 hardcoded T_pass=True (per §3.6.1)
- Compute-then-compare for every substantive FP
- PARTIAL_compute budget: 1 max per cycle
- DEC budget: 3 max per cycle

### §4.2 — Anti-Lakatos discipline
- NO modification of pre-registered falsifiers post Phase 1 start
- NO re-frame of FAIL as "marginal consistency"
- NO post-hoc threshold modification

### §4.3 — Falsifier independence
- Each falsifier has standalone pass/fail
- F-PSR-A failure does NOT auto-fail B/C (could still bound something)
- F-PSR-B failure does NOT auto-fail A (derivation might be correct, just rules out TGP)
- F-PSR-C failure indicates sek08a approximation breakdown, NOT cycle invalidation

---

## §5 — LEGITIMATE inheritance (LOCKED TGP infrastructure)

| Source | Element | Status | Use in this cycle |
|--------|---------|--------|-------------------|
| sek08a (concept paper) | thm:einstein-emergence | concept paper postulate | foundation for U-expansion |
| sek08a (concept paper) | Δg_00 = -U³/6 + O(U⁴) | concept paper derivation | starting equation for F-PSR-A |
| sek08c (concept paper) | V_M911 potential | concept paper postulate | vacuum stability, IR limit |
| γ-3' Phase 3 | E_P/ℓ_P calibration | LOCKED PASS | scale fixing |
| γ-5 Phase 3 | G_eff = c³·ℓ_P²/ℏ | LOCKED PASS (gravity-synthesis PASS) | Newton matching |
| γ-7 Phase 1 | q = √(4πG)·m | LOCKED PASS | source coupling |
| Appendix E eq. 353 | m_sp ~ ℏH_0/c² | concept paper postulate | IR cutoff |

**Critical:** all inheritance is from **LOCKED PASS** items or **concept paper postulates**. NO inheritance from F8 FAIL verdicts.

---

## §6 — FORBIDDEN moves (anti-Lakatos)

| Move | Why forbidden |
|------|---------------|
| Cite γ-3/3'/5/7 FAILs as motivation | Post-hoc rescue logic |
| Cite F8_FORENSIC_2026-05-24.md as positive evidence | Forensic is information record, not motivation |
| Cite E1/E2 explorations as predictions | Explorations are scouting, not pre-registered tests |
| Lower threshold below σ_obs | Threshold inflation |
| Re-frame F-PSR-B FAIL as "marginal PASS" | Verdict tampering |
| Add new falsifiers post Phase 1 start | Anti-Lakatos LOCK violation |
| Use factor-10 threshold from γ-7 | Wrong category inheritance |

---

## §7 — Constants classification (§3.6.13 THIRD practical application)

Constants entering Phase 1 derivation:

| Constant | Category | Source | Status |
|----------|----------|--------|--------|
| G | Calibrated | γ-5 Phase 3 LOCKED | inherited |
| c | Derived (γ-5) | ℓ_P·ω_P | inherited |
| ℓ_P | Calibrated | γ-3' Phase 3 LOCKED | inherited |
| ℏ | Input (deferred γ-6) | external | input |
| m_NS (pulsar masses) | External observation | Damour-Taylor B1913+16; Kramer J0737-3039 | observational |
| a (orbit semi-major) | External observation | binary pulsar timing | observational |
| e (eccentricity) | External observation | binary pulsar timing | observational |
| P_b (orbital period) | External observation | binary pulsar timing | observational |

**No new fundamental constants introduced.** All inputs are either:
1. TGP-locked (G, c, ℓ_P from prior cycles)
2. Standard physics input (ℏ, deferred γ-6)
3. External observational (NS masses, orbital parameters)

---

## §8 — Risk register

| ID | Risk | Severity | Mitigation |
|----|------|----------|------------|
| R-PSR-1 | O(U³) ~ 10⁻¹⁸ for typical NS-NS — below observational precision | HIGH | Honest declaration in F-PSR-B; cycle outcome may be "consistent but not strong test" |
| R-PSR-2 | sek08a §3840 formula is for static spherical, binary pulsar is dynamic + 2-body | MEDIUM | Phase 1: extract 2-body O(U³) from sek08a thm:einstein-emergence-linearized (general U(t,x)) |
| R-PSR-3 | Existing PPN tests already exclude U² deviations → does O(U³) hide in PPN-residual systematics? | MEDIUM | Phase 1: compute residual AFTER subtracting GR PPN to 2PN, then check |
| R-PSR-4 | Pulsar mass uncertainties propagate to U³ → uncertainty cube | MEDIUM | Phase 1: error propagation analysis |
| R-PSR-5 | User claim "10⁻³ rad/orbita" not yet verified from sek08a derivation | MEDIUM | Phase 1: explicit magnitude derivation from sek08a — claim from sesja #8 conversation, must be verified independently |
| R-PSR-6 | sek08c V_M911 dependence on ψ → ψ_pulsar ≠ 1 may modify O(U³) coefficients | LOW | Phase 1: confirm ψ_eq for relevant U-regime |
| R-PSR-7 | BH shadow O(U³) prediction overlaps PR-012 scope | LOW | Already declared OUT OF SCOPE in §1.2 |

---

## §9 — Phase plan

| Phase | FP target | Duration | Output |
|-------|-----------|----------|--------|
| 0 | scaffold + falsifiers + reading | this | Phase0_balance.md + LOCK |
| 1 | O(U³) derivation from sek08a + magnitude verification (F-PSR-A) | 1 sesja | Phase1_sympy.py + Phase1_derivation.md |
| 2 | Comparison to observational bounds B1913+16 + J0737-3039 (F-PSR-B + C) | 1 sesja | Phase2_sympy.py + Phase2_results.md |
| FINAL | aggregate verdict + (if PASS) PR-### LOCK entry creation | 0.5 sesja | Phase_FINAL_close.md + PRE_REGISTERED_FALSIFIERS.md PR-### entry |

**Total estimate:** 2-3 sesji.

---

## §10 — Decision budget

| Budget | Cap | Used | Remaining |
|--------|-----|------|-----------|
| DEC (substantive decision) | 3 | 0 | 3 |
| PARTIAL_compute | 1 | 0 | 1 |
| Hardcoded T_pass=True | 0 | 0 | 0 |

---

## §11 — Anti-Lakatos verification (Phase 0 check)

| Item | Status |
|------|--------|
| F8 FAILs cited as motivation? | NO ✓ |
| F8_FORENSIC cited as positive evidence? | NO (only as scope discipline reference) ✓ |
| E1/E2 explorations cited as predictions? | NO ✓ |
| Threshold inherited from γ-7? | NO (uses σ_obs) ✓ |
| Cycle named γ-8? | NO (op-PSR-orbital-drift) ✓ |
| Pre-registered falsifiers BEFORE derivation? | YES (this document) ✓ |
| Standalone fail modes declared? | YES (F-PSR-A/B/C independent) ✓ |
| Independent observable? | YES (pulsar timing, not cosmological) ✓ |

**Anti-Lakatos status:** COMPLIANT.

---

## §12 — Phase 0 LOCK requirements

To LOCK Phase 0 and proceed to Phase 1, user must:

1. Read this document
2. Confirm: scope, falsifiers, thresholds, methodology, inheritance, forbidden moves
3. Authorize: "działaj" or equivalent
4. At LOCK: PR-### number assigned + entry added to PRE_REGISTERED_FALSIFIERS.md

**Awaiting authorization.**

---

## §13 — Status summary

| Field | Value |
|-------|-------|
| Phase 0 status | DRAFT (this document) |
| Authorization | PENDING user review |
| F8 independence | DECLARED ✓ |
| Anti-Lakatos | COMPLIANT ✓ |
| Inheritance | LEGITIMATE only ✓ |
| Phase 1 prerequisites | sek08a + sek08c + CALIBRATION_PROTOCOL §3.6.13 |
| Expected duration | 2-3 sesji |
| Related queued cycles | A (LAM-vacuum), D (G-substrate), C (EMT-emergent-time) |
