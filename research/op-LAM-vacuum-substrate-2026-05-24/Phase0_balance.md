---
title: "Phase 0 — balance sheet + pre-registration (op-LAM-vacuum-substrate)"
type: phase_balance
status: PHASE0_LOCKED
phase: 0
cycle: op-LAM-vacuum-substrate-2026-05-24
parent_motivation: "Appendix E eq. 207 + eq. 353 + eq. 365 + sek08c V_M911 — concept paper Phi-phonon DE candidate already pre-predicted"
created_date: 2026-05-24
locked_date: 2026-05-25
locked_by: "User 2026-05-25: 'działaj Phase 1' → Phase 0 LOCK + Phase 1 execution authorized"
PR_reserved: "PR-018 reserved for Phase FINAL closure"
methodology_binding: "CALIBRATION_PROTOCOL.md §3.6 BINDING (incl. §3.6.13 THIRD practical application)"
anti_lakatos_lock: PRESERVED
independent_of_F8_cycles: YES (different mechanism category: vacuum stress-energy, NOT kinematic/clumping/c-variation)
tests_F8_indirectly: YES (Λ_eff IS Ω_DE) — but pre-registered as standalone mechanism test
expected_duration: "3-5 sesji"
---

# Phase 0 — op-LAM-vacuum-substrate pre-registration

## §1 — Scope declaration

### §1.1 — Primary objective

Derive TGP-native effective cosmological constant Λ_eff from Phi-substrate vacuum stress-energy mechanism (Appendix E eq. 207 with V_M911 substrate potential + IR-regulated quantum loop contribution). Test against observed Λ_obs at pre-registered factor-10 threshold for magnitude; sign analysis; equation-of-state w_DE consistency.

### §1.2 — Mechanism category (this cycle ONLY)

**This cycle tests:** vacuum stress-energy contribution to T_μν^(Φ) from Phi-phonon ground state + IR-regulated 1-loop quantum corrections (Appendix E §first-iteration extension).

**OTHER F8-related categories OUT OF SCOPE for this cycle:**
- Kinematic mechanisms (γ-3 R(t)=c·t, γ-3' E_P/ℓ_P, γ-5 quasi-static — all LITERAL_FAIL, not cited as motivation)
- Geometric mechanisms (γ-7 mass-clumping pair-overlap — HALT-B, not cited)
- Emergent time effects (op-EMT-emergent-time DEFERRED separate cycle)
- Independent γ derivation (op-G-substrate-derivation QUEUED separate cycle — affects this cycle's prediction independence status, see §3.5)
- NS surface O(U³) (PR-017 op-PSR-orbital-drift CLOSED B+)

### §1.3 — Critical scope distinction: structural consistency vs. independent prediction

**This cycle, in current form, tests STRUCTURAL CONSISTENCY:**
- γ is currently CALIBRATED to observation via Appendix E eq. 353 (m_sp ~ ℏH_0/c²)
- Setting Λ_eff = γ/12 with γ = H_0²/c² gives consistency check, NOT independent prediction
- Outcome: factor-25 envelope per [[../../meta/F8_FORENSIC_2026-05-24.md]] §6 (Λ_TGP/Λ_obs ≈ 1/25)

**For this to become TRUE INDEPENDENT PREDICTION:**
- Cycle D (op-G-substrate-derivation) must succeed in deriving γ from {ℓ_P, c, ℏ, V_M911 parameters} without H_0 input
- D's outcome is PREREQUISITE for upgrading this cycle's claim_status beyond "structural consistency"

**This cycle declares (explicit):**
- PRIMARY scope: structural consistency check + sign + w_DE + loop corrections (independent of D's outcome)
- SECONDARY scope: IF D delivers independent γ during this cycle's execution timeline, upgrade interpretation
- ANTI-LAKATOS: do NOT cite F8_FORENSIC envelope factor-25 as positive prediction — envelope is informational only; threshold pre-registered independently here

### §1.4 — Out-of-scope (anti-Lakatos)

❌ NOT testing kinematic F8 mechanisms (different category)
❌ NOT inheriting factor-10 threshold from γ-7 (declared independently §3.2 below)
❌ NOT continuation of γ-3/γ-3'/γ-5/γ-7
❌ NOT citing F8 FAILs as motivation
❌ NOT citing F8_FORENSIC envelope factor-25 as predicted (envelope is information record only)
❌ NOT citing E1/E2 explorations as positive evidence
❌ NOT modifying γ-7 HALT-B or other F8 verdicts based on this cycle's outcome
✓ Independent mechanism (vacuum stress-energy)
✓ Independent observable (Λ_eff = Ω_DE in Friedmann)
✓ Own pre-registered factor-10 threshold (declared here, NOT inherited)
✓ Standalone failure modes

---

## §2 — Mandatory reading (Phase 1 start prerequisites)

Before Phase 1 sympy work:

1. **core/formalizm/dodatekE_kwantyzacja.tex**
   - eq. 97-209: TGP propagator + first iteration loop correction
   - eq. 207: Λ_eff = γ/12 (classical vacuum + IR cutoff formula)
   - eq. 353: m_sp = √γ·(ℏ_0 c_0/l_sp) ~ ℏ_0 H_0/c_0 (Phi-phonon mass)
   - eq. 365: explicit DE candidate prediction (Phi-phonon coincidence)
   - §first-iteration §170-210: quantum loop computation with IR cutoff Λ_UV^eff = √γ

2. **core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex**
   - V_M911(ψ) = -γ·ψ²·(4-3ψ)²/12 (substrate potential)
   - Critical post-2026-05-09 update: M9.1'' (4-3ψ)/ψ FORM FALSIFIED 5σ by GWTC-3
   - V_M911 derived from (4-3ψ)/ψ ansatz; per sek08c body "matematycznie pozostaje VALID"
   - **R1 #18 awareness**: sek08a §3840 gauge ambiguity (from B cycle Phase 1)

3. **core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex**
   - thm:einstein-emergence §3493: G_μν^eff = κ T_μν^(Φ) FRW exact + O(U²) general
   - §3461-3478: Friedmann emergent + Bianchi propagation
   - §3447: g_00^eff = κ T_00^(Φ) in FRW

4. **core/sek02_pole/sek02_pole.tex**
   - eq. N-preview: N[Φ] = (α/Φ_0)(∇Φ)²/Φ + β·Φ²/Φ_0² - γ·Φ³/Φ_0³
   - α=2 (from action variation); β=γ (vacuum condition)
   - K2 ontology: matter = Phi configuration

5. **meta/CALIBRATION_PROTOCOL.md**
   - §3.6.1-§3.6.14 BINDING
   - Especially §3.6.13 THIRD practical application (constants classification)
   - §3.6.11 PARTIAL_compute / PARTIAL_concept_mismatch (anticipate possible use)

6. **meta/F8_FORENSIC_2026-05-24.md**
   - §6 envelope (informational baseline, NOT prediction)
   - §6.1 critical caveat (γ = H_0²/c² is calibration, not independent)
   - §9 anti-Lakatos discipline for any future F8-related cycle

7. **meta/PRE_REGISTERED_FALSIFIERS.md**
   - PR-010 LOCKED α ∈ [-0.832, 0.832] S07 polynomial family
   - PR-017 LOCKED 2026-05-24 cycle B closure precedent

8. **External observational references** (LITERATURE_ANCHORED):
   - Planck 2018: Λ_obs = 1.10×10⁻⁵² m⁻² (Ω_Λ = 0.685 from Planck+BAO+SN)
   - w_DE constraints: w ≈ -1.03 ± 0.03 (DES+Planck+SN combined)

---

## §3 — Pre-registered falsifiers

### F-LAM-A — Sign of Λ_eff

**Hypothesis:** TGP vacuum stress-energy from V_M911(ψ_eq) + 1-loop correction → effective cosmological constant Λ_eff with definite sign.

**Pre-registered acceptance criteria:**
- **PASS**: Λ_eff > 0 (positive → drives cosmic acceleration, matches DE observation sign)
- **FAIL_SIGN**: Λ_eff < 0 (would predict cosmic deceleration, opposite of DE observation)
- **FAIL_ZERO**: Λ_eff = 0 (no vacuum energy; not a DE candidate)

**Computation route (Phase 1):**
1. Evaluate V_M911(ψ=1) classical vacuum (sek08c v3.0 form when materialized; otherwise canonical M9.1''-derived form with caveat)
2. Add 1-loop quantum correction with IR cutoff Λ_UV^eff = √γ (Appendix E §first-iteration)
3. Map to cosmological constant via sek08a G^eff_00 = κ T_00^(Φ) emergent Einstein
4. Determine sign

### F-LAM-B — Magnitude (PRIMARY OBSERVABLE)

**Hypothesis:** Λ_eff_TGP matches Λ_obs within factor 10 (own pre-registered threshold).

**Pre-registered acceptance criteria (LOCKED, immutable):**
- **PASS**: 0.1 ≤ Λ_eff_TGP / Λ_obs ≤ 10 (factor-10 threshold)
- **FAIL_HIGH**: Λ_eff_TGP / Λ_obs > 10 (TGP over-predicts; e.g., classical QFT vacuum energy at 10¹²⁰)
- **FAIL_LOW**: Λ_eff_TGP / Λ_obs < 0.1 (TGP under-predicts; e.g., the factor-25 envelope is in this regime — anticipate possible verdict)

**Computation route (Phase 1):**
1. Compute Λ_eff_TGP from V_M911 + 1-loop with explicit γ value (γ = H_0²/c² as PRIMARY input; structural consistency check)
2. Compute Λ_obs from 3·Ω_Λ·H_0²/c² (Planck 2018)
3. Compare ratio
4. **CRITICAL**: pre-registered factor-10 threshold is INDEPENDENT of γ-7's factor-10 — declared here with own justification (factor 10 = standard order-of-magnitude tolerance for natural emergence; same threshold used for many TGP cycles)

**Threshold justification** (PHASE 0 LOCK, anti-Lakatos):
- Factor 10 is appropriate for natural-emergence (no fine-tuning) prediction
- NOT inherited from γ-7 verdict structure (this cycle's mechanism is different)
- Comparable to standard cosmology "1 OOM" expectation for dimensional analysis
- DECISION RULE IMMUTABLE: do NOT loosen to factor 100 if factor 10 fails; honest verdict required

**Anticipated outcome (informational, NOT pre-registered):**
Per F8_FORENSIC §6 envelope, Λ_TGP/Λ_obs ≈ 1/25 (factor 25 OFF). This is in **FAIL_LOW** territory.

Phase 1 MUST derive ab initio, NOT inherit envelope. If derivation independently confirms factor-25 (or similar), verdict is FAIL_LOW. If Phase 3 quantum loop corrections close the gap to factor-10 ratio, verdict is PASS. Either outcome accepted honestly.

### F-LAM-C — Equation of state w_DE

**Hypothesis:** TGP phonon-vacuum mechanism predicts w_DE = -1 + O(δ) for some small δ.

**Pre-registered acceptance criteria:**
- **PASS**: |w_DE_TGP - (-1)| ≤ 0.05 (observational 2σ: DES+Planck+SN gives w = -1.03 ± 0.03)
- **FAIL**: w_DE_TGP outside observational 2σ bounds
- **PARTIAL_CONCEPT**: w_DE_TGP derivation incomplete (time-evolution requires Phase 3 quantum corrections; could declare PARTIAL_compute)

**Computation route (Phase 2):**
1. Compute time-evolution of Λ_eff(ψ_eq(t)) through cosmological epochs
2. Check if Λ_eff is approximately time-independent (→ w_DE = -1) or evolves slowly
3. Compare to observational w_DE bounds

**Reference (sek08a):** remark §10287 mentions "w_DE = -1 + O(10⁻⁹)" — Phase 2 to verify this applies to phonon-vacuum specifically.

### F-LAM-D — Quantum loop correction effectiveness

**Hypothesis:** 1-loop quantum correction (Appendix E §first-iteration extended) modifies Λ_eff by O(γ) factor, potentially closing factor-25 envelope gap to factor-10 PASS regime.

**Pre-registered acceptance criteria:**
- **PASS**: Loop correction brings Λ_eff_TGP / Λ_obs to within factor-10 (closing the gap)
- **FAIL_PRESERVES**: Loop correction preserves factor-25 (or worse) discrepancy
- **PARTIAL_compute**: Loop computation requires resources beyond cycle scope (e.g., 2-loop required for accuracy)

**Computation route (Phase 3):**
1. Implement Appendix E first-iteration loop formula (eq. 156-167)
2. Use IR cutoff Λ_UV^eff = √γ per eq. 207
3. Compute δΛ^(1) correction
4. Verify magnitude vs leading γ/12 term
5. Updated Λ_eff_TGP = γ/12 + δΛ^(1) compared to Λ_obs

---

## §4 — Methodology (BINDING)

### §4.1 — Strict cycle 1/2/7
- 0 hardcoded T_pass=True (per §3.6.1)
- Compute-then-compare for every substantive FP
- PARTIAL_compute budget: 1 max per cycle
- PARTIAL_concept_mismatch: only if formalism gap discovered honestly (declare R1)
- DEC budget: 3 max per cycle

### §4.2 — Anti-Lakatos discipline
- NO modification of pre-registered falsifiers post Phase 1 start
- NO re-framing of FAIL_LOW as "marginal consistency"
- NO post-hoc threshold loosening to factor-100 if factor-10 fails
- NO citing F8 cycles as motivation post-hoc
- NO using envelope factor-25 as predicted

### §4.3 — Falsifier independence
- Each falsifier has standalone pass/fail
- F-LAM-A (sign) failure does NOT auto-fail B/C/D
- F-LAM-B (magnitude) failure does NOT auto-fail A/C/D
- F-LAM-D (loop correction) failure does NOT auto-fail others
- Cycle aggregate verdict synthesizes individual outcomes honestly

### §4.4 — R1 #18 awareness (gauge ambiguity)
- Any sek08a §3840-style Δg formulas used must declare gauge explicitly
- Prefer gauge-invariant observables (Λ_eff is gauge-invariant in cosmological FRW context)
- If gauge ambiguity surfaces, declare R1 flag explicitly

---

## §5 — LEGITIMATE inheritance (LOCKED TGP infrastructure)

| Source | Element | Status | Use in this cycle |
|--------|---------|--------|-------------------|
| Appendix E §kwantyzacja | Eq. 97-104 KG propagator | concept paper postulate | foundation for vacuum loop |
| Appendix E eq. 207 | Λ_eff = γ/12 formula | concept paper derivation (with IR cutoff caveat) | starting equation F-LAM-B |
| Appendix E eq. 353 | m_sp ~ ℏH_0/c² | concept paper postulate | γ-H_0 mapping (structural consistency, NOT prediction) |
| Appendix E eq. 365 | Phi-phonon DE candidate prediction | concept paper postulate | primary motivation |
| Appendix E §first-iteration | δΛ^(1) loop correction sketch | concept paper sketch | F-LAM-D extension |
| sek08c V_M911 | -γ·ψ²·(4-3ψ)²/12 | concept paper (M9.1''-derived; M9.1'' falsified but V_M911 form may survive) | F-LAM-A sign derivation |
| sek08a thm:einstein-emergence | G^eff = κ T^(Φ) FRW exact | concept paper LOCKED | Friedmann mapping |
| sek02_pole eq. N-preview | α=2, β=γ vacuum | concept paper LOCKED | non-linear corrections (if needed Phase 3) |
| γ-3' Phase 3 | ℓ_P, E_P, ω_P calibration | LOCKED PASS | scale fixing |
| γ-5 Phase 3 | G_eff = c³·ℓ_P²/ℏ | LOCKED PASS | gravitational scale |
| γ-7 Phase 1 | q = √(4πG)·m | LOCKED PASS | source coupling (if needed) |

**Critical:** ALL inheritance is from **LOCKED PASS** items or **concept paper postulates**. NO inheritance from F8 FAIL verdicts.

---

## §6 — FORBIDDEN moves (anti-Lakatos)

| Move | Why forbidden |
|------|---------------|
| Cite γ-3/3'/5/7 F8 FAILs as motivation | Post-hoc rescue logic |
| Cite F8_FORENSIC envelope factor-25 as "prediction" | Envelope is information record only |
| Cite E1/E2 explorations as positive evidence | Explorations are scouting, not pre-registered tests |
| Loosen factor-10 threshold to factor-100 if magnitude FAIL_LOW | Threshold inflation |
| Re-frame F-LAM-B FAIL_LOW as "marginal PASS" | Verdict tampering |
| Add new falsifiers post Phase 1 start | Anti-Lakatos LOCK violation |
| Frame as "γ-8" or continuation of cosmology cycles | Naming convention violation per F8_FORENSIC §9 |
| Modify γ-7 HALT-B or F8 verdicts based on Λ_eff result | Verdict tampering across cycles |
| Use γ derived from H_0 then claim "predicts H_0" | Circular reasoning; must disclose as structural consistency |
| Cite cycle D outcome as already-resolved when D is QUEUED | Don't claim D's verdict before D executes |

---

## §7 — Constants classification (§3.6.13 THIRD practical application)

Constants entering Phase 1-3 derivations:

| # | Constant | Category | Source | Notes |
|---|----------|----------|--------|-------|
| 1 | γ | **CALIBRATED** | Appendix E eq. 353 (m_sp ~ ℏH_0/c²) | Currently calibrated to H_0; D-cycle prerequisite for independent derivation |
| 2 | H_0 | **EXTERNAL_OBSERVATION** | Planck 2018: 67.7 km/s/Mpc | Input to γ-calibration |
| 3 | c | **DERIVED** | γ-5 Phase 3 c = ℓ_P·ω_P LOCKED | inherited |
| 4 | ℓ_P | **CALIBRATED** | γ-3' Phase 3 LOCKED | inherited |
| 5 | ℏ | **INPUT** | external (γ-6 deferred for derivation) | input |
| 6 | G | **DERIVED** | γ-5 Phase 3 G_eff = c³·ℓ_P²/ℏ | inherited |
| 7 | Φ_0 | **CALIBRATED** | sek02 footnote g_0 = 0.87 | dimensionless |
| 8 | α (S07 polynomial) | **EXTERNAL_BOUND** | PR-010 LOCKED ∈ [-0.832, 0.832] | not core to this cycle |
| 9 | Λ_UV^eff (IR cutoff) | **DERIVED** | Appendix E eq. 207 Λ_UV^eff = √γ | inherited |
| 10 | ω_P | **DERIVED** | γ-3' Phase 3 | inherited |
| 11 | E_P | **DERIVED** | γ-3' Phase 3 E_P = ℏω_P | inherited |
| 12 | Λ_obs (DE) | **EXTERNAL_OBSERVATION** | Planck 2018 Ω_Λ = 0.685 → Λ_obs ≈ 1.10×10⁻⁵² m⁻² | observational target |
| 13 | w_DE | **EXTERNAL_OBSERVATION** | DES+Planck+SN w ≈ -1.03 ± 0.03 | F-LAM-C target |
| 14 | β (sek02 N[Φ]) | **DERIVED** | sek02 vacuum condition β = γ | inherited |
| 15 | M_sun | **EXTERNAL_CONSTANT** | 1.989×10³⁰ kg | not directly used here |
| 16 | β_q (S07 quadratic) | **EXTERNAL_BOUND** | PR-010 LOCKED ∈ [-0.4, 0.4] | not core to this cycle |

**Total: 16 constants classified.** No new fundamental constants introduced. All inputs are either:
1. TGP-locked (from prior cycles)
2. Concept-paper postulate (Appendix E, sek08c)
3. External observational (LITERATURE_ANCHORED)
4. Currently calibrated (γ — pending D cycle for derivation)

---

## §8 — Risk register

| ID | Risk | Severity | Mitigation |
|----|------|----------|------------|
| R-LAM-1 | γ calibration circularity: γ ← H_0; Λ_TGP ← γ; "predicts" H_0 — circular | HIGH | Disclose as structural consistency in Phase 0 §1.3; honest classification §6 forbidden move |
| R-LAM-2 | M9.1'' (4-3ψ)/ψ FALSIFIED 5σ; V_M911 derived from M9.1''; may need replacement | HIGH | sek08c body still mathematically VALID per pre-falsification; use V_M911 form with caveat; if S07 alternative materialized, switch |
| R-LAM-3 | Sign of Λ_eff convention-dependent (V_M911 vs +V or -V) | MEDIUM | Phase 1 explicit convention check via sek08a G^eff_00 = κ T_00^(Φ) Friedmann mapping |
| R-LAM-4 | Quantum loop UV divergences regulated by IR cutoff Λ_UV^eff = √γ — non-standard QFT regularization | MEDIUM | Appendix E §first-iteration argues this is consistent (IR cutoff motivated by m_sp scale); Phase 3 to verify finite result |
| R-LAM-5 | R1 #18 gauge ambiguity (sek08a §3840) — does Λ definition depend on gauge? | LOW | Λ in FRW context is gauge-invariant (defined via Friedmann constraint); declare explicitly Phase 1 |
| R-LAM-6 | w_DE prediction may depend on detailed time evolution of ψ_eq(t) | MEDIUM | Phase 2 explicit time-evolution; if PARTIAL_compute needed (Phase 3 chain inflation epoch tracking), declare honestly |
| R-LAM-7 | Cycle anticipated outcome: factor-25 FAIL_LOW per envelope; honest disclosure required | LOW | Anticipated outcome NOT pre-determined; Phase 1 ab initio computation may differ; do NOT cite envelope as predicted |
| R-LAM-8 | D cycle (independent γ) is QUEUED — this cycle's "prediction" status depends on D outcome | MEDIUM | Pre-register PRIMARY = structural consistency; SECONDARY = upgrade if D delivers; don't claim D's verdict |
| R-LAM-9 | PARTIAL_compute likely needed for full 1-loop computation | LOW | 1/1 budget reserved; honest declaration |
| R-LAM-10 | Concept paper Appendix E factor 1/12 may have different convention (γ/12 vs γ/24 vs 2γ etc.) — sympy verification needed | MEDIUM | Phase 1 explicit V_M911 → Λ_eff derivation from action principle |

---

## §9 — Phase plan

| Phase | FP target | Duration | Output |
|-------|-----------|----------|--------|
| 0 | scaffold + falsifiers + reading | this | Phase0_balance.md + LOCK |
| 1 | F-LAM-A sign + F-LAM-B magnitude (leading order) | 1-2 sesji | Phase1_sympy.py + Phase1_derivation.md |
| 2 | F-LAM-C equation of state w_DE | 1 sesja | Phase2_sympy.py + Phase2_results.md |
| 3 | F-LAM-D quantum loop correction extension | 1-2 sesji | Phase3_sympy.py + Phase3_results.md |
| FINAL | aggregate verdict + claim_status + PR-### LOCK entry | 0.5 sesja | Phase_FINAL_close.md + PRE_REGISTERED_FALSIFIERS.md PR entry |

**Total estimate:** 3-5 sesji.

**Decision points between phases:**
- Post Phase 1: if F-LAM-B confirms factor-25 FAIL_LOW direction, user may decide to:
  - Continue to Phase 2-3 (verify F-LAM-C/D honestly)
  - HALT (declare cycle outcome at Phase 1; document loop correction work as future scope)
  - Pivot to D cycle (if independent γ might change interpretation)
- Post Phase 2: w_DE result informs whether cycle is "consistent with DE-like behavior" regardless of magnitude

---

## §10 — Decision budget

| Budget | Cap | Used | Remaining |
|--------|-----|------|-----------|
| DEC (substantive decision) | 3 | 0 | 3 |
| PARTIAL_compute | 1 | 0 | 1 |
| PARTIAL_concept_mismatch | unrestricted (declare R1) | 0 | unlimited |
| Hardcoded T_pass=True | 0 | 0 | 0 |

---

## §11 — Anti-Lakatos verification (Phase 0 check)

| Item | Status |
|------|--------|
| F8 FAILs cited as motivation? | NO ✓ |
| F8_FORENSIC envelope cited as predicted? | NO (cited only as informational baseline §6) ✓ |
| E1/E2 explorations cited as predictions? | NO ✓ |
| Threshold inherited from γ-7? | NO (independent factor-10 declared §3.2 with own justification) ✓ |
| Cycle named γ-8? | NO (op-LAM-vacuum-substrate; reflects mechanism) ✓ |
| Pre-registered falsifiers BEFORE derivation? | YES (this document) ✓ |
| Standalone fail modes declared? | YES (F-LAM-A/B/C/D each independent) ✓ |
| Independent mechanism from F8 FAIL cycles? | YES (vacuum stress-energy vs kinematic/clumping) ✓ |
| Structural consistency vs prediction status disclosed? | YES (§1.3 explicit) ✓ |
| Forbidden moves register comprehensive? | YES (§6 explicit) ✓ |
| Anticipated FAIL_LOW outcome disclosed honestly? | YES (§3.2 anticipated outcome note) ✓ |

**Anti-Lakatos status:** COMPLIANT ✓

---

## §12 — Phase 0 LOCK status

**LOCKED 2026-05-25** by user authorization "działaj Phase 1".

- Scope (§1): CONFIRMED
- Falsifiers (§3): F-LAM-A/B/C/D LOCKED (immutable post-lock)
- Threshold §3.2: factor-10 LOCKED (independent declaration, NOT inherited from γ-7)
- Methodology (§4): cycle 1/2/7 strict; DEC budget 3; PARTIAL_compute 1
- Inheritance (§5): LEGITIMATE concept-paper + LOCKED prior cycles only
- Forbidden moves (§6): registered; ANTI-LAKATOS COMPLIANT
- Risk register (§8): 10 items registered
- Phase plan (§9): 5 phases (0/1/2/3/FINAL)
- PR-018 reserved for Phase FINAL entry

**Phase 1 authorized. Execute F-LAM-A (sign) + F-LAM-B (magnitude leading order).**

---

## §13 — Status summary

| Field | Value |
|-------|-------|
| Phase 0 status | LOCKED 2026-05-25 |
| Authorization | "działaj Phase 1" 2026-05-25 ✓ |
| F8 independence | DECLARED (different mechanism: vacuum stress-energy) ✓ |
| Structural consistency vs independent prediction | DECLARED §1.3 (current: structural; upgrade conditional on D) ✓ |
| Anti-Lakatos | COMPLIANT ✓ |
| Inheritance | LEGITIMATE only ✓ |
| Phase 1 prerequisites | Appendix E + sek08c + sek08a + CALIBRATION_PROTOCOL §3.6.13 |
| Expected duration | 3-5 sesji |
| Anticipated outcome (informational) | factor-25 FAIL_LOW envelope; Phase 1 must derive ab initio |
| Related cycles | B CLOSED-RESOLVED (PR-017); D QUEUED (prerequisite for independent prediction); C DEFERRED |
| Threshold | factor-10 (declared independently, NOT inherited) |
