---
title: "Phase 0 — balance sheet + pre-registration (op-G-substrate-derivation)"
type: phase_balance
status: PHASE0_LOCKED
phase: 0
cycle: op-G-substrate-derivation-2026-05-24
parent_motivation: "Appendix E eq. 352-355: γ currently calibrated to Hubble scale via l_sp = 1/√γ ≈ R_H (m_sp ~ ℏH_0/c). Test whether γ can be DERIVED from {ℓ_P, c, ℏ, V_M911 Lagrangian parameters, Φ_0} WITHOUT H_0 input. Outcome determines whether cycle A (LAM) Λ_eff = γ/12 can be promoted from STRUCTURAL_CONSISTENCY (C+) to INDEPENDENT_PREDICTION."
created_date: 2026-05-24
locked_date: 2026-06-01
locked_by: "User 2026-06-01: 'Ok działaj z cyklem D' → cycle D Phase 0 activation authorized"
PR_reserved: "PR-019 reserved for Phase FINAL closure"
methodology_binding: "CALIBRATION_PROTOCOL.md §3.6 BINDING (incl. §3.6.13 FOURTH practical application)"
anti_lakatos_lock: PRESERVED
independent_of_F8_cycles: YES (different mechanism category: foundational scale derivation, NOT F8 mechanism test)
tests_F8_indirectly: NO (this cycle does NOT test any F8 mechanism; F8 implication is downstream consequence of D's result, not D's primary objective)
prerequisite_for: "op-LAM-vacuum-substrate upgrade beyond STRUCTURAL_PARTIAL (C+)"
expected_duration: "2-3 sesji"
---

# Phase 0 — op-G-substrate-derivation pre-registration

## §1 — Scope declaration

### §1.1 — Primary objective

Test whether the Phi-substrate curvature parameter γ — appearing in TGP as

- effective phonon mass: m_sp² = γ (Appendix E eq. 104, 161, 324)
- vacuum cosmological constant: Λ_eff = γ/12 (Appendix E eq. 207)
- Yukawa range: l_sp = 1/√γ (Appendix E eq. 355)

— can be **derived** from a TGP-internal combination of {ℓ_P, c_0, ℏ_0, V_M911 Lagrangian parameters, Φ_0} **without** using H_0 (or any cosmological observable) as input.

The cycle has a **binary structural verdict** (derivation exists / does not exist), independent of any numerical magnitude question. If derivation exists, F-G-B/C/D test its numerical adequacy.

### §1.2 — Mechanism category (this cycle ONLY)

**This cycle tests:** existence of a closed-form γ formula expressible in TGP-internal fundamentals, plus numerical comparison of derived γ against the currently-calibrated value γ_cal ≡ H_0²/c_0².

**OTHER categories OUT OF SCOPE for this cycle:**
- F8 (DE acceleration) mechanism tests — kinematic (γ-3, γ-3′), c-variation (γ-5), clumping (γ-7), vacuum stress-energy (op-LAM-vacuum-substrate) — all separate cycles, **not cited as motivation**
- NS surface O(U³) drift (op-PSR-orbital-drift, PR-017 CLOSED-RESOLVED B+) — different observable category
- Emergent time formalism (op-EMT-emergent-time DEFERRED multi-cycle program)
- S07 alternative f(ψ) derivation (separate gravity-sector cycle)
- Λ_eff magnitude prediction — that is op-LAM-vacuum-substrate scope (PR-018 LOCKED STRUCTURAL_PARTIAL)

### §1.3 — Critical scope distinction: foundation vs. consequence

**This cycle, in current form, tests FOUNDATIONAL DERIVATION of γ:**
- Primary observable is **existence** of derivation formula γ = F(ℓ_P, c_0, ℏ_0, …) with NO H_0 input
- Secondary observable is **numerical match** γ_derived vs γ_cal = H_0²/c_0² ≈ 5.37×10⁻⁵³ m⁻²
- Cycle delivers a verdict on γ's epistemological status, NOT on F8

**Downstream consequence (informational, NOT pre-registered as this cycle's verdict):**
- If F-G-A PASS + F-G-B PASS → op-LAM-vacuum-substrate (cycle A) Λ_eff = γ/12 upgrades from STRUCTURAL_CONSISTENCY (C+) to INDEPENDENT_PREDICTION (status reclassification pending separate cycle re-assessment)
- If F-G-A FAIL (no derivation found) → γ remains OBSERVATIONAL_ANCHOR (γ); cycle A stays C+ permanently; concept paper Appendix E formalism unchanged

**This cycle declares (explicit):**
- PRIMARY scope: γ derivability test (independent of cycle A's existing C+ status)
- SECONDARY scope: numerical match if derivation exists
- TERTIARY scope: implications for H_0 (predicted vs observed) if derivation exists
- ANTI-LAKATOS: **NO motivation by F8 FAILs**. The motivation is a gap in concept paper Appendix E (eq. 353 reads "γ ~ H_0²/c_0²" — this is a calibration, not a derivation, and Appendix E itself flags it as "coincidence problem", hyp:coincidence).

### §1.4 — Out-of-scope (anti-Lakatos)

❌ NOT testing F8 mechanisms (different category)
❌ NOT citing F8 FAILs (γ-3/γ-3′/γ-5/γ-7) as motivation
❌ NOT citing F8_FORENSIC envelope factor-25 as predicted outcome
❌ NOT framed as "γ-8" or continuation of cosmology cycles
❌ NOT modifying cycle A (PR-018) verdict based on this cycle's outcome (cycle A LOCKED; its upgrade requires a separate reassessment cycle, not this one)
❌ NOT modifying γ-7 HALT-B
❌ NOT introducing new fundamental constants (cycle tests EXISTING γ derivability)
❌ NOT post-hoc selection from multi-candidate route set without pre-declared selection criterion (per §3.6 anti-pattern κ.1 multi-path warning)
✓ Independent mechanism category (foundational scale derivation)
✓ Standalone failure mode (γ may NOT be derivable; that's a valid PASS for this cycle's audit)
✓ Pre-declared route enumeration (A-E §3 below) BEFORE Phase 1 execution
✓ Binary structural verdict on F-G-A (exists / does not exist), not a numerical threshold

---

## §2 — Mandatory reading (Phase 1 start prerequisites)

Before Phase 1 sympy work:

1. **core/formalizm/dodatekE_kwantyzacja.tex**
   - eq. 97-104: TGP propagator with m_sp² = γ (foundational definition)
   - eq. 153-167: 1-loop ⟨Φ⟩^(1) with Λ_UV = ℓ_P⁻¹ (Planck cutoff scheme)
   - eq. 193-211: IR cutoff scheme Λ_UV^eff = √γ + O15 open problem note
   - eq. 207: Λ_eff = γ/12 (cosmological constant formula)
   - eq. 224-270: thm:natural-cutoff — substrate lattice a_sub = ℓ_P (cutoff is ABSOLUTE, Φ-independent)
   - **eq. 304-309**: "Podstawiając m_sp² = γ ~ H_0²/c_0²" — **THE CALIBRATION TO BE CHALLENGED**
   - eq. 339-367: phonon dispersion + l_sp = 1/√γ ≈ R_H (Hubble radius identification)
   - eq. 352-355: m_sp = √γ · ℏ_0 c_0 / l_sp ~ ℏ_0 H_0 / c_0 (the postulated relation)
   - **rem:naturalness §307-332**: Appendix E's own framing of γ as "self-consistent loop, not fine-tuning" — acknowledges status of derivation gap
   - **prob:kwantyzacja §405-430 (O15)**: open problem statement; this cycle directly addresses O15 in its "γ from Lagrangian" facet

2. **core/sek02_pole/sek02_pole.tex**
   - N[Φ] = (α/Φ_0)(∇Φ)²/Φ + β·Φ²/Φ_0² − γ·Φ³/Φ_0³
   - α = 2 derived (Phase 2 universal mass formula)
   - β = γ (vacuum condition; β and γ identified)
   - **Critical**: the γ in N[Φ] cubic coupling is the same γ as in m_sp² = V''(1)
   - Φ_0 ≈ 25 (dimensionless per g_0^e ≈ 0.87 calibration, sek04)

3. **core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex**
   - V_M911(ψ) = −γ·ψ²·(4−3ψ)²/12 (substrate potential; γ as overall coefficient)
   - Critical post-2026-05-09: M9.1″ (4−3ψ)/ψ FORM FALSIFIED 5σ; but V coefficient structure (γ as overall scale) survives in any (4−3ψ)/ψ-class potential
   - R1 #18: gauge ambiguity (declare gauge for any Δg_xx-style formula)

4. **core/sek04_stale/**
   - Φ_0 calibration via g_0^e ≈ 0.87 (lepton paper anchor)
   - Status of Φ_0: dimensionless (per Φ_0 = ⟨ŝ²⟩/(reference scale))

5. **meta/CALIBRATION_PROTOCOL.md**
   - §3.6.1-§3.6.14 BINDING
   - **Especially §3.6.13 FOURTH practical application (constants classification)** — γ classification is THE central question of this cycle
   - §3.6.11 PARTIAL_compute / PARTIAL_concept_mismatch (anticipate possible use)
   - §3.6.5 multi-path selection warning (κ.1 anti-pattern) — applies to routes A-E §3 below

6. **meta/F8_FORENSIC_2026-05-24.md**
   - §6.1 critical caveat: "γ = H_0²/c² is calibration, not derivation" — this cycle's direct mandate
   - §6 informational only (NOT predicted outcome basis)
   - §9 anti-Lakatos discipline for any γ-related cycle

7. **meta/PRE_REGISTERED_FALSIFIERS.md**
   - PR-018 LOCKED 2026-05-25: cycle A STRUCTURAL_PARTIAL C+ (target this cycle's upgrade-enabling result, but DO NOT cite as motivation)

8. **research/op-LAM-vacuum-substrate-2026-05-24/Phase_FINAL_close.md**
   - cycle A closure (read for context on Λ_eff structure, NOT as motivation)
   - **Note**: cycle A's FAIL_LOW factor 21.4 result is INFORMATIONAL ONLY for this cycle; D's outcome must be derived independently

9. **External**:
   - Standard QFT dimensional analysis (mass/length scales, Wilsonian RG)
   - Coincidence problem literature (e.g., Weinberg 1989, Padmanabhan 2003) — for honest context, not as evidence

---

## §3 — Pre-registered derivation routes (LOCKED, immutable post-LOCK)

Five candidate routes for γ derivation, pre-enumerated to prevent post-hoc selection (§3.6.5 anti-pattern κ.1 mitigation):

### Route A — Planck-scale natural (UV)

**Hypothesis:** γ as natural Lagrangian coupling at UV cutoff ℓ_P
**Formula attempt:** γ_A = c₁ · ℓ_P⁻² for some dimensionless coefficient c₁ = O(1)
**Anticipated outcome (informational):** γ_A ≈ ℓ_P⁻² ≈ 10⁷⁰ m⁻² vs γ_cal ≈ 10⁻⁵² m⁻² → factor 10¹²² off (the classical cosmological constant problem). FAIL_HIGH.
**Status:** route TESTABLE; expected FAIL; documents that pure UV-natural γ is incompatible with TGP framework (consistent with Appendix E rem:naturalness self-acknowledgment).

### Route B — Dimensional Φ_0 + ℓ_P combinations

**Hypothesis:** γ derivable from {ℓ_P, c_0, ℏ_0, Φ_0} via dimensional analysis with Φ_0 entering as suppressing factor
**Formula attempt:** γ_B = ℓ_P⁻² · f(Φ_0) where f is a dimensionless function (e.g., Φ_0⁻ⁿ for some n)
**Anticipated outcome (informational):** Φ_0 ≈ 25 cannot bridge 122 OOM gap via any power law; would require Φ_0⁻⁸⁵ ≈ 10⁻¹²², physically unmotivated. Likely FAIL_HIGH or FAIL_UNMOTIVATED.
**Status:** route TESTABLE; serves to rule out simple dimensional cures.

### Route C — RG dimensional transmutation (IR fixed point with anomalous dimension)

**Hypothesis:** γ at IR (cosmological scale) differs from γ at UV (Planck) by anomalous dimension via Wilsonian RG flow in the Φ⁴-class theory (sek02 N[Φ]); dimensional transmutation à la QCD Λ_QCD = M_Pl·exp(−c/g²)
**Formula attempt:** γ_C = ℓ_P⁻² · exp(−ξ · Φ_0 / g_eff) or similar exponential suppression from RG flow
**Anticipated outcome (informational):** UNKNOWN. This is the genuinely open route. For Φ_0 ≈ 25, exp(−Φ_0) ≈ 10⁻¹¹ — still 100 OOM short of 10⁻¹²² gap. Would need exp(−Φ_0³) or similar. Compatibility with TGP RG flow needs ab initio Phase 1 derivation.
**Status:** route TESTABLE; this is the MOST LIKELY route to deliver something nontrivial; may yield PARTIAL_concept_mismatch if formal Wilson-RG machinery for Φ⁴-class TGP is not yet in concept paper.

### Route D — Geometric / topological self-consistency

**Hypothesis:** γ determined by self-consistent matching l_sp = 1/√γ to some TGP-internal geometric scale (e.g., topological soliton size, causal horizon, substrate periodicity)
**Formula attempt:** γ_D from condition l_sp = L_internal where L_internal is computed from TGP-internal action
**CRITICAL CIRCULARITY RISK:** If L_internal = R_H by construction, route is circular (returns to calibration). Must declare and audit.
**Anticipated outcome:** Either FAIL (no non-circular L_internal identifiable) or PARTIAL_concept_mismatch (concept paper lacks the relevant geometric definition)
**Status:** route TESTABLE WITH MANDATORY CIRCULARITY CHECK; if any geometric scale invoked, must trace its origin and verify it does not implicitly use H_0.

### Route E — Action-principle internal relations

**Hypothesis:** γ is fixed by structural relations within V_M911 and N[Φ] action coefficients (e.g., γ = α · β / Φ_0² for some structural identity)
**Formula attempt:** γ_E = f(α, β, …) where α and β are the other N[Φ] coefficients
**Known constraint:** sek02 vacuum condition gives β = γ (not γ = f(β)); α = 2 fixed separately. So γ enters as overall scale, not derivable from action-internal ratios alone.
**Anticipated outcome:** FAIL (γ is the overall potential scale; no internal action ratio fixes it). Documents that V_M911 + N[Φ] form does NOT determine γ.
**Status:** route TESTABLE; expected FAIL; serves as structural completeness check.

### Multi-route selection rule (LOCKED, anti-Lakatos)

**If multiple routes pass F-G-A (deliver a closed-form γ formula):**
- Pre-declared selection criterion (LOCKED here): **prefer route requiring fewest external inputs**, ranked
  1. Routes using only {ℓ_P, c_0, ℏ_0} → highest preference (pure Planck-fundamental)
  2. Routes using {ℓ_P, c_0, ℏ_0, Φ_0} → next preference (Φ_0 is TGP-internal calibration)
  3. Routes using {ℓ_P, c_0, ℏ_0, Φ_0, additional Lagrangian coefficients} → lowest preference
- If two routes tie under above ranking, **declare R1 flag for ambiguity**; do NOT cherry-pick by closer-to-observation numerical match (κ.1 anti-pattern explicit ban)

---

## §4 — Pre-registered falsifiers

### F-G-A — Existence of derivation (PRIMARY, structural)

**Hypothesis:** γ admits a closed-form symbolic expression in {ℓ_P, c_0, ℏ_0} or {ℓ_P, c_0, ℏ_0, Φ_0} or {ℓ_P, c_0, ℏ_0, Φ_0, V_M911 coefficients} that does NOT depend on H_0 or any cosmological observable.

**Pre-registered acceptance criteria:**
- **PASS_PURE**: γ = F(ℓ_P, c_0, ℏ_0) — pure Planck-fundamental, no further inputs
- **PASS_WITH_PHI0**: γ = F(ℓ_P, c_0, ℏ_0, Φ_0)
- **PASS_WITH_LAGRANGIAN**: γ = F(ℓ_P, c_0, ℏ_0, Φ_0, {α, β, additional dimensionless coefficients})
- **FAIL_NO_DERIVATION**: no closed-form expression identifiable across routes A-E (γ remains OBSERVATIONAL_ANCHOR; cycle CLOSED-RESOLVED with HONEST_NEGATIVE verdict)
- **FAIL_CIRCULAR**: route D returns L_internal = R_H by construction (or any route silently uses H_0 in its formula)

**Computation route (Phase 1):**
- Execute routes A, B, C, D, E in sympy
- For each: write candidate formula; verify dimensional consistency; verify no H_0 leakage (audit)
- Aggregate to verdict per acceptance criteria above
- Apply multi-route selection rule §3 if multiple PASS

**HONEST_NEGATIVE clause:** A FAIL_NO_DERIVATION verdict is a **valid cycle PASS for audit**. It honestly establishes that γ is a fundamental free parameter in current TGP scope, that concept paper Appendix E's calibration m_sp ~ ℏH_0/c is justified as an empirical input, and that cycle A's STRUCTURAL_CONSISTENCY (C+) is the correct epistemological classification (not factor-25 FAIL_LOW per se).

### F-G-B — Numerical match (SECONDARY, conditional on F-G-A PASS)

**Hypothesis:** If γ derived (F-G-A PASS), the numerical value matches γ_cal ≡ H_0²/c_0² within order-of-magnitude tolerance.

**Pre-registered acceptance criteria (LOCKED):**
- **PASS**: 0.1 ≤ γ_derived / γ_cal ≤ 10 (factor-10 threshold)
- **FAIL_HIGH**: γ_derived / γ_cal > 10 (over-prediction; e.g., route A Planck UV ≈ 10¹²² OOM)
- **FAIL_LOW**: γ_derived / γ_cal < 0.1 (under-prediction)
- **NOT_APPLICABLE**: triggered if F-G-A returns FAIL_NO_DERIVATION

**Threshold justification** (PHASE 0 LOCK, anti-Lakatos):
- Factor 10 = standard OOM tolerance for natural-emergence prediction
- NOT inherited from γ-7 verdict structure (different cycle, different mechanism category)
- NOT inherited from cycle A factor-10 (cycle A's was for Λ_eff magnitude; this is for γ value)
- DECISION RULE IMMUTABLE: do NOT loosen to factor 100 if factor 10 fails; honest verdict required

**Computation route (Phase 2):**
- Numerical substitution: ℓ_P = 1.616×10⁻³⁵ m, c_0 = 2.998×10⁸ m/s, ℏ_0 = 1.054×10⁻³⁴ J·s, Φ_0 ≈ 25 (sek04 calibration), additional Lagrangian coefficients if needed
- Compute γ_derived from F-G-A formula
- Compute γ_cal = (H_0/c_0)² with H_0 = 67.7 km/s/Mpc (Planck 2018) ≈ 5.37×10⁻⁵³ m⁻²
- Compare ratio

### F-G-C — Consistency with Appendix E eq. 353 (CROSS-CHECK)

**Hypothesis:** If F-G-A PASS, the derived γ remains consistent with Appendix E eq. 353 m_sp = √γ · ℏ_0 c_0 / l_sp at the appropriate scale.

**Pre-registered acceptance criteria:**
- **PASS**: m_sp_derived ≈ m_sp_AppE within factor 10 (i.e., factor 10 in mass, factor 100 in γ; weaker than F-G-B)
- **FAIL**: order-of-magnitude inconsistency with Appendix E's own scaling
- **NOT_APPLICABLE**: F-G-A FAIL

**Purpose:** consistency check that derived γ embeds correctly in concept paper formalism (does not invalidate Appendix E by orthogonal numerical conflict).

### F-G-D — Implication for H_0 (TERTIARY, true-prediction status)

**Hypothesis:** If F-G-A PASS AND F-G-B PASS, then γ_derived inverts to predict H_0:
H_0_predicted = c_0 · √(γ_derived) (per Appendix E l_sp = 1/√γ ≈ R_H identification — but now γ is INPUT and H_0 is OUTPUT, reversed direction)

**Pre-registered acceptance criteria:**
- **PASS_HUBBLE**: 0.1 ≤ H_0_predicted / H_0_obs ≤ 10 (factor-10 threshold consistent with F-G-B)
- **FAIL**: H_0_predicted outside factor-10 of H_0_obs
- **NOT_APPLICABLE**: F-G-A FAIL or F-G-B FAIL

**Significance:** PASS_HUBBLE would mean TGP predicts H_0 from {ℓ_P, c_0, ℏ_0, Φ_0} alone — a true cosmological prediction. This is the maximum possible payoff of the cycle.

**Honest disclosure (anti-Lakatos):** F-G-D PASS does NOT prove TGP; it shows γ has a TGP-native derivation that is numerically consistent with observed H_0. Tension with H_0 itself (Hubble tension: 67 km/s/Mpc vs 73 km/s/Mpc) is OUT OF SCOPE for this cycle (factor-10 tolerance covers both).

---

## §5 — Methodology (BINDING)

### §5.1 — Strict cycle 1/2/7
- 0 hardcoded T_pass=True (per §3.6.1)
- Compute-then-compare for every substantive FP
- PARTIAL_compute budget: 1 max per cycle
- PARTIAL_concept_mismatch: only if formalism gap discovered honestly (declare R1)
- DEC budget: 3 max per cycle

### §5.2 — Anti-Lakatos discipline
- NO modification of pre-registered falsifiers post Phase 1 start
- NO re-framing of FAIL_HIGH (route A Planck UV) as "not a true failure"
- NO post-hoc selection between routes A-E by closer-to-H_0 numerical match (§3.6.5 κ.1 anti-pattern; selection rule pre-declared §3)
- NO citing F8 cycles as motivation post-hoc
- NO using factor-25 envelope (F8_FORENSIC §6) as predicted result

### §5.3 — Falsifier independence
- Each falsifier has standalone pass/fail
- F-G-A (structural) failure → F-G-B/C/D NOT_APPLICABLE (cycle still produces HONEST_NEGATIVE verdict; not auto-FAIL)
- F-G-B (numerical) failure → F-G-D NOT_APPLICABLE; F-G-C still tested
- Aggregate verdict synthesizes individual outcomes honestly

### §5.4 — R1 #18 awareness (gauge ambiguity)
- Not directly relevant (γ is gauge-invariant scalar Lagrangian coupling)
- However: any Phase 3 cosmological mapping (F-G-D) must declare cosmological frame explicitly (FRW comoving vs proper); R1 #18 lesson — declare gauge / frame BEFORE computation

### §5.5 — H_0 circularity audit (CRITICAL)
- Every Phase 1 candidate formula MUST be audited explicitly: does it use H_0 (directly or indirectly through R_H, t_H = 1/H_0, ρ_crit, etc.)?
- Any implicit H_0 usage → FAIL_CIRCULAR for that route
- Audit is performed at sympy-symbolic level (substitute H_0 → 0; check if formula degenerates)

---

## §6 — LEGITIMATE inheritance (LOCKED TGP infrastructure)

| Source | Element | Status | Use in this cycle |
|--------|---------|--------|-------------------|
| Appendix E §kwantyzacja | KG propagator with m_sp² = γ | concept paper LOCKED | foundational definition |
| Appendix E eq. 207 | Λ_eff = γ/12 | concept paper LOCKED | structural relation (γ as overall coupling) |
| Appendix E eq. 304-309 | "γ ~ H_0²/c_0²" calibration | concept paper postulate | **TARGET TO CHALLENGE** (not inherited as fact) |
| Appendix E eq. 352-355 | m_sp = √γ · ℏ_0 c_0 / l_sp ~ ℏ_0 H_0 / c_0 | concept paper postulate | reference relation for F-G-C consistency check |
| Appendix E thm:natural-cutoff | a_sub = ℓ_P; cutoff absolute | concept paper LOCKED | route A foundation |
| sek02 N[Φ] α = 2, β = γ | action coefficients | LOCKED | route E foundation (rules out internal action ratios) |
| sek04 Φ_0 ≈ 25 | dimensionless | LOCKED | routes B, C, D, E input |
| ℓ_P, c_0, ℏ_0 numerical values | CODATA / γ-3′ Phase 3 | LOCKED | numerical substitution |
| γ-5 Phase 3 G_eff = c_0³ · ℓ_P² / ℏ_0 | LOCKED | inheritable scale relation if route C invokes G |
| γ-7 Phase 1 q = √(4πG)·m | LOCKED | available for source-term arguments if needed (likely not used in this cycle) |
| cycle A PR-018 Λ_eff = γ/12 final form | LOCKED | reference for F-G-D upgrade implication (NOT for motivation) |

**Critical:** ALL inheritance from **LOCKED PASS** items or **concept paper postulates**. The specific Appendix E line "γ ~ H_0²/c_0²" is inherited as TARGET TO TEST, not as fact.

---

## §7 — FORBIDDEN moves (anti-Lakatos)

| Move | Why forbidden |
|------|---------------|
| Cite γ-3/3′/5/7 F8 FAILs as motivation | Post-hoc F8 rescue logic |
| Cite cycle A FAIL_LOW factor-25 as predicted | Cycle A is CLOSED-RESOLVED; D is upstream, not derivative |
| Cite F8_FORENSIC envelope factor-25 as positive evidence | Envelope is information record per F8_FORENSIC §6 explicit framing |
| Select route post-hoc by closer-to-H_0 numerical match | κ.1 anti-pattern (§3.6.5); selection rule pre-declared §3 |
| Use H_0 in any candidate formula (directly or via R_H, t_H, ρ_crit, …) | F-G-A FAIL_CIRCULAR; H_0 audit §5.5 mandatory |
| Loosen factor-10 threshold (F-G-B/D) to factor-100 if magnitude fails | Threshold inflation |
| Re-frame F-G-A FAIL_NO_DERIVATION as "implicit derivation exists" | Verdict tampering; HONEST_NEGATIVE is a valid PASS for audit |
| Add new routes post Phase 1 start | Anti-Lakatos LOCK violation (§5.2) |
| Frame as "γ-8" or continuation of cosmology cycles | Naming convention violation per F8_FORENSIC §9 |
| Modify γ-7 HALT-B or cycle A PR-018 based on D's outcome | Verdict tampering across cycles; cycle A upgrade requires separate reassessment cycle, not this one |
| Cite cycle A's Phase 1 ab initio computation (ratio 0.0406) as predicted by D | Cycle A is independent; D delivers γ-status verdict, not Λ_eff comparison |
| Introduce new fundamental constants to "fix" derivation | Cycle tests EXISTING γ derivability; introducing new constants is moving the goalposts |
| Promote γ classification (γ → β) without explicit derivation passing all routes A-E audit | §3.6.13 BINDING: classification changes require derivation evidence |

---

## §8 — Constants classification (§3.6.13 FOURTH practical application)

Constants entering Phase 1-3 derivations:

| # | Constant | Pre-cycle category | Source | Notes / cycle target |
|---|----------|---------------------|--------|----------------------|
| 1 | **γ** | **(γ) OBSERVATIONAL_ANCHOR** (current) — calibrated to H_0²/c_0² per Appendix E eq. 353 | Appendix E | **CYCLE TARGET**: F-G-A tests whether reclassify to (α) TGP_FUNDAMENTAL |
| 2 | ℓ_P | (α) TGP_FUNDAMENTAL | thm:lP (Appendix E eq. 232): ℓ_P = √(ℏ_0 G_0/c_0³) | reference scale |
| 3 | c_0 | (α) TGP_FUNDAMENTAL (reference value; vs (β) EMERGENT_FROM_PHI for local c) | sek04 reference | per §3.6.13 example: c is (β) EMERGENT_FROM_PHI, but c_0 reference value is (α) |
| 4 | ℏ_0 | (α) TGP_FUNDAMENTAL (reference) | sek04 reference | input |
| 5 | Φ_0 | (α) TGP_FUNDAMENTAL (calibrated dimensionless via g_0^e ≈ 0.87) | sek04, sek02 | candidate input for route B/C |
| 6 | α (sek02 N[Φ]) | (α) TGP_FUNDAMENTAL (= 2, derived Phase 2 universal mass formula) | sek02 LOCKED | route E input |
| 7 | β (sek02 N[Φ]) | (α) TGP_FUNDAMENTAL (= γ by vacuum condition) | sek02 LOCKED | route E input; redundant with γ |
| 8 | G_0 (Newton reference) | (γ) OBSERVATIONAL_ANCHOR (or derived via γ-5 Phase 3 G_eff = c³ℓ_P²/ℏ_0) | external CODATA / γ-5 | input to ℓ_P definition |
| 9 | H_0 | (γ) OBSERVATIONAL_ANCHOR | Planck 2018: 67.7 km/s/Mpc | **EXTERNAL TARGET** for F-G-D inversion; **FORBIDDEN as INPUT** to any F-G-A candidate formula |
| 10 | R_H (Hubble radius) | derived: c_0/H_0 | external | **FORBIDDEN as INPUT** (would inject H_0 via back door); flag for route D circularity check |
| 11 | m_sp | derived from γ: m_sp = √γ | per definition | F-G-C cross-check observable |
| 12 | l_sp | derived from γ: l_sp = 1/√γ | per definition | reference scale |
| 13 | Λ_eff | derived: γ/12 | per eq. 207 | downstream consequence (NOT this cycle's verdict) |
| 14 | Λ_obs | (γ) OBSERVATIONAL_ANCHOR | Planck 2018 | cycle A target (not this cycle) |

**Total: 14 constants classified.** Critical pre-cycle classification: γ is currently (γ) OBSERVATIONAL_ANCHOR per Appendix E eq. 353 calibration. **F-G-A PASS would reclassify γ to (α) TGP_FUNDAMENTAL.** F-G-A FAIL would confirm (γ) classification as definitive.

**No new fundamental constants introduced** (per §7 forbidden moves).

---

## §9 — Risk register

| ID | Risk | Severity | Mitigation |
|----|------|----------|------------|
| R-G-1 | All routes A-E may FAIL — cycle delivers HONEST_NEGATIVE only | LOW | Pre-disclosed §4 F-G-A FAIL_NO_DERIVATION; HONEST_NEGATIVE is valid audit PASS; cycle's value is epistemological clarification |
| R-G-2 | Route C (RG / dimensional transmutation) may exceed current TGP formalism scope (no Wilson-RG machinery for Φ⁴-class in concept paper) | HIGH | Declare PARTIAL_concept_mismatch + R1 flag for "Wilson-RG of Φ⁴-class TGP — open program" |
| R-G-3 | Route D circularity (geometric scale L_internal = R_H by silent construction) | HIGH | §5.5 mandatory H_0 audit; any geometric scale traced to non-cosmological origin or rejected |
| R-G-4 | Multi-route ambiguity (e.g., routes C and D both deliver formulas with different numerical match) | MEDIUM | §3 selection rule pre-declared (fewest inputs preferred); R1 flag if true tie under that rule |
| R-G-5 | Concept paper Appendix E formalism may be insufficient (γ defined operationally only via eq. 353) | MEDIUM | Phase 1 audit Appendix E: is γ defined intrinsically (V''(1)) or extrinsically (via m_sp ~ ℏH_0/c)? Outcome documented |
| R-G-6 | F-G-D inversion may yield H_0 prediction within Hubble tension band (67-73 km/s/Mpc) — false positive | LOW | Factor-10 tolerance F-G-D is broad enough to cover tension; honest disclosure §4 |
| R-G-7 | Sympy symbolic algebra may not capture nonperturbative effects (instantons, condensates) | MEDIUM | Phase 1 explicit assumption: perturbative + power-counting routes only; nonperturbative routes declared R1 + future scope |
| R-G-8 | Hardcoded T_pass=True risk for "yes derivation found" without first-principles | LOW | §5.1 strict 0 hardcoded; every formula tested ab initio in sympy |
| R-G-9 | Naming temptation "γ-8" (continuation framing) | LOW | §7 explicit forbidden; cycle is "op-G-substrate-derivation", not "γ-N" |
| R-G-10 | Cycle A upgrade temptation (autoplay D PASS → cycle A re-promote) | MEDIUM | §7 explicit forbidden; cycle A reassessment requires separate cycle |

---

## §10 — Phase plan

| Phase | FP target | Duration | Output |
|-------|-----------|----------|--------|
| 0 | scaffold + falsifiers + reading + route enumeration | this | Phase0_balance.md + LOCK |
| 1 | F-G-A routes A-E sympy attempts + H_0 audit per route | 1-2 sesji | Phase1_sympy.py + Phase1_derivation.md |
| 2 | F-G-B numerical substitution (if F-G-A PASS) + F-G-C consistency check | 0.5-1 sesja | Phase2_sympy.py + Phase2_results.md |
| 3 | F-G-D H_0 inversion check (if F-G-A PASS + F-G-B PASS) | 0.5 sesja | Phase3_sympy.py + Phase3_results.md (may be merged with Phase 2 if scope allows) |
| FINAL | aggregate verdict + claim_status + PR-019 LOCK entry + cycle A upgrade implication note (informational) | 0.5 sesja | Phase_FINAL_close.md + PRE_REGISTERED_FALSIFIERS.md PR entry |

**Total estimate:** 2-3 sesji.

**Decision points between phases:**
- Post Phase 1: if F-G-A returns FAIL_NO_DERIVATION across all 5 routes, cycle goes directly to FINAL (HONEST_NEGATIVE verdict; F-G-B/C/D NOT_APPLICABLE).
- Post Phase 1: if F-G-A returns PASS_PURE or PASS_WITH_PHI0, proceed to Phase 2-3 standard.
- Post Phase 1: if F-G-A returns PASS_WITH_LAGRANGIAN with multi-route ambiguity, apply §3 selection rule; if R1-tied, declare ambiguity and proceed to Phase 2 for the preferred route.
- Post Phase 2: if F-G-B FAIL_HIGH (route A Planck UV expected), document as expected outcome; cycle still has F-G-A PASS structural significance.

---

## §11 — Decision budget

| Budget | Cap | Used | Remaining |
|--------|-----|------|-----------|
| DEC (substantive decision) | 3 | 0 | 3 |
| PARTIAL_compute | 1 | 0 | 1 |
| PARTIAL_concept_mismatch | unrestricted (declare R1) | 0 | unlimited |
| Hardcoded T_pass=True | 0 | 0 | 0 |

---

## §12 — Anti-Lakatos verification (Phase 0 check)

| Item | Status |
|------|--------|
| F8 FAILs cited as motivation? | NO — motivation is Appendix E gap §1.3 ✓ |
| F8_FORENSIC envelope cited as predicted? | NO — only as informational baseline ✓ |
| Cycle A FAIL_LOW cited as predicted? | NO — cycle A is downstream consequence, not motivation ✓ |
| Route selected post-hoc by closer-to-observation match? | NO — selection rule pre-declared §3 (fewest inputs) ✓ |
| Threshold inherited from γ-7 or cycle A? | NO — factor-10 declared independently §4 with own justification ✓ |
| Cycle named γ-8? | NO — op-G-substrate-derivation reflects mechanism category ✓ |
| Pre-registered falsifiers BEFORE derivation? | YES — this document ✓ |
| Standalone fail modes declared? | YES — F-G-A/B/C/D each independent ✓ |
| Independent mechanism from F8 FAIL cycles? | YES — foundational scale derivation vs F8 mechanism test ✓ |
| H_0 circularity audit mandated? | YES — §5.5 explicit, every route audited ✓ |
| Honest disclosure of FAIL_NO_DERIVATION as valid PASS? | YES — §1.3 + §4 F-G-A explicit ✓ |
| Forbidden moves register comprehensive? | YES — §7, 12 entries ✓ |
| Cycle A upgrade auto-claimed? | NO — §1.3 explicit: cycle A upgrade requires separate reassessment cycle ✓ |
| New fundamental constants introduced? | NO — §7 explicit forbidden ✓ |

**Anti-Lakatos status:** COMPLIANT ✓

---

## §13 — Phase 0 LOCK status

**LOCKED 2026-06-01** by user authorization "Ok działaj z cyklem D".

- Scope (§1): CONFIRMED
- Routes (§3): A/B/C/D/E LOCKED with multi-route selection rule
- Falsifiers (§4): F-G-A/B/C/D LOCKED (immutable post-lock)
- Threshold §4: factor-10 (F-G-B, F-G-D) LOCKED independently
- Methodology (§5): cycle 1/2/7 strict; DEC budget 3; PARTIAL_compute 1; H_0 audit mandatory
- Inheritance (§6): LEGITIMATE concept-paper + LOCKED prior cycles only
- Forbidden moves (§7): registered; ANTI-LAKATOS COMPLIANT
- Constants (§8): §3.6.13 FOURTH practical application, 14 constants classified
- Risk register (§9): 10 items registered
- Phase plan (§10): 5 phases (0/1/2/3/FINAL)
- PR-019 reserved for Phase FINAL entry

**Phase 1 authorized upon user "działaj Phase 1" trigger.** Execute F-G-A routes A-E in sympy with mandatory H_0 audit per §5.5.

---

## §14 — Status summary

| Field | Value |
|-------|-------|
| Phase 0 status | LOCKED 2026-06-01 |
| Authorization | "Ok działaj z cyklem D" 2026-06-01 ✓ |
| F8 independence | DECLARED (different mechanism category: foundational scale derivation) ✓ |
| Foundation vs consequence distinction | DECLARED §1.3 (this cycle = γ derivability; cycle A upgrade = separate cycle) ✓ |
| Anti-Lakatos | COMPLIANT ✓ |
| Inheritance | LEGITIMATE only ✓ |
| Phase 1 prerequisites | Appendix E eq. 304-355 + sek02 + sek04 + CALIBRATION_PROTOCOL §3.6.13 |
| Expected duration | 2-3 sesji |
| Anticipated outcome distribution (informational, NOT predicted) | Route A: FAIL_HIGH (10¹²² OOM); Route B: likely FAIL_HIGH; Route C: open (most likely vehicle for nontrivial result); Route D: likely FAIL_CIRCULAR; Route E: likely FAIL (γ as overall scale). Aggregate likely: F-G-A FAIL_NO_DERIVATION with R1 flag for "Wilson-RG of Φ⁴-class TGP — open program". |
| Related cycles | B CLOSED-RESOLVED (PR-017); A CLOSED-RESOLVED (PR-018 STRUCTURAL_PARTIAL C+); C DEFERRED; this cycle D ACTIVE |
| Threshold | factor-10 (F-G-B, F-G-D) declared independently |
| Cycle A upgrade implication | informational only; explicit forbidden as motivation; separate reassessment cycle required for any cycle A status change |
