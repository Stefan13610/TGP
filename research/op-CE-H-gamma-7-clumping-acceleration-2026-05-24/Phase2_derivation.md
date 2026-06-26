---
title: "Phase 2 — V_eff(t) formal derivation + dimensional reconciliation (γ-7 v2)"
type: phase_derivation
status: LOCKED_PENDING_USER_REVIEW
phase: 2
parent_cycle: op-CE-H-gamma-7-clumping-acceleration-2026-05-24
parent_handoff: meta/HANDOFF_GAMMA_7_2026-05-24.md
parent_phase1: Phase1_derivation.md
preregistration_version: v2 (field-based; HANDOFF §11.3 final formula CORRECTED Phase 2 per derivation)
created_date: 2026-05-24
authorization: "User 'tak działaj z fazą 2' 2026-05-24 (post Phase 1 LOCK)"
substantive_fp_total: 10
substantive_fp_pass: 9
hardcoded_T_pass_count: 0
partial_compute_count: 0
partial_compute_cumulative: "1/1 (Phase 1 T_P1_10; Phase 2 0)"
dec_budget_used: "1/3 (DEC 1 Phase 1; Phase 2 NIE additional DEC)"
substantive_findings:
  - "R1 #14 dimensional reconciliation RESOLVED: V_eff = ∫⟨Φ⟩²/v² dV literal volume integral as primary"
  - "V_eff(t) corrected primary equation derived; HANDOFF v2 §11.3 final formula superseded"
  - "V_baseline(t) defined precisely (uniform N-source distribution reference)"
  - "Energy conservation: V_eff = geometric measure, NIE energy quantity"
  - "v_phi convention LOCKED: Appendix E v²≈25 Planck units → SI v²≈3.03×10⁴⁵ J/m"
  - "F-γ-7-B preview FAIL direction CONFIRMED (ratio ≈ 10⁻³ outside factor 10)"
falsifier_status_update:
  F-γ-7-A: "STRUCTURALLY_VERIFIED + DIMENSIONALLY_RECONCILED (Phase 2)"
  F-γ-7-B: "PREVIEW_FAIL direction CONFIRMED Phase 2; Phase 5 full test pending (ξ_clump dependence)"
  F-γ-7-C: "DEFERRED to Phase 3 (ξ_clump derivation)"
  F-γ-7-D: "DEFERRED to Phase 3+4 (timing)"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — V_eff(t) formal derivation + dimensional reconciliation (γ-7 v2)

**Status:** LOCKED_PENDING_USER_REVIEW.
**Authorization:** User "tak działaj z fazą 2" 2026-05-24 (post Phase 1 LOCK).
**Methodology:** Strict cycle 1/2/7 (0 hardcoded T_pass=True); compute-then-compare; PARTIAL_compute budget preserved 1/1 (Phase 1 used; Phase 2 NIE additional).

---

## §1 — Phase 2 scope (per Phase 1 §8.1 extended)

**Per BINDING contract + Phase 1 findings:**
1. Resolve R1 #14 dimensional reconciliation (action-derived vs literal volume V_eff)
2. Define V_baseline(t) precisely (Q4 finalize)
3. Derive V_eff(t) full equation z time-dependence
4. Energy conservation check (R5 mitigation)
5. v_phi normalization convention finalization (Q9)
6. F-γ-7-B preview update z corrected V_eff formula

**Key Phase 2 disposition:** Phase 1 identified ambiguity in HANDOFF v2 §11.3 between two physical observables (literal volume integral vs action-derived). Phase 2 chooses primary measure z TGP-native justification, NIE z post-hoc fitting.

---

## §2 — Substantive FP results (sympy executed)

**Script:** [[Phase2_sympy.py]] (executed; full output preserved).

| FP ID | Description | Verdict | Notes |
|-------|-------------|---------|-------|
| T_P2_1 | Diagonal self-energy ∫(δΦ_j)² d³x = q_j²/(8π μ_sp) | **PASS** ✓ | Spherical integration; UV-finite (Yukawa pole at r=0 regularized by Appendix E ℓ_P) |
| T_P2_2 | Off-diagonal pair-overlap ∫δΦ_i δΦ_j d³x = q_i q_j exp(-μ_sp r_ij)/(8π μ_sp) | **PASS** ✓ | Phase 1 T_P1_4 reconfirmed via ∂/∂μ² method |
| T_P2_3 | Linear single-source ∫δΦ_j d³x = q_j/μ_sp² (POSITION-INDEPENDENT) | **PASS** ✓ | Independent of x_j; cancels in V_eff - V_baseline |
| T_P2_4 | V_eff(t) primary equation z literal volume integral | **PASS** ✓ | Corrected formula derived |
| T_P2_5 | Dimensional consistency m³ check | **PASS** ✓ | R1 #14 RESOLVED |
| T_P2_6 | Reformulation z ξ_clump enhancement factor | **PASS** ✓ | Phase 3 hook prepared |
| T_P2_7 | Energy conservation analysis | **PASS** ✓ | V_eff = geometric, NIE energy; Newton-Yukawa consistent z γ-5 |
| T_P2_8 | v_phi convention finalization (Q9) | **PASS** ✓ | Appendix E v²≈25 Planck → SI v²≈3.03×10⁴⁵ J/m LOCKED |
| T_P2_9 | F-γ-7-B numerical preview update | **FAIL direction CONFIRMED** | Ratio ≈ 1.1×10⁻³ outside factor 10 |
| T_P2_10 | γ-5 Phase 3 cross-check | **PASS** ✓ | Yukawa form + G_eff preserved; γ-7 V_eff = separate observable |

**Summary:**
- **9/10 substantive FP PASS**
- **0 hardcoded T_pass=True** (strict cycle 1/2/7 preserved) ✓
- **0 PARTIAL_compute** in Phase 2 (cumulative 1/1 from Phase 1 T_P1_10) ✓
- **DEC 0 additional** (Phase 2 NIE used DEC; cumulative 1/3)

---

## §3 — R1 #14 dimensional reconciliation (RESOLVED)

### §3.1 — The ambiguity (Phase 1 finding)

HANDOFF v2 §11.3 defined V_eff = ∫⟨Φ⟩²/v² dV (literal volume integral) **AND** quoted final formula z 1/(4π r_ij) pair-overlap. Phase 1 identified these as **inconsistent**: literal integral gives 1/(8π μ_sp) form, NIE 1/(4π r_ij) form.

Two physical interpretations possible:
- **(I) Literal volume**: V_eff = ∫⟨Φ⟩²/v² dV — pair-overlap = q_i q_j exp(-μ_sp r_ij)/(8π μ_sp). Dimensionally **m³** ✓.
- **(II) Action-derived**: V_eff = -2 S_field/v² — pair-overlap = q_i q_j exp(-μ_sp r_ij)/(4π r_ij). Dimensionally **m** (length, NIE volume) ✗.

### §3.2 — Phase 2 resolution: Literal volume integral is PRIMARY

**Decision: Interpretation (I) literal volume integral chosen as primary V_eff measure.**

**Rationale (TGP-native, NIE post-hoc):**

1. **K2 ontology compliance** (concept paper §1.1):
   > "Masa i przestrzeń są dwiema stronami tej samej konfiguracji. Brak masy = brak struktury = brak rozróżnialnej przestrzeni."
   
   V_eff = ∫⟨Φ⟩²/v² dV directly measures "amount of Phi-structure × volume" — most TGP-native interpretation.

2. **Dimensional consistency** (R1 #14):
   - Literal volume integral form has units m³ correctly z standard SI.
   - Action-derived form gives length per pair (would require ad-hoc area prefactor to give volume).

3. **Vacuum reference clean**:
   - Uniform Φ = v: V_eff = ∫(v²/v²)dV = V_geometric (Euclidean volume).
   - Sources contribute via (δΦ/v + (δΦ/v)²) corrections — standard fluctuation expansion.

4. **HANDOFF v2 §11.3 SUPERSEDED** for final formula. HANDOFF retained for:
   - V_eff = ∫⟨Φ⟩²/v² dV definition (preserved correctly)
   - Yukawa Green's function δΦ_j = q_j exp(-μ_sp r)/(4π r) (preserved correctly)
   - Algebraic slip in pair-overlap derivation (CORRECTED Phase 2)

### §3.3 — Anti-Lakatos disposition

**This refinement is LEGITIMATE per CALIBRATION_PROTOCOL §0.3:**
- ✓ Derivation-based correction (NIE post-hoc rescue)
- ✓ Worsens F-γ-7-B (ratio remains ~10⁻³); NIE attempt to escape FAIL
- ✓ HANDOFF v2 definition (∫⟨Φ⟩²/v² dV) PRESERVED; only derivation route corrected
- ✓ γ-5 Phase 3 LOCKED preserved unchanged (γ-5 used action-derived for gravitational force; γ-7 V_eff = different observable)
- ✓ Pre-registration v2 thresholds UNCHANGED (factor 10 around Ω_DE; z_onset [0.3, 1.0])

**Phase 2 documents this as substantive Phase 2 finding, NIE retroactive pre-registration modification.**

---

## §4 — V_eff(t) primary equation (Phase 2 CORRECTED final form)

### §4.1 — Full derivation

Per literal volume integral interpretation, expand ⟨Φ⟩² around vacuum v + Σ_j δΦ_j:

$$\frac{\langle\Phi\rangle^2}{v^2} = 1 + \frac{2}{v}\sum_j \delta\Phi_j + \frac{1}{v^2}\left(\sum_j \delta\Phi_j\right)^2$$

Integration over all space:

$$\int \frac{\langle\Phi\rangle^2}{v^2}\, dV = V_{\text{geom}} + \frac{2}{v}\sum_j \int \delta\Phi_j\, dV + \frac{1}{v^2}\int\left(\sum_j \delta\Phi_j\right)^2 dV$$

**Component evaluation (T_P2_1 through T_P2_3):**

| Term | Value | Position-dependence |
|------|-------|---------------------|
| V_geom | ∫dV (Euclidean) | None |
| Linear (2/v)·∫δΦ_j dV | (2/v)·(q_j/μ_sp²) | None (boundary at infinity → q_j/μ_sp² independent of x_j) |
| Quadratic diagonal (1/v²)·∫δΦ_j² dV | q_j²/(8π μ_sp v²) | None (self-energy depends only on q_j) |
| Quadratic off-diag (1/v²)·Σ_{i≠j}∫δΦ_i δΦ_j dV | q_i q_j exp(-μ_sp r_ij)/(8π μ_sp v²) | **YES (depends on r_ij)** |

### §4.2 — V_baseline(t) definition (Q4 RESOLVED)

**V_baseline(t) ≡ ∫⟨Φ⟩²/v² dV evaluated for UNIFORM N-source distribution at time t.**

For uniform distribution at time t, all of:
- Linear terms: same (position-independent)
- Diagonal terms: same (q_j² values fixed)
- Off-diagonal: uses ⟨exp(-μ_sp r_ij)⟩_uniform (average over uniform pair PDF)

So V_eff - V_baseline retains only the off-diagonal CLUMPING contribution:

$$\boxed{V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{1}{8\pi \mu_{sp} v_\phi^2}\sum_{i\ne j} q_i q_j \left[e^{-\mu_{sp} r_{ij}(t)} - \left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_{\text{uniform}}\right]}$$

**For identical sources (q_i = q ≡ q_proton) + N >> 1:**

$$V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{N^2 q^2}{8\pi \mu_{sp} v_\phi^2}\left[\left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_t - \left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_{\text{uniform}}\right]$$

### §4.3 — ξ_clump reformulation (Phase 3 hook)

Per HANDOFF §11.4 spirit, define ξ_clump(t) as enhancement factor:

$$\left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_t = \left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_{\text{uniform}} \cdot [1 + \xi_{\text{clump}}(t)]$$

Then:

$$\boxed{V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{N^2 q^2 \left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_{\text{uniform}}}{8\pi \mu_{sp} v_\phi^2} \cdot \xi_{\text{clump}}(t)}$$

**This is the γ-7 v2 Phase 2 LOCKED primary equation.**

**Properties:**
- ξ_clump = 0 → V_eff = V_baseline (no clumping enhancement)
- ξ_clump(t) > 0 → V_eff > V_baseline (clumping → more effective space)
- Time-dependence enters ONLY through ξ_clump(t) — Phase 3 task

### §4.4 — Substituting Candidate B q² = 4π G_eff m²

$$\boxed{V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{G_{\text{eff}} M^2 \left\langle e^{-\mu_{sp} r_{ij}}\right\rangle_{\text{uniform}}}{2 \mu_{sp} v_\phi^2} \cdot \xi_{\text{clump}}(t)}$$

with M = N · m_proton ≈ M_universe.

**This expresses V_eff(t) entirely in terms of:**
- γ-5 LOCKED G_eff (Newton's gravitational constant; OBSERVATIONAL_ANCHOR)
- M_universe (OBSERVATIONAL_ANCHOR)
- μ_sp = H_0/c (Appendix E + OBSERVATIONAL_ANCHOR for H_0)
- v_phi² (TGP_FUNDAMENTAL z Appendix E rem. naturalness; LOCKED Phase 2 T_P2_8)
- ⟨exp⟩_uniform (geometric factor from uniform pair PDF)
- ξ_clump(t) (PHASE 3 task — TGP-native source dynamics derivation)

---

## §5 — Energy conservation (R5 mitigation, T_P2_7)

**Concern:** Does V_eff(t) growth imply "creation of space ex nihilo" violating energy conservation?

**Phase 2 resolution:**

### §5.1 — V_eff is GEOMETRIC, NIE energy

V_eff measures **integrated Phi-density × volume**, a geometric/configurational quantity. It is NIE an energy reservoir.

### §5.2 — Energy of multi-source configuration

Per Phase 1 T_P1_3 + γ-5 Phase 3 LOCKED:
- Pair-interaction energy E_pair(r_ij) = -q_i q_j exp(-μ_sp r_ij)/(4π r_ij) (attractive Yukawa convention)
- Substituting q² = 4π G_eff m²: E_pair = -G_eff m_i m_j exp(-μ_sp r_ij)/r_ij (Newton-Yukawa)

### §5.3 — Clumping dynamics

When sources approach (clumping):
- r_ij decreases → E_pair more negative → gravitational binding energy grows
- Released energy → kinetic energy of structure formation → eventually radiated/thermalized
- This is **standard gravitational dynamics** (Newton-Yukawa modified by Hubble screening)

V_eff growth from clumping:
- ⟨exp(-μ_sp r_ij)⟩ increases (more close pairs) → V_eff geometric grows
- **BUT total energy of Phi field DECREASES** (configuration becomes more bound)

**No energy creation ex nihilo.** Universe expansion driven by V_eff growth is **geometric consequence of substrate Phi configuration** — analogous to how GR expansion is geometric, NIE energetic.

### §5.4 — Verdict

Energy conservation: **PRESERVED** structurally.

**R5 RISK MITIGATED.** V_eff = geometric substrate measure; clumping is gravitationally bound (energy released as kinetic/thermal), NIE "free energy from clumping".

---

## §6 — v_phi convention LOCKED (Q9, T_P2_8)

### §6.1 — Appendix E natural units

Appendix E rem. naturalness LOCKED:
> "λ_UV² · γ ∼ 10⁻¹²² (Planck scale); próżniowym Φ_0 ≈ 25"

Convention: Planck units (ℏ = c = ℓ_P = 1).

### §6.2 — Canonical scalar field 3+1D

For canonical scalar field w (3+1)D Lagrangian density [L] = E/L³:
- [Φ²·μ²] = [E/L³] → [Φ²] = [E·L]/[L²·μ²] = [E·L⁻¹]·[1] = [E/L]
- [Φ] = √[E/L]
- [v_phi²] = [Φ²] = [E/L] = J/m in SI

### §6.3 — Conversion to SI

Planck-[J/m] = E_P/ℓ_P:
- E_P = √(ℏc⁵/G) ≈ 1.96 × 10⁹ J
- ℓ_P = √(ℏG/c³) ≈ 1.62 × 10⁻³⁵ m
- **Planck-[J/m] ≈ 1.21 × 10⁴⁴ J/m**

Therefore:

$$\boxed{v_\phi^2 \approx 25 \cdot \frac{E_P}{\ell_P} \approx 3.03 \times 10^{45} \text{ J/m}}$$

**Q9 RESOLVED.** v_phi² LOCKED for γ-7 numerical calculations.

### §6.4 — Caveat (alternative conventions)

If Φ defined differently (e.g., dimensionless via field redefinition Φ' = Φ/M with mass scale M), v² shifts by appropriate power of M. Current LOCKED convention assumes canonical 3+1D scalar (most standard).

Potential future R1 flag (#16 candidate): convention sensitivity dla cosmological field measurements.

---

## §7 — F-γ-7-B numerical preview update (T_P2_9 — FAIL direction CONFIRMED)

### §7.1 — Updated formula z corrected V_eff

Using Phase 2 corrected primary equation (§4.4):

$$V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{G_{\text{eff}} M^2 \left\langle e^{-\mu_{sp} r}\right\rangle_{\text{uniform}} \xi_{\text{clump}}(t)}{2 \mu_{sp} v_\phi^2}$$

### §7.2 — F-γ-7-B test: set V_eff = Ω_DE · V_univ

$$v_\phi^2_{\text{required}} = \frac{G_{\text{eff}} M^2 \left\langle e^{-\mu_{sp} r}\right\rangle_{\text{uniform}} \xi_{\text{clump}}}{2 \mu_{sp} \Omega_{DE} V_{\text{univ}}}$$

**Numerical (T_P2_9):**

| Quantity | Value |
|----------|-------|
| G_eff | 6.674 × 10⁻¹¹ m³/(kg·s²) |
| M_univ | ~10⁵³ kg |
| μ_sp = H_0/c | 7.7 × 10⁻²⁷ m⁻¹ |
| V_univ (Hubble) | ~9.3 × 10⁷⁸ m³ |
| Ω_DE | 0.7 |
| ⟨exp⟩_uniform | ~0.5 (geometric estimate) |
| ξ_clump | ~1 (order-unity for non-linear epoch; PROVISIONAL) |

**Result:**
- v²_required ≈ **3.33 × 10⁴² J/m**
- v²_natural (Appendix E) ≈ **3.03 × 10⁴⁵ J/m**
- **Ratio ≈ 1.10 × 10⁻³** (log₁₀ ≈ -2.96)

### §7.3 — Verdict: OUTSIDE factor 10 → F-γ-7-B preview FAIL direction CONFIRMED

**Phase 2 dimensional reconciliation does NIE close the gap.**

The Phase 1 preview FAIL direction is **robust** against Phase 2 corrections.

### §7.4 — Honest disposition (anti-Lakatos)

**Per HANDOFF prompt sesja #8 + Forbidden move #19:**
> "γ-7 musi rozwiązać F8 albo declare honest FAIL (HALT-B) — to jest third attempt. Anti-Lakatos: NIE pivot post-hoc do yet another mechanism."

**NIE post-hoc adjustments allowed:**
- ❌ NIE adjust v_phi² to match (LOCKED z Appendix E)
- ❌ NIE postulate larger N via dark matter (forbidden move #19; needs TGP-native dark matter derivation w γ-8+)
- ❌ NIE invent additional enhancement factors
- ❌ NIE change V_eff measure definition (Phase 2 LOCKED literal volume integral)

**Legitimate Phase 3+ pathways (within current scope):**

1. **ξ_clump(t) >> 1 in non-linear regime?** — Phase 3 derives. Even if ξ_clump ~ 10², ratio shifts only to 10⁻¹ — still outside factor 10.

2. **More accurate ⟨exp⟩_uniform?** — Could be 0.1 to 1.0 range (geometric factor). Shifts ratio by factor 2-5 at most.

3. **Effective M via clustering hierarchy?** — Phase 3 may show large clusters effectively contribute more via composite Phi-charge. But this is speculative and requires explicit derivation.

**Most likely Phase 5 verdict:** F-γ-7-B FAIL.

**Anti-Lakatos disposition:** Document honestly; NIE pivot to γ-8 z new mechanism. If F-γ-7-B FAIL confirmed Phase 5, declare HALT-B per Phase 0 §13.3 acknowledgment.

### §7.5 — IMPORTANT — NIE pre-empt FAIL declaration

**T_P2_9 is PREVIEW (continuation of Phase 1 T_P1_10 PARTIAL_compute):**
- Order-of-magnitude only
- Uses provisional ξ_clump ~ 1 (Phase 3 will derive properly)
- Uses approximate ⟨exp⟩_uniform ~ 0.5 (Phase 3 may refine)

**Phase 5 does definitive F-γ-7-B test** z Phase 3 derived ξ_clump(t).

**Phase 2 disposition:** Preview FAIL direction CONFIRMED, but Phase 5 verdict awaits Phase 3 ξ_clump derivation. NIE auto-FAIL declaration here.

---

## §8 — γ-5 Phase 3 cross-check (T_P2_10, anti-Lakatos verification)

### §8.1 — Yukawa form preserved

Both γ-5 (gravitational force) and γ-7 (V_eff measure) use SAME Yukawa Green's function:
$$\delta\Phi_j(r) = \frac{q_j e^{-\mu_{sp} r}}{4\pi r}$$

Same Appendix E eq. 101 KG propagator foundation.

### §8.2 — Distinct physical observables

| Observable | Formula | Physical meaning |
|------------|---------|------------------|
| γ-5 E_pair (action-derived) | q_i q_j exp(-μ_sp r_ij)/(4π r_ij) | Pair-interaction ENERGY (gravitational binding) |
| γ-7 V_eff (literal volume) | ∫⟨Φ⟩²/v² dV − V_baseline | Effective substrate VOLUME |

These are **different functionals of the same Phi configuration** — NIE contradiction.

### §8.3 — γ-5 LOCKED preserved

- F = -q²/(4π r²) → Newton form (γ-5 Phase 3 LOCKED)
- G_eff = c³·ℓ_P²/ℏ identification (γ-5 Phase 3 LOCKED)
- F-γ-5-A PASS_CALIBRATED, F-γ-5-B PASS_MARGINAL preserved

**γ-5 LOCKED status PRESERVED unchanged.** γ-7 V_eff is separate observable z same underlying field theory.

---

## §9 — Anti-Lakatos verification (Phase 2)

| Check | Status |
|-------|--------|
| γ-3 + γ-3' + γ-5 B+ verdicts modified? | NO ✓ (all LOCKED preserved) |
| F-γ-7-A/B/C/D pre-registered thresholds modified? | NO ✓ (factor 10 + z [0.3,1.0] LOCKED) |
| Pre-registration v2 definition (V_eff = ∫⟨Φ⟩²/v² dV) modified? | NO ✓ (literal interpretation matches pre-reg definition) |
| HANDOFF v2 §11.3 final formula corrected based on derivation? | YES (legitimate Phase 2 finding per §0.3; documented transparently) |
| Forbidden move #19 (postulate q to match Ω_DE)? | NO ✓ (q derived z γ-5 Newton matching) |
| Forbidden move #20 (mean-field aggregate)? | NO ✓ (v2 field-based ACTIVE) |
| Cycle 1/2/7 violated (hardcoded T_pass)? | NO ✓ (0/10 Phase 2; cumulative 0/20) |
| F-γ-7-B preview FAIL direction → post-hoc rescue attempted? | NO ✓ (documented honestly; legitimate Phase 3+ pathways flagged) |
| Pre-emptive FAIL declaration before Phase 5? | NO ✓ (Phase 3 ξ_clump derivation still pending) |
| Tautology paths? | NO ✓ (V_eff measure dimensionally well-defined; pair-overlap derivation z first principles) |
| PARTIAL_compute budget? | 1/1 cumulative (Phase 1 T_P1_10; Phase 2 0) ✓ |
| DEC budget? | 1/3 cumulative; Phase 2 0 additional ✓ |

**Anti-Lakatos LOCK PRESERVED across γ-3 + γ-3' + γ-5 + γ-7 v2 Phase 1-2.**

---

## §10 — R1 flag updates (Phase 2)

### §10.1 — Updates to Phase 1 R1 candidates

**R1 #13 (HANDOFF v2 §11.3 algebraic slip):** RESOLVED Phase 2 — literal volume integral chosen; HANDOFF formula superseded. **Disposition: CLOSED via Phase 2 §3.**

**R1 #14 (V_eff dimensional consistency):** RESOLVED Phase 2 — literal volume integral gives m³ dimensionally clean. **Disposition: CLOSED via Phase 2 §3 + §5.**

**R1 #15 (F-γ-7-B preview FAIL direction):** CONFIRMED Phase 2 — ratio remains ~10⁻³ outside factor 10. **Disposition: STAYS HIGH; Phase 5 definitive (post Phase 3 ξ_clump derivation).**

### §10.2 — NEW R1 from Phase 2

**R1 #16 (NEW, sesja #8 Phase 2):** v_phi convention sensitivity dla cosmological measurements.
- **Pattern:** F-γ-7-B numerical depends on convention assumption (canonical scalar vs alternative field redefinitions).
- **Severity:** LOW (current convention is standard); MEDIUM if alternative conventions are explored.
- **Disposition:** CANDIDATE; documented Phase 2 §6.4 caveat. NIE blocking γ-7 progression.

---

## §11 — Open questions update (post Phase 2)

| Q | Phase 1 status | Phase 2 update |
|---|----------------|----------------|
| Q1: Soliton charge integral well-defined? | Confirmed structurally | Unchanged |
| Q2: TGP gravitational instability | Deferred Phase 3 | Unchanged |
| Q3: TGP virialization analog | Deferred Phase 3 | Unchanged |
| Q4: V_baseline definition | Phase 2 task | **RESOLVED** § 4.2 (uniform N-source reference at time t) |
| Q5: Integration measure | Confirmed PROPER | Unchanged |
| Q6: Self-consistency iteration depth | First iteration | Unchanged |
| Q7: w_eff(t) from V_eff(t) | Deferred Phase 4-5 | Unchanged |
| Q8: Dimensional reconciliation | OPEN | **RESOLVED** §3 (literal volume integral) |
| Q9: v_phi normalization | OPEN | **RESOLVED** §6 (Appendix E v²≈25 Planck → SI ≈ 3.03×10⁴⁵ J/m) |
| Q10: Dark matter sources | NIE included | UNCHANGED — outside γ-7 scope (γ-8+) |

**3 questions resolved in Phase 2** (Q4, Q8, Q9). All remaining QC deferred to appropriate phases.

---

## §12 — F-γ-7-A v2 verdict update

**Pre-registration v2 (LOCKED 2026-05-24):**
> "V_eff(t) MUST be derivable jako functional of ⟨Φ⟩(x,t) z multi-source Yukawa configuration."

**Phase 2 update:** **STRUCTURALLY_VERIFIED + DIMENSIONALLY_RECONCILED**
- Derivation z KG propagator (Phase 1 T_P1_1)
- Multi-source Yukawa configuration (Phase 1 T_P1_2)
- Pair-overlap term explicit (Phase 1 T_P1_4 + Phase 2 T_P2_2)
- V_eff form derived (Phase 2 T_P2_4) — dimensionally m³ ✓
- V_baseline(t) defined precisely (Phase 2 T_P2_4)
- NIE mean-field aggregate ✓ (forbidden move #20 preserved)
- q derived z TGP fundamentals via Candidate A + B (Phase 1 §3.4)

**F-γ-7-A v2 PASS structurally + dimensionally** — full structural deliverable achieved.

---

## §13 — Phase 3 prep notes

### §13.1 — Phase 3 inheritance from Phase 2

Phase 2 LOCKED primary equation ready dla Phase 3:

$$V_{\text{eff}}(t) - V_{\text{baseline}}(t) = \frac{N^2 q^2 \langle e^{-\mu_{sp} r}\rangle_{\text{uniform}}}{8\pi \mu_{sp} v_\phi^2} \cdot \xi_{\text{clump}}(t)$$

**Phase 3 task (DEC 2):** Derive ξ_clump(t) growth z TGP-native source dynamics (NIE Press-Schechter).

### §13.2 — Phase 3 scope per Phase 0 §5.2 + HANDOFF §11.7

**Candidates pre-registered:**
- **Candidate A PRIMARY**: Linear gravitational instability z γ-5 1/r force (analog δ̈ + 2Hδ̇ - 4πG_eff ρ̄ δ = 0)
- **Candidate B NON-LINEAR**: Virialization via local Phi-saturation (γ-5 Phase 2 c(n_critical)=0)
- **Candidate C EXTENSION**: S_creation coupling (concept paper EQ-5)
- **Candidate D FORBIDDEN**: Press-Schechter borrowing (forbidden move #18)

### §13.3 — Phase 3 elevated importance

Phase 3 is CRITICAL for γ-7 cycle outcome:
- ξ_clump(t) determines F-γ-7-B numerical (currently preview ratio ~10⁻³)
- ξ_clump(t) determines F-γ-7-C acceleration condition (d²V_eff/dt² > 0)
- ξ_clump(t) determines F-γ-7-D timing (z_onset)

**Phase 3 outcomes will determine γ-7 final verdict** (A+/A/A-/B+/HALT-B per README §8).

---

## §14 — Cross-references

- Phase 1: [[Phase1_derivation.md]] (Phase 1 substantive findings carried forward)
- Phase 0: [[Phase0_balance.md]] (§3.6.13 + DEC pre-allocation)
- HANDOFF: [[../../meta/HANDOFF_GAMMA_7_2026-05-24.md]] §11 (v2 pre-registration)
- README: [[README.md]] §10 (v2 BINDING contract)
- γ-5 Phase 3: [[../op-CE-H-gamma-5-c-interpretation-2026-05-24/Phase_FINAL_close.md]] §3.3 (G_eff identification + Yukawa LOCKED)
- Appendix E: [[../../core/formalizm/dodatekE_kwantyzacja.tex]] (KG propagator eq. 101 + phonon dispersion eq. 350-353 + dark energy REMARK eq. 365)
- Concept paper: [[../../meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md]] §1.1 K2 + §3.2 Lagrangian + §3.3 D2 soliton
- Sympy script: [[Phase2_sympy.py]]
- Falsifiers: [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §15 (F-γ-7 v2)

---

## §15 — Phase 2 LOCK declaration

**Phase 2 LOCKED PENDING USER REVIEW.**

### §15.1 — Summary statistics

- ✅ **10/10 substantive FP executed** in Phase 2
- ✅ **9/10 PASS** (T_P2_1 through T_P2_8 + T_P2_10)
- ⚠ **1/10 FAIL direction confirmed** (T_P2_9 F-γ-7-B preview — NIE binding; Phase 5 definitive)
- ✅ **0 hardcoded T_pass=True** in Phase 2 (strict cycle 1/2/7 preserved)
- ✅ **0 additional PARTIAL_compute** in Phase 2 (cumulative 1/1 from Phase 1)
- ✅ **0 additional DEC** in Phase 2 (cumulative 1/3)
- ✅ **Anti-Lakatos LOCK preserved**

### §15.2 — Substantive Phase 2 deliverables

1. **V_eff(t) primary equation LOCKED** (corrected, dimensionally clean):
   V_eff(t) - V_baseline(t) = (N² q² ⟨exp⟩_uniform · ξ_clump(t))/(8π μ_sp v²)

2. **R1 #14 dimensional reconciliation RESOLVED** via literal volume integral choice

3. **V_baseline(t) defined precisely** (Q4 RESOLVED)

4. **Energy conservation: V_eff = geometric measure** (R5 mitigation; Newton-Yukawa consistent z γ-5)

5. **v_phi convention LOCKED** (Q9 RESOLVED; Appendix E v²≈25 Planck → SI ≈ 3.03×10⁴⁵ J/m)

6. **F-γ-7-B preview FAIL direction CONFIRMED** (ratio ~10⁻³ outside factor 10; honest disposition)

7. **γ-5 LOCKED preserved** (Yukawa form + G_eff unchanged; γ-7 V_eff = separate observable)

### §15.3 — F-γ-7-A v2 verdict

**STRUCTURALLY_VERIFIED + DIMENSIONALLY_RECONCILED** — Phase 2 closes F-γ-7-A v2 question.

### §15.4 — F-γ-7-B preview disposition

**FAIL direction CONFIRMED at current understanding** — but **NIE pre-empt FAIL declaration**. Phase 5 definitive test pending Phase 3 ξ_clump(t) derivation.

### §15.5 — Phase 3 authorization gate

**Phase 3 (ξ_clump(t) TGP-native structure formation theory derivation) AWAITS EXPLICIT USER AUTHORIZATION.**

Phase 3 scope per HANDOFF §11.7 + README §10.5 + Phase 0 §5.2:
1. Linear gravitational instability z γ-5 1/r force (Candidate A)
2. Non-linear virialization via local Phi-saturation (Candidate B)
3. ξ_clump(t) growth equation
4. Connection to ⟨exp(-μ_sp r_ij)⟩_t evolution

**Estimated:** 1-2 sesji (HIGH uncertainty per Phase 0 §13 R1 CRITICAL).
**DEC 2 will be used** for growth model selection.

### §15.6 — Honest expectation

Per Phase 2 §7.4 disposition + anti-Lakatos:
- F-γ-7-B preview FAIL direction is ROBUST under Phase 2 reconciliation
- Phase 3 ξ_clump derivation may close gap partially OR confirm FAIL
- **Most likely γ-7 outcome: B+ z F8 still UNEXPLAINED, OR HALT-B z honest declaration**
- **NIE pivot do new mechanism post-Phase-5 FAIL** (anti-Lakatos LOCK)

---

**END OF PHASE 2 — V_eff(t) formal derivation LOCKED 2026-05-24 (PENDING USER REVIEW)**

**F-γ-7-A v2: STRUCTURALLY_VERIFIED + DIMENSIONALLY_RECONCILED.**

**F-γ-7-B v2: PREVIEW_FAIL direction CONFIRMED (Phase 5 definitive pending).**

**Anti-Lakatos LOCK PRESERVED across γ-3 + γ-3' + γ-5 + γ-7 v2 Phase 1-2.**

**Forbidden moves #18-20 ENFORCED.**
