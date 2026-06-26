---
title: "Phase 3+4+5 batch plan — gravity synthesis + F8 re-test + Schwarzschild R_s/δt/t"
type: batch_phase_plan
status: LOCKED
phases: [3, 4, 5]
parent_cycle: op-CE-H-gamma-5-c-interpretation-2026-05-24
authorization_date: 2026-05-24
authorization_quote: "Phase 3+4+5 batch"
---

# Phase 3+4+5 batch plan

**Status:** LOCKED 2026-05-24. **Authorization:** User "Phase 3+4+5 batch" 2026-05-24.

---

## §1 — Phase 3: Gravity-as-configuration-constraint (HANDOFF §3.8)

### §1.1 — Objectives

1. Formal synthesis of Q1 (multi-source frustration → c saturates globally) and Q3 (local density → c drops locally) — HANDOFF §3.8
2. Derive gravity emergence from forced Phi-gradient overlap (NIE postulate metric)
3. Compute pair-source Phi overlap integral
4. Identify Newtonian/GR analog in far-field limit
5. Derive G_eff in terms of TGP parameters

### §1.2 — Derivation strategy

**Per HANDOFF §3.8 user explicit:**
> "globalnie cząstki chcą się oderwać, ale grawitacja wynika z nakładania się gradientów"

**Mathematical synthesis:**

For two solitons at positions x_1, x_2 with separation r:
- Each soliton's Phi field decays Yukawa-like: Φ_i(x) ~ q·exp(-m_σ|x-x_i|)/(4π|x-x_i|)
- Overlap integral: E_overlap(r) = ∫ ∇Φ_1 · ∇Φ_2 dV
- In massless limit (long-range, m_σ → 0): E_overlap ∝ -q²/r (Coulomb-like)

**Configuration constraint interpretation:**
- E_overlap > 0 → sources can spread (repulsive at short range — Q1 frustration)
- ∇E_overlap as function of r → effective force F = -∇E
- Far-field 1/r potential → 1/r² force → **Newtonian gravity analog**

**Per anti-Lakatos forbidden move #15:** Derivation MUST come from configuration counting, NIE from postulating g_μν.

**Substantive FP tests:**
- T_P3_1: Yukawa Green's function form verification
- T_P3_2: Pair overlap integral computation (massless limit → 1/r)
- T_P3_3: Effective force F = -dE/dr → 1/r² scaling
- T_P3_4: G_eff derivation from TGP parameters
- T_P3_5: §3.8 Q1 + Q3 reconciliation verification (globally repulsive, locally attractive emerges from same overlap)

---

## §2 — Phase 4: F8 re-test under c(N(t)) framework

### §2.1 — Objectives

1. Apply Phase 1 c(N) saturating form to N(t) cosmological evolution
2. Compute R(t) integral under c(N(t)) variation
3. Verify F8 thresholds: w_eff ∈ [-1.2, -0.8], ä > 0
4. Pre-derivation honest: predicted FAIL_LITERAL per Phase 1 saturation finding

### §2.2 — Derivation strategy

**Cosmological N(t) growth (per γ-3 Phase 2 + frontier creation):**
- Frontier area A(t) = 4π·R(t)² ∝ t² (linear R = c·t)
- Creation rate S_creation = 3Hv (per γ-3 (EQ-5))
- N(t) ∝ ∫ A(t') dt' ∝ ∫ t'² dt' ∝ t³

**c(N(t)) evaluation:**
- Phase 1: c(N) saturates to c_0 by N≈11
- For N(t) ∝ t³, even at very early times N(t) >> 11
- c(t) ≈ c_0 throughout observable epoch
- R(t) = ∫ c(t') dt' ≈ c_0·t (linear, same as γ-3)
- ä = 0; w_eff = -1/3

**Substantive FP tests:**
- T_P4_1: N(t) cosmological growth: cube-time scaling
- T_P4_2: c(N(t)) evaluation at multiple epochs (recombination, matter era, today)
- T_P4_3: R(t) integral computation — confirm linear (same as γ-3)
- T_P4_4: w_eff computation = -1/3 (NIE in [-1.2, -0.8])
- T_P4_5: F8 verdict FAIL_LITERAL confirmed (compute-then-compare)

---

## §3 — Phase 5: Schwarzschild R_s + gravitational time dilation

### §3.1 — Objectives

1. Compute R_s from Phase 2 c(n_local) → critical density condition
2. Test F-γ-5-A: R_s within factor 2 of 2GM/c² for {M_⊙, 1.4 M_⊙ NS, M_⊕}
3. Compute Earth surface δt/t from Phase 3 cumulative-potential approach
4. Test F-γ-5-B: δt/t within factor 2 of 7×10⁻¹⁰

### §3.2 — Derivation strategy

**Path A — Local density (naive):**
R_s³ = 3M / (4π·n_critical·m_eff) → R_s ∝ M^(1/3)
**Expected FAIL F-γ-5-A** (GR has R_s ∝ M linear)

**Path B — Cumulative potential (Phase 3 derived):**
- Use Phase 3 ⟨Φ⟩(r) Yukawa far-field ~ M/r
- δc/c_0 at r ∝ ⟨Φ⟩(r)/v ∝ M/(r·v·n_critical)
- δt/t ≈ δc/c_0
- R_s: condition c_eff = 0 → M/(R_s·v·n_critical) = 1 → R_s ∝ M ✓ (linear, GR-like!)

**Path B should match GR scaling.** Prefactor comparison with G/c² gives:
1/(v·n_critical) = G/c² ⇒ n_critical·v = c²/G

For n_critical = 1/ℓ_P³: v = c²·ℓ_P³/G = c²·(ℏG/c³)·G/G = ... let me compute:
- ℓ_P² = ℏG/c³, so ℓ_P³ = ℓ_P · ℏG/c³
- v = c²·ℓ_P·ℏG/(c³·G) = ℏ·ℓ_P/c

In Planck units (ℏ=c=ℓ_P=1): v = 1. So v ~ Planck mass M_P·c² (energy scale).

**Substantive FP tests:**
- T_P5_1: R_s Path A (local density mean field) for M_⊙, 1.4 M_⊙, M_⊕
- T_P5_2: R_s Path B (cumulative potential) — derive linear M scaling
- T_P5_3: G_eff from TGP params (v·n_critical = c²/G derivation)
- T_P5_4: F-γ-5-A verdict for each benchmark mass (factor 2 test)
- T_P5_5: Earth δt/t numerical evaluation
- T_P5_6: F-γ-5-B verdict (factor 2 around 7×10⁻¹⁰)

---

## §4 — DEC budget status

**Phase 3:** No DEC needed (uses Phase 2 c(n_local) + standard Phi field analysis).

**Phase 4:** No DEC needed (uses Phase 1 c(N) + γ-3 cosmological growth).

**Phase 5:** May invoke DEC 3 if n_critical scaling needs revision (Path A vs Path B reconciliation).

**Budget total after batch: 2-3 / 3 used.**

---

## §5 — Anti-Lakatos discipline (batch-wide)

- ✅ γ-3 + γ-3' verdicts STAY (NIE modified)
- ✅ Phase 1 + Phase 2 verdicts STAY (NIE modified)
- ✅ All F-γ-5 thresholds inherited unchanged
- ✅ F8 expected FAIL — predicted Phase 1 honestly; NIE retroactive rescue
- ✅ F-γ-5-A uncertain — Path A vs Path B comparison transparent
- ✅ §3.6.13 BINDING applied throughout

---

## §6 — Pre-derivation honest disposition (summary)

**Phase 3:** Expected SUCCESS — formal pair overlap gives 1/r potential → Newtonian gravity analog → G_eff identification.

**Phase 4:** Expected **FAIL_LITERAL F8** — c(N) saturates too fast to drive acceleration.

**Phase 5:**
- F-γ-5-A R_s: Path B (cumulative) expected PASS factor 2; Path A (naive) expected FAIL. Phase 5 will test both honestly.
- F-γ-5-B δt/t Earth: depends on Path B prefactor; expected PASS if v·n_critical = c²/G derivation works dimensionally.

**Aggregate γ-5 verdict prediction (pre-derivation, BEFORE batch execution):**
- F-γ-5-C ✓ PASS (Phase 1)
- F-γ-5-D ✓ PASS (Phase 2)
- Gravity synthesis: PASS expected (Phase 3)
- F8 re-test: FAIL expected (Phase 4)
- F-γ-5-A R_s: PASS via Path B expected (Phase 5)
- F-γ-5-B δt/t: PASS via Path B expected (Phase 5)

**Predicted aggregate: A- or B+ (mixed — F-γ-5 GR predictions PASS but F8 still FAIL)**

---

**END OF PHASE 3+4+5 BATCH PLAN — LOCKED 2026-05-24**
