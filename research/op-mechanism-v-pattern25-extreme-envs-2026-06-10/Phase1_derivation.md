---
title: "Phase 1 — F-P25-A: TGP-native near-horizon source derivation — PARTIAL_SOURCE_NS_ONLY"
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-mechanism-v-pattern25-extreme-envs-2026-06-10
created_date: 2026-06-10
authorization: "User 2026-06-10: 'działaj z phase 1' → F-P25-A execution"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "15/15 PASS; 0 hardcoded T_pass=True"
falsifier_resolved: "F-P25-A = PARTIAL_SOURCE_NS_ONLY (pre-registered Phase 0 §3 criterion, verbatim)"
anti_lakatos_lock: PRESERVED
methodology_binding: "CALIBRATION_PROTOCOL §3.6 BINDING; TGP_NATIVE_COMPUTATIONAL_PATTERNS BINDING; m_Φ category per foundations §3.5.6.3 rule: m_Φ_intrinsic used for AMBIENT vacuum screening (explicit); ⟨Φ⟩_local shift is the DERIVED object"
---

# Phase 1 — F-P25-A: source derivation + scaling class

## §0 — Verdict at a glance

**F-P25-A = PARTIAL_SOURCE_NS_ONLY** (pre-registered criterion, Phase 0 §3, applied verbatim):

| Element | Result |
|---------|--------|
| BH-BH exterior native source | **= 0 identically** (ρ_matter = 0; no-hair analog at level-0) → **BH-BH branch NEGATIVE at the gate** |
| NS-NS near-contact source | well-defined: −q·ρ_matter (LOCKED M9.2/T3 convention) → continues to F-P25-B |
| Scaling class | **(S-ρ) density-type FORCED** by Branch A screening (regime selector σ̃·m̃ ≈ 2.1·10³⁹ ≫ 1) |
| (S-κ) compactness channel | **STRUCTURALLY EXCLUDED** — it is exactly the massless-limit (unscreened) response; under Branch A carries exp(−2.1·10³⁹) |
| NS-NS magnitude preview | δψ ≈ 1.46·10⁻⁷⁹ (~77 orders below factor-10 PARTIAL band) — **preview; F-P25-B verdict recorded in Phase 2** |
| Sympy | **15/15 PASS**; 0 hardcoded; regression gate FP7 rel dev 0.000 |

**The cycle's pre-declared bimodality (Phase 0 §1.2) is hereby RESOLVED on the negative branch
of the decision structure:** the TGP-native source under Branch A scales with LOCAL DENSITY,
not compactness. The "genuinely open" question closed cleanly at the derivation gate.

## §1 — Derivation chain (LOCKED-machinery-only; no new ingredients)

### §1.1 — Source form (inherited verbatim)

Level-0 LOCKED field equation with source (T3 Phase 2 §1.1 + M9.2 convention):

```
ψ'' + (2/r̃)ψ' + 2(ψ')²/ψ + W(ψ) = −q·ρ̃_matter ,   W(ψ) = (1/3)ψ(8−18ψ+9ψ²)
```

with physical coupling q = 8πG/c² (M9.2 setup, LOCKED). Re-verified in sympy: W(2/3) = 0;
m̃² = −W'(2/3) = 4/3; ψ_± = (6 ± 2√3)/9; δψ_critical = ψ_+ − 2/3 = **2√3/9 EXACT** ≈ 0.3849
(FP1-FP3). The level-0 machinery contains exactly ONE source channel: matter density ρ_matter.
No curvature coupling exists at level-0 (that would be candidate (c) framework extension —
out of scope per Phase 0 §1.3, forbidden to introduce here).

### §1.2 — Response kernel and the two regimes

Linearized around vacuum, the response is the Yukawa convolution
δψ(x) = q·∫ G(|x−x'|)·ρ̃(x') d³x', G(r) = exp(−m̃r)/(4πr), ∫G d³x = 1/m̃² (FP4). Two regimes:

- **R_source ≫ λ_C = 1/m̃ (local-screening regime):** δψ(x) = q·ρ̃(x)/m̃² = (3/4)·ρ̃(x) —
  the field tracks LOCAL density; it cannot integrate over the source (**S-ρ density-type**).
- **R_source ≪ λ_C (unscreened regime):** δψ(r) = q·M̃/(4πr̃) — the field integrates the
  whole mass; response ∝ compactness (**S-κ compactness-type**).

**Validation of the kernel (FP6):** the exact linear Yukawa-Gaussian central response at
M = 0.01, σ = 1 gives δψ = 1.907·10⁻⁴ vs LOCKED T3 Phase 2 full-nonlinear BVP 1.91·10⁻⁴ —
**ratio 1.00**. The kernel derivation is quantitatively exact against LOCKED numerics.

### §1.3 — Branch A regime selector (the decisive step)

Under Branch A (γ ~ M_Pl², IMMUTABLE per γ-cascade + PR-019): λ_C ~ Planck length. For ANY
astrophysical source (≥ km scale): σ̃·m̃ ≈ 2.14·10³⁹ ≫ 1 (FP9). **The local-screening regime
is FORCED — scaling class (S-ρ) is not a choice but a structural consequence.**

### §1.4 — (S-κ) exclusion + audit of the foundations "δψ ~ 0.3+" estimate

The compactness channel is exactly the **massless limit** of the kernel. Symbolically (FP10):

```
δψ(r)|_(m̃→0) = q·M/(4πr) = 2GM/(c²r)   →   δψ(R_s) = 1 EXACT
```

So IF the field were unscreened, near-horizon activation would be automatic (1 > 0.385) —
this identifies the provenance of the foundations §3.5.6 "extreme environments (δψ ~ 0.3+)"
figure: it is the **unscreened (massless-limit) intuition**. Under Branch A the channel carries
the Yukawa factor exp(−m̃·R̃_s) = exp(−2.11·10³⁹) ≈ 10^(−9.2·10³⁸) (FP11) and is
**STRUCTURALLY EXCLUDED**. Per Phase 0 §6 #15 the "0.3+" figure was treated as a test target,
NOT an input — **audit result: the estimate does NOT survive Branch A screening.**

### §1.5 — Source inventory per pre-declared class

**(i) BH-BH exterior:** ρ_matter = 0 identically → native source ≡ 0 (FP12). The gradient
self-source (∂δψ)² is O(δψ²) and vanishes on a zero perturbation; curvature coupling is absent
from level-0. **No-hair analog: BH-BH branch NEGATIVE at the gate.**

**(ii) NS-NS near-contact:** ρ_matter well-defined. Branch A density unit = Planck density
(m_Φ⁴ ≈ 5·10⁹⁶ kg/m³); massive-NS central density ~10¹⁸ kg/m³ → ρ̃ ≈ 1.9·10⁻⁷⁹ →
δψ ≈ (3/4)·ρ̃ ≈ **1.46·10⁻⁷⁹** (FP13; preview — F-P25-B verdict belongs to Phase 2).
Shortfall: ~77 orders vs the factor-10 PARTIAL band (0.0385); ~78 orders vs δψ_critical.

### §1.6 — Self-consistency (no screening bootstrap)

W''(2/3) = 0 EXACT (the vacuum is a symmetric point of W') — the leading correction to the
screening mass is SECOND order: δm̃²/m̃² = 6.75·δψ² ≈ 1.4·10⁻¹⁵⁷ (FP14). No positive-feedback
loop ("response lowers mass → larger response"): linearization is self-consistent to ~157 orders.
The T3 nonlinear amplification (×1.35 at δψ ~ 0.4) is real but operates only NEAR the threshold —
it cannot bridge 77 orders.

## §2 — Regression gate + INFORMATIONAL documentation flag

**FP7 (mandatory weak-field regression):** pipeline reproduces the LOCKED T3 Phase 3 executed
output δψ(typical LIGO, Branch A) = 6.833·10⁻⁸¹ with rel dev 0.000 (identical constants,
identical formula). Gate PASS — near-horizon claims recordable.

**FP8 (INFORMATIONAL; no predecessor modification):** T3 `Phase3_results.md` §0 quotes
δψ ≈ 1.74·10⁻⁷⁹ while the LOCKED executed output `Phase3_dimensional.txt` line 78 gives
**6.83·10⁻⁸¹** (factor 25.5 transcription slip in the markdown; the also-quoted "10⁻¹⁰⁴" is an
earlier rough figure). Verdict-irrelevant (all variants are ≥ 77 orders below threshold).
Analog of the S07-INT "documentation observation" — flagged for a future doc-cleanup pass;
predecessor verdict PRESERVED per Phase 0 §4.5.

## §3 — Anti-Lakatos verification (Phase 1)

| Check | Status |
|-------|--------|
| Branch A → B/C switching? | NO — Branch A used throughout; FP10-FP11 massless limit is an AUDIT of (S-κ), not a branch switch ✓ |
| Pattern 2.5 cited as realization evidence? | NO — §1.4 audits and REJECTS the unscreened estimate under Branch A ✓ |
| Thresholds loosened? | NO — 2√3/9 EXACT + factor-10 band, both in comparison gates only ✓ |
| Threshold constants leaked into source forms? | NO — FP15 contamination check clean ✓ |
| Post-hoc source classes? | NO — (i)/(ii) only ✓ |
| Source engineered toward activation? | NO — derivation FORCED (S-ρ); the activation-friendly (S-κ) was EXCLUDED, not adopted ✓ |
| Predecessor verdicts modified? | NO — mPhi "mechanism (iii) FAILS typical LIGO" untouched; FP8 informational only ✓ |
| m_Φ category rule (§3.5.6.3) honored? | YES — m_Φ_intrinsic explicitly for ambient screening; ⟨Φ⟩_local shift derived; self-consistency FP14 ✓ |
| Hardcoded T_pass=True? | 0/15 ✓ |
| New fundamental constants? | 0 ✓ |
| Circularity audit performed? | YES — FP15 (ρ→0, M→0 degeneration + threshold independence) ✓ |
| Honest-negative left reachable? | YES — and the BH-BH branch in fact RETURNED negative at the gate ✓ |

**Anti-Lakatos status: COMPLIANT ✓ (12/12).**

## §4 — Decision budget (Phase 1)

| Budget | Cap | Used |
|--------|-----|------|
| DEC | 3 | 0 |
| PARTIAL_compute | 1 | 0 |
| PARTIAL_concept_mismatch | R1 | 0 |
| Hardcoded T_pass=True | 0 | 0 |

## §5 — Implications + handoff to Phase 2

**F-P25-A = PARTIAL_SOURCE_NS_ONLY.** Per Phase 0 §3 + §9 decision points:

1. **BH-BH branch: NEGATIVE at the gate** (native source ≡ 0). Recorded; no BVP needed.
2. **NS-NS branch continues to F-P25-B (Phase 2)** — but the Phase 1 derivation strongly
   pre-conditions the outcome: (S-ρ) forced + ρ̃_NS ~ 10⁻⁷⁹ means the BVP scan is expected to
   return **FAIL_NEGATIVE by ~77 orders**. Honest framing: Phase 2 is now a **verification
   formality** (regression-gated BVP confirmation), not an open test.
3. **Recommended Phase 2 scope (condensed):** single focused BVP run at NS-NS near-contact
   parameters + the weak-field regression gate, recording F-P25-B formally; estimated 0.5 sesji
   (vs 1-2 in Phase 0 plan). Alternative per user choice: accept the analytical FP13 preview as
   sufficient and go directly to Phase FINAL with F-P25-B = FAIL_NEGATIVE recorded from the
   Phase 1 kernel computation (the kernel was validated to ratio 1.00 against LOCKED nonlinear
   numerics, FP6) + F-P25-C = NOT_APPLICABLE + F-P25-D = NEGATIVE.
4. **F-P25-C: anticipated NOT_APPLICABLE** (conditional gate; no activation to propagate).
5. **F-P25-D anticipated: NEGATIVE** — P6 R5 confirmed for extreme environments; mechanism v
   routes to candidate (c) framework extension (multi-cycle program per enumeration F-MECH-V-D).
   NO PR-021 append.

**Awaiting user decision:** "działaj Phase 2" (condensed BVP verification, recommended for
numerical closure symmetry with T3) **or** "Phase FINAL" (direct closure on Phase 1 analytics).
