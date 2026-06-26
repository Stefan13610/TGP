---
title: "Phase 2 — F-P25-B: condensed BVP verification (NS-NS) — FAIL_NEGATIVE"
type: phase_result
status: PHASE2_COMPLETE
phase: 2
cycle: op-mechanism-v-pattern25-extreme-envs-2026-06-10
created_date: 2026-06-11
authorization: "User 2026-06-11: 'ok działaj z P25' → recommended path (Phase 2 → FINAL)"
sympy_script: "[[./Phase2_bvp.py]]"
sympy_output: "[[./Phase2_bvp.txt]]"
sympy: "9/9 PASS; 0 hardcoded T_pass=True; all BVP runs converged (rms ~1e-11)"
falsifier_resolved: "F-P25-B = FAIL_NEGATIVE (NS-NS branch; pre-registered Phase 0 §3 thresholds applied mechanically); BH-BH branch NEGATIVE at gate (Phase 1)"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — F-P25-B: condensed BVP verification (NS-NS branch)

## §0 — Verdict at a glance

**F-P25-B = FAIL_NEGATIVE** (thresholds LOCKED Phase 0 §3, applied mechanically):

| Element | Result |
|---|---|
| Regression gate (mandatory) | **PASS** — T3 Phase 3 LOCKED .txt δψ = 6.833×10⁻⁸¹, rel dev 0.0004 |
| Nonlinear BVP anchor (M=0.01, σ=1) | **1.907×10⁻⁴** vs LOCKED T3 Phase 2 1.91×10⁻⁴ (rel dev 0.2%; rms 1.8×10⁻¹¹) |
| Amplitude ladder ×100 | log-log slope **1.00001**; ratio spread 0.005% — linear regime exact |
| Local formula validation (σ·m̃ ≈ 11.5) | central δψ vs (3/4)ρ̃(0): rel dev 2.2% — **S-ρ formula numerically confirmed by full nonlinear BVP** |
| NS-NS δψ_max (×2 contact bound) | **2.92×10⁻⁷⁹** |
| vs PARTIAL band 0.0385 | shortfall **77.1 orders of magnitude** |
| Numerical honesty clause | all runs converged; no tachyonic breakdown (far from T3 M_critical ≈ 15.80) |

## §1 — What was computed (numerical symmetry with T3, as recommended)

1. **EOM re-verification** (FP1): W = −V′/γ EXACT; m̃² = 4/3; δψ_critical = 2√3/9 EXACT.
2. **Mandatory regression gate** (FP2): T3 Phase 3 Branch A typical-LIGO recipe reproduced to 0.04%.
3. **Full nonlinear BVP** (T3 Phase 2 template: ψ″ + 2ψ′/r̃ + 2(ψ′)²/ψ + W(ψ) = −ρ̃; BC ψ′(0)=0,
   ψ(∞)=2/3; Pattern 2.1, NOT linearized): anchor at M=0.01/σ=1 matches LOCKED 1.91×10⁻⁴ (FP3).
4. **Amplitude ladder** M ∈ [10⁻⁴, 10⁻²] (FP4): response exactly linear (slope 1.00001) →
   linear extrapolation to NS amplitude is licensed (with Phase 1 FP14 no-bootstrap:
   W″(2/3)=0 EXACT, δm̃²/m̃² ~ 10⁻¹⁵⁷).
5. **Local-regime direct validation** (FP5): wide source (σ·m̃ ≈ 11.5) — full nonlinear BVP central
   response matches δψ = (3/4)ρ̃(0) to 2.2%. **This is the regime NS-NS sits in** (regime selector
   R_NS·m̃ ≈ 7×10³⁸ under Branch A) — the formula used for the verdict is BVP-validated, not assumed.
6. **NS-NS formal evaluation** (FP6): ρ_NS = 10¹⁸ kg/m³ (massive-NS central; Phase 1 FP13 verbatim);
   Branch A density unit = Planck density ≈ 5×10⁹⁶ kg/m³ → ρ̃ = 1.95×10⁻⁷⁹ →
   δψ_single = 1.46×10⁻⁷⁹; near-contact superposition bound ×2 → **δψ_max = 2.92×10⁻⁷⁹**.
7. **Degeneration audit** (FP7): ρ→0 ⇒ δψ→0; thresholds appear ONLY as comparison targets.

## §2 — Verdict mechanics

Pre-registered (Phase 0 §3, immutable): PASS_REALIZED ≥ 0.385 / PARTIAL_APPROACH ≥ 0.0385 /
FAIL_NEGATIVE < 0.0385. Computed δψ_max = 2.92×10⁻⁷⁹ → **FAIL_NEGATIVE**, shortfall 77.1 orders.
No interpretation discretion exists or was used.

**Why the verdict is robust:** to reach even the PARTIAL band, local density would need to be
~10⁷⁷× Planck-density-units higher, i.e. ρ ~ 10⁹⁵ kg/m³ — within ~1.5 OOM of the Planck density
itself. No astrophysical matter-bearing source approaches this; the gap is not closable by
modeling refinements (geometry, EOS, contact dynamics give O(1) factors).

## §3 — Anti-Lakatos verification (Phase 2)

- ✓ 0/9 hardcoded T_pass; verdict computed from LOCKED thresholds
- ✓ Branch A IMMUTABLE (no B/C switching despite NEGATIVE direction)
- ✓ No threshold loosening; no post-hoc source classes
- ✓ Numerical honesty clause honored: initial solver non-convergence (tol=10⁻¹⁰, uniform mesh)
  was FIXED (tol=10⁻⁸, geometric+linear mesh, r_min=10⁻³) rather than gate-loosened; final runs
  rms ~10⁻¹¹, all converged — fix documented here per §3 clause
- ✓ Predecessor verdicts PRESERVED (T3, mPhi-verification "mechanism (iii) FAILS typical LIGO", P6 R5)
- ✓ Phase 1 PARTIAL_SOURCE_NS_ONLY unchanged; this phase records the NS-NS formal verdict it deferred

## §4 — Continuation (per Phase 0 §1.2 decision structure)

- F-P25-C: **NOT_APPLICABLE** (conditional gate; F-P25-B = FAIL_NEGATIVE) — recorded at FINAL
- F-P25-D: **NEGATIVE** (mechanical: A PARTIAL_SOURCE_NS_ONLY + B FAIL_NEGATIVE both branches)
- → Phase FINAL closure; NO PR-021 (forbidden move); mechanism v routes to candidate (c)
  framework extension per enumeration handoff
