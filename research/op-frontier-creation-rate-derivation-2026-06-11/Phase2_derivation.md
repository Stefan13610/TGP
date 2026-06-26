---
title: "Phase 2 — F-FCR-C target #2: bulk transport DERIVED — C-form selected; B1 prediction 10^3.25 (PASS_BAND, 0.25 dex from observed)"
type: phase_result
status: PHASE2_COMPLETE
phase: 2
cycle: op-frontier-creation-rate-derivation-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'ok działaj z #2'"
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
sympy: "8/8 PASS; 0 hardcoded T_pass=True"
falsifier_resolved: "F-FCR-C: PARTIAL_concept_mismatch → C-DERIVED_CONDITIONAL (assumptions A-i/A-ii/A-iii declared; caveat C-2); F-FCR-D remains STRUCTURAL_CONDITIONAL; NO PR-022"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — bulk transport of frontier-created matter: derivation

## §0 — Verdict at a glance

**F-FCR-C upgraded: PARTIAL_concept_mismatch → C-DERIVED_CONDITIONAL.** The bulk equation form
is no longer a 3-way proxy menu — it is **derived** (under declared assumptions) and it is NONE
of the original C2a/b/c forms:

```
C-DERIVED:  δ″ + (4/3τ)δ′ − (ε_G/τ²)δ = 0        [H_m = 2/(3t); source-free bulk]
```

| Cell | p₊ | log₁₀G (τ: 1/1091 → 1) | Verdict |
|---|---|---|---|
| **B1 (ε = 3/2 EXACT)** | **(√55−1)/6 = 1.06937 EXACT** | **3.249** | **PASS_BAND** ⭐ (0.25 dex from observed 3.0) |
| B2 (ε = 0.465) | 0.535 | 1.626 | PARTIAL_BAND |

## §1 — Derivation chain (each step forced; assumptions declared, not hidden)

1. **A-i (concept paper LOCKED-claim, §6):** bulk spontaneous creation BLOCKED ("Bulk saturated,
   blocked — E2 property") ⇒ bulk continuity is **source-free exactly**. This EXCLUDES the
   uniform-in-bulk creation pictures behind C2a/C2b (and moots C2c's δ_S question).
2. **A-ii (homogeneity, imposed as consistency requirement; concept §6 CMB row):**
   ρ(x,t) = ρ̄(t) ∝ t⁻² (Phase 1 B-skeleton) in the bulk.
3. **FP1:** source-free continuity + homogeneity ⇒ **∇·u = −ρ̄̇/ρ̄ = 2/t** (forced).
4. **A-iii (isotropy):** u = f(t)·x ⇒ f = 2/(3t) ⇒ **matter-flow scale factor a_m ∝ t^(2/3)**
   (FP2) — DISTINCT from the space frontier R = c·t. Matter lags the frontier; the gap is filled
   by frontier creation (qualitatively consistent with A-i). Photon-side epoch mapping stays
   γ-3-native (1+z = t₀/t) — photons ride generated space; matter rides the derived flow.
5. **FP4:** perturbations around the derived flow (standard continuity+Euler+Poisson, NO source
   terms — they are excluded by A-i, not assumed away) ⇒ C-DERIVED form above.
   **Validation: EdS limit (ε = 2/3) returns p = 2/3 EXACT** — the machinery reproduces the
   textbook Einstein–de Sitter growth law as a special case.
6. **FP5-6:** B1 ⇒ p = (√55−1)/6 EXACT; log₁₀G = 3.2485 (numeric cross-check rel dev 3×10⁻¹⁴).

## §2 — Honest caveats (carried into aggregate)

- **C-2 (FP3):** background Euler residual Δ_bulk = |3ε−2|/4 = 0.625 (B1) — bounded,
  time-independent (per R17 lemma: no-runaway class) — but its balance requires an O(1)·H_m²·x
  substrate force. TGP-native in spirit (expansion is substrate-driven, γ-3) yet NOT derived from
  the action here.
- **A-ii** is imposed (consistency with CMB isotropy), not derived from frontier dynamics.
- **B1** remains STRUCTURAL_POSTULATE (zero-energy condition — missing piece #1, unchanged).
- δ growth target G_obs = 10³ never used in any derivation (FP7 guard; free symbols = {ε}).

## §3 — F-FCR-D aggregate (post Phase 2): STRUCTURAL_CONDITIONAL (unchanged class)

**NO PR-022** (forbidden move #5 — B1 postulate + C-2/A-ii caveats outstanding). The prediction
is now **one derivation away** from PR-022: the zero-energy condition (missing piece #1). If
derived, the chain closes: parameter-free log₁₀G = [(√55−1)/6]·log₁₀(1091) = 3.249 vs observed
3.0 as a genuine pre-registered prediction.

**Updated missing-piece list:** (1) zero-energy condition ⭐ NEXT; (2) ~~bulk transport~~ **DONE
(THIS phase, conditional)**; (3) A2 Φ→matter bridge (subsumable under #1 energetics); (+) C-2
substrate-balance + A-ii derivation (both natural corollaries of a frontier-energetics cycle).

## §4 — Anti-Lakatos (Phase 2): COMPLIANT ✓

0/8 hardcoded; derived form supersedes proxy menu BY DERIVATION (not selection-by-result —
the chain A-i→FP1→FP2 contains no reference to growth values; FP7 circularity guard);
EdS validation independent; predecessor verdicts (incl. Phase 1 matrix as state-of-record)
PRESERVED; bands LOCKED Phase 0 unchanged; PR-022 withheld.
