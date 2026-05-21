---
title: "Phase 0 — Balance red-giant tension analysis"
date: 2026-05-17
parent: "[[./README.md]]"
type: phase-zero-balance
phase: 0
status: 🟢 ACTIVE
sympy_substance_ratio: "6 FP / 1 LIT / 1 DEC = 75% FP"
hardcoded_T_pass: 0
---

# Phase 0 — Balance red-giant tension cyklu

## §1 — Inputs

| Element | Source | Status |
|---|---|---|
| μ_ν^TGP central = 3.55·10⁻¹² μ_B | Cycle 3 spinor B | LIVE prediction |
| m_X = 60 MeV (anchor) | L06 | NUMERICAL ANCHOR |
| Suppression heuristic (L_kink/λ_C)² | Cycle 3 T6 | LIVE placeholder |
| Capozzi-Raffelt 2020: μ_ν < 1.2·10⁻¹² μ_B (TRGB 2σ) | arXiv:2007.03694 | LIT bound |
| Raffelt 1990: μ_ν < 3·10⁻¹² μ_B (classical) | Phys Rep 198, 1 | LIT bound (conservative) |
| Viaux+2013 (M5): μ_ν < 4.5·10⁻¹² μ_B (95% CL) | A&A 558 A12 | LIT bound |

## §2 — Outputs

**Phase 1:**
- Best red-giant bound z literature (T1)
- Statistical 1σ/2σ interpretation (T2)
- TGP prediction propagated z m_X uncertainty (T3-T4)
- Suppression form sensitivity (T5)
- Tension assessment in σ units (T6)
- Falsifiability re-assessment (T7)

**Phase FINAL:**
- Verdict: TENSION REAL / MARGINAL / NONE
- Implications for cycle 3 prediction status

## §3 — Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| R1 Stellar model systematics | medium | Use Capozzi-Raffelt 2σ (already includes systematics) |
| R2 Suppression placeholder | medium | T5 sensitivity scan n ∈ [1,3] |
| R3 m_X anchor 1.7x range | medium | T4 sensitivity scan |
| R4 Interpretation choices | low | Honest CI propagation |

## §4 — 8/8 gate

- [x] G1-G8 all checked ✓

**Verdict:** 🟢 **8/8 PASS.**

---

**Sign-off:** Claudian @ 2026-05-17 (4th cycle sesji).
