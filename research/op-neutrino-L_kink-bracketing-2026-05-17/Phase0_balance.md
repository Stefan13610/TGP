---
title: "Phase 0 — Balance L_kink bracketing"
date: 2026-05-17
parent: "[[./README.md]]"
type: phase-zero-balance
phase: 0
status: 🟢 ACTIVE
sympy_substance_ratio: "6 FP / 1 LIT / 1 DEC = 75% FP"
hardcoded_T_pass: 0
---

# Phase 0 — Balance L_kink bracketing cyklu

## §1 — Inputs

| Element | Source | Status |
|---|---|---|
| β-task source S | Cycle 1 (2026-05-17) | LIVE LOCK |
| Spinor channel μ_spinor | Cycle 2 RP² (2026-05-17) | LIVE |
| m_X ≈ 60 MeV substrate anchor | L06 [[../op-L06-axion-mass-derivation-2026-05-16/]] | NUMERICAL ANCHOR (B+ partial) |
| m_eff = c_M·A_tail²·g_0^(e²/2) | why_n3 PHASE2 | LIVE |
| g_0_ν ≈ 0.22, A_tail_ν ≈ 6.1·10⁻⁴ | Exploration playground | LIVE calibration |
| m_ν = 0.1 eV (light ν typical) | PDG 2024 ν₃ scale | LIVE phenomenological |
| XENONnT < 6.3·10⁻¹² μ_B (2022) | PRL 129.161805 | LIT bound |
| SM Dirac ≈ 3·10⁻²⁰ μ_B | Marciano-Sanda 1977 | LIT reference |

## §2 — Outputs

**Phase 1:**
- L_kink range bracket z 3 scenarios
- μ_ν^TGP range dla scalar + spinor channels
- Falsifiability window check
- Honest classification: bracketing (NOT first-principles)

**Phase FINAL:**
- B+ partial verdict per pre-registered tree
- Empirical commitments + handoff

## §3 — Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| R1 m_X anchor status | medium | Honestly inherited; documented w T2 LIT |
| R2 m_eff vs m_ν observable | low | Use m_ν = 0.1 eV PDG (observable) |
| R3 Loop-conversion δθ→μ | medium | Heuristic Larmor mapping; quantitative deferred |
| R4 Range overlap | low (cosmetic) | Acceptable; honest uncertainty propagation |

## §4 — 8/8 gate

- [x] G1 contract written ✓
- [x] G2 pre-registered falsification ✓
- [x] G3 native-first methodology ✓
- [x] G4 pre-flight read ✓
- [x] G5 substance 75% FP, 0 hardcoded ✓
- [x] G6 risks (4) ✓
- [x] G7 decision tree explicit ✓
- [x] G8 downstream impact ✓

**Verdict:** 🟢 **8/8 PASS.**

---

**Sign-off:** Claudian @ 2026-05-17 (3rd cycle sesji).
