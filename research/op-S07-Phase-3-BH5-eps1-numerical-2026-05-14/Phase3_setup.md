---
title: "Phase 3 Setup — numerical projections + family discriminability matrix + 4-way M9.1'' anchor"
date: 2026-05-14
type: phase-setup
parent: "[[./README.md]]"
phase: 3
phase_focus: "Numerical projections per family at fiducials + LIGO-O5/CE/ngEHT/LISA σ family discriminability + cross-cycle 4-way anchor"
status: 🟢 ACTIVE — Phase 3 sympy in progress (2026-05-14 sesja P3-numerical)
predecessors:
  - "[[./Phase1_results.md]] (BH5 KEY DERIVATION)"
  - "[[./Phase2_results.md]] (ε.1 KEY DERIVATION + cross-channel ratio invariant)"
sympy_target: "10 tests; ≥7 FP (≥70%); ≤3 LIT; 0 hardcoded"
substance_target: "Family discriminability matrix per detector; 4-way anchor matrix at α=-4 (BH5 + ε.1 + S07-reset α/3 + emergent-metric c_0·κ_σ=4/3); cross-channel coupled bound; LISA 2035+ projection"
tags:
  - phase-3
  - numerical-projections
  - family-discriminability-matrix
  - 4-way-anchor-matrix
  - cross-cycle-coherence
---

# Phase 3 Setup — numerical projections + family discriminability matrix

> **Phase scope:** Numerical projections per S07 family at fiducial parameter values + LIGO-O5/CE/ngEHT/LISA σ thresholds → family discriminability matrix + cross-channel coupled bound + 4-way M9.1'' anchor consistency matrix (BH5 + ε.1 + S07-reset α/3 + emergent-metric c_0·κ_σ=4/3).

## §0 — Phase 3 ASK-RULE pre-flight

### §0.1 — Trigger A → ✅ PASS

Numerical projections z fiducial α/β_q values (inherited z S07-reset PR-010 recovery region + Phase 2 pre-bounded β_q range); detector sensitivity envelopes z literature LIGO-O5/CE/ngEHT/LISA. **Brak gap'u.**

### §0.2 — Trigger B → ✅ PASS

- Phase 1+2 derivations LIVE (kategoria a)
- LIGO-O5/CE/ngEHT/LISA literature (kategoria n/a — observational sensitivity envelopes)

### §0.3 — Trigger C → ✅ PASS (no new derivation in standard physics tools)

Phase 3 takes Phase 1+2 SYMBOLIC formulas + numerical evaluation at fiducials. Brak nowych dérivations w standard physics tools. ASK-RULE Trigger C zaadresowany w Phase 1+2 §0.3 form-meaning split — preserved throughout Phase 3.

### §0.4 — Trigger D → ✅ PASS (all LOCKs preserved through Phase 1+2)

### §0.5 — Phase 3 sympy substance plan (10 tests)

| Test | Type | Substance question | Classification |
|---|---|---|---|
| T1 | FP | BH5 family discriminability LIGO-O5: per-pair Δω/σ_BH5_O5 (poly-vs-quad, poly-vs-trans, quad-vs-trans) at fiducials α/β_q∈{-σ, 0, +σ} | **FIRST_PRINCIPLES** |
| T2 | FP | BH5 family discriminability Cosmic Explorer: same matrix z CE σ_BH5_CE = 0.025% (factor ×20 improvement vs LIGO-O5 stack) | **FIRST_PRINCIPLES** |
| T3 | FP | ε.1 family discriminability ngEHT: per-pair Δε_ph²/σ_ε_ngEHT z 10-SMBH stack σ_stack ≈ 6.3% (Phase 2 T11 inheritance) | **FIRST_PRINCIPLES** |
| T4 | FP | Cross-channel coupled BH5+ε.1 bound: simultaneous fit reduces α uncertainty via cross-channel ratio invariant (T12 Phase 2) | **FIRST_PRINCIPLES** |
| T5 | FP | LISA 2035+ EMRI ringdown ~0.1% σ_BH5_LISA projection: family discriminability ×4 vs LIGO-O5 stack | **FIRST_PRINCIPLES** |
| T6 | FP | **4-way M9.1'' anchor matrix at α=-4** (effective): BH5 trans channel + ε.1 trans channel + S07-reset α/3 ppE + emergent-metric c_0·κ_σ=4/3 simultaneously consistent | **FIRST_PRINCIPLES** |
| T7 | FP | Recovery region α∈[-0.832, 0.832] PASSING test per channel: each fiducial α value within recovery PASSES current GWTC-3 1σ for ppE channel; family-channel fiducials pre-bound β_q | **FIRST_PRINCIPLES** |
| T8 | LIT | LIGO-O5 PSD-coupled SNR_BH5 at fiducial M_BH=30 M_⊙ (GW150914-class): SNR=20 single-event per BH5 LIVE T3.2 | **LITERATURE_ANCHORED** |
| T9 | LIT | Photon orbit numerical at Sgr A* M=4.1·10⁶ M_⊙ θ_shadow=50 μas (per Phase 2 T10 ngEHT inheritance) | **LITERATURE_ANCHORED** |
| T10 | FP | Pre-observational discriminability sigmas per family pair: minimum signal threshold per detector per channel za 5σ family discrimination | **FIRST_PRINCIPLES** |

**Target:** ≥7 FP (70%) — wzór: T1-T7 + T10 = 8 FP guaranteed; T8+T9 = LIT (2). FP% ≈ 80% expected.

**DEC structural declarations (separate, NOT counted in 10/10):**
- DEC-5: Numerical fiducials są post-PE values (z S07-reset PR-010 recovery), NIE first-principles
- DEC-6: Pre-observational discrimination ceiling A− (full A reserved dla actual detection data)

## §1 — Risk register Phase 3

| # | Risk | Probability | Mitigation |
|---|---|---|---|
| **P3.R1** | Numerical fiducials may shift signal magnitudes outside detector sensitivity per family | LOW (recovery region α∈[-0.832, 0.832] inherited; signals scale linearly/quadratically) | Per-fiducial check w T1-T5 explicit |
| **P3.R2** | 4-way anchor matrix at α=-4 may NOT consistently match (4-way over-constrained) | LOW (each anchor independently derived in different cycles; convergence non-trivial test) | Symbolic computation per anchor; verify each independently; report tension if any |
| **P3.R3** | LISA EMRI projection extrapolation z conceptual paper | MEDIUM (LISA timeline ~2035+; speculative numerics) | Honest LIT classification; explicit annotation forecasts |

## §2 — Phase 3 entry gate confirmation

**All gates from Phase 2 closure PASS:**
1. ✅ Cumulative cycle 24/24 PASS, 20 FP (83.3%)
2. ✅ Three-layer L1/L2/L3 explicit (Phase 1+2 results)
3. ✅ Cross-cycle inheritance preserved 9/9
4. ✅ Anti-Lakatos PR-012 6/6 sub-checks PASS
5. ✅ ASK-RULE 4/4 Triggers PASS
6. ✅ Cross-channel ratio invariant T12 SUBSTANTIVELY NOVEL discriminator
7. ✅ User authorization "Authorize Phase 3 numerical + Phase FINAL combined same session (Opcja A heroic)"

**Phase 3 sympy entry: AUTHORIZED.**

## §3 — Sign-off + next deliverables

**Phase 3 setup committed:** 2026-05-14 sesja P3-numerical (Claudian).

**Status:** 🟢 ACTIVE — Phase 3 sympy in progress.

**Next deliverables (Phase 3 + Phase FINAL combined per Opcja A):**
- `Phase3_sympy.py` — 10 tests
- `Phase3_sympy.txt` — output
- `Phase3_results.md` — three-layer + family discriminability matrix + 4-way anchor + verdict draft
- `Phase_FINAL_close.md` — closure ceremony A−
- README.md update (closed-resolved; claim_status A−)
- PR-012 → LOCKED-PENDING-DATA
- STATE.md final closure record
- PREDICTIONS_REGISTRY entry

**Cross-references:** [[./Phase1_results.md]] + [[./Phase2_results.md]] + [[./README.md]]
