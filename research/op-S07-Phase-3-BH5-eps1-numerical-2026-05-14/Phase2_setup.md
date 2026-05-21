---
title: "Phase 2 Setup — ε.1 photon ring symbolic family marker mapping"
date: 2026-05-14
type: phase-setup
parent: "[[./README.md]]"
phase: 2
phase_focus: "ε.1 photon ring quadrant shift δε_ph²/ε_ph²_GR symbolic mapping per S07 family"
status: 🟢 ACTIVE — Phase 2 sympy in progress (2026-05-14 sesja P2-eps1)
predecessors:
  - "[[./Phase1_results.md]] (12/12 PASS; BH5 KEY DERIVATION δω/ω = κ_geom·d²f/dψ²/2·Δψ²)"
  - "[[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] §3.4 (ε.1 quadrant channel coefficients {0, β_q/9, α²/18})"
  - "[[../op-eht/]] (M9.1'' +14.6% photon ring shadow shift data point)"
sympy_target: "12 tests; ≥9 FP (≥75%); ≤3 LIT; 0 hardcoded `T_pass = True`"
substance_target: "Symbolic derivation δε_ph²/ε_ph²_GR (quad channel) = κ_ε · d²f/dψ²(ψ_0)/2 z κ_ε=1/9 (S07-reset Phase 2 inheritance) + cross-channel BH5/ε.1 ratio invariant test for trans family"
tags:
  - phase-2
  - phase-setup
  - eps1-photon-ring
  - quad-channel-family-distinct
  - cross-channel-ratio-invariant
---

# Phase 2 Setup — ε.1 photon ring symbolic family marker mapping

> **Phase scope:** Symbolic derivation ε.1 photon ring quadrant shift `δε_ph²/ε_ph²_GR` (quadratic family-discriminator channel) per S07 family {polynomial, quadratic, transcendental}, mapping family marker `d²f/dψ²(ψ_0)` na obserwabile via form-meaning analog Cunha-Herdeiro photon ring formulas. Critical Phase 2 deliverable: cross-channel BH5/ε.1 ratio invariant test for transcendental family.

## §0 — Phase 2 ASK-RULE pre-flight (per TGP_NATIVE_COMPUTATIONAL_PATTERNS §1.1)

### §0.1 — Trigger A (No TGP analogy?) → ✅ PASS

S07-reset Phase 2 ε.1 quadrant channel coefficients {0, β_q/9, α²/18} są LIVE TGP-native LOCK per family marker derivation. Pattern 2.2 form-meaning analog dla photon orbit modification z g_eff[Φ_eq(BH)] established Phase 1.

### §0.2 — Trigger B (Predecessor analogy outdated?) → ✅ PASS z conditional citation

- **op-S07-reset Phase 2 (LIVE):** ε.1 quad channel coefficients {0, β_q/9, α²/18} z derivation w polynomial/quadratic/transcendental basis (kategoria a — LIVE LOCK)
- **op-eht (M9.1'' +14.6% photon ring shadow shift):** observational data point (kategoria not-applicable — observational)
- **op-eps-photon-ring/Phase3 E2-E6 LIVE:** ε_ph² = 23²/137² = 529/18769 fine-structure decomposition (kategoria b — F4 chain anchor; conditional citation w T7 cross-check)

### §0.3 — Trigger C (Just reproducing standard physics?) → ⚠️ HIGH RISK identified + MITIGATED

**Standard physics default (anti-pattern):** Cunha-Herdeiro 2018 photon ring formulas — modified BH metric (Schwarzschild + scalar hair) → photon orbit equation Solver → r_ph + b_crit shift z scalar coupling.

**TGP-meaning per Pattern 2.2:** photon (test particle) moves on g_eff[Φ_eq(BH)] geodesic. Photon ring r_ph determined by extremum of effective potential w g_eff. Family marker d²f/dψ²(ψ_0) modifies SECOND-ORDER curvature near photon orbit. Per Pattern 2.5: κ_ε(BH-environment) ≠ κ_cosmological.

**Key form-meaning split (per Phase 1 §0.3 analog):**

| BD-form (Cunha-Herdeiro) | TGP-meaning (this cycle) |
|---|---|
| `δε_ph²/ε_ph²_GR = κ_CH · g_scalar` (scalar exchange at r_ph) | `δε_ph²/ε_ph²_GR = κ_ε · d²f/dψ²(ψ_0)/2` (g_eff curvature at photon orbit) |
| `g_scalar` = scalar coupling Lagrangian vertex | `d²f/dψ²(ψ_0)/2` from S07 family marker (LIVE Phase 2) |
| `κ_CH` = analytic factor z scalar field BC at photon ring | `κ_ε = 1/9` from photon orbit geometry r_ph = 3M (no Δψ² factor; single radius sampling) |

**Critical scope clarification:** Phase 2 derives QUADRATIC family-discriminator channel ONLY. The TOTAL photon ring shift (including LINEAR channel α/3 shared across families) is derived separately in S07-reset Phase 2 / parent emergent-metric cycles. M9.1'' anchor +14.6% z op-eht is dominated by LINEAR channel; QUADRATIC contribution z this Phase = 4/9 (M9.1'' specific) — represents family-discriminator small-add component, NOT total shift.

### §0.4 — Trigger D (Inheriting LOCK suspect?) → ✅ PASS

- S07-reset Phase 2 ε.1 quad coefficients (clean LOCK; kategoria a)
- κ_ε = 1/9 (S07-reset Phase 2 derivation inheritance; kategoria a)
- ε_ph² = 23²/137² F4 chain (op-eps-photon-ring LIVE; kategoria b z conditional citation w T7)
- ngEHT σ_θ ~5 μas (literature; kategoria not-applicable)

**Brak suspect inheritance.**

### §0.5 — Phase 2 sympy substance plan (per README §0.5 + Phase0_balance §4)

| Test | Type | Substance question | Classification |
|---|---|---|---|
| T1 | FP | Schwarzschild photon orbit baseline: ε_ph²_GR = 1/27 (r_ph=3M, b_crit=3√3·M); symbolic from Schwarzschild geodesic equation | **FIRST_PRINCIPLES** |
| T2 | FP | f(ψ) Taylor expansion at ψ_ph per family (analog Phase 1 T2 dla photon orbit) | **FIRST_PRINCIPLES** |
| T3 | FP | δε_ph²/ε_ph²_GR (quad channel) = κ_ε · d²f/dψ²(ψ_0)/2 z κ_ε=1/9 inherited derivation; linear in d²f | **FIRST_PRINCIPLES** |
| T4 | FP | Polynomial channel: d²f/dψ²=0 → δε_ph²/ε_ph²_GR (quad) = 0 EXACT (GR limit recovered for quad channel) | **FIRST_PRINCIPLES** |
| T5 | FP | Quadratic channel: d²f/dψ²=2β_q → δε_ph²/ε_ph²_GR = (1/9)·β_q (matches S07-reset Phase 2 derivation) | **FIRST_PRINCIPLES** |
| T6 | FP | Transcendental channel: d²f/dψ²=α² → δε_ph²/ε_ph²_GR = α²/18 (matches S07-reset Phase 2 derivation) | **FIRST_PRINCIPLES** |
| T7 | FP | M9.1'' anchor: d²f_M911/dψ²(1)=8 → quad channel = 4/9 ≈ 44.4%; clarify NIE total +14.6% (linear-dominated); quad-channel-only family discriminator | **FIRST_PRINCIPLES** |
| T8 | FP | c_0·κ_σ=4/3 LOCK preservation (Path 2 anchor); photon ring quad channel independent at leading order | **FIRST_PRINCIPLES** |
| T9 | FP | ASK-RULE Trigger C resolution: κ_ε=1/9 NIE BD coupling (specific photon-orbit geometric factor); symbolic distinguishability | **FIRST_PRINCIPLES** |
| T10 | LIT | ngEHT angular resolution σ_θ ~5 μas; σ_ε_ph² translation per source | **LITERATURE_ANCHORED** |
| T11 | LIT | 10-SMBH stack discriminability (Sgr A*, M87*, NGC1277, etc.) per BH4 LIVE source | **LITERATURE_ANCHORED** |
| T12 | FP | **Cross-channel ratio invariant**: BH5(trans)/ε.1(trans) = κ_geom·9·(Δψ_ringdown)² INDEPENDENT of α — geometric only | **FIRST_PRINCIPLES** |

**Target:** ≥9 FP (75%) — strong: T1-T9 + T12 = 10 FP guaranteed; T10+T11 = LIT.

**DEC structural declarations (separate, NOT counted in 12/12):**
- DEC-3: S05 single-Φ axiom preserved across ε.1 channel (no Φ_2)
- DEC-4: ax:metric-coupling preserved (universal g_eff coupling at photon orbit)

**Phase 2 KEY NEW DERIVATION (vs Phase 1):**

T12 cross-channel ratio invariant for trans family:
```
δω_QNM/ω_GR (BH5, trans family)     κ_geom · α²/2 · (Δψ_ringdown)²
─────────────────────────────── = ─────────────────────────────── = 9 · κ_geom · (Δψ_ringdown)²
δε_ph²/ε_ph²_GR (ε.1, trans family)            α²/18
```

**Result: ratio depends ONLY on geometric factors (κ_geom, Δψ_ringdown), NIE on α family parameter.** This is a PRE-OBSERVATIONAL test — measuring both BH5 and ε.1 simultaneously gives κ_geom·(Δψ)² extraction independent of family-marker amplitude. **Substantively novel cross-channel discriminator.**

## §1 — Risk register Phase-specific

| # | Risk | Probability | Mitigation |
|---|---|---|---|
| **P2.R1** | T7 M9.1'' anchor consistency for ε.1 quad channel: 4/9 ≈ 44% NIE matches op-eht +14.6% (different observable; quad-only vs total) | MEDIUM (interpretive risk; need explicit annotation) | T7 explicit annotation: quad channel derived ONLY; total shadow shift dominated by linear α/3 channel (S07-reset Phase 2 / parent cycles) |
| **P2.R2** | T3 derivation κ_ε=1/9 may need geometric re-derivation if S07-reset Phase 2 inheritance is non-trivial | LOW (S07-reset Phase 2 verified) | Inherit κ_ε=1/9 z S07-reset Phase 2 LIVE LOCK; explicit cite w T3 doc |
| **P2.R3** | T12 cross-channel ratio test may be tautological if both BH5 and ε.1 channels share α² factor | LOW (substantive: ratio cancels α², leaves pure geometric) | T12 explicit symbolic: `(κ_geom·α²/2·Δψ²) / (α²/18) = 9·κ_geom·Δψ²` — α² CANCELS, ratio = pure geometric |
| **P2.R4** | LIT classification creep: T10 + T11 over-budget (≤3 LIT target) | LOW | T7 borderline FP/LIT — deliberately counted FP via symbolic computation `sp.diff(f_M911, ψ, 2).subs(ψ, 1) == 8` (NIE pure LIT lookup) |

## §2 — Phase 2 entry gate confirmation

**All gates from Phase 1 closure PASS:**
1. ✅ Phase 1 sympy 12/12 PASS, 10 FP (83.3% > 75% binding)
2. ✅ Three-layer L1/L2/L3 explicit (Phase1_results §2)
3. ✅ Cross-cycle inheritance preserved 7/7
4. ✅ Anti-Lakatos PR-012 6/6 sub-checks PASS
5. ✅ ASK-RULE 4/4 Triggers PASS
6. ✅ User authorization "Authorize Phase 2 ε.1 same session (Opcja A continuation)" 2026-05-14 sesja P2-eps1
7. ✅ WIP slot 1/5 remains occupied (cycle ACTIVE)

**Phase 2 sympy entry: AUTHORIZED.**

## §3 — Sign-off + next deliverables

**Phase 2 setup committed:** 2026-05-14 sesja P2-eps1 (Claudian).

**Status:** 🟢 ACTIVE — Phase 2 sympy in progress.

**Next deliverables:**
- `Phase2_sympy.py` — 12 tests target ≥9 FP / ≤3 LIT / 0 hardcoded
- `Phase2_sympy.txt` — output saved (PYTHONIOENCODING=utf-8)
- `Phase2_results.md` — three-layer L1/L2/L3 + per-family table + cross-channel ratio invariant + cumulative cycle metrics

**Cross-references:**
- [[./README.md]] (BINDING contract; PR-012 LOCKED-PHASE-1-COMPLETE)
- [[./Phase1_results.md]] (BH5 KEY DERIVATION; cross-cycle 7/7 PASSED)
- [[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] §3.4 (ε.1 quadrant coefficients inheritance)
- [[../op-eht/]] (M9.1'' +14.6% data point; total shadow shift, NIE quad-only)
- [[../op-eps-photon-ring/Phase3_results.md]] E2-E6 (ε_ph² fine structure inheritance reference)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1+§2.2+§4 (ASK-RULE + Pattern 2.2)
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-012 LOCKED-PHASE-1-COMPLETE → LOCKED-PHASE-2-COMPLETE
