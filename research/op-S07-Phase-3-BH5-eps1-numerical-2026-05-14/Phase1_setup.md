---
title: "Phase 1 Setup — BH5 QNM symbolic family marker mapping"
date: 2026-05-14
type: phase-setup
parent: "[[./README.md]]"
phase: 1
phase_focus: "BH5 QNM ringdown frequency shift δω_QNM/ω_GR symbolic mapping per S07 family"
status: 🟢 ACTIVE — Phase 1 sympy in progress (2026-05-14 sesja P1-bh5)
predecessors:
  - "[[./README.md]] (BINDING contract; PR-012 LOCKED)"
  - "[[./Phase0_balance.md]] (6/6 P-gate scope-PASS)"
  - "[[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] (family marker source)"
  - "[[../op-bh-alpha-threshold/Phase3_results.md]] T3.2 (BH5 anchor δf/f ~8-16% at α(ψ_ringdown=1.20)=0.1608)"
sympy_target: "12 tests; ≥9 FP (≥75%); ≤3 LIT; 0 hardcoded `T_pass = True`"
substance_target: "First-principles symbolic derivation δω_QNM/ω_GR = κ_QNM·d²f/dψ²(ψ_0)/2 per family + M9.1'' anchor consistency check"
tags:
  - phase-1
  - phase-setup
  - bh5-qnm-symbolic
  - family-marker-mapping
  - form-meaning-split
---

# Phase 1 Setup — BH5 QNM symbolic family marker mapping

> **Phase scope:** Symbolic derivation BH5 QNM ringdown frequency shift `δω_QNM/ω_GR` per S07 alternative f(ψ) family {polynomial, quadratic, transcendental}, mapping family marker `d²f/dψ²(ψ_0) = {0, 2β_q, α²}` na obserwabile via form-meaning analog Berti-Cardoso QNM perturbation spectroscopy.

## §0 — Phase 1 ASK-RULE pre-flight (per TGP_NATIVE_COMPUTATIONAL_PATTERNS §1.1)

### §0.1 — Trigger A (No TGP analogy?) → ✅ PASS

S07-reset Phase 2 family marker `d²f/dψ²(ψ_0)` jest TGP-native LIVE LOCK (kategoria a per M9_RESTRUCTURE §1.4). Pattern 2.2 form-meaning analog dla Newton/momentum-flux ↔ QNM/photon-orbit-perturbation już istnieje konceptualnie. **Brak gap'u.**

### §0.2 — Trigger B (Predecessor analogy outdated?) → ✅ PASS z conditional citation

- **op-S07-reset Phase FINAL A− (2026-05-13):** family marker source LIVE TGP-native LOCK (kategoria a)
- **op-bh-alpha-threshold/Phase3 T3.2 (BH5 LIVE):** używa `α(ψ) = α₀·(ψ-1)²` form (essentially second-order Taylor expansion around ψ=1) — kategoria (b) per M9_RESTRUCTURE §1.4 (M9.1''-specific anchor parameter; conditional na specific α₀=4.02 from Phase 2 of that cycle). Inheritance preserved jako anchor-consistency-check w T7; new derivation rebases na S07 family marker abstraction.
- **op-emergent-metric Phase 4 LOCKs:** c_0·κ_σ=4/3 (kategoria a) preserved; Path 2 anchor preserved (kategoria b)

### §0.3 — Trigger C (Just reproducing standard physics?) → ⚠️ HIGH RISK identified + MITIGATED

**Standard physics default that I would use (anti-pattern):** Berti-Cardoso 2009 QNM perturbation spectroscopy formulas:
```
ω_QNM_modified = ω_GR + δω(ζ)
δω(ζ) = perturbation z scalar coupling ζ na BH metric (Schwarzschild + scalar hair)
```
This is BD-mode framing — couples scalar field directly to graviton perturbation w sense effective Lagrangian.

**TGP-meaning per Pattern 2.2 (form-meaning split per README §0.1):**

W TGP, QNM modification pochodzi z structural change w `g_eff[Φ_eq(BH)]` jako funkcjonał of full nonlinear Φ-EOM solution Φ_eq(r) wokół BH. **Test particle (graviton perturbation) moves on geodezykach `g_eff` — NIE couples bezpośrednio do Φ-quanta.** Per Pattern 2.5: m_Φ_observable jest BH-environment-specific.

Family marker `d²f/dψ²(ψ_0)` controls **second-order Taylor expansion** of `f(ψ)` (and hence `g_eff[Φ_eq(BH)]`) around ψ_0 = ψ_horyzont. QNM shift jest sensitive do tej drugiej pochodnej (curvature near horizon).

**Mapping (per Pattern 2.2 form-meaning):**

| BD-form (Berti-Cardoso) | TGP-meaning (this cycle) |
|---|---|
| `δω_QNM/ω_GR = κ_BC · g_scalar_coupling` (scalar exchange) | `δω_QNM/ω_GR = κ_QNM · d²f/dψ²(ψ_0)/2` (g_eff[Φ_eq(BH)] structural curvature) |
| `g_scalar_coupling` from Lagrangian vertex | `d²f/dψ²(ψ_0)/2` from S07 f(ψ) family expansion |
| `κ_BC` = 1-loop integral over BH spacetime | `κ_QNM = (Δψ_ringdown)² · geometric_factor` from r_horizon → r_ringdown |

**ASK-RULE Trigger C resolution:** form preserved (Berti-Cardoso template), meaning replaced (TGP substrate g_eff[Φ_eq(BH)]). Każdy Phase 1 sympy test cite za §0.1+§0.3 dla anti-BD-drift documentation.

### §0.4 — Trigger D (Inheriting LOCK suspect?) → ✅ PASS

- S07-reset PR-010 LOCKED (clean LOCK; family marker LIVE; kategoria a)
- emergent-metric Phase 4 c_0·κ_σ=4/3 (clean LOCK; kategoria a)
- BH5 LIVE α(ψ_ringdown)=0.1608 (M9.1''-specific anchor; kategoria b z conditional citation w T7)
- LIGO-O5 PSD curve (literature, kategoria not-applicable — observational data)

**Brak suspect inheritance.**

### §0.5 — Phase 1 sympy substance plan (per README §0.5 + Phase0_balance §4)

| Test | Type | Substance question | Classification |
|---|---|---|---|
| T1 | FP | Schwarzschild QNM symbolic baseline (Berti-Cardoso form template; ω_QNM(M_BH) functional form) | **FIRST_PRINCIPLES** |
| T2 | FP | g_eff[Φ_eq(BH)] second-order Taylor expansion w r_horyzont; form-meaning per §0.3 | **FIRST_PRINCIPLES** |
| T3 | FP | δω_QNM/ω_GR derivation z perturbation g_eff^Schwarzschild + family marker; analytical κ_QNM(ψ_ringdown) | **FIRST_PRINCIPLES** |
| T4 | FP | Polynomial channel: d²f/dψ² = 0 → δω_QNM/ω_GR = 0 EXACT (GR limit recovered) | **FIRST_PRINCIPLES** |
| T5 | FP | Quadratic channel: d²f/dψ² = 2β_q → δω_QNM/ω_GR = κ_QNM·β_q | **FIRST_PRINCIPLES** |
| T6 | FP | Transcendental channel: d²f/dψ² = α² → δω_QNM/ω_GR = κ_QNM·α²/2 | **FIRST_PRINCIPLES** |
| T7 | FP/LIT | M9.1'' anchor consistency: α(ψ_ringdown=1.20)=0.1608 → δf/f within 8-16% LIVE inheritance | **FIRST_PRINCIPLES** + LIT cite |
| T8 | FP | c_0·κ_σ=4/3 LOCK preservation w κ_QNM derivation (Path 2 anchor) | **FIRST_PRINCIPLES** |
| T9 | FP | ASK-RULE Trigger C symbolic check: TGP κ_QNM functional dependence ≠ standard scalar coupling | **FIRST_PRINCIPLES** |
| T10 | LIT | LIGO O5 A+ PSD curve at 250 Hz (Hild+2010 ET-D / LIGO design) | **LITERATURE_ANCHORED** |
| T11 | LIT | Cosmic Explorer ~2030 stacked 5σ projection (CE concept paper) | **LITERATURE_ANCHORED** |
| T12 | FP | Pattern 2.5 environment-dependent: κ_QNM(BH-env) ≠ κ_cosmological (symbolic distinguishability) | **FIRST_PRINCIPLES** |

**Target:** ≥9 FP (75%) — wzór: T1-T6 + T8 + T9 + T12 = 9 FP guaranteed; T7 borderline FP/LIT depending on classification of M9.1'' anchor cite.

**DEC structural declarations (separate, NOT counted in 12/12 PASS):**
- DEC-1: S05 single-Φ axiom preserved across BH5 channel (no Φ_2, hidden field)
- DEC-2: ax:metric-coupling preserved (universal g_eff coupling, NIE Φ-quantum exchange)

## §1 — Risk register Phase-specific

| # | Risk | Probability | Mitigation |
|---|---|---|---|
| **P1.R1** | T3 symbolic derivation κ_QNM may not converge to clean closed form | LOW (geometric factor + Taylor expansion are standard symbolic operations) | Use sympy series expansion + simplify; if no closed form, parametrize as κ_QNM(Δψ²) + numerical evaluation w T7 |
| **P1.R2** | T7 M9.1'' anchor consistency may NOT match 8-16% range | MEDIUM (depends on specific κ_QNM derivation; existing α(ψ)=4.02·(ψ-1)² embeds α₀=4.02 specific to M9.1'' Phase 2 ansatz, NIE my new family marker abstraction) | Honest annotation jeśli range mismatch: declare TENSION; defer reconciliation to Phase 3 numerical or cross-cycle audit. Substance ceiling A− preserved. |
| **P1.R3** | T9 Trigger C symbolic check may be tautological (T_pass = True artifact) | LOW (symbolic check checks κ_QNM functional dependence on BH-radius geometric factor) | Explicit symbolic test: κ_QNM ≠ const (depends on Δψ²); failure mode = const independent of Δψ would indicate BD-style coupling |
| **P1.R4** | LIT classification creep: T7 + T10 + T11 over-budget (target ≤3 LIT) | LOW-MEDIUM | T7 borderline; if LIT, bump T10-T11 careful |

## §2 — Phase 1 entry gate confirmation

**All gates from Phase0_balance §6 PASS:**
1. ✅ README.md BINDING contract validator PASS (verified 2026-05-14 sesja P0)
2. ✅ Phase0_balance.md scope-PASS 6/6 P-gate
3. ✅ PR-012 LOCKED-PENDING-PHASE-1
4. ✅ User authorization "Authorize Phase 1 BH5 same session (Opcja A)" 2026-05-14 sesja P1
5. ✅ folder_status updated parking → active
6. ✅ WIP slot 1/5 occupied

**Phase 1 sympy entry: AUTHORIZED.**

## §3 — Sign-off + next deliverables

**Phase 1 setup committed:** 2026-05-14 sesja P1-bh5 (Claudian).

**Status:** 🟢 ACTIVE — Phase 1 sympy in progress.

**Next deliverables:**
- `Phase1_sympy.py` — 12 tests target ≥9 FP / ≤3 LIT / 0 hardcoded
- `Phase1_sympy.txt` — sympy output saved (PYTHONIOENCODING=utf-8 dla Windows console)
- `Phase1_results.md` — three-layer L1/L2/L3 + per-family table + verdict draft + cross-cycle anchor consistency

**Cross-references:**
- [[./README.md]] (BINDING contract; PR-012)
- [[./Phase0_balance.md]] (6/6 P-gate scope-PASS)
- [[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] §3.4 (family marker source)
- [[../op-bh-alpha-threshold/Phase3_results.md]] T3.2 (BH5 anchor δf/f 8-16%)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1+§2.2+§4 (ASK-RULE + Pattern 2.2 + form-meaning)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 (three-layer L1/L2/L3)
- [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 (Bucket reading)
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-012 LOCKED
