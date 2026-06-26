---
title: "Phase 4 plan — Ω_m, CMB, BBN compatibility (F5, F6, F7)"
type: phase_plan
status: PRE_REGISTERED_LOCKED
phase: 4
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
pre_registration_date: 2026-05-23
---

# Phase 4 plan — F5 (Ω_m), F6 (CMB), F7 (BBN) compatibility

**Status:** PRE_REGISTERED_LOCKED 2026-05-23.
**Methodology note:** F5-F7 są **strukturalnie trudniejsze** od F4 (PRIMARY KILLER). F4 miał czyste geometric derivation; F5-F7 wymagają full early-universe TGP-native model + nucleosynthesis machinery. Pre-registered: **PARTIAL OR DEFERRED status expected** dla części F5-F7.

---

## §1 — Pre-registered FP

**0 hardcoded T_pass=True.** Każdy FP = compute-then-compare.

### T_P4_1 — Ω_m structural argument (F5 SECONDARY KILLER)

**Pre-registered:** F5 = Ω_m,critical ∈ [0.155, 0.62] factor 2 (target 0.31).

**TGP-native concept (paper §5):** "Equilibrium density (E2 stability condition)".

**Pre-registered test:** Compute TGP-native characteristic densities:
- Soliton volume scale: V_sol ~ (1/m_σ)³
- Saturation density: n_sat ~ m_σ³
- Observed baryon density: n_obs ~ 10⁻⁷ /cm³

**Pre-registered (literal):**
- PASS: TGP-native derivation of Ω_m ∈ [0.155, 0.62]
- PARTIAL: structural argument provided; numerical match NIE direct
- FAIL: TGP-native predicts Ω_m outside factor 2 → F5 FAIL

**Anti-Lakatos:** Honest disposition pre-acknowledged — F5 may PARTIAL.

### T_P4_2 — CMB blackbody shape (F6 HARD CONSTRAINT)

**Pre-registered:** F6 = blackbody deviation < 10⁻⁴.

**TGP-native concept (paper §5):** "Relict thermal sygnatura saturacji E2".

**Pre-registered test:**
- Thermal equilibrium at frontier saturation → blackbody spectrum (generic)
- Shape match: ✓ trivially if thermal equilibrium achieved
- Temperature T = 2.725 K: NIE pre-registered as derived; observed value

**Pre-registered (literal):**
- PASS shape: blackbody generic for thermal equilibrium ✓
- PARTIAL temperature: T = 2.725 K accepted as observational input (analogous to t_universe)
- FAIL: TGP-native predicts non-blackbody OR fundamental thermodynamic inconsistency

### T_P4_3 — BBN compatibility (F7 HARD CONSTRAINT)

**Pre-registered:** F7 = D/H, ⁴He/H, ⁷Li/H w standard uncertainty.

**TGP-native concept:** Nucleosynthesis happens in early universe; depends on:
- Expansion rate history H(t) at z ~ 10⁹
- Temperature history T(t) at z ~ 10⁹
- Neutron-proton ratio at freezeout

**TGP-native disposition:**
- H(t) = 1/t (geometric, late-time) — extrapolation to early universe NIE necessarily valid
- TGP early-universe model NIE yet developed
- Pre-registered as **DEFERRED** for full TGP early-universe cycle (likely Poziom γ-4 or Poziom δ)

**Pre-registered (literal):**
- PASS: TGP-native BBN calculation matches standard ratios
- DEFERRED: explicitly outside Phase 4 scope (requires multi-cycle effort)
- FAIL: TGP-native predicts incompatible BBN

### T_P4_4 — Early-universe expansion model honest disposition

**Pre-registered:** Document explicitly that R(t) = c·t formula (Phase 2) is **late-time approximation**, NIE valid universally.

**Issues:**
- Early universe: matter-dominated → radiation-dominated requires DIFFERENT R(t) form
- TGP-native early model NIE yet developed (Phase 0 §4.3 (c) "Late-time: Focus on current epoch (z << 1)")
- F6 + F7 strictly require early universe → PARTIAL or DEFERRED

**Pre-registered (literal):**
- PASS: honest declaration of late-time limit; F6 + F7 deferred z explicit acknowledgment
- FAIL: false claim of universal R = c·t validity → Lakatos hazard

---

## §2 — Anti-Lakatos disposition

**Critical (after Phase 3 PASS_TARGET):**

1. NIE inflate Phase 4 results beyond honest status. F4 PASS does NOT imply F5-F7 PASS automatically.

2. NIE claim BBN/CMB without derivation. Acknowledge limits of TGP-native model at current development stage.

3. Mark DEFERRED clearly (NOT FAIL). DEFERRED = future work, NOT falsification.

4. **Risk R1 flag:** "F4 geometric easy, F5-F7 model-development hard" pattern may signal **methodology gap** — does TGP-native cosmology have enough internal structure to test F5-F7 in current form, or do these tests require pre-Bing-Bang style model development?

---

## §3 — DEC tracking

- DEC 2 NIE expended (no numerical method needed)
- DEC 3 reserved

---

## §4 — Computational plan (Phase4_sympy.py)

**Tool:** sympy symbolic.

**Sections:**
1. T_P4_1: Ω_m structural argument z m_σ scale
2. T_P4_2: CMB blackbody generic from thermal equilibrium
3. T_P4_3: BBN explicitly DEFERRED with structural rationale
4. T_P4_4: Late-time vs early-universe disposition declaration

---

## §5 — Honest expected outcome

**Pre-registered Phase 4 disposition:**

- T_P4_1 (Ω_m): **PARTIAL** likely — structural argument provided; numerical match limited
- T_P4_2 (CMB shape): **PASS** generic; temperature observational input
- T_P4_3 (BBN): **DEFERRED** — outside Phase 4 scope (future cycle)
- T_P4_4 (early-universe): **PASS** structural — late-time limit honest acknowledgment

**Expected verdict:** **Mixed (PARTIAL + DEFERRED)** dla F5-F7.

**Phase FINAL implications:** γ-3 claim_status likely **A** (not A+) jeśli F4 PASS + F5-F7 PARTIAL/DEFERRED.

---

**Phase 4 plan LOCKED 2026-05-23. Honest disposition pre-acknowledged. Ready dla Phase4_sympy.py.**
