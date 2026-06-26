---
title: "Phase 6 plan — F9 (NULL CONSISTENCY) + F-γ-4 (SPECULATIVE confinement)"
type: phase_plan
status: PRE_REGISTERED_LOCKED
phase: 6
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
pre_registration_date: 2026-05-23
---

# Phase 6 plan — F9 + F-γ-4

**Status:** PRE_REGISTERED_LOCKED 2026-05-23.

---

## §1 — Pre-registered FP

**0 hardcoded T_pass=True.**

### T_P6_1 — F9 null consistency structural argument

**Hypothesis:** TGP-native E2 bulk saturated (Φ = v at VEV) → no driving force for spontaneous local creation.

**Pre-registered (literal):**
- Symbolic argument: V'_TGP(Φ = v) = 0 (potential minimum)
- d⟨Φ⟩/dt = 0 in bulk (already at VEV)
- No spontaneous creation source in bulk
- PASS: derivation consistent

### T_P6_2 — F9 observational match

**Pre-registered:** Observed = no spontaneous proton/quark creation in any lab.

**Pre-registered (literal):**
- PASS: TGP predicts null + observation = null → match
- FAIL: TGP predicts local creation OR observation contradicts (neither expected)

### T_P6_3 — F-γ-4 D_critical mapping to QCD T_c

**Pre-registered:** TGP-native characteristic scale (m_σ from γ-1 retry) compared do QCD T_c ≈ 150 MeV.

**Pre-registered (literal):**
- m_σ = v√(2λ); v ~ 200 MeV (Λ_QCD); √(2λ) ~ O(1) → m_σ ~ 200-300 MeV
- T_c QCD ~ 150-170 MeV
- Ratio m_σ / T_c ~ 1-2 (well within factor 10)
- PASS if ratio < 10

### T_P6_4 — F-γ-4 SPECULATIVE verdict

**Pre-registered:** F-γ-4 was labeled SPECULATIVE w concept paper §7. Speculative test może być PARTIAL acceptable.

**Pre-registered (literal):**
- PASS: m_σ / T_c < 10 + structural mapping reasonable
- PARTIAL: speculative, but ballpark agreement
- FAIL: ratio > 10 or mapping conceptually wrong

---

## §2 — Anti-Lakatos disposition

- F9: auto-PASS expected (null prediction matches null observation)
- F-γ-4: order-of-magnitude PASS expected; honestly speculative

---

## §3 — DEC tracking

- DEC 2, 3 still reserved (likely unused)

---

## §4 — Computational plan (Phase6_sympy.py)

**Sections:**
1. T_P6_1: V'(v) = 0 sympy verification
2. T_P6_2: F9 null match declaration
3. T_P6_3: m_σ / T_c ratio computation
4. T_P6_4: F-γ-4 verdict

---

**Phase 6 plan LOCKED 2026-05-23.**
