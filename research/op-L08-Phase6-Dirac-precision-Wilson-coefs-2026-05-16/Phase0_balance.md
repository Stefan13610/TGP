---
title: "Phase 0 — Balance sheet: Wilson coefs Φ-dependent corrections + a_e lab-scale"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-balance
phase: 0
gate_status: 8/8 ☑ PASS
tags:
  - balance-sheet
  - phase0
  - Wilson-coefs
  - a_e-prediction
---

# Phase 0 — Balance sheet (8/8 ☑ gate)

## ☑ Check 1 — External inputs declared

| Input | Source |
|---|---|
| PDG 2024 a_e = 0.00115965218073(28) | Particle Data Group |
| α = 1/137.0360 (CODATA) | CODATA |
| α/(2π) ≈ 0.0011614098 | Schwinger 1948 (QED leading order) |
| QED 5-loop a_e theory | Aoyama et al. 2020 (a_e^QED ≈ 0.001159652181643) |
| Bjorken-Drell §15.3 propagator | textbook |

**LIVE inheritance:**
- L08-Dirac S_F^TGP(p)|_{Φ_0} = S_F^Dirac canonical
- L01-N1 prefactor α/(3π) ≈ 7.74·10⁻⁴ (1-loop QED structure)
- B9 MICROSCOPE η_TGP_lab = 1.32·10⁻²⁶ baseline
- L05 m_obs LIVE; k_obs(α=2,d=3)=3
- Q2 F1 Φ_eq(today) = H_0 anchor

**Status:** ☑ External inputs + LIVE inheritance enumerated.

## ☑ Check 2 — Structural axioms LOCKED

- **S05** single-Φ axiom (FOUNDATIONS §1) — Phase 1 T11 DEC
- **S04 + B9** universal coupling g_eff[Φ] (MICROSCOPE 6/6 PASS, 8.3·10¹⁰× safe)
- **ax:c-ax:G-ax:ℏ** (sek04_stale) — coupling rules dla Φ-dependence
- **L01** ρ ≡ -T^μ_μ/c_0² formal definition

**Status:** ☑ All LOCKED axioms enumerated.

## ☑ Check 3 — Derived outputs

1. **Wilson coef expansion** S_F^TGP[Φ] = S_F|_{Φ_0} + ζ_1(δΦ/Φ_0) S^(1) + O((δΦ)²)
2. **ζ_1 derived** z g_eff[Φ] background structure (NIE free parameter)
3. **Lab-scale (δΦ/Φ_0)_lab** estimate z B9 MICROSCOPE inheritance
4. **a_e^TGP correction magnitude** prediction
5. **PDG comparison** vs 0.28 ppb experimental precision

**Status:** ☑ Operationalized outputs.

## ☑ Check 4 — Tautology test

- ❌ NIE zakładamy Wilson coefs a priori — derive z g_eff[Φ] structure
- ❌ NIE zakładamy lab-scale (δΦ/Φ_0) — inherit z B9 baseline (independent input)
- ❌ NIE zakładamy a_e^TGP = a_e^QED — compute correction explicitly
- ✅ PDG comparison jest **test**, NIE assumption

**Status:** ☑ NO tautology.

## ☑ Check 5 — Falsifiability

| Sub-test | Falsifies if... |
|---|---|
| T1-T3 Wilson expansion | NIE applicable structure → HALT-B |
| T4-T5 ζ_i derivation | ζ_i undetermined z g_eff → HALT-B |
| T6 lab-scale δΦ | NIE compatible z B9 baseline → revisit inheritance |
| T7-T8 a_e correction | magnitude > 1 ppb → tension z PDG |
| T9-T10 PDG comparison | > 0.28 ppb (PDG precision) → falsification at lab |
| T11 DEC S05+B9 | violation → critical, cycle abandon |

**Pre-registered rule §0.2 BINDING. No post-hoc weakening.**

**Status:** ☑ Explicit gates.

## ☑ Check 6 — Independent path

≥3 sources:
1. **L08-Dirac** S_F^TGP construction (A− 2026-05-16)
2. **L01-N1** Wilson coef framework reference (A− 2026-05-11)
3. **B9 MICROSCOPE** η_lab baseline (LOCKED 2026-05-01)
4. **PDG 2024** a_e measurement (external LIT)
5. **QED multi-loop** literature (Aoyama et al. 2020)

**Status:** ☑ 5 sources independent.

## ☑ Check 7 — Risk assessment

| R# | Risk | Mitigation |
|---|---|---|
| **R1** | Wilson coef expansion may break in non-perturbative regime | Lab-scale (δΦ/Φ_0)~10⁻²⁶ ⇒ perturbative dominate; documented |
| **R2** | ζ_i derivation może być ambiguous (gauge choice) | Use universal coupling form g_eff[Φ] = (Φ/Φ_0)^k canonical — fix gauge via S04 |
| **R3** | a_e^TGP correction may be unconstrained beyond order-of-magnitude | Use B9 inheritance + L01-N1 prefactor structure dla scale estimate |
| **R4** | Multi-loop QED comparison not in scope | Phase 2 deferred; this cycle: leading-order estimate only |
| **R5** | Lab-scale δΦ/Φ_0 not directly measurable | Use B9 universal coupling bound as proxy (8.3·10¹⁰× safe) |

**Status:** ☑ 5 R-flags z mitigation.

## ☑ Check 8 — Honesty test (0 hardcoded)

Phase 1 sympy commitment: **0 hardcoded `T_pass = True`**.

Pattern:
```python
T_pass = (sp.simplify(expr_lhs - expr_rhs) == 0) AND (some condition)
```

Substance target Phase 1: 8 FP + 2 LIT + 1 DEC = 11 tests total; ~73% FP.

**Status:** ☑ Phase 6 ABSOLUTE BINDING enforced.

---

## Verdict

| Check | Status |
|---|---|
| 1-8 | ☑ ALL PASS |

🟢 **8/8 ☑ PASS — Phase 1 ENABLED.**

## Phase 1 sympy plan (11 tests)

| # | Test | Type |
|---|---|---|
| T1 | S_F^TGP[Φ] expansion zeroth-order = S_F^Dirac (L08-Dirac inheritance) | FP |
| T2 | First-order Wilson coef expansion δS_F = ζ_1 (δΦ/Φ_0) S^(1) | FP |
| T3 | S^(1) explicit form z g_eff[Φ] = (Φ/Φ_0)^k_geo perturbation | FP |
| T4 | ζ_1 derived z m_obs[Φ] = m_obs|_{Φ_0} · (Φ/Φ_0)^k_mass expansion | FP |
| T5 | k_mass = e²(1-α/4)·(1/2) z why_n3 + L05 (α=2) reconciliation | FP |
| T6 | Lab-scale δΦ/Φ_0|_lab ≤ B9 baseline ~10⁻²⁶ inheritance | LIT |
| T7 | a_e^QED leading order = α/(2π) (LIT reference) | LIT |
| T8 | a_e^TGP correction = ζ_1 · (δΦ/Φ_0)_lab · prefactor estimate | FP |
| T9 | a_e^TGP correction magnitude vs 0.28 ppb PDG precision | FP |
| T10 | Free-field limit consistency: lim_{Φ→Φ_0} correction → 0 | FP |
| T11 | DECLARATIVE S05 + B9 + universal coupling preservation | DEC |

**Substance target:** 8 FP + 2 LIT + 1 DEC = 11 tests; 72.7% FP.

## Cross-references

- [[./README.md]] §0 BINDING contract
- [[../op-L08-Phase6-Dirac-propagator-2026-05-16/Phase_FINAL_close.md]] (predecessor inheritance)
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase_FINAL_close.md]] (Wilson framework reference)
