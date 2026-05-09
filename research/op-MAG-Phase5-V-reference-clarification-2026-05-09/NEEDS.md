---
title: "NEEDS — op-MAG-Phase5-V-reference-clarification-2026-05-09"
date: 2026-05-09
type: needs-list
status: COMPLETE
parent: "[[./README.md]]"
---

# NEEDS — op-MAG-Phase5-V-reference-clarification

## CRITICAL needs (B1-B3)

### B1: Sympy verification λ_4 V_orig vs V_M9.1''
**Priority:** CRITICAL
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** Compare V''''(at vacuum)/4 dla obu potencjałów:
- V_orig: V''''(Phi_0)|β=γ = 6γ/Phi_0² → λ_4 = 3γ/(2 Phi_0²)
- V_M9.1'': V''''(ψ_eq=2/3) = ? (compute z sympy)

**Resolution path:** sympy compute, compare, report sign + magnitude.

### B2: Impact assessment on m_Mach formula
**Priority:** CRITICAL
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** Czy m_Mach = λ_4·⟨δΦ²_bg⟩·∫δΦ²_sol zachowuje znak i scale
pod V_orig → V_M9.1'' canonical?

**Resolution path:** sympy: substitute λ_4_V_M9.1'' do m_Mach, check sign.

### B3: Re-interpretation feasibility (Path A)
**Priority:** HIGH
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** Czy Path A (re-interpretation Phi_0_Phase5 = (2/3)·Phi_0_V_M911)
naturally rozwiązuje gap, czy Path B (re-derivation) jest unavoidable?

**Resolution path:** sympy: assume Phi_0_Phase5 = (2/3)·Phi_0_V_M911, check
czy expansion around (2/3)·Phi_0_V_M911 = Phi_0_Phase5 reproduces V_orig structure.

## SUPPORTING needs (B4-B5)

### B4: Cross-impact w innych Phase 5 quantitative results
**Priority:** MEDIUM
**Phase:** 1
**Status:** OPEN

**Description:** Phase 5 quantitative claim m_e=511 keV przy Phi_0=v_EW.
Jeśli λ_4 zmienia znak/magnitude, czy v_EW scenario nadal działa?

### B5: Recommendations dla follow-up cycles
**Priority:** HIGH (deliverable)
**Phase:** 1 (results)
**Status:** OPEN

**Description:** Po Phase 1, jasne recommendation:
- Path A (re-interpretation) — minor update do Phase 5 documentation
- Path B (re-derivation) — spawn `op-Phase5-V-M911-rederivation`
- Path C (multi-vacuum) — koordynacja z `op-multi-vacuum-identification`

## Priority matrix

| ID | Title | Priority | Phase | Status |
|---|---|---|---|---|
| B1 | λ_4 sign comparison | CRITICAL | 1 | OPEN |
| B2 | m_Mach formula impact | CRITICAL | 1 | OPEN |
| B3 | Path A feasibility | HIGH | 1 | OPEN |
| B4 | v_EW quantitative impact | MEDIUM | 1 | OPEN |
| B5 | Recommendations | HIGH | 1 | OPEN |

## Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]]
- [[../op-V-canonical-consistency-audit-2026-05-09/]]
