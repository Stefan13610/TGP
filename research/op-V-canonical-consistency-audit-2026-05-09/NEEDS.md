---
title: "NEEDS — op-V-canonical-consistency-audit-2026-05-09"
date: 2026-05-09
type: needs-list
status: COMPLETE_FOR_AUDIT
parent: "[[./README.md]]"
tags:
  - needs
  - phase0
  - V-canonical-audit
---

# NEEDS — op-V-canonical-consistency-audit-2026-05-09

> **Note:** scope post-inventory **dramatycznie zawężony** — G.0 P33 audit
> (2026-05-02) już objął większość framework. Realny scope to **3 residual
> gaps**: T-Λ + Phase 5 MAG + op-Phi-vacuum-scale acknowledgment.

## CRITICAL needs (A1-A2)

### A1: T-Λ closure ρ_vac numerical re-verification z V_M9.1'' canonical
**Priority:** CRITICAL — **flagship gap**
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** T-Λ closed 2026-04-26 (PRZED G.0 V_M9.1'' lock 2026-05-02).
Formula T-Λ:
```
ρ_vac,TGP = V(Phi_eq) = γ·Phi_eq²/12 = M_Pl²·H_0²/12
```

Z V_orig vacuum: V_eq = -γ/12 (przy psi_eq = 1, β=γ)
Z V_M9.1'' canonical: V_min = -4γ/27 (przy psi_eq = 2/3)

**Key question:** czy z V_M9.1'' canonical numerical match z observation
(ratio ~ 1.020) jest zachowany?

**Resolution path:**
1. Sympy: compute |V_M9.1''(psi=2/3)|·Phi_0² i compare to ρ_vac,obs
2. Z Phi_eq=(2/3)·Phi_0 i Phi_eq=H_0 → Phi_0=(3/2)·H_0
3. ρ_vac,V_M9.1'' = (4γ/27)·Phi_0² = (4·M_Pl²/27)·(3H_0/2)² = M_Pl²·H_0²/3
4. Compare to obs 2.518e-11 eV⁴

**Pre-cycle expectation:** ratio LIKELY 4× off (M_Pl²·H_0²/3 vs M_Pl²·H_0²/12).
Jeśli tak, T-Λ wymaga re-derivation lub re-interpretacja "Phi_0" w T-Λ
oznacza coś innego niż V_M9.1'' parameter.

### A2: Phase 5 MAG Phi_0 reference clarification
**Priority:** CRITICAL — **wymaga przed jakąkolwiek follow-up cyklu**
**Phase:** 1 (audit)
**Status:** OPEN

**Description:** Phase 5 MAG formula:
```
m_Mach = (3γq²)/(16π Phi_0² m_C) · ⟨δΦ²_bg⟩
```

Phi_0 nie ma explicit V-source citation. Ambiguity:
- (a) V_orig Phi_0 (Phi_eq = Phi_0)
- (b) V_M9.1'' Phi_0 (Phi_eq = (2/3)·Phi_0)
- (c) Independent parameter

**Resolution path:**
1. Read Phase 5 derivation [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]]
2. Check git history / closure date — pre-G.0 czy post-G.0?
3. Sympy: czy formula invariant pod V change, czy zmienia coefficient?

**Falsifier:** jeśli Phase 5 explicit cytuje V_orig, formula może wymagać re-derivation.

## HIGH needs (A3-A4)

### A3: op-Phi-vacuum-scale Phase 1 acknowledgment
**Priority:** HIGH — **already partially done w Phase 1.6 user iteration**
**Phase:** 1 (cross-cycle update)
**Status:** PARTIAL ✅

**Description:** op-Phi-vacuum-scale Phase 1 reconnaissance (2026-05-09)
cytowała deprecated V_orig w T13. Phase 1.6 user iteration explicit acknowledged.

**Resolution path:** post-Phase 1 niniejszego audit cyklu, dodać formal
cross-reference do op-Phi-vacuum-scale z aktualnym findings (T-Λ status).

### A4: Cross-cycle classification table
**Priority:** HIGH (deliverable D1)
**Phase:** 1 (documentation)
**Status:** OPEN

**Description:** Produce explicit table per cykl:

| Cykl | Closure date | V usage | Status (pod V_M9.1'') | Action required |
|---|---|---|---|---|

Z minimum 8 cykli: T-Λ, MAG Phase 5, MAG-Lorentz, UV.3, particle_sector,
op-Phi-decomposition-photon, op-Phi-vacuum-scale Phase 1, sek08a master.

**Resolution path:** kompilacja z G.0 P33 + niniejszy Phase 1 sympy.

## SUPPORTING needs (A5)

### A5: Recommendations dla follow-up cykli
**Priority:** MEDIUM (deliverable D4)
**Phase:** 1 (results)
**Status:** OPEN

**Description:** Po Phase 1 audit, jasne recommendations:
- Cykle wymagające re-derivation (BLOCKER)
- Cykle gdzie V_orig był valid approximation (NO ACTION)
- Cykle independent of V (NO ACTION)
- Recommended order dla follow-up cycles (op-multi-vacuum-identification, op-EWSB-from-substrate)

## Priority matrix

| ID | Title | Priority | Phase | Status | Resolution |
|---|---|---|---|---|---|
| A1 | T-Λ ρ_vac z V_M9.1'' | CRITICAL | 1 | OPEN | Sympy verification |
| A2 | Phase 5 Phi_0 reference | CRITICAL | 1 | OPEN | Audit + sympy |
| A3 | op-Phi-vacuum-scale ack | HIGH | 1 | ✅ PARTIAL | Cross-ref update |
| A4 | Cross-cycle table | HIGH | 1 | OPEN | Compilation |
| A5 | Follow-up recommendations | MEDIUM | 1 | OPEN | Post-audit synthesis |

## Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[../op-g0-r3-from-canonical-projection/P33_audit_results.md]]
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]]
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]]
- [[../op-Phi-vacuum-scale-2026-05-09/Phase1_reconnaissance_results.md]]
