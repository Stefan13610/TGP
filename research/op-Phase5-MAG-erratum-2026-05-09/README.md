---
title: "op-Phase5-MAG-erratum — fix internal inconsistency w Phase 5 gamma identification"
date: 2026-05-09
type: erratum-cycle
status: PHASE0_PHASE1_IN_PROGRESS
folder_status: closed-resolved
classification: ERRATUM_LIGHTWEIGHT
parent: "[[../op-Phi-vacuum-scale-2026-05-09/Phase2_results.md]]"
related_cycles:
  - "[[../op-MAG-resonance-formalization-2026-05-09/]]"
  - "[[../op-Phi-vacuum-scale-2026-05-09/]]"
tgp_owner: research/op-Phase5-MAG-erratum-2026-05-09
tags:
  - erratum-cycle
  - phase5-MAG
  - gamma-identification
  - m_C-correction
---

# op-Phase5-MAG-erratum-2026-05-09

## Geneza

Cykl spawned 2026-05-09 jako **erratum** dla [[../op-MAG-resonance-formalization-2026-05-09/]]
Phase 5 Mach inertia derivation, na podstawie odkrycia w
[[../op-Phi-vacuum-scale-2026-05-09/Phase2_results.md]] §2.

**KRYTYCZNE FINDING z Phase 2:**

Phase 5 sympy ([[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]]
**linie 197-200**) zakłada **simultaneously**:

```python
# In TGP: V''(Phi_0) = -2 beta + 3 gamma = m_C^2
# Assuming beta << gamma (typical), m_C^2 ~ 3 gamma, so gamma ~ m_C^2/3
```

**Internal inconsistency:** Phase 5 zakłada:
- (a) **β = γ** (V_orig vacuum condition, V'(Phi_0) = 0, Phi_eq = Phi_0)
- (b) **β << γ** (żeby uzyskać m_C² ≈ 3γ)

Te dwa założenia są **wzajemnie sprzeczne**.

## Cel cyklu (lightweight erratum)

1. **Confirm** correct identification: **m_C² = γ** EXACTLY (z β=γ vacuum)
2. **Verify** że corrected formula daje consistent results we wszystkich
   Phi_0 scenariuszach (Phi_0 = H_0, v_EW, M_Pl)
3. **Apply erratum** do Phase 5 source files (results.md + sympy)
4. **Document** correction trail (transparency)

## Hipoteza centralna H1

**H1:** Z corrected m_C² = γ (β=γ vacuum):
- m_C = √γ identification depends on what γ is set to
- **Z γ = M_Pl² (T-Λ canonical):** m_C = M_Pl, ratio sqrt(⟨δΦ²_bg⟩)/Phi_0 ~ 10⁻¹⁰ uniwersalny
- Hierarchia v_EW/H_0 NIE jest forced — Phi_0 wolny parameter

**Już potwierdzone w Phase 2 sympy (12/13 PASS).** Niniejszy cykl:
- formal erratum documentation
- apply edits do Phase 5 source files

## Zakres zmian

### Pliki do update:

1. **[[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]]**
   - Add ERRATUM 2026-05-09 sekcja
   - Document corrected γ identification
   - Update quantitative results table

2. **[[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]]**
   - Add comment in lines 197-200 flagging inconsistency
   - NIE rewriting old code — preserve historical record
   - Direct user do new corrected sympy (this cycle)

3. **[[../op-MAG-resonance-formalization-2026-05-09/README.md]]** (jeśli applicable)
   - Mark Phase 5 z erratum status

## Plan szkic

### Phase 0: Balance (DONE w niniejszym README)
### Phase 1: Sympy verification corrected formula + erratum application
- Re-run corrected Phase 5 calculation z m_C = M_Pl
- Verify all Phi_0 scenariusze działają
- Apply edits do Phase 5 source files
- Phase 1 results

## Time budget

Lightweight: ~30 min cycle.

## Cross-references

- [[../op-Phi-vacuum-scale-2026-05-09/Phase2_results.md]] §2 — discovery
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] — to update
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]] — to update

## Status

**SCOPED. Phase 1 in progress.**
