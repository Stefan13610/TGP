---
title: "op-MAG-Phase5-V-reference-clarification — explicit clarification V_orig (DEPRECATED) cytowane w Phase 5 MAG"
date: 2026-05-09
type: audit-clarification-cycle
status: PHASE0_PHASE1_IN_PROGRESS
folder_status: closed-resolved
classification: AUDIT_LIGHTWEIGHT_CLARIFICATION
parent: "[[../op-V-canonical-consistency-audit-2026-05-09/Phase1_audit_results.md]]"
related_cycles:
  - "[[../op-V-canonical-consistency-audit-2026-05-09/]]"
  - "[[../op-MAG-resonance-formalization-2026-05-09/]]"
  - "[[../op-Phi-vacuum-scale-2026-05-09/]]"
tgp_owner: research/op-MAG-Phase5-V-reference-clarification-2026-05-09
tags:
  - audit-cycle
  - V-reference-clarification
  - phase5-MAG
  - mach-inertia
  - V-orig-cited
---

# op-MAG-Phase5-V-reference-clarification-2026-05-09

## Geneza

Cykl spawned 2026-05-09 jako residual gap resolver z
[[../op-V-canonical-consistency-audit-2026-05-09/Phase1_audit_results.md]] §1.2.

**KRYTYCZNE odkrycie:** Phase 5 MAG Mach inertia derivation
([[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] linia 38)
**explicitnie cytuje DEPRECATED V_orig:**

```
V(Φ) = -β Φ³/(3 Φ_0) + γ Φ⁴/(4 Φ_0²)
```

To jest formula **V_orig (DEPRECATED 2026-05-02)** w sek08a (linie 95-110).
Phase 5 closed 2026-05-09 (**PO** V_M9.1'' canonical lock), ale **nadal używa
V_orig**. To NIE jest tylko ambiguity — to realny residual gap.

## Cel cyklu (lightweight)

1. **Document** Phase 5 explicit V_orig usage (already verified)
2. **Sympy verify** structural impact: czy m_Mach formula zmienia się qualitatively
   pod V_orig→V_M9.1'' canonical?
3. **Clarification recommendation:** czy Phase 5 wymaga full re-derivation (heavy)
   czy wystarczy re-interpretation Phi_0 reference (light)?

## Hipoteza centralna H1

**H1:** Phase 5 m_Mach formula zmienia się **qualitatively** pod V_M9.1''
canonical, bo λ_4 (quartic coupling przy V''''/4) zmienia znak:
- V_orig λ_4 = +3γ/(2 Phi_0²) (positive)
- V_M9.1'' λ_4 = -9γ/2 (NEGATIVE — przy ψ=2/3 vacuum)

Negative λ_4 implikuje m_Mach zmienia znak — quantitative match m_e=511 keV
w V_orig framework **NIE jest preserved** w V_M9.1'' canonical.

**Falsifier:** jeśli λ_4 sign jest invariant pod V change (np. specyficzny inne
expansion point), Phase 5 jest invariant.

## Three resolution paths

**Path A (lightweight):** Re-interpretation: Phi_0_Phase5 = (2/3)·Phi_0_V_M911.
Phase 5 derivation jest valid w V_orig konwencji ale "Phi_0" tam to NIE
V_M9.1'' Phi_0 parameter — to V_orig vacuum value (= V_M9.1'' (2/3)·Phi_0).
Wtedy v_EW=246 GeV w Phase 5 odpowiada Phi_0_V_M911 = 369 GeV.

**Path B (medium):** Re-derivation Phase 5 around V_M9.1'' true minimum (ψ=2/3).
Wymaga recompute V'', V''', V'''' przy ψ=2/3, redo path-integral expansion.
Może dać różne formuła m_Mach.

**Path C (heavy):** Full Phase 5 re-derivation z V_M9.1'' i wszystkimi vacuum
configurations rozważanymi (ψ=2/3 minimum, ψ=4/3 horizon, ψ=1 reference).
To może być **inny cykl** (multi-vacuum identification).

## Plan szkic Phase 0-1

### Phase 0: Balance sheet (DONE w niniejszym README)
- Identify Phase 5 explicit V_orig citation
- Inventory: czy inne cykle MAG też cytują V_orig?
- Hipoteza H1 + 3 paths

### Phase 1: Audit sympy + verdict
- T1-T8 sympy: V''''(Phi_0) w V_orig vs V_M9.1''
- λ_4 sign comparison
- m_Mach formula structural impact
- Path A/B/C verdict

## Probability assessment

| Outcome | Prob |
|---|---|
| Path A clean (re-interpretation OK) | 50-60% |
| Path B needed (re-derivation) | 30-40% |
| Path C unavoidable (multi-vacuum spawn) | 10-20% |

## Time budget

Lightweight cycle: ~1 session.

## Cross-references

- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] — Phase 5 source
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]] — sympy verification
- [[../op-V-canonical-consistency-audit-2026-05-09/]] — parent audit
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] — V_orig deprecated source

## Status

**SCOPED. Phase 0-1 in progress.**
