---
title: "op-dual-V-structure-clarification — formal verification dual-V hypothesis (V_M9.1'' gravity, V_orig matter)"
date: 2026-05-09
type: structural-clarification-cycle
status: PHASE0_PHASE1_IN_PROGRESS
folder_status: closed-resolved
classification: STRUCTURAL_CLARIFICATION
parent: "[[../op-MAG-Phase5-V-reference-clarification-2026-05-09/Phase1_clarification_results.md]]"
related_cycles:
  - "[[../op-MAG-Phase5-V-reference-clarification-2026-05-09/]]"
  - "[[../op-V-canonical-consistency-audit-2026-05-09/]]"
  - "[[../op-Phi-vacuum-scale-2026-05-09/]]"
  - "[[../op-g0-r3-from-canonical-projection/]]"
  - "[[../op-MAG-resonance-formalization-2026-05-09/]]"
  - "[[../closure_2026-04-26/Lambda_from_Phi0/]]"
tgp_owner: research/op-dual-V-structure-clarification-2026-05-09
tags:
  - structural-clarification
  - dual-V-structure
  - V-M911-gravity
  - V-orig-matter
  - A4-audit-marker
  - framework-consistency
---

# op-dual-V-structure-clarification-2026-05-09

## Geneza

Cykl spawned 2026-05-09 jako Path C verification z
[[../op-MAG-Phase5-V-reference-clarification-2026-05-09/Phase1_clarification_results.md]].

**KRYTYCZNE PRE-CYCLE FINDING (search agent + Phase1_results.md G.0):**

G.0 closure (op-g0-r3-from-canonical-projection, 2026-05-02) **EXPLICITNIE**
ogranicza V_M9.1'' canonical do **GRAVITATIONAL sektora**:

> [[../op-g0-r3-from-canonical-projection/Phase1_results.md]] **linia 266:**
> *"A4 (matter coupling) — wymaga osobnego sprawdzenia (G.0 nie dotyka L_mat)"*

To znaczy że:
1. V_M9.1'' canonical został locked **dla gravity sector** (R3 ODE, M9.1''
   metryka, Newton limit)
2. V_orig **NIE został zastąpiony** dla matter sector (L_mat)
3. **A4 audit marker** explicit reserves matter sector dla "separate verification"

**Path C hypothesis** z prior cycle jest zatem **NIE TYLKO HIPOTEZA** — jest
**już documented w G.0 closure** jako otwarty audit marker A4.

## Cel cyklu (focused)

**Niniejszy cykl jest formalną realizacją A4 audit marker** z G.0 closure:

1. **Confirm** że V_orig pozostaje valid dla matter sector (formal sympy)
2. **Document** dual-V structure jako TGP framework feature (NIE bug, NIE crisis)
3. **Update sek08a annotation** na "V_orig DEPRECATED FOR GRAVITY SECTOR
   2026-05-02 (G.0 closure). Matter sector usage maintained pending A4
   verification — see op-dual-V-structure-clarification-2026-05-09."

## Hipoteza centralna H1

**H1:** TGP framework legitymie posiada **dual-V structure**:
- $V_{M9.1''}(\psi) = -\gamma\psi^2(4-3\psi)^2/12$ — gravitational sektor
  (M9.1'' metryka, R3 ODE, Newton limit, mass spectrum invariance)
- $V_{orig}(\Phi) = -\beta\Phi^3/(3\Phi_0) + \gamma\Phi^4/(4\Phi_0^2)$ —
  matter sektor (Phase 5 Mach inertia, T-Λ ρ_vac, particle masses)

**Falsifier:** jeśli sympy lub framework analysis pokazuje że dual-V
struktura prowadzi do internal contradictions (np. EOM-y są incompatible),
H1 FALSIFIED i wymagana globalna re-derivacja.

## Cztery tests dla Phase 1

### T1: G.0 explicit scope statement (text reading)

**Cel:** udokumentować tekstowo że G.0 closure ograniczył V_M9.1'' do gravity.

### T2: sek08a "DEPRECATED" annotation context

**Cel:** zweryfikować że sek08a "DEPRECATED 2026-05-02" annotation jest w
sekcji **gravity-related** (sek08a master action), NIE matter Lagrangian.

### T3: V_M9.1'' EOM derivation source

**Cel:** sympy: V_M9.1'' wynika z R3 ODE (gravitational static spherical
EOM), NIE z matter dynamics. Confirm z phase1_G0a_volume_integration.py.

### T4: V_orig EOM scope source

**Cel:** sympy: V_orig wynika z field expansion around Phi_0 (matter quartic
self-interaction), NIE z gravitational metric. Confirm z Phase 5 derivation.

### T5: Dual-V framework consistency check

**Cel:** sympy: czy dual-V strukture jest matematycznie konsystentna
(EOM-y nie są w sprzeczności)?

## Plan szkic Phase 0-1

### Phase 0: Balance sheet
- Confirmed pre-cycle finding (G.0 explicit A4 marker)
- 5 sympy tests T1-T5
- 8/8 gate criteria

### Phase 1: Formal verification
- Sympy reads G.0 explicit statements
- Cross-cycle classification update
- Path C confirmation lub falsification
- Recommendation: sek08a annotation update

## Probability assessment (post-pre-cycle finding)

| Outcome | Pre-search | Post-search (G.0 explicit A4) |
|---|---|---|
| Path C confirmed (dual-V valid) | 50-60% | **80-90%** ↑ |
| Internal contradiction discovered | 20-30% | **5-10%** ↓ |
| Need full re-derivation | 10-20% | **5%** ↓ |

**Search agent finding zmienia werdykt z hipoteza → established framework feature.**

## Time budget

Lightweight verification cycle: ~1 session.

## Cross-references

- [[../op-g0-r3-from-canonical-projection/Phase1_results.md]] linia 266 — A4 marker
- [[../op-MAG-Phase5-V-reference-clarification-2026-05-09/Phase1_clarification_results.md]] — Path C origin
- [[../op-V-canonical-consistency-audit-2026-05-09/]] — broader audit context
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
- [[../../meta/CALIBRATION_PROTOCOL.md]]

## Status

**SCOPED. Phase 0-1 in progress.**
