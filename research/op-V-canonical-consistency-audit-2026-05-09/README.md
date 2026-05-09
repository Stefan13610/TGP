---
title: "op-V-canonical-consistency-audit — audit V_orig (DEPRECATED) vs V_M9.1'' (canonical) we wszystkich cyklach TGP"
date: 2026-05-09
type: audit-cycle
status: PHASE0_PHASE1_IN_PROGRESS
folder_status: closed-resolved
classification: AUDIT_DOCUMENTATION_CYCLE
parent: "[[../op-Phi-vacuum-scale-2026-05-09/Phase1_reconnaissance_results.md]]"
related_cycles:
  - "[[../op-Phi-vacuum-scale-2026-05-09/]]"
  - "[[../closure_2026-04-26/Lambda_from_Phi0/]]"
  - "[[../op-MAG-resonance-formalization-2026-05-09/]]"
  - "[[../op-uv3-phi0-renormalization/]]"
  - "[[../particle_sector_closure/]]"
tgp_owner: research/op-V-canonical-consistency-audit-2026-05-09
tags:
  - audit-cycle
  - V-canonical
  - V-M911
  - V-orig-deprecated
  - framework-consistency
  - blocker-resolution
---

# op-V-canonical-consistency-audit-2026-05-09

## Geneza

Cykl spawned 2026-05-09 jako **BLOCKER resolution** dla
[[../op-Phi-vacuum-scale-2026-05-09/]] **P11** (NEEDS.md).

Phase 1.6 user iteration ([[../op-Phi-vacuum-scale-2026-05-09/Phase1_6_strong_field_canonical_sympy.py]],
15/15 PASS) odkryła krytyczne pominięcie:

**V_orig(Φ) = (β/3)Φ³ - (γ/4)Φ⁴ jest DEPRECATED 2026-05-02** (sek08a linie
95-110). Canonical TGP używa:

$$V_{M9.1''}(\psi) = -\frac{\gamma}{12}\psi^2(4-3\psi)^2$$

W canonical V_M9.1'':
- **NIE MA** wolnego parametru β (tylko γ overall scale)
- **Multi-vacuum strukture**: ψ ∈ {0, 2/3, 4/3}
- Inne fundamentalne wartości (V_min = -4γ/27 zamiast V_min = -γ/12)

Phase 1 reconnaissance op-Phi-vacuum-scale **cytowała deprecated V_orig**, co
oznacza że **inne cykle mogą również używać deprecated formuły**. To **blokuje**
wszystkie follow-up cycles (op-multi-vacuum-identification, op-EWSB-from-substrate)
do czasu, gdy framework consistency zostanie zaudytowana.

## Cel cyklu

**Audit każdego cyklu TGP** używającego V(Φ) potencjału dla:
1. **Klasyfikacja**: który cykl używa V_orig (deprecated) vs V_M9.1'' (canonical)?
2. **Impact assessment**: jakie wnioski cyklu zmieniają się jeśli przepiszemy z V_orig → V_M9.1''?
3. **Recommendation**: czy cykl wymaga formal update / re-derivation, czy też V_orig był approximation valid w specyficznym regime?

**Nie jest celem:** re-derivacja każdego cyklu. Tylko **identyfikacja** rozbieżności i flagi do dalszej pracy.

## Hipoteza centralna H1

**H1:** Większość cykli TGP closed pre-2026-05-02 używała V_orig (gdy
formuła była canonical). Po deprecacji (2026-05-02), nowe cykle MUSZĄ używać
V_M9.1''. **Niektóre wnioski mogą wymagać aktualizacji** (np. T-Λ V_min = -γ/12
vs canonical -4γ/27 zmienia liczbowy współczynnik o 16/9 ≈ 1.78).

**Falsifier:** jeśli wszystkie audit-cykle są strukturalnie OK z V_M9.1'',
niezależnie od oryginalnej formuły, wnioskujemy że V_orig był legalna
approximation w pewnym regime — **brak BLOCKER**, framework consistent.

## Zaplanowany scope (audit list)

Cykle używające V(Φ) potencjału (do audit):

| Cykl | Status | V_orig czy V_M9.1''? | Impact |
|---|---|---|---|
| **closure_2026-04-26/Lambda_from_Phi0** (T-Λ) | CLOSED 2026-04-26 | TBD (audit Phase 1) | wpływ na ρ_vac numeric |
| **op-MAG-resonance-formalization-2026-05-09 Phase 5** | CLOSED 2026-05-09 | TBD (audit Phase 1) | wpływ na m_e Mach |
| **op-uv3-phi0-renormalization** | CLOSED 2026-05-04 | używa P(g) NIE V(Phi) | likely OK |
| **particle_sector_closure** (P4) | CLOSED 2026-04-21 | używa A_tail NIE V(Phi) | likely OK |
| **op-MAG-Lorentz-A-mu-coupling** | CLOSED-FALSIFIED 2026-05-04 | TBD | impact na dygresje |
| **op-Phi-decomposition-photon** | CLOSED 2026-05-07 | TBD | niski impact |
| **sek08a master document** | LIVE | używa V_orig (deprecated) Z V_M9.1'' addendum | core, fundamental |
| **op-Phi-vacuum-scale-2026-05-09 Phase 1** | RECONNAISSANCE COMPLETE | używał V_orig (FLAGGED) | acknowledged, needs P11 audit |

## Plan szkicowy Phase 0-N

### Phase 0: Balance sheet
- Identify all V(Φ) usage in TGP framework
- Establish V_M9.1'' canonical reference (sek08a addendum 2026-05-02)
- 8/8 gate criteria
- NEEDS list (A1-A5)

### Phase 1: Audit reconnaissance
- Search agent: find all `V_orig`, `V_M9.1''`, `(β/3)Φ³`, `-γψ²(4-3ψ)²/12` patterns
- Per cycle classification table
- Sympy verification: V_orig vs V_M9.1'' specific quantities (ρ_vac, V_min, ψ_eq)
- Impact assessment per cycle

### Phase 2 (jeśli potrzebne): Targeted re-derivations
- Tylko cykle z critical impact (np. T-Λ rho_vac)
- NIE wszystkie cykle (out of scope)

### Phase 3 (final): Recommendations
- Lista cykli wymagających update
- Lista cykli OK as-is (V_orig był valid approximation)
- Updated cross-cycle map

## Probability assessment

| Outcome | Prob |
|---|---|
| Wszystkie cykle OK (V_orig valid jako approximation) | 30% |
| Większość OK, 1-2 wymaga update (np. T-Λ ρ_vac coefficient) | 50% |
| Większość wymaga update — framework crisis | 15% |
| Inconclusive — wymaga deeper review | 5% |

## Time budget

Audit cycle: ~1 session (similar do Phase 0-1 reconnaissance op-Phi-vacuum-scale).

## CALIBRATION_PROTOCOL compliance

Audit cycle, NIE derivation cycle. Anti-patterns 1-6 NIE applicable (brak
multi-candidate fit). Główne wymagania:
- Honest reporting: jeśli cykl używa V_orig, FLAG it explicit
- NIE forced "OK" jeśli rzeczywiście inconsistency
- Sympy verification dla każdej zidentyfikowanej rozbieżności

## Cross-references

- [[../op-Phi-vacuum-scale-2026-05-09/Phase1_6_strong_field_canonical_sympy.py]] — discovery
- [[../op-Phi-vacuum-scale-2026-05-09/NEEDS.md]] — P11 BLOCKER source
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] — V_orig deprecated linie 95-110
- [[../../meta/CALIBRATION_PROTOCOL.md]] — audit binding rules

## Status

**SCOPED. Phase 0-1 in progress.**
