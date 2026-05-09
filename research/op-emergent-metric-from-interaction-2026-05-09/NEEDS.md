---
title: "NEEDS — op-emergent-metric-from-interaction"
date: 2026-05-09
type: needs-list
status: 🟢 OPEN
parent: "[[./README.md]]"
tags:
  - needs
  - emergent-metric
---

# NEEDS — op-emergent-metric-from-interaction

## Six requirements (P1-P6 z README)

| # | Wymaganie | Sub-needs |
|---|-----------|-----------|
| **P1** | Formal definition g_eff = G[{Φ_i}] zgodne z S05 i §5.1 | N1, N2, N3 |
| **P2** | 1PN reproduction γ=β=1 z derivation, NIE postulat | N4, N5 |
| **P3** | 2.5PN β_ppE alternatywa do -15/4 z gradient cross-terms | N6, N7, N8 |
| **P4** | M9.2 Lenz back-reakcja sympy → m_inertial obliczalna | N9, N10 |
| **P5** | Cross-consistency z 3 SU(2) paths z SPIN cycle | N11 |
| **P6** | Falsifiability w |δφ̂_4| ≲ 0.18 1σ (GWTC-3) | N12 |

## Detailed needs

### N1: Wielociałowe rozszerzenie akcji TGP
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Akcja `S_TGP` z FOUNDATIONS §3 jest jednociałowa
(jedno pole Φ z source ρ). Cykl wymaga rozszerzenia na konfigurację
{Φ_i, ρ_i} wielu źródeł, gdzie każde Φ_i jest *lokalnym profilem*
**tego samego** pola Φ wokół źródła i. Formalizm musi:
- Zachować jednoznaczność (jedno Φ globalne, dekomponowane na Φ_i + tło Φ̄)
- Identyfikować cross-terms `∂_μΦ_i · ∂_νΦ_j` jako strukturalnie różne
  od jednoźródłowego `(∂Φ_total)²`
- Mieć linear superposition w słabopolowym limicie (sprawdzenie)

**Resolution path:**
- Adaptacja standard "many-body field theory" z osciloscylator/plasma physics
- Φ(x) = Φ̄(x) + Σ_i δΦ_i(x; x_i) + radiation
- Cross-terms emerge naturally z (∇Φ)² expansion

**Sympy:** verification superposition + cross-terms identification

### N2: Identyfikacja roli σ_ab jako gradient strain composite
**Priority:** HIGH
**Phase:** 1
**Status:** OPEN

**Description:** σ_ab istnieje od OP-7 T2 (2026-04-25) jako poziom 0 obiekt
(`K_ab − (1/3)δ_ab Tr K`). Cykl powinien aktywować σ_ab jako *embrion*
tensorowej struktury g_eff. Formalnie:
- Pokazać, że σ_ab transformuje się jak tensor pod współrzędnymi (level 0
  promotion to level 1 tensor)
- Sprawdzić, czy σ_ab w 2-source przypadku ma cross-terms naturalnie
- Identyfikować część g_eff zależną od σ_ab eksplicite

**Sympy:** transformation properties σ_ab pod boost/rotation

### N3: Demarkacja od BD/Horndeski (R1)
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Test bezsource'owy: w limicie {Φ_i} → ∅ (lub Φ_i → const
trywialnie), g_eff musi redukować się do flat-Minkowski (z Φ̄ background
at most z tłem cosmologicznym). To różni TGP od BD, gdzie g_μν ma
niezależną dynamikę vacuum.

**Resolution path:**
- Variational analysis: pokazać że NIE ma niezależnego równania pola dla g_eff
- δS/δg_eff jest *trywialnie zerowe* (g_eff = funkcjonał, nie zmienna)
- Dynamics tylko w Φ-EOM

**Sympy:** verification że g_eff variation gives no new EOM

### N4: 1PN expansion w schemacie wielociałowym
**Priority:** CRITICAL
**Phase:** 2
**Status:** OPEN

**Description:** Wyprowadzić γ_PPN, β_PPN dla 1 dominującego źródła w
nowym formalizmie. Cel: γ=β=1 jako konsekwencja struktury G[{Φ_i}],
nie postulat per-form.

**Resolution path:**
- Slabopolowe rozwinięcie δΦ_i wokół Φ̄
- Identyfikacja g_eff^00 → -(1+2U+...), g_eff^ij → δ_ij(1-2U+...)
- Match na PPN parametry standardowe

**Sympy:** PPN expansion automated

### N5: Solar system constraint check
**Priority:** HIGH
**Phase:** 2
**Status:** OPEN

**Description:** Mercury perihelion, Cassini Shapiro delay, lunar
laser ranging — wszystkie wymagają |γ-1|, |β-1| ≲ 10⁻⁵. Cykl musi
być z tym zgodny (już w 1PN level).

**Resolution path:** numerical check from N4 results.

### N6: 2-source case formalization
**Priority:** CRITICAL
**Phase:** 3
**Status:** OPEN

**Description:** Najprostszy nontrivial case: dwie masy m_1, m_2 w binary
inspiral. Wyprowadzić g_eff^μν w obecności obu Φ_i, identyfikować
gradient cross-terms.

**Resolution path:**
- Standardowy 2-body PN expansion (Blanchet review)
- TGP-natywne cross-terms identification
- Sympy verification

**Sympy:** pełen 2-body 2.5PN expansion

### N7: Effective phase modification
**Priority:** CRITICAL
**Phase:** 3
**Status:** OPEN

**Description:** Z 2-source g_eff wyprowadzić waveform phase modification
δφ(f) dla GW signal. Mapping na ppE framework (β_ppE).

**Resolution path:**
- Stationary phase approximation (SPA), z G_SPA correction (op-ppE-mapping
  Phase 1.5: G_SPA=48 sympy-exact dla structural-modification theories)
- δφ(f) = β_ppE · u^b, b=-1 dla 0.5PN

**Sympy:** SPA derivation + ppE mapping

### N8: β_ppE^new vs jednoźródłowe -15/4
**Priority:** CRITICAL
**Phase:** 3
**Status:** OPEN

**Description:** Pokazać, że β_ppE^new (z gradient cross-terms 2-source)
jest *strukturalnie różne* od jednoźródłowego -15/4. Konkretnie:
- Czy β_ppE^new zawiera *additional* terms, nieobecne w 1-source case?
- Czy te terms mają konkretną wartość liczbową?

**Resolution path:**
- Side-by-side comparison: 1-source ppE (= -15/4 z M9.1'') vs 2-source ppE
  (z cross-terms)
- Difference identification

**Sympy:** explicit comparison

### N9: Lenz back-reaction formalization
**Priority:** CRITICAL
**Phase:** 5
**Status:** OPEN

**Description:** §6 FOUNDATIONS deklaruje, że m_inertial = współczynnik
back-reakcji Φ-pola na zmianę konfiguracji równowagi. To NIGDY nie było
formalnie wyprowadzone z sympy.

**Resolution path:**
- Statyczne źródło: Φ_eq[ρ static] = równowaga
- Poruszające się źródło: Φ_eq[ρ(x-vt)] = przesunięty profil
- Przyspieszające źródło: Φ_eq[ρ(x-x(t))] niejednoznaczne, back-reakcja
- Linearyzacja: F_back = -m_inertial · a, gdzie m_inertial = ∫(...) [Φ_eq]

**Sympy:** linearizacja back-reakcji, identification m_inertial

### N10: Zasada równoważności automatycznie
**Priority:** HIGH
**Phase:** 5
**Status:** OPEN

**Description:** m_inertial (z N9) i m_grav (z 1PN N4) muszą być równe
*automatycznie* z S05 (to samo q i ρ). Cykl musi to pokazać explicite.

**Resolution path:** porównanie expression z N9 i N4.

**Sympy:** equality verification

### N11: Cross-consistency z SU(2) paths
**Priority:** CRITICAL (programmatic)
**Phase:** 6
**Status:** OPEN

**Description:** Mechanizm interakcji generujący g_eff (gradient cross-terms)
musi być zgodny z mechanizmem generującym SU(2) z SPIN cycle (3 paths).
Konkretnie:
- Path A (N18, dynamic-equilibrium bifurkacja) — czy wymagany dla g_eff też?
- Path B (N21, horizon multipole z M9.1'') — w niniejszym cyklu g_eff nie ma
  ψ=4/3 horyzontu strict, więc B może wymagać re-examinacji
- Path C (N19, external embedding) — naturalnie zgodne z "tensor z interakcji"

**Resolution path:** structural analysis, nie sympy.

### N12: GWTC-3 falsifier check
**Priority:** CRITICAL (hard test)
**Phase:** 4
**Status:** OPEN

**Description:** β_ppE^new (z N8) skonwertowane na δφ̂_4 (TIGER framework).
Required: |δφ̂_4| ≲ 0.18 (1σ).

**Resolution path:**
- Conversion β_ppE → δφ̂_4 (op-GWTC3-reanalysis Phase 1 mapping)
- Numerical check vs GWTC-3 posterior

**Numerical:** posterior comparison.

### N13: c_GW = c check
**Priority:** HIGH (hard test)
**Phase:** 4
**Status:** OPEN

**Description:** GW170817 wymaga c_GW = c do 10⁻¹⁵. Emergentna metryka
NIE może wprowadzić scalar mode propagującego z c' ≠ c.

**Resolution path:**
- Identyfikacja propagating modes w nowym g_eff
- Sprawdzenie speed (z dispersion relation linearyzacji)

**Sympy:** dispersion relation derivation

### N14: LIGO scalar mode bound check
**Priority:** HIGH
**Phase:** 4
**Status:** OPEN

**Description:** LIGO scalar polarization < few %. Emergentna metryka
może mieć additional polarizations — sprawdzić, czy są w bound.

**Resolution path:**
- Polarization decomposition w nowym g_eff
- Amplitude estimate scalar component vs tensor

**Numerical:** check vs LIGO bounds.

## Status summary

| Need | Phase | Priority | Status |
|------|-------|----------|--------|
| N1 (action many-body) | 1 | CRITICAL | ✅ RESOLVED (Phase 1, 16/16) |
| N2 (σ_ab activation) | 1 | HIGH | ✅ RESOLVED (Phase 1, 16/16) |
| N3 (BD demarcation) | 1 | CRITICAL | ✅ RESOLVED (Phase 1, 16/16) |
| N4 (1PN derivation) | 2 | CRITICAL | ✅ RESOLVED (Phase 2, 7/7) |
| N5 (solar system) | 2 | HIGH | ✅ RESOLVED (Phase 2, 7/7) |
| N6 (2-source formal) | 3 | CRITICAL | ✅ RESOLVED (Phase 3, 5/5) |
| N7 (phase modification) | 3 | CRITICAL | ✅ RESOLVED (Phase 3, 5/5) |
| N8 (β_ppE^new) | 3 | CRITICAL | ✅ RESOLVED (Phase 3, 5/5) |
| N9 (Lenz formal) | 5 | CRITICAL | ✅ RESOLVED (Phase 5, 10/10) |
| N10 (equiv principle) | 5 | HIGH | ✅ RESOLVED (Phase 5, 10/10) |
| N11 (cross-SU(2)) | 6 | CRITICAL | ✅ RESOLVED (Phase 6, 11/11) |
| N12 (GWTC-3 check) | 4 | CRITICAL | ✅ RESOLVED (Phase 4, 8/8) |
| N13 (c_GW=c) | 4 | HIGH | ✅ RESOLVED (Phase 4, 8/8) |
| N14 (LIGO scalar) | 4 | HIGH | DEFERRED (R5 risk, multi-session) |

**14 needs total. RESOLVED: 13. DEFERRED: 1 (N14 R5 risk, multi-session).**

**Cumulative sympy: 57/57 PASS (100%) across Phase 1-6.**

## Dependencies

```
N1 → N2 → N3 (Phase 1, sequential)
N3 → N4 → N5 (Phase 2)
N4 → N6 → N7 → N8 (Phase 3)
N8 → N12, N13, N14 (Phase 4 falsifier)
N1, N4 → N9 → N10 (Phase 5 momentum)
all → N11 (Phase 6 cross-consistency)
```

Phase 4 jest hard gate. Failure tam → STRUCTURAL_NO_GO (cycle terminates honestly).
