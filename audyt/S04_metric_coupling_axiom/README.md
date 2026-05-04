---
title: "S04 — ax:metric-coupling vs L_mat = -(q/Φ_0)·φ·ρ"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S04_metric_coupling_axiom
tags:
  - audit
  - axiom
  - metric-coupling
  - L_mat
  - critical
  - architectural
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../L01_rho_operational]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["A4-formal-derivation", "fifth-force-test"]
  depends_on: []
  impacts:
    - "[[../L01_rho_operational]]"
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.4, §N"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# S04 — `ax:metric-coupling` ⊥ `L_mat = -(q/Φ_0)·φ·ρ`

## Klasa: STRUKTURALNA SPRZECZNOŚĆ AKSJOMATU • Priorytet: P1

## Diagnoza

W `core/` współistnieją:

### Aksjomat metric-coupling

`sek08_formalizm.tex` lin. 11115–11132 (`ax:metric-coupling`):

> Materia sprzęga z Φ **wyłącznie przez metrykę efektywną**, nie ma
> bezpośredniego dostępu do Φ.

### Lagrangean materii

`sek08a_akcja_zunifikowana.tex` lin. 74–83 (`eq:L-mat-unified`):

```
L_mat = -(q/Φ_0)·φ·ρ
```

To zawiera **bezpośredni czynnik φ poza metryką** — sprzęganie typu
dilatonowego (Brans–Dicke).

### Co to znaczy fizycznie

Po pełnej akcji z volume element (zakładając formę I):

```
S_mat = ∫d⁴x √(-g_eff)·L_mat = -(q·c_0/Φ_0)·∫d⁴x·φ²·ρ
```

Pojawia się jawny czynnik `φ²·ρ`. **Aksjomat dawałby brak piątej siły**
(equivalence principle automatycznie spełniony przez czyste sprzęganie
metryczne), ale `φ·ρ` w `L_mat` generuje **siłę typu Brans–Dicke**:

```
∇_μ(φρ) → dφ/dr · ρ = piąta siła ∝ ∇φ
```

Brans–Dicke parameter `ω_BD` jest mocno ograniczone:
- Cassini (Shapiro delay): ω_BD > 40 000
- Eöt-Wash, MICROSCOPE: equivalence principle do 10⁻¹³–10⁻¹⁵

**Brakuje weryfikacji**, że TGP nie generuje siły piątego typu w
zakresie obserwowalnym.

## „Rozwiązanie" audytu (Option-2)

§ N decyzja audytowa:

> Aksjomat `ax:metric-coupling` zachowuje pierwszeństwo. Czynnik `φ` w
> `eq:L-mat-unified` interpretowany jest jako **derived consequence**
> struktury volume element `√(-g_eff) ∝ φ` plus ontologicznej roli `Φ`
> jako konstytuującego `ρ` (`ax:zrodlo`); `ρ` traktowane jest jako
> `T^μ_μ/c_0²` pochodne ze stress-energy tensor materii w kanonicznym
> sprzęganiu metrycznym, NIE jako independent dilaton coupling
> Brans–Dicke type.

### Co to faktycznie znaczy

To jest **re-interpretacja, nie derywacja**. Konkretnie brakuje:

1. **Formalnej derywacji** `φ·ρ = T^μ_μ/c_0² · √(-g_eff)/c_0` z
   kanonicznego matter Lagrangianu `L_mat[ψ_m, g_eff]`.
2. **Sprawdzenia**, że dla materii nierelatywistycznej (`T^μ_μ ≠ ρ_rest`
   jednoznacznie) ta interpretacja jest spójna.
3. **Eksperymentalnego testu Eöt-Wash/MICROSCOPE**, że `dφ/dr · ρ`
   nie generuje obserwowalnego naruszenia equivalence principle dla
   ciał o różnym składzie chemicznym.
4. **Porównania z bound** Brans–Dicke ω_BD > 40 000.

Audit § N.4 sam acknowledges:

| element | status |
|---|---|
| `eq:L-mat-unified` formuła | ZACHOWANA (φ·ρ retained) |
| **Interpretacja** φ·ρ coupling | REFORMA: derived from √(-g_eff) ∝ φ + ax:zrodlo |
| Brans–Dicke fifth force | WYKLUCZONE (explicit Option-2 decyzja) |

„Wykluczone explicit decyzją" to **postulat**, nie wyprowadzenie.

## Wpływ na predykcje

- **MICROSCOPE 10⁻¹⁵ WEP** (M9.2): TGP twierdzi pass, ale tylko jeśli
  φ·ρ jest faktycznie equivalent do czysto metrycznego sprzęgania.
  Audit § B.9 wymienia jako pending: „Test kompozycji (dwa źródła o
  różnym q lub niejednorodnym ρ)".
- **Cassini Shapiro** (γ=1.0000, |γ-1| < 2.3·10⁻⁵): TGP twierdzi pass,
  ale `dφ/dr · ρ` może modyfikować dyspersję sygnału w sposób zależny
  od mocy Slońca (ω_BD analog).
- **Test piątej siły z LLR** (Lunar Laser Ranging): η < 4.4·10⁻⁴, TGP
  twierdzi η=0 ale przez czyste sprzęganie metryczne, nie z Option-2
  re-interpretation.

## Status w `meta/AUDYT_TGP_2026-05-01.md`

§A.4 — **CRITICAL**, oznaczony „CLOSED architecturally" (§ N.5)
poprzez Option-2 decyzję. Faktyczny status: **CLOSED-decision-only**.
Formalna derywacja, niezależny test Eöt-Wash, i porównanie z ω_BD nie
zostały wykonane.

## Rekomendacja

Otworzyć dedykowany cykl `op-A4-formal-derivation/` z trzema fazami:

### Phase 1 — formal kowariantna derywacja

Pokazać, że dla czysto metrycznego sprzęgania `L_mat[ψ_m, g_eff]` z
g_eff = (Φ/Φ_0)·δ:

```
∫d⁴x √(-g_eff) L_mat[ψ_m, g_eff] → ?·∫d⁴x φ²·ρ
```

Wymaga: explicit `T^μ_ν[ψ_m]`, kontrakcja z `g_eff^μν`, śledzenie
rozkładu na `ρ` w slot `L_mat` w sek08a.

### Phase 2 — test piątej siły TGP

Numeryczny test: dla solitona o danym `q`, gradient `dφ/dr` w
otoczeniu pola materii, siła piątego typu vs MICROSCOPE bound.

### Phase 3 — porównanie z Brans–Dicke

Mapowanie TGP `q/Φ_0` na `ω_BD`. Sprawdzenie, czy `q/Φ_0 ≪ 1/ω_BD`
(40 000)⁻¹ ≈ 2.5·10⁻⁵.

**Estymata:** 4–6 tygodni.

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `sek08a` lin. 74–83 | dodać formal derivation φ·ρ z T^μ_μ |
| `sek08_formalizm.tex` lin. 11115–11132 | sync z derived form |
| `TGP_FOUNDATIONS.md` §4 (warstwy 3a/3b) | dodać explicit Option-2 |
| nowy: `research/op-A4-formal-derivation/` | nowy cykl |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §I.S4
- [[../PRIORITY_MATRIX.md]] klaster B
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.4, §N
- [[../L01_rho_operational]] (wspólna ontologia ρ)
- [[../S05_tensor_sector_singleField]] (równoległy aksjomatyczny problem)
