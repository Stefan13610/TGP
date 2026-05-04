---
title: "L02 — β/γ kolizja semantyczna: WF FP vs faza złamana"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L02_beta_gamma_semantics
tags:
  - audit
  - ontology
  - beta-gamma
  - WF-fixed-point
  - broken-symmetry
  - semantic-collision
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
  - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["beta-gamma-renotation"]
  depends_on: []
  impacts: []
  source_of_status:
    - "[[../../TGP_FOUNDATIONS.md]] §7"
    - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LK-1d"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L02 — β/γ kolizja semantyczna

## Klasa: LUKA ONTOLOGICZNA (semantyczna) • Priorytet: P3

## Diagnoza

Symbol `β/γ` jest używany w TGP w **dwóch różnych znaczeniach**, które
dają sprzeczne wyniki — i nigdzie w rdzeniu nie ma explicit
różnicowania.

### Znaczenie 1: stosunek wykładników krytycznych w punkcie stałym Wilsona–Fishera

[[../../TGP_FOUNDATIONS.md]] §7 — cykle M3–M8 (kwiecień 2026):

> W jednoskładnikowym skalarnym Z₂ na poziomie LPA β/γ w punkcie stałym
> WF jest **ujemne (≈ −0.3 do −0.5)** w obu schematach RG (MK-RG, NPRG/Wetterich).

To są **wykładniki krytyczne** w sensie statystyczno-fizycznym:
- β = wykładnik order-parameter (`m ~ |T-T_c|^β`)
- γ = wykładnik susceptibility (`χ ~ |T-T_c|^{-γ}`)

W 3D Ising universality class: β ≈ 0.326, γ ≈ 1.237 ⇒ β/γ ≈ 0.264 (nie −0.3
— LPA daje rough estimate, dokładne wartości z bootstrap).

Audit w FOUNDATIONS broni się: *„TGP nie żyje na T_c, żyje w fazie złamanej"*.

### Znaczenie 2: stosunek współczynników GL w fazie złamanej (warunek próżni)

[[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LK-1d:

> β/γ → 1 zweryfikowane numerycznie — **6/6 PASS**, R = 0.88 ± 0.38
> (1 w 0.3σ), skewness ≈ 0.005

To są **współczynniki potencjału GL** w sek08a:

```
V(φ) = (β/3)·φ³ - (γ/4)·φ⁴
```

Warunek próżni `V'(1) = 0` daje `β = γ` ⇒ `β/γ = 1`. To jest *algebraiczna
tożsamość*, nie zmierzona własność — wymóg vacuum stability.

## Sprzeczność semantyczna

W jednym artykule TGP czytelnik może natknąć się na:

> *„β/γ → 1 (warunek próżni, LK-1d 6/6 PASS)"*

i

> *„β/γ ≈ −0.3 (LPA WF FP, M3–M8 sfalsyfikowane)"*

bez explicit różnicowania, że to **dwie różne wielkości pod tą samą
nazwą**. Każdy zewnętrzny audyt (referee, reviewer) zatrzyma się na tym
i zażąda wyjaśnienia.

## Konsekwencja

To nie jest fizyczna sprzeczność — to **notation collision**. Ale
rozprzestrzeniona przez 25+ plików `core/sek*` i `research/op*`, plus
dziesiątki skryptów. Bez renotacji manuskrypt jest niespójny.

## Obecna obrona w FOUNDATIONS §7

> **Wnioski reinterpretacyjne (pod nową ramą "GR w limicie"):**
> - M3–M8 odpowiadało na pytanie o **uniwersalność klasy krytycznej**
>   (wykładniki w punkcie stałym WF), nie na pytanie o **fenomenologię
>   TGP** (klasyczna dynamika pola średniego w fazie złamanej, na
>   skończonej skali korelacji).
> - TGP fizycznie nie żyje na T = T_c (faza krytyczna), tylko w **fazie
>   złamanej** (`m₀² < 0`, `v² = |m₀²|/λ₀`), z **klasyczną dynamiką
>   Φ-EOM** plus poprawki z fluktuacji.
> - Wartość `β/γ` w punkcie stałym WF **nie jest** liczbą wiążącą dla
>   fenomenologii TGP. **Wiążące jest:** liczbowa zgodność rozwiązań
>   Φ-EOM z fenomenologią GR w odpowiednim regimie.

To poprawnie tłumaczy *epistemicznie*, że są to dwie różne wielkości,
ale nie wprowadza explicit notational distinction w plikach LaTeX.

## Status w audycie

Nie wymieniony explicit jako CRITICAL/HIGH/MEDIUM. To jest semantyczna
luka, którą sam audit nie wychwytuje, bo skupia się na fizycznych
sprzecznościach.

## Rekomendacja: renotacja

Krótki editorial cycle (1 tydzień, mechaniczny):

1. **Wprowadzić jawne oznaczenia:**
   - `(β/γ)_WF` = stosunek wykładników krytycznych w punkcie stałym WF
   - `(β/γ)_GL` lub `(β/γ)_BS` = stosunek współczynników GL w fazie złamanej
   (BS = broken-symmetry)

2. **Globalny grep + replace** w `core/`, `research/`, `scripts/`:
   - Wszystkie miejsca cytujące M3–M8 → `(β/γ)_WF`
   - Wszystkie miejsca cytujące LK-1d, vacuum condition `β=γ` → `(β/γ)_GL`

3. **Dodać sekcję notation w `dodatekA_notacja.tex`** explicit
   różnicującą oba.

4. **Update [[../../TGP_FOUNDATIONS.md]] §7** z explicit notation.

## Pliki dotknięte (estymata)

- `core/sek08*.tex` — ~10 lokalizacji
- `core/sek10_N0_wyprowadzenie/` — ~5 lokalizacji
- `axioms/notacja/dodatekA_notacja.tex` — dodać sekcję
- `research/op1-op2-op4/M{3..8}_*.md` — ~20 lokalizacji
- `research/op-newton-momentum/M9*.md` — ~10 lokalizacji
- `scripts/lk1d_beta_gamma_ratio.c` — komentarze + naming
- `TGP_FOUNDATIONS.md` §7 — explicit renotacja

**Estymata:** 1 tydzień (głównie mechaniczne grep+replace + sanity check).

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L2
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../TGP_FOUNDATIONS.md]] §7
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LK-1d
