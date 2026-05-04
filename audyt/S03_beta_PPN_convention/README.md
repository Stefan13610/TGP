---
title: "S03 — β_PPN: 1/2 (metric) vs 1 (master formula)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S03_beta_PPN_convention
tags:
  - audit
  - PPN
  - beta
  - critical
  - convention-tension
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../S01_metric_four_forms]]"
  - "[[../S02_volume_element_M9]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L2
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["PPN-convention-lock", "B8-Lagrangean-derivation"]
  depends_on:
    - "[[../S01_metric_four_forms]]"
  impacts: []
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.3, §M.2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# S03 — β_PPN: czysto metryczne 1/2 vs master formula 1

## Klasa: STRUKTURALNA SPRZECZNOŚĆ KONWENCJI • Priorytet: P1

## Diagnoza

Bezpośrednia ekspansja metryki M9.1'' wokół ψ=1+ε:

```
f(ψ) = (4-3ψ)/ψ
f(1+ε) = 1 - 4ε + 4ε² - 4ε³ + …
```

Newton matching `-4ε = -2U` daje `ε = U/2`. Wstawienie:

```
g_tt/(-c²) = 1 - 2U + U² - U³/2 + …
```

PPN wymaga `g_tt = -(1 - 2U + 2β·U²)`, więc:

```
β_PPN_metric = 1/2
```

Wartość `β=1` (zgodna z testem Mercury 1±10⁻⁴) osiąga się **tylko** przez
„master formula" z kinetic correction:

```
β_PPN = f''(1)/f'(1)² + 2c₂/f'(1) = 1/2 + 1/2 = 1
```

gdzie `c₂ = -1` pochodzi z **dokładnie tej samej** nieliniowości α=2,
która w formie potęgowej (I) dawała `β=4` (a została stamtąd
sfalsyfikowana jako wykluczająca). To jest **algebraiczne dopasowanie
metryki+kinetyki**, nie niezależne potwierdzenie.

Dodatkowo: linearyzacja w `sek08c.tex` lin. 183–191 (czysto metryczna,
forma (I) → β=2) wciąż jest w pliku — **logicznie sprzeczna z M9.1''
β=1**.

## Czego brakuje

Niezależne potwierdzenie c₂ = -1 musiałoby pochodzić z:

1. **Lagrangianu**, z którego forma f(ψ) = (4-3ψ)/ψ wynika wariacyjnie
   (B8 audytu). Aktualnie M9.1'' jest *postulatem* na podstawie
   „potrójnej motywacji" (E1/E2/E3), które audit § B.8 ocenia jako
   *cyrkularne* — sprowadzają się do jednego ad-hoc postulatu f(4/3)=0.
2. **Niezależnego eksperymentu** mierzącego c₂ inaczej niż przez β_PPN
   (np. z dispersji GW lub wyższych rzędów PPN).

Audit § A.3 sam acknowledges: master formula jest *„algebraic dopasowanie
metryki+kinetyki"*, nie niezależne potwierdzenie.

## Status w `meta/AUDYT_TGP_2026-05-01.md`

§A.3 — **CRITICAL**, oznaczony „CLOSED konwencyjnie" (§ M.6) z wyborem
master formula jako obowiązującej (§ M.2). Audit wprost mówi:

> czysto metryczna konwencja daje β=1/2 dla M9.1'' (tension z Mercury),
> master formula daje β=1 ✓; wybór master jest *algebraicznym
> dopasowaniem* (acknowledged), nie niezależną derivation.

To znaczy: status `β_PPN=1` w README/FOUNDATIONS jest **zależny od
wyboru konwencji**. W zewnętrznym audycie (referee, reviewer artykułu)
ta arbitralność będzie pierwszą rzeczą do podważenia.

## Powiązany problem: linearyzacja w sek08c lin. 183–191

`sek08c.tex` lin. 183–191 stosuje czysto metryczną PPN dla formy (I)
`g_tt = -c²/(1+2U)` → **β=2**, nie β=1. To jest sprzeczne z M9.1''
β=1 (master formula). Audit § A.1 i § B.10 oba wskazują na ten konflikt.

## Rekomendacja

Trzy ścieżki, w kolejności preferencji:

### Ścieżka A: Lagrangean derivation (B8 cycle)

Otworzyć dedykowany cykl `op-M911-lagrangean-derivation/` z celem
**wyprowadzenia f(ψ) = (4-3ψ)/ψ z akcji `S[Φ, g]` przez wariację**.
Jeśli się uda, c₂ = -1 będzie konsekwencją struktury Lagrangianu, a nie
post-hoc dopasowaniem.

**Estymata:** 6–10 tygodni (otwarty problem fizyczny).

### Ścieżka B: explicit acknowledgment w manuskrypcie

Dodać do README, FOUNDATIONS, i sek08c jawne stwierdzenie:

> β_PPN = 1 dla M9.1'' uzyskane przez master formula z kinetic
> correction `c₂ = -1` z α=2 nieliniowości — to algebraiczne
> dopasowanie, nie niezależne potwierdzenie. Niezależny test wymaga
> wyprowadzenia formy f(ψ) z Lagrangianu (cykl B8 pending).

To utrzymuje uczciwość epistemiczną, ale nie rozwiązuje sprzeczności
strukturalnej.

### Ścieżka C: alternatywna konwencja PPN

Przyjąć czysto metryczną konwencję → β=1/2 → 0.5σ napięcie z Mercury
(pomiar β=1.000 ± 0.001) jest **5·10²σ** odchyleniem. To by wymagało
falsyfikacji M9.1''.

→ Ścieżka C **wycofuje M9.1''** jako kanoniczną. Niewskazana, ale
formalnie spójniejsza niż obecna sytuacja.

**Decyzja:** Ścieżka A (preferowana) lub Ścieżka B (interim).

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `sek08c.tex` lin. 183–191 | przepisać dla M9.1'' lub usunąć |
| `sek08c.tex` (nagłówek) | dodać explicit konwencja master formula |
| `README.md` | dodać caveat „β=1 z master formula, B8-pending" |
| `TGP_FOUNDATIONS.md` lin. 64–69 | dodać caveat |
| nowy: `research/op-M911-lagrangean-derivation/` | nowy cykl |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §I.S3
- [[../PRIORITY_MATRIX.md]] klaster A
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.3, §B.8 (cyrkularność M9.1''
  „potrójnej motywacji"), §B.10 (β=2 forma minimalna), §M.2
- [[../S01_metric_four_forms]] (zależność wstępna)
- [[../S02_volume_element_M9]] (zależność równoległa)
