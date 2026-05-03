<!--
FINDINGS.template.md
Co JEST GOTOWE w tym folderze. Eksportowalne wyniki — numeryczne, formalne,
strukturalne — które inny folder może skonsumować.

Reguły:
  1. Każdy finding musi cytować WERYFIKOWALNE źródło (`source:` =
     plik wewnątrz folderu z PASS-em / wartością / formułą / sekcją).
  2. NIE wpisujemy "FULL CONVERGENCE" / "DERIVED" bez liczbowego PASS-a
     z konkretnego pliku. Reguła "anty-overclaim" twarda.
  3. Format wartości: jeśli to liczba — z jednostką; jeśli to formuła —
     LaTeX inline; jeśli to status strukturalny — werdykt + plik.
  4. Każdy finding ma `consumers:` — listę folderów, które MOGĄ go użyć.
     Pusta lista jest OK — broadcast w RESEARCH_BUS.md jest osobno.
  5. Synchronizacja: po dopisaniu nowego finding'a — wpis w
     [[meta/research/RESEARCH_BUS.md]] (broadcast).
-->

---
title: "FINDINGS — <nazwa folderu>"
date: <YYYY-MM-DD>
parent: "[[README.md]]"
type: findings
tgp_owner: research/<nazwa-folderu>
tags:
  - findings
  - <topic-tags>
---

# FINDINGS — <nazwa folderu>

> Eksportowalne wyniki tego folderu. Każdy item musi mieć cytowane źródło
> (`source:`). Wpisy bez source są zabronione (anty-overclaim).

## Wyniki numeryczne

| ID | Wynik | Wartość / formuła | Source | Consumers |
|----|-------|-------------------|--------|-----------|
| F1 | <co policzono> | <wartość ± błąd, jednostka> | <`<plik>:<linia>` lub `<plik>` z PASS X/Y> | [`research/<konsument>`] |
| F2 | … | … | … | … |

## Wyniki strukturalne / analityczne

| ID | Statement | Forma | Source | Consumers |
|----|-----------|-------|--------|-----------|
| F3 | <lemma / zamknięcie strukturalne> | <LaTeX / opis> | <plik> | [`research/<konsument>`] |

## Status PASS-ów per faza (jeśli folder jest 3-fazowy)

| Faza | Tests | Pass-rate | Werdykt | Plik wynikowy |
|------|-------|-----------|---------|---------------|
| Phase 1 | <X/Y> | <%> | <CLOSED / OPEN / STRUCTURAL HINT> | `Phase1_results.md` |
| Phase 2 | <X/Y> | <%> | … | `Phase2_results.md` |
| Phase 3 | <X/Y> | <%> | … | `Phase3_results.md` |

## Falsyfikatory ustawione przez ten folder

| ID | Predykcja | Kryterium falsyfikacji | Eksperyment / timeline |
|----|-----------|-------------------------|------------------------|
| FX1 | <co przewiduje TGP> | <co go obali> | <eksperyment + horyzont> |

## Pre-existing flag

<!--
  Jeśli `pre_existing_findings: true` w README — to znaczy, że ten plik
  ISTNIAŁ przed Sesją 5. Agent NIE nadpisał. Zachowujemy oryginalną treść,
  agent dopisuje tylko strukturalny header (frontmatter + sekcję
  "Pre-existing flag"). Pełna treść poniżej tej notki.
-->

## Cross-references

- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych wyników
- [[meta/research/FOLDER_STATUS_INDEX.md]] — pozycja folderu w globalnej mapie
